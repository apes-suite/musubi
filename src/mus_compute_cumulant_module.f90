! Copyright (c) 2021-2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF SIEGEN “AS IS” AND ANY EXPRESS
! OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
! OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL UNIVERSITY OF SIEGEN OR CONTRIBUTORS BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!> This module contains functions for Cumulant and Cascaded for D3Q27 stencil.
!! author: Gregorio Gerardo Spinelli
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011-2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013-2015, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF SIEGEN “AS IS” AND ANY EXPRESS
! OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
! OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL UNIVERSITY OF SIEGEN OR CONTRIBUTORS BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
module mus_compute_cumulant_module
  use iso_c_binding,            only: c_f_pointer

  use env_module,               only: rk
  use tem_varSys_module,        only: tem_varSys_type
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit
  use tem_debug_module,         only: dbgUnit
  use tem_float_module,         only: operator(.fne.)
  use tem_timer_module,         only: tem_starttimer, tem_stoptimer
  use tem_param_module,         only: div1_3, div2_3, div1_9, div1_2,     &
    &                                 div1_27, div1_6

  use mus_directions_module,    only: qN00, q0N0, q00N, q100, q010, q001, &
    &                                 q0NN, q0N1, q01N, q011, qN0N, q10N, &
    &                                 qN01, q101, qNN0, qN10, q1N0, q110, &
    &                                 qNNN, qNN1, qN1N, qN11, q1NN, q1N1, &
    &                                 q11N, q111
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_varSys_module,        only: mus_varSys_data_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_param_module,         only: mus_param_type
  use mus_derVarPos_module,     only: mus_derVarPos_type
  use mus_gradData_module,      only: mus_gradData_type
  use mus_scheme_type_module,   only: mus_scheme_type
  !use mus_timer_module,         only: mus_timerHandles

  implicit none
  private

  public :: weights_from_layout
  public :: weights_ijg
  public :: weights_ibg
  public :: weights_abg
  public :: central_moment
  public :: central_moment_split
  public :: cm_to_pdf
  public :: cumulant_d3q27
  public :: kum_ij_g
  public :: kum_i_bg
  public :: kum_abg
  public :: kpc_i_bg
  public :: kpc_ij_g
  public :: kpc_ijk
  public :: cumulant_d3q27_extended_generic
  public :: cumulant_d3q27_extended_fast
  public :: cascaded_d3q27

  integer, parameter :: QQ = 27  !< number of pdf directions
  integer, parameter :: q000 = 27


contains

! **************************************************************************** !
  !> allocate weights from D3Q27 ordered disposition to cumulant disposition
  pure function weights_from_layout( layout_weight ) result ( w )
    !> layout weights
    real(kind=rk), intent(in) :: layout_weight(27)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: w(-1:1, -1:1, -1:1)
    ! --------------------------------------------------------------------------

    !allocate weights from D3Q27 ordered disposition to cumulant disposition
    w(-1, 0, 0) = layout_weight( qN00 )
    w( 0,-1, 0) = layout_weight( q0N0 )
    w( 0, 0,-1) = layout_weight( q00N )
    w( 1, 0, 0) = layout_weight( q100 )
    w( 0, 1, 0) = layout_weight( q010 )
    w( 0, 0, 1) = layout_weight( q001 )
    w( 0,-1,-1) = layout_weight( q0NN )
    w( 0,-1, 1) = layout_weight( q0N1 )
    w( 0, 1,-1) = layout_weight( q01N )
    w( 0, 1, 1) = layout_weight( q011 )
    w(-1, 0,-1) = layout_weight( qN0N )
    w( 1, 0,-1) = layout_weight( q10N )
    w(-1, 0, 1) = layout_weight( qN01 )
    w( 1, 0, 1) = layout_weight( q101 )
    w(-1,-1, 0) = layout_weight( qNN0 )
    w(-1, 1, 0) = layout_weight( qN10 )
    w( 1,-1, 0) = layout_weight( q1N0 )
    w( 1, 1, 0) = layout_weight( q110 )
    w(-1,-1,-1) = layout_weight( qNNN )
    w(-1,-1, 1) = layout_weight( qNN1 )
    w(-1, 1,-1) = layout_weight( qN1N )
    w(-1, 1, 1) = layout_weight( qN11 )
    w( 1,-1,-1) = layout_weight( q1NN )
    w( 1,-1, 1) = layout_weight( q1N1 )
    w( 1, 1,-1) = layout_weight( q11N )
    w( 1, 1, 1) = layout_weight( q111 )
    w( 0, 0, 0) = layout_weight( q000 )

  end function
! **************************************************************************** !

! **************************************************************************** !
  !> Calculating central moment weights
  !! This follows equation 3 in cumulent paper (Geier .et al 2017)
  pure function weights_ijg( ii, jj, ig, w ) result ( weight )
    !> indeces of the pdf
    integer, intent(in) :: ii, jj
    !> order gamma of moments
    integer, intent(in) :: ig
    !> cumulant ordered weights
    real(kind=rk), intent(in) :: w(-1:1, -1:1, -1:1)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: weight
    integer :: kk
    real(kind=rk), parameter :: ii_rk(-1:1) = [ -1._rk, 0._rk, 1._rk ]
    ! --------------------------------------------------------------------------
    weight = 0.0_rk

    do kk = -1, 1
      weight = weight + w(ii, jj, kk) * ( ii_rk(kk) ** ig )
    end do

  end function
! **************************************************************************** !

! **************************************************************************** !
  !> Calculating central moment weights
  !! This follows equation 3 in cumulent paper (Geier .et al 2017)
  pure function weights_ibg( ii, ib, ig, w_ij_g ) result ( weight )
    !> partial weights_ij_g
    real(kind=rk), intent(in) :: w_ij_g(-1:1, -1:1, 0:2)
    !> indeces of the pdf
    integer, intent(in) :: ii
    !> order gamma of moments
    integer, intent(in) :: ib, ig
    ! --------------------------------------------------------------------------
    real(kind=rk) :: weight
    integer :: jj
    real(kind=rk), parameter :: ii_rk(-1:1) = [ -1._rk, 0._rk, 1._rk ]
    ! --------------------------------------------------------------------------
    weight = 0.0_rk

    do jj = -1, 1
      weight = weight + w_ij_g( ii, jj, ig ) * ( ii_rk(jj) ** ib )
    end do

  end function
! **************************************************************************** !

! **************************************************************************** !
  !> Calculating central moment weights
  !! This follows equation 5 in cumulent paper (Geier .et al 2017)
  pure function weights_abg( ia, ib, ig, w_i_bg ) result ( weight )
    !> partial weights_i_bg
    real(kind=rk), intent(in) :: w_i_bg(-1:1, 0:2, 0:2)
    !> order gamma of moments
    integer, intent(in) :: ia, ib, ig
    ! --------------------------------------------------------------------------
    real(kind=rk) :: weight
    integer :: ii
    real(kind=rk), parameter :: ii_rk(-1:1) = [ -1._rk, 0._rk, 1._rk ]
    ! --------------------------------------------------------------------------
    weight = 0.0_rk

    do ii = -1, 1
      weight = weight + w_i_bg( ii, ib, ig ) * ( ii_rk(ii) ** ia )
    end do

  end function
! **************************************************************************** !

! ****************************************************************************** !
  !> Calculating central moment by spliting among directions.
  !! This follows equations 43, 44, 45 in cumulent paper (Geier .et al 2015)
  !! We first do x direction for better performance.
  pure function central_moment_split( f, a, b, g, ux, uy, uz ) result ( kappa )
    !> PDF
    real(kind=rk), intent(in) :: f(-1:1, -1:1, -1:1)
    !> order of central moments
    integer, intent(in) :: a, b, g
    real(kind=rk), intent(in) :: ux, uy, uz
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: ka(-1:1, -1:1), kb(-1:1)
    real(kind=rk) :: kappa
    integer :: ii, jj, kk
    real(kind=rk), parameter :: ii_rk(-1:1) = [ -1._rk, 0._rk, 1._rk ]
    ! ---------------------------------------------------------------------------

    ka = 0.0_rk
    kb = 0.0_rk
    kappa = 0.0_rk

    do kk = -1, 1
      do jj = -1, 1
        do ii = -1, 1
          ka(jj,kk) = ka(jj,kk) + f( ii, jj, kk ) * ( ( ii_rk(ii) - ux ) ** a )
        end do
        kb(kk) = kb(kk) + ka(jj,kk) * ( ( ii_rk(jj) - uy ) ** b )
      end do
      kappa = kappa + kb(kk) * ( ( ii_rk(kk) - uz ) ** g )
    end do

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !> Calculating central moment.
  !! This follows equations 21 in cumulent paper ( Geier .et al 2015 )
  pure function central_moment( f, a, b, g, ux, uy, uz ) result ( kappa )
    !> PDF
    real(kind=rk), intent(in) :: f(-1:1, -1:1, -1:1)
    !> order of central moments
    integer, intent(in) :: a, b, g
    real(kind=rk), intent(in) :: ux, uy, uz
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: kappa
    integer :: ii, jj, kk
    real(kind=rk), parameter :: ii_rk(-1:1) = [ -1._rk, 0._rk, 1._rk ]
    ! ---------------------------------------------------------------------------

    kappa = 0.0_rk

    do kk = -1, 1
      do jj = -1, 1
        do ii = -1, 1
          kappa = kappa + f( ii, jj, kk ) * ( ( ii_rk(ii) - ux ) ** a )  &
            &                             * ( ( ii_rk(jj) - uy ) ** b )  &
            &                             * ( ( ii_rk(kk) - uz ) ** g )
        end do
      end do
    end do

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !> No comment yet!
  !!
  !! TODO Add comment!
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine cumulant_d3q27( fieldProp, inState, outState, auxField, &
    &                        neigh, nElems, nSolve, level, layout,   &
    &                        params, varSys, derVarPos               )
    ! -------------------------------------------------------------------- !
    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> current layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> number of elements in state Array
    integer, intent(in) :: nElems
    !> input  pdf vector
    real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
    !> output pdf vector
    real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
    !> Auxiliary field computed from pre-collision state
    !! Is updated with correct velocity field for multicomponent models
    real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
    !> connectivity vector
    integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
    !> number of elements solved in kernel
    integer, intent(in) :: nSolve
    !> current level
    integer,intent(in) :: level
    !> global parameters
    type(mus_param_type),intent(in) :: params
    !> position of derived quantities in varsys for all fields
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, ii, jj, kk, nScalars
    real(kind=rk) :: ux, uy, uz, rho, inv_rho, omega
    real(kind=rk) :: dxu, dyv, dzw, AA, BB, CC, com_omega
    ! PDF
    real(kind=rk) :: f(-1:1,-1:1,-1:1)
    ! k = central moment, c = cumulant
    real(kind=rk) :: k(0:2,0:2,0:2), c(0:2,0:2,0:2)
    integer :: dens_pos, vel_pos(3), elemOff
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

    nodeloop: do iElem = 1, nSolve

      ! perform streaming step
      f(-1, 0, 0) = inState(neigh (( qn00-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0,-1, 0) = inState(neigh (( q0n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 0,-1) = inState(neigh (( q00n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 0, 0) = inState(neigh (( q100-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 1, 0) = inState(neigh (( q010-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 0, 1) = inState(neigh (( q001-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0,-1,-1) = inState(neigh (( q0nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0,-1, 1) = inState(neigh (( q0n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 1,-1) = inState(neigh (( q01n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 1, 1) = inState(neigh (( q011-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 0,-1) = inState(neigh (( qn0n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 0,-1) = inState(neigh (( q10n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 0, 1) = inState(neigh (( qn01-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 0, 1) = inState(neigh (( q101-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1,-1, 0) = inState(neigh (( qnn0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 1, 0) = inState(neigh (( qn10-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1,-1, 0) = inState(neigh (( q1n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 1, 0) = inState(neigh (( q110-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1,-1,-1) = inState(neigh (( qnnn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1,-1, 1) = inState(neigh (( qnn1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 1,-1) = inState(neigh (( qn1n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 1, 1) = inState(neigh (( qn11-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1,-1,-1) = inState(neigh (( q1nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1,-1, 1) = inState(neigh (( q1n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 1,-1) = inState(neigh (( q11n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 1, 1) = inState(neigh (( q111-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 0, 0) = inState(neigh (( q000-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      ! access rho and velocity from auxField-----------------------------------
      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      inv_rho = 1.0_rk / rho
      ! local x-, y- and z-velocity
      ux = auxField(elemOff + vel_pos(1))
      uy = auxField(elemOff + vel_pos(2))
      uz = auxField(elemOff + vel_pos(3))
      ! ------------------------------------------------------------------------

      ! Central moments ---------------------------------------------------------
      ! eq 43 - 45
      do kk = 0,2
        do jj = 0,2
          do ii = 0,2
            k( ii, jj, kk ) = central_moment_split( f, ii, jj, kk, ux, uy, uz )
          end do
        end do
      end do
      ! Central moments ---------------------------------------------------------


      ! this omega refers to the w1 in the paper
      ! other relaxation parameters are assumed to be 1, thus omitted
      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)

      ! Cumulants and Collision -------------------------------------------------
      ! Not sure how to calcuate the post collision value for 0th and 1st order
      c(0,0,0) = k(0,0,0)
      c(1,0,0) = k(1,0,0)
      c(0,1,0) = k(0,1,0)
      c(0,0,1) = k(0,0,1)
      ! Relaxation for 2nd order terms
      com_omega = ( 1.0_rk - omega )
      c(1,1,0) = com_omega * k(1,1,0)
      c(1,0,1) = com_omega * k(1,0,1)
      c(0,1,1) = com_omega * k(0,1,1)

      ! eq 58 - 60
      Dxu = -omega * 0.5_rk * inv_rho * ( 2.0_rk * k(2,0,0) - k(0,2,0) - k(0,0,2))&
        &   - 0.5_rk * inv_rho * ( k(2,0,0) + k(0,2,0) + k(0,0,2) - k(0,0,0) )
      Dyv = Dxu + 1.5_rk * omega * inv_rho * ( k(2,0,0) - k(0,2,0) )
      Dzw = Dxu + 1.5_rk * omega * inv_rho * ( k(2,0,0) - k(0,0,2) )

      ! eq 61 - 63
      ! first calculate the rights parts
      AA = com_omega * ( k(2,0,0) - k(0,2,0) ) &
        & - 3.0_rk * rho * ( 1.0_rk - omega*0.5_rk ) * ( ux*ux*Dxu - uy*uy*Dyv )
      BB = com_omega * ( k(2,0,0) - k(0,0,2) ) &
        & - 3.0_rk * rho * ( 1.0_rk - omega*0.5_rk ) * ( ux*ux*Dxu - uz*uz*Dzw )
      CC = k(0,0,0) - 1.5_rk * rho * ( ux*ux*Dxu + uy*uy*Dyv + uz*uz*Dzw )

      ! then solve the three moments
      c(2,0,0) = (AA+BB+CC) * div1_3
      c(0,2,0) = (BB+CC) * div1_3 - AA * div2_3
      c(0,0,2) = (AA+CC) * div1_3 - BB * div2_3

      ! for 3rd order or higher, as relaxation == 1, they become simply 0
      c(2,1,0) = 0.0_rk ! 3rd order
      c(2,0,1) = 0.0_rk
      c(1,2,0) = 0.0_rk
      c(0,2,1) = 0.0_rk
      c(1,0,2) = 0.0_rk
      c(0,1,2) = 0.0_rk
      c(1,1,1) = 0.0_rk

      c(2,1,1) = 0.0_rk ! 4th order
      c(1,2,1) = 0.0_rk
      c(1,1,2) = 0.0_rk
      c(2,2,0) = 0.0_rk
      c(2,0,2) = 0.0_rk
      c(0,2,2) = 0.0_rk

      c(1,2,2) = 0.0_rk ! 5th order
      c(2,1,2) = 0.0_rk
      c(2,2,1) = 0.0_rk
      c(2,2,2) = 0.0_rk ! 6th order
      ! Collision --------------------------------------------------------------

      ! Back to central moment -------------------------------------------------
      k = c
      k(1,0,0) = -k(1,0,0)
      k(0,1,0) = -k(0,1,0)
      k(0,0,1) = -k(0,0,1)
      ! get 4th and 6th order from 2nd order
      k(2,1,1) = ( k(2,0,0)*k(0,1,1) + 2.0_rk*k(1,1,0)*k(1,0,1) ) * inv_rho
      k(1,2,1) = ( k(0,2,0)*k(1,0,1) + 2.0_rk*k(1,1,0)*k(0,1,1) ) * inv_rho
      k(1,1,2) = ( k(0,0,2)*k(1,1,0) + 2.0_rk*k(1,0,1)*k(0,1,1) ) * inv_rho
      k(2,2,0) = ( k(2,0,0)*k(0,2,0) + 2.0_rk*k(1,1,0)*k(1,1,0) ) * inv_rho
      k(2,0,2) = ( k(2,0,0)*k(0,0,2) + 2.0_rk*k(1,0,1)*k(1,0,1) ) * inv_rho
      k(0,2,2) = ( k(0,2,0)*k(0,0,2) + 2.0_rk*k(0,1,1)*k(0,1,1) ) * inv_rho
      ! for k 222 half ot the formula was missing!!!! but it has all 0 entries!!!
      k(2,2,2) = - ( 16.0_rk*k(1,1,0)*k(1,0,1)*k(0,1,1) &
        &          + 4.0_rk*(k(1,0,1)**2*k(0,2,0) + k(0,1,1)**2*k(2,0,0) + k(1,1,0)**2*k(0,0,2) ) &
        &          + 2.0_rk*k(2,0,0)*k(0,2,0)*k(0,0,2) ) * inv_rho * inv_rho &
        &        + ( k(2,0,0)*k(0,2,2) + k(0,2,0)*k(2,0,2) + k(0,0,2)*k(2,2,0)  &
        &          + 4.0_rk * ( k(0,1,1)*k(2,1,1) + k(1,0,1)*k(1,2,1) + k(1,1,0)*k(1,1,2) ) &
        &          ) * inv_rho
        ! the previous 3 lines are zeros
      ! Back to central moment -------------------------------------------------

      ! Back to PDF ------------------------------------------------------------
      f = cm_to_pdf( k, ux, uy, uz )
      ! Back to PDF ------------------------------------------------------------

      ! write to state array
      outState( (ielem-1)*nscalars+ q000+(1-1)*qq) = f( 0, 0, 0)
      outState( (ielem-1)*nscalars+ qn00+(1-1)*qq) = f(-1, 0, 0)
      outState( (ielem-1)*nscalars+ q100+(1-1)*qq) = f( 1, 0, 0)
      outState( (ielem-1)*nscalars+ q0n0+(1-1)*qq) = f( 0,-1, 0)
      outState( (ielem-1)*nscalars+ q010+(1-1)*qq) = f( 0, 1, 0)
      outState( (ielem-1)*nscalars+ q00n+(1-1)*qq) = f( 0, 0,-1)
      outState( (ielem-1)*nscalars+ q001+(1-1)*qq) = f( 0, 0, 1)
      outState( (ielem-1)*nscalars+ q0nn+(1-1)*qq) = f( 0,-1,-1)
      outState( (ielem-1)*nscalars+ q01n+(1-1)*qq) = f( 0, 1,-1)
      outState( (ielem-1)*nscalars+ q0n1+(1-1)*qq) = f( 0,-1, 1)
      outState( (ielem-1)*nscalars+ q011+(1-1)*qq) = f( 0, 1, 1)
      outState( (ielem-1)*nscalars+ qn0n+(1-1)*qq) = f(-1, 0,-1)
      outState( (ielem-1)*nscalars+ q10n+(1-1)*qq) = f( 1, 0,-1)
      outState( (ielem-1)*nscalars+ qn01+(1-1)*qq) = f(-1, 0, 1)
      outState( (ielem-1)*nscalars+ q101+(1-1)*qq) = f( 1, 0, 1)
      outState( (ielem-1)*nscalars+ qnn0+(1-1)*qq) = f(-1,-1, 0)
      outState( (ielem-1)*nscalars+ q1n0+(1-1)*qq) = f( 1,-1, 0)
      outState( (ielem-1)*nscalars+ qn10+(1-1)*qq) = f(-1, 1, 0)
      outState( (ielem-1)*nscalars+ q110+(1-1)*qq) = f( 1, 1, 0)
      outState( (ielem-1)*nscalars+ q1nn+(1-1)*qq) = f( 1,-1,-1)
      outState( (ielem-1)*nscalars+ q11n+(1-1)*qq) = f( 1, 1,-1)
      outState( (ielem-1)*nscalars+ q1n1+(1-1)*qq) = f( 1,-1, 1)
      outState( (ielem-1)*nscalars+ q111+(1-1)*qq) = f( 1, 1, 1)
      outState( (ielem-1)*nscalars+ qnnn+(1-1)*qq) = f(-1,-1,-1)
      outState( (ielem-1)*nscalars+ qn1n+(1-1)*qq) = f(-1, 1,-1)
      outState( (ielem-1)*nscalars+ qnn1+(1-1)*qq) = f(-1,-1, 1)
      outState( (ielem-1)*nscalars+ qn11+(1-1)*qq) = f(-1, 1, 1)

    end do nodeloop

  end subroutine cumulant_d3q27
! ****************************************************************************** !

! ****************************************************************************** !
  !> Calculating central moment
  !! This follows equations 6-8 in cumulent paper (Geier .et al 2017)
  pure function kum_ij_g( f, ii, jj, gg, uz, w_ij_g ) result ( kappa )
    !> PDF
    real(kind=rk), intent(in) :: f(-1:1, -1:1, -1:1)
    !> indeces of the pdf
    integer, intent(in) :: ii, jj
    !> order gamma of moments
    integer, intent(in) :: gg
    !> z-component of velocity
    real(kind=rk), intent(in) :: uz
    !> partial weights_ij_g
    real(kind=rk), intent(in) :: w_ij_g(-1:1, -1:1, 0:2)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: kappa, k_0
    ! ---------------------------------------------------------------------------

    k_0 = f(ii,jj,1) + f(ii,jj,-1) + f(ii,jj,0)
    if (gg == 0) then
      kappa = k_0
    else if (gg == 1) then
      kappa = f(ii,jj,1) - f(ii,jj,-1) - uz * (k_0 + w_ij_g( ii, jj, 0 ))
    else
      kappa = f(ii,jj,1) + f(ii,jj,-1) - 2._rk * uz * (f(ii,jj,1) - f(ii,jj,-1)) &
      & + uz**2 * (k_0 + w_ij_g( ii, jj, 0 ))
    endif

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !> Calculating central moment
  !! This follows equations 9-11 in cumulent paper (Geier .et al 2017)
  pure function kum_i_bg( ii, bb, gg, uy, w_i_bg, k_ij_g ) result ( kappa )
    !> indeces of the pdf
    integer, intent(in) :: ii
    !> order gamma of moments
    integer, intent(in) :: bb, gg
    !> y-, z-component of velocity
    real(kind=rk), intent(in) :: uy
    !> partial weights_i_bg
    real(kind=rk), intent(in) :: w_i_bg(-1:1, 0:2, 0:2)
    !> partial cumulants_ij_g
    real(kind=rk), intent(in) :: k_ij_g(-1:1, -1:1, 0:2)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: kappa, k_0
    ! ---------------------------------------------------------------------------

    k_0 = k_ij_g( ii, 1, gg ) + k_ij_g( ii, -1, gg ) + k_ij_g( ii, 0, gg )
    if (bb == 0) then
      kappa = k_0
    else if (bb == 1) then
      kappa = k_ij_g( ii, 1, gg ) - k_ij_g( ii, -1, gg ) - uy * ( k_0 + w_i_bg( ii, 0, gg ) )
    else
      kappa = k_ij_g( ii, 1, gg ) + k_ij_g( ii, -1, gg ) - 2._rk * uy * ( k_ij_g( ii, 1, gg ) &
        & - k_ij_g( ii, -1, gg ) ) + uy**2 * ( k_0 + w_i_bg( ii, 0, gg ) )
    endif

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !> Calculating central moment
  !! This follows equations 12-14 in cumulent paper (Geier .et al 2017)
  pure function kum_abg( aa, bb, gg, ux, w_abg, k_i_bg ) result ( kappa )
    !> order gamma of moments
    integer, intent(in) :: aa, bb, gg
    !> x-, y-, z-component of velocity
    real(kind=rk), intent(in) :: ux
    !> partial weights_abg
    real(kind=rk), intent(in) :: w_abg(0:2, 0:2, 0:2)
    !> partial cumulants_i_bg
    real(kind=rk), intent(in) :: k_i_bg(-1:1, 0:2, 0:2)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: kappa, k_0
    ! ---------------------------------------------------------------------------

    k_0 = k_i_bg( 1, bb, gg) + k_i_bg( -1, bb, gg) + k_i_bg( 0, bb, gg)
    if (aa == 0) then
      kappa = k_0
    else if (aa == 1) then
      kappa = k_i_bg( 1, bb, gg) - k_i_bg( -1, bb, gg) - ux * (k_0 + w_abg( 0, bb, gg ))
    else
      kappa = k_i_bg( 1, bb, gg) + k_i_bg( -1, bb, gg) - 2._rk * ux * (k_i_bg( 1, bb, gg) &
        & - k_i_bg( -1, bb, gg) ) + ux**2 * (k_0 + w_abg( 0, bb, gg ))
    endif

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !> Back to central moment
  !! This follows equations 57-59 in cumulent paper (Geier .et al 2017)
  pure function kpc_i_bg( k, ii, bb, gg, ux, w_abg ) result ( kappa )
    !> cumulant
    real(kind=rk), intent(in) :: k(0:2, 0:2, 0:2)
    !> indeces of the pdf
    integer, intent(in) :: ii
    !> order gamma of moments
    integer, intent(in) :: bb, gg
    !> x-component of velocity
    real(kind=rk), intent(in) :: ux
    !> partial weights_abg
    real(kind=rk), intent(in) :: w_abg(0:2, 0:2, 0:2)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: kappa
    ! ---------------------------------------------------------------------------

    if (ii == 0) then
      kappa = k(0,bb,gg) * (1._rk - ux**2) - 2._rk * ux * k(1,bb,gg) - k(2,bb,gg) &
      & - w_abg(0, bb, gg) * ux**2
    else if (ii == -1) then
      kappa = ( ( k(0,bb,gg) + w_abg(0, bb, gg) ) * ux * (ux - 1._rk) &
        & + k(1,bb,gg) * (2._rk * ux - 1._rk) + k(2,bb,gg) ) * div1_2
    else
      kappa = ( ( k(0,bb,gg) + w_abg(0, bb, gg) ) * ux * (ux + 1._rk) &
        & + k(1,bb,gg) * (2._rk * ux + 1._rk) + k(2,bb,gg) ) * div1_2
    endif

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !> Back to central moment
  !! This follows equations 60-62 in cumulent paper (Geier .et al 2017)
  pure function kpc_ij_g( ii, jj, gg, uy, w_i_bg, k_i_bg ) result ( kappa )
    !> indeces of the pdf
    integer, intent(in) :: ii, jj
    !> order gamma of moments
    integer, intent(in) :: gg
    !> y-component of velocity
    real(kind=rk), intent(in) :: uy
    !> partial weights_i_bg
    real(kind=rk), intent(in) :: w_i_bg(-1:1, 0:2, 0:2)
    !> partial cumulant_i_bg
    real(kind=rk), intent(in) :: k_i_bg(-1:1, 0:2, 0:2)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: kappa
    ! ---------------------------------------------------------------------------

    if (jj == 0) then
      kappa = k_i_bg(ii, 0, gg) * (1._rk - uy**2) - 2._rk * uy * k_i_bg(ii, 1, gg) &
        & - k_i_bg(ii, 2, gg) - w_i_bg(ii, 0, gg) * uy**2
    else if (jj == -1) then
      kappa = ( ( k_i_bg(ii, 0, gg) + w_i_bg(ii, 0, gg) ) * uy * (uy - 1._rk) &
        & + k_i_bg(ii, 1, gg) * (2._rk * uy - 1._rk) + k_i_bg(ii, 2, gg) ) * div1_2
    else
      kappa = ( ( k_i_bg(ii, 0, gg) + w_i_bg(ii, 0, gg) ) * uy * (uy + 1._rk) &
        & + k_i_bg(ii, 1, gg) * (2._rk * uy + 1._rk) + k_i_bg(ii, 2, gg) ) * div1_2
    endif

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !> Back to central moment
  !! This follows equations 63-65 in cumulent paper (Geier .et al 2017)
  pure function kpc_ijk( ii, jj, kk, uz, w_ij_g, k_ij_g ) result ( kappa )
    !> indeces of the pdf
    integer, intent(in) :: ii, jj, kk
    !> x-, y-, z-component of velocity
    real(kind=rk), intent(in) :: uz
    !> partial weights_ij_g
    real(kind=rk), intent(in) :: w_ij_g(-1:1, -1:1, 0:2)
    !> partial cumulants_ij_g
    real(kind=rk), intent(in) :: k_ij_g(-1:1, -1:1, 0:2)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: kappa
    ! ---------------------------------------------------------------------------

    if (kk == 0) then
      kappa = k_ij_g(ii, jj, 0) * (1._rk - uz**2) - 2._rk * uz * k_ij_g(ii, jj, 1) &
        & - k_ij_g(ii, jj, 2) - w_ij_g(ii, jj, 0) * uz**2
    else if (kk == -1) then
      kappa = ( ( k_ij_g(ii, jj, 0) + w_ij_g(ii, jj, 0) ) * uz * (uz - 1._rk) + k_ij_g(ii, jj, 1) &
        & * (2._rk * uz - 1._rk) + k_ij_g(ii, jj, 2) ) * div1_2
    else
      kappa = ( ( k_ij_g(ii, jj, 0) + w_ij_g(ii, jj, 0) ) * uz * (uz + 1._rk) + k_ij_g(ii, jj, 1) &
        & * (2._rk * uz + 1._rk) + k_ij_g(ii, jj, 2) ) * div1_2
    endif

  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !> Cumulant kernel based on Geier2017 and optimized.
  !! Just omega(2) is given in input. omega(2)=-1 means omega2=omegaBulk.
  !! Limiters read from input. lim(N)=10^10 means unlimited.
  !! lim(N) is for omega(N+2). Just omega(3:5) are limited as in the paper.
  !! omega(6:10) = 1._rk
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine cumulant_d3q27_extended_fast( fieldProp, inState, outState, auxField, &
    &                        neigh, nElems, nSolve, level, layout,   &
    &                        params, varSys, derVarPos               )
    ! -------------------------------------------------------------------- !
    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> current layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> number of elements in state Array
    integer, intent(in) :: nElems
    !> input  pdf vector
    real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
    !> output pdf vector
    real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
    !> Auxiliary field computed from pre-collision state
    !! Is updated with correct velocity field for multicomponent models
    real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
    !> connectivity vector
    integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
    !> number of elements solved in kernel
    integer, intent(in) :: nSolve
    !> current level
    integer,intent(in) :: level
    !> global parameters
    type(mus_param_type),intent(in) :: params
    !> position of derived quantities in varsys for all fields
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, ii, jj, kk, nScalars
    real(kind=rk) :: ux, uy, uz, rho, inv_rho, inv_rho_sq, delta_rho, A, B
    real(kind=rk) :: dxu, dyv, dzw, Dxvyu, Dxwzu, Dywzv, AA, BB, CC, par_omega1
    ! PDF
    real(kind=rk) :: f(-1:1,-1:1,-1:1), w(-1:1,-1:1,-1:1), w_i_bg(-1:1,0:2,0:2)
    real(kind=rk) :: w_abg(0:2,0:2,0:2), w_ij_g(-1:1,-1:1,0:2)
    real(kind=rk) :: k_i_bg(-1:1,0:2,0:2), k_ij_g(-1:1,-1:1,0:2)
    ! k = central moment, c = cumulant
    real(kind=rk) :: k(0:2,0:2,0:2), c(0:2,0:2,0:2), omega(10), omega_diff
    integer :: dens_pos, vel_pos(3), elemOff
    !> gradient data
    type(mus_gradData_type), pointer :: gradData
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: inv_omega1, omega_lim, omega_lim_vec(3) !, gradXXU(3,1)
    real(kind=rk) :: comp_omega1, par_omega1_2
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

    ! access gradData
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( 1 )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme
    gradData => scheme%gradData(level)

    omega = fieldProp(1)%fluid%omega_Cum
    omega_lim_vec = fieldProp(1)%fluid%omega_Lim
    if (omega(2) < 0._rk) then
      ! this is omegabulk
      omega(2) = fieldProp(1)%fluid%omegaBulkLvl(level)
    endif
    ! weights for cumulant to pdfs transformation and viceversa
    ! GGS tried to store this in the fluid module, but evaluating
    ! them here results in a faster kernel execution.
    !w      = fieldProp(1)%fluid%w_cumulant
    !w_ij_g = fieldProp(1)%fluid%w_ij_g
    !w_i_bg = fieldProp(1)%fluid%w_i_bg
    !w_abg  = fieldProp(1)%fluid%w_abg
    w = weights_from_layout( layout%weight )
    do kk = 0,2
      do jj = -1,1
        do ii = -1,1
          w_ij_g( ii, jj, kk ) = weights_ijg( ii, jj, kk, w )
        end do
      end do
    end do
    do kk = 0,2
      do jj = 0,2
        do ii = -1,1
          w_i_bg( ii, jj, kk ) = weights_ibg( ii, jj, kk, w_ij_g )
        end do
      end do
    end do
    do kk = 0,2
      do jj = 0,2
        do ii = 0,2
          w_abg( ii, jj, kk ) = weights_abg( ii, jj, kk, w_i_bg )
        end do
      end do
    end do

    ! Test this routine with the following settings.
    ! It should behave like the non extended version
    !omega = 1._rk
    !A = 0._rk
    !B = 0._rk

    !call tem_starttimer(timerHandle = mus_timerHandles%inner_loop)

    nodeloop: do iElem = 1, nSolve

      ! perform streaming step
      f(-1, 0, 0) = inState(neigh (( qn00-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(-1, 0, 0) !layout%weight( qN00 )
      f( 0,-1, 0) = inState(neigh (( q0n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(0,-1, 0) !layout%weight( q0N0 )
      f( 0, 0,-1) = inState(neigh (( q00n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(0, 0,-1) !layout%weight( q00N )
      f( 1, 0, 0) = inState(neigh (( q100-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(1, 0, 0) !layout%weight( q100 )
      f( 0, 1, 0) = inState(neigh (( q010-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(0, 1, 0) !layout%weight( q010 )
      f( 0, 0, 1) = inState(neigh (( q001-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(0, 0, 1) !layout%weight( q001 )
      f( 0,-1,-1) = inState(neigh (( q0nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(0,-1,-1) !layout%weight( q0NN )
      f( 0,-1, 1) = inState(neigh (( q0n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(0,-1, 1) !layout%weight( q0N1 )
      f( 0, 1,-1) = inState(neigh (( q01n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(0, 1,-1) !layout%weight( q01N )
      f( 0, 1, 1) = inState(neigh (( q011-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(0, 1, 1) !layout%weight( q011 )
      f(-1, 0,-1) = inState(neigh (( qn0n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(-1, 0,-1) !layout%weight( qN0N )
      f( 1, 0,-1) = inState(neigh (( q10n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(1, 0,-1) !layout%weight( q10N )
      f(-1, 0, 1) = inState(neigh (( qn01-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(-1, 0, 1) !layout%weight( qN01 )
      f( 1, 0, 1) = inState(neigh (( q101-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(1, 0, 1) !layout%weight( q101 )
      f(-1,-1, 0) = inState(neigh (( qnn0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(-1,-1, 0) !layout%weight( qNN0 )
      f(-1, 1, 0) = inState(neigh (( qn10-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(-1, 1, 0) !layout%weight( qN10 )
      f( 1,-1, 0) = inState(neigh (( q1n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(1,-1, 0) !layout%weight( q1N0 )
      f( 1, 1, 0) = inState(neigh (( q110-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(1, 1, 0) !layout%weight( q110 )
      f(-1,-1,-1) = inState(neigh (( qnnn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(-1,-1,-1) !layout%weight( qNNN )
      f(-1,-1, 1) = inState(neigh (( qnn1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(-1,-1, 1) !layout%weight( qNN1 )
      f(-1, 1,-1) = inState(neigh (( qn1n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(-1, 1,-1) !layout%weight( qN1N )
      f(-1, 1, 1) = inState(neigh (( qn11-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(-1, 1, 1) !layout%weight( qN11 )
      f( 1,-1,-1) = inState(neigh (( q1nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(1,-1,-1) !layout%weight( q1NN )
      f( 1,-1, 1) = inState(neigh (( q1n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(1,-1, 1) !layout%weight( q1N1 )
      f( 1, 1,-1) = inState(neigh (( q11n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(1, 1,-1) !layout%weight( q11N )
      f( 1, 1, 1) = inState(neigh (( q111-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(1, 1, 1) !layout%weight( q111 )
      f( 0, 0, 0) = inState(neigh (( q000-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - w(0, 0, 0) !layout%weight( q000 )

      ! access rho and velocity from auxField-----------------------------------
      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      inv_rho = 1.0_rk / rho
      inv_rho_sq = inv_rho**2
      delta_rho = rho - 1._rk
      ! local x-, y- and z-velocity
      ux = auxField(elemOff + vel_pos(1))
      uy = auxField(elemOff + vel_pos(2))
      uz = auxField(elemOff + vel_pos(3))
      ! ------------------------------------------------------------------------

      ! Central moments ---------------------------------------------------------
      ! eq 6 - 14
      !do kk = 0,2
      !  do jj = 0,2
      !    do ii = 0,2
      !      if (ii + jj + kk < 4) then
      !        k( ii, jj, kk ) = kum_abg( ii, jj, kk, ux, uy, uz, w_ij_g, w_i_bg, w_abg )
      !      endif
      !    end do
      !  end do
      !end do
      !call tem_starttimer(timerHandle = mus_timerHandles%kappa_abg)
      do kk = 0,2
        do jj = -1, 1
          do ii = -1, 1
            k_ij_g(ii, jj, kk) = kum_ij_g(f, ii, jj, kk, uz, w_ij_g)
          end do
        end do
      end do
      do kk = 0,2
        do jj = 0,2
          do ii = -1, 1
            k_i_bg(ii, jj, kk) = kum_i_bg(ii, jj, kk, uy, w_i_bg, k_ij_g)
          end do
        end do
      end do
      ! 0-th order
      k( 0, 0, 0 ) = kum_abg( 0, 0, 0, ux, w_abg, k_i_bg )
      ! 1-st order
      k( 1, 0, 0 ) = kum_abg( 1, 0, 0, ux, w_abg, k_i_bg )
      k( 0, 1, 0 ) = kum_abg( 0, 1, 0, ux, w_abg, k_i_bg )
      k( 0, 0, 1 ) = kum_abg( 0, 0, 1, ux, w_abg, k_i_bg )
      !2-nd order
      k( 1, 1, 0 ) = kum_abg( 1, 1, 0, ux, w_abg, k_i_bg )
      k( 0, 1, 1 ) = kum_abg( 0, 1, 1, ux, w_abg, k_i_bg )
      k( 1, 0, 1 ) = kum_abg( 1, 0, 1, ux, w_abg, k_i_bg )
      k( 2, 0, 0 ) = kum_abg( 2, 0, 0, ux, w_abg, k_i_bg )
      k( 0, 2, 0 ) = kum_abg( 0, 2, 0, ux, w_abg, k_i_bg )
      k( 0, 0, 2 ) = kum_abg( 0, 0, 2, ux, w_abg, k_i_bg )
      !3-rd order
      k( 1, 1, 1 ) = kum_abg( 1, 1, 1, ux, w_abg, k_i_bg )
      k( 2, 1, 0 ) = kum_abg( 2, 1, 0, ux, w_abg, k_i_bg )
      k( 0, 2, 1 ) = kum_abg( 0, 2, 1, ux, w_abg, k_i_bg )
      k( 1, 0, 2 ) = kum_abg( 1, 0, 2, ux, w_abg, k_i_bg )
      k( 0, 1, 2 ) = kum_abg( 0, 1, 2, ux, w_abg, k_i_bg )
      k( 2, 0, 1 ) = kum_abg( 2, 0, 1, ux, w_abg, k_i_bg )
      k( 1, 2, 0 ) = kum_abg( 1, 2, 0, ux, w_abg, k_i_bg )
      !call tem_stoptimer(timerHandle = mus_timerHandles%kappa_abg)
      c = k

      ! 4-th order eq 20 - 21
      !c(2,1,1) = k(2,1,1) - ( (k(2,0,0) + div1_3) * k(0,1,1) &
      !  & + 2._rk * k(1,1,0) * k(1,0,1) ) * inv_rho
      !c(1,2,1) = k(1,2,1) - ( (k(0,2,0) + div1_3) * k(1,0,1) &
      !  & + 2._rk * k(0,1,1) * k(1,1,0) ) * inv_rho
      !c(1,1,2) = k(1,1,2) - ( (k(0,0,2) + div1_3) * k(1,1,0) &
      !  & + 2._rk * k(1,0,1) * k(0,1,1) ) * inv_rho

      !c(2,2,0) = k(2,2,0) - ( ( (k(2,0,0) * k(0,2,0) &
      !  & + 2._rk * k(1,1,0)**2) + ( k(2,0,0) + k(0,2,0) ) * div1_3 ) &
      !  & * inv_rho - delta_rho * inv_rho * div1_9 )
      !c(2,0,2) = k(2,0,2) - ( ( (k(0,0,2) * k(2,0,0) &
      !  & + 2._rk * k(1,0,1)**2) + ( k(0,0,2) + k(2,0,0) ) * div1_3 ) &
      !  & * inv_rho - delta_rho * inv_rho * div1_9 )
      !c(0,2,2) = k(0,2,2) - ( ( (k(0,2,0) * k(0,0,2) &
      !  & + 2._rk * k(0,1,1)**2) + ( k(0,2,0) + k(0,0,2) ) * div1_3 ) &
      !  & * inv_rho - delta_rho * inv_rho * div1_9 )

      ! 5-th order eq 22
      !c(1,2,2) = k(1,2,2) - ( ( k(0,0,2) * k(1,2,0) + k(0,2,0) * k(1,0,2) &
      !  & + 4._rk * k(0,1,1) * k(1,1,1) + 2._rk * ( k(1,0,1) * k(0,2,1) &
      !  & + k(1,1,0) * k(0,1,2) ) ) + (k(1,2,0) + k(1,0,2)) * div1_3 ) &
      !  & * inv_rho
      !c(2,1,2) = k(2,1,2) - ( ( k(2,0,0) * k(0,1,2) + k(0,0,2) * k(2,1,0) &
      !  & + 4._rk * k(1,0,1) * k(1,1,1) + 2._rk * ( k(1,1,0) * k(1,0,2) &
      !  & + k(0,1,1) * k(2,0,1) ) ) + (k(0,1,2) + k(2,1,0)) * div1_3 ) &
      !  & * inv_rho
      !c(2,2,1) = k(2,2,1) - ( ( k(0,2,0) * k(2,0,1) + k(2,0,0) * k(0,2,1) &
      !  & + 4._rk * k(1,1,0) * k(1,1,1) + 2._rk * ( k(0,1,1) * k(2,1,0) &
      !  & + k(1,0,1) * k(1,2,0) ) ) + (k(2,0,1) + k(0,2,1)) * div1_3 ) &
      !  & * inv_rho

      ! 6-th order eq 23
      !c(2,2,2) = k(2,2,2) - (4._rk * k(1,1,1)**2 + k(2,0,0) * k(0,2,2) &
      !  & + k(0,2,0) * k(2,0,2) + k(0,0,2) * k(2,2,0) + 4._rk * (k(0,1,1) &
      !  & * k(2,1,1) + k(1,0,1) * k(1,2,1) + k(1,1,0) * k(1,1,2) ) &
      !  & + 2._rk * (k(1,2,0) * k(1,0,2) + k(2,1,0) * k(0,1,2) + k(2,0,1) &
      !  & * k(0,2,1) ) ) * inv_rho + (16._rk * k(1,1,0) * k(1,0,1) * k(0,1,1) &
      !  & + 4._rk * (k(1,0,1)**2 * k(0,2,0) + k(0,1,1)**2 * k(2,0,0) &
      !  & + k(1,1,0)**2 * k(0,0,2)) + 2._rk * k(2,0,0) * k(0,2,0) * k(0,0,2) ) &
      !  & * inv_rho_sq - (3._rk * (k(0,2,2) + k(2,0,2) + k(2,2,0)) + (k(2,0,0) &
      !  & + k(0,2,0) + k(0,0,2) ) ) * div1_9 * inv_rho + 2._rk * (2._rk &
      !  & * ( k(1,0,1)**2 + k(0,1,1)**2 + k(1,1,0)**2) + (k(0,0,2) * k(0,2,0) &
      !  & + k(0,0,2) * k(2,0,0) + k(0,2,0) * k(2,0,0)) + (k(0,0,2) + k(0,2,0) &
      !  & + k(2,0,0) ) * div1_3 ) * div1_3 * inv_rho_sq &
      !  & + delta_rho * (delta_rho - 1._rk) * div1_27 * inv_rho_sq
      ! Central moments ---------------------------------------------------------

      ! this omegas
      omega(1) = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      inv_omega1 = 1._rk / omega(1)
      ! eq. 111 - 113
      omega(3) = (8._rk * (omega(1) - 2._rk) * (omega(2) * (3._rk * omega(1) - 1._rk) &
        & - 5._rk * omega(1) ) ) / (8._rk * (5._rk - 2._rk * omega(1)) * omega(1) &
        & + omega(2) * (8._rk + omega(1) * (9._rk * omega(1) - 26._rk) ) )
      omega(4) = (8._rk * (omega(1) - 2._rk) * ( omega(1) + omega(2) * (3._rk * omega(1) &
        & - 7._rk) ) ) / ( omega(2) * (56._rk - 42._rk * omega(1) + 9._rk * omega(1)**2 ) &
        & - 8._rk + omega(1) )
      omega(5) = ( 24._rk * ( omega(1) - 2._rk ) * ( 4._rk * omega(1)**2 + omega(1) &
        & * omega(2) * ( 18._rk - 13._rk * omega(1) ) + omega(2)**2 * ( 2._rk + omega(1) &
        & * ( 6._rk * omega(1) - 11._rk ) ) ) ) / ( 16._rk * omega(1)**2 * ( omega(1) - 6._rk ) &
        & - 2._rk * omega(1) * omega(2) * ( 216._rk + 5._rk * omega(1) * (9._rk * omega(1) &
        & - 46._rk ) ) + omega(2)**2 * ( omega(1) * ( 3._rk * omega(1) - 10._rk ) * (15._rk &
        & * omega(1) - 28._rk ) - 48._rk ) )


      ! eq 114 - 115
      A = ( 4._rk * omega(1)**2 + 2._rk * omega(1) * omega(2) * ( omega(1) - 6._rk ) &
        & + omega(2)**2 * ( omega(1) * ( 10._rk - 3._rk * omega(1) ) - 4._rk ) ) &
        & / ( ( omega(1) - omega(2) ) * ( omega(2) * ( 2._rk + 3._rk * omega(1) ) &
        & - 8._rk * omega(1) ) )
      B = ( 4._rk * omega(1) * omega(2) * ( 9._rk * omega(1) - 16._rk ) - 4._rk * omega(1)**2 &
        & - 2._rk * omega(2)**2 * ( 2._rk + 9._rk * omega(1) * ( omega(1) - 2._rk ) ) ) &
        & / ( 3._rk * ( omega(1) - omega(2) ) * ( omega(2) * ( 2._rk + 3._rk * omega(1) ) &
        & - 8._rk * omega(1) ) )


      ! switch to non parametrized Cumulant when omega 3, 4, 5 are not within
      ! the limits 0 < omega < 2
      if (      omega(3) <= 0._rk .or. omega(3) >= 2._rk &
        &  .or. omega(4) <= 0._rk .or. omega(4) >= 2._rk &
        &  .or. omega(5) <= 0._rk .or. omega(5) >= 2._rk ) then
        omega(2:10) = 1._rk
        A = 0._rk
        B = 0._rk
      endif

      ! Post collision Cumulants -------------------------------------------------
      !k = c
      ! 2-nd order terms eq 24 - 26
      comp_omega1 = ( 1.0_rk - omega(1) )
      c(1,1,0) = comp_omega1 * k(1,1,0)
      c(1,0,1) = comp_omega1 * k(1,0,1)
      c(0,1,1) = comp_omega1 * k(0,1,1)

      ! eq 27 - 29
      Dxu = - omega(1) * 0.5_rk * inv_rho * ( 2.0_rk * k(2,0,0) - k(0,2,0) - k(0,0,2)) &
        &   - omega(2) * 0.5_rk * inv_rho * ( k(2,0,0) + k(0,2,0) + k(0,0,2) - k(0,0,0) )
      Dyv = Dxu + 3.0_rk * omega(1) * 0.5_rk * inv_rho * ( k(2,0,0) - k(0,2,0) )
      Dzw = Dxu + 3.0_rk * omega(1) * 0.5_rk * inv_rho * ( k(2,0,0) - k(0,0,2) )

      ! eq 33 - 35 NOT modified according to B.1 - B.3, this is not stable for small viscosities
      ! evaluate second order second derivatives of velocities
      !gradXXU(:,1:1) = scheme%Grad%XXU_ptr(       &
      !     &   auxField     = auxField,           &
      !     &   gradData     = gradData,           &
      !     &   velPos       = vel_pos,            &
      !     &   nAuxScalars  = varSys%nAuxScalars, &
      !     &   nDims        = 3,                  &
      !     &   nSolve       = 1,                  &
      !     &   elemOffset   = iElem - 1           )

      ! first calculate the rights parts
      !par_omega1 = rho * omega(1) * ( 2._rk * ( inv_omega1 - div1_2 )**2 - div1_6 )
      par_omega1_2 = 3.0_rk * rho * ( 1.0_rk - omega(1) * 0.5_rk )
      AA = comp_omega1 * ( k(2,0,0) - k(0,2,0) ) - par_omega1_2 * ( ux*ux*Dxu - uy*uy*Dyv ) !&
        !& + par_omega1 * ( Dxu**2 + ux * gradXXU(1,1) - Dyv**2 - uy * gradXXU(2,1) )
      BB = comp_omega1 * ( k(2,0,0) - k(0,0,2) ) - par_omega1_2 * ( ux*ux*Dxu - uz*uz*Dzw ) !&
        !& + par_omega1 * ( Dxu**2 + ux * gradXXU(1,1) - Dzw**2 - uz * gradXXU(3,1) )
      CC = k(0,0,0) * omega(2) + (1.0_rk - omega(2)) * (k(2,0,0) + k(0,2,0) + k(0,0,2)) &
        & - 3._rk * rho * (1._rk - omega(2) * 0.5_rk) * ( ux*ux*Dxu + uy*uy*Dyv + uz*uz*Dzw ) !&
        !& + rho * ( 6._rk - 3._rk * ( omega(1) + omega(2) ) + omega(1) * omega(2) ) &
        !& * div1_3 * inv_omega1 * ( Dxu**2 + ux * gradXXU(1,1) + Dyv**2 + uy * gradXXU(2,1) &
        !& + Dzw**2 + uz * gradXXU(3,1) )

      ! then solve the three moments
      c(2,0,0) = (AA+BB+CC) * div1_3
      c(0,2,0) = (BB+CC) * div1_3 - AA * div2_3
      c(0,0,2) = (AA+CC) * div1_3 - BB * div2_3

      ! for 3rd order eq 36 - 42 + limiter eq 116-123
      omega_diff = abs(k(1,2,0) + k(1,0,2))
      omega_lim = ( (1._rk - omega(3)) * omega_diff ) &
        & / ( rho * omega_lim_vec(1) + omega_diff )
      AA = (1._rk - (omega(3) + omega_lim) ) * (k(1,2,0) + k(1,0,2))
      omega_diff = abs(k(1,2,0) - k(1,0,2))
      omega_lim = ( (1._rk - omega(4)) * omega_diff ) &
        & / ( rho * omega_lim_vec(2) + omega_diff )
      BB = (1._rk - (omega(4) + omega_lim) ) * (k(1,2,0) - k(1,0,2))
      c(1,2,0) = (AA + BB) * div1_2
      c(1,0,2) = (AA - BB) * div1_2

      omega_diff = abs(k(2,1,0) + k(0,1,2))
      omega_lim = ( (1._rk - omega(3)) * omega_diff ) &
        & / ( rho * omega_lim_vec(1) + omega_diff )
      AA = (1._rk - (omega(3) + omega_lim) ) * (k(2,1,0) + k(0,1,2))
      omega_diff = abs(k(2,1,0) - k(0,1,2))
      omega_lim = ( (1._rk - omega(4)) * omega_diff ) &
        & / ( rho * omega_lim_vec(2) + omega_diff )
      BB = (1._rk - (omega(4) + omega_lim) ) * (k(2,1,0) - k(0,1,2))
      c(2,1,0) = (AA + BB) * div1_2
      c(0,1,2) = (AA - BB) * div1_2

      omega_diff = abs(k(2,0,1) + k(0,2,1))
      omega_lim = ( (1._rk - omega(3)) * omega_diff ) &
        & / ( rho * omega_lim_vec(1) + omega_diff )
      AA = (1._rk - (omega(3) + omega_lim) ) * (k(2,0,1) + k(0,2,1))
      omega_diff = abs(k(2,0,1) - k(0,2,1))
      omega_lim = ( (1._rk - omega(4)) * omega_diff ) &
        & / ( rho * omega_lim_vec(2) + omega_diff )
      BB = (1._rk - (omega(4) + omega_lim) ) * (k(2,0,1) - k(0,2,1))
      c(2,0,1) = (AA + BB) * div1_2
      c(0,2,1) = (AA - BB) * div1_2

      omega_diff = abs(k(1,1,1))
      omega_lim = ( (1._rk - omega(5)) * omega_diff ) &
        & / ( rho * omega_lim_vec(3) + omega_diff )
      c(1,1,1) = (1._rk - (omega(5) + omega_lim) ) * k(1,1,1)

      ! 4th order
      ! eq 30 - 32
      Dxvyu = - 3._rk * omega(1) * inv_rho * k(1,1,0)
      Dxwzu = - 3._rk * omega(1) * inv_rho * k(1,0,1)
      Dywzv = - 3._rk * omega(1) * inv_rho * k(0,1,1)

      ! eq 43 - 45
      par_omega1 = (inv_omega1 - div1_2) * div1_3
      AA = 2._rk * par_omega1 * A * rho * & ! omega(6) *
        & (Dxu - 2._rk * Dyv + Dzw) !+ (1._rk - omega(6)) * (k(2,2,0) &
        !& - 2._rk * k(2,0,2) + k(0,2,2) )
      BB = 2._rk * par_omega1 * A * rho * & ! omega(6) *
        & (Dxu + Dyv - 2._rk * Dzw) !+ (1._rk - omega(6)) * (k(2,2,0) &
        !& + k(2,0,2) - 2._rk * k(0,2,2) )
      CC = - 4._rk * par_omega1 * A * rho * & ! omega(7) *
        & (Dxu + Dyv + Dzw) !+ (1._rk - omega(7)) * (k(2,2,0) + k(2,0,2) + k(0,2,2) )
      c(2,2,0) = (AA + BB + CC) * div1_3
      c(2,0,2) = (CC - AA) * div1_3
      c(0,2,2) = (CC - BB) * div1_3

      ! eq 46 - 48
      c(2,1,1) = - par_omega1 * B * rho * Dywzv !  * omega(8) &
        ! & + (1._rk - omega(8)) * k(2,1,1)
      c(1,2,1) = - par_omega1 * B * rho * Dxwzu !  * omega(8) &
        ! & + (1._rk - omega(8)) * k(1,2,1)
      c(1,1,2) = - par_omega1 * B * rho * Dxvyu !  * omega(8) &
        ! & + (1._rk - omega(8)) * k(1,1,2)

      ! 5th order
      ! eq 49 - 51
      c(2,2,1) = 0._rk !(1._rk - omega(9)) * k(2,2,1)
      c(2,1,2) = 0._rk !(1._rk - omega(9)) * k(2,1,2)
      c(1,2,2) = 0._rk !(1._rk - omega(9)) * k(1,2,2)

      ! 6th order eq 52
      c(2,2,2) = 0._rk !(1._rk - omega(10)) * k(2,2,2)
      ! Collision --------------------------------------------------------------

      ! Back to central moment -------------------------------------------------
      k = c
      ! 1st order with inverted sign
      k(1,0,0) = -c(1,0,0)
      k(0,1,0) = -c(0,1,0)
      k(0,0,1) = -c(0,0,1)
      ! 4th order eq 53-54
      k(2,1,1) = c(2,1,1) + ( (k(2,0,0) + div1_3) * k(0,1,1) + 2.0_rk * k(1,1,0) &
        & * k(1,0,1) ) * inv_rho
      k(1,2,1) = c(1,2,1) + ( (k(0,2,0) + div1_3) * k(1,0,1) + 2.0_rk * k(0,1,1) &
        & * k(1,1,0) ) * inv_rho
      k(1,1,2) = c(1,1,2) + ( (k(0,0,2) + div1_3) * k(1,1,0) + 2.0_rk * k(1,0,1) &
        & * k(0,1,1) ) * inv_rho

      k(2,2,0) = c(2,2,0) + ( ( ( k(2,0,0) * k(0,2,0) + 2.0_rk * k(1,1,0)**2 ) &
        & + (k(2,0,0) + k(0,2,0) ) * div1_3 ) * inv_rho - (delta_rho * inv_rho * div1_9) )
      k(0,2,2) = c(0,2,2) + ( ( ( k(0,2,0) * k(0,0,2) + 2.0_rk * k(0,1,1)**2 ) &
        & + (k(0,2,0) + k(0,0,2) ) * div1_3 ) * inv_rho - (delta_rho * inv_rho * div1_9) )
      k(2,0,2) = c(2,0,2) + ( ( ( k(0,0,2) * k(2,0,0) + 2.0_rk * k(1,0,1)**2 ) &
        & + (k(0,0,2) + k(2,0,0) ) * div1_3 ) * inv_rho - (delta_rho * inv_rho * div1_9) )

      ! 5th order eq 55
      k(1,2,2) = ( ( k(0,0,2) * k(1,2,0) + k(0,2,0) * k(1,0,2) + 4._rk * k(0,1,1) &
        & * k(1,1,1) + 2._rk * (k(1,0,1) * k(0,2,1) + k(1,1,0) * k(0,1,2) ) ) &
        & + ( k(1,2,0) + k(1,0,2) ) * div1_3 ) * inv_rho
      k(2,1,2) = ( ( k(2,0,0) * k(0,1,2) + k(0,0,2) * k(2,1,0) + 4._rk * k(1,0,1) &
        & * k(1,1,1) + 2._rk * (k(1,1,0) * k(1,0,2) + k(0,1,1) * k(2,0,1) ) ) &
        & + ( k(0,1,2) + k(2,1,0) ) * div1_3 ) * inv_rho
      k(2,2,1) = ( ( k(0,2,0) * k(2,0,1) + k(2,0,0) * k(0,2,1) + 4._rk * k(1,1,0) &
        & * k(1,1,1) + 2._rk * (k(0,1,1) * k(2,1,0) + k(1,0,1) * k(1,2,0) ) ) &
        & + ( k(2,0,1) + k(0,2,1) ) * div1_3 ) * inv_rho

      ! 6th order eq 56
      k(2,2,2) = ( 4._rk * k(1,1,1)**2 + k(2,0,0) * k(0,2,2) + k(0,2,0) * k(2,0,2) &
        & + k(0,0,2) * k(2,2,0) + 4.0_rk * ( k(0,1,1) * k(2,1,1) + k(1,0,1) * k(1,2,1) &
        & + k(1,1,0) * k(1,1,2) ) + 2.0_rk * ( k(1,2,0) * k(1,0,2) + k(2,1,0) * k(0,1,2) &
        & + k(2,0,1) * k(0,2,1) ) ) * inv_rho - ( 16.0_rk * k(1,1,0) * k(1,0,1) * k(0,1,1) &
        & + 4.0_rk * ( k(1,0,1)**2 * k(0,2,0) + k(0,1,1)**2 * k(2,0,0) + k(1,1,0)**2 * k(0,0,2) ) &
        & + 2.0_rk * k(2,0,0) * k(0,2,0) * k(0,0,2) ) * inv_rho_sq + ( 3._rk * ( &
        & k(0,2,2) + k(2,0,2) + k(2,2,0) ) + ( k(2,0,0) + k(0,2,0) + k(0,0,2) ) ) * div1_9 &
        & * inv_rho - 2._rk * ( 2._rk * ( k(1,0,1)**2 + k(0,1,1)**2 + k(1,1,0)**2 ) + ( k(0,0,2) &
        & * k(0,2,0) + k(0,0,2) * k(2,0,0) + k(0,2,0) * k(2,0,0) ) + ( k(0,0,2) + k(0,2,0) &
        & + k(2,0,0) ) * div1_3 ) * div1_3 * inv_rho_sq - delta_rho * (delta_rho - 1._rk) * div1_27 &
        & * inv_rho_sq

      ! Back to central moment -------------------------------------------------

      ! Back to PDF ------------------------------------------------------------
      !call tem_starttimer(timerHandle = mus_timerHandles%inv_kappa_abg)
      !f = cm_to_pdf_extended( k, ux, uy, uz, w_ij_g, w_i_bg, w_abg )
      do kk = 0,2
        do jj = 0,2
          do ii = -1, 1
            k_i_bg(ii, jj, kk) = kpc_i_bg(k, ii, jj, kk, ux, w_abg)
          end do
        end do
      end do
      do kk = 0,2
        do jj = -1, 1
          do ii = -1, 1
            k_ij_g(ii, jj, kk) = kpc_ij_g(ii, jj, kk, uy, w_i_bg, k_i_bg)
          end do
        end do
      end do
      do kk = -1,1
        do jj = -1,1
          do ii = -1, 1
            f(ii, jj, kk) = kpc_ijk( ii, jj, kk, uz, w_ij_g, k_ij_g )
          end do
        end do
      end do
      !call tem_stoptimer(timerHandle = mus_timerHandles%inv_kappa_abg)
      ! Back to PDF ------------------------------------------------------------

      ! write to state array
      outState( (ielem-1)*nscalars+ q000+(1-1)*qq) = f( 0, 0, 0) &
        & + w(0, 0, 0) !layout%weight( q000 )
      outState( (ielem-1)*nscalars+ qn00+(1-1)*qq) = f(-1, 0, 0) &
        & + w(-1, 0, 0) !layout%weight( qN00 )
      outState( (ielem-1)*nscalars+ q100+(1-1)*qq) = f( 1, 0, 0) &
        & + w(1, 0, 0) !layout%weight( q100 )
      outState( (ielem-1)*nscalars+ q0n0+(1-1)*qq) = f( 0,-1, 0) &
        & + w(0, -1, 0) !layout%weight( q0N0 )
      outState( (ielem-1)*nscalars+ q010+(1-1)*qq) = f( 0, 1, 0) &
        & + w(0, 1, 0) !layout%weight( q010 )
      outState( (ielem-1)*nscalars+ q00n+(1-1)*qq) = f( 0, 0,-1) &
        & + w(0, 0, -1) !layout%weight( q00N )
      outState( (ielem-1)*nscalars+ q001+(1-1)*qq) = f( 0, 0, 1) &
        & + w(0, 0, 1) !layout%weight( q001 )
      outState( (ielem-1)*nscalars+ q0nn+(1-1)*qq) = f( 0,-1,-1) &
        & + w(0, -1, -1) !layout%weight( q0NN )
      outState( (ielem-1)*nscalars+ q01n+(1-1)*qq) = f( 0, 1,-1) &
        & + w(0, 1, -1) !layout%weight( q01N )
      outState( (ielem-1)*nscalars+ q0n1+(1-1)*qq) = f( 0,-1, 1) &
        & + w(0, -1, 1) !layout%weight( q0N1 )
      outState( (ielem-1)*nscalars+ q011+(1-1)*qq) = f( 0, 1, 1) &
        & + w(0, 1, 1) !layout%weight( q011 )
      outState( (ielem-1)*nscalars+ qn0n+(1-1)*qq) = f(-1, 0,-1) &
        & + w(-1, 0,-1) !layout%weight( qN0N )
      outState( (ielem-1)*nscalars+ q10n+(1-1)*qq) = f( 1, 0,-1) &
        & + w(1, 0,-1) !layout%weight( q10N )
      outState( (ielem-1)*nscalars+ qn01+(1-1)*qq) = f(-1, 0, 1) &
        & + w(-1, 0, 1) !layout%weight( qN01 )
      outState( (ielem-1)*nscalars+ q101+(1-1)*qq) = f( 1, 0, 1) &
        & + w(1, 0, 1) !layout%weight( q101 )
      outState( (ielem-1)*nscalars+ qnn0+(1-1)*qq) = f(-1,-1, 0) &
        & + w(-1,-1, 0) !layout%weight( qNN0 )
      outState( (ielem-1)*nscalars+ q1n0+(1-1)*qq) = f( 1,-1, 0) &
        & + w(1,-1, 0) !layout%weight( q1N0 )
      outState( (ielem-1)*nscalars+ qn10+(1-1)*qq) = f(-1, 1, 0) &
        & + w(-1, 1, 0) !layout%weight( qN10 )
      outState( (ielem-1)*nscalars+ q110+(1-1)*qq) = f( 1, 1, 0) &
        & + w(1, 1, 0) !layout%weight( q110 )
      outState( (ielem-1)*nscalars+ q1nn+(1-1)*qq) = f( 1,-1,-1) &
        & + w(1,-1,-1) !layout%weight( q1NN )
      outState( (ielem-1)*nscalars+ q11n+(1-1)*qq) = f( 1, 1,-1) &
        & + w(1, 1,-1) !layout%weight( q11N )
      outState( (ielem-1)*nscalars+ q1n1+(1-1)*qq) = f( 1,-1, 1) &
        & + w(1,-1, 1) !layout%weight( q1N1 )
      outState( (ielem-1)*nscalars+ q111+(1-1)*qq) = f( 1, 1, 1) &
        & + w(1, 1, 1) !layout%weight( q111 )
      outState( (ielem-1)*nscalars+ qnnn+(1-1)*qq) = f(-1,-1,-1) &
        & + w(-1,-1,-1) !layout%weight( qNNN )
      outState( (ielem-1)*nscalars+ qn1n+(1-1)*qq) = f(-1, 1,-1) &
        & + w(-1, 1,-1) !layout%weight( qN1N )
      outState( (ielem-1)*nscalars+ qnn1+(1-1)*qq) = f(-1,-1, 1) &
        & + w(-1,-1, 1) !layout%weight( qNN1 )
      outState( (ielem-1)*nscalars+ qn11+(1-1)*qq) = f(-1, 1, 1) &
        & + w(-1, 1, 1) !layout%weight( qN11 )

    end do nodeloop

   !call tem_stoptimer(timerHandle = mus_timerHandles%inner_loop)

  end subroutine cumulant_d3q27_extended_fast
! ****************************************************************************** !

! ****************************************************************************** !
  !> Cumulant based on Geier2017 paper. With all modifications and limiters.
  !! All omegas are given in input. When omega(N) = -1 means adjusted in
  !! this routine. omega(2)=-1 means omega2=omegaBulk. Limiters read from
  !! input. lim(N)=10^10 means unlimited.
  !! lim(N) is for omega(N+2). Just omega(3:5) are limited as in the paper.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine cumulant_d3q27_extended_generic( fieldProp, inState, outState, auxField, &
    &                        neigh, nElems, nSolve, level, layout,   &
    &                        params, varSys, derVarPos               )
    ! -------------------------------------------------------------------- !
    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> current layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> number of elements in state Array
    integer, intent(in) :: nElems
    !> input  pdf vector
    real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
    !> output pdf vector
    real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
    !> Auxiliary field computed from pre-collision state
    !! Is updated with correct velocity field for multicomponent models
    real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
    !> connectivity vector
    integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
    !> number of elements solved in kernel
    integer, intent(in) :: nSolve
    !> current level
    integer,intent(in) :: level
    !> global parameters
    type(mus_param_type),intent(in) :: params
    !> position of derived quantities in varsys for all fields
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, ii, jj, kk, nScalars
    real(kind=rk) :: ux, uy, uz, rho, inv_rho, inv_rho_sq, delta_rho, A, B
    real(kind=rk) :: dxu, dyv, dzw, Dxvyu, Dxwzu, Dywzv, AA, BB, CC, par_omega1
    ! PDF
    real(kind=rk) :: f(-1:1,-1:1,-1:1), w(-1:1,-1:1,-1:1), w_i_bg(-1:1,0:2,0:2)
    real(kind=rk) :: w_abg(0:2,0:2,0:2), w_ij_g(-1:1,-1:1,0:2)
    real(kind=rk) :: k_i_bg(-1:1,0:2,0:2), k_ij_g(-1:1,-1:1,0:2)
    ! k = central moment, c = cumulant
    real(kind=rk) :: k(0:2,0:2,0:2), c(0:2,0:2,0:2), omega(10), omega_diff
    integer :: dens_pos, vel_pos(3), elemOff
    !> gradient data
    type(mus_gradData_type), pointer :: gradData
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: gradXXU(3,1), inv_omega1, omega_lim, omega_lim_vec(3)
    real(kind=rk) :: comp_omega1, par_omega1_2
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

    ! access gradData
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( 1 )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme
    gradData => scheme%gradData(level)


    omega = fieldProp(1)%fluid%omega_Cum
    omega_lim_vec = fieldProp(1)%fluid%omega_Lim
    if (omega(2) < 0._rk) then
      ! this is omegabulk
      omega(2) = fieldProp(1)%fluid%omegaBulkLvl(level)
    endif
    !w      = fieldProp(1)%fluid%w_cumulant
    !w_ij_g = fieldProp(1)%fluid%w_ij_g
    !w_i_bg = fieldProp(1)%fluid%w_i_bg
    !w_abg  = fieldProp(1)%fluid%w_abg
    w = weights_from_layout( layout%weight )
    do kk = 0,2
      do jj = -1,1
        do ii = -1,1
          w_ij_g( ii, jj, kk ) = weights_ijg( ii, jj, kk, w )
        end do
      end do
    end do
    do kk = 0,2
      do jj = 0,2
        do ii = -1,1
          w_i_bg( ii, jj, kk ) = weights_ibg( ii, jj, kk, w_ij_g )
        end do
      end do
    end do
    do kk = 0,2
      do jj = 0,2
        do ii = 0,2
          w_abg( ii, jj, kk ) = weights_abg( ii, jj, kk, w_i_bg )
        end do
      end do
    end do

    ! Test this routine with the following settings.
    ! It should behave like the non extended version
    !omega = 1._rk
    !A = 0._rk
    !B = 0._rk

    nodeloop: do iElem = 1, nSolve

      ! perform streaming step
      f(-1, 0, 0) = inState(neigh (( qn00-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( qN00 )
      f( 0,-1, 0) = inState(neigh (( q0n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q0N0 )
      f( 0, 0,-1) = inState(neigh (( q00n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q00N )
      f( 1, 0, 0) = inState(neigh (( q100-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q100 )
      f( 0, 1, 0) = inState(neigh (( q010-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q010 )
      f( 0, 0, 1) = inState(neigh (( q001-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q001 )
      f( 0,-1,-1) = inState(neigh (( q0nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q0NN )
      f( 0,-1, 1) = inState(neigh (( q0n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q0N1 )
      f( 0, 1,-1) = inState(neigh (( q01n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q01N )
      f( 0, 1, 1) = inState(neigh (( q011-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q011 )
      f(-1, 0,-1) = inState(neigh (( qn0n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( qN0N )
      f( 1, 0,-1) = inState(neigh (( q10n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q10N )
      f(-1, 0, 1) = inState(neigh (( qn01-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( qN01 )
      f( 1, 0, 1) = inState(neigh (( q101-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q101 )
      f(-1,-1, 0) = inState(neigh (( qnn0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( qNN0 )
      f(-1, 1, 0) = inState(neigh (( qn10-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( qN10 )
      f( 1,-1, 0) = inState(neigh (( q1n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q1N0 )
      f( 1, 1, 0) = inState(neigh (( q110-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q110 )
      f(-1,-1,-1) = inState(neigh (( qnnn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( qNNN )
      f(-1,-1, 1) = inState(neigh (( qnn1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( qNN1 )
      f(-1, 1,-1) = inState(neigh (( qn1n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( qN1N )
      f(-1, 1, 1) = inState(neigh (( qn11-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( qN11 )
      f( 1,-1,-1) = inState(neigh (( q1nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q1NN )
      f( 1,-1, 1) = inState(neigh (( q1n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q1N1 )
      f( 1, 1,-1) = inState(neigh (( q11n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q11N )
      f( 1, 1, 1) = inState(neigh (( q111-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q111 )
      f( 0, 0, 0) = inState(neigh (( q000-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0) &
        & - layout%weight( q000 )

      ! access rho and velocity from auxField-----------------------------------
      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      inv_rho = 1.0_rk / rho
      inv_rho_sq = inv_rho**2
      delta_rho = rho - 1._rk
      ! local x-, y- and z-velocity
      ux = auxField(elemOff + vel_pos(1))
      uy = auxField(elemOff + vel_pos(2))
      uz = auxField(elemOff + vel_pos(3))
      ! ------------------------------------------------------------------------


      ! Central moments ---------------------------------------------------------
      ! eq 6 - 14
      do kk = 0,2
        do jj = -1, 1
          do ii = -1, 1
            k_ij_g(ii, jj, kk) = kum_ij_g(f, ii, jj, kk, uz, w_ij_g)
          end do
        end do
      end do
      do kk = 0,2
        do jj = 0,2
          do ii = -1, 1
            k_i_bg(ii, jj, kk) = kum_i_bg(ii, jj, kk, uy, w_i_bg, k_ij_g)
          end do
        end do
      end do
      do kk = 0,2
        do jj = 0,2
          do ii = 0,2
            k( ii, jj, kk ) = kum_abg( ii, jj, kk, ux, w_abg, k_i_bg )
          end do
        end do
      end do

      c = k

      ! 4-th order eq 20 - 21
      c(2,1,1) = k(2,1,1) - ( (k(2,0,0) + div1_3) * k(0,1,1) &
        & + 2._rk * k(1,1,0) * k(1,0,1) ) * inv_rho
      c(1,2,1) = k(1,2,1) - ( (k(0,2,0) + div1_3) * k(1,0,1) &
        & + 2._rk * k(0,1,1) * k(1,1,0) ) * inv_rho
      c(1,1,2) = k(1,1,2) - ( (k(0,0,2) + div1_3) * k(1,1,0) &
        & + 2._rk * k(1,0,1) * k(0,1,1) ) * inv_rho

      c(2,2,0) = k(2,2,0) - ( ( (k(2,0,0) * k(0,2,0) &
        & + 2._rk * k(1,1,0)**2) + ( k(2,0,0) + k(0,2,0) ) * div1_3 ) &
        & * inv_rho - delta_rho * inv_rho * div1_9 )
      c(2,0,2) = k(2,0,2) - ( ( (k(0,0,2) * k(2,0,0) &
        & + 2._rk * k(1,0,1)**2) + ( k(0,0,2) + k(2,0,0) ) * div1_3 ) &
        & * inv_rho - delta_rho * inv_rho * div1_9 )
      c(0,2,2) = k(0,2,2) - ( ( (k(0,2,0) * k(0,0,2) &
        & + 2._rk * k(0,1,1)**2) + ( k(0,2,0) + k(0,0,2) ) * div1_3 ) &
        & * inv_rho - delta_rho * inv_rho * div1_9 )

      ! 5-th order eq 22
      c(1,2,2) = k(1,2,2) - ( ( k(0,0,2) * k(1,2,0) + k(0,2,0) * k(1,0,2) &
        & + 4._rk * k(0,1,1) * k(1,1,1) + 2._rk * ( k(1,0,1) * k(0,2,1) &
        & + k(1,1,0) * k(0,1,2) ) ) + (k(1,2,0) + k(1,0,2)) * div1_3 ) &
        & * inv_rho
      c(2,1,2) = k(2,1,2) - ( ( k(2,0,0) * k(0,1,2) + k(0,0,2) * k(2,1,0) &
        & + 4._rk * k(1,0,1) * k(1,1,1) + 2._rk * ( k(1,1,0) * k(1,0,2) &
        & + k(0,1,1) * k(2,0,1) ) ) + (k(0,1,2) + k(2,1,0)) * div1_3 ) &
        & * inv_rho
      c(2,2,1) = k(2,2,1) - ( ( k(0,2,0) * k(2,0,1) + k(2,0,0) * k(0,2,1) &
        & + 4._rk * k(1,1,0) * k(1,1,1) + 2._rk * ( k(0,1,1) * k(2,1,0) &
        & + k(1,0,1) * k(1,2,0) ) ) + (k(2,0,1) + k(0,2,1)) * div1_3 ) &
        & * inv_rho

      ! 6-th order eq 23
      c(2,2,2) = k(2,2,2) - (4._rk * k(1,1,1)**2 + k(2,0,0) * k(0,2,2) &
        & + k(0,2,0) * k(2,0,2) + k(0,0,2) * k(2,2,0) + 4._rk * (k(0,1,1) &
        & * k(2,1,1) + k(1,0,1) * k(1,2,1) + k(1,1,0) * k(1,1,2) ) &
        & + 2._rk * (k(1,2,0) * k(1,0,2) + k(2,1,0) * k(0,1,2) + k(2,0,1) &
        & * k(0,2,1) ) ) * inv_rho + (16._rk * k(1,1,0) * k(1,0,1) * k(0,1,1) &
        & + 4._rk * (k(1,0,1)**2 * k(0,2,0) + k(0,1,1)**2 * k(2,0,0) &
        & + k(1,1,0)**2 * k(0,0,2)) + 2._rk * k(2,0,0) * k(0,2,0) * k(0,0,2) ) &
        & * inv_rho_sq - (3._rk * (k(0,2,2) + k(2,0,2) + k(2,2,0)) + (k(2,0,0) &
        & + k(0,2,0) + k(0,0,2) ) ) * div1_9 * inv_rho + 2._rk * (2._rk &
        & * ( k(1,0,1)**2 + k(0,1,1)**2 + k(1,1,0)**2) + (k(0,0,2) * k(0,2,0) &
        & + k(0,0,2) * k(2,0,0) + k(0,2,0) * k(2,0,0)) + (k(0,0,2) + k(0,2,0) &
        & + k(2,0,0) ) * div1_3 ) * div1_3 * inv_rho_sq &
        & + delta_rho * (delta_rho - 1._rk) * div1_27 * inv_rho_sq
      ! Central moments ---------------------------------------------------------

      ! this omegas
      omega(1) = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      ! eq. 111 - 113
      if (omega(3) < 0._rk) then
        omega(3) = (8._rk * (omega(1) - 2._rk) * (omega(2) * (3._rk * omega(1) - 1._rk) &
          & - 5._rk * omega(1) ) ) / (8._rk * (5._rk - 2._rk * omega(1)) * omega(1) &
          & + omega(2) * (8._rk + omega(1) * (9._rk * omega(1) - 26._rk) ) )
      endif
      if (omega(4) < 0._rk) then
        omega(4) = (8._rk * (omega(1) - 2._rk) * ( omega(1) + omega(2) * (3._rk * omega(1) &
          & - 7._rk) ) ) / ( omega(2) * (56._rk - 42._rk * omega(1) + 9._rk * omega(1)**2 ) &
          & - 8._rk + omega(1) )
      endif
      if (omega(5) < 0._rk) then
        omega(5) = ( 24._rk * ( omega(1) - 2._rk ) * ( 4._rk * omega(1)**2 + omega(1) &
          & * omega(2) * ( 18._rk - 13._rk * omega(1) ) + omega(2)**2 * ( 2._rk + omega(1) &
          & * ( 6._rk * omega(1) - 11._rk ) ) ) ) / ( 16._rk * omega(1)**2 * ( omega(1) - 6._rk ) &
          & - 2._rk * omega(1) * omega(2) * ( 216._rk + 5._rk * omega(1) * (9._rk * omega(1) &
          & - 46._rk ) ) + omega(2)**2 * ( omega(1) * ( 3._rk * omega(1) - 10._rk ) * (15._rk &
          & * omega(1) - 28._rk ) - 48._rk ) )
      endif


      ! eq 114 - 115
      if (omega(5) < 0._rk) then
        A = ( 4._rk * omega(1)**2 + 2._rk * omega(1) * omega(2) * ( omega(1) - 6._rk ) &
          & + omega(2)**2 * ( omega(1) * ( 10._rk - 3._rk * omega(1) ) - 4._rk ) ) &
          & / ( ( omega(1) - omega(2) ) * ( omega(2) * ( 2._rk + 3._rk * omega(1) ) &
          & - 8._rk * omega(1) ) )
        B = ( 4._rk * omega(1) * omega(2) * ( 9._rk * omega(1) - 16._rk ) - 4._rk * omega(1)**2 &
          & - 2._rk * omega(2)**2 * ( 2._rk + 9._rk * omega(1) * ( omega(1) - 2._rk ) ) ) &
          & / ( 3._rk * ( omega(1) - omega(2) ) * ( omega(2) * ( 2._rk + 3._rk * omega(1) ) &
          & - 8._rk * omega(1) ) )
      else
        A = 0._rk
        B = 0._rk
      endif


      ! switch to non parametrized Cumulant when omega 3, 4, 5 are not within
      ! the limits 0 < omega < 2
      if (      omega(3) <= 0._rk .or. omega(3) >= 2._rk &
        &  .or. omega(4) <= 0._rk .or. omega(4) >= 2._rk &
        &  .or. omega(5) <= 0._rk .or. omega(5) >= 2._rk ) then
        omega(2:10) = 1._rk
        A = 0._rk
        B = 0._rk
      endif

      ! Post collision Cumulants -------------------------------------------------
      k = c
      ! 2-nd order terms eq 24 - 26
      comp_omega1 = ( 1.0_rk - omega(1) )
      c(1,1,0) = comp_omega1 * k(1,1,0)
      c(1,0,1) = comp_omega1 * k(1,0,1)
      c(0,1,1) = comp_omega1 * k(0,1,1)

      ! eq 27 - 29
      Dxu = - omega(1) * 0.5_rk * inv_rho * ( 2.0_rk * k(2,0,0) - k(0,2,0) - k(0,0,2)) &
        &   - omega(2) * 0.5_rk * inv_rho * ( k(2,0,0) + k(0,2,0) + k(0,0,2) - k(0,0,0) )
      Dyv = Dxu + 3.0_rk * omega(1) * 0.5_rk * inv_rho * ( k(2,0,0) - k(0,2,0) )
      Dzw = Dxu + 3.0_rk * omega(1) * 0.5_rk * inv_rho * ( k(2,0,0) - k(0,0,2) )

      ! eq 33 - 35 modified according to B.1 - B.3
      ! evaluate second order second derivatives of velocities
      gradXXU(:,1:1) = scheme%Grad%XXU_ptr(       &
           &   auxField     = auxField,           &
           &   gradData     = gradData,           &
           &   velPos       = vel_pos,            &
           &   nAuxScalars  = varSys%nAuxScalars, &
           &   nDims        = 3,                  &
           &   nSolve       = 1,                  &
           &   elemOffset   = iElem - 1           )

      ! first calculate the rights parts
      inv_omega1 = 1._rk / omega(1)
      par_omega1 = rho * omega(1) * ( 2._rk * ( inv_omega1 - div1_2 )**2 - div1_6 )
      par_omega1_2 = 3.0_rk * rho * ( 1.0_rk - omega(1) * 0.5_rk )
      AA = comp_omega1 * ( k(2,0,0) - k(0,2,0) ) &
        & - par_omega1_2 * ( ux*ux*Dxu - uy*uy*Dyv ) &
        & + par_omega1 * ( Dxu**2 + ux * gradXXU(1,1) - Dyv**2 - uy * gradXXU(2,1) )
      BB = comp_omega1 * ( k(2,0,0) - k(0,0,2) ) &
        & - par_omega1_2 * ( ux*ux*Dxu - uz*uz*Dzw ) &
        & + par_omega1 * ( Dxu**2 + ux * gradXXU(1,1) - Dzw**2 - uz * gradXXU(3,1) )
      CC = k(0,0,0) * omega(2) + (1.0_rk - omega(2)) * (k(2,0,0) + k(0,2,0) + k(0,0,2)) &
        & - 3._rk * rho * (1._rk - omega(2) * 0.5_rk) * ( ux*ux*Dxu + uy*uy*Dyv + uz*uz*Dzw ) &
        & + rho * ( 6._rk - 3._rk * ( omega(1) + omega(2) ) + omega(1) * omega(2) ) &
        & * div1_3 * inv_omega1 * ( Dxu**2 + ux * gradXXU(1,1) + Dyv**2 + uy * gradXXU(2,1) &
        & + Dzw**2 + uz * gradXXU(3,1) )

      ! then solve the three moments
      c(2,0,0) = (AA+BB+CC) * div1_3
      c(0,2,0) = (BB+CC) * div1_3 - AA * div2_3
      c(0,0,2) = (AA+CC) * div1_3 - BB * div2_3

      ! for 3rd order eq 36 - 42 + limiter eq 116-123
      omega_diff = abs(k(1,2,0) + k(1,0,2))
      omega_lim = ( (1._rk - omega(3)) * omega_diff ) &
        & / ( rho * omega_lim_vec(1) + omega_diff )
      AA = (1._rk - (omega(3) + omega_lim) ) * (k(1,2,0) + k(1,0,2))
      omega_diff = abs(k(1,2,0) - k(1,0,2))
      omega_lim = ( (1._rk - omega(4)) * omega_diff ) &
        & / ( rho * omega_lim_vec(2) + omega_diff )
      BB = (1._rk - (omega(4) + omega_lim) ) * (k(1,2,0) - k(1,0,2))
      c(1,2,0) = (AA + BB) * div1_2
      c(1,0,2) = (AA - BB) * div1_2

      omega_diff = abs(k(2,1,0) + k(0,1,2))
      omega_lim = ( (1._rk - omega(3)) * omega_diff ) &
        & / ( rho * omega_lim_vec(1) + omega_diff )
      AA = (1._rk - (omega(3) + omega_lim) ) * (k(2,1,0) + k(0,1,2))
      omega_diff = abs(k(2,1,0) - k(0,1,2))
      omega_lim = ( (1._rk - omega(4)) * omega_diff ) &
        & / ( rho * omega_lim_vec(2) + omega_diff )
      BB = (1._rk - (omega(4) + omega_lim) ) * (k(2,1,0) - k(0,1,2))
      c(2,1,0) = (AA + BB) * div1_2
      c(0,1,2) = (AA - BB) * div1_2

      omega_diff = abs(k(2,0,1) + k(0,2,1))
      omega_lim = ( (1._rk - omega(3)) * omega_diff ) &
        & / ( rho * omega_lim_vec(1) + omega_diff )
      AA = (1._rk - (omega(3) + omega_lim) ) * (k(2,0,1) + k(0,2,1))
      omega_diff = abs(k(2,0,1) - k(0,2,1))
      omega_lim = ( (1._rk - omega(4)) * omega_diff ) &
        & / ( rho * omega_lim_vec(2) + omega_diff )
      BB = (1._rk - (omega(4) + omega_lim) ) * (k(2,0,1) - k(0,2,1))
      c(2,0,1) = (AA + BB) * div1_2
      c(0,2,1) = (AA - BB) * div1_2

      omega_diff = abs(k(1,1,1))
      omega_lim = ( (1._rk - omega(5)) * omega_diff ) &
        & / ( rho * omega_lim_vec(3) + omega_diff )
      c(1,1,1) = (1._rk - (omega(5) + omega_lim) ) * k(1,1,1)

      ! 4th order
      ! eq 30 - 32
      Dxvyu = - 3._rk * omega(1) * inv_rho * k(1,1,0)
      Dxwzu = - 3._rk * omega(1) * inv_rho * k(1,0,1)
      Dywzv = - 3._rk * omega(1) * inv_rho * k(0,1,1)

      ! eq 43 - 45
      par_omega1 = (inv_omega1 - div1_2) * div1_3
      AA = 2._rk * par_omega1 * omega(6) * A * rho * (Dxu &
        & - 2._rk * Dyv + Dzw) + (1._rk - omega(6)) * (k(2,2,0) &
        & - 2._rk * k(2,0,2) + k(0,2,2) )
      BB = 2._rk * par_omega1 * omega(6) * A * rho * (Dxu &
        & + Dyv - 2._rk * Dzw) + (1._rk - omega(6)) * (k(2,2,0) &
        & + k(2,0,2) - 2._rk * k(0,2,2) )
      CC = - 4._rk * par_omega1 * omega(7) * A * rho * (Dxu &
        & + Dyv + Dzw) + (1._rk - omega(7)) * (k(2,2,0) + k(2,0,2) + k(0,2,2) )
      c(2,2,0) = (AA + BB + CC) * div1_3
      c(2,0,2) = (CC - AA) * div1_3
      c(0,2,2) = (CC - BB) * div1_3

      ! eq 46 - 48
      c(2,1,1) = - par_omega1 * omega(8) * B * rho * Dywzv &
        & + (1._rk - omega(8)) * k(2,1,1)
      c(1,2,1) = - par_omega1 * omega(8) * B * rho * Dxwzu &
        & + (1._rk - omega(8)) * k(1,2,1)
      c(1,1,2) = - par_omega1 * omega(8) * B * rho * Dxvyu &
        & + (1._rk - omega(8)) * k(1,1,2)

      ! 5th order
      ! eq 49 - 51
      c(2,2,1) = (1._rk - omega(9)) * k(2,2,1)
      c(2,1,2) = (1._rk - omega(9)) * k(2,1,2)
      c(1,2,2) = (1._rk - omega(9)) * k(1,2,2)

      ! 6th order eq 52
      c(2,2,2) = (1._rk - omega(10)) * k(2,2,2)
      ! Collision --------------------------------------------------------------

      ! Back to central moment -------------------------------------------------
      k = c
      ! 1st order with inverted sign
      k(1,0,0) = -c(1,0,0)
      k(0,1,0) = -c(0,1,0)
      k(0,0,1) = -c(0,0,1)
      ! 4th order eq 53-54
      k(2,1,1) = c(2,1,1) + ( (k(2,0,0) + div1_3) * k(0,1,1) + 2.0_rk * k(1,1,0) &
        & * k(1,0,1) ) * inv_rho
      k(1,2,1) = c(1,2,1) + ( (k(0,2,0) + div1_3) * k(1,0,1) + 2.0_rk * k(0,1,1) &
        & * k(1,1,0) ) * inv_rho
      k(1,1,2) = c(1,1,2) + ( (k(0,0,2) + div1_3) * k(1,1,0) + 2.0_rk * k(1,0,1) &
        & * k(0,1,1) ) * inv_rho

      k(2,2,0) = c(2,2,0) + ( ( ( k(2,0,0) * k(0,2,0) + 2.0_rk * k(1,1,0)**2 ) &
        & + (k(2,0,0) + k(0,2,0) ) * div1_3 ) * inv_rho - (delta_rho * inv_rho * div1_9) )
      k(0,2,2) = c(0,2,2) + ( ( ( k(0,2,0) * k(0,0,2) + 2.0_rk * k(0,1,1)**2 ) &
        & + (k(0,2,0) + k(0,0,2) ) * div1_3 ) * inv_rho - (delta_rho * inv_rho * div1_9) )
      k(2,0,2) = c(2,0,2) + ( ( ( k(0,0,2) * k(2,0,0) + 2.0_rk * k(1,0,1)**2 ) &
        & + (k(0,0,2) + k(2,0,0) ) * div1_3 ) * inv_rho - (delta_rho * inv_rho * div1_9) )

      ! 5th order eq 55
      k(1,2,2) = c(1,2,2) + ( ( k(0,0,2) * k(1,2,0) + k(0,2,0) * k(1,0,2) + 4._rk * k(0,1,1) &
        & * k(1,1,1) + 2._rk * (k(1,0,1) * k(0,2,1) + k(1,1,0) * k(0,1,2) ) ) &
        & + ( k(1,2,0) + k(1,0,2) ) * div1_3 ) * inv_rho
      k(2,1,2) = c(2,1,2) + ( ( k(2,0,0) * k(0,1,2) + k(0,0,2) * k(2,1,0) + 4._rk * k(1,0,1) &
        & * k(1,1,1) + 2._rk * (k(1,1,0) * k(1,0,2) + k(0,1,1) * k(2,0,1) ) ) &
        & + ( k(0,1,2) + k(2,1,0) ) * div1_3 ) * inv_rho
      k(2,2,1) = c(2,2,1) + ( ( k(0,2,0) * k(2,0,1) + k(2,0,0) * k(0,2,1) + 4._rk * k(1,1,0) &
        & * k(1,1,1) + 2._rk * (k(0,1,1) * k(2,1,0) + k(1,0,1) * k(1,2,0) ) ) &
        & + ( k(2,0,1) + k(0,2,1) ) * div1_3 ) * inv_rho

      ! 6th order eq 56
      k(2,2,2) = c(2,2,2) + ( 4._rk * k(1,1,1)**2 + k(2,0,0) * k(0,2,2) + k(0,2,0) * k(2,0,2) &
        & + k(0,0,2) * k(2,2,0) + 4.0_rk * ( k(0,1,1) * k(2,1,1) + k(1,0,1) * k(1,2,1) &
        & + k(1,1,0) * k(1,1,2) ) + 2.0_rk * ( k(1,2,0) * k(1,0,2) + k(2,1,0) * k(0,1,2) &
        & + k(2,0,1) * k(0,2,1) ) ) * inv_rho - ( 16.0_rk * k(1,1,0) * k(1,0,1) * k(0,1,1) &
        & + 4.0_rk * ( k(1,0,1)**2 * k(0,2,0) + k(0,1,1)**2 * k(2,0,0) + k(1,1,0)**2 * k(0,0,2) ) &
        & + 2.0_rk * k(2,0,0) * k(0,2,0) * k(0,0,2) ) * inv_rho_sq + ( 3._rk * ( &
        & k(0,2,2) + k(2,0,2) + k(2,2,0) ) + ( k(2,0,0) + k(0,2,0) + k(0,0,2) ) ) * div1_9 &
        & * inv_rho - 2._rk * ( 2._rk * ( k(1,0,1)**2 + k(0,1,1)**2 + k(1,1,0)**2 ) + ( k(0,0,2) &
        & * k(0,2,0) + k(0,0,2) * k(2,0,0) + k(0,2,0) * k(2,0,0) ) + ( k(0,0,2) + k(0,2,0) &
        & + k(2,0,0) ) * div1_3 ) * div1_3 * inv_rho_sq - delta_rho * (delta_rho - 1._rk) * div1_27 &
        & * inv_rho_sq

      ! Back to central moment -------------------------------------------------

      ! Back to PDF ------------------------------------------------------------
      !f = cm_to_pdf_extended( k, ux, uy, uz, w_ij_g, w_i_bg, w_abg )
      do kk = 0,2
        do jj = 0,2
          do ii = -1, 1
            k_i_bg(ii, jj, kk) = kpc_i_bg(k, ii, jj, kk, ux, w_abg)
          end do
        end do
      end do
      do kk = 0,2
        do jj = -1, 1
          do ii = -1, 1
            k_ij_g(ii, jj, kk) = kpc_ij_g(ii, jj, kk, uy, w_i_bg, k_i_bg)
          end do
        end do
      end do
      do kk = -1,1
        do jj = -1,1
          do ii = -1, 1
            f(ii, jj, kk) = kpc_ijk( ii, jj, kk, uz, w_ij_g, k_ij_g )
          end do
        end do
      end do
      ! Back to PDF ------------------------------------------------------------



      ! write to state array
      outState( (ielem-1)*nscalars+ q000+(1-1)*qq) = f( 0, 0, 0) &
        & + layout%weight( q000 )
      outState( (ielem-1)*nscalars+ qn00+(1-1)*qq) = f(-1, 0, 0) &
        & + layout%weight( qN00 )
      outState( (ielem-1)*nscalars+ q100+(1-1)*qq) = f( 1, 0, 0) &
        & + layout%weight( q100 )
      outState( (ielem-1)*nscalars+ q0n0+(1-1)*qq) = f( 0,-1, 0) &
        & + layout%weight( q0N0 )
      outState( (ielem-1)*nscalars+ q010+(1-1)*qq) = f( 0, 1, 0) &
        & + layout%weight( q010 )
      outState( (ielem-1)*nscalars+ q00n+(1-1)*qq) = f( 0, 0,-1) &
        & + layout%weight( q00N )
      outState( (ielem-1)*nscalars+ q001+(1-1)*qq) = f( 0, 0, 1) &
        & + layout%weight( q001 )
      outState( (ielem-1)*nscalars+ q0nn+(1-1)*qq) = f( 0,-1,-1) &
        & + layout%weight( q0NN )
      outState( (ielem-1)*nscalars+ q01n+(1-1)*qq) = f( 0, 1,-1) &
        & + layout%weight( q01N )
      outState( (ielem-1)*nscalars+ q0n1+(1-1)*qq) = f( 0,-1, 1) &
        & + layout%weight( q0N1 )
      outState( (ielem-1)*nscalars+ q011+(1-1)*qq) = f( 0, 1, 1) &
        & + layout%weight( q011 )
      outState( (ielem-1)*nscalars+ qn0n+(1-1)*qq) = f(-1, 0,-1) &
        & + layout%weight( qN0N )
      outState( (ielem-1)*nscalars+ q10n+(1-1)*qq) = f( 1, 0,-1) &
        & + layout%weight( q10N )
      outState( (ielem-1)*nscalars+ qn01+(1-1)*qq) = f(-1, 0, 1) &
        & + layout%weight( qN01 )
      outState( (ielem-1)*nscalars+ q101+(1-1)*qq) = f( 1, 0, 1) &
        & + layout%weight( q101 )
      outState( (ielem-1)*nscalars+ qnn0+(1-1)*qq) = f(-1,-1, 0) &
        & + layout%weight( qNN0 )
      outState( (ielem-1)*nscalars+ q1n0+(1-1)*qq) = f( 1,-1, 0) &
        & + layout%weight( q1N0 )
      outState( (ielem-1)*nscalars+ qn10+(1-1)*qq) = f(-1, 1, 0) &
        & + layout%weight( qN10 )
      outState( (ielem-1)*nscalars+ q110+(1-1)*qq) = f( 1, 1, 0) &
        & + layout%weight( q110 )
      outState( (ielem-1)*nscalars+ q1nn+(1-1)*qq) = f( 1,-1,-1) &
        & + layout%weight( q1NN )
      outState( (ielem-1)*nscalars+ q11n+(1-1)*qq) = f( 1, 1,-1) &
        & + layout%weight( q11N )
      outState( (ielem-1)*nscalars+ q1n1+(1-1)*qq) = f( 1,-1, 1) &
        & + layout%weight( q1N1 )
      outState( (ielem-1)*nscalars+ q111+(1-1)*qq) = f( 1, 1, 1) &
        & + layout%weight( q111 )
      outState( (ielem-1)*nscalars+ qnnn+(1-1)*qq) = f(-1,-1,-1) &
        & + layout%weight( qNNN )
      outState( (ielem-1)*nscalars+ qn1n+(1-1)*qq) = f(-1, 1,-1) &
        & + layout%weight( qN1N )
      outState( (ielem-1)*nscalars+ qnn1+(1-1)*qq) = f(-1,-1, 1) &
        & + layout%weight( qNN1 )
      outState( (ielem-1)*nscalars+ qn11+(1-1)*qq) = f(-1, 1, 1) &
        & + layout%weight( qN11 )

    end do nodeloop

  end subroutine cumulant_d3q27_extended_generic
! ****************************************************************************** !

! ****************************************************************************** !
  pure function cm_to_pdf( k, ux, uy, uz ) result ( f )
    real(kind=rk), intent(in) :: k(0:2,0:2,0:2)
    real(kind=rk), intent(in) :: ux, uy, uz

    real(kind=rk) :: f(-1:1,-1:1,-1:1)
    real(kind=rk) :: g(-1:1, 0:2, 0:2)
    real(kind=rk) :: h(-1:1,-1:1, 0:2)
    integer :: ii, jj
    real(kind=rk) :: uu, mm, nn, pp, qq, rr, u2

    uu = ux*ux
    mm = 1.0_rk - uu
    nn = uu + ux
    pp = uu - ux
    u2 = 2.0_rk * ux
    qq = u2 - 1.0_rk
    rr = u2 + 1.0_rk
    do jj = 0,2
      do ii = 0,2
        g( 0,ii,jj) =   k(0,ii,jj)*mm - k(1,ii,jj)*u2 - k(2,ii,jj)
        g(-1,ii,jj) = ( k(0,ii,jj)*pp + k(1,ii,jj)*qq + k(2,ii,jj) ) * div1_2
        g( 1,ii,jj) = ( k(0,ii,jj)*nn + k(1,ii,jj)*rr + k(2,ii,jj) ) * div1_2
      end do
    end do

    uu = uy*uy
    mm = 1.0_rk - uu
    nn = uu + uy
    pp = uu - uy
    u2 = 2.0_rk * uy
    qq = u2 - 1.0_rk
    rr = u2 + 1.0_rk
    do jj = 0,2
      do ii = -1,1
        h(ii, 0,jj) =   g(ii,0,jj)*mm - g(ii,1,jj)*u2 - g(ii,2,jj)
        h(ii,-1,jj) = ( g(ii,0,jj)*pp + g(ii,1,jj)*qq + g(ii,2,jj) ) * div1_2
        h(ii, 1,jj) = ( g(ii,0,jj)*nn + g(ii,1,jj)*rr + g(ii,2,jj) ) * div1_2
      end do
    end do

    uu = uz*uz
    mm = 1.0_rk - uu
    nn = uu + uz
    pp = uu - uz
    u2 = 2.0_rk * uz
    qq = u2 - 1.0_rk
    rr = u2 + 1.0_rk
    do jj = -1,1
      do ii = -1,1
        f(ii, jj, 0) =   h(ii,jj,0)*mm - h(ii,jj,1)*u2 - h(ii,jj,2)
        f(ii, jj,-1) = ( h(ii,jj,0)*pp + h(ii,jj,1)*qq + h(ii,jj,2) ) * div1_2
        f(ii, jj, 1) = ( h(ii,jj,0)*nn + h(ii,jj,1)*rr + h(ii,jj,2) ) * div1_2
      end do
    end do
  end function
! ****************************************************************************** !

! ****************************************************************************** !
  !> No comment yet!
  !!
  !! TODO Add coment!
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine cascaded_d3q27( fieldProp, inState, outState, auxField, &
    &                        neigh, nElems, nSolve, level, layout,   &
    &                        params, varSys, derVarPos               )
    ! -------------------------------------------------------------------- !
    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> current layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> number of elements in state Array
    integer, intent(in) :: nElems
    !> input  pdf vector
    real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
    !> output pdf vector
    real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
    !> Auxiliary field computed from pre-collision state
    !! Is updated with correct velocity field for multicomponent models
    real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
    !> connectivity vector
    integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
    !> number of elements solved in kernel
    integer, intent(in) :: nSolve
    !> current level
    integer,intent(in) :: level
    !> global parameters
    type(mus_param_type),intent(in) :: params
    !> position of derived quantities in varsys for all fields
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, ii, jj, kk, nScalars
    real(kind=rk) :: ux, uy, uz, rho, omega, inv_rho
    real(kind=rk) :: dxu, dyv, dzw, AA, BB, CC, com_omega
    real(kind=rk) ::  f(-1:1,-1:1,-1:1)
    real(kind=rk) :: ff(-1:1,-1:1,-1:1)
    ! k = central moment, c = cumulant
    real(kind=rk) :: k(0:2,0:2,0:2), c(0:2,0:2,0:2)
    integer :: dens_pos, vel_pos(3), elemOff
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

    nodeloop: do iElem = 1, nSolve

      f(-1, 0, 0) = inState(neigh (( qn00-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0,-1, 0) = inState(neigh (( q0n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 0,-1) = inState(neigh (( q00n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 0, 0) = inState(neigh (( q100-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 1, 0) = inState(neigh (( q010-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 0, 1) = inState(neigh (( q001-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0,-1,-1) = inState(neigh (( q0nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0,-1, 1) = inState(neigh (( q0n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 1,-1) = inState(neigh (( q01n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 1, 1) = inState(neigh (( q011-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 0,-1) = inState(neigh (( qn0n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 0,-1) = inState(neigh (( q10n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 0, 1) = inState(neigh (( qn01-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 0, 1) = inState(neigh (( q101-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1,-1, 0) = inState(neigh (( qnn0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 1, 0) = inState(neigh (( qn10-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1,-1, 0) = inState(neigh (( q1n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 1, 0) = inState(neigh (( q110-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1,-1,-1) = inState(neigh (( qnnn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1,-1, 1) = inState(neigh (( qnn1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 1,-1) = inState(neigh (( qn1n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f(-1, 1, 1) = inState(neigh (( qn11-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1,-1,-1) = inState(neigh (( q1nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1,-1, 1) = inState(neigh (( q1n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 1,-1) = inState(neigh (( q11n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 1, 1, 1) = inState(neigh (( q111-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f( 0, 0, 0) = inState(neigh (( q000-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      ! calculate rho and velocity ----------------------------------------------
      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      inv_rho = 1.0_rk / rho
      ! local x-, y- and z-velocity
      ux = auxField(elemOff + vel_pos(1))
      uy = auxField(elemOff + vel_pos(2))
      uz = auxField(elemOff + vel_pos(3))
      ! calculate rho and velocity ----------------------------------------------

      ! get all central moments
      ! eq 43 - 45
      do kk = 0,2
        do jj = 0,2
          do ii = 0,2
            k( ii, jj, kk ) = central_moment_split( f, ii, jj, kk, ux, uy, uz )
          end do
        end do
      end do

      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)

      ! Collision --------------------------------------------------------------
      c(0,0,0) = k(0,0,0)
      c(1,0,0) = k(1,0,0)
      c(0,1,0) = k(0,1,0)
      c(0,0,1) = k(0,0,1)
      com_omega = ( 1.0_rk - omega )
      c(1,1,0) = com_omega * k(1,1,0)
      c(1,0,1) = com_omega * k(1,0,1)
      c(0,1,1) = com_omega * k(0,1,1)

      ! eq D.1 - D.3
      Dxu = -omega * 0.5_rk * inv_rho * ( 2.0_rk * k(2,0,0) - k(0,2,0) - k(0,0,2))&
        &   - 0.5_rk * inv_rho * ( k(2,0,0) + k(0,2,0) + k(0,0,2) - rho )
      Dyv = Dxu + 1.5_rk * omega * inv_rho * ( k(2,0,0) - k(0,2,0) )
      Dzw = Dxu + 1.5_rk * omega * inv_rho * ( k(2,0,0) - k(0,0,2) )

      ! eq D.4 - D.6
      ! first calculate the rights parts
      AA = com_omega * ( k(2,0,0) - k(0,2,0) ) &
        & - 3.0_rk * rho * ( 1.0_rk - omega*0.5_rk ) * ( ux*ux*Dxu - uy*uy*Dyv )
      BB = com_omega * ( k(2,0,0) - k(0,0,2) ) &
        & - 3.0_rk * rho * ( 1.0_rk - omega*0.5_rk ) * ( ux*ux*Dxu - uz*uz*Dzw )
      CC = rho - 1.5_rk * rho * ( ux*ux*Dxu + uy*uy*Dyv + uz*uz*Dzw )

      ! then solve the three moments
      c(2,0,0) = (AA+BB+CC) * div1_3
      c(0,2,0) = (BB+CC) * div1_3 - AA * div2_3
      c(0,0,2) = (AA+CC) * div1_3 - BB * div2_3

      ! for 3rd order or higher, set to 0
      c(2,1,0) = 0.0_rk ! 3rd order
      c(2,0,1) = 0.0_rk
      c(1,2,0) = 0.0_rk
      c(0,2,1) = 0.0_rk
      c(1,0,2) = 0.0_rk
      c(0,1,2) = 0.0_rk
      c(1,1,1) = 0.0_rk

      c(2,1,1) = 0.0_rk ! 4th order
      c(1,2,1) = 0.0_rk
      c(1,1,2) = 0.0_rk
      c(2,2,0) = rho * div1_9
      c(2,0,2) = rho * div1_9
      c(0,2,2) = rho * div1_9

      c(1,2,2) = 0.0_rk ! 5th order
      c(2,1,2) = 0.0_rk
      c(2,2,1) = 0.0_rk
      c(2,2,2) = rho * div1_27 ! 6th order
      ! Collision --------------------------------------------------------------

      ! Back to PDF ------------------------------------------------------------
      ff = cm_to_pdf( c, ux, uy, uz )
      ! Back to PDF ------------------------------------------------------------

      ! write to state array
outState( (ielem-1)*qq+ q000+(1-1)*qq) = ff( 0,  0,  0)
outState( (ielem-1)*qq+ qn00+(1-1)*qq) = ff(-1,  0,  0)
outState( (ielem-1)*qq+ q100+(1-1)*qq) = ff( 1,  0,  0)
outState( (ielem-1)*qq+ q0n0+(1-1)*qq) = ff( 0, -1,  0)
outState( (ielem-1)*qq+ q010+(1-1)*qq) = ff( 0,  1,  0)
outState( (ielem-1)*qq+ q00n+(1-1)*qq) = ff( 0,  0, -1)
outState( (ielem-1)*qq+ q001+(1-1)*qq) = ff( 0,  0,  1)
outState( (ielem-1)*qq+ q0nn+(1-1)*qq) = ff( 0, -1, -1)
outState( (ielem-1)*qq+ q01n+(1-1)*qq) = ff( 0,  1, -1)
outState( (ielem-1)*qq+ q0n1+(1-1)*qq) = ff( 0, -1,  1)
outState( (ielem-1)*qq+ q011+(1-1)*qq) = ff( 0,  1,  1)
outState( (ielem-1)*qq+ qn0n+(1-1)*qq) = ff(-1,  0, -1)
outState( (ielem-1)*qq+ q10n+(1-1)*qq) = ff( 1,  0, -1)
outState( (ielem-1)*qq+ qn01+(1-1)*qq) = ff(-1,  0,  1)
outState( (ielem-1)*qq+ q101+(1-1)*qq) = ff( 1,  0,  1)
outState( (ielem-1)*qq+ qnn0+(1-1)*qq) = ff(-1, -1,  0)
outState( (ielem-1)*qq+ q1n0+(1-1)*qq) = ff( 1, -1,  0)
outState( (ielem-1)*qq+ qn10+(1-1)*qq) = ff(-1,  1,  0)
outState( (ielem-1)*qq+ q110+(1-1)*qq) = ff( 1,  1,  0)
outState( (ielem-1)*qq+ q1nn+(1-1)*qq) = ff( 1, -1, -1)
outState( (ielem-1)*qq+ q11n+(1-1)*qq) = ff( 1,  1, -1)
outState( (ielem-1)*qq+ q1n1+(1-1)*qq) = ff( 1, -1,  1)
outState( (ielem-1)*qq+ q111+(1-1)*qq) = ff( 1,  1,  1)
outState( (ielem-1)*qq+ qnnn+(1-1)*qq) = ff(-1, -1, -1)
outState( (ielem-1)*qq+ qn1n+(1-1)*qq) = ff(-1,  1, -1)
outState( (ielem-1)*qq+ qnn1+(1-1)*qq) = ff(-1, -1,  1)
outState( (ielem-1)*qq+ qn11+(1-1)*qq) = ff(-1,  1,  1)

    end do nodeloop

  end subroutine cascaded_d3q27
! ****************************************************************************** !

end module mus_compute_cumulant_module
