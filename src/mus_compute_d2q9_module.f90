! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011, 2017, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2016, 2018-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2021-2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
! Copyright (c) 2021 Tobias Horstmann <tobias.horstmann@dlr.de>
! Copyright (c) 2021 Jana Gericke <jana.gericke@dlr.de>
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
! **************************************************************************** !
!> author: Jiaxing Qi
!! Routines and parameter definitions for the standard D2Q9 model
!!
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
module mus_d2q9_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,              only: rk
  use tem_varSys_module,       only: tem_varSys_type, tem_varSys_op_type
  use tem_param_module,        only: cs2inv, div1_12, cs4inv, cs6inv, cs8inv,  &
    &                                div1_4, div1_9, div1_6, div1_18, div1_36, &
    &                                div1_2, div4_9, div2_9, div1_3, div2_3,   &
    &                                cs2, rho0
  use tem_aux_module,          only: tem_abort
  use tem_property_module,     only: prp_fluid, prp_fineGhostClosestToFluid
  use tem_construction_module, only: tem_levelDesc_type

  ! include musubi modules
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_param_module,              only: mus_param_type
  use mus_varSys_module,             only: mus_varSys_data_type
  use mus_directions_module,         only: qN0, q0N, q10, q01, qNN, qN1, q1N, &
    &                                      q11, q__W, q__S, q__E, q__N,       &
    &                                      q_SW, q_NW, q_SE, q_NE
  use mus_gradData_module,           only: mus_gradData_type
  use mus_derivedQuantities_module2, only: secondMom_2D
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_hrrInit_module,            only: HRR_Correction_d2q9, &
    &                                      getHermitepolynomials

  implicit none

  private

  public :: mus_advRel_kFluid_rBGK_vImproved_lD2Q9
  public :: mus_advRel_kFluid_rBGK_vStd_lD2Q9
  public :: mus_advRel_kFluid_rMRT_vStd_lD2Q9
  public :: mus_advRel_kFluidIncomp_rMRT_vStd_lD2Q9
  public :: mus_advRel_kFluidIncomp_rBGK_vStd_lD2Q9
  public :: bgk_Regularized_d2q9
  public :: bgk_RecursiveRegularized_d2q9
  public :: bgk_HybridRecursiveRegularized_d2q9
  public :: bgk_HybridRecursiveRegularizedCorr_d2q9
  public :: bgk_ProjectedRecursiveRegularized_d2q9
  public :: bgk_DualRelaxationTime_RR_d2q9
  public :: mus_advRel_kFluidGNS_rBGK_vStd_lD2Q9

  ! ============================================================================
  ! D2Q9 flow model
  ! ============================================================================
  !> Definition of the discrete velocity set
  integer, parameter :: QQ = 9   !< number of pdf directions
  ! General direction index
  integer, parameter :: q00 = 9  !< rest
  integer, parameter :: q__0 = 9 !< rest density is last

contains

! **************************************************************************** !
  !> Improved BGK model (with Galilean correction term)
  !! taken from Martin Geier cumulent paper 2015
  !! Geier, M., Schönherr, M., Pasquali, A., & Krafczyk, M. (2015).
  !! The cumulant lattice Boltzmann equation in three dimensions : Theory and
  !! validation. Computers and Mathematics with Applications.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type;compute]] function pointer.
  subroutine mus_advRel_kFluid_rBGK_vImproved_lD2Q9( fieldProp, inState,       &
    &                                                outState, auxField, neigh,&
    &                                                nElems, nSolve, level,    &
    &                                                layout, params, varSys,   &
    &                                                derVarPos )
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
    integer       :: iElem
    real(kind=rk) :: f(-1:1,-1:1)
    real(kind=rk) :: u,v,u2,v2
    real(kind=rk) :: rho, rho_omg
    real(kind=rk) :: inv_rho
    real(kind=rk) :: omega, cmpl_o, fac
    real(kind=rk) :: m20, m02
    real(kind=rk) :: sumX1, sumXN, sumY1, sumYN
    real(kind=rk) :: Gx, Gy
    real(kind=rk) :: x0, x1, xn, y0, y1, yn
    integer :: dens_pos, vel_pos(2), elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve
      f(-1, 0) = inState(  neigh(( qn0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f( 0,-1) = inState(  neigh(( q0n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f( 1, 0) = inState(  neigh(( q10-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f( 0, 1) = inState(  neigh(( q01-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f(-1,-1) = inState(  neigh(( qnn-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f(-1, 1) = inState(  neigh(( qn1-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f( 1,-1) = inState(  neigh(( q1n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f( 1, 1) = inState(  neigh(( q11-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f( 0, 0) = inState(  neigh(( q00-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u = auxField(elemOff + vel_pos(1))
      v = auxField(elemOff + vel_pos(2))

      ! u v square
      u2 = u*u
      v2 = v*v

      sumX1 = sum(f( 1,:))
      sumXN = sum(f(-1,:))

      sumY1 = sum(f(:, 1))
      sumYN = sum(f(:,-1))

      inv_rho = 1.0_rk / rho
      ! second moments, by equation A.7 and A.8
      m20 = ( sumX1 + sumXN ) * inv_rho
      m02 = ( sumY1 + sumYN ) * inv_rho

      ! relaxation rate
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1._rk - omega
      fac    = 4.5_rk - 2.25_rk * omega

      ! calculate Galilean correction term, by equation A.13 and A.14
      Gx = fac * u2 * ( m20 - div1_3 - u2 )
      Gy = fac * v2 * ( m02 - div1_3 - v2 )

      ! X Y components of eq
      ! by equation A.19 - A.21
      X0 = -div2_3 + u2 + Gx
      X1 = - ( X0 + 1.0_rk + u ) * 0.5_rk
      XN = X1 + u

      Y0 = -div2_3 + v2 + Gy
      Y1 = - ( Y0 + 1.0_rk + v ) * 0.5_rk
      YN = Y1 + v

      rho_omg = rho * omega
      X0 = X0 * rho_omg
      X1 = X1 * rho_omg
      XN = XN * rho_omg

      outState( ( ielem-1)* qq+ q00+( 1-1)* qq) &
        & = cmpl_o * f(0,0) + X0 * Y0
      outState( ( ielem-1)* qq+ qn0+( 1-1)* qq) &
        & = cmpl_o * f(-1,0) + XN * Y0
      outState( ( ielem-1)* qq+ q10+( 1-1)* qq) &
        & = cmpl_o * f(1,0) + X1 * Y0
      outState( ( ielem-1)* qq+ q0n+( 1-1)* qq) &
        & = cmpl_o * f(0,-1) + X0 * YN
      outState( ( ielem-1)* qq+ q01+( 1-1)* qq) &
        & = cmpl_o * f(0,1) + X0 * Y1
      outState( ( ielem-1)* qq+ qnn+( 1-1)* qq) &
        & = cmpl_o * f(-1,-1) + XN * YN
      outState( ( ielem-1)* qq+ q11+( 1-1)* qq) &
        & = cmpl_o * f(1,1) + X1 * Y1
      outState( ( ielem-1)* qq+ q1n+( 1-1)* qq) &
        & = cmpl_o * f(1,-1) + X1 * YN
      outState( ( ielem-1)* qq+ qn1+( 1-1)* qq) &
        & = cmpl_o * f(-1,1) + XN * Y1

    end do nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluid_rBGK_vImproved_lD2Q9
! **************************************************************************** !

! **************************************************************************** !
  !> No comment yet!
  !!
  !! TODO add comment
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rBGK_vStd_lD2Q9( fieldProp, inState, outState, &
    &                                           auxField, neigh, nElems,      &
    &                                           nSolve, level, layout, params,&
    &                                           varSys, derVarPos )
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
    integer :: iElem
    ! temporary distribution variables
    real(kind=rk) :: f1, f2, f3, f4, f5, f6, f7, f8, f9
    real(kind=rk) :: fEq1, fEq2, fEq3, fEq4, fEq5, fEq6, fEq7, fEq8, fEq9
    real(kind=rk) :: u_x, u_y
    real(kind=rk) :: rho, usq, ucx
    real(kind=rk) :: omega, cmpl_o
    real(kind=rk) :: c0, c1, c2, c3
    integer :: dens_pos, vel_pos(2), elemOff
    !type(mus_varSys_data_type), pointer :: fPtr
    !type( tem_levelDesc_type ), pointer :: levelDesc
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

    !call c_f_pointer( varSys%method%val( 1 )%method_data, fPtr )
    !levelDesc => fPtr%solverdata%scheme%levelDesc(level)

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve

      f1 = inState(  neigh(( 1-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f2 = inState(  neigh(( 2-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f3 = inState(  neigh(( 3-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f4 = inState(  neigh(( 4-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f5 = inState(  neigh(( 5-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f6 = inState(  neigh(( 6-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f7 = inState(  neigh(( 7-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f8 = inState(  neigh(( 8-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f9 = inState(  neigh(( 9-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)

      !! no collision for coarse ghosts and for fine ghost cells far from fluid cell
      !if ( .not. btest(levelDesc%property( iElem ), prp_fluid) ) then
      !  if ( .not. btest(levelDesc%property( iElem ), prp_fineGhostClosestToFluid) ) then
      !    outState( ( ielem-1)* 9+ 9+( 1-1)* 9) = f9
      !    outState( ( ielem-1)* 9+ 1+( 1-1)* 9) = f1
      !    outState( ( ielem-1)* 9+ 2+( 1-1)* 9) = f2
      !    outState( ( ielem-1)* 9+ 3+( 1-1)* 9) = f3
      !    outState( ( ielem-1)* 9+ 4+( 1-1)* 9) = f4
      !    outState( ( ielem-1)* 9+ 5+( 1-1)* 9) = f5
      !    outState( ( ielem-1)* 9+ 6+( 1-1)* 9) = f6
      !    outState( ( ielem-1)* 9+ 7+( 1-1)* 9) = f7
      !    outState( ( ielem-1)* 9+ 8+( 1-1)* 9) = f8
      !  CYCLE
      !  ! fine ghost cells closest to fluid cells should also skip the second collision
      !  else if ( params%iNesting( level ) == 2 ) then
      !    outState( ( ielem-1)* 9+ 9+( 1-1)* 9) = f9
      !    outState( ( ielem-1)* 9+ 1+( 1-1)* 9) = f1
      !    outState( ( ielem-1)* 9+ 2+( 1-1)* 9) = f2
      !    outState( ( ielem-1)* 9+ 3+( 1-1)* 9) = f3
      !    outState( ( ielem-1)* 9+ 4+( 1-1)* 9) = f4
      !    outState( ( ielem-1)* 9+ 5+( 1-1)* 9) = f5
      !    outState( ( ielem-1)* 9+ 6+( 1-1)* 9) = f6
      !    outState( ( ielem-1)* 9+ 7+( 1-1)* 9) = f7
      !    outState( ( ielem-1)* 9+ 8+( 1-1)* 9) = f8
      !    CYCLE
      !  endif
      !end if

      ! element offset for auxField array
      elemOff = (iElem - 1) * varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))

      ! calculate fEq
      usq = u_x * u_x + u_y * u_y
      c1 = rho - rho * usq * div1_2 * cs2inv
      feq9 = div4_9 * c1

      c0 = rho * cs2inv * cs2inv * div1_2
      c2 = rho * cs2inv * u_x
      c3 = rho * cs2inv * u_y

      feq1 = div1_9 * ( c1 - c2 + u_x * u_x * c0 )
      feq3 = feq1 + div2_9 * c2

      feq2 = div1_9 * ( c1 - c3 + u_y * u_y * c0 )
      feq4 = feq2 + div2_9 * c3

      ucx  = u_x + u_y
      feq5 = div1_36 * ( c1 - c2 - c3 + ucx * ucx * c0 )
      feq8 = feq5 + div1_18 * ( c2 + c3 )

      ucx  = u_x - u_y
      feq6 = div1_36 * ( c1 - c2 + c3 + ucx * ucx * c0 )
      feq7 = feq6 + div1_18 * ( c2 - c3 )


      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1._rk - omega

      outState( ( ielem-1)* 9+ 9+( 1-1)* 9) &
        & = cmpl_o * f9 + omega * fEq9
      outState( ( ielem-1)* 9+ 1+( 1-1)* 9) &
        & = cmpl_o * f1 + omega * fEq1
      outState( ( ielem-1)* 9+ 2+( 1-1)* 9) &
        & = cmpl_o * f2 + omega * fEq2
      outState( ( ielem-1)* 9+ 3+( 1-1)* 9) &
        & = cmpl_o * f3 + omega * fEq3
      outState( ( ielem-1)* 9+ 4+( 1-1)* 9) &
        & = cmpl_o * f4 + omega * fEq4
      outState( ( ielem-1)* 9+ 5+( 1-1)* 9) &
        & = cmpl_o * f5 + omega * fEq5
      outState( ( ielem-1)* 9+ 6+( 1-1)* 9) &
        & = cmpl_o * f6 + omega * fEq6
      outState( ( ielem-1)* 9+ 7+( 1-1)* 9) &
        & = cmpl_o * f7 + omega * fEq7
      outState( ( ielem-1)* 9+ 8+( 1-1)* 9) &
        & = cmpl_o * f8 + omega * fEq8

    end do nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluid_rBGK_vStd_lD2Q9
! **************************************************************************** !

! **************************************************************************** !
  !> Advection relaxation routine for the D2Q9 MRT model
  !! f( x+c, t+1 ) = f(x,t) - M^(-1)S( m - meq )
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rMRT_vStd_lD2Q9( fieldProp, inState, outState, &
    &                                           auxField, neigh, nElems,      &
    &                                           nSolve, level, layout, params,&
    &                                           varSys, derVarPos )
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
    integer       :: iElem
    real(kind=rk) :: f1, f2, f3, f4, f5, f6, f7, f8, f9
    real(kind=rk) :: fneq1, fneq2, fneq3, fneq4, fneq5, fneq6, fneq7, fneq8, &
      &              fneq9
    real(kind=rk) :: mom2, mom3, mom5, mom7, mom8, mom9
    real(kind=rk) :: meq2, meq3, meq8, meq9
    real(kind=rk) :: mneq2, mneq3, mneq5, mneq7, mneq8, mneq9
    real(kind=rk) :: s_mrt(9), omegaBulk
    real(kind=rk) :: rho ! local density
    real(kind=rk) :: ux, uy, u2, v2
    integer :: dens_pos, vel_pos(2), elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

    ! overwrite omegaKine indices 8 and 9 in element loop
    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)
    s_mrt = fieldProp(1)%fluid%mrtPtr( omegaKine = 1.0_rk, &
      &               omegaBulk = omegaBulk, &
      &               QQ        = 9       )

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve
      ! First load all local values into temp array
      f1 = inState(  neigh(( 1-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f2 = inState(  neigh(( 2-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f3 = inState(  neigh(( 3-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f4 = inState(  neigh(( 4-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f5 = inState(  neigh(( 5-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f6 = inState(  neigh(( 6-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f7 = inState(  neigh(( 7-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f8 = inState(  neigh(( 8-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f9 = inState(  neigh(( 9-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)

      ! element offset for auxField array
      elemOff = (iElem - 1) * varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      ux = auxField(elemOff + vel_pos(1))
      uy = auxField(elemOff + vel_pos(2))

      ! square of velocity
      u2 = ux * ux
      v2 = uy * uy
      ! meq1 =  rho ! rho
      meq2 =  -2._rk * rho + 3._rk * rho * (u2 + v2) ! e
      meq3 =  rho - 3._rk * rho * (u2 + v2) ! eps
      ! meq4 =  ux ! jx
      ! meq5 = -ux ! qx
      ! meq6 =  uy ! jy
      ! meq7 = -uy ! qy
      meq8 =  rho * (u2 - v2) ! pxx
      meq9 =  rho * ux * uy ! pxy

      ! convert pdf into moment
      mom2 = - f1 - f2 - f3 - f4 + 2.0_rk * ( f5 + f6 + f7 + f8 ) - 4.0_rk * f9
      mom3 = -2.0_rk * (f1 + f2 + f3 + f4) + f5 + f6 + f7 + f8 + 4.0_rk * f9
      mom5 = 2.0_rk * ( f1 - f3 ) - f5 - f6 + f7 + f8
      mom7 = 2.0_rk * ( f2 - f4 ) - f5 + f6 - f7 + f8
      mom8 = f1 - f2 + f3 - f4
      mom9 = f5 - f6 - f7 + f8

      ! omega
      s_mrt(8) = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      s_mrt(9) = s_mrt(8)

      ! compute neq moment
      ! mneq1 = s_mrt1 * ( mom1 - meq1 ) = 0
      mneq2 = s_mrt(2) * ( mom2 - meq2 )
      mneq3 = s_mrt(3) * ( mom3 - meq3 )
      ! mneq4 = s_mrt4 * ( mom4 - meq4 ) = 0
      mneq5 = s_mrt(5) * ( mom5 + rho*ux ) ! meq5 = -ux
      ! mneq6 = s_mrt6 * ( mom6 - meq6 ) = 0
      mneq7 = s_mrt(7) * ( mom7 + rho*uy ) ! meq7 = -uy
      mneq8 = s_mrt(8) * ( mom8 - meq8 )
      mneq9 = s_mrt(9) * ( mom9 - meq9 )

      ! compute fNeq
      ! do iDir = 1, 9
      !   fneq(iDir) = sum( MMIvD2Q9(iDir,:) * mneq(:) )
      ! end do
      ! mneq1 = mneq4 = mneq6 = 0
      fNeq1 = -div1_36 * mneq2 - div1_18 * mneq3 + div1_6 * mneq5 &
        &       + div1_4 * mneq8
      fNeq2 = -div1_36 * mneq2 - div1_18 * mneq3 + div1_6 * mneq7 &
        &       - div1_4 * mneq8
      fNeq3 = -div1_36 * mneq2 - div1_18 * mneq3 - div1_6 * mneq5 &
        &       + div1_4 * mneq8
      fNeq4 = -div1_36 * mneq2 - div1_18 * mneq3 - div1_6 * mneq7 &
        &       - div1_4 * mneq8
      fNeq5 = div1_18 * mneq2 + div1_36 * mneq3 - div1_12 * (mneq5 + mneq7) &
        &       + div1_4 * mneq9
      fNeq6 = div1_18 * mneq2 + div1_36 * mneq3 - div1_12 * (mneq5 - mneq7) &
        &       - div1_4 * mneq9
      fNeq7 = div1_18 * mneq2 + div1_36 * mneq3 + div1_12 * (mneq5 - mneq7) &
        &       - div1_4 * mneq9
      fNeq8 = div1_18 * mneq2 + div1_36 * mneq3 + div1_12 * (mneq5 + mneq7) &
        &       + div1_4 * mneq9
      fNeq9 = div1_9 * (-mneq2 + mneq3)

      outState( ( ielem-1)* 9+ 1+( 1-1)* 9) = f1 - fNeq1
      outState( ( ielem-1)* 9+ 2+( 1-1)* 9) = f2 - fNeq2
      outState( ( ielem-1)* 9+ 3+( 1-1)* 9) = f3 - fNeq3
      outState( ( ielem-1)* 9+ 4+( 1-1)* 9) = f4 - fNeq4
      outState( ( ielem-1)* 9+ 5+( 1-1)* 9) = f5 - fNeq5
      outState( ( ielem-1)* 9+ 6+( 1-1)* 9) = f6 - fNeq6
      outState( ( ielem-1)* 9+ 7+( 1-1)* 9) = f7 - fNeq7
      outState( ( ielem-1)* 9+ 8+( 1-1)* 9) = f8 - fNeq8
      outState( ( ielem-1)* 9+ 9+( 1-1)* 9) = f9 - fNeq9

    enddo nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluid_rMRT_vStd_lD2Q9
! **************************************************************************** !

! **************************************************************************** !
  !> Advection relaxation routine for the D2Q9 MRT model
  !! f( x+c, t+1 ) = f(x,t) - M^(-1)S( m - meq )
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluidIncomp_rMRT_vStd_lD2Q9( fieldProp, inState,      &
    &                                                 outState, auxField,      &
    &                                                 neigh, nElems, nSolve,   &
    &                                                 level, layout, params,   &
    &                                                 varSys, derVarPos )
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
    integer       :: iElem
    real(kind=rk) :: f1, f2, f3, f4, f5, f6, f7, f8, f9
    real(kind=rk) :: fneq1, fneq2, fneq3, fneq4, fneq5, fneq6, fneq7, fneq8, &
      &              fneq9
    real(kind=rk) :: mom2, mom3, mom5, mom7, mom8, mom9
    real(kind=rk) :: meq2, meq3, meq8, meq9
    real(kind=rk) :: mneq2, mneq3, mneq5, mneq7, mneq8, mneq9
    real(kind=rk) :: s_mrt(9), omegaBulk
    real(kind=rk) :: rho ! local density
    real(kind=rk) :: ux, uy, u2, v2
    integer :: dens_pos, vel_pos(2), elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

    ! overwrite omegaKine indices 8 and 9 in element loop
    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)
    s_mrt = fieldProp(1)%fluid%mrtPtr( omegaKine = 1.0_rk, &
      &               omegaBulk = omegaBulk, &
      &               QQ        = 9       )

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve
      ! First load all local values into temp array
      f1 = inState(  neigh(( 1-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f2 = inState(  neigh(( 2-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f3 = inState(  neigh(( 3-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f4 = inState(  neigh(( 4-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f5 = inState(  neigh(( 5-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f6 = inState(  neigh(( 6-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f7 = inState(  neigh(( 7-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f8 = inState(  neigh(( 8-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f9 = inState(  neigh(( 9-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)

      ! element offset for auxField array
      elemOff = (iElem - 1) * varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      ux = auxField(elemOff + vel_pos(1))
      uy = auxField(elemOff + vel_pos(2))

      u2 = ux * ux
      v2 = uy * uy
      ! meq1 =  rho ! rho
      meq2 =  -2._rk * rho + 3._rk * rho0 * (u2 + v2) ! e
      meq3 =  rho - 3._rk * rho0 * (u2 + v2) ! eps
      ! meq4 =  ux ! jx
      ! meq5 = -ux ! qx
      ! meq6 =  uy ! jy
      ! meq7 = -uy ! qy
      meq8 =  rho0 * (u2 - v2) ! pxx
      meq9 =  rho0 * ux * uy ! pxy

      ! convert pdf into moment
      mom2 = - f1 - f2 - f3 - f4 + 2.0_rk * (f5 + f6 + f7 + f8) - 4.0_rk * f9
      mom3 = -2.0_rk * (f1 + f2 + f3 + f4) + f5 + f6 + f7 + f8 + 4.0_rk * f9
      mom5 = 2.0_rk * (f1 - f3) - f5 - f6 + f7 + f8
      mom7 = 2.0_rk * (f2 - f4) - f5 + f6 - f7 + f8
      mom8 = f1 - f2 + f3 - f4
      mom9 = f5 - f6 - f7 + f8

      ! omega
      s_mrt(8) = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      s_mrt(9) = s_mrt(8)

      ! compute neq moment
      ! mneq1 = s_mrt1 * (mom1 - meq1) = 0
      mneq2 = s_mrt(2) * (mom2 - meq2)
      mneq3 = s_mrt(3) * (mom3 - meq3)
      ! mneq4 = s_mrt4 * (mom4 - meq4) = 0
      mneq5 = s_mrt(5) * (mom5 + rho0 * ux) ! meq5 = -ux
      ! mneq6 = s_mrt6 * (mom6 - meq6) = 0
      mneq7 = s_mrt(7) * (mom7 + rho0 * uy) ! meq7 = -uy
      mneq8 = s_mrt(8) * (mom8 - meq8)
      mneq9 = s_mrt(9) * (mom9 - meq9)

      ! compute fNeq
      ! do iDir = 1, 9
      !   fneq(iDir) = sum( MMIvD2Q9(iDir,:) * mneq(:) )
      ! end do
      ! mneq1 = mneq4 = mneq6 = 0
      fNeq1 = -div1_36 * mneq2 - div1_18 * mneq3 + div1_6 * mneq5 &
        &       + div1_4 * mneq8
      fNeq2 = -div1_36 * mneq2 - div1_18 * mneq3 + div1_6 * mneq7 &
        &       - div1_4 * mneq8
      fNeq3 = -div1_36 * mneq2 - div1_18 * mneq3 - div1_6 * mneq5 &
        &       + div1_4 * mneq8
      fNeq4 = -div1_36 * mneq2 - div1_18 * mneq3 - div1_6 * mneq7 &
        &       - div1_4 * mneq8
      fNeq5 = div1_18 * mneq2 + div1_36 * mneq3 - div1_12 * (mneq5 + mneq7) &
        &       + div1_4 * mneq9
      fNeq6 = div1_18 * mneq2 + div1_36 * mneq3 - div1_12 * (mneq5 - mneq7) &
        &       - div1_4 * mneq9
      fNeq7 = div1_18 * mneq2 + div1_36 * mneq3 + div1_12 * (mneq5 - mneq7) &
        &       - div1_4 * mneq9
      fNeq8 = div1_18 * mneq2 + div1_36 * mneq3 + div1_12 * (mneq5 + mneq7) &
        &       + div1_4 * mneq9
      fNeq9 = div1_9 * (-mneq2 + mneq3)

      outState( ( ielem-1)* 9+ 1+( 1-1)* 9) = f1 - fNeq1
      outState( ( ielem-1)* 9+ 2+( 1-1)* 9) = f2 - fNeq2
      outState( ( ielem-1)* 9+ 3+( 1-1)* 9) = f3 - fNeq3
      outState( ( ielem-1)* 9+ 4+( 1-1)* 9) = f4 - fNeq4
      outState( ( ielem-1)* 9+ 5+( 1-1)* 9) = f5 - fNeq5
      outState( ( ielem-1)* 9+ 6+( 1-1)* 9) = f6 - fNeq6
      outState( ( ielem-1)* 9+ 7+( 1-1)* 9) = f7 - fNeq7
      outState( ( ielem-1)* 9+ 8+( 1-1)* 9) = f8 - fNeq8
      outState( ( ielem-1)* 9+ 9+( 1-1)* 9) = f9 - fNeq9

    enddo nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluidIncomp_rMRT_vStd_lD2Q9
! **************************************************************************** !

! **************************************************************************** !
  !> No comment yet!
  !!
  !! TODO add comment
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluidIncomp_rBGK_vStd_lD2Q9( fieldProp, inState,      &
    &                                                 outState, auxField,      &
    &                                                 neigh, nElems, nSolve,   &
    &                                                 level, layout, params,   &
    &                                                 varSys, derVarPos )
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
    integer :: iElem
    ! temporary distribution variables
    real(kind=rk) :: f1, f2, f3, f4, f5, f6, f7, f8, f9
    real(kind=rk) :: fEq1, fEq2, fEq3, fEq4, fEq5, fEq6, fEq7, fEq8, fEq9
    real(kind=rk) :: u_x, u_y
    real(kind=rk) :: rho, usq, ucx
    real(kind=rk) :: omega, cmpl_o
    real(kind=rk) :: c0, c1, c2, c3
    integer :: dens_pos, vel_pos(2), elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve

      f1 = inState(  neigh(( 1-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f2 = inState(  neigh(( 2-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f3 = inState(  neigh(( 3-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f4 = inState(  neigh(( 4-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f5 = inState(  neigh(( 5-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f6 = inState(  neigh(( 6-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f7 = inState(  neigh(( 7-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f8 = inState(  neigh(( 8-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f9 = inState(  neigh(( 9-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)

      ! element offset for auxField array
      elemOff = (iElem - 1) * varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))

      ! calculate fEq
      usq = u_x * u_x + u_y * u_y
      c1 = rho - rho0 * usq * div1_2 * cs2inv
      feq9 = div4_9 * c1

      c0 = rho0 * cs2inv * cs2inv * div1_2
      c2 = rho0 * cs2inv * u_x
      c3 = rho0 * cs2inv * u_y

      feq1 = div1_9 * ( c1 - c2 + u_x*u_x*c0 )
      feq3 = feq1 + div2_9 * c2

      feq2 = div1_9 * ( c1 - c3 + u_y*u_y*c0 )
      feq4 = feq2 + div2_9 * c3

      ucx  = u_x+u_y
      feq5 = div1_36 * ( c1 - c2 - c3 + ucx*ucx*c0 )
      feq8 = feq5 + div1_18 * ( c2 + c3 )

      ucx  = u_x-u_y
      feq6 = div1_36 * ( c1 - c2 + c3 + ucx*ucx*c0 )
      feq7 = feq6 + div1_18 * ( c2 - c3 )


      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1._rk - omega

      outState(( ielem-1)* 9+9+( 1-1)* 9) &
        & = cmpl_o * f9 + omega * fEq9
      outState(( ielem-1)* 9+1+( 1-1)* 9) &
        & = cmpl_o * f1 + omega * fEq1
      outState(( ielem-1)* 9+2+( 1-1)* 9) &
        & = cmpl_o * f2 + omega * fEq2
      outState(( ielem-1)* 9+3+( 1-1)* 9) &
        & = cmpl_o * f3 + omega * fEq3
      outState(( ielem-1)* 9+4+( 1-1)* 9) &
        & = cmpl_o * f4 + omega * fEq4
      outState(( ielem-1)* 9+5+( 1-1)* 9) &
        & = cmpl_o * f5 + omega * fEq5
      outState(( ielem-1)* 9+6+( 1-1)* 9) &
        & = cmpl_o * f6 + omega * fEq6
      outState(( ielem-1)* 9+7+( 1-1)* 9) &
        & = cmpl_o * f7 + omega * fEq7
      outState(( ielem-1)* 9+8+( 1-1)* 9) &
        & = cmpl_o * f8 + omega * fEq8

    end do nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluidIncomp_rBGK_vStd_lD2Q9
! **************************************************************************** !

! ****************************************************************************** !
  ! > Regularizes f anf feq up to 2nd order
  !
  ! Based on:
  ! Latt, J. and Chopard, B. (2005) 'Lattice Boltzmann Method with regularized
  ! non-equilibrium distribution functions'
pure subroutine f_f_eq_regularized_2nd_ord_d2q9 ( weight, rho, u_x, u_y, feq, &
  &                                               f1, a12xx, a12yy, a12xy     )
    ! -------------------------------------------------------------------- !
    !> Weights of the stencil
    real(kind=rk), intent(in) :: weight(QQ)
    !> Density, velocity components
    real(kind=rk), intent(in) :: rho
    real(kind=rk), intent(in) :: u_x
    real(kind=rk), intent(in) :: u_y
    !> Equilibrium pdf and full pdf
    real(kind=rk), intent(out) :: feq(QQ)
    real(kind=rk), intent(out) :: f1(QQ)
    !> Coefficients of f1: a12xx, a12yy, a12xy
    real(kind=rk), intent(in) :: a12xx, a12yy, a12xy
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: u_x_sqr, u_y_sqr, u_x_u_y, f00, f01, f02, f12
    ! --------------------------------------------------------------------------

    u_x_sqr = u_x**2
    u_y_sqr = u_y**2
    u_x_u_y = u_x * u_y
    f00 = 1.0_rk

    !iDir = 1
    f01 = -cs2inv*u_x
    f02 = div1_6*cs4inv*(2._rk*u_x_sqr - u_y_sqr)
    feq(1) = weight(1) * rho * (f00 + f01 + f02)
    f12 = div1_6*cs4inv*(2._rk*a12xx - a12yy)
    f1(1) = weight(1) * f12

    !iDir = 3
    f01 = -f01
    feq(3) = weight(3) * rho * (f00 + f01 + f02)
    f1(3) = weight(3) * f12

    !iDir = 2
    f01 = -cs2inv*u_y
    f02 = div1_6*cs4inv*(-u_x_sqr + 2._rk*u_y_sqr)
    feq(2) = weight(2) * rho * (f00 + f01 + f02)
    f12 = div1_6*cs4inv*(-a12xx + 2._rk*a12yy)
    f1(2) = weight(2) * f12

    !iDir = 4
    f01 = -f01
    feq(4) = weight(4) * rho * (f00 + f01 + f02)
    f1(4) = weight(4) * f12

    !iDir = 5
    f01 = cs2inv*(-u_x - u_y)
    f02 = cs4inv*(div1_3*(u_x_sqr + u_y_sqr) + u_x_u_y)
    feq(5) = weight(5) * rho * (f00 + f01 + f02)
    f12 = cs4inv*(div1_3*(a12xx + a12yy) + a12xy)
    f1(5) = weight(5) * f12

    !iDir = 8
    f01 = -f01
    feq(8) = weight(8) * rho * (f00 + f01 + f02)
    f1(8) = weight(8) * f12

    !iDir = 6
    f01 = cs2inv*(-u_x + u_y)
    f02 = f02 - 2._rk*cs4inv*u_x_u_y
    feq(6) = weight(6) * rho * (f00 + f01 + f02)
    f12 = f12 - 2._rk*cs4inv*a12xy
    f1(6) = weight(6) * f12

    !iDir = 7
    f01 = -f01
    feq(7) = weight(7) * rho * (f00 + f01 + f02)
    f1(7) = weight(7) * f12

    !iDir = 9
    !f01 = 0._rk
    f02 = -div1_6*cs4inv*(u_x_sqr + u_y_sqr)
    feq(9) = weight(9) * rho * (f00 + f02)
    f12 = -div1_6*cs4inv*(a12xx + a12yy)
    f1(9) = weight(9) * f12

  end subroutine f_f_eq_regularized_2nd_ord_d2q9
! ****************************************************************************** !

! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
pure subroutine f_f_eq_regularized_4th_ord_d2q9 ( weight, rho, u_x, u_y, feq, &
  &  f1, a12xx, a12yy, a12xy )
    ! -------------------------------------------------------------------- !
    !> weights of the stencil
    real(kind=rk), intent(in) :: weight(QQ)
    !> density, velocity components
    real(kind=rk), intent(in) :: rho
    real(kind=rk), intent(in) :: u_x
    real(kind=rk), intent(in) :: u_y
    !> equilibrium pdf and full pdf
    real(kind=rk), intent(out) :: feq(QQ)
    real(kind=rk), intent(out) :: f1(QQ)
    !> coefficients of f1: a12xx, a12yy, a12xy
    real(kind=rk), intent(in) :: a12xx, a12yy, a12xy
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: u_x_sqr_u_y, u_y_sqr_u_x, u_x_sqr_u_y_sqr, f03, f04
    real(kind=rk) :: u_x_sqr, u_y_sqr, a13xyy, a13xxy, a14xxyy, f13, f14
    ! ---------------------------------------------------------------------------

      call f_f_eq_regularized_2nd_ord_d2q9 ( weight, rho, u_x, u_y, feq, f1, &
        &                                    a12xx, a12yy, a12xy )

      u_x_sqr = u_x**2
      u_y_sqr = u_y**2
      u_x_sqr_u_y = u_x_sqr * u_y
      u_y_sqr_u_x = u_y_sqr * u_x
      u_x_sqr_u_y_sqr = u_x_sqr * u_y_sqr

      a13xyy = 2.0_rk * u_y * a12xy + u_x * a12yy
      a13xxy = 2.0_rk * u_x * a12xy + u_y * a12xx
      a14xxyy = u_y*a13xxy + u_x_sqr*a12yy + 2.0_rk*u_x*u_y*a12xy

      !iDir = 1
      f03 = div1_6*cs6inv*(u_y_sqr_u_x)
      f04 = div1_18*cs8inv*(-u_x_sqr_u_y_sqr)
      feq(1) = feq(1) + weight(1) * rho * (f03 + f04)
      f13 = div1_6*cs6inv*(a13xyy)
      f14 = div1_18*cs8inv*(-a14xxyy)
      f1(1) = f1(1) + weight(1) * (f13 + f14)

      !iDir = 3
      f03 = -f03
      feq(3) = feq(3) + weight(3) * rho * (f03 + f04)
      f13 = -f13
      f1(3) = f1(3) + weight(3) * (f13 + f14)

      !iDir = 2
      f03 = div1_6*cs6inv*(u_x_sqr_u_y)
      feq(2) = feq(2) + weight(2) * rho * (f03 + f04)
      f13 = div1_6*cs6inv*(a13xxy)
      f1(2) = f1(2) + weight(2) * (f13 + f14)

      !iDir = 4
      f03 = -f03
      feq(4) = feq(4) + weight(4) * rho * (f03 + f04)
      f13 = -f13
      f1(4) = f1(4) + weight(4) * (f13 + f14)

      !iDir = 5
      f03 = div1_3*cs6inv*(-(u_x_sqr_u_y + u_y_sqr_u_x))
      f04 = -2._rk*f04
      feq(5) = feq(5) + weight(5) * rho * (f03 + f04)
      f13 = div1_3*cs6inv*(-(a13xxy + a13xyy))
      f14 = -2._rk*f14
      f1(5) = f1(5) + weight(5) * (f13 + f14)

      !iDir = 8
      f03 = -f03
      feq(8) = feq(8) + weight(8) * rho * (f03 + f04)
      f13 = -f13
      f1(8) = f1(8) + weight(8) * (f13 + f14)

      !iDir = 6
      f03 = div1_3*cs6inv*((u_x_sqr_u_y - u_y_sqr_u_x))
      feq(6) = feq(6) + weight(6) * rho * (f03 + f04)
      f13 = div1_3*cs6inv*((a13xxy - a13xyy))
      f1(6) = f1(6) + weight(6) * (f13 + f14)

      !iDir = 7
      f03 = -f03
      feq(7) = feq(7) + weight(7) * rho * (f03 + f04)
      f13 = -f13
      f1(7) = f1(7) + weight(7) * (f13 + f14)

      !iDir = 9
      !f03 = 0._rk
      f04 = div1_4*f04
      feq(9) = feq(9) + weight(9) * rho * (f04) !f03 +
      !f13 = 0._rk
      f14 = div1_4*f14
      f1(9) = f1(9) + weight(9) * (f14) !f13 +

  end subroutine f_f_eq_regularized_4th_ord_d2q9
! ****************************************************************************** !

! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
subroutine bgk_Regularized_d2q9 ( fieldProp, inState, outState, auxField, &
    &                                  neigh, nElems, nSolve, level, layout,   &
    &                                  params, varSys, derVarPos) !, gradData              )
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
    ! indeces
    integer :: iElem, iDir
    ! temporary distribution variables
    real(kind=rk) :: f( QQ ), SOM(3), SOM_neq(3)
    real(kind=rk) :: rho, u_x, u_y, a12xx, a12xy, a12yy
    real(kind=rk) :: omega, cmpl_o, feq(QQ), f1(QQ)
    integer :: denspos, velpos(3), elemOff, nScalars
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:2) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

    nScalars = varSys%nScalars

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve

      do iDir = 1, QQ
        f(iDir) = inState( neigh(( idir-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      enddo

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + denspos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + velpos(1))
      u_y = auxField(elemOff + velpos(2))

      ! non equilibrium second-order moments
      ! SOM_neq = SOM - SOM_eq
      SOM = secondMom_2D(layout%fStencil%cxcx, f, layout%fStencil%QQ)
      SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      SOM_neq(3) = SOM(3) - rho * u_x * u_y

      ! Hermitian coefficients
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1.0_rk - omega
      a12xx = SOM_neq(1)
      a12yy = SOM_neq(2)
      a12xy = SOM_neq(3)

      call f_f_eq_regularized_2nd_ord_d2q9( layout%weight(:), rho, u_x, u_y, feq, f1, &
        &                                   a12xx, a12yy, a12xy                       )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                   + cmpl_o * f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_Regularized_d2q9
! ****************************************************************************** !

! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
subroutine bgk_RecursiveRegularized_d2q9 ( fieldProp, inState, outState, auxField, &
    &                                  neigh, nElems, nSolve, level, layout,   &
    &                                  params, varSys, derVarPos) !, gradData              )
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
    ! indeces
    integer :: iElem, iDir
    ! temporary distribution variables
    real(kind=rk) :: f( QQ ), SOM(3), SOM_neq(3)
    real(kind=rk) :: rho, u_x, u_y, a12xx, a12xy, a12yy
    real(kind=rk) :: omega, cmpl_o, feq(QQ), f1(QQ)
    integer :: denspos, velpos(3), elemOff, nScalars
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:2) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

    nScalars = varSys%nScalars

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve

      do iDir = 1, QQ
        f(iDir) = inState( neigh(( idir-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      enddo

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + denspos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + velpos(1))
      u_y = auxField(elemOff + velpos(2))

      ! non equilibrium second-order moments
      ! SOM_neq = SOM - SOM_eq
      SOM = secondMom_2D(layout%fStencil%cxcx, f, layout%fStencil%QQ)
      SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      SOM_neq(3) = SOM(3) - rho * u_x * u_y

      ! Hermitian coefficients
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1.0_rk - omega
      a12xx = SOM_neq(1)
      a12yy = SOM_neq(2)
      a12xy = SOM_neq(3)

      call f_f_eq_regularized_4th_ord_d2q9( layout%weight(:), rho, u_x, u_y, feq, f1, &
        &                                   a12xx, a12yy, a12xy                       )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                   + cmpl_o * f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_RecursiveRegularized_d2q9
! ****************************************************************************** !


! ****************************************************************************** !
  ! based on High-order extension of the recursive regularized lattice
  ! Boltzmann method, PhD Thesis, COREIXAS 2018
  ! Correction term from but doesnt work!
  ! Solid wall and open boundary conditionsin hybrid recursive regularized
  ! latticeBoltzmann method for compressible flows, Feng 2019, PoF
subroutine bgk_ProjectedRecursiveRegularized_d2q9 ( fieldProp, inState, outState, auxField, &
    &                                  neigh, nElems, nSolve, level, layout,   &
    &                                  params, varSys, derVarPos) !, gradData              )
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
    !> gradient data
    type(mus_gradData_type), pointer :: gradData
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! indeces
    integer :: iElem, iDir
    ! temporary distribution variables
    real(kind=rk) :: f( QQ ), SR(3), gradU(2,2,1)!TODO::, tr_SR
    real(kind=rk) :: rho, u_x, u_y, a12xx, a12xy, a12yy
    real(kind=rk) :: omega, cmpl_o, feq(QQ), f1(QQ), taup
    integer :: denspos, velpos(3), elemOff, nScalars
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:2) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

    nScalars = varSys%nScalars

    ! access gradData
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( 1 )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme
    gradData => scheme%gradData(level)

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve

      do iDir = 1, QQ
        f(iDir) = inState( neigh(( idir-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      enddo

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + denspos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + velpos(1))
      u_y = auxField(elemOff + velpos(2))

      ! Stress tensor components
      gradU(:,:,1:1) = scheme%Grad%U_ptr(         &
           &   auxField     = auxField,           &
           &   gradData     = gradData,           &
           &   velPos       = velPos,             &
           &   nAuxScalars  = varSys%nAuxScalars, &
           &   nDims        = 2,                  &
           &   nSolve       = 1,                  &
           &   elemOffset   = iElem -1            )

      ! symmetric strain rate tensors
      ! transformed inro RHS of a1 FD equation
      ! the trace is needed only for energy conservation, which is not
      ! done in Musubi at the moment.
      !TODO::tr_SR = (gradU(1, 1, 1) + gradU(2, 2, 1)) ! * 2 / D = 1
      SR(1) = 2._rk * gradU(1, 1, 1)            !S_XX  !TODO:: - tr_SR
      SR(2) = 2._rk * gradU(2, 2, 1)            !S_YY  !TODO:: - tr_SR
      SR(3) = gradU(1, 2, 1) + gradU(2, 1, 1)   !S_XY

      ! Hermitian coefficients
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1.0_rk - omega
      taup = rho * cs2 / omega
      a12xx = -taup * SR(1)
      a12yy = -taup * SR(2)
      a12xy = -taup * SR(3)

      call f_f_eq_regularized_4th_ord_d2q9( layout%weight(:), rho, u_x, u_y, feq, f1, &
        &                                   a12xx, a12yy, a12xy                       )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                   + cmpl_o * f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_ProjectedRecursiveRegularized_d2q9
! **************************************************************************** !


! **************************************************************************** !
  !> Hybrid recursive regularization relaxation routine for the BGK model.
  !! based on: Feng et al., JCP 2019,
  !! "Hybrid recursive regularized thermal lattice Boltzmann model
  !! for high subsonic compressible flows"
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
subroutine bgk_HybridRecursiveRegularized_d2q9( fieldProp, inState, outState,  &
    &                                           auxField, neigh, nElems,       &
    &                                           nSolve, level, layout, params, &
    &                                           varSys, derVarPos              )
    ! -------------------------------------------------------------------- !
    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    !> Variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Vurrent layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> Number of elements in state Array
    integer, intent(in) :: nElems
    !> Input  pdf vector
    real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
    !> Output pdf vector
    real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
    !> Auxiliary field computed from pre-collision state
    !! is updated with correct velocity field for multicomponent models
    real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
    !> Connectivity vector
    integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
    !> Number of elements solved in kernel
    integer, intent(in) :: nSolve
    !> Current level
    integer,intent(in) :: level
    !> Global parameters
    type(mus_param_type),intent(in) :: params
    !> Position of derived quantities in varsys for all fields
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    !> Gradient data
    type(mus_gradData_type), pointer :: gradData
    type(mus_varSys_data_type), pointer :: fPtr
    ! Self-describing variables, loop indices, etc.
    integer :: iElem, iDir
    ! Temporary distribution variables
    real(kind=rk) :: pdfTmp(QQ)
    !symmetric strain rate tensor (SR) and it's trace (tr)
    real(kind=rk) :: SR(3)!TODO::, tr_SR
    ! Velocity gradient
    real(kind=rk) :: gradU(2,2,1)
    ! Second-order moments (SOM) and it's equilibrium
    real(kind=rk) :: SOM(3), SOM_neq(3)
    ! Density and velocity
    real(kind=rk) :: rho, u_x, u_y
    ! Hermitian coefficents
    real(kind=rk) :: a12xx, a12xy, a12yy
    ! Relaxation parameter omega and it's complementary
    real(kind=rk) :: omega, cmpl_o
    ! tau * pressure, f_eq and f_neq
    real(kind=rk) :: taup, feq(QQ), f1(QQ)
    ! sigma and it's complementary
    real(kind=rk) :: sigma, cmpl_sgm
    ! Self-describing variables
    integer :: denspos, velpos(3), elemOff, nScalars
    ! --------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    ! Get position of density and velocity in auxField to determine them later
    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:2) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

    nScalars = varSys%nScalars

    !call getHermitepolynomials( 2, QQ, layout, 4)

    ! Access data via variable method data which is a state variable
    call c_f_pointer( varSys%method%val( 1 )%method_data, fPtr )
    ! Set pointers for gradData via method data
    gradData => fPtr%solverData%scheme%gradData(level)

    ! Sigma value
    sigma = fieldProp(1)%fluid%HRR_sigma

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve

      ! Set pdf for all QQ-directions
      !> Generic fetching step:
      !! Streaming for pull
      !! Local copy for push
      do iDir = 1, QQ
        pdfTmp(iDir) = inState( neigh(( idir-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      enddo

      ! Calculate element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! Set local density
      rho = auxField(elemOff + denspos)
      ! Set local x- and y-velocity
      u_x = auxField(elemOff + velpos(1))
      u_y = auxField(elemOff + velpos(2))

      ! Calculate velocity gradient using central or forward difference, latter
      ! one if element has a boundary as neighbor
      gradU(:,:,1:1) = fPtr%solverData%scheme%Grad%U_ptr(           &
        &                        auxField     = auxField,           &
        &                        gradData     = gradData,           &
        &                        velPos       = velPos,             &
        &                        nAuxScalars  = varSys%nAuxScalars, &
        &                        nDims        = 2,                  &
        &                        nSolve       = 1,                  &
        &                        elemOffset   = iElem-1             )

      ! Calculate symmetric strain rate tensor (SR) and it's trace (tr)
      ! the trace is needed only for energy conservation, which is not
      ! done in Musubi at the moment
      !TODO::tr_SR = (gradU(1, 1, 1) + gradU(2, 2, 1)) ! * 2 / D = 1
      SR(1) = 2._rk * gradU(1, 1, 1)            !S_XX  !TODO:: - tr_SR
      SR(2) = 2._rk * gradU(2, 2, 1)            !S_YY  !TODO:: - tr_SR
      SR(3) = gradU(1, 2, 1) + gradU(2, 1, 1)   !S_XY

      ! Determine the non-equilibrium second-order moments via
      ! SOM_neq = SOM - SOM_eq
      ! First, get the second-order moments
      SOM = secondMom_2D(layout%fStencil%cxcx, pdfTmp, layout%fStencil%QQ)
      ! Second, subtract it's equilibrium: SOM_eq = rho*(cs²+u²)
      SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      SOM_neq(3) = SOM(3) - rho * u_x * u_y

      ! Hermitian coefficients
      ! Relaxation parameter
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      ! Complmentary of relaxation parameter (cmpl_o)
      cmpl_o = 1.0_rk - omega
      ! Complmentary of sigma (cmpl_sigma)
      cmpl_sgm = 1.0_rk - sigma
      ! tau * pressure
      taup = rho * cs2 / omega
      ! Hermite coefficient of second order
      a12xx = SOM_neq(1) * sigma + cmpl_sgm * (-taup * SR(1))
      a12yy = SOM_neq(2) * sigma + cmpl_sgm * (-taup * SR(2))
      a12xy = SOM_neq(3) * sigma + cmpl_sgm * (-taup * SR(3))

      ! regularization of f and feq up to 4th order
      call f_f_eq_regularized_4th_ord_d2q9( layout%weight(:), rho, u_x, u_y, &
        &                                   feq, f1, a12xx, a12yy, a12xy     )

      ! Relaxation
      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = &
          & feq(iDir) + cmpl_o * f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_HybridRecursiveRegularized_d2q9
! ****************************************************************************** !

! **************************************************************************** !
  !> Hybrid recursive regularization relaxation routine for the BGK model.
  !! based on: Feng et al., JCP 2019,
  !! "Hybrid recursive regularized thermal lattice Boltzmann model
  !! for high subsonic compressible flows"
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
subroutine bgk_HybridRecursiveRegularizedCorr_d2q9( fieldProp, inState, outState,  &
  &                                           auxField, neigh, nElems,       &
  &                                           nSolve, level, layout, params, &
  &                                           varSys, derVarPos              )
  ! -------------------------------------------------------------------- !
  !> Array of field properties (fluid or species)
  type(mus_field_prop_type), intent(in) :: fieldProp(:)
  !> Variable system definition
  type(tem_varSys_type), intent(in) :: varSys
  !> Vurrent layout
  type(mus_scheme_layout_type), intent(in) :: layout
  !> Number of elements in state Array
  integer, intent(in) :: nElems
  !> Input  pdf vector
  real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
  !> Output pdf vector
  real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
  !> Auxiliary field computed from pre-collision state
  !! is updated with correct velocity field for multicomponent models
  real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
  !> Connectivity vector
  integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
  !> Number of elements solved in kernel
  integer, intent(in) :: nSolve
  !> Current level
  integer,intent(in) :: level
  !> Global parameters
  type(mus_param_type),intent(in) :: params
  !> Position of derived quantities in varsys for all fields
  type(mus_derVarPos_type), intent(in) :: derVarPos(:)
  ! -------------------------------------------------------------------- !
  !> Gradient data
  type(mus_gradData_type), pointer :: gradData
  type(mus_varSys_data_type), pointer :: fPtr
  type(mus_scheme_type), pointer :: scheme
  ! Self-describing variables, loop indices, etc.
  integer :: iElem, iDir
  ! Temporary distribution variables
  real(kind=rk) :: pdfTmp(QQ)
  !symmetric strain rate tensor (SR) and it's trace (tr)
  real(kind=rk) :: SR(3)!TODO::, tr_SR
  ! Velocity gradient
  real(kind=rk) :: gradU(2,2,1), gradRHOU3(2,1)
  ! Second-order moments (SOM) and it's equilibrium
  real(kind=rk) :: SOM(3), SOM_neq(3)
  !real(kind=rk) :: gradRHOU3(2,1), Scorr
  ! Density and velocity
  real(kind=rk) :: rho, u_x, u_y
  ! Hermitian coefficents
  real(kind=rk) :: a12xx, a12xy, a12yy
  ! Relaxation parameter omega and it's complementary
  real(kind=rk) :: omega, cmpl_o
  ! tau * pressure, f_eq and f_neq, and Correction term
  real(kind=rk) :: taup, feq(QQ), f1(QQ), S_Corr(QQ)
  ! sigma and it's complementary
  real(kind=rk) :: sigma, cmpl_sgm
  ! Self-describing variables
  integer :: denspos, velpos(3), elemOff, nScalars, iSrc
  ! --------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

  ! Get position of density and velocity in auxField to determine them later
  denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
  velpos(1:2) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

  nScalars = varSys%nScalars

  ! Access data via variable method data which is a state variable
  call c_f_pointer( varSys%method%val( 1 )%method_data, fPtr )
  ! Set pointers for gradData via method data
  scheme => fPtr%solverData%scheme
  gradData => scheme%gradData(level)

  ! Sigma value, read from input
  ! standard value is 0.98
  sigma = fieldProp(1)%fluid%HRR_sigma

  ! allocate internalSource element array
  do iSrc = 1, scheme%field(1)%internalSource%varDict%nVals
    if ( trim(scheme%field(1)%internalSource%varDict%val(iSrc)%key) == 'hrr_correction' ) exit
  end do

  associate( HRR_Corr => scheme%field(1)%internalSource%method(iSrc)%elemLvl(Level)%HRR_Corr  )

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve

      ! Set pdf for all QQ-directions
      !> Generic fetching step:
      !! Streaming for pull
      !! Local copy for push
      do iDir = 1, QQ
        pdfTmp(iDir) = inState( neigh(( idir-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      enddo

      ! Calculate element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! Set local density
      rho = auxField(elemOff + denspos)
      ! Set local x- and y-velocity
      u_x = auxField(elemOff + velpos(1))
      u_y = auxField(elemOff + velpos(2))

      ! Calculate velocity gradient using central or forward difference, latter
      ! one if element has a boundary as neighbor
      gradU(:,:,1:1) = scheme%Grad%U_ptr( auxField     = auxField,           &
        &                                 gradData     = gradData,           &
        &                                 velPos       = velPos,             &
        &                                 nAuxScalars  = varSys%nAuxScalars, &
        &                                 nDims        = 2,                  &
        &                                 nSolve       = 1,                  &
        &                                 elemOffset   = iElem-1             )

      ! 1 = x, 2 = y, 3 = z, no xy returned
      gradRhoU3(:,1:1) = scheme%Grad%RhoU3_ptr( &
        &   auxField     = auxField,            &
        &   gradData     = gradData,            &
        &   velPos       = velpos,              &
        &   densPos      = denspos,             &
        &   nAuxScalars  = varSys%nAuxScalars,  &
        &   nDims        = 2,                   &
        &   nSolve       = 1,                   &
        &   elemOffset   = iElem-1              )

      ! Calculate correction
      call HRR_Correction_d2q9 (               &
        &    QQ        = QQ,                   &
        &    weight    = layout%weight(:),     &
        &    gradRHOU3 = gradRHOU3(:, 1),      &
        &    phi       = S_corr(:),            &
        &    dens      = HRR_Corr%dens(iElem), &
        &    vel       = HRR_Corr%vel(iElem,:) )

      ! Calculate symmetric strain rate tensor (SR) and it's trace (tr)
      ! the trace is needed only for energy conservation, which is not
      ! done in Musubi at the moment
      !TODO::tr_SR = (gradU(1, 1, 1) + gradU(2, 2, 1)) ! * 2 / D = 1
      SR(1) = 2._rk * gradU(1, 1, 1)            !S_XX  !TODO:: - tr_SR
      SR(2) = 2._rk * gradU(2, 2, 1)            !S_YY  !TODO:: - tr_SR
      SR(3) = gradU(1, 2, 1) + gradU(2, 1, 1)   !S_XY

      ! Determine the non-equilibrium second-order moments via
      ! SOM_neq = SOM - SOM_eq
      ! Apply correction
      pdfTmp(:) = pdfTmp(:) + 0.5_rk * S_corr(:)
      ! First, get the second-order moments
      SOM = secondMom_2D(layout%fStencil%cxcx, pdfTmp, layout%fStencil%QQ)
      ! Second, subtract it's equilibrium: SOM_eq = rho*(cs²+u²)
      SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      SOM_neq(3) = SOM(3) - rho * u_x * u_y

      ! Hermitian coefficients
      ! Relaxation parameter
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      ! Complmentary of relaxation parameter (cmpl_o)
      cmpl_o = 1.0_rk - omega
      ! Complmentary of sigma (cmpl_sigma)
      cmpl_sgm = 1.0_rk - sigma
      ! tau * pressure
      taup = rho * cs2 / omega
      ! Hermite coefficient of second order
      a12xx = SOM_neq(1) * sigma + cmpl_sgm * (-taup * SR(1))
      a12yy = SOM_neq(2) * sigma + cmpl_sgm * (-taup * SR(2))
      a12xy = SOM_neq(3) * sigma + cmpl_sgm * (-taup * SR(3))

      ! regularization of f and feq up to 4th order
      call f_f_eq_regularized_4th_ord_d2q9( layout%weight(:), rho, u_x, u_y, &
        &                                   feq, f1, a12xx, a12yy, a12xy     )

      ! Relaxation
      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = &
          & feq(iDir) + cmpl_o * f1(iDir) + 0.5_rk * S_corr(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end associate

end subroutine bgk_HybridRecursiveRegularizedCorr_d2q9
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine bgk_DualRelaxationTime_RR_d2q9 ( fieldProp, inState, outState, auxField, &
  &                                           neigh, nElems, nSolve, level, layout,   &
  &                                           params, varSys, derVarPos)
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
  ! indeces
  integer :: iElem, iDir
  ! temporary distribution variables
  real(kind=rk) :: f( QQ ), SOM(3), SOM_neq(3)
  real(kind=rk) :: rho, u_x, u_y, a12xx, a12xy, a12yy
  real(kind=rk) :: omega, tau, tauN, CoefTauNTau
  real(kind=rk) :: feq(QQ), f1(QQ), f_temp
  !real(kind=rk) :: sigma
  integer :: denspos, velpos(3), elemOff, nScalars
  ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:2) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)

    nScalars = varSys%nScalars

    tauN = fieldProp(1)%fluid%DRT_tauN

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve
      do iDir = 1, QQ
        f(iDir) = inState( neigh(( idir-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      enddo

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + denspos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + velpos(1))
      u_y = auxField(elemOff + velpos(2))

      ! non equilibrium second-order moments
      ! SOM_neq = SOM - SOM_eq
      SOM = secondMom_2D(layout%fStencil%cxcx, f, layout%fStencil%QQ)
      SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      SOM_neq(3) = SOM(3) - rho * u_x * u_y

      !Relaxation coefficients
      ! remains constant on uniform mesh. Does loop over elements include voxels off different size?
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      tau = 1.0_rk / omega
      CoefTauNTau = ( tau - tauN ) / ( tau * tauN )

      a12xx = SOM_neq(1)
      a12yy = SOM_neq(2)
      a12xy = SOM_neq(3)

      call f_f_eq_regularized_2nd_ord_d2q9( layout%weight(:), rho, u_x, u_y, feq, f1, &
        &                                   a12xx, a12yy, a12xy                       )

      do iDir = 1, QQ
        f_temp = f(iDir) - 1.0_rk/tauN * ( f(iDir) - feq(iDir) )
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = f_temp     &
          &                                                + CoefTauNTau * f1(iDir)
      enddo
    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_DualRelaxationTime_RR_d2q9
! ****************************************************************************** !

  !> BGK relaxation routine using equilibrium distribution for simulating the
  !! Generalized Navier Stokes (GNS) aka Volume-Averaged Navier-Stokes (VANS)
  !! equations for coupled LBM-DEM simulations
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluidGNS_rBGK_vStd_lD2Q9(     &
    &          fieldProp, inState, outState, auxField, &
    &          neigh, nElems, nSolve, level, layout,   &
    &          params, varSys, derVarPos               )
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
    integer :: iElem
    ! temporary distribution variables
    real(kind=rk) :: f1, f2, f3, f4, f5, f6, f7, f8, f9
    real(kind=rk) :: fEq1, fEq2, fEq3, fEq4, fEq5, fEq6, fEq7, fEq8, fEq9
    real(kind=rk) :: u_x, u_y
    real(kind=rk) :: rho, usq, ucx
    real(kind=rk) :: eps_inv     ! fluid volume fraction
    real(kind=rk) :: omega, cmpl_o
    real(kind=rk) :: c0, c1, c2, c3
    integer :: dens_pos, vel_pos(2), vol_frac_pos, elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)
    vol_frac_pos = varSys%method%val(derVarPos(1)%vol_frac)%auxField_varPos(1)

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve

      f1 = inState(  neigh(( 1-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f2 = inState(  neigh(( 2-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f3 = inState(  neigh(( 3-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f4 = inState(  neigh(( 4-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f5 = inState(  neigh(( 5-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f6 = inState(  neigh(( 6-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f7 = inState(  neigh(( 7-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f8 = inState(  neigh(( 8-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)
      f9 = inState(  neigh(( 9-1)* nelems+ ielem)+( 1-1)* 9+ 9*0)

      ! element offset for auxField array
      elemOff = (iElem - 1) * varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))

      ! Inverse of local fluid volume fraction
      eps_inv = 1.0_rk / auxField(elemOff + vol_frac_pos)

      ! calculate fEq
      usq = u_x * u_x + u_y * u_y
      c1 = rho - rho * usq * div1_2 * cs2inv * eps_inv
      feq9 = div4_9 * c1

      c0 = rho * cs2inv * cs2inv * div1_2 * eps_inv
      c2 = rho * cs2inv * u_x
      c3 = rho * cs2inv * u_y

      feq1 = div1_9 * ( c1 - c2 + u_x * u_x * c0 )
      feq3 = feq1 + div2_9 * c2

      feq2 = div1_9 * ( c1 - c3 + u_y * u_y * c0 )
      feq4 = feq2 + div2_9 * c3

      ucx  = u_x + u_y
      feq5 = div1_36 * ( c1 - c2 - c3 + ucx * ucx * c0 )
      feq8 = feq5 + div1_18 * ( c2 + c3 )

      ucx  = u_x - u_y
      feq6 = div1_36 * ( c1 - c2 + c3 + ucx * ucx * c0 )
      feq7 = feq6 + div1_18 * ( c2 - c3 )


      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1._rk - omega

      outState( ( ielem-1)* 9+ 9+( 1-1)* 9) &
        & = cmpl_o * f9 + omega * fEq9
      outState( ( ielem-1)* 9+ 1+( 1-1)* 9) &
        & = cmpl_o * f1 + omega * fEq1
      outState( ( ielem-1)* 9+ 2+( 1-1)* 9) &
        & = cmpl_o * f2 + omega * fEq2
      outState( ( ielem-1)* 9+ 3+( 1-1)* 9) &
        & = cmpl_o * f3 + omega * fEq3
      outState( ( ielem-1)* 9+ 4+( 1-1)* 9) &
        & = cmpl_o * f4 + omega * fEq4
      outState( ( ielem-1)* 9+ 5+( 1-1)* 9) &
        & = cmpl_o * f5 + omega * fEq5
      outState( ( ielem-1)* 9+ 6+( 1-1)* 9) &
        & = cmpl_o * f6 + omega * fEq6
      outState( ( ielem-1)* 9+ 7+( 1-1)* 9) &
        & = cmpl_o * f7 + omega * fEq7
      outState( ( ielem-1)* 9+ 8+( 1-1)* 9) &
        & = cmpl_o * f8 + omega * fEq8

    end do nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluidGNS_rBGK_vStd_lD2Q9
! **************************************************************************** !

end module mus_d2q9_module
! **************************************************************************** !
