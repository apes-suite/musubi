! Copyright (c) 2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
! **************************************************************************** !
!> author: Kannan Masilamani
!! This module provides the definition and methods for
!! MRT advection relaxation scheme for D3Q27 stencil.
!! The weighted MRT is based on the following paper
!! Abbas Fakhari, Diogo Bolster, Li-Shi Luo
!! "A weighted multiple-relaxation-time lattice Boltzmann method for multiphase
!! flows and its application to partial coalescence cascades"
!! Journal of Computational Physics, 2017
!!
!! Density (rho) and velocity (ux, uy, uz) are conserved during collision.
!!  i.e. m(1) = meq(1) --> mneq(1) = 0
!!       m(2) = meq(2) --> mneq(2) = 0
!!       m(3) = meq(3) --> mneq(3) = 0
!!       m(4) = meq(4) --> mneq(4) = 0
!!
!! Collision parameters are chosen as
!! s(1:4) = 0
!! s(5:9) = omega
!! s(10) = bulk_omega
!! s(11:27) = 1.0
module mus_mrt_d3q27_module
  use iso_c_binding,            only: c_f_pointer

  ! include treelm modules
  use env_module,               only: rk
  use tem_varSys_module,        only: tem_varSys_type, tem_varSys_op_type
  use tem_param_module,         only: rho0, cs2inv, cs4inv, t2cs2inv, t2cs4inv

  ! include musubi modules
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_param_module,         only: mus_param_type
  use mus_varSys_module,        only: mus_varSys_data_type
  use mus_derVarPos_module,     only: mus_derVarPos_type
  use mus_mrtInit_module,       only: WMMtrD3Q27, WMMIvD3Q27

  implicit none

  private

  public :: mus_advRel_kCFD_rMRT_vStdNoOpt_lD3Q27
  public :: mus_advRel_kFluid_rMRT_vStd_lD3Q27
  public :: mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q27

  !=============================================================================
  ! D3Q27 flow model
  !=============================================================================
  !> Definition of the discrete velocity set

  integer,parameter :: QQ   = 27   !< number of pdf directions
  integer,parameter :: q000 = 27    !< rest density is last

contains

! **************************************************************************** !
  !> Unoptimized explicit implementation
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kCFD_rMRT_vStdNoOpt_lD3Q27( fieldProp, inState,        &
    &                                               outState, auxField, neigh, &
    &                                               nElems, nSolve, level,     &
    &                                               layout, params, varSys,    &
    &                                               derVarPos )
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
    integer       :: iElem, iDir
    real(kind=rk) :: pdfTmp( QQ ) ! temporary local pdf values
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: vel(3)  ! local velocity
    real(kind=rk) :: fEq( QQ ) !< equilibrium distribution
    ! MRT Variables
    real(kind=rk) :: s_mrt( QQ )
    real(kind=rk) :: mneq( QQ )
    real(kind=rk) :: fneq( QQ )
    real(kind=rk) :: omegaBulk
    integer       :: dens_pos, vel_pos(3), elemOff
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)
    ! collision matrix
    s_mrt = fieldProp(1)%fluid%mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1,nSolve

      !> First load all local values into temp array
      do iDir = 1, QQ
        pdfTmp( iDir ) = inState( neigh (( idir-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      end do

      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      vel(1) = auxField(elemOff + vel_pos(1))
      vel(2) = auxField(elemOff + vel_pos(2))
      vel(3) = auxField(elemOff + vel_pos(3))

      ! Calculate the equilibrium distribution function
      fEq = layout%quantities%pdfEq_ptr(rho=rho, vel=vel, QQ=QQ)

      ! Calculate the non-equilibrium part
      fnEq(:) = pdfTmp(:) - fEq(:)

      ! convert non-equilibrium PDF into moments
      do iDir = 1, QQ
        mneq(iDir) =  sum( fnEq(:) * WMMtrD3Q27(iDir,:) )
      end do

      ! viscosity omega
      s_mrt( 5:9 ) = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)

      ! multiply it with collision matrix
      mneq(:) = s_mrt(:) * mneq(:)

      ! compute fNeq
      do iDir = 1, QQ
        fneq(iDir) = sum( WMMIvD3Q27(iDir,:) * mneq(:) )
        outState((ielem-1)*qq+idir+(1-1)*qq) = pdfTmp(iDir) - fneq(iDir)
      end do

    enddo nodeloop

  end subroutine mus_advRel_kCFD_rMRT_vStdNoOpt_lD3Q27
! **************************************************************************** !

! **************************************************************************** !
  !> Semi-optimized explicit implementation
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rMRT_vStd_lD3Q27( fieldProp, inState, outState, &
    &                                            auxField, neigh, nElems,      &
    &                                            nSolve, level, layout, params,&
    &                                            varSys, derVarPos )
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
    integer       :: iElem, iDir
    integer :: nScalars
    real(kind=rk) :: f(QQ)
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: u_x     ! local x-velocity
    real(kind=rk) :: u_y     ! local y-velocity
    real(kind=rk) :: u_z     ! local z-velocity
    real(kind=rk) :: sum_19_22, sum_23_26, sum_19_26, sum_21_22_23_24
    real(kind=rk) :: sum_20_21_24_25, sum_20_22_23_25, sum_7_10, sum_11_14
    real(kind=rk) :: sum_15_18, sum_11_18, sum_19_20_23_24, sum_19_21_23_25
    ! MRT Variables
    real(kind=rk) :: s_mrt( QQ ), meq( QQ )
    real(kind=rk) :: mneq( QQ ), mom( QQ )
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: omegaBulk
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    nScalars = varSys%nScalars

    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)
    ! overwrite omegaKine term in the element loop
    s_mrt = fieldProp(1)%fluid%mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1,nSolve

      !> First load all local values into temp array
      do iDir = 1, QQ
        f( iDir ) = inState( neigh (( idir-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      end do

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)

      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      ! Zero moment
      mom( 1) = rho

      ! First moments
      mom( 2) = rho*u_x
      mom( 3) = rho*u_y
      mom( 4) = rho*u_z

      sum_19_22 = f(19) + f(20) + f(21) + f(22)
      sum_23_26 = f(23) + f(24) + f(25) + f(26)
      sum_19_26 = sum_19_22 + sum_23_26
      sum_21_22_23_24 = f(21) + f(22) + f(23) + f(24)
      mom( 5) = f(15) - f(16) - f(17) + f(18) + sum_19_26 - 2._rk*sum_21_22_23_24
      sum_20_21_24_25 = f(20) + f(21) + f(24) + f(25)
      mom( 6) = f(7) - f(8) - f(9) + f(10) + sum_19_26 - 2._rk*sum_20_21_24_25
      sum_20_22_23_25 = f(20) + f(22) + f(23) + f(25)
      mom( 7) = f(11) - f(12) - f(13) + f(14) + sum_19_26 - 2._rk*sum_20_22_23_25

      sum_7_10 = f(7) + f(8) + f(9) + f(10)
      sum_11_14 = f(11) + f(12) + f(13) + f(14)
      sum_15_18 = f(15) + f(16) + f(17) + f(18)
      sum_11_18 = sum_11_14 + sum_15_18
      mom( 8) = 2.0_rk*( f(1) + f(4) - sum_7_10 ) - f(2) - f(3) - f(5) - f(6) &
        &       + sum_11_18
      mom( 9) = f(2) - f(3) + f(5) - f(6) - sum_11_14 + sum_15_18
      mom(10) = sum_7_10 + sum_11_18 + 2._rk*( sum_19_26 ) - f(27)

      mom(11) = 2.0_rk*( f(1) - f(4) ) - f(11) + f(12) - f(13) + f(14) - f(15) &
        &       - f(16) + f(17) + f(18) + 4._rk*( -sum_19_22 + sum_23_26 )
      sum_19_20_23_24 = f(19) + f(20) + f(23) + f(24)
      mom(12) = 2.0_rk*( f(2) - f(5) ) - f(7) - f(8) + f(9) + f(10) - f(15) + f(16) &
        &       - f(17) + f(18) + 4._rk*( sum_19_26 - 2._rk*sum_19_20_23_24 )
      sum_19_21_23_25 = f(19) + f(21) + f(23) + f(25)
      mom(13) = 2.0_rk*( f(3) - f(6) ) - f(7) + f(8) - f(9) + f(10) - f(11) - f(12) &
        &       + f(13) + f(14) + 4._rk*( sum_19_26 - 2._rk*sum_19_21_23_25 )

      mom(14) = f(11) - f(12) + f(13) - f(14) - f(15) - f(16) + f(17) + f(18)
      mom(15) = -f(7) - f(8) + f(9) + f(10) + f(15) - f(16) + f(17) - f(18)
      mom(16) = f(7) - f(8) + f(9) - f(10) - f(11) - f(12) + f(13) + f(14)
      mom(17) = -f(19) + f(20) + f(21) - f(22) + f(23) - f(24) - f(25) + f(26)

      mom(18) = -f(1) - f(2) - f(3) - f(4) - f(5) - f(6) &
        &       + 4._rk*( sum_19_26 ) + f(27)
      mom(19) = 2.0_rk*( -f(1) - f(4) ) + f(2) + f(3) + f(5) + f(6) - 4._rk*sum_7_10 &
        &       + 2._rk*( sum_11_18 )
      mom(20) = -f(2) + f(3) - f(5) + f(6) + 2._rk*( -sum_11_14 + sum_15_18 )

      mom(21) = -f(15) + f(16) + f(17) - f(18) + 2.0_rk*( sum_19_26 - 2._rk*sum_21_22_23_24 )
      mom(22) = -f(7) + f(8) + f(9) - f(10) + 2.0_rk*( sum_19_26 - 2._rk*sum_20_21_24_25 )
      mom(23) = -f(11) + f(12) + f(13) - f(14) + 2.0_rk*( sum_19_26 - 2._rk*sum_20_22_23_25 )

      ! third order pseudo-vector
      mom(24) = -f(1) + f(4) + 2._rk*( f(11) - f(12) + f(13) - f(14) + f(15) + f(16) - f(17) &
        &       -f(18) ) + 4._rk*( -sum_19_22 + sum_23_26 )
      mom(25) = -f(2) + f(5) + 2._rk*( f(7) + f(8) - f(9) - f(10) + f(15) - f(16) + f(17) &
        &       -f(18) ) + 4._rk*( sum_19_26 - 2._rk*sum_19_20_23_24 )
      mom(26) = -f(3) + f(6) + 2._rk*( f(7) - f(8) + f(9) - f(10) + f(11) + f(12) - f(13) &
        &       -f(14) ) + 4._rk*( sum_19_26 - 2._rk*sum_19_21_23_25 )

      mom(27) = 2._rk*( f(1) + f(2) + f(3) + f(4) + f(5) + f(6) ) + 4._rk*( -sum_7_10 &
        &       -sum_11_18 ) + 8._rk*( sum_19_26 ) - f(27)

      ! -------------------------------------------------------------------------
      ! Equilibrium moments
      meq(1:QQ) =  0.0_rk
      meq( 1) = rho
      meq( 2) = rho*u_x
      meq( 3) = rho*u_y
      meq( 4) = rho*u_z
      meq( 5) = meq(2) * u_y
      meq( 6) = meq(3) * u_z
      meq( 7) = meq(4) * u_x
      meq( 8) = rho * (2.0_rk*u_x*u_x - u_y*u_y - u_z*u_z)
      meq( 9) = rho * (u_y*u_y - u_z*u_z)
      meq(10) = rho * (u_x*u_x + u_y*u_y + u_z*u_z)

      ! update kinematic omega
      s_mrt(5:9) = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      ! compute neq moment
      mneq(:) = s_mrt(:) * ( mom(:) - meq(:) )

      ! compute fNeq
      do iDir = 1, QQ
        f(iDir) = f(iDir) - sum( WMMIvD3Q27(iDir,:) * mneq(:) )
      end do

      do iDir = 1, QQ
        outState(( ielem-1)* nscalars+ idir+( 1-1)* qq) = &
          & f( iDir )
      end do

    enddo nodeloop

  end subroutine mus_advRel_kFluid_rMRT_vStd_lD3Q27
! **************************************************************************** !

! **************************************************************************** !
  !> Semi-optimized explicit implementation
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q27( fieldProp, inState,     &
    &                                                  outState, auxField,     &
    &                                                  neigh, nElems, nSolve,  &
    &                                                  level, layout, params,  &
    &                                                  varSys, derVarPos )
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
    integer       :: iElem, iDir
    integer :: nScalars
    real(kind=rk) :: f(QQ)
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: u_x     ! local x-velocity
    real(kind=rk) :: u_y     ! local y-velocity
    real(kind=rk) :: u_z     ! local z-velocity
    real(kind=rk) :: sum_19_22, sum_23_26, sum_19_26, sum_21_22_23_24
    real(kind=rk) :: sum_20_21_24_25, sum_20_22_23_25, sum_7_10, sum_11_14
    real(kind=rk) :: sum_15_18, sum_11_18, sum_19_20_23_24, sum_19_21_23_25
    ! MRT Variables
    real(kind=rk) :: s_mrt( QQ ), meq( QQ )
    real(kind=rk) :: mneq( QQ ), mom( QQ )
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: omegaBulk
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    nScalars = varSys%nScalars

    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)
    ! overwrite omegaKine term in the element loop
    s_mrt = fieldProp(1)%fluid%mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1,nSolve

      !> First load all local values into temp array
      do iDir = 1, QQ
        f( iDir ) = inState( neigh (( idir-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      end do

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)

      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      ! Zero moment
      mom( 1) = rho

      ! First moments
      mom( 2) = rho*u_x
      mom( 3) = rho*u_y
      mom( 4) = rho*u_z

      sum_19_22 = f(19) + f(20) + f(21) + f(22)
      sum_23_26 = f(23) + f(24) + f(25) + f(26)
      sum_19_26 = sum_19_22 + sum_23_26
      sum_21_22_23_24 = f(21) + f(22) + f(23) + f(24)
      mom( 5) = f(15) - f(16) - f(17) + f(18) + sum_19_26 - 2._rk*sum_21_22_23_24
      sum_20_21_24_25 = f(20) + f(21) + f(24) + f(25)
      mom( 6) = f(7) - f(8) - f(9) + f(10) + sum_19_26 - 2._rk*sum_20_21_24_25
      sum_20_22_23_25 = f(20) + f(22) + f(23) + f(25)
      mom( 7) = f(11) - f(12) - f(13) + f(14) + sum_19_26 - 2._rk*sum_20_22_23_25

      sum_7_10 = f(7) + f(8) + f(9) + f(10)
      sum_11_14 = f(11) + f(12) + f(13) + f(14)
      sum_15_18 = f(15) + f(16) + f(17) + f(18)
      sum_11_18 = sum_11_14 + sum_15_18
      mom( 8) = 2.0_rk*( f(1) + f(4) - sum_7_10 ) - f(2) - f(3) - f(5) - f(6) &
        &       + sum_11_18
      mom( 9) = f(2) - f(3) + f(5) - f(6) - sum_11_14 + sum_15_18
      mom(10) = sum_7_10 + sum_11_18 + 2._rk*( sum_19_26 ) - f(27)

      mom(11) = 2.0_rk*( f(1) - f(4) ) - f(11) + f(12) - f(13) + f(14) - f(15) &
        &       - f(16) + f(17) + f(18) + 4._rk*( -sum_19_22 + sum_23_26 )
      sum_19_20_23_24 = f(19) + f(20) + f(23) + f(24)
      mom(12) = 2.0_rk*( f(2) - f(5) ) - f(7) - f(8) + f(9) + f(10) - f(15) + f(16) &
        &       - f(17) + f(18) + 4._rk*( sum_19_26 - 2._rk*sum_19_20_23_24 )
      sum_19_21_23_25 = f(19) + f(21) + f(23) + f(25)
      mom(13) = 2.0_rk*( f(3) - f(6) ) - f(7) + f(8) - f(9) + f(10) - f(11) - f(12) &
        &       + f(13) + f(14) + 4._rk*( sum_19_26 - 2._rk*sum_19_21_23_25 )

      mom(14) = f(11) - f(12) + f(13) - f(14) - f(15) - f(16) + f(17) + f(18)
      mom(15) = -f(7) - f(8) + f(9) + f(10) + f(15) - f(16) + f(17) - f(18)
      mom(16) = f(7) - f(8) + f(9) - f(10) - f(11) - f(12) + f(13) + f(14)
      mom(17) = -f(19) + f(20) + f(21) - f(22) + f(23) - f(24) - f(25) + f(26)

      mom(18) = -f(1) - f(2) - f(3) - f(4) - f(5) - f(6) &
        &       + 4._rk*( sum_19_26 ) + f(27)
      mom(19) = 2.0_rk*( -f(1) - f(4) ) + f(2) + f(3) + f(5) + f(6) - 4._rk*sum_7_10 &
        &       + 2._rk*( sum_11_18 )
      mom(20) = -f(2) + f(3) - f(5) + f(6) + 2._rk*( -sum_11_14 + sum_15_18 )

      mom(21) = -f(15) + f(16) + f(17) - f(18) + 2.0_rk*( sum_19_26 - 2._rk*sum_21_22_23_24 )
      mom(22) = -f(7) + f(8) + f(9) - f(10) + 2.0_rk*( sum_19_26 - 2._rk*sum_20_21_24_25 )
      mom(23) = -f(11) + f(12) + f(13) - f(14) + 2.0_rk*( sum_19_26 - 2._rk*sum_20_22_23_25 )

      ! third order pseudo-vector
      mom(24) = -f(1) + f(4) + 2._rk*( f(11) - f(12) + f(13) - f(14) + f(15) + f(16) - f(17) &
        &       -f(18) ) + 4._rk*( -sum_19_22 + sum_23_26 )
      mom(25) = -f(2) + f(5) + 2._rk*( f(7) + f(8) - f(9) - f(10) + f(15) - f(16) + f(17) &
        &       -f(18) ) + 4._rk*( sum_19_26 - 2._rk*sum_19_20_23_24 )
      mom(26) = -f(3) + f(6) + 2._rk*( f(7) - f(8) + f(9) - f(10) + f(11) + f(12) - f(13) &
        &       -f(14) ) + 4._rk*( sum_19_26 - 2._rk*sum_19_21_23_25 )

      mom(27) = 2._rk*( f(1) + f(2) + f(3) + f(4) + f(5) + f(6) ) + 4._rk*( -sum_7_10 &
        &       -sum_11_18 ) + 8._rk*( sum_19_26 ) - f(27)

      ! -------------------------------------------------------------------------
      ! Equilibrium moments
      meq(1:QQ) =  0.0_rk
      meq( 1) = rho
      meq( 2) = rho0*u_x
      meq( 3) = rho0*u_y
      meq( 4) = rho0*u_z
      meq( 5) = meq(2) * u_y
      meq( 6) = meq(3) * u_z
      meq( 7) = meq(4) * u_x
      meq( 8) = rho0 * (2.0_rk*u_x*u_x - u_y*u_y - u_z*u_z)
      meq( 9) = rho0 * (u_y*u_y - u_z*u_z)
      meq(10) = rho0 * (u_x*u_x + u_y*u_y + u_z*u_z)

      ! update kinematic omega
      s_mrt(5:9) = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      ! compute neq moment
      mneq(:) = s_mrt(:) * ( mom(:) - meq(:) )

      ! compute fNeq
      do iDir = 1, QQ
        f(iDir) = f(iDir) - sum( WMMIvD3Q27(iDir,:) * mneq(:) )
      end do

      do iDir = 1, QQ
        outState(( ielem-1)* nscalars+ idir+( 1-1)* qq) = &
          & f( iDir )
      end do

    enddo nodeloop

  end subroutine mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q27
! **************************************************************************** !

end module mus_mrt_d3q27_module
