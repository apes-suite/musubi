! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! ****************************************************************************** !
!> author: Jiaxing Qi
!! Routines and parameter definitions for the standard D3Q27 model
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
module mus_d3q27_module
  use iso_c_binding,            only: c_f_pointer

  ! include treelm modules
  use env_module,              only: rk
  use tem_varSys_module,       only: tem_varSys_type, tem_varSys_op_type
  use tem_param_module,        only: div1_2, div1_3, cs2inv, cs4inv, t2cs2inv,  &
    &                                rho0, cs2, cs6inv, cs8inv, t2cs4inv,       &
    &                                cs10inv, cs12inv, div2_3, div1_36, div1_6, &
    &                                div1_9, div2_9, div1_108, div1_18, div1_8, &
    &                                div1_6
  use tem_aux_module,          only: tem_abort

  ! include musubi modules
  use mus_field_prop_module,          only: mus_field_prop_type
  use mus_scheme_layout_module,       only: mus_scheme_layout_type
  use mus_param_module,               only: mus_param_type
  use mus_varSys_module,              only: mus_varSys_data_type
  use mus_derVarPos_module,           only: mus_derVarPos_type
  use mus_directions_module,          only: qN00, q0N0, q00N, q100, q010, q001, &
    &                                       q0NN, q0N1, q01N, q011, qN0N, q10N, &
    &                                       qN01, q101, qNN0, qN10, q1N0, q110, &
    &                                       qNNN, qNN1, qN1N, qN11, q1NN, q1N1, &
    &                                       q11N, q111
  use mus_gradData_module,           only: mus_gradData_type
  use mus_derivedQuantities_module2, only: secondMom_3D
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_hrrInit_module,            only: HRR_Correction_d3q27, &
    &                                      getHermitepolynomials

  implicit none

  private

  public :: mus_advRel_kFluid_rBGK_vImproved_lD3Q27
  public :: mus_advRel_kCFD_rBGK_vStd_lD3Q27
  public :: mus_advRel_kFluid_rTRT_vStd_lD3Q27
  public :: bgk_Regularized_d3q27
  public :: bgk_RecursiveRegularized_d3q27
  public :: bgk_HybridRecursiveRegularized_d3q27
  public :: bgk_HybridRecursiveRegularizedCorr_d3q27
  public :: bgk_ProjectedRecursiveRegularized_d3q27
  public :: bgk_DualRelaxationTime_RR_d3q27
  public :: mus_intp_getPdfs_D3Q27

  integer, parameter :: QQ = 27  !< number of pdf directions
  integer, parameter :: q000 = 27

contains

! ****************************************************************************** !
  !> Improved BGK model (with Galilean correction term)
  !! taken from Martin Geier cumulent paper 2015
  !! Geier, M., Schönherr, M., Pasquali, A., & Krafczyk, M. (2015).
  !! The cumulant lattice Boltzmann equation in three dimensions : Theory and
  !! validation. Computers and Mathematics with Applications.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rBGK_vImproved_lD3Q27( fieldProp, inState,      &
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
    integer :: nScalars
    real(kind=rk) :: f(-1:1,-1:1,-1:1)
    real(kind=rk) :: u, v, w , u2, v2, w2
    real(kind=rk) :: rho, rho_omg
    real(kind=rk) :: inv_rho
    real(kind=rk) :: omega, cmpl_o, fac
    real(kind=rk) :: sumX1, sumXN, X0, X1, XN
    real(kind=rk) :: sumY1, sumYN, Y0, Y1, YN
    real(kind=rk) :: sumZ1, sumZN, Z0, Z1, ZN
    real(kind=rk) :: m200, m020, m002
    real(kind=rk) :: Gx, Gy, Gz
    real(kind=rk) :: x0y0, x0y1, x0yn, x1y0, xny0, x1y1, x1yn, xny1, xnyn
    integer :: dens_pos, vel_pos(3), elemOff
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

!$omp do schedule(static)
    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
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

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u = auxField(elemOff + vel_pos(1))
      v = auxField(elemOff + vel_pos(2))
      w = auxField(elemOff + vel_pos(3))

      ! u v square
      u2 = u*u
      v2 = v*v
      w2 = w*w

      sumX1 = sum(f( 1, :, :))
      sumXN = sum(f(-1, :, :))

      sumY1 = sum(f(:, 1, :))
      sumYN = sum(f(:,-1, :))

      sumZ1 = sum(f(:, :, 1))
      sumZN = sum(f(:, :,-1))

      ! second moments, by equation A.7, A.8 and A.9
      inv_rho = 1.0_rk / rho
      m200 = (sumX1 + sumXN) * inv_rho
      m020 = (sumY1 + sumYN) * inv_rho
      m002 = (sumZ1 + sumZN) * inv_rho

      ! relaxation rate
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1._rk - omega
      fac    = 4.5_rk - 2.25_rk * omega

      !----------------------------------------------------------------
      ! calculate Galilean correction term, by equation A.13, A.14, A.15
      Gx = fac * u2 * ( m200 - div1_3 - u2 )
      Gy = fac * v2 * ( m020 - div1_3 - v2 )
      Gz = fac * w2 * ( m002 - div1_3 - w2 )
      !----------------------------------------------------------------

      ! X Y Z components of eq
      ! by equation A.19 - A.21
      X0 = -div2_3 + u2 + Gx
      X1 = - ( X0 + 1.0_rk + u ) * 0.5_rk
      XN = X1 + u

      Y0 = -div2_3 + v2 + Gy
      Y1 = - ( Y0 + 1.0_rk + v ) * 0.5_rk
      YN = Y1 + v

      Z0 = -div2_3 + w2 + Gz
      Z1 = - ( Z0 + 1.0_rk + w ) * 0.5_rk
      ZN = Z1 + w

      ! rho * omega
      rho_omg = rho * omega
      X0 = -rho_omg * X0
      X1 = -rho_omg * X1
      XN = -rho_omg * XN

! fEq000 = X0 * Y0 * Z0
! fEq00N = X0 * Y0 * ZN
! fEq001 = X0 * Y0 * Z1
X0Y0 = X0 * Y0
outState( (ielem-1)*qq+ q000+(1-1)*qq) = cmpl_o*f( 0, 0, 0) + X0Y0*Z0
outState( (ielem-1)*qq+ q00n+(1-1)*qq) = cmpl_o*f( 0, 0,-1) + X0Y0*ZN
outState( (ielem-1)*qq+ q001+(1-1)*qq) = cmpl_o*f( 0, 0, 1) + X0Y0*Z1

! fEq010 = X0 * Y1 * Z0
! fEq011 = X0 * Y1 * Z1
! fEq01N = X0 * Y1 * ZN
X0Y1 = X0 * Y1
outState( (ielem-1)*qq+ q010+(1-1)*qq) = cmpl_o*f( 0, 1, 0) + X0Y1*Z0
outState( (ielem-1)*qq+ q011+(1-1)*qq) = cmpl_o*f( 0, 1, 1) + X0Y1*Z1
outState( (ielem-1)*qq+ q01n+(1-1)*qq) = cmpl_o*f( 0, 1,-1) + X0Y1*ZN

! fEq0N0 = X0 * YN * Z0
! fEq0N1 = X0 * YN * Z1
! fEq0NN = X0 * YN * ZN
X0YN = X0 * YN
outState( (ielem-1)*qq+ q0n0+(1-1)*qq) = cmpl_o*f( 0,-1, 0) + X0YN*Z0
outState( (ielem-1)*qq+ q0n1+(1-1)*qq) = cmpl_o*f( 0,-1, 1) + X0YN*Z1
outState( (ielem-1)*qq+ q0nn+(1-1)*qq) = cmpl_o*f( 0,-1,-1) + X0YN*ZN

! fEq100 = X1 * Y0 * Z0
! fEq10N = X1 * Y0 * ZN
! fEq101 = X1 * Y0 * Z1
X1Y0 = X1 * Y0
outState( (ielem-1)*qq+ q100+(1-1)*qq) = cmpl_o*f( 1, 0, 0) + X1Y0*Z0
outState( (ielem-1)*qq+ q101+(1-1)*qq) = cmpl_o*f( 1, 0, 1) + X1Y0*Z1
outState( (ielem-1)*qq+ q10n+(1-1)*qq) = cmpl_o*f( 1, 0,-1) + X1Y0*ZN

! fEqN00 = XN * Y0 * Z0
! fEqN01 = XN * Y0 * Z1
! fEqN0N = XN * Y0 * ZN
XNY0 = XN * Y0
outState( (ielem-1)*qq+ qn00+(1-1)*qq) = cmpl_o*f(-1, 0, 0) + XNY0*Z0
outState( (ielem-1)*qq+ qn01+(1-1)*qq) = cmpl_o*f(-1, 0, 1) + XNY0*Z1
outState( (ielem-1)*qq+ qn0n+(1-1)*qq) = cmpl_o*f(-1, 0,-1) + XNY0*ZN

! fEq110 = X1 * Y1 * Z0
! fEq111 = X1 * Y1 * Z1
! fEq11N = X1 * Y1 * ZN
X1Y1 = X1 * Y1
outState( (ielem-1)*qq+ q110+(1-1)*qq) = cmpl_o*f( 1, 1, 0) + X1Y1*Z0
outState( (ielem-1)*qq+ q111+(1-1)*qq) = cmpl_o*f( 1, 1, 1) + X1Y1*Z1
outState( (ielem-1)*qq+ q11n+(1-1)*qq) = cmpl_o*f( 1, 1,-1) + X1Y1*ZN

! fEq1N0 = X1 * YN * Z0
! fEq1N1 = X1 * YN * Z1
! fEq1NN = X1 * YN * ZN
X1YN = X1 * YN
outState( (ielem-1)*qq+ q1n0+(1-1)*qq) = cmpl_o*f( 1,-1, 0) + X1YN*Z0
outState( (ielem-1)*qq+ q1n1+(1-1)*qq) = cmpl_o*f( 1,-1, 1) + X1YN*Z1
outState( (ielem-1)*qq+ q1nn+(1-1)*qq) = cmpl_o*f( 1,-1,-1) + X1YN*ZN

! fEqN10 = XN * Y1 * Z0
! fEqN11 = XN * Y1 * Z1
! fEqN1N = XN * Y1 * ZN
XNY1 = XN * Y1
outState( (ielem-1)*qq+ qn10+(1-1)*qq) = cmpl_o*f(-1, 1, 0) + XNY1*Z0
outState( (ielem-1)*qq+ qn11+(1-1)*qq) = cmpl_o*f(-1, 1, 1) + XNY1*Z1
outState( (ielem-1)*qq+ qn1n+(1-1)*qq) = cmpl_o*f(-1, 1,-1) + XNY1*ZN

! fEqNN0 = XN * YN * Z0
! fEqNN1 = XN * YN * Z1
! fEqNNN = XN * YN * ZN
XNYN = XN * YN
outState( (ielem-1)*qq+ qnn0+(1-1)*qq) = cmpl_o*f(-1,-1, 0) + XNYN*Z0
outState( (ielem-1)*qq+ qnn1+(1-1)*qq) = cmpl_o*f(-1,-1, 1) + XNYN*Z1
outState( (ielem-1)*qq+ qnnn+(1-1)*qq) = cmpl_o*f(-1,-1,-1) + XNYN*ZN

    end do nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluid_rBGK_vImproved_lD3Q27
! ****************************************************************************** !


! ****************************************************************************** !
  !> Advection relaxation routine for the D3Q27 model with BGK with
  !! standard equilibrium function
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kCFD_rBGK_vStd_lD3Q27( fieldProp, inState, outState,  &
    &                                          auxField, neigh, nElems,       &
    &                                          nSolve, level, layout, params, &
    &                                          varSys, derVarPos )
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
    integer :: iElem, nScalars
    real(kind=rk) :: pdfTmp(QQ)
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: vel(3)  ! local velocity
    real(kind=rk) :: fEq(QQ) !< equilibrium distribution
    real(kind=rk) :: omega
    integer :: dens_pos, vel_pos(3), elemOff
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars


    !$omp do schedule(static)

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector

    nodeloop: do iElem = 1,nSolve
      ! First load all local values into temp array
      ! Generic! PUSH+PULL is possible
      pdfTmp( qN00 ) = inState(neigh (( qn00-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q0N0 ) = inState(neigh (( q0n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q00N ) = inState(neigh (( q00n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q100 ) = inState(neigh (( q100-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q010 ) = inState(neigh (( q010-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q001 ) = inState(neigh (( q001-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      pdfTmp( q0NN ) = inState(neigh (( q0nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q0N1 ) = inState(neigh (( q0n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q01N ) = inState(neigh (( q01n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q011 ) = inState(neigh (( q011-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( qN0N ) = inState(neigh (( qn0n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q10N ) = inState(neigh (( q10n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( qN01 ) = inState(neigh (( qn01-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q101 ) = inState(neigh (( q101-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( qNN0 ) = inState(neigh (( qnn0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( qN10 ) = inState(neigh (( qn10-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q1N0 ) = inState(neigh (( q1n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q110 ) = inState(neigh (( q110-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      pdfTmp( qNNN ) = inState(neigh (( qnnn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( qNN1 ) = inState(neigh (( qnn1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( qN1N ) = inState(neigh (( qn1n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( qN11 ) = inState(neigh (( qn11-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q1NN ) = inState(neigh (( q1nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q1N1 ) = inState(neigh (( q1n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q11N ) = inState(neigh (( q11n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( q111 ) = inState(neigh (( q111-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      pdfTmp( q000 ) = inState(neigh (( q000-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      vel(1) = auxField(elemOff + vel_pos(1))
      vel(2) = auxField(elemOff + vel_pos(2))
      vel(3) = auxField(elemOff + vel_pos(3))


      ! Calculate the equilibrium distribution function
      fEq = layout%quantities%pdfEq_ptr(rho = rho, vel = vel, QQ = QQ)

      ! omega
      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)

      ! relaxation
      outState((ielem-1)*qq+qn00+(1-1)*qq) =                     &
        &      pdfTmp( qN00 ) - omega*(pdfTmp( qN00 ) - fEq( qN00 ))
      outState((ielem-1)*qq+q0n0+(1-1)*qq) =                     &
        &      pdfTmp( q0N0 ) - omega*(pdfTmp( q0N0 ) - fEq( q0N0 ))
      outState((ielem-1)*qq+q00n+(1-1)*qq) =                     &
        &      pdfTmp( q00N ) - omega*(pdfTmp( q00N ) - fEq( q00N ))
      outState((ielem-1)*qq+q100+(1-1)*qq) =                     &
        &      pdfTmp( q100 ) - omega*(pdfTmp( q100 ) - fEq( q100 ))
      outState((ielem-1)*qq+q010+(1-1)*qq) =                     &
        &      pdfTmp( q010 ) - omega*(pdfTmp( q010 ) - fEq( q010 ))
      outState((ielem-1)*qq+q001+(1-1)*qq) =                     &
        &      pdfTmp( q001 ) - omega*(pdfTmp( q001 ) - fEq( q001 ))
      outState((ielem-1)*qq+q0nn+(1-1)*qq)  =                     &
        &      pdfTmp( q0NN ) - omega*(pdfTmp( q0NN ) - fEq( q0NN ))
      outState((ielem-1)*qq+q0n1+(1-1)*qq)  =                     &
        &      pdfTmp( q0N1 ) - omega*(pdfTmp( q0N1 ) - fEq( q0N1 ))
      outState((ielem-1)*qq+q01n+(1-1)*qq)  =                     &
        &      pdfTmp( q01N ) - omega*(pdfTmp( q01N ) - fEq( q01N ))
      outState((ielem-1)*qq+q011+(1-1)*qq)  =                     &
        &      pdfTmp( q011 ) - omega*(pdfTmp( q011 ) - fEq( q011 ))
      outState((ielem-1)*qq+qn0n+(1-1)*qq)  =                     &
        &      pdfTmp( qN0N ) - omega*(pdfTmp( qN0N ) - fEq( qN0N ))
      outState((ielem-1)*qq+q10n+(1-1)*qq)  =                     &
        &      pdfTmp( q10N ) - omega*(pdfTmp( q10N ) - fEq( q10N ))
      outState((ielem-1)*qq+qn01+(1-1)*qq)  =                     &
        &      pdfTmp( qN01 ) - omega*(pdfTmp( qN01 ) - fEq( qN01 ))
      outState((ielem-1)*qq+q101+(1-1)*qq)  =                     &
        &      pdfTmp( q101 ) - omega*(pdfTmp( q101 ) - fEq( q101 ))
      outState((ielem-1)*qq+qnn0+(1-1)*qq)  =                     &
        &      pdfTmp( qNN0 ) - omega*(pdfTmp( qNN0 ) - fEq( qNN0 ))
      outState((ielem-1)*qq+qn10+(1-1)*qq)  =                     &
        &      pdfTmp( qN10 ) - omega*(pdfTmp( qN10 ) - fEq( qN10 ))
      outState((ielem-1)*qq+q1n0+(1-1)*qq)  =                     &
        &      pdfTmp( q1N0 ) - omega*(pdfTmp( q1N0 ) - fEq( q1N0 ))
      outState((ielem-1)*qq+q110+(1-1)*qq)  =                     &
        &      pdfTmp( q110 ) - omega*(pdfTmp( q110 ) - fEq( q110 ))

      outState((ielem-1)*qq+qnnn+(1-1)*qq)  =                     &
        &      pdfTmp( qNNN ) - omega*(pdfTmp( qNNN ) - fEq( qNNN ))
      outState((ielem-1)*qq+qnn1+(1-1)*qq)  =                     &
        &      pdfTmp( qNN1 ) - omega*(pdfTmp( qNN1 ) - fEq( qNN1 ))
      outState((ielem-1)*qq+qn1n+(1-1)*qq)  =                     &
        &      pdfTmp( qN1N ) - omega*(pdfTmp( qN1N ) - fEq( qN1N ))
      outState((ielem-1)*qq+qn11+(1-1)*qq)  =                     &
        &      pdfTmp( qN11 ) - omega*(pdfTmp( qN11 ) - fEq( qN11 ))
      outState((ielem-1)*qq+q1nn+(1-1)*qq)  =                     &
        &      pdfTmp( q1NN ) - omega*(pdfTmp( q1NN ) - fEq( q1NN ))
      outState((ielem-1)*qq+q1n1+(1-1)*qq)  =                     &
        &      pdfTmp( q1N1 ) - omega*(pdfTmp( q1N1 ) - fEq( q1N1 ))
      outState((ielem-1)*qq+q11n+(1-1)*qq)  =                     &
        &      pdfTmp( q11N ) - omega*(pdfTmp( q11N ) - fEq( q11N ))
      outState((ielem-1)*qq+q111+(1-1)*qq)  =                     &
        &      pdfTmp( q111 ) - omega*(pdfTmp( q111 ) - fEq( q111 ))

      outState((ielem-1)*qq+q000+(1-1)*qq) =                      &
        &      pdfTmp( q000 ) - omega*(pdfTmp( q000 ) - fEq( q000 ))


    enddo nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kCFD_rBGK_vStd_lD3Q27
! ****************************************************************************** !

! ****************************************************************************** !
  !> No comment yet!
  !!
  !! TODO add comment
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  ! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rTRT_vStd_lD3Q27( fieldProp, inState, outState, &
    &                                            auxField, neigh, nElems,      &
    &                                            nSolve, level, layout,        &
    &                                            params, varSys, derVarPos )
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
    integer :: nScalars
    real(kind=rk) :: f(-1:1,-1:1,-1:1)
    real(kind=rk) :: u, v, w, u2, v2, w2
    real(kind=rk) :: rho, inv_rho!, rho_omg
    real(kind=rk) :: wP, wN, p_part, n_part
    real(kind=rk) :: X0, X1, XN
    real(kind=rk) :: Y0, Y1, YN
    real(kind=rk) :: Z0, Z1, ZN
    integer :: dens_pos, vel_pos(3), elemOff
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
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
      u = auxField(elemOff + vel_pos(1))
      v = auxField(elemOff + vel_pos(2))
      w = auxField(elemOff + vel_pos(3))
      ! rho_omg = rho * omega
      ! calculate rho and velocity ----------------------------------------------

      ! u v square
      u2 = u*u
      v2 = v*v
      w2 = w*w

      ! X Y Z components of eq
      ! by equation A.19 - A.21
      X0 = -div2_3 + u2
      X1 = - ( X0 + 1.0_rk + u ) * 0.5_rk
      XN = X1 + u

      Y0 = -div2_3 + v2
      Y1 = - ( Y0 + 1.0_rk + v ) * 0.5_rk
      YN = Y1 + v

      Z0 = -div2_3 + w2
      Z1 = - ( Z0 + 1.0_rk + w ) * 0.5_rk
      ZN = Z1 + w

      ! X0 = -rho_omg * X0
      ! X1 = -rho_omg * X1
      ! XN = -rho_omg * XN

      wP = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      wN = 1.0_rk / ( 0.5_rk + fieldProp(1)%fluid%lambda &
        &                   / ( 1.0_rk/wP - 0.5_rk )     )

! fEq = -rho * X * Y * Z
! fEq000 = X0 * Y0 * Z0
outState( (ielem-1)*qq+ q000+(1-1)*qq) = (1.0_rk - wP)*f(0,0,0) &
  &                                                  - rho*wP*X0*Y0*Z0

p_part = wP * ( ( f(1,0,0) + f(-1,0,0) ) - (-rho*X1*Y0*Z0 - rho*XN*Y0*Z0) ) * div1_2
n_part = wN * ( ( f(1,0,0) - f(-1,0,0) ) - (-rho*X1*Y0*Z0 + rho*XN*Y0*Z0) ) * div1_2
outState( (ielem-1)*qq+ q100+(1-1)*qq) = f( 1,0,0) - p_part - n_part
outState( (ielem-1)*qq+ qn00+(1-1)*qq) = f(-1,0,0) - p_part + n_part

p_part = wP * ( ( f(0,1,0) + f(0,-1,0) ) - (-rho*X0*Y1*Z0 - rho*X0*YN*Z0) ) * div1_2
n_part = wN * ( ( f(0,1,0) - f(0,-1,0) ) - (-rho*X0*Y1*Z0 + rho*X0*YN*Z0) ) * div1_2
outState( (ielem-1)*qq+ q010+(1-1)*qq) = f(0, 1,0) - p_part - n_part
outState( (ielem-1)*qq+ q0n0+(1-1)*qq) = f(0,-1,0) - p_part + n_part

p_part = wP * ( ( f(0,0,1) + f(0,0,-1) ) - (-rho*X0*Y0*Z1 - rho*X0*Y0*ZN) ) * div1_2
n_part = wN * ( ( f(0,0,1) - f(0,0,-1) ) - (-rho*X0*Y0*Z1 + rho*X0*Y0*ZN) ) * div1_2
outState( (ielem-1)*qq+ q001+(1-1)*qq) = f(0,0, 1) - p_part - n_part
outState( (ielem-1)*qq+ q00n+(1-1)*qq) = f(0,0,-1) - p_part + n_part

p_part = wP * ( ( f(0,1,1) + f(0,-1,-1) ) - (-rho*X0*Y1*Z1 - rho*X0*YN*ZN) ) * div1_2
n_part = wN * ( ( f(0,1,1) - f(0,-1,-1) ) - (-rho*X0*Y1*Z1 + rho*X0*YN*ZN) ) * div1_2
outState( (ielem-1)*qq+ q011+(1-1)*qq) = f(0, 1, 1) - p_part - n_part
outState( (ielem-1)*qq+ q0nn+(1-1)*qq) = f(0,-1,-1) - p_part + n_part

p_part = wP * ( ( f(0,1,-1) + f(0,-1,1) ) - (-rho*X0*Y1*ZN - rho*X0*YN*Z1) ) * div1_2
n_part = wN * ( ( f(0,1,-1) - f(0,-1,1) ) - (-rho*X0*Y1*ZN + rho*X0*YN*Z1) ) * div1_2
outState( (ielem-1)*qq+ q01n+(1-1)*qq) = f(0, 1,-1) - p_part - n_part
outState( (ielem-1)*qq+ q0n1+(1-1)*qq) = f(0,-1, 1) - p_part + n_part

p_part = wP * ( ( f(1,0,1) + f(-1,0,-1) ) - (-rho*X1*Y0*Z1 - rho*XN*Y0*ZN) ) * div1_2
n_part = wN * ( ( f(1,0,1) - f(-1,0,-1) ) - (-rho*X1*Y0*Z1 + rho*XN*Y0*ZN) ) * div1_2
outState( (ielem-1)*qq+ q101+(1-1)*qq) = f( 1,0, 1) - p_part - n_part
outState( (ielem-1)*qq+ qn0n+(1-1)*qq) = f(-1,0,-1) - p_part + n_part

p_part = wP * ( ( f(1,0,-1) + f(-1,0,1) ) - (-rho*X1*Y0*ZN - rho*XN*Y0*Z1) ) * div1_2
n_part = wN * ( ( f(1,0,-1) - f(-1,0,1) ) - (-rho*X1*Y0*ZN + rho*XN*Y0*Z1) ) * div1_2
outState( (ielem-1)*qq+ q10n+(1-1)*qq) = f( 1,0,-1) - p_part - n_part
outState( (ielem-1)*qq+ qn01+(1-1)*qq) = f(-1,0, 1) - p_part + n_part

p_part = wP * ( ( f(1,1,0) + f(-1,-1,0) ) - (-rho*X1*Y1*Z0 - rho*XN*YN*Z0) ) * div1_2
n_part = wN * ( ( f(1,1,0) - f(-1,-1,0) ) - (-rho*X1*Y1*Z0 + rho*XN*YN*Z0) ) * div1_2
outState( (ielem-1)*qq+ q110+(1-1)*qq) = f( 1, 1,0) - p_part - n_part
outState( (ielem-1)*qq+ qnn0+(1-1)*qq) = f(-1,-1,0) - p_part + n_part

p_part = wP * ( ( f(1,-1,0) + f(-1,1,0) ) - (-rho*X1*YN*Z0 - rho*XN*Y1*Z0) ) * div1_2
n_part = wN * ( ( f(1,-1,0) - f(-1,1,0) ) - (-rho*X1*YN*Z0 + rho*XN*Y1*Z0) ) * div1_2
outState( (ielem-1)*qq+ q1n0+(1-1)*qq) = f( 1,-1,0) - p_part - n_part
outState( (ielem-1)*qq+ qn10+(1-1)*qq) = f(-1, 1,0) - p_part + n_part

p_part = wP * ( ( f(1,-1,-1) + f(-1,1,1) ) - (-rho*X1*YN*ZN - rho*XN*Y1*Z1) ) * div1_2
n_part = wN * ( ( f(1,-1,-1) - f(-1,1,1) ) - (-rho*X1*YN*ZN + rho*XN*Y1*Z1) ) * div1_2
outState( (ielem-1)*qq+ q1nn+(1-1)*qq) = f( 1,-1,-1) - p_part - n_part
outState( (ielem-1)*qq+ qn11+(1-1)*qq) = f(-1, 1, 1) - p_part + n_part

p_part = wP * ( ( f(1,1,-1) + f(-1,-1,1) ) - (-rho*X1*Y1*ZN - rho*XN*YN*Z1) ) * div1_2
n_part = wN * ( ( f(1,1,-1) - f(-1,-1,1) ) - (-rho*X1*Y1*ZN + rho*XN*YN*Z1) ) * div1_2
outState( (ielem-1)*qq+ q11n+(1-1)*qq) = f( 1, 1,-1) - p_part - n_part
outState( (ielem-1)*qq+ qnn1+(1-1)*qq) = f(-1,-1, 1) - p_part + n_part

p_part = wP * ( ( f(1,-1,1) + f(-1,1,-1) ) - (-rho*X1*YN*Z1 - rho*XN*Y1*ZN) ) * div1_2
n_part = wN * ( ( f(1,-1,1) - f(-1,1,-1) ) - (-rho*X1*YN*Z1 + rho*XN*Y1*ZN) ) * div1_2
outState( (ielem-1)*qq+ q1n1+(1-1)*qq) = f( 1,-1, 1) - p_part - n_part
outState( (ielem-1)*qq+ qn1n+(1-1)*qq) = f(-1, 1,-1) - p_part + n_part

p_part = wP * ( ( f(1,1,1) + f(-1,-1,-1) ) - (-rho*X1*Y1*Z1 - rho*XN*YN*ZN) ) * div1_2
n_part = wN * ( ( f(1,1,1) - f(-1,-1,-1) ) - (-rho*X1*Y1*Z1 + rho*XN*YN*ZN) ) * div1_2
outState( (ielem-1)*qq+ q111+(1-1)*qq) = f( 1, 1, 1) - p_part - n_part
outState( (ielem-1)*qq+ qnnn+(1-1)*qq) = f(-1,-1,-1) - p_part + n_part

    end do nodeloop

  end subroutine mus_advRel_kFluid_rTRT_vStd_lD3Q27
! ****************************************************************************** !


! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
pure subroutine f_f_eq_regularized_2nd_ord_d3q27 ( weight, rho, u_x, u_y, u_z, feq, &
  &  f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )
    ! -------------------------------------------------------------------- !
    !> weights of the stencil
    real(kind=rk), intent(in) :: weight(QQ)
    !> density, velocity components
    real(kind=rk), intent(in) :: rho
    real(kind=rk), intent(in) :: u_x
    real(kind=rk), intent(in) :: u_y
    real(kind=rk), intent(in) :: u_z
    !> equilibrium pdf and full pdf
    real(kind=rk), intent(out) :: feq(QQ)
    real(kind=rk), intent(out) :: f1(QQ)
    !> coefficients of f1: a12xx, a12yy, a12xy, etc ...
    real(kind=rk), intent(in) :: a12xx, a12yy, a12zz, a12xy, a12xz, a12yz
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: u_x_sqr, u_y_sqr, u_z_sqr, u_x_u_y, u_x_u_z, u_y_u_z
    real(kind=rk) :: f00, f01, f02, f12
    ! ---------------------------------------------------------------------------

      u_x_sqr = u_x**2
      u_y_sqr = u_y**2
      u_z_sqr = u_z**2
      u_x_u_y = u_x * u_y
      u_x_u_z = u_x * u_z
      u_y_u_z = u_y * u_z

      f00 = 1.0_rk

      !iDir = 1
      f01 = cs2inv*(-u_x)
      f02 = div1_6*cs4inv*(2._rk*u_x_sqr - (u_y_sqr + u_z_sqr))
      feq(1) = weight(1) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*(2._rk*a12xx - (a12yy + a12zz))
      f1(1) = weight(1) * f12

      !iDir = 4
      f01 = -f01
      feq(4) = weight(4) * rho * (f00 + f01 + f02)
      f1(4) = weight(4) * f12

      !iDir = 2
      f01 = cs2inv*(-u_y)
      f02 = div1_6*cs4inv*( -(u_x_sqr + u_z_sqr) + 2._rk*u_y_sqr )
      feq(2) = weight(2) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*( -(a12xx + a12zz) + 2._rk*a12yy )
      f1(2) = weight(2) * f12

      !iDir = 5
      f01 = -f01
      feq(5) = weight(5) * rho * (f00 + f01 + f02)
      f1(5) = weight(5) * f12

      !iDir = 3
      f01 = cs2inv*(-u_z)
      f02 = div1_6*cs4inv*( -(u_x_sqr + u_y_sqr) + 2._rk*u_z_sqr )
      feq(3) = weight(3) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*( -(a12xx + a12yy) + 2._rk*a12zz)
      f1(3) = weight(3) * f12

      !iDir = 6
      f01 = -f01
      feq(6) = weight(6) * rho * (f00 + f01 + f02)
      f1(6) = weight(6) * f12

      !iDir = 7
      f01 = cs2inv*(-u_y - u_z)
      f02 = div1_6*cs4inv*(-u_x_sqr + 2._rk*(u_y_sqr + u_z_sqr) + 6._rk*u_y_u_z )
      feq(7) = weight(7) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*( -a12xx + 2._rk*(a12yy + a12zz) + 6._rk*a12yz )
      f1(7) = weight(7) * f12

      !iDir = 10
      f01 = -f01
      feq(10) = weight(10) * rho * (f00 + f01 + f02)
      f1(10) = weight(10) * f12

      !iDir = 8
      f01 = cs2inv*(-u_y + u_z)
      f02 = f02 - 2._rk*cs4inv*( u_y_u_z )
      feq(8) = weight(8) * rho * (f00 + f01 + f02)
      f12 = f12 - 2._rk*cs4inv*( a12yz )
      f1(8) = weight(8) * f12

      !iDir = 9
      f01 = -f01
      feq(9) = weight(9) * rho * (f00 + f01 + f02)
      f1(9) = weight(9) * f12

      !iDir = 11
      f01 = cs2inv*(-u_x - u_z)
      f02 = div1_6*cs4inv*( 2._rk*(u_x_sqr + u_z_sqr) - u_y_sqr + 6._rk*u_x_u_z )
      feq(11) = weight(11) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*( 2._rk*(a12xx + a12zz) - a12yy + 6._rk*a12xz )
      f1(11) = weight(11) * f12

      !iDir = 14
      f01 = -f01
      feq(14) = weight(14) * rho * (f00 + f01 + f02)
      f1(14) = weight(14) * f12

      !iDir = 12
      f01 = cs2inv*(u_x - u_z)
      f02 = f02 - 2._rk*cs4inv*( u_x_u_z )
      feq(12) = weight(12) * rho * (f00 + f01 + f02)
      f12 = f12 - 2._rk*cs4inv*( a12xz )
      f1(12) = weight(12) * f12

      !iDir = 13
      f01 = -f01
      feq(13) = weight(13) * rho * (f00 + f01 + f02)
      f1(13) = weight(13) * f12

      !iDir = 15
      f01 = cs2inv*(-u_x - u_y)
      f02 = div1_6*cs4inv*( 2._rk*(u_x_sqr + u_y_sqr) - u_z_sqr + 6._rk*u_x_u_y )
      feq(15) = weight(15) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*( 2._rk*(a12xx + a12yy) - a12zz + 6.0_rk*a12xy)
      f1(15) = weight(15) * f12

      !iDir = 18
      f01 = -f01
      feq(18) = weight(18) * rho * (f00 + f01 + f02)
      f1(18) = weight(18) * f12

      !iDir = 16
      f01 = cs2inv*(-u_x + u_y)
      f02 = f02 - 2._rk*cs4inv*( u_x_u_y )
      feq(16) = weight(16) * rho * (f00 + f01 + f02)
      f12 = f12 - 2._rk*cs4inv*( a12xy )
      f1(16) = weight(16) * f12

      !iDir = 17
      f01 = -f01
      feq(17) = weight(17) * rho * (f00 + f01 + f02)
      f1(17) = weight(17) * f12

      !iDir = 19
      f01 = cs2inv*(-u_x - u_y - u_z)
      f02 = cs4inv*( div1_3 * (u_x_sqr + u_y_sqr + u_z_sqr) &
        &                      + (u_x_u_y + u_x_u_z + u_y_u_z) )
      feq(19) = weight(19) * rho * (f00 + f01 + f02)
      f12 = cs4inv*( div1_3*(a12xx + a12yy + a12zz) + (a12xy + a12xz + a12yz) )
      f1(19) = weight(19) * f12

      !iDir = 26
      f01 = -f01
      feq(26) = weight(26) * rho * (f00 + f01 + f02)
      f1(26) = weight(26) * f12

      !iDir = 22
      f01 = f01 - 2._rk*cs2inv*( u_x )
      f02 = f02 + 2._rk*cs4inv*( -u_x_u_y - u_x_u_z )
      feq(22) = weight(22) * rho * (f00 + f01 + f02)
      f12 = f12 + 2._rk*cs4inv*( -a12xy - a12xz )
      f1(22) = weight(22) * f12

      !iDir = 23
      f01 = -f01
      feq(23) = weight(23) * rho * (f00 + f01 + f02)
      f1(23) = weight(23) * f12

      !iDir = 25
      f01 = f01 + 2._rk*cs2inv*u_y
      f02 = f02 + 2._rk*cs4inv*( u_x_u_y - u_y_u_z )
      feq(25) = weight(25) * rho * (f00 + f01 + f02)
      f12 = f12 + 2._rk*cs4inv*( a12xy - a12yz )
      f1(25) = weight(25) * f12

      !iDir = 20
      f01 = -f01
      feq(20) = weight(20) * rho * (f00 + f01 + f02)
      f1(20) = weight(20) * f12

      !iDir = 24
      f01 = f01 + 2._rk*cs2inv*( u_x )
      f02 = f02 + 2._rk*cs4inv*( -u_x_u_y + u_x_u_z )
      feq(24) = weight(24) * rho * (f00 + f01 + f02)
      f12 = f12 + 2._rk*cs4inv*( -a12xy + a12xz )
      f1(24) = weight(24) * f12

      !iDir = 21
      f01 = -f01
      feq(21) = weight(21) * rho * (f00 + f01 + f02)
      f1(21) = weight(21) * f12

      !iDir = 27
      !f01 = 0._rk
      f02 = -div1_6*cs4inv*(u_x_sqr + u_y_sqr + u_z_sqr)
      feq(27) = weight(27) * rho * (f00 + f02)
      f12 = -div1_6*cs4inv*(a12xx + a12yy + a12zz)
      f1(27) = weight(27) * f12

  end subroutine f_f_eq_regularized_2nd_ord_d3q27
! ****************************************************************************** !


! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
pure subroutine f_f_eq_regularized_4th_ord_d3q27 ( weight, rho, u_x, u_y, u_z, feq, &
  &  f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )
    ! -------------------------------------------------------------------- !
    !> weights of the stencil
    real(kind=rk), intent(in) :: weight(QQ)
    !> density, velocity components
    real(kind=rk), intent(in) :: rho
    real(kind=rk), intent(in) :: u_x
    real(kind=rk), intent(in) :: u_y
    real(kind=rk), intent(in) :: u_z
    !> equilibrium pdf and full pdf
    real(kind=rk), intent(out) :: feq(QQ)
    real(kind=rk), intent(out) :: f1(QQ)
    !> coefficients of f1: a12xx, a12yy, a12xy, etc ...
    real(kind=rk), intent(in) :: a12xx, a12yy, a12zz, a12xy, a12xz, a12yz
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: u_x_sqr_u_y, u_y_sqr_u_x, f03, u_z_sqr_u_x, u_z_sqr_u_y
    real(kind=rk) :: u_x_sqr, u_y_sqr, u_z_sqr, u_x_sqr_u_z, u_y_sqr_u_z
    real(kind=rk) :: u_x_u_y, u_x_u_z, u_y_u_z, u_x_u_y_u_z, f04, f05, f06
    real(kind=rk) :: u_x_sqr_u_y_sqr, u_x_sqr_u_z_sqr, u_y_sqr_u_z_sqr
    real(kind=rk) :: u_x_sqr_u_y_sqr_u_z, u_x_sqr_u_y_u_z_sqr, u_x_u_y_sqr_u_z_sqr
    real(kind=rk) :: u_x_sqr_u_y_sqr_u_z_sqr, f13, f14, f15, f16
    real(kind=rk) :: a13xyy, a13xxy, a13xxz, a13xzz, a13yzz, a13yyz, a13xyz
    real(kind=rk) :: a14xxyy, a14xxzz, a14yyzz
    real(kind=rk) :: a15xxyzz, a15xxyyz, a15xyyzz, a16xxyyzz
    real(kind=rk) :: u_x_u_y_u_z_sqr, u_x_u_y_sqr_u_z, u_x_sqr_u_y_u_z
    real(kind=rk) :: a14xyzz, a14xyyz, a14xxyz
    ! ---------------------------------------------------------------------------

      call f_f_eq_regularized_2nd_ord_d3q27 ( weight, rho, u_x, u_y, u_z, feq, f1, &
        &    a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      u_x_sqr = u_x**2
      u_y_sqr = u_y**2
      u_z_sqr = u_z**2
      u_x_u_y = u_x * u_y
      u_x_u_z = u_x * u_z
      u_y_u_z = u_y * u_z
      u_x_sqr_u_y = u_x_sqr * u_y
      u_x_sqr_u_z = u_x_sqr * u_z
      u_y_sqr_u_x = u_y_sqr * u_x
      u_y_sqr_u_z = u_y_sqr * u_z
      u_z_sqr_u_x = u_z_sqr * u_x
      u_z_sqr_u_y = u_z_sqr * u_y
      u_x_u_y_u_z = u_x_u_y * u_z
      u_x_sqr_u_y_sqr = u_x_sqr * u_y_sqr
      u_x_sqr_u_z_sqr = u_x_sqr * u_z_sqr
      u_y_sqr_u_z_sqr = u_y_sqr * u_z_sqr
      u_x_u_y_u_z_sqr = u_x_u_y * u_z_sqr
      u_x_u_y_sqr_u_z = u_x * u_y_sqr_u_z
      u_x_sqr_u_y_u_z = u_x_sqr * u_y_u_z
      u_x_sqr_u_y_sqr_u_z = u_x_sqr_u_y_sqr * u_z
      u_x_sqr_u_y_u_z_sqr = u_x_sqr_u_z_sqr * u_y
      u_x_u_y_sqr_u_z_sqr = u_y_sqr_u_z_sqr * u_x
      u_x_sqr_u_y_sqr_u_z_sqr = u_x_sqr_u_y_sqr * u_z_sqr

      a13xxy = 2.0_rk * u_x * a12xy + u_y * a12xx
      a13xxz = 2.0_rk * u_x * a12xz + u_z * a12xx
      a13xyy = u_x * a12yy + 2.0_rk * u_y * a12xy
      a13xzz = u_x * a12zz + 2.0_rk * u_z * a12xz
      a13yzz = u_y * a12zz + 2.0_rk * u_z * a12yz
      a13yyz = 2.0_rk * u_y * a12yz + u_z * a12yy
      a13xyz = u_x * a12yz + u_y * a12xz + u_z * a12xy

      a14xxyy = u_x_sqr * a12yy + 2.0_rk*u_x_u_y*a12xy + u_y *a13xxy
      a14xxzz = u_x_sqr * a12zz + 2.0_rk*u_x_u_z*a12xz + u_z *a13xxz
      a14yyzz = u_y_sqr * a12zz + 2.0_rk*u_y_u_z*a12yz + u_z *a13yyz
      a14xyzz = u_x_u_y * a12zz + u_x_u_z * a12yz + u_y_u_z * a12xz + u_z * a13xyz
      a14xyyz = 2.0_rk * u_x_u_y * a12yz + u_y_sqr * a12xz + u_z * a13xyy
      a14xxyz = u_x_sqr * a12yz + 2.0_rk * u_x_u_y * a12xz + u_z * a13xxy

      a15xxyzz = u_x_sqr_u_y * a12zz + u_x_sqr_u_z * a12yz &
        &          + 2.0_rk * u_x_u_y_u_z * a12xz + u_z * a14xxyz
      a15xxyyz = 2.0_rk * u_x_sqr_u_y * a12yz + 2.0_rk * u_y_sqr_u_x * a12xz  &
        &          + u_z * a14xxyy
      a15xyyzz = u_y_sqr_u_x * a12zz + 2.0_rk * u_x_u_y_u_z * a12yz &
        &          + u_y_sqr_u_z * a12xz + u_z*a14xyyz

      a16xxyyzz = u_x_sqr_u_y_sqr * a12zz + 2.0_rk * u_x_sqr_u_y_u_z * a12yz &
        &          + 2.0_rk * u_x_u_y_sqr_u_z * a12xz + u_z * a15xxyyz

      !iDir = 1
      f03 = div1_6*cs6inv*( (u_y_sqr_u_x + u_z_sqr_u_x) )
      f04 = div1_36*cs8inv*( -2._rk*(u_x_sqr_u_y_sqr + u_x_sqr_u_z_sqr) &
        &     + u_y_sqr_u_z_sqr )
      f05 = div1_36*cs10inv*( -u_x_u_y_sqr_u_z_sqr )
      f06 = div1_108*cs12inv*( u_x_sqr_u_y_sqr_u_z_sqr )
      feq(1) = feq(1) + weight(1) * rho * (f03 + f04 + f05 + f06)
      f13 = div1_6*cs6inv*( (a13xyy + a13xzz) )
      f14 = div1_36*cs8inv*( -2._rk*(a14xxyy + a14xxzz) + a14yyzz)
      f15 = div1_36*cs10inv*( -a15xyyzz )
      f16 = div1_108*cs12inv*( a16xxyyzz )
      f1(1) = f1(1) + weight(1) * ( f13 + f14 + f15 + f16 )

      !iDir = 4
      f03 = -f03
      f05 = -f05
      feq(4) = feq(4) + weight(4) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(4) = f1(4) + weight(4) * ( f13 + f14 + f15 + f16 )

      !iDir = 2
      f03 = div1_6*cs6inv*( (u_x_sqr_u_y + u_z_sqr_u_y) )
      f04 = div1_36*cs8inv*( -2._rk*(u_x_sqr_u_y_sqr + u_y_sqr_u_z_sqr) &
        &     + u_x_sqr_u_z_sqr )
      f05 = div1_36*cs10inv*( -u_x_sqr_u_y_u_z_sqr )
      feq(2) = feq(2) + weight(2) * rho * (f03 + f04 + f05 + f06)
      f13 = div1_6*cs6inv*( (a13xxy + a13yzz) )
      f14 = div1_36*cs8inv*( -2._rk*(a14xxyy + a14yyzz) + a14xxzz )
      f15 = div1_36*cs10inv*( -a15xxyzz )
      f1(2) = f1(2) + weight(2) * ( f13 + f14 + f15 + f16 )

      !iDir = 5
      f03 = -f03
      f05 = -f05
      feq(5) = feq(5) + weight(5) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(5) = f1(5) + weight(5) * ( f13 + f14 + f15 + f16 )

      !iDir = 3
      f03 = div1_6*cs6inv*( (u_x_sqr_u_z + u_y_sqr_u_z) )
      f04 = div1_36*cs8inv*( u_x_sqr_u_y_sqr - 2._rk*(u_x_sqr_u_z_sqr &
        &     + u_y_sqr_u_z_sqr) )
      f05 = div1_36*cs10inv*( -u_x_sqr_u_y_sqr_u_z )
      feq(3) = feq(3) + weight(3) * rho * (f03 + f04 + f05 + f06)
      f13 = div1_6*cs6inv*( (a13xxz + a13yyz) )
      f14 = div1_36*cs8inv*( a14xxyy - 2._rk*(a14xxzz + a14yyzz) )
      f15 = div1_36*cs10inv*( -a15xxyyz )
      f1(3) = f1(3) + weight(3) * ( f13 + f14 + f15 + f16 )

      !iDir = 6
      f03 = -f03
      f05 = -f05
      feq(6) = feq(6) + weight(6) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(6) = f1(6) + weight(6) * ( f13 + f14 + f15 + f16 )

      !iDir = 7
      f03 = div1_6*cs6inv*( (u_x_sqr_u_y + u_x_sqr_u_z) &
        &      - 2._rk*(u_z_sqr_u_y + u_y_sqr_u_z) )
      f04 = div1_18*cs8inv*( -(u_x_sqr_u_y_sqr + u_x_sqr_u_z_sqr) &
        &     + 2._rk*u_y_sqr_u_z_sqr - 3._rk*u_x_sqr_u_y_u_z)
      f05 = div1_18*cs10inv*( (u_x_sqr_u_y_u_z_sqr + u_x_sqr_u_y_sqr_u_z) )
      f06 = -2._rk * f06
      feq(7) = feq(7) + weight(7) * rho * (f03 + f04 + f05 + f06)
      f13 = div1_6*cs6inv*( (a13xxy + a13xxz) - 2._rk*(a13yzz + a13yyz) )
      f14 = div1_18*cs8inv*( -(a14xxyy + a14xxzz) + 2._rk*a14yyzz &
        &                    -3._rk*a14xxyz)
      f15 = div1_18*cs10inv*( (a15xxyzz + a15xxyyz) )
      f16 = -2._rk * f16
      f1(7) = f1(7) + weight(7) * ( f13 + f14 + f15 + f16 )

      !iDir = 10
      f03 = -f03
      f05 = -f05
      feq(10) = feq(10) + weight(10) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(10) = f1(10) + weight(10) * ( f13 + f14 + f15 + f16 )

      !iDir = 8
      f03 = f03 + div1_3*cs6inv*( u_x_sqr_u_y - 2._rk*u_z_sqr_u_y )
      f04 = f04 + div1_3*cs8inv*u_x_sqr_u_y_u_z
      f05 = div1_18*cs10inv*( (u_x_sqr_u_y_u_z_sqr - u_x_sqr_u_y_sqr_u_z) )
      feq(8) = feq(8) + weight(8) * rho * (f03 + f04 + f05 + f06)
      f13 = f13 + div1_3*cs6inv* ( a13xxy - 2._rk*a13yzz )
      f14 = f14 + div1_3*cs8inv*a14xxyz
      f15 = div1_18*cs10inv*( (a15xxyzz - a15xxyyz) )
      f1(8) = f1(8) + weight(8) * ( f13 + f14 + f15 + f16 )

      !iDir = 9
      f03 = -f03
      f05 = -f05
      feq(9) = feq(9) + weight(9) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(9) = f1(9) + weight(9) * ( f13 + f14 + f15 + f16 )

      !iDir = 11
      f03 = div1_6*cs6inv*( -2._rk*(u_x_sqr_u_z + u_z_sqr_u_x) + (u_y_sqr_u_x &
        &      + u_y_sqr_u_z) )
      f04 = div1_18*cs8inv*( -(u_x_sqr_u_y_sqr + u_y_sqr_u_z_sqr) &
        &     + 2._rk*u_x_sqr_u_z_sqr - 3._rk*u_x_u_y_sqr_u_z)
      f05 = div1_18*cs10inv*( (u_x_sqr_u_y_sqr_u_z + u_x_u_y_sqr_u_z_sqr) )
      feq(11) = feq(11) + weight(11) * rho * (f03 + f04 + f05 + f06)
      f13 = div1_6*cs6inv*( -2._rk*(a13xxz + a13xzz) + (a13xyy + a13yyz) )
      f14 = div1_18*cs8inv*( -(a14xxyy + a14yyzz) + 2._rk*a14xxzz &
        &                   -3._rk*a14xyyz)
      f15 = div1_18*cs10inv*( (a15xxyyz + a15xyyzz) )
      f1(11) = f1(11) + weight(11) * ( f13 + f14 + f15 + f16 )

      !iDir = 14
      f03 = -f03
      f05 = -f05
      feq(14) = feq(14) + weight(14) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(14) = f1(14) + weight(14) * ( f13 + f14 + f15 + f16 )

      !iDir = 12
      f03 = f03 + div1_3*cs6inv*( -2._rk*u_x_sqr_u_z + u_y_sqr_u_z)
      f04 = f04 + div1_3*cs8inv*u_x_u_y_sqr_u_z
      f05 = div1_18*cs10inv*( (u_x_sqr_u_y_sqr_u_z - u_x_u_y_sqr_u_z_sqr) )
      feq(12) = feq(12) + weight(12) * rho * (f03 + f04 + f05 + f06)
      f13 = f13 + div1_3*cs6inv*( -2._rk*a13xxz + a13yyz)
      f14 = f14 + div1_3*cs8inv*a14xyyz
      f15 = div1_18*cs10inv*( (a15xxyyz - a15xyyzz) )
      f1(12) = f1(12) + weight(12) * ( f13 + f14 + f15 + f16 )

      !iDir = 13
      f03 = -f03
      f05 = -f05
      feq(13) = feq(13) + weight(13) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(13) = f1(13) + weight(13) * ( f13 + f14 + f15 + f16 )

      !iDir = 15
      f03 = div1_6*cs6inv*( -2._rk*(u_x_sqr_u_y + u_y_sqr_u_x) &
        &      + (u_z_sqr_u_x + u_z_sqr_u_y) )
      f04 = div1_18*cs8inv*( 2._rk*u_x_sqr_u_y_sqr - (u_x_sqr_u_z_sqr &
        &     + u_y_sqr_u_z_sqr) - 3._rk*u_x_u_y_u_z_sqr)
      f05 = div1_18*cs10inv*( (u_x_sqr_u_y_u_z_sqr + u_x_u_y_sqr_u_z_sqr) )
      feq(15) = feq(15) + weight(15) * rho * (f03 + f04 + f05 + f06)
      f13 = div1_6*cs6inv*( -2._rk*(a13xxy + a13xyy) + (a13xzz + a13yzz) )
      f14 = div1_18*cs8inv*( 2._rk*a14xxyy - (a14xxzz + a14yyzz) &
        &                   -3._rk*a14xyzz)
      f15 = div1_18*cs10inv*( (a15xxyzz + a15xyyzz) )
      f1(15) = f1(15) + weight(15) * ( f13 + f14 + f15 + f16 )

      !iDir = 18
      f03 = -f03
      f05 = -f05
      feq(18) = feq(18) + weight(18) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(18) = f1(18) + weight(18) * ( f13 + f14 + f15 + f16 )

      !iDir = 16
      f03 = f03 + div1_3*cs6inv*( -2._rk*u_y_sqr_u_x + u_z_sqr_u_x )
      f04 = f04 + div1_3*cs8inv*u_x_u_y_u_z_sqr
      f05 = div1_18*cs10inv*( (-u_x_sqr_u_y_u_z_sqr + u_x_u_y_sqr_u_z_sqr) )
      feq(16) = feq(16) + weight(16) * rho * (f03 + f04 + f05 + f06)
      f13 = f13 + div1_3*cs6inv*( -2._rk*a13xyy + a13xzz )
      f14 = f14 + div1_3*cs8inv*a14xyzz
      f15 = div1_18*cs10inv*( (-a15xxyzz + a15xyyzz) )
      f1(16) = f1(16) + weight(16) * ( f13 + f14 + f15 + f16 )

      !iDir = 17
      f03 = -f03
      f05 = -f05
      feq(17) = feq(17) + weight(17) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(17) = f1(17) + weight(17) * ( f13 + f14 + f15 + f16 )

      !iDir = 19
      f03 = cs6inv*( -div1_3*(u_x_sqr_u_y + u_x_sqr_u_z + u_y_sqr_u_x &
        &      + u_z_sqr_u_x + u_z_sqr_u_y + u_y_sqr_u_z) &
        &      - u_x_u_y_u_z )
      f04 = div1_9*cs8inv*( (u_x_sqr_u_y_sqr + u_x_sqr_u_z_sqr &
        &     + u_y_sqr_u_z_sqr) + 3._rk*( u_x_u_y_u_z_sqr + u_x_u_y_sqr_u_z &
        &                                  + u_x_sqr_u_y_u_z) )
      f05 = div1_9*cs10inv*( -(u_x_sqr_u_y_u_z_sqr + u_x_sqr_u_y_sqr_u_z &
        &     + u_x_u_y_sqr_u_z_sqr) )
      f06 = -2._rk*f06
      feq(19) = feq(19) + weight(19) * rho * (f03 + f04 + f05 + f06)
      f13 = cs6inv*( -div1_3*(a13xxy + a13xxz + a13xyy + a13xzz + a13yzz &
        &     + a13yyz) - a13xyz )
      f14 = div1_9*cs8inv*( (a14xxyy + a14xxzz + a14yyzz) &
        &                   + 3._rk*(a14xyzz + a14xyyz + a14xxyz) )
      f15 = div1_9*cs10inv*( -(a15xxyzz + a15xxyyz + a15xyyzz) )
      f16 = -2._rk*f16
      f1(19) = f1(19) + weight(19) * ( f13 + f14 + f15 + f16 )

      !iDir = 26
      f03 = -f03
      f05 = -f05
      feq(26) = feq(26) + weight(26) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(26) = f1(26) + weight(26) * ( f13 + f14 + f15 + f16 )

      !iDir = 22
      f03 = f03 + 2._rk*cs6inv*( -div1_3*(u_y_sqr_u_x + u_z_sqr_u_x) &
        &      - u_x_u_y_u_z )
      f04 = f04 + div2_3*cs8inv*( -(u_x_u_y_u_z_sqr + u_x_u_y_sqr_u_z) )
      f05 = f05 + div2_9*cs10inv*( -u_x_u_y_sqr_u_z_sqr )
      feq(22) = feq(22) + weight(22) * rho * (f03 + f04 + f05 + f06)
      f13 = f13 + 2._rk*cs6inv*( -div1_3*(a13xyy + a13xzz) - a13xyz )
      f14 = f14 + div2_3*cs8inv*( -(a14xyzz + a14xyyz) )
      f15 = f15 + div2_9*cs10inv*( -a15xyyzz )
      f1(22) = f1(22) + weight(22) * ( f13 + f14 + f15 + f16 )

      !iDir = 23
      f03 = -f03
      f05 = -f05
      feq(23) = feq(23) + weight(23) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(23) = f1(23) + weight(23) * ( f13 + f14 + f15 + f16 )

      !iDir = 25
      f03 = f03 + 2._rk*cs6inv*( div1_3*(u_x_sqr_u_y + u_z_sqr_u_y) - u_x_u_y_u_z )
      f04 = f04 + div2_3*cs8inv*( u_x_u_y_u_z_sqr - u_x_sqr_u_y_u_z )
      f05 = f05 + div2_9*cs10inv*( u_x_sqr_u_y_u_z_sqr )
      feq(25) = feq(25) + weight(25) * rho * (f03 + f04 + f05 + f06)
      f13 = f13 + 2._rk*cs6inv*( div1_3*(a13xxy + a13yzz) - a13xyz )
      f14 = f14 + div2_3*cs8inv*( a14xyzz - a14xxyz )
      f15 = f15 + div2_9*cs10inv*( a15xxyzz )
      f1(25) = f1(25) + weight(25) * ( f13 + f14 + f15 + f16 )

      !iDir = 20
      f03 = -f03
      f05 = -f05
      feq(20) = feq(20) + weight(20) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(20) = f1(20) + weight(20) * ( f13 + f14 + f15 + f16 )

      !iDir = 24
      f03 = f03 + 2._rk*cs6inv*( div1_3*(u_y_sqr_u_x + u_z_sqr_u_x) - u_x_u_y_u_z )
      f04 = f04 + div2_3*cs8inv*( -u_x_u_y_u_z_sqr + u_x_u_y_sqr_u_z )
      f05 = f05 + div2_9*cs10inv*( u_x_u_y_sqr_u_z_sqr )
      feq(24) = feq(24) + weight(24) * rho * (f03 + f04 + f05 + f06)
      f13 = f13 + 2._rk*cs6inv*( div1_3*(a13xyy + a13xzz) - a13xyz )
      f14 = f14 + div2_3*cs8inv*( -a14xyzz + a14xyyz )
      f15 = f15 + div2_9*cs10inv*( a15xyyzz )
      f1(24) = f1(24) + weight(24) * ( f13 + f14 + f15 + f16 )

      !iDir = 21
      f03 = -f03
      f05 = -f05
      feq(21) = feq(21) + weight(21) * rho * (f03 + f04 + f05 + f06)
      f13 = -f13
      f15 = -f15
      f1(21) = f1(21) + weight(21) * ( f13 + f14 + f15 + f16 )

      !iDir = 27
      !f03 = 0._rk
      f04 = div1_36*cs8inv*( u_x_sqr_u_y_sqr + u_x_sqr_u_z_sqr + u_y_sqr_u_z_sqr )
      !f05 = 0._rk
      f06 = -div1_8 * f06
      feq(27) = feq(27) + weight(27) * rho * (f04 + f06)
      !f13 = 0._rk
      f14 = div1_36*cs8inv*( (a14xxyy + a14xxzz + a14yyzz) )
      !f15 = 0._rk
      f16 = -div1_8 * f16
      f1(27) = f1(27) + weight(27) * ( f14 + f16 )

  end subroutine f_f_eq_regularized_4th_ord_d3q27
! ****************************************************************************** !

! **************************************************************************** !
  !> Regularized relaxation routine for the D3Q19 and 27 model with BGK.
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005

  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  !! works for both d3q19 and d3q27
  subroutine bgk_Regularized_d3q27( fieldProp, inState, outState, auxField, &
    &                          neigh, nElems, nSolve, level, layout,   &
    &                          params, varSys, derVarPos )
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
    real(kind=rk) :: f( QQ ), SOM(6), SOM_neq(6)
    real(kind=rk) :: rho, u_x, u_y, u_z, a12xx, a12xy, a12yy, a12zz, a12xz, a12yz
    real(kind=rk) :: omega, cmpl_o, feq(QQ), f1(QQ)
    integer :: denspos, velpos(3), elemOff, nScalars
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:3) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

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
      u_z = auxField(elemOff + velpos(3))

      ! non equilibrium second-order moments
      ! SOM_neq = SOM - SOM_eq
      ! SOM = Second order moments
      ! 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=xz
      SOM = secondMom_3D(layout%fStencil%cxcx, f, layout%fStencil%QQ)
      SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      SOM_neq(3) = SOM(3) - rho * (cs2 + (u_z * u_z))
      SOM_neq(4) = SOM(4) - rho * u_x * u_y
      SOM_neq(5) = SOM(5) - rho * u_y * u_z
      SOM_neq(6) = SOM(6) - rho * u_x * u_z

      ! Hermitian coefficients
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1.0_rk - omega
      a12xx = SOM_neq(1)
      a12yy = SOM_neq(2)
      a12zz = SOM_neq(3)
      a12xy = SOM_neq(4)
      a12yz = SOM_neq(5)
      a12xz = SOM_neq(6)

      call f_f_eq_regularized_2nd_ord_d3q27 ( layout%weight(:), rho, u_x, u_y, u_z, feq, f1, &
        &    a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                      + cmpl_o*f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_Regularized_d3q27
! **************************************************************************** !

! **************************************************************************** !
  !> Recursive Regularized relaxation routine for the D3Q27
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005

  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_RecursiveRegularized_d3q27( fieldProp, inState, outState, auxField, &
    &                          neigh, nElems, nSolve, level, layout,   &
    &                          params, varSys, derVarPos )
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
    real(kind=rk) :: f( QQ ), SOM(6), SOM_neq(6)
    real(kind=rk) :: rho, u_x, u_y, u_z, a12xx, a12xy, a12yy, a12zz, a12xz, a12yz
    real(kind=rk) :: omega, cmpl_o, feq(QQ), f1(QQ)
    integer :: denspos, velpos(3), elemOff, nScalars
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:3) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

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
      u_z = auxField(elemOff + velpos(3))

      ! non equilibrium second-order moments
      ! SOM_neq = SOM - SOM_eq
      ! SOM = Second order moments
      ! 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=xz
      SOM = secondMom_3D(layout%fStencil%cxcx, f, layout%fStencil%QQ)
      SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      SOM_neq(3) = SOM(3) - rho * (cs2 + (u_z * u_z))
      SOM_neq(4) = SOM(4) - rho * u_x * u_y
      SOM_neq(5) = SOM(5) - rho * u_y * u_z
      SOM_neq(6) = SOM(6) - rho * u_x * u_z

      ! Hermitian coefficients
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1.0_rk - omega
      a12xx = SOM_neq(1)
      a12yy = SOM_neq(2)
      a12zz = SOM_neq(3)
      a12xy = SOM_neq(4)
      a12yz = SOM_neq(5)
      a12xz = SOM_neq(6)

      call f_f_eq_regularized_4th_ord_d3q27 ( layout%weight(:), rho, u_x, u_y, u_z, feq, &
        &  f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                      + cmpl_o*f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_RecursiveRegularized_d3q27
! **************************************************************************** !


! **************************************************************************** !
  !> Projected Recursive Regularized relaxation routine for the D3Q27
  ! based on High-order extension of the recursive regularized lattice
  ! Boltzmann method, PhD Thesis, COREIXAS 2018
  ! the suggested basis in this thesis has an error.
  ! GGS: fixed by adding 4th order terms xxyz
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_ProjectedRecursiveRegularized_d3q27( fieldProp, inState, outState, auxField, &
    &                          neigh, nElems, nSolve, level, layout,   &
    &                          params, varSys, derVarPos )
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
    !> gradient data
    type(mus_gradData_type), pointer :: gradData
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! temporary distribution variables
    real(kind=rk) :: f( QQ ), SR(6), gradU(3,3,1)!!TODO::, tr_SR
    real(kind=rk) :: rho, u_x, u_y, u_z, a12xx, a12xy, a12yy, a12zz, a12xz, a12yz
    real(kind=rk) :: omega, taup, cmpl_o, feq(QQ), f1(QQ)
    integer :: denspos, velpos(3), elemOff, nScalars
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:3) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

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
      u_z = auxField(elemOff + velpos(3))

      ! Stress tensor components
      gradU(:,:,1:1) = scheme%Grad%U_ptr(         &
           &   auxField     = auxField,           &
           &   gradData     = gradData,           &
           &   velPos       = velPos,             &
           &   nAuxScalars  = varSys%nAuxScalars, &
           &   nDims        = 3,                  &
           &   nSolve       = 1,                  &
           &   elemOffset   = iElem -1            )

      ! symmetric strain rate tensors
      ! transformed inro RHS of a1 FD equation
      ! the trace is needed only for energy conservation, which is not
      ! done in Musubi at the moment
      !!TODO:tr_SR = div2_3 * (gradU(1, 1, 1) + gradU(2, 2, 1) + gradU(3, 3, 1))
      SR(1) = 2._rk * gradU(1, 1, 1)         !S_XX !!TODO:- tr_SR
      SR(2) = 2._rk * gradU(2, 2, 1)         !S_YY !!TODO:- tr_SR
      SR(3) = 2._rk * gradU(3, 3, 1)         !S_ZZ !!TODO:- tr_SR
      SR(4) = gradU(1, 2, 1)+gradU(2, 1, 1)  !S_XY
      SR(5) = gradU(2, 3, 1)+gradU(3, 2, 1)  !S_YZ
      SR(6) = gradU(1, 3, 1)+gradU(3, 1, 1)  !S_XZ

      ! non equilibrium second-order moments
      ! SOM_neq = SOM - SOM_eq
      ! SOM = Second order moments
      ! 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=xz
      !SOM = secondMom_3D(layout%fStencil%cxcx, f, layout%fStencil%QQ)
      !SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      !SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      !SOM_neq(3) = SOM(3) - rho * (cs2 + (u_z * u_z))
      !SOM_neq(4) = SOM(4) - rho * u_x * u_y
      !SOM_neq(5) = SOM(5) - rho * u_y * u_z
      !SOM_neq(6) = SOM(6) - rho * u_x * u_z

      ! Hermitian coefficients
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1.0_rk - omega
      taup = rho * cs2 / omega
      a12xx = -taup * SR(1)
      a12yy = -taup * SR(2)
      a12zz = -taup * SR(3)
      a12xy = -taup * SR(4)
      a12yz = -taup * SR(5)
      a12xz = -taup * SR(6)

      call f_f_eq_regularized_4th_ord_d3q27 ( layout%weight(:), rho, u_x, u_y, u_z, feq, &
        &  f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                      + cmpl_o*f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_ProjectedRecursiveRegularized_d3q27
! **************************************************************************** !

! **************************************************************************** !
  !> Projected Recursive Regularized relaxation routine for the D3Q19
  ! based on High-order extension of the recursive regularized lattice
  ! Boltzmann method, PhD Thesis, COREIXAS 2018
  ! the suggested basis in this thesis has an error.
  ! GGS: fixed by adding 4th order terms xxyz
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  !! works for both d3q19 and d3q27
  subroutine bgk_HybridRecursiveRegularized_d3q27( fieldProp, inState, outState, auxField, &
    &                          neigh, nElems, nSolve, level, layout,   &
    &                          params, varSys, derVarPos )
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
    !> gradient data
    type(mus_gradData_type), pointer :: gradData
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! temporary distribution variables
    real(kind=rk) :: f( QQ ), SR(6), gradU(3,3,1), SOM(6), SOM_neq(6)!!TODO:, tr_SR
    real(kind=rk) :: rho, u_x, u_y, u_z, a12xx, a12xy, a12yy, a12zz, a12xz, a12yz
    real(kind=rk) :: omega, taup, cmpl_o, feq(QQ), f1(QQ), sigma
    integer :: denspos, velpos(3), elemOff, nScalars
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:3) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

    !call getHermitepolynomials( 3, QQ, layout, 6)

    ! access gradData
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( 1 )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme
    gradData => scheme%gradData(level)

    ! sigma value, ideally read from input
    sigma = fieldProp(1)%fluid%HRR_sigma ! if sigma == 1 --> no HRR

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
      u_z = auxField(elemOff + velpos(3))

      ! Stress tensor components
      gradU(:,:,1:1) = scheme%Grad%U_ptr(         &
           &   auxField     = auxField,           &
           &   gradData     = gradData,           &
           &   velPos       = velPos,             &
           &   nAuxScalars  = varSys%nAuxScalars, &
           &   nDims        = 3,                  &
           &   nSolve       = 1,                  &
           &   elemOffset   = iElem -1            )

      ! symmetric strain rate tensors
      ! transformed inro RHS of a1 FD equation
      ! the trace is needed only for energy conservation, which is not
      ! done in Musubi at the moment
      !!TODO:tr_SR = div2_3 * (gradU(1, 1, 1) + gradU(2, 2, 1) + gradU(3, 3, 1))
      SR(1) = 2._rk * gradU(1, 1, 1)         !S_XX !!TODO:- tr_SR
      SR(2) = 2._rk * gradU(2, 2, 1)         !S_YY !!TODO:- tr_SR
      SR(3) = 2._rk * gradU(3, 3, 1)         !S_ZZ !!TODO:- tr_SR
      SR(4) = gradU(1, 2, 1)+gradU(2, 1, 1)  !S_XY
      SR(5) = gradU(2, 3, 1)+gradU(3, 2, 1)  !S_YZ
      SR(6) = gradU(1, 3, 1)+gradU(3, 1, 1)  !S_XZ

      ! non equilibrium second-order moments
      ! SOM_neq = SOM - SOM_eq
      ! SOM = Second order moments
      ! 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=xz
      SOM = secondMom_3D(layout%fStencil%cxcx, f, layout%fStencil%QQ)
      SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      SOM_neq(3) = SOM(3) - rho * (cs2 + (u_z * u_z))
      SOM_neq(4) = SOM(4) - rho * u_x * u_y
      SOM_neq(5) = SOM(5) - rho * u_y * u_z
      SOM_neq(6) = SOM(6) - rho * u_x * u_z

      ! Hermitian coefficients
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      cmpl_o = 1.0_rk - omega
      taup = rho * cs2 / omega
      a12xx = SOM_neq(1) * sigma + (1.0_rk - sigma) * (-taup * SR(1))
      a12yy = SOM_neq(2) * sigma + (1.0_rk - sigma) * (-taup * SR(2))
      a12zz = SOM_neq(3) * sigma + (1.0_rk - sigma) * (-taup * SR(3))
      a12xy = SOM_neq(4) * sigma + (1.0_rk - sigma) * (-taup * SR(4))
      a12yz = SOM_neq(5) * sigma + (1.0_rk - sigma) * (-taup * SR(5))
      a12xz = SOM_neq(6) * sigma + (1.0_rk - sigma) * (-taup * SR(6))

      call f_f_eq_regularized_4th_ord_d3q27 ( layout%weight(:), rho, u_x, u_y, u_z, feq, &
        &  f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                      + cmpl_o*f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_HybridRecursiveRegularized_d3q27
! **************************************************************************** !

! **************************************************************************** !
  !> Projected Recursive Regularized relaxation routine for the D3Q19
  ! based on High-order extension of the recursive regularized lattice
  ! Boltzmann method, PhD Thesis, COREIXAS 2018
  ! the suggested basis in this thesis has an error.
  ! GGS: fixed by adding 4th order terms xxyz
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  !! works for both d3q19 and d3q27
  subroutine bgk_HybridRecursiveRegularizedCorr_d3q27( fieldProp, inState, outState, auxField, &
    &                          neigh, nElems, nSolve, level, layout,   &
    &                          params, varSys, derVarPos )
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
    !> gradient data
    type(mus_gradData_type), pointer :: gradData
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! temporary distribution variables
    real(kind=rk) :: f( QQ ), SR(6), gradU(3,3,1), SOM(6), SOM_neq(6)
    real(kind=rk) :: rho, u_x, u_y, u_z, a12xx, a12xy, a12yy, a12zz, a12xz, a12yz
    real(kind=rk) :: omega, taup, cmpl_o, feq(QQ), f1(QQ), sigma
    real(kind=rk) :: gradRhoU3(3,1), S_Corr(QQ)!!TODO:, tr_SR
    integer :: denspos, velpos(3), elemOff, nScalars, iSrc
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:3) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

    ! access gradData
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( 1 )%method_data, &
      &               fPtr )
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
        u_z = auxField(elemOff + velpos(3))

        ! Stress tensor components
        gradU(:,:,1:1) = scheme%Grad%U_ptr(         &
             &   auxField     = auxField,           &
             &   gradData     = gradData,           &
             &   velPos       = velPos,             &
             &   nAuxScalars  = varSys%nAuxScalars, &
             &   nDims        = 3,                  &
             &   nSolve       = 1,                  &
             &   elemOffset   = iElem -1            )

          ! 1 = x, 2 = y, 3 = z, no xy returned
          gradRhoU3(:,1:1) = scheme%Grad%RhoU3_ptr( &
            &   auxField     = auxField,            &
            &   gradData     = gradData,            &
            &   velPos       = velpos,              &
            &   densPos      = denspos,             &
            &   nAuxScalars  = varSys%nAuxScalars,  &
            &   nDims        = 3,                   &
            &   nSolve       = 1,                   &
            &   elemOffset   = iElem-1              )

          ! Calculate correction
          call HRR_Correction_d3q27 (               &
            &    QQ         = QQ,                   &
            &    weight     = layout%weight(:),     &
            &    gradRHOU3  = gradRHOU3(:, 1),      &
            &    phi        = S_corr(:),            &
            &    dens       = HRR_Corr%dens(iElem), &
            &    vel        = HRR_Corr%vel(iElem,:) )

        ! symmetric strain rate tensors
        ! transformed inro RHS of a1 FD equation
        ! the trace is needed only for energy conservation, which is not
        ! done in Musubi at the moment
        !!TODO:tr_SR = div2_3 * (gradU(1, 1, 1) + gradU(2, 2, 1) + gradU(3, 3, 1))
        SR(1) = 2._rk * gradU(1, 1, 1)         !S_XX !!TODO:- tr_SR
        SR(2) = 2._rk * gradU(2, 2, 1)         !S_YY !!TODO:- tr_SR
        SR(3) = 2._rk * gradU(3, 3, 1)         !S_ZZ !!TODO:- tr_SR
        SR(4) = gradU(1, 2, 1)+gradU(2, 1, 1)  !S_XY
        SR(5) = gradU(2, 3, 1)+gradU(3, 2, 1)  !S_YZ
        SR(6) = gradU(1, 3, 1)+gradU(3, 1, 1)  !S_XZ

        ! non equilibrium second-order moments
        ! SOM_neq = SOM - SOM_eq
        ! Apply correction
        f(:) = f(:) + 0.5_rk * S_corr(:)
        ! SOM = Second order moments
        ! 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=xz
        SOM = secondMom_3D(layout%fStencil%cxcx, f, layout%fStencil%QQ)
        SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
        SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
        SOM_neq(3) = SOM(3) - rho * (cs2 + (u_z * u_z))
        SOM_neq(4) = SOM(4) - rho * u_x * u_y
        SOM_neq(5) = SOM(5) - rho * u_y * u_z
        SOM_neq(6) = SOM(6) - rho * u_x * u_z

        ! Hermitian coefficients
        omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
        cmpl_o = 1.0_rk - omega
        taup = rho * cs2 / omega
        a12xx = SOM_neq(1) * sigma + (1.0_rk - sigma) * (-taup * SR(1))
        a12yy = SOM_neq(2) * sigma + (1.0_rk - sigma) * (-taup * SR(2))
        a12zz = SOM_neq(3) * sigma + (1.0_rk - sigma) * (-taup * SR(3))
        a12xy = SOM_neq(4) * sigma + (1.0_rk - sigma) * (-taup * SR(4))
        a12yz = SOM_neq(5) * sigma + (1.0_rk - sigma) * (-taup * SR(5))
        a12xz = SOM_neq(6) * sigma + (1.0_rk - sigma) * (-taup * SR(6))

        call f_f_eq_regularized_4th_ord_d3q27 ( layout%weight(:), rho, u_x, u_y, u_z, feq, &
          &  f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

        do iDir = 1, QQ
          outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
            &                              + cmpl_o*f1(iDir) + 0.5_rk * S_corr(iDir)
        enddo

      enddo nodeloop
!$omp end do nowait

    end associate

  end subroutine bgk_HybridRecursiveRegularizedCorr_d3q27
! **************************************************************************** !


! **************************************************************************** !
  !> This function returns pdfs from state
  pure function mus_intp_getPdfs_D3Q27(state, Neigh, elem, nScalars, nSize) &
    & result(f)
    ! --------------------------------------------------------------------------
    !> State vector
    real(kind=rk), intent(in) :: state(:)
    integer, intent(in) :: Neigh(:)
    !> element position in state array
    integer, intent(in) :: elem
    !> number of scalars in state vector
    integer, intent(in) :: nScalars
    !> Size of state vector
    integer, intent(in) :: nSize
    !> moments
    real(kind=rk) :: f(-1:1,-1:1,-1:1)
    ! --------------------------------------------------------------------------
    integer :: QQ
    ! --------------------------------------------------------------------------
    QQ = 27
    f(-1, 0, 0) = state(( elem-1)* nscalars+ qn00+( 1-1)* qq) !- layout%weight( qN00 )
    f( 0,-1, 0) = state(( elem-1)* nscalars+ q0n0+( 1-1)* qq) !- layout%weight( q0N0 )
    f( 0, 0,-1) = state(( elem-1)* nscalars+ q00n+( 1-1)* qq) !- layout%weight( q00N )
    f( 1, 0, 0) = state(( elem-1)* nscalars+ q100+( 1-1)* qq) !- layout%weight( q100 )
    f( 0, 1, 0) = state(( elem-1)* nscalars+ q010+( 1-1)* qq) !- layout%weight( q010 )
    f( 0, 0, 1) = state(( elem-1)* nscalars+ q001+( 1-1)* qq) !- layout%weight( q001 )
    f( 0,-1,-1) = state(( elem-1)* nscalars+ q0nn+( 1-1)* qq) !- layout%weight( q0NN )
    f( 0,-1, 1) = state(( elem-1)* nscalars+ q0n1+( 1-1)* qq) !- layout%weight( q0N1 )
    f( 0, 1,-1) = state(( elem-1)* nscalars+ q01n+( 1-1)* qq) !- layout%weight( q01N )
    f( 0, 1, 1) = state(( elem-1)* nscalars+ q011+( 1-1)* qq) !- layout%weight( q011 )
    f(-1, 0,-1) = state(( elem-1)* nscalars+ qn0n+( 1-1)* qq) !- layout%weight( qN0N )
    f( 1, 0,-1) = state(( elem-1)* nscalars+ q10n+( 1-1)* qq) !- layout%weight( q10N )
    f(-1, 0, 1) = state(( elem-1)* nscalars+ qn01+( 1-1)* qq) !- layout%weight( qN01 )
    f( 1, 0, 1) = state(( elem-1)* nscalars+ q101+( 1-1)* qq) !- layout%weight( q101 )
    f(-1,-1, 0) = state(( elem-1)* nscalars+ qnn0+( 1-1)* qq) !- layout%weight( qNN0 )
    f(-1, 1, 0) = state(( elem-1)* nscalars+ qn10+( 1-1)* qq) !- layout%weight( qN10 )
    f( 1,-1, 0) = state(( elem-1)* nscalars+ q1n0+( 1-1)* qq) !- layout%weight( q1N0 )
    f( 1, 1, 0) = state(( elem-1)* nscalars+ q110+( 1-1)* qq) !- layout%weight( q110 )
    f(-1,-1,-1) = state(( elem-1)* nscalars+ qnnn+( 1-1)* qq) !- layout%weight( qNNN )
    f(-1,-1, 1) = state(( elem-1)* nscalars+ qnn1+( 1-1)* qq) !- layout%weight( qNN1 )
    f(-1, 1,-1) = state(( elem-1)* nscalars+ qn1n+( 1-1)* qq) !- layout%weight( qN1N )
    f(-1, 1, 1) = state(( elem-1)* nscalars+ qn11+( 1-1)* qq) !- layout%weight( qN11 )
    f( 1,-1,-1) = state(( elem-1)* nscalars+ q1nn+( 1-1)* qq) !- layout%weight( q1NN )
    f( 1,-1, 1) = state(( elem-1)* nscalars+ q1n1+( 1-1)* qq) !- layout%weight( q1N1 )
    f( 1, 1,-1) = state(( elem-1)* nscalars+ q11n+( 1-1)* qq) !- layout%weight( q11N )
    f( 1, 1, 1) = state(( elem-1)* nscalars+ q111+( 1-1)* qq) !- layout%weight( q111 )
    f( 0, 0, 0) = state(( elem-1)* nscalars+ q000+( 1-1)* qq) !- layout%weight( q000 )

  end function mus_intp_getPdfs_D3Q27
! **************************************************************************** !


! **************************************************************************** !
  !> Recursive Regularized relaxation routine for the D3Q27
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005

  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_DualRelaxationTime_RR_d3q27( fieldProp, inState, outState, auxField, &
    &                          neigh, nElems, nSolve, level, layout,   &
    &                          params, varSys, derVarPos )
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
    real(kind=rk) :: f( QQ ), SOM(6), SOM_neq(6)
    real(kind=rk) :: rho, u_x, u_y, u_z, a12xx, a12xy, a12yy, a12zz, a12xz, a12yz
    real(kind=rk) :: omega, tau, tauN, CoefTauNTau, feq(QQ), f_temp, f1(QQ)
    integer :: denspos, velpos(3), elemOff, nScalars
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:3) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

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
      u_z = auxField(elemOff + velpos(3))

      ! non equilibrium second-order moments
      ! SOM_neq = SOM - SOM_eq
      ! SOM = Second order moments
      ! 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=xz
      SOM = secondMom_3D(layout%fStencil%cxcx, f, layout%fStencil%QQ)
      SOM_neq(1) = SOM(1) - rho * (cs2 + (u_x * u_x))
      SOM_neq(2) = SOM(2) - rho * (cs2 + (u_y * u_y))
      SOM_neq(3) = SOM(3) - rho * (cs2 + (u_z * u_z))
      SOM_neq(4) = SOM(4) - rho * u_x * u_y
      SOM_neq(5) = SOM(5) - rho * u_y * u_z
      SOM_neq(6) = SOM(6) - rho * u_x * u_z

      ! Hermitian coefficients
      omega  = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      tau = 1.0_rk / omega
      CoefTauNTau = ( tau - tauN ) / ( tau * tauN )
      a12xx = SOM_neq(1)
      a12yy = SOM_neq(2)
      a12zz = SOM_neq(3)
      a12xy = SOM_neq(4)
      a12yz = SOM_neq(5)
      a12xz = SOM_neq(6)

      call f_f_eq_regularized_2nd_ord_d3q27 ( layout%weight(:), rho, u_x, u_y, u_z, feq, &
        &  f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        f_temp = f(iDir) - 1.0_rk/tauN * ( f(iDir) - feq(iDir) )

        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = &
          & f_temp + CoefTauNTau * f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_DualRelaxationTime_RR_d3q27
! **************************************************************************** !


end module mus_d3q27_module
! ****************************************************************************** !
