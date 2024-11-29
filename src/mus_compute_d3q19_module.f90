! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011, 2013, 2016-2017, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2016, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012, 2014-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
! **************************************************************************** !
!> author: Manuel Hasert
!! author: Kannan Masilamani
!! author: Jiaxing Qi
!! Routines and parameter definitions for the standard D3Q19 model
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
module mus_d3q19_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,              only: rk
  use tem_compileconf_module,  only: vlen
  use tem_varSys_module,       only: tem_varSys_type, tem_varSys_op_type
  use tem_param_module,        only: div1_3, div1_6, div1_36, div1_8, div2_8, &
    &                                div3_4h, cs2inv, cs4inv, t2cs4inv, cs2,  &
    &                                rho0, cs6inv, div2_3
  use tem_dyn_array_module,    only: PositionOfVal
  use tem_aux_module,          only: tem_abort
  use tem_property_module,     only: prp_fluid
  use tem_construction_module, only: tem_levelDesc_type

  ! include musubi modules
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_param_module,              only: mus_param_type
  use mus_varSys_module,             only: mus_varSys_data_type
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_directions_module,         only: qN00, q0N0, q00N, q100, q010, q001, &
    &                                      q0NN, q0N1, q01N, q011, qN0N, q10N, &
    &                                      qN01, q101, qNN0, qN10, q1N0, q110
  use mus_gradData_module,           only: mus_gradData_type
  use mus_derivedQuantities_module2, only: secondMom_3D
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_hrrInit_module,            only: HRR_Correction_d3q19,        &
    &                                      getHermitepolynomials,       &
    &                                      getHermitepolynomials_D3Q19

  implicit none

  private

  public :: mus_advRel_kFluid_rBGK_vStd_lD3Q19
  public :: mus_advRel_kFluidIncomp_rBGK_vStd_lD3Q19
  public :: mus_advRel_kFluid_rTRT_vStd_lD3Q19
  public :: mus_advRel_kFluidIncomp_rTRT_vStd_lD3Q19
  public :: mus_advRel_kFluid_rBGK_vBlock_lD3Q19
  public :: bgk_Regularized_d3q19
  public :: bgk_RecursiveRegularized_d3q19
  public :: bgk_HybridRecursiveRegularized_d3q19
  public :: bgk_ProjectedRecursiveRegularized_d3q19
  public :: bgk_HybridRecursiveRegularizedCorr_d3q19
  public :: bgk_DualRelaxationTime_RR_d3q19

  ! ============================================================================
  ! D3Q19 flow model
  ! ============================================================================
  !> Definition of the discrete velocity set

  ! integer,parameter :: block = 32
  integer,parameter :: QQ   = 19  !< number of pdf directions
  integer,parameter :: q000 = 19  !< rest density is last


contains


! **************************************************************************** !
  !> Advection relaxation routine for the D3Q19 model with BGK.
  !!
  !! \[ f_\alpha(x_i+e_{\alpha,i},t+1) =
  !! f_\alpha(x_i,t) - \omega(f_\alpha(x_i,t)-f^{eq}_{\alpha}(x_i,t)) \]
  !!
  !! The number of floating point operation in this routine is 160 roughly.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rBGK_vStd_lD3Q19( fieldProp, inState, outState, &
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
    integer :: iElem, nScalars
    real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
      &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
      &              f000
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: u_x     ! local x-velocity
    real(kind=rk) :: u_y     ! local y-velocity
    real(kind=rk) :: u_z     ! local z-velocity
    real(kind=rk) :: usq     ! square velocity
    ! derived constants
    real(kind=rk) :: usqn, usqn_o1, usqn_o2
    real(kind=rk) :: omega_2, cmpl_o, omega
    real(kind=rk) :: coeff_1, coeff_2
    real(kind=rk) :: ui1, ui3, ui10, ui11, ui12, ui13
    real(kind=rk) :: fac_1, fac_2, fac_3, fac_4, fac_9, fac_10, fac_11, fac_12,&
      &              fac_13
    real(kind=rk) :: sum1_1, sum1_2, sum2_1, sum2_2, sum3_1, sum3_2, sum4_1,   &
      &              sum4_2, sum9_1, sum9_2, sum10_1, sum10_2, sum11_1,        &
      &              sum11_2, sum12_1, sum12_2, sum13_1, sum13_2
    integer :: dens_pos, vel_pos(3), elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    ! nElems = size(neigh)/QQ
    nScalars = varSys%nScalars


!$omp do schedule(static)

!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1, nSolve

      ! First load all local values into temp array
      fN00 = inState( neigh (( qn00-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f0N0 = inState( neigh (( q0n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f00N = inState( neigh (( q00n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f100 = inState( neigh (( q100-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f010 = inState( neigh (( q010-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f001 = inState( neigh (( q001-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f0NN = inState( neigh (( q0nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f0N1 = inState( neigh (( q0n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f01N = inState( neigh (( q01n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f011 = inState( neigh (( q011-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      fN0N = inState( neigh (( qn0n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f10N = inState( neigh (( q10n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      fN01 = inState( neigh (( qn01-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f101 = inState( neigh (( q101-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      fNN0 = inState( neigh (( qnn0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      fN10 = inState( neigh (( qn10-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f1N0 = inState( neigh (( q1n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f110 = inState( neigh (( q110-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      f000 = inState( neigh (( q000-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      ! element offset for auxField array
      elemOff = (iElem-1) * varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))


      ! square velocity and derived constants
      usq  = (u_x * u_x) + (u_y * u_y) + (u_z * u_z)
      usqn = div1_36 * (1._rk - 1.5_rk * usq) * rho

      ! read the relaxation parameter omega for the current level
      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      ! pre-calculate partial collision constants
      cmpl_o = 1._rk - omega

      ! f = (1-w) * f + w * fEq
      outState(( ielem-1)* nscalars+q000+( 1-1)* qq) &
        & = f000*cmpl_o+omega*rho*(div1_3-0.5_rk*usq)

      coeff_1 = div1_8 * omega * rho

      usqn_o1 = omega * usqn

      ui1 = u_x + u_y
      fac_1 = coeff_1 * ui1
      sum1_1 = fac_1 * div3_4h
      sum1_2 = fac_1 * ui1 + usqn_o1

      outState(( ielem-1)* nscalars+q110+( 1-1)* qq) &
        & = f110 * cmpl_o + sum1_1 + sum1_2
      outState(( ielem-1)* nscalars+qnn0+( 1-1)* qq) &
        & = fNN0 * cmpl_o - sum1_1 + sum1_2

      ui3 = -u_x + u_y
      fac_3 = coeff_1 * ui3
      sum3_1 = fac_3 * div3_4h
      sum3_2 = fac_3 * ui3 + usqn_o1

      outState(( ielem-1)* nscalars+qn10+( 1-1)* qq) &
        & = fN10 * cmpl_o + sum3_1 + sum3_2
      outState(( ielem-1)* nscalars+q1n0+( 1-1)* qq) &
        & = f1N0 * cmpl_o - sum3_1 + sum3_2

      ui10 = u_x + u_z
      fac_10 = coeff_1 * ui10
      sum10_1 = fac_10 * div3_4h
      sum10_2 = fac_10 * ui10 + usqn_o1

      outState(( ielem-1)* nscalars+q101+( 1-1)* qq) &
        & = f101 * cmpl_o + sum10_1 + sum10_2
      outState(( ielem-1)* nscalars+qn0n+( 1-1)* qq) &
        & = fN0N * cmpl_o - sum10_1 + sum10_2

      ui12 = -u_x + u_z
      fac_12 = coeff_1 * ui12
      sum12_1 = fac_12 * div3_4h
      sum12_2 = fac_12 * ui12 + usqn_o1

      outState(( ielem-1)* nscalars+qn01+( 1-1)* qq) &
        & = fN01 * cmpl_o + sum12_1 + sum12_2
      outState(( ielem-1)* nscalars+q10n+( 1-1)* qq) &
        & = f10N * cmpl_o - sum12_1 + sum12_2

      ui11 = u_y + u_z
      fac_11 = coeff_1 * ui11
      sum11_1 = fac_11 * div3_4h
      sum11_2 = fac_11 * ui11 + usqn_o1

      outState(( ielem-1)* nscalars+q011+( 1-1)* qq) &
        & = f011 * cmpl_o + sum11_1 + sum11_2
      outState(( ielem-1)* nscalars+q0nn+( 1-1)* qq) &
        & = f0NN * cmpl_o - sum11_1 + sum11_2

      ui13 = -u_y + u_z
      fac_13 = coeff_1 * ui13
      sum13_1 = fac_13 * div3_4h
      sum13_2 = fac_13 * ui13 + usqn_o1

      outState(( ielem-1)* nscalars+q0n1+( 1-1)* qq) &
        & = f0N1 * cmpl_o + sum13_1 + sum13_2
      outState(( ielem-1)* nscalars+q01n+( 1-1)* qq) &
        & = f01N * cmpl_o - sum13_1 + sum13_2

      omega_2 = 2._rk * omega
      coeff_2 = div1_8 * omega_2 * rho
      usqn_o2 = omega_2 * usqn

      fac_2 = coeff_2 * u_y
      sum2_1 = fac_2 * div3_4h
      sum2_2 = fac_2 * u_y + usqn_o2

      outState(( ielem-1)* nscalars+q010+( 1-1)* qq) &
        & = f010 * cmpl_o + sum2_1 + sum2_2
      outState(( ielem-1)* nscalars+q0n0+( 1-1)* qq) &
        & = f0N0 * cmpl_o - sum2_1 + sum2_2

      fac_4 = coeff_2 * u_x
      sum4_1 = fac_4 * div3_4h
      sum4_2 = fac_4 * u_x + usqn_o2

      outState(( ielem-1)* nscalars+qn00+( 1-1)* qq) &
        & = fN00 * cmpl_o - sum4_1 + sum4_2
      outState(( ielem-1)* nscalars+q100+( 1-1)* qq) &
        & = f100 * cmpl_o + sum4_1 + sum4_2

      fac_9 = coeff_2 * u_z
      sum9_1 = fac_9 * div3_4h
      sum9_2 = fac_9 * u_z + usqn_o2

      outState(( ielem-1)* nscalars+q001+( 1-1)* qq) &
        & = f001 * cmpl_o + sum9_1 + sum9_2
      outState(( ielem-1)* nscalars+q00n+( 1-1)* qq) &
        & = f00N * cmpl_o - sum9_1 + sum9_2


    enddo nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluid_rBGK_vStd_lD3Q19
! **************************************************************************** !


! **************************************************************************** !
  !> Advection relaxation routine for the D3Q19 model with BGK for
  !! incompressible lbm model
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluidIncomp_rBGK_vStd_lD3Q19( fieldProp, inState,    &
    &                                                  outState, auxField,    &
    &                                                  neigh, nElems, nSolve, &
    &                                                  level, layout, params, &
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
    integer :: iElem
    real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
      &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
      &              f000
    real(kind=rk) :: rho, u_x, u_y, u_z, usq
    real(kind=rk) :: usqn_o1, usqn_o2
    real(kind=rk) :: cmpl_o, omega
    real(kind=rk) :: coeff_1, coeff_2
    real(kind=rk) :: ui1, ui3, ui10, ui11, ui12, ui13
    real(kind=rk) :: fac_1, fac_2, fac_3, fac_4, fac_9, &
      &              fac_10, fac_11, fac_12, fac_13
    real(kind=rk) :: sum1_1, sum1_2, sum2_1, sum2_2, sum3_1, sum3_2, sum4_1, &
      &              sum4_2, sum9_1, sum9_2, sum10_1, sum10_2, sum11_1,      &
      &              sum11_2, sum12_1, sum12_2, sum13_1, sum13_2
    integer :: dens_pos, vel_pos(3), elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

!$omp do schedule(static)

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem=1,nSolve
      ! First load all local values into temp array
      fN00 = inState( neigh((qn00-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f0N0 = inState( neigh((q0n0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f00N = inState( neigh((q00n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f100 = inState( neigh((q100-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f010 = inState( neigh((q010-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f001 = inState( neigh((q001-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f0NN = inState( neigh((q0nn-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f0N1 = inState( neigh((q0n1-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f01N = inState( neigh((q01n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f011 = inState( neigh((q011-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fN0N = inState( neigh((qn0n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f10N = inState( neigh((q10n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fN01 = inState( neigh((qn01-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f101 = inState( neigh((q101-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fNN0 = inState( neigh((qnn0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fN10 = inState( neigh((qn10-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f1N0 = inState( neigh((q1n0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f110 = inState( neigh((q110-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f000 = inState( neigh((q000-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)

      ! element offset for auxField array
      elemOff = (iElem-1) * varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      ! square velocity and derived constants
      usq = u_x * u_x + u_y * u_y + u_z * u_z

      ! read the relaxation parameter omega for the current level
      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      ! pre-calculate partial collision constants
      cmpl_o = 1._rk - omega

      ! usqn = div1_36 * (rho - 1.5_rk * usq * rho0 )
      usqn_o1 = omega * div1_36 * ( rho - 1.5_rk * usq )

      outState(( ielem-1)* qq+q000+( 1-1)* qq) &
        & = f000 * cmpl_o + 12_rk * usqn_o1

      coeff_1 = div1_8 * omega

      ui1 = u_x + u_y
      fac_1 = coeff_1 * ui1
      sum1_1 = fac_1 * div3_4h
      sum1_2 = fac_1 * ui1 + usqn_o1

      outState(( ielem-1)* qq+q110+( 1-1)* qq) &
        & = f110 * cmpl_o + sum1_1 + sum1_2
      outState(( ielem-1)* qq+qnn0+( 1-1)* qq) &
        & = fNN0 * cmpl_o - sum1_1 + sum1_2

      ui3 = -u_x + u_y
      fac_3 = coeff_1 * ui3
      sum3_1 = fac_3 * div3_4h
      sum3_2 = fac_3 * ui3 + usqn_o1

      outState(( ielem-1)* qq+qn10+( 1-1)* qq) &
        & = fN10 * cmpl_o + sum3_1 + sum3_2
      outState(( ielem-1)* qq+q1n0+( 1-1)* qq) &
        & = f1N0 * cmpl_o - sum3_1 + sum3_2

      ui10 = u_x + u_z
      fac_10 = coeff_1 * ui10
      sum10_1 = fac_10 * div3_4h
      sum10_2 = fac_10 * ui10 + usqn_o1

      outState(( ielem-1)* qq+q101+( 1-1)* qq) &
        & = f101 * cmpl_o + sum10_1 + sum10_2
      outState(( ielem-1)* qq+qn0n+( 1-1)* qq) &
        & = fN0N * cmpl_o - sum10_1 + sum10_2

      ui12 = -u_x + u_z
      fac_12 = coeff_1 * ui12
      sum12_1 = fac_12 * div3_4h
      sum12_2 = fac_12 * ui12 + usqn_o1

      outState(( ielem-1)* qq+qn01+( 1-1)* qq) &
        & = fN01 * cmpl_o + sum12_1 + sum12_2
      outState(( ielem-1)* qq+q10n+( 1-1)* qq) &
        & = f10N * cmpl_o - sum12_1 + sum12_2

      ui11 = u_y + u_z
      fac_11 = coeff_1 * ui11
      sum11_1 = fac_11 * div3_4h
      sum11_2 = fac_11 * ui11 + usqn_o1

      outState(( ielem-1)* qq+q011+( 1-1)* qq) &
        & = f011 * cmpl_o + sum11_1 + sum11_2
      outState(( ielem-1)* qq+q0nn+( 1-1)* qq) &
        & = f0NN * cmpl_o - sum11_1 + sum11_2

      ui13 = -u_y + u_z
      fac_13 = coeff_1 * ui13
      sum13_1 = fac_13 * div3_4h
      sum13_2 = fac_13 * ui13 + usqn_o1

      outState(( ielem-1)* qq+q0n1+( 1-1)* qq) &
        & = f0N1 * cmpl_o + sum13_1 + sum13_2
      outState(( ielem-1)* qq+q01n+( 1-1)* qq) &
        & = f01N * cmpl_o - sum13_1 + sum13_2

      coeff_2 = div1_8 * omega * 2.0_rk
      usqn_o2 = 2_rk * usqn_o1

      fac_2 = coeff_2 * u_y
      sum2_1 = fac_2 * div3_4h
      sum2_2 = fac_2 * u_y + usqn_o2

      outState(( ielem-1)* qq+q010+( 1-1)* qq) &
        & = f010 * cmpl_o + sum2_1 + sum2_2
      outState(( ielem-1)* qq+q0n0+( 1-1)* qq) &
        & = f0N0 * cmpl_o - sum2_1 + sum2_2

      fac_4 = coeff_2 * u_x
      sum4_1 = fac_4 * div3_4h
      sum4_2 = fac_4 * u_x + usqn_o2

      outState(( ielem-1)* qq+qn00+( 1-1)* qq) &
        & = fN00 * cmpl_o - sum4_1 + sum4_2
      outState(( ielem-1)* qq+q100+( 1-1)* qq) &
        & = f100 * cmpl_o + sum4_1 + sum4_2

      fac_9 = coeff_2 * u_z
      sum9_1 = fac_9 * div3_4h
      sum9_2 = fac_9 * u_z + usqn_o2

      outState(( ielem-1)* qq+q001+( 1-1)* qq) &
        & = f001 * cmpl_o + sum9_1 + sum9_2
      outState(( ielem-1)* qq+q00n+( 1-1)* qq) &
        & = f00N * cmpl_o - sum9_1 + sum9_2

    enddo nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluidIncomp_rBGK_vStd_lD3Q19
! **************************************************************************** !


! **************************************************************************** !
  !> Advection relaxation routine for the D3Q19 model with TRT collision
  !! operator
  !! In TRT, there are two relaxation parameters one can choose.
  !! They have a relationship, which is so-called magic number:
  !! Lambda = ( 1/omegaP - 1/2 ) * ( 1/omegaN - 1/2 )
  !! Different value of Lambda results different error:
  !! Lambda = 1/4 is the best stability for the LBE. As well, this number gives
  !! the solution for the steady-state case dependant only on the equilibirium
  !! funciton.
  !! Lambda = 1/12 removes the third-order advection error
  !! Lambda = 1/6 removes fourth-order diffusion errors
  !! Lambda = 3/16 gives exact location of bounce-back walls for the Poiseuille
  !! flow.
  !! omegaP is usually fixed by viscosity, another one is fixed through the
  !! above magic number combination.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rTRT_vStd_lD3Q19( fieldProp, inState, outState, &
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
    real(kind=rk), intent(in) :: inState(nElems * varSys%nScalars)
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
    real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
      &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
      &              f000
    real(kind=rk) :: rho
    real(kind=rk) :: u_x, u_y, u_z, usq
    real(kind=rk) :: omega, omega_h, asym_omega, asym_omega_h
    real(kind=rk) :: ui
    real(kind=rk) :: asym, sym, feq_common
    real(kind=rk) :: t1x2,t2x2,fac1,fac2
    real(kind=rk), parameter :: t1x2_0 = 1._rk/18._rk * 2._rk
    real(kind=rk), parameter :: t2x2_0 = 1._rk/36._rk * 2._rk
    integer :: dens_pos, vel_pos(3), elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

!$omp do schedule(static)

!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem=1,nSolve
      ! First load all local values into temp array
      fN00 = inState( neigh((qn00-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f0N0 = inState( neigh((q0n0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f00N = inState( neigh((q00n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f100 = inState( neigh((q100-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f010 = inState( neigh((q010-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f001 = inState( neigh((q001-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f0NN = inState( neigh((q0nn-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f0N1 = inState( neigh((q0n1-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f01N = inState( neigh((q01n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f011 = inState( neigh((q011-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fN0N = inState( neigh((qn0n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f10N = inState( neigh((q10n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fN01 = inState( neigh((qn01-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f101 = inState( neigh((q101-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fNN0 = inState( neigh((qnn0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fN10 = inState( neigh((qn10-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f1N0 = inState( neigh((q1n0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f110 = inState( neigh((q110-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f000 = inState( neigh((q000-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)

      ! element offset for auxField array
      elemOff = (iElem-1) * varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      ! square velocity and derived constants
      usq  = (u_x * u_x) + (u_y * u_y) + (u_z * u_z)
      feq_common = 1._rk - 1.5_rk * usq

      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      omega_h = 0.5_rk * omega ! half omega

      asym_omega = 1.0_rk / (0.5_rk + fieldProp(1)%fluid%lambda &
        &                     / (1.0_rk / omega - 0.5_rk)       )
      asym_omega_h = 0.5_rk * asym_omega  ! half asymmetric omega

      ! let's start the relaxation process
      outstate(( ielem-1)* qq+q000+( 1-1)* qq)          &
        & = f000 * (1.0_rk - omega) + omega * div1_3 * rho * feq_common

      t2x2 = t2x2_0 * rho
      fac2 = t2x2 * t2cs4inv !inv2csq2

      ui = u_x + u_y
      sym = omega_h * (f110 + fNN0 - fac2 * ui * ui - t2x2 * feq_common)
      asym = asym_omega_h * (f110 - fNN0 - 3._rk * t2x2 * ui)
      outstate(( ielem-1)* qq+q110+( 1-1)* qq) &
        & = f110 - sym - asym
      outstate(( ielem-1)* qq+qnn0+( 1-1)* qq) &
        & = fNN0 - sym + asym

      ui = u_x - u_y
      sym = omega_h * (f1N0 + fN10 - fac2 * ui * ui - t2x2 * feq_common)
      asym = asym_omega_h * (f1N0 - fN10 - 3._rk * t2x2 * ui)
      outstate(( ielem-1)* qq+q1n0+( 1-1)* qq) &
        & = f1N0 - sym - asym
      outstate(( ielem-1)* qq+qn10+( 1-1)* qq) &
        & = fN10 - sym + asym

      ui = u_x + u_z
      sym = omega_h * (f101 + fN0N - fac2 * ui * ui - t2x2 * feq_common)
      asym = asym_omega_h * (f101 - fN0N - 3._rk * t2x2 * ui)
      outstate(( ielem-1)* qq+q101+( 1-1)* qq) &
        & = f101 - sym - asym
      outstate(( ielem-1)* qq+qn0n+( 1-1)* qq) &
        & = fN0N - sym + asym

      ui = u_x - u_z
      sym = omega_h * (f10N + fN01 - fac2 * ui * ui - t2x2 * feq_common)
      asym = asym_omega_h * (f10N - fN01 - 3._rk * t2x2 * ui)
      outstate(( ielem-1)* qq+q10n+( 1-1)* qq) &
        & = f10N - sym - asym
      outstate(( ielem-1)* qq+qn01+( 1-1)* qq) &
        & = fN01 - sym + asym

      ui = u_y + u_z
      sym = omega_h * (f011 + f0NN - fac2 * ui * ui - t2x2 * feq_common)
      asym = asym_omega_h * (f011 - f0NN - 3._rk * t2x2 * ui)
      outstate(( ielem-1)* qq+q011+( 1-1)* qq) &
        & = f011 - sym - asym
      outstate(( ielem-1)* qq+q0nn+( 1-1)* qq) &
        & = f0NN - sym + asym

      ui = u_y - u_z
      sym = omega_h * (f01N + f0N1 - fac2 * ui * ui - t2x2 * feq_common)
      asym = asym_omega_h * (f01N - f0N1 - 3._rk * t2x2 * ui)
      outstate(( ielem-1)* qq+q01n+( 1-1)* qq) &
        & = f01N - sym - asym
      outstate(( ielem-1)* qq+q0n1+( 1-1)* qq) &
        & = f0N1 - sym + asym

      t1x2 = t1x2_0 * rho
      fac1 = t1x2 * t2cs4inv !inv2csq2

      sym = omega_h * (f100 + fN00 - fac1 * u_x * u_x - t1x2 * feq_common)
      asym = asym_omega_h * (f100 - fN00 - 3._rk * t1x2 * u_x)
      outstate(( ielem-1)* qq+q100+( 1-1)* qq) &
        & = f100 - sym - asym
      outstate(( ielem-1)* qq+qn00+( 1-1)* qq) &
        & = fN00 - sym + asym

      sym = omega_h * (f010 + f0N0 - fac1 * u_y * u_y - t1x2 * feq_common)
      asym = asym_omega_h * (f010 - f0N0 - 3._rk * t1x2 * u_y)
      outstate(( ielem-1)* qq+q010+( 1-1)* qq) &
        & = f010 - sym - asym
      outstate(( ielem-1)* qq+q0n0+( 1-1)* qq) &
        & = f0N0 - sym + asym

      sym = omega_h * (f001 + f00N - fac1 * u_z * u_z - t1x2 * feq_common)
      asym = asym_omega_h * (f001 - f00N - 3._rk*t1x2*u_z)
      outstate(( ielem-1)* qq+q001+( 1-1)* qq) &
        & = f001 - sym - asym
      outstate(( ielem-1)* qq+q00n+( 1-1)* qq) &
        & = f00N - sym + asym
    enddo nodeloop
!$omp end do


  end subroutine mus_advRel_kFluid_rTRT_vStd_lD3Q19
! **************************************************************************** !


! **************************************************************************** !
  !> Advection relaxation routine for the D3Q19 model with TRT collision
  !! operator.
  !!
  !! Incompressible model
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluidIncomp_rTRT_vStd_lD3Q19( fieldProp, inState,    &
    &                                                  outState, auxField,    &
    &                                                  neigh, nElems, nSolve, &
    &                                                  level, layout, params, &
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
    real(kind=rk), intent(in) :: inState(nElems * varSys%nScalars)
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
    integer :: iElem ! element counter
    real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
      &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
      &              f000
    real(kind=rk) :: rho, u_x, u_y, u_z, usq
    real(kind=rk) :: omega, omega_h, asym_omega, asym_omega_h
    real(kind=rk) :: ui
    real(kind=rk) :: asym, sym, feq_common, t1_feq, t2_feq
    real(kind=rk), parameter :: t1x2 = 1._rk/ 9._rk
    real(kind=rk), parameter :: t2x2 = 1._rk/18._rk
    real(kind=rk), parameter :: fac1 = t1x2*t2cs4inv !inv2csq2
    real(kind=rk), parameter :: fac2 = t2x2*t2cs4inv !inv2csq2
    integer :: dens_pos, vel_pos(3), elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

!$omp do schedule(static)

!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1, nSolve
      ! First load all local values into temp array
      fN00 = inState( neigh((qn00-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f0N0 = inState( neigh((q0n0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f00N = inState( neigh((q00n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f100 = inState( neigh((q100-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f010 = inState( neigh((q010-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f001 = inState( neigh((q001-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f0NN = inState( neigh((q0nn-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f0N1 = inState( neigh((q0n1-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f01N = inState( neigh((q01n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f011 = inState( neigh((q011-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fN0N = inState( neigh((qn0n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f10N = inState( neigh((q10n-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fN01 = inState( neigh((qn01-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f101 = inState( neigh((q101-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fNN0 = inState( neigh((qnn0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      fN10 = inState( neigh((qn10-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f1N0 = inState( neigh((q1n0-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f110 = inState( neigh((q110-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      f000 = inState( neigh((q000-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)

      ! element offset for auxField array
      elemOff = (iElem-1) * varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)
      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      ! square velocity and derived constants
      usq  = (u_x * u_x) + (u_y * u_y) + (u_z * u_z)
      feq_common = rho - 1.5_rk * usq

      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      omega_h = 0.5_rk * omega ! half omega

      asym_omega = 1.0_rk / (0.5_rk + fieldProp(1)%fluid%lambda &
        &                     / (1.0_rk/omega - 0.5_rk)         )
      asym_omega_h = 0.5_rk * asym_omega  ! half asymmetric omega

      ! let's start the relaxation process
      outstate(( ielem-1)* qq+q000+( 1-1)* qq) &
        & = f000 * (1._rk - omega) + omega * div1_3 * feq_common

      ! t2x2 = 1._rk/18._rk
      t2_feq = t2x2 * feq_common

      ui = u_x + u_y
      sym = omega_h * (f110 + fNN0 - fac2 * ui * ui - t2_feq)
      asym = asym_omega_h * (f110 - fNN0 - div1_6 * ui)
      outstate(( ielem-1)* qq+q110+( 1-1)* qq) &
        & = f110 - sym - asym
      outstate(( ielem-1)* qq+qnn0+( 1-1)* qq) &
        & = fNN0 - sym + asym

      ui = u_x - u_y
      sym = omega_h * (f1N0 + fN10 - fac2 * ui * ui - t2_feq)
      asym = asym_omega_h * (f1N0 - fN10 - div1_6 * ui)
      outstate(( ielem-1)* qq+q1n0+( 1-1)* qq) &
        & = f1N0 - sym - asym
      outstate(( ielem-1)* qq+qn10+( 1-1)* qq) &
        & = fN10 - sym + asym

      ui = u_x + u_z
      sym = omega_h * (f101 + fN0N - fac2 * ui * ui - t2_feq)
      asym = asym_omega_h * (f101 - fN0N - div1_6 * ui)
      outstate(( ielem-1)* qq+q101+( 1-1)* qq) &
        & = f101 - sym - asym
      outstate(( ielem-1)* qq+qn0n+( 1-1)* qq) &
        & = fN0N - sym + asym

      ui = u_x - u_z
      sym = omega_h * (f10N + fN01 - fac2 * ui * ui-t2_feq)
      asym = asym_omega_h * (f10N - fN01 - div1_6 * ui)
      outstate(( ielem-1)* qq+q10n+( 1-1)* qq) &
        & = f10N - sym - asym
      outstate(( ielem-1)* qq+qn01+( 1-1)* qq) &
        & = fN01 - sym + asym

      ui = u_y + u_z
      sym = omega_h * (f011 + f0NN - fac2 * ui * ui - t2_feq)
      asym = asym_omega_h * (f011 - f0NN - div1_6 * ui)
      outstate(( ielem-1)* qq+q011+( 1-1)* qq) &
        & = f011 - sym - asym
      outstate(( ielem-1)* qq+q0nn+( 1-1)* qq) &
        & = f0NN - sym + asym

      ui = u_y - u_z
      sym = omega_h * (f01N + f0N1 - fac2 * ui * ui - t2_feq)
      asym = asym_omega_h * (f01N - f0N1 - div1_6 * ui)
      outstate(( ielem-1)* qq+q01n+( 1-1)* qq) &
        & = f01N - sym - asym
      outstate(( ielem-1)* qq+q0n1+( 1-1)* qq) &
        & = f0N1 - sym + asym

      ! t1x2 = 1._rk/ 9._rk
      t1_feq = t1x2*feq_common

      ! ui   = u_y
      sym = omega_h * (f010 + f0N0 - fac1 * u_y * u_y - t1_feq)
      asym = asym_omega_h * (f010 - f0N0 - div1_3 * u_y)
      outstate(( ielem-1)* qq+q010+( 1-1)* qq) &
        & = f010 - sym - asym
      outstate(( ielem-1)* qq+q0n0+( 1-1)* qq) &
        & = f0N0 - sym + asym

      ! ui   = u_x
      sym = omega_h * (f100 + fN00 - fac1 * u_x * u_x - t1_feq)
      asym = asym_omega_h * (f100 - fN00 - div1_3 * u_x)
      outstate(( ielem-1)* qq+q100+( 1-1)* qq) &
        & = f100 - sym - asym
      outstate(( ielem-1)* qq+qn00+( 1-1)* qq) &
        & = fN00 - sym + asym

      ! ui   = u_z
      sym = omega_h * (f001 + f00N - fac1 * u_z * u_z - t1_feq)
      asym = asym_omega_h * (f001 - f00N - div1_3 * u_z)
      outstate(( ielem-1)* qq+q001+( 1-1)* qq) &
        & = f001 - sym - asym
      outstate(( ielem-1)* qq+q00n+( 1-1)* qq) &
        & = f00N - sym + asym
    enddo nodeloop

  end subroutine mus_advRel_kFluidIncomp_rTRT_vStd_lD3Q19
! **************************************************************************** !


! **************************************************************************** !
  !> No comment yet!
  !!
  !! TODO add comment
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rBGK_vBlock_lD3Q19( fieldProp, inState,         &
    &                                              outState, auxField, neigh,  &
    &                                              nElems, nSolve, level,      &
    &                                              layout, params, varSys,     &
    &                                              derVarPos )
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
    integer :: iElem, nScalars, minElem, maxElem
    real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
      &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
      &              f000
    real(kind=rk) :: rho     !< local density
    real(kind=rk) :: u_x     !< local x-velocity
    real(kind=rk) :: u_y     !< local y-velocity
    real(kind=rk) :: u_z     !< local z-velocity
    real(kind=rk) :: usq     !< square velocity
    ! derived constants
    real(kind=rk) :: usqn, usqn_o1, usqn_o2
    real(kind=rk) :: cmpl_o, omega
    real(kind=rk) :: coeff_1, coeff_2
    real(kind=rk) :: ui1, ui3, ui10, ui11, ui12, ui13
    real(kind=rk) :: fac_1, fac_2, fac_3, fac_4, fac_9, fac_10, fac_11, fac_12,&
      &              fac_13
    real(kind=rk) :: sum1_1, sum1_2, sum2_1, sum2_2, sum3_1, sum3_2, sum4_1,   &
      &              sum4_2, sum9_1, sum9_2, sum10_1, sum10_2, sum11_1,        &
      &              sum11_2, sum12_1, sum12_2, sum13_1, sum13_2
    integer :: dens_pos, vel_pos(3), elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

!$omp do schedule(static)
    do minElem = 1, nSolve, vlen

      maxElem = min( minElem + vlen - 1, nSolve )

      !NEC$ ivdep
      !NEC$ shortloop
      do iElem = minElem, maxElem
        ! First load all local values into temp array
        fN00 = inState( neigh((qn00-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f0N0 = inState( neigh((q0n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f00N = inState( neigh((q00n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f100 = inState( neigh((q100-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f010 = inState( neigh((q010-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f001 = inState( neigh((q001-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f0NN = inState( neigh((q0nn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f0N1 = inState( neigh((q0n1-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f01N = inState( neigh((q01n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f011 = inState( neigh((q011-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        fN0N = inState( neigh((qn0n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f10N = inState( neigh((q10n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        fN01 = inState( neigh((qn01-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f101 = inState( neigh((q101-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        fNN0 = inState( neigh((qnn0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        fN10 = inState( neigh((qn10-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f1N0 = inState( neigh((q1n0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f110 = inState( neigh((q110-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
        f000 = inState( neigh((q000-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

        ! element offset for auxField array
        elemOff = (iElem-1) * varSys%nAuxScalars
        ! local density
        rho = auxField(elemOff + dens_pos)
        ! local x-, y- and z-velocity
        u_x = auxField(elemOff + vel_pos(1))
        u_y = auxField(elemOff + vel_pos(2))
        u_z = auxField(elemOff + vel_pos(3))

        ! square velocity and derived constants
        usq = (u_x * u_x) + (u_y * u_y) + (u_z * u_z)
        usqn = div1_36 * (1._rk - 1.5_rk * usq) * rho

        ! read the relaxation parameter omega for the current level
        omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
        ! pre-calculate partial collision constants
        cmpl_o = 1._rk - omega

        ! f = (1-w) * f + w * fEq
        outState(( ielem-1)* nscalars+q000+( 1-1)* qq) &
          & = f000 * cmpl_o + omega * rho * (div1_3 - 0.5_rk * usq)

        coeff_1 = div1_8 * omega * rho

        usqn_o1 = omega * usqn

        ui1 = u_x + u_y
        fac_1 = coeff_1 * ui1
        sum1_1 = fac_1 * div3_4h
        sum1_2 = fac_1 * ui1 + usqn_o1

        outState(( ielem-1)* nscalars+q110+( 1-1)* qq) &
          & = f110 * cmpl_o + sum1_1 + sum1_2
        outState(( ielem-1)* nscalars+qnn0+( 1-1)* qq) &
          & = fNN0 * cmpl_o - sum1_1 + sum1_2

        ui3 = -u_x + u_y
        fac_3 = coeff_1 * ui3
        sum3_1 = fac_3 * div3_4h
        sum3_2 = fac_3 * ui3 + usqn_o1

        outState(( ielem-1)* nscalars+qn10+( 1-1)* qq) &
          & = fN10 * cmpl_o + sum3_1 + sum3_2
        outState(( ielem-1)* nscalars+q1n0+( 1-1)* qq) &
          & = f1N0 * cmpl_o - sum3_1 + sum3_2

        ui10 =  u_x + u_z
        fac_10 = coeff_1 * ui10
        sum10_1 = fac_10 * div3_4h
        sum10_2 = fac_10 * ui10 + usqn_o1

        outState(( ielem-1)* nscalars+q101+( 1-1)* qq) &
          & = f101 * cmpl_o + sum10_1 + sum10_2
        outState(( ielem-1)* nscalars+qn0n+( 1-1)* qq) &
          & = fN0N * cmpl_o - sum10_1 + sum10_2

        ui12 = -u_x + u_z
        fac_12 = coeff_1 * ui12
        sum12_1 = fac_12 * div3_4h
        sum12_2 = fac_12 * ui12 + usqn_o1

        outState(( ielem-1)* nscalars+qn01+( 1-1)* qq) &
          & = fN01 * cmpl_o + sum12_1 + sum12_2
        outState(( ielem-1)* nscalars+q10n+( 1-1)* qq) &
          & = f10N * cmpl_o - sum12_1 + sum12_2

        ui11 =  u_y + u_z
        fac_11 = coeff_1 * ui11
        sum11_1 = fac_11 * div3_4h
        sum11_2 = fac_11 * ui11 + usqn_o1

        outState(( ielem-1)* nscalars+q011+( 1-1)* qq) &
          & = f011 * cmpl_o + sum11_1 + sum11_2
        outState(( ielem-1)* nscalars+q0nn+( 1-1)* qq) &
          & = f0NN * cmpl_o - sum11_1 + sum11_2

        ui13 = -u_y + u_z
        fac_13 = coeff_1 * ui13
        sum13_1 = fac_13 * div3_4h
        sum13_2 = fac_13 * ui13 + usqn_o1

        outState(( ielem-1)* nscalars+q0n1+( 1-1)* qq) &
          & = f0N1 * cmpl_o + sum13_1 + sum13_2
        outState(( ielem-1)* nscalars+q01n+( 1-1)* qq) &
          & = f01N * cmpl_o - sum13_1 + sum13_2

        coeff_2 = div1_8 * omega * 2.0_rk * rho
        usqn_o2 = omega * 2.0_rk * usqn

        fac_2 = coeff_2 * u_y
        sum2_1 = fac_2 * div3_4h
        sum2_2 = fac_2 * u_y + usqn_o2

        outState(( ielem-1)* nscalars+q010+( 1-1)* qq) &
          & = f010 * cmpl_o + sum2_1 + sum2_2
        outState(( ielem-1)* nscalars+q0n0+( 1-1)* qq) &
          & = f0N0 * cmpl_o - sum2_1 + sum2_2

        fac_4 = coeff_2 * u_x
        sum4_1 = fac_4 * div3_4h
        sum4_2 = fac_4 * u_x + usqn_o2

        outState(( ielem-1)* nscalars+qn00+( 1-1)* qq) &
          & = fN00 * cmpl_o - sum4_1 + sum4_2
        outState(( ielem-1)* nscalars+q100+( 1-1)* qq) &
          & = f100 * cmpl_o + sum4_1 + sum4_2

        fac_9 = coeff_2 * u_z
        sum9_1 = fac_9 * div3_4h
        sum9_2 = fac_9 * u_z + usqn_o2

        outState(( ielem-1)* nscalars+q001+( 1-1)* qq) &
          & = f001 * cmpl_o + sum9_1 + sum9_2
        outState(( ielem-1)* nscalars+q00n+( 1-1)* qq) &
          & = f00N * cmpl_o - sum9_1 + sum9_2

      end do
    end do
!$omp end do nowait

  end subroutine mus_advRel_kFluid_rBGK_vBlock_lD3Q19
! **************************************************************************** !


! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
pure subroutine f_f_eq_regularized_2nd_ord_d3q19 ( weight, rho, u_x, u_y, u_z, feq, &
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
      f01 = -cs2inv*u_x
      f02 = div1_6*cs4inv*(2._rk*u_x_sqr - (u_y_sqr + u_z_sqr))
      feq(1) = weight(1) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*(2._rk*a12xx - (a12yy + a12zz))
      f1(1) = weight(1) * f12

      !iDir = 4
      f01 = -f01
      feq(4) = weight(4) * rho * (f00 + f01 + f02)
      f1(4) = weight(4) * f12

      !iDir = 2
      f01 = -cs2inv*u_y
      f02 = div1_6*cs4inv*(-(u_x_sqr + u_z_sqr) + 2._rk*u_y_sqr)
      feq(2) = weight(2) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*(-(a12xx + a12zz) + 2._rk*a12yy)
      f1(2) = weight(2) * f12

      !iDir = 5
      f01 = -f01
      feq(5) = weight(5) * rho * (f00 + f01 + f02)
      f1(5) = weight(5) * f12

      !iDir = 3
      f01 = -cs2inv*u_z
      f02 = div1_6*cs4inv*(-(u_x_sqr + u_y_sqr) + 2._rk*u_z_sqr)
      feq(3) = weight(3) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*(-(a12xx + a12yy) + 2._rk*a12zz)
      f1(3) = weight(3) * f12

      !iDir = 6
      f01 = -f01
      feq(6) = weight(6) * rho * (f00 + f01 + f02)
      f1(6) = weight(6) * f12

      !iDir = 7
      f01 = cs2inv*(-u_y - u_z)
      f02 = f02 + 0.5_rk*cs4inv*(u_y_sqr + 2._rk*u_y_u_z)
      feq(7) = weight(7) * rho * (f00 + f01 + f02)
      f12 = f12 + 0.5_rk*cs4inv*(a12yy + 2.0_rk*a12yz)
      f1(7) = weight(7) * f12

      !iDir = 10
      f01 = -f01
      feq(10) = weight(10) * rho * (f00 + f01 + f02)
      f1(10) = weight(10) * f12

      !iDir = 8
      f01 = cs2inv*(-u_y + u_z)
      f02 = f02 - 2._rk*cs4inv*(u_y_u_z)
      feq(8) = weight(8) * rho * (f00 + f01 + f02)
      f12 = f12 - 2._rk*cs4inv*(a12yz)
      f1(8) = weight(8) * f12

      !iDir = 9
      f01 = -f01
      feq(9) = weight(9) * rho * (f00 + f01 + f02)
      f1(9) = weight(9) * f12

      !iDir = 11
      f01 = cs2inv*(-u_x - u_z)
      f02 = div1_6*cs4inv*(2._rk*(u_x_sqr + u_z_sqr) - u_y_sqr + 6._rk*u_x_u_z)
      feq(11) = weight(11) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*(2._rk*(a12xx + a12zz) - a12yy + 6.0_rk*a12xz)
      f1(11) = weight(11) * f12

      !iDir = 14
      f01 = -f01
      feq(14) = weight(14) * rho * (f00 + f01 + f02)
      f1(14) = weight(14) * f12

      !iDir = 12
      f01 = cs2inv*(u_x - u_z)
      f02 = f02 - 2._rk*cs4inv*(u_x_u_z)
      feq(12) = weight(12) * rho * (f00 + f01 + f02)
      f12 = f12 - 2._rk*cs4inv*(a12xz)
      f1(12) = weight(12) * f12

      !iDir = 13
      f01 = -f01
      feq(13) = weight(13) * rho * (f00 + f01 + f02)
      f1(13) = weight(13) * f12

      !iDir = 15
      f01 = cs2inv*(-u_x - u_y)
      f02 = div1_6*cs4inv*(2._rk*(u_x_sqr + u_y_sqr) - u_z_sqr + 6._rk*u_x_u_y)
      feq(15) = weight(15) * rho * (f00 + f01 + f02)
      f12 = div1_6*cs4inv*(2._rk*(a12xx + a12yy) - a12zz + 6.0_rk*a12xy)
      f1(15) = weight(15) * f12

      !iDir = 18
      f01 = -f01
      feq(18) = weight(18) * rho * (f00 + f01 + f02)
      f1(18) = weight(18) * f12

      !iDir = 16
      f01 = cs2inv*(-u_x + u_y)
      f02 = f02 - 2._rk*cs4inv*(u_x_u_y)
      feq(16) = weight(16) * rho * (f00 + f01 + f02)
      f12 = f12 - 2._rk*cs4inv*(a12xy)
      f1(16) = weight(16) * f12

      !iDir = 17
      f01 = -f01
      feq(17) = weight(17) * rho * (f00 + f01 + f02)
      f1(17) = weight(17) * f12

      !iDir = 19
      !f01 = 0._rk
      f02 = -div1_6*cs4inv*(u_x_sqr + u_y_sqr + u_z_sqr)
      feq(19) = weight(19) * rho * (f00 + f02)
      f12 = -div1_6*cs4inv*(a12xx + a12yy + a12zz)
      f1(19) = weight(19) * f12

  end subroutine f_f_eq_regularized_2nd_ord_d3q19
! ****************************************************************************** !

! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
  ! with correction from Guo et al., "An efficient lattice Boltzmann method for
  ! compressible aerodynamics on D3Q19 lattice", JCP, 2020
pure subroutine f_f_eq_regularized_4th_ord_d3q19 ( weight, rho, u_x, u_y, u_z, &
&  feq, f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )
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
  real(kind=rk) :: a13xyy, a13xxy, a13xxz, a13xzz, a13yzz, a13yyz, f13
  ! ---------------------------------------------------------------------------

    call f_f_eq_regularized_2nd_ord_d3q19 ( weight, rho, u_x, u_y, u_z, feq, f1, &
      &    a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

    u_x_sqr = u_x**2
    u_y_sqr = u_y**2
    u_z_sqr = u_z**2
    u_x_sqr_u_y = u_x_sqr * u_y
    u_x_sqr_u_z = u_x_sqr * u_z
    u_y_sqr_u_x = u_y_sqr * u_x
    u_y_sqr_u_z = u_y_sqr * u_z
    u_z_sqr_u_x = u_z_sqr * u_x
    u_z_sqr_u_y = u_z_sqr * u_y

    a13xxy = 2.0_rk * u_x * a12xy + u_y * a12xx
    a13xxz = 2.0_rk * u_x * a12xz + u_z * a12xx
    a13xyy = u_x * a12yy + 2.0_rk * u_y * a12xy
    a13xzz = u_x * a12zz + 2.0_rk * u_z * a12xz
    a13yzz = u_y * a12zz + 2.0_rk * u_z * a12yz
    a13yyz = 2.0_rk * u_y * a12yz + u_z * a12yy

    !iDir = 1
    f03 = div1_3*cs6inv*( (u_z_sqr_u_x + u_y_sqr_u_x) )
    feq(1) = feq(1) + weight(1) * rho * (f03)
    f13 = div1_3*cs6inv*( (a13xzz + a13xyy) )
    f1(1) = f1(1) + weight(1) * f13

    !iDir = 4
    f03 = -f03
    feq(4) = feq(4) + weight(4) * rho * (f03)
    f13 = -f13
    f1(4) = f1(4) + weight(4) * f13

    !iDir = 2
    f03 = div1_3*cs6inv*( (u_x_sqr_u_y + u_z_sqr_u_y) )
    feq(2) = feq(2) + weight(2) * rho * (f03)
    f13 = div1_3*cs6inv*( (a13xxy + a13yzz) )
    f1(2) = f1(2) + weight(2) * f13

    !iDir = 5
    f03 = -f03
    feq(5) = feq(5) + weight(5) * rho * (f03)
    f13 = -f13
    f1(5) = f1(5) + weight(5) * f13

    !iDir = 3
    f03 = div1_3*cs6inv*( (u_x_sqr_u_z + u_y_sqr_u_z))
    feq(3) = feq(3) + weight(3) * rho * (f03)
    f13 = div1_3*cs6inv*( (a13xxz + a13yyz) )
    f1(3) = f1(3) + weight(3) * f13

    !iDir = 6
    f03 = -f03
    feq(6) = feq(6) + weight(6) * rho * (f03)
    f13 = -f13
    f1(6) = f1(6) + weight(6) * f13

    !iDir = 7
    !f03 = div1_6*cs6inv*( -(u_x_sqr_u_y + u_z_sqr_u_y) + (u_x_sqr_u_y - u_z_sqr_u_y) &
    !  &      - (u_y_sqr_u_z + u_x_sqr_u_z) - (u_y_sqr_u_z - u_x_sqr_u_z) )
    f03 = div1_3*cs6inv*( - (u_z_sqr_u_y + u_y_sqr_u_z) )
    feq(7) = feq(7) + weight(7) * rho * (f03)
    f13 = div1_3*cs6inv*( - (a13yzz + a13yyz) )
    f1(7) = f1(7) + weight(7) * f13

    !iDir = 10
    f03 = -f03
    feq(10) = feq(10) + weight(10) * rho * (f03)
    f13 = -f13
    f1(10) = f1(10) + weight(10) * f13

    !iDir = 8
    !f03 = div1_6*cs6inv*( -(u_x_sqr_u_y + u_z_sqr_u_y) + (u_x_sqr_u_y - u_z_sqr_u_y) &
    !  &      + (u_y_sqr_u_z + u_x_sqr_u_z) + (u_y_sqr_u_z - u_x_sqr_u_z) )
    f03 = div1_3*cs6inv*( (-u_z_sqr_u_y + u_y_sqr_u_z) )
    feq(8) = feq(8) + weight(8) * rho * (f03)
    f13 = div1_3*cs6inv*( (-a13yzz + a13yyz) )
    f1(8) = f1(8) + weight(8) * f13

    !iDir = 9
    f03 = -f03
    feq(9) = feq(9) + weight(9) * rho * (f03)
    f13 = -f13
    f1(9) = f1(9) + weight(9) * f13

    !iDir = 11
    !f03 = div1_6*cs6inv*( -(u_z_sqr_u_x + u_y_sqr_u_x) - (u_z_sqr_u_x - u_y_sqr_u_x) &
    !  &      - (u_y_sqr_u_z + u_x_sqr_u_z) + (u_y_sqr_u_z - u_x_sqr_u_z) )
    f03 = div1_3*cs6inv*( -(u_z_sqr_u_x + u_x_sqr_u_z ) )
    feq(11) = feq(11) + weight(11) * rho * (f03)
    f13 = div1_3*cs6inv*( -(a13xxz + a13xzz) )
    f1(11) = f1(11) + weight(11) * f13

    !iDir = 14
    f03 = -f03
    feq(14) = feq(14) + weight(14) * rho * (f03)
    f13 = -f13
    f1(14) = f1(14) + weight(14) * f13

    !iDir = 12
    !f03 = div1_6*cs6inv*( (u_z_sqr_u_x + u_y_sqr_u_x) + (u_z_sqr_u_x - u_y_sqr_u_x) &
    !  &      - (u_y_sqr_u_z + u_x_sqr_u_z) + (u_y_sqr_u_z - u_x_sqr_u_z) )
    f03 = div1_3*cs6inv*( (u_z_sqr_u_x - u_x_sqr_u_z) )
    feq(12) = feq(12) + weight(12) * rho * (f03)
    f13 =  div1_3*cs6inv*( (-a13xxz + a13xzz) )
    f1(12) = f1(12) + weight(12) * ( f13 )

    !iDir = 13
    f03 = -f03
    feq(13) = feq(13) + weight(13) * rho * (f03)
    f13 = -f13
    f1(13) = f1(13) + weight(13) * f13

    !iDir = 15
    !f03 = div1_6*cs6inv*( -(u_x_sqr_u_y + u_z_sqr_u_y) - (u_x_sqr_u_y - u_z_sqr_u_y) &
    !  &      - (u_z_sqr_u_x + u_y_sqr_u_x) + (u_z_sqr_u_x - u_y_sqr_u_x) )
    f03 = div1_3*cs6inv*( -(u_x_sqr_u_y + u_y_sqr_u_x) )
    feq(15) = feq(15) + weight(15) * rho * (f03)
    f13 = div1_3*cs6inv*( -(a13xxy + a13xyy) )
    f1(15) = f1(15) + weight(15) * f13

    !iDir = 18
    f03 = -f03
    feq(18) = feq(18) + weight(18) * rho * (f03)
    f13 = -f13
    f1(18) = f1(18) + weight(18) * f13

    !iDir = 16
    !f03 = div1_6*cs6inv*( (u_x_sqr_u_y + u_z_sqr_u_y) + (u_x_sqr_u_y - u_z_sqr_u_y) &
    !  &      - (u_z_sqr_u_x + u_y_sqr_u_x) + (u_z_sqr_u_x - u_y_sqr_u_x) )
    f03 = div1_3*cs6inv*( (u_x_sqr_u_y - u_y_sqr_u_x) )
    feq(16) = feq(16) + weight(16) * rho * (f03)
    f13 = div1_3*cs6inv*( (a13xxy - a13xyy) )
    f1(16) = f1(16) + weight(16) * (f13)

    !iDir = 17
    f03 = -f03
    feq(17) = feq(17) + weight(17) * rho * (f03)
    f13 = -f13
    f1(17) = f1(17) + weight(17) * ( f13 )

    !iDir = 19
    !f03 = 0._rk
    !feq(19) = feq(19) + 0._rk
    !f13 = 0._rk
    !f1(19) = f1(19) + 0._rk
end subroutine f_f_eq_regularized_4th_ord_d3q19
! ****************************************************************************** !

! **************************************************************************** !
  !> Regularized relaxation routine for the D3Q19 and 27 model with BGK.
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005

  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  !! works for both d3q19 and d3q27
  subroutine bgk_Regularized_d3q19( fieldProp, inState, outState, auxField, &
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

      call f_f_eq_regularized_2nd_ord_d3q19( layout%weight(:), rho, u_x, u_y, u_z, feq, &
        &                                   f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                      + cmpl_o*f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_Regularized_d3q19
! **************************************************************************** !


! **************************************************************************** !
  !> Recursive Regularized relaxation routine for the D3Q19
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005

  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  !! works for both d3q19 and d3q27
  subroutine bgk_RecursiveRegularized_d3q19( fieldProp, inState, outState, auxField, &
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

      call f_f_eq_regularized_4th_ord_d3q19( layout%weight(:), rho, u_x, u_y, u_z, &
        &                        feq, f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                      + cmpl_o*f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_RecursiveRegularized_d3q19
! **************************************************************************** !



! **************************************************************************** !
  !> Projected Recursive Regularized relaxation routine for the D3Q19
  ! based on High-order extension of the recursive regularized lattice
  ! Boltzmann method, PhD Thesis, COREIXAS 2018
  ! comemnted terms are not well approximated by d3q19 discretization model

  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  !! works for both d3q19 and d3q27
  subroutine bgk_ProjectedRecursiveRegularized_d3q19( fieldProp, inState, outState, auxField, &
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
    real(kind=rk) :: f( QQ ), SR(6), gradU(3,3,1)!!TODO:, tr_SR
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

      call f_f_eq_regularized_4th_ord_d3q19( layout%weight(:), rho, u_x, u_y, u_z, &
        &                        feq, f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                      + cmpl_o*f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_ProjectedRecursiveRegularized_d3q19
! **************************************************************************** !

! **************************************************************************** !
  !> Projected Recursive Regularized relaxation routine for the D3Q19
  ! based on High-order extension of the recursive regularized lattice
  ! Boltzmann method, PhD Thesis, COREIXAS 2018
  ! comemnted terms are not well approximated by d3q19 discretization model

  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  !! works for both d3q19 and d3q27
  subroutine bgk_HybridRecursiveRegularized_d3q19( fieldProp, inState, outState, auxField, &
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
    real(kind=rk) :: omega, taup, cmpl_o, feq(QQ), f1(QQ)
    real(kind=rk) :: sigma
    integer :: denspos, velpos(3), elemOff, nScalars
    ! ---------------------------------------------------------------------------

!cdir nodep
!ibm* novector
!dir$ novector

    denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    velpos(1:3) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    nScalars = varSys%nScalars

    !call getHermitepolynomials_D3Q19( layout )

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

      call f_f_eq_regularized_4th_ord_d3q19( layout%weight(:), rho, u_x, u_y, u_z, &
        &                        feq, f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
          &                                                      + cmpl_o*f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_HybridRecursiveRegularized_d3q19
! **************************************************************************** !

! **************************************************************************** !
  !> Projected Recursive Regularized relaxation routine for the D3Q19
  ! based on High-order extension of the recursive regularized lattice
  ! Boltzmann method, PhD Thesis, COREIXAS 2018
  ! commented terms are not well approximated by d3q19 discretization model

  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  !! works for both d3q19 and d3q27
  subroutine bgk_HybridRecursiveRegularizedCorr_d3q19( fieldProp, inState, outState, auxField, &
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
    real(kind=rk) :: omega, taup, cmpl_o, feq(QQ), f1(QQ)
    real(kind=rk) :: sigma, gradRhoU3(3,1), S_Corr(QQ), gradRhoUVZ(3,1)
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
        gradU(:,:,1:1) = scheme%Grad%U_ptr(      &
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

        ! 1 = x, 2 = y, 3 = z, no xy returned
        gradRhoUVZ(:,1:1) = scheme%Grad%RhoUVZ_ptr( &
          &   auxField     = auxField,              &
          &   gradData     = gradData,              &
          &   velPos       = velpos,                &
          &   densPos      = denspos,               &
          &   nAuxScalars  = varSys%nAuxScalars,    &
          &   nDims        = 3,                     &
          &   nSolve       = 1,                     &
          &   elemOffset   = iElem-1                )

        ! Calculate correction
        call HRR_Correction_d3q19 (               &
          &    QQ         = QQ,                   &
          &    weight     = layout%weight(:),     &
          &    gradRHOU3  = gradRHOU3(:, 1),      &
          &    gradRHOUVZ = gradRHOUVZ(:, 1),     &
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

        call f_f_eq_regularized_4th_ord_d3q19( layout%weight(:), rho, u_x, u_y, u_z, &
          &                        feq, f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

        do iDir = 1, QQ
          outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = feq(iDir) &
            &                              + cmpl_o*f1(iDir) + 0.5_rk * S_corr(iDir)
        enddo

      enddo nodeloop
!$omp end do nowait

    end associate

  end subroutine bgk_HybridRecursiveRegularizedCorr_d3q19
! **************************************************************************** !

! **************************************************************************** !
  !> Recursive Regularized relaxation routine for the D3Q19
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005

  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  !! works for both d3q19 and d3q27
  subroutine bgk_DualRelaxationTime_RR_d3q19( fieldProp, inState, outState, auxField, &
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
    real(kind=rk) :: f( QQ ), SOM(6), SOM_neq(6), feq(QQ), f_temp, f1(QQ)
    real(kind=rk) :: rho, u_x, u_y, u_z, a12xx, a12xy, a12yy, a12zz, a12xz, a12yz
    real(kind=rk) :: omega, tau, tauN, CoefTauNTau
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

      call f_f_eq_regularized_2nd_ord_d3q19( layout%weight(:), rho, u_x, u_y, u_z, feq, &
        &                                   f1, a12xx, a12yy, a12zz, a12xy, a12xz, a12yz )

      do iDir = 1, QQ
        f_temp = f(iDir) - 1.0_rk/tauN * ( f(iDir) - feq(iDir) )

        outState( ( ielem-1)* nscalars+ idir+(1-1)* qq) = &
          & f_temp + CoefTauNTau * f1(iDir)
      enddo

    enddo nodeloop
!$omp end do nowait

  end subroutine bgk_DualRelaxationTime_RR_d3q19
! **************************************************************************** !


end module mus_d3q19_module
! **************************************************************************** !
