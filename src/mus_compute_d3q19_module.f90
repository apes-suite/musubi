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
! Copyright (c) 2025 Mengyu Wang <m.wang-2@utwente.nl>
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
  use mus_species_module, only: Dxx, Dyy, Dzz, Dxy, Dxz, Dyz

  implicit none

  private

  public :: mus_advRel_kFluid_rBGK_vStd_lD3Q19
  public :: bgk_advRel_d3q19_GNS
  public :: mus_advRel_kFluidIncomp_rBGK_vStd_lD3Q19
  public :: mus_advRel_kFluidIncompGNS_rBGK_vStd_lD3Q19
  public :: mus_advRel_kFluid_rTRT_vStd_lD3Q19
  public :: mus_advRel_kFluidIncomp_rTRT_vStd_lD3Q19
  public :: mus_advRel_kFluid_rBGK_vBlock_lD3Q19
  public :: bgk_Regularized_d3q19
  public :: bgk_RecursiveRegularized_d3q19
  public :: bgk_HybridRecursiveRegularized_d3q19
  public :: bgk_ProjectedRecursiveRegularized_d3q19
  public :: bgk_HybridRecursiveRegularizedCorr_d3q19
  public :: bgk_DualRelaxationTime_RR_d3q19

  ! passive scalar routines (anisotropic diffusion)
  public :: mus_advRel_kPS_rBGK_vEmodel_lD3Q19
  public :: mus_advRel_kPS_rTRT_vEmodel_lD3Q19
  public :: mus_advRel_kPS_rBGK_vEmodelCorr_lD3Q19
  public :: mus_advRel_kPS_rTRT_vEmodelCorr_lD3Q19
  public :: mus_advRel_kPS_rTRT_vLmodel_lD3Q19
  public :: mus_advRel_kPS_rMRT_vEmodelCorr_lD3Q19

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
  !> Advection relaxation routine for the D3Q19 model with BGK based on the !!
  !! equilibrium distribution function for the generalized Navier Stokes equations
  !! (GNS) aka Volume Averaged Navier-Stokes !! equations (VANS).
  !! feq definition from: Z. Guo and T. S. Zhao, “Lattice Boltzmann model for
  !! incompressible flows through porous media,” Phys. Rev. E, vol. 66, no. 3, p.
  !! 036304, Sep. 2002, doi: 10.1103/PhysRevE.66.036304.
  !!
  !! \[ f_\alpha(x_i+e_{\alpha,i},t+1) =
  !! f_\alpha(x_i,t) - \omega(f_\alpha(x_i,t)-f^{eq}_{\alpha}(x_i,t)) \]
  !!
  !! The number of floating point operation in this routine is 160 roughly.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_advRel_d3q19_GNS( fieldProp, inState, outState, auxField, &
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
integer :: iElem, nScalars
real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
  &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
  &              f000
real(kind=rk) :: rho     ! local density
real(kind=rk) :: u_x     ! local x-velocity
real(kind=rk) :: u_y     ! local y-velocity
real(kind=rk) :: u_z     ! local z-velocity
real(kind=rk) :: eps_f   ! local fluid volume fraction
real(kind=rk) :: eps_f_inv   ! 1 divided by local fluid volume fraction
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
integer :: dens_pos, vel_pos(3), vol_frac_pos, elemOff
! -------------------------------------------------------------------- !
dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
vol_frac_pos = varSys%method%val(derVarPos(1)%vol_frac)%auxField_varPos(1)

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

  ! local fluid volume fraction
  eps_f = auxField(elemOff + vol_frac_pos)
  eps_f_inv = 1/eps_f

  ! square velocity and derived constants
  usq  = (u_x * u_x) + (u_y * u_y) + (u_z * u_z)

  ! Modified for Generalized Navier-Stokes!
  usqn = div1_36 * (1._rk - 1.5_rk * usq * eps_f_inv) * rho
  ! usqn = div1_36 * (1._rk - 1.5_rk * usq) * rho

  ! read the relaxation parameter omega for the current level
  omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
  ! pre-calculate partial collision constants
  cmpl_o = 1._rk - omega

  ! f = (1-w) * f + w * fEq

  ! --- Compute f000: PDF at rest position --- !
  outState(( ielem-1)* nscalars+q000+( 1-1)* qq) &
    & = f000*cmpl_o+omega*rho*(div1_3-0.5_rk*usq*eps_f_inv)

  ! --- Compute f110 and fNN0 --- !
  coeff_1 = div1_8 * omega * rho

  usqn_o1 = omega * usqn

  ui1 = u_x + u_y
  fac_1 = coeff_1 * ui1
  sum1_1 = fac_1 * div3_4h
  sum1_2 = fac_1 * ui1 * eps_f_inv + usqn_o1    ! GNS
  ! sum1_2 = fac_1 * ui1 + usqn_o1

  outState(( ielem-1)* nscalars+q110+( 1-1)* qq) &
    & = f110 * cmpl_o + sum1_1 + sum1_2
  outState(( ielem-1)* nscalars+qnn0+( 1-1)* qq) &
    & = fNN0 * cmpl_o - sum1_1 + sum1_2

  ! --- Compute fN10 and f1N0 --- !
  ui3 = -u_x + u_y
  fac_3 = coeff_1 * ui3
  sum3_1 = fac_3 * div3_4h
  sum3_2 = fac_3 * ui3 * eps_f_inv + usqn_o1   ! GNS
  ! sum3_2 = fac_3 * ui3 + usqn_o1

  outState(( ielem-1)* nscalars+qn10+( 1-1)* qq) &
    & = fN10 * cmpl_o + sum3_1 + sum3_2
  outState(( ielem-1)* nscalars+q1n0+( 1-1)* qq) &
    & = f1N0 * cmpl_o - sum3_1 + sum3_2

  ! --- Compute f101 and fN0N --- !
  ui10 = u_x + u_z
  fac_10 = coeff_1 * ui10
  sum10_1 = fac_10 * div3_4h
  sum10_2 = fac_10 * ui10 * eps_f_inv + usqn_o1   ! GNS
  ! sum10_2 = fac_10 * ui10 + usqn_o1

  outState(( ielem-1)* nscalars+q101+( 1-1)* qq) &
    & = f101 * cmpl_o + sum10_1 + sum10_2
  outState(( ielem-1)* nscalars+qn0n+( 1-1)* qq) &
    & = fN0N * cmpl_o - sum10_1 + sum10_2

  ! --- Compute fN01 and f10N --- !
  ui12 = -u_x + u_z
  fac_12 = coeff_1 * ui12
  sum12_1 = fac_12 * div3_4h
  sum12_2 = fac_12 * ui12 * eps_f_inv + usqn_o1    ! GNS
  ! sum12_2 = fac_12 * ui12 + usqn_o1

  outState(( ielem-1)* nscalars+qn01+( 1-1)* qq) &
    & = fN01 * cmpl_o + sum12_1 + sum12_2
  outState(( ielem-1)* nscalars+q10n+( 1-1)* qq) &
    & = f10N * cmpl_o - sum12_1 + sum12_2

  ! --- Compute f011 and f0NN --- !
  ui11 = u_y + u_z
  fac_11 = coeff_1 * ui11
  sum11_1 = fac_11 * div3_4h
  sum11_2 = fac_11 * ui11 * eps_f_inv + usqn_o1 ! GNS
  ! sum11_2 = fac_11 * ui11 + usqn_o1

  outState(( ielem-1)* nscalars+q011+( 1-1)* qq) &
    & = f011 * cmpl_o + sum11_1 + sum11_2
  outState(( ielem-1)* nscalars+q0nn+( 1-1)* qq) &
    & = f0NN * cmpl_o - sum11_1 + sum11_2

  ! --- Compute f0N1 and f01N --- !
  ui13 = -u_y + u_z
  fac_13 = coeff_1 * ui13
  sum13_1 = fac_13 * div3_4h
  sum13_2 = fac_13 * ui13 * eps_f_inv + usqn_o1 ! GNS
  ! sum13_2 = fac_13 * ui13 + usqn_o1

  outState(( ielem-1)* nscalars+q0n1+( 1-1)* qq) &
    & = f0N1 * cmpl_o + sum13_1 + sum13_2
  outState(( ielem-1)* nscalars+q01n+( 1-1)* qq) &
    & = f01N * cmpl_o - sum13_1 + sum13_2

  ! --- Compute f010 and f0N0 --- !
  omega_2 = 2._rk * omega
  coeff_2 = div1_8 * omega_2 * rho
  usqn_o2 = omega_2 * usqn

  fac_2 = coeff_2 * u_y
  sum2_1 = fac_2 * div3_4h
  sum2_2 = fac_2 * u_y * eps_f_inv + usqn_o2

  outState(( ielem-1)* nscalars+q010+( 1-1)* qq) &
    & = f010 * cmpl_o + sum2_1 + sum2_2
  outState(( ielem-1)* nscalars+q0n0+( 1-1)* qq) &
    & = f0N0 * cmpl_o - sum2_1 + sum2_2

  ! --- Compute fN00 and f100 --- !
  fac_4 = coeff_2 * u_x
  sum4_1 = fac_4 * div3_4h
  sum4_2 = fac_4 * u_x * eps_f_inv + usqn_o2

  outState(( ielem-1)* nscalars+qn00+( 1-1)* qq) &
    & = fN00 * cmpl_o - sum4_1 + sum4_2
  outState(( ielem-1)* nscalars+q100+( 1-1)* qq) &
    & = f100 * cmpl_o + sum4_1 + sum4_2

  ! --- Compute f00N and f001 --- !
  fac_9 = coeff_2 * u_z
  sum9_1 = fac_9 * div3_4h
  sum9_2 = fac_9 * u_z * eps_f_inv + usqn_o2

  outState(( ielem-1)* nscalars+q001+( 1-1)* qq) &
    & = f001 * cmpl_o + sum9_1 + sum9_2
  outState(( ielem-1)* nscalars+q00n+( 1-1)* qq) &
    & = f00N * cmpl_o - sum9_1 + sum9_2

enddo nodeloop
!$omp end do nowait

end subroutine bgk_advRel_d3q19_GNS


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
  !> Advection relaxation routine for the D3Q19 model with BGK based on the !!
  !! equilibrium distribution function for the generalized Navier Stokes equations
  !! (GNS) aka Volume Averaged Navier-Stokes !! equations (VANS).
  !! feq definition from: Z. Guo and T. S. Zhao, “Lattice Boltzmann model for
  !! incompressible flows through porous media,” Phys. Rev. E, vol. 66, no. 3, p.
  !! 036304, Sep. 2002, doi: 10.1103/PhysRevE.66.036304.
  !! Incompressible version
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluidIncompGNS_rBGK_vStd_lD3Q19( &
    &          fieldProp, inState, outState, auxField,    &
    &          neigh, nElems, nSolve, level, layout,      &
    &          params, varSys, derVarPos                  )
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
    real(kind=rk) :: eps_f, eps_f_inv
    real(kind=rk) :: usqn_o1, usqn_o2
    real(kind=rk) :: cmpl_o, omega
    real(kind=rk) :: coeff_1, coeff_2
    real(kind=rk) :: ui1, ui3, ui10, ui11, ui12, ui13
    real(kind=rk) :: fac_1, fac_2, fac_3, fac_4, fac_9, &
      &              fac_10, fac_11, fac_12, fac_13
    real(kind=rk) :: sum1_1, sum1_2, sum2_1, sum2_2, sum3_1, sum3_2, sum4_1, &
      &              sum4_2, sum9_1, sum9_2, sum10_1, sum10_2, sum11_1,      &
      &              sum11_2, sum12_1, sum12_2, sum13_1, sum13_2
    integer :: dens_pos, vel_pos(3), vol_frac_pos, elemOff
    ! -------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    vol_frac_pos = varSys%method%val(derVarPos(1)%vol_frac)%auxField_varPos(1)

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

      ! local fluid volume fraction
      eps_f = auxField(elemOff + vol_frac_pos)
      eps_f_inv = 1/eps_f

      ! square velocity and derived constants
      usq = u_x * u_x + u_y * u_y + u_z * u_z

      ! read the relaxation parameter omega for the current level
      omega = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      ! pre-calculate partial collision constants
      cmpl_o = 1._rk - omega

      ! usqn = div1_36 * (rho - 1.5_rk * usq * rho0 )
      usqn_o1 = omega * div1_36 * ( rho - 1.5d0 * usq * eps_f_inv )

      outState(( ielem-1)* qq+q000+( 1-1)* qq) &
        & = f000 * cmpl_o + 12d0 * usqn_o1

      coeff_1 = div1_8 * omega

      ui1 = u_x + u_y
      fac_1 = coeff_1 * ui1
      sum1_1 = fac_1 * div3_4h
      sum1_2 = fac_1 * ui1 * eps_f_inv + usqn_o1

      outState(( ielem-1)* qq+q110+( 1-1)* qq) &
        & = f110 * cmpl_o + sum1_1 + sum1_2
      outState(( ielem-1)* qq+qnn0+( 1-1)* qq) &
        & = fNN0 * cmpl_o - sum1_1 + sum1_2

      ui3 = -u_x + u_y
      fac_3 = coeff_1 * ui3
      sum3_1 = fac_3 * div3_4h
      sum3_2 = fac_3 * ui3 * eps_f_inv + usqn_o1

      outState(( ielem-1)* qq+qn10+( 1-1)* qq) &
        & = fN10 * cmpl_o + sum3_1 + sum3_2
      outState(( ielem-1)* qq+q1n0+( 1-1)* qq) &
        & = f1N0 * cmpl_o - sum3_1 + sum3_2

      ui10 = u_x + u_z
      fac_10 = coeff_1 * ui10
      sum10_1 = fac_10 * div3_4h
      sum10_2 = fac_10 * ui10 * eps_f_inv + usqn_o1

      outState(( ielem-1)* qq+q101+( 1-1)* qq) &
        & = f101 * cmpl_o + sum10_1 + sum10_2
      outState(( ielem-1)* qq+qn0n+( 1-1)* qq) &
        & = fN0N * cmpl_o - sum10_1 + sum10_2

      ui12 = -u_x + u_z
      fac_12 = coeff_1 * ui12
      sum12_1 = fac_12 * div3_4h
      sum12_2 = fac_12 * ui12 * eps_f_inv + usqn_o1

      outState(( ielem-1)* qq+qn01+( 1-1)* qq) &
        & = fN01 * cmpl_o + sum12_1 + sum12_2
      outState(( ielem-1)* qq+q10n+( 1-1)* qq) &
        & = f10N * cmpl_o - sum12_1 + sum12_2

      ui11 = u_y + u_z
      fac_11 = coeff_1 * ui11
      sum11_1 = fac_11 * div3_4h
      sum11_2 = fac_11 * ui11 * eps_f_inv + usqn_o1

      outState(( ielem-1)* qq+q011+( 1-1)* qq) &
        & = f011 * cmpl_o + sum11_1 + sum11_2
      outState(( ielem-1)* qq+q0nn+( 1-1)* qq) &
        & = f0NN * cmpl_o - sum11_1 + sum11_2

      ui13 = -u_y + u_z
      fac_13 = coeff_1 * ui13
      sum13_1 = fac_13 * div3_4h
      sum13_2 = fac_13 * ui13 * eps_f_inv + usqn_o1

      outState(( ielem-1)* qq+q0n1+( 1-1)* qq) &
        & = f0N1 * cmpl_o + sum13_1 + sum13_2
      outState(( ielem-1)* qq+q01n+( 1-1)* qq) &
        & = f01N * cmpl_o - sum13_1 + sum13_2

      coeff_2 = div1_8 * omega * 2.0_rk
      usqn_o2 = 2d0 * usqn_o1

      fac_2 = coeff_2 * u_y
      sum2_1 = fac_2 * div3_4h
      sum2_2 = fac_2 * u_y * eps_f_inv + usqn_o2

      outState(( ielem-1)* qq+q010+( 1-1)* qq) &
        & = f010 * cmpl_o + sum2_1 + sum2_2
      outState(( ielem-1)* qq+q0n0+( 1-1)* qq) &
        & = f0N0 * cmpl_o - sum2_1 + sum2_2

      fac_4 = coeff_2 * u_x
      sum4_1 = fac_4 * div3_4h
      sum4_2 = fac_4 * u_x * eps_f_inv + usqn_o2

      outState(( ielem-1)* qq+qn00+( 1-1)* qq) &
        & = fN00 * cmpl_o - sum4_1 + sum4_2
      outState(( ielem-1)* qq+q100+( 1-1)* qq) &
        & = f100 * cmpl_o + sum4_1 + sum4_2

      fac_9 = coeff_2 * u_z
      sum9_1 = fac_9 * div3_4h
      sum9_2 = fac_9 * u_z * eps_f_inv + usqn_o2

      outState(( ielem-1)* qq+q001+( 1-1)* qq) &
        & = f001 * cmpl_o + sum9_1 + sum9_2
      outState(( ielem-1)* qq+q00n+( 1-1)* qq) &
        & = f00N * cmpl_o - sum9_1 + sum9_2

    enddo nodeloop
!$omp end do nowait

  end subroutine mus_advRel_kFluidIncompGNS_rBGK_vStd_lD3Q19
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

! **************************************************************************** !
  !> Advection relaxation routine based on the E-Model for the D3Q19
  !! Lattice Boltzmann model for the generic advection and anisotropic-dispersion
  !! equation
  !!
  !!   Irina Ginzburg (2005), "Equilibrium-type and link-type lattice Boltzmann
  !!   models for generic advection and anisotropic-dispersion equation",
  !!   Advances in Water Resources, Volume 28, Issue 11
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kPS_rBGK_vEmodel_lD3Q19( fieldProp, inState, outState,   &
    &                            auxField, neigh, nElems, nSolve, level, layout, &
    &                            params, varSys, derVarPos                       )
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
    integer :: iElem, iDir
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
    real(kind=rk) :: rho, feq
    real(kind=rk) :: d_omega, nu_q
    real(kind=rk) :: transVel( nSolve*3 ) ! velocity from the transport field
    real(kind=rk) :: uc, cqx2, cqy2, cqz2, cq2
    real(kind=rk) :: a_e, a_xx, a_ww, a_xy, a_xz, a_yz
    real(kind=rk) :: P_e, P_xx, P_ww, P_xy, P_xz, P_yz
    integer :: vel_varPos ! position of transport velocity variable in varSys
    real(kind=rk) :: inv_vel, u_fluid(3)
    ! --------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(derVarPos(1)%pdf)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    ! passive scalar has only one transport Variable
    vel_varPos = scheme%transVar%method(1)%data_varPos
    ! Get velocity field
    call varSys%method%val(vel_varPos)%get_valOfIndex( &
      & varSys  = varSys,                              &
      & time    = params%general%simControl%now,       &
      & iLevel  = level,                               &
      & idx     = scheme%transVar%method(1)            &
      &           %pntIndex%indexLvl(level)            &
      &           %val(1:nSolve),                      &
      & nVals   = nSolve,                              &
      & res     = transVel                             )

    ! convert physical velocity into LB velocity
    inv_vel = 1.0_rk / params%physics%fac( level )%vel
    transVel = transVel * inv_vel

    ! for isotropic diffusion factor, it turns out to be 1st order bgk
    d_omega = fieldProp(1)%species%omega
    ! reciprocal of the free parameter nu, i.e. nu_q = 1/nu
    nu_q = cs2inv / (1.0_rk / d_omega - 0.5_rk)
    a_e = ( fieldProp(1)%species%diff_tensor(Dxx)  &
      &    + fieldProp(1)%species%diff_tensor(Dyy) &
      &    + fieldProp(1)%species%diff_tensor(Dzz) ) * nu_q / 3.0_rk - 1.0_rk
    a_xx = 2._rk / 3._rk * nu_q * ( fieldProp(1)%species%diff_tensor(Dxx)                &
      &                            - 0.5_rk * ( fieldProp(1)%species%diff_tensor(Dyy)    &
      &                                        + fieldProp(1)%species%diff_tensor(Dzz) ) &
      &                           )
    a_ww = 0.5_rk * nu_q * ( fieldProp(1)%species%diff_tensor(Dyy)  &
      &                     - fieldProp(1)%species%diff_tensor(Dzz) )

    ! attention: D_ab has a multiplier of 2 in ADE
    a_xy = fieldProp(1)%species%diff_tensor(Dxy) * nu_q
    a_xz = fieldProp(1)%species%diff_tensor(Dxz) * nu_q
    a_yz = fieldProp(1)%species%diff_tensor(Dyz) * nu_q

    elemloop: do iElem = 1, nSolve
      ! x-, y- and z-velocity from transport field
      u_fluid = transVel( (iElem-1)*3+1 : iElem*3 )

      do iDir = 1, layout%fStencil%QQ
        pdfTmp( iDir ) = instate( neigh( (idir - 1) * nelems + ielem ) )
      end do
      rho = sum( pdfTmp )

      do iDir = 1, layout%fStencil%QQ
        ! compute c_i * u
        uc = real(layout%fStencil%cxDir(1, iDir), kind=rk ) * u_fluid(1)  &
          & + real(layout%fStencil%cxDir(2, iDir), kind=rk ) * u_fluid(2) &
          & + real(layout%fStencil%cxDir(3, iDir), kind=rk ) * u_fluid(3)

        cqx2 = real( layout%fStencil%cxDir(1, iDir), kind=rk ) &
          &   * real( layout%fStencil%cxDir(1, iDir), kind=rk )
        cqy2 = real( layout%fStencil%cxDir(2, iDir), kind=rk ) &
          &   * real( layout%fStencil%cxDir(2, iDir), kind=rk )
        cqz2 = real( layout%fStencil%cxDir(3, iDir), kind=rk ) &
          &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )
        cq2 = cqx2 + cqy2 + cqz2

        P_e = (19 * cq2 - 30) / 42
        P_xx = (3 * cqx2 - cq2) / 12
        P_ww = (cqy2 - cqz2) / 6
        P_xy = (real( layout%fStencil%cxDir(1, iDir), kind=rk )    &
          &   * real( layout%fStencil%cxDir(2, iDir), kind=rk )) / 4
        P_xz = (real( layout%fStencil%cxDir(1, iDir), kind=rk )    &
          &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )) / 4
        P_yz = (real( layout%fStencil%cxDir(2, iDir), kind=rk )    &
          &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )) / 4

        feq = rho * ( layout%weight( iDir )               &   ! e_0
          &         + layout%weight( iDir ) * cs2inv * uc &   ! C_\alpha
          &         + cs2 * (a_e * P_e + a_xy * P_xy      &
          &                 + a_xz * P_xz + a_yz * P_yz   &
          &                 + a_xx * P_xx + a_ww * P_ww)  )

        outstate( (ielem - 1) * varsys%nscalars + idir ) = &
          &                pdfTmp(iDir) + d_omega * ( feq - pdfTmp(iDir) )

      end do

    end do elemloop

    end subroutine mus_advRel_kPS_rBGK_vEmodel_lD3Q19
  ! ****************************************************************************** !


! **************************************************************************** !
  !> Advection relaxation routine based on the E-Model for the D3Q19
  !! Lattice Boltzmann model for the generic advection and anisotropic-dispersion
  !! equation
  !!
  !!   Irina Ginzburg (2005), "Equilibrium-type and link-type lattice Boltzmann
  !!   models for generic advection and anisotropic-dispersion equation",
  !!   Advances in Water Resources, Volume 28, Issue 11
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kPS_rTRT_vEmodel_lD3Q19( fieldProp, inState, outState,   &
    &                            auxField, neigh, nElems, nSolve, level, layout, &
    &                            params, varSys, derVarPos               )
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
      integer :: iElem, iDir, invDir
      type(mus_varSys_data_type), pointer :: fPtr
      type(mus_scheme_type), pointer :: scheme
      real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
      real(kind=rk) :: rho, feqPlus, feqMinus, fPlus, fMinus
      real(kind=rk) :: d_omega, aux_omega, nu_q, d_lambda
      real(kind=rk) :: transVel( nSolve*3 ) ! velocity from the transport field
      real(kind=rk) :: uc, cqx2, cqy2, cqz2, cq2
      real(kind=rk) :: a_e, a_xx, a_ww, a_xy, a_xz, a_yz
      real(kind=rk) :: P_e, P_xx, P_ww, P_xy, P_xz, P_yz
      integer :: vel_varPos ! position of transport velocity variable in varSys
      real(kind=rk) :: inv_vel, u_fluid(3)
      ! --------------------------------------------------------------------------
      ! access scheme via 1st variable method data which is a state variable
      call C_F_POINTER( varSys%method%val(derVarPos(1)%pdf)%method_Data, fPtr )
      scheme => fPtr%solverData%scheme

      ! passive scalar has only one transport Variable
      vel_varPos = scheme%transVar%method(1)%data_varPos
      ! Get velocity field
      call varSys%method%val(vel_varPos)%get_valOfIndex( &
        & varSys  = varSys,                              &
        & time    = params%general%simControl%now,       &
        & iLevel  = level,                               &
        & idx     = scheme%transVar%method(1)            &
        &           %pntIndex%indexLvl(level)            &
        &           %val(1:nSolve),                      &
        & nVals   = nSolve,                              &
        & res     = transVel                             )

      ! convert physical velocity into LB velocity
      inv_vel = 1.0_rk / params%physics%fac( level )%vel
      transVel = transVel * inv_vel

      d_omega = fieldProp(1)%species%omega
      d_lambda = fieldProp(1)%species%lambda
      aux_omega = 1.0_rk / (d_lambda / (1.0_rk / d_omega - 0.5_rk) + 0.5_rk)

      ! for isotropic diffusion factor, it turns out to be 1st order trt
      ! reciprocal of the free parameter nu, i.e. nu_q = 1/nu
      nu_q = cs2inv / (1.0_rk / d_omega - 0.5_rk)
      a_e = ( fieldProp(1)%species%diff_tensor(Dxx)   &
        &    + fieldProp(1)%species%diff_tensor(Dyy)  &
        &    + fieldProp(1)%species%diff_tensor(Dzz) ) * nu_q / 3.0_rk - 1.0_rk
      a_xx = 2._rk / 3._rk * nu_q * ( fieldProp(1)%species%diff_tensor(Dxx)               &
        &                            - 0.5_rk * (fieldProp(1)%species%diff_tensor(Dyy)    &
        &                                        + fieldProp(1)%species%diff_tensor(Dzz)) &
        &                           )
      a_ww = 0.5_rk * nu_q * ( fieldProp(1)%species%diff_tensor(Dyy) &
        &                     - fieldProp(1)%species%diff_tensor(Dzz) )

      ! attention: D_ab has a multiplier of 2 in ADE
      a_xy = fieldProp(1)%species%diff_tensor(Dxy) * nu_q
      a_xz = fieldProp(1)%species%diff_tensor(Dxz) * nu_q
      a_yz = fieldProp(1)%species%diff_tensor(Dyz) * nu_q


      elemloop: do iElem = 1, nSolve
        ! x-, y- and z-velocity from transport field
        u_fluid = transVel( (iElem-1)*3+1 : iElem*3 )

        do iDir = 1, layout%fStencil%QQ
          pdfTmp( iDir ) = instate( neigh( (idir - 1) * nelems + ielem ) )
        end do
        rho = sum( pdfTmp )

        do iDir = 1, layout%fStencil%QQ
          ! compute c_i * u
          uc = real( layout%fStencil%cxDir(1, iDir), kind=rk ) * u_fluid(1)  &
            & + real( layout%fStencil%cxDir(2, iDir), kind=rk ) * u_fluid(2) &
            & + real( layout%fStencil%cxDir(3, iDir), kind=rk ) * u_fluid(3)

          cqx2 = real( layout%fStencil%cxDir(1, iDir), kind=rk ) &
            &   * real( layout%fStencil%cxDir(1, iDir), kind=rk )
          cqy2 = real( layout%fStencil%cxDir(2, iDir), kind=rk ) &
            &   * real( layout%fStencil%cxDir(2, iDir), kind=rk )
          cqz2 = real( layout%fStencil%cxDir(3, iDir), kind=rk ) &
            &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )
          cq2 = cqx2 + cqy2 + cqz2

          P_e = (19 * cq2 - 30) / 42
          P_xx = (3 * cqx2 - cq2) / 12
          P_ww = (cqy2 - cqz2) / 6
          P_xy = (real( layout%fStencil%cxDir(1, iDir), kind=rk )    &
            &   * real( layout%fStencil%cxDir(2, iDir), kind=rk )) / 4
          P_xz = (real( layout%fStencil%cxDir(1, iDir), kind=rk )    &
            &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )) / 4
          P_yz = (real( layout%fStencil%cxDir(2, iDir), kind=rk )    &
            &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )) / 4

          ! compute the equilibrium (fi_eq = weight_i * rho * ( 1+c_i*u / cs^2))
          feqPlus = rho * ( layout%weight( iDir )                & ! e_0
            &              + cs2 * (a_e * P_e + a_xy * P_xy      &
            &                       + a_xz * P_xz + a_yz * P_yz  &
            &                       + a_xx * P_xx + a_ww * P_ww) &
            &              )

          feqMinus = rho * layout%weight( iDir ) * 3._rk * uc

          invDir = layout%fStencil%cxDirInv(iDir)
          fPlus = 0.5_rk * (pdfTmp(iDir) + pdfTmp(invDir))
          fMinus = 0.5_rk * (pdfTmp(iDir) - pdfTmp(invDir))

          outstate( (ielem - 1) * varsys%nscalars + idir ) =              &
            &                pdfTmp(iDir) + d_omega * (feqMinus - fMinus) &
            &                 + aux_omega * (feqPlus - fPlus)
        end do

      end do elemloop

      end subroutine mus_advRel_kPS_rTRT_vEmodel_lD3Q19
    ! ****************************************************************************** !



! **************************************************************************** !
  !> Advection relaxation routine based on the E-Model for the D3Q19
  !! Lattice Boltzmann model for the generic advection and anisotropic-dispersion
  !! equation. The E-Model is corrected to account for the effect of numerical
  !! diffusion with the nonlinear equilibrium correction
  !! \[ f_i^{eq,cor} = f_i^{eq} + \frac{w_i}{2 c_s^4} (c_{i\alpha} c_{i\beta}
  !!      - c_s^2 \delta_{\alpha\beta}) u_\alpha u_\beta \]
  !!
  !!   Irina Ginzburg (2005), "Equilibrium-type and link-type lattice Boltzmann
  !!   models for generic advection and anisotropic-dispersion equation",
  !!   Advances in Water Resources, Volume 28, Issue 11
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kPS_rBGK_vEmodelCorr_lD3Q19( fieldProp, inState, outState, &
    &                            auxField, neigh, nElems, nSolve, level, layout,   &
    &                            params, varSys, derVarPos               )
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
    integer :: iElem, iDir
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
    real(kind=rk) :: rho, feq
    real(kind=rk) :: d_omega, nu_q
    real(kind=rk) :: transVel( nSolve*3 ) ! velocity from the transport field
    real(kind=rk) :: uc, cqx2, cqy2, cqz2, cq2, usq
    real(kind=rk) :: a_e, a_xx, a_ww, a_xy, a_xz, a_yz
    real(kind=rk) :: P_e, P_xx, P_ww, P_xy, P_xz, P_yz
    integer :: vel_varPos ! position of transport velocity variable in varSys
    real(kind=rk) :: inv_vel, u_fluid(3)
    ! --------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(derVarPos(1)%pdf)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    ! passive scalar has only one transport Variable
    vel_varPos = scheme%transVar%method(1)%data_varPos
    ! Get velocity field
    call varSys%method%val(vel_varPos)%get_valOfIndex( &
      & varSys  = varSys,                              &
      & time    = params%general%simControl%now,       &
      & iLevel  = level,                               &
      & idx     = scheme%transVar%method(1)            &
      &           %pntIndex%indexLvl(level)            &
      &           %val(1:nSolve),                      &
      & nVals   = nSolve,                              &
      & res     = transVel                             )

    ! convert physical velocity into LB velocity
    inv_vel = 1.0_rk / params%physics%fac( level )%vel
    transVel = transVel * inv_vel

    ! for isotropic diffusion factor, it turns out to be 1st order bgk
    d_omega = fieldProp(1)%species%omega
    ! reciprocal of the free parameter nu, i.e. nu_q = 1/nu
    nu_q = cs2inv / (1.0_rk / d_omega - 0.5_rk)
    a_e = ( fieldProp(1)%species%diff_tensor(Dxx)  &
      &    + fieldProp(1)%species%diff_tensor(Dyy) &
      &    + fieldProp(1)%species%diff_tensor(Dzz) ) * nu_q / 3.0_rk - 1.0_rk
    a_xx = 2._rk / 3._rk * nu_q * ( fieldProp(1)%species%diff_tensor(Dxx)               &
      &                            - 0.5_rk * (fieldProp(1)%species%diff_tensor(Dyy)    &
      &                                        + fieldProp(1)%species%diff_tensor(Dzz)) &
      &                           )
    a_ww = 0.5_rk * nu_q * ( fieldProp(1)%species%diff_tensor(Dyy)  &
      &                     - fieldProp(1)%species%diff_tensor(Dzz) )

    ! attention: D_ab has a multiplier of 2 in ADE
    a_xy = fieldProp(1)%species%diff_tensor(Dxy) * nu_q
    a_xz = fieldProp(1)%species%diff_tensor(Dxz) * nu_q
    a_yz = fieldProp(1)%species%diff_tensor(Dyz) * nu_q

    elemloop: do iElem = 1, nSolve
      ! x-, y- and z-velocity from transport field
      u_fluid = transVel( (iElem-1)*3+1 : iElem*3 )

      do iDir = 1, layout%fStencil%QQ
        pdfTmp( iDir ) = instate( neigh( (idir - 1) * nelems + ielem ) )
      end do
      rho = sum( pdfTmp )

      do iDir = 1, layout%fStencil%QQ
        ! compute c_i * u
        uc = real( layout%fStencil%cxDir(1, iDir), kind=rk ) * u_fluid(1)  &
          & + real( layout%fStencil%cxDir(2, iDir), kind=rk ) * u_fluid(2) &
          & + real( layout%fStencil%cxDir(3, iDir), kind=rk ) * u_fluid(3)

        usq = u_fluid(1)*u_fluid(1) + u_fluid(2)*u_fluid(2) + u_fluid(3)*u_fluid(3)

        cqx2 = real( layout%fStencil%cxDir(1, iDir), kind=rk ) &
          &   * real( layout%fStencil%cxDir(1, iDir), kind=rk )
        cqy2 = real( layout%fStencil%cxDir(2, iDir), kind=rk ) &
          &   * real( layout%fStencil%cxDir(2, iDir), kind=rk )
        cqz2 = real( layout%fStencil%cxDir(3, iDir), kind=rk ) &
          &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )
        cq2 = cqx2 + cqy2 + cqz2

        P_e = (19 * cq2 - 30) / 42
        P_xx = (3 * cqx2 - cq2) / 12
        P_ww = (cqy2 - cqz2) / 6
        P_xy = (real( layout%fStencil%cxDir(1, iDir), kind=rk )    &
          &   * real( layout%fStencil%cxDir(2, iDir), kind=rk )) / 4
        P_xz = (real( layout%fStencil%cxDir(1, iDir), kind=rk )    &
          &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )) / 4
        P_yz = (real( layout%fStencil%cxDir(2, iDir), kind=rk )    &
          &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )) / 4

        feq = rho * layout%weight( iDir )               &
          &   * ( 1 + cs2inv * uc                       & ! e_0 + C_\alpha
          &     + cs4inv * uc * uc * 0.5_rk             &
          &     - usq * 0.5_rk * cs2inv )               & ! Correction of numerical diffusion
          &   + rho * cs2 * ( a_e * P_e + a_xy * P_xy   &
          &                 + a_xz * P_xz + a_yz * P_yz &
          &                 + a_xx * P_xx + a_ww * P_ww )

        outstate( (ielem - 1) * varsys%nscalars + idir ) = &
          &                pdfTmp(iDir) + d_omega * ( feq - pdfTmp(iDir) )

      end do

    end do elemloop

    end subroutine mus_advRel_kPS_rBGK_vEmodelCorr_lD3Q19
  ! ****************************************************************************** !


! **************************************************************************** !
  !> Advection relaxation routine based on the E-Model for the D3Q19
  !! Lattice Boltzmann model for the generic advection and anisotropic-dispersion
  !! equation. The E-Model is corrected to account for the effect of numerical
  !! diffusion with the nonlinear equilibrium correction
  !! \[ f_i^{eq,cor} = f_i^{eq} + \frac{w_i}{2 c_s^4} (c_{i\alpha} c_{i\beta}
  !!      - c_s^2 \delta_{\alpha\beta}) u_\alpha u_\beta \]
  !!
  !!   Irina Ginzburg (2005), "Equilibrium-type and link-type lattice Boltzmann
  !!   models for generic advection and anisotropic-dispersion equation",
  !!   Advances in Water Resources, Volume 28, Issue 11
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kPS_rTRT_vEmodelCorr_lD3Q19( fieldProp, inState, outState, &
    &                            auxField, neigh, nElems, nSolve, level, layout,   &
    &                            params, varSys, derVarPos               )
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
      integer :: iElem, iDir, invDir
      type(mus_varSys_data_type), pointer :: fPtr
      type(mus_scheme_type), pointer :: scheme
      real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
      real(kind=rk) :: rho, feqPlus, feqMinus, fPlus, fMinus
      real(kind=rk) :: d_omega, aux_omega, nu_q, d_lambda
      real(kind=rk) :: transVel( nSolve*3 ) ! velocity from the transport field
      real(kind=rk) :: uc, cqx2, cqy2, cqz2, cq2, usq
      real(kind=rk) :: a_e, a_xx, a_ww, a_xy, a_xz, a_yz
      real(kind=rk) :: P_e, P_xx, P_ww, P_xy, P_xz, P_yz
      integer :: vel_varPos ! position of transport velocity variable in varSys
      real(kind=rk) :: inv_vel, u_fluid(3)
      ! --------------------------------------------------------------------------
      ! access scheme via 1st variable method data which is a state variable
      call C_F_POINTER( varSys%method%val(derVarPos(1)%pdf)%method_Data, fPtr )
      scheme => fPtr%solverData%scheme

      ! passive scalar has only one transport Variable
      vel_varPos = scheme%transVar%method(1)%data_varPos
      ! Get velocity field
      call varSys%method%val(vel_varPos)%get_valOfIndex( &
        & varSys  = varSys,                              &
        & time    = params%general%simControl%now,       &
        & iLevel  = level,                               &
        & idx     = scheme%transVar%method(1)            &
        &           %pntIndex%indexLvl(level)            &
        &           %val(1:nSolve),                      &
        & nVals   = nSolve,                              &
        & res     = transVel                             )

      ! convert physical velocity into LB velocity
      inv_vel = 1.0_rk / params%physics%fac( level )%vel
      transVel = transVel * inv_vel

      d_omega = fieldProp(1)%species%omega
      d_lambda = fieldProp(1)%species%lambda
      aux_omega = 1.0_rk / (d_lambda / (1.0_rk / d_omega - 0.5_rk) + 0.5_rk)

      ! for isotropic diffusion factor, it turns out to be 1st order trt
      ! reciprocal of the free parameter nu, i.e. nu_q = 1/nu
      nu_q = cs2inv / (1.0_rk / d_omega - 0.5_rk)
      a_e = ( fieldProp(1)%species%diff_tensor(Dxx)  &
        &    + fieldProp(1)%species%diff_tensor(Dyy) &
        &    + fieldProp(1)%species%diff_tensor(Dzz) ) * nu_q / 3.0_rk - 1.0_rk
      a_xx = 2._rk / 3._rk * nu_q * ( fieldProp(1)%species%diff_tensor(Dxx)               &
        &                            - 0.5_rk * (fieldProp(1)%species%diff_tensor(Dyy)    &
        &                                        + fieldProp(1)%species%diff_tensor(Dzz)) &
        &                           )
      a_ww = 0.5_rk * nu_q * (fieldProp(1)%species%diff_tensor(Dyy)  &
        &                     - fieldProp(1)%species%diff_tensor(Dzz))

      ! attention: D_ab has a multiplier of 2 in ADE
      a_xy = fieldProp(1)%species%diff_tensor(Dxy) * nu_q
      a_xz = fieldProp(1)%species%diff_tensor(Dxz) * nu_q
      a_yz = fieldProp(1)%species%diff_tensor(Dyz) * nu_q


      elemloop: do iElem = 1, nSolve
        ! x-, y- and z-velocity from transport field
        u_fluid = transVel( (iElem-1)*3+1 : iElem*3 )

        do iDir = 1, layout%fStencil%QQ
          pdfTmp( iDir ) = instate( neigh( (idir - 1) * nelems + ielem ) )
        end do
        rho = sum( pdfTmp )

        do iDir = 1, layout%fStencil%QQ
          ! compute c_i * u
          uc = real( layout%fStencil%cxDir(1, iDir), kind=rk ) * u_fluid(1)  &
            & + real( layout%fStencil%cxDir(2, iDir), kind=rk ) * u_fluid(2) &
            & + real( layout%fStencil%cxDir(3, iDir), kind=rk ) * u_fluid(3)

          usq = u_fluid(1)*u_fluid(1) + u_fluid(2)*u_fluid(2) + u_fluid(3)*u_fluid(3)

          cqx2 = real( layout%fStencil%cxDir(1, iDir), kind=rk ) &
            &   * real( layout%fStencil%cxDir(1, iDir), kind=rk )
          cqy2 = real( layout%fStencil%cxDir(2, iDir), kind=rk ) &
            &   * real( layout%fStencil%cxDir(2, iDir), kind=rk )
          cqz2 = real( layout%fStencil%cxDir(3, iDir), kind=rk ) &
            &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )
          cq2 = cqx2 + cqy2 + cqz2

          P_e = (19 * cq2 - 30) / 42
          P_xx = (3 * cqx2 - cq2) / 12
          P_ww = (cqy2 - cqz2) / 6
          P_xy = (real( layout%fStencil%cxDir(1, iDir), kind=rk )    &
            &   * real( layout%fStencil%cxDir(2, iDir), kind=rk )) / 4
          P_xz = (real( layout%fStencil%cxDir(1, iDir), kind=rk )    &
            &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )) / 4
          P_yz = (real( layout%fStencil%cxDir(2, iDir), kind=rk )    &
            &   * real( layout%fStencil%cxDir(3, iDir), kind=rk )) / 4

          ! compute the equilibrium (fi_eq = weight_i * rho * ( 1+c_i*u / cs^2))
          feqPlus = rho * layout%weight( iDir ) * ( 1                        & ! e_0
            &                                    + cs4inv * uc * uc * 0.5_rk &
            &                                    - usq * 0.5_rk * cs2inv )   & ! Correction
            &       + rho * cs2 * ( a_e * P_e + a_xy * P_xy                  &
            &                      + a_xz * P_xz + a_yz * P_yz               &
            &                      + a_xx * P_xx + a_ww * P_ww               )

          feqMinus = rho * layout%weight( iDir ) * 3._rk * uc

          invDir = layout%fStencil%cxDirInv(iDir)
          fPlus = 0.5_rk * (pdfTmp(iDir) + pdfTmp(invDir))
          fMinus = 0.5_rk * (pdfTmp(iDir) - pdfTmp(invDir))

          outstate( (ielem - 1) * varsys%nscalars + idir ) =              &
            &                pdfTmp(iDir) + d_omega * (feqMinus - fMinus) &
            &                 + aux_omega * (feqPlus - fPlus)

        end do

      end do elemloop

      end subroutine mus_advRel_kPS_rTRT_vEmodelCorr_lD3Q19
    ! ****************************************************************************** !


! **************************************************************************** !
  !> Advection relaxation routine based on the L-Model for the D3Q19
  !! Lattice Boltzmann model for the generic advection and anisotropic-dispersion
  !! equation.
  !!
  !!   Irina Ginzburg (2005), "Equilibrium-type and link-type lattice Boltzmann
  !!   models for generic advection and anisotropic-dispersion equation",
  !!   Advances in Water Resources, Volume 28, Issue 11
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kPS_rTRT_vLmodel_lD3Q19( fieldProp, inState, outState,   &
    &                            auxField, neigh, nElems, nSolve, level, layout, &
    &                            params, varSys, derVarPos               )
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
    integer :: iElem, iDir
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
    real(kind=rk) :: rho, sxy, sxz, syz
    real(kind=rk) :: omega(9), aux_omega
    real(kind=rk) :: transVel( nSolve*3 ) ! velocity from the transport field
    integer :: vel_varPos ! position of transport velocity variable in varSys
    real(kind=rk) :: inv_vel, u_fluid(3)
    ! symmetric and antisymmetric pairs in L-Model
    integer, parameter :: p = 1
    integer, parameter :: m = 2
    real(kind=rk) :: pm(2, 9)
    ! --------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(derVarPos(1)%pdf)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    ! passive scalar has only one transport Variable
    vel_varPos = scheme%transVar%method(1)%data_varPos
    ! Get velocity field
    call varSys%method%val(vel_varPos)%get_valOfIndex( &
      & varSys  = varSys,                              &
      & time    = params%general%simControl%now,       &
      & iLevel  = level,                               &
      & idx     = scheme%transVar%method(1)            &
      &           %pntIndex%indexLvl(level)            &
      &           %val(1:nSolve),                      &
      & nVals   = nSolve,                              &
      & res     = transVel                             )

    ! convert physical velocity into LB velocity
    inv_vel = 1.0_rk / params%physics%fac( level )%vel
    transVel = transVel * inv_vel

    ! define the 3 free parameters
    sxy = ( abs(fieldProp(1)%species%diff_tensor(Dxy))      &
      &      + min( fieldProp(1)%species%diff_tensor(Dxx),  &
      &             fieldProp(1)%species%diff_tensor(Dyy) ) &
      &    ) / 2 ! free parameter for xy direction

    sxz = ( abs(fieldProp(1)%species%diff_tensor(Dxz))      &
      &      + min( fieldProp(1)%species%diff_tensor(Dxx),  &
      &             fieldProp(1)%species%diff_tensor(Dzz) ) &
      &    ) / 2 ! free parameter for xz direction

    syz = ( abs(fieldProp(1)%species%diff_tensor(Dyz))      &
      &      + min( fieldProp(1)%species%diff_tensor(Dyy),  &
      &             fieldProp(1)%species%diff_tensor(Dzz) ) &
      &    ) / 2 ! free parameter for xz direction

    ! \Lambda_xx, yy, zz with t_1 = 1/6
    omega(1) = (fieldProp(1)%species%diff_tensor(Dxx) - sxy - sxz) * 9.0_rk
    omega(2) = (fieldProp(1)%species%diff_tensor(Dyy) - sxy - syz) * 9.0_rk
    omega(3) = (fieldProp(1)%species%diff_tensor(Dzz) - sxz - syz) * 9.0_rk
    ! \Lambda_\pn xy, \pn xz, \pn yz with t_2 = 1/12
    omega(4) = 9.0_rk * (sxy + fieldProp(1)%species%diff_tensor(Dxy))
    omega(5) = 9.0_rk * (sxy - fieldProp(1)%species%diff_tensor(Dxy))
    omega(6) = 9.0_rk * (sxz + fieldProp(1)%species%diff_tensor(Dxz))
    omega(7) = 9.0_rk * (sxz - fieldProp(1)%species%diff_tensor(Dxz))
    omega(8) = 9.0_rk * (syz + fieldProp(1)%species%diff_tensor(Dyz))
    omega(9) = 9.0_rk * (syz - fieldProp(1)%species%diff_tensor(Dyz))
    ! get omega from \Lambda
    omega = 1.0_rk / (omega + 0.5_rk)
    aux_omega = 1.0_rk ! free parameter

    elemloop: do iElem = 1, nSolve
      ! x-, y- and z-velocity from transport field
      u_fluid = transVel( (iElem-1)*3+1 : iElem*3 )

      do iDir = 1, layout%fStencil%QQ
        pdfTmp( iDir ) = instate( neigh( (idir - 1) * nelems + ielem ) )
      end do
      rho = sum( pdfTmp )

      !> compute the link-wise omega
      !! the directions of omegas are
      !!  1     !< west             x-
      !!  2     !< south            y-
      !!  3     !< bottom           z-
      !!  4     !< east             x+
      !!  5     !< north            y+
      !!  6     !< top              z+
      !!  7     !<                  z-,y-
      !!  8     !<                  z+,y-
      !!  9     !<                  z-,y+
      !!  10    !<                  z+,y+
      !!  11    !<                  x-,z-
      !!  12    !<                  x+,z-
      !!  13    !<                  x-,z+
      !!  14    !<                  x+,z+
      !!  15    !<                  y-,x-
      !!  16    !<                  y+,x-
      !!  17    !<                  y-,x+
      !!  18    !<                  y+,x+
      !!  19    !< rest density is last
      !!
      pm(p, 1) = aux_omega * (rho * layout%weight(4) - 0.5_rk  &
        &       * (pdfTmp(4) + pdfTmp(1)))
      pm(m, 1) = omega(1) * (rho * layout%weight(4) * 3.0_rk   &
        &       * u_fluid(1) - 0.5_rk * (pdfTmp(4) - pdfTmp(1)))

      pm(p, 2) = aux_omega * (rho * layout%weight(5) - 0.5_rk  &
        &       * (pdfTmp(5) + pdfTmp(2)))
      pm(m, 2) = omega(2) * (rho * layout%weight(5) * 3.0_rk   &
        &       * u_fluid(2) - 0.5_rk * (pdfTmp(5) - pdfTmp(2)))

      pm(p, 3) = aux_omega * (rho * layout%weight(6) - 0.5_rk  &
        &       * (pdfTmp(6) + pdfTmp(3)))
      pm(m, 3) = omega(3) * (rho * layout%weight(6) * 3.0_rk   &
        &       * u_fluid(3) - 0.5_rk * (pdfTmp(6) - pdfTmp(3)))

      pm(p, 4) = aux_omega * (rho * layout%weight(18) - 0.5_rk &
        &       * (pdfTmp(18) + pdfTmp(15)))
      pm(m, 4) = omega(4) * (rho * layout%weight(18) * 3.0_rk  &
        &       * (u_fluid(1) + u_fluid(2)) - 0.5_rk * (pdfTmp(18) - pdfTmp(15)))

      pm(p, 5) = aux_omega * (rho * layout%weight(17) - 0.5_rk &
        &       * (pdfTmp(17) + pdfTmp(16)))
      pm(m, 5) = omega(5) * (rho * layout%weight(17) * 3.0_rk  &
        &       * (u_fluid(1) - u_fluid(2)) - 0.5_rk * (pdfTmp(17) - pdfTmp(16)))

      pm(p, 6) = aux_omega * (rho * layout%weight(14) - 0.5_rk &
        &       * (pdfTmp(14) + pdfTmp(11)))
      pm(m, 6) = omega(6) * (rho * layout%weight(14) * 3.0_rk  &
        &       * (u_fluid(1) + u_fluid(3)) - 0.5_rk * (pdfTmp(14) - pdfTmp(11)))

      pm(p, 7) = aux_omega * (rho * layout%weight(12) - 0.5_rk &
        &       * (pdfTmp(12) + pdfTmp(13)))
      pm(m, 7) = omega(7) * (rho * layout%weight(12) * 3.0_rk  &
        &       * (u_fluid(1) - u_fluid(3)) - 0.5_rk * (pdfTmp(12) - pdfTmp(13)))

      pm(p, 8) = aux_omega * (rho * layout%weight(10) - 0.5_rk &
        &       * (pdfTmp(10) + pdfTmp(7)))
      pm(m, 8) = omega(8) * (rho * layout%weight(10) * 3.0_rk  &
        &       * (u_fluid(2) + u_fluid(3)) - 0.5_rk * (pdfTmp(10) - pdfTmp(7)))

      pm(p, 9) = aux_omega * (rho * layout%weight(9) - 0.5_rk  &
        &       * (pdfTmp(9) + pdfTmp(8)))
      pm(m, 9) = omega(9) * (rho * layout%weight(9) * 3.0_rk   &
        &       * (u_fluid(2) - u_fluid(3)) - 0.5_rk * (pdfTmp(9) - pdfTmp(8)))


      outstate(( ielem-1)* varsys%nscalars+4) = pdfTmp(4) + pm(p, 1) + pm(m, 1)
      outstate(( ielem-1)* varsys%nscalars+1) = pdfTmp(1) + pm(p, 1) - pm(m, 1)

      outstate(( ielem-1)* varsys%nscalars+5) = pdfTmp(5) + pm(p, 2) + pm(m, 2)
      outstate(( ielem-1)* varsys%nscalars+2) = pdfTmp(2) + pm(p, 2) - pm(m, 2)

      outstate(( ielem-1)* varsys%nscalars+6) = pdfTmp(6) + pm(p, 3) + pm(m, 3)
      outstate(( ielem-1)* varsys%nscalars+3) = pdfTmp(3) + pm(p, 3) - pm(m, 3)

      outstate(( ielem-1)* varsys%nscalars+18) = pdfTmp(18) + pm(p, 4) + pm(m, 4)
      outstate(( ielem-1)* varsys%nscalars+15) = pdfTmp(15) + pm(p, 4) - pm(m, 4)

      outstate(( ielem-1)* varsys%nscalars+17) = pdfTmp(17) + pm(p, 5) + pm(m, 5)
      outstate(( ielem-1)* varsys%nscalars+16) = pdfTmp(16) + pm(p, 5) - pm(m, 5)

      outstate(( ielem-1)* varsys%nscalars+14) = pdfTmp(14) + pm(p, 6) + pm(m, 6)
      outstate(( ielem-1)* varsys%nscalars+11) = pdfTmp(11) + pm(p, 6) - pm(m, 6)

      outstate(( ielem-1)* varsys%nscalars+12) = pdfTmp(12) + pm(p, 7) + pm(m, 7)
      outstate(( ielem-1)* varsys%nscalars+13) = pdfTmp(13) + pm(p, 7) - pm(m, 7)

      outstate(( ielem-1)* varsys%nscalars+10) = pdfTmp(10) + pm(p, 8) + pm(m, 8)
      outstate(( ielem-1)* varsys%nscalars+7) = pdfTmp(7) + pm(p, 8) - pm(m, 8)

      outstate(( ielem-1)* varsys%nscalars+9) = pdfTmp(9) + pm(p, 9) + pm(m, 9)
      outstate(( ielem-1)* varsys%nscalars+8) = pdfTmp(8) + pm(p, 9) - pm(m, 9)

      outstate(( ielem-1)* varsys%nscalars+19) = pdfTmp(19) + aux_omega &
        &                                   * (rho * layout%weight(19) - pdfTmp(19))

    end do elemloop

    end subroutine mus_advRel_kPS_rTRT_vLmodel_lD3Q19
    ! ****************************************************************************** !


! **************************************************************************** !
  !> Advection relaxation routine based on the E-Model for the D3Q19
  !! Lattice Boltzmann model for the generic advection and anisotropic-dispersion
  !! equation. The E-Model is corrected to account for the effect of numerical
  !! diffusion with the nonlinear equilibrium correction
  !! \[ f_i^{eq,cor} = f_i^{eq} + \frac{w_i}{2 c_s^4} (c_{i\alpha} c_{i\beta}
  !!      - c_s^2 \delta_{\alpha\beta}) u_\alpha u_\beta \]
  !!
  !!   Irina Ginzburg (2005), "Equilibrium-type and link-type lattice Boltzmann
  !!   models for generic advection and anisotropic-dispersion equation",
  !!   Advances in Water Resources, Volume 28, Issue 11
  !!
  !! The MRT model is based on:
  !!
  !!   Dominique d’Humières. “Multiple–relaxation–time lattice Boltzmann models
  !!   in three dimensions”. In: Philosophical Transactions of the Royal Society
  !!   A: Mathematical, Physical and Engineering Sciences 360.1792 (2002),
  !!   pp. 437–451. doi: 10.1098/rsta.2001.0955
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kPS_rMRT_vEmodelCorr_lD3Q19( fieldProp, inState, outState, &
    &                            auxField, neigh, nElems, nSolve, level, layout,   &
    &                            params, varSys, derVarPos               )
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
      integer :: iElem, iDir
      type(mus_varSys_data_type), pointer :: fPtr
      type(mus_scheme_type), pointer :: scheme
      real(kind=rk) :: rho
      real(kind=rk) :: d_omega, nu_q, d_lambda, aux_omega
      real(kind=rk) :: transVel( nSolve*3 ) ! velocity from the transport field
      real(kind=rk) :: usq
      real(kind=rk) :: a_e, a_xx, a_ww, a_xy, a_xz, a_yz
      integer :: vel_varPos ! position of transport velocity variable in varSys
      real(kind=rk) :: inv_vel, u_fluid(3)
      real(kind=rk) :: s_1, s_2, s_4, s_9, s_10, s_13, s_16 ! relaxation parameters
      real(kind=rk) :: s_mrt( layout%fStencil%QQ ) ! relaxation parameters for MRT
      real(kind=rk) :: meq( layout%fStencil%QQ ) ! equilibrium moments
      real(kind=rk) :: mneq( layout%fStencil%QQ ) ! non-equilibrium moments of the pdf
      real(kind=rk) :: sum_2x, sum_2y, sum_2z ! sums of the 2nd order moments
      real(kind=rk) :: sum_4x, sum_4y, sum_4z ! sums of the 4th order moments
      real(kind=rk) :: sum_6, sum_12
      real(kind=rk) :: t2_1, t2_2, t4_1, t4_2, t4_3, t6_1, t6_2, t8_1, t8_2, &
                        t8_3, t8_4, t8_5, t8_6, t10_1, t10_2, t12_1, t12_2
      real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
        &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
        &              f000 ! local pdf values
      ! logical :: switch_neg, switch_print
      ! --------------------------------------------------------------------------
      ! access scheme via 1st variable method data which is a state variable
      call C_F_POINTER( varSys%method%val(derVarPos(1)%pdf)%method_Data, fPtr )
      scheme => fPtr%solverData%scheme
      ! switch_neg = .true.
      ! passive scalar has only one transport Variable
      vel_varPos = scheme%transVar%method(1)%data_varPos
      ! Get velocity field
      call varSys%method%val(vel_varPos)%get_valOfIndex( &
        & varSys  = varSys,                              &
        & time    = params%general%simControl%now,       &
        & iLevel  = level,                               &
        & idx     = scheme%transVar%method(1)            &
        &           %pntIndex%indexLvl(level)            &
        &           %val(1:nSolve),                      &
        & nVals   = nSolve,                              &
        & res     = transVel                             )

      ! convert physical velocity into LB velocity
      inv_vel = 1.0_rk / params%physics%fac( level )%vel
      transVel = transVel * inv_vel

      d_omega = fieldProp(1)%species%omega

      ! reciprocal of the free parameter nu, i.e. nu_q = 1/nu
      nu_q = cs2inv / (1.0_rk / d_omega - 0.5_rk)
      a_e = ( fieldProp(1)%species%diff_tensor(Dxx)  &
        &    + fieldProp(1)%species%diff_tensor(Dyy) &
        &    + fieldProp(1)%species%diff_tensor(Dzz) ) * nu_q / 3.0_rk - 1.0_rk
      a_xx = 2._rk / 3._rk * nu_q * ( fieldProp(1)%species%diff_tensor(Dxx)               &
        &                            - 0.5_rk * (fieldProp(1)%species%diff_tensor(Dyy)    &
        &                                        + fieldProp(1)%species%diff_tensor(Dzz)) &
        &                           )
      a_ww = 0.5_rk * nu_q * ( fieldProp(1)%species%diff_tensor(Dyy) &
        &                     - fieldProp(1)%species%diff_tensor(Dzz) )

      ! attention: D_ab has a multiplier of 2 in ADE
      a_xy = fieldProp(1)%species%diff_tensor(Dxy) * nu_q
      a_xz = fieldProp(1)%species%diff_tensor(Dxz) * nu_q
      a_yz = fieldProp(1)%species%diff_tensor(Dyz) * nu_q

      d_lambda = 0.25_rk
      aux_omega = 1.0_rk / (d_lambda / (1.0_rk / d_omega - 0.5_rk) + 0.5_rk)
      s_1 = 1.2_rk
      s_2 = aux_omega
      s_4 = d_omega
      s_9 = aux_omega
      s_10 = s_2
      s_13 = s_9
      s_16 = d_omega

      s_mrt(1) = 0.0_rk
      s_mrt(2) = s_1 / 2394._rk
      s_mrt(3) = s_2 / 252._rk
      s_mrt(4) = d_omega / 10._rk
      s_mrt(5) = s_4 / 40._rk
      s_mrt(6) = d_omega / 10._rk
      s_mrt(7) = s_4 / 40._rk
      s_mrt(8) = d_omega / 10._rk
      s_mrt(9) = s_4 / 40._rk
      s_mrt(10) = s_9 / 36._rk
      s_mrt(11) = s_10 / 72._rk
      s_mrt(12) = s_9 / 12._rk
      s_mrt(13) = s_10 / 24._rk
      s_mrt(14) = s_13 / 4._rk
      s_mrt(15) = s_13 / 4._rk
      s_mrt(16) = s_13 / 4._rk
      s_mrt(17) = s_16 / 8._rk
      s_mrt(18) = s_16 / 8._rk
      s_mrt(19) = s_16 / 8._rk

      elemloop: do iElem = 1, nSolve

        !> First load all local values into temp array
        fN00 = inState(neigh (( qn00-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f0N0 = inState(neigh (( q0n0-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f00N = inState(neigh (( q00n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f100 = inState(neigh (( q100-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f010 = inState(neigh (( q010-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f001 = inState(neigh (( q001-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f0NN = inState(neigh (( q0nn-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f0N1 = inState(neigh (( q0n1-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f01N = inState(neigh (( q01n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f011 = inState(neigh (( q011-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        fN0N = inState(neigh (( qn0n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f10N = inState(neigh (( q10n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        fN01 = inState(neigh (( qn01-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f101 = inState(neigh (( q101-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        fNN0 = inState(neigh (( qnn0-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        fN10 = inState(neigh (( qn10-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f1N0 = inState(neigh (( q1n0-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f110 = inState(neigh (( q110-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
        f000 = inState(neigh (( q000-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)

        ! x-, y- and z-velocity from transport field
        u_fluid = transVel( (iElem-1)*3+1 : iElem*3 )
        usq = u_fluid(1) * u_fluid(1) + u_fluid(2) * u_fluid(2) + u_fluid(3) * u_fluid(3)

        rho = 0.0_rk
        do iDir = 1, layout%fStencil%QQ
          rho = rho + instate( neigh( (idir - 1) * nelems + ielem ) )
        end do

        ! compute the equilibrium moments of the pdf
        meq = 0._rk
        meq(2) = rho * (19._rk * (usq + a_e) - 11._rk)
        meq(3) = 3._rk * rho - 11._rk / 2._rk * rho * usq
        meq(4) = rho * u_fluid(1)
        meq(5) = -2._rk / 3._rk * meq(4)
        meq(6) = rho * u_fluid(2)
        meq(7) = -2._rk / 3._rk * meq(6)
        meq(8) = rho * u_fluid(3)
        meq(9) = -2._rk / 3._rk * meq(8)
        meq(10) = rho * (2._rk * u_fluid(1)*u_fluid(1) - u_fluid(2)*u_fluid(2) - &
          & u_fluid(3)*u_fluid(3) + a_xx)
        meq(11) = -meq(10) / 2._rk
        meq(12) = rho * (u_fluid(2)*u_fluid(2) - u_fluid(3)*u_fluid(3) + &
          & a_ww * 2._rk / 3._rk)
        meq(13) = -meq(12) / 2._rk
        meq(14) = rho * (u_fluid(1)*u_fluid(2) + &
          & a_xy / 3._rk)
        meq(15) = rho * (u_fluid(2)*u_fluid(3) + &
          & a_yz / 3._rk)
        meq(16) = rho * (u_fluid(1)*u_fluid(3) + &
          & a_xz / 3._rk)

        ! compute the non-equilbrium moments of the pdf
        mneq = 0._rk
        sum_2x = fN00 + f100
        sum_2y = f0N0 + f010
        sum_2z = f00N + f001
        sum_6 =  sum_2x + sum_2y + sum_2z
        sum_4x = f0NN + f0N1 + f01N + f011
        sum_4y = fN0N + f10N + fN01 + f101
        sum_4z = fNN0 + fN10 + f1N0 + f110
        sum_12 = sum_4x + sum_4y + sum_4z
        mneq(2) = -30._rk * f000 - 11._rk * sum_6 + 8._rk * sum_12 - meq(2)
        mneq(3) = 12._rk * f000 - 4._rk * sum_6 + sum_12 - meq(3)
        t4_1 = f100 - fN00
        t4_2 = f10N - fN0N + f101 - fN01 + &
          & f110 - fN10 + f1N0 - fNN0
        mneq(4) = t4_1 + t4_2 - meq(4)
        mneq(5) = -4.0_rk * t4_1 + t4_2 - meq(5)
        t6_1 = f010 - f0N0
        t6_2 = f01N - f0NN + f011 - f0N1 + &
          & f110 - f1N0 + fN10 - fNN0
        mneq(6) = t6_1 + t6_2 - meq(6)
        mneq(7) = -4.0_rk * t6_1 + t6_2 - meq(7)
        t8_1 = f001 - f00N
        t8_2 = f0N1 - f0NN + f011 - f01N + &
          & f101 - f10N + fN01 - fN0N
        mneq(8) = t8_1 + t8_2 - meq(8)
        mneq(9) = -4.0_rk * t8_1 + t8_2 - meq(9)
        t10_1 = 2.0_rk * sum_2x - sum_2y - sum_2z
        t10_2 = sum_4y + sum_4z - 2.0_rk * sum_4x
        mneq(10) = t10_1 + t10_2 - meq(10)
        mneq(11) = -2.0_rk * t10_1 + t10_2 - meq(11)
        t12_1 = sum_2y - sum_2z
        t12_2 = sum_4z - sum_4y
        mneq(12) = t12_1 + t12_2 - meq(12)
        mneq(13) = -2.0_rk * t12_1 + t12_2 - meq(13)
        mneq(14) = f110 + fNN0 - f1N0 - fN10 - meq(14)
        mneq(15) = f011 + f0NN - f01N - f0N1 - meq(15)
        mneq(16) = fN0N + f101 - f10N - fN01 - meq(16)
        mneq(17) = f110 - fNN0 + f1N0 - fN10  + &
          & fN0N - f101 - f10N + fN01 - meq(17)
        mneq(18) = -f110 + fNN0 + f1N0 - fN10  + &
          & f011 - f0NN + f01N - f0N1 - meq(18)
        mneq(19) = -fN0N + f101 - f10N + fN01 - &
          & f011 + f0NN + f01N - f0N1 - meq(19)

        mneq = -mneq * s_mrt

        outstate( (ielem-1)*qq+ q000+(1-1)*qq) = &
          & f000 - 30.0_rk*mneq(2) + 12.0_rk*mneq(3)


        t2_1 = 2.0_rk*mneq(10) - 4.0_rk*mneq(11) - &
          & 11.0_rk * mneq(2) - 4.0_rk * mneq(3)
        t2_2 = mneq(4) - 4.0_rk*mneq(5)
        outstate( (ielem-1)*qq+  q100+(1-1)*qq) = &
          & f100 + t2_1 + t2_2
        outstate( (ielem-1)*qq+  qn00+(1-1)*qq) = &
          & fN00 + t2_1 - t2_2


        t4_1 = - 11.0_rk*mneq(2) - 4.0_rk*mneq(3) - &
          & mneq(10) + 2.0_rk*mneq(11)
        t4_2 = mneq(6) - 4.0_rk*mneq(7)
        t4_3 = mneq(12) - 2.0_rk*mneq(13)
        outstate( (ielem-1)*qq+  q010+(1-1)*qq) = &
          & f010 + t4_1 + t4_2 + t4_3
        outstate( (ielem-1)*qq+  q0n0+(1-1)*qq) = &
          & f0N0 + t4_1 - t4_2 + t4_3
        t6_2 = mneq(8) - 4.0_rk*mneq(9)
        outstate( (ielem-1)*qq+  q001+(1-1)*qq) = &
          & f001 + t4_1 + t6_2 - t4_3
        outstate( (ielem-1)*qq+  q00n+(1-1)*qq) = &
          & f00N + t4_1 - t6_2 - t4_3


        t8_1 = 8.0_rk*mneq(2) + mneq(3)
        t8_2 = mneq(4) + mneq(5)
        t8_3 = mneq(6) + mneq(7)
        t8_4 = mneq(10) + mneq(11)
        t8_5 = mneq(12) + mneq(13)
        t8_6 = mneq(8) + mneq(9)
        outstate( (ielem-1)*qq+ q110+(1-1)*qq) = &
          & f110 + t8_1 + t8_2 + t8_3 + t8_4 + t8_5 + mneq(14) + mneq(17) - mneq(18)
        outstate( (ielem-1)*qq+ qn10+(1-1)*qq) = &
          & fN10 + t8_1 - t8_2 + t8_3 + t8_4 + t8_5 - mneq(14) - mneq(17) - mneq(18)
        outstate( (ielem-1)*qq+ q1n0+(1-1)*qq) = &
          & f1N0 + t8_1 + t8_2 - t8_3 + t8_4 + t8_5 - mneq(14) + mneq(17) + mneq(18)
        outstate( (ielem-1)*qq+ qnn0+(1-1)*qq) = &
          &fNN0 + t8_1 - t8_2 - t8_3 + t8_4 + t8_5 + mneq(14) - mneq(17) + mneq(18)
        outstate( (ielem-1)*qq+ q101+(1-1)*qq) = &
          & f101 + t8_1 + t8_2 + t8_6 + t8_4 - t8_5 + mneq(16) - mneq(17) + mneq(19)
        outstate( (ielem-1)*qq+ qn01+(1-1)*qq) = &
          & fN01 + t8_1 - t8_2 + t8_6 + t8_4 - t8_5 - mneq(16) + mneq(17) + mneq(19)
        outstate( (ielem-1)*qq+ q10n+(1-1)*qq) = &
          & f10N + t8_1 + t8_2 - t8_6 + t8_4 - t8_5 - mneq(16) - mneq(17) - mneq(19)
        outstate( (ielem-1)*qq+ qn0n+(1-1)*qq) = &
          & fN0N + t8_1 - t8_2 - t8_6 + t8_4 - t8_5 + mneq(16) + mneq(17) - mneq(19)
        outstate( (ielem-1)*qq+ q011+(1-1)*qq) = &
          & f011 + t8_1 + t8_3 + t8_6 - 2.0_rk * t8_4 + mneq(15) + mneq(18) - mneq(19)
        outstate( (ielem-1)*qq+  q0n1+(1-1)*qq) = &
          & f0N1 + t8_1 - t8_3 + t8_6 - 2.0_rk * t8_4 - mneq(15) - mneq(18) - mneq(19)
        outstate( (ielem-1)*qq+  q01n+(1-1)*qq) = &
          & f01N + t8_1 + t8_3 - t8_6 - 2.0_rk * t8_4 - mneq(15) + mneq(18) + mneq(19)
        outstate( (ielem-1)*qq+  q0nn+(1-1)*qq) = &
          & f0NN + t8_1 - t8_3 - t8_6 - 2.0_rk * t8_4 + mneq(15) - mneq(18) + mneq(19)

      end do elemloop

      end subroutine mus_advRel_kPS_rMRT_vEmodelCorr_lD3Q19
    ! ****************************************************************************** !


end module mus_d3q19_module
! **************************************************************************** !
