! Copyright (c) 2012, 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> author: Simon Zimny
!! author: Kartik Jain
!! This module contains the type definition(s) required in the scheme routines
!! In addition depend type and condition type for geometry increase routine
!! are defined. Compute kernel definition is also defined in this module
!!
module mus_scheme_type_module

  ! include treelm modules
  use env_module,                only: rk
  use tem_tracking_module,       only: tem_tracking_type
  use tem_varSys_module,         only: tem_varSys_type
  use tem_variable_module,       only: tem_variable_type
  use tem_varMap_module,         only: tem_varMap_type,          &
    &                                  tem_possible_variable_type
  use tem_spacetime_fun_module,  only: tem_st_fun_linkedList_type
  use tem_construction_module,   only: tem_levelDesc_type

  ! include musubi modules
  use mus_bc_header_module,          only: glob_boundary_type
  use mus_source_type_module,        only: mus_source_type
  use mus_field_module,              only: mus_field_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_pdf_module,                only: pdf_data_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_param_module,              only: mus_param_type
  use mus_mixture_module,            only: mus_mixture_type
  use mus_interpolate_header_module, only: mus_interpolation_type
  use mus_transport_var_module,      only: mus_transport_var_type
  use mus_nernstPlanck_module,       only: mus_nernstPlanck_type
  use mus_auxField_module,           only: mus_auxFieldVar_type, &
    &                                      mus_proc_calcAuxField
  use mus_gradData_module,           only: mus_gradData_type, mus_Grad_type

  implicit none
  private

  public :: mus_scheme_type
  public :: kernel
  public :: array2D_type

  type array2D_type
    !! To allow Intel AVX SIMD streaming store instructions,
    !! the array must be aligned at 32 bytes
    real(kind=rk), allocatable, dimension(:,:) :: val
    !dir$ attributes align : 32 :: val
  end type array2D_type

  !> Datatype containing all information on the scheme.
  !!
  !! The mus_scheme_type contains of all information that are needed
  !! to run a simulation (including informations on the: fluid, boundary
  !! conditions, levelDescriptor, state vector, layout, diffusion info,
  !! tracking).
  type mus_scheme_type

    !> Interpolation description for each scheme to do its own interpolation
    type(mus_interpolation_type) :: intp

    !> contains mixture information for multispecies
    type(mus_mixture_type) :: mixture

    !. Contains information for nernst-planck
    type(mus_nernstPlanck_type) :: nernstPlanck

    !> number of fields in the current scheme
    integer :: nFields = 0
    !> array of field type for each field
    type(mus_field_type), allocatable :: field(:)

    !> array of boundary types contains elems of each boundary
    type(glob_boundary_type), allocatable :: globBC(:)

    !> global source applied to all fields
    type(mus_source_type) :: globSrc

    !> possible source variables depends on scheme kind
    type(tem_possible_variable_type) :: poss_srcVar

    !> transport variables
    type(mus_transport_var_type) :: transVar

    !> possible transport variables depends on scheme kind
    !! This variables might be used in compute kernel
    type(tem_possible_variable_type) :: poss_transVar

    !> identifier of the scheme
    type(mus_scheme_header_type) :: header

    type(tem_levelDesc_type), allocatable :: levelDesc(:)
    !> pdf_data_types for every level
    !! size: minLevel:maxLevel
    type(pdf_data_type), allocatable :: pdf(:)

    !> Data vector containing the pdf state
    !! allocated in routine: mus_construct
    !! size: minLevel:maxLevel
    type( array2D_type ), allocatable :: state(:)

    !> the scheme representation used in this scheme
    type(mus_scheme_layout_type) :: layout

    !> function pointer to compute kernel
    procedure( kernel ), pointer, nopass :: compute => null()

    !> Contains trackingControl, config and instances
    type( tem_tracking_type ) :: track

    !> Position of reduction transient variable in varSys
    type(tem_varMap_type) :: redTransVarMap

    !> store position of derived variable each field and total field
    !! in the global system
    type(mus_derVarPos_type), allocatable :: derVarPos(:)

    !> global variable system definition
    type(tem_varSys_type) :: varSys

    !> Variables defined in the lua file
    type(tem_variable_type), allocatable :: luaVar(:)

    !> state variable position in the global varSys
    type(tem_varMap_type) :: stateVarMap

    !> contains spacetime functions defined for lua variables
    type(tem_st_fun_linkedList_type) :: st_funList

    !> Used in mus_harvesting to check whether variables loaded from
    !! restart file has pdf variable
    logical :: readVarIsPdf

    !> stores auxField variable values and function pointer to compute
    !! auxiliary field
    !! Size: minlevel:maxLevel
    type(mus_auxFieldVar_type), allocatable :: auxField(:)

    !> Contains direct neighbor position in the state and
    !! finite difference coefficients to compute gradient
    type(mus_gradData_type), allocatable :: gradData(:)

    !> Function pointer to evaluate auxilary variable
    procedure(mus_proc_calcAuxField), pointer, nopass :: calcAuxField => null()

    !> Contains the different pointers to calculate the gradients
    type(mus_Grad_type) :: Grad

  end type mus_scheme_type

  !> What does the kernel interface look like?
  !! Every kernel's argument list must correspond to this one.
  !!
  !! Adhere to the below naming convection for the kernel names
  !!
  !! mus_advRel_k<kind>_r<relaxation>_v<variant>_l<layout>
  !!
  !! Examples:
  !!   mus_advRel_kFluid_rBGK_vStd_lD3Q19
  !!   mus_advRel_kFluidIncomp_rBGK_vStd_lD3Q19
  !!   mus_advRel_kFluid_rBGK_vHRR_lD3Q19
  !!
  !! For non-specific implementation leave names out.
  !! So we do not use generic keyword in the kernel names anymore.
  !! Examples:
  !!   mus_advRel_kFluid_rBGK_vStd_l
  !!   mus_advRel_kMsLiquid_rBGK_vStd_l
  !!
  abstract interface
    !> The common subroutine interface for compute kernels. All kernels have to
    !! implement this interface in order to be callable via
    !! mus_scheme_type%compute function pointer.
    subroutine kernel(fieldProp, inState, outState, auxField, neigh, nElems, &
      &               nSolve, level, layout, params, varSys, derVarPos       )
      import :: rk, mus_field_prop_type, mus_scheme_layout_type,  &
        &       mus_scheme_type, tem_varSys_type, mus_param_type, &
        &       mus_derVarPos_type
      ! ---------------------------------------------------------------- !
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
    end subroutine kernel
  end interface

end module mus_scheme_type_module
! ****************************************************************************** !
