! Copyright (c) 2013, 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
! ************************************************************************** !
!> author: Kannan Masilamani
!! author: Jiaxing Qi
!! This module provides the MUSUBI specific functions for calculating
!! macroscopic quantities from the state variables for incompressible LBM
!! models.\n
!! Notice that only those quantities that related to density should have a
!! formula which differs from normal LBM model.
!!
!! The depending common interface between MUSUBI and ATELES is defined in the
!! tem_derived_module. The functionality for accessing a variable from the state
!! and evaluating a lua function are also provided in the tem_derived module.
!!
!! Do not use get_Element or get_Point routines to update the state !
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
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015, 2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!--------------------------------------------
!    A O S - Array of structures layout new
!-------------------------------------------
! Access to get_point value output
! Access to get_element value output
module mus_derQuanIncomp_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use tem_param_module,         only: div1_2, div1_3, div1_54, div1_9, div3_4, &
    &                                 div1_36, div3_4h,                        &
    &                                 sqrt3, cs2inv, cs2, t2cs4inv, t2cs2inv,  &
    &                                 cs4inv, rho0, rho0Inv, q000
  use env_module,               only: rk, long_k, labelLen
  use tem_float_module,         only: operator(.feq.), operator(.fge.), &
    &                                 operator(.fle.)
  use tem_varSys_module,        only: tem_varSys_type, tem_varSys_op_type,     &
    &                                 tem_varSys_append_derVar,                &
    &                                 tem_varSys_proc_point,                   &
    &                                 tem_varSys_proc_element,                 &
    &                                 tem_varSys_proc_setParams,               &
    &                                 tem_varSys_proc_getParams,               &
    &                                 tem_varSys_proc_setupIndices,            &
    &                                 tem_varSys_proc_getValOfIndex,           &
    &                                 tem_varSys_getPoint_dummy,               &
    &                                 tem_varSys_getElement_dummy,             &
    &                                 tem_varSys_setupIndices_dummy,           &
    &                                 tem_varSys_getValOfIndex_dummy,          &
    &                                 tem_varSys_setParams_dummy,              &
    &                                 tem_varSys_getParams_dummy
  use tem_variable_module,      only: tem_variable_type
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_topology_module,      only: tem_levelOf
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_subTree_type_module,  only: tem_subTree_type, tem_treeIDfrom_subTree
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit
  use tem_operation_var_module, only: tem_evalMag_forElement,      &
    &                                 tem_evalMag_forPoint,        &
    &                                 tem_evalMag_fromIndex,       &
    &                                 tem_evalDiff_forElement,     &
    &                                 tem_evalDiff_forPoint,       &
    &                                 tem_evalDiff_fromIndex,      &
    &                                 tem_opVar_setupIndices,      &
    &                                 tem_get_new_varSys_data_ptr, &
    &                                 tem_opVar_setParams,         &
    &                                 tem_opVar_getParams
  use tem_debug_module,         only: dbgUnit
  use tem_grow_array_module,    only: grw_labelarray_type, append


  ! include musubi modules
  use mus_source_type_module,        only: mus_source_op_type
  use mus_pdf_module,                only: pdf_data_type
  use mus_varSys_module,             only: mus_varSys_data_type,               &
    &                                      mus_varSys_solverData_type,         &
    &                                      mus_get_new_solver_ptr,             &
    &                                      mus_deriveVar_ForPoint,             &
    &                                      mus_generic_varFromPDF_fromIndex,   &
    &                                      mus_generic_fromPDF_forElement,     &
    &                                      mus_derive_fromPDF
  use mus_stateVar_module,           only: mus_accessVar_setupIndices,         &
    &                                      mus_stateVar_Fetch_fromIndex,       &
    &                                      mus_stateVar_Fetch_now_fromIndex,   &
    &                                      mus_access_stateFetch_ForElement,   &
    &                                      mus_access_stateFetch_now_ForElement
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_derivedQuantities_module2, only: secondMom_3D
  use mus_derQuan_module,            only: deriveDensity,                      &
    &                                      deriveDensity_fromIndex,            &
    &                                      derivePressure,                     &
    &                                      derivePressure_forPoint,            &
    &                                      derivePressure_fromIndex,           &
    &                                      deriveShearStress,                  &
    &                                      deriveShearMag,                     &
    &                                      deriveWSS2D,                        &
    &                                      deriveWSS3D,                        &
    &                                      deriveTemp,                         &
    &                                      deriveShearRate,                    &
    &                                      deriveBndForce,                     &
    &                                      deriveNonEquil,                     &
    &                                      deriveNonEquil_fromIndex,           &
    &                                      deriveKE, deriveKE_forPoint,        &
    &                                      deriveKe_fromIndex,                 &
    &                                      deriveStrainRate,                   &
    &                                      deriveStrainRate_fromIndex,         &
    &                                      deriveMomentum,                     &
    &                                      deriveMomentum_forPoint,            &
    &                                      deriveMomentum_fromIndex,           &
    &                                      deriveEquil, deriveEquil_fromIndex
  use mus_operation_var_module,      only: mus_opVar_setupIndices,         &
    &                                      mus_opVar_gradU_forElement,     &
    &                                      mus_opVar_vorticity_forElement, &
    &                                      mus_opVar_QCriterion_forElement
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_physics_module,            only: mus_convertFac_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none

  private

  public :: mus_append_derVar_fluidIncomp

  ! source variables
  public :: derive_absorbLayerIncomp
  public :: applySrc_absorbLayerIncomp

contains


  ! ************************************************************************ !
  !> subroutine to add derive variables for incompressible LBM
  !! (schemekind = 'fluid_incompressible') to the varsys.
  subroutine mus_append_derVar_fluidIncomp( varSys, solverData, schemeHeader, &
    &                                       stencil, fldLabel, derVarName     )
    ! -------------------------------------------------------------------- !
    !> global variable system
    type(tem_varSys_type), intent(inout)  :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in)  :: schemeHeader

    !> compute stencil defintion
    type(tem_stencilHeader_type), intent(in)       :: stencil

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel

    !> array of derive physical variables
    type(grw_labelarray_type), intent(inout) :: derVarName
    ! -------------------------------------------------------------------- !
    ! number of derive variables
    integer :: nDerVars, iVar, nComponents, addedPos, iIn
    logical :: wasAdded
    character(len=labelLen), allocatable ::  input_varname(:)
    character(len=labelLen)  ::  varName
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    type(c_ptr) :: method_data
    character(len=labelLen), allocatable :: derVarName_loc(:)
    ! -------------------------------------------------------------------- !
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    nDerVars = 20
    allocate(derVarName_loc(nDerVars))
    derVarName_loc = [ 'fetch_pdf         ', 'fetch_pdf_now     ', &
      &                'pressure          ', 'equilibrium       ', &
      &                'non_equilibrium   ', 'kinetic_energy    ', &
      &                'shear_stress      ', 'strain_rate       ', &
      &                'shear_rate        ', 'wss               ', &
      &                'momentum          ', 'vel_mag           ', &
      &                'bnd_force         ', 'fetch_pdf         ', &
      &                'shear_mag         ', 'temperature       ', &
      &                'grad_velocity     ', 'vorticity         ', &
      &                'q_criterion       ', 'pressure_deviation'  ]

    do iVar = 1, nDerVars
      call append(derVarName, derVarName_loc(iVar))
      ! set default pointers, overwrite if neccessary
      get_element => tem_varSys_getElement_dummy
      get_point => mus_deriveVar_ForPoint
      setup_indices => mus_opVar_setupIndices
      get_valOfIndex => tem_varSys_getvalOfIndex_dummy
      set_params => tem_varSys_setParams_dummy
      get_params => tem_varSys_getParams_dummy
      method_data  = mus_get_new_solver_ptr(solverData)

      select case(trim(adjustl(derVarName_loc(iVar))))
      case ('fetch_pdf')
        get_element => mus_access_stateFetch_ForElement
        get_valOfIndex => mus_stateVar_Fetch_fromIndex
        setup_indices => mus_accessVar_setupIndices
        nComponents = stencil%QQ
        allocate(input_varname(1))
        input_varname(1) = 'pdf'

      case ('fetch_pdf_now')
        get_element => mus_access_stateFetch_now_ForElement
        get_valOfIndex => mus_stateVar_Fetch_now_fromIndex
        setup_indices => mus_accessVar_setupIndices
        nComponents = stencil%QQ
        allocate(input_varname(1))
        input_varname(1) = 'pdf'

      case ('pressure')
        get_element => derivePressure
        get_point => derivePressure_forPoint
        get_valOfIndex => derivePressure_fromIndex
        nComponents = 1
        allocate(input_varname(1))
        input_varname(1) = 'density'

      case ('bnd_force')
        get_element => deriveBndForce
        nComponents = 3
        allocate(input_varname(1))
        input_varname(1) = 'pdf'

      case ('equilibrium')
        get_element => deriveEquil
        get_valOfIndex => deriveEquil_fromIndex
        nComponents = stencil%QQ
        allocate(input_varname(1))
        input_varname(1) = 'pdf'

      case ('non_equilibrium')
        get_element => deriveNonEquil
        get_valOfIndex => deriveNonEquil_fromIndex
        nComponents = stencil%QQ
        allocate(input_varname(1))
        input_varname(1) = 'fetch_pdf_now'

      case ('kinetic_energy')
        get_element => deriveKe
        get_point => deriveKe_forPoint
        get_ValOfIndex => deriveKE_fromIndex
        nComponents = 1
        allocate(input_varname(2))
        ! because we use the same routine of compressible and that one uses density!
        input_varname(1) = 'density'
        input_varname(2) = 'velocity'

      case ('temperature')
        get_element => deriveTemp
        nComponents = 1
        allocate(input_varname(0))

      case ('shear_stress')
        nComponents = 6
        get_element => deriveShearStress
        allocate(input_varname(1))
        input_varname(1) = 'non_equilibrium'

      case ('strain_rate')
        nComponents = 6
        get_element => deriveStrainRate
        get_ValOfIndex => deriveStrainRate_fromIndex
        allocate(input_varname(1))
        input_varname(1) = 'fetch_pdf_now'

      case ('shear_rate')
        get_element => deriveShearRate
        nComponents = 1
        allocate(input_varname(1))
        input_varname(1) = 'strain_rate'

      case ('wss')
        nComponents = 1
        allocate(input_varname(1))
        input_varname(1) = 'shear_stress'
        if (stencil%nDims == 2) then
          get_element => deriveWSS2D
        else if (stencil%nDims == 3) then
          get_element => deriveWSS3D
        else
          write(logUnit(1),*) 'WARNING: WSS does not support 1D'
        end if

      case ('momentum')
        get_element => deriveMomentum
        get_point => deriveMomentum_forPoint
        get_valOfIndex => deriveMomentum_fromIndex
        nComponents = 3
        allocate(input_varname(2))
        ! because we use the same function of compressible which uses density!
        input_varname(1) = 'density'
        input_varname(2) = 'velocity'

      case ('grad_velocity')
        get_element => mus_opVar_gradU_forElement
        nComponents = 9
        allocate(input_varname(1))
        input_varname(1) = 'velocity'

      case ('vorticity')
        get_element => mus_opVar_vorticity_forElement
        nComponents = 3
        allocate(input_varname(1))
        input_varname(1) = 'velocity'

      case ('q_criterion')
        get_element => mus_opVar_QCriterion_forElement
        nComponents = 1
        allocate(input_varname(1))
        input_varname(1) = 'velocity'

      case ('vel_mag')
        get_element => tem_evalMag_forElement
        get_point => tem_evalMag_forPoint
        get_valOfIndex => tem_evalMag_fromIndex
        setup_indices => tem_opVar_setupIndices
        set_params => tem_opVar_setParams
        get_params => tem_opVar_getParams
        method_data = tem_get_new_varSys_data_ptr(method_data)
        nComponents = 1
        allocate(input_varname(1))
        input_varname(1) = 'velocity'

      case ('shear_mag')
        get_element => deriveShearMag
        nComponents = 1
        allocate(input_varname(1))
        input_varname(1) = 'shear_stress'

      case ('pressure_deviation')
        get_element => tem_evalDiff_forElement
        get_point => tem_evalDiff_forPoint
        get_valOfIndex => tem_evalDiff_fromIndex
        setup_indices => tem_opVar_setupIndices
        set_params => tem_opVar_setParams
        get_params => tem_opVar_getParams
        method_data = tem_get_new_varSys_data_ptr(method_data)
        nComponents = 1
        allocate(input_varname(2))
        input_varname(1) = 'pressure'
        input_varname(2) = 'pressure_reference'

      case default
        write(logUnit(1),*) 'WARNING: Unknown variable: '//&
          &                 trim(derVarName_loc(iVar))
        cycle !go to next variable
      end select

      ! update variable names with field label
      varname = trim(fldLabel)//trim(adjustl(derVarName_loc(iVar)))
      do iIn = 1, size(input_varname)
        input_varname(iIn) = trim(fldLabel)//trim(input_varname(iIn))
      end do

      ! append variable to varSys
      call tem_varSys_append_derVar( me             = varSys,         &
        &                            varName        = trim(varname),  &
        &                            nComponents    = nComponents,    &
        &                            input_varname  = input_varname,  &
        &                            method_data    = method_data,    &
        &                            get_point      = get_point,      &
        &                            get_element    = get_element,    &
        &                            set_params     = set_params,     &
        &                            get_params     = get_params,     &
        &                            setup_indices  = setup_indices,  &
        &                            get_valOfIndex = get_valOfIndex, &
        &                            pos            = addedPos,       &
        &                            wasAdded       = wasAdded        )

      if (wasAdded) then
        write(logUnit(10),*) ' Appended variable: '//trim(varname)
      else if (addedpos < 1) then
        write(logUnit(1),*) 'Error: variable '//trim(varname)// &
          &                 ' is not added to variable system'
      end if

      deallocate(input_varname)
    end do

  end subroutine mus_append_derVar_fluidIncomp
  ! ************************************************************************ !


! **************************************************************************** !
!        Subroutines with common interface for the function pointers           !
! **************************************************************************** !


! **************************************************************************** !
!       Subroutines with common interface for the function pointers            !
!                              getValOfIndex                                   !
! **************************************************************************** !

! **************************************************************************** !
!                           Calculation routines                               !
! **************************************************************************** !

! **************************************************************************** !
   !> Derive absorb layer variable defined as a source term.
  recursive subroutine derive_absorbLayerIncomp(fun, varsys, elempos, time, &
    &                                           tree, nElems, nDofs, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !

    call tem_abort('Not implemented')

  end subroutine derive_absorbLayerIncomp
! **************************************************************************** !


! **************************************************************************** !
  !> Update state with source variable "absorb_layer".
  !! absorb_layer is used to absorb the flow and gradually reduce the flow
  !! quantities like pressure and velocity to a fixed value for incompressible
  !! model. It is based on:
  !! Xu, H., & Sagaut, P. (2013). Analysis of the absorbing layers for the
  !! weakly-compressible lattice Boltzmann methods. Journal of Computational
  !! Physics, 245(x), 14–42.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_absorbLayerIncomp( fun, inState, outState, neigh,      &
    &                                    auxField, nPdfSize, iLevel, varSys, &
    &                                    time, phyConvFac, derVarPos         )
    ! -------------------------------------------------------------------- !
    !> Description of method to apply source terms
    class(mus_source_op_type), intent(in) :: fun

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    call tem_abort('Error: Absorb layer is not yet implemented')
  end subroutine applySrc_absorbLayerIncomp
! **************************************************************************** !


end module mus_derQuanIncomp_module
! **************************************************************************** !
