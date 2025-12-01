! Copyright (c) 2013, 2016, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015-2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Philipp Otte <otte@mathcces.rwth-aachen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! ****************************************************************************** !
!> author: Kannan Masilamani
!! author: Jiaxing Qi
!! This module provides the MUSUBI specific functions for calculating
!! macroscopic quantities from the state variables.
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
module mus_derQuan_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use tem_param_module,         only: div1_2, div1_3, div1_54, div1_9, div3_4, &
    &                                 sqrt3, cs2inv, cs2, t2cs2inv, t2cs4inv,  &
    &                                 cs4inv, q000, rho0, csInv
  use env_module,               only: rk, long_k, labelLen
  use tem_float_module,         only: operator(.feq.), operator(.fge.), &
    &                                 operator(.fle.)
  use tem_variable_module,      only: tem_variable_type
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_topology_module,      only: tem_levelOf
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_logging_module,       only: logUnit
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
  use tem_aux_module,           only: tem_abort
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
  use tem_property_module,      only: prp_hasBnd, prp_hasQval
  use tem_tools_module,         only: tem_PositionInSorted
  use tem_debug_module,         only: dbgUnit
  use tem_grow_array_module,    only: grw_labelarray_type, append
  use tem_dyn_array_module,     only: PositionOfVal

  ! include musubi modules
  use mus_source_type_module,        only: mus_source_op_type
  use mus_pdf_module,                only: pdf_data_type
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_derivedQuantities_module2, only: secondMom_3D, getShearRate
  use mus_varSys_module,             only: mus_varSys_data_type,             &
    &                                      mus_varSys_solverData_type,       &
    &                                      mus_get_new_solver_ptr,           &
    &                                      mus_deriveVar_forPoint,           &
    &                                      mus_generic_varFromPDF_fromIndex, &
    &                                      mus_generic_fromPDF_forElement,   &
    &                                      mus_derive_fromPDF
  use mus_stateVar_module,           only: mus_accessVar_setupIndices,         &
    &                                      mus_stateVar_Fetch_fromIndex,       &
    &                                      mus_stateVar_Fetch_now_fromIndex,   &
    &                                      mus_access_stateFetch_ForElement,   &
    &                                      mus_access_stateFetch_now_ForElement
  use mus_operation_var_module,      only: mus_opVar_setupIndices,         &
    &                                      mus_opVar_gradU_forElement,     &
    &                                      mus_opVar_vorticity_forElement, &
    &                                      mus_opVar_QCriterion_forElement
  use mus_stateVar_module,           only: mus_accessVar_setupIndices
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_physics_module,            only: mus_convertFac_type
  use mus_hrrInit_module,            only: HRR_Correction_d2q9,  &
    &                                      HRR_Correction_d3q19, &
    &                                      HRR_Correction_d3q27
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none

  private

  public :: mus_append_derVar_fluid

  public :: derivePressure, derivePressure_fromIndex
  public :: derivePressure_forPoint
  public :: deriveDensity
  public :: deriveDensity_fromIndex
  public :: deriveShearStress
  public :: deriveShearMag
  public :: deriveWSS3D
  public :: deriveWSS2D
  public :: deriveTemp
  public :: deriveShearRate
  public :: deriveBndForce
  public :: deriveBndMoment
  ! equilbrium from macro uses different interface defined in
  ! mus_variable_module
  public :: deriveEquil_FromMacro
  public :: deriveRho_fromState
  public :: deriveVel_fromState
  public :: deriveVel_FromPreColState
  public :: deriveEq_fromState
  public :: deriveEquil_fromAux
  public :: deriveAux_fromState

  ! derive functions common also to incompressible fluids
  public :: deriveNonEquil
  public :: deriveNonEquil_fromIndex
  public :: deriveKE
  public :: deriveKE_forPoint
  public :: deriveKe_fromIndex
  public :: deriveStrainRate
  public :: deriveStrainRate_fromIndex
  public :: deriveMomentum
  public :: deriveMomentum_forPoint
  public :: deriveMomentum_fromIndex
  public :: deriveEquil
  public :: deriveEquil_fromIndex

  ! source variable
  public :: derive_absorbLayer
  public :: derive_force_MRT
  public :: derive_force1stOrd
  public :: derive_HRRCorrection_d2q9
  public :: derive_HRRCorrection_d3q19
  public :: derive_HRRCorrection_d3q27
  public :: derive_brinkmanForce
  public :: derive_brinkmanForce_TRT

  ! Apply source add source term to state
  public :: applySrc_absorbLayer
  public :: applySrc_absorbLayer_MRT
  public :: applySrc_absorbLayerDyn
  public :: applySrc_absorbLayerDyn_MRT
  public :: applySrc_force
  public :: applySrc_force_GNS
  public :: applySrc_force_MRT
  public :: applySrc_force_MRT_d2q9
  public :: applySrc_force_MRT_d3q19
  public :: applySrc_force_MRT_d3q27
  public :: applySrc_turbChanForce
  public :: applySrc_turbChanForce_MRT
  public :: applySrc_turbChanForce_MRT_d2q9
  public :: applySrc_turbChanForce_MRT_d3q19
  public :: applySrc_turbChanForce_MRT_d3q27
  public :: applySrc_force1stOrd
  public :: applySrc_brinkmanForce
  public :: applySrc_brinkmanForce_TRT

contains

  ! **************************************************************************** !
  !> subroutine to add derive variables for weakly compressible LBM
  !! (schemekind = 'fluid') to the varsys.
  subroutine mus_append_derVar_fluid( varSys, solverData, schemeHeader, &
    &                                 stencil, fldLabel, derVarName     )
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)  :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in)  :: schemeHeader

    !> compute stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel

    !> array of derive physical variables
    type(grw_labelarray_type), intent(inout) :: derVarName
    ! ---------------------------------------------------------------------------
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
    ! ---------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    nDerVars = 21
    allocate(derVarName_loc(nDerVars))
    derVarName_loc    = [ 'fetch_pdf         ', 'fetch_pdf_now     ', &
      &                   'pressure          ', 'equilibrium       ', &
      &                   'non_equilibrium   ', 'kinetic_energy    ', &
      &                   'shear_stress      ', 'strain_rate       ', &
      &                   'shear_rate        ', 'wss               ', &
      &                   'momentum          ', 'vel_mag           ', &
      &                   'bnd_force         ', 'moment            ', &
      &                   'shear_mag         ', 'temperature       ', &
      &                   'grad_velocity     ', 'vorticity         ', &
      &                   'q_criterion       ', 'mach              ', &
      &                   'pressure_deviation' ]

    do iVar = 1, nDerVars

      call append(derVarName, derVarName_loc(iVar))

      ! set default pointers, overwrite if neccessary
      get_element => tem_varSys_getElement_dummy
      get_point => mus_deriveVar_forPoint
      setup_indices => mus_opVar_setupIndices
      get_valOfIndex => tem_varSys_getValOfIndex_dummy
      method_data  = mus_get_new_solver_ptr(solverData)
      set_params => tem_varSys_setParams_dummy
      get_params => tem_varSys_getParams_dummy

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

      case ('mach')
        get_element => deriveMachNr
        get_point => deriveMachNr_forPoint
        get_valOfIndex => deriveMachNr_fromIndex
        nComponents = 1
        allocate(input_varname(1))
        input_varname(1) = 'vel_mag'

      case ('bnd_force')
        get_element => deriveBndForce
        nComponents = 3
        allocate(input_varname(1))
        input_varname(1) = 'pdf'

      case ('bnd_moment')
        get_element => deriveBndMoment
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
        get_element => deriveKE
        get_point => deriveKE_forPoint
        get_valOfIndex => deriveKe_fromIndex
        nComponents = 1
        allocate(input_varname(2))
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
        get_valOfIndex => deriveStrainRate_fromIndex
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

      case ('moment') !ONLY FOR 2D
        get_element => deriveMoment
        nComponents = 9
        allocate(input_varname(1))
        input_varname(1) = 'pdf'

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
      call tem_varSys_append_derVar(  me             = varSys,         &
        &                             varName        = trim(varname),  &
        &                             nComponents    = nComponents,    &
        &                             input_varname  = input_varname,  &
        &                             method_data    = method_data,    &
        &                             get_point      = get_point,      &
        &                             get_element    = get_element,    &
        &                             set_params     = set_params,     &
        &                             get_params     = get_params,     &
        &                             setup_indices  = setup_indices,  &
        &                             get_valOfIndex = get_valOfIndex, &
        &                             pos            = addedPos,       &
        &                             wasAdded       = wasAdded        )

      if (wasAdded) then
        write(logUnit(10),*) ' Appended variable:'//trim(varname)
      else if (addedpos < 1) then
        write(logUnit(1),*) 'Error: variable '//trim(varname)// &
          &                 ' is not added to variable system'
      end if

      deallocate(input_varname)
    end do

  end subroutine mus_append_derVar_fluid
  ! **************************************************************************** !


! ****************************************************************************** !
!       Subroutines with common interface for the function pointers            !
! ****************************************************************************** !

! ****************************************************************************** !
  !> Initiates the calculation of density
  !! This routine sets the function Pointer for density calcualtion and calls
  !! the generice get Element from PDF routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveDensity(fun, varsys, elempos, time, tree, nElems, &
    &                                nDofs, res                                )
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer :: fnCalcPtr
    ! -------------------------------------------------------------------- !

    fnCalcPtr => mus_derivedensity

    call mus_generic_fromPDF_forElement( &
      &  fun       = fun,                &
      &  varSys    = varSys,             &
      &  elempos   = elempos,            &
      &  tree      = tree,               &
      &  time      = time,               &
      &  nVals     = nElems,             &
      &  fnCalcPtr = fnCalcPtr,          &
      &  nDofs     = nDofs,              &
      &  res       = res                 )

  end subroutine deriveDensity
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the temperature of a given set of elements (sum up all links).
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveTemp(fun, varsys, elempos, time, tree, nElems, &
    &                             nDofs, res                                )
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

    res( 1:nElems ) = div1_3

  end subroutine deriveTemp
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the pressure of a given set of elements (sum up all links).
  !!
  !! Pressure calculation according to the isentropic equation of state for
  !! the LBM \( p = \rho c_s^2 \)
  !! with the calculation of density as in deriveDensity
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine derivePressure(fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res                )
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
    integer :: dens_pos
    ! -------------------------------------------------------------------- !
    ! position of density in glob system
    dens_pos = fun%input_varPos(1)
    ! get density variable from auxField
    call varSys%method%val(dens_pos)%get_element( varSys  = varSys,  &
      &                                           elemPos = elemPos, &
      &                                           time    = time,    &
      &                                           tree    = tree,    &
      &                                           nElems  = nElems,  &
      &                                           nDofs   = nDofs,   &
      &                                           res     = res      )

    ! convert density to pressure
    res(1:nElems) = res(1:nElems) * cs2

  end subroutine derivePressure
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the mach number of a given set of elements (sum up all links).
  !!
  !! Mach number calculation according to
  !! \( Ma = vel / c_s \)
  !! with the calculation of velicty as in deriveVelMag
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveMachNr( fun, varsys, elempos, time, tree, &
    &                                nElems, nDofs, res                )
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
    integer :: vel_mag_pos
    ! -------------------------------------------------------------------- !
    ! position of vel_mag in glob system
    vel_mag_pos = fun%input_varPos(1)
    ! get density variable from auxField
    call varSys%method%val(vel_mag_pos)%get_element( varSys  = varSys,  &
      &                                              elemPos = elemPos, &
      &                                              time    = time,    &
      &                                              tree    = tree,    &
      &                                              nElems  = nElems,  &
      &                                              nDofs   = nDofs,   &
      &                                              res     = res      )

    ! compute mach nr = vel/cs
    res(1:nElems) = res(1:nElems) * csInv

  end subroutine deriveMachNr
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the force on the boundary of a given set of elements
  !!
  !! Force calculation in LBM:
  !! \[ F = \sum_\alpha{\alpha \ne 0} \vec{e}_{\alphabar} (2f_\alphabar
  !!       + 6\omega_\alpha \vec{e}_\alpha \rho \vec{u_b})  \]
  !! \n where, \n
  !! \alpha = direction towards boundary
  !! \alphabar = direction opposite to boundary
  !! \vec{u_b} = velocity of the obstacle cell.
  !! For stationary object, this should be zero.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveBndForce(fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: iElem, posInBndID
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    ! compute bnd_force only if shape is defined as boundary kind
    do iElem = 1, nElems
      ! compute force only if current element has boundary property
      if ( btest( tree%elemPropertyBits( elemPos(iElem) ), &
        &         prp_hasBnd )  ) then

        ! position of current element in boundary_ID list
        posInBndID = fPtr%solverData%geometry%posInBndID( elempos(iElem) )
        if (posInBndID > 0) then
          res( (iElem-1)*fun%nComponents+1: iElem*fun%nComponents ) &
            & = fPtr%solverData%geometry%bndForce(posInBndID, :)
        end if ! posInBndID > 0

      end if ! element has_boundary
    end do ! iElem

  end subroutine deriveBndForce
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the moment on the boundary of a given set of elements
  !!
  !! Force calculation in LBM:
  !! \[ F = \sum_\alpha{\alpha \ne 0} \vec{e}_{\alphabar} (2f_\alphabar
  !!       + 6\omega_\alpha \vec{e}_\alpha \rho \vec{u_b})  \]
  !! \n where, \n
  !! \alpha = direction towards boundary
  !! \alphabar = direction opposite to boundary
  !! \vec{u_b} = velocity of the obstacle cell.
  !! Moment is computed from force using cross_product(arm, force).
  !! arm is computed from boundary element bary - reference point.
  !! For stationary object, this should be zero.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveBndMoment(fun, varsys, elempos, time, tree, &
    &                                  nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: iElem, posInBndID
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    ! compute bnd_force only if shape is defined as boundary kind
    do iElem = 1, nElems
      ! compute force only if current element has boundary property
      if ( btest( tree%elemPropertyBits( elemPos(iElem) ), &
        &         prp_hasBnd )  ) then

        ! position of current element in boundary_ID list
        posInBndID = fPtr%solverData%geometry%posInBndID( elempos(iElem) )
        if (posInBndID > 0) then
          res( (iElem-1)*fun%nComponents+1: iElem*fun%nComponents ) &
            & = fPtr%solverData%geometry%bndMoment(posInBndID, :)
        end if ! posInBndID > 0

      end if ! element has_boundary
    end do ! iElem

  end subroutine deriveBndMoment
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate momentum from density and velocity stored in auxField
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveMomentum(fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res    )
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
    !> Function pointer to perform specific operation.
    integer :: statePos, iElem, iLevel, elemOff
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos, vel_pos(3)
    real(kind=rk) :: rho, momentum(3), vel(3)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer, &
      &        auxField => fPtr%solverData%scheme%auxField,           &
      &        quantities => fPtr%solverData%scheme%layout%quantities )

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars
        ! position of density in auxField array
        dens_pos = varSys%method%val( fun%input_varPos(1) ) &
          &                     %auxField_varPos(1)

        ! position of velocity in auxField array
        vel_pos = varSys%method%val( fun%input_varPos(2) ) &
          &                     %auxField_varPos(1:3)

        ! density
        rho = auxField(iLevel)%val( elemOff + dens_pos )

        ! velocity
        vel = auxField(iLevel)%val(elemOff + vel_pos(1) : elemOff + vel_pos(3))
        ! compute and store momentum
        momentum = quantities%momentum_from_vel_dens_ptr(vel=vel, dens=rho)
        res((iElem-1)*3 + 1) = momentum(1)
        res((iElem-1)*3 + 2) = momentum(2)
        res((iElem-1)*3 + 3) = momentum(3)

      end do ! iElem
    end associate

  end subroutine deriveMomentum
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initiates the calculation of equlibrium
  !! This routine sets the function Pointer for equlibrium calcualtion and calls
  !! the generice get Element from PDF routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveEquil(fun, varsys, elempos, time, tree, nElems, &
    &                              nDofs, res                                )
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer :: fnCalcPtr
    ! -------------------------------------------------------------------- !
    fnCalcPtr => mus_deriveEquil

    call mus_generic_fromPDF_forElement( &
      &  fun       = fun,                &
      &  varSys    = varSys,             &
      &  elempos   = elempos,            &
      &  tree      = tree,               &
      &  time      = time,               &
      &  nVals     = nElems,             &
      &  fnCalcPtr = fnCalcPtr,          &
      &  nDofs     = nDofs,              &
      &  res       = res                 )

  end subroutine deriveEquil
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initiates the calculation of NonEquil
  !! This routine sets the function Pointer for NonEquil calcualtion and calls
  !! the generice get Element from PDF routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveNonEquil(fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res                )
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer :: fnCalcPtr
    ! -------------------------------------------------------------------- !
    fnCalcPtr => mus_deriveNonEquil

    call mus_generic_fromPDF_forElement( &
      &  fun       = fun,                &
      &  varSys    = varSys,             &
      &  elempos   = elempos,            &
      &  tree      = tree,               &
      &  time      = time,               &
      &  nVals     = nElems,             &
      &  fnCalcPtr = fnCalcPtr,          &
      &  nDofs     = nDofs,              &
      &  res       = res                 )

  end subroutine deriveNonEquil
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initiates the calculation of StrainRate
  !! This routine sets the function Pointer for StrainRate calcualtion and calls
  !! the generice get Element from PDF routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveStrainRate(fun, varsys, elempos, time, tree, &
    &                                   nElems, nDofs, res                )
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer :: fnCalcPtr
    ! -------------------------------------------------------------------- !
    fnCalcPtr => mus_deriveStrainRate

    call mus_generic_fromPDF_forElement( &
      &  fun       = fun,                &
      &  varSys    = varSys,             &
      &  elempos   = elempos,            &
      &  tree      = tree,               &
      &  time      = time,               &
      &  nVals     = nElems,             &
      &  fnCalcPtr = fnCalcPtr,          &
      &  nDofs     = nDofs,              &
      &  res       = res                 )
  end subroutine deriveStrainRate
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Jiaxing Qi
  !! Calculate the shear rate
  !!
  !! The Shear Rate is defined as
  !! \[
  !!  \dot{\gamma} = \sqrt{ 2D_{II} }
  !! \]
  !! where \( D_{II} \) is the second invariant of the strain rate tensor and
  !! defined as
  !! \[
  !!    D_{II} = \sum^{l}_{\alpha,\beta=l} S_{\alpha\beta} S_{\alpha\beta}
  !! \]
  !! where \( S_{\alpha\beta} \) is the strain rate tensor.
  !!
  recursive subroutine deriveShearRate(fun, varsys, elempos, time, tree, &
    &                                  nElems, nDofs, res                )
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
    integer :: posStrain, iElem
    real(kind=rk) :: strain(6*nElems)
    ! -------------------------------------------------------------------- !

    posStrain = fun%input_varPos(1)

    call varSys%method%val(posStrain)%get_element( varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            nElems  = nElems,  &
      &                                            nDofs   = nDofs,   &
      &                                            res     = strain   )

    do iElem = 1, nElems

      res( iElem ) = getShearRate( strain( (iElem-1)*6+1:iElem*6 ) )

    end do ! iElem

  end subroutine deriveShearRate
! ****************************************************************************** !

! ****************************************************************************** !
  !> author: Jiaxing Qi
  !! Calculate the deviatoric shear stress for Newtonian fluid
  !! (exclude pressure) (no mixtures).\n
  !! Shear Stress depends on variable: nonEquilibirium
  !!
  !! The formula is
  !! \[
  !!  \tau_{\alpha \beta}=
  !!    -(1-\frac{\omega}{2}) \sum_{i} f^{neq}_{i} c_{i\alpha} c_{i\beta}
  !! \]
  !! where \( \tau_{\alpha \beta}\) is the stress
  !! in the \(\beta\)-direction on a face normal to the \(\alpha\)-axis,\n
  !! \( f^{neq}_i = f_i - f^{eq}_i\) is the non-equilibirium density.\n
  !! For more information, please refer to:\n
  !! Krueger T, Varnik F, Raabe D. Shear stress in lattice Boltzmann
  !! simulations. Physical Review E. 2009;79(4):1-14.\n
  !! Thus this variable is dependent on another variable: nonEquilibirium.
  !!
  !! For multi-level mesh, Omega on finer level needs to be adjusted in order to
  !! get the correct shearstress calculation.\n
  !! First, we defines c as the dx ratio between finer and coarse level.\n
  !! \( c={ \Delta dx }_{ c }/{ \Delta dx }_{ f } \)
  !! Then the viscosity on the different levels must satisfy:\n
  !! \( \frac { { \nu  }_{ f } }{ { \nu  }_{ c } } =c \)
  !! This constrain leads to a relationship of omega on different levels:\n
  !! \( {\omega}_f = \frac {1}{ {\lambda}(\frac{1}{{\omega}_c}-0.5)+0.5 } \)
  !! For more information, please refer to:\n
  !! Manuel H, Harald K, Joerg B, Sabine R. Aeroacoustic validation of the
  !! lattice boltzmann method on non-uniform grids. ECCOMAS 2012
  !!
  recursive subroutine deriveShearStress(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: posNonEq, nCompsNonEq, iLevel, iElem
    real(kind=rk) :: omega
    real(kind=rk), allocatable :: nonEq(:)
    real(kind=rk), allocatable :: tau(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    posNonEq = fun%input_varPos(1)
    nCompsNonEq = varSys%method%val( posNonEq )%nComponents
    allocate( nonEq( nCompsNonEq * nElems ) )
    allocate(   tau( fun%nComponents ) )

    ! calculate nonEq for all elements
    call varSys%method%val(posNonEq)%get_element(  varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            nElems  = nElems,  &
      &                                            nDofs   = nDofs,   &
      &                                            res     = nonEq    )

    do iElem = 1, nElems

      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
      ! get the correct omega value
      omega = scheme%field(1)%fieldProp%fluid%viscKine%omLvl(iLevel)%val(iElem)

      ! compute second moment of nonEq
      ! here the 3D function works also for 2D, maybe in debug mode this
      ! will fail. cxcx is allocated according to the stencil:
      !   - 1 entry for 1D
      !   - 3 entries for 2D
      !   - 6 entries for 3D
      tau(:) = secondMom_3D( cxcx = scheme%layout%fStencil%cxcx(:,:), &
        &                    f    = nonEq( (iElem-1)*nCompsNonEq+1    &
        &                                 :iElem*nCompsNonEq ),       &
        &                    QQ   = scheme%layout%fStencil%QQ         )

      res( (iElem-1)*fun%nComponents+1: iElem*fun%nComponents ) = &
        &                   tau(:) * ( 0.5_rk * omega - 1._rk )

    end do ! iElem

    deallocate( nonEq )
    deallocate( tau )

  end subroutine deriveShearStress
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Jiaxing Qi
  !! Calculate the wall shear stress (WSS) of a given element with the given
  !! input
  !!
  !! The wall shear stress is shear stress exerted on the wall. Since it is well
  !! known that shear stress gets its maximum value on the wall, here we can
  !! directly calculate the principle stress of the stress tensor instead of
  !! multiplying the stress tensor by normal vector of the plane.
  !! The principle stresses are just the eigenvalues of the stress tensor. To
  !! get those eigenvalues, we need solve the characteristic equation:\n
  !! \( {\tau}^3 + a_{2}{\tau}^2 + a_1{\tau} + a_0 = 0 \)
  !! where \( a_2 = \tau_{ii} \)
  !! \( a_1 = (\tau_{ii}\tau_{jj} - \tau_{ij}\tau_{ji}) / 2 \)
  !! \( a_0 = det(\tau_{ij}) \). Here Einstein notation is used.
  !!
  recursive subroutine deriveWSS3D(fun, varsys, elempos, time, tree, nElems, &
    &                              nDofs, res                                )
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
    ! local variables
    integer :: iElem, shear_pos
    real(kind=rk) :: shear(6*nElems) ! shear stress of all
    ! -------------------------------------------------------------------- !

    ! wall shear stress is depend on shear stress. get its position in globSys
    shear_pos = fun%input_varPos(1)

    ! first compute shear stress, save it to array tau
    ! the shear tensor has the following order: x, y, z, xy, xz, yz
    call varSys%method%val(shear_pos)%get_element( varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            nElems  = nElems,  &
      &                                            nDofs   = nDofs,   &
      &                                            res     = shear    )

    do iElem = 1, nElems

      res( iElem ) = getWSS( shear( (iElem-1)*6+1 : iElem*6 ) )

    end do ! iElem

  end subroutine deriveWSS3D
! ****************************************************************************** !

! ****************************************************************************** !
  recursive subroutine deriveWSS2D(fun, varsys, elempos, time, tree, nElems, &
    &                              nDofs, res                                )
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
    ! local variables
    integer :: iElem, shear_pos
    real(kind=rk) :: a0, a1
    real(kind=rk) :: tau(6*nElems) ! shear stress of all
    ! -------------------------------------------------------------------- !

    ! wall shear stress is depend on shear stress. get its position in globSys
    shear_pos = fun%input_varPos(1)

    ! first compute shear stress, save it to array tau
    ! the shear tensor has the following order: x, y, z, xy, xz, yz
    call varSys%method%val(shear_pos)%get_element( varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            nElems  = nElems,  &
      &                                            nDofs   = nDofs,   &
      &                                            res     = tau      )


    do iElem = 1, nElems

      a1 = - tau((iElem-1)*6+1) - tau((iElem-1)*6+2)
      a0 =   tau((iElem-1)*6+1) * tau((iElem-1)*6+2) &
        &  - tau((iElem-1)*6+3) * tau((iElem-1)*6+3)

      res(iElem) = sqrt( a1 * a1 - 4._rk * a0 ) * div1_2

    end do ! iElem

  end subroutine deriveWSS2D
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the shear stress magnitude of a given element number with the
  !! given
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveShearMag(fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res                )
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
    integer :: shear_pos, iElem
    real(kind=rk) :: shear(6*nElems) ! shear stress of all
    real(kind=rk) :: tau(6)
    ! -------------------------------------------------------------------- !

    shear_pos    = fun%input_varPos(1)

    ! first compute shear stress, save it to array tau
    ! the shear tensor has the following order: x, y, z, xy, xz, yz
    call varSys%method%val(shear_pos)%get_element( varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            nElems  = nElems,  &
      &                                            nDofs   = nDofs,   &
      &                                            res     = shear    )


    do iElem = 1, nElems

      tau = shear( (iElem-1)*6+1 : iElem*6 )

      ! Von mises-criterion (see wikipedia article on von Mises yield criterion)
      ! wss = (tau_xx-tau_yy)**2 + (tau_yy-tau_zz)**2 + (tau_xx-tau_zz)**2  &
      ! &   + 6._rk*( tau_xy*tau_xy + tau_yz*tau_yz + tau_xz*tau_xz)
      ! wss = sqrt(0.5_rk*wss )
      res( iElem ) = sqrt( 0.5_rk*(   (tau(1)-tau(2))**2 &
        &                           + (tau(2)-tau(3))**2 &
        &                           + (tau(1)-tau(3))**2 &
        &                           + 6._rk * (   tau(4)*tau(4) &
        &                                       + tau(5)*tau(5) &
        &                                       + tau(6)*tau(6) )))

    end do ! iElem

  end subroutine deriveShearMag
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate kinetic energy from density and velocity in auxField
  !! This routine sets the function Pointer for kinetic energy calcualtion and
  !! calls the generice get Element from PDF routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveKE(fun, varsys, elempos, time, tree, nElems, &
    &                           nDofs, res                    )
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
    !> Function pointer to perform specific operation.
    integer :: statePos, iElem, iLevel, elemOff
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos, vel_pos(3)
    real(kind=rk) :: rho, vel(3)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer, &
      &        auxField => fPtr%solverData%scheme%auxField,           &
      &        quantities => fPtr%solverData%scheme%layout%quantities )

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars
        ! position of density in auxField array
        dens_pos = varSys%method%val( fun%input_varPos(1) ) &
          &                     %auxField_varPos(1)

        ! position of velocity in auxField array
        vel_pos = varSys%method%val( fun%input_varPos(2) ) &
          &                     %auxField_varPos(1:3)

        ! density
        rho = auxField(iLevel)%val( elemOff + dens_pos )

        ! velocity
        vel(1) = auxField(iLevel)%val(elemOff + vel_pos(1))
        vel(2) = auxField(iLevel)%val(elemOff + vel_pos(2))
        vel(3) = auxField(iLevel)%val(elemOff + vel_pos(3))

        ! compute and store kinetic energy
        res(iElem) = quantities%kineticEnergy_from_vel_dens_ptr(vel=vel, dens=rho)
      end do ! iElem
    end associate

  end subroutine deriveKE
! ****************************************************************************** !

! ****************************************************************************** !
   !> Derive absorb layer variable defined as a source term.
  recursive subroutine derive_absorbLayer(fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res                )
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

  end subroutine derive_absorbLayer
! ****************************************************************************** !


! ****************************************************************************** !
   !> Derive external force variable defined as a source term.
   !! It evaluates spacetime function defined in lua file for force variable
   !! and convert it to state value which is to be added to the state
  recursive subroutine derive_force_MRT(fun, varsys, elempos, time, tree,     &
    &                                   nElems, nDofs, res                    )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: forceField(nElems*3)
    real(kind=rk) :: G(3), velocity(3), ucx, uMinusCX(3)
    integer :: iElem, iDir, QQ, nScalars, posInTotal, elemOff
    integer :: vel_pos(3), iLevel
    real(kind=rk) :: omegaKine, omegaBulk, discForce
    real(kind=rk) :: forceTerm(27), momForce(27), s_mrt(27)
    real(kind=rk) :: mInvXOmega(27,27)
    ! -------------------------------------------------------------------- !
    !call tem_abort('Not implemented')

    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%input_varPos(1) )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme
    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%input_varPos(2))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = forceField                                 )

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ

    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    do iElem = 1, nElems
      ! get iLevel for element
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem ) ) )
      posInTotal = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )
      omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity(1) = scheme%auxField(iLevel)%val(elemOff + vel_pos(1))
      velocity(2) = scheme%auxField(iLevel)%val(elemOff + vel_pos(2))
      velocity(3) = scheme%auxField(iLevel)%val(elemOff + vel_pos(3))

      ! convert physical to lattice.
      ! force field on current element
      ! For incompressible model: this forceField should be divided by rho0.
      ! Since rho0 =1, this term is also valid for incompressible model
      G = forceField((iElem-1)*3+1 : iElem*3) &
        & / fPtr%solverData%physics%fac(iLevel)%body_force

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)
      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      s_mrt(1:QQ) = scheme%field(1)%fieldProp%fluid                              &
        &           %mrtPtr(omegaKine=omegaKine, omegaBulk=omegaBulk, QQ=QQ)

      ! M^-1 * (I-0.5 S)
      s_mrt(1:QQ) = 1.0_rk - 0.5_rk * s_mrt(1:QQ)
      do iDir = 1, QQ
        mInvXOmega(1:QQ,iDir) = scheme%layout%moment%toPDF%A(1:QQ,iDir) &
          &                   * s_mrt(iDir)
      end do

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )
        uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - velocity

        forceTerm(iDir) = scheme%layout%weight(iDir)             &
          &       * dot_product( uMinusCx * cs2inv               &
          &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
          &       * cs4inv, G )
      end do

      ! Force moments: M * F
      !do iDir = 1, QQ
      !  momForce(iDir) = sum(scheme%layout%moment%toMoments%A(iDir,:) * forceTerm)
      !end do
      momForce(1:QQ) = matmul( scheme%layout%moment%toMoments%A(1:QQ, 1:QQ), &
        &                forceTerm(1:QQ) )

      !discForce = matmul( omegaTerm, forceTerm )
      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = dot_product(mInvXOmega(iDir,1:QQ),  momForce(1:QQ))
        res((iElem-1)*fun%nComponents + iDir) = discForce
      end do

    end do !iElem

  end subroutine derive_force_MRT
! ****************************************************************************** !


! ************************************************************************** !
   !> Derive external force variable defined as a source term.
   !! It evaluates spacetime function defined in lua file for force variable
   !! and convert it to state value which is to be added to the state
   !! @todo KM: Not use we need seperate force for incompressible fluid model
  recursive subroutine derive_force1stOrd(fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res                )
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

  end subroutine derive_force1stOrd
! ************************************************************************** !


! ****************************************************************************** !
   !> Derive external force variable defined as a source term.
   !! It evaluates spacetime function defined in lua file for force variable
   !! and convert it to state value which is to be added to the state
  recursive subroutine derive_HRRCorrection_d2q9(fun, varsys, elempos, time, &
    &                                            tree, nElems, nDofs, res    )
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
    integer :: statePos, iElem, iLevel, QQ, iVar, HRRVar
    integer :: denspos, velpos(3)
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: gradRhoU3(2,1), S_Corr(9), dens, vel(2)
    ! -------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, fPtr )

    associate ( source => fPtr%solverData%scheme%field(1)%internalSource, &
      &         derVarPos => fPtr%solverData%scheme%derVarPos             )

      HRRVar = -1
      ! res is always AOS layout
      !write(*,*) 'Source nVals = ', source%varDict%nVals
      do iVar = 1, source%varDict%nVals
        !write(*,*) 'Source varDictKey = ', trim(source%varDict%val(iVar)%key)
        if (trim(source%varDict%val(iVar)%key) == 'hrr_correction') then
          HRRVar = iVar
        end if
      end do
      if (HRRVar == -1) then
        call tem_abort ('Source variable HRR_Correction not found')
      end if

      ! Get position of density and velocity in auxField to determine them later
      denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
      velpos(1:2) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:2)
    end associate

    res = 0.0_rk
    QQ = fun%nComponents
    do iElem = 1, nElems
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      statePos = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

      associate ( gradData => fPtr%solverData%scheme%gradData(iLevel),        &
        &         auxfield => fPtr%solverData%scheme%auxField(iLevel)%val(:), &
        &         weight   => fPtr%solverData%scheme%layout%weight,           &
        &         Grad     => fPtr%solverData%scheme%Grad                     )
        ! 1 = x, 2 = y, 3 = z, no xy returned
        gradRhoU3(:,1:1) = Grad%RhoU3_ptr(       &
          &   auxField     = auxField,           &
          &   gradData     = gradData,           &
          &   velPos       = velpos,             &
          &   densPos      = denspos,            &
          &   nAuxScalars  = varSys%nAuxScalars, &
          &   nDims        = 2,                  &
          &   nSolve       = 1,                  &
          &   elemOffset   = statePos-1          )

        ! Calculate correction
        call HRR_Correction_d2q9 (                    &
          &    QQ        = QQ,                        &
          &    weight    = weight,                    &
          &    gradRHOU3 = gradRHOU3(:, 1),           &
          &    phi       = S_corr(:),                 &
          &    dens      = dens,                      &
          &    vel       = vel(:)                     )

        res( (iElem-1)*QQ+1:(iElem-1)*QQ+QQ ) = S_Corr(:)

      end associate

    end do !iElem


  end subroutine derive_HRRCorrection_d2q9
! ****************************************************************************** !


! ****************************************************************************** !
   !> Derive external force variable defined as a source term.
   !! It evaluates spacetime function defined in lua file for force variable
   !! and convert it to state value which is to be added to the state
  recursive subroutine derive_HRRCorrection_d3q19(fun, varsys, elempos, time, &
    &                                             tree, nElems, nDofs, res    )
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
    integer :: statePos, iElem, iLevel, QQ, iVar, HRRVar
    integer :: denspos, velpos(3)
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: gradRhoU3(3,1), gradRhoUVZ(3,1)
    real(kind=rk) :: S_Corr(19), dens, vel(3)
    ! -------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, fPtr )

    associate ( source => fPtr%solverData%scheme%field(1)%internalSource, &
      &         derVarPos => fPtr%solverData%scheme%derVarPos             )

      HRRVar = -1
      ! res is always AOS layout
      !write(*,*) 'Source nVals = ', source%varDict%nVals
      do iVar = 1, source%varDict%nVals
        !write(*,*) 'Source varDictKey = ', trim(source%varDict%val(iVar)%key)
        if (trim(source%varDict%val(iVar)%key) == 'hrr_correction') then
          HRRVar = iVar
        end if
      end do
      if (HRRVar == -1) then
        call tem_abort ('Source variable HRR_Correction not found')
      end if

      ! Get position of density and velocity in auxField to determine them later
      denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
      velpos(1:3) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    end associate

    res = 0.0_rk
    QQ = fun%nComponents
    do iElem = 1, nElems
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      statePos = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

      associate ( gradData => fPtr%solverData%scheme%gradData(iLevel),        &
        &         auxfield => fPtr%solverData%scheme%auxField(iLevel)%val(:), &
        &         weight   => fPtr%solverData%scheme%layout%weight,           &
        &         Grad     => fPtr%solverData%scheme%Grad                     )
        ! 1 = x, 2 = y, 3 = z, no xy returned
        gradRhoU3(:,1:1) = Grad%RhoU3_ptr(       &
          &   auxField     = auxField,           &
          &   gradData     = gradData,           &
          &   velPos       = velpos,             &
          &   densPos      = denspos,            &
          &   nAuxScalars  = varSys%nAuxScalars, &
          &   nDims        = 3,                  &
          &   nSolve       = 1,                  &
          &   elemOffset   = statePos-1          )

        ! 1 = x, 2 = y, 3 = z, no xy returned
        gradRhoUVZ(:,1:1) = Grad%RhoUVZ_ptr(     &
          &   auxField     = auxField,           &
          &   gradData     = gradData,           &
          &   velPos       = velpos,             &
          &   densPos      = denspos,            &
          &   nAuxScalars  = varSys%nAuxScalars, &
          &   nDims        = 3,                  &
          &   nSolve       = 1,                  &
          &   elemOffset   = iElem-1             )

        ! Calculate correction
        call HRR_Correction_d3q19 (            &
          &    QQ         = QQ,                &
          &    weight     = weight,            &
          &    gradRHOU3  = gradRHOU3(:, 1),   &
          &    gradRHOUVZ = gradRHOUVZ(:, 1),  &
          &    phi        = S_corr(:),         &
          &    dens       = dens,              &
          &    vel        = vel(:)             )

        res( (iElem-1)*QQ+1:(iElem-1)*QQ+QQ ) = S_Corr(:)

      end associate

    end do !iElem

  end subroutine derive_HRRCorrection_d3q19
! ****************************************************************************** !

! ****************************************************************************** !
   !> Derive external force variable defined as a source term.
   !! It evaluates spacetime function defined in lua file for force variable
   !! and convert it to state value which is to be added to the state
  recursive subroutine derive_HRRCorrection_d3q27(fun, varsys, elempos, time, &
    &                                             tree, nElems, nDofs, res    )
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
    integer :: statePos, iElem, iLevel, QQ, iVar, HRRVar
    integer :: denspos, velpos(3)
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: gradRhoU3(3,1), S_Corr(27), dens, vel(3)
    ! -------------------------------------------------------------------- !

    call C_F_POINTER( fun%method_Data, fPtr )

    associate ( source => fPtr%solverData%scheme%field(1)%internalSource, &
      &         derVarPos => fPtr%solverData%scheme%derVarPos             )

      HRRVar = -1
      ! res is always AOS layout
      !write(*,*) 'Source nVals = ', source%varDict%nVals
      do iVar = 1, source%varDict%nVals
        !write(*,*) 'Source varDictKey = ', trim(source%varDict%val(iVar)%key)
        if (trim(source%varDict%val(iVar)%key) == 'hrr_correction') then
          HRRVar = iVar
        end if
      end do
      if (HRRVar == -1) then
        call tem_abort ('Source variable HRR_Correction not found')
      end if

      ! Get position of density and velocity in auxField to determine them later
      denspos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
      velpos(1:3) = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    end associate

    res = 0.0_rk
    QQ = fun%nComponents
    do iElem = 1, nElems
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      statePos = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

      associate ( gradData => fPtr%solverData%scheme%gradData(iLevel),        &
        &         auxfield => fPtr%solverData%scheme%auxField(iLevel)%val(:), &
        &         weight   => fPtr%solverData%scheme%layout%weight,           &
        &         Grad     => fPtr%solverData%scheme%Grad                     )
        ! 1 = x, 2 = y, 3 = z, no xy returned
        gradRhoU3(:,1:1) = Grad%RhoU3_ptr(       &
          &   auxField     = auxField,           &
          &   gradData     = gradData,           &
          &   velPos       = velpos,             &
          &   densPos      = denspos,            &
          &   nAuxScalars  = varSys%nAuxScalars, &
          &   nDims        = 3,                  &
          &   nSolve       = 1,                  &
          &   elemOffset   = statePos-1          )

        ! Calculate correction
        call HRR_Correction_d3q27 (                   &
          &    QQ        = QQ,                        &
          &    weight    = weight,                    &
          &    gradRHOU3 = gradRHOU3(:, 1),           &
          &    phi       = S_corr(:),                 &
          &    dens      = dens,                      &
          &    vel       = vel(:)                     )

        res( (iElem-1)*QQ+1:(iElem-1)*QQ+QQ ) = S_Corr(:)

      end associate

    end do !iElem

  end subroutine derive_HRRCorrection_d3q27
! ****************************************************************************** !

! ****************************************************************************** !
   !> Derive the Brinkman force variable defined as a source term.
   !!
   !! It multiplies the Brinkman coefficient variable as defined by the user
   !! with the negative velocity field to obtain the force term.
   !! This force term is then converted to a state value that is to be added to the
   !! state.
   !!
   !! Reference:
   !! 1) Zhaoli Guo and T. S. Zhao. “Lattice Boltzmann Model for Incompressible
   !!    Flows through Porous Media”. In: Physical Review E 66.3 (2002), p. 036304.
   !!    doi: 10.1103/PhysRevE.66.036304.
   !! 2) Irina Ginzburg. “Consistent lattice Boltzmann schemes for the Brinkman
   !!    model of porous flow and infinite Chapman-Enskog expansion”. In: Phys.
   !!    Rev. E 77 (6 June 2008), p. 066704. doi: 10.1103/PhysRevE.77.066704.
  recursive subroutine derive_brinkmanForce(fun, varsys, elempos, time, tree, &
    &                                        nElems, nDofs, res               )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: bCoeffField(nElems), bCoeffTerm
    real(kind=rk) :: velocity(3), ucx
    integer :: iElem, iDir, QQ, nScalars, posInTotal, elemOff
    integer :: vel_pos(3), iLevel
    real(kind=rk) :: omegaKine
    ! -------------------------------------------------------------------- !
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%input_varPos(1) )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Obtain the Brinkman coefficient from the respective variable
    ! defined in the configuration for all requested elements.
    call varSys%method%val(fun%input_varPos(2))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = bCoeffField                                )

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ

    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    do iElem = 1, nElems
      ! get iLevel for element
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
      posInTotal = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )

      ! element offset
      elemoff = (posInTotal - 1) * varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity = scheme%auxField(iLevel)%val(elemOff + vel_pos)

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)

      bCoeffTerm = bCoeffField(iElem) &
        &           / fPtr%solverData%physics%fac(iLevel)%sourceCoeff

      ! Brinkman force term:
      ! B_i = w_i \vec{e}_i/cs2 \cdot \vec{B}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )

        res( (iElem-1) * fun%nComponents + iDir ) = (omegaKine / 2.0_rk - 1.0_rk)     &
          &                                        * scheme%layout%weight(iDir)       &
          &                                        * cs2inv * bCoeffField(iElem) * ucx * rho0

      end do

    end do !iElem

  end subroutine derive_brinkmanForce
! ****************************************************************************** !


! ****************************************************************************** !
   !> Derive Brinkman force variable defined as a source term with TRT operator.
   !! It evaluates spacetime function defined in lua file for Brinkman coefficient
   !! variable and multiplies it with negative velocity to obtain the force term.
   !! This force term is then converted to state value which is to be added to the
   !! state with respect to TRT operator.
   !! Reference:
   !! 1) Zhaoli Guo and T. S. Zhao. “Lattice Boltzmann Model for Incompressible
   !!    Flows through Porous Media”. In: Physical Review E 66.3 (2002), p. 036304.
   !!    doi: 10.1103/PhysRevE.66.036304.
   !! 2) Irina Ginzburg. “Consistent lattice Boltzmann schemes for the Brinkman
   !!    model of porous flow and infinite Chapman-Enskog expansion”. In: Phys.
   !!    Rev. E 77 (6 June 2008), p. 066704. doi: 10.1103/PhysRevE.77.066704.
  recursive subroutine derive_brinkmanForce_TRT(fun, varsys, elempos, time, tree, &
    &                                            nElems, nDofs, res               )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: bCoeffField(nElems), bCoeffTerm
    real(kind=rk) :: velocity(3), ucx
    integer :: iElem, iDir, QQ, nScalars, posInTotal, elemOff
    integer :: vel_pos(3), iLevel
    real(kind=rk) :: omegaKine, omegaMinus
    ! -------------------------------------------------------------------- !
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%input_varPos(1) )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme
    ! Get Brinkman coefficient which is refered in config file either
    ! its spacetime variable or operation variable
    call varSys%method%val(fun%input_varPos(2))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = bCoeffField                                )

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ

    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    do iElem = 1, nElems
      ! get iLevel for element
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
      posInTotal = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity = scheme%auxField(iLevel)%val(elemOff + vel_pos)

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)
      omegaMinus = 1.0_rk / ( scheme%field(1)%fieldProp%fluid%lambda  &
        &                    / (1.0_rk / omegaKine - 0.5_rk) + 0.5_rk )

      bCoeffTerm = bCoeffField(iElem) &
        &           / fPtr%solverData%physics%fac(iLevel)%sourceCoeff

      ! force term:
      ! B_i = w_i \vec{e}_i/cs2 \cdot \vec{B}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )

        res( (iElem - 1) * fun%nComponents + iDir ) = (omegaMinus / 2.0_rk - 1.0_rk)    &
          &                                          * scheme%layout%weight(iDir)       &
          &                                          * cs2inv * bCoeffField(iElem) * ucx * rho0

      end do

    end do !iElem

  end subroutine derive_brinkmanForce_TRT
! ****************************************************************************** !



! **************************************************************************** !
  !> Update state with source variable "absorb_layer".
  !! absorb_layer is used to absorb the flow and gradually reduce the flow
  !! quantities like pressure and velocity to a fixed value.
  !! It is based on:
  !! Xu, H., & Sagaut, P. (2013). Analysis of the absorbing layers for the
  !! weakly-compressible lattice Boltzmann methods. Journal of Computational
  !! Physics, 245(x), 14–42.
  !! Jacob, J.; Sagaut, P. (2019): Solid wall and open boundary conditions in
  !! hybrid recursive regularized lattice Boltzmann method for compressible
  !! flows. In Physics of Fluids 31 (12), p. 126103.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_absorbLayer( fun, inState, outState, neigh, auxField,    &
    &                              nPdfSize, iLevel, varSys, time, phyConvFac, &
    &                              derVarPos                                   )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: spongeField(fun%elemLvl(iLevel)%nElems*5)
    real(kind=rk) :: dens, vel(3), ucx, uMinusCX(3), sponge_velTerm
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: sigma, dens_ref, vel_ref(3), sponge_vel(3), sponge_dens
    real(kind=rk) :: inv_rho_phy, inv_vel_phy
    real(kind=rk) :: omega, omega_fac
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    !call tem_abort('Error: Absorb layer is not yet implemented')
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = spongeField                               )

    inv_rho_phy = 1.0_rk / fPtr%solverData%physics%fac(iLevel)%press * cs2inv
    inv_vel_phy = 1.0_rk / fPtr%solverData%physics%fac(iLevel)%vel

    ! target pressure and velocity in lattice unit
    dens_ref = fun%absLayer%config%target_pressure * inv_rho_phy
    vel_ref(1:3) = fun%absLayer%config%target_velocity(1:3) * inv_vel_phy

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      dens = auxField(elemOff + dens_pos)
      ! obtain velocity from auxField
      vel(1) = auxField(elemOff + vel_pos(1))
      vel(2) = auxField(elemOff + vel_pos(2))
      vel(3) = auxField(elemOff + vel_pos(3))

      ! get the correct omega value
      omega = scheme%field(1)%fieldProp%fluid%viscKine              &
        &                              %omLvl(iLevel)%val(posInTotal)
      omega_fac = 1.0_rk - omega * 0.5_rk

      ! SpongeField contains: spongeStrength
      sigma = spongeField(iElem)

      sponge_vel(:) = -sigma*(vel - vel_ref)
      sponge_dens = -sigma*(dens - dens_ref)

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                vel )
        uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - vel

        sponge_velTerm = dot_product( uMinusCx * cs2inv          &
          &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
          &       * cs4inv, sponge_vel )

        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos)        &
          & + omega_fac * scheme%layout%weight( iDir ) &
          & * (sponge_dens + sponge_velTerm)

      end do

    end do !iElem

  end subroutine applySrc_absorbLayer
! **************************************************************************** !

! **************************************************************************** !
  !> Update state with source variable "absorb_layer".
  !! absorb_layer is used to absorb the flow and gradually reduce the flow
  !! quantities like pressure and velocity to a fixed value.
  !! It is based on:
  !! Xu, H., & Sagaut, P. (2013). Analysis of the absorbing layers for the
  !! weakly-compressible lattice Boltzmann methods. Journal of Computational
  !! Physics, 245(x), 14–42.
  !! Jacob, J.; Sagaut, P. (2019): Solid wall and open boundary conditions in
  !! hybrid recursive regularized lattice Boltzmann method for compressible
  !! flows. In Physics of Fluids 31 (12), p. 126103.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_absorbLayer_MRT( fun, inState, outState, neigh, auxField,&
    &                              nPdfSize, iLevel, varSys, time, phyConvFac, &
    &                              derVarPos                                   )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: spongeField(fun%elemLvl(iLevel)%nElems*5)
    real(kind=rk) :: dens, vel(3), ucx, uMinusCX(3)
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: sigma, dens_ref, vel_ref(3), sponge_vel(3), sponge_dens
    real(kind=rk) :: inv_rho_phy, inv_vel_phy
    real(kind=rk) :: omegaKine, omegaBulk, discForce, sponge_velTerm
    real(kind=rk) :: momForce(27), s_mrt(27)
    real(kind=rk) :: mInvXOmega(27,27), sponge_Term(27)
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    !call tem_abort('Error: Absorb layer is not yet implemented')
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = spongeField                               )

    inv_rho_phy = 1.0_rk / fPtr%solverData%physics%fac(iLevel)%press * cs2inv
    inv_vel_phy = 1.0_rk / fPtr%solverData%physics%fac(iLevel)%vel

    ! target pressure and velocity in lattice unit
    dens_ref = fun%absLayer%config%target_pressure * inv_rho_phy
    vel_ref(1:3) = fun%absLayer%config%target_velocity(1:3) * inv_vel_phy

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars

    ! Position of velocity variable in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)
    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      dens = auxField(elemOff + dens_pos)
      ! obtain velocity from auxField
      vel(1) = auxField(elemOff + vel_pos(1))
      vel(2) = auxField(elemOff + vel_pos(2))
      vel(3) = auxField(elemOff + vel_pos(3))

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)

      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      s_mrt(1:QQ) = scheme%field(1)%fieldProp%fluid                              &
        &           %mrtPtr(omegaKine=omegaKine, omegaBulk=omegaBulk, QQ=QQ)

      ! M^-1 * (I-0.5 S)
      s_mrt(1:QQ) = 1.0_rk - 0.5_rk * s_mrt(1:QQ)
      do iDir = 1, QQ
        mInvXOmega(1:QQ,iDir) = scheme%layout%moment%toPDF%A(1:QQ,iDir) &
          &                   * s_mrt(iDir)
      end do

      ! SpongeField contains: spongeStrength
      sigma = spongeField(iElem)

      ! Sponge factor for density and velocity field
      sponge_dens = -sigma*(dens - dens_ref)
      sponge_vel(:) = -sigma*(vel - vel_ref)

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), vel )
        uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - vel

        sponge_velTerm = dot_product( uMinusCx * cs2inv          &
          &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
          &       * cs4inv, sponge_vel )

        sponge_Term(iDir) = scheme%layout%weight(iDir)     &
          &               * ( sponge_dens  + sponge_velTerm )
      end do

      ! Force moments: M * F
      momForce(1:QQ) = matmul( scheme%layout%moment%toMoments%A(1:QQ, 1:QQ), &
        &                sponge_Term(1:QQ) )

      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = dot_product(mInvXOmega(iDir,1:QQ),  momForce(1:QQ))
        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos) + discForce
      end do

    end do !iElem

  end subroutine applySrc_absorbLayer_MRT
! **************************************************************************** !

! **************************************************************************** !
  !> Update state with source variable "absorb_layer".
  !! absorb_layer is used to absorb the flow and gradually reduce the flow
  !! quantities like pressure and velocity to a fixed value.
  !! It is based on:
  !! Xu, H., & Sagaut, P. (2013). Analysis of the absorbing layers for the
  !! weakly-compressible lattice Boltzmann methods. Journal of Computational
  !! Physics, 245(x), 14–42.
  !! Jacob, J.; Sagaut, P. (2019): Solid wall and open boundary conditions in
  !! hybrid recursive regularized lattice Boltzmann method for compressible
  !! flows. In Physics of Fluids 31 (12), p. 126103.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_absorbLayerDyn( fun, inState, outState, neigh, auxField, &
    &                              nPdfSize, iLevel, varSys, time, phyConvFac, &
    &                              derVarPos                                   )
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
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: spongeField(fun%elemLvl(iLevel)%nElems*5)
    real(kind=rk) :: dens, vel(3), ucx, uMinusCX(3), sponge_velTerm
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: sigma, sponge_vel(3), sponge_dens
    real(kind=rk) :: omega, omega_fac
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    !call tem_abort('Error: Absorb layer is not yet implemented')
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = spongeField                               )

    ! constant parameter
    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    associate( layout => fPtr%solverData%scheme%layout,                  &
      &        fluid => fPtr%solverData%scheme%field(1)%fieldProp%fluid, &
      &        dynAvg => fun%elemLvl(iLevel)%dynAvg                      )
      QQ = layout%fStencil%QQ
!$omp parallel do schedule(static), private( posInTotal, vel, ucx, uMinusCx, omega )
      do iElem = 1, nElems
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

        ! element offset
        elemoff = (posInTotal-1)*varSys%nAuxScalars
        dens = auxField(elemOff + dens_pos)
        ! obtain velocity from auxField
        vel(1) = auxField(elemOff + vel_pos(1))
        vel(2) = auxField(elemOff + vel_pos(2))
        vel(3) = auxField(elemOff + vel_pos(3))

        ! get the correct omega value
        omega = fluid%viscKine%omLvl(iLevel)%val(posInTotal)
        omega_fac = 1.0_rk - omega * 0.5_rk

        ! SpongeField contains: spongeStrength
        sigma = spongeField(iElem)

        sponge_dens = -sigma*(dens - dynAvg%dens(iElem))
        sponge_vel(1) = -sigma*(vel(1) - dynAvg%velX(iElem))
        sponge_vel(2) = -sigma*(vel(2) - dynAvg%velY(iElem))
        sponge_vel(3) = -sigma*(vel(3) - dynAvg%velZ(iElem))

        ! force term:
        ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
        !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
        do iDir = 1, QQ
          ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), vel )
          uMinusCx = layout%fStencil%cxDirRK(:, iDir) - vel

          sponge_velTerm = dot_product( uMinusCx * cs2inv   &
            &       + ucx * layout%fStencil%cxDirRK(:,iDir) &
            &       * cs4inv, sponge_vel )

          ! position in state array
          statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
          ! update outstate
          outState(statePos) = outState(statePos) &
            & + omega_fac * layout%weight( iDir ) &
            & * (sponge_dens + sponge_velTerm)

        end do

      end do !iElem
!$omp end parallel do
    end associate

  end subroutine applySrc_absorbLayerDyn
! **************************************************************************** !

! **************************************************************************** !
  !> Update state with source variable "absorb_layer".
  !! absorb_layer is used to absorb the flow and gradually reduce the flow
  !! quantities like pressure and velocity to a fixed value.
  !! It is based on:
  !! Xu, H., & Sagaut, P. (2013). Analysis of the absorbing layers for the
  !! weakly-compressible lattice Boltzmann methods. Journal of Computational
  !! Physics, 245(x), 14–42.
  !! Jacob, J.; Sagaut, P. (2019): Solid wall and open boundary conditions in
  !! hybrid recursive regularized lattice Boltzmann method for compressible
  !! flows. In Physics of Fluids 31 (12), p. 126103.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_absorbLayerDyn_MRT( fun, inState, outState, neigh,      &
    &                                     auxField, nPdfSize, iLevel, varSys, &
    &                                     time, phyConvFac, derVarPos         )
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
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: spongeField(fun%elemLvl(iLevel)%nElems*5)
    real(kind=rk) :: dens, vel(3), ucx, uMinusCX(3)
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: sigma, sponge_vel(3), sponge_dens
    real(kind=rk) :: omegaKine, omegaBulk, discForce, sponge_velTerm
    real(kind=rk) :: momForce(27), s_mrt(27)
    real(kind=rk) :: mInvXOmega(27,27), sponge_Term(27)
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    !call tem_abort('Error: Absorb layer is not yet implemented')
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = spongeField                               )

    ! constant parameter
    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    associate( layout => fPtr%solverData%scheme%layout,                  &
      &        fluid => fPtr%solverData%scheme%field(1)%fieldProp%fluid, &
      &        dynAvg => fun%elemLvl(iLevel)%dynAvg                      )
      QQ = layout%fStencil%QQ
      omegaBulk = fluid%omegaBulkLvl(iLevel)

      do iElem = 1, nElems
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

        ! element offset
        elemoff = (posInTotal-1)*varSys%nAuxScalars
        dens = auxField(elemOff + dens_pos)
        ! obtain velocity from auxField
        vel(1) = auxField(elemOff + vel_pos(1))
        vel(2) = auxField(elemOff + vel_pos(2))
        vel(3) = auxField(elemOff + vel_pos(3))

        ! get the correct omega value
        omegaKine = fluid%viscKine%omLvl(iLevel)%val(posInTotal)

        ! MRT omegas
        ! overwrite omegaKine term in the element loop
        s_mrt(1:QQ) = fluid%mrtPtr(omegaKine=omegaKine, omegaBulk=omegaBulk, QQ=QQ)

        ! M^-1 * (I-0.5 S)
        s_mrt(1:QQ) = 1.0_rk - 0.5_rk * s_mrt(1:QQ)
        do iDir = 1, QQ
          mInvXOmega(1:QQ,iDir) = layout%moment%toPDF%A(1:QQ,iDir) &
            &                   * s_mrt(iDir)
        end do

        ! SpongeField contains: spongeStrength
        sigma = spongeField(iElem)

        ! Sponge factor for density and velocity field
        sponge_dens = -sigma*(dens - dynAvg%dens(iElem))
        sponge_vel(1) = -sigma*(vel(1) - dynAvg%velX(iElem))
        sponge_vel(2) = -sigma*(vel(2) - dynAvg%velY(iElem))
        sponge_vel(3) = -sigma*(vel(3) - dynAvg%velZ(iElem))

        ! force term:
        ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
        !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
        do iDir = 1, QQ
          ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), vel )
          uMinusCx = layout%fStencil%cxDirRK(:, iDir) - vel

          sponge_velTerm = dot_product( uMinusCx * cs2inv   &
            &       + ucx * layout%fStencil%cxDirRK(:,iDir) &
            &       * cs4inv, sponge_vel )

          sponge_Term(iDir) = layout%weight(iDir)             &
            &               * ( sponge_dens  + sponge_velTerm )
        end do

        ! Force moments: M * F
        momForce(1:QQ) = matmul( layout%moment%toMoments%A(1:QQ, 1:QQ), &
          &                sponge_Term(1:QQ) )

        do iDir = 1, QQ
          ! discrete force
          ! \bar{F} =  M^-1 (I-S/2) M F
          discForce = dot_product(mInvXOmega(iDir,1:QQ),  momForce(1:QQ))
          ! position in state array
          statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
          ! update outstate
          outState(statePos) = outState(statePos) + discForce
        end do

      end do !iElem

    end associate

  end subroutine applySrc_absorbLayerDyn_MRT
! **************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "force".
  !! Force term used here is from:
  !! "Discrete lattice effects on the forcing term in the lattice Boltzmann
  !!  method", Zhaoli Guo, Chugung Zheng and Baochang Shi.
  !! In the paper, use force term is referred as Method 2 as:
  !! \[ F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !! Force must be defined as body force per unit volume
  !! KM: If this force formula is used then velocity needs to be
  !! computed as u = \sum c_i f_i + \vec{F}/2
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_force( fun, inState, outState, neigh, auxField,    &
    &                        nPdfSize, iLevel, varSys, time, phyConvFac, &
    &                        derVarPos                                   )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: G(3), velocity(3), ucx, uMinusCX(3), forceTerm
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omega, omega_fac
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / fPtr%solverData%physics%fac(iLevel)%body_force

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))
      velocity(3) = auxField(elemOff + vel_pos(3))

      ! force field on current element
      ! For incompressible model: this forceField should be divided by rho0.
      ! Since rho0 =1, this term is also valid for incompressible model
      G = forceField((iElem-1)*3+1 : iElem*3)

      ! get the correct omega value
      omega = scheme%field(1)%fieldProp%fluid%viscKine              &
        &                              %omLvl(iLevel)%val(posInTotal)
      omega_fac = 1.0_rk - omega * 0.5_rk

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )
        uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - velocity

        forceTerm = dot_product( uMinusCx * cs2inv               &
          &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
          &       * cs4inv, G )

        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos)                         &
          & + omega_fac * scheme%layout%weight( iDir ) * forceTerm

      end do

    end do !iElem

  end subroutine applySrc_force
! ****************************************************************************** !
  ! ****************************************************************************** !
  !> Update state with source variable "force" for Generalized Navier-Stokes equations.
  !! The implementation is taken from "Lattice Boltzmann model for incompressible flows
  !! through porous media" by Z. Guo and T.S. Zhao (2002) Physical Review E.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_force_GNS( fun, inState, outState, neigh, auxField,    &
    &                            nPdfSize, iLevel, varSys, time, phyConvFac, &
    &                            derVarPos                                   )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: G(3), velocity(3), ucx, uMinusCX(3), forceTerm
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), vol_frac_pos
    real(kind=rk) :: eps_f, eps_f_inv          ! Inverse of local fluid volume fraction for GNS
    real(kind=rk) :: omega, omega_fac
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / fPtr%solverData%physics%fac(iLevel)%body_force

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    vol_frac_pos = varSys%method%val(derVarPos(1)%vol_frac)%auxField_varPos(1)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))
      velocity(3) = auxField(elemOff + vel_pos(3))
      eps_f       = auxField(elemOff + vol_frac_pos)
      eps_f_inv   = 1.0_rk / eps_f

      ! force field on current element
      ! For incompressible model: this forceField should be divided by rho0.
      ! Since rho0 =1, this term is also valid for incompressible model
      G = forceField((iElem-1)*3+1 : iElem*3)

      ! get the correct omega value
      omega = scheme%field(1)%fieldProp%fluid%viscKine              &
        &                              %omLvl(iLevel)%val(posInTotal)
      omega_fac = 1.0_rk - omega * 0.5_rk

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )
        uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - velocity * eps_f_inv

        forceTerm = dot_product( uMinusCx * cs2inv               &
          &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
          &       * eps_f_inv * cs4inv, G*eps_f )

        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos)                         &
          & + omega_fac * scheme%layout%weight( iDir ) * forceTerm

      end do

    end do !iElem

  end subroutine applySrc_force_GNS
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "force" for MRT collision model.
  !! Force term used here is from:
  !! Chai, Z., & Zhao, T. (2012). Effect of the forcing term in the
  !! multiple-relaxation-time lattice Boltzmann equation on the shear stress
  !! or the strain rate tensor. Physical Review E, 86(1), 1–11.
  !! Force term for MRT is
  !! \[ \bar{F} = M^-1 (I-0.5 S) M F' \] and
  !! \[ F'_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !!
  !! \vec{F} is the force that must be defined as body force per unit volume
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_force_MRT( fun, inState, outState, neigh, auxField,    &
    &                            nPdfSize, iLevel, varSys, time, phyConvFac, &
    &                            derVarPos                                   )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: G(3), velocity(3), ucx, uMinusCX(3)
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omegaKine, omegaBulk, discForce
    real(kind=rk) :: forceTerm(27), momForce(27), s_mrt(27)
    real(kind=rk) :: mInvXOmega(27,27)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / fPtr%solverData%physics%fac(iLevel)%body_force

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ

    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))
      velocity(3) = auxField(elemOff + vel_pos(3))

      ! force field on current element
      ! For incompressible model: this forceField should be divided by rho0.
      ! Since rho0 =1, this term is also valid for incompressible model
      G = forceField((iElem-1)*3+1 : iElem*3)

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)
      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      s_mrt(1:QQ) = scheme%field(1)%fieldProp%fluid                              &
        &           %mrtPtr(omegaKine=omegaKine, omegaBulk=omegaBulk, QQ=QQ)

      ! M^-1 * (I-0.5 S)
      s_mrt(1:QQ) = 1.0_rk - 0.5_rk * s_mrt(1:QQ)
      do iDir = 1, QQ
        mInvXOmega(1:QQ,iDir) = scheme%layout%moment%toPDF%A(1:QQ,iDir) &
          &                   * s_mrt(iDir)
      end do

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )
        uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - velocity

        forceTerm(iDir) = scheme%layout%weight(iDir)             &
          &       * dot_product( uMinusCx * cs2inv               &
          &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
          &       * cs4inv, G )
      end do

      ! Force moments: M * F
      !do iDir = 1, QQ
      !  momForce(iDir) = sum(scheme%layout%moment%toMoments%A(iDir,:) * forceTerm)
      !end do
      momForce(1:QQ) = matmul( scheme%layout%moment%toMoments%A(1:QQ, 1:QQ), &
        &                forceTerm(1:QQ) )

      !discForce = matmul( omegaTerm, forceTerm )
      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = dot_product(mInvXOmega(iDir,1:QQ),  momForce(1:QQ))
        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos) + discForce
      end do

    end do !iElem

  end subroutine applySrc_force_MRT
! ****************************************************************************** !
! ****************************************************************************** !
  !> Update state with source variable "force" for MRT collision model.
  !! Force term used here is from:
  !! Chai, Z., & Zhao, T. (2012). Effect of the forcing term in the
  !! multiple-relaxation-time lattice Boltzmann equation on the shear stress
  !! or the strain rate tensor. Physical Review E, 86(1), 1–11.
  !! Force term for MRT is
  !! \[ \bar{F} = M^-1 (I-0.5 S) M F' \] and
  !! \[ F'_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !!
  !! \vec{F} is the force that must be defined as body force per unit volume
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_force_MRT_d3q27( fun, inState, outState, neigh,      &
    &                            auxField, nPdfSize, iLevel, varSys, time, &
    &                            phyConvFac,derVarPos                      )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: F(3), velocity(3)
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omegaBulk, discForce
    real(kind=rk) :: momForce(27), s_mrt(27)
    real(kind=rk) :: mInvXOmega(27,27)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / fPtr%solverData%physics%fac(iLevel)%body_force

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ

    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)
    ! MRT omegas
    s_mrt(1:QQ) = scheme%field(1)%fieldProp%fluid                              &
      &           %mrtPtr(omegaKine=1._rk, omegaBulk=omegaBulk, QQ=QQ)
    ! M^-1 * (I-0.5 S)
    s_mrt(2:4) = 1.0_rk - 0.5_rk * s_mrt(2:4)
    s_mrt(10) = 1.0_rk - 0.5_rk * s_mrt(10)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))
      velocity(3) = auxField(elemOff + vel_pos(3))

      ! force field on current element
      ! For incompressible model: this forceField should be divided by rho0.
      ! Since rho0 =1, this term is also valid for incompressible model
      F = forceField((iElem-1)*3+1 : iElem*3)

      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      ! get the correct omega value
      s_mrt(5:9) = scheme%field(1)%fieldProp%fluid%viscKine          &
      &                              %omLvl(iLevel)%val(posInTotal)

      ! M^-1 * (I-0.5 S)
      s_mrt(5:9) = 1.0_rk - 0.5_rk * s_mrt(5:9)
      do iDir = 2, 10
        mInvXOmega(1:QQ,iDir) = scheme%layout%moment%toPDF%A(1:QQ,iDir) &
          &                   * s_mrt(iDir)
      end do

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      ! Force moments: M * F
      momForce(1:QQ) = 0._rk
      momForce(2:4) = F(1:3)
      momForce(5) = F(1) * velocity(2) + F(2) * velocity(1)
      momForce(6) = F(2) * velocity(3) + F(3) * velocity(2)
      momForce(7) = F(1) * velocity(3) + F(3) * velocity(1)
      momForce(8) = -2._rk * ( F(2) * velocity(2) - 2._rk * F(1) * velocity(1) &
        &                      + F(3) * velocity(3) )
      momForce(9) = 2._rk * ( F(2) * velocity(2) - F(3) * velocity(3) )
      momForce(10) = 2._rk * ( F(1) * velocity(1) + F(2) * velocity(2) &
        &                      + F(3) * velocity(3) )

      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = dot_product(mInvXOmega(iDir,2:10),  momForce(2:10))
        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos) + discForce
      end do

    end do !iElem

  end subroutine applySrc_force_MRT_d3q27
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update state with source variable "force" for d3q19 MRT collision model.
  !! Force term used here is from:
  !! Chai, Z., & Zhao, T. (2012). Effect of the forcing term in the
  !! multiple-relaxation-time lattice Boltzmann equation on the shear stress
  !! or the strain rate tensor. Physical Review E, 86(1), 1–11.
  !! Force term for MRT is
  !! \[ \bar{F} = M^-1 (I-0.5 S) M F' \] and
  !! \[ F'_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !!
  !! Force must be defined as body force per unit volume
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_force_MRT_d3q19( fun, inState, outState, neigh,      &
    &                                  auxField, nPdfSize, iLevel, varSys, &
    &                                  time, phyConvFac, derVarPos         )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: F(3), velocity(3)
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omegaKine, omegaBulk, discForce
    real(kind=rk) :: momForce(19), s_mrt(19)
    real(kind=rk) :: mInvXOmega(19,19)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / fPtr%solverData%physics%fac(iLevel)%body_force

    ! constant parameter
    QQ = 19

    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)
    ! MRT omegas
    ! overwrite omegaKine term in the element loop
    ! KM: For incompressible model: omegaBulk is unused in mrtPtr
    s_mrt = scheme%field(1)%fieldProp%fluid                           &
      &           %mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    ! F = M^-1 (I-0.5 S) M F
    ! (I-0.5 S) - omega for force term
    s_mrt = 1.0_rk - 0.5_rk * s_mrt

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))
      velocity(3) = auxField(elemOff + vel_pos(3))

      ! force field on current element
      ! For incompressible model: this forceField should be divided by rho0.
      ! Since rho0 =1, this term is also valid for incompressible model
      F = forceField((iElem-1)*3+1 : iElem*3)

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)
      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      s_mrt(10) = 1.0_rk - 0.5_rk * omegaKine
      s_mrt(12) = s_mrt(10)
      s_mrt(14) = s_mrt(10)
      s_mrt(15) = s_mrt(10)
      s_mrt(16) = s_mrt(10)

      ! M^-1 (1-0.5 S)
      do iDir = 1, QQ
        mInvXOmega(:,iDir) = scheme%layout%moment%toPDF%A(:,iDir) * s_mrt(iDir)
      end do

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      ! Force moments: M * F
      momForce(1:QQ) = 0._rk
      momForce(2) = 2._rk * ( F(1) * velocity(1) + F(2) * velocity(2) &
        &                     + F(3) * velocity(3) )
      momForce(4) = F(1)
      momForce(6) = F(2)
      momForce(8) = F(3)
      momForce(10) = -2._rk * ( F(2) * velocity(2) - 2._rk * F(1) * velocity(1) &
        &                       + F(3) * velocity(3) )
      momForce(12) = 2._rk * ( F(2) * velocity(2) - F(3) * velocity(3) )
      momForce(14) = F(1) * velocity(2) + F(2) * velocity(1)
      momForce(15) = F(2) * velocity(3) + F(3) * velocity(2)
      momForce(16) = F(1) * velocity(3) + F(3) * velocity(1)

      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = sum(mInvXOmega(iDir,1:QQ) * momForce(1:QQ))
        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos) + discForce
      end do

    end do !iElem

  end subroutine applySrc_force_MRT_d3q19
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update state with source variable "force" for d3q19 MRT collision model.
  !! Force term used here is from:
  !! Chai, Z., & Zhao, T. (2012). Effect of the forcing term in the
  !! multiple-relaxation-time lattice Boltzmann equation on the shear stress
  !! or the strain rate tensor. Physical Review E, 86(1), 1–11.
  !! Force term for MRT is
  !! \[ \bar{F} = M^-1 (I-0.5 S) M F' \] and
  !! \[ F'_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
  !!       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F} \]
  !!
  !! Force must be defined as body force per unit volume
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_force_MRT_d2q9( fun, inState, outState, neigh,      &
    &                                  auxField, nPdfSize, iLevel, varSys, &
    &                                  time, phyConvFac, derVarPos         )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: F(3), velocity(3)
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omegaKine, omegaBulk, discForce
    real(kind=rk) :: momForce(9), s_mrt(9)
    real(kind=rk) :: mInvXOmega(9,9)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / fPtr%solverData%physics%fac(iLevel)%body_force

    ! constant parameter
    QQ = 9

    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)
    ! MRT omegas
    ! overwrite omegaKine term in the element loop
    ! KM: For incompressible model: omegaBulk is unused in mrtPtr
    s_mrt = scheme%field(1)%fieldProp%fluid                           &
      &           %mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    ! F = M^-1 (I-0.5 S) M F
    ! (I-0.5 S) - omega for force term
    s_mrt = 1.0_rk - 0.5_rk * s_mrt

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))

      ! force field on current element
      ! For incompressible model: this forceField should be divided by rho0.
      ! Since rho0 =1, this term is also valid for incompressible model
      F = forceField((iElem-1)*3+1 : iElem*3)

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)
      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      s_mrt(8:9) = 1.0_rk - 0.5_rk * omegaKine

      ! M^-1 (1-0.5 S)
      do iDir = 1, QQ
        mInvXOmega(:,iDir) = scheme%layout%moment%toPDF%A(:,iDir) * s_mrt(iDir)
      end do

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      ! Force moments: M * F
      momForce(1) = 0._rk
      momForce(2) = 6._rk * ( F(1) * velocity(1) + F(2) * velocity(2) )
      momForce(3) = -momForce(2)
      momForce(4) = F(1)
      momForce(5) = -F(1)
      momForce(6) = F(2)
      momForce(7) = -F(2)
      momForce(8) = 2._rk * ( F(1) * velocity(1) - F(2) * velocity(2) )
      momForce(9) = F(1) * velocity(2) + F(2) * velocity(1)

      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = sum(mInvXOmega(iDir,1:QQ) * momForce(1:QQ))
        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos) + discForce
      end do

    end do !iElem

  end subroutine applySrc_force_MRT_d2q9
! ****************************************************************************** !

! ************************************************************************** !
  !> Update state with source variable "force_1stOrd"
  !! Force term used here is from:
  !! "A D3Q27 multiple-relaxation-time lattice Boltzmann method for
  !! turbulent flows", K. Suga, Y. Kuwata, K. Takashima, R. Chikasue
  !!
  !! \[ F_i = w_i/c_s^2 ( \vec{e}_i \cdot \vec{F} ) \]
  !! Force must be defined as body force per unit volume
  !! This force term can be applied for both compressible and incompressible
  !! LBM models
  !!
  !! Similar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_force1stOrd( fun, inState, outState, neigh, auxField,    &
    &                              nPdfSize, iLevel, varSys, time, phyConvFac, &
    &                              derVarPos                                   )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: FF_elem(3)
    real(kind=rk) :: forceTerm
    integer :: nElems, iElem, iDir, QQ, nScalars
    integer :: posInTotal
    ! ---------------------------------------------------------------------- !
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

!write(dbgUnit(1),*) 'ApplySrc_force1stOrdIncomp'
!    do iElem = 1, nElems
!      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)
!      FF_elem = forceField((iElem-1)*3+1 : iElem*3)
!write(dbgUnit(1),*) 'treeID ', scheme%levelDesc(iLevel)%total(posInTotal) &
!  & , 'force ', FF_elem
!    end do

    ! convert physical to lattice
    forceField = forceField / fPtr%solverData%physics%fac(iLevel)%body_force

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars

    do iElem = 1, nElems
      ! to access level wise state array
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! force field on current element
      FF_elem = forceField((iElem-1)*3+1 : iElem*3)

      ! force term:
      ! F_i = w_i/cs2 (\vec{e}_i \cdot \vec{F}
      do iDir = 1, QQ

        forceTerm = dot_product( scheme%layout%fStencil    &
          &                      %cxDirRK(:,iDir), FF_elem )

        outState( (posintotal-1)*nscalars+idir+(1-1)*qq ) &
          & = outState(                                                  &
          & (posintotal-1)*nscalars+idir+(1-1)*qq )       &
          & + scheme%layout%weight( iDir ) * cs2inv * forceTerm

      end do

    end do !iElem
  end subroutine applySrc_force1stOrd
! ************************************************************************** !


! ****************************************************************************** !
  !> Update state with source variable "turb_channel_force" for BGK.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_turbChanForce( fun, inState, outState, neigh, auxField,  &
    &                                nPdfSize, iLevel, varSys, time,           &
    &                                phyConvFac, derVarPos                     )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: velocity(3), ucx, uMinusCX(3), forceTerm, dens
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: omega, omega_fac
    real(kind=rk) :: forceDynL(3)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Convert Dynamic force in m/s^2 to lattice
    forceDynL = fun%turbChanForce%forceDyn              &
      &       / fPtr%solverData%physics%fac(iLevel)%accel

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars
    ! Position of density and velocity variable in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain density from auxField
      dens = auxField(elemOff + dens_pos)
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))
      velocity(3) = auxField(elemOff + vel_pos(3))

      ! get the correct omega value
      omega = scheme%field(1)%fieldProp%fluid%viscKine              &
        &                              %omLvl(iLevel)%val(posInTotal)
      omega_fac = 1.0_rk - omega * 0.5_rk

      ! force term:
      ! F_i = w_i*(1-omega/2)*rho*( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{g}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )
        uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - velocity

        forceTerm = dot_product( uMinusCx * cs2inv               &
          &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
          &       * cs4inv, forceDynL )

        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos)                         &
          & + omega_fac * scheme%layout%weight( iDir ) * dens * forceTerm

      end do

    end do !iElem

  end subroutine applySrc_turbChanForce
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "force" for generic MRT collision model
  !! for turb_channel_force. It uses velocityX average in bulk to adapt the
  !! driving force for turbulent channel.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_turbChanForce_MRT( fun, inState, outState, neigh,  &
    &                                    auxField, nPdfSize, iLevel,     &
    &                                    varSys, time, phyConvFac,       &
    &                                    derVarPos                       )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: velocity(3), ucx, uMinusCX(3), dens
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: omegaKine, omegaBulk, discForce
    real(kind=rk) :: forceTerm(27), momForce(27), s_mrt(27)
    real(kind=rk) :: mInvXOmega(27,27)
    real(kind=rk) :: forceDynL(3)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! convert physical to lattice, force is defined in acceleration m/s^2
    ! Convert Dynamic force in m/s^2 to lattice
    forceDynL = fun%turbChanForce%forceDyn              &
      &       / fPtr%solverData%physics%fac(iLevel)%accel

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ

    nScalars = varSys%nScalars
    ! Position of density and velocity variable in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain density from auxField
      dens = auxField(elemOff + dens_pos)
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))
      velocity(3) = auxField(elemOff + vel_pos(3))

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)
      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      s_mrt(1:QQ) = scheme%field(1)%fieldProp%fluid                              &
        &           %mrtPtr(omegaKine=omegaKine, omegaBulk=omegaBulk, QQ=QQ)

      ! M^-1 * (I-0.5 S)
      s_mrt(1:QQ) = 1.0_rk - 0.5_rk * s_mrt(1:QQ)
      do iDir = 1, QQ
        mInvXOmega(1:QQ,iDir) = scheme%layout%moment%toPDF%A(1:QQ,iDir) &
          &                   * s_mrt(iDir)
      end do

      ! force term:
      ! F_i = w_i*(1-omega/2)*rho*( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )
        uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - velocity

        forceTerm(iDir) = scheme%layout%weight(iDir) * dens      &
          &       * dot_product( uMinusCx * cs2inv               &
          &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
          &       * cs4inv, forceDynL )
      end do

      ! Force moments: M * F
      !do iDir = 1, QQ
      !  momForce(iDir) = sum(scheme%layout%moment%toMoments%A(iDir,:) * forceTerm)
      !end do
      momForce(1:QQ) = matmul(scheme%layout%moment%toMoments%A(1:QQ,1:QQ), &
        &               forceTerm(1:QQ))

      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = sum(mInvXOmega(iDir,1:QQ) * momForce(1:QQ))
        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos) + discForce
      end do

    end do !iElem

  end subroutine applySrc_turbChanForce_MRT
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "force" for generic MRT collision model
  !! for turb_channel_force. It uses velocityX average in bulk to adapt the
  !! driving force for turbulent channel.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_turbChanForce_MRT_d3q27( fun, inState, outState, neigh, &
    &                                    auxField, nPdfSize, iLevel,          &
    &                                    varSys, time, phyConvFac,            &
    &                                    derVarPos                            )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: F(3), velocity(3), dens
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: omegaBulk, discForce
    real(kind=rk) :: momForce(27), s_mrt(27)
    real(kind=rk) :: mInvXOmega(27,27)
    real(kind=rk) :: forceDynL(3)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! convert physical to lattice, force is defined in acceleration m/s^2
    ! Convert Dynamic force in m/s^2 to lattice
    forceDynL = fun%turbChanForce%forceDyn              &
      &       / fPtr%solverData%physics%fac(iLevel)%accel

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ

    nScalars = varSys%nScalars
    ! Position of density and velocity variable in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)
    ! MRT omegas
    s_mrt(1:QQ) = scheme%field(1)%fieldProp%fluid                              &
      &           %mrtPtr(omegaKine=1._rk, omegaBulk=omegaBulk, QQ=QQ)
    ! M^-1 * (I-0.5 S)
    s_mrt(2:4) = 1.0_rk - 0.5_rk * s_mrt(2:4)
    s_mrt(10) = 1.0_rk - 0.5_rk * s_mrt(10)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain density from auxField
      dens = auxField(elemOff + dens_pos)
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))
      velocity(3) = auxField(elemOff + vel_pos(3))

      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      ! get the correct omega value
      s_mrt(5:9) = scheme%field(1)%fieldProp%fluid%viscKine          &
      &                              %omLvl(iLevel)%val(posInTotal)

      ! M^-1 * (I-0.5 S)
      s_mrt(5:9) = 1.0_rk - 0.5_rk * s_mrt(5:9)
      do iDir = 2, 10
        mInvXOmega(1:QQ,iDir) = scheme%layout%moment%toPDF%A(1:QQ,iDir) &
          &                   * s_mrt(iDir)
      end do

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      ! Force moments: M * F
      F = dens * forceDynL
      momForce(1:QQ) = 0._rk
      momForce(2:4) = F(1:3)
      momForce(5) = F(1) * velocity(2) + F(2) * velocity(1)
      momForce(6) = F(2) * velocity(3) + F(3) * velocity(2)
      momForce(7) = F(1) * velocity(3) + F(3) * velocity(1)
      momForce(8) = -2._rk * ( F(2) * velocity(2) - 2._rk * F(1) * velocity(1) &
        &                      + F(3) * velocity(3) )
      momForce(9) = 2._rk * ( F(2) * velocity(2) - F(3) * velocity(3) )
      momForce(10) = 2._rk * ( F(1) * velocity(1) + F(2) * velocity(2) &
        &                      + F(3) * velocity(3) )

      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = sum(mInvXOmega(iDir,2:10) * momForce(2:10))
        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos) + discForce
      end do

    end do !iElem

  end subroutine applySrc_turbChanForce_MRT_d3q27
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "force" for d3q19 MRT collision model
  !! for turb_channel_force. It uses velocityX average in bulk to adapt the
  !! driving force for turbulent channel.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_turbChanForce_MRT_d3q19( fun, inState, outState, neigh,  &
    &                                          auxField, nPdfSize, iLevel,     &
    &                                          varSys, time, phyConvFac,       &
    &                                          derVarPos                       )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: F(3), velocity(3), dens
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: omegaKine, omegaBulk, discForce
    real(kind=rk) :: momForce(19), s_mrt(19)
    real(kind=rk) :: mInvXOmega(19,19)
    real(kind=rk) :: forceDynL(3)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! convert physical to lattice, force is defined in acceleration m/s^2
    ! Convert Dynamic force in m/s^2 to lattice
    forceDynL = fun%turbChanForce%forceDyn              &
      &       / fPtr%solverData%physics%fac(iLevel)%accel

    ! constant parameter
    QQ = 19

    nScalars = varSys%nScalars
    ! Position of density and velocity variable in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)
    ! MRT omegas
    ! overwrite omegaKine term in the element loop
    ! KM: For incompressible model: omegaBulk is unused in mrtPtr
    s_mrt = scheme%field(1)%fieldProp%fluid                           &
      &           %mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    ! F = M^-1 (I-0.5 S) M F
    ! (I-0.5 S) - omega for force term
    s_mrt = 1.0_rk - 0.5_rk * s_mrt

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain density from auxField
      dens = auxField(elemOff + dens_pos)
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))
      velocity(3) = auxField(elemOff + vel_pos(3))

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)
      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      s_mrt(10) = 1.0_rk - 0.5_rk * omegaKine
      s_mrt(12) = s_mrt(10)
      s_mrt(14) = s_mrt(10)
      s_mrt(15) = s_mrt(10)
      s_mrt(16) = s_mrt(10)

      ! M^-1 (1-0.5 S)
      do iDir = 1, QQ
        mInvXOmega(:,iDir) = scheme%layout%moment%toPDF%A(:,iDir) * s_mrt(iDir)
      end do

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      ! Force moments: M * F
      F = dens * forceDynL
      momForce(1:QQ) = 0._rk
      momForce(2) = 2._rk * ( F(1) * velocity(1) + F(2) * velocity(2) &
        &                     + F(3) * velocity(3) )
      momForce(4) = F(1)
      momForce(6) = F(2)
      momForce(8) = F(3)
      momForce(10) = -2._rk * ( F(2) * velocity(2) - 2._rk * F(1) * velocity(1) &
        &                       + F(3) * velocity(3) )
      momForce(12) = 2._rk * ( F(2) * velocity(2) - F(3) * velocity(3) )
      momForce(14) = F(1) * velocity(2) + F(2) * velocity(1)
      momForce(15) = F(2) * velocity(3) + F(3) * velocity(2)
      momForce(16) = F(1) * velocity(3) + F(3) * velocity(1)

      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = sum(mInvXOmega(iDir,1:QQ) * momForce(1:QQ))
        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos) + discForce
      end do

    end do !iElem

  end subroutine applySrc_turbChanForce_MRT_d3q19
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "force" for d3q19 MRT collision model
  !! for turb_channel_force. It uses velocityX average in bulk to adapt the
  !! driving force for turbulent channel.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_turbChanForce_MRT_d2q9( fun, inState, outState, neigh,  &
    &                                          auxField, nPdfSize, iLevel,     &
    &                                          varSys, time, phyConvFac,       &
    &                                          derVarPos                       )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: F(3), velocity(3), dens
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: omegaKine, omegaBulk, discForce
    real(kind=rk) :: momForce(9), s_mrt(9)
    real(kind=rk) :: mInvXOmega(9,9)
    real(kind=rk) :: forceDynL(3)
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! convert physical to lattice, force is defined in acceleration m/s^2
    ! Convert Dynamic force in m/s^2 to lattice
    forceDynL = fun%turbChanForce%forceDyn              &
      &       / fPtr%solverData%physics%fac(iLevel)%accel

    ! constant parameter
    QQ = 9

    nScalars = varSys%nScalars
    ! Position of density and velocity variable in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = scheme%field(1)%fieldProp%fluid%omegaBulkLvl(iLevel)
    ! MRT omegas
    ! overwrite omegaKine term in the element loop
    ! KM: For incompressible model: omegaBulk is unused in mrtPtr
    s_mrt = scheme%field(1)%fieldProp%fluid                           &
      &           %mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    ! F = M^-1 (I-0.5 S) M F
    ! (I-0.5 S) - omega for force term
    s_mrt = 1.0_rk - 0.5_rk * s_mrt

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! obtain density from auxField
      dens = auxField(elemOff + dens_pos)
      ! obtain velocity from auxField
      velocity(1) = auxField(elemOff + vel_pos(1))
      velocity(2) = auxField(elemOff + vel_pos(2))

      ! get the correct omega value
      omegaKine = scheme%field(1)%fieldProp%fluid%viscKine          &
        &                              %omLvl(iLevel)%val(posInTotal)
      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      s_mrt(8:9) = 1.0_rk - 0.5_rk * omegaKine

      ! M^-1 (1-0.5 S)
      do iDir = 1, QQ
        mInvXOmega(:,iDir) = scheme%layout%moment%toPDF%A(:,iDir) * s_mrt(iDir)
      end do

      ! force term:
      ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
      !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}
      ! Force moments: M * F
      F = dens * forceDynL
      momForce(1) = 0._rk
      momForce(2) = 6._rk * ( F(1) * velocity(1) + F(2) * velocity(2) )
      momForce(3) = -momForce(2)
      momForce(4) = F(1)
      momForce(5) = -F(1)
      momForce(6) = F(2)
      momForce(7) = -F(2)
      momForce(8) = 2._rk * ( F(1) * velocity(1) - F(2) * velocity(2) )
      momForce(9) = F(1) * velocity(2) + F(2) * velocity(1)

      do iDir = 1, QQ
        ! discrete force
        ! \bar{F} =  M^-1 (I-S/2) M F
        discForce = sum(mInvXOmega(iDir,1:QQ) * momForce(1:QQ))
        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos) + discForce
      end do

    end do !iElem

  end subroutine applySrc_turbChanForce_MRT_d2q9
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable Brinkman force obtained from Brinkman
  !! coefficient.
  !!
  !! Reference:
  !! 1) Zhaoli Guo and T. S. Zhao. “Lattice Boltzmann Model for Incompressible
  !!    Flows through Porous Media”. In: Physical Review E 66.3 (2002), p. 036304.
  !!    doi: 10.1103/PhysRevE.66.036304.
  !! 2) Irina Ginzburg. “Consistent lattice Boltzmann schemes for the Brinkman
  !!    model of porous flow and infinite Chapman-Enskog expansion”. In: Phys.
  !!    Rev. E 77 (6 June 2008), p. 066704. doi: 10.1103/PhysRevE.77.066704.
  !!
  !! Similar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_brinkmanForce( fun, inState, outState, neigh, auxField, &
    &                                 nPdfSize, iLevel, varSys, time,         &
    &                                 phyConvFac, derVarPos                   )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: bCoeffField(fun%elemLvl(iLevel)%nElems)
    real(kind=rk) :: velocity(3), ucx
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omega, omega_fac
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = bCoeffField                               )

    bCoeffField = bCoeffField / fPtr%solverData%physics%fac(iLevel)%sourceCoeff

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal - 1) * varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity = auxField(elemOff + vel_pos)

      ! get the correct omega value
      omega = scheme%field(1)%fieldProp%fluid%viscKine              &
        &                              %omLvl(iLevel)%val(posInTotal)
      omega_fac = 1.0_rk - omega * 0.5_rk

      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )

        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos)                                &
          &                   - omega_fac * scheme%layout%weight( iDir ) * ucx &
          &                   * cs2inv * bCoeffField(iElem) * rho0

      end do

    end do !iElem

  end subroutine applySrc_brinkmanForce
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable Brinkman force obtained from Brinkman
  !! coefficient. The force update is for TRT collision model.
  !!
  !! Reference:
  !! 1) Zhaoli Guo and T. S. Zhao. “Lattice Boltzmann Model for Incompressible
  !!    Flows through Porous Media”. In: Physical Review E 66.3 (2002), p. 036304.
  !!    doi: 10.1103/PhysRevE.66.036304.
  !! 2) Irina Ginzburg. “Consistent lattice Boltzmann schemes for the Brinkman
  !!    model of porous flow and infinite Chapman-Enskog expansion”. In: Phys.
  !!    Rev. E 77 (6 June 2008), p. 066704. doi: 10.1103/PhysRevE.77.066704.
  !!
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_var_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_brinkmanForce_TRT( fun, inState, outState, neigh,       &
    &                                     auxField, nPdfSize, iLevel, varSys, &
    &                                     time, phyConvFac, derVarPos         )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: bCoeffField(fun%elemLvl(iLevel)%nElems)
    real(kind=rk) :: velocity(3), ucx
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos, elemOff
    integer :: vel_pos(3)
    real(kind=rk) :: omega, omega_fac, omegaMinus
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = bCoeffField                               )

    bCoeffField = bCoeffField / fPtr%solverData%physics%fac(iLevel)%sourceCoeff

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars
    ! Position of velocity variable in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! element offset
      elemoff = (posInTotal - 1) * varSys%nAuxScalars
      ! obtain velocity from auxField
      velocity = auxField(elemOff + vel_pos)

      ! get the correct omega value
      omega = scheme%field(1)%fieldProp%fluid%viscKine              &
        &                              %omLvl(iLevel)%val(posInTotal)
      omegaMinus = 1.0_rk / ( scheme%field(1)%fieldProp%fluid%lambda &
        &                    / (1.0_rk / omega - 0.5_rk) + 0.5_rk )
      omega_fac = 1.0_rk - omegaMinus * 0.5_rk

      do iDir = 1, QQ
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
          &                velocity )

        ! position in state array
        statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq
        ! update outstate
        outState(statePos) = outState(statePos)                          &
          &                   - omega_fac * scheme%layout%weight(iDir)   &
          &                   * ucx * cs2inv * bCoeffField(iElem) * rho0

      end do

    end do !iElem

  end subroutine applySrc_brinkmanForce_TRT
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate wss from shear stress (tau)
  !! tau: x, y, z, xy, yz, xz
  pure function getWSS( tau ) result( wss )
    ! ---------------------------------------------------------------------------
    real(kind=rk), intent(in) :: tau(6)
    real(kind=rk) :: a0, a1, a2
    real(kind=rk) :: q_stress, q_stress_cube, q_term, r_stress, d_stress
    real(kind=rk) :: theta_stress, cos_term, shear_1, shear_2, shear_3
    real(kind=rk) :: wss
    ! ---------------------------------------------------------------------------

      ! tau_x + tau_y + tau_z
      a2 = -1._rk * ( tau(1) + tau(2) + tau(3) )

      ! tau_x * tau_y + tau_y * tau_z + tau_x * tau_z - tau_xy^2 - tau_xz^2 - tau_yz^2
      a1 =    tau(1) * tau(2) + tau(2) * tau(3) + tau(1) * tau(3)              &
        &   - tau(4) * tau(4) - tau(5) * tau(5) - tau(6) * tau(6)

      a0 = -1._rk * (   tau(1) * tau(2) * tau(3)                               &
        &             + 2._rk * tau(4) * tau(5) * tau(6)                       &
        &             - tau(1) * tau(6) * tau(6)                               &
        &             - tau(2) * tau(5) * tau(5)                               &
        &             - tau(3) * tau(4) * tau(4) )

      q_stress = div1_9 * ( 3._rk * a1 - a2 * a2 )
      q_stress_cube = q_stress * q_stress * q_stress
      if ( q_stress < 0._rk ) then
        q_term = 2._rk * sqrt( -q_stress )
      else
        q_term = 0._rk
      end if
      r_stress = div1_54 * ( 9._rk*a2*a1 - 27._rk*a0 - 2._rk*a2*a2*a2 )
      d_stress = q_stress_cube + r_stress * r_stress

      if ( d_stress < 0._rk ) then
        ! the solutions are real and unequal
        theta_stress = div1_3 * acos( r_stress / sqrt( -q_stress_cube ))

        cos_term = div3_4 * q_term * cos( theta_stress )

        ! three values of maximum shear stress
        shear_1 = abs( sqrt3 * div1_2 * q_term * sin( theta_stress ))

        shear_2 = abs( cos_term - shear_1*div1_2 )

        shear_3 = abs( cos_term + shear_1*div1_2 )

        ! choose the max shear as Wall Shear Stress(wss)
        wss = max(shear_1, shear_2, shear_3)
      else
        wss = 0._rk
      end if
  end function getWSS
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine computes equilbrium from density and velocity
  !! This must comply with interface in mus_variable_module derive_FromMacro
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromMacro]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:equilFromMacro]] function pointer.
  subroutine deriveEquil_FromMacro( density, velocity, iField, nElems, varSys, &
    &                               layout, res                                )
    ! -------------------------------------------------------------------- !
    !> Array of density.
    !! Single species: dens_1, dens_2 .. dens_n
    !! multi-species: dens_1_sp1, dens_1_sp2, dens_2_sp1, dens_2_sp2 ...
    !!                dens_n_sp1, dens_n_sp2
    real(kind=rk), intent(in) :: density(:)

    !> Array of velocity.
    !! Size: dimension 1: n*nFields. dimension 2: 3 (nComp)
    !! 1st dimension arrangement for multi-species is same as density
    real(kind=rk), intent(in) :: velocity(:, :)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n*nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: fEq(layout%fStencil%QQ), vel(3)
    integer :: QQ, iElem
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    do iElem = 1, nElems
      vel = velocity(:, iElem)

      fEq = layout%quantities%pdfEq_ptr( rho = density(iElem),    &
        &                                vel = vel,               &
        &                                QQ = QQ                  )

      res( (iElem-1)*QQ+1: iElem*QQ ) = fEq
    end do
  end subroutine deriveEquil_FromMacro
! ****************************************************************************** !

! ************************************************************************** !
  !> This routine computes equilbrium from auxField
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_equilFromAux]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:equilFromAux]] function pointer.
  subroutine deriveEquil_fromAux( derVarPos, auxField, iField, nElems, varSys, &
    &                             layout, fEq                                  )
    ! -------------------------------------------------------------------- !
    !> Position of derive variable in variable system
    class(mus_derVarPos_type), intent(in) :: derVarPos
    !> Array of auxField.
    !! Single species: dens_1, vel_1, dens_2, vel_2, .. dens_n, vel_n
    !! multi-species: dens_1_sp1, vel_1_spc1, dens_1_sp2, vel_1_spc2,
    !!                dens_2_sp1, vel_2_spc2, dens_2_sp2, vel_2_spc2 ...
    !!                dens_n_sp1, vel_n_sp1, dens_n_sp2, vel_n_spc2
    !! Access: (iElem-1)*nAuxScalars + auxField_varPos
    real(kind=rk), intent(in) :: auxField(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n*QQ of res
    real(kind=rk), intent(out) :: fEq(:)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rho, vel(3)
    integer :: QQ, iElem, elemOff
    integer :: dens_pos, vel_pos(3)
    ! ---------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)

    QQ = layout%fStencil%QQ
    !NEC$ ivdep
    do iElem = 1, nElems
      ! element offset
      elemoff = (iElem-1)*varSys%nAuxScalars
      ! density
      rho = auxField(elemOff + dens_pos)
      ! velocity
      vel(1) = auxField(elemOff + vel_pos(1))
      vel(2) = auxField(elemOff + vel_pos(2))
      vel(3) = auxField(elemOff + vel_pos(3))

      ! calculate equilibrium density
      fEq( (iElem-1)*QQ + 1 : iElem*QQ ) = layout%quantities%pdfEq_ptr( rho = rho, &
        &                                                       vel = vel, &
        &                                                       QQ = QQ    )

    end do
  end subroutine deriveEquil_fromAux
! ************************************************************************** !

! ****************************************************************************** !
  !> This routine computes density from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! mus_derVarPos_type:equilFromState, mus_derVarPos_type:momFromState,
  !! mus_derVarPos_type:velocitiesFromState, and
  !! mus_derVarPos_type:momentaFromState function pointers.
  subroutine deriveRho_FromState( state, iField, nElems, varSys, layout, res )
    ! -------------------------------------------------------------------- !
    !> Array of state
    !! n * layout%fStencil%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iDir
    ! ---------------------------------------------------------------------------

    res = 0.0_rk
    do iElem = 1, nElems
      do iDir = 1, layout%fStencil%QQ
        res( iElem ) = res( iElem ) + state( iDir+(iElem-1)*varSys%nScalars )
      end do
    end do

  end subroutine deriveRho_FromState
! ****************************************************************************** !

! ****************************************************************************** !
  !> This routine computes velocity from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]], [[mus_derVarPos_type:momFromState],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState] function pointers.
  subroutine deriveVel_FromState( state, iField, nElems, varSys, layout, res )
    ! -------------------------------------------------------------------- !
    !> Array of state
    !! n * layout%fStencil%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iDir
    real(kind=rk) :: rho, pdf( layout%fStencil%QQ ), vel(3)
    ! ---------------------------------------------------------------------------

    do iElem = 1, nElems
      do iDir = 1, layout%fStencil%QQ
        pdf(iDir) = state( iDir+(iElem-1)*varSys%nScalars )
      end do
      rho = sum( pdf )
      vel = layout%quantities%vel_from_pdf_ptr(pdf = pdf, dens = rho)
      res( (iElem-1)*3+1 ) = vel(1)
      res( (iElem-1)*3+2 ) = vel(2)
      res( (iElem-1)*3+3 ) = vel(3)
    end do

  end subroutine deriveVel_FromState
! ****************************************************************************** !

! ****************************************************************************** !
  !> This routine computes velocity from pre collision state array using Fetch
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromPreColState]] in derived/[[mus_derVarPos_module]].f90 in order
  !! to be callable via [[mus_derVarPos_type:velFromPreColState]] function
  !! pointer.
  subroutine deriveVel_FromPreColState( state, neigh, iField, nSize, nElems, &
                                        varSys, layout, res                  )
    ! -------------------------------------------------------------------- !
    !> Array of state
    !! n * layout%fStencil%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> connectivity array
    integer, intent(in) :: neigh(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements in state array
    integer, intent(in) :: nSize

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iDir
    real(kind=rk) :: rho, pdf( layout%fStencil%QQ ), vel(3)
    integer :: QQ, nScalars, nDims
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nScalars = varSys%nScalars
    nDims = layout%fStencil%nDims

    do iElem = 1, nElems
      do iDir = 1, QQ
        pdf(iDir) = state(                                          &
          &  neigh((idir-1)* nsize+ ielem)+( ifield-1)* qq+ nscalars*0)
      end do
      ! density
      rho = sum( pdf )
      ! velocity
      vel = layout%quantities%vel_from_pdf_ptr(pdf = pdf, dens = rho)

      ! return velocity field according on stencil dimensions
      res( (iElem-1)*nDims+1:iElem*nDims) = vel(1:nDims)
    end do

  end subroutine deriveVel_FromPreColState
! ****************************************************************************** !

! **************************************************************************** !
  !> This routine computes auxField from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_auxFromState]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:auxFieldFromState]] function pointer.
  subroutine deriveAux_fromState( derVarPos, state, neigh, iField, nElems, &
    &                             nSize, iLevel, stencil, varSys, auxField, quantities )
    ! -------------------------------------------------------------------- !
    !> Position of derive variable in variable system
    class(mus_derVarPos_type), intent(in) :: derVarPos
    !> Array of state
    !! n * layout%stencil(1)%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> connectivity vector
    integer, intent(in) :: neigh(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> number of elements in state array
    integer, intent(in) :: nSize

    !> current level
    integer, intent(in) :: iLevel

    !> stencil header contains discrete velocity vectors
    type(tem_stencilHeader_type), intent(in) :: stencil

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities

    !> Output of this routine
    !! Size: nElems*nAuxScalars
    real(kind=rk), intent(inout) :: auxField(:)
    ! -------------------------------------------------------------------- !
    integer :: dens_pos, vel_pos(3)
    integer :: iElem, iDir, elemOff
    real(kind=rk) :: pdf( stencil%QQ ), rho, vel(3)
    ! ------------------------------------------------------------------------ !
    dens_pos = varSys%method%val(derVarPos%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)
    !NEC$ ivdep
    do iElem = 1, nElems
      !NEC$ shortloop
      do iDir = 1, stencil%QQ
        pdf(iDir) = state(                                                     &
          & ( ielem-1)* varsys%nscalars+idir+( 1-1)* stencil%qq)
      end do
      ! element offset
      elemoff = (iElem-1)*varSys%nAuxScalars

      ! density
      rho = sum( pdf )
      auxField(elemOff+dens_pos) = rho

      ! velocity
      vel = quantities%vel_from_pdf_ptr(pdf = pdf, dens = rho, cxDirRK = stencil%cxDirRK)
      auxField(elemOff+vel_pos(1)) = vel(1)
      auxField(elemOff+vel_pos(2)) = vel(2)
      auxField(elemOff+vel_pos(3)) = vel(3)
    end do

  end subroutine deriveAux_fromState
! **************************************************************************** !


! ****************************************************************************** !
  !> This routine computes velocity from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveEq_FromState( state, iField, nElems, varSys, layout, res )
    ! -------------------------------------------------------------------- !
    !> Array of state
    !! n * layout%fStencil%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: QQ, iElem, iDir
    real(kind=rk) :: pdf( layout%fStencil%QQ ), fEq( layout%fStencil%QQ )
    real(kind=rk) :: rho, vel(3)
    ! ---------------------------------------------------------------------------

    QQ = layout%fStencil%QQ

    do iElem = 1, nElems

      do iDir = 1, QQ
        pdf(iDir) = state( iDir+(iElem-1)*varSys%nScalars )
      end do

      rho = sum( pdf )
      vel = layout%quantities%vel_from_pdf_ptr(pdf = pdf, dens = rho)

      ! calculate equilibrium density
      fEq = layout%quantities%pdfEq_ptr( rho = rho, &
        &                                vel = vel, &
        &                                QQ = QQ    )
      res( (iElem-1)*QQ+1 : iElem*QQ ) = fEq

    end do

  end subroutine deriveEq_FromState
! ****************************************************************************** !



! ****************************************************************************** !
  !> Initiates the calculation of moment for 2D
  !! This routine sets the function Pointer for moment for 2D calcualtion and
  !! calls the generice get Element from PDF routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  !> For 2D only!
  recursive subroutine deriveMoment(fun, varsys, elempos, time, tree, nElems, &
    &                               nDofs, res                                )
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer :: fnCalcPtr
    ! -------------------------------------------------------------------- !
    fnCalcPtr => mus_deriveMoment

    call mus_generic_fromPDF_forElement( &
      &  fun       = fun,                &
      &  varSys    = varSys,             &
      &  elempos   = elempos,            &
      &  tree      = tree,               &
      &  time      = time,               &
      &  nVals     = nElems,             &
      &  fnCalcPtr = fnCalcPtr,          &
      &  nDofs     = nDofs,              &
      &  res       = res                 )

  end subroutine deriveMoment
! ****************************************************************************** !

  ! ************************************************************************* !
  !         Subroutines with common interface for values from point           !
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> Calculates pressure for given set of points.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_point]].
  !!
  recursive subroutine derivePressure_forPoint(fun, varsys, point, time, &
    &                                          tree, nPnts, res          )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: dens_pos
    ! ---------------------------------------------------------------------------
    ! position of density in glob system
    dens_pos = fun%input_varPos(1)
    ! get density variable from auxField
    call varSys%method%val( dens_pos )%get_point( &
      &    varSys = varSys,                       &
      &    point  = point,                        &
      &    time   = time,                         &
      &    tree   = tree,                         &
      &    nPnts  = nPnts,                        &
      &    res    = res                           )

    ! convert density to pressure
    res(1:nPnts) = res(1:nPnts) * cs2

  end subroutine derivePressure_forPoint
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> Calculates Mach nr for given set of points.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_point]].
  !!
  recursive subroutine deriveMachNr_forPoint(fun, varsys, point, time, &
    &                                        tree, nPnts, res          )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: vel_mag_pos
    ! ---------------------------------------------------------------------------
    ! position of vel_mag in glob system
    vel_mag_pos = fun%input_varPos(1)
    ! get vel_mag variable
    call varSys%method%val( vel_mag_pos )%get_point( &
      &    varSys = varSys,                          &
      &    point  = point,                           &
      &    time   = time,                            &
      &    tree   = tree,                            &
      &    nPnts  = nPnts,                           &
      &    res    = res                              )

    ! compute Ma = vel/cs
    res(1:nPnts) = res(1:nPnts) * csInv

  end subroutine deriveMachNr_forPoint
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> Calculates Kinetic energy from density and velocity for given set of points.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_point]].
  !!
  recursive subroutine deriveKE_forPoint(fun, varsys, point, time,   &
    &                                    tree, nPnts, res )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: dens_pos, vel_pos, iPnt
    real(kind=rk) :: mass_dens(nPnts), vel_res(3*nPnts), vel(3)
    !> Function pointer to perform specific operation.
    type(mus_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !
    ! ---------------------------------------------------------------------------
    ! get mass density values for IDX
    dens_pos = fun%input_varPos(1)
    ! get density variable from auxField
    call varSys%method%val( dens_pos )%get_point( &
      &    varSys = varSys,                       &
      &    point  = point,                        &
      &    time   = time,                         &
      &    tree   = tree,                         &
      &    nPnts  = nPnts,                        &
      &    res    = mass_dens                     )

    ! position of velocity in glob system
    vel_pos = fun%input_varPos(2)
    ! get velocity variable from auxField
    call varSys%method%val( vel_pos )%get_point( &
      &    varSys = varSys,                      &
      &    point  = point,                       &
      &    time   = time,                        &
      &    tree   = tree,                        &
      &    nPnts  = nPnts,                       &
      &    res    = vel_res                      )


    call C_F_POINTER( fun%method_Data, fPtr )
    associate( quantities => fPtr%solverData%scheme%layout%quantities )
      ! compute kinetic energy
      do iPnt = 1, nPnts
        vel(1) = vel_res((iPnt-1)*3 + 1)
        vel(2) = vel_res((iPnt-1)*3 + 2)
        vel(3) = vel_res((iPnt-1)*3 + 3)
        res(iPnt) = quantities%kineticEnergy_from_vel_dens_ptr(vel=vel, dens=mass_dens(iPnt))
      end do
    end associate
  end subroutine deriveKE_forPoint
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> Calculates momentum from density and velocity for given set of points.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_point]].
  !!
  recursive subroutine deriveMomentum_forPoint(fun, varsys, point, time,    &
    &                                          tree, nPnts, res )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: dens_pos, vel_pos, iPnt
    real(kind=rk) :: mass_dens(nPnts), vel(3)
    !> Function pointer to perform specific operation.
    type(mus_varSys_data_type), pointer :: fPtr
    ! ---------------------------------------------------------------------------
    ! get mass density values for IDX
    dens_pos = fun%input_varPos(1)
    ! get density variable from auxField
    call varSys%method%val( dens_pos )%get_point( &
      &    varSys = varSys,                      &
      &    point  = point,                       &
      &    time   = time,                        &
      &    tree   = tree,                        &
      &    nPnts  = nPnts,                       &
      &    res    = mass_dens                    )

    ! position of velocity in glob system
    vel_pos = fun%input_varPos(2)
    ! get velocity variable from auxField
    call varSys%method%val( vel_pos )%get_point( &
      &    varSys = varSys,                      &
      &    point  = point,                       &
      &    time   = time,                        &
      &    tree   = tree,                        &
      &    nPnts  = nPnts,                       &
      &    res    = res                          )

    call C_F_POINTER( fun%method_Data, fPtr )
    associate( quantities => fPtr%solverData%scheme%layout%quantities )
      ! convert velocity to momentum
      do iPnt = 1, nPnts
        vel = res( (iPnt-1)*3 + 1 : iPnt*3 )
        res( (iPnt-1)*3 + 1 : iPnt*3 ) = quantities%momentum_from_vel_dens_ptr(vel=vel,dens=mass_dens(iPnt))
      end do
    end associate

  end subroutine deriveMomentum_forPoint
  ! ************************************************************************* !

  ! ************************************************************************* !
  !         Subroutines with common interface for values from index           !
  ! ************************************************************************* !

! ****************************************************************************** !
  !> Initiates the calculation of density.
  !! This routine sets the function Pointer for density calcualtion and calls
  !! the generice get Value of Index routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveDensity_fromIndex( fun, varSys, time, iLevel,  &
     &                                          idx, idxLen, nVals, res )
      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in)     :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)       :: time

      !> Level on which values are requested
      integer, intent(in)                   :: iLevel

      !> Index of points in the growing array and variable val array to
      !! return.
      !! Size: most times nVals, if contiguous arrays are used it depends
      !! on the number of first indices
      integer, intent(in)                   :: idx(:)

      !> With idx as start index in contiguous memory,
      !! idxLength defines length of each contiguous memory
      !! Size: dependes on number of first index for contiguous array,
      !! but the sum of all idxLen is equal to nVals
      integer, optional, intent(in)         :: idxLen(:)

      !> Number of values to obtain for this variable (vectorized access).
      integer, intent(in)                   :: nVals

      !> Resulting values for the requested variable.
      !!
      !! Dimension: n requested entries x nComponents of this variable
      !! Access: (iElem-1)*fun%nComponents + iComp
      real(kind=rk), intent(out)            :: res(:)
    ! ---------------------------------------------------------------------------
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer  :: fnCalcPtr
    ! ---------------------------------------------------------------------------

    fnCalcPtr => mus_derivedensity

    call mus_generic_varFromPDF_fromIndex( &
      &  fun       = fun,                  &
      &  varSys    = varSys,               &
      &  time      = time,                 &
      &  iLevel    = iLevel,               &
      &  idx       = idx,                  &
      &  nVals     = nVals,                &
      &  fnCalcPtr = fnCalcPtr,            &
      &  res       = res                   )

  end subroutine deriveDensity_fromIndex
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initiates the calculation of pressure.
  !! This routine sets the function Pointer for pressure calcualtion and calls
  !! the generice get Value of Index routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine derivePressure_fromIndex( fun, varSys, time, iLevel, &
     &                                           idx, idxLen, nVals, res )
      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in)     :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)       :: time

      !> Level on which values are requested
      integer, intent(in)                   :: iLevel

      !> Index of points in the growing array and variable val array to
      !! return.
      !! Size: most times nVals, if contiguous arrays are used it depends
      !! on the number of first indices
      integer, intent(in)                   :: idx(:)

      !> With idx as start index in contiguous memory,
      !! idxLength defines length of each contiguous memory
      !! Size: dependes on number of first index for contiguous array,
      !! but the sum of all idxLen is equal to nVals
      integer, optional, intent(in)         :: idxLen(:)

      !> Number of values to obtain for this variable (vectorized access).
      integer, intent(in)                   :: nVals

      !> Resulting values for the requested variable.
      !!
      !! Dimension: n requested entries x nComponents of this variable
      !! Access: (iElem-1)*fun%nComponents + iComp
      real(kind=rk), intent(out)            :: res(:)
    ! ---------------------------------------------------------------------------
    integer :: dens_pos
    type(mus_varSys_data_type), pointer     :: fPtr
    ! ---------------------------------------------------------------------------
    !convert pointer from C to Fotran
    call C_F_POINTER( fun%method_Data, fPtr )
    ! position of density in glob system
    dens_pos = fun%input_varPos(1)
    ! get density variable from auxField
    call varSys%method%val( dens_pos )%get_ValOfIndex( &
      &     varSys = varSys,                           &
      &     time   = time,                             &
      &     iLevel = iLevel,                           &
      &     idx    = fPtr%opData%input_pntIndex(1)     &
      &              %indexLvl(iLevel)%val( idx(:) ),  &
      &     nVals  = nVals,                            &
      &     res    = res                               )

    ! convert density to pressure
    res(1:nVals) = res(1:nVals) * cs2

  end subroutine derivePressure_fromIndex
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initiates the calculation of mach number.
  !! This routine sets the function Pointer for pressure calcualtion and calls
  !! the generice get Value of Index routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveMachNr_fromIndex( fun, varSys, time, iLevel, &
     &                                         idx, idxLen, nVals, res )
      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in)     :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)       :: time

      !> Level on which values are requested
      integer, intent(in)                   :: iLevel

      !> Index of points in the growing array and variable val array to
      !! return.
      !! Size: most times nVals, if contiguous arrays are used it depends
      !! on the number of first indices
      integer, intent(in)                   :: idx(:)

      !> With idx as start index in contiguous memory,
      !! idxLength defines length of each contiguous memory
      !! Size: dependes on number of first index for contiguous array,
      !! but the sum of all idxLen is equal to nVals
      integer, optional, intent(in)         :: idxLen(:)

      !> Number of values to obtain for this variable (vectorized access).
      integer, intent(in)                   :: nVals

      !> Resulting values for the requested variable.
      !!
      !! Dimension: n requested entries x nComponents of this variable
      !! Access: (iElem-1)*fun%nComponents + iComp
      real(kind=rk), intent(out)            :: res(:)
    ! ---------------------------------------------------------------------------
    integer :: vel_mag_pos
    type(mus_varSys_data_type), pointer     :: fPtr
    ! ---------------------------------------------------------------------------
    !convert pointer from C to Fotran
    call C_F_POINTER( fun%method_Data, fPtr )
    ! position of density in glob system
    vel_mag_pos = fun%input_varPos(1)
    ! get density variable from auxField
    call varSys%method%val( vel_mag_pos )%get_ValOfIndex( &
      &        varSys = varSys,                           &
      &        time   = time,                             &
      &        iLevel = iLevel,                           &
      &        idx    = fPtr%opData%input_pntIndex(1)     &
      &                 %indexLvl(iLevel)%val( idx(:) ),  &
      &        nVals  = nVals,                            &
      &        res    = res                               )

    ! convert density to pressure
    res(1:nVals) = res(1:nVals) * csInv

  end subroutine deriveMachNr_fromIndex
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initiates the calculation of equilibrium.
  !! This routine sets the function Pointer for equilibrium calcualtion and calls
  !! the generice get Value of Index routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveEquil_fromIndex( fun, varSys, time, iLevel, idx, &
     &                                        idxLen, nVals, res )
      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in)     :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)       :: time

      !> Level on which values are requested
      integer, intent(in)                   :: iLevel

      !> Index of points in the growing array and variable val array to
      !! return.
      !! Size: most times nVals, if contiguous arrays are used it depends
      !! on the number of first indices
      integer, intent(in)                   :: idx(:)

      !> With idx as start index in contiguous memory,
      !! idxLength defines length of each contiguous memory
      !! Size: dependes on number of first index for contiguous array,
      !! but the sum of all idxLen is equal to nVals
      integer, optional, intent(in)         :: idxLen(:)

      !> Number of values to obtain for this variable (vectorized access).
      integer, intent(in)                   :: nVals

      !> Resulting values for the requested variable.
      !!
      !! Dimension: n requested entries x nComponents of this variable
      !! Access: (iElem-1)*fun%nComponents + iComp
      real(kind=rk), intent(out)            :: res(:)
    ! ---------------------------------------------------------------------------
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer  :: fnCalcPtr
    ! ---------------------------------------------------------------------------

    fnCalcPtr => mus_deriveEquil

    call mus_generic_varFromPDF_fromIndex( &
      &  fun       = fun,                  &
      &  varSys    = varSys,               &
      &  time      = time,                 &
      &  iLevel    = iLevel,               &
      &  idx       = idx,                  &
      &  nVals     = nVals,                &
      &  fnCalcPtr = fnCalcPtr,            &
      &  res       = res                   )

  end subroutine deriveEquil_fromIndex
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initiates the calculation of non_equilibrium.
  !! This routine sets the function Pointer for non_equilibrium calcualtion and
  !! calls the generice get Value of Index routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveNonEquil_fromIndex( fun, varSys, time, iLevel, idx, &
     &                                        idxLen, nVals, res )
      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in)     :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)       :: time

      !> Level on which values are requested
      integer, intent(in)                   :: iLevel

      !> Index of points in the growing array and variable val array to
      !! return.
      !! Size: most times nVals, if contiguous arrays are used it depends
      !! on the number of first indices
      integer, intent(in)                   :: idx(:)

      !> With idx as start index in contiguous memory,
      !! idxLength defines length of each contiguous memory
      !! Size: dependes on number of first index for contiguous array,
      !! but the sum of all idxLen is equal to nVals
      integer, optional, intent(in)         :: idxLen(:)

      !> Number of values to obtain for this variable (vectorized access).
      integer, intent(in)                   :: nVals

      !> Resulting values for the requested variable.
      !!
      !! Dimension: n requested entries x nComponents of this variable
      !! Access: (iElem-1)*fun%nComponents + iComp
      real(kind=rk), intent(out)            :: res(:)
    ! ---------------------------------------------------------------------------
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer  :: fnCalcPtr
    ! ---------------------------------------------------------------------------

    fnCalcPtr => mus_deriveNonEquil

    call mus_generic_varFromPDF_fromIndex( &
      &  fun       = fun,                  &
      &  varSys    = varSys,               &
      &  time      = time,                 &
      &  iLevel    = iLevel,               &
      &  idx       = idx,                  &
      &  nVals     = nVals,                &
      &  fnCalcPtr = fnCalcPtr,            &
      &  res       = res                   )

  end subroutine deriveNonEquil_fromIndex
! ****************************************************************************** !

! ****************************************************************************** !
  !> Initiates the calculation of kinetic_energy.
  !! This routine sets the function Pointer for kinetic_energy calcualtion and
  !! calls the generice get Value of Index routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveKe_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                      idxLen, nVals, res  )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in)     :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)       :: time

    !> Level on which values are requested
    integer, intent(in)                   :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: most times nVals, if contiguous arrays are used it depends
    !! on the number of first indices
    integer, intent(in)                   :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: dependes on number of first index for contiguous array,
    !! but the sum of all idxLen is equal to nVals
    integer, optional, intent(in)         :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in)                   :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out)            :: res(:)
    ! ---------------------------------------------------------------------------
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos, vel_pos, iVal
    real(kind=rk) :: mass_dens(nVals), vel_res(3*nVals), vel(3)
    ! ---------------------------------------------------------------------------


    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk

    ! get mass density values for IDX
    dens_pos = fun%input_varPos(1)
    call varSys%method%val( dens_pos )%get_ValOfIndex( &
      &     varSys = varSys,                           &
      &     time   = time,                             &
      &     iLevel = iLevel,                           &
      &     idx    = fPtr%opData%input_pntIndex(1)     &
      &              %indexLvl(iLevel)%val( idx(:) ),  &
      &     nVals  = nVals,                            &
      &     res    = mass_dens                         )

    ! get velocity values for IDX
    vel_pos = fun%input_varPos(2)
    call varSys%method%val( vel_pos )%get_ValOfIndex( &
      &     varSys = varSys,                          &
      &     time   = time,                            &
      &     iLevel = iLevel,                          &
      &     idx    = fPtr%opData%input_pntIndex(2)    &
      &              %indexLvl(iLevel)%val( idx(:) ), &
      &     nVals  = nVals,                           &
      &     res    = vel_res(:)                       )

    associate( quantities => fPtr%solverData%scheme%layout%quantities )
      ! compute kinetic energy
      do iVal = 1, nVals
        vel(1) = vel_res((iVal-1)*3 + 1)
        vel(2) = vel_res((iVal-1)*3 + 2)
        vel(3) = vel_res((iVal-1)*3 + 3)
        res(iVal) = quantities%kineticEnergy_from_vel_dens_ptr(vel=vel, dens=mass_dens(iVal))
      end do
    end associate

  end subroutine deriveKe_fromIndex
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initiates the calculation of kinetic_energy.
  !! This routine sets the function Pointer for kinetic_energy calcualtion and
  !! calls the generice get Value of Index routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveStrainRate_fromIndex( fun, varSys, time, iLevel,  &
     &                                                 idx, idxLen, nVals, res )
      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in)     :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)       :: time

      !> Level on which values are requested
      integer, intent(in)                   :: iLevel

      !> Index of points in the growing array and variable val array to
      !! return.
      !! Size: most times nVals, if contiguous arrays are used it depends
      !! on the number of first indices
      integer, intent(in)                   :: idx(:)

      !> With idx as start index in contiguous memory,
      !! idxLength defines length of each contiguous memory
      !! Size: dependes on number of first index for contiguous array,
      !! but the sum of all idxLen is equal to nVals
      integer, optional, intent(in)         :: idxLen(:)

      !> Number of values to obtain for this variable (vectorized access).
      integer, intent(in)                   :: nVals

      !> Resulting values for the requested variable.
      !!
      !! Dimension: n requested entries x nComponents of this variable
      !! Access: (iElem-1)*fun%nComponents + iComp
      real(kind=rk), intent(out)            :: res(:)
    ! ---------------------------------------------------------------------------
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer  :: fnCalcPtr
    ! ---------------------------------------------------------------------------

    fnCalcPtr => mus_deriveStrainRate

    call mus_generic_varFromPDF_fromIndex( &
      &  fun       = fun,                  &
      &  varSys    = varSys,               &
      &  time      = time,                 &
      &  iLevel    = iLevel,               &
      &  idx       = idx,                  &
      &  nVals     = nVals,                &
      &  fnCalcPtr = fnCalcPtr,            &
      &  res       = res                   )

  end subroutine deriveStrainRate_fromIndex
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate Momentum from density and velocity in auxField.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveMomentum_fromIndex( fun, varSys, time, iLevel,          &
    &                                            idx, idxLen, nVals, res )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in)     :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)       :: time

    !> Level on which values are requested
    integer, intent(in)                   :: iLevel

    !> Index
    !! return.
    !! Size: most times nVals, if contiguous arrays are used it depends
    !! on the number of first indices
    integer, intent(in)                   :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: dependes on number of first index for contiguous array,
    !! but the sum of all idxLen is equal to nVals
    integer, optional, intent(in)         :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in)                   :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out)            :: res(:)
    ! ---------------------------------------------------------------------------
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos, vel_pos, iVal
    real(kind=rk) :: mass_dens(nVals), vel(3)
    ! ---------------------------------------------------------------------------


    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk

    ! get mass density values for IDX
    dens_pos = fun%input_varPos(1)
    call varSys%method%val( dens_pos )%get_ValOfIndex( &
      &     varSys = varSys,                           &
      &     time   = time,                             &
      &     iLevel = iLevel,                           &
      &     idx    = fPtr%opData%input_pntIndex(1)     &
      &              %indexLvl(iLevel)%val( idx(:) ),  &
      &     nVals  = nVals,                            &
      &     res    = mass_dens                         )

    ! get velocity values for IDX
    vel_pos = fun%input_varPos(2)
    call varSys%method%val( vel_pos )%get_ValOfIndex( &
      &     varSys = varSys,                          &
      &     time   = time,                            &
      &     iLevel = iLevel,                          &
      &     idx    = fPtr%opData%input_pntIndex(2)    &
      &              %indexLvl(iLevel)%val( idx(:) ), &
      &     nVals  = nVals,                           &
      &     res    = res(:)                           )

    associate( quantities => fPtr%solverData%scheme%layout%quantities )
      ! convert velocity to momentum
      do iVal = 1, nVals
        vel = res( (iVal-1)*3 + 1 : iVal*3 )
        res( (iVal-1)*3 + 1 : iVal*3 ) = quantities%momentum_from_vel_dens_ptr(vel=vel,dens=mass_dens(iVal))
      end do
    end associate

  end subroutine deriveMomentum_fromIndex
! ****************************************************************************** !


  ! ************************************************************************* !
  !          Subroutines with common interface called by fnCalcPtr            !
  ! ************************************************************************* !

! ****************************************************************************** !
  !> Calculate the density of a given set of elements (sum up all links).
  !! This routine is used to compute density for all scheme kinds
  !! For multispecies, it can compute both species density and mixture density
  !!
  !! The interface has to comply to the abstract interface
  !! [[mus_varSys_module:mus_derive_fromPDF]].
  !!
  recursive subroutine mus_derivedensity(fun, varsys, stencil, iLevel, &
    &                                      posInState, pdf, res, nVals )
    !> description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varsys_op_type), intent(in)     :: fun
    !> the variable system to obtain the variable from.
    type(tem_varsys_type), intent(in)         :: varsys
    !> fluid stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil
    !> current Level
    integer, intent(in)                       :: iLevel
    !> Position of element in levelwise state array
    integer, intent(in)                       :: posInState(:)
    !> pdf array
    real(kind=rk), intent(in)                 :: pdf(:)
    !> results
    real(kind=rk), intent(out)                :: res(:)
    !> nVals to get
    integer, intent(in)                       :: nVals
    ! ---------------------------------------------------------------------------
    real(kind=rk), allocatable                :: tmpPDF(:)
    real(kind=rk)                             :: dens
    integer                                   :: iDep, iVal, iComp
    integer                                   :: pdfPos, nCompPDF
    ! ---------------------------------------------------------------------------
    pdfPos = fun%input_varPos(1)
    nCompPDF = varSys%method%val(pdfPos)%nComponents
    allocate( tmpPDF( nCompPDF ) )
    res = 0.0_rk

    do iVal = 1, nVals
      dens = 0.0_rk
      do iDep = 1, fun%nInputs
        tmpPDF = pdf(                                &
          &  (iDep-1)*nVals + (iVal-1)*nCompPDF + 1  &
          &  :                                       &
          &  (iDep-1)*nVals + iVal*nCompPDF          )
        do iComp = 1, nCompPDF
          dens = dens + tmpPDF( iComp )
        end do !iComp
      end do !iDep
      res( iVal ) = dens
    end do !iVal
    deallocate( tmpPDF )
  end subroutine mus_deriveDensity
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the equlibrium of given elements with the given input state
  !! array.
  !!
  !! The interface has to comply to the abstract interface
  !! [[mus_varSys_module:mus_derive_fromPDF]].
  !!
  recursive subroutine mus_deriveEquil(fun, varsys, stencil, iLevel, &
    &                                  posInState, pdf, res, nVals )
    !> description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varsys_op_type), intent(in)     :: fun
    !> the variable system to obtain the variable from.
    type(tem_varsys_type), intent(in)         :: varsys
    !> fluid stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil
    !> current Level
    integer, intent(in)                       :: iLevel
    !> Position of element in levelwise state array
    integer, intent(in)                       :: posInState(:)
    !> pdf array
    real(kind=rk), intent(in)                 :: pdf(:)
    !> results
    real(kind=rk), intent(out)                :: res(:)
    !> nVals to get
    integer, intent(in)                       :: nVals
    ! ---------------------------------------------------------------------------
    type(mus_varSys_data_type), pointer       :: fPtr
    type(mus_scheme_type), pointer            :: scheme
    real(kind=rk), allocatable                :: tmpPDF(:)
    real(kind=rk), allocatable                :: fEq(:)
    real(kind=rk)                             :: dens, vel(3)
    integer                                   :: pdfPos, nCompsPDF, iVal
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), pointer :: quantities
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    quantities => scheme%layout%quantities
    pdfPos = fun%input_varPos(1)
    nCompsPDF = varSys%method%val( pdfPos )%nComponents
    allocate( tmpPDF( nCompsPDF ) )
    allocate( fEq( fun%nComponents ) )
    res = 0.0_rk

    do iVal = 1 , nVals
      tmpPDF = pdf( (iVal-1)*nCompsPDF+1 : iVal*nCompsPDF )
      ! computes density and velocity
      dens   = sum(tmpPDF)
      vel = quantities%vel_from_pdf_ptr(pdf = tmpPDF, dens = dens)

      fEq = quantities%pdfEq_ptr( rho = dens,    &
        &                         vel = vel,     &
        &                         QQ = nCompsPDF )

      res( (iVal-1)*fun%nComponents+1: iVal*fun%nComponents ) = fEq
    end do !iVal
    deallocate( tmpPDF )
    deallocate( fEq )

  end subroutine mus_deriveEquil
! ****************************************************************************** !


! ****************************************************************************** !
  !> For 2D only!
  !!
  !! The interface has to comply to the abstract interface
  !! [[mus_varSys_module:mus_derive_fromPDF]].
  !!
  recursive subroutine mus_deriveMoment(fun, varsys, stencil, iLevel, &
    &                                     posInState, pdf, res, nVals )
    !> description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varsys_op_type), intent(in)     :: fun
    !> the variable system to obtain the variable from.
    type(tem_varsys_type), intent(in)         :: varsys
    !> fluid stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil
    !> current Level
    integer, intent(in)                       :: iLevel
    !> Position of element in levelwise state array
    integer, intent(in)                       :: posInState(:)
    !> pdf array
    real(kind=rk), intent(in)                 :: pdf(:)
    !> results
    real(kind=rk), intent(out)                :: res(:)
    !> nVals to get
    integer, intent(in)                       :: nVals
    ! ---------------------------------------------------------------------------
    real(kind=rk), allocatable                :: tmpPDF(:)
    integer                                   :: pdfPos, nCompsPDF
    integer                                   :: iVal
    ! ---------------------------------------------------------------------------
    pdfPos = fun%input_varPos(1)
    nCompsPDF = varSys%method%val( pdfPos )%nComponents
    allocate( tmpPDF( nCompsPDF ) )
    res = 0.0_rk

    do iVal = 1, nVals
      tmpPDF = pdf( (iVal-1)*nCompsPDF+1 : iVal*nCompsPDF )

      res( (iVal-1)*9+1 ) = sum( tmpPDF ) - 1.0_rk
      res( (iVal-1)*9+2 ) = - sum(tmpPDF(1:4)) + 2.0_rk &
        &                   * sum(tmpPDF(5:8)) - 4.0_rk &
        &                   * tmpPDF(9) + 2.0_rk
      res( (iVal-1)*9+3 ) = - 2.0_rk * sum(tmpPDF(1:4)) &
        &                   + sum(tmpPDF(5:8)) + 4.0_rk &
        &                   * tmpPDF(9) - 1.0_rk
      res( (iVal-1)*9+4 ) = tmpPDF(3) + tmpPDF(7) + tmpPDF(8) &
        &                   - tmpPDF(1) - tmpPDF(5) - tmpPDF(6)
      res( (iVal-1)*9+6 ) = tmpPDF(4) + tmpPDF(6) + tmpPDF(8) &
        &                   - tmpPDF(2) - tmpPDF(5) - tmpPDF(7)
      res( (iVal-1)*9+5 ) = 2.0_rk * ( tmpPDF(1) - tmpPDf(3) ) - tmpPDF(5) &
        &                   - tmpPDF(6) + tmpPDF(7) + tmpPDF(8)
      res( (iVal-1)*9+7 ) = 2.0_rk * ( tmpPDF(2) - tmpPDF(4) ) - tmpPDf(5) &
        &                   + tmpPDF(6) - tmpPDF(7) + tmpPDF(8)
      res( (iVal-1)*9+8 ) = tmpPDF(1) - tmpPDF(2) + tmpPDF(3) - tmpPDF(4)
      res( (iVal-1)*9+9 ) = tmpPDF(5) - tmpPDF(6) - tmpPDF(7) + tmpPDF(8)
    end do ! iVal
    deallocate( tmpPDF )

  end subroutine mus_deriveMoment
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the Non-Equlibrium
  !!
  recursive subroutine mus_deriveNonEquil(fun, varsys, stencil, iLevel, &
    &                                     posInState, pdf, res, nVals )
    !> description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varsys_op_type), intent(in)     :: fun
    !> the variable system to obtain the variable from.
    type(tem_varsys_type), intent(in)         :: varsys
    !> fluid stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil
    !> current Level
    integer, intent(in)                       :: iLevel
    !> Position of element in levelwise state array
    integer, intent(in)                       :: posInState(:)
    !> pdf array
    real(kind=rk), intent(in)                 :: pdf(:)
    !> results
    real(kind=rk), intent(out)                :: res(:)
    !> nVals to get
    integer, intent(in)                       :: nVals
    ! ---------------------------------------------------------------------------
    type(mus_varSys_data_type), pointer       :: fPtr
    type(mus_scheme_type), pointer            :: scheme
    real(kind=rk), allocatable                :: tmpPDF(:)
    real(kind=rk), allocatable                :: fEq(:)
    real(kind=rk)                             :: dens, vel(3)
    integer                                   :: iVal
    integer                                   :: pdfPos, nCompsPDF
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), pointer :: quantities
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    quantities => scheme%layout%quantities

    pdfPos = fun%input_varPos(1)
    nCompsPDF = varSys%method%val( pdfPos )%nComponents
    allocate( fEq( fun%nComponents ) )
    allocate( tmpPDF( nCompsPDF ) )
    res = 0.0_rk

    do iVal = 1 , nVals
      tmpPDF = pdf( (iVal-1)*nCompsPDF+1 : iVal*nCompsPDF )
      ! computes density and velocity
      dens = sum(tmpPDF)
      vel = quantities%vel_from_pdf_ptr(pdf = tmpPDF, dens = dens)

      ! computes equilibrium
      fEq = quantities%pdfEq_ptr( rho = dens,      &
        &                         vel = vel,       &
        &                         QQ = nCompsPDF   )

      res( (iVal-1)*fun%nComponents+1 : ival*fun%nComponents ) = tmpPDF - fEq
    end do !iVal
    deallocate( tmpPDF )
    deallocate( fEq )

  end subroutine mus_deriveNonEquil
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Jiaxing Qi
  !! Calculate the strain rate ( or rate of strain, or rate of deformation)
  !!
  !! The equation is:
  !! \[
  !!  \tau_{\alpha \beta}=
  !!    -\frac{3\omega}{2\rho} \sum_{i} f^{neq}_{i} c_{i\alpha} c_{i\beta}
  !! \]
  !! where \( \tau_{\alpha \beta}\) is the stress
  !! in the \(\beta\)-direction on a face normal to the \(\alpha\)-axis,\n
  !! \( f^{neq}_i = f_i - f^{eq}_i\) is the non-equilibrium pdf.\n
  !! For more information, please refer to: equation 45 in\n
  !! Krueger T, Varnik F, Raabe D. Shear stress in lattice Boltzmann
  !! simulations. Physical Review E. 2009;79(4):1-14.\n
  !!
  !! For multi-level mesh, Omega on finer level needs to be adjusted in order to
  !! get the correct shearstress calculation.\n
  !! First, we defines c as the dx ratio between finer and coarse level.\n
  !! \( c={ \Delta dx }_{ c }/{ \Delta dx }_{ f } \)
  !! Then the viscosity on the different levels must satisfy:\n
  !! \( \frac { { \nu  }_{ f } }{ { \nu  }_{ c } } =c \)
  !! This constrain leads to a relationship of omega on different levels:\n
  !! \( {\omega}_f = \frac {1}{ {\lambda}(\frac{1}{{\omega}_c}-0.5)+0.5 } \)
  !! For more information, please refer to:\n
  !! Manuel H, Harald K, Joerg B, Sabine R. Aeroacoustic validation of the
  !! lattice boltzmann method on non-uniform grids. ECCOMAS 2012
  !!
  recursive subroutine mus_deriveStrainRate(fun, varsys, stencil, iLevel, &
      &                                     posInState, pdf, res, nVals )
    !> description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varsys_op_type), intent(in)     :: fun
    !> the variable system to obtain the variable from.
    type(tem_varsys_type), intent(in)         :: varsys
    !> fluid stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil
    !> current level
    integer, intent(in) :: iLevel
    !> Position of element in levelwise state array
    integer, intent(in) :: posInState(:)
    !> pdf array
    real(kind=rk), intent(in)                 :: pdf(:)
    !> results
    real(kind=rk), intent(out)                :: res(:)
    !> nVals to get
    integer, intent(in)                       :: nVals
    ! ---------------------------------------------------------------------------
    type(mus_varSys_data_type), pointer       :: fPtr
    type(mus_scheme_type), pointer            :: scheme
    real(kind=rk), allocatable                :: tmpPDF(:)
    real(kind=rk), allocatable                :: nonEq(:)
    real(kind=rk), allocatable                :: tau(:)
    real(kind=rk), allocatable                :: fEq(:)
    real(kind=rk)                             :: dens, vel(3), omega
    integer                                   :: pdfPos, nCompsPDF, iVal
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), pointer :: quantities
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    quantities => scheme%layout%quantities

    pdfPos = fun%input_varPos(1)
    nCompsPDF = varSys%method%val( pdfPos )%nComponents

    allocate( tmpPDF( nCompsPDF ) )
    allocate( nonEq( nCompsPDF ) )
    allocate( fEq( nCompsPDF ) )
    allocate( tau( fun%nComponents ) )

    do iVal = 1, nVals
      tmpPDF = pdf( (iVal-1)*nCompsPDF+1 : iVal*nCompsPDF )
      ! computes density and velocity
      dens   = sum(tmpPDF)
      vel = quantities%vel_from_pdf_ptr(pdf = tmpPDF, dens = dens)

      ! computes equilibrium
      fEq = quantities%pdfEq_ptr( rho = dens,       &
        &                         vel = vel,        &
        &                         QQ = nCompsPDF    )

      ! Non-Eq
      nonEq = tmpPDF - fEq

      ! get the correct omega value
      omega = scheme%field(1)%fieldProp%fluid%viscKine%omLvl(iLevel) &
        &                                    %val(posInState(iVal))

      ! compute shear stress
      ! here the 3D function works also for 2D, maybe in debug mode this
      ! will fail. cxcx is allocated according to the stencil:
      !   - 1 entry for 1D
      !   - 3 entries for 2D
      !   - 6 entries for 3D
      tau(:) = secondMom_3D( scheme%layout%fStencil%cxcx(:,:), &
        &                    nonEq(:), &
        &                    scheme%layout%fStencil%QQ )

      ! S_{ij} = - (omega/(2*rho*cs2))* sum_k(c_{ki} c{kj} f^neq_k)
      res( (iVal-1)*fun%nComponents + 1 : iVal*fun%nComponents ) &
        &  = tau(:) * (-1.5_rk) * omega * quantities%rho0Inv_ptr(dens=dens)
    end do !iVal

    deallocate( nonEq )
    deallocate( tau )
    deallocate( fEq )
    deallocate( tmpPDF )

  end subroutine mus_deriveStrainRate
! ****************************************************************************** !

end module mus_derQuan_module
