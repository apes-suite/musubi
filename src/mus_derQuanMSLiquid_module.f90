! Copyright (c) 2013-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!! This module provides the MUSUBI specific functions for calculating
!! macroscopic quantities from the state variables for multispecies liquid model.
!!
!! The depending common interface between MUSUBI and ATELES is defined in the
!! [[tem_derived_module]]. The functionality for accessing a variable from the
!! state
!! and evaluating a lua function are also provided in the
!! [[tem_derived_module]].
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
module mus_derQuanMSLiquid_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,               only: rk, long_k, globalMaxLevels, labelLen
  use tem_param_module,         only: div1_2, div1_3, div1_54, div1_9, div3_4, &
    &                                 sqrt3, cs2inv, cs2, t2cs2inv, t2cs4inv,  &
    &                                 cs4inv
  use tem_aux_module,           only: tem_abort
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
  use tem_logging_module,       only: logUnit
  use tem_topology_module,      only: tem_levelOf
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_math_module,          only: invert_matrix
  use tem_spatial_module,       only: tem_spatial_for
  use tem_subTree_type_module,  only: tem_subTree_type, tem_treeIDfrom_subTree
  use tem_debug_module,         only: dbgUnit
  use tem_operation_var_module, only: tem_evalMag_forElement,                &
    &                                 tem_evalMag_forPoint,                  &
    &                                 tem_evalMag_fromIndex,                 &
    &                                 tem_opVar_setupIndices,                &
    &                                 tem_get_new_varSys_data_ptr,           &
    &                                 tem_evalAdd_forElement,                &
    &                                 tem_evalAdd_fromIndex,                 &
    &                                 tem_opVar_setParams, tem_opVar_getParams
  use tem_spacetime_fun_module, only: tem_spacetime_for,       &
    &                                 tem_st_fun_listElem_type
  use tem_tools_module,         only: tem_PositionInSorted
  use tem_dyn_array_module,     only: PositionOfVal
  use tem_grow_array_module,    only: grw_labelarray_type, append

  ! include musubi modules
  use mus_source_type_module,    only: mus_source_op_type
  use mus_pdf_module,            only: pdf_data_type
  use mus_scheme_header_module,  only: mus_scheme_header_type
  use mus_scheme_type_module,    only: mus_scheme_type
  use mus_eNRTL_module,          only: mus_calc_thermFactor, &
    &                                  mus_calc_MS_DiffMatrix
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_varSys_module,         only: mus_varSys_data_type,        &
    &                                  mus_varSys_solverData_type,  &
    &                                  mus_get_new_solver_ptr,      &
    &                                  mus_deriveVar_ForPoint
  use mus_stateVar_module,       only: mus_accessVar_setupIndices
  use mus_operation_var_module,  only: mus_opVar_setupIndices
  use mus_derVarPos_module,      only: mus_derVarPos_type
  use mus_physics_module,            only: mus_convertFac_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  ! include aotus modules
  use aotus_module, only: flu_State

  implicit none

  private

  public :: mus_append_derVar_MSLiquid
  public :: mus_append_derMixVar_MS
  public :: deriveMoleDensityMS
  public :: deriveMoleDensityMS_fromIndex
  public :: derivePressureMS
  public :: deriveMoleFluxMS
  public :: deriveMoleFluxMS_fromIndex
  public :: deriveMoleFracMS
  public :: deriveMassFracMS
  public :: deriveVelocityMS, deriveVelocityMS_fromIndex

  public :: deriveEquilMSLiquid_FromMacro
  public :: deriveEqMSLiquid_FromState
  public :: deriveVelMSLiquid_FromState
  public :: deriveMomMSLiquid_FromState
  public :: deriveMomentaMSLiquid_FromState
  public :: deriveVelocitiesMSLiquid_FromState
  public :: deriveAuxMSLiquid_fromState
  public :: deriveAuxMSLiquid_fromState_WTDF
  public :: deriveEquilMSLiquid_fromAux

  public :: momentumFromMacroLSE, momentumFromMacroLSE_WTDF
  ! routines which uses thermodynamic factor and variable diffusivities
  ! public :: deriveVelocityWTDF_MSLiquid_FromState
  !@todo TO IMPLEMENT
!  public :: deriveEquilWTDF_MSLiquid_FromState

  ! source variables
  public :: applySrc_electricMSLiquid_2ndOrd
  public :: applySrc_electricMSLiquid_2ndOrd_WTDF
  public :: applySrc_electricMSLiquid_1stOrd
  public :: applySrc_electricMSLiquid_1stOrd_WTDF
  public :: applySrc_forceMSLiquid_2ndOrd
  public :: applySrc_forceMSLiquid_2ndOrd_WTDF
  public :: applySrc_forceMSLiquid_1stOrd
  public :: applySrc_forceMSLiquid_1stOrd_WTDF

contains

  ! ************************************************************************ !
  !> subroutine to add derive variables for multispecies-liquid
  !! (schemekind = 'multispecies_liquid') to the varsys.
  !! @todo KM: set proper velocity and equilbrium pointer for bgk_thermodynfac
  subroutine mus_append_derVar_MSLiquid( varSys, solverData, schemeHeader, &
    &                                    stencil, nFields, fldLabel,       &
    &                                    derVarName )
    ! -------------------------------------------------------------------- !
    !> global variable system
    type(tem_varSys_type), intent(inout)  :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in)  :: schemeHeader

    !> compute stencil defintion
    type(tem_stencilHeader_type), intent(in)       :: stencil

    !> number of fields
    integer, intent(in)                      :: nFields

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel(:)

    !> array of derive physical variables
    type(grw_labelarray_type), intent(inout) :: derVarName
    ! -------------------------------------------------------------------- !
    ! number of derive variables
    integer :: nDerVars, iVar, nComponents, addedPos, iIn
    integer :: iField
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

    nDerVars = 9
    allocate(derVarName_loc(nDerVars))
    derVarName_loc    = [ 'pressure          ', 'velocity          ', &
      &                   'equilibrium       ', 'equilibrium_vel   ', &
      &                   'mole_density      ', 'mole_fraction     ', &
      &                   'mass_fraction     ', 'mole_flux         ', &
      &                   'vel_mag           '   ]

    ! append local derVarName to growing array.
    ! should be done only once for all fields
    do iVar = 1, nDerVars
      call append(derVarName, derVarName_loc(iVar))
    end do

    ! add variable for each field and add mixture variable later
    do iField = 1, nFields
      do iVar = 1, nDerVars

        ! set default pointers, overwrite if neccessary
        get_element => tem_varSys_getElement_dummy
        get_point => mus_deriveVar_ForPoint
        setup_indices => mus_opVar_setupIndices
        get_valOfIndex => tem_varSys_getValOfIndex_dummy
        method_data  = mus_get_new_solver_ptr(solverData)
        set_params => tem_varSys_setParams_dummy
        get_params => tem_varSys_getParams_dummy

        select case(trim(adjustl(derVarName_loc(iVar))))
        case ('mole_density')
          get_element => deriveMoleDensityMS
          get_valOfIndex => deriveMoleDensityMS_fromIndex
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'density'

        case ('mole_fraction')
          get_element => deriveMoleFracMS
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'density'

        case ('mass_fraction')
          get_element => deriveMassFracMS
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'density'

        case ('velocity')
          get_element => deriveVelocityMS
          get_valOfIndex => deriveVelocityMS_fromIndex
          nComponents = 3
          allocate(input_varname(2))
          input_varname(1) = 'density'
          input_varname(2) = 'momentum'

        case ('mole_flux')
          ! momentum / molecularWeight
          get_element => deriveMoleFluxMS
          get_valOfIndex => deriveMoleFluxMS_fromIndex
          nComponents = 3
          allocate(input_varname(1))
          input_varname(1) = 'momentum'

        case ('pressure')
          get_element => derivePressureMS
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'density'

        case ('equilibrium_vel')
          if (trim(schemeHeader%relaxation) == 'bgk_withthermodynfac' .or. &
            & trim(schemeHeader%relaxation) == 'mrt_withthermodynfac') then
            get_element => deriveEquilVelWTDF_MSLiquid
          else
            get_element => deriveEquilVelMSLiquid
          end if
          nComponents = 3
          allocate(input_varname(2))
          input_varname(1) = 'density'
          input_varname(2) = 'momentum'

        case ('equilibrium')
          if (trim(schemeHeader%relaxation) == 'bgk_withthermodynfac' .or. &
            & trim(schemeHeader%relaxation) == 'mrt_withthermodynfac') then
            get_element => deriveEquilWTDF_MSLiquid
          else
            get_element => deriveEquilMSLiquid
          end if
          nComponents = stencil%QQ
          allocate(input_varname(2))
          input_varname(1) = 'density'
          input_varname(2) = 'momentum'

        case ('vel_mag')
          get_element => tem_evalMag_forElement
          get_point => tem_evalMag_forPoint
          get_valOfIndex => tem_evalMag_fromIndex
          setup_indices => tem_opVar_setupIndices
          method_data = tem_get_new_varSys_data_ptr(method_data)
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'velocity'

        case default
          write(logUnit(1),*) 'WARNING: Unknown variable: '//&
            &                 trim(derVarName_loc(iVar))
          cycle !go to next variable
        end select

        ! update variable names with field label
        varname = trim(fldLabel(iField))//trim(adjustl(derVarName_loc(iVar)))
        do iIn = 1, size(input_varname)
          input_varname(iIn) = trim(fldLabel(iField))//trim(input_varname(iIn))
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
      end do !iVar
    end do !iField

    ! append mixture variable for every derived species variable
    call mus_append_derMixVar_MS( varSys     = varSys,     &
      &                           solverData = solverData, &
      &                           nFields    = nFields,    &
      &                           fldLabel   = fldLabel,   &
      &                           derVarName = derVarName  )

    ! append liquid mixture variables
    call mus_append_derLiquidMixVar( varSys     = varSys,     &
      &                              solverData = solverData, &
      &                              nFields    = nFields,    &
      &                              fldLabel   = fldLabel,   &
      &                              derVarName = derVarName  )

  end subroutine mus_append_derVar_MSLiquid
  ! ************************************************************************ !

  ! ************************************************************************** !
  !> Append mixture variables for multicomponent models
  subroutine mus_append_derMixVar_MS( varSys, solverData, nFields, fldLabel, &
    &                                 derVarName )
    ! -------------------------------------------------------------------- !
    !> global variable system
    type(tem_varSys_type), intent(inout)                 :: varSys
    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData
    !> number of fields
    integer, intent(in)                                  :: nFields
    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)                         :: fldLabel(:)
    !> array of derive physical variable
    type(grw_labelarray_type), intent(in)                :: derVarName
    ! -------------------------------------------------------------------- !
    ! number of derive variables
    integer :: iVar, nComponents, addedPos
    integer, allocatable :: inPos(:)
    integer :: iField
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
    ! -------------------------------------------------------------------- !

    ! append mixture variable for every derived species variable
    do iVar = 1, derVarName%nVals
      ! set pointers for mixture. mixture quantity is obtained by summing up
      ! species quantities
      get_element => tem_evalAdd_forElement
      get_point => mus_deriveVar_ForPoint
      setup_indices => tem_opVar_setupIndices
      get_valOfIndex => tem_evalAdd_fromIndex
      method_data = tem_get_new_varSys_data_ptr(method_data)
      set_params => tem_opVar_setParams
      get_params => tem_opVar_getParams

      varname = trim(adjustl(derVarName%val(iVar)))

      ! overwrite get_element and get_valOfIndex for certain mixture variable
      select case (trim(varname))
      case ('velocity')
        get_element => deriveMixVelMS
        get_valOfIndex => deriveMixVelMS_fromIndex
        method_data  = mus_get_new_solver_ptr(solverData)
        allocate(input_varname(nFields*2))
        allocate(inPos(nFields*2))
        do iField = 1, nFields
          input_varname(iField) = trim(fldlabel(iField))//'density'
          ! position of input variable in varSys
          inPos(iField) = PositionofVal(varSys%varname, input_varname(iField))
        end do
        do iField = nFields+1, nFields*2
          input_varname(iField) = trim(fldlabel(iField-nFields))//'momentum'
          ! position of input variable in varSys
          inPos(iField) = PositionofVal(varSys%varname, input_varname(iField))
        end do

        ! nComponents same as species momentum/velocity
        nComponents = 3

      case ('vel_mag')
        get_element => tem_evalMag_forElement
        get_point => tem_evalMag_forPoint
        get_valOfIndex => tem_evalMag_fromIndex
        setup_indices => tem_opVar_setupIndices
        method_data = tem_get_new_varSys_data_ptr(method_data)
        nComponents = 1
        allocate(input_varname(1))
        allocate(inPos(1))
        input_varname(1) = 'velocity'
        inPos(1) = PositionofVal(varSys%varname, input_varname(1))

      case ( 'density', 'pressure', 'mole_density', 'mole_fraction', &
        &    'mass_fraction', 'momentum', 'mole_flux'                )
        allocate(input_varname(nFields))
        allocate(inPos(nFields))
        do iField = 1, nFields
          input_varname(iField) = trim(fldlabel(iField))//trim(varname)
          ! position of input variable in varSys
          inPos(iField) = PositionofVal(varSys%varname, input_varname(iField))
        end do

        ! get nComponents from dependent variable
        if (inPos(1) > 0) then
          nComponents = varSys%method%val(inPos(1))%nComponents
        else
          write(logUnit(1),*) 'Input variable for mixture variable ' &
            &                //trim(varname)//' not found in varSys'
          call tem_abort()
        end if
      case default
        ! Not supported as a mixture variable
        cycle
      end select


      ! input variable not found in varSys. Goto next variable
      if (any(inPos <= 0)) then
        deallocate(input_varname)
        deallocate(inPos)
        cycle
      end if

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
        write(logUnit(10),*) ' Appended variable: '//trim(varname)
      else if (addedpos < 1) then
        write(logUnit(1),*) 'Error: variable '//trim(varname)// &
          &                 ' is not added to variable system'
      end if

      deallocate(input_varname)
      deallocate(inPos)
    end do !iVar

  end subroutine mus_append_derMixVar_MS
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Append mixture variables for multicomponent liquid models
  subroutine mus_append_derLiquidMixVar(varSys, solverData, nFields, &
    &                                   fldLabel, derVarName)
    ! -------------------------------------------------------------------- !
    !> global variable system
    type(tem_varSys_type), intent(inout)                 :: varSys
    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData
    !> number of fields
    integer, intent(in)                                  :: nFields
    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)                         :: fldLabel(:)
    !> array of derive physical variable
    type(grw_labelarray_type), intent(in)                :: derVarName
    ! -------------------------------------------------------------------- !
    ! number of derive variables
    integer :: iVar, nComponents, addedPos, iIn, nInputs
    integer, allocatable :: inPos(:)
    integer :: iField
    logical :: wasAdded
    ! input varname with field label
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
    integer :: nDerVars
    character(len=labelLen), allocatable :: derVarName_loc(:)
    ! -------------------------------------------------------------------- !

    nDerVars = 3
    allocate(derVarName_loc(nDerVars))
    derVarName_loc    = [ 'kinematic_pressure', 'charge_density    ', &
      &                   'current_density   '                        ]

    ! append mixture variable dedicated to liquid model
    do iVar = 1, nDerVars
      ! append local derVarName to growing array.
      call append(derVarName, derVarName_loc(iVar))
      ! set pointers for mixture. mixture quantity is obtained by summing up
      ! species quantities
      get_element => tem_varSys_getElement_dummy
      get_point => mus_deriveVar_ForPoint
      get_valOfIndex => tem_varSys_getValOfIndex_dummy
      setup_indices => mus_opVar_setupIndices
      method_data  = mus_get_new_solver_ptr(solverData)
      set_params => tem_opVar_setParams
      get_params => tem_opVar_getParams

      varname = trim(adjustl(derVarName_loc(iVar)))

      select case (trim(varname))
      case ('kinematic_pressure')
        get_element => deriveKinePressMSLiquid
        nComponents = 1
        nInputs = 1
        allocate(input_varname(nInputs))
        input_varname(1) = 'pressure'

      case ('charge_density')
        get_element => deriveChargeDensity
        get_valOfIndex => deriveChargeDensity_fromIndex
        nComponents = 1
        nInputs = nFields
        allocate(input_varname(nInputs))
        do iField = 1, nFields
          input_varname(iField) = trim(fldlabel(iField))//'density'
        end do

      case ('current_density')
        get_element => deriveCurrentDensity
        get_valOfIndex => deriveCurrentDensity_fromIndex
        nComponents = 3
        nInputs = nFields
        allocate(input_varname(nInputs))
        do iField = 1, nFields
          input_varname(iField) = trim(fldlabel(iField))//'momentum'
        end do

      case default
        write(logUnit(1),*) 'WARNING: Unknown variable: '//trim(varname)
        cycle !go to next variable
      end select

      ! position of input variable in varSys
      allocate(inPos(nInputs))
      do iIn = 1, nInputs
        inPos(iIn) = PositionofVal(varSys%varname, input_varname(iIn))
      end do

      ! input variable not found in varSys. Goto next variable
      if (any(inPos <= 0)) then
        deallocate(input_varname)
        deallocate(inPos)
        cycle
      end if

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
        write(logUnit(10),*) ' Appended variable: '//trim(varname)
      else if (addedpos < 1) then
        write(logUnit(1),*) 'Error: variable '//trim(varname)// &
          &                 ' is not added to variable system'
      end if

      deallocate(input_varname)
      deallocate(inPos)

    end do !iVar

  end subroutine mus_append_derLiquidMixVar
  ! ************************************************************************** !

  ! ************************************************************************ !
  !      Subroutines with common interface for the function pointers         !
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Calculate the number density of a given element for single species
  !! from the density stored in auxField array.
  !! Mixture number density is computed by summing species number density
  !! using tem_evalAdd_forElement
  !! Number density = density/molecular weight
  !! mixture number density = sum(number_density)
  recursive subroutine deriveMoleDensityMS(fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res                )
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
    integer :: statePos, iElem, iLevel, dens_pos, elemOff
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: iField, depField
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer, &
      &        auxField => fPtr%solverData%scheme%auxField,           &
      &        field => fPtr%solverData%scheme%field                  )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%moleDensity ) &
          & depField = iField
      end do

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state and auxField array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars

        ! position of density field in auxField array
        dens_pos = varSys%method%val( fun%input_varPos(1) ) &
          &                     %auxField_varPos(1)

        ! number density of species: mass_dens / MolWeight
        res(iElem) = auxField(iLevel)%val( elemOff + dens_pos ) &
          &        * field( depField )%fieldProp%species%molWeightInv

      end do ! iElem
    end associate

  end subroutine deriveMoleDensityMS
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate charge density of a given element for mixture.
  !! Charge density is mixture quantity to it returns same value for all
  !! species
  !! charge_density = Faraday * \sum_i z_i*density_i/molecularWeight_i
  recursive subroutine deriveChargeDensity(fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res                )
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
    integer :: statePos, iElem, iLevel, iField, dens_pos, elemOff
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: mass_dens, charge_dens
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                            &
      &        levelPointer => fPtr%solverData%geometry%levelPointer,       &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state and auxField array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars

        charge_dens = 0.0_rk
        do iField = 1, scheme%nFields
          ! position of density field in auxField array
          dens_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens = auxField(iLevel)%val( elemOff + dens_pos )

          ! charge density
          charge_dens = charge_dens + mass_dens * species(iField)%molWeightInv &
            &         * species(iField)%chargeNr
        end do !iField

        res(iElem) = charge_dens * scheme%mixture%faradayLB

      end do ! iElem
    end associate

  end subroutine deriveChargeDensity
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Current density is computed from species momentum stored in auxField array.
  !! Current density, J = charge_density*velocity = \rho_e * u
  !! = \sum_k z_k F p_k / M_k
  !! where p_k is the species momentum
  recursive subroutine deriveCurrentDensity(fun, varsys, elempos, time, tree, &
    &                                       nElems, nDofs, res                )
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
    integer :: statePos, iElem, iComp, iLevel, iField, mom_pos(3)
    type(mus_varSys_data_type), pointer :: fPtr
    ! current density
    real(kind=rk) :: curr_dens(3)
    ! species momentum
    real(kind=rk) :: momentum(3)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                            &
      &        levelPointer => fPtr%solverData%geometry%levelPointer,       &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        curr_dens = 0.0_rk
        do iField = 1, scheme%nFields
          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1)

          ! species momentum
          do iComp = 1, 3
            momentum(iComp) = auxField(iLevel)%val(                            &
              &               (statePos-1)*varSys%nAuxScalars + mom_pos(iComp) )
          end do

          curr_dens = curr_dens + momentum(:) * species(iField)%molWeightInv &
            &       * species(iField)%chargeNr
        end do !iField

        ! copy the results to the res
        res( (iElem-1)*3 + 1 : iElem*3) = curr_dens * scheme%mixture%faradayLB

      end do ! iElem
    end associate

  end subroutine deriveCurrentDensity
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate mixture kinematic pressure.
  !! This routine requires initial total mole density which is defined in
  !! mixture table.
  !! Formula to compute kinematic pressure
  !! \( p = c^2_s (\sum_k \rho_k \phi_k - min_l (m_l) n_0)/\rho_0 \)
  !! here, \( \rho_k \) - species density, \\
  !! \( \phi_k \) - species molecular weight ratio, \\
  !! \( n_0 \) - reference mixture number density,\\
  !! \( \rho_0 \) - reference density.
  !! In tracking,
  !!```lua
  !! variable = {{"kinematic_pressure"}}
  !!```
  recursive subroutine deriveKinePressMSLiquid(fun, varsys, elempos, time, &
    &                                          tree, nElems, nDofs, res    )
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
    integer :: press_pos
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: min_molWeight, const_press
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    associate( mixture => fPtr%solverData%scheme%mixture,                   &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! position of mixture pressue in glob system
      press_pos = fun%input_varPos(1)

      ! derive dependent variable, mixture pressue
      call varSys%method%val(press_pos)%get_element( varSys  = varSys,  &
        &                                            elemPos = elemPos, &
        &                                            time    = time,    &
        &                                            tree    = tree,    &
        &                                            nElems  = nElems,  &
        &                                            nDofs   = nDofs,   &
        &                                            res     = res      )

      min_molWeight = minval(species(:)%molWeight)

      const_press = min_molWeight * mixture%moleDens0 * cs2 / mixture%rho0

      ! p = cs2 (sum_k( rho_k*phi_k ) - min_l (m_l) n_0 ) / rho_0
      res = ( res / mixture%rho0 - const_press )
    end associate

  end subroutine deriveKinePressMSLiquid
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate species  pressure both for gas and liquid model
  !! In case of gas mixture, it is partial pressure where as in
  !! liquid mixture this is not valid.
  !! However, it is used to compute mixture pressure and then the
  !! kinematic_pressure from the mixture pressure.
  !! Formula to calculate pressure:
  !! \( p_k = c^2_s ( \rho_k \phi_k ) \)
  !! here, \( \rho_k \) - species density, \\
  !! \( \phi_k \) - species molecular weight ratio, \\
  recursive subroutine derivePressureMS(fun, varsys, elempos, time, tree, &
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
    integer :: statePos, iElem, iLevel, dens_Pos, iField, depField, elemOff
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: mass_dens
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer, &
      &        auxField => fPtr%solverData%scheme%auxField,           &
      &        field => fPtr%solverData%scheme%field                  )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%pressure ) &
          & depField = iField
      end do

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state and auxField array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars

        ! position of density field in auxField array
        dens_pos = varSys%method%val( fun%input_varPos(1) ) &
          &                     %auxField_varPos(1)

        ! mass density of species
        mass_dens = auxField(iLevel)%val( elemOff + dens_pos )

        ! pressue p_k = cs2 * rho_k * phi_k
        ! multiply factor cs2 after dep since for mixture pressue
        ! p = cs2 * sum_k(rho_k*phi_k)
        res(iElem) = cs2 * mass_dens                                &
          &        * field( depField )%fieldProp%species%molWeigRatio

      end do ! iElem
    end associate

  end subroutine derivePressureMS
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Calculate the velocity of a given element for single species.
  !! from the momentum and density stored in auxField array for liquid mixture.
  !! auxField was updated with momentum of untransformed PDF which was computed
  !! by solving LSE in compute kernel.
  recursive subroutine deriveVelocityMS(fun, varsys, elempos, time, tree, &
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
    integer :: statePos, iElem, iComp, iLevel, elemOff
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos, mom_pos(3)
    ! mass density of species
    real(kind=rk) :: mass_dens
    ! species equation
    real(kind=rk) :: momentum(3)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer, &
      &        auxField => fPtr%solverData%scheme%auxField            )

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars
        ! position of current field density in auxField array
        dens_pos = varSys%method%val( fun%input_varPos(1) ) &
          &                     %auxField_varPos(1)

        ! position of current field momentum in auxField array
        mom_pos = varSys%method%val( fun%input_varPos(2) ) &
          &                     %auxField_varPos(1:3)

        ! mass density of species
        mass_dens = auxField(iLevel)%val( elemOff + dens_pos )

        ! species momentum
        do iComp = 1, 3
          momentum(iComp) = auxField(iLevel)%val( elemOff + mom_pos(iComp) )
        end do

        ! compute and store velocity
        res( (iElem-1)*3 + 1 : iElem*3) = momentum / mass_dens

      end do ! iElem
    end associate

  end subroutine deriveVelocityMS
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Equilibrium velocity from density and momentum in auxField.
  !!
  recursive subroutine deriveEquilVelMSLiquid(fun, varsys, elempos, time, &
    &                                         tree, nElems, nDofs, res    )
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
    integer :: statePos, iElem, iComp, iLevel, iField, depField, elemOff
    integer :: dens_pos, mom_pos(3)
    type(mus_varSys_data_type), pointer :: fPtr
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !momentum from auxField
    real(kind=rk) :: momentum( 3 )
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    real(kind=rk) :: eqVel(3)
    !mixture info
    real(kind=rk) :: paramBInv
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                            &
      &        levelPointer => fPtr%solverData%geometry%levelPointer,       &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%equilibriumVel ) &
          & depField = iField
      end do

      paramBInv = 1.0_rk / scheme%mixture%paramB

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1)

          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1:3)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! species momentum
          do iComp = 1, 3
            momentum(iComp) = auxField(iLevel)%val( elemOff + mom_pos(iComp) )
          end do

          ! species velocity
          vel(:, iField) = momentum / mass_dens(iField)

          ! number density
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv

        end do !iField

        ! total number density Inv
        totNum_densInv = 1.0_rk/sum(num_dens(:))

        ! molefraction
        moleFraction(:) = num_dens(:)*totNum_densInv

        eqVel = equilVelFromMacro( iField       = depField,             &
          &                        moleFraction = moleFraction,         &
          &                        velocity     = vel,                  &
          &                        nFields      = scheme%nFields,       &
          &                        paramBInv    = paramBInv,            &
          &                        phi          = species(depField)     &
          &                                       %molWeigRatio,        &
          &                        resi_coeff   = species(depField)     &
          &                                       %resi_coeff(:)        )

        ! copy the results to the res
        res( (iElem-1)*3 + 1 : iElem*3) = eqVel

      end do ! iElem
    end associate

  end subroutine deriveEquilVelMSLiquid
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Equilibrium velocity from density and momentum in auxField
  !! with thermodynamic factor
  !!
  recursive subroutine deriveEquilVelWTDF_MSLiquid(fun, varsys, elempos, time, &
    &                                              tree, nElems, nDofs, res    )
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
    integer :: statePos, iElem, iComp, iLevel, iField, depField, elemOff
    integer :: dens_pos, mom_pos(3)
    type(mus_varSys_data_type), pointer :: fPtr
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !momentum from auxField
    real(kind=rk) :: momentum( 3 )
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    real(kind=rk) :: eqVel(3)
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: &
      & resi_coeff, thermodynamic_fac, &
      & inv_thermodyn_fac, diff_coeff
    !mixture info
    real(kind=rk) :: paramBInv
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                            &
      &        levelPointer => fPtr%solverData%geometry%levelPointer,       &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%equilibriumVel ) &
          & depField = iField
      end do

      paramBInv = 1.0_rk / mixture%paramB

      do iField = 1, scheme%nFields
        ! diffusivity coefficients
        diff_coeff(iField, :) = species(iField)%diff_coeff(:)
      end do

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1)

          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1:3)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! species momentum
          do iComp = 1, 3
            momentum(iComp) = auxField(iLevel)%val( elemOff + mom_pos(iComp) )
          end do

          ! species velocity
          vel(:, iField) = momentum / mass_dens(iField)

          ! number density
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv

        end do !iField

        ! total number density Inv
        totNum_densInv = 1.0_rk/sum(num_dens(:))

        ! molefraction
        moleFraction(:) = num_dens(:)*totNum_densInv

        ! MS-Diff coeff matrix from C++ code
        call mus_calc_MS_DiffMatrix( nFields   = scheme%nFields,             &
          &                          temp      = mixture%temp0,              &
          &                          press     = mixture%atm_press,          &
          &                          mole_dens = num_dens*physics%moleDens0, &
          &                          D_ij_out  = diff_coeff                  )

        ! Convert to lattice unit
        resi_coeff = physics%fac(iLevel)%diffusivity/diff_coeff

        ! Thermodynamic factor from C++ code
        call mus_calc_thermFactor( nFields       = scheme%nFields,    &
          &                        temp          = mixture%temp0,     &
          &                        press         = mixture%atm_press, &
          &                        mole_frac     = moleFraction,      &
          &                        therm_factors = thermodynamic_fac  )

        ! invert thermodynamic factor
        inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

        ! compute equilibrium velocity with thermodynamic factor
        eqVel = equilVelFromMacroWTDF( iField            = depField,          &
          &                            mass_dens         = mass_dens,         &
          &                            moleFraction      = moleFraction,      &
          &                            velocity          = vel,               &
          &                            nFields           = scheme%nFields,    &
          &                            inv_thermodyn_fac = inv_thermodyn_fac, &
          &                            paramBInv         = paramBInv,         &
          &                            phi               = species(:)         &
          &                                                %molWeigRatio,     &
          &                            resi_coeff        = resi_coeff         )

        ! copy the results to the res
        res( (iElem-1)*3 + 1 : iElem*3) = eqVel

      end do ! iElem
    end associate

  end subroutine deriveEquilVelWTDF_MSLiquid
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Equilibrium from density and momentum stored in auxField
  recursive subroutine deriveEquilMSLiquid(fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res                )
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
    integer :: statePos, iElem, iComp, iLevel, iField, depField, elemOff
    integer :: dens_pos, mom_pos(3)
    type(mus_varSys_data_type), pointer :: fPtr
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !momentum from auxField
    real(kind=rk) :: momentum( 3 )
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    !mixture info
    real(kind=rk) :: paramBInv
    real(kind=rk) :: fEq(fun%nComponents)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                            &
      &        levelPointer => fPtr%solverData%geometry%levelPointer,       &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%equilibrium ) &
          & depField = iField
      end do

      paramBInv = 1.0_rk / scheme%mixture%paramB

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1)

          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1:3)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! species momentum
          do iComp = 1, 3
            momentum(iComp) = auxField(iLevel)%val( elemOff + mom_pos(iComp) )
          end do

          ! species velocity
          vel(:, iField) = momentum / mass_dens(iField)

          ! number density
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv

        end do !iField

        ! total number density Inv
        totNum_densInv = 1.0_rk/sum(num_dens(:))

        ! molefraction
        moleFraction(:) = num_dens(:)*totNum_densInv

        ! compute equilibrium from macroscopic quantities
        fEq = equilFromMacro( iField       = depField,          &
          &                   mass_dens    = mass_dens,         &
          &                   moleFraction = moleFraction,      &
          &                   velocity     = vel,               &
          &                   layout       = scheme%layout,     &
          &                   nFields      = scheme%nFields,    &
          &                   paramBInv    = paramBInv,         &
          &                   phi          = species(depField)  &
          &                                  %molWeigRatio,     &
          &                   resi_coeff   = species(depField)  &
          &                                  %resi_coeff(:),    &
          &                   theta_eq     = mixture%theta_eq   )

        ! copy the results to the res
        res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents) = fEq

      end do ! iElem
    end associate

   end subroutine deriveEquilMSLiquid
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Equilibrium from density and momentum in auxField with thermodynamic factor
  recursive subroutine deriveEquilWTDF_MSLiquid(fun, varsys, elempos, time, &
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
    integer :: statePos, iElem, iComp, iLevel, iField, depField, elemOff
    integer :: dens_pos, mom_pos(3)
    type(mus_varSys_data_type), pointer :: fPtr
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !momentum from auxField
    real(kind=rk) :: momentum( 3 )
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: &
      & resi_coeff, thermodynamic_fac, &
      & inv_thermodyn_fac, diff_coeff
    !mixture info
    real(kind=rk) :: paramBInv
    real(kind=rk) :: fEq(fun%nComponents)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                            &
      &        levelPointer => fPtr%solverData%geometry%levelPointer,       &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%equilibrium ) &
          & depField = iField
      end do

      paramBInv = 1.0_rk / mixture%paramB

      do iField = 1, scheme%nFields
        ! diffusivity coefficients
        diff_coeff(iField, :) = species(iField)%diff_coeff(:)
      end do

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1)

          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( fun%input_varPos(iField) ) &
            &                     %auxField_varPos(1:3)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! species momentum
          do iComp = 1, 3
            momentum(iComp) = auxField(iLevel)%val( elemOff + mom_pos(iComp) )
          end do

          ! species velocity
          vel(:, iField) = momentum / mass_dens(iField)

          ! number density
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv

        end do !iField

        ! total number density Inv
        totNum_densInv = 1.0_rk/sum(num_dens(:))

        ! molefraction
        moleFraction(:) = num_dens(:)*totNum_densInv

        ! MS-Diff coeff matrix from C++ code
        call mus_calc_MS_DiffMatrix( nFields   = scheme%nFields,             &
          &                          temp      = mixture%temp0,              &
          &                          press     = mixture%atm_press,          &
          &                          mole_dens = num_dens*physics%moleDens0, &
          &                          D_ij_out  = diff_coeff                  )

        ! Convert to lattice unit
        resi_coeff = physics%fac(iLevel)%diffusivity/diff_coeff

        ! Thermodynamic factor from C++ code
        call mus_calc_thermFactor( nFields       = scheme%nFields,    &
          &                        temp          = mixture%temp0,     &
          &                        press         = mixture%atm_press, &
          &                        mole_frac     = moleFraction,      &
          &                        therm_factors = thermodynamic_fac  )

        ! invert thermodynamic factor
        inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

        ! compute equilibrium from macroscopic quantities
        fEq = equilFromMacroWTDF( iField            = depField,          &
          &                       mass_dens         = mass_dens,         &
          &                       moleFraction      = moleFraction,      &
          &                       velocity          = vel,               &
          &                       inv_thermodyn_fac = inv_thermodyn_fac, &
          &                       layout            = scheme%layout,     &
          &                       nFields           = scheme%nFields,    &
          &                       paramBInv         = paramBInv,         &
          &                       phi               = species(:)         &
          &                                           %molWeigRatio,     &
          &                       resi_coeff        = resi_coeff,        &
          &                       theta_eq          = mixture%theta_eq   )

        ! copy the results to the res
        res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents) = fEq

      end do ! iElem
    end associate

   end subroutine deriveEquilWTDF_MSLiquid
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> mole fraction from density stored in auxField
  recursive subroutine deriveMoleFracMS(fun, varsys, elempos, time, tree, &
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
    integer :: statePos, iElem, iLevel, depField, iField, elemOff
    integer :: dens_pos
    type(mus_varSys_data_type), pointer :: fPtr
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer, &
      &        auxField => fPtr%solverData%scheme%auxField,           &
      &        field => fPtr%solverData%scheme%field                  )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%moleFrac ) &
          & depField = iField
      end do

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state and auxField array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars
        do iField = 1, scheme%nFields
          ! position of density field in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! number density of species: mass density / MolWeight
          num_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos ) &
            &              * field( iField )%fieldProp%species%molWeightInv
        end do !iField

        ! Mole fraction for current field
        res(iElem) = num_dens(depField) / sum( num_dens(:) )

      end do ! iElem
    end associate

  end subroutine deriveMoleFracMS
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> mass fraction from density stored in auxField
  recursive subroutine deriveMassFracMS(fun, varsys, elempos, time, tree, &
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
    integer :: statePos, iElem, iLevel, depField, iField, elemOff
    integer :: dens_pos
    type(mus_varSys_data_type), pointer :: fPtr
    !density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer, &
      &        auxField => fPtr%solverData%scheme%auxField,           &
      &        field => fPtr%solverData%scheme%field                  )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%massFrac ) &
          & depField = iField
      end do

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state and auxField array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars
        do iField = 1, scheme%nFields
          ! position of density field in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )
        end do !iField

        ! mass fraction for current field
        res(iElem) = mass_dens(depField) / sum( mass_dens(:) )
      end do ! iElem
    end associate

  end subroutine deriveMassFracMS
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Compute mole flux from momentum stored in auxField.
  !! mole flux = numDens_i*velocity_i = momentum / molWeight
  recursive subroutine deriveMoleFluxMS(fun, varsys, elempos, time, tree, &
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
    integer :: statePos, iElem, iComp, iLevel, iField, depField
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: mom_pos(3)
    ! species equation
    real(kind=rk) :: momentum(3)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer, &
      &        auxField => fPtr%solverData%scheme%auxField,           &
      &        field => fPtr%solverData%scheme%field                  )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%moleFlux ) &
          & depField = iField
      end do

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! position of current field momentum in auxField array
        mom_pos = varSys%method%val( fun%input_varPos(1) ) &
          &                     %auxField_varPos(1:3)

        ! species momentum
        do iComp = 1, 3
          momentum(iComp) = auxField(iLevel)%val(                            &
            &               (statePos-1)*varSys%nAuxScalars + mom_pos(iComp) )
        end do

        ! mole flux = momentum / molWeight
        res( (iElem-1)*3 + 1 : iElem*3) = momentum                  &
          &        * field( depField )%fieldProp%species%molWeightInv

      end do ! iElem
    end associate

  end subroutine deriveMoleFluxMS
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Calculate mixture velocity of a given element
  !! from the momentum and density stored in auxField array for liquid mixture.
  !! auxField was updated with momentum of untransformed PDF which was computed
  !! by solving LSE in compute kernel.
  recursive subroutine deriveMixVelMS(fun, varsys, elempos, time, tree, &
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
    integer :: statePos, iElem, iLevel, elemOff, iField
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos, mom_pos(3)
    ! species equation
    real(kind=rk) :: vel(3,varSys%nStateVars), mixVel(3), inv_rho
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, massFrac
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer, &
      &        auxField => fPtr%solverData%scheme%auxField            )

      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! element offset for auxField
        elemoff = (statePos-1)*varSys%nAuxScalars

        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( scheme%derVarPos(iField)%momentum ) &
            &                     %auxField_varPos(1:3)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! species momentum
          inv_rho = 1.0_rk / mass_dens(iField)
          vel(1, iField) = auxField(iLevel)%val( elemOff + mom_pos(1) )*inv_rho
          vel(2, iField) = auxField(iLevel)%val( elemOff + mom_pos(2) )*inv_rho
          vel(3, iField) = auxField(iLevel)%val( elemOff + mom_pos(3) )*inv_rho
        end do

        ! mass fraction
        massFrac = mass_dens / sum(mass_dens)

        ! mixture velocity: \sum_k y_k v_k
        mixVel(1) = dot_product(massFrac, vel(1, :))
        mixVel(2) = dot_product(massFrac, vel(2, :))
        mixVel(3) = dot_product(massFrac, vel(3, :))

        ! compute and store velocity
        res( (iElem-1)*3 + 1 : iElem*3) = mixVel(1:3)

      end do ! iElem
    end associate

  end subroutine deriveMixVelMS
  ! ************************************************************************ !


  ! ************************************************************************* !
  !         Subroutines with common interface for values from index           !
  ! ************************************************************************* !

  ! ************************************************************************** !
  !> Calculate mole density from species concentration for getValOfIndex
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveMoleDensityMS_fromIndex(fun, varSys, time,   &
    &                                                iLevel, idx, idxLen, &
    &                                                nVals, res           )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos
    integer :: iField, depField
    ! -------------------------------------------------------------------- !


    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                            &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%moleDensity ) &
          & depField = iField
      end do

      dens_pos = fun%input_varPos(1)
      ! get mass density values for IDX for iField
      call varSys%method%val( dens_pos )%get_ValOfIndex( &
        &     varSys = varSys,                           &
        &     time   = time,                             &
        &     iLevel = iLevel,                           &
        &     idx    = fPtr%opData%input_pntIndex(1)     &
        &              %indexLvl(iLevel)%val( idx(:) ),  &
        &     nVals  = nVals,                            &
        &     res    = res(:)                            )

      ! convert mass density to mole density
      res(1:nVals) = res(1:nVals) * species(depField)%molWeightInv

    end associate

  end subroutine deriveMoleDensityMS_fromIndex
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate velocity from species density and momentum in auxField
  !! for getValOfIndex
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveVelocityMS_fromIndex(fun, varSys, time, iLevel, &
    &                                             idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos, mom_pos, iVal
    real(kind=rk) :: mass_dens(nVals)
    ! -------------------------------------------------------------------- !


    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk

    ! get mass density values for IDX for iField
    dens_pos = fun%input_varPos(1)
    call varSys%method%val( dens_pos )%get_ValOfIndex( &
      &     varSys = varSys,                           &
      &     time   = time,                             &
      &     iLevel = iLevel,                           &
      &     idx    = fPtr%opData%input_pntIndex(1)     &
      &              %indexLvl(iLevel)%val( idx(:) ),  &
      &     nVals  = nVals,                            &
      &     res    = mass_dens                         )

    ! get species momentum values for IDX for iField
    mom_pos = fun%input_varPos(2)
    call varSys%method%val( mom_pos )%get_ValOfIndex( &
      &     varSys = varSys,                          &
      &     time   = time,                            &
      &     iLevel = iLevel,                          &
      &     idx    = fPtr%opData%input_pntIndex(2)    &
      &              %indexLvl(iLevel)%val( idx(:) ), &
      &     nVals  = nVals,                           &
      &     res    = res(:)                           )

    ! convert momentum to velocity
    do iVal = 1, nVals
      res( (iVal-1)*3 + 1 : iVal*3 ) =  res( (iVal-1)*3 + 1 : iVal*3 ) &
        &                            / mass_dens(iVal)
    end do

  end subroutine deriveVelocityMS_fromIndex
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculate mole flux from species momentum for getValOfIndex
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveMoleFluxMS_fromIndex(fun, varSys, time, iLevel, &
    &                                             idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: mom_pos
    integer :: iField, depField
    ! -------------------------------------------------------------------- !


    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( scheme => fPtr%solverData%scheme,                            &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! find dependent field of current variable by comparing the position
      ! of this variable against derVarPos of this variable for all fields in
      ! scheme
      depField = 0
      do iField = 1, scheme%nFields
        if ( fun%myPos == scheme%derVarpos(iField)%moleDensity ) &
          & depField = iField
      end do

      mom_pos = fun%input_varPos(1)
      ! get mass density values for IDX for iField
      call varSys%method%val( mom_pos )%get_ValOfIndex(  &
        &     varSys = varSys,                           &
        &     time   = time,                             &
        &     iLevel = iLevel,                           &
        &     idx    = fPtr%opData%input_pntIndex(1)     &
        &              %indexLvl(iLevel)%val( idx(:) ),  &
        &     nVals  = nVals,                            &
        &     res    = res(:)                            )

      ! convert mass flux to mole flux
      res(:) = res(:) * species(depField)%molWeightInv

    end associate

  end subroutine deriveMoleFluxMS_fromIndex
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate charge density from species concentration
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveChargeDensity_fromIndex(fun, varSys, time,   &
    &                                                iLevel, idx, idxLen, &
    &                                                nVals, res           )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos
    integer :: iField
    real(kind=rk) :: mass_dens(nVals)
    ! -------------------------------------------------------------------- !


    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( mixture => fPtr%solverData%scheme%mixture,                   &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      do iField = 1, fun%nInputs
        dens_pos = fun%input_varPos(iField)
        ! get mass density values for IDX for iField
        call varSys%method%val( dens_pos )%get_ValOfIndex(  &
          &     varSys = varSys,                            &
          &     time   = time,                              &
          &     iLevel = iLevel,                            &
          &     idx    = fPtr%opData%input_pntIndex(iField) &
          &              %indexLvl(iLevel)%val( idx(:) ),   &
          &     nVals  = nVals,                             &
          &     res    = mass_dens(:)                       )

        ! compute charge density
        res(1:nVals) = res(1:nVals) + mass_dens(1:nVals) &
          &          * species(iField)%molWeightInv      &
          &          * species(iField)%chargeNr
      end do

      ! Multiply faraday constant
      res(1:nVals) = res(1:nVals) * mixture%faradayLB

    end associate

  end subroutine deriveChargeDensity_fromIndex
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate current density from species momentum
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveCurrentDensity_fromIndex(fun, varSys, time,   &
    &                                                 iLevel, idx, idxLen, &
    &                                                 nVals, res           )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: mom_pos
    integer :: iField, iVal
    real(kind=rk) :: momentum(nVals*3)
    ! -------------------------------------------------------------------- !


    call C_F_POINTER( fun%method_Data, fPtr )
    res = 0.0_rk
    associate( mixture => fPtr%solverData%scheme%mixture,                   &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      do iField = 1, fun%nInputs
        mom_pos = fun%input_varPos(iField)
        ! get mass density values for IDX for iField
        call varSys%method%val( mom_pos )%get_ValOfIndex(   &
          &     varSys = varSys,                            &
          &     time   = time,                              &
          &     iLevel = iLevel,                            &
          &     idx    = fPtr%opData%input_pntIndex(iField) &
          &              %indexLvl(iLevel)%val( idx(:) ),   &
          &     nVals  = nVals,                             &
          &     res    = momentum(:)                        )

        ! compute current density
        do iVal = 1, nVals
          res(:) = res(:) + momentum(:) * species(iField)%molWeightInv &
            &    * species(iField)%chargeNr
        end do ! iVal
      end do ! iField

      ! Multiply faraday constant
      res(:) = res(:) * mixture%faradayLB

    end associate

  end subroutine deriveCurrentDensity_fromIndex
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate mixture velocity from density from species momentum
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveMixVelMS_fromIndex(fun, varSys, time, iLevel, &
    &                                           idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !

    call tem_abort('deriveMixVelMS_fromIndex is not implemented yet!')

  end subroutine deriveMixVelMS_fromIndex
  ! ************************************************************************** !

  ! ************************************************************************* !
  !         Subroutines with common interface for apply source                !
  ! ************************************************************************* !

  ! ************************************************************************ !
  !> Update state with source variable "electric_field" with 2nd order force
  !! integration.
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_electricMSLiquid_2ndOrd( fun, inState, outState, neigh, &
    &                                          auxField, nPdfSize, iLevel,    &
    &                                          varSys, time, phyConvFac,      &
    &                                          derVarPos                      )
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
    real(kind=rk) :: electricField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: EF_elem(3)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, depField, nScalars, QQ, nInputStates
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: massFrac( varSys%nStateVars )
    real(kind=rk) :: charge_dens, diffForce_cs2inv, diffForce_cs2inv_sqr
    real(kind=rk) :: omegaTerm, mixVel(3), inv_rho, ucx, uMinusCX(3)
    real(kind=rk), dimension(varSys%nStateVars) :: chargeTerm
    real(kind=rk) :: minMolWeight, forceTerm
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce, vel
    integer :: statePos, dens_pos, mom_pos(3), elemOff
    ! -------------------------------------------------------------------- !
!write(dbgUnit(1),*) 'source variable: ', trim(varSys%varname%val(fun%srcTerm_varPos))
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    associate( scheme => fPtr%solverData%scheme,                            &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! Get electrical force which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = electricField                             )

      ! convert physical to lattice
      electricField = electricField * physics%coulomb0 &
        &           / physics%fac(iLevel)%force

      ! minimum molecular weight
      minMolWeight = minval(species(:)%molWeight)

      ! constant term to multiply forcing term
      diffForce_cs2inv = minMolWeight / ( mixture%gasConst_R_LB &
        &           * mixture%temp0LB )
      diffForce_cs2inv_sqr = diffForce_cs2inv * diffForce_cs2inv

      ! omega term to multiply forceTerm
      omegaTerm = 1.0_rk                                                  &
        &       / ( 1.0_rk + mixture%relaxLvl(iLevel)%omega_diff * 0.5_rk )

      ! number of pdf states this source depends on
      ! last input is spacetime function so it is neglected
      nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      ! update source for each element
      do iElem = 1, nElems

        ! to access level wise state array
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)
        ! element offset for auxField
        elemoff = (posInTotal-1)*varSys%nAuxScalars

        ! get mass density from auxField
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! chargeTerm for each species: \rho_k z_k Faraday / M_k
          chargeTerm(iField) = mass_dens(iField)            &
            &                * species(iField)%molWeightInv &
            &                * species(iField)%chargeNr     &
            &                * mixture%faradayLB

          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( scheme%derVarPos(iField)%momentum ) &
            &                    %auxField_varPos(1:3)

          inv_rho = 1.0_rk / mass_dens(iField)
          ! species velocity
          vel(1, iField) = auxField(iLevel)%val(elemOff + mom_pos(1)) * inv_rho
          vel(2, iField) = auxField(iLevel)%val(elemOff + mom_pos(2)) * inv_rho
          vel(3, iField) = auxField(iLevel)%val(elemOff + mom_pos(3)) * inv_rho

        end do !iField

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! Mixture velocity
        mixVel(1) = dot_product(massFrac(:), vel(1,:))
        mixVel(2) = dot_product(massFrac(:), vel(2,:))
        mixVel(3) = dot_product(massFrac(:), vel(3,:))

        ! compute charge density: \sum_k \rho_k z_k Faraday / M_k
        charge_dens = sum(chargeTerm)

        ! electric field for current element
        EF_elem = electricField((iElem-1)*3+1 : iElem*3)
        ! compute electrical migrating force each species
        ! F_k = (\rho_k z_k/M_k - y_k \sum_l \rho_l z_l / M_l) F E / (RT)
        ! Above term is multiplied by minMolWeight which comes from lattice
        ! force term
        do iField = 1, scheme%nFields
          spcForce(:, iField) = EF_elem(:) * diffForce_cs2inv * cs2 &
            & * (chargeTerm(iField) - massFrac(iField) * charge_dens )
        end do

        ! compute external forcing term
        ! d^m_k = w_m*c_m*( min_a(m_a)*F_k  )
        ! F_k is diffusive forcing term
        !
        ! Update souce depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)

          do iDir = 1, QQ
            ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
              &                mixVel(:) )
            uMinusCX = scheme%layout%fStencil%cxDirRK(:, iDir) - mixVel(:)

            forceTerm = dot_product( uMinusCx * cs2inv               &
              &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
              &       * cs4inv, spcForce(:, depField)                )
            !forceTerm = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
            !  &                      spcForce(:, depField) )

            statePos                                                        &
              & = (posintotal-1)*nscalars+idir+(depfield-1)*qq

            outState( statePos ) = outState( statePos )              &
              & + omegaTerm * scheme%layout%weight( iDir ) * forceTerm
          end do ! iDir

        end do !iField
      end do !iElem
    end associate

  end subroutine applySrc_electricMSLiquid_2ndOrd
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Update state with source variable "electric_field"
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !! @todo species electricField
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_electricMSLiquid_1stOrd( fun, inState, outState, neigh, &
    &                                          auxField, nPdfSize, iLevel,    &
    &                                          varSys, time, phyConvFac,      &
    &                                          derVarPos                      )
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
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: electricField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: EF_elem(3)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, depField, nScalars, QQ, nInputStates
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: massFrac( varSys%nStateVars )
    real(kind=rk) :: charge_dens, diffForce_cs2inv
    real(kind=rk), dimension(varSys%nStateVars) :: chargeTerm
    real(kind=rk) :: minMolWeight, forceTerm
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce
    integer :: dens_pos, elemOff
    ! -------------------------------------------------------------------- !
!write(dbgUnit(1),*) 'source variable: ', trim(varSys%varname%val(fun%srcTerm_varPos))
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    associate( scheme => fPtr%solverData%scheme,                            &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! Get electrical force which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = electricField                             )

      ! convert physical to lattice
      electricField = electricField * physics%coulomb0 &
        &           / physics%fac(iLevel)%force

      ! minimum molecular weight
      minMolWeight = minval(species(:)%molWeight)

      ! constant term to multiply forcing term
      diffForce_cs2inv = minMolWeight / ( mixture%gasConst_R_LB &
        &           * mixture%temp0LB )

      ! number of pdf states this source depends on
      ! last input is spacetime function so it is neglected
      nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      ! update source for each element
      do iElem = 1, nElems

        ! to access level wise state array
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)
        ! element offset for auxField
        elemoff = (posInTotal-1)*varSys%nAuxScalars

        ! get mass density from auxField
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! chargeTerm for each species: \rho_k z_k Faraday / M_k
          chargeTerm(iField) = mass_dens(iField)            &
            &                * species(iField)%molWeightInv &
            &                * species(iField)%chargeNr     &
            &                * mixture%faradayLB
        end do !iField

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! compute charge density: \sum_k \rho_k z_k Faraday / M_k
        charge_dens = sum(chargeTerm)

        ! electric field for current element
        EF_elem = electricField((iElem-1)*3+1 : iElem*3)

        ! compute electrical migrating force each species
        ! F_k = (\rho_k z_k/M_k - y_k \sum_l \rho_l z_l / M_l) F E / (RT)
        ! Above term is multiplied by minMolWeight which comes from lattice
        ! force term
        do iField = 1, scheme%nFields
          spcForce(:, iField) = EF_elem(:) * diffForce_cs2inv        &
            & * (chargeTerm(iField) - massFrac(iField) * charge_dens )
        end do

        ! compute external forcing term
        ! d^m_k = w_m*c_m*( min_a(m_a)*F_k  )
        ! F_k is diffusive forcing term
        !
        ! Update souce depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)

          do iDir = 1, QQ
            forceTerm = scheme%layout%fStencil%cxDirRK( 1, iDir ) &
              &         * spcForce(1, depField)                   &
              &       + scheme%layout%fStencil%cxDirRK( 2, iDir ) &
              &         * spcForce(2, depField)                   &
              &       + scheme%layout%fStencil%cxDirRK( 3, iDir ) &
              &         * spcForce(3, depField)

            outState(                                                         &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & = outState(                                                   &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & + scheme%layout%weight( iDir ) * forceTerm
          end do ! iDir

        end do !iField
      end do !iElem
    end associate

  end subroutine applySrc_electricMSLiquid_1stOrd
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Update state with source variable "electric_field" with 2nd order
  !! integration of force term in LBE with thermodynamic factor
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_electricMSLiquid_2ndOrd_WTDF( fun, inState, outState,    &
    &                                               neigh, auxField, nPdfSize, &
    &                                               iLevel, varSys, time,      &
    &                                               phyConvFac, derVarPos      )
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
    real(kind=rk) :: electricField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: EF_elem(3)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, iField_2, depField, nScalars, QQ, nInputStates
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: massFrac( varSys%nStateVars )
    real(kind=rk) :: num_dens( varSys%nStateVars )
    real(kind=rk) :: moleFrac( varSys%nStateVars )
    real(kind=rk) :: charge_dens, diffForce_cs2inv, diffForce_cs2inv_sqr
    real(kind=rk) :: omegaTerm, mixVel(3), inv_rho, ucx, uMinusCX(3)
    real(kind=rk), dimension(varSys%nStateVars) :: chargeTerm
    real(kind=rk) :: minMolWeight, forceTerm
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce, vel, &
      & spcForce_WTDF
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: &
      & thermodynamic_fac, inv_thermodyn_fac
    integer :: statePos, dens_pos, mom_pos(3), elemOff
    ! -------------------------------------------------------------------- !
!write(dbgUnit(1),*) 'source variable: ', trim(varSys%varname%val(fun%srcTerm_varPos))
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    associate( scheme => fPtr%solverData%scheme,                            &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! Get electrical force which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = electricField                             )

      ! convert physical to lattice
      electricField = electricField * physics%coulomb0 &
        &           / physics%fac(iLevel)%force

      ! minimum molecular weight
      minMolWeight = minval(species(:)%molWeight)

      ! constant term to multiply forcing term
      diffForce_cs2inv = minMolWeight / ( mixture%gasConst_R_LB &
        &           * mixture%temp0LB )
      diffForce_cs2inv_sqr = diffForce_cs2inv * diffForce_cs2inv

      ! omega term to multiply forceTerm
      omegaTerm = 1.0_rk/(1.0_rk + mixture%relaxLvl(iLevel)%omega_diff * 0.5_rk)

      ! number of pdf states this source depends on
      ! last input is spacetime function so it is neglected
      nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      ! update source for each element
      do iElem = 1, nElems

        ! to access level wise state array
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)
        ! element offset for auxField
        elemoff = (posInTotal-1)*varSys%nAuxScalars

        ! get mass density from auxField
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! chargeTerm for each species: \rho_k z_k Faraday / M_k
          chargeTerm(iField) = mass_dens(iField)            &
            &                * species(iField)%molWeightInv &
            &                * species(iField)%chargeNr     &
            &                * mixture%faradayLB

          ! number density of species
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv

          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( scheme%derVarPos(iField)%momentum ) &
            &                    %auxField_varPos(1:3)

          inv_rho = 1.0_rk / mass_dens(iField)
          ! species velocity
          vel(1, iField) = auxField(iLevel)%val(elemOff + mom_pos(1)) * inv_rho
          vel(2, iField) = auxField(iLevel)%val(elemOff + mom_pos(2)) * inv_rho
          vel(3, iField) = auxField(iLevel)%val(elemOff + mom_pos(3)) * inv_rho

        end do !iField

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! Mixture velocity
        mixVel(1) = dot_product(massFrac(:), vel(1,:))
        mixVel(2) = dot_product(massFrac(:), vel(2,:))
        mixVel(3) = dot_product(massFrac(:), vel(3,:))

        ! compute charge density: \sum_k \rho_k z_k Faraday / M_k
        charge_dens = sum(chargeTerm)

        ! electric field for current element
        EF_elem = electricField((iElem-1)*3+1 : iElem*3)
        ! compute electrical migrating force each species
        ! F_k = (\rho_k z_k/M_k - y_k \sum_l \rho_l z_l / M_l) F E / (RT)
        ! Above term is multiplied by minMolWeight which comes from lattice
        ! force term
        do iField = 1, scheme%nFields
          spcForce(:, iField) = EF_elem(:)                           &
            & * (chargeTerm(iField) - massFrac(iField) * charge_dens )
        end do

        ! mole fraction
        moleFrac(:) = num_dens(:)/sum(num_dens)

        ! Thermodynamic factor from C++ code
        call mus_calc_thermFactor( nFields       = scheme%nFields,    &
          &                        temp          = mixture%temp0,     &
          &                        press         = mixture%atm_press, &
          &                        mole_frac     = moleFrac,          &
          &                        therm_factors = thermodynamic_fac  )

        ! invert thermodynamic factor
        inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

        ! compute external forcing term
        ! d^m_k = w_m*c_m*( \sum_l \gamma^{-1}_{k,l} min_a(m_a)*F_k  )
        ! F_k is diffusive forcing term
        spcForce_WTDF = 0.0_rk
        do iField = 1, scheme%nFields
          do iField_2 = 1, scheme%nFields
            spcForce_WTDF(:, iField ) = spcForce_WTDF(:, iField)           &
              &                      + inv_thermodyn_fac(iField, iField_2) &
              &                      * spcForce(:, iField_2)
          end do
        end do

        ! compute external forcing term
        ! d^m_k = w_m*c_m*( min_a(m_a)*F_k  )
        ! F_k is diffusive forcing term
        !
        ! Update souce depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)

          do iDir = 1, QQ
            ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
              &                mixVel )
            uMinusCX = scheme%layout%fStencil%cxDirRK(:, iDir) - mixVel

            forceTerm = dot_product( uMinusCx * diffForce_cs2inv         &
              &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir)     &
              &       * diffForce_cs2inv_sqr, spcForce_WTDF(:, depField) )

            statePos                                                        &
              & = (posintotal-1)*nscalars+idir+(depfield-1)*qq

            outState( statePos ) = outState( statePos )              &
              & + omegaTerm * scheme%layout%weight( iDir ) * forceTerm
          end do ! iDir

        end do !iField
      end do !iElem
    end associate

  end subroutine applySrc_electricMSLiquid_2ndOrd_WTDF
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Update state with source variable "electric_field" with thermodynamic
  !! factor.
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !! @todo species electricField
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_electricMSLiquid_1stOrd_WTDF( fun, inState, outState,    &
    &                                               neigh, auxField, nPdfSize, &
    &                                               iLevel, varSys, time,      &
    &                                               phyConvFac, derVarPos      )
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
    real(kind=rk) :: electricField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: EF_elem(3)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, iField_2, depField, nScalars, QQ, nInputStates
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: num_dens( varSys%nStateVars )
    real(kind=rk) :: massFrac( varSys%nStateVars )
    real(kind=rk) :: moleFrac( varSys%nStateVars )
    real(kind=rk) :: charge_dens, diffForce_cs2inv
    real(kind=rk), dimension(varSys%nStateVars) :: chargeTerm
    real(kind=rk) :: minMolWeight, forceTerm
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce, spcForce_WTDF
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: &
      & thermodynamic_fac, inv_thermodyn_fac
    integer :: dens_pos, elemOff
    ! -------------------------------------------------------------------- !
!write(dbgUnit(1),*) 'source variable: ', trim(varSys%varname%val(fun%srcTerm_varPos))
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    associate( scheme => fPtr%solverData%scheme,                            &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! Get electrical force which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = electricField                             )

      ! convert physical to lattice
      electricField = electricField * physics%coulomb0 &
        &           / physics%fac(iLevel)%force

      ! minimum molecular weight
      minMolWeight = minval(species(:)%molWeight)

      ! constant term to multiply forcing term
      diffForce_cs2inv = minMolWeight / ( mixture%gasConst_R_LB &
        &           * mixture%temp0LB )

      ! number of pdf states this source depends on
      ! last input is spacetime function so it is neglected
      nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      ! update source for each element
      do iElem = 1, nElems

        ! to access level wise state array
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)
        ! element offset for auxField
        elemoff = (posInTotal-1)*varSys%nAuxScalars

        ! get mass density from auxField
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! number density of species
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv

          ! chargeTerm for each species: \rho_k z_k Faraday / M_k
          chargeTerm(iField) = num_dens(iField) * species(iField)%chargeNr     &
            &                * mixture%faradayLB
        end do !iField

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! compute charge density: \sum_k \rho_k z_k Faraday / M_k
        charge_dens = sum(chargeTerm)

        ! electric field for current element
        EF_elem = electricField((iElem-1)*3+1 : iElem*3)

        ! compute electrical migrating force each species
        ! F_k = (\rho_k z_k/M_k - y_k \sum_l \rho_l z_l / M_l) F E / (RT)
        ! Above term is multiplied by minMolWeight which comes from lattice
        ! force term
        do iField = 1, scheme%nFields
          spcForce(:, iField) = EF_elem(:) * diffForce_cs2inv        &
            & * (chargeTerm(iField) - massFrac(iField) * charge_dens )
        end do

        ! mole fraction
        moleFrac(:) = num_dens(:)/sum(num_dens)

        ! Thermodynamic factor from C++ code
        call mus_calc_thermFactor( nFields       = scheme%nFields,    &
          &                        temp          = mixture%temp0,     &
          &                        press         = mixture%atm_press, &
          &                        mole_frac     = moleFrac,          &
          &                        therm_factors = thermodynamic_fac  )

        ! invert thermodynamic factor
        inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

        ! compute external forcing term
        ! d^m_k = w_m*c_m*( \sum_l \gamma^{-1}_{k,l} min_a(m_a)*F_k  )
        ! F_k is diffusive forcing term
        spcForce_WTDF = 0.0_rk
        do iField = 1, scheme%nFields
          do iField_2 = 1, scheme%nFields
            spcForce_WTDF(:, iField ) = spcForce_WTDF(:, iField)           &
              &                      + inv_thermodyn_fac(iField, iField_2) &
              &                      * spcForce(:, iField_2)
          end do
        end do

        ! Update souce depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)

          do iDir = 1, QQ
            forceTerm = scheme%layout%fStencil%cxDirRK( 1, iDir ) &
              &         * spcForce_WTDF(1, depField)              &
              &       + scheme%layout%fStencil%cxDirRK( 2, iDir ) &
              &         * spcForce_WTDF(2, depField)              &
              &       + scheme%layout%fStencil%cxDirRK( 3, iDir ) &
              &         * spcForce_WTDF(3, depField)

            outState(                                                         &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & = outState(                                                   &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & + scheme%layout%weight( iDir ) * forceTerm
          end do ! iDir

        end do !iField
      end do !iElem
    end associate

  end subroutine applySrc_electricMSLiquid_1stOrd_WTDF
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Update state with source variable "force" with 2nd order integration
  !! of force in lattice Boltzmann equation.
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_forceMSLiquid_2ndOrd( fun, inState, outState, neigh, &
    &                                       auxField, nPdfSize, iLevel,    &
    &                                       varSys, time, phyConvFac,      &
    &                                       derVarPos                      )
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
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, depField, nScalars, QQ, nInputStates
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: massFrac( varSys%nStateVars )
    real(kind=rk) :: forceTerm, force_elem(3), ucx, uMinusCX(3)
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce, vel
    real(kind=rk) :: omegaTerm, mixVel(3), inv_rho
    integer :: statePos, dens_pos, mom_pos(3), elemOff
    ! -------------------------------------------------------------------- !
    !write(*,*) 'source variable: ', trim(varSys%varname%val(fun%srcTerm_varPos))
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    associate( scheme => fPtr%solverData%scheme,                            &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! Get body force which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = forceField                                )

      ! convert physical to lattice
      forceField = forceField / physics%fac(iLevel)%body_force

      ! omega term to multiply forceTerm
      omegaTerm = 1.0_rk/(1.0_rk + mixture%relaxLvl(iLevel)%omega_kine * 0.5_rk)

      ! number of pdf states this source depends on
      ! last input is spacetime function so it is neglected
      nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      ! update source for each element
      do iElem = 1, nElems

        ! to access level wise state array
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

        ! element offset for auxField
        elemoff = (posInTotal-1)*varSys%nAuxScalars

        ! get mass density from auxField
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( scheme%derVarPos(iField)%momentum ) &
            &                    %auxField_varPos(1:3)

          inv_rho = 1.0_rk / mass_dens(iField)
          ! species velocity
          vel(1, iField) = auxField(iLevel)%val(elemOff + mom_pos(1)) * inv_rho
          vel(2, iField) = auxField(iLevel)%val(elemOff + mom_pos(2)) * inv_rho
          vel(3, iField) = auxField(iLevel)%val(elemOff + mom_pos(3)) * inv_rho

        end do !iField

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! Mixture velocity
        mixVel(1) = dot_product(massFrac(:), vel(1,:))
        mixVel(2) = dot_product(massFrac(:), vel(2,:))
        mixVel(3) = dot_product(massFrac(:), vel(3,:))

        ! force field for current element
        force_elem = forceField((iElem-1)*3+1 : iElem*3)

        ! compute external force for each species
        ! F_k = y_k F, F - body force per unit volume of form
        ! \rho g or \rho_e E.
        ! Above term is multiplied by cs2inv which comes from lattice
        ! force term
        do iField = 1, scheme%nFields
          spcForce(:, iField) =  massFrac(iField) * force_elem(:)
        end do

        ! compute external forcing term
        ! d^m_k = w_m*c_m*( F_k / cs2  )
        ! F_k is the external forcing term
        !
        ! Update souce depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)

          do iDir = 1, QQ
            ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
              &                mixVel )
            uMinusCX = scheme%layout%fStencil%cxDirRK(:, iDir) - mixVel

            forceTerm = dot_product( uMinusCx * cs2inv               &
              &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
              &       * cs4inv, spcForce(:, depField)                )

            statePos                                                        &
              & = (posintotal-1)*nscalars+idir+(depfield-1)*qq
            outState(statePos) = outState(statePos)                &
              & + omegaTerm * scheme%layout%weight(iDir) * forceTerm
          end do ! iDir

        end do !iField
      end do !iElem
    end associate

  end subroutine applySrc_forceMSLiquid_2ndOrd
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Update state with source variable "force" with 1st order integration
  !! of force in lattice Boltzmann equation.
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_forceMSLiquid_1stOrd( fun, inState, outState, neigh, &
    &                                       auxField, nPdfSize, iLevel,    &
    &                                       varSys, time, phyConvFac,      &
    &                                       derVarPos                      )
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
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, depField, nScalars, QQ, nInputStates
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: massFrac( varSys%nStateVars )
    real(kind=rk) :: forceTerm, force_elem(3)
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce
    integer :: dens_pos, elemOff
    ! -------------------------------------------------------------------- !
    !write(*,*) 'source variable: ', trim(varSys%varname%val(fun%srcTerm_varPos))
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    associate( scheme => fPtr%solverData%scheme,                            &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! Get body force which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = forceField                                )

      ! convert physical to lattice
      forceField = forceField / physics%fac(iLevel)%body_force

      ! number of pdf states this source depends on
      ! last input is spacetime function so it is neglected
      nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      ! update source for each element
      do iElem = 1, nElems

        ! to access level wise state array
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

        ! element offset for auxField
        elemoff = (posInTotal-1)*varSys%nAuxScalars

        ! get mass density from auxField
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

        end do !iField

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! force field for current element
        force_elem = forceField((iElem-1)*3+1 : iElem*3)

        ! compute external force for each species
        ! F_k = y_k F, F - body force per unit volume of form
        ! \rho g or \rho_e E.
        ! Above term is multiplied by cs2inv which comes from lattice
        ! force term
        do iField = 1, scheme%nFields
          spcForce(:, iField) = cs2inv * massFrac(iField) * force_elem(:)
        end do

        ! compute external forcing term
        ! d^m_k = w_m*c_m*( F_k / cs2  )
        ! F_k is the external forcing term
        !
        ! Update souce depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)

          do iDir = 1, QQ
            forceTerm =  dot_product(                               &
              &          scheme%layout%fStencil%cxDirRK( :, iDir ), &
              &          spcForce(:, depField) )

            outState(                                                         &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & = outState(                                                   &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & + scheme%layout%weight( iDir ) * forceTerm
          end do ! iDir

        end do !iField
      end do !iElem
    end associate

  end subroutine applySrc_forceMSLiquid_1stOrd
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Update state with source variable "force" with 2nd order integration
  !! of force in lattice Boltzmann equation with thermodynamic factor.
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_forceMSLiquid_2ndOrd_WTDF( fun, inState, outState,    &
    &                                            neigh, auxField, nPdfSize, &
    &                                            iLevel, varSys, time,      &
    &                                            phyConvFac, derVarPos      )
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
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, iField_2, depField, nScalars, QQ, nInputStates
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: num_dens( varSys%nStateVars )
    real(kind=rk) :: massFrac( varSys%nStateVars )
    real(kind=rk) :: moleFrac( varSys%nStateVars )
    real(kind=rk) :: forceTerm, force_elem(3), ucx, uMinusCX(3)
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce, vel, &
      & spcForce_WTDF
    real(kind=rk) :: omegaTerm, mixVel(3), inv_rho
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: &
      & thermodynamic_fac, inv_thermodyn_fac
    integer :: statePos, dens_pos, mom_pos(3), elemOff
    ! -------------------------------------------------------------------- !
    !write(*,*) 'source variable: ', trim(varSys%varname%val(fun%srcTerm_varPos))
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    associate( scheme => fPtr%solverData%scheme,                            &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! Get body force which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = forceField                                )

      ! convert physical to lattice
      forceField = forceField / physics%fac(iLevel)%body_force

      ! omega term to multiply forceTerm
      omegaTerm = 1.0_rk/(1.0_rk + mixture%relaxLvl(iLevel)%omega_kine * 0.5_rk)

      ! number of pdf states this source depends on
      ! last input is spacetime function so it is neglected
      nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      ! update source for each element
      do iElem = 1, nElems

        ! to access level wise state array
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

        ! element offset for auxField
        elemoff = (posInTotal-1)*varSys%nAuxScalars

        ! get mass density from auxField
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! position of current field momentum in auxField array
          mom_pos = varSys%method%val( scheme%derVarPos(iField)%momentum ) &
            &                    %auxField_varPos(1:3)

          inv_rho = 1.0_rk / mass_dens(iField)
          ! species velocity
          vel(1, iField) = auxField(iLevel)%val(elemOff + mom_pos(1)) * inv_rho
          vel(2, iField) = auxField(iLevel)%val(elemOff + mom_pos(2)) * inv_rho
          vel(3, iField) = auxField(iLevel)%val(elemOff + mom_pos(3)) * inv_rho

          ! number density of species
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv
        end do !iField

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! Mixture velocity
        mixVel(1) = dot_product(massFrac(:), vel(1,:))
        mixVel(2) = dot_product(massFrac(:), vel(2,:))
        mixVel(3) = dot_product(massFrac(:), vel(3,:))

        ! force field for current element
        force_elem = forceField((iElem-1)*3+1 : iElem*3)

        ! compute external force for each species
        ! F_k = y_k F, F - body force per unit volume of form
        ! \rho g or \rho_e E.
        ! Above term is multiplied by cs2inv which comes from lattice
        ! force term
        do iField = 1, scheme%nFields
          spcForce(:, iField) =  massFrac(iField) * force_elem(:)
        end do

        ! mole fraction
        moleFrac(:) = num_dens(:)/sum(num_dens)

        ! Thermodynamic factor from C++ code
        call mus_calc_thermFactor( nFields       = scheme%nFields,    &
          &                        temp          = mixture%temp0,     &
          &                        press         = mixture%atm_press, &
          &                        mole_frac     = moleFrac,          &
          &                        therm_factors = thermodynamic_fac  )

        ! invert thermodynamic factor
        inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

        ! compute external forcing term
        ! d^m_k = w_m*c_m*(  \gamma^{-1}_{k,l} F_k / cs2  )
        ! F_k is diffusive forcing term
        spcForce_WTDF = 0.0_rk
        do iField = 1, scheme%nFields
          do iField_2 = 1, scheme%nFields
            spcForce_WTDF(:, iField ) = spcForce_WTDF(:, iField)           &
              &                      + inv_thermodyn_fac(iField, iField_2) &
              &                      * spcForce(:, iField_2)
          end do
        end do

        ! compute external forcing term
        ! d^m_k = w_m*c_m*( F_k / cs2  )
        ! F_k is the external forcing term
        !
        ! Update souce depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)

          do iDir = 1, QQ
            ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
              &                mixVel )
            uMinusCX = scheme%layout%fStencil%cxDirRK(:, iDir) - mixVel

            forceTerm = dot_product( uMinusCx * cs2inv               &
              &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
              &       * cs4inv, spcForce_WTDF(:, depField)           )

            statePos                                                        &
              & = (posintotal-1)*nscalars+idir+(depfield-1)*qq
            outState(statePos) = outState(statePos)                &
              & + omegaTerm * scheme%layout%weight(iDir) * forceTerm
          end do ! iDir

        end do !iField
      end do !iElem
    end associate

  end subroutine applySrc_forceMSLiquid_2ndOrd_WTDF
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Update state with source variable "force" with thermodynamic factor
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_forceMSLiquid_1stOrd_WTDF( fun, inState, outState,    &
    &                                            neigh, auxField, nPdfSize, &
    &                                            iLevel, varSys, time,      &
    &                                            phyConvFac, derVarPos      )
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
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, iField_2, depField, nScalars, QQ, nInputStates
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: num_dens( varSys%nStateVars )
    real(kind=rk) :: massFrac( varSys%nStateVars )
    real(kind=rk) :: moleFrac( varSys%nStateVars )
    real(kind=rk) :: forceTerm, force_elem(3)
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce, spcForce_WTDF
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: &
      & thermodynamic_fac, inv_thermodyn_fac
    integer :: dens_pos, elemOff
    ! -------------------------------------------------------------------- !
    !write(*,*) 'source variable: ', trim(varSys%varname%val(fun%srcTerm_varPos))
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    associate( scheme => fPtr%solverData%scheme,                            &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

      ! Get body force which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = forceField                                )

      ! convert physical to lattice
      forceField = forceField / physics%fac(iLevel)%body_force

      ! number of pdf states this source depends on
      ! last input is spacetime function so it is neglected
      nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      ! update source for each element
      do iElem = 1, nElems

        ! to access level wise state array
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

        ! element offset for auxField
        elemoff = (posInTotal-1)*varSys%nAuxScalars

        ! get mass density from auxField
        do iField = 1, scheme%nFields
          ! position of current field density in auxField array
          dens_pos = varSys%method%val( scheme%derVarPos(iField)%density ) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField(iLevel)%val( elemOff + dens_pos )

          ! number density of species
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv
        end do !iField

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! mole fraction
        moleFrac(:) = num_dens(:)/sum(num_dens)

        ! Thermodynamic factor from C++ code
        call mus_calc_thermFactor( nFields       = scheme%nFields,    &
          &                        temp          = mixture%temp0,     &
          &                        press         = mixture%atm_press, &
          &                        mole_frac     = moleFrac,          &
          &                        therm_factors = thermodynamic_fac  )

        ! invert thermodynamic factor
        inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

        ! force field for current element
        force_elem = forceField((iElem-1)*3+1 : iElem*3)

        ! compute external force for each species
        ! F_k = y_k F, F - body force per unit volume of form
        ! \rho g or \rho_e E.
        ! Above term is multiplied by cs2inv which comes from lattice
        ! force term
        do iField = 1, scheme%nFields
          spcForce(:, iField) = cs2inv * massFrac(iField) * force_elem(:)
        end do

        ! compute external forcing term
        ! d^m_k = w_m*c_m*(  \gamma^{-1}_{k,l} F_k / cs2  )
        ! F_k is diffusive forcing term
        spcForce_WTDF = 0.0_rk
        do iField = 1, scheme%nFields
          do iField_2 = 1, scheme%nFields
            spcForce_WTDF(:, iField ) = spcForce_WTDF(:, iField)           &
              &                      + inv_thermodyn_fac(iField, iField_2) &
              &                      * spcForce(:, iField_2)
          end do
        end do

        ! Update souce depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)

          do iDir = 1, QQ
            forceTerm =  dot_product(                               &
              &          scheme%layout%fStencil%cxDirRK( :, iDir ), &
              &          spcForce_WTDF(:, depField) )

            outState(                                                         &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & = outState(                                                   &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & + scheme%layout%weight( iDir ) * forceTerm
          end do ! iDir

        end do !iField
      end do !iElem
    end associate

  end subroutine applySrc_forceMSLiquid_1stOrd_WTDF
  ! ************************************************************************ !


  ! ************************************************************************* !
  !         Subroutines with common interface for deriveFromMacro,            !
  !         deriveFromState and deriveFromAux                                 !
  ! ************************************************************************* !

  ! ************************************************************************ !
  !> This routine computes equilbrium from density and velocity
  !! This must comply with mus_variable_module%derive_FromMacro
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromMacro]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:equilFromMacro]] function pointer.
  subroutine deriveEquilMSLiquid_FromMacro( density, velocity, iField, nElems, &
    &                                       varSys, layout, res                )
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
    real(kind=rk) :: fEq(layout%fStencil%QQ)
    integer :: QQ, iElem, iFld, nFields, offset
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: phi, resi_coeff(varSys%nStateVars), paramBInv, theta_eq
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    QQ = layout%fStencil%QQ
    nFields = scheme%nFields
    ! molecular weight ratio
    phi = scheme%field(iField)%fieldProp%species%molWeigRatio
    ! resistivity coefficients
    resi_coeff(:) =                                     &
      & scheme%field(iField)%fieldProp%species%resi_coeff(:)

    paramBInv = 1.0_rk / scheme%mixture%paramB

    theta_eq = scheme%mixture%theta_eq

    do iElem = 1, nElems
      offset = (iElem-1)*nFields
      ! get species density and velocity of iElem
      do ifld = 1, nFields
        mass_dens(ifld) = density( offset + ifld )
        vel(:, ifld) = velocity(:,offset+ifld)
      end do
      num_dens(:) = mass_dens(:) * scheme%field(:)%fieldProp%species &
        &                                         %molWeightInv
      ! number density
      totNum_densInv = 1.0_rk/sum(num_dens)
      ! molefraction
      moleFraction = num_dens * totNum_densInv

      ! compute equilibrium from macroscopic quantities
      fEq = equilFromMacro( iField       = iField,       &
        &                   mass_dens    = mass_dens,    &
        &                   moleFraction = moleFraction, &
        &                   velocity     = vel,          &
        &                   layout       = layout,       &
        &                   nFields      = nFields,      &
        &                   paramBInv    = paramBInv,    &
        &                   phi          = phi,          &
        &                   resi_coeff   = resi_coeff,   &
        &                   theta_eq     = theta_eq      )


      res( (iElem-1)*QQ+1: iElem*QQ ) = fEq
    end do
  end subroutine deriveEquilMSLiquid_FromMacro
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine computes velocity from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveVelMSLiquid_FromState( state, iField, nElems, varSys, &
    &                                     layout, res                    )
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
    integer :: iElem, iComp, iFld
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos, QQ, nFields
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    real(kind=rk) :: vel( 3 )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    real(kind=rk) :: omega_diff, omega_fac, paramBInv
    integer :: stateVarMap(varSys%nStateVars)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    QQ = layout%fStencil%QQ
    omega_diff = scheme%mixture%omega_diff
    omega_fac = omega_diff * 0.5_rk

    paramBInv = 1.0_rk / scheme%mixture%paramB

    do iElem = 1, nElems
      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iFld))%state_varPos(iComp)

          tmpPDF(iComp) = state(                       &
            ! position of this state variable in the state array
            & ( ielem-1)* varsys%nscalars+varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! number density of species
        num_dens(iFld) = mass_dens(iFld)                      &
          & * scheme%field(iFld)%fieldProp%species%molWeightInv

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF *        &
            & scheme%layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total number density Inv
      totNum_densInv = 1.0_rk/sum(num_dens(:))

      ! molefraction
      moleFraction(:) = num_dens(:)*totNum_densInv

      ! momentum of all species
      momentum = momentumFromMacroLSE( moleFraction  = moleFraction,   &
        &                              first_moments = first_moments,  &
        &                              nFields       = nFields,        &
        &                              phi           = phi,            &
        &                              resi_coeff    = resi_coeff,     &
        &                              omega_fac     = omega_fac,      &
        &                              paramBInv     = paramBInv       )

      ! find required species velocity
      vel = momentum(:, iField) / mass_dens(iField)

      ! copy the results to the res
      res( (iElem-1)*3 + 1 : iElem*3) = vel
    end do

  end subroutine deriveVelMSLiquid_FromState
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine computes momentum from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveMomMSLiquid_FromState( state, iField, nElems, varSys, &
    &                                     layout, res                    )
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
    integer :: iElem, iComp, iFld
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos, QQ, nFields
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    real(kind=rk) :: omega_diff, omega_fac, paramBInv
    integer :: stateVarMap(varSys%nStateVars)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    QQ = layout%fStencil%QQ
    omega_diff = scheme%mixture%omega_diff
    omega_fac = omega_diff * 0.5_rk

    paramBInv = 1.0_rk / scheme%mixture%paramB

    do iElem = 1, nElems
      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iFld))%state_varPos(iComp)

          tmpPDF(iComp) = state(                       &
            ! position of this state variable in the state array
            & ( ielem-1)* varsys%nscalars+varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! number density of species
        num_dens(iFld) = mass_dens(iFld)                      &
          & * scheme%field(iFld)%fieldProp%species%molWeightInv

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF *   &
            & scheme%layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total number density Inv
      totNum_densInv = 1.0_rk/sum(num_dens(:))

      ! molefraction
      moleFraction(:) = num_dens(:)*totNum_densInv

      ! momentum of all species
      momentum = momentumFromMacroLSE( moleFraction  = moleFraction,   &
        &                              first_moments = first_moments,  &
        &                              nFields       = nFields,        &
        &                              phi           = phi,            &
        &                              resi_coeff    = resi_coeff,     &
        &                              omega_fac     = omega_fac,      &
        &                              paramBInv     = paramBInv       )

      ! copy the results to the res
      res( (iElem-1)*3 + 1 : iElem*3) = momentum(:, iField)
    end do

  end subroutine deriveMomMSLiquid_FromState
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine computes velocities of all species from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveVelocitiesMSLiquid_FromState( state, iField, nElems, &
    &                                            varSys, layout, res    )
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
    integer :: iElem, iComp, iFld
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos, QQ, nFields
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    real(kind=rk) :: vel( 3 )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    real(kind=rk) :: omega_diff, omega_fac, paramBInv
    integer :: stateVarMap(varSys%nStateVars)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    QQ = layout%fStencil%QQ
    omega_diff = scheme%mixture%omega_diff
    omega_fac = omega_diff * 0.5_rk

    paramBInv = 1.0_rk / scheme%mixture%paramB
    res = 0.0_rk

    do iElem = 1, nElems
      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iFld))%state_varPos(iComp)

          tmpPDF(iComp) = state(                       &
            ! position of this state variable in the state array
            & ( ielem-1)* varsys%nscalars+varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! number density of species
        num_dens(iFld) = mass_dens(iFld)                      &
          & * scheme%field(iFld)%fieldProp%species%molWeightInv

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF *   &
            & scheme%layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total number density Inv
      totNum_densInv = 1.0_rk/sum(num_dens(:))

      ! molefraction
      moleFraction(:) = num_dens(:)*totNum_densInv

      ! momentum of all species
      momentum = momentumFromMacroLSE( moleFraction  = moleFraction,   &
        &                              first_moments = first_moments,  &
        &                              nFields       = nFields,        &
        &                              phi           = phi,            &
        &                              resi_coeff    = resi_coeff,     &
        &                              omega_fac     = omega_fac,      &
        &                              paramBInv     = paramBInv       )

     ! copy the results to the res
      do iFld = 1, nFields
        ! find required species velocity
        vel = momentum(:, iFld)/mass_dens(iFld)
        res( (iElem-1)*nFields*3 + (iFld-1)*3 + 1) = vel(1)
        res( (iElem-1)*nFields*3 + (iFld-1)*3 + 2) = vel(2)
        res( (iElem-1)*nFields*3 + (iFld-1)*3 + 3) = vel(3)
      end do
    end do

  end subroutine deriveVelocitiesMSLiquid_FromState
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine computes momentum of all species from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveMomentaMSLiquid_FromState( state, iField, nElems, varSys, &
    &                                         layout, res                    )
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
    integer :: iElem, iComp, iFld
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos, QQ, nFields
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    real(kind=rk) :: omega_diff, omega_fac, paramBInv
    integer :: stateVarMap(varSys%nStateVars)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    QQ = layout%fStencil%QQ
    omega_diff = scheme%mixture%omega_diff
    omega_fac = omega_diff * 0.5_rk

    paramBInv = 1.0_rk / scheme%mixture%paramB

    do iElem = 1, nElems
      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(iFld)%state_varPos(iComp)

          tmpPDF(iComp) = state(                       &
            ! position of this state variable in the state array
            & ( ielem-1)* varsys%nscalars+ varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! number density of species
        num_dens(iFld) = mass_dens(iFld)                      &
          & * scheme%field(iFld)%fieldProp%species%molWeightInv

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF *   &
            & scheme%layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total number density Inv
      totNum_densInv = 1.0_rk/sum(num_dens(:))

      ! molefraction
      moleFraction(:) = num_dens(:)*totNum_densInv

      ! momentum of all species
      momentum = momentumFromMacroLSE( moleFraction  = moleFraction,   &
        &                              first_moments = first_moments,  &
        &                              nFields       = nFields,        &
        &                              phi           = phi,            &
        &                              resi_coeff    = resi_coeff,     &
        &                              omega_fac     = omega_fac,      &
        &                              paramBInv     = paramBInv       )

     ! copy the results to the res
      do iFld = 1, nFields
        res( (iElem-1)*nFields*3 + (iFld-1)*3 + 1) = momentum(1, iFld)
        res( (iElem-1)*nFields*3 + (iFld-1)*3 + 2) = momentum(2, iFld)
        res( (iElem-1)*nFields*3 + (iFld-1)*3 + 3) = momentum(3, iFld)
      end do

    end do

  end subroutine deriveMomentaMSLiquid_FromState
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine computes equilibrium from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveEqMSLiquid_FromState( state, iField, nElems, varSys, &
    &                                    layout, res                    )
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
    integer :: iElem, iComp, QQ, iFld, nFields
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    real(kind=rk) :: omega_diff, paramBInv, theta_eq, omega_fac
    real(kind=rk) :: fEq(layout%fStencil%QQ)
    integer :: stateVarMap(varSys%nStateVars)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = layout%fStencil%QQ
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    omega_diff = scheme%mixture%omega_diff
    omega_fac = omega_diff * 0.5_rk

    paramBInv = 1.0_rk / scheme%mixture%paramB
    theta_eq = scheme%mixture%theta_eq

    do iElem = 1, nElems
      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iFld))%state_varPos(iComp)

          tmpPDF(iComp) = state(                       &
            ! position of this state variable in the state array
            & ( ielem-1)* varsys%nscalars+varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! number density of species
        num_dens(iFld) = mass_dens(iFld)                      &
          & * scheme%field(iFld)%fieldProp%species%molWeightInv

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF *   &
            & scheme%layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total number density Inv
      totNum_densInv = 1.0_rk/sum(num_dens(:))

      ! molefraction
      moleFraction(:) = num_dens(:)*totNum_densInv

      ! momentum of all species
      momentum = momentumFromMacroLSE( moleFraction  = moleFraction,   &
        &                              first_moments = first_moments,  &
        &                              nFields       = nFields,        &
        &                              phi           = phi,            &
        &                              resi_coeff    = resi_coeff,     &
        &                              omega_fac     = omega_fac,      &
        &                              paramBInv     = paramBInv       )

      !velocity of all species
      do iFld = 1, nFields
        vel( :, iFld) = momentum( :, iFld) / mass_dens(iFld)
      end do

      ! compute equilibrium from macroscopic quantities
      fEq = equilFromMacro( iField       = iField,                &
        &                   mass_dens    = mass_dens,             &
        &                   moleFraction = moleFraction,          &
        &                   velocity     = vel,                   &
        &                   layout       = scheme%layout,         &
        &                   nFields      = nFields,               &
        &                   paramBInv    = paramBInv,             &
        &                   phi          = phi(iField),           &
        &                   resi_coeff   = resi_coeff(iField, :), &
        &                   theta_eq     = theta_eq               )

      ! copy the results to the res
      res( (iElem-1)*QQ + 1 : iElem*QQ) = fEq

    end do !iElem
  end subroutine deriveEqMSLiquid_FromState
  ! ************************************************************************ !

! **************************************************************************** !
  !> This routine computes auxField 'density and velocity' of given field
  !! from state array.
  !! velocity of original PDF is computed in this routine by solving LSE.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_auxFromState]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:auxFieldFromState]] function pointer.
  subroutine deriveAuxMSLiquid_fromState( derVarPos, state, neigh, iField, &
    &                                     nElems, nSize, iLevel, stencil,  &
    &                                     varSys, auxField, quantities     )
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
    integer :: iElem, iDir, iFld, nFields, elemOff, dens_pos, mom_pos(3)
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    real(kind=rk) :: omega_fac, paramBInv
    real(kind=rk) :: pdf( stencil%QQ )
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! ------------------------------------------------------------------------ !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields

    do iFld = 1, nFields
      ! species properties
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    omega_fac = scheme%mixture%omega_diff * 0.5_rk

    paramBInv = 1.0_rk / scheme%mixture%paramB

    !NEC$ ivdep
    do iElem = 1, nElems
      !NEC$ shortloop
      do iFld = 1, nFields
        do iDir = 1, stencil%QQ
          pdf(iDir) = state(                                           &
& ( ielem-1)* varsys%nscalars+idir+( ifld-1)* stencil%qq )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum(pdf)

        ! number density of species
        num_dens(iFld) = mass_dens(iFld)                      &
          & * scheme%field(iFld)%fieldProp%species%molWeightInv

        ! momentum of all species
        first_moments(1, iFld) = sum( pdf * stencil%cxDirRK(1, :) )
        first_moments(2, iFld) = sum( pdf * stencil%cxDirRK(2, :) )
        first_moments(3, iFld) = sum( pdf * stencil%cxDirRK(3, :) )
      end do !iFld

      ! total number density Inv
      totNum_densInv = 1.0_rk/sum(num_dens(:))

      ! molefraction
      moleFraction(:) = num_dens(:)*totNum_densInv

      ! momentum of all species
      momentum = momentumFromMacroLSE( moleFraction  = moleFraction,   &
        &                              first_moments = first_moments,  &
        &                              nFields       = nFields,        &
        &                              phi           = phi,            &
        &                              resi_coeff    = resi_coeff,     &
        &                              omega_fac     = omega_fac,      &
        &                              paramBInv     = paramBInv       )

      ! position of density and momentum of current field in auxField array
      dens_pos = varSys%method%val(derVarPos%density)%auxField_varPos(1)
      mom_pos = varSys%method%val(derVarPos%momentum)%auxField_varPos(:)

      ! element offset for auxField
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! store field density
      auxField(elemOff + dens_pos) = mass_dens(iField)
      ! store field momentum
      auxField(elemOff + mom_pos(1)) = momentum(1, iField)
      auxField(elemOff + mom_pos(2)) = momentum(2, iField)
      auxField(elemOff + mom_pos(3)) = momentum(3, iField)

    end do !iElem

  end subroutine deriveAuxMSLiquid_fromState
! **************************************************************************** !

! **************************************************************************** !
  !> This routine computes auxField 'density and velocity' of given field
  !! from state array with thermodynamic factot.
  !! velocity of original PDF is computed in this routine by solving LSE.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_auxFromState]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:auxFieldFromState]] function pointer.
  subroutine deriveAuxMSLiquid_fromState_WTDF( derVarPos, state, neigh, iField,&
    &                                          nElems, nSize, iLevel,          &
    &                                          stencil, varSys, auxField, quantities )
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
    integer :: iElem, iDir, iFld, nFields, elemOff, dens_pos, mom_pos(3)
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mole fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv, press, temp, phy_moleDens_fac
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: phi
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: &
      & resi_coeff, thermodynamic_fac, &
      & inv_thermodyn_fac, diff_coeff
    !mixture info
    real(kind=rk) :: omega_fac, paramBInv
    real(kind=rk) :: pdf( stencil%QQ )
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! ------------------------------------------------------------------------ !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields

    do iFld = 1, nFields
      ! species properties
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    phy_moleDens_fac = fPtr%solverData%physics%moleDens0
    omega_fac = scheme%mixture%omega_diff * 0.5_rk

    paramBInv = 1.0_rk / scheme%mixture%paramB

    ! temperature
    temp = scheme%mixture%temp0
    ! atmospheric pressure
    press = scheme%mixture%atm_press

    !NEC$ ivdep
    do iElem = 1, nElems
      !NEC$ shortloop
      do iFld = 1, nFields
        do iDir = 1, stencil%QQ
          pdf(iDir) = state(                                           &
& ( ielem-1)* varsys%nscalars+idir+( ifld-1)* stencil%qq )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum(pdf)

        ! number density of species
        num_dens(iFld) = mass_dens(iFld)                      &
          & * scheme%field(iFld)%fieldProp%species%molWeightInv

        ! momentum of all species
        first_moments(1, iFld) = sum( pdf * stencil%cxDirRK(1, :) )
        first_moments(2, iFld) = sum( pdf * stencil%cxDirRK(2, :) )
        first_moments(3, iFld) = sum( pdf * stencil%cxDirRK(3, :) )
      end do !iFld

      ! total number density Inv
      totNum_densInv = 1.0_rk/sum(num_dens(:))

      ! molefraction
      moleFraction(:) = num_dens(:)*totNum_densInv

      ! MS-Diff coeff matrix from C++ code
      call mus_calc_MS_DiffMatrix( nFields, temp, press,                 &
        &                          num_dens*phy_moleDens_fac, diff_coeff )

      ! Convert to lattice unit
      resi_coeff = fPtr%solverData%physics%fac(iLevel)%diffusivity/diff_coeff

      ! Thermodynamic factor from C++ code
      call mus_calc_thermFactor( nFields, temp, press, moleFraction, &
        &                        thermodynamic_fac )

      inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

      ! momentum of all species
      momentum = momentumFromMacroLSE_WTDF(                                 &
        &                            moleFraction      = moleFraction,      &
        &                            first_moments     = first_moments,     &
        &                            nFields           = nFields,           &
        &                            inv_thermodyn_fac = inv_thermodyn_fac, &
        &                            phi               = phi,               &
        &                            resi_coeff        = resi_coeff,        &
        &                            omega_fac         = omega_fac,         &
        &                            paramBInv         = paramBInv          )

      ! position of density and momentum of current field in auxField array
      dens_pos = varSys%method%val(derVarPos%density)%auxField_varPos(1)
      mom_pos = varSys%method%val(derVarPos%momentum)%auxField_varPos(:)

      ! element offset for auxField
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! store field density
      auxField(elemOff + dens_pos) = mass_dens(iField)
      ! store field momentum
      auxField(elemOff + mom_pos(1)) = momentum(1, iField)
      auxField(elemOff + mom_pos(2)) = momentum(2, iField)
      auxField(elemOff + mom_pos(3)) = momentum(3, iField)

    end do !iElem

  end subroutine deriveAuxMSLiquid_fromState_WTDF
! **************************************************************************** !

  ! ************************************************************************** !
  !> This routine computes equilbrium from auxField
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_equilFromAux]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:equilFromAux]] function pointer.
  subroutine deriveEquilMSLiquid_fromAux( derVarPos, auxField, iField, nElems, &
    &                                     varSys, layout, fEq                  )
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
    integer :: iElem, QQ, iFld, nFields, elemOff, dens_pos, mom_pos(3)
    real(kind=rk) :: fEq_loc(layout%fStencil%QQ)
    real(kind=rk) :: phi, resi_coeff(varSys%nStateVars), paramBInv, theta_eq
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: num_dens( varSys%nStateVars )
    !mass fraction
    real(kind=rk) :: moleFraction( varSys%nStateVars )
    real(kind=rk) :: totNum_densInv
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    type(mus_scheme_type), pointer :: scheme
    type(mus_varSys_data_type), pointer :: fPtr
    ! ------------------------------------------------------------------------ !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = layout%fStencil%QQ

    ! molecular weight ratio
    phi = scheme%field(iField)%fieldProp%species%molWeigRatio
    ! resistivity coefficients
    resi_coeff(:) =                                     &
      & scheme%field(iField)%fieldProp%species%resi_coeff(:)

    paramBInv = 1.0_rk / scheme%mixture%paramB

    theta_eq = scheme%mixture%theta_eq

    !NEC$ ivdep
    do iElem = 1, nElems
      ! element offset
      elemOff = (iElem-1)*varSys%nAuxScalars

      ! get all field density and velocity to compute equilibrium velocity
      !NEC$ shortloop
      do iFld = 1, nFields
        ! field density
        dens_pos = varSys%method%val(scheme%derVarPos(iFld)%density) &
          &                     %auxField_varPos(1)
        mass_dens(iFld) = auxField(elemOff + dens_pos)

        ! field velocity
        mom_pos = varSys%method%val(scheme%derVarPos(iFld)%momentum) &
          &                    %auxField_varPos(:)
        vel(1, iFld) = auxField(elemOff+mom_pos(1)) / mass_dens(iFld)
        vel(2, iFld) = auxField(elemOff+mom_pos(2)) / mass_dens(iFld)
        vel(3, iFld) = auxField(elemOff+mom_pos(3)) / mass_dens(iFld)
      end do !iFld

      ! number density
      num_dens(:) = mass_dens(:) * scheme%field(:)%fieldProp%species &
        &                                         %molWeightInv
      ! number density
      totNum_densInv = 1.0_rk/sum(num_dens)
      ! molefraction
      moleFraction = num_dens * totNum_densInv

      ! compute equilibrium from macroscopic quantities
      fEq_loc = equilFromMacro( iField       = iField,       &
        &                       mass_dens    = mass_dens,    &
        &                       moleFraction = moleFraction, &
        &                       velocity     = vel,          &
        &                       layout       = layout,       &
        &                       nFields      = nFields,      &
        &                       paramBInv    = paramBInv,    &
        &                       phi          = phi,          &
        &                       resi_coeff   = resi_coeff,   &
        &                       theta_eq     = theta_eq      )

      fEq( (iElem-1)*QQ+1: iElem*QQ ) = fEq_loc
    end do
  end subroutine deriveEquilMSLiquid_fromAux
  ! ************************************************************************** !

  ! ************************************************************************* !
  !         Pure functions used in this module                                !
  ! ************************************************************************* !

  ! ************************************************************************ !
  !> Equlibrium velocity from macro
  pure function equilVelFromMacro( iField, moleFraction, velocity, nFields,    &
    &                              paramBInv, phi, resi_coeff ) result( eqVel )
    ! ---------------------------------------------------------------------------
    !> current field
    integer, intent(in) :: iField
    !> number of species
    integer, intent(in) :: nFields
    !> mole fraction of all species
    real(kind=rk), intent(in) :: moleFraction(nFields)
    !> velocity of all species
    real(kind=rk), intent(in) :: velocity(3,nFields)
    !> free parameter B
    real(kind=rk), intent(in) :: paramBInv
    !> molecular weight ratio of iField
    real(kind=rk), intent(in) :: phi
    !> resistivity coefficients of iField
    real(kind=rk), intent(in) :: resi_coeff(nFields)
    !> return equilibrium velocity
    real(kind=rk) :: eqVel(3)
    ! ---------------------------------------------------------------------------
    integer :: ifld
    ! ---------------------------------------------------------------------------
    eqVel(:) = velocity(:,iField)
    do ifld = 1, nFields
      eqVel(:) = eqVel(:) + resi_coeff(ifld)*phi*molefraction(ifld)  &
        &      * ( velocity(:,ifld) - velocity(:,iField) ) * paramBInv
    end do

  end function equilVelFromMacro
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Equlibrium velocity from macro with thermodynamic factor
  pure function equilVelFromMacroWTDF( iField, mass_dens, moleFraction,        &
    &                                  velocity, nFields, inv_thermodyn_fac,   &
    &                                  paramBInv, phi, resi_coeff )            &
    &                                  result( eqVel )
    ! -------------------------------------------------------------------- !
    !> current field
    integer, intent(in) :: iField
    !> number of species
    integer, intent(in) :: nFields
    !> mass density of all species
    real(kind=rk), intent(in) :: mass_dens(nFields)
    !> mole fraction of all species
    real(kind=rk), intent(in) :: moleFraction(nFields)
    !> velocity of all species
    real(kind=rk), intent(in) :: velocity(3,nFields)
    !> inverse of thermodynamic factor
    real(kind=rk), intent(in) :: inv_thermodyn_fac(nFields, nFields)
    !> free parameter B
    real(kind=rk), intent(in) :: paramBInv
    !> molecular weight ratio
    real(kind=rk), intent(in) :: phi(nFields)
    !> resistivity coefficients
    real(kind=rk), intent(in) :: resi_coeff(nFields, nFields)
    !> return equilibrium velocity
    real(kind=rk) :: eqVel(3)
    ! -------------------------------------------------------------------- !
    integer :: iField_2, iField_3
    ! -------------------------------------------------------------------- !
    ! computation equilibrium vel with thermodynamic factor is done on momentum
    ! space
    eqVel(:) = mass_dens(iField)*velocity( :, iField )
    do iField_2 = 1, nFields
      do iField_3 = 1, nFields
        eqVel(:) = eqVel(:) + inv_thermodyn_fac(iField, iField_2)             &
          &                 * mass_dens(iField_2)                             &
          &                 * resi_coeff( iField_2, iField_3 ) * phi(iField_2)&
          &                 * moleFraction(iField_3)                          &
          &                 * (velocity(:, iField_3) - velocity(:,iField_2))  &
          &                 * paramBInv
      end do
    end do
    eqVel = eqVel / mass_dens(iField)

  end function equilVelFromMacroWTDF
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> derive untransformed pdf velocity of species by solving system of
  !! equations of nSpecies
  function momentumFromMacroLSE( moleFraction, first_Moments, nFields, phi,    &
    &                            resi_coeff, omega_fac, paramBInv )            &
    &                            result(momentum)
    ! -------------------------------------------------------------------- !
    !> number of species
    integer, intent(in) :: nFields
    !> molefraction  of all species
    real(kind=rk), intent(in) :: moleFraction(nFields)
    !> momentum from transformed pdf of all species
    real(kind=rk), intent(in) :: first_Moments(3, nFields)
    !> free parameter B
    real(kind=rk), intent(in) :: paramBInv
    !> molecular weight ratio of all species
    real(kind=rk), intent(in) :: phi(nFields)
    !> resistivity coefficients
    real(kind=rk), intent(in) :: resi_coeff(nFields, nFields)
    !> relaxation parameter, omega_diff*0.5_rk
    real(kind=rk), intent(in) :: omega_fac
    !> return actual momentum
    real(kind=rk) :: momentum(3, nFields)
    ! -------------------------------------------------------------------- !
    integer :: ifield, ifieldDia, ifieldNonDia, icomp
    real(kind=rk) :: matrixA( nFields, nFields ), invA( nFields, nFields )
    ! -------------------------------------------------------------------- !
    ! build up the equation system for momentum
    matrixA = 0.0_rk
    do iField = 1, nFields
      ! set diagonal part
      matrixA( iField, iField ) = 1.0_rk
      do iFieldDia = 1, nFields
        matrixA( iField, iField ) = matrixA( iField, iField )     &
          & + omega_fac * resi_coeff( iField, iFieldDia )         &
          & * phi( iField ) * moleFraction( iFieldDia ) * paramBInv
      end do
      ! set nonDiagonal
      do iFieldNonDia = 1, nFields
        matrixA( iField, iFieldNonDia ) = matrixA( iField, iFieldNonDia ) &
          & - omega_fac * resi_coeff( iField, iFieldNonDia )              &
          & * phi( iFieldNonDia ) * moleFraction( iField ) * paramBInv
      end do
    end do

    invA = invert_matrix( matrixA )

    ! momentum of all species
    momentum = 0.0_rk
    do iComp = 1, 3
      momentum( iComp, : ) = matmul( invA, first_Moments( iComp, : ) )
    end do

  end function momentumFromMacroLSE
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> derive untransformed pdf velocity of species by solving system of
  !! equations of nSpecies
  function momentumFromMacroLSE_WTDF( moleFraction, first_Moments, nFields,    &
    &                                 inv_thermodyn_fac, phi, resi_coeff,      &
    &                                 omega_fac, paramBInv ) result(momentum)
    ! -------------------------------------------------------------------- !
    !> number of species
    integer, intent(in) :: nFields
    !> molefraction  of all species
    real(kind=rk), intent(in) :: moleFraction(nFields)
    !> momentum from transformed pdf of all species
    real(kind=rk), intent(in) :: first_Moments(3, nFields)
    !> inverse of thermodynamic factor
    real(kind=rk), intent(in) :: inv_thermodyn_fac(nFields, nFields)
    !> free parameter B
    real(kind=rk), intent(in) :: paramBInv
    !> molecular weight ratio of all species
    real(kind=rk), intent(in) :: phi(nFields)
    !> resistivity coefficients
    real(kind=rk), intent(in) :: resi_coeff(nFields, nFields)
    !> relaxation parameter, omega_diff*0.5_rk
    real(kind=rk), intent(in) :: omega_fac
    !> return actual momentum
    real(kind=rk) :: momentum(3, nFields)
    ! -------------------------------------------------------------------- !
    integer :: ifield, ifield_2, ifield_3, icomp
    real(kind=rk) :: matrixA( nFields, nFields ), invA( nFields, nFields )
    ! -------------------------------------------------------------------- !
    matrixA = 0.0_rk
    !build up matrix to solver LSE for actual velocity
    do iField = 1, nFields
      !set diagonal part
      matrixA(iField, iField) = 1.0_rk
      do iField_2 = 1, nFields
        do iField_3 = 1, nFields
          matrixA(iField, iField_2) = matrixA(iField, iField_2) + omega_fac  &
            &                    * inv_thermodyn_fac(iField, iField_2)       &
            &                    * resi_coeff(iField_2, iField_3)            &
            &                    * phi(iField_2) * moleFraction(iField_3)    &
            &                    * paramBInv
        end do
      end do
      !set non-diagonal part
      do iField_2 = 1, nFields
        do iField_3 = 1, nFields
          matrixA(iField, iField_3) = matrixA(iField, iField_3) - omega_fac  &
            &                    * inv_thermodyn_fac(iField, iField_2)       &
            &                    * resi_coeff(iField_2, iField_3)            &
            &                    * phi(iField_3) * moleFraction(iField_2)    &
            &                    * paramBInv
        end do
      end do
    end do

    invA = invert_matrix( matrixA )

    ! momentum of all species
    momentum = 0.0_rk
    do iComp = 1, 3
      momentum( iComp, : ) = matmul( invA, first_Moments( iComp, : ) )
    end do

  end function momentumFromMacroLSE_WTDF
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> derive equilibrium from macro
  pure function equilFromMacro( iField, mass_dens, moleFraction, velocity,     &
    &                           layout, nFields, phi, paramBInv, resi_coeff,   &
    &                           theta_eq )  result(fEq)
    ! -------------------------------------------------------------------- !
    !> current field
    integer, intent(in) :: iField
    !> number of species
    integer, intent(in) :: nFields
    !> mass density  of all species
    real(kind=rk), intent(in) :: mass_dens(nFields)
    !> molefraction  of all species
    real(kind=rk), intent(in) :: moleFraction(nFields)
    !> velocity of all species
    real(kind=rk), intent(in) :: velocity(3, nFields)
    !> scheme layout contains stencil definition and lattice weight
    type(mus_scheme_layout_type), intent(in) :: layout
    !> molecular weight ratio of iField
    real(kind=rk), intent(in) :: phi
    !> free parameter B
    real(kind=rk), intent(in) :: paramBInv
    !> resistivity coefficients
    real(kind=rk), intent(in) :: resi_coeff(nFields)
    !> parameter to tune mixture velocity in equilibrium quadratic term
    real(kind=rk), intent(in) :: theta_eq
    !> return equilibrium
    real(kind=rk) :: fEq(layout%fStencil%QQ)
    ! -------------------------------------------------------------------- !
    integer :: iDir, QQ
    real(kind=rk) :: totMass_densInv
    real(kind=rk) :: ucx, usq, ucxQuadTerm, velAvg(3), velQuadTerm(3), eqVel(3)
    !> Inverse of lattice weight ar restPosition
    real(kind=rk) :: weight0Inv
    ! -------------------------------------------------------------------- !
    QQ = layout%fStencil%QQ
    weight0Inv = 1.0_rk / layout%weight(layout%fStencil%restPosition)

    ! compute equilibrium velocity
    eqVel = equilVelFromMacro( iField       = iField,       &
      &                        moleFraction = moleFraction, &
      &                        velocity     = velocity,     &
      &                        nFields      = nFields,      &
      &                        paramBInv    = paramBInv,    &
      &                        phi          = phi,          &
      &                        resi_coeff   = resi_coeff    )

    ! total mass density inverse
    totMass_densInv = 1.0_rk/sum(mass_dens)

    ! mass averaged mixture velocity
    velAvg(1) = dot_product( mass_dens(:), velocity(1, :) )*totMass_densInv
    velAvg(2) = dot_product( mass_dens(:), velocity(2, :) )*totMass_densInv
    velAvg(3) = dot_product( mass_dens(:), velocity(3, :) )*totMass_densInv

    ! velocity in quadratic term of equilibrium
    velQuadTerm(:) = theta_eq*velAvg(:) + (1.0_rk-theta_eq)*eqVel(:)

    ! Calculate the square of velocity
    usq = dot_product( velQuadTerm,velQuadTerm ) * t2cs2inv

    do iDir = 1, QQ
      ! Velocity times lattice unit velocity
      ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), eqVel )

      ucxQuadTerm = dot_product( layout%fStencil%cxDirRK(:, iDir), &
        &                        velQuadTerm )

      ! calculate equilibrium
      fEq(iDir) = layout%weight( iDir ) * mass_dens(iField)  &
        &         * ( phi + ucx * cs2inv                     &
        &         + ucxQuadTerm * ucxQuadTerm * t2cs4inv     &
        &         - usq                                      )
    end do ! iDir

    fEq(layout%fStencil%restPosition) =                 &
      & layout%weight( layout%fStencil%restPosition )   &
      & * mass_dens(iField)                               &
      & * ( weight0Inv + (1.0_rk - weight0Inv)*phi - usq  )

  end function equilFromMacro
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> derive equilibrium from macro
  pure function equilFromMacroWTDF( iField, mass_dens, moleFraction, velocity, &
    &                               inv_thermodyn_fac, layout, nFields, phi,   &
    &                               paramBInv, resi_coeff, theta_eq )          &
    &                               result(fEq)
    ! -------------------------------------------------------------------- !
    !> current field
    integer, intent(in) :: iField
    !> number of species
    integer, intent(in) :: nFields
    !> mass density  of all species
    real(kind=rk), intent(in) :: mass_dens(nFields)
    !> molefraction  of all species
    real(kind=rk), intent(in) :: moleFraction(nFields)
    !> velocity of all species
    real(kind=rk), intent(in) :: velocity(3, nFields)
    !> inverse of thermodynamic factor
    real(kind=rk), intent(in) :: inv_thermodyn_fac(nFields, nFields)
    !> scheme layout contains stencil definition and lattice weight
    type(mus_scheme_layout_type), intent(in) :: layout
    !> molecular weight ratio of iField
    real(kind=rk), intent(in) :: phi(nFields)
    !> free parameter B
    real(kind=rk), intent(in) :: paramBInv
    !> resistivity coefficients
    real(kind=rk), intent(in) :: resi_coeff(nFields, nFields)
    !> parameter to tune mixture velocity in equilibrium quadratic term
    real(kind=rk), intent(in) :: theta_eq
    !> return equilibrium
    real(kind=rk) :: fEq(layout%fStencil%QQ)
    ! -------------------------------------------------------------------- !
    integer :: iDir, QQ
    real(kind=rk) :: totMass_densInv
    real(kind=rk) :: ucx, usq, ucxQuadTerm, velAvg(3), velQuadTerm(3), eqVel(3)
    real(kind=rk) :: weight0Inv
    ! -------------------------------------------------------------------- !
    QQ = layout%fStencil%QQ
    weight0Inv = 1.0_rk / layout%weight(layout%fStencil%restPosition)

    ! compute equilibrium velocity
    eqVel = equilVelFromMacroWTDF( iField            = iField,            &
      &                            mass_dens         = mass_dens,         &
      &                            moleFraction      = moleFraction,      &
      &                            velocity          = velocity,          &
      &                            nFields           = nFields,           &
      &                            inv_thermodyn_fac = inv_thermodyn_fac, &
      &                            paramBInv         = paramBInv,         &
      &                            phi               = phi,               &
      &                            resi_coeff        = resi_coeff         )

    ! total mass density inverse
    totMass_densInv = 1.0_rk/sum(mass_dens)

    ! mass averaged mixture velocity
    velAvg(1) = dot_product( mass_dens(:), velocity(1, :) )*totMass_densInv
    velAvg(2) = dot_product( mass_dens(:), velocity(2, :) )*totMass_densInv
    velAvg(3) = dot_product( mass_dens(:), velocity(3, :) )*totMass_densInv

    ! velocity in quadratic term of equilibrium
    velQuadTerm(:) = theta_eq*velAvg(:) + (1.0_rk-theta_eq)*eqVel(:)

    ! Calculate the square of velocity
    usq = dot_product( velQuadTerm,velQuadTerm ) * t2cs2inv

    do iDir = 1, QQ
      ! Velocity times lattice unit velocity
      ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), eqVel )

      ucxQuadTerm = dot_product( layout%fStencil%cxDirRK(:, iDir), &
        &                        velQuadTerm )

      ! calculate equilibrium
      fEq(iDir) = layout%weight( iDir ) * mass_dens(iField) &
        &         * ( phi(iField) + ucx * cs2inv            &
        &         + ucxQuadTerm * ucxQuadTerm * t2cs4inv    &
        &         - usq                                     )
    end do ! iDir

    fEq(layout%fStencil%restPosition) =                           &
      & layout%weight( layout%fStencil%restPosition )             &
      & * mass_dens(iField)                                       &
      & * ( weight0Inv + (1.0_rk - weight0Inv)*phi(iField) - usq  )

  end function equilFromMacroWTDF
  ! ************************************************************************ !

end module mus_derQuanMSLiquid_module
! **************************************************************************** !
