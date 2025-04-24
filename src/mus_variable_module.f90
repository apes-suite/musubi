! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Philipp Otte <otte@mathcces.rwth-aachen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2017, 2020 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
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
!> This module provides all possible 'pure variables' (= no prefixes) for a
!!  given kind of simulation.
!!
!! IMPORTANT NOTE: When you create a new variable, do not forget to add this
!! variable
!! name into [[mus_append_auxField]], [[mus_store_derVarPos]] and also in
!! [[mus_store_bcVarPos]] routine.
!!
module mus_variable_module

  use iso_c_binding,            only: c_loc
  ! include treelm modules
  use env_module,               only: rk, long_k, labelLen
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit
  use tem_debug_module,         only: dbgUnit
  use tem_varSys_module,        only: tem_varSys_type,                     &
    &                                 tem_varSys_append_stateVar,          &
    &                                 tem_varSys_append_derVar,            &
    &                                 tem_varSys_append_auxFieldVar,       &
    &                                 tem_varSys_proc_point,               &
    &                                 tem_varSys_proc_element,             &
    &                                 tem_varSys_proc_setParams,           &
    &                                 tem_varSys_proc_getParams,           &
    &                                 tem_varSys_proc_setupIndices,        &
    &                                 tem_varSys_proc_getValOfIndex,       &
    &                                 tem_varSys_solverData_evalElem_type, &
    &                                 tem_varSys_dump,                     &
    &                                 tem_varSys_getPoint_dummy,           &
    &                                 tem_varSys_getElement_dummy,         &
    &                                 tem_varSys_setupIndices_dummy,       &
    &                                 tem_varSys_getValOfIndex_dummy,      &
    &                                 tem_varSys_setParams_dummy,          &
    &                                 tem_varSys_getParams_dummy
  use tem_varMap_module,        only: tem_varMap_type,         &
    &                                 tem_possible_variable_type
  use tem_variable_module,      only: tem_variable_type
  use tem_dyn_array_module,     only: dyn_labelArray_type, init, append,       &
    &                                 PositionOfVal
  use tem_grow_array_module,    only: grw_labelarray_type, init, append, &
    &                                 truncate
  use tem_derived_module,       only: tem_varSys_append_luaVar
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_meshInfo_module,      only: tem_varSys_append_meshInfoVar
  use tem_spacetime_fun_module, only: tem_st_fun_linkedList_type

  ! include musubi modules
  use mus_scheme_header_module,   only: mus_scheme_header_type
  use mus_varSys_module,          only: mus_varSys_solverData_type,  &
    &                                   mus_get_new_solver_ptr,      &
    &                                   mus_deriveVar_ForPoint,      &
    &                                   mus_generic_varFromPDF_fromIndex,&
    &                                   mus_set_stFun_getElement
  use mus_stateVar_module,        only: mus_access_state_ForElement, &
    &                                   mus_stateVar_forPoint,       &
    &                                   mus_accessVar_setupIndices,  &
    &                                   mus_stateVar_fromIndex
  use mus_source_type_module,     only: mus_source_type,            &
    &                                   mus_applySrc_dummy,         &
    &                                   mus_updateSrcVar_dummy
  use mus_source_type_module,     only: mus_addSrcToAuxField_dummy
  use mus_source_var_module,      only: mus_updateSrcVar_dynSponFld
  use mus_source_var_turbChanForce_module, only: mus_updateSrcVar_turbChanForce
  use mus_field_module,           only: mus_field_type
  use mus_scheme_layout_module,   only: mus_scheme_layout_type
  use mus_derVarPos_module,       only: mus_derVarPos_type,            &
    &                                   mus_derive_FromMacro_dummy,    &
    &                                   mus_derive_FromState_dummy,    &
    &                                   mus_derive_FromPreColState_dummy
  use mus_derQuan_module,         only: mus_append_derVar_fluid,          &
    &                                   deriveEquil_FromMacro,            &
    &                                   deriveVel_FromState,              &
    &                                   deriveVel_FromPreColState,        &
    &                                   deriveEq_FromState,               &
    &                                   deriveAux_fromState,              &
    &                                   deriveEquil_fromAux,              &
    &                                   derive_absorbLayer,               &
    &                                   derive_force_MRT,                 &
    &                                   derive_force1stOrd,               &
    &                                   derive_HRRCorrection_d2q9,        &
    &                                   derive_HRRCorrection_d3q19,       &
    &                                   derive_HRRCorrection_d3q27,       &
    &                                   applySrc_absorbLayer,             &
    &                                   applySrc_absorbLayer_MRT,         &
    &                                   applySrc_absorbLayerDyn,          &
    &                                   applySrc_absorbLayerDyn_MRT,      &
    &                                   applySrc_force,                   &
    &                                   applySrc_force_GNS,               &
    &                                   applySrc_force_MRT,               &
    &                                   applySrc_force_MRT_d2q9,          &
    &                                   applySrc_force_MRT_d3q19,         &
    &                                   applySrc_force_MRT_d3q27,         &
    &                                   applySrc_turbChanForce,           &
    &                                   applySrc_turbChanForce_MRT,       &
    &                                   applySrc_turbChanForce_MRT_d2q9,  &
    &                                   applySrc_turbChanForce_MRT_d3q19, &
    &                                   applySrc_turbChanForce_MRT_d3q27, &
    &                                   applySrc_force1stOrd
  use mus_derQuanIncomp_module,   only: mus_append_derVar_fluidIncomp,     &
    &                                   derive_absorbLayerIncomp,          &
    &                                   applySrc_absorbLayerIncomp
  use mus_derQuanPS_module,       only: mus_append_derVar_lbmPS, &
    &                                   deriveEquilPS_FromMacro, &
    &                                   deriveEquilPS2ndOrder_FromMacro,   &
    &                                   derive_equalInjectionPS, &
    &                                   deriveAuxPS_fromState,   &
    &                                   deriveEquilPS_fromAux,   &
    &                                   derive_injectionPS,      &
    &                                   applySrc_injectionPS,    &
    &                                   applySrc_equalInjectionPS
  use mus_derQuanMSGas_module,    only: mus_append_derVar_MSGas,         &
    &                                   deriveAuxMSGas_fromState,        &
    &                                   deriveEquilMSGas_fromAux,        &
    &                                   deriveEquilMSGas_FromMacro,      &
    &                                   deriveVelMSGas_FromState,        &
    &                                   deriveMomMSGas_FromState,        &
    &                                   deriveVelocitiesMSGas_FromState, &
    &                                   deriveMomentaMSGas_FromState,    &
    &                                   deriveEqMSGas_FromState
  use mus_derQuanMSLiquid_module, only: mus_append_derVar_MSLiquid,          &
    &                                   mus_append_derMixVar_MS,             &
    &                                   deriveEquilMSLiquid_FromMacro,       &
    &                                   deriveVelMSLiquid_FromState,         &
    &                                   deriveMomMSLiquid_FromState,         &
    &                                   deriveVelocitiesMSLiquid_FromState,  &
    &                                   deriveMomentaMSLiquid_FromState,     &
    &                                   deriveEqMSLiquid_FromState,          &
    &                                   deriveAuxMSLiquid_fromState,         &
    &                                   deriveAuxMSLiquid_fromState_WTDF,    &
    &                                   deriveEquilMSLiquid_fromAux,         &
    &                                   applySrc_electricMSLiquid_1stOrd,      &
    &                                   applySrc_forceMSLiquid_1stOrd,         &
    &                                   applySrc_electricMSLiquid_1stOrd_WTDF, &
    &                                   applySrc_forceMSLiquid_1stOrd_WTDF,    &
    &                                   applySrc_electricMSLiquid_2ndOrd,      &
    &                                   applySrc_forceMSLiquid_2ndOrd,         &
    &                                   applySrc_electricMSLiquid_2ndOrd_WTDF, &
    &                                   applySrc_forceMSLiquid_2ndOrd_WTDF
  use mus_derQuanPoisson_module,   only: mus_append_derVar_poisson,     &
    &                                    applySrc_chargeDensity_2ndOrd, &
    &                                    applySrc_chargeDensity_1stOrd, &
    &                                    deriveSrc_chargeDensity,       &
    &                                    deriveAuxPoisson_fromState,    &
    &                                    deriveEquilPoisson_fromAux
  use mus_derQuanNernstPlanck_module,only: applySrc_electricFieldNP, &
    &                                      deriveAuxNP_fromState,    &
    &                                      deriveEquilNP_fromAux
  use mus_derQuanPhysics_module,   only: mus_append_derVar_physics
  use mus_derQuanIsothermAcEq_module, only: &
    &   mus_append_derVar_isotherm_acEq,    &
    &   deriveEquil_FromMacro_IsothermAcEq, &
    &   deriveEq_FromState_IsothermAcEq,    &
    &   deriveEquilIsoThermAcEq_fromAux,    &
    &   deriveVelocity_FromState_IsothermAcEq
  use mus_operation_var_module,       only: mus_opVar_setupIndices, &
    &                                   mus_set_opVar_getElement
  use mus_auxFieldVar_module,         only: mus_addForceToAuxField_fluid,        &
    &                                       mus_addForceToAuxField_fluid_GNS,    &
    &                                       mus_addForceToAuxField_fluidIncomp,  &
    &                                       mus_addForceToAuxField_MSL,          &
    &                                       mus_addForceToAuxField_MSL_WTDF,     &
    &                                       mus_addElectricToAuxField_MSL,       &
    &                                       mus_addElectricToAuxField_MSL_WTDF,  &
    &                                       mus_addSrcToAuxField_poisson,        &
    &                                       mus_addSponFldToAuxField_fluid,      &
    &                                       mus_addDynSponFldToAuxField_fluid,   &
    &                                       mus_access_auxFieldVar_forElement,   &
    &                                       mus_auxFieldVar_forPoint,            &
    &                                       mus_auxFieldVar_fromIndex,           &
    &                                       mus_addTurbChanForceToAuxField_fluid,&
    &                                       mus_addHRRCorrToAuxField_fluid_2D,   &
    &                                       mus_addHRRCorrToAuxField_fluid_3D
  use mus_turbulence_var_module,      only: mus_append_turbVar
  use mus_material_var_module,        only: mus_append_materialVar
  use mus_bc_var_module,              only: mus_append_bcVar

  implicit none

  private

  public :: mus_build_varSys
  public :: mus_append_stateVar
  public :: mus_append_auxField
  public :: mus_store_derVarPos
  public :: mus_append_readVarAsStateVar

contains

  ! **************************************************************************** !
  !> Build global variable system for Musubi
  subroutine mus_build_varSys( varSys, solverData, schemeHeader, stencil,  &
    &                          nFields, derVarPos, luaVar, field, globSrc, &
    &                          poss_srcVar, st_funList )
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout) :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in) :: schemeHeader

    !> Compute stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil

    !> number of fields
    integer, intent(in)                      :: nFields

    !> store position of each variable for each field and mixture
    !! size: nFields+1
    type(mus_derVarPos_type), allocatable, intent(out)  :: derVarPos(:)

    !> additional variable defined in the lua file.
    !! Function pointer for this variables depends on its varType.
    type(tem_variable_type), allocatable, intent(in)    :: luaVar(:)

    !> Field contains sources and boundary infos
    !KM: Passed complete field array as work around for GNU compiler
    !bug
    type(mus_field_type), intent(inout) :: field(:)
    !type(mus_source_type), intent(inout) :: fldSrc(:)

    !> global source
    type(mus_source_type), intent(inout) :: globSrc

    !> possible source variables
    type(tem_possible_variable_type), intent(in)  :: poss_srcVar

    !> contains spacetime functions of all variables
    type(tem_st_fun_linkedList_type), intent(out) :: st_funList
    ! ---------------------------------------------------------------------------
    integer :: nVars, iField, iWave, iDerVP
    ! array of derive variable names depends on scheme kind
    type(grw_labelarray_type) ::  derVarName
    type(tem_varSys_solverData_evalElem_type) :: solverData_evalElem
    ! ---------------------------------------------------------------------------
    write(logUnit(1),*) 'Building variable system for scheme '&
      &                 //'kind: '//trim(schemeHeader%kind)

    ! nFields > 1 only scheme kind multispecies
    if (nFields > 1) then
      allocate(derVarPos(nFields+1))
    else
      allocate(derVarPos(nFields))
    end if

    ! initialize derVarname list
    call init(me = derVarName, length=nFields)

    ! assign default to dummy pointers
    do iDerVP = 1, size(derVarPos)
      derVarPos(iDerVP)%equilFromMacro => mus_derive_FromMacro_dummy
      derVarPos(iDerVP)%velFromState => mus_derive_FromState_dummy
      derVarPos(iDerVP)%velFromPreColState => mus_derive_FromPreColState_dummy
      derVarPos(iDerVP)%momFromState => mus_derive_FromState_dummy
      derVarPos(iDerVP)%equilFromState => mus_derive_FromState_dummy
      derVarPos(iDerVP)%velocitiesFromState => mus_derive_FromState_dummy
      derVarPos(iDerVP)%MomentaFromState => mus_derive_FromState_dummy
    end do

    ! append auxField variable depending on scheme kinds
    call mus_append_auxField( varSys       = varSys,         &
      &                       solverData   = solverData,     &
      &                       schemeHeader = schemeHeader,   &
      &                       nFields      = nFields,        &
      &                       fldLabel     = field(:)%label, &
      &                       derVarName   = derVarName      )

    ! do append variables until all variables with dependent variable
    ! are added recursively. With this variables with dependent variable
    ! can be appended in arbitrary order
    nVars = 0
    iWave = 0
    do
      if (nVars == varSys%varname%nVals) EXIT
      iWave = iWave + 1
      nVars = varSys%varname%nVals
      write(logUnit(5),*) 'Current append variable wave loop: ', iWave
      write(logUnit(5),*) 'Append derive variables to varSys'
      ! append derive vars depends on scheme kind
      select case ( trim(schemeHeader%kind) )

      case ( 'fluid', 'fluid_GNS' )
        ! append derived variables
        call mus_append_derVar_fluid( varSys       = varSys,       &
          &                         solverData   = solverData,     &
          &                         schemeHeader = schemeHeader,   &
          &                         stencil      = stencil,        &
          &                         fldLabel     = field(1)%label, &
          &                         derVarname   = derVarname      )
        derVarPos(1)%equilFromMacro => deriveEquil_FromMacro
        derVarPos(1)%velFromState   => deriveVel_fromState
        derVarPos(1)%velFromPreColState => deriveVel_FromPreColState
        derVarPos(1)%equilFromState => deriveEq_fromState
        derVarPos(1)%equilFromAux => deriveEquil_fromAux
        derVarPos(1)%auxFieldFromState => deriveAux_fromState

        ! append turbulence variable if turbulence is active
        if (field(1)%fieldProp%fluid%turbulence%active) then
          call mus_append_turbVar(                                      &
            & varSys       = varSys,                                    &
            & solverData   = solverData,                                &
            & derVarName   = derVarName,                                &
            & turbConfig   = field(1)%fieldProp%fluid%turbulence%config )
        end if

      case ( 'fluid_incompressible', 'fluid_incompressible_GNS' )
        ! append derived variables
        call mus_append_derVar_fluidIncomp( varSys       = varSys,       &
          &                               solverData   = solverData,     &
          &                               schemeHeader = schemeHeader,   &
          &                               stencil      = stencil,        &
          &                               fldLabel     = field(1)%label, &
          &                               derVarname   = derVarname      )
        derVarPos(1)%equilFromState => deriveEq_FromState
        derVarPos(1)%velFromPreColState => deriveVel_FromPreColState
        derVarPos(1)%velFromState   => deriveVel_fromState
        derVarPos(1)%equilFromMacro => deriveEquil_FromMacro
        derVarPos(1)%auxFieldFromState => deriveAux_fromState
        derVarPos(1)%equilFromAux => deriveEquil_fromAux

        ! append turbulence variable if turbulence is active
        if (field(1)%fieldProp%fluid%turbulence%active) then
          call mus_append_turbVar(                                      &
            & varSys       = varSys,                                    &
            & solverData   = solverData,                                &
            & derVarName   = derVarName,                                &
            & turbConfig   = field(1)%fieldProp%fluid%turbulence%config )
        end if

      case ( 'passive_scalar' )
        call mus_append_derVar_lbmPS( varSys       = varSys,         &
          &                           solverData   = solverData,     &
          &                           fldLabel     = field(1)%label, &
          &                           derVarname   = derVarname      )
        select case (trim(schemeHeader%relaxHeader%variant))
        case ('first')
          derVarPos(1)%equilFromMacro => deriveEquilPS_FromMacro
        case ('second')
          derVarPos(1)%equilFromMacro => deriveEquilPS2ndOrder_FromMacro
        case default
          derVarPos(1)%equilFromMacro => deriveEquilPS2ndOrder_FromMacro
        end select
        derVarPos(1)%auxFieldFromState => deriveAuxPS_fromState
        derVarPos(1)%equilFromAux => deriveEquilPS_fromAux
      case ('poisson', 'poisson_boltzmann_linear', &
        &   'poisson_boltzmann_nonlinear')
        call mus_append_derVar_poisson( varSys       = varSys,           &
          &                             solverData   = solverData,       &
          &                             fldLabel     = field(1)%label,   &
          &                             stencil      = stencil,          &
          &                             derVarname   = derVarname,       &
          &                             schemeKind   = schemeHeader%kind )
        derVarPos(1)%auxFieldFromState => deriveAuxPoisson_fromState
        derVarPos(1)%equilFromAux => deriveEquilPoisson_fromAux
      case ('nernst_planck')
        ! Nernst_planck model has only mole_density as derived variable
        ! and this variable is append to varSys as auxField
        derVarPos(1)%auxFieldFromState => deriveAuxNP_fromState
        derVarPos(1)%equilFromAux => deriveEquilNP_fromAux
      case ( 'multispecies_gas' )
        call mus_append_derVar_MSGas( varSys     = varSys,         &
          &                           solverData = solverData,     &
          &                           stencil    = stencil,        &
          &                           nFields    = nFields,        &
          &                           fldLabel   = field(:)%label, &
          &                           derVarname = derVarname      )

        do iField = 1, nFields
          derVarPos(iField)%equilFromMacro => deriveEquilMSGas_FromMacro
          derVarPos(iField)%velFromState   => deriveVelMSGas_fromState
          derVarPos(iField)%momFromState   => deriveMomMSGas_fromState
          derVarPos(iField)%equilFromState => deriveEqMSGas_fromState
          derVarPos(iField)%momentaFromState  => deriveMomentaMSGas_fromState
          derVarPos(iField)%velocitiesFromState                             &
            &                              => deriveVelocitiesMSGas_fromState
          derVarPos(iField)%auxFieldFromState => deriveAuxMSGas_fromState
          derVarPos(iField)%equilFromAux => deriveEquilMSGas_fromAux
        end do
      case ( 'multispecies_liquid' )
        call mus_append_derVar_MSLiquid( varSys       = varSys,         &
          &                              solverData   = solverData,     &
          &                              schemeHeader = schemeHeader,   &
          &                              stencil      = stencil,        &
          &                              nFields      = nFields,        &
          &                              fldLabel     = field(:)%label, &
          &                              derVarname   = derVarname      )

        do iField = 1, nFields
          derVarPos(iField)%equilFromMacro => deriveEquilMSLiquid_FromMacro
          derVarPos(iField)%velFromState   => deriveVelMSLiquid_fromState
          derVarPos(iField)%momFromState   => deriveMomMSLiquid_fromState
          derVarPos(iField)%momentaFromState  => deriveMomentaMSLiquid_fromState
          derVarPos(iField)%velocitiesFromState  &
            & => deriveVelocitiesMSLiquid_fromState
          derVarPos(iField)%equilFromState => deriveEqMSLiquid_fromState
          derVarPos(iField)%equilFromAux => deriveEquilMSLiquid_fromAux
          select case (trim(schemeHeader%relaxation))
          case('bgk_withthermodynfac', 'mrt_withthermodynfac')
            derVarPos(iField)%auxFieldFromState => deriveAuxMSLiquid_fromState_WTDF
          case default
            derVarPos(iField)%auxFieldFromState => deriveAuxMSLiquid_fromState
          end select
        end do
      case ( 'isotherm_acEq' )
        call mus_append_derVar_isotherm_acEq( varSys       = varSys,         &
          &                                   solverData   = solverData,     &
          &                                   schemeHeader = schemeHeader,   &
          &                                   stencil      = stencil,        &
          &                                   fldLabel     = field(1)%label, &
          &                                   derVarname   = derVarname      )
        derVarPos(1)%equilFromMacro => deriveEquil_FromMacro_IsothermAcEq
        derVarPos(1)%velFromState   => deriveVelocity_FromState_IsothermAcEq
        derVarPos(1)%equilFromState => deriveEq_FromState_IsothermAcEq
        derVarPos(1)%auxFieldFromState => deriveAux_fromState
        derVarPos(1)%equilFromAux => deriveEquilIsoThermAcEq_fromAux

      case default
        write(logUnit(1),*) ' The selecited scheme kind is unknown '//        &
          &                 trim(schemeHeader%kind)
        call tem_abort()
      end select

      ! append boundary variables
      call mus_append_bcVar( varSys     = varSys,     &
        &                    solverData = solverData, &
        &                    derVarName = derVarName, &
        &                    nFields    = nFields,    &
        &                    field      = field,      &
        &                    stencil    = stencil     )

      ! append material variable
      call mus_append_materialVar( varSys       = varSys,       &
        &                          solverData   = solverData,   &
        &                          schemeHeader = schemeHeader, &
        &                          derVarName   = derVarName    )

      ! append physical variable
      ! get the list of physical variable labels from scheme specific
      ! mus_append_derVar_<schemeKind> routine
      call mus_append_dervar_physics( derVarname = derVarname,    &
        &                             varSys     = varSys,        &
        &                             solverData = solverData,    &
        &                             nFields    = nFields,       &
        &                             fldLabel   = field(:)%label )

      ! append extra variables defined in lua file
      solverData_evalElem%solver_bundle = c_loc(solverData)
      solverData_evalElem%stFun_setter => mus_set_stfun_getElement
      solverData_evalElem%opVar_setter => mus_set_opVar_getElement
      if (allocated(luaVar)) then
        call tem_varSys_append_luaVar(                     &
          & luaVar                   = luaVar,             &
          & varSys                   = varSys,             &
          & st_funList               = st_funList,         &
          & solverData_evalElem      = solverData_evalElem )
      end if

      ! append mesh info variables
      call tem_varSys_append_meshInfoVar( varSys = varSys )

      ! append field source variable to varSys and store position in varSys in
      ! mus_source_op_type
      do iField = 1, nFields
        write(logUnit(10),*) 'Append field source: iField ', iField
        call mus_append_sourceVar( me           = field(iField)%source, &
          &                        solverData   = solverData,           &
          &                        schemeHeader = schemeHeader,         &
          &                        varSys       = varSys,               &
          &                        nFields      = nFields,              &
          &                        stencil      = stencil,              &
          &                        poss_srcVar  = poss_srcVar,          &
          &                        fldLabel     = field(iField)%label   )

        write(logUnit(10),*) 'Append internal field source: iField ', iField
        call mus_append_sourceVar( me           = field(iField)%internalSource, &
          &                        solverData   = solverData,                   &
          &                        schemeHeader = schemeHeader,                 &
          &                        varSys       = varSys,                       &
          &                        nFields      = nFields,                      &
          &                        stencil      = stencil,                      &
          &                        poss_srcVar  = poss_srcVar,                  &
          &                        fldLabel     = field(iField)%label           )
      end do

      ! append global source variable to varSys and store position in varSys in
      ! mus_source_op_type
      write(logUnit(10),*) 'Append global source'
      call mus_append_sourceVar( me           = globSrc,      &
        &                        solverData   = solverData,   &
        &                        schemeHeader = schemeHeader, &
        &                        varSys       = varSys,       &
        &                        nFields      = nFields,      &
        &                        stencil      = stencil,      &
        &                        poss_srcVar  = poss_srcVar   )

    end do
    write(logUnit(1),*) 'Done appending variables to varSys'

    call tem_varSys_dump( varSys, dbgUnit(10) )

    ! store derVarPos
    call mus_store_derVarPos( derVarPos  = derVarPos,     &
      &                       derVarname = derVarname,    &
      &                       varSys     = varSys,        &
      &                       nFields    = nFields,       &
      &                       fldLabel   = field(:)%label )

    ! Store boundary variable position
    call mus_store_bcVarPos( field   = field,   &
      &                      nFields = nFields, &
      &                      varSys  = varSys   )

  end subroutine mus_build_varSys
  ! **************************************************************************** !

  ! **************************************************************************** !
  !> Append variable read from restart file as state variables
  subroutine mus_append_readVarAsStateVar( varSys, readVarIsPdf, read_varSys,  &
    &                                      stateVarMap, solverData, nFields,   &
    &                                      fldLabel )
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)      :: varSys

    !> Is true if read_varSys has pdf variable
    logical, intent(out)                      :: readVarIsPdf

    !> Variable system loaded from restart header file
    type(tem_varSys_type), intent(in)         :: read_varSys

    !> Store position of state variable in global varSys
    type(tem_varMap_type), intent(out)        :: stateVarMap

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> number of fields
    integer, intent(in)                       :: nFields

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel(:)
    ! ---------------------------------------------------------------------------
    character(len=labelLen) :: varname
    integer :: nComponents
    integer :: iVar, addedPos
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    ! ---------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    write(logUnit(1),*) 'Appending state variables read from restart file'

    ! Determine variables read from restart file are pdf variable
    ! or derived variable.
    ! Derived variable requires different treatment in mus_construction
    ! to create communication buffer and interpolation routines.
    !
    ! Input to mus_harvesting can have either only pdf or derived variable
    ! not combination of both
    call check_varSys_forPdfVar( readVarIsPdf = readVarIsPdf, &
      &                          varSys       = read_varSys,  &
      &                          nFields      = nFields,      &
      &                          fldLabel     = fldLabel      )

    ! get all state variable using access_state
    get_element => mus_access_state_ForElement
    get_point => mus_stateVar_forPoint
    setup_indices => mus_accessVar_setupIndices
    get_valOfIndex => mus_stateVar_fromIndex

    ! initialize state varMap
    call init(stateVarMap%varName)
    call init(stateVarMap%varPos)

    do iVar = 1, read_varSys%varname%nVals
      varname = trim(read_varSys%varname%val(iVar))
      nComponents = read_varSys%method%val(iVar)%nComponents
      call tem_varSys_append_stateVar(                          &
        &  me             = varSys,                             &
        &  varName        = varname,                            &
        &  nComponents    = nComponents,                        &
        &  method_data    = mus_get_new_solver_ptr(solverData), &
        &  get_point      = get_point,                          &
        &  get_element    = get_element,                        &
        &  set_params     = set_params,                         &
        &  get_params     = get_params,                         &
        &  setup_indices  = setup_indices,                      &
        &  get_valOfIndex = get_valOfIndex,                     &
        &  pos            = addedPos,                           &
        &  wasAdded       = wasAdded                            )

      if (wasAdded) then
        write(logUnit(10),*) 'Appended state variable: '//trim(varname)
        call append( me = stateVarMap%varPos, val = addedPos )
        call append( me  = stateVarMap%varname, val = varname )
      else
        write(logUnit(1),*) 'Error: State variable '//trim(varname)// &
          &                 ' is not added to variable system'
        call tem_abort()
      end if
    end do

    stateVarMap%nScalars = varSys%nScalars
    call truncate(me = stateVarMap%varPos)
    call truncate(me = stateVarMap%varname)

    ! debug output
    write(logUnit(10),*) '  nStateVars in varSys: ', varSys%nStateVars
    write(logUnit(10),*) '  nScalars in varSys: ', varSys%nScalars

  end subroutine mus_append_readVarAsStateVar
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> Append state variable depends on the scheme kind
  !!
  !! allocate function pointers, and append pdf to scheme%varSys
  subroutine mus_append_stateVar( varSys, stateVarMap, solverData,         &
    &                             schemeHeader, stencil, nFields, fldLabel )
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)      :: varSys

    !> Store position of state variable in global varSys
    type(tem_varMap_type), intent(out)        :: stateVarMap

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in)  :: schemeHeader

    !> compute stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil

    !> number of fields
    integer, intent(in)                       :: nFields

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel(:)
    ! ---------------------------------------------------------------------------
    integer :: iField, addedPos
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    character(len=labelLen) :: varname
    ! ---------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    write(logUnit(1),*) 'Appending state variables '

    ! check there no variables added to varSys before state variable
    if (varSys%varName%nVals /= 0) then
      write(logUnit(1),*) 'Error: Found variables before state variables'
      call tem_abort()
    end if

    ! get all state variable using access_state
    get_element => mus_access_state_ForElement
    get_point => mus_stateVar_forPoint
    setup_indices => mus_accessVar_setupIndices
    get_valOfIndex => mus_stateVar_fromIndex

    ! initialize state varMap
    call init(stateVarMap%varName)
    call init(stateVarMap%varPos)

    ! append pdf for each field
    do iField = 1, nFields
      write(varname,'(a)') trim( fldLabel( iField ) )//'pdf'
      call tem_varSys_append_stateVar(                          &
        &  me             = varSys,                             &
        &  varName        = varname,                            &
        &  nComponents    = stencil%QQ,                         &
        &  method_data    = mus_get_new_solver_ptr(solverData), &
        &  get_point      = get_point,                          &
        &  get_element    = get_element,                        &
        &  set_params     = set_params,                         &
        &  get_params     = get_params,                         &
        &  setup_indices  = setup_indices,                      &
        &  get_valOfIndex = get_valOfIndex,                     &
        &  pos            = addedPos,                           &
        &  wasAdded       = wasAdded                            )

      if (wasAdded) then
        write(logUnit(10),*) 'Appended state variable: '//trim(varname)
        call append( me = stateVarMap%varPos, val = addedPos )
        call append( me = stateVarMap%varname, val = varname )
      else
        write(logUnit(1),*) 'Error: State variable '//trim(varname)// &
          &                 ' is not added to variable system'
        call tem_abort()
      end if
    end do

    stateVarMap%nScalars = varSys%nScalars
    call truncate(me = stateVarMap%varPos)
    call truncate(me = stateVarMap%varname)

    ! debug output
    write(logUnit(10),"(A,I0)") '  nStateVars in varSys: ', varSys%nStateVars
    write(logUnit(10),"(A,I0)") '  nScalars   in varSys: ', varSys%nScalars

  end subroutine mus_append_stateVar
  ! **************************************************************************** !


  ! *************************************************************************** !
  !> Append auxiliary variables which are computed from state and stored
  !! in auxField array using calcAuxField function
  subroutine mus_append_auxField(varSys, solverData, schemeHeader, nFields, &
    &                            fldLabel, derVarname)
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)      :: varSys
    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData
    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in)  :: schemeHeader
    !> number of fields
    integer, intent(in)                       :: nFields
    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel(:)
    !> array of derive physical variables
    type(grw_labelarray_type), intent(inout) :: derVarName
    ! ---------------------------------------------------------------------------
    integer :: iVar, iField, addedPos, nComponents, nDerVars
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    character(len=labelLen) :: varname
    character(len=labelLen), allocatable :: derVarName_loc(:)
    ! ---------------------------------------------------------------------------
    ! Initialize variables
    nDerVars = 0

    write(logUnit(1),*) 'Appending auxiliary variables '
    select case (trim(schemeHeader%kind))
    case ('fluid', 'fluid_incompressible')
      ! auxiliary variable for incompressible is same as compressible model
      ! append density and velocity as auxField variables
      nDerVars = 2
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'density ', 'velocity'  ]
    case ('fluid_GNS', 'fluid_incompressible_GNS') 
      ! auxiliary variable for simulating generalized Navier-Stokes
      ! equation for flows through porous media
      ! append density, velocity and fluid volume fraction as auxField variables
      nDerVars = 3
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'density ', 'velocity', 'vol_frac' ]
    case ('passive_scalar')
      ! append density as auxField variable
      nDerVars = 1
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'density' ]
    case ('poisson', 'poisson_boltzmann_linear', &
      &   'poisson_boltzmann_nonlinear')
      ! append potential as auxField variable
      nDerVars = 1
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'potential' ]
    case ('nernst_planck')
      ! append mole density as auxField variable
      nDerVars = 1
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'mole_density' ]
    case ('multispecies_gas', 'multispecies_liquid')
      ! append density and velocity of species as auxField variable
      ! mixture density and velocity are appended later are derive variable
      nDerVars = 2
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'density ', 'momentum' ]
    case ('isotherm_acEq')
      ! auxiliary variable for isotherm_acEq is same as incompressible model
      ! append density and velocity as auxField variables
      nDerVars = 2
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'density ', 'velocity' ]
    case default
      write(logUnit(1),*) ' The selected scheme kind is unknown '//        &
        &                 trim(schemeHeader%kind)
      call tem_abort()
    end select

    ! get all auxField variable uses same access routines. auxField_varPos
    ! is used to access exact variable
    get_element => mus_access_auxFieldVar_forElement
    get_point => mus_auxFieldVar_forPoint
    setup_indices => mus_accessVar_setupIndices
    get_valOfIndex => mus_auxFieldVar_fromIndex
    set_params => tem_varSys_setParams_dummy
    get_params => tem_varSys_getParams_dummy

    ! append dervarnames to growing array to create physics variable
    do iVar = 1, nDerVars
      call append(derVarName, derVarName_loc(iVar))
    end do

    do iField = 1, nFields
      do iVar = 1, nDerVars
        select case(trim(adjustl(derVarName_loc(iVar))))
        case ('density', 'mole_density', 'potential', 'vol_frac')
          nComponents = 1
        case ('velocity', 'momentum')
          nComponents = 3
        !case ('grad_velocity')
        !  nComponents = 9
        case default
          write(logUnit(1),*) 'WARNING: Unknown variable: '//&
            &                 trim(derVarName_loc(iVar))
          cycle !go to next variable
        end select

        write(varname,'(a)') trim(fldLabel(iField)) &
          &               //trim(adjustl(derVarName_loc(ivar)))
        call tem_varSys_append_auxFieldVar(                       &
          &  me             = varSys,                             &
          &  varName        = varname,                            &
          &  nComponents    = nComponents,                        &
          &  method_data    = mus_get_new_solver_ptr(solverData), &
          &  get_point      = get_point,                          &
          &  get_element    = get_element,                        &
          &  set_params     = set_params,                         &
          &  get_params     = get_params,                         &
          &  setup_indices  = setup_indices,                      &
          &  get_valOfIndex = get_valOfIndex,                     &
          &  pos            = addedPos,                           &
          &  wasAdded       = wasAdded                            )

        if (wasAdded) then
          write(logUnit(10),*) 'Appended auxField variable: '//trim(varname)
        else if (addedPos < 1) then
          write(logUnit(1),*) 'Error: AuxField variable '//trim(varname)// &
            &                 ' is not added to variable system'
          call tem_abort()
        end if
      end do !iVar
    end do !iField

   ! debug output
    write(logUnit(10),"(A,I0)") '  nAuxVars in varSys: ', varSys%nAuxVars
    write(logUnit(10),"(A,I0)") '  nAuxScalars in varSys: ', varSys%nAuxScalars

  end subroutine mus_append_auxField
  ! *************************************************************************** !

  ! *************************************************************************** !
  !> Build a variable system of all possible source terms for the given
  !! schemeKind
  subroutine mus_append_sourceVar( me, solverData, schemeHeader, varSys,   &
    &                              nFields, stencil, poss_srcVar, fldLabel )
    ! --------------------------------------------------------------------------
    !> Contains source function pointer,
    !! source variable definition from lua and
    !! mapping of source variable in global varSys
    type(mus_source_type), intent(inout)          :: me

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in)     :: solverData

    !> Identifier of the scheme
    type(mus_scheme_header_type), intent(in)      :: schemeHeader

    !> Global variable system
    type(tem_varSys_type), intent(inout)          :: varSys

    !> number of fields
    integer, intent(in)                           :: nFields

    !> compute stencil defintion
    type(tem_stencilHeader_type), intent(in)      :: stencil

    !> possible source variables
    type(tem_possible_variable_type), intent(in)  :: poss_srcVar

    !> array of field label prefix required only for field source.
    !! If not present, it is assumed as global source
    character(len=*), optional, intent(in)        :: fldLabel
    ! --------------------------------------------------------------------------
    logical :: wasAdded
    character(len=labelLen), allocatable ::  input_varname(:)
    character(len=labelLen)  ::  varName, fldLabel_loc
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    integer :: iSrc, dataVar_InInVar, nSrcVars, nComponents
    integer :: addedPos, iField
    integer :: nComp_defined, nComp_expected, data_varPos, possSrc_varPos
    ! --------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    get_element => tem_varSys_getElement_dummy
    get_point => mus_deriveVar_ForPoint
    setup_indices => mus_opVar_setupIndices
    get_valOfIndex => tem_varSys_getValOfIndex_dummy

    nSrcVars = me%varDict%nVals

    ! do nothing if there no source variable defined in lua file
    if (nSrcVars == 0) return

    ! if present then it is field source else global source.
    ! If field source then inputs for source variable are
    ! field pdf and source spacetime function.
    ! Else global source then inputs for source variable are
    ! all fields pdf and source spacetime function.
    if (present(fldLabel)) then
      fldLabel_loc = trim(fldLabel)
      allocate(input_varname(2))
      input_varname(1) = trim(fldLabel)//'pdf'
      dataVar_InInVar = 2
      ! use stencil%QQ since source variable return value which is
      ! added to state as source term so both state and source should
      ! have same nComponents
      nComponents = stencil%QQ
    else
      fldLabel_loc = ''
      allocate(input_varname(nFields+1))
      do iField = 1, nFields
        input_varname(iField) = trim(varSys%varname%val(iField))
      end do
      dataVar_InInVar = nFields + 1
      ! mixture source: return source value for all species
      nComponents = stencil%QQ * nFields
    end if

    srcLoop: do iSrc = 1, nSrcVars
      input_varname(dataVar_InInVar) = trim(me%varDict%val(iSrc)%value)
      ! get actual source variable name from temSource varname
      ! which stores the name as defined in possible sources
      varname = trim(me%varDict%val(iSrc)%key)
      me%method(iSrc)%varname = trim(varname)

      ! Check number of components expected by possible variable is
      ! same as number of components defined by space time function variable
      data_varPos = PositionOfVal( me  = varSys%varname,                      &
        &                          val = trim(input_varname(dataVar_InInVar)) )
      if (data_varPos > 0) then
        ! position of variable name in possible variable list
        possSrc_varPos = PositionOfVal( me  = poss_srcVar%varname, &
          &                             val = trim(varname)        )

        nComp_defined = varSys%method%val(data_varPos)%nComponents
        nComp_expected = poss_srcVar%nComponents%val(possSrc_varPos)

        if ( nComp_defined /= nComp_expected ) then
          write(logUnit(1),'(a)') 'Error: Appending source variable'
          write(logUnit(1),'(a,i0)') 'nComponent of defined variable: "'//     &
            & trim(input_varname(dataVar_InInVar))//'"= ', nComp_defined
          write(logUnit(1),'(a,i0)') '/= nComponent of expected variable: "'// &
            & trim(varname)//'"= ', nComp_expected
          call tem_abort()
        end if
      else
        ! user st_fun variable is not found in varSys
        write(logUnit(1),*) 'Error: User defined space-time function variable'
        write(logUnit(1),*) '"'//trim(input_varName(dataVar_InInVar))// &
          &                 '" not found in varSys'
        call tem_abort()
      end if

      ! set to default to dummy routine which does noting to auxField
      me%method(iSrc)%applySrc => mus_applySrc_dummy
      me%method(iSrc)%addSrcToAuxField => mus_addSrcToAuxField_dummy
      me%method(iSrc)%updateSourceVar => mus_updateSrcVar_dummy

      ! choose appropriate function pointer
      ! select get_element, applySrc and addSrcToAuxField according to scheme
      ! kind
      select case (trim(schemeHeader%kind))
      case ('fluid', 'fluid_incompressible')
        select case (trim(varname))
        case ('force')
          ! select pointer according to order
          if (me%method(iSrc)%order == 2) then
            select case (trim(schemeHeader%relaxation))
            case ('mrt')
              ! get_element and applySrc are same of fluid and
              ! fluid_incompressible
              get_element => derive_force_MRT
              select case (trim(schemeHeader%layout))
              case ('d2q9')
                me%method(iSrc)%applySrc => applySrc_force_MRT_d2q9
              case ('d3q19')
                me%method(iSrc)%applySrc => applySrc_force_MRT_d3q19
              case ('d3q27')
                me%method(iSrc)%applySrc => applySrc_force_MRT_d3q27
              case default
                me%method(iSrc)%applySrc => applySrc_force_MRT
              end select
            case default
              me%method(iSrc)%applySrc => applySrc_force
            end select

            ! select addSrcToAuxField according to scheme kind
            select case (trim(schemeHeader%kind))
            case ('fluid')
              me%method(iSrc)%addSrcToAuxField => mus_addForceToAuxField_fluid
            case ('fluid_incompressible')
              me%method(iSrc)%addSrcToAuxField                             &
                &                      => mus_addForceToAuxField_fluidIncomp
            end select
          else ! 1st order
            ! addSrcToAuxField is not required for force1stOrd 1st order
            get_element => derive_force1stOrd
            me%method(iSrc)%applySrc => applySrc_force1stOrd
          end if

        case ('turb_channel_force_accel')
          ! Select apply force according to scheme kind.
          select case (trim(schemeHeader%kind))
          case ('fluid')
            select case (trim(schemeHeader%relaxation))
            case ('mrt')
              select case (trim(schemeHeader%layout))
              case ('d2q9')
                me%method(iSrc)%applySrc => applySrc_turbChanForce_MRT_d2q9
              case ('d3q19')
                me%method(iSrc)%applySrc => applySrc_turbChanForce_MRT_d3q19
              case ('d3q27')
                me%method(iSrc)%applySrc => applySrc_turbChanForce_MRT_d3q27
              case default
                me%method(iSrc)%applySrc => applySrc_turbChanForce_MRT
              end select
            case default
              me%method(iSrc)%applySrc => applySrc_turbChanForce
            end select
          case ('fluid_incompressible')
            call tem_abort('TurbChannel only implemented for fluid scheme')
          end select

          ! Add force acceleration to velocity field.
          ! It is same for both fluid and fluid_incompressible.
          me%method(iSrc)%addSrcToAuxField                      &
            &             => mus_addTurbChanForceToAuxField_fluid
          me%method(iSrc)%updateSourceVar => mus_updateSrcVar_turbChanForce

        case ('absorb_layer', 'absorb_layer_inlet', 'absorb_layer_outlet')
          ! Absorb layer is dependent of collision but independent on
          ! scheme kind. So implemented seperate routine for MRT.
          ! select addSrcToAuxField according to scheme relaxation.
          ! Use time average quantities if pressure or velocity is defined
          ! as dynamic.
          if (me%method(iSrc)%absLayer%config%isPressDyn &
            & .or. me%method(iSrc)%absLayer%config%isVelDyn) then
            ! \todo KM: 20210301 Implement seperate routine for
            ! absorb_layer_inlet and absorb_layer_outlet when target_velocity
            ! and target_pressure respectively are defined as stFun.
            me%method(iSrc)%addSrcToAuxField                   &
              &             => mus_addDynSponFldToAuxField_fluid

            me%method(iSrc)%updateSourceVar => mus_updateSrcVar_dynSponFld

            select case (trim(schemeHeader%relaxation))
            case ('mrt')
              ! KM: \todo 25012021 Implement optimized routine for d3q19
              me%method(iSrc)%applySrc => applySrc_absorbLayerDyn_MRT
            case default
              me%method(iSrc)%applySrc => applySrc_absorbLayerDyn
            end select
          else
            me%method(iSrc)%addSrcToAuxField => mus_addSponFldToAuxField_fluid
            select case (trim(schemeHeader%relaxation))
            case ('mrt')
              ! KM: \todo 25012021 Implement optimized routine for d3q19
              me%method(iSrc)%applySrc => applySrc_absorbLayer_MRT
            case default
              me%method(iSrc)%applySrc => applySrc_absorbLayer
            end select
          end if

        case ('hrr_correction')
          if (trim(schemeHeader%kind) == 'fluid') then
            ! get_element and applySrc are same of fluid and fluid_incompressible
            ! select pointer according to order
            select case (trim(schemeHeader%relaxation))
            case ('hrr_bgk_corrected', 'prr_bgk_corrected', 'rr_bgk_corrected')
              me%method(iSrc)%applySrc => mus_applySrc_dummy
            case default
              call tem_abort('HRR Correction not supported for '  &
              &            //trim(schemeHeader%relaxation)        )
            end select

            ! select addSrcToAuxField  and get_element according to scheme kind
            select case (trim(schemeHeader%layout))
            case ('d2q9')
              me%method(iSrc)%addSrcToAuxField => mus_addHRRCorrToAuxField_fluid_2D
              get_element => derive_HRRCorrection_d2q9
            case ('d3q19')
              me%method(iSrc)%addSrcToAuxField => mus_addHRRCorrToAuxField_fluid_3D
              get_element => derive_HRRCorrection_d3q19
            case ('d3q27')
              me%method(iSrc)%addSrcToAuxField => mus_addHRRCorrToAuxField_fluid_3D
              get_element => derive_HRRCorrection_d3q27
            case default
              call tem_abort('HRR Correction not supported for '  &
              &            //trim(schemeHeader%layout)        )
            end select
          else
            call tem_abort('HRR Correction not supported for '  &
            &            //trim(schemeHeader%kind)              )
          end if

        case default
          call tem_abort('Unknown source variable for ' &
            &            //trim(schemeHeader%kind)      )
        end select

      case('fluid_GNS', 'fluid_incompressible_GNS')
        select case (trim(varname))
        case ('force')
          ! select pointer according to order
          if (me%method(iSrc)%order == 2) then
            select case (trim(schemeHeader%relaxation))
            case ('bgk')
              if (trim(schemeHeader%layout) == 'd3q19' &
                & .OR. trim(schemeHeader%layout) == 'd2q9') then
                me%method(iSrc)%applySrc => applySrc_force_GNS
              else
                write(logUnit(1),*) 'Error loading body force routine ', &
                & 'layout ', trim(schemeHeader%layout),' not supported for fluid_GNS scheme'
                call tem_abort() 
              end if
            case default
              write(logUnit(1),*) 'Error loading body force routine ', &
              & 'relaxation ', trim(schemeHeader%relaxation),' not supported for fluid_GNS scheme'
              call tem_abort() 
            end select ! relaxation kind

            ! Assign addSrcToAuxField which is used to incorporate the force term
            ! in the macroscopic velocity
            me%method(iSrc)%addSrcToAuxField => mus_addForceToAuxField_fluid_GNS
          else ! 1st order
            ! addSrcToAuxField is not required for force1stOrd 1st order
            write(logUnit(1),*) 'Error loading body force routine ', &
            & 'first order forcing terms are not supported for fluid_GNS!'
            call tem_abort() 
          end if ! force order == 2
        case default
          call tem_abort('Unknown source variable for ' &
            &            //trim(schemeHeader%kind)      )
        end select ! varname

      case ('multispecies_liquid')
        select case (trim(varname))
        case ('force')
          ! select pointer according to order
          if (me%method(iSrc)%order == 2) then
            select case (trim(schemeHeader%relaxation))
            case ('bgk_withthermodynfac', 'mrt_withthermodynfac')
              me%method(iSrc)%addSrcToAuxField                          &
                &                      => mus_addForceToAuxField_MSL_WTDF
              me%method(iSrc)%applySrc => applySrc_forceMSLiquid_2ndOrd_WTDF
            case default
              me%method(iSrc)%addSrcToAuxField => mus_addForceToAuxField_MSL
              me%method(iSrc)%applySrc => applySrc_forceMSLiquid_2ndOrd
            end select
          else ! 1st order
            select case (trim(schemeHeader%relaxation))
            case ('bgk_withthermodynfac', 'mrt_withthermodynfac')
              me%method(iSrc)%applySrc => applySrc_forceMSLiquid_1stOrd_WTDF
            case default
              me%method(iSrc)%applySrc => applySrc_forceMSLiquid_1stOrd
            end select
          end if

        case ('electric_field')
          ! select pointer according to order
          if (me%method(iSrc)%order == 2) then
            select case (trim(schemeHeader%relaxation))
            case ('bgk_withthermodynfac', 'mrt_withthermodynfac')
              me%method(iSrc)%addSrcToAuxField                             &
                &                      => mus_addElectricToAuxField_MSL_WTDF
              me%method(iSrc)%applySrc => applySrc_electricMSLiquid_2ndOrd_WTDF
            case default
              me%method(iSrc)%addSrcToAuxField => mus_addElectricToAuxField_MSL
              me%method(iSrc)%applySrc => applySrc_electricMSLiquid_2ndOrd
            end select
          else ! 1st order
            select case (trim(schemeHeader%relaxation))
            case ('bgk_withthermodynfac', 'mrt_withthermodynfac')
              me%method(iSrc)%applySrc => applySrc_electricMSLiquid_1stOrd_WTDF
            case default
              me%method(iSrc)%applySrc => applySrc_electricMSLiquid_1stOrd
            end select
          end if
        case default
          call tem_abort('Unknown source variable for ' &
            &            //trim(schemeHeader%kind)      )
        end select

      case ('nernst_planck')
        select case (trim(varname))
        case ('electric_field')
          me%method(iSrc)%applySrc => applySrc_electricFieldNP
        case default
          call tem_abort('Unknown source variable for ' &
            &            //trim(schemeHeader%kind)      )
        end select

      case ('passive_scalar')
        select case (trim(varname))
        case ('injection')
          if ( solverData%scheme%transVar%varDict%nVals > 0) then
            if ( trim(solverData%scheme%transVar%varDict%val(1)%key) &
              & /= 'transport_velocity') then
              write(logUnit(1),*) 'Error: transport_velocity variable required ' &
                &              // 'for injection source is not defined'
              call tem_abort()
            end if
          else
            write(logUnit(1),*) 'Error: injection source requires '&
              &               //'transport_velocity variable'
            call tem_abort()
          end if
          get_element => derive_injectionPS
          me%method(iSrc)%applySrc => applySrc_injectionPS
        case ('equal_injection')
          get_element => derive_equalInjectionPS
          me%method(iSrc)%applySrc => applySrc_equalInjectionPS
        case default
          call tem_abort('Unknown source variable for ' &
            &            //trim(schemeHeader%kind)      )
        end select

      case ('poisson')
        select case (trim(varname))
        case ('charge_density')
          ! select pointer according to order
          if (me%method(iSrc)%order == 2) then
            me%method(iSrc)%addSrcToAuxField => mus_addSrcToAuxField_poisson
            get_element => deriveSrc_chargeDensity
            me%method(iSrc)%applySrc => applySrc_chargeDensity_2ndOrd
          else ! 1st order
            get_element => deriveSrc_chargeDensity
            me%method(iSrc)%applySrc => applySrc_chargeDensity_1stOrd
          end if
        case default
          call tem_abort( 'Unknown source variable for ' &
            &             //trim(schemeHeader%kind)      )
        end select
      case default
        write(logUnit(1),*) 'ERROR: Scheme kind: '//trim(schemeHeader%kind) &
          &                 // 'does not support source variable: '         &
          &                 // trim(varname)
        call tem_abort()
      end select

      varname = trim(fldLabel)//'src_'//trim(varname)

      ! append variable to varSys
      call tem_varSys_append_derVar(                              &
        &    me             = varSys,                             &
        &    varName        = trim(varname),                      &
        &    nComponents    = nComponents,                        &
        &    input_varname  = input_varname,                      &
        &    method_data    = mus_get_new_solver_ptr(solverData), &
        &    get_point      = get_point,                          &
        &    get_element    = get_element,                        &
        &    set_params     = set_params,                         &
        &    get_params     = get_params,                         &
        &    setup_indices  = setup_indices,                      &
        &    get_valOfIndex = get_valOfIndex,                     &
        &    pos            = addedPos,                           &
        &    wasAdded       = wasAdded                            )

      if (wasAdded) then
        write(logUnit(10),*) ' Appended variable:'//trim(varname)
        me%method(iSrc)%srcTerm_varPos = addedPos
        me%method(iSrc)%data_varPos = varSys%method%val(addedPos) &
          &                                 %input_varPos(dataVar_InInVar)
      else if (addedpos < 1) then
        write(logUnit(1),*) 'Error: variable '//trim(varname)// &
          &                 ' is not added to variable system'
        call tem_abort()
      end if
    end do srcLoop

    ! debug output
    call tem_varSys_dump( varSys, dbgUnit(10), me%method(:)%srcTerm_varPos )

  end subroutine mus_append_sourceVar
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> Store the position of each variable in the global system in the derVarPos
  !! This function is also called in Harvester.
  subroutine mus_store_derVarPos( derVarPos, derVarname, varSys, nFields,      &
    &                             fldLabel )
    ! ---------------------------------------------------------------------------
    !> Position of derived variables
    type(mus_derVarPos_type), intent(inout) :: derVarPos(:)

    !> array of derive physical variables
    type(grw_labelarray_type), intent(in) :: derVarName

    !> global variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> number of fields
    integer, intent(in)                       :: nFields

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in) :: fldLabel(:)
    ! ------------------------------------------------------------------------
    integer :: iField, iVar, nFields_loc, varPos, nDerVars
    character(len=labelLen)  ::  varName
    ! ------------------------------------------------------------------------
    nDerVars = derVarname%nVals

    if (nFields > 1) then
      nFields_loc = nFields + 1 ! nSpecies + 1 mixture
    else
      nFields_loc = 1
    end if

    write(logUnit(1),*) 'Storing state variable position '
    do iField = 1, nFields
      ! store pdf pos
      varname = trim(fldLabel(iField))//'pdf'
      varPos = PositionOfVal( me  = varSys%varname, val = trim(varname) )
      derVarPos(iField)%pdf = varPos

      ! store omega pos
      varname = trim(fldLabel(iField))//'omega'
      varPos = PositionOfVal( me  = varSys%varname, val = trim(varname) )
      derVarPos(iField)%omega = varPos
    end do

    write(logUnit(1),*) 'Storing derive variable position '
    do iField = 1, nFields_loc

      do iVar = 1, nDerVars
        if (iField > nFields) then ! mixture
          varname = trim(adjustl(derVarName%val(iVar)))
        else
          varname = trim(fldLabel(iField))//trim(adjustl(derVarName%val(iVar)))
        end if

        varPos = PositionOfVal( me  = varSys%varname, val = trim(varname) )

        select case(trim(adjustl(derVarName%val(iVar))))
        case ('pdf')
          derVarPos(iField)%pdf = varPos
        case ('fetch_pdf')
          derVarPos(iField)%fetch_pdf = varPos
        case ('density')
          derVarPos(iField)%density = varPos
        case ('mole_density')
          derVarPos(iField)%moleDensity = varPos
        case ('pressure')
          derVarPos(iField)%pressure = varPos
        case ('kinematic_pressure')
          derVarPos(iField)%kinePress = varPos
        case ('velocity')
          derVarPos(iField)%velocity = varPos
        !case ('grad_velocity')
        !  derVarPos(iField)%grad_velocity = varPos
        case ('vel_mag')
          derVarPos(iField)%velMag = varPos
        case ('momentum')
          derVarPos(iField)%momentum = varPos
        case ('shear_stress')
          derVarPos(iField)%shearStress = varPos
        case ('wss' )
          derVarPos(iField)%wss = varPos
        case ('shear_mag')
          derVarPos(iField)%shearMag = varPos
        case ('strain_rate')
          derVarPos(iField)%strainRate = varPos
        case ('shear_rate')
          derVarPos(iField)%shearRate = varPos
        case ('kinetic_energy')
          derVarPos(iField)%kineticEnergy = varPos
        case ('temperature')
          derVarPos(iField)%temperature = varPos
        case ('mole_fraction')
          derVarPos(iField)%moleFrac = varPos
        case ('mass_fraction')
          derVarPos(iField)%massFrac = varPos
        case ('mole_flux')
          derVarPos(iField)%moleFlux = varPos
        case ('equilibrium')
          derVarPos(iField)%equilibrium = varPos
        case ('non_equilibrium')
          derVarPos(iField)%nonEquilibrium = varPos
        case ('equilibrium_vel')
          derVarPos(iField)%equilibriumVel = varPos
        case ('potential')
          derVarPos(iField)%potential = varPos
        case ('vol_frac')
          derVarPos(iField)%vol_frac = varPos
        case default
          write(logUnit(10),*) 'WARNING: Unknown variable: '//trim(varname)
          !write(*,*) 'derVarName ', derVarName%val(iVar)
          cycle !go to next variable
        end select
      enddo
    end do

  end subroutine mus_store_derVarPos
! ****************************************************************************** !


  ! *************************************************************************** !
  !> Store the position of each boundary variable in the global varSys
  !! in the field%bc%varPos%<variable>.
  !! This routine also checks if boundary variable defined in config file
  !! has same number of components as expected.
  subroutine mus_store_bcVarPos( field, nFields, varSys )
    ! ---------------------------------------------------------------------------
    !> Field containing boundary infos
    type(mus_field_type), intent(inout) :: field(:)

    !> number of fields
    integer, intent(in)                 :: nFields

    !> global variable system
    type(tem_varSys_type), intent(in)   :: varSys
    ! ------------------------------------------------------------------------
    integer :: iField, iBC, iVar, defVar_pos
    character(len=labelLen) :: def_varName, bc_varName
    integer :: nComp_defined, nComp_expected
    ! ------------------------------------------------------------------------
    write(logUnit(10),*) 'Storing boundary variable position in varSys'

    do iField = 1, nFields
      do iBC = 1, size(field(iField)%bc)
        do iVar = 1, field(iField)%bc(iBC)%varDict%nVals
          ! Variable name loaded from boundary table
          def_varName = trim(field(iField)%bc(iBC)%varDict%val(iVar)%value)
          ! position of boundary variable in varSys
          defVar_pos = PositionOfVal( me  = varSys%varName, &
            &                       val = trim(def_varName) )
          ! continue only if this variable exist in varSys
          if (defVar_pos>0) then
            nComp_defined = varSys%method%val(defVar_pos)%nComponents

            ! check number of components defined for the variable is
            ! same as expected, if so store the position if boundary type
            bc_varName = trim(field(iField)%bc(iBC)%varDict%val(iVar)%key)
            select case(trim(bc_varName))
            case ('velocity')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %velocity%nComponents
              field(iField)%bc(iBC)%bc_states%velocity%varPos = defVar_pos
            case ('pdf')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %pdf%nComponents
              field(iField)%bc(iBC)%bc_states%pdf%varPos = defVar_pos
            case ('pressure')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %pressure%nComponents
              field(iField)%bc(iBC)%bc_states%pressure%varPos = defVar_pos
            case ('mass_flowrate')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %massFlowRate%nComponents
              field(iField)%bc(iBC)%bc_states%massFlowRate%varPos = defVar_pos
            case ('mole_fraction')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %moleFrac%nComponents
              field(iField)%bc(iBC)%bc_states%moleFrac%varPos = defVar_pos
            case ('mole_density')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %moleDens%nComponents
              field(iField)%bc(iBC)%bc_states%moleDens%varPos = defVar_pos
            case ('mole_flux')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %moleFlux%nComponents
              field(iField)%bc(iBC)%bc_states%moleFlux%varPos = defVar_pos
            case ('mole_diff_flux')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %moleDiff_flux%nComponents
              field(iField)%bc(iBC)%bc_states%moleDiff_flux%varPos = defVar_pos
            case ('potential')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %potential%nComponents
              field(iField)%bc(iBC)%bc_states%potential%varPos = defVar_pos
            case ('surface_charge_density')
              nComp_expected = field(iField)%bc(iBC)%bc_states  &
                &                           %surChargeDens%nComponents
              field(iField)%bc(iBC)%bc_states%surChargeDens%varPos = defVar_pos
            case default
              write(logUnit(1),*) 'Error: Unknown boundary variable: "'// &
                & trim(bc_varName)//'"'
              call tem_abort()
            end select

            if (nComp_defined /= nComp_expected) then
              write(logUnit(1),*) 'Error: Storing boundary variable position'
              write(logUnit(1),'(a,i0)') 'nComponent of defined variable: "'// &
                & trim(def_varName)//'"= ', nComp_defined
              write(logUnit(1),'(a,i0)') '/= nComponent of expected '// &
                & 'variable: "'//trim(bc_varName)//'"= ', nComp_expected
              call tem_abort()
            end if

          else
            write(logUnit(1),*) 'Error: User defined space-time function '
            write(logUnit(1),*) 'or reference to varName for boundary variable'
            write(logUnit(1),*) '"'//trim(def_varName)//'" not found in varSys'
            call tem_abort()
          end if !defVar_pos

        end do !iVar
      end do !iiBC
    end do !iField

  end subroutine mus_store_bcVarPos
  ! **************************************************************************** !


  ! **************************************************************************** !
  !> This function runs over all variable loaded from restart file check if
  !! variables loaded are pdf variable or derive variable
  !!
  !! Variable read from restart file can have variables other than pdf but pdf
  !! must be the 1st variable if not stateVarIsPdf will be false
  subroutine check_varSys_forPdfVar( readVarIsPdf, varSys, nFields, fldLabel )
    ! ---------------------------------------------------------------------------
    !> return true if variable read from restart file has pdf
    logical, intent(out) :: readVarIsPdf

    !> variable system loaded from restart file
    type(tem_varSys_type), intent(in) :: varSys

    !> number of fields
    integer, intent(in)               :: nFields

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)      :: fldLabel(:)
    ! ---------------------------------------------------------------------------
    integer :: iField
    character(len=labelLen) :: buffer
    ! ---------------------------------------------------------------------------

    readVarIsPdf = .true.
    do iField = 1, nFields
      write(buffer,'(a)') trim( fldLabel( iField ) )//'pdf'
      if ( trim(buffer) /= varSys%varName%val(iField) ) then
         readVarIsPdf = .false.
         write(logUnit(1),*) 'WARNING: Variable loaded from restart file does'
         write(logUnit(1),*) '         not have pdf as 1st variables.'
         write(logUnit(1),*) '         Following are deactivated:'
         write(logUnit(1),*) '           * Connectivity neigh array'
         write(logUnit(1),*) '           * Interpolation of pdf for multilevel'
         exit
      end if
    end do

  end subroutine check_varSys_forPdfVar
! ****************************************************************************** !


end module mus_variable_module
! ****************************************************************************** !
