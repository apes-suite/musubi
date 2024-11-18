! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2013, 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
! Copyright (c) 2022-2023 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
! ******************************************************************************!
!> author: Simon Zimny
!! author: Kannan Masilamani
!! Definition of the datatypes for the scheme implementation.
!!
!! For a detailed description see
!! [Scheme Implementation](../page/features/scheme.html).
!!
module mus_scheme_module

  ! include treelm modules
  use env_module,            only: rk, newUnit, LabelLen, PathLen
  use tem_aux_module,        only: tem_checkLabel, tem_abort
  use treelmesh_module,      only: treelmesh_type
  ! use tem_bc_prop_module,    only: tem_bc_prop_type
  use tem_tools_module,      only: tem_horizontalSpacer
  use tem_tracking_module,   only: tem_load_tracking, tem_tracking_config_type
  use tem_variable_module,   only: tem_variable_load
  use tem_stencil_module,    only: tem_stencil_map_ToTreelmDef,                &
    &                              tem_treelmDef_map_toStencil,                &
    &                              init, append, destroy
  use tem_logging_module,    only: logUnit
  use tem_dyn_array_module,  only: append, destroy
  use tem_varSys_module,     only: tem_varSys_init
  use tem_debug_module,      only: dbgUnit
  use tem_spacetime_fun_module, only: tem_destroy_subTree_of_st_funList

  ! include musubi modules
  use mus_param_module,         only: mus_param_type
  use mus_scheme_layout_module, only: mus_load_newLayout, &
    &                                 mus_define_layout,  &
    &                                 mus_init_layout
  use mus_scheme_header_module, only: mus_load_scheme_header, &
    &                                 mus_scheme_header_out
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_moments_module,       only: mus_init_moments
  use mus_field_module,         only: mus_load_fields, mus_fields_out, &
    &                                 mus_load_fieldBaseInfos, mus_field_cleanup
  use mus_mixture_module,       only: mus_mixture_out
  use mus_variable_module,      only: mus_build_varSys,          &
    &                                 mus_append_stateVar,       &
    &                                 mus_append_readVarAsStateVar
  use mus_source_type_module,   only: mus_source_cleanup,  &
    &                                 mus_load_source_var, &
    &                                 mus_create_poss_srcVar
  use mus_varSys_module,        only: mus_varSys_solverData_type, &
    &                                 mus_init_varSys_solverData
  use mus_interpolate_header_module, only: mus_load_interpolate,               &
    &                                      mus_set_nSources,                   &
    &                                      mus_interpolate_out
  use mus_geom_module,           only: mus_geom_type
  use mus_transport_var_module,  only: mus_create_poss_transVar, &
    &                                 mus_load_transport_var
  use mus_auxFieldVar_module,    only: mus_assign_calcAuxField_ptr
  use mus_gradData_module,       only: mus_assign_GradCalculation_ptr
  use mus_scheme_derived_quantities_module, only: mus_assign_derived_functions_ptr

  ! include aotus modules
  use aotus_module,     only: flu_state, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close,                 &
    &                         aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_open_table,                &
    &                         aot_out_close_table, aot_out_open, aot_out_close

  implicit none

  private

  public :: mus_load_scheme
  public :: mus_scheme_out
  public :: mus_init_scheme
  public :: mus_scheme_cleanup

contains

! ****************************************************************************** !
  !> load single scheme defined in lua file with or without scheme handle
  !!
  !! This routines checks whether schemes table
  !!```lua
  !! boundary_conditions = {...},
  !!```
  !! is defined. If yes then it will
  !! load schemes table and set mus_scheme_type.
  !! If no special scheme table is defined, default tables are loaded from the
  !! root level of the lua file
  !! fluid, boundary conditions and initial condiitions
  !!```lua
  !! fluid = { ... }
  !! boundary_conditions = {...}
  !! initial_conditions = {...}
  !!```
  !!
  !!
  subroutine mus_load_scheme( me, solverData, geometry, conf, params, &
    &                         parent, isMusHvs                        )
    ! ---------------------------------------------------------------------------
    !> scheme type
    type( mus_scheme_type ), target, intent(inout) :: me
    !> contains pointer to scheme, physics types
    type( mus_varSys_solverData_type ), target :: solverData
    !> geometry information like tree and boundary
    type( mus_geom_type ), intent(in), target    :: geometry
    !> global parameter type
    type( mus_param_type ), target, intent(inout) :: params
    !> flu state
    type( flu_State ) :: conf
    !> parent handle if scheme table is defined
    integer, intent(in), optional :: parent
    !> Logic to not to load tracking and variable table if this routine
    !! is called from mus_hvs_config_load.
    !! Default is False
    logical, optional, intent(in) :: isMusHvs
    ! ---------------------------------------------------------------------------
    integer, allocatable :: varErr(:)
    logical :: isMusHvs_loc
    ! ---------------------------------------------------------------------------
    if (present(isMusHvs)) then
      isMusHvs_loc = isMusHvs
    else
      isMusHvs_loc = .false.
    end if

    ! Set solverData pointer for variable method data
    call mus_init_varSys_solverData( me        = solverData,     &
      &                              scheme    = me,             &
      &                              physics   = params%physics, &
      &                              geometry  = geometry        )

    ! load scheme identifier table
    call mus_load_scheme_header( me      = me%header,     &
      &                          conf    = conf,          &
      &                          parent  = parent,        &
      &                          scaling = params%scaling )

    if (geometry%tree%global%minLevel /= geometry%tree%global%maxLevel) then
      ! Ensure to have the proper scaling available for multi-level meshes!
      if ('acoustic' == trim(params%scaling)) then
        select case(trim(me%header%kind))
        case ('fluid', 'fluid_incompressible', 'isotherm_acEq')
          write(logUnit(1),*) 'Using ' // trim(params%scaling) &
            &                 // ' scaling in multi-level mesh!'
        case default
            write(logUnit(1),*) 'Scaling is set to acoustic'
            write(logUnit(1),*) 'But ' // trim(me%header%kind) // ' requires' &
              &                 // ' diffusive scaling.'
            write(logUnit(1),*) 'Sorry, can not proceed on a multi-level mesh!'
            call tem_abort('ERROR: diffusive scaling not available for'    &
              &            // ' multi-level with ' // trim(me%header%kind) )
        end select
      else
        write(logUnit(1),*) 'Multi-level with diffusive scaling not' &
          &                 // ' implemented!'
        write(logUnit(1),*) 'Sorry, can not proceed on a multi-level mesh!'
        call tem_abort('ERROR: diffusive scaling not available for' &
          &            // ' multi-level'                            )
      end if
    end if

    ! initialize scheme layout here to append state variable to varSys
    ! call mus_init_layout( layout   = me%layout )

    ! create fStencil
    if (trim(me%header%layout) == 'new_stencil') then
      call tem_horizontalSpacer(fUnit = logUnit(1))
      write(logUnit(1),*) 'Reading the new layout...'
      call mus_load_newLayout( me            = me%layout, &
        &                      parent_handle = parent,    &
        &                      conf          = conf       )
    else
      ! create fStencil
      call mus_define_layout( layout      = me%layout,           &
        &                     stencilName = me%header%layout,    &
        &                     nElems      = geometry%tree%nElems )
    end if

    ! Allocate field array and load field lables.
    ! it must be loaded before to append state variables
    call mus_load_fieldBaseInfos( me      = me%field,   &
      &                           nFields = me%nFields, &
      &                           parent  = parent,     &
      &                           conf    = conf        )

    ! Variable system must be initialized here so annouymous source and
    ! boundary variables can be appended to variable system during
    ! loading
    write(logUnit(1),*) 'Initializing global variable system for scheme '&
      &                 //'kind: '//trim(me%header%kind)
    call tem_varSys_init( me         = me%varSys,            &
      &                   systemName = trim(me%header%kind), &
      &                   length     = 8                     )

    ! if load scheme is called from harvesting then append variable loaded
    ! for restart file as state variables
    if (isMusHvs_loc) then
      ! Append variables loaded from restart header as state variables
      ! with get_element pointing to access_state
      call mus_append_readVarAsStateVar(          &
        & varSys       = me%varSys,               &
        & readVarIsPdf = me%readVarIsPdf,         &
        & read_varSys  = params%general%restart   &
        &                %header%varSys,          &
        & stateVarMap  = me%stateVarMap,          &
        & solverData   = solverData,              &
        & nFields      = me%nFields,              &
        & fldLabel     = me%field(:)%label        )
    else
      me%readVarIsPdf = .true.

      ! append state variable depends on scheme kind
      call mus_append_stateVar( varSys       = me%varSys,          &
        &                       stateVarMap  = me%stateVarMap,     &
        &                       solverData   = solverData,         &
        &                       schemeHeader = me%header,          &
        &                       stencil      = me%layout%fstencil, &
        &                       nFields      = me%nFields,         &
        &                       fldLabel     = me%field(:)%label   )
    end if

    ! load interpolation parameters
    call mus_load_interpolate( me     = me%intp%config, &
      &                        parent = parent,         &
      &                        conf   = conf            )

    ! create possible source variables depends on scheme kind
    call mus_create_poss_srcVar( poss_srcVar  = me%poss_srcVar, &
      &                          schemeHeader = me%header       )

    ! create possible transport variables depends on scheme kind
    call mus_create_poss_transVar( poss_transVar = me%poss_transVar, &
      &                            schemeHeader  = me%header         )

    ! load fields from parent, including possible sources for each field
    call mus_load_fields( me           = me%field,                      &
      &                   varSys       = me%varSys,                     &
      &                   nFields      = me%nFields,                    &
      &                   mixture      = me%mixture,                    &
      &                   nernstPlanck = me%nernstPlanck,               &
      &                   bc_prop      = geometry%boundary,             &
      &                   conf         = conf,                          &
      &                   minLevel     = geometry%tree%global%minLevel, &
      &                   maxLevel     = geometry%tree%global%maxLevel, &
      &                   parent       = parent,                        &
      &                   schemeHeader = me%header,                     &
      &                   poss_srcVar  = me%poss_srcVar,                &
      &                   physics      = params%physics,                &
      &                   scaling      = params%scaling,                &
      &                   layout       = me%layout,                     &
      &                   isMusHvs     = isMusHvs                       )


    ! load tracking and variable from conf only if it is not called from
    ! mus_hvs_config_load
    if (.not. isMusHvs_loc) then

      ! load source variables
      call mus_load_source_var( me       = me%globSrc,         &
        &                       possVars = me%poss_srcVar,     &
        &                       conf     = conf,               &
        &                       parent   = parent,             &
        &                       key      = 'glob_source',      &
        &                       varSys   = me%varSys           )

      ! load transport variables
      call mus_load_transport_var( me           = me%transVar,      &
        &                          possVars     = me%poss_transVar, &
        &                          conf         = conf,             &
        &                          parent       = parent,           &
        &                          varSys       = me%varSys,        &
        &                          schemeHeader = me%header         )

    end if

    ! add analytical function loaded from lua file
    ! and in the end add error variable without dependency
    call tem_variable_load( me     = me%luaVar,  &
      &                     conf   = conf,       &
      &                     parent = parent,     &
      &                     key    = 'variable', &
      &                     vError = varErr      )

    ! load tracking
    call tem_load_tracking( me     = me%track, &
      &                     conf   = conf,     &
      &                     parent = parent    )

    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_load_scheme
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initialize single scheme stencil and variable system
  subroutine mus_init_scheme( me, tree, solverData )
    ! ---------------------------------------------------------------------------
    !> single scheme to initialize
    type(mus_scheme_type), intent(inout) :: me
    !> global treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> contains pointer to scheme, physics types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------
    write(logUnit(1),*) 'Initializing the scheme ...'

    write(logUnit(1),*) 'Map stencil to treelm definition and vice versa ...'
    ! map the stencil offsets to the treelm definitions
    ! required to map boundary ID direction defined by treeLM directions
    ! from fluid stencil directions
    call tem_stencil_map_toTreelmDef( me%layout%fStencil )
    ! map the stencil offsets to the treelm definitions
    call tem_treelmDef_map_toStencil( me%layout%fStencil )

    ! append fluid stencil as 1st stencil in growing array of stencil
    call append( me  = me%layout%grwStencil,  &
      &          val = me%layout%fStencil     )

    ! append fluid stencil label
    call append( me  = me%layout%stencil_labels, &
      &          val = me%layout%fStencil%label, &
      &          pos = me%layout%fStencil_pos    )

    ! initialize the moments matrix
    call mus_init_moments( me           = me%layout%moment,         &
      &                    QQ           = me%layout%fStencil%QQ,    &
      &                    cxDir        = me%layout%fStencil%cxDir, &
      &                    label        = me%layout%fStencil%label, &
      &                    schemeHeader = me%header                 )

    if ( tree%global%minLevel /= tree%global%maxLevel ) then
      ! set number of sources required by interpolation
      call mus_set_nSources( me     = me%intp,                  &
        &                    nDims  = me%layout%fStencil%nDims, &
        &                    QQ     = me%layout%fStencil%QQ,    &
        &                    layout = me%header%layout          )
    end if

    ! build variable system only if read var from restart is PDF
    if (me%readVarIsPDF) then
      call mus_build_varSys( varSys       = me%varSys,          &
        &                    solverData   = solverData,         &
        &                    schemeHeader = me%header,          &
        &                    stencil      = me%layout%fStencil, &
        &                    nFields      = me%nFields,         &
        &                    derVarPos    = me%derVarPos,       &
        &                    luaVar       = me%luaVar,          &
        &                    field        = me%field(:),        &
        &                    globSrc      = me%globSrc,         &
        &                    poss_srcVar  = me%poss_srcVar,     &
        &                    st_funList   = me%st_funList       )
    end if

    ! assign function pointer to compute auxFieldVar
    call mus_assign_calcAuxField_ptr(me%header, me%calcAuxField)

    ! Initialize gradients pointers
    me%Grad = mus_assign_GradCalculation_ptr(label = me%layout%fStencil%label)

    ! Initialize quantities type in the layout class
    me%layout%quantities = mus_assign_derived_functions_ptr( &
      & label_stencil = me%layout%fStencil%label, &
      & label_fluid = me%header%kind              )

    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_init_scheme
! ****************************************************************************** !


! ****************************************************************************** !
  !> Dump single scheme info into restart solver specific conf to dump
  !! solver specific information in restart header file
  subroutine mus_scheme_out( me, conf )
    ! ---------------------------------------------------------------------------
    !> schemes to dump to restart header file
    type(mus_scheme_type), intent(in) :: me
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), optional, intent(inout) :: conf
    ! ---------------------------------------------------------------------------
    ! Dump scheme identify table
    call mus_scheme_header_out( me%header, conf )

    ! \todo KM Dump: 1. stencil
    !                2. Interpolation
    select case ( trim( me%header%layout ) )
    case ( 'new_stencil' )
      write(logUnit(1),*) 'WARNING: New stencil is not dumped in solver '//&
        &                 'specific unit. So cannot use mus_harvesting'
    end select

    ! Dump interpolation
    call mus_interpolate_out( me   = me%intp, &
      &                         conf = conf     )

    ! Dump mixture info for multi-scheme
    call mus_mixture_out( me           = me%mixture, &
      &                   conf         = conf,       &
      &                   schemeHeader = me%header   )

    ! Dump field info
    call mus_fields_out(me%field, conf, me%header)

  end subroutine mus_scheme_out
! ****************************************************************************** !

! ****************************************************************************** !
  !> This subroutine acts as a destructor for the construct routine
  subroutine mus_scheme_cleanup( me, minLevel, maxLevel, nBCs )
    ! ---------------------------------------------------------------------------
    !> scheme information including fluid, boundary and flow information
    type( mus_scheme_type ), intent(inout) :: me
    !> minlevel
    integer, intent(in) :: minLevel
    !> maxlevel
    integer, intent(in) :: maxLevel
    !> Number of boundary conditions
    integer, intent(in) :: nBCs
    ! ---------------------------------------------------------------------------

    write(dbgUnit(1),*) "Enter mus_scheme_cleanup"

    ! ... deallocate the globBC array
    if( allocated( me%globBC ))then
      deallocate( me%globBC )
    end if

    deallocate(me%pdf)
    deallocate(me%state)
    if ( me%track%control%active ) then
      deallocate ( me%track%instance )
      allocate( me%track%instance( me%track%control%nDefined ) )
    end if

    ! ... deallocate the level descriptor array
    if( allocated( me%levelDesc ))then
      deallocate( me%levelDesc )
    end if

    ! ... loop over all stencils and ...
    if ( allocated( me%layout%stencil ) ) then
      deallocate( me%layout%stencil )
    end if
    call destroy(me%layout%grwStencil)
    call destroy(me%layout%stencil_labels)

    ! ... deallcoate the auxField
    if( allocated(me%auxField) ) then
      deallocate( me%auxField )
    end if

    ! ... deallocate gradData
    if ( allocated(me%gradData) ) then
      deallocate(me%gradData)
    end if

    ! ... deallocate viscosity and omega array in fluid type
    ! and arrays allocated in field BC type
    ! KM: DO NOT DESTROY FIELD ARRAY AS IT CONTAINS ALL CONFIG INFO
    call mus_field_cleanup(me%field, me%header, minLevel, maxLevel, nBCs, &
      &                    me%nFields)

    ! ... destroy glob source info
    call mus_source_cleanup(me%globSrc)

    ! ... destroy space-time function in stFun_list
    call tem_destroy_subTree_of_st_funList( me = me%st_funList )

    write(dbgUnit(1),*) "Done!"

  end subroutine mus_scheme_cleanup
! ****************************************************************************** !

end module mus_scheme_module
! ****************************************************************************** !


! ****************************************************************************** !
!> \page scheme_implementation Scheme Implementation
!! The concept of schemes provides the user with a bigger flexibility. It is now
!! possible to run multiple simulations with different layouts on the same mesh.
!! This is the basis for LBM passive scalar transport simulations in MUSUBI.
!! In general the scheme implementation is able to run on old and new lua files.
!!
!! \subsection scheme_definition Scheme Definition
!!
!! The scheme is part of the individual solvers and contains all information about
!! initial conditions, boundary conditions, the \ref mus_scheme_layout_module
!! "Scheme Layout", fluid and flow properties.
!!
!! \subsection scheme_usage Usage
!!
!! To start a flow simulation one can either use the 'old' lua files or add the
!! table species with the following quantities:
!!
!! - the scheme identifier which include label, kind, relaxation, layout
!! - the initial conditions
!! - the boundary conditions
!! - the fluid quantities
!!
!!```lua
!! -- Sample Flow Scheme Definition
!!identify = {
!!   kind = 'fluid',
!!   relaxation = 'bgk',
!!   layout = 'd3q19',
!! } -- identify table
!!
!!-- Initial condition
!!initial_condition = {
!!  density = 1.0,
!!  velocityX = 0.0,
!!  velocityY = 0.0,
!!  velocityZ = 0.0,
!!  Sxx = 0.0, Syy = 0.0, Szz = 0.0, Sxy = 0.0, Syz = 0.0, Sxz = 0.0,
!!}
!!
!!-- Boundary conditions
!!boundary_condition = {
!!{ label = 'wall',
!!  kind = 'velocity_bounceback',
!!  velocity = {0.03, 0.0, 0.0 }
!!},
!!{ label = 'wall',
!!  kind = 'wall',
!!} -- boundary table
!!
!!fluid = {
!!   omega = 1.8,
!!   rho0 = 1.0
!!} -- fluid table
!!```
!!
!! Full example:
