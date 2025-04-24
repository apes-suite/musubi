! Copyright (c) 2012-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012, 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Philipp Otte <otte@mathcces.rwth-aachen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
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
! **************************************************************************** !
!> author: Kannan Masilamani
!! author: Simon Zimny
!! This module contains information about all fields like fluid,
!! species, temperature etc. This field type will be used for
!! multispecies and passive scalar transport.
!!
!! [mus_field_prop]: @ref mus_field_prop_module "mus_field_prop_module"
!!
module mus_field_module

  ! include treelm modules
  use env_module,               only: labelLen, rk
  use tem_param_module,         only: cs2inv, cs2, div1_2
  use tem_math_module,          only: invert_matrix
  use tem_time_module,          only: tem_time_type
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_debug_module,         only: dbgUnit
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_aux_module,           only: tem_abort, tem_CheckLabel
  use tem_ini_condition_module, only: tem_ini_condition_type, tem_load_ic
  use tem_varSys_module,        only: tem_varSys_type
  use tem_varMap_module,        only: tem_possible_variable_type
  use tem_restart_module,       only: tem_restart_type
  use tem_logging_module,       only: logUnit
  use tem_temporal_module,      only: tem_temporal_for
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_stringKeyValuePair_module, only: init, truncate

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_NonExistent
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_open_table,                &
    &                         aot_out_close_table, aot_out_val

  ! include musubi modules
  use mus_bc_header_module,     only: boundary_type, mus_load_bc,            &
    &                                 glob_boundary_type,                    &
    &                                 check_solid_in_bc, rearrange_bc_elems, &
    &                                 mus_fieldBC_cleanup
  use mus_scheme_header_module, only: mus_scheme_header_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_field_prop_module,    only: mus_field_prop_type, mus_field_prop_out, &
    &                                 mus_load_field_prop
  use mus_fluid_module,         only: mus_fluid_cleanup
  use mus_species_module,       only: compute_molWeightRatio,                  &
    &                                 compute_bulkViscOmega
  use mus_physics_module,       only: mus_physics_type
  use mus_mixture_module,       only: mus_mixture_type, mus_load_mixture
  use mus_source_type_module,   only: mus_source_type,     &
    &                                 mus_load_source_var, &
    &                                 mus_source_cleanup
  use mus_source_var_module,    only: mus_add_internal_source_var
  use mus_nernstPlanck_module,  only: mus_nernstPlanck_type, &
    &                                 mus_load_nernstPlanck

  implicit none

  private

  public :: mus_field_type
  public :: mus_load_fields, mus_fields_out
  public :: mus_load_fieldBaseInfos
  public :: setParameters_multispecies
  public :: remove_solid_in_bc
  public :: mus_check_allWall
  public :: mus_field_getSymmetricBCs
  public :: mus_field_cleanup

  !> This type contains all information on fields with ic and bc
  !! Example fields: fluid, species etc.
  !! Each field contains one initial condition and array of
  !! boundary conditions
  !!
  type mus_field_type
    !> field label. Should be unique for each field
    character(len=labelLen) :: label

    !> physics parameters (fluid and species) for field
    type(mus_field_prop_type) :: fieldProp

    !> array of field boundary types for each field
    !! size: #BCs in the boundary_condition table
    !! allocated in mus_load_bc
    type(boundary_type), allocatable :: bc(:)

    !> initialization case, one initial condition for each field
    type(tem_ini_condition_type) :: ic

    !> field source applied only to current field
    type(mus_source_type) :: source

    !> field internal source applied only to current field
    !! not needed input from musubi.lua
    !! used only for HRR correction at the moment
    type(mus_source_type) :: internalSource

    !> An instance of restart type
    type(tem_restart_type) :: restart
  end type mus_field_type

  !> Interface for dumping a single field or a set of fields in a file in lua
  !! format.
  !!
  interface mus_fields_out
    module procedure mus_fields_out_vec
    module procedure mus_field_out_scal
  end interface mus_fields_out


contains


! **************************************************************************** !
  !> Subroutine to load the field table from the lua configuration file.
  !!
  !! If field table is not defined than load bc, ic, fluid, species from scheme
  !! table. If scheme table is not defined than load field variables from
  !! config parent.
  !!
  subroutine mus_load_fields( me, varSys, nFields, mixture, nernstPlanck,  &
    &                         bc_prop, conf, parent, minLevel, maxLevel,   &
    &                         schemeHeader, poss_srcVar, physics, scaling, &
    &                         layout, isMusHvs                             )
    ! --------------------------------------------------------------------------
    !> array of field type
    type(mus_field_type), intent(inout) :: me(:)
    !> Global variable system required to append annoymous source and
    !! boundary variables
    type(tem_varSys_type), intent(inout) :: varSys
    integer, intent(in) :: nFields !< number of fields defined in lua file
    !> contains mixture information
    type(mus_mixture_type), intent(out) :: mixture
    !> contains solvent information
    type(mus_nernstPlanck_type), intent(out) :: nernstPlanck
    !> boundary data from mesh
    type(tem_bc_prop_type), intent(in):: bc_prop
    !> flu state
    type(flu_State), intent(inout) :: conf
    !> global pdf info
    integer, intent(in) :: minLevel, maxLevel
    !> identifier of the scheme
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    !> parent handle if scheme table is defined
    integer, intent(in), optional :: parent
    !> possible source variables
    type(tem_possible_variable_type), intent(in) :: poss_srcVar
    !> physics type to convert physics to lattice unit or vice versa
    type( mus_physics_type ), intent(in) :: physics
    !> scaling type
    character(len=labelLen), intent(in) :: scaling
    !> fluid stencil info
    type(mus_scheme_layout_type), intent(in) :: layout
    !> Logic to not to load tracking and variable table if this routine
    !! is called from mus_hvs_config_load.
    !! Default is False
    logical, optional, intent(in) :: isMusHvs
    ! --------------------------------------------------------------------------
    ! counter variables
    integer :: iField
    ! aotus handles
    integer :: field_handle, field_sub_handle
    ! --------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*) 'Loading the fields ...'
    write(logUnit(1),*) 'Number of fields to load: ', nFields
    call aot_table_open( L       = conf,         &
      &                  parent  = parent,       &
      &                  thandle = field_handle, &
      &                  key     = 'field'       )

    if ( field_handle == 0 ) then
      if (present(parent)) then
        ! field table is not defined load field variable from scheme
        write(logUnit(1),*) 'No field table defined.' &
          &              // 'Loading field data from scheme'

        call mus_load_field_single( me           = me(1),        &
          &                         varSys       = varSys,       &
          &                         nFields      = nFields,      &
          &                         bc_prop      = bc_prop,      &
          &                         conf         = conf,         &
          &                         minLevel     = minLevel,     &
          &                         maxLevel     = maxLevel,     &
          &                         parent       = parent,       &
          &                         poss_srcVar  = poss_srcVar,  &
          &                         physics      = physics,      &
          &                         schemeHeader = schemeHeader, &
          &                         scaling      = scaling,      &
          &                         layout       = layout,       &
          &                         isMusHvs     = isMusHvs      )

      ! scheme table is not defined load field variable from config
      else
        write(logUnit(1),*) 'No field table and scheme table are defined.' &
          &              // 'Loading field data from config'

        call mus_load_field_single( me           = me(1),        &
          &                         varSys       = varSys,       &
          &                         nFields      = nFields,      &
          &                         bc_prop      = bc_prop,      &
          &                         conf         = conf,         &
          &                         minLevel     = minLevel,     &
          &                         maxLevel     = maxLevel,     &
          &                         poss_srcVar  = poss_srcVar,  &
          &                         physics      = physics,      &
          &                         schemeHeader = schemeHeader, &
          &                         scaling      = scaling,      &
          &                         layout       = layout,       &
          &                         isMusHvs     = isMusHvs      )
      endif ! present parent?
    else ! check whether field is a single table or multiple table

      call aot_table_open( L       = conf,             &
        &                  parent  = field_handle,     &
        &                  thandle = field_sub_handle, &
        &                  pos     = 1                 )
      if ( field_sub_handle == 0 ) then
        ! field is a single table
        call aot_table_close( L=conf, thandle=field_sub_handle )
        write(logUnit(1),*) 'Field is a single table'

        call mus_load_field_single( me           = me(1),        &
          &                         varSys       = varSys,       &
          &                         nFields      = nFields,      &
          &                         bc_prop      = bc_prop,      &
          &                         conf         = conf,         &
          &                         parent       = field_handle, &
          &                         minLevel     = minLevel,     &
          &                         maxLevel     = maxLevel,     &
          &                         poss_srcVar  = poss_srcVar,  &
          &                         physics      = physics,      &
          &                         schemeHeader = schemeHeader, &
          &                         scaling      = scaling,      &
          &                         layout       = layout,       &
          &                         isMusHvs     = isMusHvs      )

      else ! field is multiple table
        call aot_table_close( L=conf, thandle=field_sub_handle )
        write(logUnit(1),*) 'Field is a multiple table'
        do iField = 1, nFields
          call aot_table_open( L       = conf,             &
            &                  parent  = field_handle,     &
            &                  thandle = field_sub_handle, &
            &                  pos     = iField            )

          write(logUnit(1),*) 'Properties for Field ', iField

          call tem_horizontalSpacer(fUnit = logUnit(1))
          call mus_load_field_single( me           = me(iField),       &
            &                         varSys       = varSys,           &
            &                         nFields      = nFields,          &
            &                         bc_prop      = bc_prop,          &
            &                         conf         = conf,             &
            &                         minLevel     = minLevel,         &
            &                         maxLevel     = maxLevel,         &
            &                         parent       = field_sub_handle, &
            &                         poss_srcVar  = poss_srcVar,      &
            &                         physics      = physics,          &
            &                         schemeHeader = schemeHeader,     &
            &                         scaling      = scaling,          &
            &                         layout       = layout,           &
            &                         isMusHvs     = isMusHvs          )

          call aot_table_close( L=conf, thandle=field_sub_handle )
        end do ! number of fields
      end if ! if field table is present
    end if

    call aot_table_close( L=conf, thandle=field_handle )

    !load scheme specific field
    select case(trim(schemeHeader%kind))
    case('fluid', 'fluid_incompressible','fluid_GNS', &
      & 'fluid_incompressible_GNS', 'passive_scalar', 'isotherm_acEq', &
      & 'poisson', 'poisson_boltzmann_linear', 'poisson_boltzmann_nonlinear' )
      !nFields must be = 1 scheme other than multispecies and nernst_planck
      if( nFields > 1 ) then
        write(logUnit(1),*) ' ERROR: Number of fields defined for '
        write(logUnit(1),*) '        '//trim(schemeHeader%kind)//' scheme: >1. '
        write(logUnit(1),*) '        Specify only one field. Aborting'
        call tem_abort()
      endif
    case('multispecies_gas','multispecies_liquid')
      !nFields must be > 1 for multispecies
      if( nFields == 1 ) then
        write(logUnit(1),*) ' ERROR: Number of fields(species) defined for '
        write(logUnit(1),*) '        multispecies scheme = 1. '
        write(logUnit(1),*) '        Specify nFields > 1. Aborting'
        call tem_abort()
      endif

      !compute molecular weight ratios incase of multispecies simulation
      call compute_molWeightRatio(                              &
        &  molWeights    = me(:)%fieldProp%species%molWeight,   &
        &  molWeigRatios = me(:)%fieldProp%species%molWeigRatio )

      !load mixture table and set bulk viscosity relaxation time
      ! for each species which is depends on molecular weight ratios
      call mus_load_mixture( me           = mixture,      &
        &                    conf         = conf,         &
        &                    parent       = parent,       &
        &                    schemeHeader = schemeHeader, &
        &                    minLevel     = minLevel,     &
        &                    maxLevel     = maxLevel,     &
        &                    physics      = physics,      &
        &                    nFields      = nFields       )

      !compute bulk viscosity omega for each field
      write(logUnit(1),*)'  Bulk omega for each species:'
      do iField = 1,nFields
        write(logUnit(1),*) '   species ', iField
        call compute_bulkViscOmega( species = me(iField)%fieldProp%species,    &
          &                         bulkvisc = mixture%bulk_viscosityLB,       &
          &                         bulkviscLvl = mixture%relaxLvl(:)%bulkVisc,&
          &                         minLevel = minLevel, &
          &                         maxLevel = maxLevel )
      end do
    case('nernst_planck')
      write(logUnit(1),"(A)") ' Loading properties for nersnt_planck:'
      call mus_load_nernstPlanck(  me         = nernstPlanck, &
        &                          conf       = conf,         &
        &                          parent     = parent,       &
        &                          physics    = physics       )
    case default
      write(logUnit(1),*) 'The selected scheme kind is unknown '//           &
        &             trim(schemeHeader%kind)
      call tem_abort()
    end select

  end subroutine mus_load_fields
! **************************************************************************** !


! **************************************************************************** !
  !> load a single field table
  !! In includes:
  !!   load field property
  !!   load source variables
  !!   load boundary defination
  !!   load immersed boundary method
  !!   load initial condition defination and its property
  !!
  subroutine mus_load_field_single( me, varSys, nFields, bc_prop, conf,      &
    &                               parent, minLevel, maxLevel, poss_srcVar, &
    &                               physics, schemeHeader, scaling, layout,  &
    &                               isMusHvs                                 )
    ! --------------------------------------------------------------------------
    !> field type
    type(mus_field_type), intent(inout) :: me
    !> Global variable system required to append annoymous source and
    !! boundary variables
    type(tem_varSys_type), intent(inout) :: varSys
    !> number of fields defined in lua file
    integer, intent(in) :: nFields
    !> boundary data from mesh
    type(tem_bc_prop_type), intent(in):: bc_prop
    !> flu state
    type(flu_State), intent(inout) :: conf
    !> parent handle if scheme table is defined
    integer, intent(in), optional :: parent
    !> global pdf info
    integer, intent(in) :: minLevel, maxLevel
    !> possible source variables
    type(tem_possible_variable_type), intent(in) :: poss_srcVar
    !> physics type to convert physics to lattice unit or vice versa
    type(mus_physics_type), intent(in) :: physics
    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    !> scaling type
    character(len=labelLen), intent(in) :: scaling
    !> fluid stencil info
    type(mus_scheme_layout_type), intent(in) :: layout
    !> Logic to not to load tracking and variable table if this routine
    !! is called from mus_hvs_config_load.
    !! Default is False
    logical, optional, intent(in) :: isMusHvs
    ! --------------------------------------------------------------------------
    character(len=labelLen), allocatable :: ic_states(:)
    logical :: isMusHvs_loc
    integer :: IC_nVars
    !> Error code of loading ic variables
    integer, allocatable :: iError(:)
    ! --------------------------------------------------------------------------
    if (present(isMusHvs)) then
      isMusHvs_loc = isMusHvs
    else
      isMusHvs_loc = .false.
    end if

    write(logUnit(1),"(A)") ' Loading field: '//trim(me%label)
    ! load field properties
    call mus_load_field_prop( me           = me%fieldProp, &
      &                       nFields      = nFields,      &
      &                       conf         = conf,         &
      &                       minLevel     = minLevel,     &
      &                       parent       = parent,       &
      &                       physics      = physics,      &
      &                       cs_lattice   = layout%cs,    &
      &                       schemeHeader = schemeHeader  )

    ! Load Boundary, initial condition only for solver
    if (.not. isMusHvs_loc) then

      ! Read initial condition
      call mus_set_ic_states( schemeHeader%kind, ic_states, IC_nVars )
      allocate( iError( size(ic_states) ) )

      call tem_load_ic( me        = me%ic,               &
        &               conf      = conf,                &
        &               parent    = parent,              &
        &               key       = 'initial_condition', &
        &               ErrCode   = iError,              &
        &               StateName = ic_states            )

      if ( any( iError(1:IC_nVars) /= 0 ) ) then
        write(logUnit(1), "(A)") " Not all the KEY variables in initial " &
          &                      //"condition are well defined!"
        call tem_abort()
      end if

      deallocate( ic_states )
      deallocate( iError    )
    end if

    ! load boundary
    call mus_load_bc( me      = me%bc,          &
      &               bc_prop = bc_prop,        &
      &               conf    = conf,           &
      &               parent  = parent,         &
      &               varSys  = varSys,         &
      &               stencil = layout%fStencil )

    ! load the field specific source variables
    call mus_load_source_var( me       = me%source,       &
      &                       possVars = poss_srcVar,     &
      &                       conf     = conf,            &
      &                       parent   = parent,          &
      &                       key      = 'source',        &
      &                       varSys   = varSys           )

    call mus_add_internal_source_var( me           = me%internalSource, &
      &                               possVars     = poss_srcVar,       &
      &                               varSys       = varSys,            &
      &                               schemeHeader = schemeHeader       )

    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_load_field_single
! **************************************************************************** !


  ! ************************************************************************ !
  !> This routine returns nFields and field labels from config file.
  !! It is required to initialize variable system.
  !! labels are loaded only if field table is present else default
  !! is set to empty string.
  subroutine mus_load_fieldBaseInfos( me, nFields, parent, conf )
    ! --------------------------------------------------------------------------
    !> array of field type
    type( mus_field_type ), allocatable, intent(out) :: me(:)
    !> number of fields defined in lua file
    integer, intent(out) :: nFields
    !> parent handle if scheme table is defined
    integer, intent(in), optional :: parent
    !> flu state
    type(flu_State), intent(inout) :: conf
    ! --------------------------------------------------------------------------
    integer :: iError, iField
    integer :: field_handle, field_sub_handle
    ! --------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*) 'Loading the field base info ...'

    call aot_table_open( L       = conf,         &
      &                  parent  = parent,       &
      &                  thandle = field_handle, &
      &                  key     = 'field'       )

    ! allocate the fields
    nFields = 1
    allocate(me(nFields))
    me(nFields)%label = ''

    ! load label only if field table is present
    if ( field_handle /= 0 ) then
      ! check whether field is a single table or multiple table
      call aot_table_open( L       = conf,             &
        &                  parent  = field_handle,     &
        &                  thandle = field_sub_handle, &
        &                  pos     = 1                 )

      if ( field_sub_handle == 0 ) then
        ! field is a single table
        call aot_get_val( L       = conf,         &
           &              thandle = field_handle, &
           &              key     = 'label',      &
           &              val     = me(1)%label,  &
           &              default = '',           &
           &              ErrCode = iError        )

        if( trim(me(1)%label) /= '')             &
          & me(1)%label = trim( me(1)%label )//'_'
        write(logUnit(1),*) 'Field label: ', trim(me(1)%label)

      else ! field is multiple table
        call aot_table_close( L=conf, thandle=field_sub_handle )
        nFields = aot_table_length( L=conf, thandle=field_handle )
        write(logUnit(1),*) 'Number of fields: ', nFields

        ! reallocate the fields
        deallocate( me )
        allocate( me( nFields ))

        do iField = 1, nFields
          call aot_table_open( L       = conf,             &
            &                  parent  = field_handle,     &
            &                  thandle = field_sub_handle, &
            &                  pos     = iField            )

          call aot_get_val( L       = conf,             &
            &               thandle = field_sub_handle, &
            &               key     = 'label',          &
            &               val     = me(iField)%label, &
            &               default = '',               &
            &               ErrCode = iError            )

          if (btest(iError, aoterr_NonExistent)) then
            write(logUnit(1),*) 'Error: field label is not specified'
            write(logUnit(1),*) 'field label is neccessary when nFields>1'
            call tem_abort()
          end if

          ! add underscore here to differentiate state variables and derived
          ! variables of each field
          if( trim(me(iField)%label) /= '' )                          &
            &          me(iField)%label = trim( me(iField)%label )//'_'

          ! Check, if the names of the fields are unique. If not, add an
          ! increasing number
          call tem_checkLabel( label = me(:)%label, nLabels = iField )

          call aot_table_close( L=conf, thandle=field_sub_handle )
          write(logUnit(1),*) iField, 'Field label: ', trim(me(iField)%label)
        end do !number of fields

      end if ! single field
    else
      write(logUnit(1),*) 'No field table defined.'
      write(logUnit(1),*) 'Assuming single field and label = "".'
    end if
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_load_fieldBaseInfos
  ! ************************************************************************ !


! **************************************************************************** !
  !> write array of fields into a lua file
  !!
  subroutine mus_fields_out_vec( me, conf, schemeHeader )
    ! --------------------------------------------------------------------------
    !> array of field type
    type( mus_field_type ), intent(in) :: me(:)
    !> aotus out type
    type( aot_out_type ), intent(inout) :: conf
    !> identifier of the scheme
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------
    integer :: iField ! counter variable
    ! --------------------------------------------------------------------------

    call aot_out_val( put_conf = conf,                                         &
      &               val      = size(me),                                     &
      &               vname    = 'nFields' )

    if (size(me) > 1) then
      call aot_out_open_table(put_conf = conf, tname = 'field')
      do iField = 1, size(me)
        call mus_field_out_scal( me           = me( iField ), &
          &                      conf         = conf,         &
          &                      schemeHeader = schemeHeader, &
          &                      level        = 1             )
      end do
      call aot_out_close_table(put_conf = conf)
    else
      call mus_field_out_scal( me           = me( 1 ),      &
        &                      conf         = conf,         &
        &                      schemeHeader = schemeHeader, &
        &                      level        = 0             )
    end if
  end subroutine mus_fields_out_vec
! **************************************************************************** !


! **************************************************************************** !
  !> write single field into a lua file
  !!
  subroutine mus_field_out_scal( me, conf, schemeHeader, level )
    ! --------------------------------------------------------------------------
    !> single field type
    type( mus_field_type ), intent(in) :: me
    !> aotus out type
    type( aot_out_type ), intent(inout) :: conf
    !> identifier of the scheme
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    !> To dump field with or without key
    integer, optional, intent(in) :: level
    ! --------------------------------------------------------------------------
    integer :: level_loc, pos
    ! --------------------------------------------------------------------------
    if (present(level)) then
      level_loc = level
    else
      level_loc = 0
    end if

    if( level_loc == 0) then
      call aot_out_open_table( put_conf = conf, tname = 'field' )
    else
      call aot_out_open_table( put_conf = conf )
    end if

    if (len(trim(me%label)) > 0) then
      ! Remove underscore from field name before dumping
      pos = INDEX(trim(me%label),'_', back=.true.)
      call aot_out_val( put_conf = conf,                   &
        &               vname    = 'label',                &
        &               val      = trim(me%label(1:pos-1)) )
    end if

    call mus_field_prop_out( me           = me%fieldProp, &
      &                      conf         = conf,         &
      &                      schemeHeader = schemeHeader  )
    call aot_out_close_table(put_conf = conf)

  end subroutine mus_field_out_scal
! **************************************************************************** !


! **************************************************************************** !
  !> Set ic states labels by scheme kind
  !!
  subroutine mus_set_ic_states( scheme_kind, ic_states, IC_nVars )
    ! --------------------------------------------------------------------------
    character(len=labelLen),  intent(in) :: scheme_kind
    character(len=labelLen), allocatable :: ic_states(:)
    !> Number of initial condition variables required to initialize state
    integer, intent(out) :: IC_nVars
    ! --------------------------------------------------------------------------

    ! load scheme specific field
    select case( trim(scheme_kind) )
    case('fluid', 'fluid_incompressible', 'fluid_GNS', 'fluid_incompressible_GNS')
      IC_nVars = 4
      allocate(ic_states(10))
      ic_states = (/ 'pressure ', 'velocityX', 'velocityY', 'velocityZ', &
        &            'Sxx      ', 'Syy      ', 'Szz      ', 'Sxy      ', &
        &            'Syz      ', 'Sxz      ' /)
    case('passive_scalar')
      IC_nVars = 4
      allocate(ic_states(4))
      ic_states = (/ 'pressure ', 'velocityX', 'velocityY', 'velocityZ' /)
    case('multispecies_gas')
      IC_nVars = 4
      allocate( ic_states(4) )
      ic_states = (/ 'pressure ', 'velocityX', 'velocityY', &
        &            'velocityZ' /)
    case('multispecies_liquid', 'nernst_planck')
      IC_nVars = 4
      allocate( ic_states(4) )
      ic_states = (/ 'mole_fraction', 'velocityX    ', 'velocityY    ', &
        &            'velocityZ    ' /)
    case('isotherm_acEq')
      IC_nVars = 4
      allocate(ic_states(10))
      ic_states = (/ 'pressure ', 'velocityX', 'velocityY', 'velocityZ', &
        &            'Sxx      ', 'Syy      ', 'Szz      ', 'Sxy      ', &
        &            'Syz      ', 'Sxz      ' /)
    case('poisson', 'poisson_boltzmann_linear', 'poisson_boltzmann_nonlinear')
      IC_nVars = 1
      allocate(ic_states(1))
      ic_states = (/ 'potential' /)
    case default
      write(logUnit(1),*) 'The selected scheme kind model is not '&
         &              //'supported to initialize ic state names: '&
         &                 //trim( scheme_kind )
      call tem_abort()
    end select

  end subroutine mus_set_ic_states
! **************************************************************************** !

! **************************************************************************** !
  !> Set parameters for multispecies
  subroutine setParameters_multispecies( field, nFields, mixture, header, &
    &                                    layout, iLevel, tNow )
    ! --------------------------------------------------------------------------
    integer,                        intent(in) :: nFields
    type( mus_field_type ),      intent(inout) :: field(nFields)
    type( mus_scheme_layout_type ), intent(in) :: layout
    type( mus_scheme_header_type ), intent(in) :: header
    type( mus_mixture_type ),    intent(inout) :: mixture
    integer,                        intent(in) :: iLevel
    !> solver general info
    type( tem_time_type ),       intent(in) :: tNow
    ! --------------------------------------------------------------------------
    real(kind=rk) :: omega_diff ! diffusive relaxation parameter
    real(kind=rk) :: omega_kine ! kinematic shear relaxation parameter
    real(kind=rk) :: omega_bulk ! bulk relaxation parameter
    real(kind=rk) :: fac ! ramping multiplication factor
    integer :: iField, iDir
    real(kind=rk), dimension( layout%fStencil%QQ, layout%fStencil%QQ ) :: &
      &                                       identity, tmpMatrix, tmpMatrixInv
    ! --------------------------------------------------------------------------
!write(dbgUnit,*) 'setting omega paramter for multispecies'
    ! get the temporal multiplication factor
    fac = tem_temporal_for( temporal = mixture%omega_ramping, &
      &                     time     = tNow )
    !! Relaxation parameter for each level
    !! kine_viscosity = dt*c^2*(1/omega - 0.5)
    !!    =>   omega = dt / (viscosity/c^2 + dt/2)
    !! omega for each level is stored at fluid%omLvl
    omega_diff  =  fac * mixture%relaxLvl( iLevel )%paramB * cs2
!write(dbgUnit,*) 'omega_diff ', omega_diff
    mixture%relaxLvl( iLevel )%omega_diff = omega_diff


    ! if relaxation is mrt set kinematic and bulk viscosity omega
    if( header%relaxation(1:3) == 'mrt' ) then
      do iField = 1, nFields
        if( .not. allocated( field( iField )%fieldProp%species          &
          &          %mrt( iLevel )%s_mrt))                             &
          & allocate( field( iField )%fieldProp%species                 &
          &          %mrt( iLevel )%s_mrt( layout%fStencil%QQ,          &
          &                                layout%fStencil%QQ))

        if( .not. allocated( field( iField )%fieldProp%species            &
          &                                 %mrt( iLevel )%omegaMoments)) &
          & allocate( field( iField )%fieldProp%species                   &
          &          %mrt( iLevel )%omegaMoments( layout%fStencil%QQ,     &
          &                                       layout%fStencil%QQ)     )

        if( .not. allocated( field( iField )%fieldProp%species             &
          &                                 %mrt( iLevel )%omegaMomForce)) &
          & allocate( field( iField )%fieldProp%species                    &
          &          %mrt( iLevel )%omegaMomForce( layout%fStencil%QQ,     &
          &                                        layout%fStencil%QQ)     )

        field( iField )%fieldProp%species%mrt( iLevel )       &
          &                              %omegaMoments = 0.0_rk
        field( iField )%fieldProp%species%mrt( iLevel )       &
          &                              %omegaMomForce = 0.0_rk

        tmpMatrix    = 0.0_rk
        tmpMatrixInv = 0.0_rk
        Identity     = 0.0_rk
        do iDir = 1, layout%fStencil%QQ
          Identity(iDir, iDir) = 1.0_rk
        end do

        ! viscous omega
        omega_kine = mixture%relaxLvl( iLevel )%omega_kine * fac
        ! bulk omega
        omega_bulk = field(iField)%fieldProp%species%omBulkLvl(iLevel)  &
          &        * fac
!write(*,*) 'omega_diff ', omega_diff
!write(*,*) 'omega_kine ', omega_kine
!write(*,*) 'omega_bulk ', omega_bulk
        ! Set the relaxation parameters (only the non-relaxing modes)
        ! \ref Multiple-relaxation-time lattice Boltzmann scheme for
        !      homogeneous mixture flows with external force- Pietro Asinari
        ! \ref Jens ICCMES paper
        ! initialize all entries in relaxation matrix
        field( iField )%fieldProp%species%mrt( iLevel )%s_mrt(:,:) = 0.0_rk
        select case(trim(header%layout))
        case('d2q9' )
!write(*,*) 'omega_diaCenter ', (omega_bulk + omega_kine )/2.0_rk
!write(*,*) 'omega_nondia ', (omega_bulk - omega_kine )/2.0_rk
          ! diffusivity omegas
          do iDir=2,3
            field( iField )%fieldProp%species%mrt( iLevel )             &
              & %s_mrt( iDir, iDir ) = omega_diff
          end do

          do iDir=4,5
            field( iField )%fieldProp%species%mrt( iLevel )             &
              & %s_mrt( iDir, iDir ) = ( omega_bulk + omega_kine ) * div1_2
          end do

          field( iField )%fieldProp%species%mrt( iLevel )               &
              & %s_mrt( 5, 4 ) = ( omega_bulk - omega_kine ) * div1_2

          field( iField )%fieldProp%species%mrt( iLevel )               &
              & %s_mrt( 4, 5 ) = ( omega_bulk - omega_kine ) * div1_2

          field( iField )%fieldProp%species%mrt( iLevel )               &
            & %s_mrt( 6, 6 ) = omega_kine

          do iDir=7,9
            field( iField )%fieldProp%species%mrt( iLevel )             &
              & %s_mrt( iDir, iDir ) = mixture%omega_hom
          end do

        case('d3q19')
!write(*,*) 'omega_diaCenter ', (omega_bulk + 2.0_rk * omega_kine )/3.0_rk
!write(*,*) 'omega_nondia ', (omega_bulk - omega_kine )/3.0_rk
          do iDir=2,4
            field( iField )%fieldProp%species%mrt( iLevel )             &
              & %s_mrt( iDir, iDir ) = omega_diff
          end do

          do iDir=5,7
            field( iField )%fieldProp%species%mrt( iLevel )                 &
              & %s_mrt( iDir, iDir ) = ( omega_bulk + 2.0_rk * omega_kine ) &
              &                        / 3.0_rk
          end do

          field( iField )%fieldProp%species%mrt( iLevel )%s_mrt( 5, 6 ) &
            & = ( omega_bulk - omega_kine ) / 3.0_rk
          field( iField )%fieldProp%species%mrt( iLevel )%s_mrt( 5, 7 ) &
            & = ( omega_bulk - omega_kine ) / 3.0_rk

          field( iField )%fieldProp%species%mrt( iLevel )%s_mrt( 6, 5 ) &
            & = ( omega_bulk - omega_kine ) / 3.0_rk
          field( iField )%fieldProp%species%mrt( iLevel )%s_mrt( 6, 7 ) &
            & = ( omega_bulk - omega_kine ) / 3.0_rk

          field( iField )%fieldProp%species%mrt( iLevel )%s_mrt( 7, 5 ) &
            & = ( omega_bulk - omega_kine ) / 3.0_rk
          field( iField )%fieldProp%species%mrt( iLevel )%s_mrt( 7, 6 ) &
            & = ( omega_bulk - omega_kine ) / 3.0_rk

          do iDir=8,10
            field( iField )%fieldProp%species%mrt( iLevel )             &
              & %s_mrt( iDir, iDir ) = omega_kine
          end do

          do iDir=11,19
            field( iField )%fieldProp%species%mrt( iLevel )             &
              & %s_mrt( iDir, iDir ) = mixture%omega_hom
          end do
        end select

        !omega moments = (M^-1 RelaxMat M)*(I+(M^-1 RelMat M)/2)^-1
        tmpMatrix = matmul( matmul( layout%moment%toPdf%A,              &
          & field( iField )%fieldProp%species%mrt( iLevel )%s_mrt ),    &
          & layout%moment%toMoments%A )

        tmpMatrixInv = invert_matrix( Identity + tmpMatrix * 0.5_rk )

        field( iField )%fieldProp%species%mrt( iLevel )%omegaMoments =  &
          & matmul(tmpMatrix, tmpMatrixInv)

        field( iField )%fieldProp%species%mrt( iLevel )      &
          &                      %omegaMomForce = tmpMatrixInv

      end do ! iField
    end if ! mrt?

  end subroutine setParameters_multispecies
! **************************************************************************** !

! **************************************************************************** !
!> First check count number of valid elements (non-solid) in each BC.
!! Then rearrange BC elements list so it contains only valid elements.
!! Update fields%bc%elemLvl%stencilPos fields%bc%elemLvl%posInNghElems
!! accordingly.
!!
  subroutine remove_solid_in_bc( minLevel, maxLevel, nBCs, nFields, &
    &                            levelPointer, levelDesc, globBC, fields )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: minLevel, maxLevel, nBCs, nFields
    !> Level pointer
    integer, intent(in) :: levelPointer(:)
    !> Level Descriptor
    type( tem_levelDesc_type ), intent(in) :: levelDesc(minLevel:maxLevel)
    !>
    type( glob_boundary_type ) :: globBC( nBCs )
    !>
    type( mus_field_type ) :: fields( nFields )
    ! --------------------------------------------------------------------------
    !> number of valid BC elements
    integer :: nValid(minLevel:maxLevel)
    !> positions of only valid elements (non-solid) in bc_elems_type
    integer, allocatable :: posInBCElem(:,:)
    integer :: iBC, iField, iLevel, iElem
    logical :: allWall
    ! --------------------------------------------------------------------------

    write(dbgUnit(3),*) 'Get into remove_solid_in_bc routine'

    do iBC = 1, nBCs

      !! \todo: do this for both SBB and LIBB?
      allWall = .true.
      do iField = 1, nFields
        if( fields(iField)%bc(iBC)%BC_kind(1:4) /= 'wall' ) then
        ! if( trim(fields(iField)%bc(iBC)%BC_kind) /= 'wall' ) then
          allWall = .false.
          exit
        end if
      end do

      if ( .not. allWall ) then

        write(logUnit(3),"(A)") 'Try to remove solid from BC: ' &
          &                     // trim(globBC(iBC)%label)

        ! ----------------------------------------------------------------------
        allocate( posInBCElem( maxval(globBC(iBC)%nElems), minLevel:maxLevel ) )

        ! check number of valid elements in each BC
        call check_solid_in_bc( minLevel, maxLevel, levelPointer, levelDesc,   &
          &                     globBC(iBC)%nElems, globBC(iBC)%elemLvl,       &
          &                     nValid, posInBCElem )

        ! remove solid in globBC%elemLvl
        call rearrange_bc_elems( minLevel, maxLevel, nValid, posInBCElem, &
          &                      globBC(iBC)%nElems, globBC(iBC)%elemLvl )

        globBC(iBC)%nElems(:) = nValid(:)

        do iField = 1, nFields
          ! re-arrange stencilPos, posInNghElems only if this BC requires
          ! neighbors
          if ( fields(iField)%bc(iBC)%nNeighs > 0 ) then
            do iLevel = minLevel, maxLevel

              ! ----------------------------------------------------------------
              do iElem = 1, nValid( iLevel )
                fields(iField)%bc(iBC)%elemLvl(iLevel)%stencilPos( iElem ) &
                  & = fields(iField)%bc(iBC)%elemLvl(iLevel) &
                  &     %stencilPos( posInBCElem(iElem, iLevel) )

                fields(iField)%bc(iBC)%elemLvl(iLevel)%posInNghElems( iElem ) &
                  & = fields(iField)%bc(iBC)%elemLvl(iLevel) &
                  &     %posInNghElems( posInBCElem(iElem, iLevel) )
              end do ! iElem = 1, nValid
              ! ----------------------------------------------------------------

            end do ! iLevel
          end if ! nNeighs > 0
        end do ! iField

        deallocate( posInBCElem )

      end if ! not allWall

    end do ! iBC

    write(dbgUnit(3),*) 'Leave out of remove_solid_in_bc routine'
    write(dbgUnit(3),*) ''

  end subroutine remove_solid_in_bc
! **************************************************************************** !

! **************************************************************************** !
  !> Check if a BC is wall or symmetry for all fields
  pure function mus_check_allWall( nFields, fields, iBC ) result ( allWall )
    ! --------------------------------------------------------------------------
    integer,                intent(in) :: iBC, nFields
    type( mus_field_type ), intent(in) :: fields( nFields )
    logical :: allWall
    ! --------------------------------------------------------------------------
    integer :: iField
    ! --------------------------------------------------------------------------

    allWall = .true.
    do iField = 1, nFields
      if( trim(fields(iField)%bc(iBC)%BC_kind) /= 'wall' .or. &
        & trim(fields(iField)%bc(iBC)%BC_kind) /= 'symmetry' ) then
        allWall = .false.
        exit
      end if ! not wall
    end do  ! iField

  end function mus_check_allWall
! **************************************************************************** !

! **************************************************************************** !
  !> This routine checks for the existence of symmetric boundaries and
  !! returns the boundary IDs which are defined as symmetry
  subroutine mus_field_getSymmetricBCs(symmetricBCs, nSymBCs, nBCs, nFields, &
    & field)
    ! --------------------------------------------------------------------------
    !> number of boundary conditions
    integer, intent(in) :: nBCs
    !> Symmetric boundary ids
    integer, intent(out) :: symmetricBCs(nBCs)
    !> Number of symmetric boundary conditions
    integer, intent(out) :: nSymBCs
    !> number of fields
    integer, intent(in) :: nFields
    !> all fields to access their boundary definitions
    type(mus_field_type), intent(in) :: field(nFields)
    ! --------------------------------------------------------------------------
    integer :: iField, iBC
    logical :: isSymmetry(nFields)
    ! --------------------------------------------------------------------------
    nSymBCs = 0
    do iBC = 1, nBCs
      isSymmetry = .false.
      do iField = 1, nFields
        if (trim(field(iField)%bc(iBC)%BC_kind) == 'symmetry') then
          isSymmetry(iField) = .true.
        end if
        if ( any(isSymmetry) ) then
          if (.not. all(isSymmetry) ) then
            call tem_abort('Error: Not all fields boundary label:'&
              &          //trim(field(iField)%bc(iBC)%label)//' are symmetry!')
          end if
          nSymBCs = nSymBCs + 1
          symmetricBCs(nSymBCs) = iBC
        end if
      end do
    end do

  end subroutine mus_field_getSymmetricBCs
! **************************************************************************** !

  ! ************************************************************************** !
  !> This routines act as a destructor for field type.
  !! Only allocatable arrays which are allocated in mus_construct routine
  !! are deallocated.
  !! KM: DO NOT DESTROY FIELD ARRAY AS IT CONTAINS ALL CONFIG INFO
  subroutine mus_field_cleanup(me, schemeHeader, minLevel, maxLevel, nBCs, &
    &                          nFields )
    ! --------------------------------------------------------------------------
    !> single field type
    type(mus_field_type), intent(inout) :: me(:)
    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    !> minlevel
    integer, intent(in) :: minLevel
    !> maxlevel
    integer, intent(in) :: maxLevel
    !> Number of boundary conditions
    integer, intent(in) :: nBCs
    !> Number of fields
    integer, intent(in) :: nFields
    ! --------------------------------------------------------------------------
    integer :: iFld
    ! --------------------------------------------------------------------------
    write(dbgUnit(1),"(A)") 'Enter mus_field_cleanup'
    ! Cleanup allocatable arrays in field properties
    select case( trim(schemeHeader%kind) )
    case('fluid', 'fluid_incompressible', 'isotherm_acEq')
      call mus_fluid_cleanup(me(1)%fieldProp%fluid)
    end select

    ! Cleanup field source
    write(dbgUnit(1),"(A)") 'Enter mus_field_cleanup'
    do iFld = 1, nFields
      call mus_source_cleanup(me(iFld)%source)
    end do

    ! Cleanup field boundary
    do iFld = 1, nFields
      write(dbgUnit(3),*) "Cleanup BC: iField=",iFld
      call mus_fieldBC_cleanup( me       = me(iFld)%bc, &
        &                       nBCs     = nBCs,        &
        &                       minLevel = minLevel,    &
        &                       maxLevel = maxLevel     )
    end do

  end subroutine mus_field_cleanup
  ! ************************************************************************** !

end module mus_field_module
! **************************************************************************** !
