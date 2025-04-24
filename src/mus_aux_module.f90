! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2011-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2011-2013,2020-2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!! Auxiliary functions to maintain the simulation
!!
module mus_aux_module

  ! include musubi modules
  use mus_abortCriteria_module, only: mus_abortCriteria_type
  use mus_param_module,              only: mus_param_type, mus_latticeUnit_type
  use mus_geom_module,               only: mus_geom_type
  use mus_tools_module,              only: perform_checks, &
    &                                      check_streaming_layout
  use mus_interpolate_module,        only: mus_init_interpolate
  use mus_source_module,             only: mus_init_sourceTerms
  use mus_transport_var_module,      only: mus_init_transport_var
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_geomIncr_module,           only: mus_geomIncr, mus_init_geomIncr
  use mus_interpolate_tools_module,  only: debug_dependencies, dump_intpLists
  use mus_time_module,               only: mus_time_homogenize,                &
    &                                      mus_timeControl_homogenize
  use mus_IBM_module,                only: mus_IBM_setParentIDs, &
    &                                      mus_IBM_globType
  use mus_fluid_module,              only: mus_init_fluid
  use mus_gradData_module,           only: mus_init_gradData
  use mus_relaxationParam_module,    only: mus_update_relaxParamKine,       &
    &                                      mus_update_viscKine,             &
    &                                      mus_update_relaxParamFromViscSTfun
  use mus_turbulence_module,         only: mus_turb_updateViscOfTurbWall
  use mus_field_module,              only: setParameters_multispecies
  use mus_tracking_module,           only: mus_init_tracker
  use mus_restart_module,            only: mus_writeRestart
  use mus_timer_module,              only: mus_timerHandles
  use mus_physics_module,            only: mus_physics_type
  use mus_ppInfo_module,             only: mus_print_ppInfo
  use mus_bndForce_module,           only: mus_init_bndForce, mus_calcBndForce

  ! include treelm modules
  use env_module,               only: rk, PathLen, pathSep, long_k
  use treelmesh_module,         only: treelmesh_type
  use tem_aux_module,           only: tem_print_execInfo, utc_date_string, &
    &                                 tem_abort
  use tem_debug_module,         only: dbgUnit
  use tem_restart_module,       only: tem_init_restart
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_tracking_module,      only: tem_tracker
  use tem_convergence_module,   only: tem_convergence_check, &
    &                                 tem_init_convergence
  use tem_timeControl_module,   only: tem_timeControl_check, &
    &                                 tem_timeControl_update
  use tem_simControl_module,    only: tem_simControl_syncUpdate
  use tem_time_module,          only: tem_time_dump, tem_time_type
  use tem_debug_module,         only: main_debug
  use tem_solveHead_module,     only: tem_solveHead_type
  use tem_logging_module,       only: logUnit
  use tem_global_module,        only: tem_global_type
  use tem_depend_module,        only: tem_init_depend
  use tem_spacetime_fun_module, only: tem_create_subTree_of_st_funList
  use tem_dyn_array_module,     only: dyn_intArray_type
  use tem_timer_module,         only: tem_getTimerVal
  use tem_general_module,       only: tem_general_type
  use tem_operation_var_module, only: tem_opVar_reduction_transient_update


  implicit none
  private

  public :: check_flow_status
  public :: mus_init_aux
  public :: mus_update_relaxParams
  public :: mus_banner
  public :: mus_dumpData


contains


  ! ------------------------------------------------------------------------ !
  !> This routine performs several tasks: geometry increment, time updating,
  !! tracking, density checking, restart
  subroutine check_flow_status( scheme, geometry, general, physics, mus_aborts,&
    &                           restart_triggered )
    ! -------------------------------------------------------------------- !
    !> containers for the different schemes
    type(mus_scheme_type), target, intent(inout) :: scheme
    !> geometry infomation
    type(mus_geom_type),intent(inout) :: geometry
    !> Global parameters
    type(tem_general_type), intent(inout)  :: general
    !> physics conversion tyoe
    type(mus_physics_type), intent(in) :: physics
    type(mus_abortCriteria_type), intent(in) :: mus_aborts
    !> Indication whether a restart output was triggered
    logical, intent(inout) :: restart_triggered
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    ! Force on boundary elements are computed from post-collision so calculate
    ! here after apply source and before tracking output
    write(logUnit(10), "(A)") 'Calculated force on boundary elements'
    call mus_calcBndForce( bndForce   = geometry%bndForce,             &
      &                    bndMoment  = geometry%bndMoment,            &
      &                    posInBndID = geometry%posInBndID,           &
      &                    nBCs       = geometry%boundary%nBCtypes,    &
      &                    field      = scheme%field,                  &
      &                    globBC     = scheme%globBC,                 &
      &                    minLevel   = geometry%tree%global%minLevel, &
      &                    maxLevel   = geometry%tree%global%maxLevel, &
      &                    state      = scheme%state,                  &
      &                    pdf        = scheme%pdf,                    &
      &                    levelDesc  = scheme%levelDesc,              &
      &                    layout     = scheme%layout,                 &
      &                    varSys     = scheme%varSys,                 &
      &                    physics    = physics                        )

    ! Call the geometry increment module which  performs solidification
    ! or liquification based on certain criteria defined by the user
    if( geometry%dynamicGeom ) then
      call mus_geomIncr(                              &
        &    geometry    = geometry, scheme = scheme, &
        &    commPattern = general%commPattern,       &
        &    general     = general                    )
    end if

    ! Perform run-time checks to ensure we did not encounter some unphysical
    ! state.
    call perform_checks(                               &
      &    scheme     = scheme,                        &
      &    minLevel   = geometry%tree%global%minLevel, &
      &    maxLevel   = geometry%tree%global%maxLevel, &
      &    general    = general,                       &
      &    mus_aborts = mus_aborts,                    &
      &    initCheck  = .false.                        )

    ! check for convergence only if abortCriteria is set for steady_state
    if( general%simControl%abortCriteria%steady_state) then
      ! check for steady state convergence based on convergence criteria
      ! defined in convergence table
      call tem_convergence_check(                      &
        &    me     = general%simControl%abortCriteria &
        &                     %convergence,            &
        &    time   = general%simControl%now,          &
        &    status = general%simControl%status,       &
        &    varSys = scheme%varSys,                   &
        &    tree   = geometry%tree                    )
    end if

    ! Update time reduction variables if there are any
    if (scheme%redTransVarMap%varPos%nVals > 0) then
      call tem_opVar_reduction_transient_update(                          &
        &    redTransVarPos = scheme%redTransVarMap%varPos                &
        &                     %val(1:scheme%redTransVarMap%varPos%nVals), &
        &    varSys         = scheme%varSys,                              &
        &    tree           = geometry%tree,                              &
        &    time           = general%simControl%now                      )
    end if

    ! in musubi, advance time step separtely in control_routine so no need
    ! to pass dt to syncUpdate
    ! This sync update: check for stop file and time
    ! control interval trigger, communicate status bits and update timeControl
    call tem_simControl_syncUpdate( me   = general%simControl, &
      &                             proc = general%proc        )

    ! Dump tracking and restart if they are active
    call mus_dumpData(                               &
      &    scheme            = scheme,               &
      &    tree              = geometry%tree,        &
      &    restart_triggered = restart_triggered,    &
      &    general           = general,              &
      &    levelPointer      = geometry%levelPointer )

  end subroutine check_flow_status
  ! ------------------------------------------------------------------------ !

  ! ------------------------------------------------------------------------ !
  !> Initialize musubi solverHead and print musubi banner to screen
  subroutine mus_banner( solver )
    ! -------------------------------------------------------------------- !
    !> solver definition
    type(tem_solveHead_type), intent(in) :: solver
    ! -------------------------------------------------------------------- !
    character(len=26) :: dat_string
    ! -------------------------------------------------------------------- !
    write(logUnit(1),"(A)")''
    write(logUnit(1),"(A)")" .___  ___.  __    __       _______. __    __  .______    __  "
    write(logUnit(1),"(A)")" |   \/   | |  |  |  |     /       ||  |  |  | |   _  \  |  | "
    write(logUnit(1),"(A)")" |  \  /  | |  |  |  |    |   (----`|  |  |  | |  |_)  | |  | "
    write(logUnit(1),"(A)")" |  |\/|  | |  |  |  |     \   \    |  |  |  | |   _  <  |  | "
    write(logUnit(1),"(A)")" |  |  |  | |  `--'  | .----)   |   |  `--'  | |  |_)  | |  | "
    write(logUnit(1),"(A)")" |__|  |__|  \______/  |_______/     \______/  |______/  |_"  &
      &              //trim(solver%version)
    write(logUnit(1),"(A)") ''
    write(logUnit(1),"(A)") " (C) 2012 German Research School"&
      &                     //" for Simulation Sciences"
    write(logUnit(1),"(A)") " (C) 2013-2020 University Siegen"
    write(logUnit(1),"(A)") " (C) 2021 German Aerospace Center (DLR) -"
    write(logUnit(1),"(A)") "          Institute of Software Methods" &
      &                     //" for Product Virtualization"
    write(logUnit(1),"(A)") " "
    call tem_print_execInfo()
    call mus_print_ppInfo()
    write(logUnit(1),"(A)") ''
    dat_string = utc_date_string()
    write(logUnit(1),"(A)") "Run at: "//trim(dat_string)
    write(logUnit(1),"(A)") ''

  end subroutine mus_banner
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Init auxiliary features such as interpolation boundaries, restart and
  !! the tracker
  subroutine mus_init_aux( scheme, geometry, params)
    ! -------------------------------------------------------------------- !
    !> containers for the different schemes
    type(mus_scheme_type), target, intent(inout) :: scheme
    !> geometry information
    type(mus_geom_type), intent(inout)           :: geometry
    !> global parameters
    type(mus_param_type),intent(inout)           :: params
    ! -------------------------------------------------------------------- !
    integer :: iConv
    integer :: iLevel, minLevel, maxLevel
    ! -------------------------------------------------------------------- !
    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel

    ! When the restart is read from separate restart table, we start the
    ! simulation from the time step given in restart file. In the case,
    ! when restart is read from initial condition table, the simulation start
    ! time step is taken as the one defined in configuration file
    call mus_time_homogenize(                              &
      &    me          = params%general%simControl%now,    &
      &    dt          = params%physics%dtLvl( maxLevel ), &
      &    readRestart = params%general%restart%controller &
      &                        %readRestart                )
    call mus_timeControl_homogenize(                       &
      &    me     = params%general%simControl%timeControl, &
      &    dt     = params%physics%dtLvl( maxLevel ),      &
      &    reqInt = params%reqInterval                     )

    !> initialize fluid type which contains relaxation parameter
    !! and function pointers to get mrt paramter and nonEqScaling factor
    !! for interpolation
    select case( trim(scheme%header%kind) )
    case('fluid', 'fluid_incompressible', 'fluid_GNS', &
        & 'fluid_incompressible_GNS', 'isotherm_acEq'  )
      if (scheme%nFields > 1) then
        call tem_abort('chosen scheme kind supports only one field')
      end if
      ! initialize fluid viscosity relaxation paramters
      call mus_init_fluid(                                   &
        &    me           = scheme%field(1)%fieldProp%fluid, &
        &    physics      = params%physics,                  &
        &    schemeHeader = scheme%header,                   &
        &    minLevel     = minLevel,                        &
        &    maxLevel     = maxLevel,                        &
        &    levelDesc    = scheme%levelDesc,                &
        &    pdf          = scheme%pdf,                      &
        &    stencil      = scheme%layout%fStencil,          &
        &    general      = params%general,                  &
        &    tNow         = params%general%simControl%now    )
    end select

    ! Initialize gradient data. Required for LES turbulent and evaluating
    ! gradient of a variable
    allocate(scheme%gradData(minLevel:maxLevel))
    do iLevel = minLevel, maxLevel
      call mus_init_gradData(                             &
        &    me        = scheme%gradData(iLevel),         &
        &    neigh     = scheme%pdf(iLevel)%neigh,        &
        !&    levelDesc = scheme%levelDesc(iLevel),        &
        &    stencil   = scheme%layout%fStencil,          &
        &    nSize     = scheme%pdf(iLevel)%nSize,        &
        &    nSolve    = scheme%pdf(iLevel)%nElems_solve, &
        &    nScalars  = scheme%varSys%nScalars           )
    end do

    ! create subTree for all spacetime function in the linked list of
    ! spacetime function
    call tem_create_subTree_of_st_funList(  &
      &    me      = scheme%st_funList,     &
      &    tree    = geometry%tree,         &
      &    bc_prop = geometry%boundary,     &
      &    stencil = scheme%layout%fStencil )

    ! initialize the source terms for all fields and global source
    call mus_init_sourceTerms(                        &
      &    field        = scheme%field,               &
      &    nFields      = scheme%nFields,             &
      &    globSrc      = scheme%globSrc,             &
      &    varSys       = scheme%varSys,              &
      &    tree         = geometry%tree,              &
      &    bc_prop      = geometry%boundary,          &
      &    nElems_solve = scheme%pdf(:)%nElems_solve, &
      &    levelDesc    = scheme%levelDesc,           &
      &    stencil      = scheme%layout%fStencil      )

    ! initialize transport variables like velocity for passive scalar
    call mus_init_transport_var(                      &
      &    me           = scheme%transVar,            &
      &    varSys       = scheme%varSys,              &
      &    tree         = geometry%tree,              &
      &    nElems_solve = scheme%pdf(:)%nElems_solve, &
      &    levelDesc    = scheme%levelDesc            )

    ! verify correct settings for the streaming layout
    call check_streaming_layout( minLevel, maxLevel )

    ! dynamic load balance time control homogenize
    if ( params%general%balance%dynamic ) then
      call mus_timeControl_homogenize(                     &
        &     me     = params%general%balance%timeControl, &
        &     dt     = params%physics%dtLvl( maxLevel ),   &
        &     reqInt = params%reqInterval                  )
    endif

    if ( ( params%general%restart%controller%writeRestart  .or. &
      &   params%general%restart%controller%readRestart ) ) then

      ! initialize the restart
      write(logUnit(1),*) 'Initializing restart...'

      call tem_init_restart(                        &
        &    me           = params%general%restart, &
        &    solver       = params%general%solver,  &
        &    varMap       = scheme%stateVarMap,     &
        &    tree         = geometry%tree           )

      call mus_timeControl_homogenize(                               &
        &    me     = params%general%restart%controller%timeControl, &
        &    dt     = params%physics%dtLvl( maxLevel ),              &
        &    reqInt = params%reqInterval                             )
    end if


    ! initialize tracking objects.
    call mus_init_tracker(       &
      &    scheme    = scheme,   &
      &    geometry  = geometry, &
      &    params    = params    )

    ! convergence objects
    if ( params%general%simControl%abortCriteria%steady_state ) then
      write(logUnit(1),*) 'Initializing convergence...'

      do iConv = 1, size( params%general%simControl%abortCriteria%convergence)
        call mus_timeControl_homogenize(                          &
          &    me = params%general%simControl%abortCriteria       &
          &               %convergence(iConv)%header%timeControl, &
          &    dt     = params%physics%dtLvl( maxLevel ),         &
          &    reqInt = params%reqInterval                        )
      end do

      call tem_init_convergence(                         &
        &    me       = params%general%simControl        &
        &                     %abortCriteria%convergence,&
        &    tree     = geometry%tree,                   &
        &    bc_prop  = geometry%boundary,               &
        &    stencil  = scheme%layout%fStencil,          &
        &    globProc = params%general%proc,             &
        &    varSys   = scheme%varSys                    )
    end if

    if( minLevel /= maxlevel ) then
      write(logUnit(1),"(A)") 'Initializing interpolation...'
      ! initialize the interpolation
      call mus_init_interpolate(                      &
        &    intp         = scheme%intp,              &
        &    levelDesc    = scheme%levelDesc,         &
        &    schemeHeader = scheme%header,            &
        &    stencil      = scheme%layout%fStencil,   &
        &    minLevel     = minLevel,                 &
        &    maxLevel     = maxLevel,                 &
        &    fieldProp    = scheme%field(:)%fieldProp )
    end if
    if( main_debug%debugDependencies ) then
      call debug_dependencies( scheme%intp, scheme%levelDesc,          &
        &                      geometry%tree, params%general%proc%rank )
      call dump_intpLists( minLevel, maxLevel, scheme%intp%config%order, &
        &                  scheme%levelDesc, params%general%proc%rank    )
    end if

    ! Boundary force calculation is valid only for single field schemes
    ! like fluid and fluid_incompressible so initialize only if nFields = 1
    if (geometry%boundary%nBCtypes > 0 .and. scheme%nFields==1) then
      call mus_init_bndForce(bndForce     = geometry%bndForce,  &
        &                    bndMoment    = geometry%bndMoment, &
        &                    bc_prop      = geometry%boundary,  &
        &                    schemeHeader = scheme%header,      &
        &                    bc           = scheme%field(1)%bc  )
    end if

    ! initialize the surface data for the immersed boundary method
    write(logUnit(1),*)'Initializing IBM surface data'
    call mus_IBM_setParentIDs(                 &
      &    nIBMs     = geometry%globIBM%nIBMs, &
      &    me        = geometry%globIBM%ibm,   &
      &    levelDesc = scheme%levelDesc,       &
      &    tree      = geometry%tree           )

    if( geometry%dynamicGeom ) then
      write(logUnit(1),*) 'Initializing geomIncr ...'
      call mus_init_geomIncr(                           &
        &    me     = geometry%geomIncr,                &
        &    varSys = scheme%varSys,                    &
        &    dt     = params%physics%dtLvl( maxLevel ), &
        &    reqInt = params%reqInterval                )
      write(logUnit(1),*) 'Done initializing geomIncr.'
    end if ! Time for geom incr

    call tem_horizontalSpacer( after = 1, fUnit = logUnit(1) )

  end subroutine mus_init_aux
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Set relaxation parameters for MRT
  !!
  subroutine mus_update_relaxParams( scheme, iLevel, tNow, physics, lattice, &
    &                                nBCs )
    ! -------------------------------------------------------------------- !
    !> scheme type
    type(mus_scheme_type), intent(inout) :: scheme
    !> level
    integer, intent(in)            :: iLevel
    !> global parameters
    type(tem_time_type),intent(in) :: tNow
    !> contains factors to convert physical to lattice unit and vice versa
    type(mus_physics_type), intent(in) :: physics
    !> contains lattice dx and dt
    type(mus_latticeUnit_type), intent(in) :: lattice
    !> Number of boundaries
    integer, intent(in) :: nBCs
    ! -------------------------------------------------------------------- !
    integer :: iBnd
    ! -------------------------------------------------------------------- !

    select case(trim(scheme%header%kind))
    case('fluid', 'fluid_incompressible')
      ! Update kinematic viscosity from STfun and calculate turbulent viscosity
      ! from velocity gradient or nonEqPDF
      call mus_update_viscKine(                                              &
        &    viscKine          = scheme%field(1)%fieldProp%fluid%viscKine,   &
        &    state             = scheme%state(iLevel)%val(:,                 &
        &                               scheme%pdf(iLevel)%nNow),            &
        &    neigh             = scheme%pdf(iLevel)%neigh,                   &
        &    auxField          = scheme%auxField(iLevel)%val,                &
        &    gradData          = scheme%gradData(iLevel),                    &
        &    nSize             = scheme%pdf(iLevel)%nSize,                   &
        &    nFluids           = scheme%pdf(iLevel)%nElems_fluid,            &
        &    nGhostFromCoarser = scheme%pdf(iLevel)%nElems_ghostFromCoarser, &
        &    nGhostFromFiner   = scheme%pdf(iLevel)%nElems_ghostFromFiner,   &
        &    nHalo             = scheme%pdf(iLevel)%nElems_halo,             &
        &    baryOfTotal       = scheme%levelDesc(iLevel)%baryOfTotal,       &
        &    varSys            = scheme%varSys,                              &
        &    iLevel            = iLevel,                                     &
        &    layout            = scheme%layout,                              &
        &    tNow              = tnow,                                       &
        &    convFac           = physics%fac(iLevel),                        &
        &    dxL               = lattice%dxLvl(iLevel),                      &
        &    dtL               = lattice%dtLvl(iLevel),                      &
        &    derVarPos         = scheme%derVarPos(1),                        &
        &    turb              = scheme%field(1)%fieldProp%fluid%turbulence, &
        &    nNwtn             = scheme%field(1)%fieldProp%fluid%nNwtn,      &
        &    Grad              = scheme%Grad                                 )

      ! Update turbulent viscosity for boundary elements of turbulent wall
      ! only for smagorinsky
      if (scheme%field(1)%fieldProp%fluid%turbulence%active .and.       &
        & trim(scheme%field(1)%fieldProp%fluid%turbulence%config%model) &
        & == 'smagorinsky') then
        do iBnd = 1, nBCs
          select case (trim(scheme%field(1)%bc(iBnd)%BC_kind))
          case ( 'turbulent_wall', 'turbulent_wall_noneq_expol', &
            & 'turbulent_wall_eq' )
            call mus_turb_updateViscOfTurbWall(                              &
              &    turbData     = scheme%field(1)%fieldProp%fluid%turbulence &
              &                     %dataOnLvl(iLevel),                      &
              &    viscTurbWall = scheme%field(1)%bc(iBnd)%turbwallFunc      &
              &                         %dataOnLvl(iLevel)%tVisc(:),         &
              &    nElems_bnd   = scheme%globBC(iBnd)%nElems(iLevel),        &
              &    elemPos      = scheme%globBC(iBnd)%elemLvl(iLevel)%elem   &
              &                         %val(:)                              )
          end select
        end do
      end if

      ! update fluid relaxation parameter omega kine
      call mus_update_relaxParamKine(                                   &
        &    viscKine     = scheme%field(1)%fieldProp%fluid%viscKine,   &
        &    turb         = scheme%field(1)%fieldProp%fluid%turbulence, &
        &    nSolve       = scheme%pdf(iLevel)%nElems_solve,            &
        &    iLevel       = iLevel                                      )

    case('multispecies_gas','multispecies_liquid')
      call setParameters_multispecies( &
        &    field   = scheme%field,   &
        &    nFields = scheme%nFields, &
        &    mixture = scheme%mixture, &
        &    header  = scheme%header,  &
        &    layout  = scheme%layout,  &
        &    iLevel  = iLevel,         &
        &    tNow    = tNow            )
    end select

  end subroutine mus_update_relaxParams
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine dumps tracking and restart when timeControl is active
  subroutine mus_dumpData( scheme, tree, levelPointer, general, &
    &                      restart_triggered                    )
    ! -------------------------------------------------------------------- !
    !> scheme type
    type( mus_scheme_type ), intent(inout) :: scheme
    !> Treelmesh data
    type(treelmesh_type), intent(in)       :: tree
    !> Level Pointer
    integer, intent(in) :: levelPointer(:)
    !> Global parameters
    type( tem_general_type ), intent(inout)  :: general
    !>
    logical, intent(inout) :: restart_triggered
    ! -------------------------------------------------------------------- !

    if (scheme%track%control%active) then
      call tem_tracker(                       &
        &    track      = scheme%track,       &
        &    simControl = general%simControl, &
        &    tree       = tree,               &
        &    varSys     = scheme%varSys       )
    end if

    if (general%restart%controller%writeRestart) then
      ! check time control and update it if triggered
      call tem_timeControl_check(                             &
        &    me     = general%restart%controller%timeControl, &
        &    now    = general%simControl%now,                 &
        &    comm   = general%proc%comm,                      &
        &    triggered = restart_triggered                    )

      if (restart_triggered) then
        if ( general%restart%lastWritten%iter &
           & == general%simControl%now%iter   ) then
          write(logUnit(1),*) ' Restart file already written for current' &
            &                 // ' point in time!'
        else
          write(logUnit(1),*) ' Writing restart at:'
          call tem_time_dump( general%simControl%now, logUnit(1) )
          call mus_writeRestart( levelPointer = levelPointer,             &
            &                    restart      = general%restart,          &
            &                    scheme       = scheme,                   &
            &                    tree         = tree,                     &
            &                    timing       = general%simControl%now,   &
            &                    timerHandle  = mus_timerHandles%wRestart )
          write(logUnit(1),*) 'Done writing'
        end if
      end if
    end if

  end subroutine mus_dumpData
  ! ------------------------------------------------------------------------ !

end module mus_aux_module
! **************************************************************************** !
