! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2011-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2012, 2021, 2025 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2025 Tristan Vlogman <t.g.vlogman@utwente.nl>
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
!> In this module, the control structure for computing each time step is
!! established.
!!
!! Various different control structures can be defined and set by function
!! pointers.
!! The current main routine is recursive_multilevel which recursively computes
!! the new time step for all the schemes defined on all levels, starting from
!! the coarsest.
!!
module mus_control_module

  ! include treelm modules
  use env_module,         only: rk, labelLen
  use tem_timer_module,   only: tem_startTimer, tem_stopTimer
  use tem_general_module, only: tem_general_type
  use tem_debug_module,   only: main_debug, dbgUnit
  use tem_time_module,    only: tem_time_advance, tem_time_dump
  use tem_logging_module, only: logUnit
  use tem_varSys_module,  only: tem_varSys_type
  use tem_stencil_module, only: tem_stencilHeader_type

  ! include musubi modules
  use mus_pdf_module,                only: pdf_data_type, mus_swap_now_next
  use mus_param_module,              only: mus_param_type
  use mus_geom_module,               only: mus_geom_type
  use mus_bc_general_module,         only: set_boundary
  use mus_aux_module,                only: mus_update_relaxParams
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_field_module,              only: mus_field_type
  use mus_source_module,             only: mus_apply_sourceTerms, &
    &                                      mus_update_sourceVars
  use mus_debug_module,              only: dump_debug_Info
  use mus_IBM_module,                only: mus_inamuro_IBM, mus_buildBuffIBM
  use mus_timer_module,              only: mus_timerHandles, nStages
  use mus_physics_module,            only: mus_convertFac_type
  use mus_auxField_module,           only: mus_calcAuxFieldAndExchange, &
    &                                      mus_intpAuxFieldCoarserAndExchange
  use mus_derVarPos_module,          only: mus_derVarPos_type

  ! include particle musubi modules
  use mus_particle_module,           only: mus_particle_group_type
  use mus_particle_DPS_module,       only: mus_particles_updateFluidVolumeFraction
  use mus_particle_logging_module,   only: mus_particles_logdata_DPS, &
    &                                      mus_particles_logdata_MEM
  use mus_particle_creator_module,   only: check_and_create_new_particles_DPS, &
    &                                      check_and_create_new_particles_MEM, &
    &                                      particle_creator

  implicit none

  private

  public :: mus_control_type
  public :: mus_init_control

  !> Datatype containing mapping of control routines to function pointers
  type mus_control_type
    type(mus_scheme_type),         pointer :: scheme => null()
    type(mus_geom_type),           pointer :: geometry => null()
    type(mus_param_type),          pointer :: params => null()
    type(mus_particle_group_type), pointer :: particleGroup => null()

    logical :: DPS_do_volfract = .false.
    logical :: DPS_do_advance = .true.
    integer :: curlvl = 0

    procedure(computation), pointer :: do_computation => null()
    procedure(update_particles_if), pointer :: check_particles => null()
    procedure(update_particles_if), pointer :: advance_particles => null()
  end type mus_control_type

  abstract interface
    !> Interface describes the main control routine which does computation
    !! set boundary and check flow status
    subroutine computation( me, iLevel)
      import :: mus_control_type
      !> self control type
      class(mus_control_type) :: me
      !> Level counter variable
      integer, intent(in) :: iLevel
    end subroutine computation

    subroutine update_particles_if(me)
      import :: mus_control_type
      !> self control type
      class(mus_control_type) :: me
    end subroutine update_particles_if
  end interface

  integer, save :: iStage = 0
  logical, save :: running = .false.


contains


  ! ------------------------------------------------------------------------ !
  !> This routines sets the function pointer to main control routine
  !!
  !! control routine is chosen based on the type of simulation.
  !! like single level, multi-level because multilevel requires
  !! recursive routine
  !!
  !! in the lua file you can now select
  !! `control_routine = '...'`
  !! where possible values are currently
  !!
  !! - `benchmark`: strongly reduced control routine with only single level,
  !!                no sources, etc.
  !! mainly for benchmarking
  !!
  !! - `multiLevel`: full multilevel, multiLevel routine
  !! - if nothing is given, the full multilevel, multiLevel routine is chosen
  subroutine mus_init_control( controlRoutine, me, minLevel, &
    &                          maxLevel, particle_kind       )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: controlRoutine
    !> contains function pointer to point control routine
    type(mus_control_type), intent(inout) :: me
    integer, intent(in) :: minLevel, maxLevel
    !> string containing kind of solid particles for coupled LBM-DEM simulations
    !! Can be fully resolved 'MEM' or unresolved 'DPS', 'DPS_twoway' 
    !! or 'DPS_oneway'
    character(len=labelLen), intent(in) :: particle_kind
    ! -------------------------------------------------------------------- !

    ! Select according to special need
    if ( trim( controlRoutine ) == 'benchmark' ) then

      me%do_computation => do_benchmark
      write(logUnit(5),"(A)") "Select benchmark control routine."

    else if ( trim( controlRoutine ) == 'debug' ) then

      ! me%do_computation => do_recursive_debug
      write(logUnit(5),"(A)") "Select debug recursive control routine."

    else

      ! No special requirement from user, then select by single or multi levels
      if ( minLevel == maxLevel ) then
        me%do_computation => do_fast_singleLevel
        write(logUnit(5),"(A)") "Select fast single level control routine."

        ! Choose control routine for coupled LBM-DEM simulations of particulate 
        ! flows.
        select case( trim( particle_kind) )
        case ('MEM')
          me%check_particles => check_particles_MEM
          me%advance_particles => advance_particles_MEM
          write(logUnit(5),"(A)") "Select fast single level control routine for ", &
            & " fully coupled MEM particles"
        case ('DPS')
          ! TV: use alternate compute routine for particles!
          me%DPS_do_VolFract = .true.
          me%DPS_do_advance = .true.
          me%check_particles => check_particles_DPS
          me%advance_particles => advance_particles_DPS
          write(logUnit(5),"(A)") "Select fast single level control routine for ", &
            & " fully coupled DPS particles using the Generalized Navier-Stokes equations"
        case ('DPS_twoway')
          me%DPS_do_VolFract = .false.
          me%DPS_do_advance = .true.
          me%check_particles => check_particles_DPS
          me%advance_particles => advance_particles_DPS
          write(logUnit(5),"(A)") "Select fast single level control routine ", &
            & "for two-way coupled DPS particles (neglecting the effect of local ", &
            & "fluid volume fraction )."
        case ('DPS_oneway')
          me%DPS_do_VolFract = .false.
          me%DPS_do_advance = .false.
          me%check_particles => check_particles_DPS
          me%advance_particles => advance_particles_DPS
          write(logUnit(5),"(A)") "Select fast single level control routine ", &
            & "for one-way coupled DPS particles."
        case default
          write(logUnit(5),"(A)") "Particulate flow simulations disabled!"
        end select
      else
        me%do_computation => do_recursive_multiLevel
        write(logUnit(5),"(A)") "Select recursive multi level control routine."
      end if

    end if

  end subroutine mus_init_control
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Main control routine: Update the time step for all levels.
  !! Main steps:
  !!   * if iLevel < maxLevel do recursive at iLevel+1
  !!   * do BC at iLevel
  !!   * do auxField calculation at iLevel
  !!   * do compute kernel at iLevel
  !!   * do apply source at iLevel
  !!   * do do_IntpFinerAndExchange at iLevel if iLevel < maxLevel
  !!     * intp My Coarser ghost (iLevel) from Finer (iLevel+1)
  !!     * do exchange bufferFromFiner at iLevel
  !!   * exchange buffer at iLevel
  !!   * exchange bufferFromCoarser at iLevel if iLevel > minLevel
  !!   * do do_intpCoarserAndExchange at iLevel if iLevel < maxLevel
  !!     * intp Finer Ghost (iLevel+1) from my coarser (iLevel)
  !!     * exchange bufferFromCoarser at iLevel+1
  !!
  recursive subroutine do_recursive_multiLevel(me, iLevel)
    ! -------------------------------------------------------------------- !
    !> self control type
    class(mus_control_type) :: me
    !> the current level
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: now, next, iNestingLoop
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------------
    write(logUnit(10), "(A,I0)") 'Entering do_computation with iLevel: ', iLevel

    ! Update auxField dependent source fields before adding source term to state
    ! and auxField such that both auxField and apply_source uses same source in
    ! one multilevel cycle.
    write(logUnit(10), "(A)") 'Update source variables which depend on auxField'
    call mus_update_sourceVars( nFields    = me%scheme%nFields,              &
      &                         field      = me%scheme%field,                &
      &                         globSrc    = me%scheme%globSrc,              &
      &                         varSys     = me%scheme%varSys,               &
      &                         iLevel     = iLevel,                         &
      &                         auxField   = me%scheme%auxField(iLevel)%val, &
      &                         phyConvFac = me%params%physics%fac(iLevel),  &
      &                         derVarPos  = me%scheme%derVarPos             )

    ! when not on finest level, go to next level
    if( iLevel < me%geometry%tree%global%maxLevel ) then
      ! Perform the number of nested time steps on the finer level L+1
      ! according to scaling type :
      !   diffusive: nNesting = 4
      !   acoustic:  nNesting = 2
      do iNestingLoop = 1, me%params%nNesting
        write(logUnit(10), "(A,I0)") 'Nesting loop ', iNestingloop
        call me%do_computation( iLevel+1 )
      end do
    end if

    write(logUnit(10), "(A,I0)") 'Compute on level: ', iLevel
    call start_stageTimer()

    ! update the time counters. MH: checked. Please dont move
    ! Increasing with the smallest time step (maxLevel)
    ! KM: time is advanced here since new time is required to update sources.
    if( iLevel == me%geometry%tree%global%maxLevel ) then
      write(logUnit(10), "(A)") 'Advance time t+dt_maxLevel'
      call tem_time_advance( me     = me%params%general%simControl%now,    &
        &                    sim_dt = me%params%physics%dtLvl( me%geometry &
        &                                   %tree%global%maxLevel )        )
    endif

    write(logUnit(10), "(A)") 'Set boundary condition'
    ! set boundary for each field in current me%scheme
    call set_boundary( field       = me%scheme%field,                  &
      &                pdf         = me%scheme%pdf(iLevel),            &
      &                state       = me%scheme%state(iLevel)%val,      &
      &                levelDesc   = me%scheme%levelDesc(iLevel),      &
      &                tree        = me%geometry%tree,                 &
      &                iLevel      = iLevel,                           &
      &                nBCs        = me%geometry%boundary%nBCtypes,    &
      &                params      = me%params,                        &
      &                layout      = me%scheme%layout,                 &
      &                physics     = me%params%physics,                &
      &                varSys      = me%scheme%varSys,                 &
      &                mixture     = me%scheme%mixture,                &
      &                derVarPos   = me%scheme%derVarPos,              &
      &                globBC      = me%scheme%globBC                  )
    ! -------------------------------------------------------------------------


    write(logUnit(10), "(A)") 'Swap now and next'
    ! swap double buffer index for current level
    call mus_swap_now_next( me%scheme%pdf( iLevel ) )
    now  = me%scheme%pdf(iLevel)%nNow
    next = me%scheme%pdf(iLevel)%nNext

    ! --------------------------------------------------------------------------
    ! Compute auxField from pre-collision state for fluid and ghostFromCoarser
    ! and exchange them if turbulence is active
    call tem_startTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    write(logUnit(10), "(A)") 'Calculate auxField'
    call mus_calcAuxFieldAndExchange(                                &
      & auxField          = me%scheme%auxField(iLevel),              &
      & calcAuxField      = me%scheme%calcAuxField,                  &
      & state             = me%scheme%state(iLevel)%val(:, now),     &
      & pdfData           = me%scheme%pdf(iLevel),                   &
      & nFields           = me%scheme%nFields,                       &
      & field             = me%scheme%field(:),                      &
      & globSrc           = me%scheme%globSrc,                       &
      & stencil           = me%scheme%layout%fStencil,               &
      & varSys            = me%scheme%varSys,                        &
      & derVarPos         = me%scheme%derVarPos,                     &
      & general           = me%params%general,                       &
      & phyConvFac        = me%params%physics%fac(iLevel),           &
      & iLevel            = iLevel,                                  &
      & minLevel          = me%geometry%tree%global%minLevel,        &
      & schemeHeader      = me%scheme%header,                        &
      & quantities        = me%scheme%layout%quantities              )

    if (iLevel < me%geometry%tree%global%maxLevel) then
      write(logUnit(10), "(A)") 'Interpolate and exchange auxField in ' &
        &                     //'ghostFromFiner'
      call mus_intpAuxFieldCoarserAndExchange(        &
        & intp        = me%scheme%intp,               &
        & tAuxField   = me%scheme%auxField(iLevel),   &
        & sAuxField   = me%scheme%auxField(iLevel+1), &
        & tLevelDesc  = me%scheme%levelDesc(iLevel),  &
        & stencil     = me%scheme%layout%fStencil,    &
        & iLevel      = iLevel,                       &
        & nAuxScalars = me%scheme%varSys%nAuxScalars, &
        & general     = me%params%general             )
    end if

    call tem_stopTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    ! --------------------------------------------------------------------------

    write(logUnit(10), "(A)") 'Update relaxparams'
    ! --------------------------------------------------------------------------
    ! Update parameters, relaxation time .etc
    call tem_startTimer( timerHandle =  mus_timerHandles%relax(iLevel) )
    call mus_update_relaxParams( scheme  = me%scheme,                        &
      &                          iLevel  = iLevel,                           &
      &                          tNow    = me%params%general%simControl%now, &
      &                          physics = me%params%physics,                &
      &                          lattice = me%params%lattice,                &
      &                          nBCs    = me%geometry%boundary%nBCtypes     )
    call tem_stopTimer( timerHandle =  mus_timerHandles%relax(iLevel) )
    ! --------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    write(logUnit(10), "(A)") 'Stream and collide'
    ! Compute current scheme of current level
    call tem_startTimer( timerHandle =  mus_timerHandles%compute(iLevel) )

!$omp parallel
    call me%scheme%compute(                               &
      &  fieldProp = me%scheme%field(:)%fieldProp,        &
      &  inState   = me%scheme%state(iLevel)%val(:,Now),  &
      &  outState  = me%scheme%state(iLevel)%val(:,Next), &
      &  auxField  = me%scheme%auxField(ilevel)%val(:),   &
      &  neigh     = me%scheme%pdf(iLevel)%neigh(:),      &
      &  nElems    = me%scheme%pdf(iLevel)%nSize,         &
      &  nSolve    = me%scheme%pdf(iLevel)%nElems_solve,  &
      &  level     = iLevel,                              &
      &  layout    = me%scheme%layout,                    &
      &  params    = me%params,                           &
      &  derVarPos = me%scheme%derVarPos,                 &
      &  varSys    = me%scheme%varSys                     )
!$omp end parallel

    call tem_stopTimer( timerHandle =  mus_timerHandles%compute(iLevel) )
    ! -------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    write(logUnit(10), "(A)") 'Apply source'
    call mus_apply_sourceTerms( field      = me%scheme%field(:),                &
      &                         nFields    = me%scheme%nFields,                 &
      &                         globSrc    = me%scheme%globSrc,                 &
      &                         pdf        = me%scheme%pdf(iLevel),             &
      &                         varSys     = me%scheme%varSys,                  &
      &                         iLevel     = iLevel,                            &
      &                         time       = me%params%general%simControl%now,  &
      &                         state      = me%scheme%state(iLevel)%val,       &
      &                         auxField   = me%scheme%auxField(iLevel)%val(:), &
      &                         derVarPos  = me%scheme%derVarPos(:),            &
      &                         phyConvFac = me%params%physics%fac(iLevel)      )

    ! --------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    write(logUnit(10), "(A)") 'Communicate fluids'
    ! Communicate the halo elements of each scheme on current level
    call tem_startTimer( timerHandle =  mus_timerHandles%comm(iLevel) )
    ! communicate halo elements for Next
    call me%params%general%commPattern%exchange_real(            &
      &  send         = me%scheme%levelDesc(iLevel)%sendbuffer,  &
      &  recv         = me%scheme%levelDesc(iLevel)%recvbuffer,  &
      &  state        = me%scheme%state(iLevel)%val(:,Next),     &
      &  message_flag = iLevel,                                  &
      &  comm         = me%params%general%proc%comm              )

    ! communicate turbulent viscosity, required for interpolation
    if (trim(me%scheme%header%kind) == 'fluid' .or. &
      & trim(me%scheme%header%kind) == 'fluid_incompressible') then
      if (me%scheme%field(1)%fieldProp%fluid%turbulence%active) then
        call me%params%general%commPattern%exchange_real(                &
          & recv         = me%scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%recvbuffer,          &
          & send         = me%scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%sendbuffer,          &
          & state        = me%scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%visc(:),             &
          & message_flag = iLevel+100,                                   &
          & comm         = me%params%general%proc%comm                   )
      end if
    end if
    call tem_stopTimer( timerHandle =  mus_timerHandles%comm(iLevel) )
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! communicate my finer ghost element from coarse level
    if ( iLevel > me%geometry%tree%global%minLevel ) then
      write(logUnit(10), "(A)") 'Communicate ghostFromCoarser'
      call tem_startTimer( timerHandle                                &
        &                  = mus_timerHandles%commFromCoarser(iLevel) )

      call me%params%general%commPattern%exchange_real(                &
        & send    = me%scheme%levelDesc(iLevel)%sendbufferFromCoarser, &
        & recv    = me%scheme%levelDesc(iLevel)%recvbufferFromCoarser, &
        & state   = me%scheme%state(iLevel)%val(:, Next),              &
        & message_flag = iLevel,                                       &
        & comm    = me%params%general%proc%comm                        )

      ! communicate turbulent viscosity, required for interpolation
      if (trim(me%scheme%header%kind) == 'fluid' .or. &
        & trim(me%scheme%header%kind) == 'fluid_incompressible') then

        if (me%scheme%field(1)%fieldProp%fluid%turbulence%active) then
          call me%params%general%commPattern%exchange_real(                  &
            & recv         = me%scheme%field(1)%fieldProp%fluid%turbulence   &
            &                      %dataOnLvl(iLevel)%recvBufferFromCoarser, &
            & send         = me%scheme%field(1)%fieldProp%fluid%turbulence   &
            &                      %dataOnLvl(iLevel)%sendBufferFromCoarser, &
            & state        = me%scheme%field(1)%fieldProp%fluid%turbulence   &
            &                      %dataOnLvl(iLevel)%visc(:),               &
            & message_flag = iLevel+200,                                     &
            & comm         = me%params%general%proc%comm                     )
        end if

      end if

      call tem_stopTimer( timerHandle                                 &
        &                 = mus_timerHandles%commFromCoarser(iLevel)  )
    end if
    ! --------------------------------------------------------------------------

    call stop_stageTimer()

    ! Interpolate ghost elements
    if( iLevel < me%geometry%tree%global%maxLevel ) then
      ! Fill my coarser element (L) from finer (L+1)
      call do_intpFinerAndExchange( me%scheme, me%params, iLevel )

      ! Interpolate the ghost elements on the finer level(L+1) with data provided
      ! from current level(L).
      call do_intpCoarserAndExchange( me%scheme, me%params, iLevel )
    end if ! if not on finest level
   ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    if( iLevel == me%geometry%tree%global%minLevel ) then
      iStage = 0
      running = .false.
    end if
    ! --------------------------------------------------------------------------

  end subroutine do_recursive_multiLevel
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Control routine for an optimized workflow with reduced functionality.
  !!
  !! No sources, no multilevel, no multiLevel.
  !! Use for benchmarking
  !!
  subroutine do_fast_singleLevel( me, iLevel )
    ! -------------------------------------------------------------------- !
    !> self control type
    !! dummy variable in this routine, required by interface
    class(mus_control_type) :: me
    !> Level counter variable
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: now, next
    ! -------------------------------------------------------------------- !

    me%curlvl = iLevel

    ! Update auxField dependent source fields before adding source term to state
    ! and auxField such that both auxField and apply_source uses same source.
    call mus_update_sourceVars( nFields    = me%scheme%nFields,              &
      &                         field      = me%scheme%field,                &
      &                         globSrc    = me%scheme%globSrc,              &
      &                         varSys     = me%scheme%varSys,               &
      &                         iLevel     = iLevel,                         &
      &                         auxField   = me%scheme%auxField(iLevel)%val, &
      &                         phyConvFac = me%params%physics%fac(iLevel),  &
      &                         derVarPos  = me%scheme%derVarPos             )

    ! -------------------------------------------------------------------------
    ! Increasing with the smallest time step (maxLevel)
    ! KM: time is advanced here since new time is required to update sources and BCs
    call tem_time_advance( me = me%params%general%simControl%now,   &
      &                    sim_dt = me%params%physics%dtLvl(iLevel ))

    ! --------------------------------------------------------------------------
    !set boundary for each field in current me%scheme
    call set_boundary( field       = me%scheme%field,                  &
      &                pdf         = me%scheme%pdf(iLevel),            &
      &                state       = me%scheme%state(iLevel)%val,      &
      &                levelDesc   = me%scheme%levelDesc(iLevel),      &
      &                tree        = me%geometry%tree,                 &
      &                iLevel      = iLevel,                           &
      &                nBCs        = me%geometry%boundary%nBCtypes,    &
      &                params      = me%params,                        &
      &                layout      = me%scheme%layout,                 &
      &                physics     = me%params%physics,                &
      &                varSys      = me%scheme%varSys,                 &
      &                mixture     = me%scheme%mixture,                &
      &                derVarPos   = me%scheme%derVarPos,              &
      &                globBC      = me%scheme%globBC                  )
    ! --------------------------------------------------------------------------

    ! swap double buffer index for current level
    call mus_swap_now_next( me%scheme%pdf( iLevel ) )
    now  = me%scheme%pdf(iLevel)%nNow
    next = me%scheme%pdf(iLevel)%nNext

    ! --------------------------------------------------------------------------
    ! Compute auxField from pre-collision state for fluid and ghostFromCoarser
    ! and exchange them if turbulence is active
    call tem_startTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    call mus_calcAuxFieldAndExchange(                            &
      & auxField          = me%scheme%auxField(iLevel),          &
      & calcAuxField      = me%scheme%calcAuxField,              &
      & state             = me%scheme%state(iLevel)%val(:, now), &
      & pdfData           = me%scheme%pdf(iLevel),               &
      & nFields           = me%scheme%nFields,                   &
      & field             = me%scheme%field(:),                  &
      & globSrc           = me%scheme%globSrc,                   &
      & stencil           = me%scheme%layout%fStencil,           &
      & varSys            = me%scheme%varSys,                    &
      & derVarPos         = me%scheme%derVarPos,                 &
      & general           = me%params%general,                   &
      & phyConvFac        = me%params%physics%fac(iLevel),       &
      & iLevel            = iLevel,                              &
      & minLevel          = me%geometry%tree%global%minLevel,    &
      & schemeHeader      = me%scheme%header,                    &
      & quantities        = me%scheme%layout%quantities          )
    call tem_stopTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    ! --------------------------------------------------------------------------

    if (me%particleGroup%nParticles > 0) then
      call me%check_particles()
    end if

    ! --------------------------------------------------------------------------
    ! Update parameters, relaxation time .etc
    call tem_startTimer( timerHandle =  mus_timerHandles%relax(iLevel) )
    call mus_update_relaxParams( scheme  = me%scheme,                        &
      &                          iLevel  = iLevel,                           &
      &                          tNow    = me%params%general%simControl%now, &
      &                          physics = me%params%physics,                &
      &                          lattice = me%params%lattice,                &
      &                          nBCs    = me%geometry%boundary%nBCtypes     )
    call tem_stopTimer( timerHandle =  mus_timerHandles%relax(iLevel) )
    ! --------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Compute current scheme of current level
    call tem_startTimer( timerHandle =  mus_timerHandles%compute(iLevel) )

!$omp parallel
    call me%scheme%compute(                                        &
      &           fieldProp = me%scheme%field(:)%fieldProp,        &
      &           inState   = me%scheme%state(iLevel)%val(:,Now),  &
      &           outState  = me%scheme%state(iLevel)%val(:,Next), &
      &           auxField  = me%scheme%auxField(ilevel)%val(:),   &
      &           neigh     = me%scheme%pdf(iLevel)%neigh(:),      &
      &           nElems    = me%scheme%pdf(iLevel)%nSize,         &
      &           nSolve    = me%scheme%pdf(iLevel)%nElems_solve,  &
      &           level     = iLevel,                              &
      &           layout    = me%scheme%layout,                    &
      &           params    = me%params,                           &
      &           derVarPos = me%scheme%derVarPos,                 &
      &           varSys    = me%scheme%varSys                     )
!$omp end parallel

    call tem_stopTimer( timerHandle =  mus_timerHandles%compute(iLevel) )
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    call mus_apply_sourceTerms( field      = me%scheme%field(:),                &
      &                         nFields    = me%scheme%nFields,                 &
      &                         globSrc    = me%scheme%globSrc,                 &
      &                         pdf        = me%scheme%pdf(iLevel),             &
      &                         varSys     = me%scheme%varSys,                  &
      &                         iLevel     = iLevel,                            &
      &                         time       = me%params%general%simControl%now,  &
      &                         state      = me%scheme%state(iLevel)%val,       &
      &                         auxField   = me%scheme%auxField(iLevel)%val(:), &
      &                         derVarPos  = me%scheme%derVarPos(:),            &
      &                         phyConvFac = me%params%physics%fac(iLevel)      )
    ! -------------------------------------------------------------------------

    if (me%particleGroup%nParticles > 0) then
      call me%advance_particles()
    end if
    ! Communicate the halo elements of each scheme on current level
    ! KM: Communicate post-collision before set_boundary because nonEq_expol
    ! BC depends on post-collision from neighbor at next time step
    call tem_startTimer( timerHandle =  mus_timerHandles%comm(iLevel) )
    call me%params%general%commPattern%exchange_real(             &
      &    send    = me%scheme%levelDesc(iLevel)%sendbuffer,      &
      &    recv    = me%scheme%levelDesc(iLevel)%recvbuffer,      &
      &    state   = me%scheme%state(iLevel)%val(:,Next),         &
      &    message_flag   = iLevel,                               &
      &    comm    = me%params%general%proc%comm                  )

    ! communicate turbulent viscosity, required for interpolation
    if (trim(me%scheme%header%kind) == 'fluid' .or. &
      & trim(me%scheme%header%kind) == 'fluid_incompressible') then
      if (me%scheme%field(1)%fieldProp%fluid%turbulence%active) then
        call me%params%general%commPattern%exchange_real(                &
          & recv         = me%scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%recvbuffer,          &
          & send         = me%scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%sendbuffer,          &
          & state        = me%scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%visc(:),             &
          & message_flag = iLevel+100,                                   &
          & comm         = me%params%general%proc%comm                   )
      end if
    end if
    call tem_stopTimer( timerHandle =  mus_timerHandles%comm(iLevel) )

    ! ... check if at least one of the IBMs is active
    if ( me%geometry%globIBM%nIBMs > 0 ) then
      call mus_buildBuffIBM(                                 &
        &       me          = me%geometry%globIBM%IBM,       &
        &       commPattern = me%params%general%commPattern, &
        &       globTree    = me%geometry%tree,              &
        &       params      = me%params,                     &
        &       layout      = me%scheme%layout,              &
        &       levelDesc   = me%scheme%levelDesc(iLevel),   &
        &       iLevel      = iLevel                         )
    end if

    ! update the immersed boundaries if available
    ! ... and over the schemes
    ! ... check if at least one of the IBMs is active
    if( me%geometry%globIBM%nIBMs > 0 )then
      call mus_inamuro_IBM(                                       &
        &      me          = me%geometry%globIBM%IBM,             &
        &      commPattern = me%params%general%commPattern,       &
        &      globTree    = me%geometry%tree,                    &
        &      general     = me%params%general,                   &
        &      pdf         = me%scheme%pdf(iLevel),               &
        &      layout      = me%scheme%layout,                    &
        &      levelDesc   = me%scheme%levelDesc(iLevel),         &
        &      globSys     = me%scheme%varSys,                    &
        &      stateVarMap = me%scheme%stateVarMap%varPos%val(:), &
        &      convFac     = me%params%physics%fac(iLevel),       &
        &      iField      = 1,                                   &
        &      state       = me%scheme%state(iLevel)%val,         &
        &      iLevel      = iLevel                               )
    end if
    ! -------------------------------------------------------------------------

  end subroutine do_fast_singleLevel
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine do_benchmark(me,  iLevel)
    ! -------------------------------------------------------------------- !
    !> self control type
    !! dummy variable in this routine, required by interface
    class(mus_control_type) :: me
    !> Level counter variable
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: now, next
    ! -------------------------------------------------------------------- !

    ! Update auxField dependent source fields before adding source term to state
    ! and auxField such that both auxField and apply_source uses same source.
    call mus_update_sourceVars( nFields    = me%scheme%nFields,              &
      &                         field      = me%scheme%field,                &
      &                         globSrc    = me%scheme%globSrc,              &
      &                         varSys     = me%scheme%varSys,               &
      &                         iLevel     = iLevel,                         &
      &                         auxField   = me%scheme%auxField(iLevel)%val, &
      &                         phyConvFac = me%params%physics%fac(iLevel),  &
      &                         derVarPos  = me%scheme%derVarPos             )

    ! Increasing with the smallest time step (maxLevel)
    call tem_time_advance( me = me%params%general%simControl%now,   &
      &                    sim_dt = me%params%physics%dtLvl(iLevel ))

    ! --------------------------------------------------------------------------
    !set boundary for each field in current me%scheme
    call set_boundary( field       = me%scheme%field,               &
      &                pdf         = me%scheme%pdf(iLevel),         &
      &                state       = me%scheme%state(iLevel)%val,   &
      &                levelDesc   = me%scheme%levelDesc(iLevel),   &
      &                tree        = me%geometry%tree,              &
      &                iLevel      = iLevel,                        &
      &                nBCs        = me%geometry%boundary%nBCtypes, &
      &                params      = me%params,                     &
      &                layout      = me%scheme%layout,              &
      &                physics     = me%params%physics,             &
      &                varSys      = me%scheme%varSys,              &
      &                mixture     = me%scheme%mixture,             &
      &                derVarPos   = me%scheme%derVarPos,           &
      &                globBC      = me%scheme%globBC               )
    ! --------------------------------------------------------------------------

    ! swap double buffer index for current level
    call mus_swap_now_next( me%scheme%pdf( iLevel ) )
    now  = me%scheme%pdf(iLevel)%nNow
    next = me%scheme%pdf(iLevel)%nNext

    ! --------------------------------------------------------------------------
    ! Compute auxField from pre-collision state for fluid and ghostFromCoarser
    ! and exchange them if turbulence is active
    call tem_startTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    call mus_calcAuxFieldAndExchange(                                &
      & auxField          = me%scheme%auxField(iLevel),              &
      & calcAuxField      = me%scheme%calcAuxField,                  &
      & state             = me%scheme%state(iLevel)%val(:, now),     &
      & pdfData           = me%scheme%pdf(iLevel),                   &
      & nFields           = me%scheme%nFields,                       &
      & field             = me%scheme%field(:),                      &
      & globSrc           = me%scheme%globSrc,                       &
      & stencil           = me%scheme%layout%fStencil,               &
      & varSys            = me%scheme%varSys,                        &
      & derVarPos         = me%scheme%derVarPos,                     &
      & general           = me%params%general,                       &
      & phyConvFac        = me%params%physics%fac(iLevel),           &
      & iLevel            = iLevel,                                  &
      & minLevel          = me%geometry%tree%global%minLevel,        &
      & schemeHeader      = me%scheme%header,                        &
      & quantities        = me%scheme%layout%quantities              )
    call tem_stopTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! Update parameters, relaxation time .etc
    call mus_update_relaxParams( scheme  = me%scheme,                        &
      &                          iLevel  = iLevel,                           &
      &                          tNow    = me%params%general%simControl%now, &
      &                          physics = me%params%physics,                &
      &                          lattice = me%params%lattice,                &
      &                          nBCs    = me%geometry%boundary%nBCtypes     )
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! Compute current scheme of current level
    call tem_startTimer( timerHandle = mus_timerHandles%compute(iLevel) )

!$omp parallel
    call me%scheme%compute(                                        &
      &           fieldProp = me%scheme%field(:)%fieldProp,        &
      &           inState   = me%scheme%state(iLevel)%val(:,Now),  &
      &           outState  = me%scheme%state(iLevel)%val(:,Next), &
      &           auxField  = me%scheme%auxField(ilevel)%val(:),   &
      &           neigh     = me%scheme%pdf(iLevel)%neigh(:),      &
      &           nElems    = me%scheme%pdf(iLevel)%nSize,         &
      &           nSolve    = me%scheme%pdf(iLevel)%nElems_solve,  &
      &           level     = iLevel,                              &
      &           layout    = me%scheme%layout,                    &
      &           params    = me%params,                           &
      &           derVarPos = me%scheme%derVarPos,                 &
      &           varSys    = me%scheme%varSys                     )
!$omp end parallel

    call tem_stopTimer( timerHandle = mus_timerHandles%compute(iLevel) )
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    call mus_apply_sourceTerms( field      = me%scheme%field(:),                &
      &                         nFields    = me%scheme%nFields,                 &
      &                         globSrc    = me%scheme%globSrc,                 &
      &                         pdf        = me%scheme%pdf(iLevel),             &
      &                         varSys     = me%scheme%varSys,                  &
      &                         iLevel     = iLevel,                            &
      &                         time       = me%params%general%simControl%now,  &
      &                         state      = me%scheme%state(iLevel)%val,       &
      &                         auxField   = me%scheme%auxField(iLevel)%val(:), &
      &                         derVarPos  = me%scheme%derVarPos(:),            &
      &                         phyConvFac = me%params%physics%fac(iLevel)      )
    ! --------------------------------------------------------------------------

    ! Communicate the halo elements of each me%scheme on current level
    call tem_startTimer( timerHandle = mus_timerHandles%comm(iLevel) )
    call me%params%general%commPattern%exchange_real(                  &
      &         send    = me%scheme%levelDesc(iLevel)%sendbuffer,      &
      &         recv    = me%scheme%levelDesc(iLevel)%recvbuffer,      &
      &         state   = me%scheme%state(iLevel)%val(:,Next),         &
      &         message_flag   = iLevel,                               &
      &         comm    = me%params%general%proc%comm                  )

    ! communicate turbulent viscosity, required for interpolation
    if (trim(me%scheme%header%kind) == 'fluid' .or. &
      & trim(me%scheme%header%kind) == 'fluid_incompressible') then
      if (me%scheme%field(1)%fieldProp%fluid%turbulence%active) then
        call me%params%general%commPattern%exchange_real(                &
          & recv         = me%scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%recvbuffer,          &
          & send         = me%scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%sendbuffer,          &
          & state        = me%scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%visc(:),             &
          & message_flag = iLevel+100,                                   &
          & comm         = me%params%general%proc%comm                   )
      end if
    end if
    call tem_stopTimer( timerHandle = mus_timerHandles%comm(iLevel) )
   ! -------------------------------------------------------------------------

  end subroutine do_benchmark
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine does:
  !! 1. interpolate my coarse ghost element (iLevel) from finer level (iLevel+1)
  !! 2. exchange the data of my coarse ghost elements between process
  subroutine do_intpFinerAndExchange( scheme, params, iLevel )
    ! -------------------------------------------------------------------- !
    !> containers for the different schemes
    type( mus_scheme_type ), intent(inout), target  :: scheme
    !> global parameters
    type( mus_param_type ), intent(inout)  :: params
    integer, intent(in) :: iLevel     !< Level counter variable
    ! -------------------------------------------------------------------- !
    integer :: thisLevelNext, nextLevelNext
    ! -------------------------------------------------------------------- !

    ! Interpolate my coarse ghost elements from the finer ones
    thisLevelNext = scheme%pdf(iLevel)%nNext
    nextLevelNext = scheme%pdf(iLevel+1)%nNext

    call start_stageTimer()
    ! write(logUnit(5), "(A,I0)") 'intpFromFiner on level: ', iLevel
    ! -------------------------------------------------------------------------
    call tem_startTimer( timerHandle =  mus_timerHandles%intpFromFiner(iLevel) )
    ! -------------------------------------------------------------------------
    write(logUnit(10), "(A)") 'Interpolate ghostFromFiner'
!$omp parallel
    call scheme%intp%fillMineFromFiner%do_intp(                         &
      &     fieldProp   = scheme%field(:)%fieldProp,                    &
      &     sState      = scheme%state(iLevel+1)%val(:,nextLevelNext),  &
      &     sNeigh      = scheme%pdf(iLevel+1)%Neigh,                   &
      &     snSize      = scheme%pdf(iLevel+1)%nSize,                   &
      &     sAuxField   = scheme%auxField(iLevel+1)%val(:),             &
      &     tnSize      = scheme%pdf(iLevel)%nSize,                     &
      &     tState      = scheme%state(iLevel)%val(:,thisLevelNext),    &
      &     tNeigh      = scheme%pdf(iLevel)%Neigh,                     &
      &     tLevelDesc  = scheme%levelDesc(iLevel),                     &
      &     level       = iLevel,                                       &
      &     nTargets    = scheme%levelDesc(iLevel)%intpFromFiner%nVals, &
      &     targetList  = scheme%levelDesc(iLevel)%intpFromFiner%val,   &
      &     layout      = scheme%layout,                                &
      &     physics     = params%physics,                               &
      &     varSys      = scheme%varSys,                                &
      &     derVarPos   = scheme%derVarPos(:),                          &
      &     time        = params%general%simControl%now                 )
!$omp end parallel
    call tem_stopTimer( timerHandle = mus_timerHandles%intpFromFiner(iLevel) )
    ! -------------------------------------------------------------------------

    ! if ( main_debug%checkEachAlgorithmicStep ) then
    !   buffer=' after fill Mine FromFiner       '
    !   call dump_debug_info( buffer, scheme, params, iLevel, 1, &
    !     &                   pdf=scheme%pdf(iLevel)             )
    ! end if
    call stop_stageTimer()
    write(logUnit(10), "(A)") 'Communicate ghostFromFiner'
    ! -------------------------------------------------------------------------
    ! Exchange the coarse ghost elements
    call tem_startTimer( timerHandle = mus_timerHandles%commFromFiner(iLevel) )
    call params%general%commPattern%exchange_real(                      &
      &    send         = scheme%levelDesc(iLevel)%sendbufferFromFiner, &
      &    recv         = scheme%levelDesc(iLevel)%recvbufferFromFiner, &
      &    state        = scheme%state(iLevel)%val(:,thisLevelNext),    &
      &    message_flag = iLevel,                                       &
      &    comm         = params%general%proc%comm                      )

    if (trim(scheme%header%kind) == 'fluid' .or. &
      & trim(scheme%header%kind) == 'fluid_incompressible') then
      if (scheme%field(1)%fieldProp%fluid%turbulence%active) then
        ! communicate turbulent viscosity, required for interpolation
        call params%general%commPattern%exchange_real(                    &
          &  recv         = scheme%field(1)%fieldProp%fluid%turbulence    &
          &                       %dataOnLvl(iLevel)%recvBufferFromFiner, &
          &  send         = scheme%field(1)%fieldProp%fluid%turbulence    &
          &                       %dataOnLvl(iLevel)%sendBufferFromFiner, &
          &  state        = scheme%field(1)%fieldProp%fluid%turbulence    &
          &                       %dataOnLvl(iLevel)%visc(:),             &
          &  message_flag = iLevel+200,                                   &
          &  comm         = params%general%proc%comm                      )
      end if
    end if

    call tem_stopTimer( timerHandle = mus_timerHandles%commFromFiner(iLevel) )
    ! -------------------------------------------------------------------------

    ! if ( main_debug%checkEachAlgorithmicStep ) then
    !   buffer = ' after exch Mine FromFiner'
    !   call dump_debug_info( buffer, scheme, params, iLevel, 1, &
    !     &                   pdf=scheme%pdf(iLevel)             )
    ! end if

  end subroutine do_intpFinerandExchange
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine utilizes fluid elements on my level (L) to fill finer
  !! ghost elements on next level (L+1).
  !! Then it exchanges the datas of finer ghost elements (L+1) between process.
  subroutine do_intpCoarserAndExchange( scheme, params, iLevel )
    ! -------------------------------------------------------------------- !
    !> containers for the different schemes
    type( mus_scheme_type ), target, intent(inout)  :: scheme
    !> global parameters
    type( mus_param_type ), intent(in) :: params
    !> Level counter variable
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: thisLevelNext, nextLevelNext
    integer :: iOrder
    ! -------------------------------------------------------------------- !

    thisLevelNext = scheme%pdf(iLevel)%nNext
    nextLevelNext = scheme%pdf(iLevel+1)%nNext

    call start_stageTimer()
    ! write(logUnit(5), "(A,I0)") 'intpFromCoarser on level: ', iLevel+1
    ! Fill finer ghost elements (L+1) by low order interpolation
    call tem_startTimer(                                            &
      &    timerHandle = mus_timerHandles%intpFromCoarser(iLevel+1) )

    ! Fill finer ghost elements (L+1) by starting from lower order
    ! to higher order interpolation
    write(logUnit(10), "(A)") 'Interpolate ghostFromCoarser iL+1'
!$omp parallel
    do iOrder = 0, scheme%intp%config%order
      call scheme%intp%fillFinerFromMe(iOrder)%do_intp(              &
        & fieldProp   = scheme%field(:)%fieldProp,                   &
        & sState      = scheme%state(iLevel)%val(:,thisLevelNext),   &
        & sNeigh      = scheme%pdf(iLevel)%Neigh,                    &
        & snSize      = scheme%pdf(iLevel)%nSize,                    &
        & sAuxField   = scheme%auxField(iLevel)%val(:),              &
        & tnSize      = scheme%pdf(iLevel+1)%nSize,                  &
        & tState      = scheme%state(iLevel+1)%val(:,nextLevelNext), &
        & tNeigh      = scheme%pdf(iLevel+1)%Neigh,                  &
        & tLevelDesc  = scheme%levelDesc(iLevel+1),                  &
        & level       = iLevel,                                      &
        & nTargets    = scheme%levelDesc(iLevel+1)                   &
        &                     %intpFromCoarser(iOrder)%nVals,        &
        & targetList  = scheme%levelDesc(iLevel+1)                   &
        &                     %intpFromCoarser(iOrder)%Val,          &
        & layout      = scheme%layout,                               &
        & physics     = params%physics,                              &
        & varSys      = scheme%varSys,                               &
        & derVarPos   = scheme%derVarPos(:),                         &
        & time        = params%general%simControl%now                )
     end do
!$omp end parallel
      call tem_stopTimer( timerHandle                                  &
        &                 = mus_timerHandles%intpFromCoarser(iLevel+1) )

    ! Debug output
    ! if ( main_debug%checkEachAlgorithmicStep ) then
    !   buffer = ' after fillFinerFromMe'
    !   call dump_debug_info( buffer, scheme, params, iLevel, 1, &
    !     &                   pdf = scheme%pdf(iLevel))
    ! end if

    call stop_stageTimer()
    ! Exchange the fine ghost elements
    write(logUnit(10), "(A)") 'Communicate ghostFromCoarser iL+1'
    call tem_startTimer(                                            &
      &    timerHandle = mus_timerHandles%commFromCoarser(iLevel+1) )
    call params%general%commPattern%exchange_real(                   &
      &    send = scheme%levelDesc(iLevel+1)%SendBufferFromCoarser,  &
      &    recv = scheme%levelDesc(iLevel+1)%RecvBufferFromCoarser,  &
      &    state   = scheme%state(iLevel+1)%val(:,nextLevelNext),    &
      &    message_flag   = iLevel,                                  &
      &    comm    = params%general%proc%comm                             )

    if (trim(scheme%header%kind) == 'fluid' .or. &
      & trim(scheme%header%kind) == 'fluid_incompressible') then
      if (scheme%field(1)%fieldProp%fluid%turbulence%active) then
        ! communicate turbulent viscosity, required for interpolation
        call params%general%commPattern%exchange_real(                       &
          & recv         = scheme%field(1)%fieldProp%fluid%turbulence        &
          &                      %dataOnLvl(iLevel+1)%recvBufferFromCoarser, &
          & send         = scheme%field(1)%fieldProp%fluid%turbulence        &
          &                      %dataOnLvl(iLevel+1)%sendBufferFromCoarser, &
          & state        = scheme%field(1)%fieldProp%fluid%turbulence        &
          &                      %dataOnLvl(iLevel+1)%visc(:),               &
          & message_flag = iLevel+200,                                       &
          & comm         = params%general%proc%comm                          )
      end if
    end if
    call tem_stopTimer(                                             &
      &    timerHandle = mus_timerHandles%commFromCoarser(iLevel+1) )

    ! Debug output
    ! if ( main_debug%checkEachAlgorithmicStep ) then
    !   buffer = '  after exchFinerFromMe'
    !   call dump_debug_info( buffer, scheme, params, iLevel, 1, &
    !     &                   pdf = scheme%pdf(iLevel)           )
    ! end if

  end subroutine do_intpCoarserandExchange
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine check_particles_MEM(me)
    ! -------------------------------------------------------------------- !
    !> self control type
    class(mus_control_type) :: me
    ! -------------------------------------------------------------------- !
    call check_and_create_new_particles_MEM(                      &
      & particle_creator = particle_creator,                      &
      & iter             = me%params%general%simControl%now%iter, &
      & particleGroup    = me%particleGroup,                      &
      & scheme           = me%scheme,                             &
      & geometry         = me%geometry,                           &
      & params           = me%params,                             &
      & myRank           = me%params%general%proc%rank            )
  end subroutine check_particles_MEM

  subroutine advance_particles_MEM(me)
    ! -------------------------------------------------------------------- !
    !> self control type
    class(mus_control_type) :: me
    ! -------------------------------------------------------------------- !
    ! Update particles
    call me%particleGroup%moveParticles( &
      &    scheme   = me%scheme,         &
      &    geometry = me%geometry,       &
      &    params   = me%params          )

    call mus_particles_logdata_MEM( particleGroup = me%particleGroup, &
      &                             params        = me%params         )

    call me%particleGroup%mapParticles( &
      &    scheme   = me%scheme,        &
      &    geometry = me%geometry,      &
      &    params   = me%params         )

    call me%particleGroup%applyHydrodynamicForces( &
      &    scheme   = me%scheme,                   &
      &    geometry = me%geometry,                 &
      &    params   = me%params                    )

    call me%particleGroup%transferMomentumToFluid( &
      &    scheme   = me%scheme,                   &
      &    geometry = me%geometry,                 &
      &    params   = me%params                    )

  end subroutine advance_particles_MEM

  subroutine check_particles_DPS(me)
    ! -------------------------------------------------------------------- !
    !> self control type
    class(mus_control_type) :: me
    ! -------------------------------------------------------------------- !
    ! Check if new particles should be created at this time step
    call check_and_create_new_particles_DPS(                         &
      &    particle_creator = particle_creator,                      &
      &    iter             = me%params%general%simControl%now%iter, &
      &    particleGroup    = me%particleGroup,                      &
      &    scheme           = me%scheme,                             &
      &    geometry         = me%geometry,                           &
      &    params           = me%params,                             &
      &    myRank           = me%params%general%proc%rank            )

    if (me%DPS_do_VolFract) then
      ! Update the local fluid volume fraction in auxField
      call mus_particles_updateFluidVolumeFraction(                & 
        &    particleGroup = me%particleGroup,                     &
        &    scheme        = me%scheme,                            &
        &    geometry      = me%geometry,                          &
        &    params        = me%params,                            &
        &    nElems        = me%scheme%pdf(me%curlvl)%nElems_local )
    end if

    ! Update positions and velocities of particles
    call me%particleGroup%moveParticles( scheme        = me%scheme,   &
      &                                  geometry      = me%geometry, &
      &                                  params        = me%params    )   
                                  
    call mus_particles_logdata_DPS( particleGroup = me%particleGroup, &
      &                             params        = me%params         )
  end subroutine check_particles_DPS

  subroutine advance_particles_DPS(me)
    ! -------------------------------------------------------------------- !
    !> self control type
    class(mus_control_type) :: me
    ! -------------------------------------------------------------------- !
    if (me%DPS_do_advance) then
      call me%particleGroup%transferMomentumToFluid( &
        &    scheme   = me%scheme,                   &
        &    geometry = me%geometry,                 &
        &    params   = me%params                    )
    end if
  end subroutine advance_particles_DPS
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine start_stageTimer( )

    if ( .not. running ) then
      iStage = mod( iStage, nStages ) + 1
      call tem_startTimer( timerHandle =  mus_timerHandles%stage(iStage) )
      running = .true.
    end if

  end subroutine start_stageTimer
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine stop_stageTimer( )

    if ( running ) then
      call tem_stopTimer( timerHandle =  mus_timerHandles%stage(iStage) )
      running = .false.
    end if

  end subroutine stop_stageTimer
  ! ------------------------------------------------------------------------ !

end module mus_control_module
! ***************************************************************************** !
