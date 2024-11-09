! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2011-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2012, 2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
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

  implicit none

  private

  public :: mus_control_type
  public :: mus_init_control

  !> Datatype containing mapping of control routines to function pointers
  type mus_control_type
    procedure( computation ), pointer :: do_computation => null()
  end type mus_control_type

  abstract interface
    !> Interface describes the main control routine which does computation
    !! set boundary and check flow status
    subroutine computation( me, scheme, geometry, params, iLevel)
      import :: mus_scheme_type, mus_geom_type, mus_param_type, &
        &       mus_control_type
      !> self control type
      class( mus_control_type ) :: me
      !> container for the scheme
      type( mus_scheme_type ), intent(inout) :: scheme
      !> geometry infomation
      type( mus_geom_type ), intent(inout) :: geometry
      !> global parameters
      type( mus_param_type ), intent(inout) :: params
      !> Level counter variable
      integer, intent(in) :: iLevel
    end subroutine computation
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
  subroutine mus_init_control( controlRoutine, me, minLevel, maxLevel )
    ! -------------------------------------------------------------------- !
    character(len=labelLen), intent(in) :: controlRoutine
    !> contains function pointer to point control routine
    type( mus_control_type ), intent(out) :: me
    integer, intent(in) :: minLevel, maxLevel
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
  recursive subroutine do_recursive_multiLevel( me, scheme, geometry, params, &
    &                                           iLevel)
    ! -------------------------------------------------------------------- !
    !> self control type
    class( mus_control_type ) :: me
    !> container for the scheme
    type( mus_scheme_type ), intent(inout) :: scheme
    !> geometry infomation
    type( mus_geom_type ), intent(inout)      :: geometry
    !> global parameters
    type( mus_param_type ), intent(inout)  :: params
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
    call mus_update_sourceVars( nFields    = scheme%nFields,              &
      &                         field      = scheme%field,                &
      &                         globSrc    = scheme%globSrc,              &
      &                         varSys     = scheme%varSys,               &
      &                         iLevel     = iLevel,                      &
      &                         auxField   = scheme%auxField(iLevel)%val, &
      &                         phyConvFac = params%physics%fac(iLevel),  &
      &                         derVarPos  = scheme%derVarPos             )

    ! when not on finest level, go to next level
    if( iLevel < geometry%tree%global%maxLevel ) then
      ! Perform the number of nested time steps on the finer level L+1
      ! according to scaling type :
      !   diffusive: nNesting = 4
      !   acoustic:  nNesting = 2
      do iNestingLoop = 1, params%nNesting
        write(logUnit(10), "(A,I0)") 'Nesting loop ', iNestingloop
        call me%do_computation( scheme, geometry, params, iLevel+1 )
      end do
    end if

    write(logUnit(10), "(A,I0)") 'Compute on level: ', iLevel
    call start_stageTimer()

    ! update the time counters. MH: checked. Please dont move
    ! Increasing with the smallest time step (maxLevel)
    ! KM: time is advanced here since new time is required to update sources.
    if( iLevel == geometry%tree%global%maxLevel ) then
      write(logUnit(10), "(A)") 'Advance time t+dt_maxLevel'
      call tem_time_advance( me     = params%general%simControl%now, &
        &                    sim_dt = params%physics%dtLvl( geometry &
        &                                   %tree%global%maxLevel )  )
    endif

    write(logUnit(10), "(A)") 'Set boundary condition'
    ! set boundary for each field in current scheme
    call set_boundary( field       = scheme%field,                  &
      &                pdf         = scheme%pdf(iLevel),            &
      &                state       = scheme%state(iLevel)%val,      &
      &                levelDesc   = scheme%levelDesc(iLevel),      &
      &                tree        = geometry%tree,                 &
      &                iLevel      = iLevel,                        &
      &                nBCs        = geometry%boundary%nBCtypes,    &
      &                params      = params,                        &
      &                layout      = scheme%layout,                 &
      &                physics     = params%physics,                &
      &                varSys      = scheme%varSys,                 &
      &                mixture     = scheme%mixture,                &
      &                derVarPos   = scheme%derVarPos,              &
      &                globBC      = scheme%globBC                  )
    ! -------------------------------------------------------------------------


    write(logUnit(10), "(A)") 'Swap now and next'
    ! swap double buffer index for current level
    call mus_swap_now_next( scheme%pdf( iLevel ) )
    now  = scheme%pdf(iLevel)%nNow
    next = scheme%pdf(iLevel)%nNext

    ! --------------------------------------------------------------------------
    ! Compute auxField from pre-collision state for fluid and ghostFromCoarser
    ! and exchange them if turbulence is active
    call tem_startTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    write(logUnit(10), "(A)") 'Calculate auxField'
    call mus_calcAuxFieldAndExchange(                             &
      & auxField          = scheme%auxField(iLevel),              &
      & calcAuxField      = scheme%calcAuxField,                  &
      & state             = scheme%state(iLevel)%val(:, now),     &
      & pdfData           = scheme%pdf(iLevel),                   &
      & nFields           = scheme%nFields,                       &
      & field             = scheme%field(:),                      &
      & globSrc           = scheme%globSrc,                       &
      & stencil           = scheme%layout%fStencil,               &
      & varSys            = scheme%varSys,                        &
      & derVarPos         = scheme%derVarPos,                     &
      & general           = params%general,                       &
      & phyConvFac        = params%physics%fac(iLevel),           &
      & iLevel            = iLevel,                               &
      & minLevel          = geometry%tree%global%minLevel,        &
      & schemeHeader      = scheme%header,                        &
      & quantities        = scheme%layout%quantities              )

    if (iLevel < geometry%tree%global%maxLevel) then
      write(logUnit(10), "(A)") 'Interpolate and exchange auxField in ' &
        &                     //'ghostFromFiner'
      call mus_intpAuxFieldCoarserAndExchange(     &
        & intp        = scheme%intp,               &
        & tAuxField   = scheme%auxField(iLevel),   &
        & sAuxField   = scheme%auxField(iLevel+1), &
        & tLevelDesc  = scheme%levelDesc(iLevel),  &
        & stencil     = scheme%layout%fStencil,    &
        & iLevel      = iLevel,                    &
        & nAuxScalars = scheme%varSys%nAuxScalars, &
        & general     = params%general             )
    end if

    call tem_stopTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    ! --------------------------------------------------------------------------

    write(logUnit(10), "(A)") 'Update relaxparams'
    ! --------------------------------------------------------------------------
    ! Update parameters, relaxation time .etc
    call tem_startTimer( timerHandle =  mus_timerHandles%relax(iLevel) )
    call mus_update_relaxParams( scheme  = scheme,                        &
      &                          iLevel  = iLevel,                        &
      &                          tNow    = params%general%simControl%now, &
      &                          physics = params%physics,                &
      &                          lattice = params%lattice,                &
      &                          nBCs    = geometry%boundary%nBCtypes     )
    call tem_stopTimer( timerHandle =  mus_timerHandles%relax(iLevel) )
    ! --------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    write(logUnit(10), "(A)") 'Stream and collide'
    ! Compute current scheme of current level
    call tem_startTimer( timerHandle =  mus_timerHandles%compute(iLevel) )

!$omp parallel
    call scheme%compute(                               &
      &  fieldProp = scheme%field(:)%fieldProp,        &
      &  inState   = scheme%state(iLevel)%val(:,Now),  &
      &  outState  = scheme%state(iLevel)%val(:,Next), &
      &  auxField  = scheme%auxField(ilevel)%val(:),   &
      &  neigh     = scheme%pdf(iLevel)%neigh(:),      &
      &  nElems    = scheme%pdf(iLevel)%nSize,         &
      &  nSolve    = scheme%pdf(iLevel)%nElems_solve,  &
      &  level     = iLevel,                           &
      &  layout    = scheme%layout,                    &
      &  params    = params,                           &
      &  derVarPos = scheme%derVarPos,                 &
      &  varSys    = scheme%varSys                     )
!$omp end parallel

    call tem_stopTimer( timerHandle =  mus_timerHandles%compute(iLevel) )
    ! -------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    write(logUnit(10), "(A)") 'Apply source'
    call mus_apply_sourceTerms( field      = scheme%field(:),                &
      &                         nFields    = scheme%nFields,                 &
      &                         globSrc    = scheme%globSrc,                 &
      &                         pdf        = scheme%pdf(iLevel),             &
      &                         varSys     = scheme%varSys,                  &
      &                         iLevel     = iLevel,                         &
      &                         time       = params%general%simControl%now,  &
      &                         state      = scheme%state(iLevel)%val,       &
      &                         auxField   = scheme%auxField(iLevel)%val(:), &
      &                         derVarPos  = scheme%derVarPos(:),            &
      &                         phyConvFac = params%physics%fac(iLevel)      )

    ! --------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    write(logUnit(10), "(A)") 'Communicate fluids'
    ! Communicate the halo elements of each scheme on current level
    call tem_startTimer( timerHandle =  mus_timerHandles%comm(iLevel) )
    ! communicate halo elements for Next
    call params%general%commPattern%exchange_real(            &
      &  send         = scheme%levelDesc(iLevel)%sendbuffer,  &
      &  recv         = scheme%levelDesc(iLevel)%recvbuffer,  &
      &  state        = scheme%state(iLevel)%val(:,Next),     &
      &  message_flag = iLevel,                               &
      &  comm         = params%general%proc%comm              )

    ! communicate turbulent viscosity, required for interpolation
    if (trim(scheme%header%kind) == 'fluid' .or. &
      & trim(scheme%header%kind) == 'fluid_incompressible') then
      if (scheme%field(1)%fieldProp%fluid%turbulence%active) then
        call params%general%commPattern%exchange_real(                &
          & recv         = scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%recvbuffer,       &
          & send         = scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%sendbuffer,       &
          & state        = scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%visc(:),          &
          & message_flag = iLevel+100,                               &
          & comm         = params%general%proc%comm                  )
      end if
    end if
    call tem_stopTimer( timerHandle =  mus_timerHandles%comm(iLevel) )
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! communicate my finer ghost element from coarse level
    if ( iLevel > geometry%tree%global%minLevel ) then
      write(logUnit(10), "(A)") 'Communicate ghostFromCoarser'
      call tem_startTimer( timerHandle                                &
        &                  = mus_timerHandles%commFromCoarser(iLevel) )

      call params%general%commPattern%exchange_real(                  &
        & send    = scheme%levelDesc(iLevel)%sendbufferFromCoarser,   &
        & recv    = scheme%levelDesc(iLevel)%recvbufferFromCoarser,   &
        & state   = scheme%state(iLevel)%val(:, Next),                &
        & message_flag = iLevel,                                      &
        & comm    = params%general%proc%comm             )

      ! communicate turbulent viscosity, required for interpolation
      if (trim(scheme%header%kind) == 'fluid' .or. &
        & trim(scheme%header%kind) == 'fluid_incompressible') then

        if (scheme%field(1)%fieldProp%fluid%turbulence%active) then
          call params%general%commPattern%exchange_real(                     &
            & recv         = scheme%field(1)%fieldProp%fluid%turbulence      &
            &                      %dataOnLvl(iLevel)%recvBufferFromCoarser, &
            & send         = scheme%field(1)%fieldProp%fluid%turbulence      &
            &                      %dataOnLvl(iLevel)%sendBufferFromCoarser, &
            & state        = scheme%field(1)%fieldProp%fluid%turbulence      &
            &                      %dataOnLvl(iLevel)%visc(:),               &
            & message_flag = iLevel+200,                                     &
            & comm         = params%general%proc%comm                        )
        end if

      end if

      call tem_stopTimer( timerHandle                                 &
        &                 = mus_timerHandles%commFromCoarser(iLevel)  )
    end if
    ! --------------------------------------------------------------------------

    call stop_stageTimer()

    ! Interpolate ghost elements
    if( iLevel < geometry%tree%global%maxLevel ) then
      ! Fill my coarser element (L) from finer (L+1)
      call do_intpFinerAndExchange( scheme, params, iLevel )

      ! Interpolate the ghost elements on the finer level(L+1) with data provided
      ! from current level(L).
      call do_intpCoarserAndExchange( scheme, params, iLevel )
    end if ! if not on finest level
   ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    if( iLevel == geometry%tree%global%minLevel ) then
      iStage = 0
      running = .false.
    end if
    ! --------------------------------------------------------------------------

  end subroutine do_recursive_multiLevel
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
  !> Control routine for an optimized workflow with reduced functionality.
  !!
  !! No sources, no multilevel, no multiLevel.
  !! Use for benchmarking
  !!
  subroutine do_fast_singleLevel( me, scheme, geometry, params, iLevel )
    ! -------------------------------------------------------------------- !
    !> self control type
    !! dummy variable in this routine, required by interface
    class( mus_control_type ) :: me
    !> container for the scheme
    type( mus_scheme_type ), intent(inout)  :: scheme
    !> geometry infomation
    type( mus_geom_type ), intent(inout)    :: geometry
    !> global parameters
    type( mus_param_type ), intent(inout)   :: params
    !> Level counter variable
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: now, next
    ! -------------------------------------------------------------------- !

    ! Update auxField dependent source fields before adding source term to state
    ! and auxField such that both auxField and apply_source uses same source.
    call mus_update_sourceVars( nFields    = scheme%nFields,              &
      &                         field      = scheme%field,                &
      &                         globSrc    = scheme%globSrc,              &
      &                         varSys     = scheme%varSys,               &
      &                         iLevel     = iLevel,                      &
      &                         auxField   = scheme%auxField(iLevel)%val, &
      &                         phyConvFac = params%physics%fac(iLevel),  &
      &                         derVarPos  = scheme%derVarPos             )

    ! -------------------------------------------------------------------------
    ! Increasing with the smallest time step (maxLevel)
    ! KM: time is advanced here since new time is required to update sources and BCs
    call tem_time_advance( me = params%general%simControl%now,   &
      &                    sim_dt = params%physics%dtLvl(iLevel ))

    ! --------------------------------------------------------------------------
    !set boundary for each field in current scheme
    call set_boundary( field       = scheme%field,                  &
      &                pdf         = scheme%pdf(iLevel),            &
      &                state       = scheme%state(iLevel)%val,      &
      &                levelDesc   = scheme%levelDesc(iLevel),      &
      &                tree        = geometry%tree,                 &
      &                iLevel      = iLevel,                        &
      &                nBCs        = geometry%boundary%nBCtypes,    &
      &                params      = params,                        &
      &                layout      = scheme%layout,                 &
      &                physics     = params%physics,                &
      &                varSys      = scheme%varSys,                 &
      &                mixture     = scheme%mixture,                &
      &                derVarPos   = scheme%derVarPos,              &
      &                globBC      = scheme%globBC                  )
    ! --------------------------------------------------------------------------

    ! swap double buffer index for current level
    call mus_swap_now_next( scheme%pdf( iLevel ) )
    now  = scheme%pdf(iLevel)%nNow
    next = scheme%pdf(iLevel)%nNext

    ! --------------------------------------------------------------------------
    ! Compute auxField from pre-collision state for fluid and ghostFromCoarser
    ! and exchange them if turbulence is active
    call tem_startTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    call mus_calcAuxFieldAndExchange(                         &
      & auxField          = scheme%auxField(iLevel),          &
      & calcAuxField      = scheme%calcAuxField,              &
      & state             = scheme%state(iLevel)%val(:, now), &
      & pdfData           = scheme%pdf(iLevel),               &
      & nFields           = scheme%nFields,                   &
      & field             = scheme%field(:),                  &
      & globSrc           = scheme%globSrc,                   &
      & stencil           = scheme%layout%fStencil,           &
      & varSys            = scheme%varSys,                    &
      & derVarPos         = scheme%derVarPos,                 &
      & general           = params%general,                   &
      & phyConvFac        = params%physics%fac(iLevel),       &
      & iLevel            = iLevel,                           &
      & minLevel          = geometry%tree%global%minLevel,    &
      & schemeHeader      = scheme%header,                        &
      & quantities        = scheme%layout%quantities              )
    call tem_stopTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! Update parameters, relaxation time .etc
    call tem_startTimer( timerHandle =  mus_timerHandles%relax(iLevel) )
    call mus_update_relaxParams( scheme  = scheme,                        &
      &                          iLevel  = iLevel,                        &
      &                          tNow    = params%general%simControl%now, &
      &                          physics = params%physics,                &
      &                          lattice = params%lattice,                &
      &                          nBCs    = geometry%boundary%nBCtypes     )
    call tem_stopTimer( timerHandle =  mus_timerHandles%relax(iLevel) )
    ! --------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Compute current scheme of current level
    call tem_startTimer( timerHandle =  mus_timerHandles%compute(iLevel) )

!$omp parallel
    call scheme%compute(                                        &
      &           fieldProp = scheme%field(:)%fieldProp,        &
      &           inState   = scheme%state(iLevel)%val(:,Now),  &
      &           outState  = scheme%state(iLevel)%val(:,Next), &
      &           auxField  = scheme%auxField(ilevel)%val(:),   &
      &           neigh     = scheme%pdf(iLevel)%neigh(:),      &
      &           nElems    = scheme%pdf(iLevel)%nSize,         &
      &           nSolve    = scheme%pdf(iLevel)%nElems_solve,  &
      &           level     = iLevel,                           &
      &           layout    = scheme%layout,                    &
      &           params    = params,                           &
      &           derVarPos = scheme%derVarPos,                 &
      &           varSys    = scheme%varSys                     )
!$omp end parallel

    call tem_stopTimer( timerHandle =  mus_timerHandles%compute(iLevel) )
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    call mus_apply_sourceTerms( field      = scheme%field(:),                &
      &                         nFields    = scheme%nFields,                 &
      &                         globSrc    = scheme%globSrc,                 &
      &                         pdf        = scheme%pdf(iLevel),             &
      &                         varSys     = scheme%varSys,                  &
      &                         iLevel     = iLevel,                         &
      &                         time       = params%general%simControl%now,  &
      &                         state      = scheme%state(iLevel)%val,       &
      &                         auxField   = scheme%auxField(iLevel)%val(:), &
      &                         derVarPos  = scheme%derVarPos(:),            &
      &                         phyConvFac = params%physics%fac(iLevel)      )
    ! -------------------------------------------------------------------------

    ! Communicate the halo elements of each scheme on current level
    ! KM: Communicate post-collision before set_boundary because nonEq_expol
    ! BC depends on post-collision from neighbor at next time step
    call tem_startTimer( timerHandle =  mus_timerHandles%comm(iLevel) )
    call params%general%commPattern%exchange_real(             &
      &    send    = scheme%levelDesc(iLevel)%sendbuffer,      &
      &    recv    = scheme%levelDesc(iLevel)%recvbuffer,      &
      &    state   = scheme%state(iLevel)%val(:,Next),         &
      &    message_flag   = iLevel,                            &
      &    comm    = params%general%proc%comm                  )

    ! communicate turbulent viscosity, required for interpolation
    if (trim(scheme%header%kind) == 'fluid' .or. &
      & trim(scheme%header%kind) == 'fluid_incompressible') then
      if (scheme%field(1)%fieldProp%fluid%turbulence%active) then
        call params%general%commPattern%exchange_real(                &
          & recv         = scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%recvbuffer,       &
          & send         = scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%sendbuffer,       &
          & state        = scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%visc(:),          &
          & message_flag = iLevel+100,                               &
          & comm         = params%general%proc%comm                  )
      end if
    end if
    call tem_stopTimer( timerHandle =  mus_timerHandles%comm(iLevel) )

    ! ... check if at least one of the IBMs is active
    if ( geometry%globIBM%nIBMs > 0 ) then
      call mus_buildBuffIBM(                              &
        &       me          = geometry%globIBM%IBM,       &
        &       commPattern = params%general%commPattern, &
        &       globTree    = geometry%tree,              &
        &       params      = params,                     &
        &       layout      = scheme%layout,              &
        &       levelDesc   = scheme%levelDesc(iLevel),   &
        &       iLevel      = iLevel                      )
    end if

    ! update the immersed boundaries if available
    ! ... and over the schemes
    ! ... check if at least one of the IBMs is active
    if( geometry%globIBM%nIBMs > 0 )then
      call mus_inamuro_IBM(                                    &
        &      me          = geometry%globIBM%IBM,             &
        &      commPattern = params%general%commPattern,       &
        &      globTree    = geometry%tree,                    &
        &      general     = params%general,                   &
        &      pdf         = scheme%pdf(iLevel),               &
        &      layout      = scheme%layout,                    &
        &      levelDesc   = scheme%levelDesc(iLevel),         &
        &      globSys     = scheme%varSys,                    &
        &      stateVarMap = scheme%stateVarMap%varPos%val(:), &
        &      convFac     = params%physics%fac(iLevel),       &
        &      iField      = 1,                                &
        &      state       = scheme%state(iLevel)%val,         &
        &      iLevel      = iLevel                            )
    end if
    ! -------------------------------------------------------------------------

  end subroutine do_fast_singleLevel
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine do_benchmark( me, scheme, geometry, params, iLevel )
    ! -------------------------------------------------------------------- !
    !> self control type
    !! dummy variable in this routine, required by interface
    class( mus_control_type ) :: me
    !> containers for the different schemes
    type( mus_scheme_type ), intent(inout)  :: scheme
    !> geometry infomation
    type( mus_geom_type ), intent(inout)    :: geometry
    !> global parameters
    type( mus_param_type ), intent(inout)   :: params
    !> Level counter variable
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: now, next
    ! -------------------------------------------------------------------- !

    ! Update auxField dependent source fields before adding source term to state
    ! and auxField such that both auxField and apply_source uses same source.
    call mus_update_sourceVars( nFields    = scheme%nFields,              &
      &                         field      = scheme%field,                &
      &                         globSrc    = scheme%globSrc,              &
      &                         varSys     = scheme%varSys,               &
      &                         iLevel     = iLevel,                      &
      &                         auxField   = scheme%auxField(iLevel)%val, &
      &                         phyConvFac = params%physics%fac(iLevel),  &
      &                         derVarPos  = scheme%derVarPos             )

    ! Increasing with the smallest time step (maxLevel)
    call tem_time_advance( me = params%general%simControl%now,   &
      &                    sim_dt = params%physics%dtLvl(iLevel ))

    ! --------------------------------------------------------------------------
    !set boundary for each field in current scheme
    call set_boundary( field       = scheme%field,               &
      &                pdf         = scheme%pdf(iLevel),         &
      &                state       = scheme%state(iLevel)%val,   &
      &                levelDesc   = scheme%levelDesc(iLevel),   &
      &                tree        = geometry%tree,              &
      &                iLevel      = iLevel,                     &
      &                nBCs        = geometry%boundary%nBCtypes, &
      &                params      = params,                     &
      &                layout      = scheme%layout,              &
      &                physics     = params%physics,             &
      &                varSys      = scheme%varSys,              &
      &                mixture     = scheme%mixture,             &
      &                derVarPos   = scheme%derVarPos,           &
      &                globBC      = scheme%globBC               )
    ! --------------------------------------------------------------------------

    ! swap double buffer index for current level
    call mus_swap_now_next( scheme%pdf( iLevel ) )
    now  = scheme%pdf(iLevel)%nNow
    next = scheme%pdf(iLevel)%nNext

    ! --------------------------------------------------------------------------
    ! Compute auxField from pre-collision state for fluid and ghostFromCoarser
    ! and exchange them if turbulence is active
    call tem_startTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    call mus_calcAuxFieldAndExchange(                             &
      & auxField          = scheme%auxField(iLevel),              &
      & calcAuxField      = scheme%calcAuxField,                  &
      & state             = scheme%state(iLevel)%val(:, now),     &
      & pdfData           = scheme%pdf(iLevel),                   &
      & nFields           = scheme%nFields,                       &
      & field             = scheme%field(:),                      &
      & globSrc           = scheme%globSrc,                       &
      & stencil           = scheme%layout%fStencil,               &
      & varSys            = scheme%varSys,                        &
      & derVarPos         = scheme%derVarPos,                     &
      & general           = params%general,                       &
      & phyConvFac        = params%physics%fac(iLevel),           &
      & iLevel            = iLevel,                               &
      & minLevel          = geometry%tree%global%minLevel,        &
      & schemeHeader      = scheme%header,                        &
      & quantities        = scheme%layout%quantities              )
    call tem_stopTimer( timerHandle =  mus_timerHandles%aux(iLevel) )
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! Update parameters, relaxation time .etc
    call mus_update_relaxParams( scheme  = scheme,                        &
      &                          iLevel  = iLevel,                        &
      &                          tNow    = params%general%simControl%now, &
      &                          physics = params%physics,                &
      &                          lattice = params%lattice,                &
      &                          nBCs    = geometry%boundary%nBCtypes     )
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! Compute current scheme of current level
    call tem_startTimer( timerHandle = mus_timerHandles%compute(iLevel) )

!$omp parallel
    call scheme%compute(                                        &
      &           fieldProp = scheme%field(:)%fieldProp,        &
      &           inState   = scheme%state(iLevel)%val(:,Now),  &
      &           outState  = scheme%state(iLevel)%val(:,Next), &
      &           auxField  = scheme%auxField(ilevel)%val(:),   &
      &           neigh     = scheme%pdf(iLevel)%neigh(:),      &
      &           nElems    = scheme%pdf(iLevel)%nSize,         &
      &           nSolve    = scheme%pdf(iLevel)%nElems_solve,  &
      &           level     = iLevel,                           &
      &           layout    = scheme%layout,                    &
      &           params    = params,                           &
      &           derVarPos = scheme%derVarPos,                 &
      &           varSys    = scheme%varSys                     )
!$omp end parallel

    call tem_stopTimer( timerHandle = mus_timerHandles%compute(iLevel) )
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    call mus_apply_sourceTerms( field      = scheme%field(:),                &
      &                         nFields    = scheme%nFields,                 &
      &                         globSrc    = scheme%globSrc,                 &
      &                         pdf        = scheme%pdf(iLevel),             &
      &                         varSys     = scheme%varSys,                  &
      &                         iLevel     = iLevel,                         &
      &                         time       = params%general%simControl%now,  &
      &                         state      = scheme%state(iLevel)%val,       &
      &                         auxField   = scheme%auxField(iLevel)%val(:), &
      &                         derVarPos  = scheme%derVarPos(:),            &
      &                         phyConvFac = params%physics%fac(iLevel)      )
    ! --------------------------------------------------------------------------

    ! Communicate the halo elements of each scheme on current level
    call tem_startTimer( timerHandle = mus_timerHandles%comm(iLevel) )
    call params%general%commPattern%exchange_real(                  &
      &         send    = scheme%levelDesc(iLevel)%sendbuffer,      &
      &         recv    = scheme%levelDesc(iLevel)%recvbuffer,      &
      &         state   = scheme%state(iLevel)%val(:,Next),         &
      &         message_flag   = iLevel,                            &
      &         comm    = params%general%proc%comm                  )

    ! communicate turbulent viscosity, required for interpolation
    if (trim(scheme%header%kind) == 'fluid' .or. &
      & trim(scheme%header%kind) == 'fluid_incompressible') then
      if (scheme%field(1)%fieldProp%fluid%turbulence%active) then
        call params%general%commPattern%exchange_real(                &
          & recv         = scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%recvbuffer,       &
          & send         = scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%sendbuffer,       &
          & state        = scheme%field(1)%fieldProp%fluid%turbulence &
          &                      %dataOnLvl(iLevel)%visc(:),          &
          & message_flag = iLevel+100,                               &
          & comm         = params%general%proc%comm                  )
      end if
    end if
    call tem_stopTimer( timerHandle = mus_timerHandles%comm(iLevel) )
   ! -------------------------------------------------------------------------

  end subroutine do_benchmark
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
