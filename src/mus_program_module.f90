! Copyright (c) 2013-2018, 2020-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2020-2021 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> author: Kannan Masilamani
!! This module contains data types and routines used by musubi main program
!! i.e unifying routines for apesmate
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
module mus_program_module

  ! include treelm modules
  use mpi
  use env_module,                   only: rk, long_k, newunit, rk_mpi, &
    &                                     pathLen, OutLen, labelLen
  use treelmesh_module,             only: treelmesh_type, free_treelmesh
  use tem_global_module,            only: tem_global_mesh_free
  use tem_comm_env_module,          only: tem_comm_env_type
  use tem_comm_module,              only: tem_comm_destroy
  use tem_timer_module,             only: tem_startTimer, tem_stopTimer
  use tem_status_module,            only: tem_status_run_end, tem_status_dump, &
    &                                     tem_status_run_terminate
  use tem_simControl_module,        only: tem_simControl_clearStat
  use tem_tracking_module,          only: tem_tracking_finalize
  use tem_restart_module,           only: tem_restart_finalize
  use tem_logging_module,           only: logUnit, tem_logging_init
  use tem_debug_module,             only: tem_reportStatus, tem_debug_type
  use tem_timeControl_module,       only: tem_timeControl_check, &
    &                                     tem_timeControl_dump
  use tem_tools_module,             only: tem_horizontalSpacer
  use tem_operation_var_module,     only: tem_opVar_reduction_transient_init
  use tem_time_module,              only: tem_time_dump
  use tem_adaptation_config_module, only: tem_adapt_type
  use tem_aux_module,               only: tem_abort
  use tem_grow_array_module,        only: destroy
  use tem_tracking_module,          only: tem_tracker

  ! include musubi modules
  ! RESPECT THE ORDER !!!!!!!!
  use mus_pdf_module,               only: pdf_data_type
  use mus_scheme_type_module,       only: mus_scheme_type
  use mus_param_module,             only: mus_param_type
  use mus_geom_module,              only: mus_geom_type
  use mus_tools_module,             only: check_density, mus_perf_measure, &
    &                                     mus_bc_timing, check_potential
  use mus_flow_module,              only: mus_init_flow
  use mus_dynLoadBal_module,        only: mus_perform_dynLoadBal
  ! ESPECIALLY OF THE FOLLOWING MODULES
  use mus_construction_module,      only: mus_construct
  use mus_control_module,           only: mus_control_type, mus_init_control
  use mus_aux_module,               only: mus_init_aux, check_flow_status
  use mus_restart_module,           only: mus_writeRestart
  use mus_IBM_module,               only: mus_finishIBM
  use mus_varSys_module,            only: mus_varSys_solverData_type
  use mus_mesh_adaptation_module,   only: mus_adapt_refine
  use mus_IBM_module,               only: mus_IBM_globType, mus_init_IBM
  use mus_scheme_module,            only: mus_init_scheme
  use mus_bc_general_module,        only: mus_init_boundary
  use mus_timer_module,             only: mus_timerHandles,     &
    &                                      mus_reset_mainTimer, &
    &                                      get_mainLoopTime,    &
    &                                      get_stageRatio, nStages
  use mus_weights_module,           only: mus_getWeights, mus_dumpWeights

  ! include aotus nmodules
  use aotus_module, only: close_config

  ! include particle modules
  use mus_particle_module,        only: mus_particle_group_type, &
    &                                   mus_particles_initialize
  use mus_particle_config_module, only: mus_load_particles, &
    &                                   mus_finalize_particleGroup
  use mus_particle_DPS_module,    only: mus_particles_initFluidVolumeFraction

  implicit none

  private

  public :: mus_initialize
  public :: mus_solve
  public :: mus_finalize


contains


! **************************************************************************** !
  !> This routine load musubi configuration file and initialize construction
  !! flow, auxilary and main control routines
  subroutine mus_initialize( scheme, geometry, params, particleGroup, control, &
    &                        solverData                                        )
    ! -------------------------------------------------------------------------
    !> scheme type
    type( mus_scheme_type ),         intent(inout) :: scheme
    !> Treelmesh data
    type( mus_geom_type ),           intent(inout) :: geometry
    !> Global parameters
    type( mus_param_type ),          intent(inout) :: params
    !> Particle parameters
    type( mus_particle_group_type ), intent(inout) :: particleGroup
    !> control routine
    type( mus_control_type ),        intent(inout) :: control
    !> contains pointer to scheme, physics types.
    !! passed to init_Scheme to build varSys
    type( mus_varSys_solverData_type ), target :: solverData
    ! ------------------------------------------------------------------------!
    integer :: iLevel, maxLevel
    ! ------------------------------------------------------------------------!
    maxLevel = geometry%tree%global%maxLevel


    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*) 'Initializing musubi ...'

    ! Initialize schemes: stencil, interpolation nSources and variable system
    call mus_init_scheme( me          = scheme,        &
      &                   tree        = geometry%tree, &
      &                   solverData  = solverData     )

    ! ... initialize the immersed boundary stencils
    call mus_init_IBM( me       = geometry%globIBM, &
      &                globTree = geometry%tree     )

    ! construct levelDescriptor, connectivity array and boundary elements
    call mus_construct( scheme    = scheme,   &
      &                 geometry  = geometry, &
      &                 params    = params    )

    ! Init auxiliary features such as interpolation, boundaries, restart
    ! and the tracker
    call mus_init_aux( scheme    = scheme,      &
      &                geometry  = geometry,    &
      &                params    = params       )

    !> Initialize flow field depends on read restart or initial condition
    call mus_init_flow( scheme       = scheme,                &
      &                 tree         = geometry%tree,         &
      &                 levelPointer = geometry%levelPointer, &
      &                 physics      = params%physics,        &
      &                 scaling      = params%scaling,        &
      &                 general      = params%general         )

    ! load particles
    call mus_load_particles(                              &
      &    particleGroup = particleGroup,                 &
      &    particle_kind = trim(params%particle_kind),    &
      &    conf          = params%general%solver%conf(1), &
      &    chunkSize     = 100,                           &
      &    scheme        = scheme,                        &
      &    geometry      = geometry,                      &
      &    myRank        = params%general%proc%rank       )

    do iLevel = geometry%tree%global%minLevel, geometry%tree%global%maxLevel
      ! Initialize also nNow since fNeq calculation requires it
      scheme%state(iLevel)%val(:, scheme%pdf(iLevel)%nNow)       &
        & = scheme%state(iLevel)%val(:, scheme%pdf(iLevel)%nNext)
    end do

    if ( trim(scheme%header%kind) == 'fluid_GNS' &
      & .OR. trim(scheme%header%kind) == 'fluid_incompressible_GNS' ) then
      call mus_particles_initFluidVolumeFraction(                            &
        &                       scheme   = scheme,                           &
        &                       geometry = geometry,                         &
        &                       nElems   = scheme%pdf(maxLevel)%nElems_local )
    end if

    ! initialize each field boundary by looping over the field inside
    ! init_boundary_field
    call mus_init_boundary( field        = scheme%field,     &
      &                     pdf          = scheme%pdf,       &
      &                     state        = scheme%state,     &
      &                     auxField     = scheme%auxField,  &
      &                     tree         = geometry%tree,    &
      &                     leveldesc    = scheme%levelDesc, &
      &                     layout       = scheme%layout,    &
      &                     schemeHeader = scheme%header,    &
      &                     varSys       = scheme%varSys,    &
      &                     derVarPos    = scheme%derVarPos, &
      &                     globBC       = scheme%globBC,    &
      &                     bc_prop      = geometry%boundary )

    call tem_horizontalSpacer(fUnit = logUnit(1))

    ! Initialize time reduction operation variable and set last index
    ! with initial value.
    call tem_opVar_reduction_transient_init(                &
      &      varSys         = scheme%varSys,                &
      &      tree           = geometry%tree,                &
      &      redTransVarMap = scheme%redTransVarMap,        &
      &      time           = params%general%simControl%now )

    write(logUnit(1),*) "Starting Musubi MAIN loop with time control:"
    call tem_timeControl_dump(params%general%simControl%timeControl,logUnit(1))

    if ( trim(params%particle_kind) /= 'none' ) then
      call mus_particles_initialize(        &
        &    particleGroup = particleGroup, &
        &    scheme        = scheme,        &
        &    geometry      = geometry,      &
        &    params        = params         )
    end if

    ! choose main control routine function pointer
    call mus_init_control( params%controlRoutine, control, &
      &                    geometry%tree%global%minLevel,  &
      &                    geometry%tree%global%maxLevel,  &
      &                    params%particle_kind            )

    call tem_horizontalSpacer(fUnit = logUnit(1))
  end subroutine mus_initialize
! **************************************************************************** !


! **************************************************************************** !
  !> This routine does the main musubi computation loop
  subroutine mus_solve( scheme, geometry, params, control, solverData, adapt )
    ! ------------------------------------------------------------------------!
    !> scheme type
    type( mus_scheme_type ), intent(inout) :: scheme
    !> Treelmesh data
    type( mus_geom_type ), intent(inout) :: geometry
    !> Global parameters
    type( mus_param_type ), intent(inout) :: params
    !> control routine
    type( mus_control_type ), intent(inout) :: control
    !> contains pointer to scheme, physics types
    type( mus_varSys_solverData_type ), target :: solverData
    !> mesh adaptation
    type(tem_adapt_type), intent(inout) :: adapt
    ! ------------------------------------------------------------------------!
    integer :: minLevel, maxLevel
    ! ------------------------------------------------------------------------!

    call tem_startTimer( timerHandle =  mus_timerHandles%mainLoop )
    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel

    ! --------------------------------------------------------------------------
    ! Start main loop
    mainloop: do

      ! check if some action has to be taken based on timeControl:
      ! tracking, global reduction operations, restart.
      ! check flow status including perform checks
      ! only at minLevel to complete one multilevel cycle.
      call check_flow_status( scheme   = scheme,                    &
        &                     general  = params%general,            &
        &                     physics  = params%physics,            &
        &                  mus_aborts  = params%mus_aborts,         &
        &           restart_triggered  = params%restart_triggered,  &
        &                    geometry  = geometry                   )

      if ( tem_status_run_end(params%general%simControl%status) .or. &
        & tem_status_run_terminate(params%general%simControl%status) ) then
        exit mainLoop
      end if

      ! clear status bit
      call tem_simControl_clearStat( me = params%general%simControl )

      ! Iterate through the elements in a level-wise fashion
      ! and advance the time steps
      call control%do_computation( iLevel = minLevel  )

      call do_balance()

    end do mainloop
    write(logUnit(6),"(A)") 'Finished main loop.'
    ! Finish main loop
    ! --------------------------------------------------------------------------

    call tem_stopTimer( timerHandle = mus_timerHandles%mainLoop )

  contains

    subroutine do_balance()
      real(kind=rk) :: total_density
      logical :: balance_triggered

      ! Perform the dynamic load balancing if interval has been reached
      if ( params%general%balance%dynamic ) then
        call tem_timeControl_check( me   = params%general%balance%timeControl, &
          &                         now  = params%general%simControl%now,      &
          &                         comm = params%general%proc%comm,           &
          &                         triggered = balance_triggered )

        if( balance_triggered ) then

          ! Stop main loop timer
          call tem_stopTimer( timerHandle = mus_timerHandles%mainLoop )

          ! -------------------------------------------------------------------
          select case (trim(scheme%header%kind))
          case ('poisson', 'poisson_boltzmann_linear', &
            &   'poisson_boltzmann_nonlinear'          )
            ! check final total potential
            call check_potential( scheme, minlevel, maxLevel, params%general, &
              &                   total_density )
          case default
            call check_density( scheme, minLevel, maxLevel, params%general, &
              &                 total_density                               )
          end select

          ! Measure imbalance and dump timing to disk
          call mus_perf_measure( &
            &       totalDens   = total_density,        &
            &       domSize     = geometry%tree%global%nElems, &
            &       minLevel    = minLevel,             &
            &       maxLevel    = maxLevel,             &
            &       nElems      = scheme%pdf(minLevel:maxLevel)%nElems_fluid, &
            &       scaleFactor = params%scaleFactor,   &
            &       general     = params%general )


          ! Only dump level timing for two levels mesh
          if ( params%dump_level_timing .and. ((maxLevel-minLevel)==1) ) then
            call dump_level_timing(                          &
              & minLevel   = minLevel,                       &
              & maxLevel   = maxLevel,                       &
              & pdf        = scheme%pdf(:),                  &
              & domSize    = geometry%tree%global%nElems,    &
              & timingFile = params%general%timingFile,      &
              & revision   = params%general%solver%revision, &
              & simName    = params%general%solver%simName,  &
              & proc       = params%general%proc             )
          end if
          ! -------------------------------------------------------------------

          call mus_reset_mainTimer()

          call tem_startTimer( timerHandle =  mus_timerHandles%balance )
          call mus_perform_dynLoadBal( scheme     = scheme,             &
            &                          params     = params,             &
            &                          geometry   = geometry,           &
            &                          solverData = solverData          )
          call tem_stopTimer( timerHandle =  mus_timerHandles%balance )

          if ( geometry%globIBM%nIBMs > 0 ) then
            call mus_finishIBM( me      = geometry%globIBM, &
              &                 params  = params,           &
              &                 useTime = .true.            )
          end if

          ! Restart mainloop timer
          call tem_startTimer( timerHandle =  mus_timerHandles%mainLoop )
          ! update min iter
          params%general%simControl%timeControl%min%iter = &
            & params%general%simControl%now%iter

        end if ! balance triggered
      end if !dynamic

    end subroutine do_balance

  end subroutine mus_solve
! **************************************************************************** !


! **************************************************************************** !
  !> Do final check on check on total density,
  !! Close auxiliary stuff such as restart and the tracker,
  !! finalize treelm, dump timing and finialize mpi with fin_env
  !!
  subroutine mus_finalize(scheme, params, particleGroup, tree, levelPointer, &
    &                     nBCs, globIBM)
    ! ------------------------------------------------------------------------!
    !> scheme type
    type(mus_scheme_type), intent(inout) :: scheme
    !> Global parameters
    type(mus_param_type), intent(inout) :: params
    !> Particle data
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> geometry infomation
    type(treelmesh_type),intent(inout) :: tree
    !> global information
    integer, intent(in) :: levelPointer(:)
    !> Number of BC
    integer, intent(in) :: nBCs
    !> global IBM datatype incl. array of IBM datatypes
    type(mus_IBM_globType), intent(inout) :: globIBM
    ! ------------------------------------------------------------------------!
    integer       :: minLevel, maxLevel, ii, iLevel
    real(kind=rk) :: total_density
    character(len=labelLen) :: bc_labels(nBCs)
    ! weights for fluid elements
    real(kind=rk), allocatable :: weights(:)
    ! ------------------------------------------------------------------------!

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel

    ! dump status information on how simulation was terminated
    if ( tem_status_run_terminate(params%general%simControl%status) ) then
      write(logUnit(1),'(/A)') 'Abnormal termination with the following status:'
      call tem_status_dump( me      = params%general%simControl%status, &
        &                   outUnit = logUnit(1)                        )
    else
      write(logUnit(1),'(/A)') 'SUCCESSFUL run!'
      write(logUnit(1),*) 'Status of the simulation now:'
      call tem_status_dump( me      = params%general%simControl%status, &
        &                   outUnit = logUnit(1)                        )
    end if


    ! Destroy communication buffers
    do iLevel = minLevel,maxLevel
      !send Buffer
      call tem_comm_destroy( scheme%levelDesc(iLevel)%sendbuffer, &
        &                    params%general%commPattern           )
      call tem_comm_destroy( scheme%levelDesc(iLevel)%sendBufferFromCoarser, &
        &                    params%general%commPattern                      )
      call tem_comm_destroy( scheme%levelDesc(iLevel)%sendBufferFromFiner, &
        &                    params%general%commPattern                    )
      !recv Buffer
      call tem_comm_destroy( scheme%levelDesc(iLevel)%recvbuffer, &
        &                    params%general%commPattern           )
      call tem_comm_destroy( scheme%levelDesc(iLevel)%recvbufferFromCoarser, &
        &                    params%general%commPattern                      )
      call tem_comm_destroy( scheme%levelDesc(iLevel)%recvbufferFromFiner, &
        &                    params%general%commPattern                    )
    end do

    select case (trim(scheme%header%kind))
    case ('poisson', 'poisson_boltzmann_linear', 'poisson_boltzmann_nonlinear')
      ! check final total potential
      call check_potential( scheme, minlevel, maxLevel, params%general, &
        &                   total_density )
    case default
      ! check final total density
      call check_density( scheme, minLevel, maxLevel, params%general, &
        &                 total_density                               )
    end select

    ! write restart file when it is not triggered by restart timecontrol
    ! before exiting mainloop.
    ! restart will be dumped even if timeControl table is not defined
    ! for restart and only writeRestart is active
    if ( params%general%restart%controller%writeRestart &
      &  .and. .not. params%restart_triggered           ) then
      write(logUnit(1),*) 'Writing restart at final time:'
      call tem_time_dump( params%general%simControl%now, logUnit(1) )
      call mus_writeRestart( levelPointer = levelPointer,                  &
        &                    restart      = params%general%restart,        &
        &                    scheme       = scheme,                        &
        &                    tree         = tree,                          &
        &                    timing       = params%general%simControl%now, &
        &                    timerHandle  = mus_timerHandles%wRestart      )
      write(logUnit(1),*) 'Done writing'
    end if

    ! close musubi main configuration file
    call close_config( L = params%general%solver%conf(1) )

    if( params%general%restart%controller%writeRestart .or.  &
      & params%general%restart%controller%readRestart ) then
      call tem_restart_finalize( params%general%restart )
    end if

    call tem_tracking_finalize( scheme%track )

    ! -------------------------------------------------------------------
    ! Dump performance measurement only for successful runs
    if (.not. tem_status_run_terminate(params%general%simControl%status) ) then
      ! Measure imbalance and dump timing to disk
      call mus_perf_measure( totalDens   = total_density,                 &
        &                    domSize     = tree%global%nElems,            &
        &                    minLevel    = minLevel,                      &
        &                    maxLevel    = maxLevel,                      &
        &                    nElems      = scheme%pdf(minLevel:maxLevel)  &
        &                                        %nElems_fluid,           &
        &                    scaleFactor = params%scaleFactor,            &
        &                    general     = params%general                 )

      ! -------------------------------------------------------------------
      ! Only dump level timing for two levels mesh
      if ( params%dump_level_timing .and. ((maxLevel-minLevel)==1) ) then
        call dump_level_timing( minLevel   = minLevel,                       &
          &                     maxLevel   = maxLevel,                       &
          &                     pdf        = scheme%pdf(:),                  &
          &                     timingFile = params%general%timingFile,      &
          &                     domSize    = tree%global%nElems,             &
          &                     revision   = params%general%solver%revision, &
          &                     simName    = params%general%solver%simName,  &
          &                     proc       = params%general%proc             )
      end if

      ! Dump weights if write_weights is not empty
      if ( trim(tree%write_weights) /= '') then
        allocate( weights( tree%nElems ) )
        ! Calculate weights according to compute, intp and bc routines
        call mus_getWeights( weights   = weights,          &
          &                  tree      = tree,             &
          &                  minLevel  = minLevel,         &
          &                  maxLevel  = maxLevel,         &
          &                  levelDesc = scheme%levelDesc, &
          &                  nBCs      = nBCs,             &
          &                  globBC    = scheme%globBC     )
        ! Dump weights
        call mus_dumpWeights( tree     = tree,              &
          &                   weights  = weights,           &
          &                   basename = tree%write_weights )
        deallocate( weights )
      end if

      do ii = 1, nBCs
        bc_labels(ii) = scheme%globBC(ii)%label
      end do
      call mus_BC_timing( nBCs, bc_labels, params%general%proc%comm )

      ! dump IBM timings
      if( globIBM%nIBMs > 0 )then
        call mus_finishIBM( me      = globIBM, &
          &                 params  = params,  &
          &                 useTime = .false.  )
      end if
    end if ! if successful run

    !free mesh
    call tem_global_mesh_free(tree%global)

    ! De-allocate particleGroup object
    if ( params%particle_kind /= 'none' ) then
      call mus_finalize_particleGroup( particleGroup )
    end if


  end subroutine mus_finalize
! **************************************************************************** !


! **************************************************************************** !
  subroutine dump_level_timing( minLevel, maxLevel, pdf, timingFile, proc, &
      &                         DomSize, revision, simName )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: minLevel, maxLevel
    type( pdf_data_type ), intent(in) :: pdf( minLevel:maxLevel )
    character(len=*), intent(in) :: timingFile
    type( tem_comm_env_type ), intent(in) :: proc
    integer(kind=long_k), intent(in) :: DomSize
    character(len=*), intent(in) :: revision
    character(len=*), intent(in) :: simName
    ! --------------------------------------------------------------------------
    character(len=OutLen) :: header
    character(len=OutLen) :: output
    character(len=PathLen) :: filename
    integer :: fileunit, iLevel, ii
    !> MPI related variables
    integer :: iErr, msgsize, fileno, ErrErr
    real(kind=rk) :: tMainLoop, stageRatio(4), stageAve(4), stageMax(4)
    integer :: ioStatus( mpi_status_size )
    character(len=100) :: IOError
    integer :: resultlen = 100
    ! --------------------------------------------------------------------------

    tMainLoop = get_mainLoopTime()
    do ii = 1, 4
      stageRatio(ii) = get_stageRatio( ii )
    end do
    call mpi_reduce( stageRatio, stageMax, 4, rk_mpi, mpi_max, 0, proc%comm, &
      &              iErr                                                    )
    call mpi_reduce( stageRatio, stageAve, 4, rk_mpi, mpi_sum, 0, proc%comm, &
      &              iErr                                                    )
    stageAve = stageAve / dble( proc%comm_size )

    write(filename, "(A,I6.6,A)") 'level_'//trim(timingFile)
    if ( proc%isRoot ) then

      fileunit = newunit()
      open(unit=fileunit,file=trim(filename),position='append')

      ! Write Header ---------------------------------------------------
      !                                      12345678901234567890
      header = ''
      write(header,"(A,I0,A)") "# number of processes: ", proc%comm_size, &
        &                      new_line('A')
      write(header,"(2A,I0,A)") trim(header), "# DomSize: ", DomSize, &
        &                       new_line('A')
      write(header,"(2A,A,A )") trim(header), "# Revision: ", trim(revision), &
        &                       new_line('A')
      write(header,"(2A,2A)")   trim(header), "# SimName: ", trim(simName), &
        &                       new_line('A')
      write(header,"(2A,F8.2,A)") trim(header), "# MainLoopTime: ", tMainLoop, &
        &                         new_line('A')

      ! write min and max values
      write(header,"(2A,66X,4F12.1,A)") trim(header), "# Stage Ave:", &
        &                               stageAve(1:4), new_line("A")
      write(header,"(2A,66X,4F12.1,A)") trim(header), "# Stage Max:", &
        &                               stageMax(1:4), new_line('A')
      write(header,"(2A,62X,4F12.1,A)") trim(header), "# Ave / Max (%):",   &
        &                               stageAve(1:4)/stageMax(1:4)*100_rk, &
        &                               new_line('A')

      ! Write column information
      write(header,"(A,A6)")  trim(header), "# Rank"
      ! write elements information
      do iLevel = minLevel, maxLevel
        write(header,"(A,A10,I2.2)") trim(header), "nSolve_L", iLevel
        write(header,"(A,A10,I2.2)") trim(header), "nGFC_L", iLevel
        write(header,"(A,A10,I2.2)") trim(header), "nGFF_L", iLevel
      end do
      ! write all stages
      do ii = 1, 4
        write(header,"(A,A10,I2.2)")  trim(header), "stage_", ii
      end do

      write(fileunit,'(a)') trim(header)
      close(fileunit)
    end if
    ! Write Header ---------------------------------------------------

    ! Make sure root finishes writing header
    call MPI_BARRIER(proc%comm, IERR)

    ! Write Data ---------------------------------------------------
    output = ''
    ! Write data into string
    write(output,"(A,I6)") trim(output), proc%rank
    do iLevel = minLevel, maxLevel
      write(output,"(A,I12)") trim(output), pdf(iLevel)%nElems_solve
      write(output,"(A,I12)") trim(output), pdf(iLevel)%nElems_ghostFromCoarser
      write(output,"(A,I12)") trim(output), pdf(iLevel)%nElems_ghostFromFiner
    end do
    ! write all stages
    do ii = 1, 4
      write(output,"(A,F12.1)")  trim(output), get_stageRatio( ii )
    end do
    write(output,"(A,A)")  trim(output), new_line('A')

    ! Write string into file
    msgsize = len_trim( output )

    call MPI_File_open( proc%comm, trim(filename),                           &
      &                 ior(MPI_MODE_APPEND,MPI_MODE_WRONLY), MPI_INFO_NULL, &
      &                 fileno, ierr                                         )
    if (iErr /= MPI_SUCCESS) then
      call MPI_ERROR_STRING( iErr, IOError, resultlen, ErrErr )
      write(logUnit(0),*) 'Read File Open MPI WARNING: '// trim( IOError )
      call tem_abort()
    else
      call MPI_File_write_ordered(fileno, trim(output), msgsize, &
        &                         MPI_CHARACTER, iostatus, ierr  )
      call MPI_File_close(fileno, ierr)
    end if
    ! Write Data ---------------------------------------------------


  end subroutine dump_level_timing
! **************************************************************************** !

end module mus_program_module
! **************************************************************************** !
