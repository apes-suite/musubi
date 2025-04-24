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
!> mus_particle_comm_module contains the routines used for communication of
!! particle data between processes

module mus_particle_comm_module

  use mpi

  use env_module,                       only: rk, rk_mpi, long_k

  use tem_aux_module,                   only: tem_abort
  use tem_dyn_array_module,             only: dyn_intArray_type, append, init
  use tem_topology_module,              only: tem_IdOfCoord
  use tem_timer_module,                 only: tem_startTimer, tem_stopTimer
  use tem_logging_module,               only: logUnit
  use mus_geom_module,                  only: mus_geom_type
  use mus_scheme_type_module,           only: mus_scheme_type
  use mus_param_module,                 only: mus_param_type

  use mus_particle_type_module,         only: mus_particle_MEM_type,     &
    &                                         mus_particle_DPS_type,     &
    &                                         mus_particle_group_type,   &
    &                                         printParticleGroup,        &
    &                                         printParticleGroup2_DPS,   &
    &                                         printPIDlist,              &
    &                                         sortPosofval_particle_MEM, &
    &                                         sortPosofval_particle_DPS, &
    &                                         append_da_particle_MEM,    &
    &                                         append_da_particle_DPS
  use mus_particle_aux_module,          only: mus_particles_global_errorcheck
  use mus_particle_comm_type_module,    only:       &
    &   init_mus_particles_comm_type,               &
    &   mus_particles_comm_init_vectorbuffer,       &
    &   mus_particles_comm_init_wallbuffer,         &
    &   mus_particles_comm_init_statebuffer,        &
    &   mus_particles_comm_init_posbuffer,          &
    &   mus_particles_comm_init_IDbuffer,           &
    &   mus_particles_comm_init_particlebuffer,     &
    &   mus_particles_initForceContributionMPItype, &
    &   mus_particles_initParticleStateMPItype,     &
    &   mus_particles_initPositionUpdateMPItype,    &
    &   mus_particles_initWallPosMPItype,           &
    &   mus_particles_initParticleInfoMPItype,      &
    &   print_particles_comm,                       &
    &   print_particles_pIDvectorbuffer,            &
    &   find_particle_comm_procs4,                  &
    &   mus_particles_communication_type,           &
    &   init_mus_particles_comm_type,               &
    &   mus_positionUpdate_type,                    &
    &   mus_wallPos_type,                           &
    &   mus_pIDvector_type,                         &
    &   mus_particleState_type,                     &
    &   mus_particleInfo_type
  use mus_particle_logging_type_module, only: mus_particle_logging_type, &
    &   pgDebugLog
  use mus_particle_timer_module,        only: mus_particle_timerHandles

  implicit none

  public :: exchangeParticlesToRemove
  public :: exchangeForces
  public :: exchangePositions
  public :: exchangeVelocities
  public :: exchangeParticleStates
  public :: exchangeNewParticles_MEM


contains


  subroutine mus_particles_initialize_communication( particleGroup, scheme, &
    &                                                geometry, params       )
    !> Array of particles
    type(mus_particle_group_type), target :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    ! --------------------------------------------------- !
    real(kind=rk) :: dpad, dx
    integer :: nProcs, maxParticles
    type(dyn_intarray_type) :: particleCommProcs
    integer :: lev

    ! --------------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    dx = params%physics%dxLvl(lev)

      ! Init derived MPI datatypes
    call mus_particles_initForceContributionMPItype()
    call mus_particles_initPositionUpdateMPItype()
    call mus_particles_initParticleStateMPItype()
    call mus_particles_initParticleInfoMPItype()
    call mus_particles_initWallPosMPItype()

    call init(me = particleCommProcs, length = 1)

    ! Find processes to communicate particle data with
    ! NOTE: this is not necessary for DPS simulations - there the procs to
    ! communicate with will be the same as for the fluid since particles are
    ! smaller than a lattice site.
    select case( trim(params%particle_kind) )
      case( 'MEM' )
        if( params%particle_kind /= 'none') then
          dpad = particleGroup%halo_distance
          call find_particle_comm_procs4( prunedProcs    = particleCommProcs,        &
            &                             scheme         = scheme,                   &
            &                             geometry       = geometry,                 &
            &                             myRank         = params%general%proc%rank, &
            &                             dpad           = dpad                      )
          nProcs = particleCommProcs%nvals
          open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
            write(pgDebugLog%lu,*) 'init_particles_communication: MEM'
            write(pgDebugLog%lu,*) particleCommProcs%val(1:nProcs)
          close( pgDebugLog%lu )
        else
          dpad = 0
          nProcs = 0
        end if
      case( 'DPS', 'DPS_unittest', 'DPS_twoway' )
        if( params%particle_kind /= 'none') then
          ! For DPS the procs to communicate particles with are the same as the
          ! procs in send + recv buffer for fluid elements
          dpad = particleGroup%halo_distance
          call find_particle_comm_procs4( prunedProcs = particleCommProcs,        &
                                        & scheme      = scheme,                   &
                                        & geometry    = geometry,                 &
                                        & myRank      = params%general%proc%rank, &
                                        & dpad        = dpad                      )
          nProcs = particleCommProcs%nvals
          open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
            write(pgDebugLog%lu,*) 'init_particles_communication: DPS'
            write(pgDebugLog%lu,*) particleCommProcs%val(1:nProcs)
          close( pgDebugLog%lu )
        else
          dpad = 0
          nProcs = 0
        end if
      case( 'DPS_oneway' )
        if( params%particle_kind /= 'none') then
          dpad = dx
          call find_particle_comm_procs4( prunedProcs = particleCommProcs,        &
                                        & scheme      = scheme,                   &
                                        & geometry    = geometry,                 &
                                        & myRank      = params%general%proc%rank, &
                                        & dpad        = dpad                      )
          nProcs = particleCommProcs%nvals
          open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
            write(pgDebugLog%lu,*) 'init_particles_communication: DPS_oneway'
            write(pgDebugLog%lu,*) particleCommProcs%val(1:nProcs)
          close( pgDebugLog%lu )
        else
          dpad = 0
          nProcs = 0
        end if
      case default
        write(logUnit(1),*) 'ERROR mus_particles_initialize_communication: ', &
          &                 'unknown particle scheme kind'
    end select


      ! maximum number of particles that can be in force buffers for
    ! communication routines.
    maxParticles = particleGroup%particleBufferSize

    write(logUnit(1),*) "Allocating particle buffers with size: ", maxParticles

    ! Initialize send and recv communication types
    call init_mus_particles_comm_type( me     = particleGroup%send,             &
                                     & nProcs = nProcs,                         &
                                     & proc   = particleCommProcs%val(1:nProcs) )

    call init_mus_particles_comm_type( me     = particleGroup%recv,             &
                                     & nProcs = nProcs,                         &
                                     & proc   = particleCommProcs%val(1:nProcs) )


    call mus_particles_comm_init_buffers( particleGroup = particleGroup, &
                                        & nProcs = nProcs,               &
                                        & maxParticles = maxParticles    )


  end subroutine mus_particles_initialize_communication

  subroutine mus_particles_comm_init_buffers(particleGroup, nProcs, maxParticles)
    !> Array of particles
    type(mus_particle_group_type), target :: particleGroup
    !> Number of processes to communicate with
    integer, intent(in) :: nProcs
    !> Maximum number of particles that can be held in communication buffer
    integer, intent(in) :: maxParticles
    ! --------------------------------------------------- !
    integer :: iproc
    ! --------------------------------------------------- !
      ! Initialize the buffers
    ! NOTE: should this maybe just be part of init_mus_particles_comm_type or some larger routine?
    do iproc = 1, nProcs
      ! Force buffer
      call mus_particles_comm_init_vectorbuffer(       &
        &    me = particleGroup%send%buf_force(iproc), &
        &    maxParticles = maxParticles               )

      ! Position buffer
      call mus_particles_comm_init_posbuffer(        &
        &    me = particleGroup%send%buf_pos(iproc), &
        &    maxParticles = maxParticles             )

      ! Velocity buffer
      call mus_particles_comm_init_vectorbuffer(     &
        &    me = particleGroup%send%buf_vec(iproc), &
        &    maxParticles = maxParticles             )

      ! Wall position buffer
      call mus_particles_comm_init_wallbuffer(        &
        &    me = particleGroup%send%buf_wall(iproc), &
        &    maxParticles = maxParticles              )

      ! Buffer for sending particle state (position + velocity + coordOfOrigin)
      call mus_particles_comm_init_statebuffer(        &
        &    me = particleGroup%send%buf_state(iproc), &
        &    maxParticles = maxParticles               )

      ! Buffer for sending entire particles
      call mus_particles_comm_init_particlebuffer(        &
        &    me = particleGroup%send%buf_particle(iproc), &
        &    maxParticles = maxParticles                  )

      ! Buffer for sending particle IDs
      call mus_particles_comm_init_IDbuffer(          &
        &    me = particleGroup%send%buf_kill(iproc), &
        &    maxParticles = maxParticles              )
    end do

    do iproc = 1, nProcs
      ! Force buffer
      call mus_particles_comm_init_vectorbuffer(       &
        &    me = particleGroup%recv%buf_force(iproc), &
        &    maxParticles = maxParticles               )

      ! Position buffer
      call mus_particles_comm_init_posbuffer(        &
        &    me = particleGroup%recv%buf_pos(iproc), &
        &    maxParticles = maxParticles             )
      ! Velocity buffer
      call mus_particles_comm_init_vectorbuffer(      &
        &     me = particleGroup%recv%buf_vec(iproc), &
        &     maxParticles = maxParticles             )
      ! Wall position buffer
      call mus_particles_comm_init_wallbuffer(        &
        &    me = particleGroup%recv%buf_wall(iproc), &
        &    maxParticles = maxParticles              )

      ! Buffer for sending particle state (position + velocity + coordOfOrigin)
      call mus_particles_comm_init_statebuffer(        &
        &    me = particleGroup%recv%buf_state(iproc), &
        &    maxParticles = maxParticles               )


      ! Buffer for receiving entire particles
      call mus_particles_comm_init_particlebuffer(        &
        &    me = particleGroup%recv%buf_particle(iproc), &
        &    maxParticles = maxParticles                  )

      ! Buffer for sending particle IDs
      call mus_particles_comm_init_IDbuffer(          &
        &    me = particleGroup%recv%buf_kill(iproc), &
        &    maxParticles = maxParticles              )
    end do


  end subroutine mus_particles_comm_init_buffers

  !> If a particle is removed from the global domain (e.g. after colliding with an open boundary)
  !! then exchangeParticlesToRemove will send messages informing all other processes
  !! that this particle should be removed.
  subroutine exchangeParticlesToRemove(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    integer :: recvPID
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force and ID buffers for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_force(iproc)%nParticles = 0
      send%buf_kill(iproc)%nParticles = 0
    end do ! iproc

    ! ---- FILL KILL BUFFERS ---- !
    do iParticle = 1, this%particles_MEM%nvals
      if( this%particles_MEM%val(iParticle)%removeParticle_global ) then
        ! If particle should be removed from global domain, find the
        ! other procs on which particle existsOnProc to let them know
        do iproc = 1, send%nProcs
          if( this%particles_MEM%val(iParticle)%existsOnProc(iproc) ) then
            ! Check if buffer is large enough for the new contribution
            if( send%buf_kill(iproc)%nParticles+1                   &
              &               >= send%buf_kill(iproc)%maxParticles  ) then
              call tem_abort('Error exchangeParticlesToRemove: ID buffer too small')
            end if

            ! Increment number of particle contributions in buffer
            send%buf_kill(iproc)%nParticles = send%buf_kill(iproc)%nParticles + 1

            ! Add particle ID of particle-to-be-removed to end of buffer
            send%buf_kill(iproc)%val( send%buf_kill(iproc)%nParticles ) &
              & = this%particles_MEM%val(iParticle)%particleID

          end if ! existsOnProc
        end do ! iproc
      end if ! removeParticle_global
    end do ! iParticle

    ! Now the messages for all the processes should be complete.

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !

    ! --- START RECEIVE OPERATIONS --- !
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes

      ! start receiving particle IDs of particles to be removed/killed
      call mpi_irecv(                               &
       &      recv%buf_kill( iproc )%val,           & ! me
       &      recv%buf_kill( iproc )%maxParticles,  & ! max particles me can receive
       &      MPI_INTEGER,                          & ! data type
       &      recv%proc(iproc),                     & ! source process
       &      message_flag,                         & ! flag
       &      comm,                                 & ! communicator
       &      recv%rqhandle2(iproc),                & ! handle
       &      ierr )                                  ! error status

    end do

    ! ----- START SENDING OPERATIONS ------ !
    do iproc = 1, send%nprocs
      ! start sending particle IDs of particles to be removed/killed
      call mpi_isend(                                   &
       &      send%buf_kill( iproc )%val,               & ! buffer
       &      send%buf_kill( iproc )%nParticles,        & ! count
       &      MPI_INTEGER,                              & ! data type
       &      send%proc(iproc),                         & ! target
       &      message_flag,                             & ! integer tag
       &      comm,                                     & ! communicator
       &      send%rqhandle2( iproc ),                  & ! handle
       &      ierr )                                      ! error status

    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      ! First wait on the particle IDs
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle2,            & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! ---- READ THE PARTICLES-TO-BE-REMOVED ---- !
    do iproc = 1, recv%nprocs
      ! Get the number of particle IDs received from this proc
      call mpi_get_count( status   = recv_status(:,iproc),   & ! mpi status
                        & datatype = MPI_INTEGER,             & ! data type
                        & count    = recvcount,               & ! count
                        & ierror   = ierr                     )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit
      end do ! posOfProc

      ! Update the number of items recv buffer
      recv%buf_kill(posOfProc)%nParticles = recvcount

      ! Loop over received particle IDs
      do iBuff = 1, recv%buf_kill(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        recvPID = recv%buf_kill(posOfProc)%val(iBuff)
        pos = sortPosofval_particle_MEM( me  = this%particles_MEM, &
          &                              pID = recvPID             )
        if (pos == 0) then
          ! particle has already been deleted on this process so move on
          cycle
        else
          ! Convert to the index of corresponding particle in pIDlist
          pos = this%particles_MEM%pIDsort(pos)

          ! Indicate that this particle should be removed upon exiting this routine
          this%particles_MEM%val(pos)%removeParticle_global = .TRUE.

          open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
          write(pgDebugLog%lu,*) 'exchangeParticlesToRemove: removing particle with pID ', &
            & recvPID, ' from particle group on rank ', myRank, ' at request of proc ', proc
          close( pgDebugLog%lu )

        end if ! pos == 0

      end do ! iBuff
    end do ! iproc

    ! wait for send buffers to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,    & ! count
        &              send%rqhandle2, & ! request handles
        &              send_status,    & ! mpi status
        &              ierr            ) ! error status
    end if

  end subroutine exchangeParticlesToRemove

  !> If a particle is removed from the global domain (e.g. after colliding with an open boundary)
  !! then exchangeParticlesToRemove will send messages informing all other processes
  !! that this particle should be removed.
  subroutine exchangeParticlesToRemove_DPS(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    integer :: recvPID
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force and ID buffers for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_force(iproc)%nParticles = 0
      send%buf_kill(iproc)%nParticles = 0
    end do ! iproc

    ! ---- FILL KILL BUFFERS ---- !
    do iParticle = 1, this%particles_DPS%nvals
      if( this%particles_DPS%val(iParticle)%removeParticle_global ) then
        ! If particle should be removed from global domain, find the
        ! other procs on which particle existsOnProc to let them know
        do iproc = 1, send%nProcs
          if( this%particles_DPS%val(iParticle)%existsOnProc(iproc) ) then
            ! Check if buffer is large enough for the new contribution
            if( send%buf_kill(iproc)%nParticles+1                   &
              &               >= send%buf_kill(iproc)%maxParticles  ) then
              call tem_abort('Error exchangeParticlesToRemove: ID buffer too small')
            end if

            ! Increment number of particle contributions in buffer
            send%buf_kill(iproc)%nParticles = send%buf_kill(iproc)%nParticles + 1

            ! Add particle ID of particle-to-be-removed to end of buffer
            send%buf_kill(iproc)%val( send%buf_kill(iproc)%nParticles ) &
              & = this%particles_DPS%val(iParticle)%particleID

          end if ! existsOnProc
        end do ! iproc
      end if ! removeParticle_global
    end do ! iParticle

    ! Now the messages for all the processes should be complete.

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !

    ! --- START RECEIVE OPERATIONS --- !
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes

      ! start receiving particle IDs of particles to be removed/killed
      call mpi_irecv(                               &
       &      recv%buf_kill( iproc )%val,           & ! me
       &      recv%buf_kill( iproc )%maxParticles,  & ! max particles me can receive
       &      MPI_INTEGER,                          & ! data type
       &      recv%proc(iproc),                     & ! source process
       &      message_flag,                         & ! flag
       &      comm,                                 & ! communicator
       &      recv%rqhandle2(iproc),                & ! handle
       &      ierr )                                  ! error status

    end do

    ! ----- START SENDING OPERATIONS ------ !
    do iproc = 1, send%nprocs
      ! start sending particle IDs of particles to be removed/killed
      call mpi_isend(                                   &
       &      send%buf_kill( iproc )%val,               & ! buffer
       &      send%buf_kill( iproc )%nParticles,        & ! count
       &      MPI_INTEGER,                              & ! data type
       &      send%proc(iproc),                         & ! target
       &      message_flag,                             & ! integer tag
       &      comm,                                     & ! communicator
       &      send%rqhandle2( iproc ),                  & ! handle
       &      ierr )                                      ! error status

    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      ! First wait on the particle IDs
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle2,            & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! ---- READ THE PARTICLES-TO-BE-REMOVED ---- !
    do iproc = 1, recv%nprocs
      ! Get the number of particle IDs received from this proc
      call mpi_get_count( status   = recv_status(:,iproc), & ! mpi status
                        & datatype = MPI_INTEGER,          & ! data type
                        & count    = recvcount,            & ! count
                        & ierror   = ierr                  )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit
      end do ! posOfProc

      ! Update the number of items recv buffer
      recv%buf_kill(posOfProc)%nParticles = recvcount

      ! Loop over received particle IDs
      do iBuff = 1, recv%buf_kill(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        recvPID = recv%buf_kill(posOfProc)%val(iBuff)
        pos = sortPosofval_particle_DPS( me  = this%particles_DPS,                      &
                                       & pID = recvPID )
        if(pos == 0) then
          ! particle has already been deleted on this process so move on
          cycle
        else
          ! Convert to the index of corresponding particle in pIDlist
          pos = this%particles_DPS%pIDsort(pos)

          ! Indicate that this particle should be removed upon exiting this routine
          this%particles_DPS%val(pos)%removeParticle_global = .TRUE.

          open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
          write(pgDebugLog%lu,*) 'exchangeParticlesToRemove: removing particle with pID ', &
            & recvPID, ' from particle group on rank ', myRank, ' at request of proc ', proc
          close( pgDebugLog%lu )

        end if ! pos == 0

      end do ! iBuff
    end do ! iproc

    ! wait for send buffers to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,    & ! count
        &              send%rqhandle2, & ! request handles
        &              send_status,    & ! mpi status
        &              ierr )            ! error status
    end if

  end subroutine exchangeParticlesToRemove_DPS

  !> exchangeForces send force contributions from each process on which a particle
  !! exists to the particle owner. The total hydrodynamic force is then computed by
  !! the particle owner as the sum of all the force contributions (which physically are
  !! the surface forces from the part of the particle surface handled by each process).
  !! The owner then updates particle velocity in a subsequent call to a different routine.
  !! For momentum-exchange method only!
  subroutine exchangeForces(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    class(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    integer :: recvPID
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_force(iproc)%nParticles = 0
    end do ! iproc

    do iParticle = 1, this%particles_MEM%nvals
      if(this%particles_MEM%val(iParticle)%owner == myRank) cycle

      ! Add force contribution for this particle to the message for its owner
      ! First find the owner by looping over neighbor procs
      do iproc = 1, send%nProcs
        if(this%particles_MEM%val(iParticle)%owner == send%proc(iproc))then
          ! We found the owner of this particle, now add the
          ! force contribution to the message for this process

          ! Check if force buffer is large enough for the new contribution
          if( send%buf_force(iproc)%nParticles+1                   &
            &               >= send%buf_force(iproc)%maxParticles  ) then
            call tem_abort('Error sendForces: force buffer too small, aborting')
          end if

          ! Increment number of particle contributions in force buffer
          send%buf_force(iproc)%nParticles = send%buf_force(iproc)%nParticles + 1

          ! Add force contribution to end of the buffer
          ! send%buf_force(iproc)%val( send%buf_force(iproc)%nParticles )%V &
          !   & = this%particles%val(iParticle)%F
          send%buf_force(iproc)%val( send%buf_force(iproc)%nParticles )%V &
            & = this%particles_MEM%val(iParticle)%Fbuff( this%particles_MEM%val(iParticle)%Fnow,1:6)

          ! Add the ID of the particle to which this force contribution pertains
          send%buf_force(iproc)%val( send%buf_force(iproc)%nParticles )%pID &
            & = this%particles_MEM%val(iParticle)%particleID

          ! If I was the previous owner of this particle (but am not owner now)
          ! Send over Fbuff(Flast) which contains the TOTAL force on the particle in
          ! the last time step. We send the NEGATIVE particle ID along to indicate that
          ! this is not a force contribution but a total force
          if( this%particles_MEM%val(iParticle)%previousOwner == myRank ) then
            ! Check if force buffer is large enough for the new contribution
            if( send%buf_force(iproc)%nParticles+1                   &
              &               >= send%buf_force(iproc)%maxParticles  ) then
              call tem_abort('Error sendForces: force buffer too small, aborting')
            end if

            ! Increment number of particle contributions in force buffer
            send%buf_force(iproc)%nParticles = send%buf_force(iproc)%nParticles + 1

            ! Add force contribution to end of the buffer
            send%buf_force(iproc)%val( send%buf_force(iproc)%nParticles )%V &
              & = this%particles_MEM%val(iParticle)%Fbuff( this%particles_MEM%val(iParticle)%Flast,1:6)

            ! Add the ID of the particle to which this force contribution pertains
            send%buf_force(iproc)%val( send%buf_force(iproc)%nParticles )%pID &
              & = -this%particles_MEM%val(iParticle)%particleID


          end if ! I am the previous owner
        end if ! iproc is the current owner of this particle
      end do ! iproc
    end do ! iParticle


    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                               &
       &      recv%buf_force( iproc )%val,          & ! me
       &      recv%buf_force( iproc )%maxParticles, & ! max particles me can receive
       &      mus_pIDvector_type,                   & ! data type
       &      recv%proc(iproc),                     & ! source process
       &      message_flag,                         & ! flag
       &      comm,                                 & ! communicator
       &      recv%rqhandle(iproc),                 & ! handle
       &      ierr )                                  ! error status
    end do

    !  start the sending of force contributions
    do iproc = 1, send%nprocs
      call mpi_isend(                             &
       &      send%buf_force( iproc )%val,        & ! buffer
       &      send%buf_force( iproc )%nParticles, & ! count
       &      mus_pIDvector_type,                 & ! data type
       &      send%proc(iproc),                   & ! target
       &      message_flag,                       & ! integer tag
       &      comm,                               & ! communicator
       &      send%rqhandle( iproc ),             & ! handle
       &      ierr                                ) ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,   & ! count
        &              recv%rqhandle, & ! request handles
        &              recv_status,   & ! mpi status
        &              ierr           ) ! error status
    end if

    ! Now force contributions from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc),   & ! mpi status
                        & datatype = mus_pIDvector_type,   & ! data type
                        & count = recvcount,               & ! count
                        & ierror = ierr                    )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_force(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_force(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        recvPID = abs( recv%buf_force(posOfProc)%val(iBuff)%pID )
        pos = sortPosofval_particle_MEM( me  = this%particles_MEM, &
          &                              pID = recvPID             )
        if(pos == 0) then
          write(logUnit(1),*) '---------------------------'
          write(logUnit(1),*) &
            & 'WARNING exchangeForces: cannot find particle with ID', recv%buf_force(iProc)%val(iBuff)%pID
          write(logUnit(1),*) 'in particleGroup on rank = ', myRank
          write(logUnit(1),*) '--- RECV buff: ---'
          call print_particles_comm(recv)
          write(logUnit(1),*) '--- Rank ', myRank, ' particle group: ---'
          call printParticleGroup(this, logUnit(1))
          write(logUnit(1),*) '--- PIDlists ---'
          call printPIDlist(this)
          write(logUnit(1),*) '---------------------------'

          call tem_abort('ABORTING')
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_MEM%pIDsort(pos)

        ! If the attached PID is > 0 it indicates this is a force contribution which
        ! must be ADDED to the total force.
        if( recv%buf_force(posOfProc)%val(iBuff)%pID > 0 ) then
          this%particles_MEM%val(pos)%Fbuff( this%particles_MEM%val(pos)%Fnow,1:6) &
            & = this%particles_MEM%val(pos)%Fbuff( this%particles_MEM%val(pos)%Fnow,1:6)  &
            &                      + recv%buf_force(posOfProc)%val(iBuff)%V(1:6)
        else
          ! Negative PID indicates that this is the complete force from the last
          ! time step. We store this in Fbuff(Flast,:) for use in averaging of F later
          write(logUnit(1),*) 'Got complete force from last owner!'
          write(logUnit(1),*) 'From proc ', recv%proc(iproc), ' PID = ', recvPID
          this%particles_MEM%val(pos)%Fbuff( this%particles_MEM%val(pos)%Flast,1:6) &
            & = recv%buf_force(posOfProc)%val(iBuff)%V(1:6)

        end if

      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr           ) ! error status
    end if

  end subroutine exchangeForces

  subroutine resetForceBuffers(this)
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    ! -------------------------------------------!
    integer :: iproc
    ! -------------------------------------------!
    ! Reset send buffer
    do iproc = 1, this%send%nProcs
      this%send%buf_force(iproc)%nParticles = 0
    end do

    ! Reset recv buffer
    do iproc = 1, this%recv%nProcs
      this%recv%buf_force(iproc)%nParticles = 0
    end do

  end subroutine resetForceBuffers

  !> addCollisionForceToBuffer adds DEM collision force Fcoll to the force buffer sent to
  !! particle with id particleID on rank recvRankIndex.
  subroutine addCollisionForceToBuffer( particleGroup, recvRankIndex, Fcoll, particleID )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Index of receiving rank in particleGroup%send%buf_force(recvRankIndex)
    integer :: recvRankIndex
    !> Collision force to add to buffer
    real(kind=rk), intent(in) :: Fcoll(3)
    !> ID of particle to which collision force should be applied
    integer, intent(in) :: particleID
    ! -------------------------------------------!
    integer :: nParticles
    ! -------------------------------------------!

    if( particleGroup%send%buf_force(recvRankIndex)%nParticles + 1 &
      & >= particleGroup%send%buf_force(recvRankIndex)%maxParticles ) then

      call tem_abort('Error addCollisionForceToBuffer: force buffer too small, aborting')
    end if

    ! Add collision force to buffer for owner of other particle
    particleGroup%send%buf_force(recvRankIndex)%nParticles &
                & = particleGroup%send%buf_force(recvRankIndex)%nParticles + 1

    nParticles = particleGroup%send%buf_force(recvRankIndex)%nParticles

    particleGroup%send%buf_force(recvRankIndex)%val(nParticles)%V(1:3) = Fcoll(1:3)

    particleGroup%send%buf_force(recvRankIndex)%val(nParticles)%V(4:6) = 0.0_rk

    particleGroup%send%buf_force(recvRankIndex)%val(nParticles)%pID = particleID

  end subroutine addCollisionForceToBuffer

  !> DEM_exchangeForces sends DEM collision force contributions from collisions that
  !! cannot be resolved by the receiving process. The send force buffers in this case must
  !! already be filled in the DEM collision-handling routines!
  subroutine DEM_exchangeForces_DPS(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iBuff
    integer :: proc, posOfProc, pos, recvcount
    integer :: recvPID
    ! logical :: err_flag_local, err_flag_global    ! for debugging only
    ! -------------------------------------------!
    ! err_flag_local = .FALSE.
    ! err_flag_global = .FALSE.

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                               &
       &      recv%buf_force( iproc )%val,          & ! me
       &      recv%buf_force( iproc )%maxParticles, & ! max particles me can receive
       &      mus_pIDvector_type,                   & ! data type
       &      recv%proc(iproc),                     & ! source process
       &      message_flag,                         & ! flag
       &      comm,                                 & ! communicator
       &      recv%rqhandle(iproc),                 & ! handle
       &      ierr )                                  ! error status
    end do

    !  start the sending of force contributions
    do iproc = 1, send%nprocs
      call mpi_isend(                                   &
       &      send%buf_force( iproc )%val,              & ! buffer
       &      send%buf_force( iproc )%nParticles,       & ! count
       &      mus_pIDvector_type,                       & ! data type
       &      send%proc(iproc),                         & ! target
       &      message_flag,                             & ! integer tag
       &      comm,                                     & ! communicator
       &      send%rqhandle( iproc ),                   & ! handle
       &      ierr )                                      ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)

    ! Now force contributions from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc),   & ! mpi status
                        & datatype = mus_pIDvector_type,   & ! data type
                        & count = recvcount,               & ! count
                        & ierror = ierr                    )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_force(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_force(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        recvPID = abs( recv%buf_force(posOfProc)%val(iBuff)%pID )
        pos = sortPosofval_particle_DPS( me  = this%particles_DPS,  &
                                        & pID = recvPID             )
        if(pos == 0) then
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_DPS%pIDsort(pos)

        ! Add this DEM force contribution to the DEM force for the particle on this process
        this%particles_DPS%val(pos)%F_DEM( this%particles_DPS%val(pos)%F_DEM_next,1:6) &
          & = this%particles_DPS%val(pos)%F_DEM( this%particles_DPS%val(pos)%F_DEM_next,1:6)  &
          &                      + recv%buf_force(posOfProc)%val(iBuff)%V(1:6)
      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr           ) ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)

    ! Check for global errors
    ! call mus_particles_global_errorcheck(err_flag_global, MPI_COMM_WORLD)
    ! if(err_flag_global) then
    !   ! If ANY process encountered an error, print the particleGroup to
    !   ! the processes pgDebugLog.
    !   open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
    !     write(pgDebugLog%lu,*) 'DEM_exchangeForcesDPS encountered an error on a process'
    !     write(pgDebugLog%lu,*) 'local error flag: ', err_flag_local
    !     call printParticleGroup2_DPS(this, pgDebugLog%lu, myRank, 0)
    !   close( pgDebugLog%lu )
    !   call tem_abort("ABORTING")
    ! end if

  end subroutine DEM_exchangeForces_DPS


  !> DEM_exchangeForces sends DEM collision force contributions from collisions that
  !! cannot be resolved by the receiving process. The send force buffers in this case must
  !! already be filled in the DEM collision-handling routines!
  subroutine DEM_exchangeForces(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iBuff
    integer :: proc, posOfProc, pos, recvcount
    integer :: recvPID
    ! -------------------------------------------!

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                               &
       &      recv%buf_force( iproc )%val,          & ! me
       &      recv%buf_force( iproc )%maxParticles, & ! max particles me can receive
       &      mus_pIDvector_type,                   & ! data type
       &      recv%proc(iproc),                     & ! source process
       &      message_flag,                         & ! flag
       &      comm,                                 & ! communicator
       &      recv%rqhandle(iproc),                 & ! handle
       &      ierr )                                  ! error status
    end do

    !  start the sending of force contributions
    do iproc = 1, send%nprocs
      call mpi_isend(                                   &
       &      send%buf_force( iproc )%val,              & ! buffer
       &      send%buf_force( iproc )%nParticles,       & ! count
       &      mus_pIDvector_type,                       & ! data type
       &      send%proc(iproc),                         & ! target
       &      message_flag,                             & ! integer tag
       &      comm,                                     & ! communicator
       &      send%rqhandle( iproc ),                   & ! handle
       &      ierr )                                      ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! Now force contributions from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc),   & ! mpi status
                        & datatype = mus_pIDvector_type,   & ! data type
                        & count = recvcount,               & ! count
                        & ierror = ierr                    )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_force(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_force(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        recvPID = abs( recv%buf_force(posOfProc)%val(iBuff)%pID )
        pos = sortPosofval_particle_MEM( me  = this%particles_MEM,                      &
                                   & pID = recvPID )
        if (pos == 0) then
          write(logUnit(1),*) '---------------------------'
          write(logUnit(1),*) &
            & 'WARNING DEM_exchangeForces: cannot find particle with ID', recv%buf_force(iProc)%val(iBuff)%pID
          write(logUnit(1),*) 'in particleGroup on rank = ', myRank
          write(logUnit(1),*) '--- RECV buff: ---'
          call print_particles_comm(recv)
          write(logUnit(1),*) '--- Rank ', myRank, ' particle group: ---'
          call printParticleGroup(this, logUnit(1))
          write(logUnit(1),*) '--- PIDlists ---'
          call printPIDlist(this)
          write(logUnit(1),*) '---------------------------'

          call tem_abort('ABORTING')
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_MEM%pIDsort(pos)

        ! Add this DEM force contribution to the DEM force for the particle on this process
        this%particles_MEM%val(pos)%F_DEM( this%particles_MEM%val(pos)%F_DEM_next,1:6) &
          & = this%particles_MEM%val(pos)%F_DEM( this%particles_MEM%val(pos)%F_DEM_next,1:6)  &
          &                      + recv%buf_force(posOfProc)%val(iBuff)%V(1:6)

        ! write(logUnit(1),*) 'Proc ', myRank, 'received DEM force contribution', &
        !   & recv%buf_force(posOfProc)%val(iBuff)%V(1:6)




      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine DEM_exchangeForces

  !> DEM_exchangeWallPositions sends wall position sums close to a particle as detected on each process
  !! to the particle owner, which averages the results to get one average wall position
  !! Upon termination of this routine particle%rwall for all particles which I own should be set
  !! correctly.
  subroutine DEM_exchangeWallPositions(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iBuff, iParticle
    integer :: proc, posOfProc, pos, recvcount
    integer :: recvPID
    real(kind=rk) :: avgWallPos
    ! -------------------------------------------!
    avgWallPos = 0.0_rk

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                               &
       &      recv%buf_wall( iproc )%val,          & ! me
       &      recv%buf_wall( iproc )%maxParticles, & ! max particles me can receive
       &      mus_wallPos_type,                   & ! data type
       &      recv%proc(iproc),                     & ! source process
       &      message_flag,                         & ! flag
       &      comm,                                 & ! communicator
       &      recv%rqhandle(iproc),                 & ! handle
       &      ierr )                                  ! error status
    end do

    !  start the sending of force contributions
    do iproc = 1, send%nprocs
      call mpi_isend(                                   &
       &      send%buf_wall( iproc )%val,              & ! buffer
       &      send%buf_wall( iproc )%nParticles,       & ! count
       &      mus_wallPos_type,                       & ! data type
       &      send%proc(iproc),                         & ! target
       &      message_flag,                             & ! integer tag
       &      comm,                                     & ! communicator
       &      send%rqhandle( iproc ),                   & ! handle
       &      ierr )                                      ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! Now force contributions from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc),   & ! mpi status
                        & datatype = mus_wallPos_type,   & ! data type
                        & count = recvcount,               & ! count
                        & ierror = ierr                    )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_wall(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_wall(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        recvPID = abs( recv%buf_wall(posOfProc)%val(iBuff)%I(2) )
        pos = sortPosofval_particle_MEM( me  = this%particles_MEM,                      &
                                   & pID = recvPID )
        if (pos == 0) then
          write(logUnit(1),*) '---------------------------'
          write(logUnit(1),*) &
            & 'WARNING DEM_exchangeWallPositions: cannot find particle with ID', recv%buf_wall(iProc)%val(iBuff)%I(2)
          write(logUnit(1),*) 'in particleGroup on rank = ', myRank
          write(logUnit(1),*) '--- RECV buff: ---'
          call print_particles_comm(recv)
          write(logUnit(1),*) '--- Rank ', myRank, ' particle group: ---'
          call printParticleGroup(this, logUnit(1))
          write(logUnit(1),*) '--- PIDlists ---'
          call printPIDlist(this)
          write(logUnit(1),*) '---------------------------'

          call tem_abort('ABORTING')
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_MEM%pIDsort(pos)

        ! Increment number of elements in total wallPosSum
        this%particles_MEM%val(pos)%nWallPos = this%particles_MEM%val(pos)%nWallPos + recv%buf_wall(posOfProc)%val(iBuff)%I(1)

        ! Add the wallPosSum received from other process to the total wallPosSum
        this%particles_MEM%val(pos)%avgWallPos = this%particles_MEM%val(pos)%avgWallPos + recv%buf_wall(posOfProc)%val(iBuff)%V(1:3)

      end do ! iBuff
    end do statusLoop

    ! Right now averageWallPos is set to the sum of all the local wallPosSums for all particles this process owns.
    ! Now we divide by the total (global) amount of elements in this sum to set the averageWallPos.
    do iParticle = 1, this%particles_MEM%nvals
      if( this%particles_MEM%val(iParticle)%owner == myRank    &
        & .AND. this%particles_MEM%val(iParticle)%nWallPos > 0 ) then
        this%particles_MEM%val(iParticle)%avgWallPos &
          & = this%particles_MEM%val(iParticle)%avgWallPos / this%particles_MEM%val(iParticle)%nWallPos

        ! Compute vector pointing from particle to wall
        this%particles_MEM%val(iParticle)%rwall = this%particles_MEM%val(iParticle)%avgWallPos &
          &                                     - this%particles_MEM%val(iParticle)%pos(1:3)
      end if
    end do

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine DEM_exchangeWallPositions

  !> DEM_exchangeWallPositions sends wall position sums close to a particle as detected on each process
  !! to the particle owner, which averages the results to get one average wall position
  !! Upon termination of this routine particle%rwall for all particles which I own should be set
  !! correctly.
  subroutine DEM_exchangeWallPositions_DPS(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iBuff, iParticle
    integer :: proc, posOfProc, pos, recvcount
    integer :: recvPID
    real(kind=rk) :: avgWallPos
    ! -------------------------------------------!
    avgWallPos = 0.0_rk

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                               &
       &      recv%buf_wall( iproc )%val,          & ! me
       &      recv%buf_wall( iproc )%maxParticles, & ! max particles me can receive
       &      mus_wallPos_type,                   & ! data type
       &      recv%proc(iproc),                     & ! source process
       &      message_flag,                         & ! flag
       &      comm,                                 & ! communicator
       &      recv%rqhandle(iproc),                 & ! handle
       &      ierr )                                  ! error status
    end do

    !  start the sending of force contributions
    do iproc = 1, send%nprocs
      call mpi_isend(                                   &
       &      send%buf_wall( iproc )%val,              & ! buffer
       &      send%buf_wall( iproc )%nParticles,       & ! count
       &      mus_wallPos_type,                       & ! data type
       &      send%proc(iproc),                         & ! target
       &      message_flag,                             & ! integer tag
       &      comm,                                     & ! communicator
       &      send%rqhandle( iproc ),                   & ! handle
       &      ierr )                                      ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! Now force contributions from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc),   & ! mpi status
                        & datatype = mus_wallPos_type,   & ! data type
                        & count = recvcount,               & ! count
                        & ierror = ierr                    )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_wall(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_wall(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        recvPID = abs( recv%buf_wall(posOfProc)%val(iBuff)%I(2) )
        pos = sortPosofval_particle_DPS( me  = this%particles_DPS,                      &
                                   & pID = recvPID )
        if(pos == 0) then
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_DPS%pIDsort(pos)

        ! Increment number of elements in total wallPosSum
        this%particles_DPS%val(pos)%nWallPos = this%particles_DPS%val(pos)%nWallPos + recv%buf_wall(posOfProc)%val(iBuff)%I(1)

        ! Add the wallPosSum received from other process to the total wallPosSum
        this%particles_DPS%val(pos)%avgWallPos = this%particles_DPS%val(pos)%avgWallPos + recv%buf_wall(posOfProc)%val(iBuff)%V(1:3)
        write(logUnit(1),*) 'Added wallPosSum from proc ', recv%proc(posOfProc)

      end do ! iBuff
    end do statusLoop

    ! Right now averageWallPos is set to the sum of all the local wallPosSums for all particles this process owns.
    ! Now we divide by the total (global) amount of elements in this sum to set the averageWallPos.
    do iParticle = 1, this%particles_DPS%nvals
      if( this%particles_DPS%val(iParticle)%owner == myRank    &
        & .AND. this%particles_DPS%val(iParticle)%nWallPos > 0 ) then

        this%particles_DPS%val(iParticle)%avgWallPos &
          & = this%particles_DPS%val(iParticle)%avgWallPos / this%particles_DPS%val(iParticle)%nWallPos

        ! ! Compute vector pointing from particle to wall
        this%particles_DPS%val(iParticle)%rwall = this%particles_DPS%val(iParticle)%avgWallPos &
          &                                     - this%particles_DPS%val(iParticle)%pos(1:3)
      end if
    end do

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine DEM_exchangeWallPositions_DPS

  !> exchangePositions exchanges particles continuous positions after they
  !! have been updated by the particle owner. Each process sends position updates
  !! for particles they own and receives updates for particles which they do not
  !! own, but that do exist in their particleGroup.
  subroutine exchangePositions(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_pos(iproc)%nParticles = 0
    end do ! iproc

    do iParticle = 1, this%particles_MEM%nvals
      ! If I do not own this particle, I should not send velocity updates
      if(this%particles_MEM%val(iParticle)%owner /= myRank) cycle

      ! If I DO own this particle, check the existsOnProc mask tells us
      ! which processes need velocity updates for this particle
      do iproc = 1, send%nProcs
        if( this%particles_MEM%val(iParticle)%existsOnProc(iproc) ) then

          ! Check if buffer is large enough for the new contribution
          if( send%buf_pos(iproc)%nParticles+1                   &
            &               >= send%buf_pos(iproc)%maxParticles  ) then
            call tem_abort('Error exchangePositions: buffer too small, aborting')
          end if

          ! Increment number of particle contributions in force buffer
          send%buf_pos(iproc)%nParticles = send%buf_pos(iproc)%nParticles + 1

          ! Add velocity update to end of the buffer
          send%buf_pos(iproc)%val( send%buf_pos(iproc)%nParticles )%V &
            & = this%particles_MEM%val(iParticle)%pos

          ! Add the ID of the particle to which this velocity update pertains
          send%buf_pos(iproc)%val( send%buf_pos(iproc)%nParticles )%I(1:4) &
            & = this%particles_MEM%val(iParticle)%coordOfOrigin(1:4)

          ! Add the ID of the particle to which this velocity update pertains
          send%buf_pos(iproc)%val( send%buf_pos(iproc)%nParticles )%I(5) &
            & = this%particles_MEM%val(iParticle)%particleID

        end if ! existsOnProc

      end do ! iproc
    end do ! iParticle

  !  write(logUnit(1),*) '--- EXCHANGE VELOCTIES RANK 1 SENDBUFFER ---'
  !  call print_particles_comm( send )

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                             &
       &      recv%buf_pos( iproc )%val,          & ! me
       &      recv%buf_pos( iproc )%maxParticles, & ! max particles me can receive
       &      mus_positionUpdate_type,                 & ! data type
       &      recv%proc(iproc),                   & ! source process
       &      message_flag,                       & ! flag
       &      comm,                               & ! communicator
       &      recv%rqhandle(iproc),               & ! handle
       &      ierr                                ) ! error status
    end do

    !  start the sending of velocity updates
    do iproc = 1, send%nprocs
      call mpi_isend(                            &
       &      send%buf_pos( iproc )%val,         & ! buffer
       &      send%buf_pos( iproc )%nParticles,  & ! count
       &      mus_positionUpdate_type,                & ! data type
       &      send%proc(iproc),                  & ! target
       &      message_flag,                      & ! integer tag
       &      comm,                              & ! communicator
       &      send%rqhandle( iproc ),            & ! handle
       &      ierr                               ) ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! Now force contributions from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc), & ! mpi status
                        & datatype = mus_positionUpdate_type, & ! data type
                        & count = recvcount,             & ! count
                        & ierror = ierr                  )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_pos(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_pos(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        pos = sortPosofval_particle_MEM( me  = this%particles_MEM,                      &
                                   & pID = recv%buf_pos(posOfProc)%val(iBuff)%I(5) )
        if (pos == 0) then
          write(logUnit(1),*) '---------------------------'
          write(logUnit(1),*) &
            & 'WARNING exchangePositions: cannot find particle with ID', &
            & recv%buf_pos(iProc)%val(iBuff)%I(5)

          write(logUnit(1),*) 'in particleGroup on rank = ', myRank
          write(logUnit(1),*) '--- RECV buff: ---'
          call print_particles_comm(recv)
          write(logUnit(1),*) '--- Rank ', myRank, ' particle group: ---'
          call printParticleGroup(this, logUnit(1))
          write(logUnit(1),*) '--- PIDlists ---'
          call printPIDlist(this)
          write(logUnit(1),*) '---------------------------'

          call tem_abort('ABORTING')
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_MEM%pIDsort(pos)

        this%particles_MEM%val(pos)%pos(1:6) =                               &
          &                      recv%buf_pos(posOfProc)%val(iBuff)%V(1:6)

        this%particles_MEM%val(pos)%coordOfOrigin(1:4) =                               &
          &                      recv%buf_pos(posOfProc)%val(iBuff)%I(1:4)


      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine exchangePositions


  !> exchangePositions exchanges particles continuous positions after they
  !! have been updated by the particle owner. Each process sends position updates
  !! for particles they own and receives updates for particles which they do not
  !! own, but that do exist in their particleGroup.
  subroutine exchangePositions_DPS(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_pos(iproc)%nParticles = 0
    end do ! iproc

    do iParticle = 1, this%particles_DPS%nvals
      ! If I do not own this particle, I should not send velocity updates
      if(this%particles_DPS%val(iParticle)%owner /= myRank) cycle

      ! If I DO own this particle, check the existsOnProc mask tells us
      ! which processes need velocity updates for this particle
      do iproc = 1, send%nProcs
        if( this%particles_DPS%val(iParticle)%existsOnProc(iproc) ) then

          ! Check if buffer is large enough for the new contribution
          if( send%buf_pos(iproc)%nParticles+1                   &
            &               >= send%buf_pos(iproc)%maxParticles  ) then
            call tem_abort('Error exchangePositions: buffer too small, aborting')
          end if

          ! Increment number of particle contributions in force buffer
          send%buf_pos(iproc)%nParticles = send%buf_pos(iproc)%nParticles + 1

          ! Add velocity update to end of the buffer
          send%buf_pos(iproc)%val( send%buf_pos(iproc)%nParticles )%V &
            & = this%particles_DPS%val(iParticle)%pos

          ! Add the ID of the particle to which this velocity update pertains
          send%buf_pos(iproc)%val( send%buf_pos(iproc)%nParticles )%I(1:4) &
            & = this%particles_DPS%val(iParticle)%coordOfOrigin(1:4)

          ! Add the ID of the particle to which this velocity update pertains
          send%buf_pos(iproc)%val( send%buf_pos(iproc)%nParticles )%I(5) &
            & = this%particles_DPS%val(iParticle)%particleID

        end if ! existsOnProc

      end do ! iproc
    end do ! iParticle

  !  write(logUnit(1),*) '--- EXCHANGE VELOCTIES RANK 1 SENDBUFFER ---'
  !  call print_particles_comm( send )

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                             &
       &      recv%buf_pos( iproc )%val,          & ! me
       &      recv%buf_pos( iproc )%maxParticles, & ! max particles me can receive
       &      mus_positionUpdate_type,                 & ! data type
       &      recv%proc(iproc),                   & ! source process
       &      message_flag,                       & ! flag
       &      comm,                               & ! communicator
       &      recv%rqhandle(iproc),               & ! handle
       &      ierr                                ) ! error status
    end do

    !  start the sending of velocity updates
    do iproc = 1, send%nprocs
      call mpi_isend(                            &
       &      send%buf_pos( iproc )%val,         & ! buffer
       &      send%buf_pos( iproc )%nParticles,  & ! count
       &      mus_positionUpdate_type,                & ! data type
       &      send%proc(iproc),                  & ! target
       &      message_flag,                      & ! integer tag
       &      comm,                              & ! communicator
       &      send%rqhandle( iproc ),            & ! handle
       &      ierr                               ) ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc), & ! mpi status
                        & datatype = mus_positionUpdate_type, & ! data type
                        & count = recvcount,             & ! count
                        & ierror = ierr                  )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_pos(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_pos(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        pos = sortPosofval_particle_DPS( me  = this%particles_DPS,                 &
                                   & pID = recv%buf_pos(posOfProc)%val(iBuff)%I(5) )
        if (pos == 0) then
          ! write(logUnit(1),*) '---------------------------'
          ! write(logUnit(1),*) &
          !   & 'WARNING exchangePositions_DPS: cannot find particle with ID', &
          !   & recv%buf_pos(iProc)%val(iBuff)%I(5)

          ! write(logUnit(1),*) 'in particleGroup on rank = ', myRank
          ! write(logUnit(1),*) '--- RECV buff: ---'
          ! call print_particles_comm(recv)
          ! write(logUnit(1),*) '---------------------------'
          ! call tem_abort('ABORTING')
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_DPS%pIDsort(pos)

        this%particles_DPS%val(pos)%pos(1:6) =                               &
          &                      recv%buf_pos(posOfProc)%val(iBuff)%V(1:6)

        this%particles_DPS%val(pos)%coordOfOrigin(1:4) =                               &
          &                      recv%buf_pos(posOfProc)%val(iBuff)%I(1:4)


      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)

  end subroutine exchangePositions_DPS

  !> exchangeVelocities exchanges particles continuous velocities after they
  !! have been updated by the particle owner. Each process sends velocity updates
  !! for particles they own and receives updates for particles which they do not
  !! own, but that do exist in their particleGroup.
  subroutine exchangeVelocities(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_vec(iproc)%nParticles = 0
    end do ! iproc

    do iParticle = 1, this%particles_MEM%nvals
      ! If I do not own this particle, I should not send velocity updates
      if(this%particles_MEM%val(iParticle)%owner /= myRank) cycle

      ! If I DO own this particle, check the existsOnProc mask tells us
      ! which processes need velocity updates for this particle
      do iproc = 1, send%nProcs
        if( this%particles_MEM%val(iParticle)%existsOnProc(iproc) ) then

          ! Check if buffer is large enough for the new contribution
          if( send%buf_vec(iproc)%nParticles+1                   &
            &               >= send%buf_vec(iproc)%maxParticles  ) then
            call tem_abort('Error exchangeVelocities: buffer too small, aborting')
          end if

          ! Increment number of particle contributions in force buffer
          send%buf_vec(iproc)%nParticles = send%buf_vec(iproc)%nParticles + 1

          ! Add velocity update to end of the buffer
          send%buf_vec(iproc)%val( send%buf_vec(iproc)%nParticles )%V &
            & = this%particles_MEM%val(iParticle)%vel

          ! Add the ID of the particle to which this velocity update pertains
          send%buf_vec(iproc)%val( send%buf_vec(iproc)%nParticles )%pID &
            & = this%particles_MEM%val(iParticle)%particleID

        end if ! existsOnProc

      end do ! iproc
    end do ! iParticle

  !  write(logUnit(1),*) '--- EXCHANGE VELOCTIES RANK 1 SENDBUFFER ---'
  !  call print_particles_comm( send )

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                             &
       &      recv%buf_vec( iproc )%val,          & ! me
       &      recv%buf_vec( iproc )%maxParticles, & ! max particles me can receive
       &      mus_pIDvector_type,                 & ! data type
       &      recv%proc(iproc),                   & ! source process
       &      message_flag,                       & ! flag
       &      comm,                               & ! communicator
       &      recv%rqhandle(iproc),               & ! handle
       &      ierr                                ) ! error status
    end do

    !  start the sending of velocity updates
    do iproc = 1, send%nprocs
      call mpi_isend(                            &
       &      send%buf_vec( iproc )%val,         & ! buffer
       &      send%buf_vec( iproc )%nParticles,  & ! count
       &      mus_pIDvector_type,                & ! data type
       &      send%proc(iproc),                  & ! target
       &      message_flag,                      & ! integer tag
       &      comm,                              & ! communicator
       &      send%rqhandle( iproc ),            & ! handle
       &      ierr                               ) ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! Now force contributions from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc), & ! mpi status
                        & datatype = mus_pIDvector_type, & ! data type
                        & count = recvcount,             & ! count
                        & ierror = ierr                  )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_vec(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_vec(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        pos = sortPosofval_particle_MEM( me  = this%particles_MEM,                      &
                                   & pID = recv%buf_vec(posOfProc)%val(iBuff)%pID )
        if (pos == 0) then
          write(logUnit(1),*) '---------------------------'
          write(logUnit(1),*) &
            & 'WARNING exchangeVelocities: cannot find particle with ID', &
            & recv%buf_vec(iProc)%val(iBuff)%pID

          write(logUnit(1),*) 'in particleGroup on rank = ', myRank
          write(logUnit(1),*) '--- RECV buff: ---'
          call print_particles_comm(recv)
          write(logUnit(1),*) '--- Rank ', myRank, ' particle group: ---'
          call printParticleGroup(this, logUnit(1))
          write(logUnit(1),*) '--- PIDlists ---'
          call printPIDlist(this)
          write(logUnit(1),*) '---------------------------'

          call tem_abort('ABORTING')
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_MEM%pIDsort(pos)

        this%particles_MEM%val(pos)%vel(1:6) =                               &
          &                      recv%buf_vec(posOfProc)%val(iBuff)%V(1:6)


      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine exchangeVelocities

  !> exchangeVelocities exchanges particles continuous velocities after they
  !! have been updated by the particle owner. Each process sends velocity updates
  !! for particles they own and receives updates for particles which they do not
  !! own, but that do exist in their particleGroup.
  subroutine exchangeVelocities_DPS(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_vec(iproc)%nParticles = 0
    end do ! iproc

    do iParticle = 1, this%particles_DPS%nvals
      ! If I do not own this particle, I should not send velocity updates
      if(this%particles_DPS%val(iParticle)%owner /= myRank) cycle

      ! If I DO own this particle, check the existsOnProc mask tells us
      ! which processes need velocity updates for this particle
      do iproc = 1, send%nProcs
        if( this%particles_DPS%val(iParticle)%existsOnProc(iproc) ) then

          ! Check if buffer is large enough for the new contribution
          if( send%buf_vec(iproc)%nParticles+1                   &
            &               >= send%buf_vec(iproc)%maxParticles  ) then
            call tem_abort('Error exchangeVelocities: buffer too small, aborting')
          end if

          ! Increment number of particle contributions in force buffer
          send%buf_vec(iproc)%nParticles = send%buf_vec(iproc)%nParticles + 1

          ! Add velocity update to end of the buffer
          send%buf_vec(iproc)%val( send%buf_vec(iproc)%nParticles )%V &
            & = this%particles_DPS%val(iParticle)%vel

          ! Add the ID of the particle to which this velocity update pertains
          send%buf_vec(iproc)%val( send%buf_vec(iproc)%nParticles )%pID &
            & = this%particles_DPS%val(iParticle)%particleID

        end if ! existsOnProc

      end do ! iproc
    end do ! iParticle

  !  write(logUnit(1),*) '--- EXCHANGE VELOCTIES RANK 1 SENDBUFFER ---'
  !  call print_particles_comm( send )

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                             &
       &      recv%buf_vec( iproc )%val,          & ! me
       &      recv%buf_vec( iproc )%maxParticles, & ! max particles me can receive
       &      mus_pIDvector_type,                 & ! data type
       &      recv%proc(iproc),                   & ! source process
       &      message_flag,                       & ! flag
       &      comm,                               & ! communicator
       &      recv%rqhandle(iproc),               & ! handle
       &      ierr                                ) ! error status
    end do

    !  start the sending of velocity updates
    do iproc = 1, send%nprocs
      call mpi_isend(                            &
       &      send%buf_vec( iproc )%val,         & ! buffer
       &      send%buf_vec( iproc )%nParticles,  & ! count
       &      mus_pIDvector_type,                & ! data type
       &      send%proc(iproc),                  & ! target
       &      message_flag,                      & ! integer tag
       &      comm,                              & ! communicator
       &      send%rqhandle( iproc ),            & ! handle
       &      ierr                               ) ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)

    ! Now force contributions from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc), & ! mpi status
                        & datatype = mus_pIDvector_type, & ! data type
                        & count = recvcount,             & ! count
                        & ierror = ierr                  )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_vec(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_vec(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        pos = sortPosofval_particle_DPS( me  = this%particles_DPS,                    &
                                       & pID = recv%buf_vec(posOfProc)%val(iBuff)%pID )
        if(pos == 0) then
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_DPS%pIDsort(pos)

        this%particles_DPS%val(pos)%vel(1:6) =                               &
          &                      recv%buf_vec(posOfProc)%val(iBuff)%V(1:6)


      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)

  end subroutine exchangeVelocities_DPS

  !> exchangeHydroForces_DPS exchanges particles hydrodynamic forces after they
  !! have been updated by the particle owner. Each process sends force updates
  !! for particles they own and receives updates for particles which they do not
  !! own, but that do exist in their particleGroup.
  subroutine exchangeHydroForces_DPS(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_force(iproc)%nParticles = 0
    end do ! iproc

    do iParticle = 1, this%particles_DPS%nvals
      ! If I do not own this particle, I should not send hydro force updates
      if(this%particles_DPS%val(iParticle)%owner /= myRank) cycle

      ! If I DO own this particle, check the existsOnProc mask tells us
      ! which processes need hydro force updates for this particle
      do iproc = 1, send%nProcs
        if( this%particles_DPS%val(iParticle)%existsOnProc(iproc) ) then

          ! Check if buffer is large enough for the new contribution
          if( send%buf_force(iproc)%nParticles+1                   &
            &               >= send%buf_force(iproc)%maxParticles  ) then
            call tem_abort('Error exchangeHydroForces: buffer too small, aborting')
          end if

          ! Increment number of particle contributions in force buffer
          send%buf_force(iproc)%nParticles = send%buf_force(iproc)%nParticles + 1

          ! Add hydro force update to end of the buffer
          send%buf_force(iproc)%val( send%buf_force(iproc)%nParticles )%V &
            & = this%particles_DPS%val(iParticle)%F

          ! Add the ID of the particle to which this force update pertains
          send%buf_force(iproc)%val( send%buf_force(iproc)%nParticles )%pID &
            & = this%particles_DPS%val(iParticle)%particleID

        end if ! existsOnProc

      end do ! iproc
    end do ! iParticle

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                             &
       &      recv%buf_force( iproc )%val,          & ! me
       &      recv%buf_force( iproc )%maxParticles, & ! max particles me can receive
       &      mus_pIDvector_type,                 & ! data type
       &      recv%proc(iproc),                   & ! source process
       &      message_flag,                       & ! flag
       &      comm,                               & ! communicator
       &      recv%rqhandle(iproc),               & ! handle
       &      ierr                                ) ! error status
    end do

    !  start the sending of velocity updates
    do iproc = 1, send%nprocs
      call mpi_isend(                            &
       &      send%buf_force( iproc )%val,         & ! buffer
       &      send%buf_force( iproc )%nParticles,  & ! count
       &      mus_pIDvector_type,                & ! data type
       &      send%proc(iproc),                  & ! target
       &      message_flag,                      & ! integer tag
       &      comm,                              & ! communicator
       &      send%rqhandle( iproc ),            & ! handle
       &      ierr                               ) ! error status
    end do ! iproc

    ! wait for receive buffer to be ready

    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)

    ! Now force contributions from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc), & ! mpi status
                        & datatype = mus_pIDvector_type, & ! data type
                        & count = recvcount,             & ! count
                        & ierror = ierr                  )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_force(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_force(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        pos = sortPosofval_particle_DPS( me  = this%particles_DPS,                    &
                                       & pID = recv%buf_force(posOfProc)%val(iBuff)%pID )
        if (pos == 0) then
          ! write(logUnit(1),*) '---------------------------'
          ! write(logUnit(1),*) &
          !   & 'WARNING exchangeVelocities: cannot find particle with ID', &
          !   & recv%buf_force(iProc)%val(iBuff)%pID

          ! write(logUnit(1),*) 'in particleGroup on rank = ', myRank
          ! write(logUnit(1),*) '--- RECV buff: ---'
          ! call print_particles_comm(recv)
          ! write(logUnit(1),*) '---------------------------'

          ! call tem_abort('ABORTING')
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_DPS%pIDsort(pos)

        this%particles_DPS%val(pos)%F(1:6) =                               &
          &                      recv%buf_force(posOfProc)%val(iBuff)%V(1:6)


      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)

  end subroutine exchangeHydroForces_DPS

  !> exchangeParticleStates is used to exchange particle (continuous) position
  !! and velocity as well as the integer coordOfOrigin. This is necessary for
  !! (for example) the collision handling routines
  subroutine exchangeParticleStates(this, send, recv, comm, myRank, message_flag)
    !> particleGroup of this process
    class(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_state(iproc)%nParticles = 0
    end do ! iproc

    do iParticle = 1, this%particles_MEM%nvals
      ! I should send state updates if...
      ! * I own the particle
      ! * I do not own the particle, but just resolved a collision involving
      !   the particle
      ! if(this%particles%val(iParticle)%owner /= myRank &
      !  & .AND. (.NOT. this%particles%val(iParticle)%hasCollided) ) cycle
      if(this%particles_MEM%val(iParticle)%owner /= myRank ) cycle

      ! If I DO own this particle, check the existsOnProc mask tells us
      ! which processes need state updates for this particle
      ! If removeFromProc = TRUE we also need to send one last update so
      ! that the ExclusionList updated and particle removed on that proc
      do iproc = 1, send%nProcs
        ! if( this%particles%val(iParticle)%existsOnProc(iproc) &
        !   & .OR. this%particles%val(iParticle)%removeFromProc(iproc) ) then
        if( this%particles_MEM%val(iParticle)%existsOnProc(iproc) ) then

          ! Check if buffer is large enough for the new contribution
          if( send%buf_state(iproc)%nParticles+1                   &
            &               >= send%buf_state(iproc)%maxParticles  ) then
            call tem_abort('Error exchangeParticleStates: state buffer too small, aborting')
          end if

          ! Increment number of particle contributions in force buffer
          send%buf_state(iproc)%nParticles = send%buf_state(iproc)%nParticles + 1

          ! Add state update to end of the buffer
          send%buf_state(iproc)%val( send%buf_state(iproc)%nParticles )%V(1:6) &
            & = this%particles_MEM%val(iParticle)%pos(1:6)

          send%buf_state(iproc)%val( send%buf_state(iproc)%nParticles )%V(7:12) &
            & = this%particles_MEM%val(iParticle)%vel(1:6)

          send%buf_state(iproc)%val( send%buf_state(iproc)%nParticles )%I(1:4) &
            & = this%particles_MEM%val(iParticle)%coordOfOrigin(1:4)

          ! We send the particle ID along at the end of the integer array!
          send%buf_state(iproc)%val( send%buf_state(iproc)%nParticles )%I(5) &
            & = this%particles_MEM%val(iParticle)%particleID

        end if ! existsOnProc

      end do ! iproc
    end do ! iParticle

    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                             &
       &      recv%buf_state( iproc )%val,          & ! me
       &      recv%buf_state( iproc )%maxParticles, & ! max particles me can receive
       &      mus_particleState_type,                 & ! data type
       &      recv%proc(iproc),                   & ! source process
       &      message_flag,                       & ! flag
       &      comm,                               & ! communicator
       &      recv%rqhandle(iproc),               & ! handle
       &      ierr                                ) ! error status
    end do

    !  start the sending of velocity updates
    do iproc = 1, send%nprocs
      call mpi_isend(                            &
       &      send%buf_state( iproc )%val,         & ! buffer
       &      send%buf_state( iproc )%nParticles,  & ! count
       &      mus_particleState_type,                & ! data type
       &      send%proc(iproc),                  & ! target
       &      message_flag,                      & ! integer tag
       &      comm,                              & ! communicator
       &      send%rqhandle( iproc ),            & ! handle
       &      ierr                               ) ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! Now particle state updates from other processes can be
    ! added to particles this process owns:

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of contributions actually received

      call mpi_get_count( status = recv_status(:,iproc),     & ! mpi status
                        & datatype = mus_particleState_type, & ! data type
                        & count = recvcount,                 & ! count
                        & ierror = ierr                      )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of force contributions in recv buffer
      recv%buf_state(posOfProc)%nParticles = recvcount

      ! Loop over received particles
      do iBuff = 1, recv%buf_state(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        pos = sortPosofval_particle_MEM( me  = this%particles_MEM,                      &
                                   & pID = recv%buf_state(posOfProc)%val(iBuff)%I(5) )
        if (pos == 0) then
          write(logUnit(1),*) '---------------------------'
          write(logUnit(1),*) &
            & 'WARNING exchangeParticleStates: cannot find particle with ID', &
            & recv%buf_state(iProc)%val(iBuff)%I(5)

            write(logUnit(1),*) 'received from proc ', proc
            write(logUnit(1),*) 'in particleGroup on rank = ', myRank
            write(logUnit(1),*) '--- RECV buff: ---'
            call print_particles_comm(recv)
            write(logUnit(1),*) '--- Rank ', myRank, ' particle group: ---'
            call printParticleGroup(this, logUnit(1))
            write(logUnit(1),*) '--- PIDlists ---'
            call printPIDlist(this)
            write(logUnit(1),*) '---------------------------'

            call tem_abort('ABORTING')
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_MEM%pIDsort(pos)

        ! Set the local particle data based on the update we received
        this%particles_MEM%val(pos)%pos(1:6) =                               &
          &                      recv%buf_state(posOfProc)%val(iBuff)%V(1:6)

        this%particles_MEM%val(pos)%vel(1:6) =                               &
          &                      recv%buf_state(posOfProc)%val(iBuff)%V(7:12)

        this%particles_MEM%val(pos)%coordOfOrigin(1:4) =                     &
          &                      recv%buf_state(posOfProc)%val(iBuff)%I(1:4)


      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine exchangeParticleStates

  !> exchangeNewParticles sends all the data needed to initialize a particle on
  !! a process. This is needed when a particle travels from one process to another.
  !! The data is used to add the particle to the receiving process particleGroup.
  !! This only creates a continuous representation of the particle on the receiving
  !! process. A subsequent call to a different routine then maps this to the discrete
  !! representation in the form of elements in the exclusionList.
  subroutine exchangeNewParticles_MEM(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    class(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle
    integer :: proc, posOfProc, recvcount

    type(mus_particle_MEM_type) :: tmpParticle
    logical :: wasAdded
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_particle(iproc)%nParticles = 0
    end do ! iproc

    ! Loop over all particles in particleGroup on this process
    do iParticle = 1, this%particles_MEM%nvals
      do iproc = 1, send%nProcs
        ! If particle is new on neighbor proc, add it to the send buffer
        if( this%particles_MEM%val(iParticle)%addToProc(iproc) )then

          ! Check if buffer is large enough for the new contribution
          if( send%buf_particle(iproc)%nParticles+1                   &
            &               >= send%buf_particle(iproc)%maxParticles  ) then
            call tem_abort('Error exchangeNewParticles: buffer too small, aborting')
          end if

          ! Increment number of particle contributions in buffer
          send%buf_particle(iproc)%nParticles = send%buf_particle(iproc)%nParticles + 1

          ! Add particle info to end of the buffer
          ! Floating point data:
          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(1:6) &
            & = this%particles_MEM%val(iParticle)%pos(1:6)

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(7:12) &
            & = this%particles_MEM%val(iParticle)%vel(1:6)

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(13:18) &
            & = this%particles_MEM%val(iParticle)%Fbuff(1,1:6)

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(19:24) &
            & = this%particles_MEM%val(iParticle)%Fbuff(2,1:6)

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(25) &
            & = this%particles_MEM%val(iParticle)%radius

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(26) &
            & = this%particles_MEM%val(iParticle)%mass

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(27:32) &
            & = this%particles_MEM%val(iParticle)%Fext(1:6)

          ! Integer data:
          ! coordinate of particle origin
          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%I(1:4) &
            & = this%particles_MEM%val(iParticle)%coordOfOrigin(1:4)

          ! Current position in forcebuffer = Fnow is the 5th element
          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%I(5) &
            & = this%particles_MEM%val(iParticle)%Fnow

          ! ID of the particle to which this force contribution pertains
          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%I(6) &
            & = this%particles_MEM%val(iParticle)%particleID

        end if

      end do ! iproc
    end do ! iParticle


    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                               &
       &      recv%buf_particle( iproc )%val,          & ! me
       &      recv%buf_particle( iproc )%maxParticles, & ! max particles me can receive
       &      mus_particleInfo_type,                   & ! data type
       &      recv%proc(iproc),                     & ! source process
       &      message_flag,                         & ! flag
       &      comm,                                 & ! communicator
       &      recv%rqhandle(iproc),                 & ! handle
       &      ierr )                                  ! error status
    end do

    !  start the sending of force contributions
    do iproc = 1, send%nprocs
      call mpi_isend(                                   &
       &      send%buf_particle( iproc )%val,              & ! buffer
       &      send%buf_particle( iproc )%nParticles,       & ! count
       &      mus_particleInfo_type,                       & ! data type
       &      send%proc(iproc),                         & ! target
       &      message_flag,                             & ! integer tag
       &      comm,                                     & ! communicator
       &      send%rqhandle( iproc ),                   & ! handle
       &      ierr )                                      ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of new particles  actually received

      call mpi_get_count( status = recv_status(:,iproc),   & ! mpi status
                        & datatype = mus_particleInfo_type,   & ! data type
                        & count = recvcount,               & ! count
                        & ierror = ierr                    )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of received particles in recv buffer
      recv%buf_particle(posOfProc)%nParticles = recvcount
    end do statusLoop

    do iproc = 1, recv%nprocs
      do iParticle = 1, recv%buf_particle(iproc)%nParticles
        tmpParticle%coordOfOrigin(1:4) = recv%buf_particle(iproc)%val(iParticle)%I(1:4)
        tmpParticle%Fnow = recv%buf_particle(iproc)%val(iParticle)%I(5)

        if(tmpParticle%Fnow == 1) then
          tmpParticle%Flast = 2
        else
          tmpParticle%Flast = 1
        end if

        tmpParticle%particleID = recv%buf_particle(iproc)%val(iParticle)%I(6)

        tmpParticle%pos(1:6) = recv%buf_particle(iproc)%val(iParticle)%V(1:6)
        tmpParticle%vel(1:6) = recv%buf_particle(iproc)%val(iParticle)%V(7:12)
        tmpParticle%Fbuff(1,1:6) = recv%buf_particle(iproc)%val(iParticle)%V(13:18)
        tmpParticle%Fbuff(2,1:6) = recv%buf_particle(iproc)%val(iParticle)%V(19:24)
        tmpParticle%F(1:6) = 0.5*( tmpParticle%Fbuff(1,1:6) + tmpParticle%Fbuff(2,1:6) )
        tmpParticle%radius = recv%buf_particle(iproc)%val(iParticle)%V(25)
        tmpParticle%mass = recv%buf_particle(iproc)%val(iParticle)%V(26)
        tmpParticle%Fext(1:6) = recv%buf_particle(iproc)%val(iParticle)%V(27:32)
        tmpParticle%rotInertia = 0.4 * tmpParticle%mass * tmpParticle%radius**2

        wasAdded = .FALSE.
        ! NOTE: change length to expand with to something sensible later
        call append_da_particle_MEM( me       = this%particles_MEM, &
                               & particle = tmpParticle,    &
                               & length   = 1,              &
                               & wasAdded = wasAdded        )

        if(wasAdded) then
          this%particles_MEM%val( this%particles_MEM%nvals )%newForMe = .TRUE.
        else
          ! It could be that multiple processes send a new particle to me.
          ! This is OK as append_da_particleMEM checks for uniqueness.
          ! In this case wasAdded will be false.
          ! write(logUnit(1),*) 'WARNING: exchangeNewParticles: rejected new particle'
          ! write(logUnit(1),*) 'Particle already in my rank ', myRank, ' particlegroup'
        end if

      end do ! iParticle
    end do ! iproc

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine exchangeNewParticles_MEM
  !> exchangeNewParticles sends all the data needed to initialize a particle on
  !! a process. This is needed when a particle travels from one process to another.
  !! The data is used to add the particle to the receiving process particleGroup.
  !! This only creates a continuous representation of the particle on the receiving
  !! process. A subsequent call to a different routine then maps this to the discrete
  !! representation in the form of elements in the exclusionList.
  subroutine exchangeNewParticles_DPS(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    class(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle
    integer :: proc, posOfProc, recvcount

    type(mus_particle_DPS_type) :: tmpParticle
    logical :: wasAdded
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_particle(iproc)%nParticles = 0
    end do ! iproc

    ! Loop over all particles in particleGroup on this process
    do iParticle = 1, this%particles_DPS%nvals
      do iproc = 1, send%nProcs
        ! If particle is new on neighbor proc, add it to the send buffer
        if( this%particles_DPS%val(iParticle)%addToProc(iproc) )then

          ! Check if buffer is large enough for the new contribution
          if( send%buf_particle(iproc)%nParticles+1                   &
            &               >= send%buf_particle(iproc)%maxParticles  ) then
            call tem_abort('Error exchangeNewParticles: buffer too small, aborting')
          end if

          ! Increment number of particle contributions in buffer
          send%buf_particle(iproc)%nParticles = send%buf_particle(iproc)%nParticles + 1

          ! Add particle info to end of the buffer
          ! Floating point data:
          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(1:6) &
            & = this%particles_DPS%val(iParticle)%pos(1:6)

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(7:12) &
            & = this%particles_DPS%val(iParticle)%vel(1:6)

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(13:18) &
            & = this%particles_DPS%val(iParticle)%F(1:6)

          ! For DPS particles we do not average over two time steps
          ! To use the same buffer as for MEM particles we fill the unneccesary spots
          ! with zeros.
          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(19:24) &
            & = 0.0_rk

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(25) &
            & = this%particles_DPS%val(iParticle)%radius

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(26) &
            & = this%particles_DPS%val(iParticle)%mass

          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%V(27:32) &
            & = this%particles_DPS%val(iParticle)%Fext(1:6)

          ! Integer data:
          ! coordinate of particle origin
          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%I(1:4) &
            & = this%particles_DPS%val(iParticle)%coordOfOrigin(1:4)

          ! Again we don't average over two time steps for DPS particles so this
          ! part of the buffer is irrelevant.
          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%I(5) &
            & = 0

          ! ID of the particle
          send%buf_particle(iproc)%val( send%buf_particle(iproc)%nParticles )%I(6) &
            & = this%particles_DPS%val(iParticle)%particleID

        end if

      end do ! iproc
    end do ! iParticle


    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                               &
       &      recv%buf_particle( iproc )%val,          & ! me
       &      recv%buf_particle( iproc )%maxParticles, & ! max particles me can receive
       &      mus_particleInfo_type,                   & ! data type
       &      recv%proc(iproc),                     & ! source process
       &      message_flag,                         & ! flag
       &      comm,                                 & ! communicator
       &      recv%rqhandle(iproc),                 & ! handle
       &      ierr )                                  ! error status
    end do

    !  start the sending of force contributions
    do iproc = 1, send%nprocs
      call mpi_isend(                                   &
       &      send%buf_particle( iproc )%val,              & ! buffer
       &      send%buf_particle( iproc )%nParticles,       & ! count
       &      mus_particleInfo_type,                       & ! data type
       &      send%proc(iproc),                         & ! target
       &      message_flag,                             & ! integer tag
       &      comm,                                     & ! communicator
       &      send%rqhandle( iproc ),                   & ! handle
       &      ierr )                                      ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if

    ! for each proc in MPI status (= a neighboring process)
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of new particles  actually received

      call mpi_get_count( status = recv_status(:,iproc),   & ! mpi status
                        & datatype = mus_particleInfo_type,   & ! data type
                        & count = recvcount,               & ! count
                        & ierror = ierr                    )

      ! Get the actual process this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      procLoop: do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit procLoop
      end do procLoop

      ! Update the number of received particles in recv buffer
      recv%buf_particle(posOfProc)%nParticles = recvcount
    end do statusLoop

    do iproc = 1, recv%nprocs
      do iParticle = 1, recv%buf_particle(iproc)%nParticles
        tmpParticle%coordOfOrigin(1:4) = recv%buf_particle(iproc)%val(iParticle)%I(1:4)
        tmpParticle%particleID = recv%buf_particle(iproc)%val(iParticle)%I(6)

        tmpParticle%pos(1:6) = recv%buf_particle(iproc)%val(iParticle)%V(1:6)
        tmpParticle%vel(1:6) = recv%buf_particle(iproc)%val(iParticle)%V(7:12)
        tmpParticle%F(1:6) = recv%buf_particle(iproc)%val(iParticle)%V(13:18)
        tmpParticle%radius = recv%buf_particle(iproc)%val(iParticle)%V(25)
        tmpParticle%mass = recv%buf_particle(iproc)%val(iParticle)%V(26)
        tmpParticle%Fext(1:6) = recv%buf_particle(iproc)%val(iParticle)%V(27:32)
        tmpParticle%rotInertia = 0.4 * tmpParticle%mass * tmpParticle%radius**2

        wasAdded = .FALSE.
        ! NOTE: change length to expand with to something sensible later
        call append_da_particle_DPS( me       = this%particles_DPS, &
                                   & particle = tmpParticle,        &
                                   & length   = 1,                  &
                                   & wasAdded = wasAdded            )

        if(wasAdded) then
          this%particles_DPS%val( this%particles_DPS%nvals )%newForMe = .TRUE.
        else
          ! It could be that multiple processes send a new particle to me.
          ! This is OK as append_da_particleMEM checks for uniqueness.
          ! In this case wasAdded will be false.
          ! write(logUnit(1),*) 'WARNING: exchangeNewParticles: rejected new particle'
          ! write(logUnit(1),*) 'Particle already in my rank ', myRank, ' particlegroup'
        end if

      end do ! iParticle
    end do ! iproc

    ! wait for send buffer to be ready
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr )           ! error status
    end if

  end subroutine exchangeNewParticles_DPS

  !> exchangeMomInc_DPS send the accumulated momentum transfer FROM particles TO fluid
  !! from a particle's previous owner to its new owner.
  subroutine exchangeMomInc_DPS(this, send, recv, comm, myRank, message_flag )
    !> particleGroup of this process
    type(mus_particle_group_type), intent(inout) :: this
    !> Communication type for sending force contributions
    type(mus_particles_communication_type), intent(inout) :: send
    !> Communication type for receiving force contributions
    type(mus_particles_communication_type), intent(inout) :: recv
    !> MPI communicator
    integer, intent(in) :: comm
    !> Rank of this process
    integer, intent(in) :: myRank
    !> Flag for message (in Musubi this is just iLevel, don't think we really
    !! need this here)
    integer, intent(in) :: message_flag
    ! -------------------------------------------!
    integer :: recv_status( mpi_status_size, recv%nprocs)
    integer :: send_status( mpi_status_size, send%nprocs)
    ! integer :: status( mpi_status_size, max(send%nprocs, recv%nprocs) )
    integer :: ierr             !< error flag
    integer :: iproc
    !> starting index for a particular force contribution
    integer :: iParticle, iBuff
    integer :: proc, posOfProc, pos, recvcount
    integer :: previousOwner, newOwner
    ! -------------------------------------------!

    ! ---  1: CONSTRUCT THE MESSAGES FOR EACH PROCESS --- !
    ! Reset nParticles in the force buffer for all neighbor procs
    do iproc = 1, send%nProcs
      send%buf_vec(iproc)%nParticles = 0
    end do ! iproc


    do iParticle = 1, this%particles_DPS%nvals
      previousOwner = this%particles_DPS%val(iParticle)%previousOwner
      newOwner = this%particles_DPS%val(iParticle)%owner

      if( (previousOwner /= newOwner) .AND. (previousOwner == myRank) ) then
        ! Add momInc to the buffer for the NEW owner of the particle
        ! Find the position of that process in our send%proc array
        do iproc = 1, send%nprocs
          if( send%proc(iproc) == newOwner) exit
        end do

        ! Check if buffer is large enough for the new contribution
        if( send%buf_vec(iproc)%nParticles+1                   &
          &               >= send%buf_vec(iproc)%maxParticles  ) then
          call tem_abort('Error exchangeParticleMomInc: buffer too small, aborting')
        end if

        ! Increment number of particle contributions in buffer
        send%buf_vec(iproc)%nParticles = send%buf_vec(iproc)%nParticles + 1

        ! Add velocity update to end of the buffer
        send%buf_vec(iproc)%val( send%buf_vec(iproc)%nParticles )%V(1:3) &
          & = this%particles_DPS%val(iParticle)%Favg(1:3)

        send%buf_vec(iproc)%val( send%buf_vec(iproc)%nParticles )%V(1:3) &
          & = 0.0_rk

        ! Add the ID of the particle to which this velocity update pertains
        send%buf_vec(iproc)%val( send%buf_vec(iproc)%nParticles )%pID &
          & = this%particles_DPS%val(iParticle)%particleID

      end if
    end do ! iParticle


    ! ---  2: SENDING AND RECEIVING MESSAGES  --- !
    ! Now the messages for all the processes should be complete.
    ! Start the receiving communications
    do iproc = 1, recv%nprocs
      ! start receive communications: one for each of our neighbor processes
      call mpi_irecv(                             &
       &      recv%buf_vec( iproc )%val,          & ! me
       &      recv%buf_vec( iproc )%maxParticles, & ! max particles me can receive
       &      mus_pIDvector_type,                 & ! data type
       &      recv%proc(iproc),                   & ! source process
       &      message_flag,                       & ! flag
       &      comm,                               & ! communicator
       &      recv%rqhandle(iproc),               & ! handle
       &      ierr                                ) ! error status
    end do

    !  start the sending of velocity updates
    do iproc = 1, send%nprocs
      call mpi_isend(                            &
       &      send%buf_vec( iproc )%val,         & ! buffer
       &      send%buf_vec( iproc )%nParticles,  & ! count
       &      mus_pIDvector_type,                & ! data type
       &      send%proc(iproc),                  & ! target
       &      message_flag,                      & ! integer tag
       &      comm,                              & ! communicator
       &      send%rqhandle( iproc ),            & ! handle
       &      ierr                               ) ! error status
    end do ! iproc

    ! wait for receive buffer to be ready
    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( recv%nprocs /= 0 ) then
      call mpi_waitall(recv%nprocs,               & ! count
        &              recv%rqhandle,             & ! request handles
        &              recv_status,               & ! mpi status
        &              ierr )                       ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)

    ! Update momInc data on my process using data received from other
    ! processes

    ! Loop over processes
    statusLoop: do iproc = 1, recv%nprocs
      ! get the number of momInc updates received
      ! from this proc
      call mpi_get_count( status = recv_status(:,iproc), & ! mpi status
                        & datatype = mus_pIDvector_type, & ! data type
                        & count = recvcount,             & ! count
                        & ierror = ierr                  )

      ! Get the MPI rank of the process that this count pertains to
      proc = recv_status( MPI_source, iproc )

      ! Find the position of that process in our recv%proc array
      do posOfProc = 1, recv%nprocs
        if( recv%proc(posOfProc) == proc) exit
      end do

      ! Update the number of momInc updates in recv buffer
      recv%buf_vec(posOfProc)%nParticles = recvcount

      ! Loop over received momInc updates
      do iBuff = 1, recv%buf_vec(posOfProc)%nParticles
        ! Find position of particle with this ID in my
        ! particle array.

        ! First get the sorted position (index of particle with this ID in pIDsort)
        pos = sortPosofval_particle_DPS( me  = this%particles_DPS,                    &
                                       & pID = recv%buf_vec(posOfProc)%val(iBuff)%pID )
        if(pos == 0) then
          cycle
        end if

        ! Convert to the index of corresponding particle in pIDlist
        pos = this%particles_DPS%pIDsort(pos)

        this%particles_DPS%val(pos)%Favg(1:3) =        &
          &    recv%buf_vec(posOfProc)%val(iBuff)%V(1:3)


      end do ! iBuff
    end do statusLoop

    ! wait for send buffer to be ready
    call tem_startTimer(timerHandle = mus_particle_timerHandles%idleTimer)
    if ( send%nprocs /= 0 ) then
      call mpi_waitall(send%nprocs,   & ! count
        &              send%rqhandle, & ! request handles
        &              send_status,   & ! mpi status
        &              ierr           ) ! error status
    end if
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%idleTimer)
  end subroutine exchangeMomInc_DPS

end module mus_particle_comm_module
