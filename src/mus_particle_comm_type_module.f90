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
!> In mus_particle_comm_type_module the communication types and routines for
!! particles are defined

module mus_particle_comm_type_module

  use mpi

  use env_module,                       only: rk, rk_mpi, long_k
  use tem_aux_module,                   only: tem_abort
  use tem_logging_module,               only: logUnit
  use tem_topology_module,              only: tem_firstIdAtLevel, &
    &                                         tem_coordOfId,      &
    &                                         tem_IdOfCoord
  use tem_dyn_array_module,             only: init, append, destroy, &
    &                                         empty, dyn_intArray_type
  use tem_property_module,              only: prp_hasBnd, prp_hasRemoteNgh, &
    &                                         prp_sendHalo

  use mus_geom_module,                  only: mus_geom_type
  use mus_scheme_type_module,           only: mus_scheme_type
  use mus_particle_aux_module,          only: findPartitionOfTreeID,      &
    &                                         coordExistsOnMyRank,        &
    &                                         coordLocalOnMyRank,         &
    &                                         findPropIndex,              &
    &                                         getCartesianStencilIndices, &
    &                                         getBaryOfCoord
  use mus_particle_logging_type_module, only: mus_particle_logging_type, &
    &                                         pgDebugLog
  use mus_particle_boundary_module,     only: pgBndData, &
    &                                         getNeighborCoord

  implicit none

  integer, public :: mus_pIDvector_type
  integer, public :: mus_wallPos_type
  integer, public :: mus_particleInfo_type
  integer, public :: mus_positionUpdate_type
  integer, public :: mus_particleState_type

  ! ************************************************************************ !
  !> process-wise buffer for particle force contribution data
  !! Used to exchange particle force or velocity data with a specific process
  type mus_particles_pIDvector_type
    sequence
    !> Actual data vector. Used for:
    !! * force and moment contributions: V = (Fx, Fy, Fz, RFx, RFy, RFz)
    !! * velocity updates: V = (vx, vy, vz, RVx, RVy, RVz)
    real(kind=rk) :: V(6)
    !> Force contribution pertains to the particle with this ID
    integer :: pID
  end type mus_particles_pIDvector_type

  type mus_particles_vectorbuffer_type

    !> explicit buffer for data to be transferred
    type(mus_particles_pIDvector_type), allocatable :: val(:)

    !> number of particles to exchange, depends on number of particles at
    !! the interface of two processes.
    !! nParticles can change at every time step!
    integer :: nParticles

    !! Buffer can hold data for at most maxParticles
    !! This one we keep fixed through the duration of the program to
    !! avoid de-allocating, re-allocating memory at every time step.
    integer :: maxParticles

  end type mus_particles_vectorbuffer_type
  ! ************************************************************************ !
  ! ************************************************************************ !
  !> process-wise buffer for wall position data
  !! Used to determine the average wall location using wall position sums
  !! over all processes
  type mus_particles_wallPos_type
    sequence
    !> Actual data vector. Used for:
    !! * wall position sums ( sum of x-coordinates, sum of y, sum of z )
    real(kind=rk) :: V(3)
    !> Wall position sum pertains to the particle with this ID
    !! * I(1) = nWallPos = number of elements in wallPosSum
    !! * I(2) = particleID to which this wallPosSum pertains
    integer :: I(2)
  end type mus_particles_wallPos_type

  type mus_particles_wallbuffer_type

    !> explicit buffer for data to be transferred
    type(mus_particles_wallPos_type), allocatable :: val(:)

    !> number of particles to exchange, depends on number of particles at
    !! the interface of two processes.
    !! nParticles can change at every time step!
    integer :: nParticles

    !! Buffer can hold data for at most maxParticles
    !! This one we keep fixed through the duration of the program to
    !! avoid de-allocating, re-allocating memory at every time step.
    integer :: maxParticles

  end type mus_particles_wallbuffer_type
  ! ************************************************************************ !
  ! ************************************************************************ !
  !> process-wise buffer for particle IDs
  !! Used to inform other processes when a particle needs to be removed/killed

  type mus_particles_IDbuffer_type

    !> explicit buffer for data to be transferred
    integer, allocatable :: val(:)

    !> number of particles to exchange, depends on number of particles at
    !! the interface of two processes.
    !! nParticles can change at every time step!
    integer :: nParticles

    !! Buffer can hold data for at most maxParticles
    !! This one we keep fixed through the duration of the program to
    !! avoid de-allocating, re-allocating memory at every time step.
    integer :: maxParticles

  end type mus_particles_IDbuffer_type
  ! ************************************************************************ !

  ! ************************************************************************ !
  type mus_particles_state_type
    sequence
    !> Continuous particle data to be sent or received
    !! We gather all reals to be sent in one long vector.
    !! The order of the elements is:
    !! pos(1:6) = V(1:6)
    !! vel(1:6) = V(7:12)
    real(kind=rk) :: V(12)

    !> Also pack all integers to be sent or received in one vector
    !! coordOfOrigin = I(1:4)
    !! particleID = I(5)
    integer :: I(5)

  end type mus_particles_state_type

  type mus_particles_statebuffer_type
    type(mus_particles_state_type), allocatable :: val(:)

    !> number of particles to exchange, depends on number of particles at
    !! the interface of two processes.
    !! nParticles can change at every time step!
    integer :: nParticles

    !! Buffer can hold data for at most maxParticles
    !! This one we keep fixed through the duration of the program to
    !! avoid de-allocating, re-allocating memory at every time step.
    integer :: maxParticles

  end type mus_particles_statebuffer_type

  ! ************************************************************************ !
  ! ************************************************************************ !
  type mus_particles_positionUpdate_type
    sequence
    !> Continuous particle data to be sent or received
    !! We gather all reals to be sent in one long vector.
    !! The order of the elements is:
    !! pos(1:6) = V(1:6)
    real(kind=rk) :: V(6)

    !> Also pack all integers to be sent or received in one vector
    !! coordOfOrigin = I(1:4)
    !! particleID = I(5)
    integer :: I(5)

  end type mus_particles_positionUpdate_type

  type mus_particles_posbuffer_type
    type(mus_particles_positionUpdate_type), allocatable :: val(:)

    !> number of particles to exchange, depends on number of particles at
    !! the interface of two processes.
    !! nParticles can change at every time step!
    integer :: nParticles

    !! Buffer can hold data for at most maxParticles
    !! This one we keep fixed through the duration of the program to
    !! avoid de-allocating, re-allocating memory at every time step.
    integer :: maxParticles

  end type mus_particles_posbuffer_type

  ! ************************************************************************ !


  type mus_particles_info_type
    sequence
    !> Continuous particle data to be sent or received
    !! We gather all reals to be sent in one long vector.
    !! The order of the elements is:
    !! pos(1:6) = V(1:6)
    !! vel(1:6) = V(7:12)
    !! Fbuff(1,1:6) = V(13:18)
    !! Fbuff(2,1:6) = V(19:24)
    !! radius = V(25)
    !! mass = V(26)
    !! Fext(1:6) = V(27:32)
    real(kind=rk) :: V(32)

    !> ID of the particle that the data in V pertains to.
    !integer :: pID

    !> Also pack all integers to be sent or received in one vector
    !! coordOfOrigin = I(1:4)
    !! Fnow = I(5)
    !! particleID = I(6)
    integer :: I(6)
  end type mus_particles_info_type

  type mus_particles_infobuffer_type
    type(mus_particles_info_type), allocatable :: val(:)

    !> number of particles to exchange, depends on number of particles at
    !! the interface of two processes.
    !! nParticles can change at every time step!
    integer :: nParticles

    !! Buffer can hold data for at most maxParticles
    !! This one we keep fixed through the duration of the program to
    !! avoid de-allocating, re-allocating memory at every time step.
    integer :: maxParticles
  end type mus_particles_infobuffer_type

  ! ************************************************************************ !

  !> Description of communication data for particles
  !! Will be a member of the particleGroup data type that exists on each process
  !! particleGroup%send for sendbuffer
  !! particleGroup%recv for recvbuffer
  type mus_particles_communication_type

    integer :: nProcs=0     !< amount of partitions to send to

    !> partition MPI rank
    integer,allocatable :: proc(:)

    !> Request handle array
    integer,allocatable :: rqHandle(:)

    !> 2nd Request handle array so we can receive multiple KINDS of messages
    !! simultaneously
    integer,allocatable :: rqHandle2(:)

    !> communication buffer for force contributions
    type( mus_particles_vectorbuffer_type ), allocatable :: buf_force(:)

    !> communication buffer for position updates (inc. coordOfOrigin)
    type( mus_particles_posbuffer_type ), allocatable :: buf_pos(:)

    !> communication buffer for velocity updates
    type( mus_particles_vectorbuffer_type ), allocatable :: buf_vec(:)

    !> communication buffer for wall positions
    type( mus_particles_wallbuffer_type ), allocatable :: buf_wall(:)

    !> communication buffer for particle state (= position + velocity + origin)
    type( mus_particles_statebuffer_type ), allocatable :: buf_state(:)

    !> communication buffer for IDs of particles to be removed/killed
    type( mus_particles_IDbuffer_type ), allocatable :: buf_kill(:)

    !> communication buffer for all particle data together
    type( mus_particles_infobuffer_type ), allocatable :: buf_particle(:)
  end type mus_particles_communication_type

! ************************************************************************ !

contains
! ************************************************************************ !


  !> Allocate mus_particles_communication_type and its variables
  subroutine init_mus_particles_comm_type( me, nProcs, proc )
    ! -------------------------------------------------------------------- !
    !> Communication type to initialize
    type( mus_particles_communication_type ), intent(inout) :: me
    !> Number of neighboring processes with which to exchange particle data
    integer, intent(in) :: nProcs
    !> list of processes (ranks) to exchange with
    integer, intent(in) :: proc(nProcs)
    ! -------------------------------------------------------------------- !

    if ( allocated(me%proc) )         deallocate( me%proc )
    if ( allocated(me%rqHandle) )     deallocate( me%rqHandle )
    if ( allocated(me%rqHandle2) )    deallocate( me%rqHandle2 )
    if ( allocated(me%buf_force) )    deallocate( me%buf_force )
    if ( allocated(me%buf_pos) )      deallocate( me%buf_pos )
    if ( allocated(me%buf_vec) )      deallocate( me%buf_vec )
    if ( allocated(me%buf_wall) )     deallocate( me%buf_wall )
    if ( allocated(me%buf_state) )    deallocate( me%buf_state )
    if ( allocated(me%buf_kill) )     deallocate( me%buf_kill )
    if ( allocated(me%buf_particle) ) deallocate( me%buf_particle )

    me%nProcs = nProcs
    allocate( me%proc( nProcs ) )
    allocate( me%rqHandle( nProcs ) )
    allocate( me%rqHandle2( nProcs ) )
    allocate( me%buf_force( nProcs ) )
    allocate( me%buf_pos( nProcs ) )
    allocate( me%buf_vec( nProcs ) )
    allocate( me%buf_wall( nProcs ) )
    allocate( me%buf_state( nProcs ) )
    allocate( me%buf_kill( nProcs ) )
    allocate( me%buf_particle( nProcs ) )

    ! Set the processes to exchange with
    me%proc(1:nProcs) = proc(1:nProcs)


  end subroutine init_mus_particles_comm_type

  ! ***************************************************************************** !

  subroutine mus_particles_comm_init_vectorbuffer( me, maxParticles )
    ! -------------------------------------------------------------------- !
    type(mus_particles_vectorbuffer_type), intent(inout) :: me
    !> Number of particles to allocate space in buffer for
    !! This should be larger than the number of particles we expect
    !! to be at the interface of two processes at any time!
    integer, intent(in) :: maxParticles
    ! -------------------------------------------------------------------- !

    me%nParticles = 0
    me%maxParticles = maxParticles

    if ( allocated(me%val) ) deallocate(me%val)
    allocate( me%val(maxParticles) )

  end subroutine mus_particles_comm_init_vectorbuffer

  subroutine mus_particles_comm_init_posbuffer( me, maxParticles )
    ! -------------------------------------------------------------------- !
    type(mus_particles_posbuffer_type), intent(inout) :: me
    !> Number of particles to allocate space in buffer for
    !! This should be larger than the number of particles we expect
    !! to be at the interface of two processes at any time!
    integer, intent(in) :: maxParticles
    ! -------------------------------------------------------------------- !

    me%nParticles = 0
    me%maxParticles = maxParticles

    if ( allocated(me%val) ) deallocate(me%val)
    allocate( me%val(maxParticles) )

  end subroutine mus_particles_comm_init_posbuffer

  subroutine mus_particles_comm_init_wallbuffer( me, maxParticles )
    ! -------------------------------------------------------------------- !
    type(mus_particles_wallbuffer_type), intent(inout) :: me
    !> Number of particles to allocate space in buffer for
    !! This should be larger than the number of particles we expect
    !! to be at the interface of two processes at any time!
    integer, intent(in) :: maxParticles
    ! -------------------------------------------------------------------- !

    me%nParticles = 0
    me%maxParticles = maxParticles

    if ( allocated(me%val) ) deallocate(me%val)
    allocate( me%val(maxParticles) )

  end subroutine mus_particles_comm_init_wallbuffer

  subroutine mus_particles_comm_init_statebuffer( me, maxParticles )
    ! -------------------------------------------------------------------- !
    type(mus_particles_statebuffer_type), intent(inout) :: me
    !> Number of particles to allocate space in buffer for.
    !! This should be larger than the number of particles we expect
    !! to be at the interface of two processes at any time!
    integer, intent(in) :: maxParticles
    ! -------------------------------------------------------------------- !

    me%nParticles = 0
    me%maxParticles = maxParticles

    if ( allocated(me%val) ) deallocate(me%val)
    allocate( me%val(maxParticles) )

  end subroutine mus_particles_comm_init_statebuffer

  subroutine mus_particles_comm_init_IDbuffer( me, maxParticles )
    ! -------------------------------------------------------------------- !
    type(mus_particles_IDbuffer_type), intent(inout) :: me
    !> Number of particles to allocate space in buffer for
    !! This should be larger than the number of particles we expect
    !! to be at the interface of two processes at any time!
    integer, intent(in) :: maxParticles
    ! -------------------------------------------------------------------- !

    me%nParticles = 0
    me%maxParticles = maxParticles

    if ( allocated(me%val) ) deallocate(me%val)
    allocate( me%val(maxParticles) )

  end subroutine mus_particles_comm_init_IDbuffer

  subroutine mus_particles_comm_init_particlebuffer( me, maxParticles )
    ! -------------------------------------------------------------------- !
    type(mus_particles_infobuffer_type), intent(inout) :: me
    !> Number of particles to allocate space in buffer for
    !! This should be larger than the number of particles we expect
    !! to be at the interface of two processes at any time!
    integer, intent(in) :: maxParticles
    ! -------------------------------------------------------------------- !

    me%nParticles = 0
    me%maxParticles = maxParticles

    if ( allocated(me%val) ) deallocate(me%val)
    allocate( me%val(maxParticles) )

  end subroutine mus_particles_comm_init_particlebuffer

  ! *************************** DERIVED MPI TYPES *************************** !
  subroutine mus_particles_initForceContributionMPItype( )
    ! -------------------------------------------------------------------- !

    type(mus_particles_pIDvector_type) :: sampleForceContribution

    integer, parameter :: Nblocks = 2
    integer :: blocklengths(Nblocks), types(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: adresses(Nblocks), displacements(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: base_adress

    integer :: ierror
    ! -------------------------------------------------------------------- !

    blocklengths = (/ 6, 1 /)
    types = (/ RK_MPI, MPI_INTEGER /)
    ! Get the adresses of the sample object so we can fill the displacements array
    call MPI_get_address(sampleForceContribution, base_adress, ierror)
    call MPI_get_address(sampleForceContribution%V, adresses(1), ierror)
    call MPI_get_address(sampleForceContribution%pID, adresses(2), ierror)

    ! displacement is given by difference of subsequent adresses
    displacements(1) = MPI_Aint_diff(adresses(1), base_adress)
    displacements(2) = MPI_Aint_diff(adresses(2), base_adress)

    ! Create and commit the force contribution MPI type
    call MPI_Type_create_struct(Nblocks, blocklengths, displacements, &
                               & types, mus_pIDvector_type, ierror    )
    call MPI_Type_commit( mus_pIDvector_type, ierror)

  end subroutine mus_particles_initForceContributionMPItype


  subroutine mus_particles_initWallPosMPItype( )
    ! -------------------------------------------------------------------- !

    type(mus_particles_wallPos_type) :: sampleWallPos

    integer, parameter :: Nblocks = 2
    integer :: blocklengths(Nblocks), types(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: adresses(Nblocks), displacements(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: base_adress

    integer :: ierror
    ! -------------------------------------------------------------------- !

    blocklengths = (/ 3, 2 /)
    types = (/ RK_MPI, MPI_INTEGER /)
    ! Get the adresses of the sample object so we can fill the displacements array
    call MPI_get_address(sampleWallPos, base_adress, ierror)
    call MPI_get_address(sampleWallPos%V, adresses(1), ierror)
    call MPI_get_address(sampleWallPos%I, adresses(2), ierror)

    ! displacement is given by difference of subsequent adresses
    displacements(1) = MPI_Aint_diff(adresses(1), base_adress)
    displacements(2) = MPI_Aint_diff(adresses(2), base_adress)

    ! Create and commit the wallPos MPI type
    call MPI_Type_create_struct(Nblocks, blocklengths, displacements, &
                               & types, mus_wallPos_type, ierror    )
    call MPI_Type_commit( mus_wallPos_type, ierror)

  end subroutine mus_particles_initWallPosMPItype

  !> Routine to create and commit the MPI type used to communicate
  !! position and coordOfOrigin
  subroutine mus_particles_initPositionUpdateMPItype( )
    ! -------------------------------------------------------------------- !

    type(mus_particles_positionUpdate_type) :: sample

    integer, parameter :: Nblocks = 2
    integer :: blocklengths(Nblocks), types(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: adresses(Nblocks), displacements(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: base_adress

    integer :: ierror
    ! -------------------------------------------------------------------- !

    blocklengths = (/ 6, 5 /)
    types = (/ RK_MPI, MPI_INTEGER /)
    ! Get the adresses of the sample object so we can fill the displacements array
    call MPI_get_address(sample, base_adress, ierror)
    call MPI_get_address(sample%V, adresses(1), ierror)
    call MPI_get_address(sample%I, adresses(2), ierror)

    ! displacement is given by difference of subsequent adresses
    displacements(1) = MPI_Aint_diff(adresses(1), base_adress)
    displacements(2) = MPI_Aint_diff(adresses(2), base_adress)

    ! Create and commit the force contribution MPI type
    call MPI_Type_create_struct(Nblocks, blocklengths, displacements, &
                               & types, mus_positionUpdate_type, ierror )
    call MPI_Type_commit( mus_positionUpdate_type, ierror)

  end subroutine mus_particles_initPositionUpdateMPItype

  !> Routine to create and commit the MPI type used to communicate
  !! position, velocity and coordOfOrigin
  subroutine mus_particles_initParticleStateMPItype( )
    ! -------------------------------------------------------------------- !

    type(mus_particles_state_type) :: sample

    integer, parameter :: Nblocks = 2
    integer :: blocklengths(Nblocks), types(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: adresses(Nblocks), displacements(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: base_adress

    integer :: ierror
    ! -------------------------------------------------------------------- !

    blocklengths = (/ 12, 5 /)
    types = (/ RK_MPI, MPI_INTEGER /)
    ! Get the adresses of the sample object so we can fill the displacements array
    call MPI_get_address(sample, base_adress, ierror)
    call MPI_get_address(sample%V, adresses(1), ierror)
    call MPI_get_address(sample%I, adresses(2), ierror)

    ! displacement is given by difference of subsequent adresses
    displacements(1) = MPI_Aint_diff(adresses(1), base_adress)
    displacements(2) = MPI_Aint_diff(adresses(2), base_adress)

    ! Create and commit the force contribution MPI type
    call MPI_Type_create_struct(Nblocks, blocklengths, displacements, &
                               & types, mus_particleState_type, ierror )
    call MPI_Type_commit( mus_particleState_type, ierror)

  end subroutine mus_particles_initParticleStateMPItype



  !> Routine to create and commit the MPI type used to communicate
  !! all continuous particle data
  subroutine mus_particles_initParticleInfoMPItype( )
    ! -------------------------------------------------------------------- !

    type(mus_particles_info_type) :: sample

    integer, parameter :: Nblocks = 2
    integer :: blocklengths(Nblocks), types(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: adresses(Nblocks), displacements(Nblocks)
    integer(kind=MPI_ADDRESS_KIND) :: base_adress

    integer :: ierror
    ! -------------------------------------------------------------------- !

    blocklengths = (/ 32, 6 /)
    types = (/ RK_MPI, MPI_INTEGER /)
    ! Get the adresses of the sample object so we can fill the displacements array
    call MPI_get_address(sample, base_adress, ierror)
    call MPI_get_address(sample%V, adresses(1), ierror)
    call MPI_get_address(sample%I, adresses(2), ierror)

    ! displacement is given by difference of subsequent adresses
    displacements(1) = MPI_Aint_diff(adresses(1), base_adress)
    displacements(2) = MPI_Aint_diff(adresses(2), base_adress)

    ! Create and commit the force contribution MPI type
    call MPI_Type_create_struct(Nblocks, blocklengths, displacements, &
                               & types, mus_particleInfo_type, ierror )
    call MPI_Type_commit( mus_particleInfo_type, ierror)

  end subroutine mus_particles_initParticleInfoMPItype



  ! *************** DEBUGGING  *********************************************** !
  subroutine print_particles_comm( me )
    !> Communication type to print
    type( mus_particles_communication_type ), intent(inout) :: me

    ! -------------------------------------------------------------------- !
    character(len=1024) :: format_string
    integer :: iproc
    ! -------------------------------------------------------------------- !
    if( me%Nprocs > 0 ) then
      write(format_string, '(A,I0,A)') '(', me%Nprocs, 'I5)'
      format_string = trim(format_string)
    else
      return
    end if

    write(logUnit(1),*) 'particles_communication_type procs:'

    write(logUnit(1),'(A)', advance='no') '[ '
    write(logUnit(1),format_string, advance='no') &
      &                           ( me%proc(iproc), iproc = 1,me%nProcs )
    write(logUnit(1),'(A)') ']'

    do iproc = 1,me%nProcs
      write(logUnit(1),*) '--- Proc ', me%proc(iproc), ' ---'
      write(logUnit(1),*) 'Force: ', me%buf_force(iproc)%nParticles, ' particles'
      call print_particles_pIDvectorbuffer( me%buf_force(iproc) )
      write(logUnit(1),*) 'vec: ', me%buf_vec(iproc)%nParticles, ' particles'
      call print_particles_pIDvectorbuffer( me%buf_vec(iproc) )
      write(logUnit(1),*) 'state: ', me%buf_state(iproc)%nParticles, ' particles'
      call print_particles_statebuffer( me%buf_state(iproc) )
    end do ! iproc
  end subroutine print_particles_comm

  subroutine print_particles_pIDvectorbuffer( buff )
    !> Communication type to print
    type(mus_particles_vectorbuffer_type), intent(in) :: buff

    ! -------------------------------------------!
    integer iParticle, i
    ! -------------------------------------------!
    do iParticle = 1, buff%nParticles
      write(logUnit(1), '(A,I0,A)', advance='no') 'pID = ',  buff%val(iParticle)%pID, ' buff = ['
      write(logUnit(1), '(6E10.3)', advance = 'no') ( buff%val(iParticle)%V(i), &
                                                    & i = 1,6                      )
      write(logUnit(1), '(A)') ' ]'
    end do

  end subroutine print_particles_pIDvectorbuffer

  subroutine print_particles_statebuffer( buff )
    !> Communication type to print
    type(mus_particles_statebuffer_type), intent(in) :: buff

    ! -------------------------------------------!
    integer iParticle, i
    ! -------------------------------------------!
    do iParticle = 1, buff%nParticles
      write(logUnit(1), '(I0,A)', advance='no') buff%val(iParticle)%I(5), ' real buff = ['
      write(logUnit(1), '(6E10.3)', advance = 'no') ( buff%val(iParticle)%V(i), &
                                                    & i = 1,12                      )
      write(logUnit(1), '(A)') ' ]'
      write(logUnit(1), '(I0,A)', advance='no') buff%val(iParticle)%I(5), ' int buff = ['
      write(logUnit(1), '(6I0.5)', advance = 'no') ( buff%val(iParticle)%I(i), &
                                                    & i = 1,4                      )
      write(logUnit(1), '(A)') ' ]'
    end do

  end subroutine print_particles_statebuffer

  !> find_particle_comm_procs generates a list of all processes we may need to communicate
  !! particle data with over the duration of the simulation. We only need to communicate
  !! this data for particles on our rank which also exist on other processes. This means
  !! we will never communicate with processes whose domains are not within one particle
  !! diameter of our domain. This is the key assumption used to generate the list of procs
  !! Hence the padding distance dpad should be set to the largest particle diameter in the
  !! simulation
  subroutine find_particle_comm_procs4( prunedProcs, scheme, geometry, myRank, dpad )
    !> Dynamic array of procs which we will store the found procs in
    type(dyn_intarray_type) :: prunedProcs
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry information to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> This process's rank
    integer, intent(in) :: myRank
    !> Padding distance
    real(kind=rk) :: dpad
    ! -------------------------------------------!
    type(dyn_intarray_type) :: procs ! array of procs before pruning
    integer :: lev, Ndpad, upperBound
    integer :: searchBoxLims(6), dirs(6)
    integer :: iproc, foundProc, newPos
    integer :: nx, ny, nz
    integer :: iElem, nFluids
    integer :: elemCoord(4), coord(4)
    integer(kind=long_k) :: TIDoffset, TreeID
    logical :: wasAdded
    real(kind=rk) :: dx

    ! For debugging only
    integer :: boundingBox(6)
    real(kind=rk) :: coordX(3)
    ! -------------------------------------------!
    lev = geometry%tree%global%maxLevel
    upperBound = 2**lev
    TIDoffset = tem_firstIdAtLevel(lev)
    dx = geometry%tree%global%BoundingCubeLength / 2**lev
    Ndpad = ceiling(dpad / dx)
    nFluids = scheme%pdf(lev)%nElems_fluid

    ! Get the indices of the cartesian directions in stencil
    call getCartesianStencilIndices( scheme%layout%fStencil%cxDir, dirs )

    ! --- FOR DEBUGGING: look at the bounding box of elements we searched --- !
    ! Initialize boundingBox coordinates using the coordinates of the first halo

    elemCoord = tem_coordOfID( TreeID = scheme%levelDesc(lev)%total(1), &
                             & offset = TIDoffset                             )

    boundingBox(1:2) = elemCoord(1)
    boundingBox(3:4) = elemCoord(2)
    boundingBox(5:6) = elemCoord(3)

    ! --- /FOR DEBUGGING --- !

    ! Loop over all local fluid elements in this process
    do iElem = 1, nFluids
      ! Check if elem has prp_sendhalo
      if( btest( scheme%levelDesc(lev)%property(iElem), prp_sendHalo) ) then
        ! For elems with prp_sendhalo, we will look to the surrounding
        ! elements to find which remote processes they belong to

        ! Get the searchBoxLimits
        elemCoord = tem_coordOfID( TreeID = scheme%levelDesc(lev)%total(iElem), &
                                 & offset = TIDoffset                           )

        call getSearchBoxLimits( elemCoord = elemCoord,    &
                               & elemPos   = iElem,        &
                               & scheme    = scheme,       &
                               & geometry  = geometry,     &
                               & dirs      = dirs,         &
                               & myRank    = myRank,       &
                               & ndpad     = Ndpad,        &
                               & lims      = searchBoxLims )

        ! Now loop over these directions and identify to which neighbor
        ! procs the elements belong
        do nx = searchBoxLims(1), searchBoxLims(2)
          do ny = searchBoxLims(3), searchBoxLims(4)
            do nz = searchBoxLims(5), searchBoxLims(6)
              coord(:) = getNeighborCoord( coord        = elemCoord, &
                                         & nx           = nx,        &
                                         & ny           = ny,        &
                                         & nz           = nz,        &
                                         & boundaryData = pgBndData  )

              ! Check that the coordinates are within bounds
              if( any(coord(1:3) < 0) .OR. any(coord(1:3) > upperBound) ) then
                cycle
              else ! coordinate is within domain
                ! --- DEBUGGING --- !
                if( coord(1) < boundingBox(1) ) then
                  boundingBox(1) = coord(1)
                end if
                if( coord(2) < boundingBox(3) ) then
                  boundingBox(3) = coord(2)
                end if
                if( coord(3) < boundingBox(5) ) then
                  boundingBox(5) = coord(3)
                end if
                if( coord(1) > boundingBox(2) ) then
                  boundingBox(2) = coord(1)
                end if
                if( coord(2) > boundingBox(4) ) then
                  boundingBox(4) = coord(2)
                end if
                if( coord(3) > boundingBox(6) ) then
                  boundingBox(6) = coord(3)
                end if
                ! --- /DEBUGGING --- !

                TreeID = tem_IdOfCoord(coord = coord)

                ! Determine what process this TreeID is on
                call findPartitionOfTreeID( &
                                      & sTreeID = TreeID,                     &
                                      & geometry = geometry,                  &
                                      & nProcs = geometry%tree%global%nParts, &
                                      & outProc = foundProc,                  &
                                      & procIndex = iproc                     )

                ! Add this process to the dynamic array
                if( foundProc >= 0 .AND. foundProc /= myRank ) then
                  ! Note: append automatically only adds when element is not
                  ! already in array

                  call append( me       = procs,  &
                    &          val      = foundProc,  &
                    &          pos      = newPos,     &
                    &          wasAdded = wasAdded     )
                  if(wasAdded) then
                    ! For debugging
                    coordX = getBaryOfCoord( coord   = coord(1:3),                  &
                                           & origin  = geometry%tree%global%origin, &
                                           & dx      = dx                           )

                    open(pgDebugLog%lu, file = pgDebugLog%lfile, status = 'old', position = 'append')
                    write(pgDebugLog%lu,*) "proc", foundProc, " treeID ", TreeID
                    write(pgDebugLog%lu,'(A,3E17.9)' ) "remote bary ", coordX(1:3)

                    coordX = getBaryOfCoord( coord  = elemCoord(1:3),              &
                                           & origin = geometry%tree%global%origin, &
                                           & dx     = dx                           )

                    write(pgDebugLog%lu,*) "Elem treeID ", scheme%levelDesc(lev)%total(iElem)
                    write(pgDebugLog%lu,'(A,3E17.9)' ) "Elem bary   ", coordX(1:3)
                    write(pgDebugLog%lu,'(A,6I4)') "SearchBoxLims = ", searchBoxLims(1:6)
                    close(pgDebugLog%lu)
                  end if ! wasAdded

                end if ! foundProc > 0 and not myRank
              end if ! coord within bounds
            end do ! nz
          end do ! ny
        end do ! nx
      end if ! elem has prp_sendHalo

    end do ! iElem

    ! --- DEBUGGING --- !
    ! Convert boundingBox integer coordinates to cartesian and write to log
    coord(1) = boundingBox(1)
    coord(2) = boundingBox(3)
    coord(3) = boundingBox(5)
    coordX = getBaryOfCoord( coord = coord(1:3),                    &
                           & origin  = geometry%tree%global%origin, &
                           & dx      = dx                           )

    open(pgDebugLog%lu, file = pgDebugLog%lfile, status = 'old', position = 'append')
    write(pgDebugLog%lu,*) "find_particle_comm_procs4 boundingBox"
    write(pgDebugLog%lu,'(A,3E17.9)') "min = ", coordX(1:3)

    coord(1) = boundingBox(2)
    coord(2) = boundingBox(4)
    coord(3) = boundingBox(6)
    coordX = getBaryOfCoord( coord = coord(1:3),                    &
                           & origin  = geometry%tree%global%origin, &
                           & dx      = dx                           )
    write(pgDebugLog%lu,'(A,3E17.9)' ) "max = ", coordX(1:3)
    close(pgDebugLog%lu)
    ! --- /DEBUGGING --- !

    ! Now that we have a set of procs using local information,
    ! Prune this list by checking if the procs we talk to will also
    ! talk back.

    if( procs%nvals > 0 ) then
            call pruneParticleCommProcs( procs       = procs, &
                                       & prunedProcs = prunedProcs, &
                                       & nParts      = geometry%tree%global%nParts, &
                                       & myRank      = myRank)
    end if





  end subroutine find_particle_comm_procs4


  !> getSearchBoxLimits is used in find_particle_comm_procs to determine
  !! where the non-local elements we have to check the owner proc of are.
  subroutine getSearchBoxLimits( elemCoord, elemPos, scheme, geometry, &
                               & dirs, myRank, ndpad, lims             )
    !> Coordinate to get search box around
    integer, intent(in) :: elemCoord(4)
    !> Element position in total list
    integer, intent(in) :: elemPos
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry information to determine whether halo has boundaries
    type(mus_geom_type), intent(in) :: geometry
    !> indices of cartesian directions [-x +x -y +y -z +z] in stencil
    integer, intent(in) :: dirs(6)
    !> This procs rank
    integer, intent(in) :: myRank
    !> Length of the search box in the directions where there
    !! are no local elements
    integer, intent(in) :: ndpad
    !> lims = [xmin, xmax, ymin, ymax, zmin, zmax]
    integer, intent(out) :: lims(6)
    ! -------------------------------------------!
    integer :: coord(4)
    integer :: lev
    integer :: iDir, posInBnd
    integer(kind=long_k) :: bcID
    ! -------------------------------------------!
    lev = elemCoord(4)

    lims = (/ -ndpad, ndpad, -ndpad, ndpad, -ndpad, ndpad  /)

    ! Searchbox should not extend outside domain, so check if halo has
    ! any adjacent boundaries
    ! If so, check in which of the Cartesian dirs they are.

    ! Get position of element in boundary ID array
    if( allocated(geometry%posInBndID ) ) then
      posInBnd = geometry%posInBndID(elemPos)
      if( posInBnd > 0 ) then
        do iDir = 1, 6
          if( dirs(iDir) < 0)  then
            cycle
          end if

          bcID = geometry%boundary%boundary_ID( dirs(iDir), posInBnd )
          if(bcID /= 0) then
            ! If there is a boundary in this direction, do not search
            ! for remote processes in this direction
            lims(iDir) = 0
          end if
        end do
      end if ! posInBnd > 0
    end if ! allocated(posInBndId)

    ! We also do not have to search in directions pointing INTO
    ! the local domain of this process. So look for local elements
    ! in the stencil directions and set lims accordingly.
    ! xmin
    coord(:) = elemCoord(:) + (/ -1, 0, 0, 0/)
    if( coordLocalOnMyRank( coord, scheme, myRank ) ) then
      lims(1) = 0
    else
      lims(1) = -ndpad
    end if

    ! xmax
    coord(:) = elemCoord(:) + (/ 1, 0, 0, 0/)
    if( coordLocalOnMyRank( coord, scheme, myRank ) ) then
      lims(2) = 0
    else
      lims(2) = ndpad
    end if

    ! ymin
    coord(:) = elemCoord(:) + (/ 0, -1, 0, 0/)
    if( coordLocalOnMyRank( coord, scheme, myRank ) ) then
      lims(3) = 0
    else
      lims(3) = -ndpad
    end if

    ! ymax
    coord(:) = elemCoord(:) + (/ 0, 1, 0, 0/)
    if( coordLocalOnMyRank( coord, scheme, myRank ) ) then
      lims(4) = 0
    else
      lims(4) = ndpad
    end if

    ! zmin
    coord(:) = elemCoord(:) + (/ 0, 0, -1, 0/)
    if( coordLocalOnMyRank( coord, scheme, myRank ) ) then
      lims(5) = 0
    else
      lims(5) = -ndpad
    end if

    ! zmax
    coord(:) = elemCoord(:) + (/ 0, 0, 1, 0/)
    if( coordLocalOnMyRank( coord, scheme, myRank ) ) then
      lims(6) = 0
    else
      lims(6) = ndpad
    end if

  end subroutine getSearchBoxLimits

  !> pruneParticleCommProcs takes the initial list of processes we think
  !! we need to communicate particle data with (determined using only
  !! local data) and checks whether those processes also think they
  !! need to communicate with us. The prunedProcs array contains the
  !! final list of processes for which we know we will send and also
  !! receive messages.
  subroutine pruneParticleCommProcs(procs, prunedProcs, nParts, myRank)
    !> Initial list of particle comm procs to be pruned
    type(dyn_intarray_type) :: procs
    !> Final (pruned) list of procs
    type(dyn_intarray_type) :: prunedProcs
    !> number of processes
    integer, intent(in) :: nParts
    !> This procs rank
    integer, intent(in) :: myRank
    ! ------------------------------------- !
    integer :: iproc, ierr
    integer :: irank
    integer :: send_rqhandle(nParts)
    integer :: recv_rqhandle(nParts)
    integer :: recv_status( mpi_status_size, nParts)
    integer :: send_status( mpi_status_size, nParts)
    logical :: sendbuff(nParts)
    logical :: recvbuff(nParts)

    logical :: remoteKnowsMe
    integer :: newPos
    logical :: wasAdded
    ! ------------------------------------- !
    ! Construct the message for each process
    do iproc = 1, nParts
      irank = iproc - 1
      if ( irank == myRank .OR. any(irank == procs%val(1:procs%nvals) ) ) then
        sendbuff(iproc) = .TRUE.
      else
        sendbuff(iproc) = .FALSE.
      end if
    end do

    ! Start receive communications
    do iproc = 1, nParts
      irank = iproc - 1
      call mpi_irecv(                 &
       &      recvbuff(iproc),        & ! me
       &      1,                      & ! max count (per send proc)
       &      MPI_LOGICAL,            & ! data type
       &      irank,                  & ! source process
       &      irank,                  & ! flag
       &      MPI_COMM_WORLD,         & ! communicator
       &      recv_rqhandle(iproc),   & ! handle
       &      ierr                    ) ! error status

    end do

    ! Send the message for each proc
    do iproc = 1, nParts
      irank = iproc - 1

      call mpi_isend(                 &
       &      sendBuff(iproc),        & ! buffer
       &      1,                      & ! count
       &      MPI_LOGICAl,            & ! data type
       &      irank,                  & ! target
       &      myRank,                 & ! integer tag
       &      MPI_COMM_WORLD,         & ! communicator
       &      send_rqhandle( iproc ), & ! handle
       &      ierr                    ) ! error status

    end do

    ! wait for receive buffer to be ready
    call mpi_waitall(nParts,          & ! count
      &              recv_rqhandle,   & ! request handles
      &              recv_status,     & ! mpi status
      &              ierr             ) ! error status

    ! For each message, check if the content matches our own
    ! find_particle_comm_procs
    do iproc = 1, procs%nvals
      remoteKnowsMe = recvBuff( procs%val(iproc) + 1 )

      if( remoteKnowsMe ) then
        call append( me     = prunedProcs,       &
        &          val      = procs%val(iproc),  &
        &          pos      = newPos,            &
        &          wasAdded = wasAdded           )
      end if

    end do

    ! Wait for send buffer to be ready
    call mpi_waitall(nParts,          & ! count
      &              send_rqhandle,   & ! request handles
      &              send_status,     & ! mpi status
      &              ierr             ) ! error status

    open(pgDebugLog%lu, file = pgDebugLog%lfile, status = 'old', position = 'append')
    write(pgDebugLog%lu,*) "Original particle comm procs array"
    write(pgDebugLog%lu,*) procs%val(1:procs%nvals)
    write(pgDebugLog%lu,*) "Pruned particle comm procs array"
    write(pgDebugLog%lu,*) prunedProcs%val(1:prunedProcs%nvals)
    close(pgDebugLog%lu)

  end subroutine pruneParticleCommProcs

  !> getSearchBoxLimits is used in find_particle_comm_procs to determine
  !! where the non-local elements we have to check the owner proc of are.
  ! subroutine getSearchBoxLimits( haloCoord, scheme, myRank, ndpad, lims )
  !   !> Coordinate to get search box around
  !   integer, intent(in) :: haloCoord(4)
  !   !> Scheme for access to level descriptor
  !   type(mus_scheme_type), intent(in) :: scheme
  !   !> This procs rank
  !   integer, intent(in) :: myRank
  !   !> Length of the search box in the directions where there
  !   !! are no local elements
  !   integer, intent(in) :: ndpad
  !   !> lims = [xmin, xmax, ymin, ymax, zmin, zmax]
  !   integer, intent(out) :: lims(6)
  !   ! -------------------------------------------!
  !   integer :: coord(4)
  !   ! -------------------------------------------!
  !
  !   lims = (/ -ndpad, ndpad, -ndpad, ndpad, -ndpad, ndpad  /)
  !
  !   ! xmin
  !   coord(:) = haloCoord(:) + (/ -1, 0, 0, 0/)
  !   if( coordExistsOnMyRank( coord, scheme, myRank ) ) then
  !     lims(1) = 0
  !   else
  !     lims(1) = -ndpad
  !   end if
  !
  !   ! xmax
  !   coord(:) = haloCoord(:) + (/ 1, 0, 0, 0/)
  !   if( coordExistsOnMyRank( coord, scheme, myRank ) ) then
  !     lims(2) = 0
  !   else
  !     lims(2) = ndpad
  !   end if
  !
  !   ! ymin
  !   coord(:) = haloCoord(:) + (/ 0, -1, 0, 0/)
  !   if( coordExistsOnMyRank( coord, scheme, myRank ) ) then
  !     lims(3) = 0
  !   else
  !     lims(3) = -ndpad
  !   end if
  !
  !   ! ymax
  !   coord(:) = haloCoord(:) + (/ 0, 1, 0, 0/)
  !   if( coordExistsOnMyRank( coord, scheme, myRank ) ) then
  !     lims(4) = 0
  !   else
  !     lims(4) = ndpad
  !   end if
  !
  !   ! zmin
  !   coord(:) = haloCoord(:) + (/ 0, 0, -1, 0/)
  !   if( coordExistsOnMyRank( coord, scheme, myRank ) ) then
  !     lims(5) = 0
  !   else
  !     lims(5) = -ndpad
  !   end if
  !
  !   ! zmax
  !   coord(:) = haloCoord(:) + (/ 0, 0, 1, 0/)
  !   if( coordExistsOnMyRank( coord, scheme, myRank ) ) then
  !     lims(6) = 0
  !   else
  !     lims(6) = ndpad
  !   end if
  !
  !
  ! end subroutine getSearchBoxLimits

end module mus_particle_comm_type_module
