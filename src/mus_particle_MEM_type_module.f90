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













! --- Macros for loading particle data from lua script --- !



!> Provides the data type for the moment exchange particles
module mus_particle_MEM_type_module
  use env_module, only: rk, zeroLength, minLength

  use tem_dyn_array_module, only: dyn_intArray_type
  use tem_grow_array_module, only: grw_intarray_type
  use tem_logging_module, only: logUnit

  use mus_particle_array_module, only: maxContainerSize

  implicit none

  private

  ! Basic particle type
  type mus_particle_MEM_type
    !- Unique particleID for identification across all processes
    integer :: particleID

    !- Owner process of this particle
    !  Process is owner if coordOfOrigin is local to this process
    integer :: owner

    !> Process who was owner in last time step. We need this for the
    !! averaging of forces over two time steps
    integer :: previousOwner = -1

    !- Array matching the shape of the send%proc array
    !  the ith entry indicates whether this particle exists on
    !  send%proc(i) and hence we need to send velocity updates there
    logical, allocatable :: existsOnProc(:)

    !- Array matching the shape of the send%proc array
    !  the ith entry indicates whether this particle is new on
    !  send%proc(i) in which case we need to send over data
    !  to that proc so it can be added to its particleGroup
    logical, allocatable :: addToProc(:)

    !- Array matching the shape of the send%proc array
    !  the ith entry indicates whether this particle is new on
    !  send%proc(i) in which case we need to send over data
    !  to that proc so it can be added to its particleGroup
    logical, allocatable :: removeFromProc(:)

    !> Logical which tells us whether to initialize this particle or not
    !! is set to true only immediately after receiving this particle
    !! from a neighboring process
    logical :: newForMe = .FALSE.

    !> hasCollided tells us whether particle has just had its velocity modified
    !! in a collision and that this information needs to be sent to other processes
    logical :: hasCollided = .FALSE.

    !> removeParticle indicates that this particle needs to be removed after e.g.
    !! hitting an open boundary. This information is first sent to all other procs
    !! that know about this particle, then the particle is actually removed from
    !! the particleGroup.
    logical :: removeParticle_global = .FALSE.

    !- Radius of (spherical) particle
    real(kind=rk) :: radius

    !- mass and rotational inertia
    real(kind=rk) :: mass
    real(kind=rk) :: rotInertia

    ! radius in lattice units, rounded up
    integer :: Rn

    !- Particle origin position and translational + angular velocity (x,y,z,rx,ry,rz)
    real(kind=rk) :: vel(6)
    real(kind=rk) :: pos(6)
    real(kind=rk) :: oldPos(6)

    ! integer coordinate of desired point
    integer :: coordOfOrigin(4)
    integer :: oldCoordOfOrigin(4)


    ! Hydrodynamic force and torque acting on the particle (fx,fy,fz,mx,my,mz)
    ! F is the average of hydrodynamic force computation in current
    ! and last time step F = 0.5 * ( Fbuff(1,:) + Fbuff(2,:) )
    real(kind=rk) :: F(6) = 0.0_rk
    real(kind=rk) :: Fbuff(2,6) = 0.0_rk

    ! External force vector which is loaded from the lua script.
    ! Can only be constant as of now (for example to simulate gravity)
    real(kind=rk) :: Fext(6) = 0.0_rk

    ! Index pointing to which row of Fbuff to fill at this time step
    integer :: Fnow = 1
    ! Index pointing to which row of Fbuff contains force data of last time step
    integer :: Flast = 2

    ! Buffer for the DEM collision force and torque
    ! We need to store forces at last two times for velocity verlet integration
    real(kind=rk) :: F_DEM(2,6) = 0.0_rk

    ! Index pointing to row in Fcoll for current DEM time step (t)
    integer :: F_DEM_now = 1
    ! Index pointing to row in Fcoll at next DEM time step (t + dt_dem)
    integer :: F_DEM_next = 2

    ! For DEM treatment of particle-wall collisions we store the distance
    ! to the nearest wall and the normal vector (pointing towards the wall)
    ! of that wall
    integer :: nWallPos = 0 ! number of elements in wallPosSum
    real(kind=rk) :: avgWallPos(3) = 0.0_rk
    real(kind=rk) :: rwall(3) = 0.0_rk

    ! Logical that indicates whether this particle is close enough to a wall
    ! that we need to compute wall interactions during the subcycling loop
    logical :: interactWithWall = .FALSE.

    !- Dynamic array containing indices of the currently covered fluid elements
    !- These need to be excluded in the loop over elements in the kernel
    !- Pertains to levelDesc total list (from which kernel lists are generated)
    type(dyn_intArray_type) :: exclusionList

    !> Buffer for exclusion list used in moveParticle routine
    !  Used to determine newly uncovered fluid neighbors
    type(dyn_intArray_type) :: exclusionListBuffer

    !> Number of fluid neighbors for this particle
    integer :: NfluidNeighbors

    !> Indices in levelDesc total list of elements that
    !  need to be turned to fluid after moving particle
    type(grw_intArray_type) :: makeFluidList

  end type mus_particle_MEM_type

! \brief smart dynamic array (da) for
!
! this datatype implements a dynamic array for particles,
! which is capable of growing and adding of unique elements.
! it is available for various types of particles.
! here we deal with $tstring$.
!> dynamic array (da) type for mus_particle_mem_type
type dyn_particle_mem_array_type
  ! current number of particles in the array
  integer :: nvals = 0

  ! current container size
  integer :: containersize = 0

  ! the actual array holding the particles
  type(mus_particle_mem_type), allocatable :: val(:)

  ! array of particle id's which mirrors the actual particle array val(:)
  ! so order of the id's in this list is identical to order of val(:)
  integer, allocatable :: pidlist(:)

  ! sorted positions of the particle id list
  ! example, just like the tem dyn array types
  ! pidlist:  8, 6, 7, 4, 5
  ! pidsort:  4, 5, 2, 3, 1
  ! so we can for example traverse the pidlist in sorted order using
  ! do i = 1, nvals
  !   val = pidlist(pidsort(i))
  ! end do

  integer, allocatable :: pidsort(:)

end type

  interface allocateProcessMasks
    module procedure allocateProcessMasks_MEM
  end interface

  public :: mus_particle_MEM_type
  public :: allocateProcessMasks
  public :: dyn_particle_MEM_array_type
  public :: init_da_particle_MEM
  public :: destroy_da_particle_MEM
  public :: append_da_particle_MEM
  public :: expand_da_particle_MEM
  public :: truncate_da_particle_MEM
  public :: swap_da_particle_MEM
  public :: remove_particle_from_da_particle_MEM
  public :: sortposofval_particle_MEM


contains


! ---- Dynamic particle array methods for momentum exchange method ---- !

  subroutine init_da_particle_mem(me, length)
    type(dyn_particle_mem_array_type), intent(inout) :: me !< dynamic array to init
    integer, intent(in), optional :: length !< initial length of the container

    !----------------------------------------------!
    if (present(length)) then
      me%containersize = length
    else
      me%containersize = zerolength
    end if

    ! deallocate ...
    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%pidlist ) ) deallocate(me%pidlist)
    if( allocated( me%pidsort ) ) deallocate(me%pidsort)

    ! ... and reallocate
    allocate(me%val(me%containersize))
    allocate(me%pidlist(me%containersize))
    allocate(me%pidsort(me%containersize))
    me%nvals = 0

  end subroutine init_da_particle_mem

  ! destroy particle dynamic array
  subroutine destroy_da_particle_mem(me)
    type(dyn_particle_mem_array_type), intent(inout) :: me !< dynamic array to destroy
    !----------------------------------------------!

    me%containersize = 0
    me%nvals         = 0

    if( allocated( me%val ) ) deallocate(me%val)
    if( allocated( me%pidlist ) ) deallocate(me%pidlist)
    if( allocated( me%pidsort ) ) deallocate(me%pidsort)
  end subroutine destroy_da_particle_mem

  subroutine append_da_particle_mem(me, particle, length,  wasadded)
    ! dynamic array to append to
    type(dyn_particle_mem_array_type), intent(inout) :: me
    ! particles to append
    type(mus_particle_mem_type), intent(inout) :: particle
    ! number of elements to increase containersize with
    integer, intent(in), optional :: length
    ! logical will be set to 1 is appending was succesful, 0 if not
    logical :: wasadded !<
    !----------------------------------------------!
    integer :: foundpos
    integer :: i

    ! do a binary search on existing entries (returns closest entry next to
    ! it if not found).
    foundpos = sortposofval_particle_mem(me, particle%particleid, .true.)
    wasadded = .false.

    ! if it found the value, the position is smaller than nvals
    if (foundpos <= me%nvals) then

      ! check if particle id is already in particle array
      ! if so, do nothing
      if ( me%pidlist(me%pidsort(foundpos)) == particle%particleid ) then
        ! write(logunit(1),*) "warning append_da_particle: particle already in array "
        return
      else
        ! need to append a new value!

        if (me%nvals == huge(me%nvals)) then
           write(logunit(1),*) "reached end of integer range for dynamic particle array!"
           write(logunit(1),*) "aborting!!"
           stop
        end if

        wasadded = .true.
        if (me%nvals == me%containersize) then
          ! container is full, need to expand it
          call expand_da_particle_mem(me = me, length = length)
        end if
        me%nvals = me%nvals + 1

        ! put the new value into the last position in the
        ! array.
        me%val(me%nvals) = particle
        me%pidlist(me%nvals) = particle%particleid
        do while( foundpos < me%nvals )
          if(me%pidlist(me%pidsort(foundpos)) /= particle%particleid) then
            exit
          end if
          ! in case of multiple entries with the same value
          ! move on to the first differing entry.
          foundpos = foundpos + 1
        end do
        ! shift the sorted list of indices, to create a
        ! whole for the value to be inserted, at position
        ! foundpos.
        do i=me%nvals-1,foundpos,-1
          me%pidsort(i+1) = me%pidsort(i)
        end do
        ! put the index of the new value into the
        ! sorted list at the now freed position.
        me%pidsort(foundpos) = me%nvals

      end if

    else

      ! value to append is larger than all existing ones,
      ! just put it to the end of the list, this captures
      ! also the case of empty lists.
      ! in this case foundpos = me%nvals + 1 holds.
      wasadded = .true.
      if (foundpos > me%containersize) then
        ! expand the array, if its boundary is reached
        call expand_da_particle_mem(me = me, length = length)
      end if
      me%nvals = foundpos
      me%val(foundpos) = particle
      me%pidlist(foundpos) = particle%particleid
      me%pidsort(foundpos) = foundpos

    end if
  end subroutine append_da_particle_mem


  subroutine expand_da_particle_mem(me, length)
    !------------------------------------------------------------------------
    type(dyn_particle_mem_array_type), intent(inout) :: me !< array to resize
    !> optional length to expand the array with
    integer, intent(in), optional :: length
    !------------------------------------------------------------------------
    type(mus_particle_mem_type), allocatable :: swpval(:)
    integer, allocatable :: swpidlist(:)
    integer, allocatable :: swpidsort(:)
    integer :: explen
    !------------------------------------------------------------------------

    ! if length is present, use that, otherwise double the size
    if( present( length ) ) then
      explen = length
    else
      ! set the global minimum length, if doubling would be smaller than that
      explen = max(me%containersize, minlength)
    end if


    ! check whether the new size will exceed the max container size.
    if( (me%containersize + explen) >= maxcontainersize ) then
      ! if so, expand to the maximum size
      me%containersize = maxcontainersize
    else
      ! if not, expand to the calculated size
      me%containersize = me%containersize + explen
    end if

    ! now make a larger array and copy all the current values into it.
    ! only need to copy values, if there are actually values to append.
    if (me%nvals > 0) then
      allocate(swpval(me%containersize))
      swpval(1:me%nvals) = me%val(1:me%nvals)
      call move_alloc( swpval, me%val )

      allocate(swpidlist(me%containersize))
      swpidlist(1:me%nvals) = me%pidlist(1:me%nvals)
      call move_alloc( swpidlist, me%pidlist )

      allocate(swpidsort(me%containersize))
      swpidsort(1:me%nvals) = me%pidsort(1:me%nvals)
      call move_alloc( swpidsort, me%pidsort )

    else ! me%nvals == 0
      if( allocated(me%val) ) then
        deallocate(me%val)
        allocate(me%val(me%containersize))
      end if

      if( allocated(me%pidlist) ) then
        deallocate(me%pidlist)
        allocate(me%pidlist(me%containersize))
      end if

      if( allocated(me%pidsort) ) then
        deallocate(me%pidsort)
        allocate(me%pidsort(me%containersize))
      end if
    end if

  end subroutine expand_da_particle_mem


  !> swaps the position of two particles in particle dynamic array
  !! new position of ielem1 = old position of ielem2 and vice-versa
  !! also updates the pidlist and pidsort arrays
  subroutine swap_da_particle_mem( me, ielem1, ielem2 )
    !> particle array to operate on
    type(dyn_particle_mem_array_type), intent(inout) :: me
    !> current index of one element
    integer, intent(in) :: ielem1
    !> current index of other element
    integer, intent(in) :: ielem2
    !------------------------------------------------------------------------
    type(mus_particle_mem_type) :: tmp
    integer :: tmppid
    integer :: k, k1, k2
    !------------------------------------------------------------------------

    k1 = -1
    k2 = -1

    if(ielem1 == ielem2) then
      return
    end if
    ! 1. update particle array and pidlist
    ! copy first element to tmp
    tmp = me%val(ielem1)
    tmppid = me%pidlist(ielem1)

    ! put second element at position of first
    me%val(ielem1) = me%val(ielem2)
    me%pidlist(ielem1) = me%pidlist(ielem2)

    ! put first element at old position of second
    me%val(ielem2) = tmp
    me%pidlist(ielem2) = tmppid

    ! 2. update the sorted list
    ! find the indices of pidsort pointing to ielem1 and ielem2
    do k = 1, me%nvals
      if (me%pidsort(k) == ielem1) then
        k1 = k
        exit
      end if
    end do

    do k = 1, me%nvals
      if (me%pidsort(k) == ielem2) then
        k2 = k
        exit
      end if
    end do

    ! swap the values
    me%pidsort(k1) = ielem2
    me%pidsort(k2) = ielem1

  end subroutine swap_da_particle_mem
  ! destroy particle dynamic array
  subroutine remove_particle_from_da_particle_mem(particles, ielem)
     !> particle group to operate on
     type(dyn_particle_mem_array_type), intent(inout) :: particles
     !> current index of element to remove
     integer, intent(in) :: ielem
     !------------------------------------------------------------------------
     integer :: ks, k
     !------------------------------------------------------------------------

     ! first check if this is a positive particle id, so indicates an actual particle
     if( particles%val(ielem)%particleid < 0 ) return

     ! place element to remove at the end of the array
     call swap_da_particle_mem( me = particles,          &
                         & ielem1 = ielem,          &
                         & ielem2 = particles%nvals )

     ! set particle id to negative to indicate this particle no longer belongs
     ! to dyn_array
     particles%val( particles%nvals )%particleid = &
       &  -1*particles%val( particles%nvals )%particleid

     ! update pidsort: first do linear search to find
     ! the element of pidsort corresponding to ielem
     searchloop: do ks = 1, particles%nvals
       if( particles%pidsort(ks) == particles%nvals ) then
         ! pidsort(ks) points to the element to remove
         do k = ks, particles%nvals - 1
           ! shift all elements to the right of ks one to the left.
           particles%pidsort(k) &
             & = particles%pidsort(k + 1)
         end do
         exit searchloop
       end if
     end do searchloop

     ! after successful removal of particle, decrease nvals
     particles%nvals = particles%nvals - 1

  end subroutine remove_particle_from_da_particle_mem
  !> truncate the dynamic particle array to only fit the actual entries
  subroutine truncate_da_particle_mem( me )
    !> particle array to operate on
    type(dyn_particle_mem_array_type), intent(inout) :: me
    !------------------------------------------------------------------------
    type(mus_particle_mem_type), allocatable :: swpval(:)
    !------------------------------------------------------------------------

    if (me%nvals < me%containersize) then
      allocate(swpval(me%nvals))

      swpval = me%val(:me%nvals)

      call move_alloc(swpval, me%val)

      me%containersize = me%nvals
    end if

  end subroutine truncate_da_particle_mem

  !> return the sorted position of a value in the given dynamic array
  !!
  !! if the value was not found,
  !!  - return 0 if nextifnotfound = .false.
  !!  - return position at the end if nextifnotfound = .true.
  function sortposofval_particle_mem(me, pid, nextifnotfound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_particle_mem_array_type), intent(in) :: me !< dynamic array
    integer, intent(in) :: pid !< particle id to look for
    !> flag to indicate, if the next entry in the list should be returned,
    !! if the searched one is not found.
    logical, intent(in), optional :: nextifnotfound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< position of val in the sorted particle id list, 0 if not found
    !------------------------------------------------------------------------
    logical :: retnext
    integer :: lb, ub
    integer :: mid
    integer :: lb_val, ub_val
    integer :: mid_val
    !------------------------------------------------------------------------

    retnext = .false.
    if (present(nextifnotfound)) retnext = nextifnotfound

    lb = 1
    ub = me%nvals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retnext) pos = lb

    !> binary search on sorted list
    do while(ub >= lb)
      lb_val = me%pidlist(me%pidsort(lb))

      ! if pid is smaller than smallest val in pidlist, exit
      if (pid < lb_val) then
        if (retnext) pos = lb
        exit
      end if

      ub_val = me%pidlist(me%pidsort(ub))

      ! also if pid is greater than greatest value in list, exit
      if (pid > ub_val) then
        if (retnext) pos = ub+1
        exit
      end if

      ! safe guard against integer limit overflow
      mid = lb + (ub-lb) / 2
      mid_val = me%pidlist(me%pidsort(mid))
      if (pid == mid_val) then
        pos = mid
        exit
      end if
      if (pid > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function sortposofval_particle_mem

  ! ************************************************************************ !
  !> Routine for allocating the existsOnProc, addToProc and removeFromProc
  !! masks used to determine when particles should be sent over to new processes
  !! or which processes need to receive position, velocity updates etc.
  subroutine allocateProcessMasks_MEM( particle, nProcs )
    !> Particle to initialize
    type(mus_particle_MEM_type), intent(inout) :: particle
    !> Number of processes to communicate particle data with
    integer :: nProcs
    ! -----------------------------------------------!
    ! Allocate space for the existsOnProc mask which tells us on which other procs
    ! this particle lives at the current time step
    if( allocated(particle%existsOnProc) ) deallocate( particle%existsOnProc )
    allocate( particle%existsOnProc( nProcs ) )
    particle%existsOnProc( 1:nProcs ) = .FALSE.

    ! addToProc is used to determine whether to send over data needed to add this
    ! particle to the receiving process's particle group
    if( allocated(particle%addToProc) ) deallocate( particle%addToProc )
    allocate( particle%addToProc( nProcs ) )
    particle%addToProc( 1:nProcs ) = .FALSE.

    ! removeFromProc is used to determine whether to send over the signal that
    ! a particle needs to be removed from the receiving proc's particle group
    if( allocated(particle%removeFromProc) ) deallocate( particle%removeFromProc )
    allocate( particle%removeFromProc( nProcs ) )
    particle%removeFromProc( 1:nProcs ) = .FALSE.

  end subroutine allocateProcessMasks_MEM
  ! ************************************************************************ !

end module mus_particle_MEM_type_module
