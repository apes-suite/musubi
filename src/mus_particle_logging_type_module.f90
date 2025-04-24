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
!> Routines containing data types for logging particle data.

module mus_particle_logging_type_module

  use env_module,           only: rk, long_k, &
    &                             labelLen, newunit
  use tem_dyn_array_module, only: dyn_longArray_type
  
  
  implicit none
  
  public :: mus_particle_logging_type
  public :: pgDebugLog
  public :: init_particle_logger
  
  type mus_particle_logging_type
    !> log unit to write to particle debug log file to
    integer :: lu
    !> file name of particle debug log
    character(len=32) :: lfile
  
  end type mus_particle_logging_type
  
  !> Type containing information on fluid elements to track using 
  !! calls to dumpdata() in mus_particle_tracking_module
  !! Can track lines and planes. For lines, only length1 and dir1 
  !! need to be prescribed whereas for planes we also need 
  !! length2 and dir2
  type mus_particle_debugtracking_type
    !> file name of particle debug log
    !! time step will be appended to this name
    character(len=labelLen) :: lfile
    !> Type of tracker
    !! Can be 'line' or 'plane'
    character(len=labelLen) :: trackerType
    !> list of positions of elements (in total list) to track
    type(dyn_longArray_type) :: elemList
    !> starting coordinate of region to track
    real(kind=rk) :: xstart(3)
    !> length of line to track in first direction
    real(kind=rk) :: length1
    !> direction of line to track in second direction
    !! e.g. (/ 1, 0, 0 /) is the x-direction
    integer :: dir1(3)
    !> length of line to track in second direction
    real(kind=rk) :: length2
    !> direction of line to track in second direction
    !! e.g. (/ 1, 0, 0 /) is the x-direction
    integer :: dir2(3)
    !> Logical to indicate whether tracker is active
    logical :: active = .FALSE.
  end type mus_particle_debugtracking_type
  
  type(mus_particle_logging_type), save :: pgDebugLog
  type(mus_particle_debugtracking_type), save :: debugTracker


contains


  subroutine printDebugTrackerData(debugTracker, logUnit)
    type(mus_particle_debugtracking_type) :: debugTracker
    integer :: logUnit
    ! --------------------------------------!
    integer :: i
    ! --------------------------------------!
    write(logUnit,*) '-------- PARTICLE DEBUG TRACKER DATA -------'
    write(logUnit,*) 'lfile = ', debugTracker%lfile
    write(logUnit,*) 'trackerType = ', debugTracker%trackerType
  
    ! ----- STARTING COORDINATE ----- !
    write(logUnit,*) 'xstart = ', debugTracker%xstart
    do i = 1,3
      write(logUnit,'(E10.3)', advance='no') debugTracker%xstart(i)
    end do
    write(logUnit,'(A)')
  
    ! ----- FIRST DIRECTION ----- !
    write(logUnit,'(A,E10.3)') 'length1 = ', debugTracker%length1
  
    write(logUnit,*) 'dir1 = '
    do i = 1,3
      write(logUnit,'(I5)', advance='no') debugTracker%dir1(i)
    end do
    write(logUnit,'(A)')
  
    if( trim(debugTracker%trackerType) == 'plane') then
      ! ----- SECOND DIRECTION ----- !
      write(logUnit,'(A,E10.3)') 'length2 = ', debugTracker%length2
  
      write(logUnit,*) 'dir2 = '
      do i = 1,3
        write(logUnit,'(I5)', advance='no')  debugTracker%dir2(i)
      end do
      write(logUnit,'(A)')
    end if
      
  end subroutine printDebugTrackerData
  
  subroutine init_particle_logger( logger, myRank )
    !> logging type to initialize
    type(mus_particle_logging_type), intent(inout) :: logger
    !> rank of this process, used to construct the log file name
    integer, intent(in) :: myRank
    ! --------------------------------------!
    logical :: fileExists
    ! --------------------------------------!
    write(logger%lfile, "(A16,I0.4,A4)" ) "particleGroupLog", myRank, ".txt"
    logger%lu = newunit()
    
    inquire(file=logger%lfile, exist=fileExists)
    if(fileExists) then
      open(logger%lu, file = logger%lfile, status = 'replace')
    else
      open(logger%lu, file = logger%lfile, status = 'new')
    end if
    write(logger%lu,*) "Particle Group Log rank ", myRank
    close(logger%lu)
  
  end subroutine init_particle_logger

end module mus_particle_logging_type_module
