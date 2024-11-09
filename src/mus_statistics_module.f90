! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> In this module we collect routines  and type definitions related to
!! statistics for the runtime of the code.
!!
!!
module mus_statistics_module
  ! include treelm modules
  use mpi
  use tem_logging_module,  only: logUnit

  implicit none
  private

  public :: mus_statistics_type
  public :: mus_calc_commAmount

  !> runtime statistics
  !! this includes memory consumption, amount of data to communciate,
  !! runtimes etc.
  type mus_statistics_type

    !> number of links to communicate
    !! Subset of nLinks_total
    integer :: nLinks_comm = 0

    !> number of total links available which theoretically can be communciated
    integer :: nLinks_total = 0

    !> Mean number of neighbor processes to which one process needs to send
    integer :: nProcs_send = 0

  end type mus_statistics_type


contains


! **************************************************************************** !
  !> Calculate the number of links to be communicated
  !!
  subroutine mus_calc_commAmount( stat, comm, nProcs )
    ! --------------------------------------------------------------------------
    !> runtime statistic
    type( mus_statistics_type), intent(inout) :: stat
    !>
    integer, intent(in) :: comm, nProcs
    ! --------------------------------------------------------------------------
    integer :: iErr
    integer :: nLinks_comm, nLinks_total
    ! --------------------------------------------------------------------------

    if ( nProcs > 1 ) then
      call mpi_reduce( stat%nLinks_comm, nLinks_comm, 1, mpi_integer,   &
        &              mpi_sum, 0, comm, iErr )
      call mpi_reduce( stat%nLinks_total, nLinks_total, 1, mpi_integer, &
        &              mpi_sum, 0, comm, iErr )
      stat%nLinks_comm  = nLinks_comm
      stat%nLinks_total = nLinks_total

      if( nLinks_total > 0 ) then
        write(logUnit(3),"(A,I0)") "Amount of message per proc: ", nLinks_comm &
          &                                                        / nProcs
        write(logUnit(1),"(A,F6.2)") 'Communication amount compressed to [%]:',&
          &          real( nLinks_comm )/real( nLinks_total ) * 100.0
      endif
    endif

  end subroutine mus_calc_commAmount
! **************************************************************************** !



end module mus_statistics_module
! **************************************************************************** !
