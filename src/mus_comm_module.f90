! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011, 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
! **************************************************************************** !
!> This module provides the definition and methods for musubi-specific
!! communication.
!! It includes the wrapper functions for the actual communication, which is
!! defined in [[tem_comm_module]].
!!
module mus_comm_module

  ! include treelm modules
  use mpi
  use env_module,          only: rk
  use tem_comm_module,     only: tem_communication_type, tem_commpattern_type
  use tem_comm_env_module, only: tem_comm_env_type

  implicit none

  private

  ! public :: mus_exchange
  public :: mus_init_longBuffers


contains


! **************************************************************************** !
  !> Wrapper around the actual communication, to avoid copy-in, copy-out by the
  !! Intel compiler. (At least the intel compiler on pigeon (v12.0) seems to do
  !! copying here, if a sub-array is passed to an assumed size dummy argument.
  !! Therefore we use this wrapping with an assumed shape dummy argument, so we
  !! can pass a complete field to the actual exchange which has an assumed size
  !! argument, without copying complete state field around, just for
  !! communication. Ugly, but it doesn't seem to have an impact on performance,
  !! and right it seems to be the most suitable solution.
  !!
  ! subroutine mus_exchange(send, recv, state, pattern, level, comm)
  !   ! ------------------------------------------------------------------------
  !   !>
  !   type(tem_communication_type), intent(inout) :: send
  !   !>
  !   type(tem_communication_type), intent(inout) :: recv
  !   !>
  !   real(kind=rk), intent(inout) :: state(:)
  !   !>
  !   type(tem_commPattern_type), intent(in) :: pattern
  !   !>
  !   integer :: level
  !   !> MPI communicator
  !   integer, intent(in) :: comm
  !   ! ------------------------------------------------------------------------

  !   call pattern%exchange_real( send = send, recv = recv, state = state, &
  !     &                         message_flag = level, comm = comm        )

  ! end subroutine mus_exchange
! **************************************************************************** !


! **************************************************************************** !
  !> Copy the element position in send and recv buffer to pos array
  !! in long type buffer
  !!
  subroutine mus_init_longBuffers( comm, pattern )
    ! --------------------------------------------------------------------------
    !>
    type( tem_communication_type ), intent(inout) :: comm
    !>
    type( tem_commPattern_type ), intent(in) :: pattern
    ! --------------------------------------------------------------------------
    integer :: iProc
    ! --------------------------------------------------------------------------

    ! Allocate long buffers and initialize comm buffer
    if( allocated( comm%buf_long )) deallocate( comm%buf_long )
    allocate( comm%buf_long( comm%nProcs ) )
    do iProc = 1, comm%nProcs
      call pattern%initBuf_long( me    = comm%buf_long( iProc ),               &
        &                        pos   = comm%elemPos( iProc )%val,            &
        &                        nVals = comm%nElemsProc( iProc ))
    end do

  end subroutine mus_init_longBuffers
! **************************************************************************** !


end module mus_comm_module
! **************************************************************************** !
