! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014, 2018-2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> This module provides a data type for moment definition
!!
module mus_moments_type_module

  ! include treelm modules
  use env_module,       only: labelLen
  use tem_matrix_module,  only: tem_matrix_type

  implicit none
  private

  public :: mus_moment_type

  !> moment space definition
  type mus_moment_type
    !> is true if this type is already filled and no need
    !! to fill again after load balancing
    logical :: mom_ready = .false.

    !> transformation matrix from pdf space to moments
    type(tem_matrix_type) :: toMoments

    !> transformation matrix from moment space to pdf
    type(tem_matrix_type) :: toPDF

    !> Labels of the moments
    character(len=labelLen), allocatable :: momLabel(:)

    ! the followings are used by moment BCs
    !> position of first order moments in moments vector
    integer, allocatable :: first_moments(:)

    !> position of second order moments in moments vector
    integer, allocatable :: second_moments(:)

    !> position of third order moments in moments vector
    integer, allocatable :: third_moments(:)

    !> position of fourth order moments in moments vector
    integer, allocatable :: fourth_moments(:)

  end type mus_moment_type

end module mus_moments_type_module
! **************************************************************************** !
