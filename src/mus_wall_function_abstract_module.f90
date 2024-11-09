! Copyright (c) 2023 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> This module contains data types, function and routines for wall function
!! computations.
!!
!! author: Gregorio Gerardo Spinelli
module mus_wall_function_abstract_module
  ! include treelm modules
  use env_module,                 only: rk

  implicit none

  !> collection of properties of the wall function type
  type, abstract :: mus_wall_function_type

  contains

    !> function to get uPlus
    procedure(get_uPlus_interface), deferred, nopass :: get_uPlus

    !> function to apply the newon method
    ! F = uPlus_n - uPlus_(n+1)
    ! here we compute the derivative of uPlus_(n+1) w.r.t. u_tau
    ! this function computes the derivative of uPlus with respect to uTau
    procedure(get_d_uPlus_d_uTau_interface), deferred, nopass :: get_d_uPlus_d_uTau

  end type mus_wall_function_type


  interface
    !> function interface to get u_Plus
    pure function get_uPlus_interface( yPlus ) result( uPlus )
      import :: rk

      !> yPlus
      real(kind=rk), intent(in) :: yPlus
      !> output: uPlus
      real(kind=rk) :: uPlus

    end function get_uPlus_interface

    !> function interface to get the derivative of uPlus with respect to uTau
    pure function get_d_uPlus_d_uTau_interface( y, uTau, nu ) result( d_uPlus_d_uTau )
      import :: rk

      !> vertical distance from the wall
      real(kind=rk), intent(in) :: y
      !> uTau at iteration n
      real(kind=rk), intent(in) :: uTau
      !> dynamic viscosity
      real(kind=rk), intent(in) :: nu
      !> output: derivative of uPlus with respect to uTau
      real(kind=rk) :: d_uPlus_d_uTau

    end function get_d_uPlus_d_uTau_interface
  end interface

end module mus_wall_function_abstract_module
