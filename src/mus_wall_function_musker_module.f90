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
!! computations relative to Musker profile.
!! Haussmann, Marc; BARRETO, Alejandro CLARO; KOUYI, Gislain LIPEME;
!! Rivière, Nicolas; Nirschl, Hermann; Krause, Mathias J. (2019):
!! Large-eddy simulation coupled with wall models for turbulent channel
!! flows at high Reynolds numbers with a lattice Boltzmann method
!! — Application to Coriolis mass flowmeter. In Computers & Mathematics
!! with Applications 78 (10), pp. 3285–3302. DOI: 10.1016/j.camwa.2019.04.033.
!!
!! author: Gregorio Gerardo Spinelli
module mus_wall_function_musker_module

  ! include treelm modules
  use env_module,                        only: rk
  use mus_wall_function_abstract_module, only: mus_wall_function_type

  implicit none

  private

  public :: mus_wall_function_musker_type

  !> extend the abstract subclass mus_wall_function_type
  type, extends(mus_wall_function_type) :: mus_wall_function_musker_type

  contains

    !> function to get uPlus
    procedure, nopass :: get_uPlus

    !> function to apply the newon method
    ! F = uPlus_n - uPlus_(n+1)
    ! here we compute the derivative of uPlus_(n+1) w.r.t. u_tau
    ! this function computes the derivative of uPlus with respect to uTau
    procedure, nopass :: get_d_uPlus_d_uTau

  end type mus_wall_function_musker_type


contains


  !> function to get uPlus
  pure function get_uPlus( yPlus ) result( uPlus )

    !> yPlus
    real(kind=rk), intent(in) :: yPlus
    !> output: uPlus
    real(kind=rk) :: uPlus
    !uPlus = ( 5.424_rk * atan((2.0_rk*yPlus - 8.15_rk)/16.7_rk) &
    !  &        + log10( (yPlus + 10.6_rk)**9.6_rk               &
    !  &            / ( yPlus**2.0_rk - 8.15_rk*yPlus            &
    !  &                + 86.0_rk)**2.0_rk ) - 3.5072790194_rk   )

    ! above expression simplified with sympy in python
    uPlus = 0.434294481903252_rk * log( ( 0.0943396226415094_rk * yPlus + 1._rk )**9.6_rk /  &
      &     ( 0.0116279069767442_rk * yPlus**2 - 0.0947674418604651_rk * yPlus + 1_rk )**2 ) &
      &     + 5.424_rk * atan( 0.119760479041916_rk * yPlus - 0.488023952095808_rk )         &
      &     + 2.46666038465466_rk

  end function get_uPlus

  !> function to get the derivative of uPlus with respect to uTau
  pure function get_d_uPlus_d_uTau( y, uTau, nu ) result( d_uPlus_d_uTau )

    !> vertical distance from the wall
    real(kind=rk), intent(in) :: y
    !> uTau at iteration n
    real(kind=rk), intent(in) :: uTau
    !> dynamic viscosity
    real(kind=rk), intent(in) :: nu
    !> output: derivative of uPlus with respect to uTau
    real(kind=rk) :: d_uPlus_d_uTau
    ! ------------------------------------------------------------------------------
    real(kind=rk) :: uTau_y, yPlus, nu_sqr, uTau_y_sqr, nu_uTau_y, A, B, C, D, E
    uTau_y = uTau * y
    nu_uTau_y = uTau_y * nu
    yPlus = uTau_y / nu
    nu_sqr = nu**2
    uTau_y_sqr = uTau_y**2
    A = ( 1._rk + 0.0943396226415094_rk * yPlus )
    B = nu_sqr - 0.0947674418604651_rk * nu_uTau_y + 0.0116279069767442_rk * uTau_y_sqr
    C = ( nu - 0.245398773006135_rk * uTau_y )**2
    D = nu_sqr + 0.238167377819212_rk * C
    E = A**9.6_rk

    ! obtained as diff from sympy in python
    d_uPlus_d_uTau = 0.434294481903252_rk * y * ( 1.49571515501793_rk * nu_sqr        &
      & * E * B + D * ( nu * E * ( 0.18953488372093_rk * nu                           &
      & - 0.0465116279069767_rk * uTau_y ) + 0.905660377358491_rk * A**8.6_rk * B ) ) &
      & / ( nu * E * D * B )

  end function get_d_uPlus_d_uTau

end module mus_wall_function_musker_module
