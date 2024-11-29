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
!! computations relative to Schmitt profile.
!! Ref to following paper for Schmitt three layer equations.
!! Haussmann, M. et al. (2019) ‘Large-eddy simulation coupled with wall models
!! for turbulent channel flows at high Reynolds numbers with a lattice
!! Boltzmann method — Application to Coriolis mass flowmeter’, Computers &
!! Mathematics with Applications. Elsevier Ltd, 78(10), pp. 3285–3302.
!!
!! author: Gregorio Gerardo Spinelli
module mus_wall_function_schmitt_module

  ! include treelm modules
  use env_module,                        only: rk
  use mus_wall_function_abstract_module, only: mus_wall_function_type

  implicit none

  private

  !! Constant parameters for Schmitt's law
  real(kind=rk), parameter :: vonKA = 0.4_rk
  real(kind=rk), parameter :: sc_uLmt = 30._rk
  real(kind=rk), parameter :: sc_lLmt = 5._rk

  ! we need the limit to be public due to the function that computer the
  ! friction velocity created in mus_turb_wallFunc_module
  ! we need tot ake care of explicit and implicit part of the wall profile
  public :: sc_uLmt, sc_lLmt

  public :: mus_wall_function_schmitt_type

  ! function used to get friction velocity where the Schmitt profile
  ! has an explicit function
  public :: get_uTau_subVisousLayer
  public :: get_uTau_logLayer

  !> extend the abstract subclass mus_wall_function_type
  type, extends(mus_wall_function_type) :: mus_wall_function_schmitt_type

  contains

    !> function to get uPlus
    procedure, nopass :: get_uPlus

    !> function to apply the newton method
    ! F = uPlus_n - uPlus_(n+1)
    ! here we compute the derivative of uPlus_(n+1) w.r.t. u_tau
    ! this function computes the derivative of uPlus with respect to uTau
    procedure, nopass :: get_d_uPlus_d_uTau

  end type mus_wall_function_schmitt_type


contains


  !> function to get uPlus
  pure function get_uPlus( yPlus ) result( uPlus )

    !> yPlus
    real(kind=rk), intent(in) :: yPlus
    !> output is uPlus
    real(kind=rk) :: uPlus
    ! ------------------------------------------------------------------------------
    ! calculated with sympy
    if (yPlus >= sc_uLmt) then
      ! log layer, use powerlaw profile from werner and wengle
      uPlus = 8.3_rk * yPlus**0.142857142857143_rk
    else if ( yPlus < sc_lLmt) then
      ! viscous sublayer, use linear profile.
      uPlus = yPlus
    else ! if ( yPlus >= sc_lLmt .and. yPlus < sc_uLmt)
      ! Buffer layer, use logarithmic profile
      uPlus = 0.242384365402773_rk * ( 19.8872196713232_rk * vonKA + &
        &     ( 0.46051701859881_rk * vonKA + 3.40119738166216_rk )  &
        &       * log( yPlus ) - 5.47401601371867_rk ) / vonKA
    end if

  end function get_uPlus

  !> function to get the derivative of uPlus with respect to uTau
  pure function get_d_uPlus_d_uTau( y, uTau, nu ) result( d_uPlus_d_uTau )

    !> vertical distance from the wall
    real(kind=rk), intent(in) :: y
    !> uTau at iteration n
    real(kind=rk), intent(in) :: uTau
    !> dynamic viscosity
    real(kind=rk), intent(in) :: nu
    !> output is derivative of uPlus with respect to uTau
    real(kind=rk) :: d_uPlus_d_uTau
    ! ------------------------------------------------------------------------------

    ! obtained as diff from sympy in python
    ! Buffer layer, use logarithmic profile
    d_uPlus_d_uTau = 0.242384365402773_rk * ( 0.46051701859881_rk * vonKA &
      &               + 3.40119738166216_rk ) / ( vonKA * uTau )

  end function get_d_uPlus_d_uTau

  ! function used to get friction velocity where the Schmitt profile
  ! has an explicit function - viscous sub layer
  pure function get_uTau_subVisousLayer ( visc_div_dist, velSW ) result( uTau )

    !> dynamic viscosity divided by vertical distance from the wall
    real(kind=rk), intent(in) :: visc_div_dist
    !> velocity stream-wise parallel to wall
    real(kind=rk), intent(in) :: velSW
    !> friction velocity
    real(kind=rk) :: uTau
    ! ------------------------------------------------------------------------------

    ! obtained as diff from sympy in python
    ! Buffer layer, use logarithmic profile
    uTau = sqrt( velSW * visc_div_dist )

  end function get_uTau_subVisousLayer

  ! function used to get friction velocity where the Schmitt profile
  ! has an explicit function - log layer
  pure function get_uTau_logLayer ( visc_div_dist, velSW ) result( uTau )

    !> dynamic viscosity divided by vertical distance from the wall
    real(kind=rk), intent(in) :: visc_div_dist
    !> velocity stream-wise parallel to wall
    real(kind=rk), intent(in) :: velSW
    !> friction velocity
    real(kind=rk) :: uTau
    ! ------------------------------------------------------------------------------

    ! obtained as diff from sympy in python
    ! Buffer layer, use logarithmic profile
    uTau = 0.156966389612661_rk *                                &
      &             ( visc_div_dist**0.142857142857143_rk * velSW )**0.875_rk

  end function get_uTau_logLayer

end module mus_wall_function_schmitt_module
