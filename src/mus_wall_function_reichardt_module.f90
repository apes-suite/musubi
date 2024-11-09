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
!! computations relative to Reichardt profile.
!! Haussmann, Marc; BARRETO, Alejandro CLARO; KOUYI, Gislain LIPEME;
!! Rivière, Nicolas; Nirschl, Hermann; Krause, Mathias J. (2019):
!! Large-eddy simulation coupled with wall models for turbulent channel
!! flows at high Reynolds numbers with a lattice Boltzmann method
!! — Application to Coriolis mass flowmeter. In Computers & Mathematics
!! with Applications 78 (10), pp. 3285–3302. DOI: 10.1016/j.camwa.2019.04.033.
!!
!! The explicit power-law in terms of friction velocity given in Eq. 33 in
!! S. Wilhelm, J. Jacob, and P. Sagaut, "An explicit power-law-based wall
!! model for lattice Boltzmann method–Reynolds-averaged numerical simulations
!! of the flow around airfoils", Physics of Fluids 30, 065111 (2018)
!! https://doi.org/10.1063/1.5031764
!!
!! The model constants are chosen according to Werner and Wengle:
!! Wengle, H. and Werner, H. (1993) ‘Large-eddy Simulation of Turbulent Flow
!! Over Sharp-edged Obstacles in a Plate Channel’, (1985), pp. 192–199.
!!
!! author: Gregorio Gerardo Spinelli
module mus_wall_function_reichardt_module

  ! include treelm modules
  use env_module,                        only: rk
  use mus_wall_function_abstract_module, only: mus_wall_function_type

  implicit none

  private

  !! Constant parameters for Reichardt's law
  real(kind=rk), parameter :: vonKA = 0.4_rk
  real(kind=rk), parameter :: oneOvervonKA = 1.0_rk / vonKA
  public :: mus_wall_function_reichardt_type

  !> extend the abstract subclass mus_wall_function_type
  type, extends(mus_wall_function_type) :: mus_wall_function_reichardt_type

  contains

    !> function to get uPlus
    procedure, nopass :: get_uPlus

    !> function to apply the newon method
    ! F = uPlus_n - uPlus_(n+1)
    ! here we compute the derivative of uPlus_(n+1) w.r.t. u_tau
    ! this function computes the derivative of uPlus with respect to uTau
    procedure, nopass :: get_d_uPlus_d_uTau

  end type mus_wall_function_reichardt_type


contains


  !> function to get uPlus
  pure function get_uPlus( yPlus ) result( uPlus )

    !> yPlus
    real(kind=rk), intent(in) :: yPlus
    !> output: uPlus
    real(kind=rk) :: uPlus
    ! ------------------------------------------------------------------------------
    !uPlus = oneOvervonKA * log( 1.0_rk + vonKA*yPlus )     &
    !   &                 + 7.8 * ( 1.0_rk-exp(-yPlus/11._rk)  &
    !   &                 - (yPlus/11._rk) * exp(-yPlus/3._rk) )

    ! above expression simplified with sympy in python
    uPlus = -0.709090909090909_rk * yPlus * exp(-yPlus / 3._rk) &
      &     + log(vonKA*yPlus + 1._rk) + 7.8_rk              &
      &     - 7.8_rk * exp(-yPlus/11._rk) + oneOvervonKA

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
    real(kind=rk) :: inv_nu, yPlus, nu_sqr, y_inv_nu, A, B
    inv_nu = 1._rk / nu
    y_inv_nu = y * inv_nu
    yPlus = uTau * y_inv_nu
    nu_sqr = nu**2
    A = 0.709090909090909_rk * y_inv_nu
    B = exp( -yPlus / 3._rk)

    ! obtained as diff from sympy in python
    d_uPlus_d_uTau = vonKA * y / ( nu * ( vonKA * yPlus + 1._rk ) ) &
      & + A * ( exp( -yPlus / 11._rk) - B )                               &
      & + 0.236363636363636_rk * yPlus * y_inv_nu * B

  end function get_d_uPlus_d_uTau

end module mus_wall_function_reichardt_module
