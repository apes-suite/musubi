! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> This module contains functions for initializing Cumulant relaxation paramaters
!! author: Gregorio Gerardo Spinelli
module mus_cumulantInit_module
  
  use env_module,               only: rk
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit

  implicit none
  private

  public :: cumulant_omega_check

  integer, parameter :: QQ = 27  !< number of pdf directions
  integer, parameter :: q000 = 27


contains

! ****************************************************************************** !
!> Checking the stability regions of omegas for the parametrized Cumulant
!! Just omega(2) is given in input. omega(2)=-1 means omega2=omegaBulk.
!! Limiters read from input. lim(N)=10^10 means unlimited.
!! lim(N) is for omega(N+2). Just omega(3:5) are limited as in the paper.
!! omega(6:10) = 1._rk
subroutine cumulant_omega_check( omegaVisc, omegaBulk, omegaIn, nSolve, level )
  ! -------------------------------------------------------------------- !
  !> vector of omegas in the level
  real(kind=rk), intent(in) :: omegaVisc(:)
  !> omega bulk value in the level
  real(kind=rk), intent(in) :: omegaBulk
  !> vector of omegas as given in musubi.lua
  real(kind=rk), intent(in) :: omegaIn(:)
  !> number of elements solved in kernel
  integer, intent(in) :: nSolve
  !> current level
  integer,intent(in) :: level
  ! -------------------------------------------------------------------- !
  real(kind=rk) :: omega(10), omegaPrev
  integer       :: iElem
  ! tests to print WARNING just once for omega 1, 3, 4, 5
  logical       :: test(5)
  ! ---------------------------------------------------------------------------
  write(logUnit(1),*) 'Level: ', level

  test = .true.

  omega = omegaIn
  if (omega(2) < 0._rk) then
    ! this is omegabulk
    omega(2) = omegaBulk
  endif
  if ( omega(2) <= 0._rk .or. omega(2) >= 2._rk ) then
    write(logUnit(1),*) 'ERROR: Omega_2 is not in the admissible range: ', omega(2)
    write(logUnit(1),*) 'Change bulk_viscosity value or give omega_2 = 1_rk as input' 
    call tem_abort('Please check omega_2 value')
  endif

  omegaPrev = -100_rk
  nodeloop: do iElem = 1, nSolve

    ! this is kinetic omega
    omega(1) = omegaVisc(iElem)
    if ( omega(1) == omegaPrev ) then
      CYCLE
    else
      omegaPrev = omega(1)
    endif
    if ( omega(1) <= 0._rk .or. omega(1) >= 2._rk ) then
      write(logUnit(1),*) 'ERROR: Omega_1 is not in the admissible range: ', omega(1)
      call tem_abort('Please check omega_1 value')
    endif
    if ( omega(1) == omega(2) ) then
      call tem_Abort('omega(1) should not be equal to omega(2), please increase bulk_viscosity value')
    endif
    if ( omega(1) > 1.9995_rk .and. test(1) ) then
      test(1) = .false.
      write(logUnit(1),*) 'WARNING: Omega_1 is really close to the stability limit of 2', omega(1)
    endif

    ! eq. 111 - 113
    omega(3) = (8._rk * (omega(1) - 2._rk) * (omega(2) * (3._rk * omega(1) - 1._rk) &
      & - 5._rk * omega(1) ) ) / (8._rk * (5._rk - 2._rk * omega(1)) * omega(1) &
      & + omega(2) * (8._rk + omega(1) * (9._rk * omega(1) - 26._rk) ) )
    if ( (omega(3) <= 0._rk .or. omega(3) >= 2._rk) .and. test(3) ) then
      test(3) = .false.
      write(logUnit(1),*) 'WARNING: Omega_3 is not in the admissible range: ', omega(3)
      write(logUnit(1),*) 'A solution is to run "cumulant" as relaxation or try to set omega(2) = 1._rk'
      write(logUnit(1),*) 'or to increase the bulk_viscosity value'
      write(logUnit(1),*) 'The code internally switches to relaxation = "cumulant" for the troublesome element'
      !call tem_abort('Please check omega_3 value.')
    endif

    omega(4) = (8._rk * (omega(1) - 2._rk) * ( omega(1) + omega(2) * (3._rk * omega(1) &
      & - 7._rk) ) ) / ( omega(2) * (56._rk - 42._rk * omega(1) + 9._rk * omega(1)**2 ) &
      & - 8._rk + omega(1) )
    if ( (omega(4) <= 0._rk .or. omega(4) >= 2._rk) .and. test(4) ) then
      test(4) = .false.
      write(logUnit(1),*) 'WARNING: Omega_4 is not in the admissible range: ', omega(4)
      write(logUnit(1),*) 'A solution is to run "cumulant" as relaxation or try to set omega(2) = 1._rk'
      write(logUnit(1),*) 'or to increase the bulk_viscosity value'
      write(logUnit(1),*) 'The code internally switches to relaxation = "cumulant" for the troublesome element'
      !call tem_abort('Please check omega_4 value.')
    endif

    omega(5) = ( 24._rk * ( omega(1) - 2._rk ) * ( 4._rk * omega(1)**2 + omega(1) &
      & * omega(2) * ( 18._rk - 13._rk * omega(1) ) + omega(2)**2 * ( 2._rk + omega(1) &
      & * ( 6._rk * omega(1) - 11._rk ) ) ) ) / ( 16._rk * omega(1)**2 * ( omega(1) - 6._rk ) &
      & - 2._rk * omega(1) * omega(2) * ( 216._rk + 5._rk * omega(1) * (9._rk * omega(1) &
      & - 46._rk ) ) + omega(2)**2 * ( omega(1) * ( 3._rk * omega(1) - 10._rk ) * (15._rk &
      & * omega(1) - 28._rk ) - 48._rk ) )
    if ( (omega(5) <= 0._rk .or. omega(5) >= 2._rk) .and. test(5) ) then
      test(5) = .false.
      write(logUnit(1),*) 'WARNING: Omega_5 is not in the admissible range: ', omega(5)
      write(logUnit(1),*) 'A solution is to run "cumulant" as relaxation or try to set omega(2) = 1._rk'
      write(logUnit(1),*) 'or to increase the bulk_viscosity value'
      write(logUnit(1),*) 'The code internally switches to relaxation = "cumulant" for the troublesome element'
      !call tem_abort('Please check omega_5 value.')
    endif

  end do nodeloop

end subroutine cumulant_omega_check
! ****************************************************************************** !

end module mus_cumulantInit_module
