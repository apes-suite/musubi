! Copyright (c) 2013, 2015, 2017, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!> This module contains an interface for external C++ code to compute
!! liquid mixture property like thermodynamic factor and
!! Maxwell-Stefan Diffusivity coefficients
module mus_eNRTL_module
  use env_module,                  only: rk
  use tem_logging_module,          only: logUnit

  use, intrinsic :: iso_c_binding, only: c_char, c_int

  implicit none

  private

  public :: mus_init_eNRTL
  public :: mus_calc_thermFactor
  public :: mus_calc_MS_DiffMatrix

  interface mus_calc_thermFactor
    module procedure mus_calc_thermFactor_single
  end interface mus_calc_thermFactor

  interface mus_calc_MS_DiffMatrix
    module procedure mus_calc_MS_DiffMatrix_single
  end interface mus_calc_MS_DiffMatrix


contains


! **************************************************************************** !
  !> Dummy routine which sets thermodynamic factor matrix to diagonal matrix
  subroutine mus_calc_thermFactor_single( nFields, temp, press, mole_frac,     &
    &                                     therm_factors )
    ! -------------------------------------------------------------------------!
    !> number of fields in mixture
    integer, intent(in) :: nFields
    !> mixture temperature
    real(kind=rk), intent(in) :: temp
    !> mixture pressure
    real(kind=rk), intent(in) :: press
    !> mole fraction of all species of single element
    real(kind=rk), intent(in) :: mole_frac(nFields)
    !> thermodynamic factor matrix
    real(kind=rk), intent(out) :: therm_factors(nFields,nFields)
    ! -------------------------------------------------------------------------!
    integer :: iField
    ! -------------------------------------------------------------------------!

    therm_factors = 0.0_rk
    do iField = 1,nFields
      therm_factors(iField, iField) = 1.0_rk
    end do

    end subroutine mus_calc_thermFactor_single
! **************************************************************************** !


! **************************************************************************** !
  !> Dummy routine which sets diffusivity coeff matrix to diagonal matrix
  subroutine mus_calc_MS_DiffMatrix_single( nFields, temp, press, mole_dens,   &
    &                                     D_ij_out )
    ! -------------------------------------------------------------------------!
    !> number of fields in mixture
    integer, intent(in) :: nFields
    !> mixture temperature
    real(kind=rk), intent(in) :: temp
    !> mixture pressure
    real(kind=rk), intent(in) :: press
    !> mole density of all species of single element
    real(kind=rk), intent(in) :: mole_dens(nFields)
    !> thermodynamic factor matrix
    real(kind=rk), intent(out) :: D_ij_out(nFields,nFields)
    ! -------------------------------------------------------------------------!
    D_ij_out = 1.0_rk
  end subroutine mus_calc_MS_DiffMatrix_single
! **************************************************************************** !


! **************************************************************************** !
  !> Dummy function to init_enrtl
  function mus_init_eNRTL( filename, nFields ) result(success)
    ! -------------------------------------------------------------------------!
    character(kind=c_char), dimension(*) :: filename
    !> Number of fields defined in the property file
    integer, intent(out) :: nFields
    !> Result, indicating the status of encode
    logical :: success
    ! -------------------------------------------------------------------------!
    write(logUnit(1),*) 'Using eNRTL dummy module'
    success = .true.
    nFields = 0
  end function mus_init_eNRTL
! **************************************************************************** !

end module mus_eNRTL_module
