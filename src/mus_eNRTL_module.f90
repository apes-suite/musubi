! Copyright (c) 2013, 2015-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_float,      &
    &                                    c_int_least8_t, c_int_least64_t,      &
    &                                    c_char, c_loc, c_bool

  use env_module,                  only: rk

  implicit none

  private

  public :: mus_init_eNRTL
  public :: mus_calc_thermFactor
  public :: mus_calc_MS_DiffMatrix

  !> This function initialize eNRTL model by loading liquid mixture
  !! property from filename
  interface
    function init_enrtl_loc( filename, nSpc ) bind(c, name='init_enrtl')
    !function init_enrtl_f90( filename ) bind(c, name='_Z10init_enrtlPKc')
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: filename
      integer(kind=c_int), intent(out) :: nSpc
      !> Result, indicating the status of encode
      integer(kind=c_int) :: init_enrtl
    end function init_enrtl_loc
  end interface

  !> This routine calculates thermodynamic factor for given mole_frac
  !! of all species
  interface
    subroutine calc_therm_factor_loc( nSpc, Temp, Press, Mole_frac,            &
      & Therm_factors ) bind(c, name='calc_therm_factor_C')
      use, intrinsic :: iso_c_binding
      integer(kind=c_int), value, intent(in) :: nSpc
      real(kind=c_double), value, intent(in) :: Temp
      real(kind=c_double), value, intent(in) :: Press
      real(kind=c_double), dimension(*), intent(in) :: Mole_frac
      real(kind=c_double), dimension(*), intent(out) :: Therm_factors
    end subroutine calc_therm_factor_loc
  end interface

  !> This routine calculates Maxwell-Stefan diffusivity coeffcient Matrix
  !! for given mole_frac of all species
  interface
    subroutine calc_ms_diff_matrix_from_molefrac(nSpc, Temp, Press, Mole_frac, &
      & D_ij_out) bind(c, name='calc_ms_diff_matrix_from_molefrac_C')
      use, intrinsic :: iso_c_binding
      integer(kind=c_int), value, intent(in) :: nSpc
      real(kind=c_double), value, intent(in) :: Temp
      real(kind=c_double), value, intent(in) :: Press
      real(kind=c_double), dimension(*), intent(in) :: Mole_frac
      real(kind=c_double), dimension(*), intent(out) :: D_ij_out
    end subroutine calc_ms_diff_matrix_from_molefrac
  end interface

  !> This routine calculates Maxwell-Stefan diffusivity coeffcient Matrix
  !! for given mole_frac of all species
  interface
    subroutine calc_ms_diff_matrix_from_moledens(nSpc, Temp, Press, Mole_dens, &
      & D_ij_out) bind(c, name='calc_ms_diff_matrix_from_moledens_C')
      use, intrinsic :: iso_c_binding
      integer(kind=c_int), value, intent(in) :: nSpc
      real(kind=c_double), value, intent(in) :: Temp
      real(kind=c_double), value, intent(in) :: Press
      real(kind=c_double), dimension(*), intent(in) :: Mole_dens
      real(kind=c_double), dimension(*), intent(out) :: D_ij_out
    end subroutine calc_ms_diff_matrix_from_moledens
  end interface

  interface mus_calc_thermFactor
    module procedure mus_calc_thermFactor_single
  end interface mus_calc_thermFactor

  interface mus_calc_MS_DiffMatrix
    module procedure mus_calc_MS_DiffMatrix_single
  end interface

contains

  ! ************************************************************************** !
  !> This function loads property file using external c-function
  function mus_init_eNRTL( filename, nFields ) result(success)
    ! -------------------------------------------------------------------------!
    character(kind=c_char), dimension(*) :: filename
    !> number of fields in mixture
    integer, intent(out) :: nFields
    logical :: success
    ! -------------------------------------------------------------------------!
    integer(kind=c_int) :: nFields_c
    integer(kind=c_int) :: success_c
    ! -------------------------------------------------------------------------!

    success_c = init_eNRTL_loc( filename, nFields_c )
    if (success_c==0) then
      success = .true.
      nFields = nFields_c
    else
      success = .false.
      nFields = -1
    end if

  end function mus_init_eNRTL

  ! ************************************************************************** !
  !> This routine calculates thermodynamic factor for given mole_frac
  !! of all species for single element
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
    integer(kind=c_int) :: nFields_c
    real(kind=c_double) :: temp_c, press_c
    real(kind=c_double) :: mole_frac_c(nFields)
    real(kind=c_double):: therm_factors_c(nFields*nFields)
    integer :: iField, iField_2
    ! -------------------------------------------------------------------------!

    nFields_c = int(nFields, kind=c_int)
    temp_c = real(temp, kind=c_double)
    press_c = real(press, kind=c_double)
    do iField = 1,nFields
      mole_frac_c(iField) = real(mole_frac(iField), kind=c_double)
    end do

    call calc_therm_factor_loc( nFields_c, temp_c, press_c, mole_frac_c,       &
      &                         therm_factors_c )
    do iField = 1, nFields
      do iField_2 = 1, nFields
        therm_factors(iField, iField_2) =                                      &
          &                    therm_factors_c((iField-1)*nFields+ iField_2)
      end do
    end do
!    therm_factors = 0.0_rk
!    do iField = 1, nFields
!      therm_factors(iField, iField) = 1.0_rk
!    end do

!    write(*,*) 'Therm_factors_c ', therm_factors_c
!    write(*,*) 'Therm_factors '
!    do iField = 1, nFields
!    write(*,*) 'iField ', iField, ' factors ', therm_factors(iField,:)
!    end do
  end subroutine mus_calc_thermFactor_single

  ! ************************************************************************** !
  !> This routine calculates Diffusivity coefficients matrix for given mole_frac
  !! of all species for single element
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
    integer(kind=c_int) :: nFields_c
    real(kind=c_double) :: temp_c, press_c
    real(kind=c_double) :: mole_dens_c(nFields)
    real(kind=c_double):: D_ij_out_c(nFields*nFields)
    integer :: iField, iField_2
    ! -------------------------------------------------------------------------!

    nFields_c = int(nFields, kind=c_int)
    temp_c = real(temp, kind=c_double)
    press_c = real(press, kind=c_double)
    do iField = 1,nFields
      mole_dens_c(iField) = real(mole_dens(iField), kind=c_double)
    end do

!    write(*,*) 'nFields ', nFields
!    write(*,*) 'temp ', temp
!    write(*,*) 'press ', press
!    write(*,*) 'mole_dens ', mole_dens

    call calc_ms_diff_matrix_from_moledens( nFields_c, temp_c, press_c,        &
      & mole_dens_c, D_ij_out_c )
    do iField = 1, nFields
      do iField_2 = 1, nFields
        D_ij_out(iField, iField_2) = D_ij_out_c((iField-1)*nFields+ iField_2)
      end do
    end do
!    write(*,*) 'D_ij_out_c ', D_ij_out_c
!    write(*,*) 'D_ij_out '
!    do iField = 1, nFields
!      write(*,*) 'iField ', iField, ' factors ', D_ij_out(iField,:)
!    end do
  end subroutine mus_calc_MS_DiffMatrix_single

end module mus_eNRTL_module
