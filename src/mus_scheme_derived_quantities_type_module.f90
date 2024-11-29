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
!> This module contains data types, function and routines for gradient
!! computation.
!!
!! author: Gregorio Gerardo Spinelli
module mus_scheme_derived_quantities_module
  ! include treelm modules
  use env_module,                 only: labelLen, rk
  use tem_aux_module,             only: tem_abort
  use tem_param_module,           only: div1_3, div1_9, div1_36, div4_9, rho0, &
    &                                   div1_18, div2_27, div1_54, div1_216,   &
    &                                   div8_27, cs2inv, cs4inv, rho0Inv
  use tem_compileconf_module,     only: vlen
  use tem_logging_module,   only: logUnit

  implicit none
  private

  public :: mus_scheme_derived_quantities_type
  public :: mus_assign_derived_functions_ptr

  !> collection of properties of the scheme derived quantities type
  type mus_scheme_derived_quantities_type

    !> function pointer to get pdf equilibrium from vel and density
    procedure(get_pdfEq), nopass, pointer :: pdfEq_ptr => null()
    !> function pointer to get pdf equilibrium in a specific direction
    procedure(get_pdfEq_iDir), nopass, pointer :: pdfEq_iDir_ptr => null()
    !> function pointer to get velocities from pdf
    procedure(get_vel_from_pdf), nopass, pointer :: vel_from_pdf_ptr => null()
    !> function pointer to get velocities from pdf VECTORIZED
    procedure(get_vel_from_pdf_vectorized), nopass, pointer :: vel_from_pdf_vectorized_ptr => null()
    !> function pointer to get momentum
    procedure(get_vector_from_vel_dens), nopass, pointer :: momentum_from_vel_dens_ptr => null()
    !> function pointer to get kinetic energy
    procedure(get_scalar_from_vel_dens), nopass, pointer :: kineticEnergy_from_vel_dens_ptr => null()
    !> function pointer to get 1/rho as a mask regardless incompressibility
    procedure(get_rho0Inv), nopass, pointer :: rho0Inv_ptr => null()

  end type mus_scheme_derived_quantities_type

  abstract interface
    !> function pointer to get pdf equilibrium from vel and density
    pure function get_pdfEq( rho, vel, QQ, cxDirRK, weight ) result( fEq )
      import :: rk

      !> density
      real(kind=rk), intent(in) :: rho
      !> velocity
      real(kind=rk), intent(in) :: vel(3)
      !> size of the stencil
      integer, intent(in) :: QQ
      !> velocity streaming normal along iDir
      real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
      !> weight along iDir
      real(kind=rk), optional, intent(in) :: weight(:)
      !> output is equilibrium pdf
      real(kind=rk) :: fEq(QQ)

    end function get_pdfEq

    !> function pointer to get pdf equilibrium from vel and density along a
    !  specific direction
    pure function get_pdfEq_iDir( rho, vel, iDir, cxDirRK, weight ) result( fEq )
      import :: rk

      !> density
      real(kind=rk), intent(in) :: rho
      !> velocity
      real(kind=rk), intent(in) :: vel(3)
      !> direction of the pdf
      integer, intent(in) :: iDir
      !> velocity streaming normal along iDir
      real(kind=rk), intent(in) :: cxDirRK(3)
      !> weight along iDir
      real(kind=rk), intent(in) :: weight
      !> output is equilibrium pdf
      real(kind=rk) :: fEq

    end function get_pdfEq_iDir

    !> function pointer to get pdf equilibrium from vel and density
    pure function get_vel_from_pdf( pdf, dens, cxDirRK ) result( vel )
      import :: rk

      !> pdf
      real(kind=rk), intent(in) :: pdf(:)
      !> density
      real(kind=rk), intent(in) :: dens
      !> velocity streaming normal along iDir
      real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
      !> velocity
      real(kind=rk) :: vel(3)

    end function get_vel_from_pdf

    !> function pointer to get pdf equilibrium from vel and density VECTORIZED
    pure function get_vel_from_pdf_vectorized( pdf, dens, cxDirRK, nSolve ) result( vel )
      import :: rk, vlen

      !> pdf
      real(kind=rk), intent(in) :: pdf(:,:)
      !> density
      real(kind=rk), intent(in) :: dens(:)
      !> velocity streaming normal along iDir
      real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
      !> number of element to compute
      integer, intent(in) :: nSolve
      !> velocity
      real(kind=rk) :: vel(3,vlen)

    end function get_vel_from_pdf_vectorized

    !> function pointer to get pdf equilibrium from vel and density
    pure function get_vector_from_vel_dens( vel, dens ) result( vector )
      import :: rk

      !> velocity
      real(kind=rk), intent(in) :: vel(:)
      !> density
      real(kind=rk), intent(in) :: dens
      !> momentum
      real(kind=rk) :: vector(3)

    end function get_vector_from_vel_dens

    !> function pointer to get pdf equilibrium from vel and density
    pure function get_scalar_from_vel_dens( vel, dens ) result( scalar )
      import :: rk

      !> velocity
      real(kind=rk), intent(in) :: vel(:)
      !> density
      real(kind=rk), intent(in) :: dens
      !> momentum
      real(kind=rk) :: scalar

    end function get_scalar_from_vel_dens

    !> function pointer to get 1/rho as a mask regardless incompressibility
    pure function get_rho0Inv( dens ) result( inv_rho0 )
      import :: rk

      !> density
      real(kind=rk), intent(in) :: dens
      !> inverse of density regardless compressibility
      real(kind=rk) :: inv_rho0

    end function get_rho0Inv
  end interface

contains

! ************************************************************************** !
!> This function assigns the pointers for the respective derived function in
!  terms of stencil and fluid type
function mus_assign_derived_functions_ptr(label_stencil, label_fluid) result(getQuantities)
  ! --------------------------------------------------------------------------
  !> Scheme header information
  character(len=labelLen), intent(in) :: label_stencil
  !> Fluid label information
  character(len=labelLen), intent(in) :: label_fluid
  !> getQuantities function
  type(mus_scheme_derived_quantities_type) :: getQuantities
  ! --------------------------------------------------------------------------
  getQuantities%pdfEq_ptr => null()
  getQuantities%pdfEq_iDir_ptr => null()
  getQuantities%vel_from_pdf_ptr => null()
  getQuantities%vel_from_pdf_vectorized_ptr => null()
  getQuantities%momentum_from_vel_dens_ptr => null()
  getQuantities%kineticEnergy_from_vel_dens_ptr => null()
  getQuantities%rho0Inv_ptr => null()

  if (trim(label_fluid) == 'fluid') then
    getQuantities%pdfEq_iDir_ptr => get_pdfEq_compressible_iDir
    getQuantities%momentum_from_vel_dens_ptr => get_momentum_from_vel_dens_compressible
    getQuantities%kineticEnergy_from_vel_dens_ptr => get_kineticEnergy_from_vel_dens_compressible
    getQuantities%rho0Inv_ptr => get_rho0Inv_compressible
  else if (trim(label_fluid) == 'fluid_incompressible') then
    getQuantities%pdfEq_iDir_ptr => get_pdfEq_incompressible_iDir
    getQuantities%momentum_from_vel_dens_ptr => get_momentum_from_vel_dens_incompressible
    getQuantities%kineticEnergy_from_vel_dens_ptr => get_kineticEnergy_from_vel_dens_incompressible
    getQuantities%rho0Inv_ptr => get_rho0Inv_incompressible
  else
    write(logUnit(1),*) 'fluid type = "', trim(label_fluid), '"'
    write(logUnit(1),*) "Warning: get_pdfEq_iDir not set for fluid type"
  end if

  select case (trim(label_stencil))
  case ('d2q9')
    if (trim(label_fluid) == 'fluid') then
      getQuantities%pdfEq_ptr => get_pdfEq_d2q9
      getQuantities%vel_from_pdf_ptr => get_vel_from_pdf_d2q9
      getQuantities%vel_from_pdf_vectorized_ptr => get_vel_from_pdf_d2q9_vectorized
    else if (trim(label_fluid) == 'fluid_incompressible') then
      getQuantities%pdfEq_ptr => get_pdfEq_incomp_d2q9
      getQuantities%vel_from_pdf_ptr => get_vel_from_pdf_d2q9_incompressible
      getQuantities%vel_from_pdf_vectorized_ptr => get_vel_from_pdf_d2q9_vectorized_incompressible
    else
      write(logUnit(1),*) 'stencil label = "', trim(label_stencil), '"'
      write(logUnit(1),*) 'fluid type = "', trim(label_fluid), '"'
      write(logUnit(1),*) "Warning: get_pdfEq not set for fluid type"
    end if
  case ('d3q19')
    if (trim(label_fluid) == 'fluid') then
      getQuantities%pdfEq_ptr => get_pdfEq_d3q19
      getQuantities%vel_from_pdf_ptr => get_vel_from_pdf_d3q19
      getQuantities%vel_from_pdf_vectorized_ptr => get_vel_from_pdf_d3q19_vectorized
    elseif (trim(label_fluid) == 'fluid_incompressible') then
      getQuantities%pdfEq_ptr => get_pdfEq_incomp_d3q19
      getQuantities%vel_from_pdf_ptr => get_vel_from_pdf_d3q19_incompressible
      getQuantities%vel_from_pdf_vectorized_ptr => get_vel_from_pdf_d3q19_vectorized_incompressible
    else
      write(logUnit(1),*) 'stencil label = "', trim(label_stencil), '"'
      write(logUnit(1),*) 'fluid type = "', trim(label_fluid), '"'
      write(logUnit(1),*) "Warning: get_pdfEq not set for fluid type"
    end if
  case ('d3q27')
    if (trim(label_fluid) == 'fluid') then
      getQuantities%pdfEq_ptr => get_pdfEq_d3q27
      getQuantities%vel_from_pdf_ptr => get_vel_from_pdf_d3q27
      getQuantities%vel_from_pdf_vectorized_ptr => get_vel_from_pdf_d3q27_vectorized
    elseif (trim(label_fluid) == 'fluid_incompressible') then
      getQuantities%pdfEq_ptr => get_pdfEq_incomp_d3q27
      getQuantities%vel_from_pdf_ptr => get_vel_from_pdf_d3q27_incompressible
      getQuantities%vel_from_pdf_vectorized_ptr => get_vel_from_pdf_d3q27_vectorized_incompressible
    else
      write(logUnit(1),*) 'stencil label = "', trim(label_stencil), '"'
      write(logUnit(1),*) 'fluid type = "', trim(label_fluid), '"'
      write(logUnit(1),*) "Warning: get_pdfEq not set for fluid type"
    end if
  case default
    if (trim(label_fluid) == 'fluid') then
      getQuantities%pdfEq_ptr => get_pdfEq_compressible
      getQuantities%vel_from_pdf_ptr => get_vel_from_pdf_compressible
      getQuantities%vel_from_pdf_vectorized_ptr => get_vel_from_pdf_compressible_vectorized
    elseif (trim(label_fluid) == 'fluid_incompressible') then
      getQuantities%pdfEq_ptr => get_pdfEq_incompressible
      getQuantities%vel_from_pdf_ptr => get_vel_from_pdf_incompressible
      getQuantities%vel_from_pdf_vectorized_ptr => get_vel_from_pdf_incompressible_vectorized
    else
      write(logUnit(1),*) 'stencil label = "', trim(label_stencil), '"'
      write(logUnit(1),*) 'fluid type = "', trim(label_fluid), '"'
      write(logUnit(1),*) "Warning: get_pdfEq not set for fluid type"
    end if
  end select

end function mus_assign_derived_functions_ptr
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density along a
!  specific direction
pure function get_pdfEq_incompressible_iDir( rho, vel, iDir, cxDirRK, weight ) &
  & result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> direction of the pdf
  integer, intent(in) :: iDir
  !> velocity streaming normal along iDir
  real(kind=rk), intent(in) :: cxDirRK(3)
  !> weight along iDir
  real(kind=rk), intent(in) :: weight
  !> output is equilibrium pdf
  real(kind=rk) :: fEq
  ! --------------------------------------------------------------------------
  fEq = weight * ( rho + rho0 * ( cs2inv * sum(cxDirRK(:)*vel(:))           &
    &          + ( sum(cxDirRK(:)*vel(:)) * sum(cxDirRK(:)*vel(:)) )        &
    &          * cs4inv * 0.5_rk - sum(vel(:) * vel(:)) * 0.5_rk * cs2inv ) )

end function get_pdfEq_incompressible_iDir
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density along a
!  specific direction
pure function get_pdfEq_compressible_iDir( rho, vel, iDir, cxDirRK, weight ) &
  & result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> direction of the pdf
  integer, intent(in) :: iDir
  !> velocity streaming normal along iDir
  real(kind=rk), intent(in) :: cxDirRK(3)
  !> weight along iDir
  real(kind=rk), intent(in) :: weight
  !> output is equilibrium pdf
  real(kind=rk) :: fEq
  ! --------------------------------------------------------------------------
  fEq = weight * rho * ( 1._rk + ( cs2inv * sum(cxDirRK(:)*vel(:))        &
    &        + ( sum(cxDirRK(:)*vel(:)) * sum(cxDirRK(:)*vel(:)) )        &
    &        * cs4inv * 0.5_rk - sum(vel(:) * vel(:)) * 0.5_rk * cs2inv ) )

end function get_pdfEq_compressible_iDir
! ************************************************************************** !


! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density along a
!  specific direction
pure function get_pdfEq_incompressible( rho, vel, QQ, cxDirRK, weight ) &
  & result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> size of the pdf
  integer, intent(in) :: QQ
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> weight along iDir
  real(kind=rk), optional, intent(in) :: weight(:)
  !> output is equilibrium pdf
  real(kind=rk) :: fEq(QQ)
  ! --------------------------------------------------------------------------
  integer :: iDir
  ! --------------------------------------------------------------------------
  do iDir = 1, QQ
    fEq(iDir) = weight(iDir) * ( rho + rho0 * ( cs2inv * sum(cxDirRK(:,iDir)*vel(:)) &
      &          + ( sum(cxDirRK(:,iDir)*vel(:)) * sum(cxDirRK(:,iDir)*vel(:)) )     &
      &          * cs4inv * 0.5_rk - sum(vel(:) * vel(:)) * 0.5_rk * cs2inv ) )
  end do
end function get_pdfEq_incompressible
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density along a
!  specific direction
pure function get_pdfEq_compressible( rho, vel, QQ, cxDirRK, weight ) &
& result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> size of the pdf
  integer, intent(in) :: QQ
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> weight along iDir
  real(kind=rk), optional, intent(in) :: weight(:)
  !> output is equilibrium pdf
  real(kind=rk) :: fEq(QQ)
  ! --------------------------------------------------------------------------
  integer :: iDir
  ! --------------------------------------------------------------------------
  do iDir = 1, QQ
    fEq = weight(iDir) * rho * ( 1._rk + ( cs2inv * sum(cxDirRK(:,iDir)*vel(:)) &
      &        + ( sum(cxDirRK(:,iDir)*vel(:)) * sum(cxDirRK(:,iDir)*vel(:)) )  &
      &        * cs4inv * 0.5_rk - sum(vel(:) * vel(:)) * 0.5_rk * cs2inv ) )
  end do
end function get_pdfEq_compressible
! ************************************************************************** !

! ************************************************************************** !
!> This function computes the sigma vector necessary to get the
!  equilibrium pdf from velocity and density for d2q9 stencil.
pure function get_sigma_d2q9( vel ) result( sigma )
! --------------------------------------------------------------------------
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> output is sigma vector
  real(kind=rk) :: sigma(9)
  ! --------------------------------------------------------------------------
  sigma(9) = vel(1) + vel(2)
  sigma(8) = vel(1) - vel(2)
  sigma(7) = 3._rk * sigma(9)
  sigma(6) = 3._rk * sigma(8)
  sigma(5) = 4.5_rk * sigma(9)**2
  sigma(4) = 4.5_rk * sigma(8)**2
  sigma(3) = 4.5_rk * vel(1)**2
  sigma(2) = 4.5_rk * vel(2)**2
  sigma(1) = div1_3 * ( sigma(2) + sigma(3) )

end function get_sigma_d2q9
! ************************************************************************** !

! ************************************************************************** !
!> This function computes the equilibrium pdf from velocity
!  and density for d2q9 stencil.
pure function get_pdfEq_d2q9( rho, vel, QQ, cxDirRK, weight ) &
& result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> size of the pdf
  integer, intent(in) :: QQ
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> weight along iDir
  real(kind=rk), optional, intent(in) :: weight(:)
  !> output is equilibrium pdf
  real(kind=rk) :: fEq(QQ)
  ! --------------------------------------------------------------------------
  real(kind=rk) :: sigma(9), rho_div_9, rho_div_36
  ! --------------------------------------------------------------------------
  sigma = get_sigma_d2q9(vel)
  rho_div_9 = div1_9 * rho
  rho_div_36 = div1_36 * rho

  fEq(1) = (- rho_div_9 * ( sigma(1) + 3._rk * vel(1) - sigma(3) - 1._rk )  )
  fEq(2) = (- rho_div_9 * ( sigma(1) + 3._rk * vel(2) - sigma(2) - 1._rk ) )
  fEq(3) = (rho_div_9 * ( 3._rk * vel(1) - sigma(1) + sigma(3) + 1._rk ) )
  fEq(4) = (rho_div_9 * ( 3._rk * vel(2) - sigma(1) + sigma(2) + 1._rk ) )

  fEq(5) = (- rho_div_36 * ( sigma(1) - sigma(5) + sigma(7) - 1._rk ) )
  fEq(6) = (- rho_div_36 * ( sigma(1) - sigma(4) + sigma(6) - 1._rk ) )
  fEq(7) = (rho_div_36 * ( sigma(4) - sigma(1) + sigma(6) + 1._rk ) )
  fEq(8) = (rho_div_36 * ( sigma(5) - sigma(1) + sigma(7) + 1._rk ) )

  fEq(9) = (- div4_9 * rho * ( sigma(1) - 1._rk ) )

end function get_pdfEq_d2q9
! ************************************************************************** !


! ************************************************************************** !
!> This function computes the incompressible equilibrium pdf from velocity
!  and density for d2q9 stencil.
pure function get_pdfEq_incomp_d2q9( rho, vel, QQ, cxDirRK, weight ) &
& result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> size of the pdf
  integer, intent(in) :: QQ
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> weight along iDir
  real(kind=rk), optional, intent(in) :: weight(:)
  !> output is equilibrium pdf
  real(kind=rk) :: fEq(QQ)
  ! --------------------------------------------------------------------------
  real(kind=rk) :: sigma(9), rho_div_9, rho_div_36, rho0_div_9, rho0_div_36
  ! --------------------------------------------------------------------------
  sigma = get_sigma_d2q9(vel)
  rho_div_9 = div1_9 * rho
  rho_div_36 = div1_36 * rho
  rho0_div_9 = div1_9 * rho0
  rho0_div_36 = div1_36 * rho0

  fEq(1) = (rho_div_9 - rho0_div_9 * ( sigma(1) + 3._rk * vel(1) - sigma(3) ) )
  fEq(2) = (rho_div_9 - rho0_div_9 * ( sigma(1) + 3._rk * vel(2) - sigma(2) ) )
  fEq(3) = (rho_div_9 + rho0_div_9 * ( 3._rk * vel(1) - sigma(1) + sigma(3) ) )
  fEq(4) = (rho_div_9 + rho0_div_9 * ( 3._rk * vel(2) - sigma(1) + sigma(2) ) )

  fEq(5) = (rho_div_36 - rho0_div_36 * ( sigma(1) - sigma(5) + sigma(7) ) )
  fEq(6) = (rho_div_36 - rho0_div_36 * ( sigma(1) - sigma(4) + sigma(6) ) )
  fEq(7) = (rho_div_36 + rho0_div_36 * ( sigma(4) - sigma(1) + sigma(6) ) )
  fEq(8) = (rho_div_36 + rho0_div_36 * ( sigma(5) - sigma(1) + sigma(7) ) )

  fEq(9) = (div4_9 * rho - div4_9 * rho0 * sigma(1) )

end function get_pdfEq_incomp_d2q9
! ************************************************************************** !


! ************************************************************************** !
!> This function computes the sigma vector necessary to get the
!  equilibrium pdf from velocity and density for d3q19 stencil.
pure function get_sigma_d3q19( vel ) result( sigma )
! --------------------------------------------------------------------------
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> output is sigma vector
  real(kind=rk) :: sigma(22)
  ! --------------------------------------------------------------------------
  sigma(22) = vel(1) + vel(2)
  sigma(21) = vel(1) - vel(2)
  sigma(20) = vel(1) + vel(3)
  sigma(19) = vel(1) - vel(3)
  sigma(18) = vel(2) + vel(3)
  sigma(17) = vel(2) - vel(3)
  sigma(16) = 3._rk * sigma(22)
  sigma(15) = 3._rk * sigma(21)
  sigma(14) = 3._rk * sigma(20)
  sigma(13) = 3._rk * sigma(19)
  sigma(12) = 3._rk * sigma(18)
  sigma(11) = 3._rk * sigma(17)
  sigma(10) = 4.5_rk * sigma(22)**2
  sigma(9) = 4.5_rk * sigma(21)**2
  sigma(8) = 4.5_rk * sigma(20)**2
  sigma(7) = 4.5_rk * sigma(19)**2
  sigma(6) = 4.5_rk * sigma(18)**2
  sigma(5) = 4.5_rk * sigma(17)**2
  sigma(4) = 4.5_rk * vel(1)**2
  sigma(3) = 4.5_rk * vel(2)**2
  sigma(2) = 4.5_rk * vel(3)**2
  sigma(1) = div1_3 * ( sigma(2) + sigma(3) + sigma(4) )

end function get_sigma_d3q19
! ************************************************************************** !

! ************************************************************************** !
!> This function computes the equilibrium pdf from velocity
!  and density for d3q19 stencil.
pure function get_pdfEq_d3q19( rho, vel, QQ, cxDirRK, weight ) &
& result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> size of the pdf
  integer, intent(in) :: QQ
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> weight along iDir
  real(kind=rk), optional, intent(in) :: weight(:)
  !> output is equilibrium pdf
  real(kind=rk) :: fEq(QQ)
  ! --------------------------------------------------------------------------
  real(kind=rk) :: sigma(22), rho_div_18, rho_div_36
  ! --------------------------------------------------------------------------
  sigma = get_sigma_d3q19(vel)
  rho_div_18 = div1_18 * rho
  rho_div_36 = div1_36 * rho

  fEq(1)  = (- rho_div_18 * ( 3._rk * vel(1) - sigma(4) + sigma(1) - 1._rk ) )
  fEq(2)  = (- rho_div_18 * ( 3._rk * vel(2) - sigma(3) + sigma(1) - 1._rk ) )
  fEq(3)  = (- rho_div_18 * ( 3._rk * vel(3) - sigma(2) + sigma(1) - 1._rk ) )
  fEq(4)  = (rho_div_18 * ( 3._rk * vel(1) + sigma(4) - sigma(1) + 1._rk ) )
  fEq(5)  = (rho_div_18 * ( 3._rk * vel(2) + sigma(3) - sigma(1) + 1._rk ) )
  fEq(6)  = (rho_div_18 * ( 3._rk * vel(3) + sigma(2) - sigma(1) + 1._rk ) )

  fEq(7)  = (rho_div_36 * ( sigma(6) - sigma(12) - sigma(1) + 1._rk ) )
  fEq(8)  = (rho_div_36 * ( sigma(5) - sigma(11) - sigma(1) + 1._rk ) )
  fEq(9)  = (rho_div_36 * ( sigma(5) + sigma(11) - sigma(1) + 1._rk ) )
  fEq(10) = (rho_div_36 * ( sigma(6) + sigma(12) - sigma(1) + 1._rk ) )

  fEq(11) = (rho_div_36 * ( sigma(8) - sigma(14) - sigma(1) + 1._rk ) )
  fEq(12) = (rho_div_36 * ( sigma(7) + sigma(13) - sigma(1) + 1._rk ) )
  fEq(13) = (rho_div_36 * ( sigma(7) - sigma(13) - sigma(1) + 1._rk ) )
  fEq(14) = (rho_div_36 * ( sigma(8) + sigma(14) - sigma(1) + 1._rk ) )

  fEq(15) = (rho_div_36 * ( sigma(10) - sigma(16) - sigma(1) + 1._rk ) )
  fEq(16) = (rho_div_36 * ( sigma(9) - sigma(15) - sigma(1) + 1._rk ) )
  fEq(17) = (rho_div_36 * ( sigma(9) + sigma(15) - sigma(1) + 1._rk ) )
  fEq(18) = (rho_div_36 * ( sigma(10) + sigma(16) - sigma(1) + 1._rk ) )

  fEq(19) = (- div1_3 * rho * ( sigma(1) - 1._rk ) )

end function get_pdfEq_d3q19
! ************************************************************************** !


! ************************************************************************** !
!> This function computes the incompressible equilibrium pdf from velocity
!  and density for d3q19 stencil.
pure function get_pdfEq_incomp_d3q19( rho, vel, QQ, cxDirRK, weight ) &
& result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> size of the pdf
  integer, intent(in) :: QQ
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> weight along iDir
  real(kind=rk), optional, intent(in) :: weight(:)
  !> output is equilibrium pdf
  real(kind=rk) :: fEq(QQ)
  ! --------------------------------------------------------------------------
  real(kind=rk) :: sigma(22), rho_div_18, rho_div_36, rho0_div_18, rho0_div_36
  ! --------------------------------------------------------------------------
  sigma = get_sigma_d3q19(vel)
  rho_div_18 = div1_18 * rho
  rho_div_36 = div1_36 * rho
  rho0_div_18 = div1_18 * rho0
  rho0_div_36 = div1_36 * rho0

  fEq(1)  = (rho_div_18 - rho0_div_18 * ( 3._rk * vel(1) - sigma(4) + sigma(1) ) )
  fEq(2)  = (rho_div_18 - rho0_div_18 * ( 3._rk * vel(2) - sigma(3) + sigma(1) ) )
  fEq(3)  = (rho_div_18 - rho0_div_18 * ( 3._rk * vel(3) - sigma(2) + sigma(1) ) )
  fEq(4)  = (rho_div_18 + rho0_div_18 * ( 3._rk * vel(1) + sigma(4) - sigma(1) ) )
  fEq(5)  = (rho_div_18 + rho0_div_18 * ( 3._rk * vel(2) + sigma(3) - sigma(1) ) )
  fEq(6)  = (rho_div_18 + rho0_div_18 * ( 3._rk * vel(3) + sigma(2) - sigma(1) ) )

  fEq(7)  = (rho_div_36 + rho0_div_36 * ( sigma(6) - sigma(12) - sigma(1) ) )
  fEq(8)  = (rho_div_36 + rho0_div_36 * ( sigma(5) - sigma(11) - sigma(1) ) )
  fEq(9)  = (rho_div_36 + rho0_div_36 * ( sigma(5) + sigma(11) - sigma(1) ) )
  fEq(10) = (rho_div_36 + rho0_div_36 * ( sigma(6) + sigma(12) - sigma(1) ) )

  fEq(11) = (rho_div_36 + rho0_div_36 * ( sigma(8) - sigma(14) - sigma(1) ) )
  fEq(12) = (rho_div_36 + rho0_div_36 * ( sigma(7) + sigma(13) - sigma(1) ) )
  fEq(13) = (rho_div_36 + rho0_div_36 * ( sigma(7) - sigma(13) - sigma(1) ) )
  fEq(14) = (rho_div_36 + rho0_div_36 * ( sigma(8) + sigma(14) - sigma(1) ) )

  fEq(15) = (rho_div_36 + rho0_div_36 * ( sigma(10) - sigma(16) - sigma(1) ) )
  fEq(16) = (rho_div_36 + rho0_div_36 * ( sigma(9)  - sigma(15) - sigma(1) ) )
  fEq(17) = (rho_div_36 + rho0_div_36 * ( sigma(9)  + sigma(15) - sigma(1) ) )
  fEq(18) = (rho_div_36 + rho0_div_36 * ( sigma(10) + sigma(16) - sigma(1) ) )

  fEq(19) = (div1_3 * rho - div1_3 * rho0 * sigma(1) )

end function get_pdfEq_incomp_d3q19
! ************************************************************************** !


! ************************************************************************** !
!> This function computes the sigma vector necessary to get the
!  equilibrium pdf from velocity and density for d3q27 stencil.
pure function get_sigma_d3q27( vel ) result( sigma )
! --------------------------------------------------------------------------
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> output is sigma vector
  real(kind=rk) :: sigma(34)
  ! --------------------------------------------------------------------------
  sigma(34) = vel(1) + vel(2)
  sigma(33) = vel(1) - vel(2)
  sigma(32) = vel(1) + vel(3)
  sigma(31) = vel(1) - vel(3)
  sigma(30) = vel(2) + vel(3)
  sigma(29) = vel(2) - vel(3)
  sigma(28) = vel(1) + vel(2) + vel(3)
  sigma(27) = vel(1) + vel(2) - vel(3)
  sigma(26) = vel(1) - vel(2) + vel(3)
  sigma(25) = vel(2) - vel(1) + vel(3)
  sigma(24) = 3._rk * sigma(34)
  sigma(23) = 3._rk * sigma(33)
  sigma(22) = 3._rk * sigma(32)
  sigma(21) = 3._rk * sigma(31)
  sigma(20) = 3._rk * sigma(30)
  sigma(19) = 3._rk * sigma(29)
  sigma(18) = 3._rk * sigma(28)
  sigma(17) = 3._rk * sigma(27)
  sigma(16) = 3._rk * sigma(26)
  sigma(15) = 3._rk * sigma(25)
  sigma(14) = 4.5_rk * sigma(34)**2
  sigma(13) = 4.5_rk * sigma(33)**2
  sigma(12) = 4.5_rk * sigma(32)**2
  sigma(11) = 4.5_rk * sigma(31)**2
  sigma(10) = 4.5_rk * sigma(30)**2
  sigma(9) = 4.5_rk * sigma(29)**2
  sigma(8) = 4.5_rk * sigma(28)**2
  sigma(7) = 4.5_rk * sigma(27)**2
  sigma(6) = 4.5_rk * sigma(26)**2
  sigma(5) = 4.5_rk * sigma(25)**2
  sigma(4) = 4.5_rk * vel(1)**2
  sigma(3) = 4.5_rk * vel(2)**2
  sigma(2) = 4.5_rk * vel(3)**2
  sigma(1) = div1_3 * ( sigma(2) + sigma(3) + sigma(4) )

end function get_sigma_d3q27
! ************************************************************************** !

! ************************************************************************** !
!> This function computes the equilibrium pdf from velocity
!  and density for d3q27 stencil.
pure function get_pdfEq_d3q27( rho, vel, QQ, cxDirRK, weight ) &
& result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> size of the pdf
  integer, intent(in) :: QQ
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> weight along iDir
  real(kind=rk), optional, intent(in) :: weight(:)
  !> output is equilibrium pdf
  real(kind=rk) :: fEq(QQ)
  ! --------------------------------------------------------------------------
  real(kind=rk) :: sigma(34), rho_div2_27, rho_div_54, rho_div_216
  ! --------------------------------------------------------------------------
  sigma = get_sigma_d3q27(vel)
  rho_div2_27 = div2_27 * rho
  rho_div_54 = div1_54 * rho
  rho_div_216 = div1_216 * rho

  fEq(1)  = (- rho_div2_27 * ( 3._rk * vel(1) - sigma(4) + sigma(1) - 1._rk ) )
  fEq(2)  = (- rho_div2_27 * ( 3._rk * vel(2) - sigma(3) + sigma(1) - 1._rk ) )
  fEq(3)  = (- rho_div2_27 * ( 3._rk * vel(3) - sigma(2) + sigma(1) - 1._rk ) )
  fEq(4)  = (rho_div2_27 * ( 3._rk * vel(1) + sigma(4) - sigma(1) + 1._rk ) )
  fEq(5)  = (rho_div2_27 * ( 3._rk * vel(2) + sigma(3) - sigma(1) + 1._rk ) )
  fEq(6)  = (rho_div2_27 * ( 3._rk * vel(3) + sigma(2) - sigma(1) + 1._rk ) )

  fEq(7)  = (rho_div_54 * ( sigma(10) - sigma(20) - sigma(1) + 1._rk ) )
  fEq(8)  = (rho_div_54 * ( sigma(9) - sigma(19) - sigma(1) + 1._rk ) )
  fEq(9)  = (rho_div_54 * ( sigma(9) + sigma(19) - sigma(1) + 1._rk ) )
  fEq(10) = (rho_div_54 * ( sigma(10) + sigma(20) - sigma(1) + 1._rk ) )

  fEq(11) = (rho_div_54 * ( sigma(12) - sigma(22) - sigma(1) + 1._rk ) )
  fEq(12) = (rho_div_54 * ( sigma(11) + sigma(21) - sigma(1) + 1._rk ) )
  fEq(13) = (rho_div_54 * ( sigma(11) - sigma(21) - sigma(1) + 1._rk ) )
  fEq(14) = (rho_div_54 * ( sigma(12) + sigma(22) - sigma(1) + 1._rk ) )

  fEq(15) = (rho_div_54 * ( sigma(14) - sigma(24) - sigma(1) + 1._rk ) )
  fEq(16) = (rho_div_54 * ( sigma(13) - sigma(23) - sigma(1) + 1._rk ) )
  fEq(17) = (rho_div_54 * ( sigma(13) + sigma(23) - sigma(1) + 1._rk ) )
  fEq(18) = (rho_div_54 * ( sigma(14) + sigma(24) - sigma(1) + 1._rk ) )

  fEq(19) = (- rho_div_216 * ( sigma(18) - sigma(8) + sigma(1) - 1._rk ) )
  fEq(20) = (- rho_div_216 * ( sigma(17) - sigma(7) + sigma(1) - 1._rk ) )
  fEq(21) = (- rho_div_216 * ( sigma(16) - sigma(6) + sigma(1) - 1._rk ) )
  fEq(22) = (rho_div_216 * ( sigma(15) + sigma(5) - sigma(1) + 1._rk ) )
  fEq(23) = (- rho_div_216 * ( sigma(15) - sigma(5) + sigma(1) - 1._rk ) )
  fEq(24) = (rho_div_216 * ( sigma(16) + sigma(6) - sigma(1) + 1._rk ) )
  fEq(25) = (rho_div_216 * ( sigma(17) + sigma(7) - sigma(1) + 1._rk ) )
  fEq(26) = (rho_div_216 * ( sigma(18) + sigma(8) - sigma(1) + 1._rk ) )

  fEq(27) = (- div8_27 * rho * ( sigma(1) - 1._rk ) )

end function get_pdfEq_d3q27
! ************************************************************************** !


! ************************************************************************** !
!> This function computes the incompressible equilibrium pdf from velocity
!  and density for d3q27 stencil.
pure function get_pdfEq_incomp_d3q27( rho, vel, QQ, cxDirRK, weight ) &
& result( fEq )
  ! --------------------------------------------------------------------------
  !> density
  real(kind=rk), intent(in) :: rho
  !> velocity
  real(kind=rk), intent(in) :: vel(3)
  !> size of the pdf
  integer, intent(in) :: QQ
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> weight along iDir
  real(kind=rk), optional, intent(in) :: weight(:)
  !> output is equilibrium pdf
  real(kind=rk) :: fEq(QQ)
  ! --------------------------------------------------------------------------
  real(kind=rk) :: sigma(34), rho_div2_27, rho_div_54, rho_div_216
  real(kind=rk) :: rho0_div2_27, rho0_div_54, rho0_div_216
  ! --------------------------------------------------------------------------
  sigma = get_sigma_d3q27(vel)
  rho_div2_27 = div2_27 * rho
  rho_div_54 = div1_54 * rho
  rho_div_216 = div1_216 * rho
  rho0_div2_27 = div2_27 * rho0
  rho0_div_54 = div1_54 * rho0
  rho0_div_216 = div1_216 * rho0

  fEq(1)  = (rho_div2_27 - rho0_div2_27 * ( 3._rk * vel(1) - sigma(4) + sigma(1) ) )
  fEq(2)  = (rho_div2_27 - rho0_div2_27 * ( 3._rk * vel(2) - sigma(3) + sigma(1) ) )
  fEq(3)  = (rho_div2_27 - rho0_div2_27 * ( 3._rk * vel(3) - sigma(2) + sigma(1) ) )
  fEq(4)  = (rho_div2_27 + rho0_div2_27 * ( 3._rk * vel(1) + sigma(4) - sigma(1) ) )
  fEq(5)  = (rho_div2_27 + rho0_div2_27 * ( 3._rk * vel(2) + sigma(3) - sigma(1) ) )
  fEq(6)  = (rho_div2_27 + rho0_div2_27 * ( 3._rk * vel(3) + sigma(2) - sigma(1) ) )

  fEq(7)  = (rho_div_54 + rho0_div_54 * ( sigma(10) - sigma(20) - sigma(1) ) )
  fEq(8)  = (rho_div_54 + rho0_div_54 * ( sigma(9)  - sigma(19) - sigma(1) ) )
  fEq(9)  = (rho_div_54 + rho0_div_54 * ( sigma(9)  + sigma(19) - sigma(1) ) )
  fEq(10) = (rho_div_54 + rho0_div_54 * ( sigma(10) + sigma(20) - sigma(1) ) )

  fEq(11) = (rho_div_54 + rho0_div_54 * ( sigma(12) - sigma(22) - sigma(1) ) )
  fEq(12) = (rho_div_54 + rho0_div_54 * ( sigma(11) + sigma(21) - sigma(1) ) )
  fEq(13) = (rho_div_54 + rho0_div_54 * ( sigma(11) - sigma(21) - sigma(1) ) )
  fEq(14) = (rho_div_54 + rho0_div_54 * ( sigma(12) + sigma(22) - sigma(1) ) )

  fEq(15) = (rho_div_54 + rho0_div_54 * ( sigma(14) - sigma(24) - sigma(1) ) )
  fEq(16) = (rho_div_54 + rho0_div_54 * ( sigma(13) - sigma(23) - sigma(1) ) )
  fEq(17) = (rho_div_54 + rho0_div_54 * ( sigma(13) + sigma(23) - sigma(1) ) )
  fEq(18) = (rho_div_54 + rho0_div_54 * ( sigma(14) + sigma(24) - sigma(1) ) )

  fEq(19) = (rho_div_216 - rho0_div_216 * ( sigma(18) - sigma(8) + sigma(1) ) )
  fEq(20) = (rho_div_216 - rho0_div_216 * ( sigma(17) - sigma(7) + sigma(1) ) )
  fEq(21) = (rho_div_216 - rho0_div_216 * ( sigma(16) - sigma(6) + sigma(1) ) )
  fEq(22) = (rho_div_216 + rho0_div_216 * ( sigma(15) + sigma(5) - sigma(1) ) )
  fEq(23) = (rho_div_216 - rho0_div_216 * ( sigma(15) - sigma(5) + sigma(1) ) )
  fEq(24) = (rho_div_216 + rho0_div_216 * ( sigma(16) + sigma(6) - sigma(1) ) )
  fEq(25) = (rho_div_216 + rho0_div_216 * ( sigma(17) + sigma(7) - sigma(1) ) )
  fEq(26) = (rho_div_216 + rho0_div_216 * ( sigma(18) + sigma(8) - sigma(1) ) )

  fEq(27) = (div8_27 * rho - div8_27 * rho0 * sigma(1) )

end function get_pdfEq_incomp_d3q27
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_compressible( pdf, dens, cxDirRK ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:)
  !> density
  real(kind=rk), intent(in) :: dens
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> velocity
  real(kind=rk) :: vel(3)
  ! --------------------------------------------------------------------------
  vel(1) = sum( pdf * cxDirRK(1,:) )
  vel(2) = sum( pdf * cxDirRK(2,:) )
  vel(3) = sum( pdf * cxDirRK(3,:) )
  vel = vel / dens

end function get_vel_from_pdf_compressible
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_compressible_vectorized( pdf, dens, cxDirRK, nSolve ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:,:)
  !> density
  real(kind=rk), intent(in) :: dens(:)
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> number of element to compute
  integer, intent(in) :: nSolve
  !> velocity
  real(kind=rk) :: vel(3,vlen)
  ! --------------------------------------------------------------------------
  integer :: iSolve
  ! --------------------------------------------------------------------------
  do iSolve = 1, nSolve
    vel(:,iSolve) = get_vel_from_pdf_compressible( pdf = pdf(:,iSolve), &
      &                                            dens = dens(iSolve),  &
      &                                            cxDirRK = cxDirRK    )
  end do
end function get_vel_from_pdf_compressible_vectorized
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_incompressible( pdf, dens, cxDirRK ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:)
  !> density
  real(kind=rk), intent(in) :: dens
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> velocity
  real(kind=rk) :: vel(3)
  ! --------------------------------------------------------------------------
  vel(1) = sum( pdf * cxDirRK(1,:) )
  vel(2) = sum( pdf * cxDirRK(2,:) )
  vel(3) = sum( pdf * cxDirRK(3,:) )
  vel = vel * rho0Inv

end function get_vel_from_pdf_incompressible
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_incompressible_vectorized( pdf, dens, cxDirRK, nSolve ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:,:)
  !> density
  real(kind=rk), intent(in) :: dens(:)
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> number of element to compute
  integer, intent(in) :: nSolve
  !> velocity
  real(kind=rk) :: vel(3,vlen)
  ! --------------------------------------------------------------------------
  integer :: iSolve
  ! --------------------------------------------------------------------------
  do iSolve = 1, nSolve
    vel(:,iSolve) = get_vel_from_pdf_incompressible( pdf = pdf(:,iSolve), &
      &                                              dens = dens(iSolve),  &
      &                                              cxDirRK = cxDirRK    )
  end do
end function get_vel_from_pdf_incompressible_vectorized
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for d2q9 stencil
pure function get_vel_from_pdf_d2q9( pdf, dens, cxDirRK ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:)
  !> density
  real(kind=rk), intent(in) :: dens
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> velocity
  real(kind=rk) :: vel(3)
  ! --------------------------------------------------------------------------
  vel(1) = pdf(3) - pdf(1) - pdf(5) - pdf(6) + pdf(7) + pdf(8)
  vel(2) = pdf(4) - pdf(2) - pdf(5) + pdf(6) - pdf(7) + pdf(8)
  vel(3) = 0._rk
  vel = vel / dens

end function get_vel_from_pdf_d2q9
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_d2q9_vectorized( pdf, dens, cxDirRK, nSolve ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:,:)
  !> density
  real(kind=rk), intent(in) :: dens(:)
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> number of element to compute
  integer, intent(in) :: nSolve
  !> velocity
  real(kind=rk) :: vel(3,vlen)
  ! --------------------------------------------------------------------------
  integer :: iSolve
  ! --------------------------------------------------------------------------
  do iSolve = 1, nSolve
    vel(:,iSolve) = get_vel_from_pdf_d2q9( pdf = pdf(:,iSolve), &
      &                                    dens = dens(iSolve),  &
      &                                    cxDirRK = cxDirRK    )
  end do
end function get_vel_from_pdf_d2q9_vectorized
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for d2q9 stencil
pure function get_vel_from_pdf_d2q9_incompressible( pdf, dens, cxDirRK ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:)
  !> density
  real(kind=rk), intent(in) :: dens
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> velocity
  real(kind=rk) :: vel(3)
  ! --------------------------------------------------------------------------
  vel(1) = pdf(3) - pdf(1) - pdf(5) - pdf(6) + pdf(7) + pdf(8)
  vel(2) = pdf(4) - pdf(2) - pdf(5) + pdf(6) - pdf(7) + pdf(8)
  vel(3) = 0._rk
  vel = vel * rho0Inv

end function get_vel_from_pdf_d2q9_incompressible
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_d2q9_vectorized_incompressible( pdf, dens, cxDirRK, nSolve ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:,:)
  !> density
  real(kind=rk), intent(in) :: dens(:)
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> number of element to compute
  integer, intent(in) :: nSolve
  !> velocity
  real(kind=rk) :: vel(3,vlen)
  ! --------------------------------------------------------------------------
  integer :: iSolve
  ! --------------------------------------------------------------------------
  do iSolve = 1, nSolve
    vel(:,iSolve) = get_vel_from_pdf_d2q9_incompressible( pdf = pdf(:,iSolve), &
      &                                                   dens = dens(iSolve),  &
      &                                                   cxDirRK = cxDirRK    )
  end do
end function get_vel_from_pdf_d2q9_vectorized_incompressible
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for d3q19 stencil
pure function get_vel_from_pdf_d3q19( pdf, dens, cxDirRK ) result( vel )
! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:)
  !> density
  real(kind=rk), intent(in) :: dens
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> velocity
  real(kind=rk) :: vel(3)
  ! --------------------------------------------------------------------------
  vel(1) = pdf(4) - pdf(1) - pdf(11) + pdf(12) - pdf(13) + pdf(14) - pdf(15) &
    &      - pdf(16) + pdf(17) + pdf(18)
  vel(2) = pdf(5) - pdf(2) - pdf(7) - pdf(8) + pdf(9) + pdf(10) - pdf(15) &
    &      + pdf(16) - pdf(17) + pdf(18)
  vel(3) = pdf(6) - pdf(3) - pdf(7) + pdf(8) - pdf(9) + pdf(10) - pdf(11) &
    &      - pdf(12) + pdf(13) + pdf(14)
  vel = vel / dens

end function get_vel_from_pdf_d3q19
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_d3q19_vectorized( pdf, dens, cxDirRK, nSolve ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:,:)
  !> density
  real(kind=rk), intent(in) :: dens(:)
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> number of element to compute
  integer, intent(in) :: nSolve
  !> velocity
  real(kind=rk) :: vel(3,vlen)
  ! --------------------------------------------------------------------------
  integer :: iSolve
  ! --------------------------------------------------------------------------
  do iSolve = 1, nSolve
    vel(:,iSolve) = get_vel_from_pdf_d3q19( pdf = pdf(:,iSolve), &
      &                                    dens = dens(iSolve),  &
      &                                    cxDirRK = cxDirRK    )
  end do
end function get_vel_from_pdf_d3q19_vectorized
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for d3q19 stencil
pure function get_vel_from_pdf_d3q19_incompressible( pdf, dens, cxDirRK ) result( vel )
! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:)
  !> density
  real(kind=rk), intent(in) :: dens
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> velocity
  real(kind=rk) :: vel(3)
  ! --------------------------------------------------------------------------
  vel(1) = pdf(4) - pdf(1) - pdf(11) + pdf(12) - pdf(13) + pdf(14) - pdf(15) &
    &      - pdf(16) + pdf(17) + pdf(18)
  vel(2) = pdf(5) - pdf(2) - pdf(7) - pdf(8) + pdf(9) + pdf(10) - pdf(15) &
    &      + pdf(16) - pdf(17) + pdf(18)
  vel(3) = pdf(6) - pdf(3) - pdf(7) + pdf(8) - pdf(9) + pdf(10) - pdf(11) &
    &      - pdf(12) + pdf(13) + pdf(14)
  vel = vel * rho0Inv

end function get_vel_from_pdf_d3q19_incompressible
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_d3q19_vectorized_incompressible( pdf, dens, cxDirRK, nSolve ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:,:)
  !> density
  real(kind=rk), intent(in) :: dens(:)
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> number of element to compute
  integer, intent(in) :: nSolve
  !> velocity
  real(kind=rk) :: vel(3,vlen)
  ! --------------------------------------------------------------------------
  integer :: iSolve
  ! --------------------------------------------------------------------------
  do iSolve = 1, nSolve
    vel(:,iSolve) = get_vel_from_pdf_d3q19_incompressible( pdf = pdf(:,iSolve), &
      &                                                    dens = dens(iSolve),  &
      &                                                    cxDirRK = cxDirRK    )
  end do
end function get_vel_from_pdf_d3q19_vectorized_incompressible
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for d3q27 stencil
pure function get_vel_from_pdf_d3q27( pdf, dens, cxDirRK ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:)
  !> density
  real(kind=rk), intent(in) :: dens
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> velocity
  real(kind=rk) :: vel(3)
  ! --------------------------------------------------------------------------
  vel(1) = pdf(4) - pdf(1) - pdf(11) + pdf(12) - pdf(13) + pdf(14) - pdf(15) &
    &      - pdf(16) + pdf(17) + pdf(18) - pdf(19) - pdf(20) - pdf(21)       &
    &      - pdf(22) + pdf(23) + pdf(24) + pdf(25) + pdf(26)
  vel(2) = pdf(5) - pdf(2) - pdf(7) - pdf(8) + pdf(9) + pdf(10) - pdf(15) &
    &      + pdf(16) - pdf(17) + pdf(18) - pdf(19) - pdf(20) + pdf(21)       &
    &      + pdf(22) - pdf(23) - pdf(24) + pdf(25) + pdf(26)
  vel(3) = pdf(6) - pdf(3) - pdf(7) + pdf(8) - pdf(9) + pdf(10) - pdf(11) &
    &      - pdf(12) + pdf(13) + pdf(14) - pdf(19) + pdf(20) - pdf(21)       &
    &      + pdf(22) - pdf(23) + pdf(24) - pdf(25) + pdf(26)
  vel = vel / dens

end function get_vel_from_pdf_d3q27
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_d3q27_vectorized( pdf, dens, cxDirRK, nSolve ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:,:)
  !> density
  real(kind=rk), intent(in) :: dens(:)
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> number of element to compute
  integer, intent(in) :: nSolve
  !> velocity
  real(kind=rk) :: vel(3,vlen)
  ! --------------------------------------------------------------------------
  integer :: iSolve
  ! --------------------------------------------------------------------------
  do iSolve = 1, nSolve
    vel(:,iSolve) = get_vel_from_pdf_d3q27( pdf = pdf(:,iSolve), &
      &                                     dens = dens(iSolve),  &
      &                                     cxDirRK = cxDirRK    )
  end do
end function get_vel_from_pdf_d3q27_vectorized
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for d3q27 stencil
pure function get_vel_from_pdf_d3q27_incompressible( pdf, dens, cxDirRK ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:)
  !> density
  real(kind=rk), intent(in) :: dens
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> velocity
  real(kind=rk) :: vel(3)
  ! --------------------------------------------------------------------------
  vel(1) = pdf(4) - pdf(1) - pdf(11) + pdf(12) - pdf(13) + pdf(14) - pdf(15) &
    &      - pdf(16) + pdf(17) + pdf(18) - pdf(19) - pdf(20) - pdf(21)       &
    &      - pdf(22) + pdf(23) + pdf(24) + pdf(25) + pdf(26)
  vel(2) = pdf(5) - pdf(2) - pdf(7) - pdf(8) + pdf(9) + pdf(10) - pdf(15) &
    &      + pdf(16) - pdf(17) + pdf(18) - pdf(19) - pdf(20) + pdf(21)       &
    &      + pdf(22) - pdf(23) - pdf(24) + pdf(25) + pdf(26)
  vel(3) = pdf(6) - pdf(3) - pdf(7) + pdf(8) - pdf(9) + pdf(10) - pdf(11) &
    &      - pdf(12) + pdf(13) + pdf(14) - pdf(19) + pdf(20) - pdf(21)       &
    &      + pdf(22) - pdf(23) + pdf(24) - pdf(25) + pdf(26)
  vel = vel * rho0Inv

end function get_vel_from_pdf_d3q27_incompressible
! ************************************************************************** !

! ************************************************************************** !
!> function pointer to get pdf equilibrium from vel and density for any stencil
pure function get_vel_from_pdf_d3q27_vectorized_incompressible( pdf, dens, cxDirRK, nSolve ) result( vel )
  ! --------------------------------------------------------------------------
  !> pdf
  real(kind=rk), intent(in) :: pdf(:,:)
  !> density
  real(kind=rk), intent(in) :: dens(:)
  !> velocity streaming normal along iDir
  real(kind=rk), optional, intent(in) :: cxDirRK(:,:)
  !> number of element to compute
  integer, intent(in) :: nSolve
  !> velocity
  real(kind=rk) :: vel(3,vlen)
  ! --------------------------------------------------------------------------
  integer :: iSolve
  ! --------------------------------------------------------------------------
  do iSolve = 1, nSolve
    vel(:,iSolve) = get_vel_from_pdf_d3q27_incompressible( pdf = pdf(:,iSolve), &
      &                                                    dens = dens(iSolve), &
      &                                                    cxDirRK = cxDirRK    )
  end do
end function get_vel_from_pdf_d3q27_vectorized_incompressible
! ************************************************************************** !


! ************************************************************************** !
!> function pointer to get momentum from vel and density for any stencil
pure function get_momentum_from_vel_dens_compressible( vel, dens ) result( vector )
! --------------------------------------------------------------------------
!> velocity
real(kind=rk), intent(in) :: vel(:)
!> density
real(kind=rk), intent(in) :: dens
!> momentum
real(kind=rk) :: vector(3)
! --------------------------------------------------------------------------
vector = vel * dens

end function get_momentum_from_vel_dens_compressible
! ************************************************************************** !


! ************************************************************************** !
!> function pointer to get momentum from vel and density for any stencil
pure function get_momentum_from_vel_dens_incompressible( vel, dens ) result( vector )
! --------------------------------------------------------------------------
!> velocity
real(kind=rk), intent(in) :: vel(:)
!> density
real(kind=rk), intent(in) :: dens
!> momentum
real(kind=rk) :: vector(3)
! --------------------------------------------------------------------------
vector = vel * rho0

end function get_momentum_from_vel_dens_incompressible
! ************************************************************************** !


! ************************************************************************** !
!> function pointer to get kineticEnergy from vel and density for any stencil
pure function get_kineticEnergy_from_vel_dens_compressible( vel, dens ) result( scalar )
! --------------------------------------------------------------------------
!> velocity
real(kind=rk), intent(in) :: vel(:)
!> density
real(kind=rk), intent(in) :: dens
!> kineticEnergy
real(kind=rk) :: scalar
! --------------------------------------------------------------------------
scalar = sum( vel(:)*vel(:) ) * 0.5_rk * dens

end function get_kineticEnergy_from_vel_dens_compressible
! ************************************************************************** !


! ************************************************************************** !
!> function pointer to get kineticEnergy from vel and density for any stencil
pure function get_kineticEnergy_from_vel_dens_incompressible( vel, dens ) result( scalar )
! --------------------------------------------------------------------------
!> velocity
real(kind=rk), intent(in) :: vel(:)
!> density
real(kind=rk), intent(in) :: dens
!> kineticEnergy
real(kind=rk) :: scalar
! --------------------------------------------------------------------------
scalar = sum( vel(:)*vel(:) ) * 0.5_rk * rho0

end function get_kineticEnergy_from_vel_dens_incompressible
! ************************************************************************** !


! ************************************************************************** !
!> function pointer to get kineticEnergy from vel and density for any stencil
pure function get_rho0Inv_compressible( dens ) result( inv_rho0 )
! --------------------------------------------------------------------------
!> density
real(kind=rk), intent(in) :: dens
!> kineticEnergy
real(kind=rk) :: inv_rho0
! --------------------------------------------------------------------------
inv_rho0 = 1._rk / dens

end function get_rho0Inv_compressible
! ************************************************************************** !


! ************************************************************************** !
!> function pointer to get kineticEnergy from vel and density for any stencil
pure function get_rho0Inv_incompressible( dens ) result( inv_rho0 )
! --------------------------------------------------------------------------
!> density
real(kind=rk), intent(in) :: dens
!> kineticEnergy
real(kind=rk) :: inv_rho0
! --------------------------------------------------------------------------
inv_rho0 = rho0Inv

end function get_rho0Inv_incompressible
! ************************************************************************** !

end module mus_scheme_derived_quantities_module
