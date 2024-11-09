! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
! ****************************************************************************** !
!> This module assigns function pointer to calculate turbulent viscosity
!! according to turbulence model and scheme definition
module mus_turb_viscosity_module
  ! include treelm modules
  use env_module,          only: rk
  use tem_aux_module,      only: tem_abort

  use mus_scheme_header_module, only: mus_scheme_header_type
  use mus_turbulence_module,    only: mus_turbulence_type,        &
    &                                 mus_turbulence_config_type
  use mus_gradData_module,      only: mus_gradData_type, mus_Grad_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_WALE_module,          only: mus_turbVisc_WALE_2D, &
    &                                 mus_turbVisc_WALE_3D
  use mus_Vreman_module,        only: mus_turbVisc_Vreman_2D, &
    &                                 mus_turbVisc_Vreman_3D
  use mus_Smagorinsky_module,   only: mus_turbVisc_Smagorinsky_fromGradU3D, &
    & mus_turbVisc_Smagorinsky_fromGradU2D,             &
    & mus_turbVisc_Smagorinsky_fromGradU3D_incomp,      &
    & mus_turbVisc_Smagorinsky_fromGradU2D_incomp,      &
    & mus_turbVisc_Smagorinsky_fromPreColPDF_2D,        &
    & mus_turbVisc_Smagorinsky_fromPreColPDF_3D,        &
    & mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_2D, &
    & mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_3D

  implicit none
  private

  public :: mus_assign_turbVisc_ptr

contains

  ! ************************************************************************** !
  !> This routine assigns function pointer to compute turbulence viscosity
  !! based on turbulence model and scheme header definition
  subroutine mus_assign_turbVisc_ptr(turb, schemeHeader)
    ! --------------------------------------------------------------------------
    !> Scheme header information
    type(mus_scheme_header_type), intent(in) :: schemeHeader

    !> turbulence type
    type(mus_turbulence_type), intent(inout) :: turb
    ! --------------------------------------------------------------------------
    turb%calcVisc%fromGradU => mus_turbVisc_fromGradU_dummy
    turb%calcVisc%fromPreColPDF => mus_turbVisc_fromPreColPDF_dummy

    select case(trim(schemeHeader%layout))
    case ('d2q9')
      select case(trim(turb%config%model))
      case ('smagorinsky')
        select case (trim(schemeHeader%kind))
        case ('fluid')
          if (turb%config%compSR_fromPDF) then
            turb%calcVisc%fromPreColPDF &
              & => mus_turbVisc_Smagorinsky_fromPreColPDF_2D
          else
            turb%calcVisc%fromGradU => mus_turbVisc_Smagorinsky_fromGradU2D
          end if
        case ('fluid_incompressible')
          if (turb%config%compSR_fromPDF) then
            turb%calcVisc%fromPreColPDF                             &
              & => mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_2D
          else
            turb%calcVisc%fromGradU                        &
              & => mus_turbVisc_Smagorinsky_fromGradU2D_incomp
          end if
        case default
          call tem_abort('Error: Unknown scheme kind for turbulence '&
            &          //'viscosity ptr')
        end select

      case ('wale')
        turb%calcVisc%fromGradU => mus_turbVisc_WALE_2D
      case ('vreman')
        turb%calcVisc%fromGradU => mus_turbVisc_Vreman_2D
      end select

    case ('d3q15','d3q19','d3q27')
      select case(trim(turb%config%model))
      case ('smagorinsky')
        select case (trim(schemeHeader%kind))
        case ('fluid')
          if (turb%config%compSR_fromPDF) then
            turb%calcVisc%fromPreColPDF &
              & => mus_turbVisc_Smagorinsky_fromPreColPDF_3D
          else
            turb%calcVisc%fromGradU => mus_turbVisc_Smagorinsky_fromGradU3D
          end if
        case ('fluid_incompressible')
          if (turb%config%compSR_fromPDF) then
            turb%calcVisc%fromPreColPDF                             &
              & => mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_3D
          else
            turb%calcVisc%fromGradU                        &
              & => mus_turbVisc_Smagorinsky_fromGradU3D_incomp
          end if
        case default
          call tem_abort('Error: Unknown scheme kind for turbulence '&
            &          //'viscosity ptr')
        end select

      case ('wale')
        turb%calcVisc%fromGradU => mus_turbVisc_WALE_3D
      case ('vreman')
        turb%calcVisc%fromGradU => mus_turbVisc_Vreman_3D
      end select
    case default
      call tem_abort('Error: Unknown layout for turbulence viscosity ptr')
    end select
  end subroutine mus_assign_turbVisc_ptr
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Dummy function for turbulent viscosity from Gradu procedure
  subroutine mus_turbVisc_fromGradU_dummy(turbVisc, turbConfig, gradData, &
    & auxField, velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
    ! --------------------------------------------------------------------------
    !> output: turbulent viscosity
    real(kind=rk), intent(out) :: turbVisc(:)
    !> turbulence config contains oefficients
    type(mus_turbulence_config_type), intent(in) :: turbConfig
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Auxiliary field variable array
    real(kind=rk), intent(in) :: auxField(:)
    !> position of velocity components in auxField
    integer, intent(in) :: velPos(3)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> turbulence coefficients
    !> current level lattice element size
    real(kind=rk), intent(in) :: dxL
    !> current level lattice time step size
    real(kind=rk), intent(in) :: dtL
    !> Object that contains pointers to calculate gradients
    type(mus_Grad_type), intent(in) :: Grad
    ! --------------------------------------------------------------------------
    call tem_abort('DUMMY routine to calculate turbVisc from gradU')
  end subroutine mus_turbVisc_fromGradU_dummy
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Dummy function to compute turbulent viscosity from PDF
  subroutine mus_turbVisc_fromPreColPDF_dummy(turbVisc, turbConfig, state,  &
    & neigh, auxField, densPos, velPos, nSize, nSolve, nScalars, nAuxScalars,&
    & layout, dxL, dtL, viscKine)
    ! --------------------------------------------------------------------------
    !> output: turbulent viscosity
    real(kind=rk), intent(out) :: turbVisc(:)
    !> turbulence type is implicitly passed to access turbulence coefficients
    type(mus_turbulence_config_type), intent(in) :: turbConfig
    !> state array
    real(kind=rk), intent(in) :: state(:)
    !> neigh array to obtain precollision pdf
    integer, intent(in) :: neigh(:)
    !> Auxiliary field variable array
    real(kind=rk), intent(in) :: auxField(:)
    !> position of density in auxField
    integer, intent(in) :: densPos
    !> position of velocity components in auxField
    integer, intent(in) :: velPos(3)
    !> number of elements in state array
    integer, intent(in) :: nSize
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> number of scalars in state array
    integer, intent(in) :: nScalars
    !> number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> scheme layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> current level lattice element size
    real(kind=rk), intent(in) :: dxL
    !> current level lattice time step size
    real(kind=rk), intent(in) :: dtL
    !> Background kinematic viscosity divided by dtL
    real(kind=rk), intent(in) :: viscKine(:)
    ! --------------------------------------------------------------------------
    call tem_abort('DUMMY routine to calculate turbVisc from PreCol PDF')
  end subroutine mus_turbVisc_fromPreColPDF_dummy
  ! ************************************************************************** !

end module mus_turb_viscosity_module
