! Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!> This module contains function to compute eddy viscosity for
!! Wall-Adapting Local Eddy-Viscosity turbulence
!! model.
!! This implementation follows the LES described by Weickert et al.
!! Weickert, M., Teike, G., Schmidt, O., & Sommerfeld, M. (2010).
!! Investigation of the LES WALE turbulence model within the lattice Boltzmann
!! framework. Computers and Mathematics with Applications, 59(7), 2200–2214.
!! author: Kannan Masilamani
module mus_WALE_module
  ! treelm modules
  use env_module,                    only: rk, eps
  use tem_param_module,              only: div1_2, div1_3
  use tem_compileconf_module,        only: vlen

  ! musubi modules
  use mus_turbulence_module,         only: mus_turbulence_config_type
  use mus_gradData_module,           only: mus_gradData_type, mus_Grad_type

  implicit none
  private

  public :: mus_turbVisc_WALE_3D
  public :: mus_turbVisc_WALE_2D

contains

  ! ************************************************************************** !
  !> Calculate eddy viscosity with WALE (Wall-Adapting Local Eddy-viscosity)
  !! model
  !! \todo add reference and formula
  subroutine mus_turbVisc_WALE_3D(turbVisc, turbConfig, gradData, auxField,  &
    &                             velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
    ! --------------------------------------------------------------------------
    !> output: turbulent viscosity
    real(kind=rk), intent(out) :: turbVisc(:)
    !> Contains turbulenct coefficients
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
    !> current level lattice element size
    real(kind=rk), intent(in) :: dxL
    !> current level lattice time step size
    real(kind=rk), intent(in) :: dtL
    !> Object that contains pointers to calculate gradients
    type(mus_Grad_type), intent(in) :: Grad
    ! --------------------------------------------------------------------------
    integer :: iElem
    real(kind=rk) :: SR(6), Sd(6), oneThird_trSd, Sd_sqr, SR_sqr, OP1, OP2
    real(kind=rk) :: visc_coeff
    !> gradient of velocity
    real(kind=rk):: gradU(3,3,vlen)
    real(kind=rk):: gradU_sqr(3,3,vlen)
    integer :: ndims
    integer :: nChunks, iChunks, nChunkElems, low_bound, elempos
    ! --------------------------------------------------------------------------
    ! viscosity coeff
    visc_coeff = (turbConfig%coeff%C_w * dxL )**2
    nDims = 3
    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      gradU(:,:,1:nChunkElems) = Grad%U_ptr( auxField    = auxField,    &
        &                                    gradData    = gradData,    &
        &                                    velPos      = velPos,      &
        &                                    nAuxScalars = nAuxScalars, &
        &                                    nDims       = nDims,       &
        &                                    nSolve      = nChunkElems, &
        &                                    elemOffset  = low_bound    )

      do iElem = 1, nChunkElems
        gradU_sqr(:,:,iElem) = matmul(gradU(:,:,iElem), gradU(:,:,iElem))
      end do !iElem

      do iElem = 1, nChunkElems
        ! traceless symmetric part of the square of the velocity gradient tensor
        ! Sd_ij = 1/2(du_k/dx_i du_j/dx_k + du_k/dx_j du_i/dx_k)
        !       - 1/3\delta_ij du_k/dx_l du_l/dx_k
        oneThird_trSd = (gradU_sqr(1,1,iElem) + gradU_sqr(2,2,iElem)  &
          &             + gradU_sqr(3,3,iElem))*div1_3
        Sd(1) = gradU_sqr(1,1,iElem) - oneThird_trSd                  !XX
        Sd(2) = gradU_sqr(2,2,iElem) - oneThird_trSd                  !YY
        Sd(3) = gradU_sqr(3,3,iElem) - oneThird_trSd                  !ZZ
        Sd(4) = 0.5_rk*(gradU_sqr(1,2,iElem)+gradU_sqr(2,1,iElem))    !XY
        Sd(5) = 0.5_rk*(gradU_sqr(2,3,iElem)+gradU_sqr(3,2,iElem))    !YZ
        Sd(6) = 0.5_rk*(gradU_sqr(1,3,iElem)+gradU_sqr(3,1,iElem))    !XZ

        ! double inner product of Sd: Sd_ij Sd_ij
        Sd_sqr = Sd(1)**2 + Sd(2)**2 + Sd(3)**2            &
          &    + 2.0_rk * ( Sd(4)**2 + Sd(5)**2 + Sd(6)**2 )

        ! symmetric strain rate tensors
        SR(1) = gradU(1,1,iElem)                               !XX
        SR(2) = gradU(2,2,iElem)                               !YY
        SR(3) = gradU(3,3,iElem)                               !ZZ
        SR(4) = (gradU(1,2,iElem)+gradU(2,1,iElem))*0.5_rk     !XY
        SR(5) = (gradU(2,3,iElem)+gradU(3,2,iElem))*0.5_rk     !YZ
        SR(6) = (gradU(1,3,iElem)+gradU(3,1,iElem))*0.5_rk     !XZ

        ! double inner product of tensor
        SR_sqr = SR(1)**2 + SR(2)**2 + SR(3)**2            &
          &    + 2.0_rk * ( SR(4)**2 + SR(5)**2 + SR(6)**2 )

        ! sub-grid scale kinetic energy
        ! k_sgs = (C_w^2 * dx /C_k)^2 (OP1/OP2)^2
        ! subgrid scale eddy viscosity
        ! nu_kgs = C_k dx sqrt(k_sgs) = (C_w * dx)^2 (OP1/OP2)
        ! OP1 = (Sd_ij Sd_ij)^(3/2)
        ! OP2 = (SR_ij SR_ij)^(5/2) + (Sd_ij Sd_ij)^(5/4)
        ! Add small fraction to denominator to avoid division by zero
        OP1 = Sd_sqr**1.5_rk
        OP2 = SR_sqr**2.5_rk + Sd_sqr**1.25_rk + eps

        elemPos = low_bound + iElem
        ! turbulent viscosity
        turbVisc(elemPos) = visc_coeff * (OP1/OP2) / dtL

      end do !iElem
    end do !iChunks

  end subroutine mus_turbVisc_WALE_3D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate eddy viscosity with WALE (Wall-Adapting Local Eddy-viscosity)
  !! model
  !! \todo add reference and formula
  subroutine mus_turbVisc_WALE_2D(turbVisc, turbConfig, gradData, auxField,  &
    &                             velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
    ! --------------------------------------------------------------------------
    !> output: turbulent viscosity
    real(kind=rk), intent(out) :: turbVisc(:)
    !> Contains turbulenct coefficients
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
    !> current level lattice element size
    real(kind=rk), intent(in) :: dxL
    !> current level lattice time step size
    real(kind=rk), intent(in) :: dtL
    !> Object that contains pointers to calculate gradients
    type(mus_Grad_type), intent(in) :: Grad
    ! --------------------------------------------------------------------------
    integer :: iElem
    real(kind=rk) :: gradU(2,2,vlen)
    real(kind=rk) :: gradU_sqr(2,2,vlen)
    real(kind=rk) :: SR(3), Sd(3), onehalf_trSd, Sd_sqr, SR_sqr, OP1, OP2
    integer :: ndims
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    real(kind=rk) :: visc_coeff
    !> gradient of velocity
    ! --------------------------------------------------------------------------
    visc_coeff = (turbConfig%coeff%C_w * dxL )**2
    nDims = 2
    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      gradU(:,:,1:nChunkElems) = Grad%U_ptr( auxField    = auxField,    &
        &                                    gradData    = gradData,    &
        &                                    velPos      = velPos,      &
        &                                    nAuxScalars = nAuxScalars, &
        &                                    nDims       = nDims,       &
        &                                    nSolve      = nChunkElems, &
        &                                    elemOffset  = low_bound    )

      do iElem = 1, nChunkElems
      ! square of velocity gradient. gradU . gradU
        gradU_sqr(:,:,iElem) = matmul(gradU(:,:,iElem), gradU(:,:,iElem))
      end do !iElem

      do iElem = 1, nChunkElems
        ! traceless symmetric part of the square of the velocity gradient tensor
        ! Sd_ij = 1/2(du_k/dx_i du_j/dx_k + du_k/dx_j du_i/dx_k)
        !       - 1/3\delta_ij du_k/dx_l du_l/dx_k
        onehalf_trSd = (gradU_sqr(1,1,iElem) + gradU_sqr(2,2,iElem))*div1_2
        Sd(1) = gradU_sqr(1,1,iElem) - onehalf_trSd            !XX
        Sd(2) = gradU_sqr(2,2,iElem) - onehalf_trSd            !YY
        Sd(3) = 0.5_rk*(gradU_sqr(1,2,iElem)+gradU_sqr(2,1,iElem))    !XY

        ! double inner product of Sd: Sd_ij Sd_ij
        Sd_sqr = Sd(1)**2 + Sd(2)**2 + 2.0_rk*Sd(3)**2

        ! symmetric strain rate tensors
        SR(1) = gradU(1,1,iElem)
        SR(2) = gradU(2,2,iElem)
        SR(3) = (gradU(1,2,iElem)+gradU(2,1,iElem))*0.5_rk

        ! double inner product of tensor
        SR_sqr = SR(1)**2 + SR(2)**2 + 2.0_rk*SR(3)**2

        ! sub-grid scale kinetic energy
        ! k_sgs = (C_w^2 * dx /C_k)^2 (OP1/OP2)^2
        ! subgrid scale eddy viscosity
        ! nu_kgs = C_k dx sqrt(k_sgs) = (C_w * dx)^2 (OP1/OP2)
        ! OP1 = (Sd_ij Sd_ij)^(3/2)
        ! OP2 = (SR_ij SR_ij)^(5/2) + (Sd_ij Sd_ij)^(5/4)
        ! Add small fraction to denominator to avoid division by zero
        OP1 = Sd_sqr**1.5_rk
        OP2 = SR_sqr**2.5_rk + Sd_sqr**1.25_rk + eps

        elemPos = low_bound + iElem
        ! turbulent viscosity
        turbVisc(elemPos) = visc_coeff * (OP1/OP2) / dtL

      end do !iElem
    end do !iChunks

  end subroutine mus_turbVisc_WALE_2D
  ! ************************************************************************** !

end module mus_WALE_module
