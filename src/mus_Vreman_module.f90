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
!> This module contains function to compute eddy viscosity using
!! Vreman LES turbulence model.
!! Vreman, A. W. (2004). An eddy-viscosity subgrid-scale model for turbulent
!! shear flow: Algebraic theory and applications. Physics of Fluids, 16(10),
!! 3670–3681.
!! model.
!! author: Kannan Masilamani
module mus_Vreman_module
  ! treelm modules
  use env_module,                    only: rk
  use tem_float_module,              only: operator(.flt.)
  use tem_compileconf_module,        only: vlen

  ! musubi modules
  use mus_turbulence_module,         only: mus_turbulence_config_type
  use mus_gradData_module,           only: mus_gradData_type, mus_Grad_type

  implicit none
  private

  public :: mus_turbVisc_Vreman_3D
  public :: mus_turbVisc_Vreman_2D

contains

  ! ************************************************************************** !
  !> Calculate eddy viscosity with Vreman model for 3D stencil
  !! Fortran implementation of this model:
  !! http://www.vremanresearch.nl/Vreman_Subgridmodel_Fortran.txt
  !!
  !! $$\nu_{turb}=c_v (\Delta_x)^2 \cdot
  !!            \left( \sqrt{\frac{B_\beta}{\alpha_{ij}\alpha{ij}}} \right)$$
  !! with
  !! $$ c_v = 2.5*C_s^2$$, $$C_s$$ - Smagorinsky constant,
  !!
  !! $$ B_\beta = \beta_{11}\beta_{22} - \beta^2_{12} + \beta_{11}\beta_{33}
  !!            - \beta^2_{13} + \beta_{22}\beta_{33} - \beta^2_{23} $$,
  !!
  !! $$ \beta_{ij} = \alpha_{mi}\alpha_{mj}$$,
  !!
  !! $$ \alpha_{ij} = \frac{\partial \Bar{u}_j}{\partial x_i} $$.
  !! $$\alpha_{ij}$$ - Resolved velocity gradient.
  subroutine mus_turbVisc_Vreman_3D(turbVisc, turbConfig, gradData, auxField,  &
    &                               velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
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
    real(kind=rk) :: b_11, b_12, b_13, b_22, b_23, b_33
    real(kind=rk) :: visc_coeff
    real(kind=rk) :: bbeta
    integer :: iElem, ndims
    !> gradient of velocity
    real(kind=rk) :: gradU(3,3,vlen)
    real(kind=rk) :: gradU_sqr
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! --------------------------------------------------------------------------
    ! viscosity coeff
    visc_coeff = turbConfig%coeff%C_v * dxL**2

! TODO : AR vectorize
!      gradU = getGradU(iElem, auxField, gradData, velPos, nAuxScalars, 3)
    ndims = 3
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
        ! beta_ij = gradU_ki gradU_kj
        ! beta = gradU^T gradU
        b_11 = gradU(1,1,iElem)*gradU(1,1,iElem)             &
          &              + gradU(2,1,iElem)*gradU(2,1,iElem) &
          &              + gradU(3,1,iElem)*gradU(3,1,iElem)
        b_12 = gradU(1,1,iElem)*gradU(1,2,iElem)             &
          &              + gradU(2,1,iElem)*gradU(2,2,iElem) &
          &              + gradU(3,1,iElem)*gradU(3,2,iElem)
        b_13 = gradU(1,1,iElem)*gradU(1,3,iElem)             &
          &              + gradU(2,1,iElem)*gradU(2,3,iElem) &
          &              + gradU(3,1,iElem)*gradU(3,3,iElem)
        b_22 = gradU(1,2,iElem)*gradU(1,2,iElem)             &
          &              + gradU(2,2,iElem)*gradU(2,2,iElem) &
          &              + gradU(3,2,iElem)*gradU(3,2,iElem)
        b_23 = gradU(1,2,iElem)*gradU(1,3,iElem)             &
          &              + gradU(2,2,iElem)*gradU(2,3,iElem) &
          &              + gradU(3,2,iElem)*gradU(3,3,iElem)
        b_33 = gradU(1,3,iElem)*gradU(1,3,iElem)             &
          &              + gradU(2,3,iElem)*gradU(2,3,iElem) &
          &              + gradU(3,3,iElem)*gradU(3,3,iElem)
        ! double inner product of gradU
        gradU_sqr = gradU(1,1,iElem)**2 + gradU(1,2,iElem)**2 &
          &       + gradU(1,3,iElem)**2 + gradU(2,1,iElem)**2 &
          &       + gradU(2,2,iElem)**2 + gradU(2,3,iElem)**2 &
          &       + gradU(3,1,iElem)**2 + gradU(3,2,iElem)**2 &
          &       + gradU(3,3,iElem)**2


        ! numerator Bbeta
        bbeta=b_11*b_22-(b_12**2)+b_11*b_33-(b_13**2)+b_22*b_33-(b_23**2)

        elemPos = low_bound + iElem
        if (bbeta .flt. 1e-12_rk) then
          turbVisc(elemPos) = 0.0_rk
        else
          turbVisc(elemPos) = visc_coeff*sqrt(bbeta/gradU_sqr) / dtL
        end if
      end do !iElem
    end do !iChunks

 end subroutine mus_turbVisc_Vreman_3D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate eddy viscosity with Vreman model for 2D stencil
  !! model
  !! \todo add reference and formula
  subroutine mus_turbVisc_Vreman_2D(turbVisc, turbConfig, gradData, auxField,  &
    &                               velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
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
    integer :: iElem, ndims
    real(kind=rk) :: b_11, b_12, b_22
    real(kind=rk) :: gradU_sqr, visc_coeff
    real(kind=rk) :: bbeta
    !> gradient of velocity
    real(kind=rk) :: gradU(2,2,vlen)
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! --------------------------------------------------------------------------
    ! viscosity coeff
    visc_coeff = turbConfig%coeff%C_v * dxL**2

    !  gradU = getGradU(iElem, auxField, gradData, velPos, nAuxScalars, 2)
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

      ! beta_ij = gradU_ki gradU_kj
      ! beta = gradU^T gradU
      do iElem = 1, nChunkElems
        b_11 = gradU(1,1,iElem)*gradU(1,1,iElem) &
        &    + gradU(2,1,iElem)*gradU(2,1,iElem)
        b_12 = gradU(1,1,iElem)*gradU(1,2,iElem) &
        &    + gradU(2,1,iElem)*gradU(2,2,iElem)
        b_22 = gradU(1,2,iElem)*gradU(1,2,iElem) &
        &    + gradU(2,2,iElem)*gradU(2,2,iElem)

        ! double inner product of gradU
        gradU_sqr = gradU(1,1,iElem)**2 + gradU(1,2,iElem)**2 &
         &       + gradU(2,1,iElem)**2 + gradU(2,2,iElem)**2

        ! numerator Bbeta
        bbeta=b_11*b_22-(b_12**2)

        elemPos = low_bound + iElem
        if (bbeta .flt. 1e-12_rk) then
          turbVisc(elemPos) = 0.0_rk
        else
          turbVisc(elemPos) = visc_coeff*sqrt(bbeta/gradU_sqr) / dtL
        end if
      end do
    end do

  end subroutine mus_turbVisc_Vreman_2D
  ! ************************************************************************** !

end module mus_Vreman_module
