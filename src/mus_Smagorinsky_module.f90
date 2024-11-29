! Copyright (c) 2019-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!! Smagorinsky les turbulence model.
!! author: Kannan Masilamani
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011-2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013-2015, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
module mus_Smagorinsky_module
  ! treelm modules
  use env_module,                    only: rk
  use tem_param_module,              only: div1_2, div1_3, div2_3, &
    &                                      cs2, rho0, rho0Inv, two_sqrt2
  use tem_compileconf_module,        only: vlen

  ! musubi modules
  use mus_turbulence_module,         only: mus_turbulence_config_type
  use mus_gradData_module,           only: mus_gradData_type, mus_Grad_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_derivedQuantities_module2, only: secondMom_2D, secondMom_3D, &
    &                                      secondMom_minus_cs2_2D,     &
    &                                      secondMom_minus_cs2_3D,     &
    &                                      getEquilibrium, &
    &                                      getEquilibriumIncomp

  implicit none
  private

  public :: mus_turbVisc_Smagorinsky_fromGradU3D
  public :: mus_turbVisc_Smagorinsky_fromGradU2D
  public :: mus_turbVisc_Smagorinsky_fromGradU3D_incomp
  public :: mus_turbVisc_Smagorinsky_fromGradU2D_incomp
  public :: mus_turbVisc_Smagorinsky_fromPreColPDF_2D
  public :: mus_turbVisc_Smagorinsky_fromPreColPDF_3D
  public :: mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_2D
  public :: mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_3D

contains

  ! ************************************************************************** !
  !> Calculate eddy viscosity with smagorinsky model for compressible model
  !! using gradient of velocity
  !! Reference paper:
  !! https://link.springer.com/content/pdf/10.1007/s10494-012-9405-0.pdf?pdf=button
  !! The formula is taken from https://caefn.com/openfoam/smagorinsky-sgs-model
  !! nu_t = C_k delta sqrt(k_sgs)
  !! k_sgs = ((-b+sqrt(b^2+4ac))/2a)^2
  !! a = C_e/delta, b=2/3 tr(dev(Strain)), c = 2 C_k delta (dev(Strain):Strain)
  subroutine mus_turbVisc_Smagorinsky_fromGradU3D(turbVisc, turbConfig, &
    & gradData, auxField, velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
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
    real(kind=rk) :: SR(6), devSR(3), oneThird_trSR, devSR_SR, tr_SR
    real(kind=rk) :: sqrt_k_sgs, visc_coeff
    real(kind=rk) :: a, b, c
    !> gradient of velocity
    real(kind=rk) :: gradU(3,3,vlen)
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! --------------------------------------------------------------------------
    ! viscosity coeff
    visc_coeff = turbConfig%coeff%C_k * dxL / dtL

    ! gradU = getGradU(iElem, auxField, gradData, velPos, nAuxScalars, 3)
    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      gradU(:,:,1:nChunkElems) = Grad%U_ptr( auxField    = auxField,    &
        &                                    gradData    = gradData,    &
        &                                    velPos      = velPos,      &
        &                                    nAuxScalars = nAuxScalars, &
        &                                    nDims       = 3,           &
        &                                    nSolve      = nChunkElems, &
        &                                    elemOffset  = low_bound    )

      do iElem = 1, nChunkElems
        ! symmetric strain rate tensors
        SR(1) = gradU(1,1,iElem)                            !S_XX
        SR(2) = gradU(2,2,iElem)                            !S_YY
        SR(3) = gradU(3,3,iElem)                            !S_ZZ
        SR(4) = (gradU(1,2,iElem)+gradU(2,1,iElem))*0.5_rk  !S_XY
        SR(5) = (gradU(2,3,iElem)+gradU(3,2,iElem))*0.5_rk  !S_YZ
        SR(6) = (gradU(1,3,iElem)+gradU(3,1,iElem))*0.5_rk  !S_XZ

        ! trace of strain rate
        tr_SR = SR(1) + SR(2) + SR(3)
        ! onethird of trace of strain rate
        oneThird_trSR = tr_SR*div1_3
        ! deviatoric strain rate differs from strain rate only in diagonal terms
        devSR(1) = SR(1) - oneThird_trSR
        devSR(2) = SR(2) - oneThird_trSR
        devSR(3) = SR(3) - oneThird_trSR

        ! inner product of devSR:SR
        devSR_SR = SR(1)*devSR(1) + SR(2)*devSR(2) + SR(3)*devSR(3) &
         &      + 2.0_rk*( SR(4)**2 + SR(5)**2 + SR(6)**2)

        ! parameters to compute subgrid scale kinetic energy
        a = turbConfig%coeff%C_e / dxL
        b = tr_SR * div2_3
        c = 2.0_rk * turbConfig%coeff%C_k * dxL * devSR_SR

        ! subgrid scale kinetic energy
        sqrt_k_sgs = (-b + sqrt(b**2 + 4.0_rk*a*c))/(2.0_rk*a)

        elemPos = low_bound + iElem
        ! subgrid scale turbulent viscosity normalized to current level
        turbVisc(elemPos) = visc_coeff * sqrt_k_sgs
      end do!iELem
    end do!iChunks

  end subroutine mus_turbVisc_Smagorinsky_fromGradU3D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate eddy viscosity with smagorinsky model for compressible model
  !! using gradient of velocity for 2D layout
  subroutine mus_turbVisc_Smagorinsky_fromGradU2D(turbVisc, turbConfig, &
    & gradData, auxField, velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
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
    real(kind=rk) :: SR(3), devSR(2), oneHalf_trSR, devSR_SR, tr_SR
    real(kind=rk) :: sqrt_k_sgs, visc_coeff
    real(kind=rk) :: a, b, c
    !> gradient of velocity
    real(kind=rk) :: gradU(2,2,vlen)
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! --------------------------------------------------------------------------
    ! viscosity coeff
    visc_coeff = turbConfig%coeff%C_k * dxL / dtL

    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      gradU(:,:,1:nChunkElems) = Grad%U_ptr( auxField    = auxField,    &
        &                                    gradData    = gradData,    &
        &                                    velPos      = velPos,      &
        &                                    nAuxScalars = nAuxScalars, &
        &                                    nDims       = 2,           &
        &                                    nSolve      = nChunkElems, &
        &                                    elemOffset  = low_bound    )

      do iElem = 1, nChunkElems

        ! symmetric strain rate tensors
        SR(1) = gradU(1,1,iElem)                            !S_XX
        SR(2) = gradU(2,2,iElem)                            !S_YY
        SR(3) = (gradU(1,2,iElem)+gradU(2,1,iElem))*0.5_rk  !S_XY

        ! trace of strain rate
        tr_SR = (SR(1) + SR(2))
        ! one half of trace of strain rate
        oneHalf_trSR = tr_SR*div1_2
        ! deviatoric strain rate differs from strain rate only in diagonal terms
        devSR(1) = SR(1) - oneHalf_trSR
        devSR(2) = SR(2) - oneHalf_trSR

        ! inner product of devSR:SR
        devSR_SR = SR(1)*devSR(1) + SR(2)*devSR(2) &
          &      + 2.0_rk*( SR(3)**2 )

        ! parameters to compute subgrid scale kinetic energy
        a = turbConfig%coeff%C_e / dxL
        b = tr_SR
        c = 2.0_rk * turbConfig%coeff%C_k * dxL * devSR_SR

        ! subgrid scale kinetic energy
        sqrt_k_sgs = (-b + sqrt(b**2 + 4.0_rk*a*c))/(2.0_rk*a)

        elemPos = low_bound + iElem
        ! subgrid scale turbulent viscosity
        turbVisc(elemPos) = visc_coeff * sqrt_k_sgs
      end do
    end do

  end subroutine mus_turbVisc_Smagorinsky_fromGradU2D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculate eddy viscosity with smagorinsky model for incompressible model
  !! using gradient of velocity
  subroutine mus_turbVisc_Smagorinsky_fromGradU3D_incomp(turbVisc, turbConfig, &
    & gradData, auxField, velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
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
    real(kind=rk) :: SR(6), SR_mag_sqr
    real(kind=rk) :: k_sgs, visc_coeff
    !> gradient of velocity
    real(kind=rk) :: gradU(3,3,vlen)
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! --------------------------------------------------------------------------
    ! viscosity coeff
    visc_coeff = turbConfig%coeff%C_k * dxL / dtL

    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      gradU(:,:,1:nChunkElems) = Grad%U_ptr( auxField    = auxField,    &
        &                                    gradData    = gradData,    &
        &                                    velPos      = velPos,      &
        &                                    nAuxScalars = nAuxScalars, &
        &                                    nDims       = 3,           &
        &                                    nSolve      = nChunkElems, &
        &                                    elemOffset  = low_bound    )

      do iElem = 1, nChunkElems

        ! symmetric strain rate tensors
        SR(1) = gradU(1,1,iElem)                            !S_XX
        SR(2) = gradU(2,2,iElem)                            !S_YY
        SR(3) = gradU(3,3,iElem)                            !S_ZZ
        SR(4) = (gradU(1,2,iElem)+gradU(2,1,iElem))*0.5_rk  !S_XY
        SR(5) = (gradU(2,3,iElem)+gradU(3,2,iElem))*0.5_rk  !S_YZ
        SR(6) = (gradU(1,3,iElem)+gradU(3,1,iElem))*0.5_rk  !S_XZ

        ! magnitude of strain rate tensor, sqrt(2 S:S)
        SR_mag_sqr = 2.0_rk*( SR(1)**2 + SR(2)**2 + SR(3)**2 &
          &                  + 2.0_rk*(SR(4)**2 + SR(5)**2 + SR(6)**2))

        ! subgrid scale kinetic energy of incompressible model
        k_sgs = turbConfig%coeff%C_k * dxL**2 * SR_mag_sqr / turbConfig%coeff%C_e

        elemPos = low_bound + iElem
        ! subgrid scale turbulent viscosity
        turbVisc(elemPos) = visc_coeff * sqrt(k_sgs)
      end do
    end do

  end subroutine mus_turbVisc_Smagorinsky_fromGradU3D_incomp
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate eddy viscosity with smagorinsky model for incompressible model
  !! using gradient of velocity for 2D layout
  subroutine mus_turbVisc_Smagorinsky_fromGradU2D_incomp(turbVisc, turbConfig, &
    & gradData, auxField, velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
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
    real(kind=rk) :: SR(3), SR_mag_sqr
    real(kind=rk) :: k_sgs, visc_coeff
    !> gradient of velocity
    real(kind=rk) :: gradU(2,2,vlen)
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! --------------------------------------------------------------------------
    ! viscosity coeff
    visc_coeff = turbConfig%coeff%C_k * dxL / dtL

    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      gradU(:,:,1:nChunkElems) = Grad%U_ptr( auxField    = auxField,    &
        &                                    gradData    = gradData,    &
        &                                    velPos      = velPos,      &
        &                                    nAuxScalars = nAuxScalars, &
        &                                    nDims       = 2,           &
        &                                    nSolve      = nChunkElems, &
        &                                    elemOffset  = low_bound    )

      do iElem = 1, nChunkElems

        ! symmetric strain rate tensors
        SR(1) = gradU(1,1,iElem)                            !S_XX
        SR(2) = gradU(2,2,iElem)                            !S_YY
        SR(3) = (gradU(1,2,iElem)+gradU(2,1,iElem))*0.5_rk  !S_XY

        ! magnitude of strain rate tensor, sqrt(2 S:S)
        SR_mag_sqr = 2.0_rk*( SR(1)**2 + SR(2)**2 + 2.0_rk*(SR(3)**2) )

        ! subgrid scale kinetic energy of incompressible model
        k_sgs = turbConfig%coeff%C_k * dxL**2 * SR_mag_sqr / turbConfig%coeff%C_e

        elemPos = low_bound + iElem
        ! subgrid scale turbulent viscosity
        turbVisc(elemPos) = visc_coeff * sqrt(k_sgs)
      end do
    end do

  end subroutine mus_turbVisc_Smagorinsky_fromGradU2D_incomp
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculate eddy viscosity with smagorinsky model for compressible model
  !! using pre-collision PDF.
  !! Schneider, A. (2015). A Consistent Large Eddy Approach for Lattice
  !! Boltzmann Methods and its Application to Complex Flows.
  !! Technical University Kaiserslautern.
  subroutine mus_turbVisc_Smagorinsky_fromPreColPDF_3D(turbVisc, turbConfig, &
    & state, neigh, auxField, densPos, velPos, nSize, nSolve, nScalars,      &
    & nAuxScalars, layout, dxL, dtL, viscKine)
    ! --------------------------------------------------------------------------
    !> output: turbulent viscosity
    real(kind=rk), intent(out) :: turbVisc(:)
    !> Contains turbulenct coefficients
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
    !> Background kinematic viscosity in lattice divided by dtL
    real(kind=rk), intent(in) :: viscKine(:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: visc_coeff, rho, inv_rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(6), nEqTensMag, viscKineTerm, nEqTensTerm
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    ! viscosity coeff
    visc_coeff = (turbConfig%coeff%C_s * dxL )**2

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
       f_preCol(iDir) = state (                               &
         &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      inv_rho = 1.0_rk/rho
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = auxField( elemOff + velPos(3) )

      ! Calculate the equilibrium distribution function
      fEq(:) = getEquilibrium( rho, vel, layout)

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric deviatoric second-order tensor of
      ! nonEquilibrium part
      ! the static part cs2 I is usually neglected for weakly compressible flows
      ! however, in current implementation it is considered
      nEqTens = secondMom_minus_cs2_3D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      ! magnitude of second-order tensor
      nEqTensMag = sqrt( nEqTens(1)**2 + nEqTens(2)**2 + nEqTens(3)**2    &
        &        + 2.0_rk*(nEqTens(4)**2 + nEqTens(5)**2 + nEqTens(6)**2) )

      ! turbulent viscosity
      !nu_t = (sqrt((v+cs2 dt/2)^2+2(Cs dx)^2|nEQ|/rho) - (v+cs2 dt/2))/2
      ! viscKine is scaled to current level dt so multiply viscKine with dt
      ! to get unit consistent
      viscKineTerm = (viscKine(iElem) + div1_2 * cs2) * dtL
      nEqTensTerm = two_sqrt2 * visc_coeff * nEqTensMag * inv_rho
      turbVisc(iElem) =  ( sqrt(viscKineTerm**2 + nEqTensTerm) &
        &             - viscKineTerm ) * div1_2 / dtL
    end do

  end subroutine mus_turbVisc_Smagorinsky_fromPreColPDF_3D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculate eddy viscosity with smagorinsky model for compressible model
  !! using pre-collision PDF.
  !! Schneider, A. (2015). A Consistent Large Eddy Approach for Lattice
  !! Boltzmann Methods and its Application to Complex Flows.
  !! Technical University Kaiserslautern.
  subroutine mus_turbVisc_Smagorinsky_fromPreColPDF_2D(turbVisc, turbConfig, &
    & state, neigh, auxField, densPos, velPos, nSize, nSolve, nScalars,      &
    & nAuxScalars, layout, dxL, dtL, viscKine)
    ! --------------------------------------------------------------------------
    !> output: turbulent viscosity
    real(kind=rk), intent(out) :: turbVisc(:)
    !> Contains turbulenct coefficients
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
    !> Background kinematic viscosity in lattice divided by dtL
    real(kind=rk), intent(in) :: viscKine(:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: visc_coeff, rho, inv_rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(3), nEqTensMag, viscKineTerm, nEqTensTerm
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    ! viscosity coeff
    visc_coeff = (turbConfig%coeff%C_s * dxL )**2

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
       f_preCol(iDir) = state (                               &
         &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      inv_rho = 1.0_rk/rho
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = 0._rk

      ! Calculate the equilibrium distribution function
      fEq(:) = getEquilibrium( rho, vel, layout)

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric deviatoric second-order tensor of
      ! nonEquilibrium part
      ! the static part cs2 I is usually neglected for weakly compressible flows
      ! however, in current implementation it is considered
      nEqTens = secondMom_minus_cs2_2D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      ! magnitude of second-order tensor
      nEqTensMag = sqrt( nEqTens(1)**2 + nEqTens(2)**2 + 2.0_rk*nEqTens(3)**2 )

      ! turbulent viscosity
      !nu_t = (sqrt((v+cs2 dt/2)^2+2(Cs dx)^2|nEQ|/rho) - (v+cs2 dt/2))/2
      ! viscKine is scaled to current level dt so multiply viscKine with dt
      ! to get unit consistent
      viscKineTerm = (viscKine(iElem) + div1_2 * cs2) * dtL
      nEqTensTerm = 2.0_rk * visc_coeff * nEqTensMag * inv_rho
      turbVisc(iElem) =  ( sqrt(viscKineTerm**2 + nEqTensTerm) &
        &             - viscKineTerm ) * div1_2 / dtL
    end do

  end subroutine mus_turbVisc_Smagorinsky_fromPreColPDF_2D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculate eddy viscosity with smagorinsky model for incompressible model
  !! using pre-collision PDF
  subroutine mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_3D(turbVisc,   &
    & turbConfig, state, neigh, auxField, densPos, velPos, nSize, nSolve, &
    & nScalars, nAuxScalars, layout, dxL, dtL, viscKine)
    ! --------------------------------------------------------------------------
    !> output: turbulent viscosity
    real(kind=rk), intent(out) :: turbVisc(:)
    !> Contains turbulenct coefficients
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
    !> Background kinematic viscosity in lattice divided by dtL
    real(kind=rk), intent(in) :: viscKine(:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: visc_coeff, rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(6), nEqTensMag, viscKineTerm, nEqTensTerm
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    ! viscosity coeff
    visc_coeff = (turbConfig%coeff%C_s * dxL )**2

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
       f_preCol(iDir) = state (                               &
         &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = auxField( elemOff + velPos(3) )

      ! Calculate the equilibrium distribution function
      fEq(:) = getEquilibriumIncomp( rho, vel, layout, rho0)

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric second-order tensor of nonEquilibrium part
      nEqTens = secondMom_3D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      ! magnitude of second-order tensor
      nEqTensMag = sqrt( nEqTens(1)**2 + nEqTens(2)**2 + nEqTens(3)**2      &
        &          + 2.0_rk*(nEqTens(4)**2 + nEqTens(5)**2 + nEqTens(6)**2) )

      ! turbulent viscosity
      !nu_t = (sqrt((v+cs2 dt/2)^2+2(Cs dx)^2|nEQ|/rho) - (v+cs2 dt/2))/2
      ! viscKine is scaled to current level dt so multiply viscKine with dt
      ! to get unit consitent
      viscKineTerm = (viscKine(iElem) + div1_2 * cs2) * dtL
      ! use reference density for incompressible model
      nEqTensTerm = two_sqrt2 * visc_coeff * nEqTensMag * rho0Inv
      turbVisc(iElem) = ( sqrt(viscKineTerm**2 + nEqTensTerm) &
        &             - viscKineTerm ) * div1_2 / dtL
    end do

  end subroutine mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_3D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculate eddy viscosity with smagorinsky model for incompressible model
  !! using pre-collision PDF
  subroutine mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_2D(turbVisc,   &
    & turbConfig, state, neigh, auxField, densPos, velPos, nSize, nSolve, &
    & nScalars, nAuxScalars, layout, dxL, dtL, viscKine)
    ! --------------------------------------------------------------------------
    !> output: turbulent viscosity
    real(kind=rk), intent(out) :: turbVisc(:)
    !> Contains turbulenct coefficients
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
    !> Background kinematic viscosity in lattice divided by dtL
    real(kind=rk), intent(in) :: viscKine(:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: visc_coeff, rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(3), nEqTensMag, viscKineTerm, nEqTensTerm
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    ! viscosity coeff
    visc_coeff = (turbConfig%coeff%C_s * dxL )**2

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
       f_preCol(iDir) = state (                               &
         &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = 0._rk

      ! Calculate the equilibrium distribution function
      fEq(:) = getEquilibriumIncomp( rho, vel, layout, rho0)

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric second-order tensor of nonEquilibrium part
      nEqTens = secondMom_2D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      ! magnitude of second-order tensor
      nEqTensMag = sqrt( nEqTens(1)**2 + nEqTens(2)**2 + 2.0_rk*nEqTens(3)**2 )

      ! turbulent viscosity
      !nu_t = (sqrt((v+cs2 dt/2)^2+2(Cs dx)^2|nEQ|/rho) - (v+cs2 dt/2))/2
      ! viscKine is scaled to current level dt so multiply viscKine with dt
      ! to get unit consitent
      viscKineTerm = (viscKine(iElem) + div1_2 * cs2) * dtL
      ! use reference density for incompressible model
      nEqTensTerm = 2.0_rk * visc_coeff * nEqTensMag * rho0Inv
      turbVisc(iElem) = ( sqrt(viscKineTerm**2 + nEqTensTerm) &
        &             - viscKineTerm ) * div1_2 / dtL
    end do

  end subroutine mus_turbVisc_Smagorinsky_fromPreColPDF_incomp_2D
  ! ************************************************************************** !

end module mus_Smagorinsky_module
