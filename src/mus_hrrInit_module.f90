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
! **************************************************************************** !
!> author: Gregorio Gerardo Spinelli
!! This module provides the definitions of M and Minv for
!! MRT advection relaxation scheme for all stencils.
!!
!! The weighted MRT (D3Q27) is based on the following paper
!! Abbas Fakhari, Diogo Bolster, Li-Shi Luo
!! "A weighted multiple-relaxation-time lattice Boltzmann method for multiphase
!! flows and its application to partial coalescence cascades"
!! Journal of Computational Physics, 2017
!!
!! The MRT (D3Q19) implementation here is taken from:\n
!! J. Toelke, S. Freudiger, and M. Krafczyk,
!! "An adaptive scheme using hierarchical grids for lattice Boltzmann
!! multi-phase flow simulations," Comput. Fluids, vol. 35, pp. 820–830,
!! 2006. \n
module mus_hrrInit_module

  ! include treelm modules
  use env_module,               only: rk
  use tem_param_module,         only: div1_6, div1_3, div1_2, cs4inv, cs2
  use tem_debug_module,         only: dbgUnit
  use tem_aux_module,           only: tem_abort

  use mus_scheme_layout_module, only:mus_scheme_layout_type

  implicit none

  private

  public :: HRR_Correction_d2q9
  public :: HRR_Correction_d3q19
  public :: HRR_Correction_d3q27
  public :: getHermitepolynomials
  public :: getHermitepolynomials_D3Q19

  contains

  ! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
  pure subroutine HRR_Correction_d2q9 ( QQ, weight, gradRHOU3, phi, &
    &                                   dens, vel )
  ! -------------------------------------------------------------------- !
    !> stencil size
    integer, intent(in)        :: QQ
    !> weights of the stencil
    real(kind=rk), intent(in)  :: weight(:)
    !> gradient rho V^3
    real(kind=rk), intent(in)  :: gradRHOU3(:)
    !> correction term phi
    real(kind=rk), intent(out) :: phi(:)
    !> correction term phi, rho, vel
    real(kind=rk), intent(out) :: dens, vel(:)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: phi_temp
    ! -------------------------------------------------------------------- !
    dens = 0._rk
    vel(:) = 0._rk

    !iDir = 1
    phi_temp = div1_6 * cs4inv * ( -2._rk * gradRHOU3(1) &
      &                            + gradRHOU3(2)        )
    phi(1) = weight(1) * phi_temp
    dens = dens + phi(1)
    vel(1) = vel(1) - phi(1)
    !vel(2) = vel(2) + 0._rk

    !iDir = 3
    phi(3) = weight(3) * phi_temp
    dens = dens + phi(3)
    vel(1) = vel(1) + phi(3)
    !vel(2) = vel(2) + 0._rk

    !iDir = 2
    phi_temp = div1_6 * cs4inv * ( gradRHOU3(1)           &
      &                            - 2._rk * gradRHOU3(2) )
    phi(2) = weight(2) * phi_temp
    dens = dens + phi(2)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) - phi(2)

    !iDir = 4
    phi(4) = weight(4) * phi_temp
    dens = dens + phi(4)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) + phi(4)

    !iDir = 5
    phi_temp = div1_3 * cs4inv * ( - gradRHOU3(1) &
      &                            - gradRHOU3(2) )
    phi(5) = weight(5) * phi_temp
    dens = dens + phi(5)
    vel(1) = vel(1) - phi(5)
    vel(2) = vel(2) - phi(5)

    !iDir = 8
    phi(8) = weight(8) * phi_temp
    dens = dens + phi(8)
    vel(1) = vel(1) + phi(8)
    vel(2) = vel(2) + phi(8)

    !iDir = 6
    phi(6) = weight(6) * phi_temp
    dens = dens + phi(6)
    vel(1) = vel(1) - phi(6)
    vel(2) = vel(2) + phi(6)

    !iDir = 7
    phi(7) = weight(7) * phi_temp
    dens = dens + phi(7)
    vel(1) = vel(1) + phi(7)
    vel(2) = vel(2) - phi(7)

    !iDir = 9
    phi_temp = -div1_2 * phi_temp
    phi(9) = weight(9) * phi_temp
    dens = dens + phi(9)
    !vel(1) = vel(1) + 0._rk
    !vel(2) = vel(2) + 0._rk

  end subroutine HRR_Correction_d2q9
  ! ****************************************************************************** !

  ! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
  pure subroutine HRR_Correction_d3q19 ( QQ, weight, gradRHOU3, &
    &                                gradRHOUVZ, phi, dens, vel )
    ! -------------------------------------------------------------------- !
    !> stencil size
    integer, intent(in)        :: QQ
    !> weights of the stencil
    real(kind=rk), intent(in)  :: weight(:)
    !> gradient rho u^3, rho u v w
    real(kind=rk), intent(in)  :: gradRHOU3(:)
    real(kind=rk), intent(in)  :: gradRHOUVZ(:)
    !> correction term phi, rho, vel
    real(kind=rk), intent(out) :: phi(:), dens, vel(:)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: phi_temp
    ! -------------------------------------------------------------------- !
    dens = 0._rk
    vel(:) = 0._rk

    !iDir = 1
    phi_temp = div1_6 * cs4inv * ( -2._rk * gradRHOU3(1)         &
      &                            + gradRHOU3(2) + gradRHOU3(3) )
    phi(1) = weight(1) * phi_temp
    dens = dens + phi(1)
    vel(1) = vel(1) - phi(1)
    !vel(2) = vel(2) + 0._rk
    !vel(3) = vel(3) + 0._rk

    !iDir = 4
    phi(4) = weight(4) * phi_temp
    dens = dens + phi(4)
    vel(1) = vel(1) + phi(4)
    !vel(2) = vel(2) + 0._rk
    !vel(3) = vel(3) + 0._rk

    !iDir = 2
    phi_temp = div1_6 * cs4inv * ( -2._rk * gradRHOU3(2)         &
      &                            + gradRHOU3(1) + gradRHOU3(3) )
    phi(2) = weight(2) * phi_temp
    dens = dens + phi(2)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) - phi(2)
    !vel(3) = vel(3) + 0._rk

    !iDir = 5
    phi(5) = weight(5) * phi_temp
    dens = dens + phi(5)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) + phi(5)
    !vel(3) = vel(3) + 0._rk

    !iDir = 3
    phi_temp = div1_6 * cs4inv * ( -2._rk * gradRHOU3(3)         &
      &                            + gradRHOU3(1) + gradRHOU3(2) )
    phi(3) = weight(3) * phi_temp
    dens = dens + phi(3)
    !vel(1) = vel(1) + 0._rk
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) - phi(3)

    !iDir = 6
    phi(6) = weight(6) * phi_temp
    dens = dens + phi(6)
    !vel(1) = vel(1) + 0._rk
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) + phi(6)

    !iDir = 7
    phi_temp = phi_temp + div1_2 * cs4inv * ( -gradRHOU3(2) - gradRHOUVZ (1) )
    phi(7) = weight(7) * phi_temp
    dens = dens + phi(7)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) - phi(7)
    vel(3) = vel(3) - phi(7)

    !iDir = 10
    phi(10) = weight(10) * phi_temp
    dens = dens + phi(10)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) + phi(10)
    vel(3) = vel(3) + phi(10)

    !iDir = 8
    phi_temp = phi_temp + cs4inv * ( gradRHOUVZ(1) )
    phi(8) = weight(8) * phi_temp
    dens = dens + phi(8)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) - phi(8)
    vel(3) = vel(3) + phi(8)

    !iDir = 9
    phi(9) = weight(9) * phi_temp
    dens = dens + phi(9)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) + phi(9)
    vel(3) = vel(3) - phi(9)

    !iDir = 11
    phi_temp = div1_6 * cs4inv * (- 2._rk * ( gradRHOU3(1) + gradRHOU3(3) ) &
      &                           + gradRHOU3(2) - 3._rk * gradRHOUVZ (2)   )
    phi(11) = weight(11) * phi_temp
    dens = dens + phi(11)
    vel(1) = vel(1) - phi(11)
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) - phi(11)

    !iDir = 14
    phi(14) = weight(14) * phi_temp
    dens = dens + phi(14)
    vel(1) = vel(1) + phi(14)
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) + phi(14)

    !iDir = 12
    phi_temp = phi_temp + cs4inv * ( gradRHOUVZ(2) )
    phi(12) = weight(12) * phi_temp
    dens = dens + phi(12)
    vel(1) = vel(1) + phi(12)
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) - phi(12)

    !iDir = 13
    phi(13) = weight(13) * phi_temp
    dens = dens + phi(13)
    vel(1) = vel(1) - phi(13)
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) + phi(13)

    !iDir = 15
    phi_temp = div1_6 * cs4inv * (-2._rk * ( gradRHOU3(1) + gradRHOU3(2) ) &
      &                           + gradRHOU3(3) - 3._rk * gradRHOUVZ (3)  )
    phi(15) = weight(15) * phi_temp
    dens = dens + phi(15)
    vel(1) = vel(1) - phi(15)
    vel(2) = vel(2) - phi(15)
    !vel(3) = vel(3) + 0._rk

    !iDir = 18
    phi(18) = weight(18) * phi_temp
    dens = dens + phi(18)
    vel(1) = vel(1) + phi(18)
    vel(2) = vel(2) + phi(18)
    !vel(3) = vel(3) + 0._rk

    !iDir = 16
    phi_temp = phi_temp + cs4inv * ( gradRHOUVZ(3) )
    phi(16) = weight(16) * phi_temp
    dens = dens + phi(16)
    vel(1) = vel(1) - phi(16)
    vel(2) = vel(2) + phi(16)
    !vel(3) = vel(3) + 0._rk

    !iDir = 17
    phi(17) = weight(17) * phi_temp
    dens = dens + phi(17)
    vel(1) = vel(1) + phi(17)
    vel(2) = vel(2) - phi(17)
    !vel(3) = vel(3) + 0._rk

    !iDir = 19
    phi_temp = div1_6 * cs4inv * ( gradRHOU3(1) + gradRHOU3(2) + gradRHOU3(3) )
    phi(19) = weight(19) * phi_temp
    dens = dens + phi(19)
    !vel(1) = vel(1) + 0._rk
    !vel(2) = vel(2) + 0._rk
    !vel(3) = vel(3) + 0._rk

  end subroutine HRR_Correction_d3q19
  ! ****************************************************************************** !

  ! ****************************************************************************** !
  ! based on Lattice Boltzmann Method with regularized non-equilibrium distribution
  ! functions, Jonas Latt and Bastien Chopard 2005
  pure subroutine HRR_Correction_d3q27 ( QQ, weight, gradRHOU3, phi, &
    &                                    dens, vel )
  ! -------------------------------------------------------------------- !
    !> stencil size
    integer, intent(in)        :: QQ
    !> weights of the stencil
    real(kind=rk), intent(in)  :: weight(:)
    !> gradient rho V^3
    real(kind=rk), intent(in)  :: gradRHOU3(:)
    !> correction term phi
    real(kind=rk), intent(out) :: phi(:), dens, vel(:)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: phi_temp
    ! -------------------------------------------------------------------- !
    dens = 0._rk
    vel(:) = 0._rk

    !iDir = 1
    phi_temp = div1_6 * cs4inv * ( -2._rk * gradRHOU3(1)         &
      &                            + gradRHOU3(2) + gradRHOU3(3) )
    phi(1) = weight(1) * phi_temp
    dens = dens + phi(1)
    vel(1) = vel(1) - phi(1)
    !vel(2) = vel(2) + 0._rk
    !vel(3) = vel(3) + 0._rk

    !iDir = 4
    phi(4) = weight(4) * phi_temp
    dens = dens + phi(4)
    vel(1) = vel(1) + phi(4)
    !vel(2) = vel(2) + 0._rk
    !vel(3) = vel(3) + 0._rk

    !iDir = 2
    phi_temp = div1_6 * cs4inv * ( -2._rk * gradRHOU3(2)         &
      &                            + gradRHOU3(1) + gradRHOU3(3) )
    phi(2) = weight(2) * phi_temp
    dens = dens + phi(2)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) - phi(2)
    !vel(3) = vel(3) + 0._rk

    !iDir = 5
    phi(5) = weight(5) * phi_temp
    dens = dens + phi(5)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) + phi(5)
    !vel(3) = vel(3) + 0._rk

    !iDir = 3
    phi_temp = div1_6 * cs4inv * ( -2._rk * gradRHOU3(3)         &
      &                            + gradRHOU3(1) + gradRHOU3(2) )
    phi(3) = weight(3) * phi_temp
    dens = dens + phi(3)
    !vel(1) = vel(1) + 0._rk
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) - phi(3)

    !iDir = 6
    phi(6) = weight(6) * phi_temp
    dens = dens + phi(6)
    !vel(1) = vel(1) + 0._rk
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) + phi(6)

    !iDir = 7
    phi_temp = div1_6 * cs4inv * ( gradRHOU3(1)               &
      &                            - 2._rk * ( gradRHOU3(2)   &
      &                                      + gradRHOU3(3) ) )
    phi(7) = weight(7) * phi_temp
    dens = dens + phi(7)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) - phi(7)
    vel(3) = vel(3) - phi(7)

    !iDir = 10
    phi(10) = weight(10) * phi_temp
    dens = dens + phi(10)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) + phi(10)
    vel(3) = vel(3) + phi(10)

    !iDir = 8
    phi(8) = weight(8) * phi_temp
    dens = dens + phi(8)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) - phi(8)
    vel(3) = vel(3) + phi(8)

    !iDir = 9
    phi(9) = weight(9) * phi_temp
    dens = dens + phi(9)
    !vel(1) = vel(1) + 0._rk
    vel(2) = vel(2) + phi(9)
    vel(3) = vel(3) - phi(9)

    !iDir = 11
    phi_temp = div1_6 * cs4inv * ( gradRHOU3(2)               &
      &                            - 2._rk * ( gradRHOU3(1)   &
      &                                      + gradRHOU3(3) ) )
    phi(11) = weight(11) * phi_temp
    dens = dens + phi(11)
    vel(1) = vel(1) - phi(11)
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) - phi(11)

    !iDir = 14
    phi(14) = weight(14) * phi_temp
    dens = dens + phi(14)
    vel(1) = vel(1) + phi(14)
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) + phi(14)

    !iDir = 12
    phi(12) = weight(12) * phi_temp
    dens = dens + phi(12)
    vel(1) = vel(1) + phi(12)
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) - phi(12)

    !iDir = 13
    phi(13) = weight(13) * phi_temp
    dens = dens + phi(13)
    vel(1) = vel(1) - phi(13)
    !vel(2) = vel(2) + 0._rk
    vel(3) = vel(3) + phi(13)

    !iDir = 15
    phi_temp = div1_6 * cs4inv * ( gradRHOU3(3)               &
     &                           - 2._rk * ( gradRHOU3(1)    &
     &                                      + gradRHOU3(2) ) )
    phi(15) = weight(15) * phi_temp
    dens = dens + phi(15)
    vel(1) = vel(1) - phi(15)
    vel(2) = vel(2) - phi(15)
    !vel(3) = vel(3) + 0._rk

    !iDir = 18
    phi(18) = weight(18) * phi_temp
    dens = dens + phi(18)
    vel(1) = vel(1) + phi(18)
    vel(2) = vel(2) + phi(18)
    !vel(3) = vel(3) + 0._rk

    !iDir = 16
    phi(16) = weight(16) * phi_temp
    dens = dens + phi(16)
    vel(1) = vel(1) - phi(16)
    vel(2) = vel(2) + phi(16)
    !vel(3) = vel(3) + 0._rk

    !iDir = 17
    phi(17) = weight(17) * phi_temp
    dens = dens + phi(17)
    vel(1) = vel(1) + phi(17)
    vel(2) = vel(2) - phi(17)
    !vel(3) = vel(3) + 0._rk

    !iDir = 19
    phi_temp = div1_3 * cs4inv * ( - ( gradRHOU3(1) + gradRHOU3(2) &
      &                                 + gradRHOU3(3) )           )
    phi(19) = weight(19) * phi_temp
    dens = dens + phi(19)
    vel(1) = vel(1) - phi(19)
    vel(2) = vel(2) - phi(19)
    vel(3) = vel(3) - phi(19)

    !iDir = 26
    phi(26) = weight(26) * phi_temp
    dens = dens + phi(26)
    vel(1) = vel(1) + phi(26)
    vel(2) = vel(2) + phi(26)
    vel(3) = vel(3) + phi(26)

    !iDir = 22
    phi(22) = weight(22) * phi_temp
    dens = dens + phi(22)
    vel(1) = vel(1) - phi(22)
    vel(2) = vel(2) + phi(22)
    vel(3) = vel(3) + phi(22)

    !iDir = 23
    phi(23) = weight(23) * phi_temp
    dens = dens + phi(23)
    vel(1) = vel(1) + phi(23)
    vel(2) = vel(2) - phi(23)
    vel(3) = vel(3) - phi(23)

    !iDir = 25
    phi(25) = weight(25) * phi_temp
    dens = dens + phi(25)
    vel(1) = vel(1) + phi(25)
    vel(2) = vel(2) + phi(25)
    vel(3) = vel(3) - phi(25)

    !iDir = 20
    phi(20) = weight(20) * phi_temp
    dens = dens + phi(20)
    vel(1) = vel(1) - phi(20)
    vel(2) = vel(2) - phi(20)
    vel(3) = vel(3) + phi(20)

    !iDir = 24
    phi(24) = weight(24) * phi_temp
    dens = dens + phi(24)
    vel(1) = vel(1) + phi(24)
    vel(2) = vel(2) - phi(24)
    vel(3) = vel(3) + phi(24)

    !iDir = 21
    phi(21) = weight(21) * phi_temp
    dens = dens + phi(21)
    vel(1) = vel(1) - phi(21)
    vel(2) = vel(2) + phi(21)
    vel(3) = vel(3) - phi(21)

    !iDir = 27
    phi_temp = -div1_2 * phi_temp
    phi(27) = weight(27) * phi_temp
    dens = dens + phi(27)
    !vel(1) = vel(1) + 0._rk
    !vel(2) = vel(2) + 0._rk
    !vel(3) = vel(3) + 0._rk

end subroutine HRR_Correction_d3q27
! ****************************************************************************** !

! ************************************************************************** !
!> This function computes Hermite polinomial. It gives in output minimum
!  up to 2nd-order polynomials. It is coded for up to 6th-order polynomials
! deprecated, used only to check the optimized R, RR, HRR, PRR, DRT functions
subroutine getHermitepolynomials( nDims, QQ, layout, H_order) !, H )
  ! --------------------------------------------------------------------------
  !> number of physical dimensions
  integer, intent(in) :: nDims
  !> number of stencil streaming directions
  integer, intent(in) :: QQ
  !> current layout
  type(mus_scheme_layout_type), intent(in) :: layout
  !> maximum order of the Hermite polynomials
  integer, intent(in) :: H_order
  !> Hermite polynomials matrix
  !real(kind=rk), intent(inout) :: H(:,:)
  ! --------------------------------------------------------------------------
  integer :: iDir
  real(kind=rk) :: c_x, c_y, c_z
  real(kind=rk) :: H(QQ,QQ)
  ! --------------------------------------------------------------------------

  if (nDims == 3) then
  ! Hermite polynomials: 1=H1x, 2=H1y, 3=H1z, 4=H2xx, 5=H2yy, 6=H2zz, 7=H2xy,
  !                      8=H2xz, 9=H2yz, 10=H3xxy, 11=H3xxz, 12=H3xyy,
  !                      13=H3xzz, 14=H3yzz, 15=H3yyz,
  !    (only for D3Q27)  16=H3xyz, 17=H4xxyy, 18=H4xxzz, 19=H4yyzz, 20=H4xyzz,
  !    (only for D3Q27)  21=H4xyyz, 22=H4xxyz, 23=H5xxyzz, 24=H5xxyyz, 25=H5xyyzz,
  !    (only for D3Q27)  26=H6xxyyzz
  !                      QQ=H0
    do iDir = 1, QQ

      c_x = layout%fStencil%cxDirRK(1,iDir)
      c_y = layout%fStencil%cxDirRK(2,iDir)
      c_z = layout%fStencil%cxDirRK(3,iDir)

      ! 0th order H0 = 1
      H(iDir,QQ) = 1._rk

      ! 1st order H1n = c_n
      H(iDir, 1) = c_x
      H(iDir, 2) = c_y
      H(iDir, 3) = c_z

      ! 2nd order
      H(iDir, 4) = c_x ** 2 - cs2
      H(iDir, 5) = c_y ** 2 - cs2
      H(iDir, 6) = c_z ** 2 - cs2
      H(iDir, 7) = c_x * c_y
      H(iDir, 8) = c_x * c_z
      H(iDir, 9) = c_y * c_z

      ! 3d order
      if (H_order < 3) CYCLE
      ! recursive formula valide from 3rd order on
      ! H (n+m+l) x_n y_m z_l = H(n)x_n H(m)y_m H(l)z_l
      H(iDir, 10) = H(iDir, 4) * c_y ! H1y = c_y
      H(iDir, 11) = H(iDir, 4) * c_z
      H(iDir, 12) = c_x * H(iDir, 5)
      H(iDir, 13) = c_x * H(iDir, 6)
      H(iDir, 14) = c_y * H(iDir, 6)
      H(iDir, 15) = H(iDir, 5) * c_z


      ! Hermite polynomials: 1=H1x, 2=H1y, 3=H1z, 4=H2xx, 5=H2yy, 6=H2zz, 7=H2xy,
      !                      8=H2xz, 9=H2yz, 10=H3xxy, 11=H3xxz, 12=H3xyy, 13=H3xzz,
      !                      14=H3yzz, 15=H3yyz
      write(dbgUnit(1),*) "iDir = ", iDir
      write(dbgUnit(1),*) "  H0 = ", H(iDir,QQ)
      write(dbgUnit(1),*) "  H1x_1 = ", H(iDir,1)
      write(dbgUnit(1),*) "  H1y_2 = ", H(iDir,2)
      write(dbgUnit(1),*) "  H1z_3 = ", H(iDir,3)
      write(dbgUnit(1),*) "  H2xx_4 = ", H(iDir,4)
      write(dbgUnit(1),*) "  H2yy_5 = ", H(iDir,5)
      write(dbgUnit(1),*) "  H2zz_6 = ", H(iDir,6)
      write(dbgUnit(1),*) "  H2xy_7 = ", H(iDir,7)
      write(dbgUnit(1),*) "  H2xz_8 = ", H(iDir,8)
      write(dbgUnit(1),*) "  H2yz_9 = ", H(iDir,9)
      write(dbgUnit(1),*) "  H3xxy_10 = ", H(iDir,10)
      write(dbgUnit(1),*) "  H3xxz_11 = ", H(iDir,11)
      write(dbgUnit(1),*) "  H3xyy_12 = ", H(iDir,12)
      write(dbgUnit(1),*) "  H3xzz_13 = ", H(iDir,13)
      write(dbgUnit(1),*) "  H3yzz_14 = ", H(iDir,14)
      write(dbgUnit(1),*) "  H3yyz_15 = ", H(iDir,15)

      ! only for D3Q27
      ! 16=H3xyz, 17=H4xxyy, 18=H4xxzz, 19=H4yyzz, 20=H4xyzz, 21=H4xyyz,
      ! 22=H4xxyz, 23=H5xxyzz, 24=H5xxyyz, 25=H5xyyzz, 26=H6xxyyzz
      if (QQ == 27) then
        H(iDir, 16) = c_x * c_y * c_z
        write(dbgUnit(1),*) "  H3xyz_16 = ", H(iDir,16)

        ! 4th order
        if (H_order < 4) CYCLE
        H(iDir, 17) = H(iDir, 4) * H(iDir, 5)
        H(iDir, 18) = H(iDir, 4) * H(iDir, 6)
        H(iDir, 19) = H(iDir, 5) * H(iDir, 6)
        H(iDir, 20) = c_x * H(iDir, 14)
        H(iDir, 21) = c_x * H(iDir, 15)
        H(iDir, 22) = H(iDir,10) * c_z
        write(dbgUnit(1),*) "  H4xxyy_17 = ", H(iDir,17)
        write(dbgUnit(1),*) "  H4xxzz_18 = ", H(iDir,18)
        write(dbgUnit(1),*) "  H4yyzz_19 = ", H(iDir,19)
        write(dbgUnit(1),*) "  H4xyzz_20 = ", H(iDir,20)
        write(dbgUnit(1),*) "  H4xyyz_21 = ", H(iDir,21)
        write(dbgUnit(1),*) "  H4xxyz_22 = ", H(iDir,22)

        ! 5th order
        if (H_order < 5) CYCLE
        H(iDir, 23) = H(iDir, 18) * c_y
        H(iDir, 24) = H(iDir, 17) * c_z
        H(iDir, 25) = c_x * H(iDir, 19)
        write(dbgUnit(1),*) "  H5xxyzz_23 = ", H(iDir,23)
        write(dbgUnit(1),*) "  H5xxyyz_24 = ", H(iDir,24)
        write(dbgUnit(1),*) "  H5xyyzz_25 = ", H(iDir,25)

        ! 6th order
        if (H_order < 6) CYCLE
        H(iDir, 26) = H(iDir, 17) * H(iDir, 6)
        write(dbgUnit(1),*) "  H6xxyyzz_26 = ", H(iDir,26)
      endif
    enddo
    flush(dbgUnit(1))
    call tem_abort('done')
  else !nDims == 2
  ! Hermite polynomials: 1=H1x, 2=H1y, 3=H2xx, 4=H2yy, 5=H2xy,
  !                      6=H3xxy, 7=H3xyy, 8=H4xxyy, QQ=H0
    do iDir = 1, QQ

      c_x = layout%fStencil%cxDirRK(1,iDir)
      c_y = layout%fStencil%cxDirRK(2,iDir)

      ! 0th order H0 = 1
      H(iDir,QQ) = 1._rk

      ! 1st order H1n = c_n
      H(iDir, 1) = c_x
      H(iDir, 2) = c_y

      ! 2nd order
      H(iDir, 3) = c_x ** 2 - cs2
      H(iDir, 4) = c_y ** 2 - cs2
      H(iDir, 5) = c_x * c_y

      ! 3d order
      if (H_order < 3) CYCLE
      ! recursive formula valide from 3rd order on
      ! H (n+m+l) x_n y_m z_l = H(n)x_n H(m)y_m H(l)z_l
      H(iDir, 6) = H(iDir, 3) * c_y ! H1y = c_y
      H(iDir, 7) = c_x * H(iDir, 4)

      ! 4th order
      if (H_order < 4) CYCLE
      H(iDir, 8) = H(iDir, 3) * H(iDir, 4)

    enddo

  ! Hermite polynomials: 1=H1x, 2=H1y, 3=H2xx, 4=H2yy, 5=H2xy,
  !                      6=H3xxy, 7=H3xyy, 8=H4xxyy
    do iDir =1, QQ
      write(dbgUnit(1),*) "iDir = ", iDir
      write(dbgUnit(1),*) "  H0 = ", H(iDir,QQ)
      write(dbgUnit(1),*) "  H1x = ", H(iDir,1)
      write(dbgUnit(1),*) "  H1y = ", H(iDir,2)
      write(dbgUnit(1),*) "  H2xx = ", H(iDir,3)
      write(dbgUnit(1),*) "  H2yy = ", H(iDir,4)
      write(dbgUnit(1),*) "  H2xy = ", H(iDir,5)
      write(dbgUnit(1),*) "  H3xxy = ", H(iDir,6)
      write(dbgUnit(1),*) "  H3xyy = ", H(iDir,7)
      write(dbgUnit(1),*) "  H4xxyy = ", H(iDir,8)
    enddo
    flush(dbgUnit(1))
    call tem_abort('done')
  endif

end subroutine getHermitepolynomials
! ************************************************************************** !


! ************************************************************************** !
!> This function computes Hermite polinomial. It gives in output minimum
!  up to 2nd-order polynomials. It is coded for up to 6th-order polynomials
! deprecated, used only to check the optimized R, RR, HRR, PRR, DRT functions
subroutine getHermitepolynomials_D3Q19( layout )
  ! --------------------------------------------------------------------------
  !> current layout
  type(mus_scheme_layout_type), intent(in) :: layout
  ! --------------------------------------------------------------------------
  integer :: iDir, QQ
  real(kind=rk) :: c_x, c_y, c_z
  real(kind=rk) :: H(19,19)
  ! --------------------------------------------------------------------------

  QQ = 19

  ! Hermite polynomials: 1=H1x, 2=H1y, 3=H1z, 4=H2xx, 5=H2yy, 6=H2zz, 7=H2xy,
  !                      8=H2xz, 9=H2yz, 10=H3xxy, 11=H3xxz, 12=H3xyy,
  !                      13=H3xzz, 14=H3yzz, 15=H3yyz,
  do iDir = 1, QQ

    c_x = layout%fStencil%cxDirRK(1,iDir)
    c_y = layout%fStencil%cxDirRK(2,iDir)
    c_z = layout%fStencil%cxDirRK(3,iDir)

    ! 0th order H0 = 1
    H(iDir,QQ) = 1._rk

    ! 1st order H1n = c_n
    H(iDir, 1) = c_x
    H(iDir, 2) = c_y
    H(iDir, 3) = c_z

    ! 2nd order
    H(iDir, 4) = c_x ** 2 - cs2
    H(iDir, 5) = c_y ** 2 - cs2
    H(iDir, 6) = c_z ** 2 - cs2
    H(iDir, 7) = c_x * c_y
    H(iDir, 8) = c_x * c_z
    H(iDir, 9) = c_y * c_z

    ! 3d order
    ! recursive formula valide from 3rd order on
    ! H (n+m+l) x_n y_m z_l = H(n)x_n H(m)y_m H(l)z_l
    H(iDir, 10) = H(iDir, 4) * c_y ! H1y = c_y
    H(iDir, 11) = H(iDir, 4) * c_z
    H(iDir, 12) = c_x * H(iDir, 5)
    H(iDir, 13) = c_x * H(iDir, 6)
    H(iDir, 14) = c_y * H(iDir, 6)
    H(iDir, 15) = H(iDir, 5) * c_z

    ! Hermite polynomials: 1=H1x, 2=H1y, 3=H1z, 4=H2xx, 5=H2yy, 6=H2zz, 7=H2xy,
    !                      8=H2xz, 9=H2yz, 10=H3xxy, 11=H3xxz, 12=H3xyy, 13=H3xzz,
    !                      14=H3yzz, 15=H3yyz
    write(dbgUnit(1),*) "iDir = ", iDir
    write(dbgUnit(1),*) "  H0 = ", H(iDir,QQ)
    write(dbgUnit(1),*) "  H1x_1 = ", H(iDir,1)
    write(dbgUnit(1),*) "  H1y_2 = ", H(iDir,2)
    write(dbgUnit(1),*) "  H1z_3 = ", H(iDir,3)
    write(dbgUnit(1),*) "  H2xx_4 = ", H(iDir,4)
    write(dbgUnit(1),*) "  H2yy_5 = ", H(iDir,5)
    write(dbgUnit(1),*) "  H2zz_6 = ", H(iDir,6)
    write(dbgUnit(1),*) "  H2xy_7 = ", H(iDir,7)
    write(dbgUnit(1),*) "  H2xz_8 = ", H(iDir,8)
    write(dbgUnit(1),*) "  H2yz_9 = ", H(iDir,9)
    write(dbgUnit(1),*) "  3(H3xxy_10 + H3yzz_14) = ", 3._rk * ( H(iDir,10) + H(iDir,14) )
    write(dbgUnit(1),*) "  (H3xxy_10 - H3yzz_14) = ", ( H(iDir,10) - H(iDir,14) )
    write(dbgUnit(1),*) "  3(H3xzz_13 + H3xyy_12) = ", 3._rk * ( H(iDir,13) + H(iDir,12) )
    write(dbgUnit(1),*) "  (H3xzz_13 - H3xyy_12) = ", ( H(iDir,13) - H(iDir,12) )
    write(dbgUnit(1),*) "  3(H3yyz_15 + H3xxz_11) = ", 3._rk * ( H(iDir,15) + H(iDir,11) )
    write(dbgUnit(1),*) "  (H3yyz_15 - H3xxz_11) = ", ( H(iDir,15) - H(iDir,11) )

  enddo
  flush(dbgUnit(1))
  call tem_abort('done')

end subroutine getHermitepolynomials_D3Q19
! ************************************************************************** !

end module mus_hrrInit_module