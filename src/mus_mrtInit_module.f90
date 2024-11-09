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
module mus_mrtInit_module

  ! include treelm modules
  use env_module,               only: rk
  use tem_param_module,         only: div2_27, div2_9, div1_9, div1_18, div1_54,   &
    &                                 div1_4, div1_6, div1_36, div1_12, div1_108,  &
    &                                 div1_48, div1_72, div1_216, div8_27, div4_9, &
    &                                 div1_27, div1_24, div1_8, div1_2, div1_3,    &
    &                                 div1_16

  implicit none

  private

  public :: check_mrt_matrix_d3q19
  public :: check_mrt_matrix_d3q27

  !=============================================================================
  ! D3Q19 flow model
  !=============================================================================

  ! D3Q19 MRT pdf -> moment transformation matrix
  ! How to use:
  ! do iDir = 1, QQ
  !   moment(iDir) = sum( PDF(:) * MMtrD3Q19(iDir,:) )
  ! end do
  !  W      S     B     E     N     T     BS    TS   BN    TN    BW    BE    TW    TE    SW    NW    SE    NE     0
  real(kind=rk), dimension(19,19),parameter,public  :: MMtrD3Q19 = &
  reshape((/ &
   1._rk,   1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk, &
     & 1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk, &
   0._rk,   0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk, &
     & 1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk, -1._rk, &
  -2._rk,  -2._rk, -2._rk, -2._rk, -2._rk, -2._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk, &
    & 1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk, &
  -1._rk,   0._rk,  0._rk,  1._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk, -1._rk, &
    & 1._rk, -1._rk,  1._rk, -1._rk, -1._rk,  1._rk,  1._rk,  0._rk, &
   2._rk,   0._rk,  0._rk, -2._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk, -1._rk, &
     & 1._rk, -1._rk,  1._rk, -1._rk, -1._rk,  1._rk,  1._rk,  0._rk, &
   0._rk,  -1._rk,  0._rk,  0._rk,  1._rk,  0._rk, -1._rk, -1._rk,  1._rk,  1._rk,  0._rk, &
     & 0._rk,  0._rk,  0._rk, -1._rk,  1._rk, -1._rk,  1._rk,  0._rk, &
   0._rk,   2._rk,  0._rk,  0._rk, -2._rk,  0._rk, -1._rk, -1._rk,  1._rk,  1._rk,  0._rk, &
     & 0._rk,  0._rk,  0._rk, -1._rk,  1._rk, -1._rk,  1._rk,  0._rk, &
   0._rk,   0._rk, -1._rk,  0._rk,  0._rk,  1._rk, -1._rk,  1._rk, -1._rk,  1._rk, -1._rk, &
     &-1._rk,  1._rk,  1._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk, &
   0._rk,   0._rk,  2._rk,  0._rk,  0._rk, -2._rk, -1._rk,  1._rk, -1._rk,  1._rk, -1._rk, &
     &-1._rk,  1._rk,  1._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk, &
   2._rk,  -1._rk, -1._rk,  2._rk, -1._rk, -1._rk, -2._rk, -2._rk, -2._rk, -2._rk,  1._rk, &
     & 1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  0._rk, &
  -2._rk,   1._rk,  1._rk, -2._rk,  1._rk,  1._rk, -2._rk, -2._rk, -2._rk, -2._rk,  1._rk, &
    & 1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  0._rk, &
   0._rk,   1._rk, -1._rk,  0._rk,  1._rk, -1._rk,  0._rk,  0._rk,  0._rk,  0._rk, -1._rk, &
     &-1._rk, -1._rk, -1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  0._rk, &
   0._rk,  -1._rk,  1._rk,  0._rk, -1._rk,  1._rk,  0._rk,  0._rk,  0._rk,  0._rk, -1._rk, &
     &-1._rk, -1._rk, -1._rk,  1._rk,  1._rk,  1._rk,  1._rk,  0._rk, &
   0._rk,   0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk, &
     & 0._rk,  0._rk,  0._rk,  1._rk, -1._rk, -1._rk,  1._rk,  0._rk, &
   0._rk,   0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  1._rk, -1._rk, -1._rk,  1._rk,  0._rk, &
     & 0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk, &
   0._rk,   0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  1._rk, &
     &-1._rk, -1._rk,  1._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk, &
   0._rk,   0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  1._rk, &
     &-1._rk,  1._rk, -1._rk, -1._rk, -1._rk,  1._rk,  1._rk,  0._rk, &
   0._rk,   0._rk,  0._rk,  0._rk,  0._rk,  0._rk, -1._rk, -1._rk,  1._rk,  1._rk,  0._rk, &
     & 0._rk,  0._rk,  0._rk,  1._rk, -1._rk,  1._rk, -1._rk,  0._rk, &
   0._rk,   0._rk,  0._rk,  0._rk,  0._rk,  0._rk,  1._rk, -1._rk,  1._rk, -1._rk, -1._rk, &
     &-1._rk,  1._rk,  1._rk,  0._rk,  0._rk,  0._rk,  0._rk,  0._rk  &
  /),(/19,19/), order=(/ 2,1 /) )

  real(kind=rk), dimension(19,19),parameter,public :: MMivD3Q19 = &
  reshape((/ &
    div1_18, 0._rk, -div1_18, -div1_6, div1_6, 0._rk, 0._rk, 0._rk, 0._rk, div1_12, -div1_12, &
     0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
    div1_18, 0._rk, -div1_18, 0._rk, 0._rk, -div1_6, div1_6, 0._rk, 0._rk, -div1_24, div1_24, &
     div1_8, -div1_8, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
    div1_18, 0._rk, -div1_18, 0._rk, 0._rk, 0._rk, 0._rk, -div1_6, div1_6, -div1_24, div1_24, &
     -div1_8, div1_8, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
    div1_18, 0._rk, -div1_18, div1_6, -div1_6, 0._rk, 0._rk, 0._rk, 0._rk, div1_12, -div1_12, &
     0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
    div1_18, 0._rk, -div1_18, 0._rk, 0._rk, div1_6, -div1_6, 0._rk, 0._rk, -div1_24, div1_24, &
     div1_8, -div1_8, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
    div1_18, 0._rk, -div1_18, 0._rk, 0._rk, 0._rk, 0._rk, div1_6, -div1_6, -div1_24, div1_24, &
     -div1_8, div1_8, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
    div1_36, div1_24, div1_72, 0._rk, 0._rk, -div1_12, -div1_24, -div1_12, -div1_24, -div1_24, -div1_24, &
     0._rk, 0._rk, 0._rk, div1_4, 0._rk, 0._rk, -div1_8, div1_8, &
    div1_36, div1_24, div1_72, 0._rk, 0._rk, -div1_12, -div1_24, div1_12, div1_24, -div1_24, -div1_24, &
     0._rk, 0._rk, 0._rk, -div1_4, 0._rk, 0._rk, -div1_8, -div1_8, &
    div1_36, div1_24, div1_72, 0._rk, 0._rk, div1_12, div1_24, -div1_12, -div1_24, -div1_24, -div1_24, &
     0._rk, 0._rk, 0._rk, -div1_4, 0._rk, 0._rk, div1_8, div1_8, &
    div1_36, div1_24, div1_72, 0._rk, 0._rk, div1_12, div1_24, div1_12, div1_24, -div1_24, -div1_24, &
     0._rk, 0._rk, 0._rk, div1_4, 0._rk, 0._rk, div1_8, -div1_8, &
    div1_36, div1_24, div1_72, -div1_12, -div1_24, 0._rk, 0._rk, -div1_12, -div1_24, div1_48, div1_48, &
     -div1_16, -div1_16, 0._rk, 0._rk, div1_4, div1_8, 0._rk, -div1_8, &
    div1_36, div1_24, div1_72, div1_12, div1_24, 0._rk, 0._rk, -div1_12, -div1_24, div1_48, div1_48, &
     -div1_16, -div1_16, 0._rk, 0._rk, -div1_4, -div1_8, 0._rk, -div1_8, &
    div1_36, div1_24, div1_72, -div1_12, -div1_24, 0._rk, 0._rk, div1_12, div1_24, div1_48, div1_48, &
     -div1_16, -div1_16, 0._rk, 0._rk, -div1_4, div1_8, 0._rk, div1_8, &
    div1_36, div1_24, div1_72, div1_12, div1_24, 0._rk, 0._rk, div1_12, div1_24, div1_48, div1_48, &
     -div1_16, -div1_16, 0._rk, 0._rk, div1_4, -div1_8, 0._rk, div1_8, &
    div1_36, div1_24, div1_72, -div1_12, -div1_24, -div1_12, -div1_24, 0._rk, 0._rk, div1_48, div1_48, &
     div1_16, div1_16, div1_4, 0._rk, 0._rk, -div1_8, div1_8, 0._rk, &
    div1_36, div1_24, div1_72, -div1_12, -div1_24, div1_12, div1_24, 0._rk, 0._rk, div1_48, div1_48, &
     div1_16, div1_16, -div1_4, 0._rk, 0._rk, -div1_8, -div1_8, 0._rk, &
    div1_36, div1_24, div1_72, div1_12, div1_24, -div1_12, -div1_24, 0._rk, 0._rk, div1_48, div1_48, &
     div1_16, div1_16, -div1_4, 0._rk, 0._rk, div1_8, div1_8, 0._rk, &
    div1_36, div1_24, div1_72, div1_12, div1_24, div1_12, div1_24, 0._rk, 0._rk, div1_48, div1_48, &
     div1_16, div1_16, div1_4, 0._rk, 0._rk, div1_8, -div1_8, 0._rk, &
    div1_3, -div1_2, div1_6, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
     0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk &
  /),(/19,19/), order=(/ 2,1 /) )

  !=============================================================================
  ! D3Q27 flow model
  !=============================================================================

  ! D3Q27 WMRT pdf -> moment transformation matrix
  ! How to use:
  ! do iDir = 1, QQ
  !   moment(iDir) = sum( PDF(:) * WMMtrD3Q27(iDir,:) )
  ! end do
  !  W      S     B     E     N     T     BS    TS   BN    TN    BW    BE    TW
  !  TE    SW    NW    SE    NE   BSW    TSW   BNW  TNW   BSE   TSE   BNE   TNE  0
  real(kind=rk), dimension(27,27), parameter, public :: WMMtrD3Q27 =                                          &
&  reshape((/                                                                                                 &
&  1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk,          &
&    1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk,               &
&  -1._rk, 0._rk, 0._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -1._rk, 1._rk, -1._rk, 1._rk,       &
&    -1._rk, -1._rk, 1._rk, 1._rk, -1._rk, -1._rk, -1._rk, -1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 0._rk,         &
&  0._rk, -1._rk, 0._rk, 0._rk, 1._rk, 0._rk, -1._rk, -1._rk, 1._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk,       &
&    -1._rk, 1._rk, -1._rk, 1._rk, -1._rk, -1._rk, 1._rk, 1._rk, -1._rk, -1._rk, 1._rk, 1._rk, 0._rk,         &
&  0._rk, 0._rk, -1._rk, 0._rk, 0._rk, 1._rk, -1._rk, 1._rk, -1._rk, 1._rk, -1._rk, -1._rk, 1._rk, 1._rk,     &
&    0._rk, 0._rk, 0._rk, 0._rk, -1._rk, 1._rk, -1._rk, 1._rk, -1._rk, 1._rk, -1._rk, 1._rk, 0._rk,           &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 1._rk,   &
&    -1._rk, -1._rk, 1._rk, 1._rk, 1._rk, -1._rk, -1._rk, -1._rk, -1._rk, 1._rk, 1._rk, 0._rk,                &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 1._rk, -1._rk, -1._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk,        &
&    0._rk, 0._rk, 0._rk, 0._rk, 1._rk, -1._rk, -1._rk, 1._rk, 1._rk, -1._rk, -1._rk, 1._rk, 0._rk,           &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 1._rk, -1._rk, -1._rk, 1._rk,        &
&    0._rk, 0._rk, 0._rk, 0._rk, 1._rk, -1._rk, 1._rk, -1._rk, -1._rk, 1._rk, -1._rk, 1._rk, 0._rk,           &
&  2._rk, -1._rk, -1._rk, 2._rk, -1._rk, -1._rk, -2._rk, -2._rk, -2._rk, -2._rk, 1._rk, 1._rk, 1._rk,         &
&    1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,        &
&  0._rk, 1._rk, -1._rk, 0._rk, 1._rk, -1._rk, 0._rk, 0._rk, 0._rk, 0._rk, -1._rk, -1._rk, -1._rk, -1._rk,    &
&    1._rk, 1._rk, 1._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,               &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk, 1._rk,   &
&    1._rk, 1._rk, 1._rk, 2._rk, 2._rk, 2._rk, 2._rk, 2._rk, 2._rk, 2._rk, 2._rk, -1._rk,                     &
&  2._rk, 0._rk, 0._rk, -2._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -1._rk, 1._rk, -1._rk, 1._rk,       &
&    -1._rk, -1._rk, 1._rk, 1._rk, -4._rk, -4._rk, -4._rk, -4._rk, 4._rk, 4._rk, 4._rk, 4._rk, 0._rk,         &
&  0._rk, 2._rk, 0._rk, 0._rk, -2._rk, 0._rk, -1._rk, -1._rk, 1._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk,       &
&    -1._rk, 1._rk, -1._rk, 1._rk, -4._rk, -4._rk, 4._rk, 4._rk, -4._rk, -4._rk, 4._rk, 4._rk, 0._rk,         &
&  0._rk, 0._rk, 2._rk, 0._rk, 0._rk, -2._rk, -1._rk, 1._rk, -1._rk, 1._rk, -1._rk, -1._rk, 1._rk, 1._rk,     &
&    0._rk, 0._rk, 0._rk, 0._rk, -4._rk, 4._rk, -4._rk, 4._rk, -4._rk, 4._rk, -4._rk, 4._rk, 0._rk,           &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 1._rk, -1._rk, 1._rk, -1._rk,        &
&    -1._rk, -1._rk, 1._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,             &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -1._rk, -1._rk, 1._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk,        &
&    1._rk, -1._rk, 1._rk, -1._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,             &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 1._rk, -1._rk, 1._rk, -1._rk, -1._rk, -1._rk, 1._rk, 1._rk,      &
&    0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,               &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   &
&    0._rk, 0._rk, 0._rk, -1._rk, 1._rk, 1._rk, -1._rk, 1._rk, -1._rk, -1._rk, 1._rk, 0._rk,                  &
&  -1._rk, -1._rk, -1._rk, -1._rk, -1._rk, -1._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,    &
&    0._rk, 0._rk, 0._rk, 0._rk, 4._rk, 4._rk, 4._rk, 4._rk, 4._rk, 4._rk, 4._rk, 4._rk, 1._rk,               &
&  -2._rk, 1._rk, 1._rk, -2._rk, 1._rk, 1._rk, -4._rk, -4._rk, -4._rk, -4._rk, 2._rk, 2._rk, 2._rk, 2._rk,    &
&    2._rk, 2._rk, 2._rk, 2._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,               &
&  0._rk, -1._rk, 1._rk, 0._rk, -1._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk, -2._rk, -2._rk, -2._rk, -2._rk,    &
&    2._rk, 2._rk, 2._rk, 2._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,               &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,          &
&    -1._rk, 1._rk, 1._rk, -1._rk, 2._rk, 2._rk, -2._rk, -2._rk, -2._rk, -2._rk, 2._rk, 2._rk, 0._rk,         &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -1._rk, 1._rk, 1._rk, -1._rk, 0._rk, 0._rk, 0._rk, 0._rk,        &
&    0._rk, 0._rk, 0._rk, 0._rk, 2._rk, -2._rk, -2._rk, 2._rk, 2._rk, -2._rk, -2._rk, 2._rk, 0._rk,           &
&  0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -1._rk, 1._rk, 1._rk, -1._rk,        &
&    0._rk, 0._rk, 0._rk, 0._rk, 2._rk, -2._rk, 2._rk, -2._rk, -2._rk, 2._rk, -2._rk, 2._rk, 0._rk,           &
&  -1._rk, 0._rk, 0._rk, 1._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 2._rk, -2._rk, 2._rk, -2._rk,       &
&    2._rk, 2._rk, -2._rk, -2._rk, -4._rk, -4._rk, -4._rk, -4._rk, 4._rk, 4._rk, 4._rk, 4._rk, 0._rk,         &
&  0._rk, -1._rk, 0._rk, 0._rk, 1._rk, 0._rk, 2._rk, 2._rk, -2._rk, -2._rk, 0._rk, 0._rk, 0._rk, 0._rk,       &
&    2._rk, -2._rk, 2._rk, -2._rk, -4._rk, -4._rk, 4._rk, 4._rk, -4._rk, -4._rk, 4._rk, 4._rk, 0._rk,         &
&  0._rk, 0._rk, -1._rk, 0._rk, 0._rk, 1._rk, 2._rk, -2._rk, 2._rk, -2._rk, 2._rk, 2._rk, -2._rk, -2._rk,     &
&    0._rk, 0._rk, 0._rk, 0._rk, -4._rk, 4._rk, -4._rk, 4._rk, -4._rk, 4._rk, -4._rk, 4._rk, 0._rk,           &
&  2._rk, 2._rk, 2._rk, 2._rk, 2._rk, 2._rk, -4._rk, -4._rk, -4._rk, -4._rk, -4._rk, -4._rk, -4._rk,          &
&    -4._rk, -4._rk, -4._rk, -4._rk, -4._rk, 8._rk, 8._rk, 8._rk, 8._rk, 8._rk, 8._rk, 8._rk, 8._rk, -1._rk   &
&  /),(/27,27/), order=(/ 2,1 /) )


  ! D3Q27 MRT moment --> PDF transformation matrix
  ! How to use:
  ! do iDir = 1, QQ
  !   fneq(iDir) = sum( WMMIvD3Q27(iDir,:) * mneq(:) )
  ! end do
  real(kind=rk), dimension(27,27),parameter,public  :: WMMIvD3Q27 =                                           &
& reshape((/                                                                                                  &
& div2_27, -div2_9, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, div1_9, 0._rk, 0._rk, div1_9, 0._rk, 0._rk, 0._rk,     &
&   0._rk, 0._rk, 0._rk, -div1_18, -div1_18, 0._rk, 0._rk, 0._rk, 0._rk, -div1_18, 0._rk, 0._rk, div1_54,     &
& div2_27, 0._rk, -div2_9, 0._rk, 0._rk, 0._rk, 0._rk, -div1_18, div1_6, 0._rk, 0._rk, div1_9, 0._rk, 0._rk,  &
&   0._rk, 0._rk, 0._rk, -div1_18, div1_36, -div1_12, 0._rk, 0._rk, 0._rk, 0._rk, -div1_18, 0._rk, div1_54,   &
& div2_27, 0._rk, 0._rk, -div2_9, 0._rk, 0._rk, 0._rk, -div1_18, -div1_6, 0._rk, 0._rk, 0._rk, div1_9, 0._rk, &
&   0._rk, 0._rk, 0._rk, -div1_18, div1_36, div1_12, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -div1_18, div1_54,    &
& div2_27, div2_9, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, div1_9, 0._rk, 0._rk, -div1_9, 0._rk, 0._rk, 0._rk,     &
&   0._rk, 0._rk, 0._rk, -div1_18, -div1_18, 0._rk, 0._rk, 0._rk, 0._rk, div1_18, 0._rk, 0._rk, div1_54,      &
& div2_27, 0._rk, div2_9, 0._rk, 0._rk, 0._rk, 0._rk, -div1_18, div1_6, 0._rk, 0._rk, -div1_9, 0._rk, 0._rk,  &
&   0._rk, 0._rk, 0._rk, -div1_18, div1_36, -div1_12, 0._rk, 0._rk, 0._rk, 0._rk, div1_18, 0._rk, div1_54,    &
& div2_27, 0._rk, 0._rk, div2_9, 0._rk, 0._rk, 0._rk, -div1_18, -div1_6, 0._rk, 0._rk, 0._rk, -div1_9, 0._rk, &
&   0._rk, 0._rk, 0._rk, -div1_18, div1_36, div1_12, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, div1_18, div1_54,     &
& div1_54, 0._rk, -div1_18, -div1_18, 0._rk, div1_6, 0._rk, -div1_36, 0._rk, div1_36, 0._rk, -div1_72,        &
&   -div1_72, 0._rk, -div1_8, div1_8, 0._rk, 0._rk, -div1_36, 0._rk, 0._rk, -div1_12, 0._rk, 0._rk, div1_36,  &
&   div1_36, -div1_108,                                                                                       &
& div1_54, 0._rk, -div1_18, div1_18, 0._rk, -div1_6, 0._rk, -div1_36, 0._rk, div1_36, 0._rk, -div1_72,        &
&   div1_72, 0._rk, -div1_8, -div1_8, 0._rk, 0._rk, -div1_36, 0._rk, 0._rk, div1_12, 0._rk, 0._rk, div1_36,   &
&   -div1_36, -div1_108,                                                                                      &
& div1_54, 0._rk, div1_18, -div1_18, 0._rk, -div1_6, 0._rk, -div1_36, 0._rk, div1_36, 0._rk, div1_72,         &
&   -div1_72, 0._rk, div1_8, div1_8, 0._rk, 0._rk, -div1_36, 0._rk, 0._rk, div1_12, 0._rk, 0._rk, -div1_36,   &
&   div1_36, -div1_108,                                                                                       &
& div1_54, 0._rk, div1_18, div1_18, 0._rk, div1_6, 0._rk, -div1_36, 0._rk, div1_36, 0._rk, div1_72, div1_72,  &
&   0._rk, div1_8, -div1_8, 0._rk, 0._rk, -div1_36, 0._rk, 0._rk, -div1_12, 0._rk, 0._rk, -div1_36, -div1_36, &
&   -div1_108,                                                                                                &
& div1_54, -div1_18, 0._rk, -div1_18, 0._rk, 0._rk, div1_6, div1_72, -div1_24, div1_36, -div1_72, 0._rk,      &
&   -div1_72, div1_8, 0._rk, -div1_8, 0._rk, 0._rk, div1_72, -div1_24, 0._rk, 0._rk, -div1_12, div1_36,       &
&   0._rk, div1_36, -div1_108,                                                                                &
& div1_54, div1_18, 0._rk, -div1_18, 0._rk, 0._rk, -div1_6, div1_72, -div1_24, div1_36, div1_72, 0._rk,       &
&   -div1_72, -div1_8, 0._rk, -div1_8, 0._rk, 0._rk, div1_72, -div1_24, 0._rk, 0._rk, div1_12, -div1_36,      &
&   0._rk, div1_36, -div1_108,                                                                                &
& div1_54, -div1_18, 0._rk, div1_18, 0._rk, 0._rk, -div1_6, div1_72, -div1_24, div1_36, -div1_72, 0._rk,      &
&   div1_72, div1_8, 0._rk, div1_8, 0._rk, 0._rk, div1_72, -div1_24, 0._rk, 0._rk, div1_12, div1_36, 0._rk,   &
&   -div1_36, -div1_108,                                                                                      &
& div1_54, div1_18, 0._rk, div1_18, 0._rk, 0._rk, div1_6, div1_72, -div1_24, div1_36, div1_72, 0._rk,         &
&   div1_72, -div1_8, 0._rk, div1_8, 0._rk, 0._rk, div1_72, -div1_24, 0._rk, 0._rk, -div1_12, -div1_36,       &
&   0._rk, -div1_36, -div1_108,                                                                               &
& div1_54, -div1_18, -div1_18, 0._rk, div1_6, 0._rk, 0._rk, div1_72, div1_24, div1_36, -div1_72, -div1_72,    &
&   0._rk, -div1_8, div1_8, 0._rk, 0._rk, 0._rk, div1_72, div1_24, -div1_12, 0._rk, 0._rk, div1_36, div1_36,  &
&   0._rk, -div1_108,                                                                                         &
& div1_54, -div1_18, div1_18, 0._rk, -div1_6, 0._rk, 0._rk, div1_72, div1_24, div1_36, -div1_72, div1_72,     &
&   0._rk, -div1_8, -div1_8, 0._rk, 0._rk, 0._rk, div1_72, div1_24, div1_12, 0._rk, 0._rk, div1_36, -div1_36, &
&   0._rk, -div1_108,                                                                                         &
& div1_54, div1_18, -div1_18, 0._rk, -div1_6, 0._rk, 0._rk, div1_72, div1_24, div1_36, div1_72, -div1_72,     &
&   0._rk, div1_8, div1_8, 0._rk, 0._rk, 0._rk, div1_72, div1_24, div1_12, 0._rk, 0._rk, -div1_36, div1_36,   &
&   0._rk, -div1_108,                                                                                         &
& div1_54, div1_18, div1_18, 0._rk, div1_6, 0._rk, 0._rk, div1_72, div1_24, div1_36, div1_72, div1_72, 0._rk, &
&   div1_8, -div1_8, 0._rk, 0._rk, 0._rk, div1_72, div1_24, -div1_12, 0._rk, 0._rk, -div1_36, -div1_36,       &
&   0._rk, -div1_108,                                                                                         &
& div1_216, -div1_72, -div1_72, -div1_72, div1_24, div1_24, div1_24, 0._rk, 0._rk, div1_72, -div1_72,         &
&   -div1_72, -div1_72, 0._rk, 0._rk, 0._rk, -div1_8, div1_72, 0._rk, 0._rk, div1_24, div1_24, div1_24,       &
&   -div1_72, -div1_72, -div1_72, div1_216,                                                                   &
& div1_216, -div1_72, -div1_72, div1_72, div1_24, -div1_24, -div1_24, 0._rk, 0._rk, div1_72, -div1_72,        &
&   -div1_72, div1_72, 0._rk, 0._rk, 0._rk, div1_8, div1_72, 0._rk, 0._rk, div1_24, -div1_24, -div1_24,       &
&   -div1_72, -div1_72, div1_72, div1_216,                                                                    &
& div1_216, -div1_72, div1_72, -div1_72, -div1_24, -div1_24, div1_24, 0._rk, 0._rk, div1_72, -div1_72,        &
&   div1_72, -div1_72, 0._rk, 0._rk, 0._rk, div1_8, div1_72, 0._rk, 0._rk, -div1_24, -div1_24, div1_24,       &
&   -div1_72, div1_72, -div1_72, div1_216,                                                                    &
& div1_216, -div1_72, div1_72, div1_72, -div1_24, div1_24, -div1_24, 0._rk, 0._rk, div1_72, -div1_72,         &
&   div1_72, div1_72, 0._rk, 0._rk, 0._rk, -div1_8, div1_72, 0._rk, 0._rk, -div1_24, div1_24, -div1_24,       &
&   -div1_72, div1_72, div1_72, div1_216,                                                                     &
& div1_216, div1_72, -div1_72, -div1_72, -div1_24, div1_24, -div1_24, 0._rk, 0._rk, div1_72, div1_72,         &
&   -div1_72, -div1_72, 0._rk, 0._rk, 0._rk, div1_8, div1_72, 0._rk, 0._rk, -div1_24, div1_24, -div1_24,      &
&   div1_72, -div1_72, -div1_72, div1_216,                                                                    &
& div1_216, div1_72, -div1_72, div1_72, -div1_24, -div1_24, div1_24, 0._rk, 0._rk, div1_72, div1_72,          &
&   -div1_72, div1_72, 0._rk, 0._rk, 0._rk, -div1_8, div1_72, 0._rk, 0._rk, -div1_24, -div1_24, div1_24,      &
&   div1_72, -div1_72, div1_72, div1_216,                                                                     &
& div1_216, div1_72, div1_72, -div1_72, div1_24, -div1_24, -div1_24, 0._rk, 0._rk, div1_72, div1_72, div1_72, &
&   -div1_72, 0._rk, 0._rk, 0._rk, -div1_8, div1_72, 0._rk, 0._rk, div1_24, -div1_24, -div1_24, div1_72,      &
&   div1_72, -div1_72, div1_216,                                                                              &
& div1_216, div1_72, div1_72, div1_72, div1_24, div1_24, div1_24, 0._rk, 0._rk, div1_72, div1_72, div1_72,    &
&   div1_72, 0._rk, 0._rk, 0._rk, div1_8, div1_72, 0._rk, 0._rk, div1_24, div1_24, div1_24, div1_72, div1_72, &
&   div1_72, div1_216,                                                                                        &
& div8_27, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -div4_9, 0._rk, 0._rk, 0._rk, 0._rk,       &
&   0._rk, 0._rk, 0._rk, div2_9, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -div1_27             &
&  /),(/27,27/), order=(/ 2,1 /) )

contains

! **************************************************************************** !
  !> Unoptimized explicit implementation
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
function check_mrt_matrix_d3q19() result (test)
    logical              :: test
    ! ---------------------------------------------------------------------------
    integer       :: iDir
    real(kind=rk) :: M_Minv
    ! ---------------------------------------------------------------------------
    test = .false.
    ! check whether the M and Minv matrices are consistent.
    ! M * M_inv = I
    do iDir = 1, 19
      M_Minv = sum( MMtrD3Q19(iDir,:)*MMIvD3Q19(:,iDir) )
      write(*,*) "row = ", iDir, "; M_Sum = ", sum( MMtrD3Q19(iDir,:) ), "&
        &; MInv_Sum = ", sum( MMIvD3Q19(iDir,:) ), "&
        &; M_dot_M_inv = ", M_Minv
      if ( abs(M_Minv - 1._rk) > 1e-15 ) then
        write(*,*) 'M * M_inv = ', M_Minv, ' along direction ', iDir
        test = .true.
      endif
    end do

  end function check_mrt_matrix_d3q19
! **************************************************************************** !

! **************************************************************************** !
  !> Unoptimized explicit implementation
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  function check_mrt_matrix_d3q27() result (test)
    logical              :: test
    ! ---------------------------------------------------------------------------
    integer       :: iDir
    real(kind=rk) :: M_Minv
    ! ---------------------------------------------------------------------------
    test = .false.
    ! check whether the M and Minv matrices are consistent.
    ! M * M_inv = I
    do iDir = 1, 27
      M_Minv = sum( WMMtrD3Q27(iDir,:)*WMMIvD3Q27(:,iDir) )
      write(*,*) "row = ", iDir, "; M_Sum = ", sum( WMMtrD3Q27(iDir,:) ), "&
        &; MInv_Sum = ", sum( WMMIvD3Q27(iDir,:) ), "&
        &; M_dot_M_inv = ", M_Minv
      if ( abs(M_Minv - 1._rk) > 1e-15 ) then
        write(*,*) 'M * M_inv = ', M_Minv, ' along direction ', iDir
        test = .true.
      endif
    end do

  end function check_mrt_matrix_d3q27
! **************************************************************************** !

end module mus_mrtInit_module