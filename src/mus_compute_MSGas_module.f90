! Copyright (c) 2012-2016, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013, 2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! ****************************************************************************** !
!> This module provides the definition and methods for multispecies gas
!! bgk advection relaxation scheme.
module mus_MSGas_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,         only: rk
  use tem_varSys_module,  only: tem_varSys_type, tem_varSys_op_type
  use tem_param_module,   only: div1_3, div3_4h, div1_4, cs2inv, t2cs2inv,     &
    &                           t2cs4inv, div1_18
  use tem_debug_module,   only: dbgUnit
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit
  use tem_math_module,    only: invert_matrix

  ! include musubi modules
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_param_module,         only: mus_param_type
  use mus_varSys_module,        only: mus_varSys_data_type
  use mus_derVarPos_module,     only: mus_derVarPos_type

  implicit none

  private

  public :: bgk_advRel_MSGas_generic
  public :: bgk_advRel_d3q19f3_MSGas

  ! =======================================================================
  ! D3Q19 flow model
  ! =======================================================================
  !> Definition of the discrete velocity set

  integer, parameter :: QQ  = 19   !< number of pdf directions
  integer, parameter :: QQN = 18   !< number of neighbors

  integer,parameter :: q__W     = 1     !< west             x-
  integer,parameter :: q__S     = 2     !< south            y-
  integer,parameter :: q__B     = 3     !< bottom           z-
  integer,parameter :: q__E     = 4     !< east             x+
  integer,parameter :: q__N     = 5     !< north            y+
  integer,parameter :: q__T     = 6     !< top              z+
  integer,parameter :: q_BS     = 7     !<                  z-,y-
  integer,parameter :: q_TS     = 8     !<                  z+,y-
  integer,parameter :: q_BN     = 9     !<                  z-,y+
  integer,parameter :: q_TN     = 10    !<                  z+,y+
  integer,parameter :: q_BW     = 11    !<                  x-,z-
  integer,parameter :: q_BE     = 12    !<                  x+,z-
  integer,parameter :: q_TW     = 13    !<                  x-,z+
  integer,parameter :: q_TE     = 14    !<                  x+,z+
  integer,parameter :: q_SW     = 15    !<                  y-,x-
  integer,parameter :: q_NW     = 16    !<                  y+,x-
  integer,parameter :: q_SE     = 17    !<                  y-,x+
  integer,parameter :: q_NE     = 18    !<                  y+,x+
  integer,parameter :: q__0     = 19    !< rest density is last

  ! integer, parameter :: QQF3  = 57   !< number of pdf directions for 3 species

contains

! ******************************************************************************
  !> Optimized Advection relaxation routine for the MSGas BGK model
  !! for d3q19 layout with three species.
  !!
  !! This routine contains the implementation of semi-implicit lattice boltzmann
  !! equation using variable transformation based on the paper
  !! "Multi-species Lattice Boltzmann Model and Practical Examples. Short Course
  !! material Pietro Asinari PhD." \n
  !! Refer page: [Multispecies](../page/features/multispecies.html) for more information
  !! In the variable tranformation steps, we can skip the step 1 and step 3
  !! and evaluate only step 2 based on tranformed variable g
  !! only prerequisite is to compute feq which depends on original f not on g.
  !! feq is depend on density and velocity. Where density can be computed
  !! directly from g and velocity computed from linear system of equation
  !! given in the reference page  [Multispecies](../page/features/multispecies.html).
  !! KM: This is an non-optimized kernel
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_advRel_d3q19f3_MSGas( fieldProp, inState, outState, auxField, &
    &                                  neigh, nElems, nSolve, level, layout,   &
    &                                  params, varSys, derVarPos               )
    ! -------------------------------------------------------------------- !
    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> current layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> number of elements in state Array
    integer, intent(in) :: nElems
    !> input  pdf vector
    real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
    !> output pdf vector
    real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
    !> Auxiliary field computed from pre-collision state
    !! Is updated with correct velocity field for multicomponent models
    real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
    !> connectivity vector
    integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
    !> number of elements solved in kernel
    integer, intent(in) :: nSolve
    !> current level
    integer,intent(in) :: level
    !> global parameters
    type(mus_param_type),intent(in) :: params
    !> position of derived quantities in varsys for all fields
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    ! temporary local pdf values
    real(kind=rk) :: pdfTmp( QQ, 3 )
    integer :: iElem, nFields, nScalars, iFld
    !type(mus_varSys_data_type), pointer :: fPtr
    !type(mus_scheme_type), pointer :: scheme
    real(kind=rk), dimension(3) :: rsigma, uxsigma, uysigma, uzsigma,          &
      &                            qxsigma, qysigma, qzsigma,                  &
      &                            gqxsigma, gqysigma, gqzsigma
    real(kind=rk), dimension(3) :: uxstar, uystar, uzstar
    real(kind=rk), dimension(3) :: m_sigma, m_rsig, m_inv
    real(kind=rk) :: r, p, mm, usqr(3), ucx(6,3), theta
    real(kind=rk), dimension(3,3) :: A, B, chi
    real(kind=rk) :: zfac, yfac, a12_11, a13_11, a23fac, a21_11, a31_11, zcoeff
    real(kind=rk) :: pdivr, r1divr, r2divr, r3divr, B_inv(3), r_inv,           &
      &              a11_inv, a22fac_inv
    real(kind=rk) :: chi12, chi13, chi21, chi23, chi31, chi32
    real(kind=rk),dimension(3) :: lamfac_o, lamrho, lamrho_d, lamrho_dm,       &
      &                           lamrho_o, lamrho_om, lambda, lam_fac,        &
      &                           fac_1, fac_2, fac_3, fac_4, fac_5, fac_6,    &
      &                           fac_7, fac_8, fac_9, sum1_1, sum2_1, sum3_1, &
      &                           sum4_1, sum5_1, sum6_1, sum7_1, sum8_1,      &
      &                           sum9_1, sum1_2, sum2_2, sum3_2, sum4_2,      &
      &                           sum5_2, sum6_2, sum7_2, sum8_2, sum9_2
    real(kind=rk) :: mbb(3,3)
    integer :: dens_pos(3), mom_pos(3,3), elemOff
    ! ---------------------------------------------------------------------------
    !call tem_abort('Error: bgk_advRel_d3q19f3_MSGas need to be tested. Refered to scheme(1) beforehand.')

    ! access scheme via 1st variable method data which is a state variable
    !call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    !scheme => fPtr%solverData%scheme

    !semi-implicit lbm is recovered for theta = 0.5
    theta = 0.5_rk
    nFields = varSys%nStateVars
    nScalars = varSys%nScalars
    !molecular weights
    m_sigma = fieldProp(:)%species%molweight
    m_inv = 1.0_rk/m_sigma
    !molecular weight ratios
    m_rsig = fieldProp(:)%species%molWeigRatio
    !resistivities
    B(1,1) = fieldProp( 1 )%species%resi_coeff(1)
    B(1,2) = fieldProp( 1 )%species%resi_coeff(2)
    B(1,3) = fieldProp( 1 )%species%resi_coeff(3)

    B(2,1) = fieldProp( 2 )%species%resi_coeff(1)
    B(2,2) = fieldProp( 2 )%species%resi_coeff(2)
    B(2,3) = fieldProp( 2 )%species%resi_coeff(3)

    B(3,1) = fieldProp( 3 )%species%resi_coeff(1)
    B(3,2) = fieldProp( 3 )%species%resi_coeff(2)
    B(3,3) = fieldProp( 3 )%species%resi_coeff(3)
    !diffusivities
    B_inv(1) = fieldProp( 1 )%species%diff_coeff(1)
    B_inv(2) = fieldProp( 2 )%species%diff_coeff(2)
    B_inv(3) = fieldProp( 3 )%species%diff_coeff(3)

    mbb(1,1) = m_inv(1)*m_inv(1)*B(1,1)*B_inv(1)
    mbb(1,2) = m_inv(1)*m_inv(2)*B(1,2)*B_inv(1)
    mbb(1,3) = m_inv(1)*m_inv(3)*B(1,3)*B_inv(1)

    mbb(2,1) = m_inv(2)*m_inv(1)*B(2,1)*B_inv(2)
    mbb(2,2) = m_inv(2)*m_inv(2)*B(2,2)*B_inv(2)
    mbb(2,3) = m_inv(2)*m_inv(3)*B(2,3)*B_inv(2)

    mbb(3,1) = m_inv(3)*m_inv(1)*B(3,1)*B_inv(3)
    mbb(3,2) = m_inv(3)*m_inv(2)*B(3,2)*B_inv(3)
    mbb(3,3) = m_inv(3)*m_inv(3)*B(3,3)*B_inv(3)

    ! position of density and momentum in auxField array
    do iFld = 1, nFields
      dens_pos(iFld) = varSys%method%val(derVarPos(iFld)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:,iFld) = varSys%method%val(derVarPos(iFld)%momentum) &
        &                           %auxField_varPos(:)
    end do

    nodeloop: do iElem = 1, nSolve

      ! species 1
      pdfTmp( Q__W,1 ) = instate(  neigh (( q__w-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__S,1 ) = instate(  neigh (( q__s-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__B,1 ) = instate(  neigh (( q__b-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__E,1 ) = instate(  neigh (( q__e-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__N,1 ) = instate(  neigh (( q__n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__T,1 ) = instate(  neigh (( q__t-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,1 ) = instate(  neigh (( q_bs-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,1 ) = instate(  neigh (( q_ts-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,1 ) = instate(  neigh (( q_bn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,1 ) = instate(  neigh (( q_tn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,1 ) = instate(  neigh (( q_bw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,1 ) = instate(  neigh (( q_be-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,1 ) = instate(  neigh (( q_tw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,1 ) = instate(  neigh (( q_te-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,1 ) = instate(  neigh (( q_sw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,1 ) = instate(  neigh (( q_nw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,1 ) = instate(  neigh (( q_se-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,1 ) = instate(  neigh (( q_ne-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__0,1 ) = instate(  neigh (( q__0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      ! species 2
      pdfTmp( Q__W,2 ) = instate(  neigh (( q__w-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__S,2 ) = instate(  neigh (( q__s-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__B,2 ) = instate(  neigh (( q__b-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__E,2 ) = instate(  neigh (( q__e-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__N,2 ) = instate(  neigh (( q__n-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__T,2 ) = instate(  neigh (( q__t-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,2 ) = instate(  neigh (( q_bs-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,2 ) = instate(  neigh (( q_ts-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,2 ) = instate(  neigh (( q_bn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,2 ) = instate(  neigh (( q_tn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,2 ) = instate(  neigh (( q_bw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,2 ) = instate(  neigh (( q_be-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,2 ) = instate(  neigh (( q_tw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,2 ) = instate(  neigh (( q_te-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,2 ) = instate(  neigh (( q_sw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,2 ) = instate(  neigh (( q_nw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,2 ) = instate(  neigh (( q_se-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,2 ) = instate(  neigh (( q_ne-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__0,2 ) = instate(  neigh (( q__0-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)

      ! species 3
      pdfTmp( Q__W,3 ) = instate(  neigh (( q__w-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__S,3 ) = instate(  neigh (( q__s-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__B,3 ) = instate(  neigh (( q__b-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__E,3 ) = instate(  neigh (( q__e-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__N,3 ) = instate(  neigh (( q__n-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__T,3 ) = instate(  neigh (( q__t-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,3 ) = instate(  neigh (( q_bs-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,3 ) = instate(  neigh (( q_ts-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,3 ) = instate(  neigh (( q_bn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,3 ) = instate(  neigh (( q_tn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,3 ) = instate(  neigh (( q_bw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,3 ) = instate(  neigh (( q_be-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,3 ) = instate(  neigh (( q_tw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,3 ) = instate(  neigh (( q_te-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,3 ) = instate(  neigh (( q_sw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,3 ) = instate(  neigh (( q_nw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,3 ) = instate(  neigh (( q_se-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,3 ) = instate(  neigh (( q_ne-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__0,3 ) = instate(  neigh (( q__0-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)

      ! element offset for auxField
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! density 1
      rsigma(1) = auxField(elemOff + dens_pos(1))
      ! density 2
      rsigma(2) = auxField(elemOff + dens_pos(2))
      ! density 3
      rsigma(3) = auxField(elemOff + dens_pos(3))

      ! local x,y,z - momentum of species 1
      gqxsigma(1) = auxField(elemOff + mom_pos(1,1))
      gqysigma(1) = auxField(elemOff + mom_pos(2,1))
      gqzsigma(1) = auxField(elemOff + mom_pos(3,1))

      ! local x,y,z - momentum of species 2
      gqxsigma(2) = auxField(elemOff + mom_pos(1,2))
      gqysigma(2) = auxField(elemOff + mom_pos(2,2))
      gqzsigma(2) = auxField(elemOff + mom_pos(3,2))

      ! local x,y,z - momentum of species 3
      gqxsigma(3) = auxField(elemOff + mom_pos(1,3))
      gqysigma(3) = auxField(elemOff + mom_pos(2,3))
      gqzsigma(3) = auxField(elemOff + mom_pos(3,3))

      !total density and total pressure
      r = rsigma(1)+rsigma(2)+rsigma(3)
      r_inv = 1.0_rk/r

      p = (rsigma(1) * m_rsig(1)           &
        &  + rsigma(2) * m_rsig(2)         &
        &  + rsigma(3) * m_rsig(3)) * div1_3

      pdivr = p*r_inv
      r1divr = rsigma(1) * r_inv
      r2divr = rsigma(2) * r_inv
      r3divr = rsigma(3) * r_inv
      !relaxation time
      lambda(1) = pdivr * B(1,1)
      lambda(2) = pdivr * B(2,2)
      lambda(3) = pdivr * B(3,3)

      lam_fac(:) = lambda(:) / ( 1.0_rk + theta * lambda(:) )
      lamfac_o(:) = 1.0_rk - lam_fac(:)
      lamrho(:) = rsigma(:)*lam_fac(:)
      lamrho_d(:) = lamrho(:) * div1_4

      !Mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      mm = 1._rk / ( r1divr*m_inv(1) + r2divr*m_inv(2) + r3divr*m_inv(3) )
!      write(*,*) 'mm', mm
      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      chi = mm*mm * mbb

!      write(*,*) 'chi', chi
!      write(*,*) 'lambda', lambda
      !matrix A of the velocity lse
      A(1,1) = 1._rk + theta * (lambda(1) * (chi(1,2)*r2divr + chi(1,3)*r3divr))
      A(1,2) = - theta * lambda(1) * r1divr*chi(1,2)
      A(1,3) = - theta * lambda(1) * r1divr*chi(1,3)

      A(2,1) = - theta * lambda(2) * r2divr*chi(2,1)
      A(2,2) = 1._rk + theta * (lambda(2) * (chi(2,1)*r1divr + chi(2,3)*r3divr))
      A(2,3) = - theta * lambda(2) * r2divr*chi(2,3)

      A(3,1) = - theta * lambda(3) * r3divr*chi(3,1)
      A(3,2) = - theta * lambda(3) * r3divr*chi(3,2)
      A(3,3) = 1._rk + theta * (lambda(3) * (chi(3,1)*r1divr + chi(3,2)*r2divr))

      a11_inv = 1.0_rk/A(1,1)
      a13_11 = A(1,3)*a11_inv
      a12_11 = A(1,2)*a11_inv
      a31_11 = A(3,1)*a11_inv
      a21_11 = A(2,1)*a11_inv
      a23fac = A(2,3) - A(2,1)*a13_11
      a22fac_inv = 1.0_rk/(A(2,2) - A(2,1)*a12_11)
!
      zfac = (A(3,1)*a12_11 - A(3,2)) * a22fac_inv
      zcoeff = 1._rk / (a23fac*zfac - A(3,1)*a13_11 + A(3,3))
!      !x-momentum
!      !species 3
      yfac = gqxsigma(2) - gqxsigma(1)*a21_11
!
      qxsigma(3) = (gqxsigma(3) - a31_11*gqxsigma(1) + yfac*zfac)*zcoeff
!      write(*,*)'z',qxsigma(3)
!
!      !species 2
      qxsigma(2) = (yfac - a23fac*qxsigma(3))*a22fac_inv
!      !species 1
      qxsigma(1) = (gqxsigma(1) - A(1,2)*qxsigma(2) - A(1,3)*qxsigma(3))*a11_inv
!      write(*,*) 'qxsigma',qxsigma
!      !y-momentum
      yfac = gqysigma(2) - gqysigma(1)*a21_11
      qysigma(3) = (gqysigma(3) - a31_11*gqysigma(1) + yfac*zfac)*zcoeff
!
      qysigma(2) = (yfac - a23fac*qysigma(3))*a22fac_inv
!
      qysigma(1) = (gqysigma(1) - A(1,2)*qysigma(2) - A(1,3)*qysigma(3))*a11_inv
!      write(*,*) 'qysigma',qysigma
!      !z-momentum
      yfac = gqzsigma(2) - gqzsigma(1)*a21_11
      qzsigma(3) = (gqzsigma(3) - a31_11*gqzsigma(1) + yfac*zfac)*zcoeff
!
      qzsigma(2) = (yfac - a23fac*qzsigma(3))*a22fac_inv
!
      qzsigma(1) = (gqzsigma(1) - A(1,2)*qzsigma(2) - A(1,3)*qzsigma(3))*a11_inv
!      write(*,*) 'qzsigma',qzsigma

      uxsigma = qxsigma/rsigma
      uysigma = qysigma/rsigma
      uzsigma = qzsigma/rsigma

      ! store momentum of untransformed PDF in auxField
      ! species 1
      auxField(elemOff + mom_pos(1,1)) = qxsigma(1)
      auxField(elemOff + mom_pos(2,1)) = qysigma(1)
      auxField(elemOff + mom_pos(3,1)) = qzsigma(1)
      ! species 1
      auxField(elemOff + mom_pos(1,2)) = qxsigma(2)
      auxField(elemOff + mom_pos(2,2)) = qysigma(2)
      auxField(elemOff + mom_pos(3,2)) = qzsigma(2)
      ! species 1
      auxField(elemOff + mom_pos(1,3)) = qxsigma(3)
      auxField(elemOff + mom_pos(2,3)) = qysigma(3)
      auxField(elemOff + mom_pos(3,3)) = qzsigma(3)

      ! compute ustar to compute feq
      chi12 = chi(1,2)*r2divr
      chi13 = chi(1,3)*r3divr
      chi21 = chi(2,1)*r1divr
      chi23 = chi(2,3)*r3divr
      chi31 = chi(3,1)*r1divr
      chi32 = chi(3,2)*r2divr

      ! ux
      uxstar(1) =  uxsigma(1) + chi12*(uxsigma(2)-uxsigma(1))                  &
        &                     + chi13*(uxsigma(3)-uxsigma(1))

      uxstar(2) =  uxsigma(2) + chi21*(uxsigma(1)-uxsigma(2))                  &
        &                     + chi23*(uxsigma(3)-uxsigma(2))

      uxstar(3) =  uxsigma(3) + chi31*(uxsigma(1)-uxsigma(3))                  &
        &                     + chi32*(uxsigma(2)-uxsigma(3))

      ! uy
      uystar(1) =  uysigma(1) + chi12*(uysigma(2)-uysigma(1))                  &
        &                     + chi13*(uysigma(3)-uysigma(1))

      uystar(2) =  uysigma(2) + chi21*(uysigma(1)-uysigma(2))                  &
        &                     + chi23*(uysigma(3)-uysigma(2))

      uystar(3) =  uysigma(3) + chi31*(uysigma(1)-uysigma(3))                  &
        &                     + chi32*(uysigma(2)-uysigma(3))

      ! uz
      uzstar(1) =  uzsigma(1) + chi12*(uzsigma(2)-uzsigma(1))                  &
        &                     + chi13*(uzsigma(3)-uzsigma(1))

      uzstar(2) =  uzsigma(2) + chi21*(uzsigma(1)-uzsigma(2))                  &
        &                     + chi23*(uzsigma(3)-uzsigma(2))

      uzstar(3) =  uzsigma(3) + chi31*(uzsigma(1)-uzsigma(3))                  &
        &                     + chi32*(uzsigma(2)-uzsigma(3))

!      write(*,*) 'uxstar', uxstar
!      write(*,*) 'uystar', uystar
!      write(*,*) 'uzstar', uzstar
      !compute equilibrium and do collision
      usqr(1) = (uxstar(1) * uxstar(1)                                         &
        &     + uystar(1) * uystar(1)                                          &
        &     + uzstar(1) * uzstar(1)) * t2cs2inv

      ! equilibrium at rest
      outstate( ( ielem-1)* nscalars+ q__0+( 1-1)* qq ) =       &
              & lamfac_o(1)*pdfTmp(Q__0,1) +                                   &
              & div1_3 * lamrho(1) * (( 3._rk - 2._rk * M_rsig(1)) - usqr(1) )

      lamrho_dm(1) = lamrho(1) * (m_rsig(1) - usqr(1)) * div1_18

      !directional velocity factor

      ! species 1
      fac_1(1) = lamrho_d(1) * uxstar(1)
      sum1_1(1) = fac_1(1)*div3_4h
      sum1_2(1) = fac_1(1)*uxstar(1) + lamrho_dm(1)

      outstate( ( ielem-1)* nscalars+ q__w+( 1-1)* qq )         &
        & = lamfac_o(1) * pdfTmp(Q__W,1) - sum1_1(1) + sum1_2(1)
      outstate( ( ielem-1)* nscalars+ q__e+( 1-1)* qq )         &
        & = lamfac_o(1) * pdfTmp(Q__E,1) + sum1_1(1) + sum1_2(1)

      fac_2(1) = lamrho_d(1) * uystar(1)
      sum2_1(1) = fac_2(1)*div3_4h
      sum2_2(1) = fac_2(1)*uystar(1) + lamrho_dm(1)

      outstate( ( ielem-1)* nscalars+ q__s+( 1-1)* qq )         &
        & = lamfac_o(1) * pdfTmp(Q__S,1) - sum2_1(1) + sum2_2(1)
      outstate( ( ielem-1)* nscalars+ q__n+( 1-1)* qq )         &
        & = lamfac_o(1) * pdfTmp(Q__N,1) + sum2_1(1) + sum2_2(1)

      fac_3(1) = lamrho_d(1) * uzstar(1)
      sum3_1(1) = fac_3(1)*div3_4h
      sum3_2(1) = fac_3(1)*uzstar(1) + lamrho_dm(1)

      outstate( ( ielem-1)* nscalars+ q__b+( 1-1)* qq )         &
        & = lamfac_o(1) * pdfTmp(Q__B,1) - sum3_1(1) + sum3_2(1)
      outstate( ( ielem-1)* nscalars+ q__t+( 1-1)* qq )         &
        & = lamfac_o(1) * pdfTmp(Q__T,1) + sum3_1(1) + sum3_2(1)

      !top north / bottom south
      lamrho_o(1) = lamrho_d(1) * 0.5_rk
      lamrho_om(1) = lamrho_dm(1) * 0.5_rk
      ucx( 1, 1 ) =             + uystar(1) + uzstar(1)
      fac_4(1) = lamrho_o(1) * ucx(1,1)
      sum4_1(1) = fac_4(1)*div3_4h
      sum4_2(1) = fac_4(1)*ucx(1,1) + lamrho_om(1)

      outstate( ( ielem-1)* nscalars+ q_bs+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_BS,1) - sum4_1(1) + sum4_2(1)
      outstate( ( ielem-1)* nscalars+ q_tn+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_TN,1) + sum4_1(1) + sum4_2(1)

      !top south / bottom north
      ucx( 2, 1 ) =             - uystar(1) + uzstar(1)
      fac_5(1) = lamrho_o(1) * ucx(2,1)
      sum5_1(1) = fac_5(1)*div3_4h
      sum5_2(1) = fac_5(1)*ucx(2,1) + lamrho_om(1)

      outstate( ( ielem-1)* nscalars+ q_bn+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_BN,1) - sum5_1(1) + sum5_2(1)
      outstate( ( ielem-1)* nscalars+ q_ts+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_TS,1) + sum5_1(1) + sum5_2(1)

      !top east / bottom west
      ucx( 3, 1 ) = + uxstar(1)             + uzstar(1)
      fac_6(1) = lamrho_o(1) * ucx(3,1)
      sum6_1(1) = fac_6(1)*div3_4h
      sum6_2(1) = fac_6(1)*ucx(3,1) + lamrho_om(1)

      outstate( ( ielem-1)* nscalars+ q_bw+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_BW,1) - sum6_1(1) + sum6_2(1)
      outstate( ( ielem-1)* nscalars+ q_te+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_TE,1) + sum6_1(1) + sum6_2(1)

      !top west / bottom east
      ucx( 4, 1 ) = - uxstar(1)             + uzstar(1)
      fac_7(1) = lamrho_o(1) * ucx(4,1)
      sum7_1(1) = fac_7(1)*div3_4h
      sum7_2(1) = fac_7(1)*ucx(4,1) + lamrho_om(1)

      outstate( ( ielem-1)* nscalars+ q_be+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_BE,1) - sum7_1(1) + sum7_2(1)
      outstate( ( ielem-1)* nscalars+ q_tw+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_TW,1) + sum7_1(1) + sum7_2(1)

      !north east / south west
      ucx( 5, 1 ) = + uxstar(1) + uystar(1)
      fac_8(1) = lamrho_o(1) * ucx(5,1)
      sum8_1(1) = fac_8(1)*div3_4h
      sum8_2(1) = fac_8(1)*ucx(5,1) + lamrho_om(1)

      outstate( ( ielem-1)* nscalars+ q_sw+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_SW,1) - sum8_1(1) + sum8_2(1)
      outstate( ( ielem-1)* nscalars+ q_ne+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_NE,1) + sum8_1(1) + sum8_2(1)

      !north west / south east
      ucx( 6, 1 ) = - uxstar(1) + uystar(1)
      fac_9(1) = lamrho_o(1) * ucx(6,1)
      sum9_1(1) = fac_9(1)*div3_4h
      sum9_2(1) = fac_9(1)*ucx(6,1) + lamrho_om(1)

      outstate( ( ielem-1)* nscalars+ q_se+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_SE,1) - sum9_1(1) + sum9_2(1)
      outstate( ( ielem-1)* nscalars+ q_nw+( 1-1)* qq )         &
        & = lamfac_o(1)*pdfTmp(Q_NW,1) + sum9_1(1) + sum9_2(1)


      ! species 2
      usqr(2) = (uxstar(2) * uxstar(2)                                         &
        &     + uystar(2) * uystar(2)                                          &
        &     + uzstar(2) * uzstar(2)) * t2cs2inv

      outstate( ( ielem-1)* nscalars+ q__0+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q__0,2)                                         &
        & + div1_3 * lamrho(2) * (( 3._rk - 2._rk * M_rsig(2)) - usqr(2) )

      lamrho_dm(2) = lamrho(2) * (m_rsig(2) - usqr(2)) * div1_18
      fac_1(2) = lamrho_d(2) * uxstar(2)
      sum1_1(2) = fac_1(2)*div3_4h
      sum1_2(2) = fac_1(2)*uxstar(2) + lamrho_dm(2)

      outstate( ( ielem-1)* nscalars+ q__w+( 2-1)* qq )         &
        & = lamfac_o(2) * pdfTmp(Q__W,2) - sum1_1(2) + sum1_2(2)
      outstate( ( ielem-1)* nscalars+ q__e+( 2-1)* qq )         &
        & = lamfac_o(2) * pdfTmp(Q__E,2) + sum1_1(2) + sum1_2(2)

      fac_2(2) = lamrho_d(2) * uystar(2)
      sum2_1(2) = fac_2(2)*div3_4h
      sum2_2(2) = fac_2(2)*uystar(2) + lamrho_dm(2)

      outstate( ( ielem-1)* nscalars+ q__s+( 2-1)* qq )         &
        & = lamfac_o(2) * pdfTmp(Q__S,2) - sum2_1(2) + sum2_2(2)
      outstate( ( ielem-1)* nscalars+ q__n+( 2-1)* qq )         &
        & = lamfac_o(2) * pdfTmp(Q__N,2) + sum2_1(2) + sum2_2(2)

      fac_3(2) = lamrho_d(2) * uzstar(2)
      sum3_1(2) = fac_3(2)*div3_4h
      sum3_2(2) = fac_3(2)*uzstar(2) + lamrho_dm(2)

      outstate( ( ielem-1)* nscalars+ q__b+( 2-1)* qq )         &
        & = lamfac_o(2) * pdfTmp(Q__B,2) - sum3_1(2) + sum3_2(2)
      outstate( ( ielem-1)* nscalars+ q__t+( 2-1)* qq )         &
        & = lamfac_o(2) * pdfTmp(Q__T,2) + sum3_1(2) + sum3_2(2)

      !top north / bottom south
      lamrho_o(2) = lamrho_d(2) * 0.5_rk
      lamrho_om(2) = lamrho_dm(2) * 0.5_rk
      ucx( 1, 2 ) =             + uystar(2) + uzstar(2)
      fac_4(2) = lamrho_o(2) * ucx(1,2)
      sum4_1(2) = fac_4(2)*div3_4h
      sum4_2(2) = fac_4(2)*ucx(1,2) + lamrho_om(2)

      outstate( ( ielem-1)* nscalars+ q_bs+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_BS,2) - sum4_1(2) + sum4_2(2)
      outstate( ( ielem-1)* nscalars+ q_tn+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_TN,2) + sum4_1(2) + sum4_2(2)

      !top south / bottom north
      ucx( 2, 2 ) =             - uystar(2) + uzstar(2)
      fac_5(2) = lamrho_o(2) * ucx(2,2)
      sum5_1(2) = fac_5(2)*div3_4h
      sum5_2(2) = fac_5(2)*ucx(2,2) + lamrho_om(2)

      outstate( ( ielem-1)* nscalars+ q_bn+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_BN,2) - sum5_1(2) + sum5_2(2)
      outstate( ( ielem-1)* nscalars+ q_ts+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_TS,2) + sum5_1(2) + sum5_2(2)

      !top east / bottom west
      ucx( 3, 2 ) = + uxstar(2)             + uzstar(2)
      fac_6(2) = lamrho_o(2) * ucx(3,2)
      sum6_1(2) = fac_6(2)*div3_4h
      sum6_2(2) = fac_6(2)*ucx(3,2) + lamrho_om(2)

      outstate( ( ielem-1)* nscalars+ q_bw+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_BW,2) - sum6_1(2) + sum6_2(2)
      outstate( ( ielem-1)* nscalars+ q_te+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_TE,2) + sum6_1(2) + sum6_2(2)

      !top west / bottom east
      ucx( 4, 2 ) = - uxstar(2)             + uzstar(2)
      fac_7(2) = lamrho_o(2) * ucx(4,2)
      sum7_1(2) = fac_7(2)*div3_4h
      sum7_2(2) = fac_7(2)*ucx(4,2) + lamrho_om(2)

      outstate( ( ielem-1)* nscalars+ q_be+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_BE,2) - sum7_1(2) + sum7_2(2)
      outstate( ( ielem-1)* nscalars+ q_tw+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_TW,2) + sum7_1(2) + sum7_2(2)

      !north east / south west
      ucx( 5, 2 ) = + uxstar(2) + uystar(2)
      fac_8(2) = lamrho_o(2) * ucx(5,2)
      sum8_1(2) = fac_8(2)*div3_4h
      sum8_2(2) = fac_8(2)*ucx(5,2) + lamrho_om(2)

      outstate( ( ielem-1)* nscalars+ q_sw+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_SW,2) - sum8_1(2) + sum8_2(2)
      outstate( ( ielem-1)* nscalars+ q_ne+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_NE,2) + sum8_1(2) + sum8_2(2)

      !north west / south east
      ucx( 6, 2 ) = - uxstar(2) + uystar(2)
      fac_9(2) = lamrho_o(2) * ucx(6,2)
      sum9_1(2) = fac_9(2)*div3_4h
      sum9_2(2) = fac_9(2)*ucx(6,2) + lamrho_om(2)

      outstate( ( ielem-1)* nscalars+ q_se+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_SE,2) - sum9_1(2) + sum9_2(2)
      outstate( ( ielem-1)* nscalars+ q_nw+( 2-1)* qq )         &
        & = lamfac_o(2)*pdfTmp(Q_NW,2) + sum9_1(2) + sum9_2(2)

      ! species 3
      usqr(3) = (uxstar(3) * uxstar(3)                                         &
        &     + uystar(3) * uystar(3)                                          &
        &     + uzstar(3) * uzstar(3)) * t2cs2inv

      outstate( ( ielem-1)* nscalars+ q__0+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q__0,3)                                         &
        & + div1_3 * lamrho(3) * (( 3._rk - 2._rk * M_rsig(3)) - usqr(3) )

      lamrho_dm(3) = lamrho(3) * (m_rsig(3) - usqr(3)) * div1_18
      fac_1(3) = lamrho_d(3) * uxstar(3)
      sum1_1(3) = fac_1(3)*div3_4h
      sum1_2(3) = fac_1(3)*uxstar(3) + lamrho_dm(3)

      outstate( ( ielem-1)* nscalars+ q__w+( 3-1)* qq )         &
        & = lamfac_o(3) * pdfTmp(Q__W,3) - sum1_1(3) + sum1_2(3)
      outstate( ( ielem-1)* nscalars+ q__e+( 3-1)* qq )         &
        & = lamfac_o(3) * pdfTmp(Q__E,3) + sum1_1(3) + sum1_2(3)

      fac_2(3) = lamrho_d(3) * uystar(3)
      sum2_1(3) = fac_2(3)*div3_4h
      sum2_2(3) = fac_2(3)*uystar(3) + lamrho_dm(3)

      outstate( ( ielem-1)* nscalars+ q__s+( 3-1)* qq )         &
        & = lamfac_o(3) * pdfTmp(Q__S,3) - sum2_1(3) + sum2_2(3)
      outstate( ( ielem-1)* nscalars+ q__n+( 3-1)* qq )         &
        & = lamfac_o(3) * pdfTmp(Q__N,3) + sum2_1(3) + sum2_2(3)

      fac_3(3) = lamrho_d(3) * uzstar(3)
      sum3_1(3) = fac_3(3)*div3_4h
      sum3_2(3) = fac_3(3)*uzstar(3) + lamrho_dm(3)

      outstate( ( ielem-1)* nscalars+ q__b+( 3-1)* qq )         &
        & = lamfac_o(3) * pdfTmp(Q__B,3) - sum3_1(3) + sum3_2(3)
      outstate( ( ielem-1)* nscalars+ q__t+( 3-1)* qq )         &
        & = lamfac_o(3) * pdfTmp(Q__T,3) + sum3_1(3) + sum3_2(3)

      !top north / bottom south
      lamrho_o(3) = lamrho_d(3) * 0.5_rk
      lamrho_om(3) = lamrho_dm(3) * 0.5_rk
      ucx( 1, 3 ) =             + uystar(3) + uzstar(3)
      fac_4(3) = lamrho_o(3) * ucx(1,3)
      sum4_1(3) = fac_4(3)*div3_4h
      sum4_2(3) = fac_4(3)*ucx(1,3) + lamrho_om(3)

      outstate( ( ielem-1)* nscalars+ q_bs+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_BS,3) - sum4_1(3) + sum4_2(3)
      outstate( ( ielem-1)* nscalars+ q_tn+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_TN,3) + sum4_1(3) + sum4_2(3)

      !top south / bottom north
      ucx( 2, 3 ) =             - uystar(3) + uzstar(3)
      fac_5(3) = lamrho_o(3) * ucx(2,3)
      sum5_1(3) = fac_5(3)*div3_4h
      sum5_2(3) = fac_5(3)*ucx(2,3) + lamrho_om(3)

      outstate( ( ielem-1)* nscalars+ q_bn+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_BN,3) - sum5_1(3) + sum5_2(3)
      outstate( ( ielem-1)* nscalars+ q_ts+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_TS,3) + sum5_1(3) + sum5_2(3)

      !top east / bottom west
      ucx( 3, 3 ) = + uxstar(3)             + uzstar(3)
      fac_6(3) = lamrho_o(3) * ucx(3,3)
      sum6_1(3) = fac_6(3)*div3_4h
      sum6_2(3) = fac_6(3)*ucx(3,3) + lamrho_om(3)

      outstate( ( ielem-1)* nscalars+ q_bw+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_BW,3) - sum6_1(3) + sum6_2(3)
      outstate( ( ielem-1)* nscalars+ q_te+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_TE,3) + sum6_1(3) + sum6_2(3)

      !top west / bottom east
      ucx( 4, 3 ) = - uxstar(3)             + uzstar(3)
      fac_7(3) = lamrho_o(3) * ucx(4,3)
      sum7_1(3) = fac_7(3)*div3_4h
      sum7_2(3) = fac_7(3)*ucx(4,3) + lamrho_om(3)

      outstate( ( ielem-1)* nscalars+ q_be+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_BE,3) - sum7_1(3) + sum7_2(3)
      outstate( ( ielem-1)* nscalars+ q_tw+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_TW,3) + sum7_1(3) + sum7_2(3)

      !north east / south west
      ucx( 5, 3 ) = + uxstar(3) + uystar(3)
      fac_8(3) = lamrho_o(3) * ucx(5,3)
      sum8_1(3) = fac_8(3)*div3_4h
      sum8_2(3) = fac_8(3)*ucx(5,3) + lamrho_om(3)

      outstate( ( ielem-1)* nscalars+ q_sw+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_SW,3) - sum8_1(3) + sum8_2(3)
      outstate( ( ielem-1)* nscalars+ q_ne+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_NE,3) + sum8_1(3) + sum8_2(3)

      !north west / south east
      ucx( 6, 3 ) = - uxstar(3) + uystar(3)
      fac_9(3) = lamrho_o(3) * ucx(6,3)
      sum9_1(3) = fac_9(3)*div3_4h
      sum9_2(3) = fac_9(3)*ucx(6,3) + lamrho_om(3)

      outstate( ( ielem-1)* nscalars+ q_se+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_SE,3) - sum9_1(3) + sum9_2(3)
      outstate( ( ielem-1)* nscalars+ q_nw+( 3-1)* qq )         &
        & = lamfac_o(3)*pdfTmp(Q_NW,3) + sum9_1(3) + sum9_2(3)

!    write(*,*) 'outpdf1', outstate(varPos(:,1))
!    write(*,*) 'outpdf2', outstate(varPos(:,2))
!    write(*,*) 'outpdf3', outstate(varPos(:,3))
    end do nodeloop

  end subroutine bgk_advRel_d3q19f3_MSGas
! ****************************************************************************** !


! ******************************************************************************
  !> Unoptimized Advection relaxation routine for the multispecies BGK model
  !! for testing
  !!
  !! This routine contains the implementation of semi-implicit lattice boltzmann
  !! equation using variable transformation based on the paper
  !! "Multi-species Lattice Boltzmann Model and Practical Examples. Short Course
  !! material Pietro Asinari PhD." \n
  !! Refer page: [Multispecies](../page/features/multispecies.html) for more information
  !! In the variable tranformation steps, we can skip the step 1 and step 3
  !! and evaluate only step 2 based on tranformed variable g
  !! only prerequisite is to compute feq which depends on original f not on g.
  !! feq is depend on density and velocity. Where density can be computed
  !! directly from g and velocity computed from linear system of equation
  !! given in the reference page [Multispecies](../page/features/multispecies.html).
  !! KM: This is an non-optimized kernel
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_advRel_MSGas_generic( fieldProp, inState, outState, auxField, &
    &                                  neigh, nElems, nSolve, level, layout,   &
    &                                  params, varSys, derVarPos )
    ! -------------------------------------------------------------------- !
    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> current layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> number of elements in state Array
    integer, intent(in) :: nElems
    !> input  pdf vector
    real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
    !> output pdf vector
    real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
    !> Auxiliary field computed from pre-collision state
    !! Is updated with correct velocity field for multicomponent models
    real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
    !> connectivity vector
    integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
    !> number of elements solved in kernel
    integer, intent(in) :: nSolve
    !> current level
    integer,intent(in) :: level
    !> global parameters
    type(mus_param_type),intent(in) :: params
    !> position of derived quantities in varsys for all fields
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: pdfTmp( varSys%nScalars )  ! temporary local pdf values
    integer :: iElem, s, vs, nFields, iDir, QQ, iFld
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: vPos( layout%fStencil%QQ )
    integer :: stateVarMap(varSys%nStateVars)
    real(kind=rk), dimension(varSys%nStateVars) :: rsigma, psigma,   &
      &                             uxsigma, uysigma, uzsigma,       &
      &                             qxsigma, qysigma, qzsigma,       &
      &                             gqxsigma, gqysigma, gqzsigma
    real(kind=rk), dimension(varSys%nStateVars) :: uxstar, uystar, uzstar
    real(kind=rk), dimension(varSys%nStateVars) :: chisigma, lambda
    real(kind=rk), dimension(varSys%nStateVars) :: m_sigma, m_rsig
    real(kind=rk) :: r, p, mm, temp, usqr, ucx, feq, theta
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) &
      & :: A, B, chi, Ainv
    integer :: dens_pos(varSys%nStateVars), mom_pos(3,varSys%nStateVars)
    integer :: elemOff
    ! ---------------------------------------------------------------------------
    ! Initialize variables
    feq = 0._rk

    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    !semi-implicit lbm is recovered for theta = 0.5
    theta = 0.5_rk

    !molecular weights
    m_sigma = fieldProp(:)%species%molweight
    !molecular weight ratios
    m_rsig = fieldProp(:)%species%molWeigRatio
    !resistivities
    do s = 1, nFields
      B(s,:) = fieldProp( s )%species%resi_coeff
    enddo
!    write(dbgUnit(10),*) 'B 1', B(1,:)
!    write(dbgUnit(10),*) 'B 2', B(2,:)
!    write(dbgUnit(10),*) 'B 3', B(3,:)
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    ! KM: debugging multispecies
    !outstate = instate

    ! position of density and velocity in auxField array
    do iFld = 1, nFields
      dens_pos(iFld) = varSys%method%val(derVarPos(iFld)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:,iFld) = varSys%method%val(derVarPos(iFld)%momentum) &
        &                           %auxField_varPos(:)
    end do

    nodeloop: do iElem = 1, nSolve
!    write(dbgUnit(10),*) 'iElem', iElem

      rsigma = 0._rk
      uxsigma = 0._rk
      uysigma = 0._rk
      uzsigma = 0._rk
      gqxsigma = 0._rk
      gqysigma = 0._rk
      gqzsigma = 0._rk
      pdfTmp = 0._rk

      ! element offset for auxField
      elemOff = (iElem-1)*varSys%nAuxScalars

      do s = 1, nFields
        vPos = varSys%method%val( stateVarMap(s) )%state_varPos
        do iDir = 1, QQ
          !store all field pdf in single element to local array pdfTmp
          pdfTmp( vPos(iDir) ) =                                               &
            & instate(  neigh (( idir-1)* nelems+ ielem)+( s-1)* qq+ varsys%nscalars*0 )
        enddo
        !field density
        rsigma( s ) = auxField(elemOff + dens_pos(s))

        !field momentum (rho*u)
        gqxsigma( s ) = auxField(elemOff + mom_pos(1, s))
        gqysigma( s ) = auxField(elemOff + mom_pos(2, s))
        gqzsigma( s ) = auxField(elemOff + mom_pos(3, s))

!          write(dbgunit,*) 's',s,'iDir', iDir, layout%fStencil%cxDir(3,iDir), pdfTmp( vPos(iDir) ), gqzsigma(s)

        !partial pressure
        psigma(s) = rsigma(s)*m_rsig(s)/3.0_rk

      enddo
!      write(dgbUnit(10),*) 'gqxsigma', gqxsigma
!      write(dgbUnit(10),*) 'gqysigma', gqysigma
!      write(dgbUnit(10),*) 'gqzsigma', gqzsigma
!      write(dgbUnit(10),*) 'rsigma', rsigma
!      write(dgbUnit(10),*) 'psigma', psigma

      !total density and total pressure
      r = sum(rsigma)
      p = sum(psigma)

!      write(*,*) 'r', r
!      write(*,*) 'p', p
      !Mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      mm = 0._rk
      do s = 1, nFields
        mm = mm + (rsigma(s)/(r*m_sigma(s)))
      enddo
      mm = 1._rk/mm

!      write(*,*) 'mm', mm
      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      do s = 1, nFields
        do vs = 1, nFields
        ! \todo is B(vs,s) = B(s,vs)
          chi(s,vs) = (mm*mm/(m_sigma(s)*m_sigma(vs)))*(B(vs,s)/B(s,s))
        enddo
      enddo

!      write(dgbUnit(10),*) 'chi', chi
      !relaxation time
      do s = 1, nFields
        lambda(s) = p*B(s,s)/r
      enddo
!      write(dgbUnit(10),*) 'lambda', lambda

      !chisigma - summation in the lse for velocity
      chisigma = 0._rk
      do s = 1, nFields
        temp = 0._rk
        do vs = 1, nFields
          temp = temp + chi(s,vs) * ( rsigma(vs) / r )
        enddo
        chisigma(s) = temp
      enddo
!      write(dgbUnit(10),*) 'chisigma', chisigma

      !matrix A of the velocity lse
      do s = 1, nFields
        do vs = 1, nFields
          if ( s == vs ) then
            A(s,vs) = 1._rk + theta * lambda(s) * chisigma(s)
          else
            A(s,vs) = 0._rk
          end if
          A(s,vs) = A(s,vs) - theta * lambda(s) * ( rsigma(s) / r ) * chi(s,vs)
        enddo
      enddo

      qxsigma = 0._rk
      qysigma = 0._rk
      qzsigma = 0._rk
!      write(dgbUnit(10),*) 'A 1', A(1,:)
!      write(dgbUnit(10),*) 'A 2', A(2,:)
!      write(dgbUnit(10),*) 'A 3', A(3,:)
!      write(*,*) 'gqxsigma',gqxsigma
      Ainv = invert_matrix( A )! , scheme%nFields )
      qxsigma = matmul( Ainv, gqxsigma )
      qysigma = matmul( Ainv, gqysigma )
      qzsigma = matmul( Ainv, gqzsigma )
!      call dbg_SolveLSE(nFields,A,gqxsigma,qxsigma)
!      call dbg_SolveLSE(nFields,A,gqysigma,qysigma)
!      call dbg_SolveLSE(nFields,A,gqzsigma,qzsigma)
!      write(dgbUnit(10),*) 'A 1', A(1,:)
!      write(dgbUnit(10),*) 'A 2', A(2,:)
!      write(dgbUnit(10),*) 'A 3', A(3,:)
!      write(dgbUnit(10),*) 'qxsigma',qxsigma
!      write(dgbUnit(10),*) 'qysigma',qysigma
!      write(dgbUnit(10),*) 'qzsigma',qzsigma

      uxsigma = qxsigma/rsigma
      uysigma = qysigma/rsigma
      uzsigma = qzsigma/rsigma

      ! store momentum of untransformed PDF in auxField
      do s = 1, nFields
        auxField(elemOff + mom_pos(1,s)) = qxsigma(s)
        auxField(elemOff + mom_pos(2,s)) = qysigma(s)
        auxField(elemOff + mom_pos(3,s)) = qzsigma(s)
      end do

      !compute ustar to compute feq
      uxstar = 0._rk
      uystar = 0._rk
      uzstar = 0._rk
!      write(dgbUnit(10),*) 'uxsigma', uxsigma
!      write(dgbUnit(10),*) 'uysigma', uysigma
!      write(dgbUnit(10),*) 'uzsigma', uzsigma

      do s = 1, nFields
        uxstar(s) = uxsigma(s)
        uystar(s) = uysigma(s)
        uzstar(s) = uzsigma(s)
        do vs = 1, nFields
          uxstar(s) = uxstar(s) + chi(s,vs) * ( rsigma(vs) / r ) *             &
            & ( uxsigma(vs) - uxsigma(s) )
          uystar(s) = uystar(s) + chi(s,vs) * ( rsigma(vs) / r ) *             &
            & ( uysigma(vs) - uysigma(s) )
          uzstar(s) = uzstar(s) + chi(s,vs) * ( rsigma(vs) / r ) *             &
            & ( uzsigma(vs) - uzsigma(s) )
        enddo
      enddo

!      write(dgbUnit(10),*) 'uxstar', uxstar
!      write(dgbUnit(10),*) 'uystar', uystar
!      write(dgbUnit(10),*) 'uzstar', uzstar
      !compute equilibrium and do collision
      do s = 1, nFields
        vPos = varSys%method%val( stateVarMap(s) )%state_varPos
        usqr = 0._rk
        usqr = uxstar(s) * uxstar(s) + uystar(s) * uystar(s) + uzstar(s) * uzstar(s)
        usqr = usqr * t2cs2inv
        do iDir = 1, layout%fStencil%QQ
          ucx = dble( layout%fStencil%cxDir( 1, iDir ))*uxstar(s)            &
            & + dble( layout%fStencil%cxDir( 2, iDir ))*uystar(s)            &
            & + dble( layout%fStencil%cxDir( 3, iDir ))*uzstar(s)
          if ( iDir == layout%fStencil%QQ ) then
          ! equilibrium at rest
            select case( trim(scheme%header%layout) )
            case('d2q9')
              feq = layout%weight( iDir ) * rsigma(s) * (                      &
                & ( 9._rk - 5._rk * M_rsig(s) )/4._rk + ucx * cs2inv           &
                & + ucx * ucx * t2cs4inv - usqr )
            case('d3q19')
              feq = layout%weight( iDir ) * rsigma(s) * (                      &
                & ( 3._rk - 2._rk * M_rsig(s) ) + ucx * cs2inv                 &
                & + ucx * ucx * t2cs4inv - usqr )
            case default
              write(logUnit(1),*)'ERROR: Only cases d2q9 and d3q19 are defined '
              write(logUnit(1),*)'Aborting'
              call tem_abort()
            end select
          else
            feq = layout%weight( iDir ) * rsigma(s) *                          &
              & ( M_rsig(s) + ucx * cs2inv + ucx * ucx * t2cs4inv - usqr )
          end if

          outstate( ( ielem-1)* varsys%nscalars+ idir+( s-1)* qq )&
            & = pdfTmp( vPos(iDir) )                                           &
            & + ( lambda(s) /( 1.0_rk + theta * lambda(s) ) )                  &
            & * ( feq - pdfTmp( vPos(iDir) ) )
        enddo

!      write(*,*) 'outpdf',s, outstate(vPos(:))
      enddo
    enddo nodeloop

  end subroutine bgk_advRel_MSGas_generic
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Solve linear system of equation with gauss elimination
  !!This code is taken from [Asinari code MIXLBM.f90]
  !!(http://staff.polito.it/pietro.asinari/rome08/index.html)
  !!
!   subroutine dbg_SolveLSE(n,A,b,x)
!     ! ---------------------------------------------------------------------------
!     implicit none
!
!     integer, intent(in) :: n
!     real(kind=rk),dimension(1:n,1:n),intent(in) :: A
!     real(kind=rk),dimension(1:n),intent(inout) :: b
!     real(kind=rk),dimension(1:n),intent(out) :: x
!     ! ---------------------------------------------------------------------------
!     real(kind=rk), dimension (n) :: s
!     integer, dimension(n) :: l
!     real(kind=rk),dimension(1:n,1:n) :: A_loc
!     ! ---------------------------------------------------------------------------
!
!     A_loc = A
!     ! LU Factorization
!     call gauss(n,A_loc,l,s)
!     ! Solve linear system of equation
!     call solve(n,A_loc,l,b,x)
!
!   end subroutine dbg_SolveLSE
! ****************************************************************************** !


!! ****************************************************************************** !
!  !  SUBROUTINE: Gauss
!  !
!  !  PURPOSE:  Fattorizzazione della matrice mediante metodo di eliminazione
!  !      di Gauss
!  subroutine gauss(n,a,l,s)
!    ! ---------------------------------------------------------------------------
!    implicit none
!
!    integer, intent(in) :: n
!    real(kind=rk), dimension(1:n,1:n), intent(inout) :: a
!    real(kind=rk), dimension(n), intent(out) :: s
!    integer, dimension(n), intent(out) :: l
!    ! ---------------------------------------------------------------------------
!    integer :: i,j,k
!    real(kind=rk) :: smax, rmax, xmult, r
!    integer :: lk
!    ! ---------------------------------------------------------------------------
!
!    do  i = 1,n
!      l(i) = i
!      smax = 0.0
!      do j = 1,n
!        smax = dmax1(smax,abs(a(i,j)))
!      end do
!      s(i) = smax
!    end do
!
!    do  k = 1,n-1
!      rmax = -1.0
!      do i = k,n
!        r = abs(a(l(i),k))/s(l(i))
!        if(r <= rmax)  exit
!        j = i
!        rmax = r
!      end do
!      lk = int(l(j) )
!      l(j) = l(k)
!      l(k) = lk
!      do i = k+1,n
!        xmult = a(l(i),k)/a(lk,k)
!        do j = k+1,n
!          a(l(i),j) = a(l(i),j) - xmult*a(lk,j)
!        end do
!        a(l(i),k) = xmult
!      end do
!    end do
!
!  end subroutine gauss
!! ****************************************************************************** !
!
!
!! ****************************************************************************** !
!  !  SUBROUTINE: Solve
!  !
!  !  PURPOSE:  Risoluzione del sistema lineare
!
!  !  Ingressi / uscite
!  subroutine solve(n,a,l,b,x)
!    ! ---------------------------------------------------------------------------
!    implicit none
!
!    integer, intent(in) :: n
!    real(kind=rk), dimension (1:n,1:n), intent(in) :: a
!    integer, dimension(n), intent(in) :: l
!    real(kind=rk), dimension (1:n), intent(inout) :: b
!    real(kind=rk), dimension (n), intent(out) :: x
!    ! ---------------------------------------------------------------------------
!    integer :: k, i, j
!    real(kind=rk) :: sum
!    ! ---------------------------------------------------------------------------
!
!    do k = 1,n-1
!      do i = k+1,n
!        b(l(i)) = b(l(i)) - a(l(i),k)*b(l(k))
!      end do
!    end do
!    x(n) = b(l(n))/a(l(n),n)
!
!    do i = n-1,1,-1
!      sum = b(l(i))
!      do j = i+1,n
!        sum = sum - a(l(i),j)*x(j)
!      end do
!      x(i) = sum/a(l(i),i)
!    end do
!
!  end subroutine solve
!! ****************************************************************************** !


end module mus_MSGas_module
! ****************************************************************************** !

!>\page multispecies Multispecies
!! Multispecies approach implemented in the code is based on the paper
!! "Multi-species Lattice Boltzmann Model and Practical Examples. Short Course
!! material Pietro Asinari PhD."\n
!!
!!Equlibrium distribution function for multispecies is given in the paper as
!!i\f$  f^{\sigma(eq)}_\alpha(\rho,\mathbf{u^*_{\sigma}}) =
!! \omega_\alpha\cdot\Big[s^{\sigma}_{\alpha}+\frac{1}{c^2_s}
!! (\mathbf{e}_\alpha \cdot\mathbf{u^*_{\sigma}})+\frac{1}{2c^4_s}(\mathbf{e}_\alpha
!! \cdot\mathbf{u^*_{\sigma}})^2-\frac{1}{c^2_s}(\mathbf{u^*_{\sigma}}\cdot
!! \mathbf{u^*_{\sigma}})\Big] \f$ \n
!!
!! where,
!! \f$s^{\sigma}_0 = (9-5 \phi^\sigma)/4 \f$\n
!! \f$s^{\sigma}_{\alpha} = \phi^\sigma\f$ for \f$1\leq\alpha\leq8\f$ and
!! \f$\phi^\sigma=min_\varsigma(m^\varsigma)/m^\sigma\f$,
!! \f$ p_\sigma=\rho_\sigma\phi^\sigma/3 \f$
!! \f$m^\sigma\f$ is molecular weight for species \f$\sigma\f$.\n
!!
!! \f$u^*_{\sigma} \f$ is given as \n
!! \f$\mathbf{u}^*_{\sigma} = \mathbf{u}_\sigma + \sum_{\varsigma}
!! {\frac{m^2}{m^\sigma m^{\varsigma}}\frac{B_{\sigma \varsigma}}{B_{mm}}
!! x_\varsigma (\mathbf{u}_\varsigma-\mathbf{u}_\sigma) }  \f$ \n
!! \f$x_\sigma= \rho_\sigma/\rho\f$\n
!! Relaxation time is given as
!! \f$ \lambda_\sigma = \frac{pB_{mm}}{\rho } \f$
!!
!!\section variableTransformation Multispecies: Variable Transformation
!!
!! A semi-implicit lattice Boltzmann equation is given as,
!! \f$ f^{\sigma,+}(\mathbf{x}+\mathbf{e},t+1) = f^\sigma(\mathbf{x},t) + (1-\frac{1}{2})
!! \lambda^\sigma[f^{\sigma(eq)}-f^\sigma]
!! + \frac{1}{2} \lambda{^\sigma,+}[f^{\sigma(eq),+}-f^{\sigma,+}]  \f$\n
!! Variable transformation presented in above paper involves three steps
!! \li Step 1. Transforming f -> g
!! \f$ g^\sigma = f^\sigma-\frac{1}{2}\lambda^\sigma[f^{\sigma(eq)}-f^\sigma] \f$\n
!! \li Step 2. Stream and collide i.e g -> g^+
!! \f$ g^{\sigma,+} = g^\sigma + \frac{\lambda^\sigma}{1+\frac{1}{2}\lambda^\sigma}
!! [f^{\sigma(eq)}-g^\sigma] \f$ \n
!! \li Step 3. Back transormation to f i.e g^+ -> f^+
!! \f$ f^{\sigma,+} = \frac{g^{\sigma,+} + \frac{1}{2}\lambda^{\sigma,+}
!! f^{\sigma(eq),+}}{1+\frac{1}{2}\lambda^{\sigma,+}} \f$ \n
!! In back tranformation, to compute feq we need \f$\rho\f$ and \f$ \mathbf{u}_\sigma \f$
!! $\f$ \rho \f$ can be computed directly from g
!! \f$  \rho^+_\sigma = <g^{\sigma,+}> \f$ \n
!! where as the \f$ \mathbf{u}_\sigma \f$ computed by solving the linear
!! system of equation given below \n
!! \f$  <\mathbf{e}_i g^{\sigma,+}> = \Bigg[ 1 + \frac{1}{2} \lambda^{\sigma,+} \sum_\varsigma
!! (\chi_{\sigma \varsigma} x^+_\varsigma) \Bigg] q^+_{\sigma,i} - \frac{1}{2}
!! \lambda^{\sigma,+}x^+_\sigma \sum_\varsigma (\chi_{\sigma \varsigma}q^+_{\varsigma,i})
!!\f$ \n
!! where, \f$ \chi_{\sigma \varsigma} \f$ is,
!! \f$ \chi_{\sigma \varsigma} = \frac{m^2}{m_\sigma m_\varsigma}
!! \frac{B_{\sigma \varsigma}}{B_{\sigma \sigma}} \f$ \n
!\todo use for optimized routine
!        !West
!        u_n(1) = - uxstar(s)
!        !south
!        u_n(2) =             - uystar(s)
!        !bottom
!        u_n(3) =                         - uzstar(s)
!        !east
!        u_n(4) =   uxstar(s)
!        !north
!        u_n(5) =               uystar(s)
!        !top
!        u_n(6) =                           uzstar(s)
!        !bottom south
!        u_n(7)
!        !top south
!        u_n(7)
!        !bottom north
!        u_n(7)
!        !top north
!        u_n(7)
!        !bottom west
!        u_n(7)
!        !bo
!        u_n(7)
!        u_n(7)
!        u_n(7)
!        u_n(7)
!

