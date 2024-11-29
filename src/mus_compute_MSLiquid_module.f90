! Copyright (c) 2013-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> author: Kannan Masilamani
!! This module provides the definition and methods for MSLiquid
!! bgk advection relaxation scheme.
module mus_MSLiquid_module
  use iso_c_binding,            only: c_f_pointer
  use env_module,               only: rk, long_k, eps
  use tem_varSys_module,        only: tem_varSys_type, tem_varSys_op_type
  use tem_param_module,         only: div1_3, div1_36, div1_8, div3_4h,        &
    &                                 div1_4, div3_8, div3_2, div9_16, div3_16,&
    &                                 cs2inv, div1_6, cs4inv, cs2,             &
    &                                 t2cs2inv, t2cs4inv, div1_18
  use tem_math_module,          only: invert_matrix
  use tem_isNaN_module,         only: tem_isNaN
  use tem_debug_module,         only: dbgUnit
  use tem_spacetime_fun_module, only: tem_spacetime_for
  use tem_aux_module,           only: tem_abort

  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_param_module,         only: mus_param_type
  use mus_mixture_module,       only: mus_mixture_type
  use mus_eNRTL_module,         only: mus_calc_thermFactor,                    &
    &                                 mus_calc_MS_DiffMatrix
  use mus_varSys_module,        only: mus_varSys_data_type
  use mus_derVarPos_module,     only: mus_derVarPos_type

  implicit none

  private

  public :: mrt_advRel_MSLiquid_generic
  public :: mrt_advRel_MSLiquid_generic_WTDF
  public :: mrt_advRel_d3q19f3_MSLiquid
  public :: mrt_advRel_d3q19f3_MSLiquid_WTDF
  public :: bgk_advRel_MSLiquid_generic
  public :: bgk_advRel_MSLiquid_generic_WTDF
  public :: bgk_advRel_d3q19f3_MSLiquid
  public :: bgk_advRel_d3q19f3_MSLiquid_WTDF
  public :: bgk_forcing_advRel_MSLiquid_generic

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

  integer, parameter :: QQF3  = 57   !< number of pdf directions for 3 species

  !TG! integer,dimension(QQ,3), parameter :: varPos = &
  !TG!  & reshape(&
  !TG!  &  [q__w, q__s, q__b, q__e, q__n, q__t, q_bs, q_ts, q_bn, &
  !TG!  &   q_tn, q_bw, q_be, q_tw, q_te, q_sw, q_nw, q_se, q_ne, q__0, &
  !TG!  &   q__w+19, q__s+19, q__b+19, q__e+19, q__n+19, q__t+19, &
  !TG!  &   q_bs+19, q_ts+19, q_bn+19, q_tn+19, q_bw+19, q_be+19, &
  !TG!  &   q_tw+19, q_te+19, q_sw+19, q_nw+19, q_se+19, q_ne+19, q__0+19,&
  !TG!  &   q__w+38, q__s+38, q__b+38, q__e+38, q__n+38, q__t+38, &
  !TG!  &   q_bs+38, q_ts+38, q_bn+38, q_tn+38, q_bw+38, q_be+38, &
  !TG!  &   q_tw+38, q_te+38, q_sw+38, q_nw+38, q_se+38, q_ne+38, q__0+38] &
  !TG!  &  ,[QQ,3])
contains

! ******************************************************************************
  !> Optimized Advection relaxation routine for the MSLiquid BGK model
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
  !! given in the reference page [Multispecies](../page/features/multispecies.html).
  !! KM: This is an non-optimized kernel
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_advRel_d3q19f3_MSLiquid( fieldProp, inState, outState,    &
    &                                     auxField, neigh, nElems, nSolve, &
    &                                     level, layout, params, varSys,   &
    &                                     derVarPos                        )
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
    !local variables
    ! temporary local pdf values
    real(kind=rk) :: pdfTmp( QQ, 3 )
    integer :: iElem, nScalars, iFld
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk), dimension(3) :: rsigma, uxsigma, rsigmaInv,                 &
      & uysigma, uzsigma, qxsigma, qysigma, qzsigma, gqxsigma, gqysigma,       &
      & gqzsigma
    real(kind=rk), dimension(3) :: uxstar, uystar, uzstar
    real(kind=rk), dimension(3) :: molWeight, molWeight_inv, phi,              &
      &                            num_dens, moleFrac
    real(kind=rk) :: usqr(3), totNum_dens_inv, totMass_densInv
    real(kind=rk) :: velAvg(3)
    real(kind=rk) :: velQuadTerm_x(3), velQuadTerm_y(3), velQuadTerm_z(3)
    real(kind=rk) :: omega, omegadiv2, omega_fac, omega_o
    real(kind=rk), dimension(3,3) :: A, mbb
    real(kind=rk) :: zfac, yfac, a12_11, a13_11, &
      & a23fac, a21_11, a31_11, zcoeff
    real(kind=rk) :: theta_eq, theta_eq_spc, paramB_inv, a11_inv, a22fac_inv
    real(kind=rk) :: B_12, B_13, B_23
    real(kind=rk) :: Bratio_12, Bratio_13, Bratio_23
    real(kind=rk) :: mbbEq_12, mbbEq_13, mbbEq_21, mbbEq_23, mbbEq_31, mbbEq_32
    real(kind=rk) :: chi12, chi13, chi21, chi23, chi31, chi32
    real(kind=rk),dimension(3) :: omegaRho, omegaRho_d, omegaRho_dm, &
      & omegaRho_o, omegaRho_om, sum1_1, sum2_1, sum3_1,                       &
      & sum4_1, sum5_1, sum6_1, sum7_1, sum8_1, sum9_1, sum1_2, sum2_2, sum3_2,&
      & sum4_2, sum5_2, sum6_2, sum7_2, sum8_2, sum9_2
    integer :: dens_pos(3), mom_pos(3,3), elemOff
    ! ---------------------------------------------------------------------------

    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nScalars = varSys%nScalars
    !molecular weights
    molWeight = fieldProp(:)%species%molweight
    molWeight_inv = 1.0_rk/molWeight
    !molecular weight ratios
    phi = fieldProp(:)%species%molWeigRatio

    !omega
    omega = scheme%mixture%relaxLvl(Level)%omega_diff
    omega_fac = (1.0_rk/(1.0_rk/omega + 0.5_rk))
    omegadiv2 = omega*0.5_rk
    omega_o = 1.0_rk - omega_fac

    !free parameter B
    paramB_inv = 1.0_rk/scheme%mixture%paramB

    !equilibrium theta
    theta_eq = scheme%mixture%theta_eq
    theta_eq_spc = 1.0_rk - theta_eq

    !resistivities
    !resistivity coeff are symmetric
    !B_12 = B_21, B_13=B31, B23=B32
    B_12 = fieldProp( 1 )%species%resi_coeff(2)
    B_13 = fieldProp( 1 )%species%resi_coeff(3)
    B_23 = fieldProp( 2 )%species%resi_coeff(3)

    Bratio_12 = B_12*paramB_inv
    Bratio_13 = B_13*paramB_inv
    Bratio_23 = B_23*paramB_inv

    mbbEq_12 = Bratio_12*phi(1)
    mbbEq_13 = Bratio_13*phi(1)

    mbbEq_21 = Bratio_12*phi(2)
    mbbEq_23 = Bratio_23*phi(2)

    mbbEq_31 = Bratio_13*phi(3)
    mbbEq_32 = Bratio_23*phi(3)

    ! factor need to build matrixA to solve LSE
    ! dummy diagonal entries
    mbb(1,1) = -1.0_rk
    mbb(1,2) = mbbEq_21*omegadiv2
    mbb(1,3) = mbbEq_31*omegadiv2

    mbb(2,1) = mbbEq_12*omegadiv2
    mbb(2,2) = -1.0_rk
    mbb(2,3) = mbbEq_32*omegadiv2

    mbb(3,1) = mbbEq_13*omegadiv2
    mbb(3,2) = mbbEq_23*omegadiv2
    mbb(3,3) = -1.0_rk

    ! position of density and momentum in auxField array
    do iFld = 1, scheme%nFields
      dens_pos(iFld) = varSys%method%val(derVarPos(iFld)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:, iFld) = varSys%method%val(derVarPos(iFld)%momentum) &
        &                                %auxField_varPos(1:3)
    end do

    nodeloop: do iElem = 1, nSolve

      !species 1
      pdfTmp( Q__W,1 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__S,1 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__B,1 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__E,1 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__N,1 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__T,1 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,1 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,1 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,1 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,1 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,1 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,1 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,1 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,1 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,1 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,1 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,1 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,1 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__0,1 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      !species 2
      pdfTmp( Q__W,2 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__S,2 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__B,2 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__E,2 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__N,2 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__T,2 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,2 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,2 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,2 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,2 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,2 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,2 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,2 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,2 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,2 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,2 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,2 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,2 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__0,2 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)

      !species 3
      pdfTmp( Q__W,3 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__S,3 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__B,3 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__E,3 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__N,3 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__T,3 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,3 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,3 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,3 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,3 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,3 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,3 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,3 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,3 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,3 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,3 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,3 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,3 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__0,3 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)

      ! element offset to access auxField
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! get density and momentum from auxField
      do iFld = 1, 3
        ! density
        rsigma(iFld) = auxField(elemOff + dens_pos(iFld))
        ! momentum
        gqxsigma(iFld) = auxField(elemOff + mom_pos(1, iFld))
        gqysigma(iFld) = auxField(elemOff + mom_pos(2, iFld))
        gqzsigma(iFld) = auxField(elemOff + mom_pos(3, iFld))
      end do

      !number density
      num_dens(:) = rsigma(:)*molWeight_inv(:)

      !total number density
      totNum_dens_inv = 1.0_rk/(num_dens(1)+num_dens(2)+num_dens(3))

      !mole fraction
      moleFrac(:) = num_dens(:)*totNum_dens_inv

      !matrix A of the velocity lse
      A(1,1) = 1.0_rk + moleFrac(2)*mbb(2,1) + moleFrac(3)*mbb(3,1)
      A(1,2) = - moleFrac(1)*mbb(1,2)
      A(1,3) = - moleFrac(1)*mbb(1,3)

      A(2,1) = - moleFrac(2)*mbb(2,1)
      A(2,2) = 1.0_rk + moleFrac(1)*mbb(1,2) + moleFrac(3)*mbb(3,2)
      A(2,3) = - moleFrac(2)*mbb(2,3)

      A(3,1) = - moleFrac(3)*mbb(3,1)
      A(3,2) = - moleFrac(3)*mbb(3,2)
      A(3,3) = 1.0_rk + moleFrac(1)*mbb(1,3) + moleFrac(2)*mbb(2,3)

      a11_inv = 1.0_rk/A(1,1)
      a13_11 = A(1,3)*a11_inv
      a12_11 = A(1,2)*a11_inv
      a31_11 = A(3,1)*a11_inv
      a21_11 = A(2,1)*a11_inv
      a23fac = A(2,3) - A(2,1)*a13_11
      a22fac_inv = 1.0_rk/(A(2,2) - A(2,1)*a12_11)

      zfac = (A(3,1)*a12_11 - A(3,2)) * a22fac_inv
      zcoeff = 1._rk / (a23fac*zfac - A(3,1)*a13_11 + A(3,3))
      !x-momentum
      !species 3
      yfac = gqxsigma(2) - gqxsigma(1)*a21_11

      qxsigma(3) = (gqxsigma(3) - a31_11*gqxsigma(1) + yfac*zfac)*zcoeff

      !species 2
      qxsigma(2) = (yfac - a23fac*qxsigma(3))*a22fac_inv
      !species 1
      qxsigma(1) = (gqxsigma(1) - A(1,2)*qxsigma(2) - A(1,3)*qxsigma(3))*a11_inv
      !y-momentum
      yfac = gqysigma(2) - gqysigma(1)*a21_11
      qysigma(3) = (gqysigma(3) - a31_11*gqysigma(1) + yfac*zfac)*zcoeff

      qysigma(2) = (yfac - a23fac*qysigma(3))*a22fac_inv

      qysigma(1) = (gqysigma(1) - A(1,2)*qysigma(2) - A(1,3)*qysigma(3))*a11_inv
      !z-momentum
      yfac = gqzsigma(2) - gqzsigma(1)*a21_11
      qzsigma(3) = (gqzsigma(3) - a31_11*gqzsigma(1) + yfac*zfac)*zcoeff

      qzsigma(2) = (yfac - a23fac*qzsigma(3))*a22fac_inv

      qzsigma(1) = (gqzsigma(1) - A(1,2)*qzsigma(2) - A(1,3)*qzsigma(3))*a11_inv

      ! store momentum of untransformed PDF in auxField
      do iFld = 1, 3
        auxField(elemOff + mom_pos(1, iFld)) = qxsigma(iFld)
        auxField(elemOff + mom_pos(2, iFld)) = qysigma(iFld)
        auxField(elemOff + mom_pos(3, iFld)) = qzsigma(iFld)
      end do

      rsigmaInv(1) = 1.0_rk/rsigma(1)
      rsigmaInv(2) = 1.0_rk/rsigma(2)
      rsigmaInv(3) = 1.0_rk/rsigma(3)

      uxsigma(1) = qxsigma(1)*rsigmaInv(1)
      uxsigma(2) = qxsigma(2)*rsigmaInv(2)
      uxsigma(3) = qxsigma(3)*rsigmaInv(3)

      uysigma(1) = qysigma(1)*rsigmaInv(1)
      uysigma(2) = qysigma(2)*rsigmaInv(2)
      uysigma(3) = qysigma(3)*rsigmaInv(3)

      uzsigma(1) = qzsigma(1)*rsigmaInv(1)
      uzsigma(2) = qzsigma(2)*rsigmaInv(2)
      uzsigma(3) = qzsigma(3)*rsigmaInv(3)

!write(dbgUnit,*) 'velocity 1', uxsigma(1), uysigma(1), uzsigma(1)
!write(dbgUnit,*) 'velocity 2', uxsigma(2), uysigma(2), uzsigma(2)
!write(dbgUnit,*) 'velocity 3', uxsigma(3), uysigma(3), uzsigma(3)

      !compute equilibrium velocity to compute feq
      chi12 = mbbEq_12 * moleFrac(2)
      chi13 = mbbEq_13 * moleFrac(3)
      chi21 = mbbEq_21 * moleFrac(1)
      chi23 = mbbEq_23 * moleFrac(3)
      chi31 = mbbEq_31 * moleFrac(1)
      chi32 = mbbEq_32 * moleFrac(2)
      !ux
      uxstar(1) =  uxsigma(1) + chi12*(uxsigma(2)-uxsigma(1)) &
        &       + chi13*(uxsigma(3)-uxsigma(1))

      uxstar(2) =  uxsigma(2) + chi21*(uxsigma(1)-uxsigma(2)) &
        &       + chi23*(uxsigma(3)-uxsigma(2))

      uxstar(3) =  uxsigma(3) + chi31*(uxsigma(1)-uxsigma(3)) &
        &       + chi32*(uxsigma(2)-uxsigma(3))
      !uy
      uystar(1) =  uysigma(1) + chi12*(uysigma(2)-uysigma(1)) &
        &       + chi13*(uysigma(3)-uysigma(1))

      uystar(2) =  uysigma(2) + chi21*(uysigma(1)-uysigma(2)) &
        &       + chi23*(uysigma(3)-uysigma(2))

      uystar(3) =  uysigma(3) + chi31*(uysigma(1)-uysigma(3)) &
        &       + chi32*(uysigma(2)-uysigma(3))
      !uz
      uzstar(1) =  uzsigma(1) + chi12*(uzsigma(2)-uzsigma(1)) &
        &       + chi13*(uzsigma(3)-uzsigma(1))

      uzstar(2) =  uzsigma(2) + chi21*(uzsigma(1)-uzsigma(2)) &
        &       + chi23*(uzsigma(3)-uzsigma(2))

      uzstar(3) =  uzsigma(3) + chi31*(uzsigma(1)-uzsigma(3)) &
        &       + chi32*(uzsigma(2)-uzsigma(3))

!write(dbgUnit,*) 'eqVel 1', uxstar(1), uystar(1), uzstar(1)
!write(dbgUnit,*) 'eqVel 2', uxstar(2), uystar(2), uzstar(2)
!write(dbgUnit,*) 'eqVel 3', uxstar(3), uystar(3), uzstar(3)

      ! total mass density
      totMass_densInv = 1.0_rk/(rsigma(1)+rsigma(2)+rsigma(3))

      ! mass averaged mixture velocity three components
      velAvg(1) = ( rsigma(1) * uxsigma(1) &
        &         + rsigma(2) * uxsigma(2) &
        &         + rsigma(3) * uxsigma(3) ) * totMass_densInv * theta_eq

      velAvg(2) = ( rsigma(1) * uysigma(1) &
        &         + rsigma(2) * uysigma(2) &
        &         + rsigma(3) * uysigma(3) ) * totMass_densInv * theta_eq

      velAvg(3) = ( rsigma(1) * uzsigma(1) &
        &         + rsigma(2) * uzsigma(2) &
        &         + rsigma(3) * uzsigma(3) ) * totMass_densInv * theta_eq

      ! velocity in quadratic term of equilibrium
      velQuadTerm_x(1) = velAvg(1) + theta_eq_spc*uxstar(1)
      velQuadTerm_y(1) = velAvg(2) + theta_eq_spc*uystar(1)
      velQuadTerm_z(1) = velAvg(3) + theta_eq_spc*uzstar(1)

      !compute equilibrium and do collision
      usqr(1) = (velQuadTerm_x(1) * velQuadTerm_x(1)                           &
        &     + velQuadTerm_y(1) * velQuadTerm_y(1)                            &
        &     + velQuadTerm_z(1) * velQuadTerm_z(1)) * t2cs2inv

      omegaRho(:) = rsigma(:)*omega_fac
      omegaRho_d(:) = omegaRho(:) * div1_4

      ! species 1
      ! equilibrium at rest
      outstate( ( ielem-1)* nscalars+ q__0+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q__0,1) +                                       &
              & div1_3 * omegaRho(1) * (( 3._rk - 2._rk * phi(1)) - usqr(1) )

      omegaRho_dm(1) = omegaRho(1) * (phi(1) - usqr(1)) * div1_18

      !directional velocity factor
      sum1_1(1) = omegaRho_d(1) * uxstar(1) * div3_4h
      sum1_2(1) = omegaRho_d(1) * velQuadTerm_x(1)**2 + omegaRho_dm(1)

      outstate( ( ielem-1)* nscalars+ q__w+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__W,1) - sum1_1(1) + sum1_2(1)
      outstate( ( ielem-1)* nscalars+ q__e+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__E,1) + sum1_1(1) + sum1_2(1)

      sum2_1(1) = omegaRho_d(1) * uystar(1) * div3_4h
      sum2_2(1) = omegaRho_d(1) * velQuadTerm_y(1)**2 + omegaRho_dm(1)

      outstate( ( ielem-1)* nscalars+ q__s+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__S,1) - sum2_1(1) + sum2_2(1)
      outstate( ( ielem-1)* nscalars+ q__n+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__N,1) + sum2_1(1) + sum2_2(1)

      sum3_1(1) = omegaRho_d(1) * uzstar(1) * div3_4h
      sum3_2(1) = omegaRho_d(1) * velQuadTerm_z(1)**2 + omegaRho_dm(1)

      outstate( ( ielem-1)* nscalars+ q__b+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__B,1) - sum3_1(1) + sum3_2(1)
      outstate( ( ielem-1)* nscalars+ q__t+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__T,1) + sum3_1(1) + sum3_2(1)

      !top north / bottom south
      omegaRho_o(1) = omegaRho_d(1) * 0.5_rk
      omegaRho_om(1) = omegaRho_dm(1) * 0.5_rk
      sum4_1(1) = omegaRho_o(1) * ( uystar(1) + uzstar(1) ) *div3_4h
      sum4_2(1) = omegaRho_o(1) * ( velQuadTerm_y(1) + velQuadTerm_z(1) )**2 &
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_bs+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BS,1) - sum4_1(1) + sum4_2(1)
      outstate( ( ielem-1)* nscalars+ q_tn+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TN,1) + sum4_1(1) + sum4_2(1)

      !top south / bottom north
      sum5_1(1) = omegaRho_o(1) * ( -uystar(1) + uzstar(1) ) *div3_4h
      sum5_2(1) = omegaRho_o(1) * ( -velQuadTerm_y(1) + velQuadTerm_z(1) )**2&
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_bn+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BN,1) - sum5_1(1) + sum5_2(1)
      outstate( ( ielem-1)* nscalars+ q_ts+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TS,1) + sum5_1(1) + sum5_2(1)

      !top east / bottom west
      sum6_1(1) = omegaRho_o(1) * ( uxstar(1) + uzstar(1) ) *div3_4h
      sum6_2(1) = omegaRho_o(1) * ( velQuadTerm_x(1) + velQuadTerm_z(1) )**2 &
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_bw+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BW,1) - sum6_1(1) + sum6_2(1)
      outstate( ( ielem-1)* nscalars+ q_te+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TE,1) + sum6_1(1) + sum6_2(1)

      !top west / bottom east
      sum7_1(1) = omegaRho_o(1) * ( -uxstar(1) + uzstar(1) ) *div3_4h
      sum7_2(1) = omegaRho_o(1) * ( -velQuadTerm_x(1) + velQuadTerm_z(1) )**2&
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_be+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BE,1) - sum7_1(1) + sum7_2(1)
      outstate( ( ielem-1)* nscalars+ q_tw+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TW,1) + sum7_1(1) + sum7_2(1)

      !north east / south west
      sum8_1(1) = omegaRho_o(1) * ( uxstar(1) + uystar(1) ) *div3_4h
      sum8_2(1) = omegaRho_o(1) * ( velQuadTerm_x(1) + velQuadTerm_y(1) )**2 &
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_sw+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SW,1) - sum8_1(1) + sum8_2(1)
      outstate( ( ielem-1)* nscalars+ q_ne+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NE,1) + sum8_1(1) + sum8_2(1)

      !north west / south east
      sum9_1(1) = omegaRho_o(1) * ( -uxstar(1) + uystar(1) ) *div3_4h
      sum9_2(1) = omegaRho_o(1) * ( -velQuadTerm_x(1) + velQuadTerm_y(1) )**2&
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_se+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SE,1) - sum9_1(1) + sum9_2(1)
      outstate( ( ielem-1)* nscalars+ q_nw+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NW,1) + sum9_1(1) + sum9_2(1)

      ! species 2
      ! velocity in quadratic term of equilibrium
      velQuadTerm_x(2) = velAvg(1) + theta_eq_spc*uxstar(2)
      velQuadTerm_y(2) = velAvg(2) + theta_eq_spc*uystar(2)
      velQuadTerm_z(2) = velAvg(3) + theta_eq_spc*uzstar(2)

      usqr(2) = (velQuadTerm_x(2) * velQuadTerm_x(2)                           &
        &     + velQuadTerm_y(2) * velQuadTerm_y(2)                            &
        &     + velQuadTerm_z(2) * velQuadTerm_z(2)) * t2cs2inv

      outstate( ( ielem-1)* nscalars+ q__0+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q__0,2) + &
              & div1_3 * omegaRho(2) * (( 3._rk - 2._rk * phi(2)) - usqr(2) )

      omegaRho_dm(2) = omegaRho(2) * (phi(2) - usqr(2)) * div1_18
      !directional velocity factor
      sum1_1(2) = omegaRho_d(2) * uxstar(2) * div3_4h
      sum1_2(2) = omegaRho_d(2) * velQuadTerm_x(2)**2 + omegaRho_dm(2)

      outstate( ( ielem-1)* nscalars+ q__w+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__W,2) - sum1_1(2) + sum1_2(2)
      outstate( ( ielem-1)* nscalars+ q__e+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__E,2) + sum1_1(2) + sum1_2(2)

      sum2_1(2) = omegaRho_d(2) * uystar(2) * div3_4h
      sum2_2(2) = omegaRho_d(2) * velQuadTerm_y(2)**2 + omegaRho_dm(2)

      outstate( ( ielem-1)* nscalars+ q__s+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__S,2) - sum2_1(2) + sum2_2(2)
      outstate( ( ielem-1)* nscalars+ q__n+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__N,2) + sum2_1(2) + sum2_2(2)

      sum3_1(2) = omegaRho_d(2) * uzstar(2) * div3_4h
      sum3_2(2) = omegaRho_d(2) * velQuadTerm_z(2)**2 + omegaRho_dm(2)

      outstate( ( ielem-1)* nscalars+ q__b+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__B,2) - sum3_1(2) + sum3_2(2)
      outstate( ( ielem-1)* nscalars+ q__t+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__T,2) + sum3_1(2) + sum3_2(2)

      !top north / bottom south
      omegaRho_o(2) = omegaRho_d(2) * 0.5_rk
      omegaRho_om(2) = omegaRho_dm(2) * 0.5_rk
      sum4_1(2) = omegaRho_o(2) * ( uystar(2) + uzstar(2) ) *div3_4h
      sum4_2(2) = omegaRho_o(2) * ( velQuadTerm_y(2) + velQuadTerm_z(2) )**2 &
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_bs+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BS,2) - sum4_1(2) + sum4_2(2)
      outstate( ( ielem-1)* nscalars+ q_tn+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TN,2) + sum4_1(2) + sum4_2(2)

      !top south / bottom north
      sum5_1(2) = omegaRho_o(2) * ( -uystar(2) + uzstar(2) ) *div3_4h
      sum5_2(2) = omegaRho_o(2) * ( -velQuadTerm_y(2) + velQuadTerm_z(2) )**2&
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_bn+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BN,2) - sum5_1(2) + sum5_2(2)
      outstate( ( ielem-1)* nscalars+ q_ts+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TS,2) + sum5_1(2) + sum5_2(2)

      !top east / bottom west
      sum6_1(2) = omegaRho_o(2) * ( uxstar(2) + uzstar(2) ) *div3_4h
      sum6_2(2) = omegaRho_o(2) * ( velQuadTerm_x(2) + velQuadTerm_z(2) )**2 &
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_bw+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BW,2) - sum6_1(2) + sum6_2(2)
      outstate( ( ielem-1)* nscalars+ q_te+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TE,2) + sum6_1(2) + sum6_2(2)

      !top west / bottom east
      sum7_1(2) = omegaRho_o(2) * ( -uxstar(2) + uzstar(2) ) *div3_4h
      sum7_2(2) = omegaRho_o(2) * ( -velQuadTerm_x(2) + velQuadTerm_z(2) )**2&
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_be+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BE,2) - sum7_1(2) + sum7_2(2)
      outstate( ( ielem-1)* nscalars+ q_tw+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TW,2) + sum7_1(2) + sum7_2(2)

      !north east / south west
      sum8_1(2) = omegaRho_o(2) * ( uxstar(2) + uystar(2) ) *div3_4h
      sum8_2(2) = omegaRho_o(2) * ( velQuadTerm_x(2) + velQuadTerm_y(2) )**2 &
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_sw+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SW,2) - sum8_1(2) + sum8_2(2)
      outstate( ( ielem-1)* nscalars+ q_ne+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NE,2) + sum8_1(2) + sum8_2(2)

      !north west / south east
      sum9_1(2) = omegaRho_o(2) * ( -uxstar(2) + uystar(2) ) *div3_4h
      sum9_2(2) = omegaRho_o(2) * ( -velQuadTerm_x(2) + velQuadTerm_y(2) )**2&
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_se+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SE,2) - sum9_1(2) + sum9_2(2)
      outstate( ( ielem-1)* nscalars+ q_nw+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NW,2) + sum9_1(2) + sum9_2(2)

      ! species 3
      ! velocity in quadratic term of equilibrium
      velQuadTerm_x(3) = velAvg(1) + theta_eq_spc*uxstar(3)
      velQuadTerm_y(3) = velAvg(2) + theta_eq_spc*uystar(3)
      velQuadTerm_z(3) = velAvg(3) + theta_eq_spc*uzstar(3)

      usqr(3) = (velQuadTerm_x(3) * velQuadTerm_x(3)                           &
        &     + velQuadTerm_y(3) * velQuadTerm_y(3)                            &
        &     + velQuadTerm_z(3) * velQuadTerm_z(3)) * t2cs2inv

      outstate( ( ielem-1)* nscalars+ q__0+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q__0,3) + &
              & div1_3 * omegaRho(3) * (( 3._rk - 2._rk * phi(3)) - usqr(3) )

      omegaRho_dm(3) = omegaRho(3) * (phi(3) - usqr(3)) * div1_18
      !directional velocity factor
      sum1_1(3) = omegaRho_d(3) * uxstar(3) * div3_4h
      sum1_2(3) = omegaRho_d(3) * velQuadTerm_x(3)**2 + omegaRho_dm(3)

      outstate( ( ielem-1)* nscalars+ q__w+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__W,3) - sum1_1(3) + sum1_2(3)
      outstate( ( ielem-1)* nscalars+ q__e+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__E,3) + sum1_1(3) + sum1_2(3)

      sum2_1(3) = omegaRho_d(3) * uystar(3) * div3_4h
      sum2_2(3) = omegaRho_d(3) * velQuadTerm_y(3)**2 + omegaRho_dm(3)

      outstate( ( ielem-1)* nscalars+ q__s+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__S,3) - sum2_1(3) + sum2_2(3)
      outstate( ( ielem-1)* nscalars+ q__n+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__N,3) + sum2_1(3) + sum2_2(3)

      sum3_1(3) = omegaRho_d(3) * uzstar(3) * div3_4h
      sum3_2(3) = omegaRho_d(3) * velQuadTerm_z(3)**2 + omegaRho_dm(3)

      outstate( ( ielem-1)* nscalars+ q__b+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__B,3) - sum3_1(3) + sum3_2(3)
      outstate( ( ielem-1)* nscalars+ q__t+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__T,3) + sum3_1(3) + sum3_2(3)

      !top north / bottom south
      omegaRho_o(3) = omegaRho_d(3) * 0.5_rk
      omegaRho_om(3) = omegaRho_dm(3) * 0.5_rk
      sum4_1(3) = omegaRho_o(3) * ( uystar(3) + uzstar(3) ) *div3_4h
      sum4_2(3) = omegaRho_o(3) * ( velQuadTerm_y(3) + velQuadTerm_z(3) )**2 &
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_bs+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BS,3) - sum4_1(3) + sum4_2(3)
      outstate( ( ielem-1)* nscalars+ q_tn+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TN,3) + sum4_1(3) + sum4_2(3)

      !top south / bottom north
      sum5_1(3) = omegaRho_o(3) * ( -uystar(3) + uzstar(3) ) *div3_4h
      sum5_2(3) = omegaRho_o(3) * ( -velQuadTerm_y(3) + velQuadTerm_z(3) )**2&
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_bn+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BN,3) - sum5_1(3) + sum5_2(3)
      outstate( ( ielem-1)* nscalars+ q_ts+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TS,3) + sum5_1(3) + sum5_2(3)

      !top east / bottom west
      sum6_1(3) = omegaRho_o(3) * ( uxstar(3) + uzstar(3) ) *div3_4h
      sum6_2(3) = omegaRho_o(3) * ( velQuadTerm_x(3) + velQuadTerm_z(3) )**2 &
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_bw+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BW,3) - sum6_1(3) + sum6_2(3)
      outstate( ( ielem-1)* nscalars+ q_te+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TE,3) + sum6_1(3) + sum6_2(3)

      !top west / bottom east
      sum7_1(3) = omegaRho_o(3) * ( -uxstar(3) + uzstar(3) ) *div3_4h
      sum7_2(3) = omegaRho_o(3) * ( -velQuadTerm_x(3) + velQuadTerm_z(3) )**2&
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_be+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BE,3) - sum7_1(3) + sum7_2(3)
      outstate( ( ielem-1)* nscalars+ q_tw+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TW,3) + sum7_1(3) + sum7_2(3)

      !north east / south west
      sum8_1(3) = omegaRho_o(3) * ( uxstar(3) + uystar(3) ) *div3_4h
      sum8_2(3) = omegaRho_o(3) * ( velQuadTerm_x(3) + velQuadTerm_y(3) )**2 &
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_sw+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SW,3) - sum8_1(3) + sum8_2(3)
      outstate( ( ielem-1)* nscalars+ q_ne+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NE,3) + sum8_1(3) + sum8_2(3)

      !north west / south east
      sum9_1(3) = omegaRho_o(3) * ( -uxstar(3) + uystar(3) ) *div3_4h
      sum9_2(3) = omegaRho_o(3) * ( -velQuadTerm_x(3) + velQuadTerm_y(3) )**2&
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_se+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SE,3) - sum9_1(3) + sum9_2(3)
      outstate( ( ielem-1)* nscalars+ q_nw+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NW,3) + sum9_1(3) + sum9_2(3)

    end do nodeloop

  end subroutine bgk_advRel_d3q19f3_MSLiquid
! ******************************************************************************



! ******************************************************************************
  !> Unoptimized Advection relaxation routine for the multispecies BGK model
  !!
  !! This routine contains the implementation based on the paper
  !! "A Lattice Boltzmann Scheme for liquid mixtures - Part II: Discretization
  !! and Numerics, Jens Zudrop, Sabine Roller, Pietro Asinari. "\n
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_advRel_MSLiquid_generic( fieldProp, inState, outState,    &
    &                                     auxField, neigh, nElems, nSolve, &
    &                                     level, layout, params, varSys,   &
    &                                     derVarPos                        )
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
    !local variables
    real(kind=rk) :: pdfTmp(varSys%nScalars)  !< temporary local pdf values
    integer :: iElem, iField, iField_2, nFields, iDir, QQ, nScalars
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: vPos(layout%fStencil%QQ)
    integer, dimension(varSys%nStateVars) :: stateVarMap
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, num_dens
    real(kind=rk), dimension(varSys%nStateVars) :: moleFrac
    real(kind=rk), dimension(3, varSys%nStateVars) &
      & :: first_moments, velocity, eqVel
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk) :: totNum_dens, usqr, ucx, feq
    real(kind=rk) :: velAvg(3), velNew(3), ucxQuadTerm, theta_eq, totMassDens
    real(kind=rk) :: omega, omega_new, paramB
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) &
      & :: matA, invA, resi_coeff
    integer :: elemOff
    integer :: dens_pos(varSys%nStateVars)
    integer :: mom_pos(3, varSys%nStateVars)
    real(kind=rk) :: wRestInv, sRest(varSys%nStateVars)
    ! ---------------------------------------------------------------------------

    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    QQ     = layout%fStencil%QQ
    nFields = scheme%nFields
    nScalars = varSys%nScalars

    !molecular weights
    molWeight = fieldProp(:)%species%molweight
    !molecular weight ratios
    phi = fieldProp(:)%species%molWeigRatio
    !s_0^k - species weight for rest position
    wRestInv = 1.0_rk / layout%weight(layout%fStencil%restPosition)
    do iField = 1, nFields
      sRest(iField) = wRestInv + (1.0_rk - wRestInv) * phi(iField)
    end do

    !resistivities
    do iField = 1, nFields
      resi_coeff(iField,:) = fieldProp( iField )%species%resi_coeff
    enddo

    omega = scheme%mixture%relaxLvl(Level)%omega_diff
    omega_new = (1.0_rk/(1.0_rk/omega + 0.5_rk))

    paramB = scheme%mixture%paramB
    !equilibrium theta
    theta_eq = scheme%mixture%theta_eq

    stateVarMap = scheme%stateVarMap%varPos%val(:)

    ! position of density and momentum in auxField array
    do iField = 1, nFields
      dens_pos(iField) = varSys%method%val(derVarPos(iField)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:, iField) = varSys%method%val(derVarPos(iField)%momentum) &
        &                                  %auxField_varPos(1:3)
    end do

    nodeloop: do iElem = 1, nSolve

      mass_dens = 0._rk
      first_moments = 0.0_rk
      velocity = 0.0_rk
      eqVel = 0.0_rk
      pdfTmp = 0._rk

      ! element offset for auxField
      elemOff = (iElem-1)*varSys%nAuxScalars

      do iField = 1, nFields
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos
        ! compute field density and first moments
        do iDir = 1, layout%fStencil%QQ
          !store all field pdf in single element to local array pdfTmp
          pdfTmp( vPos(iDir) ) = instate(                               &
            &  neigh((idir-1)* nelems+ ielem)+( ifield-1)* qq+ nscalars*0 )
        enddo
        !field density
        mass_dens( iField ) = auxField(elemOff + dens_pos(iField))

        !field momentum (rho*u)
        first_moments( 1, iField ) = auxField(elemOff + mom_pos(1, iField))
        first_moments( 2, iField ) = auxField(elemOff + mom_pos(2, iField))
        first_moments( 3, iField ) = auxField(elemOff + mom_pos(3, iField))
      enddo

      !total mass density
      totmassDens = sum(mass_dens)

      ! solve linear system of equation for actual moments
      ! number density of all species
      num_dens(:) = mass_dens(:)/molWeight(:)

      !total number density
      totNum_dens = sum(num_dens(:))

      !mole fraction
      moleFrac(:) =  num_dens(:)/totNum_dens

      matA = 0.0_rk
      !build up matrix to solver LSE for actual velocity
      do iField = 1, nFields
        !set diagonal part
        matA(iField, iField) = 1.0_rk
        do iField_2 = 1, nFields
          matA(iField, iField) = matA(iField, iField) + omega * 0.5_rk         &
            &                  * resi_coeff(iField, iField_2) * phi(iField)    &
            &                  * moleFrac(iField_2) / paramB
        end do
        !set non-diagonal part
        do iField_2 = 1, nFields
          matA(iField, iField_2) = matA(iField, iField_2) - omega * 0.5_rk     &
            &                    * resi_coeff(iField, iField_2) * phi(iField_2)&
            &                    * moleFrac(iField) / paramB
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matA )

      !actual momentum of all species
      velocity(1, :) = matmul( invA, first_moments(1,:) )
      velocity(2, :) = matmul( invA, first_moments(2,:) )
      velocity(3, :) = matmul( invA, first_moments(3,:) )

      ! store momentum of untransformed PDF in auxField
      do iField = 1, nFields
        auxField(elemOff + mom_pos(1, iField)) = velocity(1, iField)
        auxField(elemOff + mom_pos(2, iField)) = velocity(2, iField)
        auxField(elemOff + mom_pos(3, iField)) = velocity(3, iField)
      end do

      ! convert momentum to velocity
      velocity(1, :) = velocity(1, :) / mass_dens(:)
      velocity(2, :) = velocity(2, :) / mass_dens(:)
      velocity(3, :) = velocity(3, :) / mass_dens(:)

      ! compute equilibrim velocity
      do iField = 1, nFields
        eqVel( :, iField ) = velocity( :, iField )
        do iField_2 = 1, nFields
          eqVel( :, iField ) = eqVel( :, iField )                              &
            &                + resi_coeff( iField, iField_2 ) * phi(iField)    &
            &                * moleFrac(iField_2)                              &
            &                * (velocity(:, iField_2) - velocity(:,iField))    &
            &                / paramB
        end do
      end do

      !compute mass averaged mixture velocity
      velAvg(1) = dot_product( mass_dens, velocity(1,:) )/totmassDens
      velAvg(2) = dot_product( mass_dens, velocity(2,:) )/totmassDens
      velAvg(3) = dot_product( mass_dens, velocity(3,:) )/totmassDens

      !compute equilibrium and do collision
      do iField = 1, nFields
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos
        !KM: TGV testcase
        velNew = theta_eq*velAvg + (1.0_rk - theta_eq)*eqVel(:, iField)

        !usqr = eqVel(1, iField) * eqVel(1, iField) &
        !  &  + eqVel(2, iField) * eqVel(2, iField) &
        !  &  + eqVel(3, iField) * eqVel(3, iField)
        usqr = dot_product( velNew, velNew )
        usqr = usqr * t2cs2inv

        do iDir = 1, QQ
          ucx = dble( layout%fStencil%cxDir( 1, iDir ))*eqVel(1, iField)     &
            & + dble( layout%fStencil%cxDir( 2, iDir ))*eqVel(2, iField)     &
            & + dble( layout%fStencil%cxDir( 3, iDir ))*eqVel(3, iField)

          ucxQuadTerm = dot_product(                                           &
            &           dble(layout%fStencil%cxDir(:, iDir)), velNew )

          !feq = layout%weight(iDir) * mass_dens(iField) * ( phi(iField)        &
          !  & + ucx * cs2inv + ucx * ucx * t2cs4inv - usqr )
          feq = layout%weight(iDir) * mass_dens(iField) * ( phi(iField)        &
            & + ucx * cs2inv + ucxQuadTerm * ucxQuadTerm * t2cs4inv - usqr )

          if ( iDir == layout%fStencil%restPosition ) then
            ! equilibrium at rest
            feq = layout%weight(iDir) * mass_dens(iField) * ( sRest(iField) &
              &                                              - usqr         )
          end if

          outstate( ( ielem-1)* nscalars+idir+( ifield-1)* qq ) &
            & = pdfTmp(vPos(iDir)) + omega_new * ( feq - pdfTmp(vPos(iDir)) )

        enddo !iDir

      enddo !iField
    enddo nodeloop

  end subroutine bgk_advRel_MSLiquid_generic
! ******************************************************************************

! ******************************************************************************
  !> Unoptimized Advection relaxation routine for the multispecies BGK model
  !! with external forcing term
  !!
  !! This routine contains the implementation based on the paper
  !! "A Lattice Boltzmann Scheme for liquid mixtures - Part II: Discretization
  !! and Numerics, Jens Zudrop, Sabine Roller, Pietro Asinari. "\n
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_forcing_advRel_MSLiquid_generic( fieldProp, inState,    &
    &                                             outState, auxField,    &
    &                                             neigh, nElems, nSolve, &
    &                                             level, layout, params, &
    &                                             varSys, derVarPos      )
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
    !local variables
    real(kind=rk) :: pdfTmp( varSys%nScalars )  !< temporary local pdf values
    integer :: iElem, iField, iField_2, nFields, iDir, QQ, nScalars
    integer :: vPos( layout%fStencil%QQ )
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer, dimension(varSys%nStateVars) :: stateVarMap
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, num_dens
    real(kind=rk), dimension(varSys%nStateVars) :: moleFrac, massFrac
    real(kind=rk), dimension(3, varSys%nStateVars) &
      & :: first_moments, velocity, eqVel, mixForce, diffForce, totForce
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi, &
      & chargeNr
    real(kind=rk) :: totNum_dens, totMass_dens, usqr, ucx, feq
    real(kind=rk) :: omega, omega_new, paramB, minMolWeight
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) &
      & :: matA, invA, resi_coeff
    real(kind=rk) :: forceTerm
    !real(kind=rk), dimension(nSolve-1):: forceState
    real(kind=rk) :: velAvg(3), velNew(3), ucxQuadTerm, theta_eq, diffForce_fac
    real(kind=rk) :: chargeTerm(varSys%nStateVars), charge_dens
    integer :: dens_pos(varSys%nStateVars)
    integer :: mom_pos(3, varSys%nStateVars)
    integer :: elemOff
    real(kind=rk) :: wRestInv, sRest(varSys%nStateVars)
    ! ------------------------------------------------------------------------
    !write(dbgUnit(1),*) 'In BGK forcing ', params%general%simControl%now%iter
    ! initialize forceState
    !forceState = -1.0_rk

    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    QQ     = layout%fStencil%QQ
    nFields = scheme%nFields
    nScalars = varSys%nScalars

    !molecular weights
    molWeight = fieldProp(:)%species%molweight
    minMolWeight = minval(molWeight)
    !write(dbgUnit(1),*) 'molWeight ', molWeight
    !write(dbgUnit(1),*) 'minMolWeight ', minMolWeight
    !specific charge of ionic species
    chargeNr = fieldProp(:)%species%chargeNr
    !molecular weight ratios
    phi = fieldProp(:)%species%molWeigRatio

    !s_0^k - species weight for rest position
    wRestInv = 1.0_rk / layout%weight(layout%fStencil%restPosition)
    do iField = 1, nFields
      sRest(iField) = wRestInv + (1.0_rk - wRestInv) * phi(iField)
    end do

    !resistivities
    do iField = 1, nFields
      resi_coeff(iField,:) = fieldProp( iField )%species%resi_coeff
    enddo

    !write(dbgUnit(1),*) 'gasCosnt ', mixture%gasConst_R_LB
    !write(dbgUnit(1),*) 'temp ', mixture%temp0LB
    !write(dbgUnit(1),*) 'molDens0 ', mixture%moleDens0LB
    omega = scheme%mixture%relaxLvl(Level)%omega_diff
    omega_new = (1.0_rk/(1.0_rk/omega + 0.5_rk))

    paramB = scheme%mixture%paramB
    !equilibrium theta
    theta_eq = scheme%mixture%theta_eq

    stateVarMap = scheme%stateVarMap%varPos%val(:)

    diffForce_fac = minMolWeight / (scheme%mixture%gasConst_R_LB &
      &           * scheme%mixture%temp0LB)

    ! position of density and velocity in auxField array
    do iField = 1, nFields
      dens_pos(iField) = varSys%method%val(derVarPos(iField)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:, iField) = varSys%method%val(derVarPos(iField)%momentum) &
        &                                  %auxField_varPos(1:3)
    end do

    ! loop over elements
    nodeloop: do iElem = 1, nSolve
      !write(dbgUnit(1),*) 'iElem ', iElem

      mass_dens = 0._rk
      first_moments = 0.0_rk
      velocity = 0.0_rk
      eqVel = 0.0_rk
      pdfTmp = 0._rk

      ! element offset to access auxField
      elemOff = (iElem-1)*varSys%nAuxScalars

      do iField = 1, nFields
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos
        ! compute field density and first moments
        do iDir = 1, QQ
          !store all field pdf in single element to local array pdfTmp
          pdfTmp( vPos(iDir) ) = instate(                               &
            &  neigh((idir-1)* nelems+ ielem)+( ifield-1)* qq+ nscalars*0 )
        enddo

        !field density
        mass_dens( iField ) = auxField(elemOff + dens_pos(iField))

        !field momentum (rho*u)
        first_moments( 1, iField ) = auxField(elemOff + mom_pos(1, iField))
        first_moments( 2, iField ) = auxField(elemOff + mom_pos(2, iField))
        first_moments( 3, iField ) = auxField(elemOff + mom_pos(3, iField))
      enddo

      !total mass density
      totmass_Dens = sum(mass_dens)

      !mass fraction
      massFrac(:) = mass_dens(:)/totMass_dens

      ! solve linear system of equation for actual moments
      ! number density of all species
      num_dens(:) = mass_dens(:)/molWeight(:)
      !write(dbgUnit(1),*) 'mass_dens ', mass_dens
      !write(dbgUnit(1),*) 'num_dens ', num_dens

      !total number density
      totNum_dens = sum(num_dens(:))

      !mole fraction
      moleFrac(:) =  num_dens(:)/totNum_dens

      ! compute external forcing term
      ! d^m_k = w_m*c_m*( n0*min_a(m_a)*F_k + (1/cs2)*(rho_k*g) )

      ! forcing term in the mixture for momentum equation
      ! mixForce = 1/cs2 * rho_k * g
      !write(dbgUnit(1),*) ' mixForce '
      do iField = 1, nFields
        mixForce(:,iField) = cs2inv * mass_dens(iField) &
          &                * scheme%mixture%gravityField(:)
        !write(dbgUnit(1),*) mixForce(:, iField)
      end do

      ! diffusive forcing term in the Maxwell-Stefan equations
      ! F_k = 1/(nRT) rho_k (z_k Faraday/m_k E)
      !     - y_k/(nRT) sum_l(rho_l  z_l Faraday/m_l E)

      ! charge density
      do iField = 1, nFields
        chargeTerm(iField) = num_dens(iField) * chargeNr(iField) &
          &                * scheme%mixture%faradayLB
      end do
      charge_dens = sum(chargeTerm)

      ! diffusive force on each species
      ! diffForce = (rho_k z_k /m_k - y_k sum_l rho_l z_l /m_l ) Faraday * E
      diffForce = 0.0_rk
      !write(dbgUnit(1),*) 'Diffusive force on species'
      do iField = 1, nFields
        diffForce(:, iField) = ( chargeTerm(iField) - massFrac(iField) &
          &                  * charge_dens ) * diffForce_fac           &
          &                  * scheme%mixture%electricField(:)
        !write(dbgUnit(1),*) diffForce(:,iField)
      end do

      ! total force in forcing term
      do iField = 1, nFields
        totForce(:, iField) = diffForce(:, iField) + mixForce(:, iField)
      end do
      !write(dbgUnit(1),*) 'totForce ', totForce

      ! Add F/2 to first_moment to recover second-order force
      ! Multiply totForce with cs2 since mixForce is divided by cs2
      ! and diffForce is divided by RT so totForce must be multiplied by
      ! cs2 to get unit for body_force N/m^3
      do iField = 1, nFields
        first_moments(:, iField) = first_moments(:, iField)         &
          &                      + cs2 * totForce(:, iField) * 0.5_rk
      end do

      matA = 0.0_rk
      !build up matrix to solver LSE for actual velocity
      do iField = 1, nFields
        !set diagonal part
        matA(iField, iField) = 1.0_rk
        do iField_2 = 1, nFields
          matA(iField, iField) = matA(iField, iField) + omega * 0.5_rk         &
            &                  * resi_coeff(iField, iField_2) * phi(iField)    &
            &                  * moleFrac(iField_2) / paramB
        end do
        !set non-diagonal part
        do iField_2 = 1, nFields
          matA(iField, iField_2) = matA(iField, iField_2) - omega * 0.5_rk     &
            &                    * resi_coeff(iField, iField_2) * phi(iField_2)&
            &                    * moleFrac(iField) / paramB
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matA )

      ! actual velocity of all species
      velocity(1, :) = matmul( invA, first_moments(1,:) )
      velocity(2, :) = matmul( invA, first_moments(2,:) )
      velocity(3, :) = matmul( invA, first_moments(3,:) )

      ! store momentum of untransformed PDF in auxField
      do iField = 1, nFields
        auxField(elemOff + mom_pos(1, iField)) = velocity(1, iField)
        auxField(elemOff + mom_pos(2, iField)) = velocity(2, iField)
        auxField(elemOff + mom_pos(3, iField)) = velocity(3, iField)
      end do

      ! convert momentum to velocity
      velocity(1, :) = velocity(1, :) / mass_dens(:)
      velocity(2, :) = velocity(2, :) / mass_dens(:)
      velocity(3, :) = velocity(3, :) / mass_dens(:)

      ! compute equilibrim velocity
      do iField = 1, nFields
        eqVel( :, iField ) = velocity( :, iField )
        do iField_2 = 1, nFields
          eqVel( :, iField ) = eqVel( :, iField )                              &
            &                + resi_coeff( iField, iField_2 ) * phi(iField)    &
            &                * moleFrac(iField_2)                              &
            &                * (velocity(:, iField_2) - velocity(:,iField))    &
            &                / paramB
        end do
      end do

      !compute mass averaged mixture velocity
      velAvg(1) = dot_product( mass_dens, velocity(1,:) ) / totmass_Dens
      velAvg(2) = dot_product( mass_dens, velocity(2,:) ) / totmass_Dens
      velAvg(3) = dot_product( mass_dens, velocity(3,:) ) / totmass_Dens

      !compute equilibrium and do collision
      do iField = 1, nFields
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos
        !KM: TGV testcase
        velNew = theta_eq*velAvg + (1.0_rk - theta_eq)*eqVel(:, iField)
        !usqr = eqVel(1, iField) * eqVel(1, iField) &
        !  &  + eqVel(2, iField) * eqVel(2, iField) &
        !  &  + eqVel(3, iField) * eqVel(3, iField)
        usqr = dot_product( velNew, velNew )
        usqr = usqr * t2cs2inv
        do iDir = 1, layout%fStencil%QQ
          !forcing terms
          forceTerm =  &
            &   dble( layout%fStencil%cxDir( 1, iDir ))*totForce(1, iField) &
            & + dble( layout%fStencil%cxDir( 2, iDir ))*totForce(2, iField) &
            & + dble( layout%fStencil%cxDir( 3, iDir ))*totForce(3, iField)

          ucx = dble( layout%fStencil%cxDir( 1, iDir ))*eqVel(1, iField)     &
            & + dble( layout%fStencil%cxDir( 2, iDir ))*eqVel(2, iField)     &
            & + dble( layout%fStencil%cxDir( 3, iDir ))*eqVel(3, iField)

          ucxQuadTerm = dble( layout%fStencil%cxDir( 1, iDir ))*velNew(1)    &
            &         + dble( layout%fStencil%cxDir( 2, iDir ))*velNew(2)    &
            &         + dble( layout%fStencil%cxDir( 3, iDir ))*velNew(3)

          feq = layout%weight(iDir) * mass_dens(iField) * ( phi(iField)        &
            & + ucx * cs2inv + ucxQuadTerm * ucxQuadTerm * t2cs4inv - usqr )

          if ( iDir == layout%fStencil%restPosition ) then
            ! equilibrium at rest
            feq = layout%weight(iDir) * mass_dens(iField) * ( sRest(iField) &
              &                                              - usqr         )
          end if

          outstate( ( ielem-1)* nscalars+idir+( ifield-1)* qq ) &
            & = pdfTmp(vPos(iDir)) + omega_new * ( feq - pdfTmp(vPos(iDir)) )  &
            ! forcing term
            & + layout%weight(iDir) * forceTerm
        enddo !iDir

      enddo !iField

    enddo nodeloop
    !write(dbgUnit(1),*)

  end subroutine bgk_forcing_advRel_MSLiquid_generic
! ******************************************************************************


! ******************************************************************************
  !> Unoptimized Advection relaxation routine for the multispecies BGK model
  !!
  !! This routine contains the implementation based on the paper
  !! "A Lattice Boltzmann Scheme for liquid mixtures - Part II: Discretization
  !! and Numerics, Jens Zudrop, Sabine Roller, Pietro Asinari. "\n
  !! MRT paper
  !! "Lattice Boltzmann liquid mixture modeling for electrodialytic engineering
  !!  applications - J. Zudrop, K. Masilamani, S. Roller and P. Asinari"\n
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mrt_advRel_MSLiquid_generic( fieldProp, inState, outState,    &
    &                                     auxField, neigh, nElems, nSolve, &
    &                                     level, layout, params, varSys,   &
    &                                     derVarPos                        )
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
    !local variables
    real(kind=rk) :: pdfTmp( varSys%nScalars )  !< temporary local pdf values
    integer :: iElem, iField, iField_2, nFields, iDir, QQ, nScalars
    integer :: vPos( layout%fStencil%QQ )
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer, dimension(varSys%nStateVars) :: stateVarMap, dens_pos
    integer :: mom_pos(3, varSys%nStateVars)
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, num_dens, &
      & moleFrac, molWeight, phi
    real(kind=rk), dimension(3, varSys%nStateVars) :: first_moments, &
      & velocity, eqVel
    real(kind=rk) :: totNum_dens, usqr, ucx, feq
    real(kind=rk) :: omega, omega_new, paramB, theta_eq
    real(kind=rk) :: velAvg(3), velNew(3), ucxQuadTerm, totMassDens
    real(kind=rk) :: omegaMoments(layout%fStencil%QQ,layout%fStencil%QQ)
    real(kind=rk) :: fneq(layout%fStencil%QQ), fneq_om(layout%fStencil%QQ)
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: matA, &
      & invA, resi_coeff
    integer :: elemOff
    real(kind=rk) :: wRestInv, sRest(varSys%nStateVars)
    ! ------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    QQ     = layout%fStencil%QQ
    nFields = scheme%nFields
    nScalars = varSys%nScalars

    !molecular weights
    molWeight = fieldProp(:)%species%molweight
    !molecular weight ratios
    phi = fieldProp(:)%species%molWeigRatio

    !s_0^k - species weight for rest position
    wRestInv = 1.0_rk / layout%weight(layout%fStencil%restPosition)
    do iField = 1, nFields
      sRest(iField) = wRestInv + (1.0_rk - wRestInv) * phi(iField)
    end do

    !resistivities
    do iField = 1, nFields
      resi_coeff(iField,:) = fieldProp( iField )%species%resi_coeff
    enddo

    omega = scheme%mixture%relaxLvl(level)%omega_diff
!    write(*,*) 'omega ', omega
    omega_new = (1.0_rk/(1.0_rk/omega + 0.5_rk))
!    write(*,*) 'omega_new ', omega_new

    paramB = scheme%mixture%paramB
    !equilibrium theta
    theta_eq = scheme%mixture%theta_eq
!    write(*,*) 'paramB ', paramB
    stateVarMap = scheme%stateVarMap%varPos%val(:)

!    write(dbgUnit,*) 'B 1', resi_coeff(1,:)
!    write(dbgUnit,*) 'B 2', resi_coeff(2,:)
!    write(dbgUnit,*) 'B 3', resi_coeff(3,:)

    ! position of density and momentum in auxField array
    do iField = 1, nFields
      dens_pos(iField) = varSys%method%val(derVarPos(iField)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:, iField) = varSys%method%val(derVarPos(iField)%momentum) &
        &                                  %auxField_varPos(1:3)
    end do

!write(dbgUnit(1),*)'time ', params%general%simControl%now%iter
    nodeloop: do iElem = 1, nSolve
!write(dbgUnit(1),*) 'iElem', iElem

      mass_dens = 0._rk
      first_moments = 0.0_rk
      velocity = 0.0_rk
      eqVel = 0.0_rk
      pdfTmp = 0._rk

      ! element offset for auxField
      elemOff = (iElem-1)*varSys%nAuxScalars

      do iField = 1, nFields
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos
        ! compute field density and first moments
        do iDir = 1, layout%fStencil%QQ
          !store all field pdf in single element to local array pdfTmp
          pdfTmp( vPos(iDir) ) = instate(                               &
            &  neigh((idir-1)* nelems+ ielem)+( ifield-1)* qq+ nscalars*0 )
        enddo
        !field density
        mass_dens( iField ) = auxField(elemOff + dens_pos(iField))

        !field momentum (rho*u)
        first_moments( 1, iField ) = auxField(elemOff + mom_pos(1, iField))
        first_moments( 2, iField ) = auxField(elemOff + mom_pos(2, iField))
        first_moments( 3, iField ) = auxField(elemOff + mom_pos(3, iField))
      enddo

      !total mass density
      totmassDens = sum(mass_dens)

      ! solve linear system of equation for actual moments
      ! number density of all species
      num_dens(:) = mass_dens(:)/molWeight(:)

      !total number density
      totNum_dens = sum(num_dens(:))

      !mole fraction
      moleFrac(:) =  num_dens(:)/totNum_dens

      matA = 0.0_rk
      !build up matrix to solver LSE for actual velocity
      do iField = 1, nFields
        !set diagonal part
        matA(iField, iField) = 1.0_rk
        do iField_2 = 1, nFields
          matA(iField, iField) = matA(iField, iField) + omega * 0.5_rk         &
            &                  * resi_coeff(iField, iField_2) * phi(iField)    &
            &                  * moleFrac(iField_2) / paramB
        end do
        !set non-diagonal part
        do iField_2 = 1, nFields
          matA(iField, iField_2) = matA(iField, iField_2) - omega * 0.5_rk     &
            &                    * resi_coeff(iField, iField_2) * phi(iField_2)&
            &                    * moleFrac(iField) / paramB
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matA )

      ! actual velocity of all species
      velocity(1, :) = matmul( invA, first_moments(1,:) )
      velocity(2, :) = matmul( invA, first_moments(2,:) )
      velocity(3, :) = matmul( invA, first_moments(3,:) )

      ! store momentum of untransformed PDF in auxField
      do iField = 1, nFields
        auxField(elemOff + mom_pos(1, iField)) = velocity(1, iField)
        auxField(elemOff + mom_pos(2, iField)) = velocity(2, iField)
        auxField(elemOff + mom_pos(3, iField)) = velocity(3, iField)
      end do

      ! convert momentum to velocity
      velocity(1, :) = velocity(1, :) / mass_dens(:)
      velocity(2, :) = velocity(2, :) / mass_dens(:)
      velocity(3, :) = velocity(3, :) / mass_dens(:)

      ! compute equilibrim velocity
      do iField = 1, nFields
        eqVel( :, iField ) = velocity( :, iField )
        do iField_2 = 1, nFields
          eqVel( :, iField ) = eqVel( :, iField )                              &
            &                + resi_coeff( iField, iField_2 ) * phi(iField)    &
            &                * moleFrac(iField_2)                              &
            &                * (velocity(:, iField_2) - velocity(:,iField))    &
            &                / paramB
        end do
      end do

      !compute mass averaged mixture velocity
      velAvg(1) = dot_product( mass_dens, velocity(1,:) )/totmassDens
      velAvg(2) = dot_product( mass_dens, velocity(2,:) )/totmassDens
      velAvg(3) = dot_product( mass_dens, velocity(3,:) )/totmassDens

      !compute equilibrium and do collision
      do iField = 1, nFields
        !omega moments = (M^-1 RelaxMat M)*(I+(M^-1 RelMat M)/2)^-1
        omegaMoments = fieldProp( iField )%species%mrt( level )%omegaMoments
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos

        !KM: TGV testcase
        velNew = theta_eq*velAvg + (1.0_rk - theta_eq)*eqVel(:, iField)

        !usqr = eqVel(1, iField) * eqVel(1, iField) &
        !  &  + eqVel(2, iField) * eqVel(2, iField) &
        !  &  + eqVel(3, iField) * eqVel(3, iField)
        usqr = dot_product( velNew, velNew )
        usqr = usqr * t2cs2inv

        do iDir = 1, layout%fStencil%QQ
          ucx = dble( layout%fStencil%cxDir( 1, iDir ))*eqVel(1, iField)     &
            & + dble( layout%fStencil%cxDir( 2, iDir ))*eqVel(2, iField)     &
            & + dble( layout%fStencil%cxDir( 3, iDir ))*eqVel(3, iField)

          ucxQuadTerm = dot_product(                                           &
            &           dble(layout%fStencil%cxDir(:, iDir)), velNew )

          feq = layout%weight(iDir) * mass_dens(iField) * ( phi(iField)        &
            & + ucx * cs2inv + ucxQuadTerm * ucxQuadTerm * t2cs4inv - usqr )

          if ( iDir == layout%fStencil%restPosition ) then
            ! equilibrium at rest
            feq = layout%weight(iDir) * mass_dens(iField) * ( sRest(iField) &
              &                                              - usqr         )
          end if
          fneq(iDir) = ( feq - pdfTmp( vPos(iDir) ) )
        enddo !iDir
!write(dbgUnit(1),*) 'omegaMoments ', omegaMoments
!write(dbgUnit(1),*) 'before fnEq ', fnEq
        fneq_om = matmul( omegaMoments, fneq )
!do iDir = 1, layout%fStencil%QQ
!  write(dbgUnit(1),*) 'iDir ', iDir, fnEq(iDir), fnEq_om(iDir)
!end do

        do iDir = 1, layout%fStencil%QQ
          outstate( ( ielem-1)* nscalars+idir+( ifield-1)* qq ) &
            & = pdfTmp( vPos(iDir) ) + fneq_om( iDir )
!          if (tem_isnan( outstate( ( varsys%nscalars-1)*neigh + vpos(idir)+( ielem-1)* nelems)) then
!            write(*,*) 'found nan', iElem, iField
        enddo !iDir

!write(dbgUnit(1),*) 'outpdf', outstate(vPos(:))
      enddo !iField
    enddo nodeloop

  end subroutine mrt_advRel_MSLiquid_generic
! ******************************************************************************

! ******************************************************************************
  !> Optimized Advection relaxation routine for the multispecies mrt model
  !! for d3q19 with 3 species
  !!
  !! This routine contains the implementation based on the paper
  !! "A Lattice Boltzmann Scheme for liquid mixtures - Part II: Discretization
  !! and Numerics, Jens Zudrop, Sabine Roller, Pietro Asinari. "\n
  !! MRT paper
  !! "Lattice Boltzmann liquid mixture modeling for electrodialytic engineering
  !!  applications - J. Zudrop, K. Masilamani, S. Roller and P. Asinari"\n
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mrt_advRel_d3q19f3_MSLiquid( fieldProp, inState, outState,    &
    &                                     auxField, neigh, nElems, nSolve, &
    &                                     level, layout, params, varSys,   &
    &                                     derVarPos                        )
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
    !local variables
    ! temporary local pdf values
    real(kind=rk) :: pdfTmp( QQ, 3 )
    integer :: iElem, iField, nScalars
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk), dimension(3) :: rsigma, uxsigma,                            &
      & uysigma, uzsigma, qxsigma, qysigma, qzsigma, gqxsigma, gqysigma,       &
      & gqzsigma
    real(kind=rk), dimension(3) :: uxstar, uystar, uzstar
    real(kind=rk), dimension(3) :: molWeight, molWeight_inv, phi,              &
      &                            num_dens, moleFrac
    real(kind=rk) :: usqr, totNum_dens_inv, totMass_dens
    real(kind=rk) :: velAvg_x, velAvg_y, velAvg_z
    real(kind=rk) :: velQuadTerm_x, velQuadTerm_y, velQuadTerm_z
    real(kind=rk) :: omega, omegadiv2, omega_fac, omega_o
    real(kind=rk), dimension(3,3) :: A, mbb
    real(kind=rk) :: zfac, yfac, a12_11, a13_11, &
      & a23fac, a21_11, a31_11, zcoeff
    real(kind=rk) :: theta_eq, theta_eq_spc, paramB_inv, a11_inv, a22fac_inv
    real(kind=rk) :: B_12, B_13, B_23
    real(kind=rk) :: Bratio_12, Bratio_13, Bratio_23
    real(kind=rk) :: mbbEq_12, mbbEq_13, mbbEq_21, mbbEq_23, mbbEq_31, mbbEq_32
    real(kind=rk) :: chi12, chi13, chi21, chi23, chi31, chi32
    real(kind=rk) :: rho_d, rho_dm, rho_o, rho_om, sum1_1, sum2_1, sum3_1,     &
      & sum4_1, sum5_1, sum6_1, sum7_1, sum8_1, sum9_1, sum1_2, sum2_2, sum3_2,&
      & sum4_2, sum5_2, sum6_2, sum7_2, sum8_2, sum9_2
    real(kind=rk) :: omegaMoments(QQ,QQ,3)
    real(kind=rk) :: fneq(QQ), fneq_om(QQ)
    integer :: dens_pos(3), mom_pos(3,3), elemOff
    ! ---------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nScalars = varSys%nScalars
    !nFields = scheme%nFields
    !molecular weights
    molWeight = fieldProp(:)%species%molweight
    molWeight_inv = 1.0_rk/molWeight
    !molecular weight ratios
    phi = fieldProp(:)%species%molWeigRatio

    !omega moments = (M^-1 RelaxMat M)*(I+(M^-1 RelMat M)/2)^-1
    omegaMoments(:,:,1) = fieldProp(1)%species%mrt( level )%omegaMoments
    omegaMoments(:,:,2) = fieldProp(2)%species%mrt( level )%omegaMoments
    omegaMoments(:,:,3) = fieldProp(3)%species%mrt( level )%omegaMoments

    !omega
    omega = scheme%mixture%relaxLvl(Level)%omega_diff
    omega_fac = (1.0_rk/(1.0_rk/omega + 0.5_rk))
    omegadiv2 = omega*0.5_rk
    omega_o = 1.0_rk - omega_fac

    !free parameter B
    paramB_inv = 1.0_rk/scheme%mixture%paramB

    !equilibrium theta
    theta_eq = scheme%mixture%theta_eq
    theta_eq_spc = 1.0_rk - theta_eq

    !resistivities
    !resistivity coeff are symmetric
    !B_12 = B_21, B_13=B31, B23=B32
    B_12 = fieldProp( 1 )%species%resi_coeff(2)
    B_13 = fieldProp( 1 )%species%resi_coeff(3)
    B_23 = fieldProp( 2 )%species%resi_coeff(3)

    Bratio_12 = B_12*paramB_inv
    Bratio_13 = B_13*paramB_inv
    Bratio_23 = B_23*paramB_inv
    mbbEq_12 = Bratio_12*phi(1)
    mbbEq_13 = Bratio_13*phi(1)

    mbbEq_21 = Bratio_12*phi(2)
    mbbEq_23 = Bratio_23*phi(2)

    mbbEq_31 = Bratio_13*phi(3)
    mbbEq_32 = Bratio_23*phi(3)

    ! factor need to build matrixA to solve LSE
    ! dummy diagonal entries
    mbb(1,1) = -1.0_rk
    mbb(1,2) = mbbEq_21*omegadiv2
    mbb(1,3) = mbbEq_31*omegadiv2

    mbb(2,1) = mbbEq_12*omegadiv2
    mbb(2,2) = -1.0_rk
    mbb(2,3) = mbbEq_32*omegadiv2

    mbb(3,1) = mbbEq_13*omegadiv2
    mbb(3,2) = mbbEq_23*omegadiv2
    mbb(3,3) = -1.0_rk

    ! position of density and momentum in auxField array
    do iField = 1, 3
      dens_pos(iField) = varSys%method%val(derVarPos(iField)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:, iField) = varSys%method%val(derVarPos(iField)%momentum) &
        &                                  %auxField_varPos(1:3)
    end do

    nodeloop: do iElem = 1, nSolve

      !species 1
      pdfTmp( Q__W,1 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__S,1 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__B,1 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__E,1 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__N,1 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__T,1 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,1 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,1 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,1 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,1 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,1 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,1 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,1 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,1 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,1 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,1 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,1 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,1 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__0,1 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      !species 2
      pdfTmp( Q__W,2 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__S,2 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__B,2 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__E,2 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__N,2 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__T,2 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,2 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,2 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,2 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,2 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,2 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,2 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,2 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,2 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,2 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,2 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,2 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,2 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__0,2 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      !species 3
      pdfTmp( Q__W,3 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__S,3 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__B,3 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__E,3 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__N,3 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__T,3 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,3 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,3 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,3 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,3 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,3 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,3 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,3 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,3 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,3 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,3 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,3 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,3 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__0,3 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)

      ! element offset to access auxField
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! get density and momentum from auxField
      do iField = 1, 3
        ! density
        rsigma(iField) = auxField(elemOff + dens_pos(iField))
        ! momentum
        gqxsigma(iField) = auxField(elemOff + mom_pos(1, iField))
        gqysigma(iField) = auxField(elemOff + mom_pos(2, iField))
        gqzsigma(iField) = auxField(elemOff + mom_pos(3, iField))
      end do

      !number density
      num_dens(:) = rsigma(:)*molWeight_inv(:)

      !total number density
      totNum_dens_inv = 1.0_rk/sum(num_dens(:))

      !mole fraction
      moleFrac(:) = num_dens(:)*totNum_dens_inv

      !matrix A of the velocity lse
      A(1,1) = 1.0_rk + moleFrac(2)*mbb(2,1) + moleFrac(3)*mbb(3,1)
      A(1,2) = - moleFrac(1)*mbb(1,2)
      A(1,3) = - moleFrac(1)*mbb(1,3)

      A(2,1) = - moleFrac(2)*mbb(2,1)
      A(2,2) = 1.0_rk + moleFrac(1)*mbb(1,2) + moleFrac(3)*mbb(3,2)
      A(2,3) = - moleFrac(2)*mbb(2,3)

      A(3,1) = - moleFrac(3)*mbb(3,1)
      A(3,2) = - moleFrac(3)*mbb(3,2)
      A(3,3) = 1.0_rk + moleFrac(1)*mbb(1,3) + moleFrac(2)*mbb(2,3)

      a11_inv = 1.0_rk/A(1,1)
      a13_11 = A(1,3)*a11_inv
      a12_11 = A(1,2)*a11_inv
      a31_11 = A(3,1)*a11_inv
      a21_11 = A(2,1)*a11_inv
      a23fac = A(2,3) - A(2,1)*a13_11
      a22fac_inv = 1.0_rk/(A(2,2) - A(2,1)*a12_11)

      zfac = (A(3,1)*a12_11 - A(3,2)) * a22fac_inv
      zcoeff = 1._rk / (a23fac*zfac - A(3,1)*a13_11 + A(3,3))
      !x-momentum
      !species 3
      yfac = gqxsigma(2) - gqxsigma(1)*a21_11

      qxsigma(3) = (gqxsigma(3) - a31_11*gqxsigma(1) + yfac*zfac)*zcoeff

      !species 2
      qxsigma(2) = (yfac - a23fac*qxsigma(3))*a22fac_inv
      !species 1
      qxsigma(1) = (gqxsigma(1) - A(1,2)*qxsigma(2) - A(1,3)*qxsigma(3))*a11_inv
      !y-momentum
      yfac = gqysigma(2) - gqysigma(1)*a21_11
      qysigma(3) = (gqysigma(3) - a31_11*gqysigma(1) + yfac*zfac)*zcoeff

      qysigma(2) = (yfac - a23fac*qysigma(3))*a22fac_inv

      qysigma(1) = (gqysigma(1) - A(1,2)*qysigma(2) - A(1,3)*qysigma(3))*a11_inv
      !z-momentum
      yfac = gqzsigma(2) - gqzsigma(1)*a21_11
      qzsigma(3) = (gqzsigma(3) - a31_11*gqzsigma(1) + yfac*zfac)*zcoeff

      qzsigma(2) = (yfac - a23fac*qzsigma(3))*a22fac_inv

      qzsigma(1) = (gqzsigma(1) - A(1,2)*qzsigma(2) - A(1,3)*qzsigma(3))*a11_inv

      ! store momentum of untransformed PDF in auxField
      do iField = 1, 3
        auxField(elemOff + mom_pos(1, iField)) = qxsigma(iField)
        auxField(elemOff + mom_pos(2, iField)) = qysigma(iField)
        auxField(elemOff + mom_pos(3, iField)) = qzsigma(iField)
      end do

      ! species velocity
      uxsigma = qxsigma/rsigma
      uysigma = qysigma/rsigma
      uzsigma = qzsigma/rsigma

      !compute equilibrium velocity to compute feq
      chi12 = mbbEq_12 * moleFrac(2)
      chi13 = mbbEq_13 * moleFrac(3)
      chi21 = mbbEq_21 * moleFrac(1)
      chi23 = mbbEq_23 * moleFrac(3)
      chi31 = mbbEq_31 * moleFrac(1)
      chi32 = mbbEq_32 * moleFrac(2)
      !ux
      uxstar(1) =  uxsigma(1) + chi12*(uxsigma(2)-uxsigma(1)) &
        & + chi13*(uxsigma(3)-uxsigma(1))

      uxstar(2) =  uxsigma(2) + chi21*(uxsigma(1)-uxsigma(2)) &
        & + chi23*(uxsigma(3)-uxsigma(2))

      uxstar(3) =  uxsigma(3) + chi31*(uxsigma(1)-uxsigma(3)) &
        & + chi32*(uxsigma(2)-uxsigma(3))
      !uy
      uystar(1) =  uysigma(1) + chi12*(uysigma(2)-uysigma(1)) &
        & + chi13*(uysigma(3)-uysigma(1))

      uystar(2) =  uysigma(2) + chi21*(uysigma(1)-uysigma(2)) &
        & + chi23*(uysigma(3)-uysigma(2))

      uystar(3) =  uysigma(3) + chi31*(uysigma(1)-uysigma(3)) &
        & + chi32*(uysigma(2)-uysigma(3))
      !uz
      uzstar(1) =  uzsigma(1) + chi12*(uzsigma(2)-uzsigma(1)) &
        & + chi13*(uzsigma(3)-uzsigma(1))

      uzstar(2) =  uzsigma(2) + chi21*(uzsigma(1)-uzsigma(2)) &
        & + chi23*(uzsigma(3)-uzsigma(2))

      uzstar(3) =  uzsigma(3) + chi31*(uzsigma(1)-uzsigma(3)) &
        & + chi32*(uzsigma(2)-uzsigma(3))

      ! total mass density
      totMass_dens = sum(rsigma)

      ! mass averaged mixture velocity
      velAvg_x = dot_product( rsigma, uxsigma ) / totMass_dens
      velAvg_y = dot_product( rsigma, uysigma ) / totMass_dens
      velAvg_z = dot_product( rsigma, uzsigma ) / totMass_dens

      !compute equilibrium and do collision
      do iField = 1, 3!nFields
        ! \todo KM: Bug remove omega factor from equilibrium
        rho_d = rsigma(iField) * div1_4

        velQuadTerm_x = theta_eq*velAvg_x + theta_eq_spc*uxstar(iField)
        velQuadTerm_y = theta_eq*velAvg_y + theta_eq_spc*uystar(iField)
        velQuadTerm_z = theta_eq*velAvg_z + theta_eq_spc*uzstar(iField)

        usqr = (velQuadTerm_x * velQuadTerm_x                                  &
          &  + velQuadTerm_y * velQuadTerm_y                                   &
          &  + velQuadTerm_z * velQuadTerm_z) * t2cs2inv

        ! at rest
        fneq(Q__0) = div1_3 * rsigma(iField) * (( 3._rk - 2._rk * phi(iField)) &
          &        - usqr ) - pdfTmp(Q__0, iField)

        rho_dm = rsigma(iField) * (phi(iField) - usqr) * div1_18
        !directional velocity factor
        sum1_1 = rho_d * uxstar(iField) * div3_4h
        sum1_2 = rho_d * velQuadTerm_x**2 + rho_dm
        fneq(Q__W) = -sum1_1 + sum1_2 - pdfTmp(Q__W,iField)
        fneq(Q__E) = sum1_1 + sum1_2 - pdfTmp(Q__E,iField)

        sum2_1 = rho_d * uystar(iField) * div3_4h
        sum2_2 = rho_d * velQuadTerm_y**2 + rho_dm
        fneq(Q__S) = -sum2_1 + sum2_2 - pdfTmp(Q__S,iField)
        fneq(Q__N) = sum2_1 + sum2_2 - pdfTmp(Q__N,iField)

        sum3_1 = rho_d * uzstar(iField) * div3_4h
        sum3_2 = rho_d * velQuadTerm_z**2 + rho_dm
        fneq(Q__B) = -sum3_1 + sum3_2 - pdfTmp(Q__B,iField)
        fneq(Q__T) = sum3_1 + sum3_2 - pdfTmp(Q__T,iField)

        rho_o = rho_d * 0.5_rk
        rho_om = rho_dm * 0.5_rk
        !top north / bottom south
        sum4_1 = rho_o * ( uystar(iField) + uzstar(iField) ) *div3_4h
        sum4_2 = rho_o * ( velQuadTerm_y + velQuadTerm_z )**2 &
        &       + rho_om
        fneq(Q_BS) = -sum4_1 + sum4_2 - pdfTmp(Q_BS,iField)
        fneq(Q_TN) = sum4_1 + sum4_2 - pdfTmp(Q_TN,iField)

        !top south / bottom north
        sum5_1 = rho_o * ( -uystar(iField) + uzstar(iField) ) *div3_4h
        sum5_2 = rho_o * ( -velQuadTerm_y + velQuadTerm_z )**2&
          &       + rho_om
        fneq(Q_BN) = -sum5_1 + sum5_2 - pdfTmp(Q_BN,iField)
        fneq(Q_TS) = sum5_1 + sum5_2 - pdfTmp(Q_TS,iField)

        !top east / bottom west
        sum6_1 = rho_o * ( uxstar(iField) + uzstar(iField) ) *div3_4h
        sum6_2 = rho_o * ( velQuadTerm_x + velQuadTerm_z )**2 &
          &       + rho_om
        fneq(Q_BW) = -sum6_1 + sum6_2 - pdfTmp(Q_BW,iField)
        fneq(Q_TE) = sum6_1 + sum6_2 - pdfTmp(Q_TE,iField)

        !top west / bottom east
        sum7_1 = rho_o * ( -uxstar(iField) + uzstar(iField) ) *div3_4h
        sum7_2 = rho_o * ( -velQuadTerm_x + velQuadTerm_z )**2&
          &       + rho_om
        fneq(Q_BE) = -sum7_1 + sum7_2 - pdfTmp(Q_BE,iField)
        fneq(Q_TW) = sum7_1 + sum7_2 - pdfTmp(Q_TW,iField)

        !north east / south west
        sum8_1 = rho_o * ( uxstar(iField) + uystar(iField) ) *div3_4h
        sum8_2 = rho_o * ( velQuadTerm_x + velQuadTerm_y )**2 &
          &       + rho_om
        fneq(Q_SW) = -sum8_1 + sum8_2 - pdfTmp(Q_SW,iField)
        fneq(Q_NE) = sum8_1 + sum8_2 - pdfTmp(Q_NE,iField)

        !north west / south east
        sum9_1 = rho_o * ( -uxstar(iField) + uystar(iField) ) *div3_4h
        sum9_2 = rho_o * ( -velQuadTerm_x + velQuadTerm_y )**2&
          &       + rho_om
        fneq(Q_SE) = -sum9_1 + sum9_2 - pdfTmp(Q_SE,iField)
        fneq(Q_NW) = sum9_1 + sum9_2 - pdfTmp(Q_NW,iField)

        ! omega * fNonEq
        fneq_om = matmul( omegaMoments(:,:, iField), fneq )

        !outstate
        outstate( ( ielem-1)* nscalars+ q__0+( ifield-1)* qq ) =       &
          & pdfTmp( Q__0, iField ) + fneq_om( Q__0 )

        outstate( ( ielem-1)* nscalars+ q__w+( ifield-1)* qq ) =       &
          & pdfTmp( Q__W, iField ) + fneq_om( Q__W )
        outstate( ( ielem-1)* nscalars+ q__e+( ifield-1)* qq ) =       &
          & pdfTmp( Q__E, iField ) + fneq_om( Q__E )
        outstate( ( ielem-1)* nscalars+ q__s+( ifield-1)* qq ) =       &
          & pdfTmp( Q__S, iField ) + fneq_om( Q__S )
        outstate( ( ielem-1)* nscalars+ q__n+( ifield-1)* qq ) =       &
          & pdfTmp( Q__N, iField ) + fneq_om( Q__N )
        outstate( ( ielem-1)* nscalars+ q__b+( ifield-1)* qq ) =       &
          & pdfTmp( Q__B, iField ) + fneq_om( Q__B )
        outstate( ( ielem-1)* nscalars+ q__t+( ifield-1)* qq ) =       &
          & pdfTmp( Q__T, iField ) + fneq_om( Q__T )
        outstate( ( ielem-1)* nscalars+ q_bs+( ifield-1)* qq ) =       &
          & pdfTmp( Q_BS, iField ) + fneq_om( Q_BS )
        outstate( ( ielem-1)* nscalars+ q_tn+( ifield-1)* qq ) =       &
          & pdfTmp( Q_TN, iField ) + fneq_om( Q_TN )
        outstate( ( ielem-1)* nscalars+ q_bn+( ifield-1)* qq ) =       &
          & pdfTmp( Q_BN, iField ) + fneq_om( Q_BN )
        outstate( ( ielem-1)* nscalars+ q_ts+( ifield-1)* qq ) =       &
          & pdfTmp( Q_TS, iField ) + fneq_om( Q_TS )
        outstate( ( ielem-1)* nscalars+ q_bw+( ifield-1)* qq ) =       &
          & pdfTmp( Q_BW, iField ) + fneq_om( Q_BW )
        outstate( ( ielem-1)* nscalars+ q_te+( ifield-1)* qq ) =       &
          & pdfTmp( Q_TE, iField ) + fneq_om( Q_TE )
        outstate( ( ielem-1)* nscalars+ q_be+( ifield-1)* qq ) =       &
          & pdfTmp( Q_BE, iField ) + fneq_om( Q_BE )
        outstate( ( ielem-1)* nscalars+ q_tw+( ifield-1)* qq ) =       &
          & pdfTmp( Q_TW, iField ) + fneq_om( Q_TW )
        outstate( ( ielem-1)* nscalars+ q_sw+( ifield-1)* qq ) =       &
          & pdfTmp( Q_SW, iField ) + fneq_om( Q_SW )
        outstate( ( ielem-1)* nscalars+ q_ne+( ifield-1)* qq ) =       &
          & pdfTmp( Q_NE, iField ) + fneq_om( Q_NE )
        outstate( ( ielem-1)* nscalars+ q_se+( ifield-1)* qq ) =       &
          & pdfTmp( Q_SE, iField ) + fneq_om( Q_SE )
        outstate( ( ielem-1)* nscalars+ q_nw+( ifield-1)* qq ) =       &
          & pdfTmp( Q_NW, iField ) + fneq_om( Q_NW )

      enddo !iField
    enddo nodeloop

  end subroutine mrt_advRel_d3q19f3_MSLiquid
! ****************************************************************************** !

! ******************************************************************************
  !> Semi-optimized Advection relaxation routine for the MSLiquid BGK model
  !! for d3q19 layout with three species with thermodynamic factor.
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
  subroutine bgk_advRel_d3q19f3_MSLiquid_WTDF( fieldProp, inState, outState,  &
    &                                          auxField, neigh, nElems,       &
    &                                          nSolve, level, layout, params, &
    &                                          varSys, derVarPos              )
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
    !local variables
    ! temporary local pdf values
    real(kind=rk) :: pdfTmp( QQ, 3 )
    integer :: iElem, nScalars, iField, iField_2, iField_3
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk), dimension(3) :: rsigma,                      &
      &                            uxsigma, uysigma, uzsigma,   &
      &                            qxsigma, qysigma, qzsigma,   &
      &                            gqxsigma, gqysigma, gqzsigma
    real(kind=rk), dimension(3) :: uxstar, uystar, uzstar
    real(kind=rk), dimension(3) :: molWeight, molWeight_inv, phi, &
      &                            num_dens, moleFrac
    real(kind=rk) :: usqr(3), totNum_dens_inv, totMass_densInv
    real(kind=rk) :: velAvg(3)
    real(kind=rk) :: velQuadTerm_x(3), velQuadTerm_y(3), velQuadTerm_z(3)
    real(kind=rk) :: omega, omegadiv2, omega_fac, omega_o
    real(kind=rk) :: theta_eq, theta_eq_spc, paramB_inv
    real(kind=rk) :: B_12, B_13, B_23
    real(kind=rk) :: Bratio_12, Bratio_13, Bratio_23
    real(kind=rk) :: mbbEq_12, mbbEq_13, mbbEq_21, mbbEq_23, mbbEq_31, mbbEq_32
    real(kind=rk) :: chi12, chi13, chi21, chi23, chi31, chi32
    real(kind=rk), dimension(3) :: omegaRho, omegaRho_d, omegaRho_dm,      &
      &                            omegaRho_o, omegaRho_om,                &
      &                            sum1_1, sum2_1, sum3_1, sum4_1, sum5_1, &
      &                            sum6_1, sum7_1, sum8_1, sum9_1,         &
      &                            sum1_2, sum2_2, sum3_2, sum4_2, sum5_2, &
      &                            sum6_2, sum7_2, sum8_2, sum9_2
    integer :: dens_pos(3), mom_pos(3,3), elemOff
    real(kind=rk), dimension(3, 3) :: matA, invA, resi_coeff, diff_coeff, &
      &                               thermodynamic_fac, inv_thermodyn_fac
    real(kind=rk) :: temp, press, phy_moleDens_fac
    ! ---------------------------------------------------------------------------

    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nScalars = varSys%nScalars
    !KM \todo check moleDens for multilevel
    phy_moleDens_fac = params%physics%moleDens0
    !molecular weights
    molWeight = fieldProp(:)%species%molweight
    molWeight_inv = 1.0_rk/molWeight
    !molecular weight ratios
    phi = fieldProp(:)%species%molWeigRatio

    !omega
    omega = scheme%mixture%relaxLvl(Level)%omega_diff
    omega_fac = (1.0_rk/(1.0_rk/omega + 0.5_rk))
    omegadiv2 = omega*0.5_rk
    omega_o = 1.0_rk - omega_fac

    !free parameter B
    paramB_inv = 1.0_rk/scheme%mixture%paramB

    !equilibrium theta
    theta_eq = scheme%mixture%theta_eq
    theta_eq_spc = 1.0_rk - theta_eq

    ! temperature
    temp = scheme%mixture%temp0
    ! atmospheric pressure
    press = scheme%mixture%atm_press

    !resistivities
    !resistivity coeff are symmetric
    !B_12 = B_21, B_13=B31, B23=B32
    B_12 = fieldProp( 1 )%species%resi_coeff(2)
    B_13 = fieldProp( 1 )%species%resi_coeff(3)
    B_23 = fieldProp( 2 )%species%resi_coeff(3)

    Bratio_12 = B_12*paramB_inv
    Bratio_13 = B_13*paramB_inv
    Bratio_23 = B_23*paramB_inv

    mbbEq_12 = Bratio_12*phi(1)
    mbbEq_13 = Bratio_13*phi(1)

    mbbEq_21 = Bratio_12*phi(2)
    mbbEq_23 = Bratio_23*phi(2)

    mbbEq_31 = Bratio_13*phi(3)
    mbbEq_32 = Bratio_23*phi(3)

    ! position of density and momentum in auxField array
    do iField = 1, scheme%nFields
      dens_pos(iField) = varSys%method%val(derVarPos(iField)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:, iField) = varSys%method%val(derVarPos(iField)%momentum) &
        &                                  %auxField_varPos(1:3)
    end do

    nodeloop: do iElem = 1, nSolve

      !species 1
      pdfTmp( Q__W,1 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__S,1 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__B,1 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__E,1 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__N,1 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__T,1 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,1 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,1 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,1 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,1 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,1 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,1 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,1 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,1 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,1 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,1 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,1 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,1 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__0,1 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

      !species 2
      pdfTmp( Q__W,2 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__S,2 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__B,2 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__E,2 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__N,2 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__T,2 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,2 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,2 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,2 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,2 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,2 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,2 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,2 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,2 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,2 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,2 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,2 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,2 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__0,2 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)

      !species 3
      pdfTmp( Q__W,3 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__S,3 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__B,3 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__E,3 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__N,3 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__T,3 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,3 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,3 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,3 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,3 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,3 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,3 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,3 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,3 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,3 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,3 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,3 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,3 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__0,3 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)

      ! element offset to access auxField
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! get density and momentum from auxField
      do iField = 1, 3
        ! density
        rsigma(iField) = auxField(elemOff + dens_pos(iField))
        ! momentum
        gqxsigma(iField) = auxField(elemOff + mom_pos(1, iField))
        gqysigma(iField) = auxField(elemOff + mom_pos(2, iField))
        gqzsigma(iField) = auxField(elemOff + mom_pos(3, iField))
      end do

      !number density
      num_dens(:) = rsigma(:)*molWeight_inv(:)

      !total number density
      totNum_dens_inv = 1.0_rk/(num_dens(1)+num_dens(2)+num_dens(3))

      !mole fraction
      moleFrac(:) = num_dens(:)*totNum_dens_inv

      ! MS-Diff coeff matrix from C++ code
      call mus_calc_MS_DiffMatrix( 3, temp, press, num_dens*phy_moleDens_fac, &
        &                          diff_coeff )

      ! Convert to lattice unit
      resi_coeff = params%physics%fac(level)%diffusivity/diff_coeff

      ! Thermodynamic factor from C++ code
      call mus_calc_thermFactor( 3, temp, press, moleFrac, thermodynamic_fac )

      inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

      matA = 0.0_rk
      !build up matrix to solver LSE for actual velocity
      do iField = 1, 3
        !set diagonal part
        matA(iField, iField) = 1.0_rk
        do iField_2 = 1, 3
          do iField_3 = 1, 3
            matA(iField, iField_2) = matA(iField, iField_2) + omega * 0.5_rk   &
              &                    * inv_thermodyn_fac(iField, iField_2)       &
              &                    * resi_coeff(iField_2, iField_3)            &
              &                    * phi(iField_2) * moleFrac(iField_3)        &
              &                    * paramB_inv
          end do
        end do
        !set non-diagonal part
        do iField_2 = 1, 3
          do iField_3 = 1, 3
            matA(iField, iField_3) = matA(iField, iField_3) - omega * 0.5_rk   &
              &                    * inv_thermodyn_fac(iField, iField_2)       &
              &                    * resi_coeff(iField_2, iField_3)            &
              &                    * phi(iField_3) * moleFrac(iField_2)        &
              &                    * paramB_inv
          end do
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matA )

      !actual momentum of all species
      qxsigma(1:3) = matmul( invA, gqxsigma(1:3) )
      qysigma(1:3) = matmul( invA, gqysigma(1:3) )
      qzsigma(1:3) = matmul( invA, gqzsigma(1:3) )

      ! store momentum of untransformed PDF in auxField
      do iField = 1, 3
        auxField(elemOff + mom_pos(1, iField)) = qxsigma(iField)
        auxField(elemOff + mom_pos(2, iField)) = qysigma(iField)
        auxField(elemOff + mom_pos(3, iField)) = qzsigma(iField)
      end do

      ! species velocity
      uxsigma = qxsigma/rsigma
      uysigma = qysigma/rsigma
      uzsigma = qzsigma/rsigma

!write(dbgUnit,*) 'velocity 1', uxsigma(1), uysigma(1), uzsigma(1)
!write(dbgUnit,*) 'velocity 2', uxsigma(2), uysigma(2), uzsigma(2)
!write(dbgUnit,*) 'velocity 3', uxsigma(3), uysigma(3), uzsigma(3)

      !compute equilibrium velocity to compute feq
      chi12 = mbbEq_12 * moleFrac(2)
      chi13 = mbbEq_13 * moleFrac(3)
      chi21 = mbbEq_21 * moleFrac(1)
      chi23 = mbbEq_23 * moleFrac(3)
      chi31 = mbbEq_31 * moleFrac(1)
      chi32 = mbbEq_32 * moleFrac(2)
      !ux
      uxstar(1) =  uxsigma(1) + chi12*(uxsigma(2)-uxsigma(1)) &
        &       + chi13*(uxsigma(3)-uxsigma(1))

      uxstar(2) =  uxsigma(2) + chi21*(uxsigma(1)-uxsigma(2)) &
        &       + chi23*(uxsigma(3)-uxsigma(2))

      uxstar(3) =  uxsigma(3) + chi31*(uxsigma(1)-uxsigma(3)) &
        &       + chi32*(uxsigma(2)-uxsigma(3))
      !uy
      uystar(1) =  uysigma(1) + chi12*(uysigma(2)-uysigma(1)) &
        &       + chi13*(uysigma(3)-uysigma(1))

      uystar(2) =  uysigma(2) + chi21*(uysigma(1)-uysigma(2)) &
        &       + chi23*(uysigma(3)-uysigma(2))

      uystar(3) =  uysigma(3) + chi31*(uysigma(1)-uysigma(3)) &
        &       + chi32*(uysigma(2)-uysigma(3))
      !uz
      uzstar(1) =  uzsigma(1) + chi12*(uzsigma(2)-uzsigma(1)) &
        &       + chi13*(uzsigma(3)-uzsigma(1))

      uzstar(2) =  uzsigma(2) + chi21*(uzsigma(1)-uzsigma(2)) &
        &       + chi23*(uzsigma(3)-uzsigma(2))

      uzstar(3) =  uzsigma(3) + chi31*(uzsigma(1)-uzsigma(3)) &
        &       + chi32*(uzsigma(2)-uzsigma(3))

!write(dbgUnit,*) 'eqVel 1', uxstar(1), uystar(1), uzstar(1)
!write(dbgUnit,*) 'eqVel 2', uxstar(2), uystar(2), uzstar(2)
!write(dbgUnit,*) 'eqVel 3', uxstar(3), uystar(3), uzstar(3)

      ! total mass density
      totMass_densInv = 1.0_rk/(rsigma(1)+rsigma(2)+rsigma(3))

      ! mass averaged mixture velocity three components
      velAvg(1) = ( rsigma(1) * uxsigma(1) &
        &         + rsigma(2) * uxsigma(2) &
        &         + rsigma(3) * uxsigma(3) ) * totMass_densInv * theta_eq

      velAvg(2) = ( rsigma(1) * uysigma(1) &
        &         + rsigma(2) * uysigma(2) &
        &         + rsigma(3) * uysigma(3) ) * totMass_densInv * theta_eq

      velAvg(3) = ( rsigma(1) * uzsigma(1) &
        &         + rsigma(2) * uzsigma(2) &
        &         + rsigma(3) * uzsigma(3) ) * totMass_densInv * theta_eq

      ! velocity in quadratic term of equilibrium
      velQuadTerm_x(1) = velAvg(1) + theta_eq_spc*uxstar(1)
      velQuadTerm_y(1) = velAvg(2) + theta_eq_spc*uystar(1)
      velQuadTerm_z(1) = velAvg(3) + theta_eq_spc*uzstar(1)

      !compute equilibrium and do collision
      usqr(1) = (velQuadTerm_x(1) * velQuadTerm_x(1)                           &
        &     + velQuadTerm_y(1) * velQuadTerm_y(1)                            &
        &     + velQuadTerm_z(1) * velQuadTerm_z(1)) * t2cs2inv

      omegaRho(:) = rsigma(:)*omega_fac
      omegaRho_d(:) = omegaRho(:) * div1_4

      ! species 1
      ! equilibrium at rest
      outstate( ( ielem-1)* nscalars+ q__0+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q__0,1) +                                       &
              & div1_3 * omegaRho(1) * (( 3._rk - 2._rk * phi(1)) - usqr(1) )

      omegaRho_dm(1) = omegaRho(1) * (phi(1) - usqr(1)) * div1_18

      !directional velocity factor
      sum1_1(1) = omegaRho_d(1) * uxstar(1) * div3_4h
      sum1_2(1) = omegaRho_d(1) * velQuadTerm_x(1)**2 + omegaRho_dm(1)

      outstate( ( ielem-1)* nscalars+ q__w+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__W,1) - sum1_1(1) + sum1_2(1)
      outstate( ( ielem-1)* nscalars+ q__e+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__E,1) + sum1_1(1) + sum1_2(1)

      sum2_1(1) = omegaRho_d(1) * uystar(1) * div3_4h
      sum2_2(1) = omegaRho_d(1) * velQuadTerm_y(1)**2 + omegaRho_dm(1)

      outstate( ( ielem-1)* nscalars+ q__s+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__S,1) - sum2_1(1) + sum2_2(1)
      outstate( ( ielem-1)* nscalars+ q__n+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__N,1) + sum2_1(1) + sum2_2(1)

      sum3_1(1) = omegaRho_d(1) * uzstar(1) * div3_4h
      sum3_2(1) = omegaRho_d(1) * velQuadTerm_z(1)**2 + omegaRho_dm(1)

      outstate( ( ielem-1)* nscalars+ q__b+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__B,1) - sum3_1(1) + sum3_2(1)
      outstate( ( ielem-1)* nscalars+ q__t+( 1-1)* qq )         &
        & = omega_o * pdfTmp(Q__T,1) + sum3_1(1) + sum3_2(1)

      !top north / bottom south
      omegaRho_o(1) = omegaRho_d(1) * 0.5_rk
      omegaRho_om(1) = omegaRho_dm(1) * 0.5_rk
      sum4_1(1) = omegaRho_o(1) * ( uystar(1) + uzstar(1) ) *div3_4h
      sum4_2(1) = omegaRho_o(1) * ( velQuadTerm_y(1) + velQuadTerm_z(1) )**2 &
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_bs+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BS,1) - sum4_1(1) + sum4_2(1)
      outstate( ( ielem-1)* nscalars+ q_tn+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TN,1) + sum4_1(1) + sum4_2(1)

      !top south / bottom north
      sum5_1(1) = omegaRho_o(1) * ( -uystar(1) + uzstar(1) ) *div3_4h
      sum5_2(1) = omegaRho_o(1) * ( -velQuadTerm_y(1) + velQuadTerm_z(1) )**2&
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_bn+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BN,1) - sum5_1(1) + sum5_2(1)
      outstate( ( ielem-1)* nscalars+ q_ts+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TS,1) + sum5_1(1) + sum5_2(1)

      !top east / bottom west
      sum6_1(1) = omegaRho_o(1) * ( uxstar(1) + uzstar(1) ) *div3_4h
      sum6_2(1) = omegaRho_o(1) * ( velQuadTerm_x(1) + velQuadTerm_z(1) )**2 &
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_bw+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BW,1) - sum6_1(1) + sum6_2(1)
      outstate( ( ielem-1)* nscalars+ q_te+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TE,1) + sum6_1(1) + sum6_2(1)

      !top west / bottom east
      sum7_1(1) = omegaRho_o(1) * ( -uxstar(1) + uzstar(1) ) *div3_4h
      sum7_2(1) = omegaRho_o(1) * ( -velQuadTerm_x(1) + velQuadTerm_z(1) )**2&
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_be+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BE,1) - sum7_1(1) + sum7_2(1)
      outstate( ( ielem-1)* nscalars+ q_tw+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TW,1) + sum7_1(1) + sum7_2(1)

      !north east / south west
      sum8_1(1) = omegaRho_o(1) * ( uxstar(1) + uystar(1) ) *div3_4h
      sum8_2(1) = omegaRho_o(1) * ( velQuadTerm_x(1) + velQuadTerm_y(1) )**2 &
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_sw+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SW,1) - sum8_1(1) + sum8_2(1)
      outstate( ( ielem-1)* nscalars+ q_ne+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NE,1) + sum8_1(1) + sum8_2(1)

      !north west / south east
      sum9_1(1) = omegaRho_o(1) * ( -uxstar(1) + uystar(1) ) *div3_4h
      sum9_2(1) = omegaRho_o(1) * ( -velQuadTerm_x(1) + velQuadTerm_y(1) )**2&
        &       + omegaRho_om(1)

      outstate( ( ielem-1)* nscalars+ q_se+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SE,1) - sum9_1(1) + sum9_2(1)
      outstate( ( ielem-1)* nscalars+ q_nw+( 1-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NW,1) + sum9_1(1) + sum9_2(1)

      ! species 2
      ! velocity in quadratic term of equilibrium
      velQuadTerm_x(2) = velAvg(1) + theta_eq_spc*uxstar(2)
      velQuadTerm_y(2) = velAvg(2) + theta_eq_spc*uystar(2)
      velQuadTerm_z(2) = velAvg(3) + theta_eq_spc*uzstar(2)

      usqr(2) = (velQuadTerm_x(2) * velQuadTerm_x(2)                           &
        &     + velQuadTerm_y(2) * velQuadTerm_y(2)                            &
        &     + velQuadTerm_z(2) * velQuadTerm_z(2)) * t2cs2inv

      outstate( ( ielem-1)* nscalars+ q__0+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q__0,2) + &
              & div1_3 * omegaRho(2) * (( 3._rk - 2._rk * phi(2)) - usqr(2) )

      omegaRho_dm(2) = omegaRho(2) * (phi(2) - usqr(2)) * div1_18
      !directional velocity factor
      sum1_1(2) = omegaRho_d(2) * uxstar(2) * div3_4h
      sum1_2(2) = omegaRho_d(2) * velQuadTerm_x(2)**2 + omegaRho_dm(2)

      outstate( ( ielem-1)* nscalars+ q__w+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__W,2) - sum1_1(2) + sum1_2(2)
      outstate( ( ielem-1)* nscalars+ q__e+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__E,2) + sum1_1(2) + sum1_2(2)

      sum2_1(2) = omegaRho_d(2) * uystar(2) * div3_4h
      sum2_2(2) = omegaRho_d(2) * velQuadTerm_y(2)**2 + omegaRho_dm(2)

      outstate( ( ielem-1)* nscalars+ q__s+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__S,2) - sum2_1(2) + sum2_2(2)
      outstate( ( ielem-1)* nscalars+ q__n+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__N,2) + sum2_1(2) + sum2_2(2)

      sum3_1(2) = omegaRho_d(2) * uzstar(2) * div3_4h
      sum3_2(2) = omegaRho_d(2) * velQuadTerm_z(2)**2 + omegaRho_dm(2)

      outstate( ( ielem-1)* nscalars+ q__b+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__B,2) - sum3_1(2) + sum3_2(2)
      outstate( ( ielem-1)* nscalars+ q__t+( 2-1)* qq )         &
        & = omega_o * pdfTmp(Q__T,2) + sum3_1(2) + sum3_2(2)

      !top north / bottom south
      omegaRho_o(2) = omegaRho_d(2) * 0.5_rk
      omegaRho_om(2) = omegaRho_dm(2) * 0.5_rk
      sum4_1(2) = omegaRho_o(2) * ( uystar(2) + uzstar(2) ) *div3_4h
      sum4_2(2) = omegaRho_o(2) * ( velQuadTerm_y(2) + velQuadTerm_z(2) )**2 &
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_bs+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BS,2) - sum4_1(2) + sum4_2(2)
      outstate( ( ielem-1)* nscalars+ q_tn+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TN,2) + sum4_1(2) + sum4_2(2)

      !top south / bottom north
      sum5_1(2) = omegaRho_o(2) * ( -uystar(2) + uzstar(2) ) *div3_4h
      sum5_2(2) = omegaRho_o(2) * ( -velQuadTerm_y(2) + velQuadTerm_z(2) )**2&
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_bn+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BN,2) - sum5_1(2) + sum5_2(2)
      outstate( ( ielem-1)* nscalars+ q_ts+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TS,2) + sum5_1(2) + sum5_2(2)

      !top east / bottom west
      sum6_1(2) = omegaRho_o(2) * ( uxstar(2) + uzstar(2) ) *div3_4h
      sum6_2(2) = omegaRho_o(2) * ( velQuadTerm_x(2) + velQuadTerm_z(2) )**2 &
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_bw+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BW,2) - sum6_1(2) + sum6_2(2)
      outstate( ( ielem-1)* nscalars+ q_te+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TE,2) + sum6_1(2) + sum6_2(2)

      !top west / bottom east
      sum7_1(2) = omegaRho_o(2) * ( -uxstar(2) + uzstar(2) ) *div3_4h
      sum7_2(2) = omegaRho_o(2) * ( -velQuadTerm_x(2) + velQuadTerm_z(2) )**2&
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_be+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BE,2) - sum7_1(2) + sum7_2(2)
      outstate( ( ielem-1)* nscalars+ q_tw+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TW,2) + sum7_1(2) + sum7_2(2)

      !north east / south west
      sum8_1(2) = omegaRho_o(2) * ( uxstar(2) + uystar(2) ) *div3_4h
      sum8_2(2) = omegaRho_o(2) * ( velQuadTerm_x(2) + velQuadTerm_y(2) )**2 &
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_sw+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SW,2) - sum8_1(2) + sum8_2(2)
      outstate( ( ielem-1)* nscalars+ q_ne+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NE,2) + sum8_1(2) + sum8_2(2)

      !north west / south east
      sum9_1(2) = omegaRho_o(2) * ( -uxstar(2) + uystar(2) ) *div3_4h
      sum9_2(2) = omegaRho_o(2) * ( -velQuadTerm_x(2) + velQuadTerm_y(2) )**2&
        &       + omegaRho_om(2)

      outstate( ( ielem-1)* nscalars+ q_se+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SE,2) - sum9_1(2) + sum9_2(2)
      outstate( ( ielem-1)* nscalars+ q_nw+( 2-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NW,2) + sum9_1(2) + sum9_2(2)

      ! species 3
      ! velocity in quadratic term of equilibrium
      velQuadTerm_x(3) = velAvg(1) + theta_eq_spc*uxstar(3)
      velQuadTerm_y(3) = velAvg(2) + theta_eq_spc*uystar(3)
      velQuadTerm_z(3) = velAvg(3) + theta_eq_spc*uzstar(3)

      usqr(3) = (velQuadTerm_x(3) * velQuadTerm_x(3)                           &
        &     + velQuadTerm_y(3) * velQuadTerm_y(3)                            &
        &     + velQuadTerm_z(3) * velQuadTerm_z(3)) * t2cs2inv

      outstate( ( ielem-1)* nscalars+ q__0+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q__0,3) + &
              & div1_3 * omegaRho(3) * (( 3._rk - 2._rk * phi(3)) - usqr(3) )

      omegaRho_dm(3) = omegaRho(3) * (phi(3) - usqr(3)) * div1_18
      !directional velocity factor
      sum1_1(3) = omegaRho_d(3) * uxstar(3) * div3_4h
      sum1_2(3) = omegaRho_d(3) * velQuadTerm_x(3)**2 + omegaRho_dm(3)

      outstate( ( ielem-1)* nscalars+ q__w+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__W,3) - sum1_1(3) + sum1_2(3)
      outstate( ( ielem-1)* nscalars+ q__e+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__E,3) + sum1_1(3) + sum1_2(3)

      sum2_1(3) = omegaRho_d(3) * uystar(3) * div3_4h
      sum2_2(3) = omegaRho_d(3) * velQuadTerm_y(3)**2 + omegaRho_dm(3)

      outstate( ( ielem-1)* nscalars+ q__s+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__S,3) - sum2_1(3) + sum2_2(3)
      outstate( ( ielem-1)* nscalars+ q__n+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__N,3) + sum2_1(3) + sum2_2(3)

      sum3_1(3) = omegaRho_d(3) * uzstar(3) * div3_4h
      sum3_2(3) = omegaRho_d(3) * velQuadTerm_z(3)**2 + omegaRho_dm(3)

      outstate( ( ielem-1)* nscalars+ q__b+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__B,3) - sum3_1(3) + sum3_2(3)
      outstate( ( ielem-1)* nscalars+ q__t+( 3-1)* qq )         &
        & = omega_o * pdfTmp(Q__T,3) + sum3_1(3) + sum3_2(3)

      !top north / bottom south
      omegaRho_o(3) = omegaRho_d(3) * 0.5_rk
      omegaRho_om(3) = omegaRho_dm(3) * 0.5_rk
      sum4_1(3) = omegaRho_o(3) * ( uystar(3) + uzstar(3) ) *div3_4h
      sum4_2(3) = omegaRho_o(3) * ( velQuadTerm_y(3) + velQuadTerm_z(3) )**2 &
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_bs+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BS,3) - sum4_1(3) + sum4_2(3)
      outstate( ( ielem-1)* nscalars+ q_tn+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TN,3) + sum4_1(3) + sum4_2(3)

      !top south / bottom north
      sum5_1(3) = omegaRho_o(3) * ( -uystar(3) + uzstar(3) ) *div3_4h
      sum5_2(3) = omegaRho_o(3) * ( -velQuadTerm_y(3) + velQuadTerm_z(3) )**2&
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_bn+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BN,3) - sum5_1(3) + sum5_2(3)
      outstate( ( ielem-1)* nscalars+ q_ts+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TS,3) + sum5_1(3) + sum5_2(3)

      !top east / bottom west
      sum6_1(3) = omegaRho_o(3) * ( uxstar(3) + uzstar(3) ) *div3_4h
      sum6_2(3) = omegaRho_o(3) * ( velQuadTerm_x(3) + velQuadTerm_z(3) )**2 &
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_bw+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BW,3) - sum6_1(3) + sum6_2(3)
      outstate( ( ielem-1)* nscalars+ q_te+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TE,3) + sum6_1(3) + sum6_2(3)

      !top west / bottom east
      sum7_1(3) = omegaRho_o(3) * ( -uxstar(3) + uzstar(3) ) *div3_4h
      sum7_2(3) = omegaRho_o(3) * ( -velQuadTerm_x(3) + velQuadTerm_z(3) )**2&
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_be+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_BE,3) - sum7_1(3) + sum7_2(3)
      outstate( ( ielem-1)* nscalars+ q_tw+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_TW,3) + sum7_1(3) + sum7_2(3)

      !north east / south west
      sum8_1(3) = omegaRho_o(3) * ( uxstar(3) + uystar(3) ) *div3_4h
      sum8_2(3) = omegaRho_o(3) * ( velQuadTerm_x(3) + velQuadTerm_y(3) )**2 &
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_sw+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SW,3) - sum8_1(3) + sum8_2(3)
      outstate( ( ielem-1)* nscalars+ q_ne+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NE,3) + sum8_1(3) + sum8_2(3)

      !north west / south east
      sum9_1(3) = omegaRho_o(3) * ( -uxstar(3) + uystar(3) ) *div3_4h
      sum9_2(3) = omegaRho_o(3) * ( -velQuadTerm_x(3) + velQuadTerm_y(3) )**2&
        &       + omegaRho_om(3)

      outstate( ( ielem-1)* nscalars+ q_se+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_SE,3) - sum9_1(3) + sum9_2(3)
      outstate( ( ielem-1)* nscalars+ q_nw+( 3-1)* qq ) =       &
              & omega_o*pdfTmp(Q_NW,3) + sum9_1(3) + sum9_2(3)

    end do nodeloop

  end subroutine bgk_advRel_d3q19f3_MSLiquid_WTDF
! ******************************************************************************



! ******************************************************************************
  !> Optimized Advection relaxation routine for the multispecies mrt model
  !! for d3q19 with 3 species with thermodynamic factor
  !!
  !! This routine contains the implementation based on the paper
  !! "A Lattice Boltzmann Scheme for liquid mixtures - Part II: Discretization
  !! and Numerics, Jens Zudrop, Sabine Roller, Pietro Asinari. "\n
  !! MRT paper
  !! "Lattice Boltzmann liquid mixture modeling for electrodialytic engineering
  !!  applications - J. Zudrop, K. Masilamani, S. Roller and P. Asinari"\n
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mrt_advRel_d3q19f3_MSLiquid_WTDF( fieldProp, inState, outState,  &
    &                                          auxField, neigh, nElems,       &
    &                                          nSolve, level, layout, params, &
    &                                          varSys, derVarPos              )
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
    !local variables
    ! temporary local pdf values
    real(kind=rk) :: pdfTmp( QQ, 3 )
    integer :: iElem, iField, nScalars, iField_2, iField_3
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk), dimension(3) :: rsigma, uxsigma,                            &
      & uysigma, uzsigma, qxsigma, qysigma, qzsigma, gqxsigma, gqysigma,       &
      & gqzsigma
    real(kind=rk), dimension(3) :: uxstar, uystar, uzstar
    real(kind=rk), dimension(3) :: molWeight, molWeight_inv, phi,              &
      &                            num_dens, moleFrac
    real(kind=rk) :: usqr, totNum_dens_inv, totMass_dens
    real(kind=rk) :: velAvg_x, velAvg_y, velAvg_z
    real(kind=rk) :: velQuadTerm_x, velQuadTerm_y, velQuadTerm_z
    real(kind=rk) :: omega, omegadiv2, omega_fac, omega_o
    real(kind=rk) :: theta_eq, theta_eq_spc, paramB_inv
    real(kind=rk) :: B_12, B_13, B_23
    real(kind=rk) :: Bratio_12, Bratio_13, Bratio_23
    real(kind=rk) :: mbbEq_12, mbbEq_13, mbbEq_21, mbbEq_23, mbbEq_31, mbbEq_32
    real(kind=rk) :: chi12, chi13, chi21, chi23, chi31, chi32
    real(kind=rk) :: rho_d, rho_dm, rho_o, rho_om, sum1_1, sum2_1, sum3_1,     &
      & sum4_1, sum5_1, sum6_1, sum7_1, sum8_1, sum9_1, sum1_2, sum2_2, sum3_2,&
      & sum4_2, sum5_2, sum6_2, sum7_2, sum8_2, sum9_2
    real(kind=rk) :: omegaMoments(QQ,QQ,3)
    real(kind=rk) :: fneq(QQ), fneq_om(QQ)
    integer :: dens_pos(3), mom_pos(3,3), elemOff
    real(kind=rk), dimension(3, 3) :: matA, invA, resi_coeff, diff_coeff, &
      & thermodynamic_fac, inv_thermodyn_fac
    real(kind=rk) :: temp, press, phy_moleDens_fac
    ! ---------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nScalars = varSys%nScalars

    !KM \todo check moleDens for multilevel
    phy_moleDens_fac = params%physics%moleDens0
    !nFields = scheme%nFields
    !molecular weights
    molWeight = fieldProp(:)%species%molweight
    molWeight_inv = 1.0_rk/molWeight
    !molecular weight ratios
    phi = fieldProp(:)%species%molWeigRatio

    !omega moments = (M^-1 RelaxMat M)*(I+(M^-1 RelMat M)/2)^-1
    omegaMoments(:,:,1) = fieldProp(1)%species%mrt( level )%omegaMoments
    omegaMoments(:,:,2) = fieldProp(2)%species%mrt( level )%omegaMoments
    omegaMoments(:,:,3) = fieldProp(3)%species%mrt( level )%omegaMoments

    !omega
    omega = scheme%mixture%relaxLvl(Level)%omega_diff
    omega_fac = (1.0_rk/(1.0_rk/omega + 0.5_rk))
    omegadiv2 = omega*0.5_rk
    omega_o = 1.0_rk - omega_fac

    !free parameter B
    paramB_inv = 1.0_rk/scheme%mixture%paramB

    !equilibrium theta
    theta_eq = scheme%mixture%theta_eq
    theta_eq_spc = 1.0_rk - theta_eq

    ! temperature
    temp = scheme%mixture%temp0
    ! atmospheric pressure
    press = scheme%mixture%atm_press

    !resistivities
    !resistivity coeff are symmetric
    !B_12 = B_21, B_13=B31, B23=B32
    B_12 = fieldProp( 1 )%species%resi_coeff(2)
    B_13 = fieldProp( 1 )%species%resi_coeff(3)
    B_23 = fieldProp( 2 )%species%resi_coeff(3)

    Bratio_12 = B_12*paramB_inv
    Bratio_13 = B_13*paramB_inv
    Bratio_23 = B_23*paramB_inv
    mbbEq_12 = Bratio_12*phi(1)
    mbbEq_13 = Bratio_13*phi(1)

    mbbEq_21 = Bratio_12*phi(2)
    mbbEq_23 = Bratio_23*phi(2)

    mbbEq_31 = Bratio_13*phi(3)
    mbbEq_32 = Bratio_23*phi(3)

    ! position of density and momentum in auxField array
    do iField = 1, 3
      dens_pos(iField) = varSys%method%val(derVarPos(iField)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:, iField) = varSys%method%val(derVarPos(iField)%momentum) &
        &                                  %auxField_varPos(1:3)
    end do

    nodeloop: do iElem = 1, nSolve

      !species 1
      pdfTmp( Q__W,1 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__S,1 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__B,1 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__E,1 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__N,1 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__T,1 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,1 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,1 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,1 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,1 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,1 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,1 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,1 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,1 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,1 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,1 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,1 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,1 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      pdfTmp( Q__0,1 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)
      !species 2
      pdfTmp( Q__W,2 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__S,2 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__B,2 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__E,2 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__N,2 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__T,2 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,2 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,2 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,2 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,2 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,2 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,2 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,2 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,2 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,2 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,2 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,2 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,2 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      pdfTmp( Q__0,2 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 2-1)* qq+ nscalars*0)
      !species 3
      pdfTmp( Q__W,3 ) = instate(neigh (( q__w-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__S,3 ) = instate(neigh (( q__s-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__B,3 ) = instate(neigh (( q__b-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__E,3 ) = instate(neigh (( q__e-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__N,3 ) = instate(neigh (( q__n-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__T,3 ) = instate(neigh (( q__t-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BS,3 ) = instate(neigh (( q_bs-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TS,3 ) = instate(neigh (( q_ts-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BN,3 ) = instate(neigh (( q_bn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TN,3 ) = instate(neigh (( q_tn-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BW,3 ) = instate(neigh (( q_bw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_BE,3 ) = instate(neigh (( q_be-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TW,3 ) = instate(neigh (( q_tw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_TE,3 ) = instate(neigh (( q_te-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SW,3 ) = instate(neigh (( q_sw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NW,3 ) = instate(neigh (( q_nw-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_SE,3 ) = instate(neigh (( q_se-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q_NE,3 ) = instate(neigh (( q_ne-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)
      pdfTmp( Q__0,3 ) = instate(neigh (( q__0-1)* nelems+ ielem)+( 3-1)* qq+ nscalars*0)

      ! element offset to access auxField
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! get density and momentum from auxField
      do iField = 1, 3
        ! density
        rsigma(iField) = auxField(elemOff + dens_pos(iField))
        ! momentum
        gqxsigma(iField) = auxField(elemOff + mom_pos(1, iField))
        gqysigma(iField) = auxField(elemOff + mom_pos(2, iField))
        gqzsigma(iField) = auxField(elemOff + mom_pos(3, iField))
      end do

      !number density
      num_dens(:) = rsigma(:)*molWeight_inv(:)

      !total number density
      totNum_dens_inv = 1.0_rk/sum(num_dens(:))

      !mole fraction
      moleFrac(:) = num_dens(:)*totNum_dens_inv

      ! MS-Diff coeff matrix from C++ code
      call mus_calc_MS_DiffMatrix( 3, temp, press, num_dens*phy_moleDens_fac, &
        &                          diff_coeff )

      ! Convert to lattice unit
      resi_coeff = params%physics%fac(level)%diffusivity/diff_coeff

      ! Thermodynamic factor from C++ code
      call mus_calc_thermFactor( 3, temp, press, moleFrac, &
        &                        thermodynamic_fac         )

      inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

      matA = 0.0_rk
      !build up matrix to solver LSE for actual velocity
      do iField = 1, 3
        !set diagonal part
        matA(iField, iField) = 1.0_rk
        do iField_2 = 1, 3
          do iField_3 = 1, 3
            matA(iField, iField_2) = matA(iField, iField_2) + omega * 0.5_rk   &
              &                    * inv_thermodyn_fac(iField, iField_2)       &
              &                    * resi_coeff(iField_2, iField_3)            &
              &                    * phi(iField_2) * moleFrac(iField_3)        &
              &                    * paramB_inv
          end do
        end do
        !set non-diagonal part
        do iField_2 = 1, 3
          do iField_3 = 1, 3
            matA(iField, iField_3) = matA(iField, iField_3) - omega * 0.5_rk   &
              &                    * inv_thermodyn_fac(iField, iField_2)       &
              &                    * resi_coeff(iField_2, iField_3)            &
              &                    * phi(iField_3) * moleFrac(iField_2)        &
              &                    * paramB_inv
          end do
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matA )

      !actual momentum of all species
      qxsigma(1:3) = matmul( invA, gqxsigma(1:3) )
      qysigma(1:3) = matmul( invA, gqysigma(1:3) )
      qzsigma(1:3) = matmul( invA, gqzsigma(1:3) )

      ! store momentum of untransformed PDF in auxField
      do iField = 1, 3
        auxField(elemOff + mom_pos(1, iField)) = qxsigma(iField)
        auxField(elemOff + mom_pos(2, iField)) = qysigma(iField)
        auxField(elemOff + mom_pos(3, iField)) = qzsigma(iField)
      end do

      ! species velocity
      uxsigma = qxsigma/rsigma
      uysigma = qysigma/rsigma
      uzsigma = qzsigma/rsigma

      !compute equilibrium velocity to compute feq
      chi12 = mbbEq_12 * moleFrac(2)
      chi13 = mbbEq_13 * moleFrac(3)
      chi21 = mbbEq_21 * moleFrac(1)
      chi23 = mbbEq_23 * moleFrac(3)
      chi31 = mbbEq_31 * moleFrac(1)
      chi32 = mbbEq_32 * moleFrac(2)
      !ux
      uxstar(1) =  uxsigma(1) + chi12*(uxsigma(2)-uxsigma(1)) &
        & + chi13*(uxsigma(3)-uxsigma(1))

      uxstar(2) =  uxsigma(2) + chi21*(uxsigma(1)-uxsigma(2)) &
        & + chi23*(uxsigma(3)-uxsigma(2))

      uxstar(3) =  uxsigma(3) + chi31*(uxsigma(1)-uxsigma(3)) &
        & + chi32*(uxsigma(2)-uxsigma(3))
      !uy
      uystar(1) =  uysigma(1) + chi12*(uysigma(2)-uysigma(1)) &
        & + chi13*(uysigma(3)-uysigma(1))

      uystar(2) =  uysigma(2) + chi21*(uysigma(1)-uysigma(2)) &
        & + chi23*(uysigma(3)-uysigma(2))

      uystar(3) =  uysigma(3) + chi31*(uysigma(1)-uysigma(3)) &
        & + chi32*(uysigma(2)-uysigma(3))
      !uz
      uzstar(1) =  uzsigma(1) + chi12*(uzsigma(2)-uzsigma(1)) &
        & + chi13*(uzsigma(3)-uzsigma(1))

      uzstar(2) =  uzsigma(2) + chi21*(uzsigma(1)-uzsigma(2)) &
        & + chi23*(uzsigma(3)-uzsigma(2))

      uzstar(3) =  uzsigma(3) + chi31*(uzsigma(1)-uzsigma(3)) &
        & + chi32*(uzsigma(2)-uzsigma(3))

      ! total mass density
      totMass_dens = sum(rsigma)

      ! mass averaged mixture velocity
      velAvg_x = dot_product( rsigma, uxsigma ) / totMass_dens
      velAvg_y = dot_product( rsigma, uysigma ) / totMass_dens
      velAvg_z = dot_product( rsigma, uzsigma ) / totMass_dens

      !compute equilibrium and do collision
      do iField = 1, 3!nFields
        ! \todo KM: Bug remove omega factor from equilibrium
        rho_d = rsigma(iField) * div1_4

        velQuadTerm_x = theta_eq*velAvg_x + theta_eq_spc*uxstar(iField)
        velQuadTerm_y = theta_eq*velAvg_y + theta_eq_spc*uystar(iField)
        velQuadTerm_z = theta_eq*velAvg_z + theta_eq_spc*uzstar(iField)

        usqr = (velQuadTerm_x * velQuadTerm_x                                  &
          &  + velQuadTerm_y * velQuadTerm_y                                   &
          &  + velQuadTerm_z * velQuadTerm_z) * t2cs2inv

        ! at rest
        fneq(Q__0) = div1_3 * rsigma(iField) * (( 3._rk - 2._rk * phi(iField)) &
          &        - usqr ) - pdfTmp(Q__0, iField)

        rho_dm = rsigma(iField) * (phi(iField) - usqr) * div1_18
        !directional velocity factor
        sum1_1 = rho_d * uxstar(iField) * div3_4h
        sum1_2 = rho_d * velQuadTerm_x**2 + rho_dm
        fneq(Q__W) = -sum1_1 + sum1_2 - pdfTmp(Q__W,iField)
        fneq(Q__E) = sum1_1 + sum1_2 - pdfTmp(Q__E,iField)

        sum2_1 = rho_d * uystar(iField) * div3_4h
        sum2_2 = rho_d * velQuadTerm_y**2 + rho_dm
        fneq(Q__S) = -sum2_1 + sum2_2 - pdfTmp(Q__S,iField)
        fneq(Q__N) = sum2_1 + sum2_2 - pdfTmp(Q__N,iField)

        sum3_1 = rho_d * uzstar(iField) * div3_4h
        sum3_2 = rho_d * velQuadTerm_z**2 + rho_dm
        fneq(Q__B) = -sum3_1 + sum3_2 - pdfTmp(Q__B,iField)
        fneq(Q__T) = sum3_1 + sum3_2 - pdfTmp(Q__T,iField)

        rho_o = rho_d * 0.5_rk
        rho_om = rho_dm * 0.5_rk
        !top north / bottom south
        sum4_1 = rho_o * ( uystar(iField) + uzstar(iField) ) *div3_4h
        sum4_2 = rho_o * ( velQuadTerm_y + velQuadTerm_z )**2 &
        &       + rho_om
        fneq(Q_BS) = -sum4_1 + sum4_2 - pdfTmp(Q_BS,iField)
        fneq(Q_TN) = sum4_1 + sum4_2 - pdfTmp(Q_TN,iField)

        !top south / bottom north
        sum5_1 = rho_o * ( -uystar(iField) + uzstar(iField) ) *div3_4h
        sum5_2 = rho_o * ( -velQuadTerm_y + velQuadTerm_z )**2&
          &       + rho_om
        fneq(Q_BN) = -sum5_1 + sum5_2 - pdfTmp(Q_BN,iField)
        fneq(Q_TS) = sum5_1 + sum5_2 - pdfTmp(Q_TS,iField)

        !top east / bottom west
        sum6_1 = rho_o * ( uxstar(iField) + uzstar(iField) ) *div3_4h
        sum6_2 = rho_o * ( velQuadTerm_x + velQuadTerm_z )**2 &
          &       + rho_om
        fneq(Q_BW) = -sum6_1 + sum6_2 - pdfTmp(Q_BW,iField)
        fneq(Q_TE) = sum6_1 + sum6_2 - pdfTmp(Q_TE,iField)

        !top west / bottom east
        sum7_1 = rho_o * ( -uxstar(iField) + uzstar(iField) ) *div3_4h
        sum7_2 = rho_o * ( -velQuadTerm_x + velQuadTerm_z )**2&
          &       + rho_om
        fneq(Q_BE) = -sum7_1 + sum7_2 - pdfTmp(Q_BE,iField)
        fneq(Q_TW) = sum7_1 + sum7_2 - pdfTmp(Q_TW,iField)

        !north east / south west
        sum8_1 = rho_o * ( uxstar(iField) + uystar(iField) ) *div3_4h
        sum8_2 = rho_o * ( velQuadTerm_x + velQuadTerm_y )**2 &
          &       + rho_om
        fneq(Q_SW) = -sum8_1 + sum8_2 - pdfTmp(Q_SW,iField)
        fneq(Q_NE) = sum8_1 + sum8_2 - pdfTmp(Q_NE,iField)

        !north west / south east
        sum9_1 = rho_o * ( -uxstar(iField) + uystar(iField) ) *div3_4h
        sum9_2 = rho_o * ( -velQuadTerm_x + velQuadTerm_y )**2&
          &       + rho_om
        fneq(Q_SE) = -sum9_1 + sum9_2 - pdfTmp(Q_SE,iField)
        fneq(Q_NW) = sum9_1 + sum9_2 - pdfTmp(Q_NW,iField)

        ! omega * fNonEq
        fneq_om = matmul( omegaMoments(:,:, iField), fneq )

        !outstate
        outstate( ( ielem-1)* nscalars+ q__0+( ifield-1)* qq ) =       &
          & pdfTmp( Q__0, iField ) + fneq_om( Q__0 )

        outstate( ( ielem-1)* nscalars+ q__w+( ifield-1)* qq ) =       &
          & pdfTmp( Q__W, iField ) + fneq_om( Q__W )
        outstate( ( ielem-1)* nscalars+ q__e+( ifield-1)* qq ) =       &
          & pdfTmp( Q__E, iField ) + fneq_om( Q__E )
        outstate( ( ielem-1)* nscalars+ q__s+( ifield-1)* qq ) =       &
          & pdfTmp( Q__S, iField ) + fneq_om( Q__S )
        outstate( ( ielem-1)* nscalars+ q__n+( ifield-1)* qq ) =       &
          & pdfTmp( Q__N, iField ) + fneq_om( Q__N )
        outstate( ( ielem-1)* nscalars+ q__b+( ifield-1)* qq ) =       &
          & pdfTmp( Q__B, iField ) + fneq_om( Q__B )
        outstate( ( ielem-1)* nscalars+ q__t+( ifield-1)* qq ) =       &
          & pdfTmp( Q__T, iField ) + fneq_om( Q__T )
        outstate( ( ielem-1)* nscalars+ q_bs+( ifield-1)* qq ) =       &
          & pdfTmp( Q_BS, iField ) + fneq_om( Q_BS )
        outstate( ( ielem-1)* nscalars+ q_tn+( ifield-1)* qq ) =       &
          & pdfTmp( Q_TN, iField ) + fneq_om( Q_TN )
        outstate( ( ielem-1)* nscalars+ q_bn+( ifield-1)* qq ) =       &
          & pdfTmp( Q_BN, iField ) + fneq_om( Q_BN )
        outstate( ( ielem-1)* nscalars+ q_ts+( ifield-1)* qq ) =       &
          & pdfTmp( Q_TS, iField ) + fneq_om( Q_TS )
        outstate( ( ielem-1)* nscalars+ q_bw+( ifield-1)* qq ) =       &
          & pdfTmp( Q_BW, iField ) + fneq_om( Q_BW )
        outstate( ( ielem-1)* nscalars+ q_te+( ifield-1)* qq ) =       &
          & pdfTmp( Q_TE, iField ) + fneq_om( Q_TE )
        outstate( ( ielem-1)* nscalars+ q_be+( ifield-1)* qq ) =       &
          & pdfTmp( Q_BE, iField ) + fneq_om( Q_BE )
        outstate( ( ielem-1)* nscalars+ q_tw+( ifield-1)* qq ) =       &
          & pdfTmp( Q_TW, iField ) + fneq_om( Q_TW )
        outstate( ( ielem-1)* nscalars+ q_sw+( ifield-1)* qq ) =       &
          & pdfTmp( Q_SW, iField ) + fneq_om( Q_SW )
        outstate( ( ielem-1)* nscalars+ q_ne+( ifield-1)* qq ) =       &
          & pdfTmp( Q_NE, iField ) + fneq_om( Q_NE )
        outstate( ( ielem-1)* nscalars+ q_se+( ifield-1)* qq ) =       &
          & pdfTmp( Q_SE, iField ) + fneq_om( Q_SE )
        outstate( ( ielem-1)* nscalars+ q_nw+( ifield-1)* qq ) =       &
          & pdfTmp( Q_NW, iField ) + fneq_om( Q_NW )

      enddo !iField
    enddo nodeloop

  end subroutine mrt_advRel_d3q19f3_MSLiquid_WTDF
! ****************************************************************************** !


! ******************************************************************************
  !> Unoptimized Advection relaxation routine for the multispecies BGK model
  !! with thermodynamic factors in Maxwell-Stefan formulation
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mrt_advRel_MSLiquid_generic_WTDF( fieldProp, inState, outState,  &
    &                                          auxField, neigh, nElems,       &
    &                                          nSolve, level, layout, params, &
    &                                          varSys, derVarPos              )
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
    !local variables
    real(kind=rk) :: pdfTmp( varSys%nScalars )  !< temporary local pdf values
    integer :: iElem, iField, iField_2, iField_3, nFields, iDir, QQ, nScalars
    integer :: vPos( layout%fStencil%QQ )
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer, dimension(varSys%nStateVars) :: stateVarMap, dens_pos
    integer :: mom_pos(3, varSys%nStateVars)
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, num_dens, &
      & moleFrac, molWeight, phi
    real(kind=rk), dimension(3, varSys%nStateVars) :: first_moments, &
      & velocity, eqVel
    real(kind=rk) :: totNum_dens, usqr, ucx, feq
    real(kind=rk) :: omega, omega_new, paramB, theta_eq
    real(kind=rk) :: velAvg(3), velNew(3), ucxQuadTerm, totMassDens
    real(kind=rk) :: omegaMoments(layout%fStencil%QQ,layout%fStencil%QQ)
    real(kind=rk) :: fneq(layout%fStencil%QQ), fneq_om(layout%fStencil%QQ)
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: matA, &
      & invA, resi_coeff, thermodynamic_fac, inv_thermodyn_fac, diff_coeff
    real(kind=rk) :: temp, press, moleDens0, phy_moleDens_fac
    integer :: elemOff
    real(kind=rk) :: wRestInv, sRest(varSys%nStateVars)
    ! ------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    QQ     = layout%fStencil%QQ
    nFields = scheme%nFields
    nScalars = varSys%nScalars

    ! constant mixture number density
    moleDens0 = scheme%mixture%moleDens0

    !KM \todo check moleDens for multilevel
    phy_moleDens_fac = params%physics%moleDens0

    !molecular weights
    molWeight = fieldProp(:)%species%molweight
    !molecular weight ratios
    phi = fieldProp(:)%species%molWeigRatio

    !s_0^k - species weight for rest position
    wRestInv = 1.0_rk / layout%weight(layout%fStencil%restPosition)
    do iField = 1, nFields
      sRest(iField) = wRestInv + (1.0_rk - wRestInv) * phi(iField)
    end do

    !resistivities
    do iField = 1, nFields
      resi_coeff(iField,:) = fieldProp( iField )%species%resi_coeff
    enddo

    omega = scheme%mixture%relaxLvl(level)%omega_diff
    omega_new = (1.0_rk/(1.0_rk/omega + 0.5_rk))

    paramB = scheme%mixture%paramB
    !equilibrium theta
    theta_eq = scheme%mixture%theta_eq

    ! temperature
    temp = scheme%mixture%temp0
    ! atmospheric pressure
    press = scheme%mixture%atm_press

    stateVarMap = scheme%stateVarMap%varPos%val(:)

    ! compute thermodynamic factor from matthias c-code for each element
    thermodynamic_fac = 0.0_rk
    do iField = 1, nFields
      thermodynamic_fac(iField, iField) = 1.0_rk
    end do

    ! position of density and momentum in auxField array
    do iField = 1, nFields
      dens_pos(iField) = varSys%method%val(derVarPos(iField)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:, iField) = varSys%method%val(derVarPos(iField)%momentum) &
        &                                  %auxField_varPos(1:3)
    end do

    nodeloop: do iElem = 1, nSolve

      mass_dens = 0._rk
      first_moments = 0.0_rk
      velocity = 0.0_rk
      eqVel = 0.0_rk
      pdfTmp = 0._rk

      ! element offset to access auxField
      elemOff = (iElem-1)*varSys%nAuxScalars

      do iField = 1, nFields
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos
        ! compute field density and first moments
        do iDir = 1, QQ
          !store all field pdf in single element to local array pdfTmp
          pdfTmp( vPos(iDir) ) =                                               &
& instate( neigh (( idir-1)* nelems+ ielem)+( ifield-1)* qq+ varsys%nscalars*0 )
        enddo
        !field density
        mass_dens( iField ) = auxField(elemOff + dens_pos(iField))

        !field momentum (rho*u)
        first_moments( 1, iField ) = auxField(elemOff + mom_pos(1, iField))
        first_moments( 2, iField ) = auxField(elemOff + mom_pos(2, iField))
        first_moments( 3, iField ) = auxField(elemOff + mom_pos(3, iField))
      enddo

      !total mass density
      totmassDens = sum(mass_dens)

      ! solve linear system of equation for actual moments
      ! number density of all species
      num_dens(:) = mass_dens(:)/molWeight(:)

      !total number density
      totNum_dens = sum(num_dens(:))

      !mole fraction
      moleFrac(:) =  num_dens(:)/totNum_dens

      ! MS-Diff coeff matrix from C++ code
      call mus_calc_MS_DiffMatrix( nFields, temp, press,                       &
        &                          num_dens*phy_moleDens_fac, diff_coeff )

      ! Convert to lattice unit
      resi_coeff = params%physics%fac(level)%diffusivity/diff_coeff

      ! Thermodynamic factor from C++ code
      call mus_calc_thermFactor( nFields, temp, press, moleFrac,               &
        &                        thermodynamic_fac )

      inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

      matA = 0.0_rk
      !build up matrix to solver LSE for actual velocity
      do iField = 1, nFields
        !set diagonal part
        matA(iField, iField) = 1.0_rk
        do iField_2 = 1, nFields
          do iField_3 = 1, nFields
            matA(iField, iField_2) = matA(iField, iField_2) + omega * 0.5_rk   &
              &                    * inv_thermodyn_fac(iField, iField_2)       &
              &                    * resi_coeff(iField_2, iField_3)            &
              &                    * phi(iField_2) * moleFrac(iField_3)        &
              &                    / paramB
          end do
        end do
        !set non-diagonal part
        do iField_2 = 1, nFields
          do iField_3 = 1, nFields
            matA(iField, iField_3) = matA(iField, iField_3) - omega * 0.5_rk   &
              &                    * inv_thermodyn_fac(iField, iField_2)       &
              &                    * resi_coeff(iField_2, iField_3)            &
              &                    * phi(iField_3) * moleFrac(iField_2)        &
              &                    / paramB
          end do
        end do
      end do


      ! invert matrix
      invA = invert_matrix( matA )

      !actual momentum of all species
      velocity(1, :) = matmul( invA, first_moments(1,:) )
      velocity(2, :) = matmul( invA, first_moments(2,:) )
      velocity(3, :) = matmul( invA, first_moments(3,:) )

      ! store momentum of untransformed PDF in auxField
      do iField = 1, nFields
        auxField(elemOff + mom_pos(1, iField)) = velocity(1, iField)
        auxField(elemOff + mom_pos(2, iField)) = velocity(2, iField)
        auxField(elemOff + mom_pos(3, iField)) = velocity(3, iField)
      end do

      ! convert momentum to velocity
      velocity(1, :) = velocity(1, :) / mass_dens(:)
      velocity(2, :) = velocity(2, :) / mass_dens(:)
      velocity(3, :) = velocity(3, :) / mass_dens(:)

      ! compute equilibrium velocity
      do iField = 1, nFields
        eqVel( :, iField ) = mass_dens(iField)*velocity( :, iField )
        do iField_2 = 1, nFields
          do iField_3 = 1, nFields
            eqVel( :, iField ) = eqVel( :, iField )                              &
              &                + inv_thermodyn_fac(iField, iField_2)             &
              &                * mass_dens(iField_2)                             &
              &                * resi_coeff( iField_2, iField_3 ) * phi(iField_2)&
              &                * moleFrac(iField_3)                              &
              &                * (velocity(:, iField_3) - velocity(:,iField_2))  &
              &                / paramB
          end do
        end do
      end do

      !compute mass averaged mixture velocity
      velAvg(1) = dot_product( mass_dens, velocity(1,:) )/totmassDens
      velAvg(2) = dot_product( mass_dens, velocity(2,:) )/totmassDens
      velAvg(3) = dot_product( mass_dens, velocity(3,:) )/totmassDens

      !compute equilibrium and do collision
      do iField = 1, nFields
        !omega moments = (M^-1 RelaxMat M)*(I+(M^-1 RelMat M)/2)^-1
        omegaMoments = fieldProp( iField )%species%mrt( level )%omegaMoments
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos

        !KM: TGV testcase
        velNew = theta_eq*velAvg                                        &
          &    + (1.0_rk - theta_eq)*eqVel(:, iField) / mass_dens(iField)

        !usqr = eqVel(1, iField) * eqVel(1, iField) &
        !  &  + eqVel(2, iField) * eqVel(2, iField) &
        !  &  + eqVel(3, iField) * eqVel(3, iField)
        usqr = dot_product( velNew, velNew ) * t2cs2inv

        do iDir = 1, QQ
          ucx = dble( layout%fStencil%cxDir( 1, iDir ))*eqVel(1, iField) &
            & + dble( layout%fStencil%cxDir( 2, iDir ))*eqVel(2, iField) &
            & + dble( layout%fStencil%cxDir( 3, iDir ))*eqVel(3, iField)

          ucxQuadTerm = dot_product(                                 &
            &           dble(layout%fStencil%cxDir(:, iDir)), velNew )

          feq = layout%weight(iDir) * ( mass_dens(iField) * ( phi(iField)    &
            & + ucxQuadTerm * ucxQuadTerm * t2cs4inv - usqr ) + ucx * cs2inv )

          if ( iDir == layout%fStencil%restPosition ) then
            ! equilibrium at rest
            feq = layout%weight(iDir) * mass_dens(iField) * ( sRest(iField) &
              &                                              - usqr         )
          end if
          fneq(iDir) = ( feq - pdfTmp( vPos(iDir) ) )
        enddo !iDir

        fneq_om = matmul( omegaMoments, fneq )

        do iDir = 1, layout%fStencil%QQ
          outstate(                                                    &
& ( ielem-1)* varsys%nscalars+ idir+( ifield-1)* qq ) = &
            & pdfTmp( vPos(iDir) ) + fneq_om( iDir )
        enddo !iDir

      enddo !iField
    enddo nodeloop

  end subroutine mrt_advRel_MSLiquid_generic_WTDF
! ******************************************************************************


! ******************************************************************************
  !> Unoptimized Advection relaxation routine for the multispecies BGK model
  !! with thermodynamic factors in Maxwell-Stefan formulation
  !!
  !! This routine contains the implementation based on the paper
  !! "A Lattice Boltzmann Scheme for liquid mixtures - Part II: Discretization
  !! and Numerics, Jens Zudrop, Sabine Roller, Pietro Asinari. "\n
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine bgk_advRel_MSLiquid_generic_WTDF( fieldProp, inState, outState,  &
    &                                          auxField, neigh, nElems,       &
    &                                          nSolve, level, layout, params, &
    &                                          varSys, derVarPos              )
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
    !local variables
    real(kind=rk) :: pdfTmp( varSys%nScalars )  !< temporary local pdf values
    integer :: iElem, iField, iField_2, iField_3, nFields, iDir, QQ, nScalars
    integer :: vPos( layout%fStencil%QQ )
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer, dimension(varSys%nStateVars) :: stateVarMap, dens_pos
    integer :: mom_pos(3, varSys%nStateVars)
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, num_dens, &
      & moleFrac, molWeight, phi
    real(kind=rk), dimension(3, varSys%nStateVars) :: first_moments, &
      & velocity, eqVel
    real(kind=rk) :: totNum_dens, usqr, ucx, feq
    real(kind=rk) :: velAvg(3), velQuad(3), ucxQuadTerm, theta_eq, totMassDens
    real(kind=rk) :: omega, omega_new, paramB
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: matA, &
      & invA, resi_coeff, thermodynamic_fac, inv_thermodyn_fac, diff_coeff
    real(kind=rk) :: temp, press, moleDens0, phy_moleDens_fac
    integer :: elemOff
    real(kind=rk) :: wRestInv, sRest(varSys%nStateVars)
    ! ---------------------------------------------------------------------------
    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    QQ     = layout%fStencil%QQ
    nFields = scheme%nFields
    nScalars = varSys%nScalars

    ! constant mixture number density
    moleDens0 = scheme%mixture%moleDens0

    !KM \todo check moleDens for multilevel
    phy_moleDens_fac = params%physics%moleDens0

    !molecular weights
    molWeight = fieldProp(:)%species%molweight
    !molecular weight ratios
    phi = fieldProp(:)%species%molWeigRatio

    !s_0^k - species weight for rest position
    wRestInv = 1.0_rk / layout%weight(layout%fStencil%restPosition)
    do iField = 1, nFields
      sRest(iField) = wRestInv + (1.0_rk - wRestInv) * phi(iField)
    end do

    !resistivities
    do iField = 1, nFields
      resi_coeff(iField,:) = fieldProp( iField )%species%resi_coeff
    enddo

    omega = scheme%mixture%relaxLvl(Level)%omega_diff
!    write(*,*) 'omega ', omega
    omega_new = (1.0_rk/(1.0_rk/omega + 0.5_rk))
!    write(*,*) 'omega_new ', omega_new

    paramB = scheme%mixture%paramB
    !equilibrium theta
    theta_eq = scheme%mixture%theta_eq

    ! temperature
    temp = scheme%mixture%temp0
    ! atmospheric pressure
    press = scheme%mixture%atm_press

!write(*,*) 'resi_coeff LB', real(resi_coeff)
!write(*,*) 'temp ', temp, ' press ', press
    ! initialize thermodynamic factor and activity coeff
    ! compute thermodynamic factor from matthias c-code for each element
    thermodynamic_fac = 0.0_rk
    do iField = 1, nFields
      thermodynamic_fac(iField, iField) = 1.0_rk
    end do
!    write(*,*) 'paramB ', paramB
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    ! position of density and momentum in auxField array
    do iField = 1, nFields
      dens_pos(iField) = varSys%method%val(derVarPos(iField)%density) &
        &                           %auxField_varPos(1)
      mom_pos(:, iField) = varSys%method%val(derVarPos(iField)%momentum) &
        &                                  %auxField_varPos(1:3)
    end do

    nodeloop: do iElem = 1, nSolve

      mass_dens = 0._rk
      first_moments = 0.0_rk
      velocity = 0.0_rk
      eqVel = 0.0_rk
      pdfTmp = 0._rk

      ! element offset to access auxField
      elemOff = (iElem-1)*varSys%nAuxScalars

      do iField = 1, nFields
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos
        ! compute field density and first moments
        do iDir = 1, layout%fStencil%QQ
          !store all field pdf in single element to local array pdfTmp
          pdfTmp( vPos(iDir) ) = instate(                               &
            &  neigh((idir-1)* nelems+ ielem)+( ifield-1)* qq+ nscalars*0 )
        enddo
        !field density
        mass_dens( iField ) = auxField(elemOff + dens_pos(iField))

        !field momentum (rho*u)
        first_moments( 1, iField ) = auxField(elemOff + mom_pos(1, iField))
        first_moments( 2, iField ) = auxField(elemOff + mom_pos(2, iField))
        first_moments( 3, iField ) = auxField(elemOff + mom_pos(3, iField))
      enddo

      !total mass density
      totmassDens = sum(mass_dens)

      ! solve linear system of equation for actual moments
      ! number density of all species
      num_dens(:) = mass_dens(:)/molWeight(:)

      !total number density
      totNum_dens = sum(num_dens(:))

      !mole fraction
      moleFrac(:) =  num_dens(:)/totNum_dens

      ! MS-Diff coeff matrix from C++ code
      call mus_calc_MS_DiffMatrix( nFields, temp, press,                       &
        &                        num_dens*phy_moleDens_fac, diff_coeff )

      ! Convert to lattice unit
      resi_coeff = params%physics%fac(level)%diffusivity/diff_coeff

      ! Thermodynamic factor from C++ code
      call mus_calc_thermFactor( nFields, temp, press, moleFrac,               &
        &                        thermodynamic_fac )

      inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

      matA = 0.0_rk
      !build up matrix to solver LSE for actual velocity
      do iField = 1, nFields
        !set diagonal part
        matA(iField, iField) = 1.0_rk
        do iField_2 = 1, nFields
          do iField_3 = 1, nFields
            matA(iField, iField_2) = matA(iField, iField_2) + omega * 0.5_rk   &
              &                    * inv_thermodyn_fac(iField, iField_2)       &
              &                    * resi_coeff(iField_2, iField_3)            &
              &                    * phi(iField_2) * moleFrac(iField_3)        &
              &                    / paramB
          end do
        end do
        !set non-diagonal part
        do iField_2 = 1, nFields
          do iField_3 = 1, nFields
            matA(iField, iField_3) = matA(iField, iField_3) - omega * 0.5_rk   &
              &                    * inv_thermodyn_fac(iField, iField_2)       &
              &                    * resi_coeff(iField_2, iField_3)            &
              &                    * phi(iField_3) * moleFrac(iField_2)        &
              &                    / paramB
          end do
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matA )

      !actual momentum of all species
      velocity(1, :) = matmul( invA, first_moments(1,:) )
      velocity(2, :) = matmul( invA, first_moments(2,:) )
      velocity(3, :) = matmul( invA, first_moments(3,:) )

      ! store momentum of untransformed PDF in auxField
      do iField = 1, nFields
        auxField(elemOff + mom_pos(1, iField)) = velocity(1, iField)
        auxField(elemOff + mom_pos(2, iField)) = velocity(2, iField)
        auxField(elemOff + mom_pos(3, iField)) = velocity(3, iField)
      end do

      ! convert momentum to velocity
      velocity(1, :) = velocity(1, :) / mass_dens(:)
      velocity(2, :) = velocity(2, :) / mass_dens(:)
      velocity(3, :) = velocity(3, :) / mass_dens(:)

      ! compute equilibrim velocity
      do iField = 1, nFields
        eqVel( :, iField ) = mass_dens(iField)*velocity( :, iField )
        do iField_2 = 1, nFields
          do iField_3 = 1, nFields
            eqVel( :, iField ) = eqVel( :, iField )                              &
              &                + inv_thermodyn_fac(iField, iField_2)             &
              &                * mass_dens(iField_2)                             &
              &                * resi_coeff( iField_2, iField_3 ) * phi(iField_2)&
              &                * moleFrac(iField_3)                              &
              &                * (velocity(:, iField_3) - velocity(:,iField_2))  &
              &                / paramB
          end do
        end do
      end do

      !compute mass averaged mixture velocity
      velAvg(1) = dot_product( mass_dens, velocity(1,:) )/totmassDens
      velAvg(2) = dot_product( mass_dens, velocity(2,:) )/totmassDens
      velAvg(3) = dot_product( mass_dens, velocity(3,:) )/totmassDens

      !compute equilibrium and do collision
      do iField = 1, nFields
        vPos = varSys%method%val( stateVarMap(iField) )%state_varPos
        !KM: TGV testcase
        velQuad = theta_eq*velAvg                                              &
          &     + (1.0_rk - theta_eq) * eqVel(:, iField) / mass_dens(iField)

        usqr = dot_product( velQuad, velQuad ) * t2cs2inv
!write(dbgUnit(1),*) 'omega_new ', omega_new
        do iDir = 1, layout%fStencil%QQ
          ucx = dot_product(                                                   &
            & dble(layout%fStencil%cxDir(:, iDir)), eqVel(:, iField) )

          ucxQuadTerm = dot_product(                                           &
            & dble(layout%fStencil%cxDir(:, iDir)), velQuad )

          !feq = layout%weight(iDir) * mass_dens(iField) * ( phi(iField)        &
          !  & + ucx * cs2inv + ucx * ucx * t2cs4inv - usqr )
          feq = layout%weight(iDir) * ( mass_dens(iField) * ( phi(iField)    &
            & + ucxQuadTerm * ucxQuadTerm * t2cs4inv - usqr ) + ucx * cs2inv )

          if ( iDir == layout%fStencil%restPosition ) then
            ! equilibrium at rest
            feq = layout%weight(iDir) * mass_dens(iField) * ( sRest(iField) &
              &                                              - usqr         )
          end if

          outstate( ( ielem-1)* nscalars+idir+( ifield-1)* qq ) &
            & = pdfTmp(vPos(iDir)) + omega_new * ( feq - pdfTmp(vPos(iDir)) )
        enddo !iDir

      enddo !iField
    enddo nodeloop

  end subroutine bgk_advRel_MSLiquid_generic_WTDF
! ******************************************************************************

end module mus_MSLiquid_module

