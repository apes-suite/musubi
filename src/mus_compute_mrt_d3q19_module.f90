! Copyright (c) 2011-2016, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011, 2017, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!> author: Jiaxing Qi
!! This module provides the definition and methods for
!! MRT advection relaxation scheme.
!! The LB equaton using MRT is
!!    f(t+dt,x+dx) = f - M^(-1) * S * ( (M*f) - m^(eq) )
!!
!! The moments m(1:19) = M * f(1:19) are labeled as
!!  m( 1) = rho
!!  m( 2) = e = rho * (ux^2 + uy^2 + uz^2)
!!  m( 3) = epsilon
!!  m( 4) = jx = rho * ux
!!  m( 5) = qx
!!  m( 6) = jy = rho * uy
!!  m( 7) = qy
!!  m( 8) = jz = rho * uz
!!  m( 9) = qz
!!  m(10) = 3 * pxx = rho * (2ux^2 - uy^2 - uz^2)
!!  m(11) = 3 * Pixx
!!  m(12) = pzz  = rho * (uy^2 - uz^2)
!!  m(13) = Piww
!!  m(14) = pxy  = rho * ux * uy
!!  m(15) = pyz  = rho * uy * uz
!!  m(16) = pzx  = rho * uz * ux
!!  m(17) = mx
!!  m(18) = my
!!  m(19) = mz
!!
!! The non-zero equilibirium moments are given by
!!  meq( 1) = rho
!!  meq( 2) = rho0 * ( ux^2 + uy^2 + uz^2 )
!!  meq( 4) = rho0 * ux
!!  meq( 6) = rho0 * uy
!!  meq( 8) = rho0 * uz
!!  meq(10) = rho0 * ( 2*ux^2 - uy^2 - uz^2 )
!!  meq(12) = rho0 * ( uy^2 - uz^2 )
!!  meq(14) = rho0 * ux * uy
!!  meq(15) = rho0 * uy * uz
!!  meq(16) = rho0 * ux * uz
!!
!! Density (rho) and velocity (ux, uy, uz) are conserved during collision.
!!  i.e. m(1) = meq(1) --> mneq(1) = 0
!!       m(4) = meq(4) --> mneq(4) = 0
!!       m(6) = meq(6) --> mneq(6) = 0
!!       m(8) = meq(8) --> mneq(8) = 0
!!
!! The collision parameters S correspondes to the omega in BGK model.
!!
!! The MRT implementation here is taken from:\n
!! J. Toelke, S. Freudiger, and M. Krafczyk,
!! "An adaptive scheme using hierarchical grids for lattice Boltzmann
!! multi-phase flow simulations," Comput. Fluids, vol. 35, pp. 820–830,
!! 2006. \n
!! Notice that the collision matrix S used in this papar corresponds to
!! -omega in BGK model, because it express the LB equation is slightly
!! different way.
!! In this paper, the following notions are used:\n
!!  s(a) = s(2)
!!  s(b) = s(3)
!!  s(c) = s(5) = s(7) = s(9)
!!  s(d) = s(11) = s(13
!!  s(e) = s(17) = s(18) = s(19)
!!  s(w) = s(10) = s(12) = s(14) = s(15) = s(16)
!! It is suggested that, for D3Q19,
!!  s(a) = s(b) = s(c) = s(d) = s(e) = max( s(w), -1.0 )
!!
!! SubGrid Stress model (SGS)
!! The implementation here is taken from:\n
!! M. Stiebler, M. Krafczyk, S. Freudiger, M. Geier
!! "Lattice Boltzmann large eddy simulation of subcritical flows around a sphere
!! on non-uniform grids", Computers and Mathematics with Applications, vol. 61
!! (2011), pp. 3475-3484
!! Equation 12:\n
!! tau_{total} = 3 * nu0 + dt * 0.5
!!               + 0.5 * ( sqrt( tau0*tau0  + 18 * Cs * Cs * dt * dt * Q ) - tau0 )
!!             = 0.5 * ( tau0 + sqrt( tau0 * tau0 + 18 * Cs * Cs * dt * dt * Q) )
!! Q = sqrt( 2.0 * sum( Pi^{neq} * Pi^{neq} ) )
!!
!! For single field LBM: QQ=nScalars
!!
module mus_mrt_d3q19_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,               only: rk
  use tem_varSys_module,        only: tem_varSys_type, tem_varSys_op_type
  use tem_param_module,         only: rho0, div1_2, div1_4, div1_8, div1_12,  &
    &                                 div1_16,  div1_24, div1_48, div1_72,    &
    &                                 cs2inv, cs4inv, t2cs2inv, t2cs4inv

  ! include musubi modules
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_param_module,         only: mus_param_type
  use mus_varSys_module,        only: mus_varSys_data_type
  use mus_derVarPos_module,     only: mus_derVarPos_type
  use mus_mrtInit_module,       only: MMtrD3Q19, MMivD3Q19

  implicit none

  private

  public :: mus_advRel_kFluid_rMRT_vStd_lD3Q19
  public :: mus_advRel_kFluid_rMRT_vStdNoOpt_lD3Q19

  public :: mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q19
  public :: mus_advRel_kFluidIncomp_rMRT_vStdNoOpt_lD3Q19

  public :: mus_advRel_kCFD_rMRT_vStdNoOpt_l

  !=============================================================================
  ! D3Q19 flow model
  !=============================================================================
  !> Definition of the discrete velocity set

  integer,parameter :: QQ   = 19   !< number of pdf directions

  integer,parameter :: qN00 = 1     !< west             x-
  integer,parameter :: q0N0 = 2     !< south            y-
  integer,parameter :: q00N = 3     !< bottom           z-
  integer,parameter :: q100 = 4     !< east             x+
  integer,parameter :: q010 = 5     !< north            y+
  integer,parameter :: q001 = 6     !< top              z+
  integer,parameter :: q0NN = 7     !<                  z-,y-
  integer,parameter :: q0N1 = 8     !<                  z+,y-
  integer,parameter :: q01N = 9     !<                  z-,y+
  integer,parameter :: q011 = 10    !<                  z+,y+
  integer,parameter :: qN0N = 11    !<                  x-,z-
  integer,parameter :: q10N = 12    !<                  x+,z-
  integer,parameter :: qN01 = 13    !<                  x-,z+
  integer,parameter :: q101 = 14    !<                  x+,z+
  integer,parameter :: qNN0 = 15    !<                  y-,x-
  integer,parameter :: qN10 = 16    !<                  y+,x-
  integer,parameter :: q1N0 = 17    !<                  y-,x+
  integer,parameter :: q110 = 18    !<                  y+,x+
  integer,parameter :: q000 = 19    !< rest density is last

contains

! ****************************************************************************** !
  !> Advection relaxation routine for the MRT model.
  !! This routine has roughly 260 FLOPS per elements.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rMRT_vStd_lD3Q19( fieldProp, inState, outState, &
    &                                            auxField, neigh, nElems,      &
    &                                            nSolve, level, layout, params,&
    &                                            varSys, derVarPos )
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
    integer       :: iElem
    real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
      &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
      &              f000
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: u_x     ! local x-velocity
    real(kind=rk) :: u_y     ! local y-velocity
    real(kind=rk) :: u_z     ! local z-velocity
    ! MRT Variables
    real(kind=rk) :: meq2, meq10, meq12
    real(kind=rk) :: m2, m6, m8, m14, m15, m16
    ! m6, m8 are temporary var
    ! mout1 are temp var
    real(kind=rk) :: mout1, mout2, mout3, mout5, mout7, mout9, mout10, mout11, &
      &              mout12, mout13, mout14, mout15, mout16, mout17, mout18,   &
      &              mout19
    real(kind=rk) :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9,   &
      &              sum10, sum11, sum12, sum13, sum14, sum15
    real(kind=rk) :: c0, c1, c2, c3, c4, c5, c6
    real(kind=rk) :: mout5_4, mout7_4, mout9_4
    real(kind=rk) :: sum_c1_c2, sub_c1_c2
    real(kind=rk) :: sum_5_17, sub_7_18, d1,d2,d3,d4
    real(kind=rk) :: sum_9_19, sub_5_17, e1,e2,e3,e4
    real(kind=rk) :: sum_7_18, sub_9_19, g1,g2,g3,g4
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: omegaKine, omegaBulk, s_mrt(QQ)
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)
    ! MRT omegas
    ! overwrite omegaKine term in the element loop
    s_mrt = fieldProp(1)%fluid%mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    s_mrt(2)  = s_mrt( 2) * div1_24
    s_mrt(3)  = s_mrt( 3) * div1_72
    s_mrt(5)  = s_mrt( 5) * div1_24
    s_mrt(7)  = s_mrt( 7) * div1_24
    s_mrt(9)  = s_mrt( 9) * div1_24
    s_mrt(17) = s_mrt(17) * div1_8
    s_mrt(18) = s_mrt(18) * div1_8
    s_mrt(19) = s_mrt(19) * div1_8

!$omp do schedule(static)
    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1, nSolve

      !> First load all local values into temp array
      fN00 = inState(neigh (( qn00-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f0N0 = inState(neigh (( q0n0-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f00N = inState(neigh (( q00n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f100 = inState(neigh (( q100-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f010 = inState(neigh (( q010-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f001 = inState(neigh (( q001-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f0NN = inState(neigh (( q0nn-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f0N1 = inState(neigh (( q0n1-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f01N = inState(neigh (( q01n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f011 = inState(neigh (( q011-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      fN0N = inState(neigh (( qn0n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f10N = inState(neigh (( q10n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      fN01 = inState(neigh (( qn01-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f101 = inState(neigh (( q101-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      fNN0 = inState(neigh (( qnn0-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      fN10 = inState(neigh (( qn10-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f1N0 = inState(neigh (( q1n0-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f110 = inState(neigh (( q110-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f000 = inState(neigh (( q000-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)

      m6     = f101 + fN0N + f10N + fN01
      m8     = f011 + f0NN + f01N + f0N1
      sum1   = f110 + fNN0 + f1N0 + fN10
      m2     =-f000 + sum1 + m6 + m8

      sum2 = f010 + f0N0
      sum3 = f001 + f00N
      sum4 = 2._rk * ( f100+ fN00 )
      sum5 = sum2 + sum3

      ! epsilon
      mout3 = ( 2._rk*(f000 - sum5) - sum4 + m2 ) * s_mrt(3)

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)

      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      omegaKine = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      s_mrt(10) = omegaKine
      s_mrt(12) = omegaKine
      s_mrt(14) = div1_4 * omegaKine
      s_mrt(15) = div1_4 * omegaKine
      s_mrt(16) = div1_4 * omegaKine

      ! Equilibrium moments
      ! non zero meq are only: 2, 10, 12, 14, 15, 16
      meq2  = rho * ( u_x*u_x + u_y*u_y + u_z*u_z )
      meq10 = rho * 3.0_rk*u_x*u_x - meq2
      meq12 = rho * ( u_y*u_y - u_z*u_z )

      mout2  = s_mrt(2) * (m2 - meq2)

      !pxy
      m14 = f110 + fNN0 - f1N0 - fN10
      mout14 = s_mrt(14) * (m14 - rho*u_x*u_y)
      !pyz
      m15 = f011 + f0NN - f01N - f0N1
      mout15 = s_mrt(15) * (m15 - rho*u_y*u_z)
      !pxz
      m16 = f101 + fN0N - f10N - fN01
      mout16 = s_mrt(16) * (m16 - rho*u_x*u_z)

      sum6 = sum1 + m6 - m8 * 2.0_rk
      sum7 = sum4 - sum5
      !3pxx
      mout10 = ( sum7 + sum6 - meq10 ) * s_mrt(10)
      !3pixx
      mout11 = ( - sum7 + sum6 ) * s_mrt(11)

      sum8 = sum1 - m6
      sum9 = sum2 - sum3
      !pww
      mout12 = ( sum8 + sum9 - meq12 ) * s_mrt(12)
      !piww
      mout13 = ( sum8 - sum9 ) * s_mrt(13)

c1 = f110 - fNN0
c2 = f1N0 - fN10
c3 = f101 - fN0N
c4 = f10N - fN01

      sum10 = c1 + c2
      sum11 = c3 + c4
      !qx
      mout5  = ( sum10 + sum11 - 2.0_rk*(f100 - fN00) ) * s_mrt(5)
      !mx
      mout17 = ( sum10 - sum11 ) * s_mrt(17)

c5 = f011 - f0NN
c6 = f01N - f0N1

      sum12 = c1 - c2
      sum13 = c5 + c6
      !qy
      mout7  = ( sum12 + sum13 - 2.0_rk*(f010 - f0N0) ) * s_mrt(7)
      !my
      mout18 = ( - sum12 + sum13 ) * s_mrt(18)

      sum14 = c3 - c4
      sum15 = c5 - c6
      !qz
      mout9  = ( sum14 + sum15 - 2.0_rk*(f001 - f00N) ) * s_mrt(9)
      !mz
      mout19 = ( sum14 - sum15 ) * s_mrt(19)

      ! Transformation back to PDF
      outstate( (ielem-1)*qq+ 19+(1-1)*qq) = f000 + 12._rk*(mout2-mout3)

! -------------------------------------------------------------------------------
      c0 = - 4._rk*mout3 + div1_12*(mout10 - mout11)
      mout5_4 = mout5 * 4.0_rk
      outstate( (ielem-1)*qq+  4+(1-1)*qq) = f100 - ( c0 - mout5_4 )
      outstate( (ielem-1)*qq+  1+(1-1)*qq) = fN00 - ( c0 + mout5_4 )
! -------------------------------------------------------------------------------

! -------------------------------------------------------------------------------
      c1 = - 4._rk*mout3 - div1_24*(mout10 - mout11)
      c2 = div1_8 *(mout12 - mout13)

      sum_c1_c2 = c1 + c2
      mout7_4 = mout7 * 4.0_rk
      outstate( (ielem-1)*qq+  5+(1-1)*qq) = f010 - (sum_c1_c2 - mout7_4)
      outstate( (ielem-1)*qq+  2+(1-1)*qq) = f0N0 - (sum_c1_c2 + mout7_4)

      sub_c1_c2 = c1 - c2
      mout9_4 = mout9 * 4.0_rk
      outstate( (ielem-1)*qq+  6+(1-1)*qq) = f001 - (sub_c1_c2 - mout9_4)
      outstate( (ielem-1)*qq+  3+(1-1)*qq) = f00N - (sub_c1_c2 + mout9_4)
! -------------------------------------------------------------------------------

! -------------------------------------------------------------------------------
      mout1 = mout2 + mout3
      c3 =   mout1  + div1_48 * (mout10 + mout11)   &
       &            + div1_16 * (mout12 + mout13)
      sum_5_17 = mout5 + mout17
      sub_7_18 = mout7 - mout18

      d1 = c3 + mout14
      d2 = sum_5_17 + sub_7_18
      outstate( (ielem-1)*qq+ 18+(1-1)*qq) = f110-(d1+d2)
      outstate( (ielem-1)*qq+ 15+(1-1)*qq) = fNN0-(d1-d2)

      d3 = c3 - mout14
      d4 = sum_5_17 - sub_7_18
      outstate( (ielem-1)*qq+ 17+(1-1)*qq) = f1N0-(d3+d4)
      outstate( (ielem-1)*qq+ 16+(1-1)*qq) = fN10-(d3-d4)
! -------------------------------------------------------------------------------

! -------------------------------------------------------------------------------
      c4 = c3 - div1_8*(mout12+mout13)

      sum_9_19 = mout9 + mout19
      sub_5_17 = mout5 - mout17

      e1 = c4 + mout16
      e2 = sum_9_19 + sub_5_17
      outstate( (ielem-1)*qq+ 14+(1-1)*qq) = f101 - ( e1 + e2 )
      outstate( (ielem-1)*qq+ 11+(1-1)*qq) = fN0N - ( e1 - e2 )

      e3 = c4 - mout16
      e4 = sum_9_19 - sub_5_17
      outstate( (ielem-1)*qq+ 12+(1-1)*qq) = f10N - ( e3 - e4 )
      outstate( (ielem-1)*qq+ 13+(1-1)*qq) = fN01 - ( e3 + e4 )
! -------------------------------------------------------------------------------

! -------------------------------------------------------------------------------
      c5 = mout1 - div1_24*(mout10+mout11)
      sum_7_18 = mout7 + mout18
      sub_9_19 = mout9 - mout19

      g1 = c5 + mout15
      g2 = sum_7_18 + sub_9_19
      outstate( (ielem-1)*qq+ 10+(1-1)*qq) = f011 - ( g1 + g2 )
      outstate( (ielem-1)*qq+  7+(1-1)*qq) = f0NN - ( g1 - g2 )

      g3 = c5 - mout15
      g4 = sum_7_18 - sub_9_19
      outstate( (ielem-1)*qq+  9+(1-1)*qq) = f01N - ( g3 + g4 )
      outstate( (ielem-1)*qq+  8+(1-1)*qq) = f0N1 - ( g3 - g4 )
! -------------------------------------------------------------------------------

    enddo nodeloop
!$omp end do

  end subroutine mus_advRel_kFluid_rMRT_vStd_lD3Q19
! ****************************************************************************** !


! ****************************************************************************** !
  !> Advection relaxation routine for the MRT model.
  !! This routine has roughly 205 FLOPS per element.
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q19( fieldProp, inState,    &
    &                                                  outState, auxField,    &
    &                                                  neigh, nElems, nSolve, &
    &                                                  level, layout, params, &
    &                                                  varSys, derVarPos )
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
    integer       :: iElem
    real(kind=rk) :: fN00, f0N0, f00N, f100, f010, f001, f0NN, f0N1, f01N, &
      &              f011, fN0N, f10N, fN01, f101, fNN0, fN10, f1N0, f110, &
      &              f000
    real(kind=rk) :: u_x     ! local x-velocity
    real(kind=rk) :: u_y     ! local y-velocity
    real(kind=rk) :: u_z     ! local z-velocity
    ! MRT Variables
    real(kind=rk) :: meq2, meq10, meq12
    real(kind=rk) :: m2, m6, m8, m14, m15, m16
    ! m6, m8 are temporary var
    ! mout1 are temp var
    real(kind=rk) :: mout1, mout2, mout3, mout5, mout7, mout9, mout10, mout11, &
      &              mout12, mout13, mout14, mout15, mout16, mout17, mout18,   &
      &              mout19
    real(kind=rk) :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9,   &
      &              sum10, sum11, sum12, sum13, sum14, sum15
    real(kind=rk) :: c0, c1, c2, c3, c4, c5, c6
    real(kind=rk) :: mout5_4, mout7_4, mout9_4
    real(kind=rk) :: sum_c1_c2, sub_c1_c2
    real(kind=rk) :: sum_5_17, sub_7_18, d1,d2,d3,d4
    real(kind=rk) :: sum_9_19, sub_5_17, e1,e2,e3,e4
    real(kind=rk) :: sum_7_18, sub_9_19, g1,g2,g3,g4
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: omegaKine, omegaBulk, s_mrt(QQ)
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)
    ! MRT omegas
    ! overwrite omegaKine term in the element loop
    s_mrt = fieldProp(1)%fluid%mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    s_mrt(2)  = s_mrt( 2) * div1_24
    s_mrt(3)  = s_mrt( 3) * div1_72
    s_mrt(5)  = s_mrt( 5) * div1_24
    s_mrt(7)  = s_mrt( 7) * div1_24
    s_mrt(9)  = s_mrt( 9) * div1_24
    s_mrt(17) = s_mrt(17) * div1_8
    s_mrt(18) = s_mrt(18) * div1_8
    s_mrt(19) = s_mrt(19) * div1_8

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1, nSolve

      !> First load all local values into temp array
      fN00 = inState(neigh (( qn00-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f0N0 = inState(neigh (( q0n0-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f00N = inState(neigh (( q00n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f100 = inState(neigh (( q100-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f010 = inState(neigh (( q010-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f001 = inState(neigh (( q001-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f0NN = inState(neigh (( q0nn-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f0N1 = inState(neigh (( q0n1-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f01N = inState(neigh (( q01n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f011 = inState(neigh (( q011-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      fN0N = inState(neigh (( qn0n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f10N = inState(neigh (( q10n-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      fN01 = inState(neigh (( qn01-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f101 = inState(neigh (( q101-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      fNN0 = inState(neigh (( qnn0-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      fN10 = inState(neigh (( qn10-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f1N0 = inState(neigh (( q1n0-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f110 = inState(neigh (( q110-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)
      f000 = inState(neigh (( q000-1)* nelems+ ielem)+( 1-1)* qq+ varsys%nscalars*0)

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars

      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      ! Equilibrium moments
      ! non zero meq are only: 2, 10, 12, 14, 15, 16
      meq2  = u_x*u_x + u_y*u_y + u_z*u_z
      meq10 = 3.0_rk*u_x*u_x - meq2
      meq12 = u_y*u_y - u_z*u_z

      omegaKine = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      s_mrt(10) = omegaKine
      s_mrt(12) = omegaKine
      s_mrt(14) = div1_4 * omegaKine
      s_mrt(15) = div1_4 * omegaKine
      s_mrt(16) = div1_4 * omegaKine

      !pxy
      m14 = f110 + fNN0 - f1N0 - fN10
      mout14 = s_mrt(14) * (m14 - u_x*u_y)
      !pyz
      m15 = f011 + f0NN - f01N - f0N1
      mout15 = s_mrt(15) * (m15 - u_y*u_z)
      !pxz
      m16 = f101 + fN0N - f10N - fN01
      mout16 = s_mrt(16) * (m16 - u_x*u_z)

      ! -------------------------------------------------------------------------
      ! now calculate mrt against the bgk
      m6 = f101 + fN0N + f10N + fN01
      m8 = f011 + f0NN + f01N + f0N1

      sum1   =   f110 + fNN0 + f1N0 + fN10
      m2     = - f000 + sum1 + m6 + m8
      mout2  =   s_mrt(2)  * (m2 - meq2)

      sum2 = f010 + f0N0
      sum3 = f001 + f00N
      sum4 = 2._rk * ( f100+ fN00 )
      sum5 = sum2 + sum3

      ! epsilon
      mout3 = ( 2._rk*(f000 - sum5) - sum4 + m2 ) * s_mrt(3)

      sum6 = sum1 + m6 - m8 * 2.0_rk
      sum7 = sum4 - sum5
      !3pxx
      mout10 = ( sum7 + sum6 - meq10 ) * s_mrt(10)
      !3pixx
      mout11 = ( - sum7 + sum6 ) * s_mrt(11)

      sum8 = sum1 - m6
      sum9 = sum2 - sum3
      !pww
      mout12 = ( sum8 + sum9 - meq12 ) * s_mrt(12)
      !piww
      mout13 = ( sum8 - sum9 ) * s_mrt(13)

c1 = f110 - fNN0
c2 = f1N0 - fN10
c3 = f101 - fN0N
c4 = f10N - fN01

      sum10 = c1 + c2
      sum11 = c3 + c4
      !qx
      mout5  = ( sum10 + sum11 - 2.0_rk*(f100 - fN00) ) * s_mrt(5)
      !mx
      mout17 = ( sum10 - sum11 ) * s_mrt(17)

c5 = f011 - f0NN
c6 = f01N - f0N1

      sum12 = c1 - c2
      sum13 = c5 + c6
      !qy
      mout7  = ( sum12 + sum13 - 2.0_rk*(f010 - f0N0) ) * s_mrt(7)
      !my
      mout18 = ( - sum12 + sum13 ) * s_mrt(18)

      sum14 = c3 - c4
      sum15 = c5 - c6
      !qz
      mout9  = ( sum14 + sum15 - 2.0_rk*(f001 - f00N) ) * s_mrt(9)
      !mz
      mout19 = ( sum14 - sum15 ) * s_mrt(19)

      ! Transformation back to PDF
      outstate( (ielem-1)*qq+ 19+(1-1)*qq) = f000 + 12._rk*(mout2-mout3)

! -------------------------------------------------------------------------------
      c0 = - 4._rk*mout3 + div1_12*(mout10 - mout11)
      mout5_4 = mout5 * 4.0_rk
      outstate( (ielem-1)*qq+  4+(1-1)*qq) = f100 - ( c0 - mout5_4 )
      outstate( (ielem-1)*qq+  1+(1-1)*qq) = fN00 - ( c0 + mout5_4 )
! -------------------------------------------------------------------------------

! -------------------------------------------------------------------------------
      c1 = - 4._rk*mout3 - div1_24*(mout10 - mout11)
      c2 = div1_8 *(mout12 - mout13)

      sum_c1_c2 = c1 + c2
      mout7_4 = mout7 * 4.0_rk
      outstate( (ielem-1)*qq+  5+(1-1)*qq) = f010 - (sum_c1_c2 - mout7_4)
      outstate( (ielem-1)*qq+  2+(1-1)*qq) = f0N0 - (sum_c1_c2 + mout7_4)

      sub_c1_c2 = c1 - c2
      mout9_4 = mout9 * 4.0_rk
      outstate( (ielem-1)*qq+  6+(1-1)*qq) = f001 - (sub_c1_c2 - mout9_4)
      outstate( (ielem-1)*qq+  3+(1-1)*qq) = f00N - (sub_c1_c2 + mout9_4)
! -------------------------------------------------------------------------------

! -------------------------------------------------------------------------------
      mout1 = mout2 + mout3
      c3 =   mout1  + div1_48 * (mout10 + mout11)   &
       &            + div1_16 * (mout12 + mout13)
      sum_5_17 = mout5 + mout17
      sub_7_18 = mout7 - mout18

      d1 = c3 + mout14
      d2 = sum_5_17 + sub_7_18
      outstate( (ielem-1)*qq+ 18+(1-1)*qq) = f110-(d1+d2)
      outstate( (ielem-1)*qq+ 15+(1-1)*qq) = fNN0-(d1-d2)

      d3 = c3 - mout14
      d4 = sum_5_17 - sub_7_18
      outstate( (ielem-1)*qq+ 17+(1-1)*qq) = f1N0-(d3+d4)
      outstate( (ielem-1)*qq+ 16+(1-1)*qq) = fN10-(d3-d4)
! -------------------------------------------------------------------------------

! -------------------------------------------------------------------------------
      c4 = c3 - div1_8*(mout12+mout13)

      sum_9_19 = mout9 + mout19
      sub_5_17 = mout5 - mout17

      e1 = c4 + mout16
      e2 = sum_9_19 + sub_5_17
      outstate( (ielem-1)*qq+ 14+(1-1)*qq) = f101 - ( e1 + e2 )
      outstate( (ielem-1)*qq+ 11+(1-1)*qq) = fN0N - ( e1 - e2 )

      e3 = c4 - mout16
      e4 = sum_9_19 - sub_5_17
      outstate( (ielem-1)*qq+ 12+(1-1)*qq) = f10N - ( e3 - e4 )
      outstate( (ielem-1)*qq+ 13+(1-1)*qq) = fN01 - ( e3 + e4 )
! -------------------------------------------------------------------------------

! -------------------------------------------------------------------------------
      c5 = mout1 - div1_24*(mout10+mout11)
      sum_7_18 = mout7 + mout18
      sub_9_19 = mout9 - mout19

      g1 = c5 + mout15
      g2 = sum_7_18 + sub_9_19
      outstate( (ielem-1)*qq+ 10+(1-1)*qq) = f011 - ( g1 + g2 )
      outstate( (ielem-1)*qq+  7+(1-1)*qq) = f0NN - ( g1 - g2 )

      g3 = c5 - mout15
      g4 = sum_7_18 - sub_9_19
      outstate( (ielem-1)*qq+  9+(1-1)*qq) = f01N - ( g3 + g4 )
      outstate( (ielem-1)*qq+  8+(1-1)*qq) = f0N1 - ( g3 - g4 )
! -------------------------------------------------------------------------------

    enddo nodeloop

  end subroutine mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q19
! ****************************************************************************** !


! ****************************************************************************** !
  !> No comment yet!
  !!
  !! TODO add comment
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluidIncomp_rMRT_vStdNoOpt_lD3Q19( fieldProp, inState,&
    &                                                       outState, auxField,&
    &                                                       neigh, nElems,     &
    &                                                       nSolve, level,     &
    &                                                       layout, params,    &
    &                                                       varSys, derVarPos  )
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
    integer       :: iElem, iDir
    real(kind=rk) :: pdfTmp( QQ ) ! temporary local pdf values
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: u_x     ! local x-velocity
    real(kind=rk) :: u_y     ! local y-velocity
    real(kind=rk) :: u_z     ! local z-velocity
    ! MRT Variables
    real(kind=rk) :: s_mrt( QQ ), meq( QQ )
    real(kind=rk) :: mneq( QQ ), mom( QQ )
    real(kind=rk) :: fneq( QQ )
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: omegaKine, omegaBulk
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)
    ! MRT omegas
    ! overwrite omegaKine term in the element loop
    s_mrt = fieldProp(1)%fluid%mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, &
      &                               QQ=QQ)

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1, nSolve

      !> First load all local values into temp array
      do iDir = 1, QQ
        pdfTmp( iDir ) = inState( neigh (( idir-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      end do

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)

      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      ! -------------------------------------------------------------------------
      ! Equilibrium moments
      meq(1:QQ) =  0.0_rk

      meq( 1) = rho
      meq( 2) = rho0 * ( u_x * u_x + u_y * u_y + u_z * u_z )
      meq( 4) = u_x*rho0
      meq( 6) = u_y*rho0
      meq( 8) = u_z*rho0
      meq(10) = rho0*(2.0_rk*u_x*u_x - ( u_y*u_y + u_z*u_z ))
      meq(12) = rho0*( u_y*u_y - u_z*u_z )
      meq(14) = rho0*u_x*u_y
      meq(15) = rho0*u_y*u_z
      meq(16) = rho0*u_x*u_z

      ! convert pdf into moment
      do iDir = 1, QQ
        mom(iDir) = sum( pdfTmp(:) * MMtrD3Q19(iDir,:) )
      end do

      ! update kinematic omega part in relaxation matrix
      omegaKine = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      s_mrt(10) = omegaKine
      s_mrt(12) = omegaKine
      s_mrt(14:16) = omegaKine

      ! compute neq moment
      mneq(:) = s_mrt(:) * ( mom(:) - meq(:) )

      ! compute fNeq
      do iDir = 1, QQ
        fneq(iDir) = sum( MMIvD3Q19(iDir,:) * mneq(:) )
        outState((ielem-1)*qq+idir+(1-1)*qq) = pdfTmp(iDir) - fneq(iDir)
      end do

    enddo nodeloop

  end subroutine mus_advRel_kFluidIncomp_rMRT_vStdNoOpt_lD3Q19
! ****************************************************************************** !

! ****************************************************************************** !
  !> No comment yet!
  !!
  !! TODO add comment
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kFluid_rMRT_vStdNoOpt_lD3Q19( fieldProp, inState,      &
    &                                                 outState, auxField,      &
    &                                                 neigh, nElems, nSolve,   &
    &                                                 level, layout, params,   &
    &                                                 varSys, derVarPos )
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
    integer       :: iElem, iDir
    real(kind=rk) :: pdfTmp( QQ ) ! temporary local pdf values
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: u_x     ! local x-velocity
    real(kind=rk) :: u_y     ! local y-velocity
    real(kind=rk) :: u_z     ! local z-velocity
    ! MRT Variables
    real(kind=rk) :: s_mrt( QQ ), meq( QQ )
    real(kind=rk) :: mneq( QQ ), mom( QQ )
    real(kind=rk) :: fneq( QQ )
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: omegaKine, omegaBulk
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)
    ! MRT omegas
    ! overwrite omegaKine term in the element loop
    s_mrt = fieldProp(1)%fluid%mrtPtr(omegaKine=1.0_rk, omegaBulk=omegaBulk, QQ=QQ)

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1,nSolve

      !> First load all local values into temp array
      do iDir = 1, QQ
        pdfTmp( iDir ) = inState( neigh (( idir-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      end do

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)

      ! local x-, y- and z-velocity
      u_x = auxField(elemOff + vel_pos(1))
      u_y = auxField(elemOff + vel_pos(2))
      u_z = auxField(elemOff + vel_pos(3))

      ! -------------------------------------------------------------------------
      ! Equilibrium moments
      meq(1:QQ) =  0.0_rk

      meq( 1) = rho
      meq( 2) = rho * ( u_x * u_x + u_y * u_y + u_z * u_z )
      meq( 4) = u_x*rho
      meq( 6) = u_y*rho
      meq( 8) = u_z*rho
      meq(10) = rho*(2.0_rk*u_x*u_x - ( u_y*u_y + u_z*u_z ))
      meq(12) = rho*( u_y*u_y - u_z*u_z )
      meq(14) = rho*u_x*u_y
      meq(15) = rho*u_y*u_z
      meq(16) = rho*u_x*u_z

      ! convert pdf into moment
      do iDir = 1, QQ
        mom(iDir) = sum( pdfTmp(:) * MMtrD3Q19(iDir,:) )
      end do

      ! update kinematic omega part in relaxation matrix
      omegaKine = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      s_mrt(10) = omegaKine
      s_mrt(12) = omegaKine
      s_mrt(14:16) = omegaKine

      ! compute neq moment
      mneq(:) = s_mrt(:) * ( mom(:) - meq(:) )

      ! compute fNeq
      do iDir = 1, QQ
        fneq(iDir) = sum( MMIvD3Q19(iDir,:) * mneq(:) )
        outState((ielem-1)*qq+idir+(1-1)*qq) = pdfTmp(iDir) - fneq(iDir)
      end do

    enddo nodeloop

  end subroutine mus_advRel_kFluid_rMRT_vStdNoOpt_lD3Q19
! ****************************************************************************** !


! **************************************************************************** !
  !> Unoptimized explicit implementation
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_advRel_kCFD_rMRT_vStdNoOpt_l( fieldProp, inState, outState, &
    &                                          auxField, neigh, nElems,      &
    &                                          nSolve, level, layout, params,&
    &                                          varSys, derVarPos )
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
    integer       :: iElem, iDir
    real(kind=rk) :: pdfTmp( layout%fStencil%QQ ) ! temporary local pdf values
    real(kind=rk) :: rho     ! local density
    real(kind=rk) :: vel(3)  ! local velocity
    real(kind=rk) :: fEq( layout%fStencil%QQ ) !< equilibrium distribution
    ! MRT Variables
    real(kind=rk) :: fneq( layout%fStencil%QQ )
    real(kind=rk) :: mneq( layout%fStencil%QQ )
    real(kind=rk) :: s_mrt( layout%fStencil%QQ )
    real(kind=rk) :: mInvXomega(layout%fStencil%QQ, layout%fStencil%QQ)
    integer :: QQ
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: omegaKine, omegaBulk, colli_term
    ! ---------------------------------------------------------------------------
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    QQ = layout%fStencil%QQ
    omegaBulk = fieldProp(1)%fluid%omegaBulkLvl(level)

    !NEC$ ivdep
!cdir nodep
!ibm* novector
!dir$ novector
    nodeloop: do iElem = 1,nSolve

      !> First load all local values into temp array
      do iDir = 1, QQ
        pdfTmp( iDir ) = inState( neigh (( idir-1)* nelems+ ielem)+( 1-1)* qq+ qq*0)
      end do

      ! element offset for auxField array
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! local density
      rho = auxField(elemOff + dens_pos)

      ! local x-, y- and z-velocity
      vel(1) = auxField(elemOff + vel_pos(1))
      vel(2) = auxField(elemOff + vel_pos(2))
      vel(3) = auxField(elemOff + vel_pos(3))

      ! -------------------------------------------------------------------------
      ! Calculate the equilibrium distribution function
      fEq = layout%quantities%pdfEq_ptr(rho=rho, vel=vel, QQ=QQ)

      ! compute fNeq
      fneq = pdfTmp - fEq

      ! convert fneq to moments neq
      mneq = matmul( layout%moment%toMoments%A(1:QQ, 1:QQ), fneq(1:QQ) )

      ! MRT omegas
      ! overwrite omegaKine term in the element loop
      omegaKine = fieldProp(1)%fluid%viscKine%omLvl(level)%val(iElem)
      s_mrt = fieldProp(1)%fluid%mrtPtr(omegaKine=omegaKine, &
        &                               omegaBulk=omegaBulk, &
        &                               QQ=QQ                )

      ! compute Minv times relaxation matrix
      do iDir = 1, QQ
        mInvXOmega(1:QQ,iDir) = layout%moment%toPDF%A(1:QQ, iDir) * s_mrt(iDir)
      end do

      do iDir = 1, QQ
        ! collision term: M^-1 S M (f-feq)
        colli_term = dot_product( mInvXOmega(iDir, 1:QQ), mneq(1:QQ) )

        outState((ielem-1)*qq+idir+(1-1)*qq) &
          & = pdfTmp(iDir) - colli_term
      end do

    enddo nodeloop

  end subroutine mus_advRel_kCFD_rMRT_vStdNoOpt_l
! **************************************************************************** !

end module mus_mrt_d3q19_module
! ****************************************************************************** !
