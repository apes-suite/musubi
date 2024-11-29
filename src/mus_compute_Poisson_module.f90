! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2017-2018, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! Make sure loglvl is defined to an integer value.
! Usually this should be defined on the command line with -Dloglvl=
! to some value.
! **************************************************************************** !
!> This module provides the definition and methods to solve
!! poisson equation.
!! author: Kannan Masilamani
!! Implementation is based on:
!! Luo, K., Wu, J., Yi, H., & Tan, H. (2016).
!! International Journal of Heat and Mass Transfer Lattice Boltzmann modelling
!! of electro-thermo-convection in a planar layer of dielectric liquid
!! subjected to unipolar injection and thermal gradient.
!! International Journal of Heat and Mass Transfer, 103, 832–846.
!!
!! Guo, Z. (2014). A Coupled Lattice Boltzmann Method to Solve Nernst
!! – Planck Model for Simulating Electro-osmotic Flows. Journal of
!! Scientific Computing, (61), 222–238. http://doi.org/10.1007/s10915-014-9820-6
!!
!! This implementation solve equation of the form
!! $$ \partial_t \phi = \nabla^2 \phi + G $$
!! $$G$$ is the source term
!! For Poisson Eq: $$G = \frac{\rho_e}{\epsilon_r \epsilon_0}$$
!! For Poisson Boltzmann Linear Eq: $$G = -k^2 \phi$$
!! For Poisson Boltzmann nonLinear Eq:
!! $$ G = \frac{1}{\epsilon_r \epsilon_0}
!!    \sum_i { e z_i c_{\infty} N_A exp(-\frac{z e }{k_b T}\phi) } $$
!! For 1:1 Electrolye solution, above equation is simplied to
!! $$G = -\frac{-2 c_{\infty} N_A z e}{\epsilon_r \epsilon_0}
!!       sinh(\frac{z e }{k_b T}\phi)  $$
module mus_compute_Poisson_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,               only: rk
  use tem_varSys_module,        only: tem_varSys_type
  use tem_geometry_module,      only: tem_baryOfId
  use tem_aux_module,           only: tem_abort
  use tem_param_module,        only: div1_36, div4_9, div1_9

  ! include musubi modules
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_param_module,         only: mus_param_type
  use mus_derVarPos_module,     only: mus_derVarPos_type

  implicit none

  private

  public :: mus_Poisson_advRel_d2q9
  public :: mus_Poisson_advRel_generic
  public :: mus_PBLinear_advRel_generic
  public :: mus_PBnonLinear_advRel_generic

  ! General direction index for 3D
  integer,parameter :: QQ_27 = 27  !< number of pdf directions
  integer,parameter :: QQ_19 = 19  !< number of pdf directions
  integer,parameter :: qN00 = 1   !< west             x-
  integer,parameter :: q0N0 = 2   !< south            y-
  integer,parameter :: q00N = 3   !< bottom           z-
  integer,parameter :: q100 = 4   !< east             x+
  integer,parameter :: q010 = 5   !< north            y+
  integer,parameter :: q001 = 6   !< top              z+
  integer,parameter :: q0NN = 7   !<                  z-,y-
  integer,parameter :: q0N1 = 8   !<                  z+,y-
  integer,parameter :: q01N = 9   !<                  z-,y+
  integer,parameter :: q011 = 10  !<                  z+,y+
  integer,parameter :: qN0N = 11  !<                  x-,z-
  integer,parameter :: q10N = 12  !<                  x+,z-
  integer,parameter :: qN01 = 13  !<                  x-,z+
  integer,parameter :: q101 = 14  !<                  x+,z+
  integer,parameter :: qNN0 = 15  !<                  y-,x-
  integer,parameter :: qN10 = 16  !<                  y+,x-
  integer,parameter :: q1N0 = 17  !<                  y-,x+
  integer,parameter :: q110 = 18  !<                  y+,x+
  integer,parameter :: q000_19 = 19 !< rest density for d3q19

  integer,parameter :: qNNN = 19  !<                  z-,y-,x-
  integer,parameter :: qNN1 = 20  !<                  z+,y-,x-
  integer,parameter :: qN1N = 21  !<                  z-,y+,x-
  integer,parameter :: qN11 = 22  !<                  z+,y+,x-
  integer,parameter :: q1NN = 23  !<                  z-,y-,x+
  integer,parameter :: q1N1 = 24  !<                  z+,y-,x+
  integer,parameter :: q11N = 25  !<                  z-,y+,x+
  integer,parameter :: q111 = 26  !<                  z+,y+,x+
  integer,parameter :: q000_27 = 27 !< rest density for d3q27

  ! General direction index for 2D
  integer, parameter :: QQ_9 = 9   !< number of pdf directions
  integer, parameter :: qN0 = 1  !< west             x-
  integer, parameter :: q0N = 2  !< south            y-
  integer, parameter :: q10 = 3  !< east             x+
  integer, parameter :: q01 = 4  !< north            y+
  integer, parameter :: qNN = 5  !<                  y-,x-
  integer, parameter :: qN1 = 6  !<                  y+,x-
  integer, parameter :: q1N = 7  !<                  y-,x+
  integer, parameter :: q11 = 8  !<                  y+,x+
  integer, parameter :: q00_9 = 9  !< rest density of d2q9

contains

  ! ************************************************************************** !
  !> Advection relaxation routine for the
  !! poisson equation with an explicit calculation of all equilibrium
  !! quantities. Slow and simple.
  !! $$ \nabla^2 \phi = - \frac{\rho_e}{\epsilon_r \epsilon_0} $$
  !! The right hand side of equation is added as a source term in
  !! mus_apply_sourceTerms routine
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_Poisson_advRel_d2q9( fieldProp, inState, outState, auxField, &
    &                                 neigh, nElems, nSolve, level, layout,   &
    &                                 params, varSys, derVarPos               )
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
    integer :: iElem  ! voxel element counter
    ! temporary distribution variables
    real(kind=rk) :: pdfTmp(QQ_9)
    real(kind=rk) :: omega, omega_fac, om_pot, fac14, fac58, fac9
    ! ---------------------------------------------------------------------------
    !>\todo KM 20170821: change omega for multilevel
    omega = fieldProp(1)%poisson%omega
    omega_fac = 1.0_rk - omega

!$omp do schedule(static)
    !NEC$ ivdep
    !DIR$ NOVECTOR
    nodeloop: do iElem = 1, nSolve
      !> Generic fetching step:
      !! Streaming for pull
      !! Local copy for push
      pdfTmp(qN0  ) = inState( neigh((qn0  -1)* nelems+ ielem)+( 1-1)* qq_9+ qq_9*0)
      pdfTmp(q0N  ) = inState( neigh((q0n  -1)* nelems+ ielem)+( 1-1)* qq_9+ qq_9*0)
      pdfTmp(q10  ) = inState( neigh((q10  -1)* nelems+ ielem)+( 1-1)* qq_9+ qq_9*0)
      pdfTmp(q01  ) = inState( neigh((q01  -1)* nelems+ ielem)+( 1-1)* qq_9+ qq_9*0)
      pdfTmp(qNN  ) = inState( neigh((qnn  -1)* nelems+ ielem)+( 1-1)* qq_9+ qq_9*0)
      pdfTmp(qN1  ) = inState( neigh((qn1  -1)* nelems+ ielem)+( 1-1)* qq_9+ qq_9*0)
      pdfTmp(q1N  ) = inState( neigh((q1n  -1)* nelems+ ielem)+( 1-1)* qq_9+ qq_9*0)
      pdfTmp(q11  ) = inState( neigh((q11  -1)* nelems+ ielem)+( 1-1)* qq_9+ qq_9*0)
      pdfTmp(q00_9) = inState( neigh((q00_9-1)* nelems+ ielem)+( 1-1)* qq_9+ qq_9*0)

      ! omega * potential
      om_pot = omega * auxField( iElem )

      fac14 = om_pot * div1_9
      outState( ( ielem-1)* qq_9 + qn0+( 1-1)* qq_9) &
        & =  omega_fac * pdfTmp(qN0) + fac14
      outState( ( ielem-1)* qq_9 + q0n+( 1-1)* qq_9) &
        & =  omega_fac * pdfTmp(q0N) + fac14
      outState( ( ielem-1)* qq_9 + q10+( 1-1)* qq_9) &
        & =  omega_fac * pdfTmp(q10) + fac14
      outState( ( ielem-1)* qq_9 + q01+( 1-1)* qq_9) &
        & =  omega_fac * pdfTmp(q01) + fac14

      fac58 = om_pot * div1_36
      outState( ( ielem-1)* qq_9 + qnn+( 1-1)* qq_9) &
        & =  omega_fac * pdfTmp(qNN) + fac58
      outState( ( ielem-1)* qq_9 + qn1+( 1-1)* qq_9) &
        & =  omega_fac * pdfTmp(qN1) + fac58
      outState( ( ielem-1)* qq_9 + q1n+( 1-1)* qq_9) &
        & =  omega_fac * pdfTmp(q1N) + fac58
      outState( ( ielem-1)* qq_9 + q11+( 1-1)* qq_9) &
        & =  omega_fac * pdfTmp(q11) + fac58

      fac9 = om_pot * div4_9
      outState( ( ielem-1)* qq_9 + q00_9+( 1-1)* qq_9) &
        & =  omega_fac * pdfTmp(q00_9) + fac9

    end do nodeloop
!$omp end do nowait

  end subroutine mus_Poisson_advRel_d2q9
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Advection relaxation routine for the
  !! poisson equation with an explicit calculation of all equilibrium
  !! quantities. Slow and simple.
  !! $$ \nabla^2 \phi = - \frac{\rho_e}{\epsilon_r \epsilon_0} $$
  !! The right hand side of equation is added as a source term in
  !! mus_apply_sourceTerms routine
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_Poisson_advRel_generic( fieldProp, inState, outState,    &
    &                                    auxField, neigh, nElems, nSolve, &
    &                                    level, layout, params, varSys,   &
    &                                    derVarPos                        )
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
    integer :: iElem, iDir                       ! voxel element counter
    integer :: QQ, nScalars
    ! temporary distribution variables
    real(kind=rk) :: pdfTmp
    real(kind=rk) :: potential
    real(kind=rk) :: eqState
    real(kind=rk) :: omega
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nScalars = varSys%nScalars

    !>\todo KM 20170821: change omega for multilevel
    omega = fieldProp(1)%poisson%omega

    nodeloop: do iElem = 1, nSolve
      ! potential
      potential = auxField( iElem )

      do iDir = 1, QQ
        !> Calculate equilibrium distribution functions fEq
        eqstate = layout%weight(iDir)*potential

        !> Generic fetching step:
        !! Streaming for pull
        !! Local copy for push
        pdfTmp = inState( neigh((idir-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

        !> Relaxation
        outState( (ielem-1)*nscalars+ idir+(1-1)*qq) &
          & = pdfTmp - omega * ( pdfTmp - eqState )
      end do ! iDir

    end do nodeloop

  end subroutine mus_Poisson_advRel_generic
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Advection relaxation routine for the linear
  !! poisson boltzmann equation with an explicit calculation of all equilibrium
  !! quantities. Slow and simple.
  !! $$ \nabla^2 \phi = k^2 \phi $$
  !! Where k^2 is inverse of debye length and in this kernel refered as
  !! RHS_coeff
  !! $$ k^2 = \sum_i \frac{ c_{\infty}z_i^2 e^2}{\epsilon_r \epsilon_0 k_b T} $$
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_PBLinear_advRel_generic( fieldProp, inState, outState,    &
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
    integer :: iElem, iDir                       ! voxel element counter
    integer :: QQ, nScalars
    ! temporary distribution variables
    real(kind=rk) :: pdfTmp
    real(kind=rk) :: potential
    real(kind=rk) :: rhs_coeff
    real(kind=rk) :: eqState, source
    real(kind=rk) :: omega, pot_diff
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nScalars = varSys%nScalars

    !>\todo KM 20170821: change omega for multilevel
    omega = fieldProp(1)%poisson%omega
    rhs_coeff = fieldProp(1)%poisson%PB%rhs_coeff
    pot_diff = fieldProp(1)%poisson%pot_diff

    nodeloop: do iElem = 1, nSolve
      ! potential
      potential = auxField( iElem )

      do iDir = 1,QQ
        !> Calculate equilibrium distribution functions fEq
        eqstate = layout%weight(iDir)*potential

        !> Calculate source term
        !KM: Negative sign is because LBE solves potential equation of form
        ! \nabla^2 \phi - k^2 \phi = 0
        source = - layout%weight(iDir)*pot_diff*rhs_coeff*potential

        !> Generic fetching step:
        !! Streaming for pull
        !! Local copy for push
        pdfTmp = inState( neigh((idir-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

        !> Relaxation
        outState( ( ielem-1)* nscalars+idir+( 1-1)* qq )  &
          & = pdfTmp - omega*( pdfTmp - eqState ) + source
      end do ! iDir

    end do nodeloop

  end subroutine mus_PBLinear_advRel_generic
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Advection relaxation routine for the nonlinear
  !! poisson boltzmann equation for electrolyte solution
  !! $$ \nabla^2 \phi = - \frac{1}{\epsilon_r \epsilon_0}
  !!    \sum_i { e z_i c_{\infty} N_A exp(-\frac{z e }{k_b T}\phi) } $$
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_PBnonLinear_advRel_generic( fieldProp, inState, outState,    &
    &                                        auxField, neigh, nElems, nSolve, &
    &                                        level, layout, params, varSys,   &
    &                                        derVarPos                        )
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
    integer :: iElem, iDir, iIon  ! voxel element counter
    integer :: QQ, nScalars
    ! temporary distribution variables
    real(kind=rk) :: pdfTmp
    real(kind=rk) :: potential
    real(kind=rk) :: eqState, source, rhs
    real(kind=rk) :: omega, pot_diff, fac, pot_fac, permit_inv
    real(kind=rk) :: charge_dens
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nScalars = varSys%nScalars

    !>\todo KM 20170821: change omega for multilevel
    omega = fieldProp(1)%poisson%omega
    pot_diff = fieldProp(1)%poisson%pot_diff
    permit_inv = 1.0_rk / fieldProp(1)%poisson%permittivity
    fac = fieldProp(1)%poisson%PB%faradayLB          &
      &     / (fieldProp(1)%poisson%PB%gasConst_R_LB &
      &     * fieldProp(1)%poisson%PB%temp )

    nodeloop: do iElem = 1, nSolve
      do iDir = 1, QQ
      end do

      ! potential
      potential = auxField( iElem )

      ! charge density
      charge_dens = 0._rk
      pot_fac = fac * potential
      associate( PB => fieldProp(1)%poisson%PB  )
        do iIon = 1, PB%nIons
          charge_dens = charge_dens + PB%moleDens0 * PB%faradayLB            &
            &         * PB%valence(iIon) * exp( - PB%valence(iIon) * pot_fac )
        end do
      end associate
      ! Charge density of symmetric 1:1 electrolyte
      !charge_dens = -2.0_rk * fieldProp(1)%poisson%PB%moleDens0 &
      !  &         * fieldProp(1)%poisson%PB%faradayLB * sinh(pot_fac*potential)
      ! RHS
      rhs = charge_dens * permit_inv

      do iDir = 1,QQ
        !> Calculate equilibrium distribution functions fEq
        eqstate = layout%weight(iDir)*potential

        !> Calculate source term
        !KM: LBE solves potential equation of form
        ! \nabla^2 \phi + source = 0
        source = layout%weight(iDir) * pot_diff * rhs

        !> Generic fetching step:
        !! Streaming for pull
        !! Local copy for push
        pdfTmp = inState( neigh((idir-1)* nelems+ ielem)+( 1-1)* qq+ nscalars*0)

        !> Relaxation
        outState( (ielem-1)*nscalars+ idir+(1-1)*qq) &
          & = pdfTmp - omega*( pdfTmp - eqState ) + source
      end do ! iDir

    end do nodeloop

  end subroutine mus_PBnonLinear_advRel_generic
  ! ************************************************************************** !

end module mus_compute_Poisson_module
! **************************************************************************** !
