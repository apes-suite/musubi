! Copyright (c) 2021-2022 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
! **************************************************************************** !
!> Boundary condition wall treatment routines
!!
!! This module contains higher order wall treatments
!! A detailed description on the implementation details are given
!! in [[tem_bc_module]].
!!
module mus_bc_fluid_turbulent_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,               only: rk, long_k
  use tem_param_module,         only: cs2inv, cs2, rho0, rho0Inv, cs4inv, &
    &                                 div1_3
  use tem_isNaN_module,         only: tem_isNaN
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_debug_module,         only: dbgUnit
  use tem_geometry_module,      only: tem_ElemSizeLevel
  use tem_property_module,      only: prp_solid
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_float_module,         only: operator(.feq.), operator(.fne.)
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_math_module,          only: cross_product3D
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit

  ! include musubi modules
  use mus_bc_header_module,       only: boundary_type, glob_boundary_type
  use mus_scheme_layout_module,   only: mus_scheme_layout_type
  use mus_field_prop_module,      only: mus_field_prop_type
  use mus_derVarPos_module,       only: mus_derVarPos_type
  use mus_param_module,           only: mus_param_type
  use mus_physics_module,         only: mus_physics_type
  use mus_mixture_module,         only: mus_mixture_type
  use mus_varSys_module,          only: mus_varSys_data_type
  use mus_relaxationParam_module, only: mus_viscosity_type
  use mus_turb_wallFunc_module,   only: mus_turb_wallFunc_type
  use mus_turbulence_module,      only: mus_turbulence_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none

  private

  public :: turbulent_wall
  public :: turbulent_wall_libb
  public :: turbulent_wall_eq
  public :: turbulent_wall_eq_curved
  public :: turbulent_wall_noneq_expol
  public :: turbulent_wall_noneq_expol_curved
  !! Constant parameters for van-driest damping function
  real(kind=rk), parameter :: vd_Aplus = 26.0_rk

contains

  ! ************************************************************************** !
  !> BC routine for turbulent wall.
  !! It uses wall model to compute velocity on the boundary node.
  !! The implementation is based on the following paper:
  !! Haussmann, Marc; Ries, Florian; Jeppener-Haltenhoff, Jonathan B.; Li,
  !! Yongxiang; Schmidt, Marius; Welch, Cooper et al. (2020): Evaluation of a
  !! Near-Wall-Modeled Large Eddy Lattice Boltzmann Method for the Analysis of
  !! Complex Flows Relevant to IC Engines. In Computation 8 (2), p. 43.
  !! DOI: 10.3390/computation8020043.
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  {
  !!    label = 'wall',
  !!    kind = 'turbulent_wall',
  !!    wall_model = 'musker',
  !!    nonlinear_solver = 'fixed_point'
  !!  }
  !!}
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine turbulent_wall( me, state, bcBuffer, globBC,                &
    &          levelDesc, tree, nSize, iLevel, sim_time, neigh, layout,  &
    &          fieldProp, varPos, nScalars, varSys, derVarPos, physics,  &
    &          iField, mixture )
    ! ------------------------------------------------------------------------ !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> contains physics conversion factors
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! ------------------------------------------------------------------------ !
    integer :: iElem, iDir, QQ, elemPos
    real(kind=rk) :: f_pre(layout%fStencil%QQ)
    real(kind=rk) :: f_neigh(layout%fStencil%QQ)
    real(kind=rk) :: dens, velSW(3)
    real(kind=rk) :: dens_neigh, vel_neigh(3)
    real(kind=rk) :: fEq_bnd(layout%fStencil%QQ), fEq_neigh(layout%fStencil%QQ)
    real(kind=rk) :: unitSW(3, globBC%nElems(iLevel))
    real(kind=rk) :: velSW_bnd(globBC%nElems(iLevel))
    ! ------------------------------------------------------------------------ !
    ! Load stencil for velocity space (DdQq) with q=QQ
    QQ = layout%fStencil%QQ

    ! Calculate stream-wise velocity component from wall function
    call calcVelSW_unitSW_velTau_tVisc(                           &
      & velSW          = velSW_bnd,                               &
      & unitSW         = unitSW,                                  &
      & turbwallFunc   = me%turbwallFunc,                         &
      & nElems         = globBC%nElems(iLevel),                   &
      & elemPos        = globBC%elemLvl(iLevel)%elem%val(:),      &
      & neighBufferPre = me%neigh(iLevel)%neighBufferPre_nNext,   &
      & viscKine       = fieldProp%fluid%viscKine,                &
      & turbulence     = fieldProp%fluid%turbulence,              &
      & stencil        = layout%fStencil,                         &
      & iLevel         = iLevel,                                  &
      & quantities     = layout%quantities                        )

    ! Velocity correction on boundary element incoming direction
    do iElem = 1, globBC%nElems(iLevel)
      ! Current element position in state array. (Used for state in step 10.)
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
      ! Compute density and velocity from pre-collision state.
      ! computeNeighBuf uses FETCH which does implicit bounce back which
      ! is valid for qVal=0.5
      f_pre = me%neigh(iLevel)%computeNeighBuf( (iElem-1)*QQ+1: iElem*QQ )
      ! Compute density and velocity from pre-collision state on neighbor
      f_neigh = me%neigh(iLevel)%neighBufferPre_nNext( 1,         &
        &                          (iElem-1)*QQ+1:(iElem-1)*QQ+QQ )
      dens = sum(f_pre)
      dens_neigh = sum(f_neigh)
      ! velocity
      vel_neigh = layout%quantities%vel_from_pdf_ptr(pdf = f_neigh, dens = dens_neigh)

      ! stream-wise velocity on boundary element from wall model
      velSW = velSW_bnd(iElem) * unitSW(:, iElem)

      ! Equilibrium on precollision density and velocity
      fEq_neigh(:) = layout%quantities%pdfEq_ptr( rho = dens_neigh, &
        &                                         vel = vel_neigh,  &
        &                                         QQ = QQ           )

      ! Equilibrium on wall model velocity
      fEq_bnd(:) = layout%quantities%pdfEq_ptr( rho = dens,  &
        &                                       vel = velSW, &
        &                                       QQ = QQ      )

      do iDir = 1, layout%fStencil%QQ
        state(  neigh(( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
          & = fEq_bnd(iDir) + f_neigh(iDir) - fEq_neigh(iDir)
      end do
    end do

    ! Calculated bndForce on boundary elements using momentum exchange method
    call calcTurbWallBndForceAndMoment(                            &
      & bndForce    = me%turbWallFunc%dataOnLvl(iLevel)%bndForce,  &
      & bndMoment   = me%turbWallFunc%dataOnLvl(iLevel)%bndMoment, &
      & momRefPnt   = me%bndMomRefPnt,                             &
      & globBC      = globBC,                                      &
      & baryOfTotal = levelDesc%baryOfTotal,                       &
      & bcBuffer    = bcBuffer,                                    &
      & state       = state,                                       &
      & nSize       = nSize,                                       &
      & neigh       = neigh,                                       &
      & layout      = layout,                                      &
      & physics     = physics,                                     &
      & nScalars    = nScalars,                                    &
      & varPos      = varPos,                                      &
      & iLevel      = iLevel,                                      &
      & iField      = iField                                       )

  end subroutine turbulent_wall
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> BC routine for turbulent wall.
  !! It uses wall model to compute velocity on the boundary node.
  !! The implementation is based on the following paper:
  !! Haussmann, Marc; Ries, Florian; Jeppener-Haltenhoff, Jonathan B.; Li,
  !! Yongxiang; Schmidt, Marius; Welch, Cooper et al. (2020): Evaluation of a
  !! Near-Wall-Modeled Large Eddy Lattice Boltzmann Method for the Analysis of
  !! Complex Flows Relevant to IC Engines. In Computation 8 (2), p. 43.
  !! DOI: 10.3390/computation8020043.
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  {
  !!    label = 'wall',
  !!    kind = 'turbulent_wall',
  !!    curved = true,
  !!    wall_model = 'musker',
  !!    nonlinear_solver = 'fixed_point'
  !!  }
  !!}
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine turbulent_wall_libb( me, state, bcBuffer, globBC,           &
    &          levelDesc, tree, nSize, iLevel, sim_time, neigh, layout,  &
    &          fieldProp, varPos, nScalars, varSys, derVarPos, physics,  &
    &          iField, mixture )
    ! ------------------------------------------------------------------------ !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> contains physics conversion factors
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! ------------------------------------------------------------------------ !
    real(kind=rk) :: fIn, fOut, fNgh
    real(kind=rk) :: cIn, cOut, cNgh
    integer :: iLink
    integer :: iElem, iDir, QQ, elemPos
    real(kind=rk) :: f_pre(layout%fStencil%QQ)
    real(kind=rk) :: f_preBuffer(layout%fStencil%QQ * globBC%nElems(iLevel))
    real(kind=rk) :: dens, velSW(3), vel(3)
    real(kind=rk) :: fEq_bnd, fEq
    real(kind=rk) :: unitSW(3, globBC%nElems(iLevel))
    real(kind=rk) :: velSW_bnd(globBC%nElems(iLevel))
    ! ------------------------------------------------------------------------ !
    ! Load stencil for velocity space (DdQq) with q=QQ
    QQ = layout%fStencil%QQ

    ! 1. First perform libb step i.e. curved boundary step
    !NEC$ ivdep
    !DIR$ ivdep
    !IBM* independent
    do iLink = 1, me%links(iLevel)%nVals

      cIn  = me%bouzidi(iLevel)% cIn( iLink )
      cOut = me%bouzidi(iLevel)%cOut( iLink )
      cNgh = me%bouzidi(iLevel)%cNgh( iLink )

      fIn  = bcBuffer( me%bouzidi(iLevel)% inPos(iLink) )
      fOut = bcBuffer( me%bouzidi(iLevel)%outPos(iLink) )
      fNgh = me%neigh(iLevel)%computeNeighBuf(me%bouzidi(iLevel)%nghPos(iLink))

      state( me%links(iLevel)%val(iLink) ) = cIn*fIn + cOut*fOut + cNgh*fNgh

    end do ! iLink

    ! 2. Calculate stream-wise velocity component from wall function
    call calcVelSW_unitSW_velTau_tVisc(                           &
      & velSW          = velSW_bnd,                               &
      & unitSW         = unitSW,                                  &
      & turbwallFunc   = me%turbwallFunc,                         &
      & nElems         = globBC%nElems(iLevel),                   &
      & elemPos        = globBC%elemLvl(iLevel)%elem%val(:),      &
      & neighBufferPre = me%neigh(iLevel)%neighBufferPre_nNext,   &
      & viscKine       = fieldProp%fluid%viscKine,                &
      & turbulence     = fieldProp%fluid%turbulence,              &
      & stencil        = layout%fStencil,                         &
      & iLevel         = iLevel,                                  &
      & quantities     = layout%quantities                        )

    ! 3. Velocity correction step
    ! Get pre-collision of boundary elementy in a buffer since it gets
    ! overwritten in the next element loop
    do iElem = 1, globBC%nElems(iLevel)
      ! Current element position in state array. (Used for state in step 10.)
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
      do iDir = 1, layout%fStencil%QQ
        f_preBuffer((iElem-1)*QQ + iDir) = state(                        &
          &  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )
      end do
    end do

    do iElem = 1, globBC%nElems(iLevel)
      ! Current element position in state array. (Used for state in step 10.)
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
      ! Compute density and velocity from pre-collision state after
      ! curved boundary step
      f_pre = f_preBuffer( (iElem-1)*QQ+1: iElem*QQ )
      dens = sum(f_pre)
      ! velocity
      vel = layout%quantities%vel_from_pdf_ptr(pdf = f_pre, dens = dens)

      ! stream-wise velocity on boundary element from wall model
      velSW = velSW_bnd(iElem) * unitSW(:, iElem)

      do iDir = 1, layout%fStencil%QQN
        if ( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Equilibrium on wall model velocity
          fEq_bnd = layout%quantities%pdfEq_iDir_ptr( rho = dens,                                &
            &                                         vel = velSW,                               &
            &                                         iDir = iDir,                               &
            &                                         cxDirRK = layout%fStencil%cxDirRK(:,iDir), &
            &                                         weight = layout%weight(iDir)               )

          fEq = layout%quantities%pdfEq_iDir_ptr( rho = dens,                                &
            &                                     vel = vel,                                 &
            &                                     iDir = iDir,                               &
            &                                     cxDirRK = layout%fStencil%cxDirRK(:,iDir), &
            &                                     weight = layout%weight(iDir)               )

          state(  neigh(( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq_bnd + f_pre(iDir) - fEq
        end if
      end do
    end do

    ! Calculated bndForce on boundary elements using momentum exchange method
    call calcTurbWallBndForceAndMoment(                            &
      & bndForce    = me%turbWallFunc%dataOnLvl(iLevel)%bndForce,  &
      & bndMoment   = me%turbWallFunc%dataOnLvl(iLevel)%bndMoment, &
      & momRefPnt   = me%bndMomRefPnt,                             &
      & globBC      = globBC,                                      &
      & baryOfTotal = levelDesc%baryOfTotal,                       &
      & bcBuffer    = bcBuffer,                                    &
      & state       = state,                                       &
      & nSize       = nSize,                                       &
      & neigh       = neigh,                                       &
      & layout      = layout,                                      &
      & physics     = physics,                                     &
      & nScalars    = nScalars,                                    &
      & varPos      = varPos,                                      &
      & iLevel      = iLevel,                                      &
      & iField      = iField                                       )

  end subroutine turbulent_wall_libb
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> BC routine for turbulent wall based on Guo's nonequilibrium extrapolation.
  !! The implementation is based on the following paper:
  !! Haussmann, M. et al. (2019) ‘Large-eddy simulation coupled with wall models
  !! for turbulent channel flows at high Reynolds numbers with a lattice
  !! Boltzmann method — Application to Coriolis mass flowmeter’, Computers &
  !! Mathematics with Applications. Elsevier Ltd, 78(10), pp. 3285–3302.
  !!
  !! It uses wall model to compute velocity on the boundary node.
  !! All directions of PDF in the boundary elements are updated with
  !! Equilibrium plus non-equilibrium.
  !! Density is computed using Zho-He approach for straight walls.
  !! "On pressure and velocity boundary conditions for the lattice Boltzmann
  !!  BGK model", Physics of Fluids 9, 1591-1598 (1997)
  !! https://doi.org/10.1063/1.869307
  !!
  !! non-equilibrium are computed from PDF on neighbor and extrapolated to
  !! boundary. This routine is used for straight wall boundaries.
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  {
  !!    label = 'wall',
  !!    kind = 'turbulent_wall_noneq_expol',
  !!    wall_model = 'musker',
  !!    nonlinear_solver = 'fixed_point'
  !!  }
  !!}
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine turbulent_wall_noneq_expol( me, state, bcBuffer, globBC,    &
    &          levelDesc, tree, nSize, iLevel, sim_time, neigh, layout,  &
    &          fieldProp, varPos, nScalars, varSys, derVarPos, physics,  &
    &          iField, mixture )
    ! ------------------------------------------------------------------------ !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> contains physics conversion factors
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! ------------------------------------------------------------------------ !
    integer :: iElem, iDir, QQ, elemPos, normalInd_inv
    real(kind=rk) :: f_pre(layout%fStencil%QQ)
    real(kind=rk) :: f_neigh(layout%fStencil%QQ)
    real(kind=rk) :: dens_neigh, velSW(3), vel_neigh(3)
    real(kind=rk) :: dens_bnd, vel_normal, normal(3)
    real(kind=rk) :: fEq(layout%fStencil%QQ), fEq_bnd(layout%fStencil%QQ)
    real(kind=rk) :: unitSW(3, globBC%nElems(iLevel))
    real(kind=rk) :: velSW_bnd(globBC%nElems(iLevel))
    ! ------------------------------------------------------------------------ !
    ! Load stencil for velocity space (DdQq) with q=QQ
    QQ = layout%fStencil%QQ

    ! Calculate stream-wise velocity component from wall function
    call calcVelSW_unitSW_velTau_tVisc(                           &
      & velSW          = velSW_bnd,                               &
      & unitSW         = unitSW,                                  &
      & turbwallFunc   = me%turbwallFunc,                         &
      & nElems         = globBC%nElems(iLevel),                   &
      & elemPos        = globBC%elemLvl(iLevel)%elem%val(:),      &
      & neighBufferPre = me%neigh(iLevel)%neighBufferPre_nNext,   &
      & viscKine       = fieldProp%fluid%viscKine,                &
      & turbulence     = fieldProp%fluid%turbulence,              &
      & stencil        = layout%fStencil,                         &
      & iLevel         = iLevel,                                  &
      & quantities     = layout%quantities                        )

    ! Velocity correction on boundary element incoming direction
    do iElem = 1, globBC%nElems(iLevel)
      ! Current element position in state array. (Used for state in step 10.)
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
      ! Compute density and velocity from pre-collision state on neighbor
      f_neigh = me%neigh(iLevel)%neighBufferPre_nNext( 1,         &
        &                          (iElem-1)*QQ+1:(iElem-1)*QQ+QQ )
      dens_neigh = sum(f_neigh)
      ! velocity
      vel_neigh = layout%quantities%vel_from_pdf_ptr(pdf = f_neigh, dens = dens_neigh)

      ! stream-wise velocity on boundary element from wall model
      velSW = velSW_bnd(iElem) * unitSW(:, iElem)

      ! calculate density at boundary for straight wall using Zho-He approach
      normalInd_inv = layout%fStencil%cxDirInv( globBC%elemLvl(iLevel) &
        &                   %normalInd%val(iElem) )
      normal = layout%fStencil%cxDirRK(:, normalInd_inv)
      vel_normal = dot_product(velSW, normal)
      dens_bnd = 0.0_rk
      f_pre = me%neigh(iLevel)%computeNeighBuf( (iElem-1)*QQ+1: iElem*QQ )
      do iDir = 1, layout%fStencil%QQ
        if ( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          dens_bnd = dens_bnd + f_pre(layout%fStencil%cxDirInv(iDir))
        else
          dens_bnd = dens_bnd + f_pre(iDir)
        end if
      end do
      dens_bnd = dens_bnd / (1.0_rk + vel_normal)

      ! Equilibrium on precollision density and velocity
      fEq(:) = layout%quantities%pdfEq_ptr( rho = dens_neigh, &
        &                                   vel = vel_neigh,  &
        &                                   QQ = QQ           )

      ! Equilibrium on wall model velocity
      fEq_bnd(:) = layout%quantities%pdfEq_ptr( rho = dens_bnd,  &
        &                                       vel = velSW,     &
        &                                       QQ = QQ          )

      do iDir = 1, layout%fStencil%QQ
        state(  neigh(( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
          & = fEq_bnd(iDir) + ( f_neigh(iDir) - fEq(iDir) )
      end do
    end do

    ! Calculated bndForce on boundary elements using momentum exchange method
    call calcTurbWallBndForceAndMoment(                            &
      & bndForce    = me%turbWallFunc%dataOnLvl(iLevel)%bndForce,  &
      & bndMoment   = me%turbWallFunc%dataOnLvl(iLevel)%bndMoment, &
      & momRefPnt   = me%bndMomRefPnt,                             &
      & globBC      = globBC,                                      &
      & baryOfTotal = levelDesc%baryOfTotal,                       &
      & bcBuffer    = bcBuffer,                                    &
      & state       = state,                                       &
      & nSize       = nSize,                                       &
      & neigh       = neigh,                                       &
      & layout      = layout,                                      &
      & physics     = physics,                                     &
      & nScalars    = nScalars,                                    &
      & varPos      = varPos,                                      &
      & iLevel      = iLevel,                                      &
      & iField      = iField                                       )

  end subroutine turbulent_wall_noneq_expol
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> BC routine for turbulent wall based on Guo's nonequilibrium extrapolation.
  !! The implementation is based on the following paper:
  !! Haussmann, M. et al. (2019) ‘Large-eddy simulation coupled with wall models
  !! for turbulent channel flows at high Reynolds numbers with a lattice
  !! Boltzmann method — Application to Coriolis mass flowmeter’, Computers &
  !! Mathematics with Applications. Elsevier Ltd, 78(10), pp. 3285–3302.
  !!
  !! It uses wall model to compute velocity on the boundary node.
  !! All directions of PDF in the boundary elements are updated with
  !! Equilibrium plus non-equilibrium. Density and non-equilibrium are commputed
  !! from PDF on neighbor and extrapolated to boundary.
  !! This routine is used for curved boundaries.
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  {
  !!    label = 'wall',
  !!    kind = 'turbulent_wall_noneq_expol',
  !!    curved = true,
  !!    wall_model = 'musker',
  !!    nonlinear_solver = 'fixed_point'
  !!  }
  !!}
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine turbulent_wall_noneq_expol_curved( me, state, bcBuffer, globBC, &
    &          levelDesc, tree, nSize, iLevel, sim_time, neigh, layout,      &
    &          fieldProp, varPos, nScalars, varSys, derVarPos, physics,      &
    &          iField, mixture )
    ! ------------------------------------------------------------------------ !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> contains physics conversion factors
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! ------------------------------------------------------------------------ !
    integer :: iElem, iDir, QQ, elemPos
    real(kind=rk) :: f_neigh(layout%fStencil%QQ)
    real(kind=rk) :: dens_neigh, velSW(3), vel_neigh(3)
    real(kind=rk) :: fEq(layout%fStencil%QQ), fEq_bnd(layout%fStencil%QQ)
    real(kind=rk) :: unitSW(3, globBC%nElems(iLevel))
    real(kind=rk) :: velSW_bnd(globBC%nElems(iLevel))
    ! ------------------------------------------------------------------------ !
    ! Load stencil for velocity space (DdQq) with q=QQ
    QQ = layout%fStencil%QQ

    ! Calculate stream-wise velocity component from wall function
    call calcVelSW_unitSW_velTau_tVisc(                           &
      & velSW          = velSW_bnd,                               &
      & unitSW         = unitSW,                                  &
      & turbwallFunc   = me%turbwallFunc,                         &
      & nElems         = globBC%nElems(iLevel),                   &
      & elemPos        = globBC%elemLvl(iLevel)%elem%val(:),      &
      & neighBufferPre = me%neigh(iLevel)%neighBufferPre_nNext,   &
      & viscKine       = fieldProp%fluid%viscKine,                &
      & turbulence     = fieldProp%fluid%turbulence,              &
      & stencil        = layout%fStencil,                         &
      & iLevel         = iLevel,                                  &
      & quantities     = layout%quantities                        )

    ! Velocity correction on boundary element incoming direction
    do iElem = 1, globBC%nElems(iLevel)
      ! Current element position in state array. (Used for state in step 10.)
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
      ! Compute density and velocity from pre-collision state on neighbor
      f_neigh = me%neigh(iLevel)%neighBufferPre_nNext( 1,         &
        &                          (iElem-1)*QQ+1:(iElem-1)*QQ+QQ )
      dens_neigh = sum(f_neigh)
      ! velocity
      vel_neigh = layout%quantities%vel_from_pdf_ptr(pdf = f_neigh, dens = dens_neigh)

      ! stream-wise velocity on boundary element from wall model
      velSW = velSW_bnd(iElem) * unitSW(:, iElem)

      ! Equilibrium on precollision density and velocity
      fEq(:) = layout%quantities%pdfEq_ptr( rho = dens_neigh, &
        &                                   vel = vel_neigh,  &
        &                                   QQ = QQ           )

      ! Equilibrium on wall model velocity
      fEq_bnd(:) = layout%quantities%pdfEq_ptr( rho = dens_neigh, &
        &                                       vel = velSW,      &
        &                                       QQ = QQ           )

      do iDir = 1, layout%fStencil%QQ
        state(  neigh(( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
          & = fEq_bnd(iDir) + ( f_neigh(iDir) - fEq(iDir) )
      end do
    end do

    ! Calculated bndForce on boundary elements using momentum exchange method
    call calcTurbWallBndForceAndMoment(                            &
      & bndForce    = me%turbWallFunc%dataOnLvl(iLevel)%bndForce,  &
      & bndMoment   = me%turbWallFunc%dataOnLvl(iLevel)%bndMoment, &
      & momRefPnt   = me%bndMomRefPnt,                             &
      & globBC      = globBC,                                      &
      & baryOfTotal = levelDesc%baryOfTotal,                       &
      & bcBuffer    = bcBuffer,                                    &
      & state       = state,                                       &
      & nSize       = nSize,                                       &
      & neigh       = neigh,                                       &
      & layout      = layout,                                      &
      & physics     = physics,                                     &
      & nScalars    = nScalars,                                    &
      & varPos      = varPos,                                      &
      & iLevel      = iLevel,                                      &
      & iField      = iField                                       )

  end subroutine turbulent_wall_noneq_expol_curved
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> BC routine for turbulent wall based on equilibrium BC.
  !! The implementation is based on the following paper:
  !! Haussmann, M. et al. (2019) ‘Large-eddy simulation coupled with wall models
  !! for turbulent channel flows at high Reynolds numbers with a lattice
  !! Boltzmann method — Application to Coriolis mass flowmeter’, Computers &
  !! Mathematics with Applications. Elsevier Ltd, 78(10), pp. 3285–3302.
  !!
  !! It uses wall model to compute velocity on the boundary node.
  !! All directions of PDF in the boundary elements are updated with
  !! Equilibrium.
  !! Density is computed using Zho-He approach for straight walls.
  !! "On pressure and velocity boundary conditions for the lattice Boltzmann
  !!  BGK model", Physics of Fluids 9, 1591-1598 (1997)
  !! https://doi.org/10.1063/1.869307
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  {
  !!    label = 'wall',
  !!    kind = 'turbulent_wall_eq',
  !!    wall_model = 'musker',
  !!    nonlinear_solver = 'fixed_point'
  !!  }
  !!}
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine turbulent_wall_eq( me, state, bcBuffer, globBC, levelDesc, tree, &
    &          nSize, iLevel, sim_time, neigh, layout, fieldProp, varPos,     &
    &          nScalars, varSys, derVarPos, physics, iField, mixture )
    ! ------------------------------------------------------------------------ !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> contains physics conversion factors
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! ------------------------------------------------------------------------ !
    integer :: iElem, iDir, QQ, elemPos, normalInd_inv
    real(kind=rk) :: f_pre(layout%fStencil%QQ)
    real(kind=rk) :: velSW(3)
    real(kind=rk) :: dens_bnd, vel_normal, normal(3)
    real(kind=rk) :: fEq_bnd(layout%fStencil%QQ)
    real(kind=rk) :: unitSW(3, globBC%nElems(iLevel))
    real(kind=rk) :: velSW_bnd(globBC%nElems(iLevel))
    ! ------------------------------------------------------------------------ !
    ! Load stencil for velocity space (DdQq) with q=QQ
    QQ = layout%fStencil%QQ

    ! Calculate stream-wise velocity component from wall function
    call calcVelSW_unitSW_velTau_tVisc(                           &
      & velSW          = velSW_bnd,                               &
      & unitSW         = unitSW,                                  &
      & turbwallFunc   = me%turbwallFunc,                         &
      & nElems         = globBC%nElems(iLevel),                   &
      & elemPos        = globBC%elemLvl(iLevel)%elem%val(:),      &
      & neighBufferPre = me%neigh(iLevel)%neighBufferPre_nNext,   &
      & viscKine       = fieldProp%fluid%viscKine,                &
      & turbulence     = fieldProp%fluid%turbulence,              &
      & stencil        = layout%fStencil,                         &
      & iLevel         = iLevel,                                  &
      & quantities     = layout%quantities                        )

    ! Velocity correction on boundary element incoming direction
    do iElem = 1, globBC%nElems(iLevel)
      ! Current element position in state array. (Used for state in step 10.)
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)

      ! stream-wise velocity on boundary element from wall model
      velSW = velSW_bnd(iElem) * unitSW(:, iElem)

      ! calculate density at boundary for straight wall using Zho-He approach
      normalInd_inv = layout%fStencil%cxDirInv( globBC%elemLvl(iLevel) &
        &                   %normalInd%val(iElem) )
      normal = layout%fStencil%cxDirRK(:, normalInd_inv)
      vel_normal = dot_product(velSW, normal)
      dens_bnd = 0.0_rk
      f_pre = me%neigh(iLevel)%computeNeighBuf( (iElem-1)*QQ+1: iElem*QQ )
      do iDir = 1, layout%fStencil%QQ
        if ( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          dens_bnd = dens_bnd + f_pre(layout%fStencil%cxDirInv(iDir))
        else
          dens_bnd = dens_bnd + f_pre(iDir)
        end if
      end do
      dens_bnd = dens_bnd / (1.0_rk + vel_normal)

      ! Equilibrium on wall model velocity
      fEq_bnd(:) = layout%quantities%pdfEq_ptr( rho = dens_bnd,    &
        &                                       vel = velSW,       &
        &                                       QQ = QQ            )

      do iDir = 1, layout%fStencil%QQ
        state(  neigh(( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
          & = fEq_bnd(iDir)
      end do
    end do

    ! Calculated bndForce on boundary elements using momentum exchange method
    call calcTurbWallBndForceAndMoment(                            &
      & bndForce    = me%turbWallFunc%dataOnLvl(iLevel)%bndForce,  &
      & bndMoment   = me%turbWallFunc%dataOnLvl(iLevel)%bndMoment, &
      & momRefPnt   = me%bndMomRefPnt,                             &
      & globBC      = globBC,                                      &
      & baryOfTotal = levelDesc%baryOfTotal,                       &
      & bcBuffer    = bcBuffer,                                    &
      & state       = state,                                       &
      & nSize       = nSize,                                       &
      & neigh       = neigh,                                       &
      & layout      = layout,                                      &
      & physics     = physics,                                     &
      & nScalars    = nScalars,                                    &
      & varPos      = varPos,                                      &
      & iLevel      = iLevel,                                      &
      & iField      = iField                                       )

  end subroutine turbulent_wall_eq
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> BC routine for turbulent wall based on equilibrium BC.
  !! The implementation is based on the following paper:
  !! Haussmann, M. et al. (2019) ‘Large-eddy simulation coupled with wall models
  !! for turbulent channel flows at high Reynolds numbers with a lattice
  !! Boltzmann method — Application to Coriolis mass flowmeter’, Computers &
  !! Mathematics with Applications. Elsevier Ltd, 78(10), pp. 3285–3302.
  !!
  !! It uses wall model to compute velocity on the boundary node.
  !! All directions of PDF in the boundary elements are updated with
  !! Equilibrium plus non-equilibrium. Density is commputed
  !! from PDF on neighbor and extrapolated to boundary.
  !! This routine is used for curved boundaries.
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  {
  !!    label = 'wall',
  !!    kind = 'turbulent_wall_eq',
  !!    curved = true,
  !!    wall_model = 'musker',
  !!    nonlinear_solver = 'fixed_point'
  !!  }
  !!}
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine turbulent_wall_eq_curved( me, state, bcBuffer, globBC,     &
    &          levelDesc, tree, nSize, iLevel, sim_time, neigh, layout, &
    &          fieldProp, varPos, nScalars, varSys, derVarPos, physics, &
    &          iField, mixture )
    ! ------------------------------------------------------------------------ !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> contains physics conversion factors
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! ------------------------------------------------------------------------ !
    integer :: iElem, iDir, QQ, elemPos
    real(kind=rk) :: f_neigh(layout%fStencil%QQ)
    real(kind=rk) :: dens_neigh, velSW(3)
    real(kind=rk) :: fEq_bnd(layout%fStencil%QQ)
    real(kind=rk) :: unitSW(3, globBC%nElems(iLevel))
    real(kind=rk) :: velSW_bnd(globBC%nElems(iLevel))
    ! ------------------------------------------------------------------------ !
    ! Load stencil for velocity space (DdQq) with q=QQ
    QQ = layout%fStencil%QQ

    ! Calculate stream-wise velocity component from wall function
    call calcVelSW_unitSW_velTau_tVisc(                           &
      & velSW          = velSW_bnd,                               &
      & unitSW         = unitSW,                                  &
      & turbwallFunc   = me%turbwallFunc,                         &
      & nElems         = globBC%nElems(iLevel),                   &
      & elemPos        = globBC%elemLvl(iLevel)%elem%val(:),      &
      & neighBufferPre = me%neigh(iLevel)%neighBufferPre_nNext,   &
      & viscKine       = fieldProp%fluid%viscKine,                &
      & turbulence     = fieldProp%fluid%turbulence,              &
      & stencil        = layout%fStencil,                         &
      & iLevel         = iLevel,                                  &
      & quantities     = layout%quantities                        )

    ! Velocity correction on boundary element incoming direction
    do iElem = 1, globBC%nElems(iLevel)
      ! Current element position in state array. (Used for state in step 10.)
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
      ! Compute density and velocity from pre-collision state on neighbor
      f_neigh = me%neigh(iLevel)%neighBufferPre_nNext( 1,         &
        &                          (iElem-1)*QQ+1:(iElem-1)*QQ+QQ )
      dens_neigh = sum(f_neigh)

      ! stream-wise velocity on boundary element from wall model
      velSW = velSW_bnd(iElem) * unitSW(:, iElem)

      ! Equilibrium on wall model velocity
      fEq_bnd(:) = layout%quantities%pdfEq_ptr( rho = dens_neigh, &
        &                                       vel = velSW,      &
        &                                       QQ = QQ           )

      do iDir = 1, layout%fStencil%QQ
        state(  neigh(( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
          & = fEq_bnd(iDir)
      end do
    end do

    ! Calculated bndForce on boundary elements using momentum exchange method
    call calcTurbWallBndForceAndMoment(                            &
      & bndForce    = me%turbWallFunc%dataOnLvl(iLevel)%bndForce,  &
      & bndMoment   = me%turbWallFunc%dataOnLvl(iLevel)%bndMoment, &
      & momRefPnt   = me%bndMomRefPnt,                             &
      & globBC      = globBC,                                      &
      & baryOfTotal = levelDesc%baryOfTotal,                       &
      & bcBuffer    = bcBuffer,                                    &
      & state       = state,                                       &
      & nSize       = nSize,                                       &
      & neigh       = neigh,                                       &
      & layout      = layout,                                      &
      & physics     = physics,                                     &
      & nScalars    = nScalars,                                    &
      & varPos      = varPos,                                      &
      & iLevel      = iLevel,                                      &
      & iField      = iField                                       )


  end subroutine turbulent_wall_eq_curved
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculation stream-wise velocity compononent from wall function and
  !! friction velocity, stream-wise unit vector and turbulent viscosity with
  !! mixing length formulation.
  subroutine calcVelSW_unitSW_velTau_tVisc(velSW, unitSW, turbwallFunc, &
    &          nElems, elemPos, neighBufferPre, viscKine, turbulence,   &
    &          stencil, iLevel, quantities)
    ! --------------------------------------------------------------------------
    !> Stream-wise velocity component from wall function
    real(kind=rk), intent(out) :: velSW(:)
    !> Stream-wise unit vector
    real(kind=rk), intent(out) :: unitSW(:,:)
    !> Turbulent wall model type contains viscosity, velTau, distToBnd and
    !! function pointers to compute velTau and velSW
    type(mus_turb_wallFunc_type), intent(inout) :: turbwallFunc
    !> Number of elements in current boundary
    integer, intent(in) :: nElems
    !> Current element position in state array. (Used to access viscosity)
    integer, intent(in) :: elemPos(:)
    !> Pre-collision from 1st and 2nd fluid neighbor
    real(kind=rk), intent(in) :: neighBufferPre(:,:)
    !> Kinematic viscosity
    type(mus_viscosity_type) :: viscKine
    !> turbulence model type
    type(mus_turbulence_type), intent(in) :: turbulence
    !> Fluid stencil
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> Current level
    integer, intent(in) :: iLevel
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! --------------------------------------------------------------------------
    integer :: iElem, QQ
    real(kind=rk) :: f_neigh(stencil%QQ)
    real(kind=rk) :: vel_neigh(3), dens_neigh, vec(3), vec_mag
    real(kind=rk) :: velSW_neigh(nElems)
    real(kind=rk) :: yPlus, gradU, velSW_neigh_2, unitNormal(3)
    ! --------------------------------------------------------------------------
    QQ = stencil%QQ
    ! Initialize unit streamwise to zero. It is set to non-zero value only
    ! if vec_mag is not zero
    unitSW = 0.0_rk
    ! Calculate stream-wise velocity component on first neighbor
    do iElem = 1, nElems
      ! Get density and velocity on second fluid element along the normal
      ! direction i.e. 1st neighbor.
      f_neigh = neighBufferPre(1, (iElem-1)*QQ+1:(iElem-1)*QQ+QQ)
      dens_neigh = sum(f_neigh)

      ! Get velocity from pdf, inverse of density and stencil as input.
      vel_neigh = quantities%vel_from_pdf_ptr( pdf = f_neigh, dens = dens_neigh )

      ! Calculate local stream-wise unit vector
      ! e' = (u - (u . n) . n) / |(u - (u . n) . n|
      !
      unitNormal = turbwallFunc%dataOnLvl(iLevel)%unitNormal(:, iElem)
      vec = vel_neigh - dot_product( vel_neigh, unitNormal ) * unitNormal
      vec_mag = sqrt(dot_product(vec, vec))

      if ( tem_isnan(vec_mag) ) then
        call tem_abort("Error: In calcVelSW_unitSW_velTau_tVisc. " &
          &          //"Stream-wise velocity mag is NaN.")
      end if
      ! Unit vector
      if (vec_mag .fne. 0.0_rk) then
        unitSW(:, iElem) = vec / vec_mag
      end if
      ! stream-wise velocity component
      velSW_neigh(iElem) = dot_product(vel_neigh, unitSW(:, iElem))
      if (velSW_neigh(iElem) < 0.0_rk) then
        call tem_abort("Error: In calcVelSW_unitSW_velTau_tVisc. " &
          &          //"Stream-wise velocity component is negative.")
      end if
    end do

    ! calculate friction velocity on neighbor element
    call turbwallFunc%calcFricVel(                                 &
      & velTau    = turbwallFunc%dataOnLvl(iLevel)%velTau,         &
      & velSW     = velSW_neigh,                                   &
      & distToBnd = turbwallFunc%dataOnLvl(iLevel)%neighDistToBnd, &
      & viscKine  = viscKine%dataOnLvl(iLevel)%val( elemPos(:) ),  &
      & nElems    = nElems                                         )

    ! calculate stream-wise velocity on boundary element
    call turbwallFunc%calcStreamWiseVel(                          &
      & velSW     = velSW,                                        &
      & velTau    = turbwallFunc%dataOnLvl(iLevel)%velTau,        &
      & distToBnd = turbwallFunc%dataOnLvl(iLevel)%distToBnd,     &
      & viscKine  = viscKine%dataOnLvl(iLevel)%val( elemPos(:) ), &
      & nElems    = nElems,                                       &
      & wall_function = turbwallFunc%wall_function                )

    ! Calculate Turbulent viscosity according to mixing length formulation
    ! with von karman constant.
    ! only if LES turbulence model is smagorinsky because turb viscosity
    ! from Vreman and WALE reduces towards the wall
    ! nu_t = (k*y)**2 * |dudy|
    if (turbulence%active .and.                       &
      & trim(turbulence%config%model) == 'smagorinsky') then
      do iElem = 1, nElems
        ! Get density and velocity on third fluid element along the normal
        ! direction i.e. 2nd neighbor.
        f_neigh = neighBufferPre(2, (iElem-1)*QQ+1:(iElem-1)*QQ+QQ)
        dens_neigh = sum(f_neigh)

        ! Get velocity from pdf, inverse of density and stencil as input.
        vel_neigh = quantities%vel_from_pdf_ptr( pdf = f_neigh, dens = dens_neigh )
        velSW_neigh_2 = dot_product(vel_neigh, unitSW(:, iElem))

        ! computed velocity gradient from second order forward difference
        gradU = abs( - 3.0_rk * velSW(iElem) + 4.0_rk * velSW_neigh(iElem) &
          & - velSW_neigh_2 ) / 2.0_rk
        turbwallFunc%dataOnLvl(iLevel)%tVisc(iElem) = ( turbWallFunc%vonKarman &
          & * turbwallFunc%dataOnLvl(iLevel)%distToBnd(iElem) )**2.0_rk        &
          & * gradU
      end do

      ! Apply van Driest damping
      ! nu_t = nu_t * (1-exp(-yplus/Aplus))**2
      if (turbWallFunc%useVanDriest) then
        do iElem = 1, nElems
          yPlus = turbwallFunc%dataOnLvl(iLevel)%distToBnd(iElem) &
            & * turbwallFunc%dataOnLvl(iLevel)%velTau(iElem)      &
            & / viscKine%dataOnLvl(iLevel)%val(elemPos(iElem))
          turbwallFunc%dataOnLvl(iLevel)%tVisc(iElem)       &
            & = turbwallFunc%dataOnLvl(iLevel)%tVisc(iElem) &
            & * (1.0_rk-exp(-yPlus/vd_Aplus))**2.0_rk
        end do
      end if
    end if
  end subroutine calcVelSW_unitSW_velTau_tVisc
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine computes bndForce on boundary elements using momentum
  !! exchange method.
  subroutine calcTurbWallBndForceAndMoment(bndForce, bndMoment, momRefPnt,     &
    &                                      globBC, baryOfTotal, bcBuffer,      &
    &                                      state, nSize, neigh, layout,        &
    &                                      physics, nScalars, varPos, iLevel,  &
    &                                      iField)
    ! ------------------------------------------------------------------------ !
    !> Boundary force
    real(kind=rk), intent(out) :: bndForce(:,:)
    !> Boundary moment
    real(kind=rk), intent(out) :: bndMoment(:,:)
    !> Reference point for moment calculation
    real(kind=rk), intent(in) :: momRefPnt(:)
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> barycenter of elements in total list
    !! Size: nElemsTotal, 3
    real(kind=rk), intent(in) :: baryOfTotal(:,:)
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> contains physics conversion factors
    type(mus_physics_type), intent(in) :: physics
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> current field
    integer, intent(in) :: iField
    ! ------------------------------------------------------------------------ !
    integer :: iElem, iDir, QQ, invDir, elemPos, posInBuffer
    real(kind=rk) :: fIn, fOut, conv2LatLen
    real(kind=rk) :: force(3), delta_x(3)
    ! ------------------------------------------------------------------------ !
    conv2LatLen = 1.0_rk / physics%fac(iLevel)%length
    ! Load stencil for velocity space (DdQq) with q=QQ
    QQ = layout%fStencil%QQ
    do iElem = 1, globBC%nElems(iLevel)
      force = 0.0_rk
      ! Current element position in state array
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
      ! position in bcBuffer
      posInBuffer = globBC%elemLvl(iLevel)%posInBcElemBuf%val( iElem )

      ! arm of the momentum
      delta_x = (baryOfTotal(elemPos, :) - momRefPnt(:)) * conv2LatLen

      do iDir = 1, layout%fStencil%QQN
        if ( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! direction towards boundary
          invDir  = layout%fStencil%cxDirInv(iDir)
          ! Incoming state before collision
          fIn = state(                                                   &
            &  neigh((idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )
          ! outgoing state after collision
          fOut = bcBuffer( (posInBuffer-1)*nScalars+varPos(invDir) )
          ! momentum exchange
          force = force + layout%fStencil%cxDirRK(:, invDir) * (fOut + fIn)
        end if
      end do
      bndForce(:, iElem) = force
      ! moment = arm x force
      bndMoment(:, iElem) = cross_product3D(delta_x, force)
    end do
  end subroutine calcTurbWallBndForceAndMoment
  ! ************************************************************************** !

end module mus_bc_fluid_turbulent_module
! *****************************************************************************!
