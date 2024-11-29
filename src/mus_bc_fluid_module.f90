! Copyright (c) 2012-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015-2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017-2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2018 Jana Gericke <jana.gericke@uni-siegen.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2021-2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> author: Manuel Hasert
!! author: Kannan Masilamani
!! Boundary condition treatment routines for fluid simulation
!! Details on implementation can be found in [[tem_bc_module]]
!!
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
module mus_bc_fluid_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,               only: rk
  use tem_param_module,         only: div1_3, div4_3, csInv, cs, cs2, cs2inv, &
    &                                 rho0, rho0Inv
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_property_module,      only: prp_hasQVal
  ! use tem_subTree_type_module,  only: tem_subTree_type
  use tem_geometry_module,      only: tem_BaryOfId,tem_ElemSizeLevel
  use tem_logging_module,       only: logUnit
  use tem_debug_module,         only: dbgUnit
  use tem_construction_module,  only: tem_levelDesc_type

  ! include musubi modules
  use mus_bc_header_module,          only: boundary_type, glob_boundary_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_param_module,              only: mus_param_type
  use mus_physics_module,            only: mus_physics_type
  use mus_mixture_module,            only: mus_mixture_type
  use mus_varSys_module,             only: mus_varSys_data_type

  implicit none

  private

  public :: velocity_eq
  public :: mfr_bounceback
  public :: mfr_eq
  public :: velocity_bounceback_incomp, velocity_bounceback
  public :: velocity_bfl, velocity_bfl_incomp
  public :: vel_neq

  public :: bc_pdf

  public :: press_neq
  public :: pressure_antiBounceBack
  public :: pressure_eq
  public :: outlet_dnt
  public :: pressure_expol
  public :: outlet_zero_prsgrd

  public :: inlet_nrbc_incomp, inlet_nrbc
  public :: outlet_nrbc_incomp, outlet_nrbc
  public :: outlet_nrbc_eq

  ! moments based BC
  public :: velocity_momentsbased
  public :: pressure_momentsbased
  public :: moments_wall

  public :: velocity_momentsbased_incomp
  public :: pressure_momentsbased_incomp

contains

! ****************************************************************************** !
  !> Inlet Velocity Equilibrium type boundary conditions for weakly
  !! compressible lbm scheme
  !!
  !! The incoming densities are set to the equilibrium distribution with
  !! macroscopic values:
  !! - density is extrapolated (0 or 1st order) from the flow domain
  !!   + 0-order: so we simply use the density of the boundary element itself
  !! - velocity is used as defined in the configuration file
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'inlet',
  !!    kind = 'velocity_eq',
  !!    velocity = 'inlet_vel',
  !!  }
  !!}
  !!variable = {
  !!  name = 'inlet_vel',
  !!  ncomponents = 3,
  !!  vartype = 'st_fun',
  !!  st_fun = {0.06, 0.0, 0.0}
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine velocity_eq( me, state, bcBuffer, globBC, levelDesc, tree, nSize, &
    &                     iLevel, sim_time, neigh, layout, fieldProp, varPos,  &
    &                     nScalars, varSys, derVarPos, physics, iField,        &
    &                     mixture                                              )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: fEq( layout%fStencil%QQ )
    real(kind=rk) :: rho(1), vel(3,1)
    real(kind=rk) :: vel_b(3*globBC%nElems(iLevel)), inv_vel
    integer :: iELem, iDir, QQ, elemPos, posInBuffer
    integer :: bcVel_pos
    ! ---------------------------------------------------------------------------

    QQ = layout%fStencil%QQ
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_b                               )

    ! convert physical velocity into LB velocity
    vel_b = vel_b * inv_vel

    ! @todo: these if condition can only be done once.
    do iElem = 1, globBC%nElems(iLevel)
      ! Calculate the density of current element
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      rho(1) = sum( bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                  (posInBuffer-1)*nScalars+varPos(QQ) ) )
      vel(:,1) = vel_b((iElem-1)*3+1:iElem*3)

      ! compute equilibrium
      call derVarPos%equilFromMacro( density  = rho,    &
        &                            velocity = vel,    &
        &                            iField   = iField, &
        &                            nElems   = 1,      &
        &                            varSys   = varSys, &
        &                            layout   = layout, &
        &                            res      = fEq     )

      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      do iDir = 1, layout%fStencil%QQN
        ! Write the values
        ! The bitmask points into the incoming direction into the flow domain,
        !  which actually we want to update
        ! * For PUSH, we write to the incoming position,
        !   as the kernel reads it from there without propagation.
        ! * For PULL, we need to write to the inverse direction, as the kernel
        !   performs a bounce back before reading it.
        !   However, this bounced back direction actually comes from the
        !   non-existent boundary element and would point into the incoming
        !   direction, so the value has to be treated and set as if it points
        !   really into the incoming direction.
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state(  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq( iDir )
        end if ! bitMask
      end do ! iDir
    end do ! iElem

  end subroutine velocity_eq
! ****************************************************************************** !

! ****************************************************************************** !
  !> Velocity Non-Equilibrium type boundary conditions from
  !! Guo, Z., & Shi, B. (2002). "Non-equilibrium extrapolation method for
  !! velocity and pressure boundary conditions in the lattice Boltzmann method."
  !! Chinese Physics, (November 2016).
  !!
  !! The incoming densities are updated as
  !! $$f_i(\vec{x}_b,t) = f^{eq}_i(\vec{x}_b,t)
  !!          + (1-\omega)(f_i(\vec{x}_f,t) - f^{eq}_i(\vec{x}_f,t)) $$
  !! macroscopic values for $$f^{eq}_i(\vec{x}_b,t)$$:
  !! - density is extrapolated (0 or 1st order) from the flow domain
  !!   + 0-order: so we simply use the density of the boundary element itself
  !! - velocity is used as defined in the configuration file
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'inlet',
  !!    kind = 'vel_neq',
  !!    velocity = 'vel',
  !!  }
  !!}
  !!variable = {
  !!  name = 'inlet_vel',
  !!  ncomponents = 3,
  !!  vartype = 'st_fun',
  !!  st_fun = {0.06, 0.0, 0.0}
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine vel_neq( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
    &                 iLevel, sim_time, neigh, layout, fieldProp, varPos,   &
    &                 nScalars, varSys, derVarPos, physics, iField, mixture )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: fEq( layout%fStencil%QQ )
    real(kind=rk) :: fEq_Bnd( layout%fStencil%QQ )
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: rho(1), vel(3,1), omega
    real(kind=rk) :: vel_b(3*globBC%nElems(iLevel)), inv_vel
    integer :: iELem, iDir, QQ, elemPos, posInBuffer
    integer :: bcVel_pos
    ! ---------------------------------------------------------------------------

    QQ = layout%fStencil%QQ
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_b                               )

    ! convert physical velocity into LB velocity
    vel_b = vel_b * inv_vel

    ! @todo: these if condition can only be done once.
    do iElem = 1, globBC%nElems(iLevel)
      ! Calculate the density of current element
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                    (posInBuffer-1)*nScalars+varPos(QQ)  )
      rho(1) = sum( fTmp )
      vel(:,1) = vel_b((iElem-1)*3+1:iElem*3)

      ! compute equilibrium on Boundary element
      call derVarPos%equilFromMacro( density  = rho,    &
        &                            velocity = vel,    &
        &                            iField   = iField, &
        &                            nElems   = 1,      &
        &                            varSys   = varSys, &
        &                            layout   = layout, &
        &                            res      = fEq_Bnd )

      ! compute equilibrium on Fluid element
      call derVarPos%equilFromState( state  = fTmp,   &
        &                            iField = iField, &
        &                            nElems = 1,      &
        &                            varSys = varSys, &
        &                            layout = layout, &
        &                            res    = fEq     )

      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! omega is STfun so get omega for current element
      omega = fieldProp%fluid%viscKine%omLvl(iLevel)%val(elemPos)

      do iDir = 1, layout%fStencil%QQN
        ! Write the values
        ! The bitmask points into the incoming direction into the flow domain,
        !  which actually we want to update
        ! * For PUSH, we write to the incoming position,
        !   as the kernel reads it from there without propagation.
        ! * For PULL, we need to write to the inverse direction, as the kernel
        !   performs a bounce back before reading it.
        !   However, this bounced back direction actually comes from the
        !   non-existent boundary element and would point into the incoming
        !   direction, so the value has to be treated and set as if it points
        !   really into the incoming direction.
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state(  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq_Bnd( iDir ) + (1.0_rk-omega)*(fTmp(iDir) - fEq(iDir))
        end if ! bitMask
      end do ! iDir
    end do ! iElem

  end subroutine vel_neq
! ****************************************************************************** !

! ****************************************************************************** !
  !> Pressure Non-Equilibrium type boundary conditions from
  !! Guo, Z., & Shi, B. (2002). "Non-equilibrium extrapolation method for
  !! velocity and pressure boundary conditions in the lattice Boltzmann method."
  !! Chinese Physics, (November 2016).
  !!
  !! The incoming densities are updated as
  !! $$f_i(\vec{x}_b,t) = f^{eq}_i(\vec{x}_b,t)
  !!          + (1-\omega)(f_i(\vec{x}_f,t) - f^{eq}_i(\vec{x}_f,t)) $$
  !! macroscopic values for $$f^{eq}_i(\vec{x}_b,t)$$:
  !! - density is defined in config file as pressure
  !! - velocity is extrapolated from fluid node
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'outlet',
  !!    kind = 'pressure_neq',
  !!    pressure = 'press',
  !!  }
  !!}
  !!variable = {
  !!  name = 'press',
  !!  ncomponents = 1,
  !!  vartype = 'st_fun',
  !!  st_fun = 0.0
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine press_neq( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
    &                   iLevel, sim_time, neigh, layout, fieldProp, varPos,   &
    &                   nScalars, varSys, derVarPos, physics, iField, mixture )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: fEq( layout%fStencil%QQ * globBC%nElems(iLevel) )
    real(kind=rk) :: fEq_Bnd( layout%fStencil%QQ * globBC%nElems(iLevel) )
    real(kind=rk) :: fTmp( layout%fStencil%QQ * globBC%nElems(iLevel) )
    real(kind=rk) :: fTmp_loc( layout%fStencil%QQ )
    real(kind=rk) :: vel(3,globBC%nElems(iLevel)), omega
    real(kind=rk) :: rhoB(globBC%nElems(iLevel)), inv_rho_phy
    integer :: iELem, iDir, QQ, elemPos, posInBuffer, offset
    integer :: bcPress_pos
    ! ---------------------------------------------------------------------------

    QQ = layout%fStencil%QQ
    inv_rho_phy = 1.0_rk / physics%fac(iLevel)%press * cs2inv

    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure variable from spacetime function
    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
      & varSys  = varSys,                               &
      & time    = sim_time,                             &
      & iLevel  = iLevel,                               &
      & idx     = me%bc_states%pressure                 &
      &           %pntIndex%indexLvl(iLevel)            &
      &           %val(1:globBC%nElems(iLevel)),        &
      & nVals   = globBC%nElems(iLevel),                &
      & res     = rhoB                                  )

    ! convert physical pressure to LB density
    rhoB =  rhoB* inv_rho_phy

    do iElem = 1, globBC%nElems(iLevel)
      ! Calculate the density of current element
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp_loc = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                  (posInBuffer-1)*nScalars+varPos(QQ)  )

      fTmp((iElem-1)*QQ+varPos(1):(iElem-1)*QQ+varPos(QQ)) = fTmp_loc

      ! compute velocity on Fluid element
      call derVarPos%velFromState( state  = fTmp_loc,     &
        &                          iField = iField,       &
        &                          nElems = 1,            &
        &                          varSys = varSys,       &
        &                          layout = layout,       &
        &                          res    = vel(:, iElem) )
    end do

    ! compute equilibrium on Boundary element
    call derVarPos%equilFromMacro( density  = rhoB,                  &
      &                            velocity = vel,                   &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq_Bnd                )

    ! compute equilibrium on Fluid element
    call derVarPos%equilFromState( state  = fTmp,                  &
      &                            iField = iField,                &
      &                            nElems = globBC%nElems(iLevel), &
      &                            varSys = varSys,                &
      &                            layout = layout,                &
      &                            res    = fEq                    )

    do iElem = 1, globBC%nElems(iLevel)
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! omega is STfun so get omega for current element
      omega = fieldProp%fluid%viscKine%omLvl(iLevel)%val(elemPos)

      do iDir = 1, layout%fStencil%QQN
        ! Write the values
        ! The bitmask points into the incoming direction into the flow domain,
        !  which actually we want to update
        ! * For PUSH, we write to the incoming position,
        !   as the kernel reads it from there without propagation.
        ! * For PULL, we need to write to the inverse direction, as the kernel
        !   performs a bounce back before reading it.
        !   However, this bounced back direction actually comes from the
        !   non-existent boundary element and would point into the incoming
        !   direction, so the value has to be treated and set as if it points
        !   really into the incoming direction.
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          offset = (iElem-1)*QQ
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state(  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq_Bnd( offset+iDir )                                         &
            & + (1.0_rk-omega)*(fTmp(offset+iDir) - fEq(offset+iDir))
        end if ! bitMask
      end do ! iDir
    end do ! iElem

  end subroutine press_neq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Inlet Velocity Equilibrium type boundary conditions with mass flow rate
  !! as input
  !!
  !! The incoming densities are set to the equilibrium distribution with
  !! macroscopic values:
  !! - density is extrapolated (0 or 1st order) from the flow domain
  !!   + 0-order: so we simply use the density of the boundary element itself
  !! - velocity is used as defined in the configuration file
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'inlet',
  !!    kind = 'mfr_eq',
  !!    massflowrate = 'mfr' }
  !!}
  !!variable = {
  !!  name = 'mfr',
  !!  ncomponents = 1,
  !!  vartype = 'st_fun',
  !!  st_fun = 0.1
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine mfr_eq( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
    &                iLevel, sim_time, neigh, layout, fieldProp, varPos,   &
    &                nScalars, varSys, derVarPos, physics, iField, mixture )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: massFlowRate(globBC%nElems(iLevel))
    real(kind=rk) :: massFlowRateToVel
    real(kind=rk) :: rho(globBC%nElems(iLevel))
    real(kind=rk) :: velocity(3,globBC%nElems(iLevel))
    integer :: iELem, iDir, iLvl, QQ, elemPos
    real(kind=rk) :: fEq( layout%fStencil%QQ*globBC%nElems(iLevel) )
    reaL(kind=rk) :: area
    integer :: posInBuffer
    integer :: bcMfr_pos
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ

    ! position of mass flow rate spacetime function variable in varSys
    bcMfr_pos = me%bc_states%massFlowRate%varPos
    ! get mass flow rate
    call varSys%method%val(bcMfr_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%massFlowRate           &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = massFlowRate                        )

    ! convertion factor to convert mass flow rate to velocity
    ! vel_lattice = massFlowRate/(density*areaOfCrossSection)
    ! area of cross section = sum(nElems(iLevel)*dx(iLevel)**2)
    area = 0.0_rk
    do iLvl = tree%global%minLevel, tree%global%maxLevel
      area = area + globBC%nElems_totalLevel(iLvl) * physics%dxLvl(iLvl)**2
    end do
    massFlowRateToVel = 1.0_rk / ( physics%rho0 * area ) &
      &               / physics%fac( iLevel )%vel

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
!write(*,*) 'velocity_P ', massFlowRateToVel, 1.0/massFlowRateToVel
!write(*,*) 'velocity_L ', massFlowRateToVel

    ! Calculate the density of current element
    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      rho(iElem) = sum(bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                        (posInBuffer-1)*nScalars+varPos(QQ) ) )
    end do

    ! compute velocity from massFlowRate
    do iElem = 1, globBC%nElems(iLevel)
      velocity(:,iElem) = massFlowRate(iElem) * massFlowRateToVel              &
        & * layout%fStencil%cxDirRK(:, globBC%elemLvl(iLevel)%                 &
        &                                                  normalInd%val(iElem))
    end do

    ! compute equilibrium
    call derVarPos%equilFromMacro( density  = rho,                   &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq                    )

    do iElem = 1, globBC%nElems(iLevel)
      ! bc element position in total
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      do iDir = 1, layout%fStencil%QQN
        ! Write the values
        ! The bitmask points into the incoming direction into the flow domain,
        !  which actually we want to update
        ! * For PUSH, we write to the incoming position,
        !   as the kernel reads it from there without propagation.
        ! * For PULL, we need to write to the inverse direction, as the kernel
        !   performs a bounce back before reading it.
        !   However, this bounced back direction actually comes from the
        !   non-existent boundary element and would point into the incoming
        !   direction, so the value has to be treated and set as if it points
        !   really into the incoming direction.
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state(  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq( (iElem-1)*QQ + iDir )
        end if ! bitMask
      end do ! iDir
    end do ! iElem

  end subroutine mfr_eq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Outlet Pressure Equilibrium type boundary conditions
  !!
  !! The incoming densities are set to the equilibrium distribution with
  !! macroscopic values used from:
  !! * Velocity: Linear extrapolation of velocity of precollision of current
  !!   time step
  !! * Pressure as defined in the configuration file
  !!
  !! Usage:\n
  !!```lua
  !!  boundary_condition = {
  !!    { label = 'outlet',
  !!       kind = 'pressure_eq',
  !!       pressure = 'pressure_out',
  !!    }
  !!  }
  !! variable = {
  !!   name = 'pressure_out',
  !!   ncomponents = 1,
  !!   vartype = 'st_fun',
  !!   st_fun = 1.0
  !! }
  !!```
  !!
  !! This is described in the paper:\n
  !! S. Izquierdo, P. Martinez-Lera, and N. Fueyo,\n
  !! Analysis of open boundary effects in unsteady lattice Boltzmann simulations
  !! Comput. Math. Applications, vol. 58, no. 5, pp. 914–921, Sep. 2009.\n
  !!
  !! The velocity is extrapolated (1st order) from the flow domain
  !! See equation (6)
  !!
  !! DISCLAIMER: This BC requires pre-collision PDF from 1st and 2nd neighbor
  !! which are not available if those neighbors are boundary halo on different
  !! process. So, if you encounter a problem with this BC, please check the
  !! domain distribution to verify whether 2 neighbors of boundary elements
  !! are in one process. The domain distribution can be changed by shifting the
  !! boundary by one element or by changing the origin of bounding cube or
  !! by using different number of MPI processes.
  !!
  !! \todo 20200410, KM: Check if this BC works if neighbufferPre is replaced
  !! by neighbufferPost
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [boundary_type:fnct]] function pointer.
  subroutine pressure_eq( me, state, bcBuffer, globBC, levelDesc, tree, nSize, &
    &                     iLevel, sim_time, neigh, layout, fieldProp, varPos,  &
    &                     nScalars, varSys, derVarPos, physics, iField,        &
    &                     mixture                                              )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: fEq( layout%fStencil%QQ*globBC%nElems(iLevel) )
    ! need to compute precollision velocity of neighbor element
    ! Density on boundary element
    real(kind=rk) :: rho( globBC%nElems(iLevel))
    ! Velocity on boundary element
    real(kind=rk) :: velocity(3,globBC%nElems(iLevel))
    ! Vel on frst neighbor element
    real(kind=rk) :: uxB_1(3*globBC%nElems(iLevel))
    ! Vel on scnd neighbor element
    real(kind=rk) :: uxB_2(3*globBC%nElems(iLevel))
    integer :: iDir, iElem, QQ, elemPos
    integer :: bcPress_pos
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ

    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure variable from spacetime function
    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
      & varSys  = varSys,                               &
      & time    = sim_time,                             &
      & iLevel  = iLevel,                               &
      & idx     = me%bc_states%pressure                 &
      &           %pntIndex%indexLvl(iLevel)            &
      &           %val(1:globBC%nElems(iLevel)),        &
      & nVals   = globBC%nElems(iLevel),                &
      & res     = rho                                   )

    ! convert physical pressure to LB density
    rho = rho / physics%fac( iLevel )%press * cs2inv

    ! KM: velocity is taken from precollision of current
    ! Calculate velocity of two neighbors
    call derVarPos%velFromState(                           &
      &     state  = me%neigh(iLevel)%neighBufferPre_nNext(1,:), &
      &     iField = iField,                               &
      &     nElems = globBC%nElems(iLevel),                &
      &     varSys = varSys,                               &
      &     layout = layout,                               &
      &     res    = uxB_1                                 )
    call derVarPos%velFromState(                           &
      &     state  = me%neigh(iLevel)%neighBufferPre_nNext(2,:), &
      &     iField = iField,                               &
      &     nElems = globBC%nElems(iLevel),                &
      &     varSys = varSys,                               &
      &     layout = layout,                               &
      &     res    = uxB_2                                 )

    ! caluclate my velocity by interpolation
    do iElem = 1, globBC%nElems(iLevel)
      velocity(:,iElem) =   1.5_rk * uxB_1((iElem-1)*3+1:iElem*3) &
        &                  - 0.5_rk * uxB_2((iElem-1)*3+1:iElem*3)
    end do

    ! compute equilibrium
    call derVarPos%equilFromMacro( density  = rho,                   &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq                    )

    ! Run over all elements in list for this boundary
    do iElem = 1, globBC%nElems( iLevel )
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      do iDir = 1, layout%fStencil%QQN
        ! the bitmask points into the incoming direction, which must
        ! be updated. As the kernel reads this value with FETCH,
        ! we store it to the FETCH position.
        ! This results for
        ! - PUSH in storing the density to where the bitmask
        !   points, i.e. incoming
        ! - PULL storing the density pointing outwards, but as
        !   this density is being bounced back before collision in
        !   the kernel, this is the correct position
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Now assign the values
state(  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) = &
          & fEq( (iElem-1)*QQ + iDir )
        end if ! bitMask
      end do ! iDir
    end do ! iElem

  end subroutine pressure_eq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Outlet Pressure do-nothing boundary is the open boundary condition for
  !! incompressible model.
  !! This BC sets reference density at boundary so pressure is not loaded
  !! config file.
  !! Here, the normal velocity is extrapolated from 1st fluid node and
  !! tangential velocity is extrapolated from 2nd fluid node in normal
  !! direction.
  !! Algorithm used in this boundary condition:
  !! fEq(1,u) and fEq(rho,u) are computed using macroscopic values from current
  !! element.
  !! In fNeq, post-collision of current time step is used for normal direction.
  !!
  !! This is taken from the paper:
  !! M. Junk and Z. Yang,
  !! Asymptotic Analysis of Lattice Boltzmann Outflow Treatments
  !! Commun. Comput. Phys., pp. 1–11, 2011.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine outlet_dnt( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
    &                    iLevel, sim_time, neigh, layout, fieldProp, varPos,   &
    &                    nScalars, varSys, derVarPos, physics, iField, mixture )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    ! fTmp uses AOS layout
    real(kind=rk) :: fTmp( nScalars*globBC%nElems(iLevel) )
    ! fEq, fEq0 use AOS layout
    real(kind=rk) :: fEq ( layout%fStencil%QQ*globBC%nElems(iLevel) )
    real(kind=rk) :: fEq0( layout%fStencil%QQ*globBC%nElems(iLevel) )
    real(kind=rk) :: fNEq
    real(kind=rk) :: rho( globBC%nElems(iLevel) ) ! my density
    real(kind=rk) :: uxB(3)  ! Velocity on boundary element
    real(kind=rk) :: uxB_n(3), uxB_n2(3), uxB_t(3), normal(3)
    integer :: normalDir
    real(kind=rk) :: velocity(3,globBC%nElems(iLevel)) ! my velocity
    real(kind=rk) :: uxF_2(globBC%nElems(iLevel)* 3) ! first  neigh velocity
    real(kind=rk) :: omega, visc, facFNeq
    integer :: iDir, invDir, iElem, QQ, elemPos, posInBuffer
    ! ---------------------------------------------------------------------------

    QQ    = layout%fStencil%QQ

    do iElem = 1, globBC%nElems(iLevel)

      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp( (iElem-1)*nScalars+1: (iElem-1)*nScalars+QQ ) &
        &       = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                   (posInBuffer-1)*nScalars+varPos(QQ)  )

      rho( iElem ) = sum( fTmp( (iElem-1)*QQ+1:iElem*QQ ) )

      ! Get my velocity
      call derVarPos%velFromState( state  = fTmp( (iElem-1)*nScalars+1:iElem*nScalars ), &
        &                          iField = iField,            &
        &                          nElems = 1,                 &
        &                          varSys = varSys,            &
        &                          layout = layout,            &
        &                          res    = velocity(:,iElem) )
    end do ! iElem


    ! compute my equilibrium using rho
    call derVarPos%equilFromMacro( density  = rho,                   &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq                    )

    ! Calculate my equilibrium distribution using rho0
    rho = 1.0_rk

    ! compute my equilibrium using rho0
    call derVarPos%equilFromMacro( density  = rho,                   &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq0                   )

    ! Get my neighbors velocity
    !KM: velocity is taken from precollision of current
    call derVarPos%velFromState(                            &
      &     state  = me%neigh(iLevel)%neighBufferPost(1,:), &
      &     iField = iField,                                &
      &     nElems = globBC%nElems(iLevel),                 &
      &     varSys = varSys,                                &
      &     layout = layout,                                &
      &     res    = uxF_2                                  )

    ! Run all elements in list for this boundary
    do iElem = 1, globBC%nElems(iLevel)
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )

      ! Extrapolate normal velocity from first neighbor
      ! and tangential velocity from second neighbor as per paper
      normalDir = globBC%elemLvl( iLevel )%normalInd%val( iElem )
      normal = layout%fStencil%cxDirRK(:, normalDir)
      ! Extract normal velocity
      uxB_n = dot_product( velocity(:, iElem), abs(normal) ) &
        &   * abs(normal)
      uxB_n2 = dot_product( uxF_2((iElem-1)*3+1 : iElem*3),   &
        &                   abs(normal)                     ) &
        &    * abs(normal)
      ! Extract tangential velocity
      uxB_t = uxF_2((iElem-1)*3+1 : iElem*3) - uxB_n2
      ! velocity on boundary
      uxB = uxB_n + uxB_t


      ! Treat all directions, but actually we only want to treat
      ! non-normal directions
      ! Equation (3.2b)
      do iDir = 1, layout%fStencil%QQN
        ! the bitmask points into the incoming direction, which must
        ! be updated. As the kernel reads this value with FETCH,
        ! we store it to the FETCH position.
        ! This results for
        ! - PUSH in storing the density to where the bitmask
        !   points, i.e. incoming
        ! - PULL storing the density pointing outwards, but as
        !   this density is being bounced back before collision in
        !   the kernel, this is the correct position
        if ( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          invDir = layout%fStencil%cxDirInv(iDir)
          ! update the incoming velocity direction
          state(  neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0)  = &
         ! This bc is stable only using FETCH with inverse direction of bitmask
         ! For PULL this means, get the outgoing one from upstream neighbor
         ! For PUSH this means, get the outgoing one from local element
         ! KM: Using FETCH gives stable result
         ! @todo test for different testcase for stability
            & state(neigh((invdir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0)&
            & + 6._rk*layout%weight( iDir )                  &
            !extrapolated velocity is used here to remove fluctuations
            & * ( layout%fStencil%cxDirRK( 1, iDir )*uxB(1)  &
            & +   layout%fStencil%cxDirRK( 2, iDir )*uxB(2)  &
            & +   layout%fStencil%cxDirRK( 3, iDir )*uxB(3)  )
        end if ! bitMask
      end do ! iDir

      ! omega is STfun so get omega for current element
      omega = fieldProp%fluid%viscKine%omLvl(iLevel)%val(elemPos)
      visc  = fieldProp%fluid%viscKine%dataOnLvl(iLevel)%val(elemPos)
      facFNeq = visc * omega - 1._rk

      ! then overwrite the normal direction with special treatment
      ! Equation (3.2a)
      iDir = globBC%elemLvl( iLevel )%normalInd%val( iElem )
      invDir = layout%fStencil%cxDirInv(iDir)
      fNEq =  fTmp( (iElem-1)*QQ + invDir ) - fEq( (iElem-1)*QQ + invDir )
      state(  neigh(( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0)  = &
        & fEq0( iDir+(iElem-1)*QQ ) - facFNeq * fNeq
    end do

  end subroutine outlet_dnt
! ****************************************************************************** !


  ! ****************************************************************************** !
  !> Outlet Pressure extrapolation boundary.
  !! qVal for this boundary must be 0.0
  !!
  !! This is taken from the paper:
  !! M. Junk and Z. Yang,
  !! Asymptotic Analysis of Lattice Boltzmann Outflow Treatments
  !! Commun. Comput. Phys., pp. 1-11, 2011.
  !!
  !! This boundary condition prescribes density = 1.0 at outlet by forcing the
  !! equlibirium part of \( f_i \) in normal direction (equation 3.13a):
  !! \[
  !!  f_i(n+1,j_0) = \left[ F^{eq}_{i}(1,u) + f^{neq}_i \right] (n,j_0)
  !! \]
  !! For non-normal direction, \( f_i \) is extrapolated from neighrbors by
  !! (equation 3.13b):
  !! \[
  !!  f_i(n+1,j_0) = 2f_i(n+1,j_0-n) - f_i(n+1,j_0-2n)
  !! \]
  !!
  !! Usage:
  !!----
  !!```lua
  !!  boundary_condition = {
  !!    { label = 'outlet',
  !!       kind = 'pressure_expol',
  !!       pressure = 'pressure_out',
  !!    }
  !!  }
  !!  variable = {
  !!    name = 'pressure_out',
  !!    ncomponents = 1,
  !!    vartype = 'st_fun',
  !!    st_fun = 1.0
  !!  }
  !!```
  !! DISCLAIMER: This BC requires pre-collision PDF from 1st and 2nd neighbor
  !! which are not available if those neighbors are boundary halo on different
  !! process. So, if you encounter a problem with this BC, please check the
  !! domain distribution to verify whether 2 neighbors of boundary elements
  !! are in one process. The domain distribution can be changed by shifting the
  !! boundary by one element or by changing the origin of bounding cube or
  !! by using different number of MPI processes.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine pressure_expol( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                        nSize, iLevel, sim_time, neigh, layout,         &
    &                        fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                        physics, iField, mixture                        )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: fPostCol
    real(kind=rk) :: fTmp_1 ! first  neighbor pdf
    real(kind=rk) :: fTmp_2 ! second neighbor pdf
    real(kind=rk) :: fEq(  layout%fStencil%QQ * globBC%nElems(iLevel) )
    real(kind=rk) :: fEq0( layout%fStencil%QQ * globBC%nElems(iLevel) )
    real(kind=rk) :: rho( globBC%nElems(iLevel) ), inv_rho_phy
    real(kind=rk) :: velocity( 3, globBC%nElems(iLevel) )
    integer :: iDir, iElem, QQ, elemPos, invDir
    integer :: iLink, bcPress_pos
    logical :: axisNormal
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos, vel_pos(3), elemOff, posInBuffer
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    dens_pos = varSys%method%val(derVarPos%density)%auxField_varPos(1)
    vel_Pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)
!write(dbgUnit(1),*) 'bclabel ', trim(me%label)
!write(dbgUnit(1),*) 'sim time ', sim_time%iter
    ! Here starts the copied part of pressure_expol
    QQ = layout%fStencil%QQ
    inv_rho_phy = 1.0_rk / physics%fac(iLevel)%press * cs2inv

    ! Get density and velocity from previous time step from auxField
    associate ( auxField => fPtr%solverData%scheme%auxField(iLevel)%val )
      do iElem = 1, globBC%nElems( iLevel )
        elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
        elemOff = (elemPos-1)*varSys%nAuxScalars
        rho( iElem ) = auxField(elemOff + dens_pos)
        velocity(1, iElem) = auxField(elemOff + vel_pos(1))
        velocity(2, iElem) = auxField(elemOff + vel_pos(2))
        velocity(3, iElem) = auxField(elemOff + vel_pos(3))
      end do
    end associate

    ! compute my equilibrium (fEq) using calculated rho
    call derVarPos%equilFromMacro( density  = rho,                   &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq                    )

    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure variable from spacetime function
    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
      & varSys  = varSys,                               &
      & time    = sim_time,                             &
      & iLevel  = iLevel,                               &
      & idx     = me%bc_states%pressure                 &
      &           %pntIndex%indexLvl(iLevel)            &
      &           %val(1:globBC%nElems(iLevel)),        &
      & nVals   = globBC%nElems(iLevel),                &
      & res     = rho                                   )

    ! convert physical pressure to LB density
    rho = rho * inv_rho_phy

    ! compute my equilibrium (fEq0)using user defined rho
    call derVarPos%equilFromMacro( density  = rho,                   &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq0                   )
    ! Here ends the copied part of pressure_expol

    ! linkwise treatment for links that have boundary
    do iLink = 1, me%links(iLevel)%nVals
      iElem = me%outletExpol(iLevel)%iElem(iLink)

      ! if any of the neighbor has property qVal then fall back to simple
      ! equilibrium boundary condition
      ! elements that has property qVal
      ! @todo: create a list for these special links.
      if( btest( levelDesc%property(                        &
        &        me%neigh(iLevel)%posInState(1, iElem)), prp_hasQVal )   &
        & .or. btest( levelDesc%property(                   &
        &        me%neigh(iLevel)%posInState(2, iElem)), prp_hasQVal ) ) &
        & then
        fTmp_1 = 0.0_rk
        fTmp_2 = 0.0_rk
        state( me%links(iLevel)%val(iLink) ) &
          &  = fEq0( me%outletExpol(iLevel)%statePos(iLink) )
!write(dbgUnit(1),*) 'fEq0 ', fEq0( me%outletExpol(iLevel)%statePos(iLink) )
      else
        fTmp_1 = me%neigh(iLevel)%neighBufferPre_nNext(1,                  &
          &                         me%outletExpol(iLevel)%statePos(iLink) )
        fTmp_2 = me%neigh(iLevel)%neighBufferPre_nNext(2,                  &
          &                         me%outletExpol(iLevel)%statePos(iLink) )
        state( me%links(iLevel)%val(iLink) ) &
          & = 1.5_rk*fTmp_1 - 0.5_rk*fTmp_2

! debug output ------------------------------------------------
! debug output ------------------------------------------------

      end if
    end do ! iLink, linkwise treatment ends here

    ! then overwrite the normal direction with special treatment
    ! Equation (3.13a)
    do iElem = 1, globBC%nElems(iLevel)
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      iDir   = globBC%elemLvl( iLevel )%normalInd%val( iElem )

      axisNormal = ( ( abs(layout%fStencil%cxDir( 1, iDir ))        &
        &            + abs(layout%fStencil%cxDir( 2, iDir ))        &
        &            + abs(layout%fStencil%cxDir( 3, iDir )) ) == 1 )

      if (axisNormal) then
        invDir = layout%fStencil%cxDirInv( iDir )

        ! get f_non_eq from postcollision
        posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
        fPostCol = bcBuffer( (posInBuffer-1)*nScalars+varPos(invDir) )

        state(  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
          & )  = &
          ! KM:we use fEq0 since this equilibrium has been computed from
          ! defined density at boundary
          & fEq0((iElem-1)*QQ+iDir) + (fPostCol                 &
          &                       -  fEq( (iElem-1)*QQ+invDir ))
      end if
    end do ! iElem

! debug output ------------------------------------------------
! debug output ------------------------------------------------

  end subroutine pressure_expol
  ! ****************************************************************************** !


  ! ****************************************************************************** !
  !> Inlet Velocity Bounce Back boundary condition with qvalues
  !!
  !! This is taken from the paper:
  !! [1] S. Izquierdo and N. Fueyo, "Characteristic  non-reflecting
  !! boundary  conditions  for  open  boundaries  in  lattice
  !!  Boltzmann  methods," Physical Review E, vol. 78, no. 46707, 2008
  !! The incoming densities \( \bar{f_{alpha}} \) are reconstructed from the
  !! already reflected densities, enforcing the given velocity me\%ubb\%velocity.
  !! [2] M. Junk and Z. Yang. One-point boundary condition for the
  !! lattice Boltzmann method, Physical Review E, vol 72, issue 8, year 2005.
  !!
  !! This boundary condition has error of 1st order pressure and 2nd order velocity
  !! only if wall is located exactly at q=1/2 else the accuracy of both pressure
  !! and velocity reduce by order 1.
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'inlet',
  !!    kind = 'velocity_bounceback',
  !!    velocity = 'inlet_vel',
  !!  }
  !!}
  !!variable = {
  !!  name = 'inlet_vel',
  !!  ncomponents = 3,
  !!  vartype = 'st_fun',
  !!  st_fun = {0.06, 0.0, 0.0}
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine velocity_bounceback_incomp( me, state, bcBuffer, globBC,         &
    &                                    levelDesc, tree, nSize, iLevel,      &
    &                                    sim_time, neigh, layout, fieldProp,  &
    &                                    varPos, nScalars, varSys, derVarPos, &
    &                                    physics, iField, mixture             )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    integer       :: iLink, iDir
    real(kind=rk) :: fOut, eqPlus, inv_vel
    real(kind=rk) :: vel_b(me%links(iLevel)%nVals*3)
    integer :: bcVel_pos, offset
    ! ---------------------------------------------------------------------------

    inv_vel = 1.0_rk / physics%fac( iLevel )%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Multiply transient velocity with spatial velocity. Spatial term includes
    ! maxvalue
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = me%links(iLevel)%nVals,             &
      & res     = vel_b                               )

    !NEC$ ivdep
    !DIR$ ivdep
    do iLink = 1, me%links(iLevel)%nVals
      offset = (iLink-1)*3

      fOut = bcBuffer( me%inletUbbQVal(iLevel)%outPos(iLink) )

      ! According to the bounce back rule,
      ! Subtract the equlibrium if its computed from alpha-(inverse of
      ! bitmask) or add the equlibrium if its computed in the direction of
      ! bitmask.
      ! In the Izquiredo paper, this is rho0, the reference density
      iDir   = me%inletUbbQVal(iLevel)%iDir(iLink)
      eqPlus = layout%weight(iDir) * 6.0_rk * rho0  &
        &    * (   layout%fStencil%cxDirRK(1, iDir) * vel_b(offset+1) * inv_vel &
        &        + layout%fStencil%cxDirRK(2, iDir) * vel_b(offset+2) * inv_vel &
        &        + layout%fStencil%cxDirRK(3, iDir) * vel_b(offset+3) * inv_vel )

      state( me%links(iLevel)%val(iLink) ) = fOut + eqPlus

    end do ! iLink

  end subroutine velocity_bounceback_incomp
! ****************************************************************************** !


  ! ****************************************************************************** !
  !> Inlet Velocity Bounce Back boundary condition with qvalues for compressible
  !! flows. It is similar to velocity_bounceback except the density is extrapolated
  !! from fluid element.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine velocity_bounceback( me, state, bcBuffer, globBC, levelDesc,      &
    &                             tree, nSize, iLevel, sim_time, neigh,        &
    &                             layout, fieldProp, varPos, nScalars, varSys, &
    &                             derVarPos, physics, iField, mixture          )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    integer       :: iLink, iDir, QQ
    real(kind=rk) :: fOut, eqPlus, inv_vel, rho
    real(kind=rk) :: vel_b(me%links(iLevel)%nVals*3)
    integer :: bcVel_pos, offset, posInBuffer
    ! ---------------------------------------------------------------------------
    QQ = layout%fstencil%QQ

    inv_vel = 1.0_rk / physics%fac( iLevel )%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Multiply transient velocity with spatial velocity. Spatial term includes
    ! maxvalue
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = me%links(iLevel)%nVals,             &
      & res     = vel_b                               )

    !NEC$ ivdep
    !DIR$ ivdep
    do iLink = 1, me%links(iLevel)%nVals
      offset = (iLink-1)*3

      fOut = bcBuffer( me%inletUbbQVal(iLevel)%outPos(iLink) )

      ! calculate density
      posInBuffer  = me%inletUbbQVal(iLevel)%posInBuffer(iLink)
      rho = sum(bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                 (posInBuffer-1)*nScalars+varPos(QQ) ))

      ! According to the bounce back rule,
      ! Subtract the equlibrium if its computed from alpha-(inverse of
      ! bitmask) or add the equlibrium if its computed in the direction of
      ! bitmask.
      ! In the Izquiredo paper, this is rho0, the reference density
      iDir   = me%inletUbbQVal(iLevel)%iDir(iLink)
      eqPlus = layout%weight(iDir) * 6.0_rk * rho  &
        &    * (   layout%fStencil%cxDirRK(1, iDir) * vel_b(offset+1) * inv_vel &
        &        + layout%fStencil%cxDirRK(2, iDir) * vel_b(offset+2) * inv_vel &
        &        + layout%fStencil%cxDirRK(3, iDir) * vel_b(offset+3) * inv_vel )

      state( me%links(iLevel)%val(iLink) ) = fOut + eqPlus

    end do ! iLink

  end subroutine velocity_bounceback
! ****************************************************************************** !


  ! ************************************************************************ !
  !> Inlet Velocity BFL rule  boundary condition
  !!
  !! This is taken from the paper:
  !! "One-point boundary condition for the lattice Boltzmann method", by
  !! M. Junk and Z. Yang published in Physical Review E 72 (2005) and leads to
  !! the chapter 'B. BFL rule'.
  !! This boundary condition shall be an improvement of the already existing
  !! velocity bounce back condition. It introduces a term called
  !! \( f_{i^{\star}}^{b} \).
  !! This one represents the particle distribution:
  !! \[ \hat{f}_{i}(n+1,\mathbf{j}) = \hat{f}_{i^{\star}}^{c}(n,\mathbf{j})
  !! + 6hf_{i}^{\star}\mathbf{c}_{i}
  !! \cdot \mathbf{\Phi}(t_{n},\mathbf{x}_{\mathbf{j}i}) \]
  !!
  !! This boundary condition has error of 1st order velocity and less 1st order
  !! pressure tested for wall is located exactly at q=1/2.
  !!
  !! Like for velocity_bounceback_qval we compute the values at the boundary
  !! points.
  !!
  !! ---------------------------------------------------------------------------
  !! Main differences to velocity_bounceback:
  !!  - implemented "link-wise", only incoming links at the boundary are updated
  !!    for each time step
  !!  - uses q-Values in pre-processing (see mus_bc_header_module)
  !!  - computes state values at the boundary nodes instead of barycentre
  !!    position
  !!  - can be an improvement of the bounce back rule for q-values ~= 0.5
  !! ---------------------------------------------------------------------------
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine velocity_bfl( me, state, bcBuffer, globBC, levelDesc, tree,      &
    &                      nSize, iLevel, sim_time, neigh, layout, fieldProp, &
    &                      varPos, nScalars, varSys, derVarPos, physics,      &
    &                      iField, mixture                                    )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    integer       :: iLink, iDir, QQ, posInBuffer
    real(kind=rk) :: weight
    real(kind=rk) :: fIn, fOut, fNgh, eqPlus, rho
    real(kind=rk) :: cIn, cOut, cNgh, cVel
    real(kind=rk) :: vel_b(me%links(iLevel)%nVals*3), inv_vel
    integer :: bcVel_pos, offset
    ! ---------------------------------------------------------------------------
    QQ = layout%fstencil%QQ

    inv_vel = 1.0_rk / physics%fac( iLevel )%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Multiply transient velocity with spatial velocity. Spatial term includes
    ! maxvalue
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = me%links(iLevel)%nVals,             &
      & res     = vel_b                               )

    do iLink = 1, me%links(iLevel)%nVals
      offset = (iLink-1)*3

      ! QValues and coefficients
      cIn  = me%inletBfl(iLevel)% cIn(iLink)
      cOut = me%inletBfl(iLevel)%cOut(iLink)
      cNgh = me%inletBfl(iLevel)%cNgh(iLink)
      cVel = me%inletBfl(iLevel)%cVel(iLink)

      ! positions of element
      fIn  = bcBuffer(me%inletBfl(iLevel)% inPos(iLink))
      fOut = bcBuffer(me%inletBfl(iLevel)%outPos(iLink))
      fNgh = me%neigh(iLevel)%computeNeighBuf(me%inletbfl(iLevel)%nghPos(iLink))

      ! local density
      posInBuffer  = me%inletBfl(iLevel)%posInBuffer(iLink)
      rho = sum(bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                 (posInBuffer-1)*nScalars+varPos(QQ) ))

      iDir   = me%inletBfl(iLevel)%iDir(iLink)
      weight = layout%weight(iDir)
      eqPlus = weight * 6._rk * rho  &
        &    * (  layout%fStencil%cxDirRK(1, iDir) * vel_b(offset+1) * inv_vel &
        &       + layout%fStencil%cxDirRK(2, iDir) * vel_b(offset+2) * inv_vel &
        &       + layout%fStencil%cxDirRK(3, iDir) * vel_b(offset+3) * inv_vel )

      ! We need to get post-collision pdf in direction
      ! alpha-, which is the inverse direction of bitmask
      ! For PULL this means, get the outgoing one, as this is the one which will be bounced back
      ! For PUSH this means, get the already bounced back pdf back, so take the incoming

      state( me%links(iLevel)%val(iLink) ) &
        &  = cIn*fIn + cOut*fOut + cNgh*fNgh + cVel*eqPlus

    end do ! iLink

  end subroutine velocity_bfl
! ****************************************************************************** !

! ****************************************************************************** !
  !> No comment yet!
  !!
  !! @TODO add comment
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine velocity_bfl_incomp( me, state, bcBuffer, globBC, levelDesc,      &
    &                             tree, nSize, iLevel, sim_time, neigh,        &
    &                             layout, fieldProp, varPos, nScalars, varSys, &
    &                             derVarPos, physics, iField, mixture          )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    integer       :: iLink, iDir
    real(kind=rk) :: weight
    real(kind=rk) :: fIn, fOut, fNgh, eqPlus
    real(kind=rk) :: cIn, cOut, cNgh, cVel
    real(kind=rk) :: vel_b(me%links(iLevel)%nVals*3), inv_vel
    integer :: bcVel_pos, offset
    ! ---------------------------------------------------------------------------

    inv_vel = 1.0_rk / physics%fac( iLevel )%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Multiply transient velocity with spatial velocity. Spatial term includes
    ! maxvalue
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = me%links(iLevel)%nVals,             &
      & res     = vel_b                               )

    do iLink = 1, me%links(iLevel)%nVals
      offset = (iLink-1)*3

      ! QValues and coefficients
      cIn  = me%inletBfl(iLevel)% cIn(iLink)
      cOut = me%inletBfl(iLevel)%cOut(iLink)
      cNgh = me%inletBfl(iLevel)%cNgh(iLink)
      cVel = me%inletBfl(iLevel)%cVel(iLink)

      ! positions of element
      fIn  = bcBuffer(me%inletBfl(iLevel)% inPos(iLink))
      fOut = bcBuffer(me%inletBfl(iLevel)%outPos(iLink))
      fNgh = me%neigh(iLevel)%computeNeighBuf(me%inletbfl(iLevel)%nghPos(iLink))

      iDir   = me%inletBfl(iLevel)%iDir(iLink)
      weight = layout%weight(iDir)
      eqPlus = weight * 6._rk * rho0  &
        &    * (  layout%fStencil%cxDirRK(1, iDir) * vel_b(offset+1) * inv_vel &
        &       + layout%fStencil%cxDirRK(2, iDir) * vel_b(offset+2) * inv_vel &
        &       + layout%fStencil%cxDirRK(3, iDir) * vel_b(offset+3) * inv_vel )

      ! We need to get post-collision pdf in direction
      ! alpha-, which is the inverse direction of bitmask
      ! For PULL this means, get the outgoing one, as this is the one which will be bounced back
      ! For PUSH this means, get the already bounced back pdf back, so take the incoming

      state( me%links(iLevel)%val(iLink) ) &
        &  = cIn*fIn + cOut*fOut + cNgh*fNgh + cVel*eqPlus

    end do ! iLink

  end subroutine velocity_bfl_incomp
! ****************************************************************************** !


  ! **************************************************************************** !
  !> Inlet Velocity Bounce Back boundary condition with mass flow rate as input
  !!
  !! This is taken from the paper:
  !! [1] S. Izquierdo and N. Fueyo, “Characteristic  non-reflecting
  !! boundary  conditions  for  open  boundaries  in  lattice
  !!  Boltzmann  methods,” Physical Review E, vol. 78, no. 46707, 2008
  !! The incoming densities \( \bar{f_{alpha}} \) are reconstructed from the
  !! already reflected densities, enforcing the given velocity
  !! me\%ubb\%velocity.
  !! [2] M. Junk and Z. Yang. One-point boundary condition for the
  !! lattice Boltzmann method,
  !! Physical Review E, vol 72, issue 8, year 2005.
  !!
  !! This boundary condition has error of 1st order pressure and 2nd order
  !! velocity only if wall is located exactly at q=1/2 else the accuracy of
  !! both pressure and velocity reduce by order 1.
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'inlet',
  !!    kind = 'mfr_bounceback',
  !!    massflowrate = 'mfr }
  !!}
  !! variable = {
  !!   name = 'mfr',
  !!   ncomponents = 1,
  !!   vartype = 'st_fun',
  !!   st_fun = 0.1
  !! }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine mfr_bounceback( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                        nSize, iLevel, sim_time, neigh, layout,         &
    &                        fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                        physics, iField, mixture                        )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: eqPlus
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: massFlowRate(globBC%nElems(iLevel))
    real(kind=rk) :: velocity(3), massFlowRateToVel
    integer :: iELem, iDir, iLvl, QQ, posInBuffer
    reaL(kind=rk) :: area
    integer :: bcMfr_pos
    ! ---------------------------------------------------------------------------

    QQ = layout%fStencil%QQ

    ! position of mass flow rate spacetime function variable in varSys
    bcMfr_pos = me%bc_states%massFlowRate%varPos
    ! get mass flow rate
    call varSys%method%val(bcMfr_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%massFlowRate           &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = massFlowRate                        )

    ! convertion factor to convert mass flow rate to velocity
    ! vel_lattice = massFlowRate/(density*nElems*dx^2)
    ! area of cross section = sum(nElems(iLevel)*dx(iLevel)**2)
    ! @todo: element area can be computed once and saved
    area = 0.0_rk
    do iLvl = tree%global%minLevel, tree%global%maxLevel
      area = area + globBC%nElems_totalLevel(iLvl) * physics%dxLvl(iLvl)**2
    end do

    massFlowRateToVel = 1.0_rk / ( physics%rho0 * area )

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
!write(*,*) 'velocity_P ', massFlowRateToVel, 1.0/massFlowRateToVel
    massFlowRateToVel = massFlowRateToVel / physics%fac( iLevel )%vel
!write(*,*) 'velocity_L ', massFlowRateToVel

    do iElem = 1, globBC%nElems( iLevel )
      ! compute velocity from massFlowRate
      velocity = massFlowRate(iElem) * massFlowRateToVel                     &
        & * layout%fStencil%cxDir(:, globBC%elemLvl(iLevel)%normalInd%val(iElem))

      !store state in temp array
      !KM:produce no checker board effect if we use pre-collision values of
      !current time step
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                    (posInBuffer-1)*nScalars+varPos(QQ)  )

      do iDir = 1, layout%fStencil%QQN
        ! Write the values
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          ! Accourding to the bounce back rule,
          ! Subtract the equlibrium if its computed from alpha-(inverse of
          ! bitmask) or add the equlibrium if its computed in the direction of
          ! bitmask.
          ! In the Izquiredo paper, this is rho0, the reference density
          eqPlus =    layout%weight( iDir ) * 6._rk * rho0             &
            &       * ( layout%fStencil%cxDir(1,iDir) * velocity(1)   &
            &       +   layout%fStencil%cxDir(2,iDir) * velocity(2)   &
            &       +   layout%fStencil%cxDir(3,iDir) * velocity(3)   )
 state( neigh((idir-1)*nsize+ globbc%elemlvl(ilevel)%elem%val(ielem))+( ifield-1)* qq+ nscalars*0 ) = &
          ! We need to get post-collision pdf in direction
          ! alpha-, which is the inverse direction of bitmask
          ! For PULL this means, get the outgoing one, as this is the one which will be bounced back
          ! For PUSH this means, get the already bounced back pdf back, so take the incoming
            & fTmp( layout%fStencil%cxDirInv( iDir ) ) + eqPlus
        end if ! bitMask
      end do ! iDir
    end do

  end subroutine mfr_bounceback
! ****************************************************************************** !


  ! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Outlet boundary conditions with zero pressure gradient.
  !!
  !! These boundary conditions use the neighbor cells density and velocity in the
  !! equilibrium function. It is not necessary to specify density at boundary in
  !! the lua configuration file \( f = f^{eq}(\rho_{b-1},u_{b-1}) \)
  !!
  !! MH: Kannan, please update them
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine outlet_zero_prsgrd( me, state, bcBuffer, globBC, levelDesc, tree, &
    &                            nSize, iLevel, sim_time, neigh, layout,       &
    &                            fieldProp, varPos, nScalars, varSys,          &
    &                            derVarPos, physics, iField, mixture           )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    ! defining local variables
    real(kind=rk) :: rho, ux(globBC%nElems(iLevel),3), usq1_3
    integer :: iElem, iDir, QQ
    real(kind=rk) :: fTmp( nScalars * globBC%nElems(iLevel) )
    integer :: posInBuffer
    ! ---------------------------------------------------------------------------

    QQ   = layout%fStencil%QQ

    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp( (iElem-1)*nScalars+1: (iElem-1)*nScalars+QQ ) &
        &       = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                   (posInBuffer-1)*nScalars+varPos(QQ)  )

      ! Get my velocity
      call derVarPos%velFromState( state  = fTmp( (iElem-1)*nScalars+1 &
        &                                        :iElem*nScalars ),    &
        &                          iField = iField,                    &
        &                          nElems = 1,                         &
        &                          varSys = varSys,                    &
        &                          layout = layout,                    &
        &                          res    = ux(iElem, :)               )
    end do

    do iElem = 1, globBC%nElems(iLevel)
      ! Calculate local density and velocity moments
      rho = sum( fTmp( (iElem-1)*QQ+1:iElem*QQ ) )

      usq1_3 =  sum( ux(iElem,:) ) * div1_3

      do iDir = 1, layout%fStencil%QQN
        ! Write the values
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          state(                                                      &
&neigh((idir-1)*nsize+globbc%elemlvl(ilevel)%elem%val(ielem))+(ifield-1)*qq+nscalars*0 &
&         )  =                                                   &
            &  layout%weight( iDir )*( 2._rk*rho + 9._rk*rho0* ( (    &
            &            layout%fStencil%cxDirRK(1,iDir)*ux(iElem,1)  &
            &         +  layout%fStencil%cxDirRK(2,iDir)*ux(iElem,2)  &
            &         +  layout%fStencil%cxDirRK(3,iDir)*ux(iElem,3)  &
            &         )**2                                                 &
            &         -  usq1_3 ) )                                        &
          ! We need to get post-collision pdf in direction
          ! alpha-, which is the inverse direction of bitmask
          ! For PULL this means, get the outgoing one, as this is the one which will be bounced back
          ! For PUSH this means, get the already bounced back pdf back, so take the incoming
            &         - fTmp( (iElem-1)*QQ + layout%fStencil%cxDirInv(iDir) )
        end if
      end do
      ! end if
    end do

  end subroutine outlet_zero_prsgrd
! ****************************************************************************** !


  ! ************************************************************************ !
  !> author: Manuel Hasert
  !! Outlet Pressure Bounce Back boundary condition
  !!
  !! An anti-bounce back is performed for the outgoing links, so that the
  !! incoming densities are set with the defined reference density while
  !! maintaining the correct non-equilibrium part.
  !! The velocity is extrapolated from two neighboring fluid elements (boundary
  !! element itself + one fluid neighbor). The position of the wall is hence
  !! between the boundary element and the (non-existing) solid element
  !! Based on the Bouzidi Boundary conditions, the current implementation is
  !! taken from the paper:
  !! S. Izquierdo and N. Fueyo, "Characteristic  non-reflecing  boundary
  !!  conditions  for  open  boundaries  in  lattice  Boltzmann
  !!  methods," Physical Review E, vol. 78, no. 46707, 2008.
  !!
  !! MH: I deactivated the correction term as it involves information from the
  !! old time step. In parallel, this can not be provided.
  !!
  !! DISCLAIMER: This BC requires pre-collision PDF from 1st and 2nd neighbor
  !! which are not available if those neighbors are boundary halo on different
  !! process. So, if you encounter a problem with this BC, please check the
  !! domain distribution to verify whether 2 neighbors of boundary elements
  !! are in one process. The domain distribution can be changed by shifting the
  !! boundary by one element or by changing the origin of bounding cube or
  !! by using different number of MPI processes.
  !!
  !! \todo 20200410, KM: Check if this BC works if neighbufferPre is replaced
  !! by neighbufferPost
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct] function pointer.
  subroutine pressure_antiBounceBack( me, state, bcBuffer, globBC, levelDesc, &
    &                                 tree, nSize, iLevel, sim_time, neigh,   &
    &                                 layout, fieldProp, varPos, nScalars,    &
    &                                 varSys, derVarPos, physics, iField,     &
    &                                 mixture                                 )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rhoDef(globBC%nElems(iLevel)) ! Density on boundary element
    real(kind=rk) :: rhoF    ! density on fluid element
    real(kind=rk) :: uxF(globBC%nElems(iLevel)*3)  ! Velocity on fluid element
    real(kind=rk) :: uxN(globBC%nElems(iLevel)*3)  ! Velocity on neighbor
    real(kind=rk) :: uxB(3)  ! Velocity on boundary surface
    real(kind=rk) :: usqB, usqF ! Squared velocity magnitude
    real(kind=rk) :: fPlusFluid ! f Plus for anti-bounceback correction
    ! Equilibrium Plus for anti-bounceback correction
    real(kind=rk) :: fEqPlus
    ! Equlibrium Plus from the old time step
    real(kind=rk) :: fEqPlusFluid
    real(kind=rk) :: inv_rho_phy
    real(kind=rk) :: fTmp( nScalars * globBC%nElems(iLevel) )
    real(kind=rk) :: omega       ! Relaxation parameter
    integer :: iDir, iElem, QQ, invDir, elemPos
    integer :: posInBuffer, bcPress_pos
    ! ---------------------------------------------------------------------------

    ! write(dbgUnit(6), *) ''
    ! write(dbgUnit(6), *) ' do pressure_antiBounceBack for bc label: '//trim(globBC%label)

    QQ = layout%fStencil%QQ
    inv_rho_phy = 1.0_rk / physics%fac(iLevel)%press * cs2inv

    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure variable from spacetime function
    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
      & varSys  = varSys,                               &
      & time    = sim_time,                             &
      & iLevel  = iLevel,                               &
      & idx     = me%bc_states%pressure                 &
      &           %pntIndex%indexLvl(iLevel)            &
      &           %val(1:globBC%nElems(iLevel)),        &
      & nVals   = globBC%nElems(iLevel),                &
      & res     = rhoDef                                )

    ! convert physical pressure into LB density
    rhoDef = rhoDef * inv_rho_phy

    do iElem = 1, globBC%nElems( iLevel )
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp( (iElem-1)*nScalars+1: (iElem-1)*nScalars+QQ ) &
        &       = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                   (posInBuffer-1)*nScalars+varPos(QQ)  )
    end do

    ! Get my velocity
    call derVarPos%velFromState( state  = fTmp,                  &
      &                          iField = iField,                &
      &                          nElems = globBC%nElems(iLevel), &
      &                          varSys = varSys,                &
      &                          layout = layout,                &
      &                          res    = uxF                    )

    ! Calculate my neighbor velocity, here Post or Pre are the same for u,rho
    call derVarPos%velFromState( state  = me%neigh(iLevel)%neighBufferPost(1,:), &
      &                          iField = iField,                &
      &                          nElems = globBC%nElems(iLevel), &
      &                          varSys = varSys,                &
      &                          layout = layout,                &
      &                          res    = uxN                    )

    do iElem = 1, globBC%nElems( iLevel )

      rhoF = sum(fTmp( (iElem-1)*nScalars+1: (iElem-1)*nScalars+QQ ))

      uxB =   1.5_rk * uxF( (iElem-1)*3+1 : iElem*3 ) &
        &   - 0.5_rk * uxN( (iElem-1)*3+1 : iElem*3 )

      usqB = uxB(1)*uxB(1) + uxB(2)*uxB(2) + uxB(3)*uxB(3)
      usqF =   uxF((iElem-1)*3+1) * uxF((iElem-1)*3+1) &
        &    + uxF((iElem-1)*3+2) * uxF((iElem-1)*3+2) &
        &    + uxF((iElem-1)*3+3) * uxF((iElem-1)*3+3)

! -----------  debug output    --------------------------------------------
! write(dbgUnit(6), *) iElem, ' treeID: ',&
!   &     levelDesc%total(globBC%elemLvl(iLevel)%elem%val(iElem)), &
!   &     'elemBuffer: ', globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
! write(dbgUnit(6), *) 'rhoDef: ', rhoDef( iElem )
! write(dbgUnit(6), *) '  rhoF: ', rhoF
! write(dbgUnit(6), *) '   uxF: ', uxF( (iElem-1)*3+1:(iElem-1)*3+3)
! write(dbgUnit(6), *) '   uxN: ', uxN( (iElem-1)*3+1:(iElem-1)*3+3)
! write(dbgUnit(6), *) '   uxB: ', uxB(1:3)
! write(dbgUnit(6), *) '  usqB: ', usqB
! write(dbgUnit(6), *) '  usqF: ', usqF
! write(dbgUnit(6), *) ' origin pdf:'
! do iDir = 1, QQ
!   write( dbgUnit(6), *) iDir, fTmp( (iElem-1)*nScalars + iDir )
! end do
! -------------------------------------------------------------------------
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! omega is STfun so get omega for current element
      omega = fieldProp%fluid%viscKine%omLvl(iLevel)%val(elemPos)

      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          ! Calc the correction term and the Dirichlet (Eq) term

          invDir = layout%fStencil%cxDirInv(iDir)

          fEqPlusFluid = layout%weight( iDir )*rhoF                  &
            &          + 4.5_rk*layout%weight( iDir )*rho0*((        &
            &            layout%fStencil%cxDirRK( 1, invDir )*uxF((iElem-1)*3+1) &
            &          + layout%fStencil%cxDirRK( 2, invDir )*uxF((iElem-1)*3+2) &
            &          + layout%fStencil%cxDirRK( 3, invDir )*uxF((iElem-1)*3+3) &
            &          )**2                                       &
            &          - div1_3*usqF)

          fEqPlus = layout%weight( iDir )*rhoDef(iElem)        &
            &     + 4.5_rk*layout%weight( iDir )*rho0*((       &
            &       layout%fStencil%cxDirRK( 1, InvDir)*uxB(1) &
            &     + layout%fStencil%cxDirRK( 2, InvDir)*uxB(2) &
            &     + layout%fStencil%cxDirRK( 3, InvDir)*uxB(3) &
            &     )**2                                         &
            &     - div1_3*usqB)

          fPlusFluid = 0.5_rk*(   fTmp( (iElem-1)*nScalars + iDir )  &
            &                   + fTmp( (iElem-1)*nScalars + invDir ) )

          ! Now assign the values
          ! Actually, we should take fPlusLAst and fEqPlusLast.
          ! But this introduces errors in parallel
          state( neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) = &
          ! antibounceback term
! & - state(  )    &
          ! We need to get post-collision pdf in direction
          ! alpha-, which is the inverse direction of bitmask
          ! For PULL this means, get the outgoing one, as this is the one which will be bounced back
          ! For PUSH this means, get the already bounced back pdf back, so take the incoming
    &   - fTmp( (iElem-1)*nScalars + invDir ) &
    &   + 2._rk*fEqPlus     &                          ! Dirichlet pressure term
    &   + (2._rk - omega )*( fPlusFluid - fEqPlusFluid ) ! Correction term

! -----------   DEBUG output     -------------------------------------------
! write(dbgUnit(6), *) '      fEqPlus: ', fEqPlus
! write(dbgUnit(6), *) ' fEqPlusFluid: ', fEqPlusFluid
! write(dbgUnit(6), *) '   fPlusFluid: ', fPlusFluid
! write(dbgUnit(6), *) '  updated pdf:'
! write( dbgUnit(6), *) 'iDir', iDir, 'invDir', invDir, state(                      &
! & neigh (( idir-1)* nsize+ globbc%elemlvl(ilevel)%elem%val( ielem ))+( ifield-1)* qq+ nscalars*0)
! --------------------------------------------------------------------------

        end if ! bitMask
      end do ! iDir
    end do ! iElem

  end subroutine pressure_antiBounceBack
! ****************************************************************************** !

! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Characteristic-based non-reflective inlet boundary conditions for
  !! incompressible flows
  !! @todo: add explaination of these steps
  !!
  !! These boundary conditions are taken from the paper:
  !!   S. Izquierdo and N. Fueyo,
  !!   "Characteristic nonreflecting boundary
  !!   conditions for open boundaries in lattice Boltzmann methods,"
  !!   Physical Review E, vol. 78, no. 46707, 2008.
  !!
  !! > Note: unstable behavior for high omega relaxation parameters.
  !! >       It is recommended to use a sponge layer with a
  !! >       \ref mus_aux_module::mus_setspatialomega "spatial omega"
  !! >       where the omega is decreased towards the boundaries to increase
  !! >       stability
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine inlet_nrbc_incomp( me, state, bcBuffer, globBC, levelDesc, &
    &                           tree, nSize, iLevel, sim_time, neigh,   &
    &                           layout, fieldProp, varPos, nScalars,    &
    &                           varSys, derVarPos, physics, iField,     &
    &                           mixture                                 )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    ! Equilibrium Plus for bounceback
    real(kind=rk) :: fEqPlus
    real(kind=rk) :: fTmp( nScalars )
    real(kind=rk) :: fTmpPre( nScalars )
    real(kind=rk) :: fTmpPre_2( nScalars )
    integer :: iDir, iElem, QQ, elemPos
    integer :: posInBuffer
    real(kind=rk) :: swap
    real(kind=rk) :: L1, L2, L3, L4, L5
    real(kind=rk) :: dudx, dvdx, dwdx, drhodx
    real(kind=rk) :: rhoN1Prev, uN1Prev(3), rhoN2Prev, uN2Prev(3)
    real(kind=rk) :: rhoLodi, uLodi(3)
    integer :: normal(3)
    real(kind=rk) :: vel_b(globBC%nElems(iLevel)*3), inv_vel, vel_bE(3)
    integer :: bcVel_pos
    ! ---------------------------------------------------------------------------

    ! write(dbgUnit(6), *) ''
    ! write(dbgUnit(6), *) ' do characteristic velocity BounceBack for ' &
    !   &                //'bc label: '//trim(globBC%label)

    QQ = layout%fStencil%QQ
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_b                               )

    ! convert physical velocity into LB velocity
    vel_b = vel_b * inv_vel

    do iElem = 1, globBC%nElems( iLevel )

      ! Dirichlet velocity of current element
      vel_bE(:) = vel_b((iElem-1)*3+1:iElem*3)

      ! my density and velocity
      rhoLodi   = me%elemLvl( iLevel )%lodi( 1,   1, iElem )
      uLodi     = me%elemLvl( iLevel )%lodi( 2:4, 1, iElem )
      ! density and velocity of my first neighbor from last iteration
      rhoN1Prev = me%elemLvl( iLevel )%lodi( 1,   2, iElem )
      uN1Prev   = me%elemLvl( iLevel )%lodi( 2:4, 2, iElem )
      ! density and velocity of my secnd neighbor from last iteration
      rhoN2Prev = me%elemLvl( iLevel )%lodi( 1,   3, iElem )
      uN2Prev   = me%elemLvl( iLevel )%lodi( 2:4, 3, iElem )
      ! Use a swapping variable to treat +x and -x directions
      ! Incoming wave amplitudes L2, L3 and L4 in Eq (9) point inwards
      ! and outgoing wave amplitude L1 points outwards.
      ! The derivatives at +x and +y direction are computed by second-order
      ! backward stencil i.e. Eq (15) and at -x and -y, the derivatives
      ! are computed by second-order forward stencil i.e. Eq (21).
      !
      ! Swap is used to switch between forward and backward stencil.
      ! It is positive for backward stencil and negative for forward stencil.
      ! NormalInd refers to incoming direction so use cxDirInv to get normal
      ! pointing away from boundary.
      normal = layout%fStencil%cxDir(:,                                   &
        &                      layout%fStencil%cxDirInv(globBC%normalInd) )
!write(dbgUnit(6),*) 'normal ', normal
      ! x-direction
      if (abs(normal(1)) == 1) then
        swap = -1.0_rk*(real(normal(1),rk))

        ! LODI Treatment
        ! calculate spatial derivatives according to (15) which will be used to
        ! compute the characteristics
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(-8._rk*rhoLodi  + 9._rk*rhoN1Prev  - rhoN2Prev)
        dudx   =  swap*div1_3*(-8._rk*uLodi(1) + 9._rk*uN1Prev(1) - uN2Prev(1))
        dvdx   =  swap*div1_3*(-8._rk*uLodi(2) + 9._rk*uN1Prev(2) - uN2Prev(2))
        dwdx   =  swap*div1_3*(-8._rk*uLodi(3) + 9._rk*uN1Prev(3) - uN2Prev(3))

        ! wave amplitudes
        ! Outgoioming wave L1 Eq (9)
        L1 = (uLodi(1) - swap * cs) &
          & * (cs2 * drhodx - swap * cs * rho0 * dudx)
        ! incoming wave amplitudes L2, L3, L4 and L5 is approached with a
        ! linear relaxation model Eq (22)
        L2 = - (vel_bE(2)-uLodi(2))
        L3 = - (vel_bE(3)-uLodi(3))
        ! L5 is valid only if cs = sqrt(kappa*p/rho)
        L4 = 0.0_rk
        L5 = L1 - 2.0_rk*rho0*swap*cs*(vel_bE(1)-uLodi(1))

        ! LODI Equations (18)
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        ! Eq. (18b) for x-velocity
        uLodi(1) = uLodi(1) - (L5-L1)/(rho0*swap*cs)*0.5_rk
        uLodi(2) = uLodi(2) - L2  ! equation (18c) for y-direction
        uLodi(3) = uLodi(3) - L3  ! equation (18c) for z-direction
      else if( abs(normal(2)) == 1 ) then
        ! Treat y-direction
        swap = -1.0_rk*(real(normal(2),rk))

        ! LODI Treatment
        ! calculate derivativies along y-direction
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(-8._rk*rhoLodi  + 9._rk*rhoN1Prev  - rhoN2Prev)
        dudx   =  swap*div1_3*(-8._rk*uLodi(1) + 9._rk*uN1Prev(1) - uN2Prev(1))
        dvdx   =  swap*div1_3*(-8._rk*uLodi(2) + 9._rk*uN1Prev(2) - uN2Prev(2))
        dwdx   =  swap*div1_3*(-8._rk*uLodi(3) + 9._rk*uN1Prev(3) - uN2Prev(3))

        ! wave amplitudes
        ! Outgoioming wave L1 Eq (9)
        L1 = (uLodi(2) - swap * cs) &
          & * (cs2 * drhodx - swap * cs * rho0 * dvdx)
        ! Outgoing wave amplitudes L2, L3, L4 and L5 are calculated from Eq (6)
        L2 = - ( vel_bE(1) - uLodi(1) )
        L3 = - ( vel_bE(3) - uLodi(3) )
        L4 = 0.0_rk
        L5 = L1 - 2.0_rk*rho0*swap*cs*(vel_bE(2)-uLodi(2))


        ! LODI Equations
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        uLodi(1) = uLodi(1) - L2
        uLodi(2) = uLodi(2) - (L5-L1)/(rho0*swap*cs)*0.5_rk
        uLodi(3) = uLodi(3) - L3
      else if( abs(normal(3)) == 1 ) then
        ! Treat z-direction
        swap = -1.0_rk*(real( normal(3),rk))

        ! LODI Treatment
        ! calculate derivatives in z-direction
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(-8._rk*rhoLodi  + 9._rk*rhoN1Prev  - rhoN2Prev)
        dudx   =  swap*div1_3*(-8._rk*uLodi(1) + 9._rk*uN1Prev(1) - uN2Prev(1))
        dvdx   =  swap*div1_3*(-8._rk*uLodi(2) + 9._rk*uN1Prev(2) - uN2Prev(2))
        dwdx   =  swap*div1_3*(-8._rk*uLodi(3) + 9._rk*uN1Prev(3) - uN2Prev(3))

        ! wave amplitudes
        ! Outgoioming wave L1 Eq (9)
        L1 = (uLodi(3) - swap * cs) &
          & * (cs2 * drhodx - swap * cs * rho0 * dwdx)
        ! Outgoing wave amplitudes L2, L3, L4 and L5 is approached with a
        ! linear relaxation model Eq (22)
        L2 = - ( vel_bE(1) - uLodi(1) )
        L3 = - ( vel_bE(2) - uLodi(2) )
        L4 = 0.0_rk
        L5 = L1 - 2.0_rk*rho0*swap*cs*(vel_bE(3)-uLodi(3))

        ! LODI Equations
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        uLodi(1) = uLodi(1) - L2
        uLodi(2) = uLodi(2) - L3
        uLodi(3) = uLodi(3) - (L5-L1)/(rho0*swap*cs)*0.5_rk
      end if


      ! Now store the previous variables for next iteration
      me%elemLvl(iLevel)%lodi(1,   1, iElem) = rhoLodi
      me%elemLvl(iLevel)%lodi(2:4, 1, iElem) = uLodi

      do iDir = 1,QQ
        fTmpPre(iDir) = me%neigh(iLevel)%neighBufferPre_nNext(1, iDir+(iElem-1)*QQ)
        fTmpPre_2(iDir) = me%neigh(iLevel)%neighBufferPre_nNext(2, iDir+(iElem-1)*QQ)
      end do

      ! 1st neighbor
      ! density
      me%elemLvl(iLevel)%lodi(1, 2, iElem) = sum( fTmpPre(1:QQ) )
      ! Calculate my neighbor velocity for next iteration
      call derVarPos%velFromState( state  = fTmpPre,             &
        &                          iField = iField,              &
        &                          nElems = 1,                   &
        &                          varSys = varSys,              &
        &                          layout = layout,              &
        &                          res    = me%elemLvl(iLevel)   &
        &                                     %lodi(2:4,2,iElem) )

      ! 2nd neighbor
      ! density
      me%elemLvl(iLevel)%lodi(1, 3, iElem) = sum( fTmpPre_2(1:QQ) )
      ! Calculate my neighbor velocity for next iteration
      call derVarPos%velFromState( state  = fTmpPre_2,           &
        &                          iField = iField,              &
        &                          nElems = 1,                   &
        &                          varSys = varSys,              &
        &                          layout = layout,              &
        &                          res    = me%elemLvl(iLevel)   &
        &                                     %lodi(2:4,3,iElem) )

      ! Position in state array
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! Position in bcBuffer
      posInBuffer = globBC%elemLvl(iLevel)%posInBcElemBuf%val( iElem )
      fTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                    (posInBuffer-1)*nScalars+varPos(QQ)  )

      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then

          ! According to the bounce back rule,
          ! Subtract the equlibrium if its computed from alpha-(inverse of
          ! bitmask) or add the equlibrium if its computed in the direction of
          ! bitmask.
          ! In the Izquiredo paper, this is rho0, the reference density
          fEqPlus = layout%weight(iDir) * 6.0_rk * rho0             &
            &     * (   layout%fStencil%cxDirRK(1, iDir) * uLodi(1) &
            &         + layout%fStencil%cxDirRK(2, iDir) * uLodi(2) &
            &         + layout%fStencil%cxDirRK(3, iDir) * uLodi(3) )

          ! Now assign the values according to bounceback rule
          state( neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fTmp(layout%fStencil%cxDirInv(iDir)) + fEqPlus
        end if ! bitMask
      end do ! iDir
    end do ! iElem


  end subroutine inlet_nrbc_incomp
! **************************************************************************** !

! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Characteristic-based non-reflective inlet boundary conditions
  !! @todo: add explaination of these steps
  !!
  !! These boundary conditions are taken from the paper:
  !!   S. Izquierdo and N. Fueyo,
  !!   "Characteristic nonreflecting boundary
  !!   conditions for open boundaries in lattice Boltzmann methods,"
  !!   Physical Review E, vol. 78, no. 46707, 2008.
  !!
  !! > Note: unstable behavior for high omega relaxation parameters.
  !! >       It is recommended to use a sponge layer with a
  !! >       \ref mus_aux_module::mus_setspatialomega "spatial omega"
  !! >       where the omega is decreased towards the boundaries to increase
  !! >       stability
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine inlet_nrbc( me, state, bcBuffer, globBC, levelDesc, &
    &                    tree, nSize, iLevel, sim_time, neigh,   &
    &                    layout, fieldProp, varPos, nScalars,    &
    &                    varSys, derVarPos, physics, iField,     &
    &                    mixture                                 )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    ! Equilibrium Plus for bounceback
    real(kind=rk) :: fEqPlus
    real(kind=rk) :: fTmp( nScalars )
    real(kind=rk) :: fTmpPre( nScalars )
    real(kind=rk) :: fTmpPre_2( nScalars )
    integer :: iDir, iElem, QQ, elemPos
    integer :: posInBuffer
    real(kind=rk) :: swap
    real(kind=rk) :: L1, L2, L3, L4, L5
    real(kind=rk) :: dudx, dvdx, dwdx, drhodx
    real(kind=rk) :: rhoN1Prev, uN1Prev(3), rhoN2Prev, uN2Prev(3)
    real(kind=rk) :: rhoLodi, uLodi(3)
    integer :: normal(3)
    real(kind=rk) :: vel_b(globBC%nElems(iLevel)*3), inv_vel, vel_bE(3)
    integer :: bcVel_pos
    ! ---------------------------------------------------------------------------

    ! write(dbgUnit(6), *) ''
    ! write(dbgUnit(6), *) ' do characteristic velocity BounceBack for ' &
    !   &                //'bc label: '//trim(globBC%label)

    QQ = layout%fStencil%QQ
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_b                               )

    ! convert physical velocity into LB velocity
    vel_b = vel_b * inv_vel

    do iElem = 1, globBC%nElems( iLevel )

      ! Dirichlet velocity of current element
      vel_bE(:) = vel_b((iElem-1)*3+1:iElem*3)

      ! my density and velocity
      rhoLodi   = me%elemLvl( iLevel )%lodi( 1,   1, iElem )
      uLodi     = me%elemLvl( iLevel )%lodi( 2:4, 1, iElem )
      ! density and velocity of my first neighbor from last iteration
      rhoN1Prev = me%elemLvl( iLevel )%lodi( 1,   2, iElem )
      uN1Prev   = me%elemLvl( iLevel )%lodi( 2:4, 2, iElem )
      ! density and velocity of my secnd neighbor from last iteration
      rhoN2Prev = me%elemLvl( iLevel )%lodi( 1,   3, iElem )
      uN2Prev   = me%elemLvl( iLevel )%lodi( 2:4, 3, iElem )
      ! Use a swapping variable to treat +x and -x directions
      ! Incoming wave amplitudes L2, L3 and L4 in Eq (9) point inwards
      ! and outgoing wave amplitude L1 points outwards.
      ! The derivatives at +x and +y direction are computed by second-order
      ! backward stencil i.e. Eq (15) and at -x and -y, the derivatives
      ! are computed by second-order forward stencil i.e. Eq (21).
      !
      ! Swap is used to switch between forward and backward stencil.
      ! It is positive for backward stencil and negative for forward stencil.
      ! NormalInd refers to incoming direction so use cxDirInv to get normal
      ! pointing away from boundary.
      normal = layout%fStencil%cxDir(:,                                   &
        &                      layout%fStencil%cxDirInv(globBC%normalInd) )
!write(dbgUnit(6),*) 'normal ', normal
      ! x-direction
      if (abs(normal(1)) == 1) then
        swap = -1.0_rk*(real(normal(1),rk))

        ! LODI Treatment
        ! calculate spatial derivatives according to (15) which will be used to
        ! compute the characteristics
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(-8._rk*rhoLodi  + 9._rk*rhoN1Prev  - rhoN2Prev)
        dudx   =  swap*div1_3*(-8._rk*uLodi(1) + 9._rk*uN1Prev(1) - uN2Prev(1))
        dvdx   =  swap*div1_3*(-8._rk*uLodi(2) + 9._rk*uN1Prev(2) - uN2Prev(2))
        dwdx   =  swap*div1_3*(-8._rk*uLodi(3) + 9._rk*uN1Prev(3) - uN2Prev(3))

        ! wave amplitudes
        ! Outgoioming wave L1 Eq (9)
        L1 = (uLodi(1) - swap * cs) &
          & * (cs2 * drhodx - swap * cs * rhoLodi * dudx)
        ! incoming wave amplitudes L2, L3, L4 and L5 is approached with a
        ! linear relaxation model Eq (22)
        L2 = - (vel_bE(2)-uLodi(2))
        L3 = - (vel_bE(3)-uLodi(3))
        ! L5 is valid only if cs = sqrt(kappa*p/rho)
        L4 = 0.0_rk
        L5 = L1 - 2.0_rk*rhoLodi*swap*cs*(vel_bE(1)-uLodi(1))

        ! LODI Equations (18)
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        ! Eq. (18b) for x-velocity
        uLodi(1) = uLodi(1) - (L5-L1)/(rhoLodi*swap*cs)*0.5_rk
        uLodi(2) = uLodi(2) - L2  ! equation (18c) for y-direction
        uLodi(3) = uLodi(3) - L3  ! equation (18c) for z-direction
      else if( abs(normal(2)) == 1 ) then
        ! Treat y-direction
        swap = -1.0_rk*(real(normal(2),rk))

        ! LODI Treatment
        ! calculate derivativies along y-direction
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(-8._rk*rhoLodi  + 9._rk*rhoN1Prev  - rhoN2Prev)
        dudx   =  swap*div1_3*(-8._rk*uLodi(1) + 9._rk*uN1Prev(1) - uN2Prev(1))
        dvdx   =  swap*div1_3*(-8._rk*uLodi(2) + 9._rk*uN1Prev(2) - uN2Prev(2))
        dwdx   =  swap*div1_3*(-8._rk*uLodi(3) + 9._rk*uN1Prev(3) - uN2Prev(3))

        ! wave amplitudes
        ! Outgoioming wave L1 Eq (9)
        L1 = (uLodi(2) - swap * cs) &
          & * (cs2 * drhodx - swap * cs * rhoLodi * dvdx)
        ! Outgoing wave amplitudes L2, L3, L4 and L5 are calculated from Eq (6)
        L2 = - ( vel_bE(1) - uLodi(1) )
        L3 = - ( vel_bE(3) - uLodi(3) )
        L4 = 0.0_rk
        L5 = L1 - 2.0_rk*rhoLodi*swap*cs*(vel_bE(2)-uLodi(2))


        ! LODI Equations
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        uLodi(1) = uLodi(1) - L2
        uLodi(2) = uLodi(2) - (L5-L1)/(rhoLodi*swap*cs)*0.5_rk
        uLodi(3) = uLodi(3) - L3
      else if( abs(normal(3)) == 1 ) then
        ! Treat z-direction
        swap = -1.0_rk*(real( normal(3),rk))

        ! LODI Treatment
        ! calculate derivatives in z-direction
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(-8._rk*rhoLodi  + 9._rk*rhoN1Prev  - rhoN2Prev)
        dudx   =  swap*div1_3*(-8._rk*uLodi(1) + 9._rk*uN1Prev(1) - uN2Prev(1))
        dvdx   =  swap*div1_3*(-8._rk*uLodi(2) + 9._rk*uN1Prev(2) - uN2Prev(2))
        dwdx   =  swap*div1_3*(-8._rk*uLodi(3) + 9._rk*uN1Prev(3) - uN2Prev(3))

        ! wave amplitudes
        ! Outgoioming wave L1 Eq (9)
        L1 = (uLodi(3) - swap * cs) &
          & * (cs2 * drhodx - swap * cs * rhoLodi * dwdx)
        ! Outgoing wave amplitudes L2, L3, L4 and L5 is approached with a
        ! linear relaxation model Eq (22)
        L2 = - ( vel_bE(1) - uLodi(1) )
        L3 = - ( vel_bE(2) - uLodi(2) )
        L4 = 0.0_rk
        L5 = L1 - 2.0_rk*rhoLodi*swap*cs*(vel_bE(3)-uLodi(3))

        ! LODI Equations
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        uLodi(1) = uLodi(1) - L2
        uLodi(2) = uLodi(2) - L3
        uLodi(3) = uLodi(3) - (L5-L1)/(rhoLodi*swap*cs)*0.5_rk
      end if

      ! Now store the previous variables for next iteration
      me%elemLvl(iLevel)%lodi(1,   1, iElem) = rhoLodi
      me%elemLvl(iLevel)%lodi(2:4, 1, iElem) = uLodi

      do iDir = 1,QQ
        fTmpPre(iDir) = me%neigh(iLevel)%neighBufferPre_nNext(1, iDir+(iElem-1)*QQ)
        fTmpPre_2(iDir) = me%neigh(iLevel)%neighBufferPre_nNext(2, iDir+(iElem-1)*QQ)
      end do

      ! 1st neighbor
      ! density
      me%elemLvl(iLevel)%lodi(1, 2, iElem) = sum( fTmpPre(1:QQ) )
      ! Calculate my neighbor velocity for next iteration
      call derVarPos%velFromState( state  = fTmpPre,             &
        &                          iField = iField,              &
        &                          nElems = 1,                   &
        &                          varSys = varSys,              &
        &                          layout = layout,              &
        &                          res    = me%elemLvl(iLevel)   &
        &                                     %lodi(2:4,2,iElem) )

      ! 2nd neighbor
      ! density
      me%elemLvl(iLevel)%lodi(1, 3, iElem) = sum( fTmpPre_2(1:QQ) )
      ! Calculate my neighbor velocity for next iteration
      call derVarPos%velFromState( state  = fTmpPre_2,           &
        &                          iField = iField,              &
        &                          nElems = 1,                   &
        &                          varSys = varSys,              &
        &                          layout = layout,              &
        &                          res    = me%elemLvl(iLevel)   &
        &                                     %lodi(2:4,3,iElem) )

      ! Position in state array
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! Position in bcBuffer
      posInBuffer = globBC%elemLvl(iLevel)%posInBcElemBuf%val( iElem )
      fTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                    (posInBuffer-1)*nScalars+varPos(QQ)  )

      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then

          ! According to the bounce back rule,
          ! Subtract the equlibrium if its computed from alpha-(inverse of
          ! bitmask) or add the equlibrium if its computed in the direction of
          ! bitmask.
          ! In the Izquiredo paper, this is rho0, the reference density
          fEqPlus = layout%weight(iDir) * 6.0_rk * rhoLodi  &
            &     * (   layout%fStencil%cxDirRK(1, iDir) * uLodi(1) &
            &         + layout%fStencil%cxDirRK(2, iDir) * uLodi(2) &
            &         + layout%fStencil%cxDirRK(3, iDir) * uLodi(3) )

          ! Now assign the values according to bounceback rule
          state( neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fTmp(layout%fStencil%cxDirInv(iDir)) + fEqPlus
        end if ! bitMask
      end do ! iDir
    end do ! iElem


  end subroutine inlet_nrbc
! **************************************************************************** !


! ****************************************************************************** !
  !> author: Kannan Masilamani for incompressible flows
  !! Characteristic-based non-reflective open boundary conditions
  !! @todo: add explaination of these steps
  !!
  !! These boundary conditions are taken from the paper:
  !!   S. Izquierdo and N. Fueyo,
  !!   "Characteristic nonreflecting boundary
  !!   conditions for open boundaries in lattice Boltzmann methods,"
  !!   Physical Review E, vol. 78, no. 46707, 2008.
  !!
  !! > Note: unstable behavior for high omega relaxation parameters.
  !! >       It is recommended to use a sponge layer with a
  !! >       \ref mus_aux_module::mus_setspatialomega "spatial omega"
  !! >       where the omega is decreased towards the boundaries to increase
  !! >       stability
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine outlet_nrbc_incomp( me, state, bcBuffer, globBC, levelDesc, &
    &                            tree, nSize, iLevel, sim_time, neigh,   &
    &                            layout, fieldProp, varPos, nScalars,    &
    &                            varSys, derVarPos, physics, iField,     &
    &                            mixture                                 )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rhoDef(globBC%nElems(iLevel)) ! Density on boundary element
    real(kind=rk) :: rhoF    ! density on fluid element
    real(kind=rk) :: uxF(globBC%nElems(iLevel)*3)  ! Velocity on fluid element
    real(kind=rk) :: usqB, usqF ! Squared velocity magnitude
    real(kind=rk) :: fPlusFluid ! f Plus for anti-bounceback correction
    ! Equilibrium Plus for anti-bounceback correction
    real(kind=rk) :: fEqPlus
    ! Equlibrium Plus from the old time step
    real(kind=rk) :: fEqPlusFluid
    real(kind=rk) :: inv_rho_phy
    real(kind=rk) :: fTmp( nScalars * globBC%nElems(iLevel) )
    real(kind=rk) :: fTmpPre( nScalars )
    real(kind=rk) :: fTmpPre_2( nScalars )
    real(kind=rk) :: omega       ! Relaxation parameter
    integer :: iDir, iElem, QQ, invDir, elemPos
    integer :: posInBuffer, bcPress_pos
    real(kind=rk) :: sigma, kappa, csMod, Ma_L, lodi_length, k_i, swap
    real(kind=rk) :: L1, L2, L3, L4, L5
    real(kind=rk) :: dudx, dvdx, dwdx, drhodx
    real(kind=rk) :: rhoN1Prev, uN1Prev(3), rhoN2Prev, uN2Prev(3)
    real(kind=rk) :: rhoLodi, uLodi(3)
    integer :: normal(3)
    ! ---------------------------------------------------------------------------

    ! write(dbgUnit(6), *) ''
    ! write(dbgUnit(6), *) ' do non-reflecting pressure_antiBounceBack for ' &
    !   &                //'bc label: '//trim(globBC%label)

    sigma = me%nrbc%sigma
    kappa = me%nrbc%kappa  ! =1 for D2Q9, =1.66666??? for D3Q19
    csMod = me%nrbc%cs_mod
    Ma_L  = me%nrbc%Ma_L
    lodi_length = me%nrbc%lodi_length
    ! Eq (17)
    k_i = sigma*(1.0_rk - Ma_L**2)*cs/lodi_length
!write(dngUnit(6),*) 'k_i ', k_i

    QQ = layout%fStencil%QQ
    inv_rho_phy = 1.0_rk / physics%fac(iLevel)%press * cs2inv

    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure variable from spacetime function
    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
      & varSys  = varSys,                               &
      & time    = sim_time,                             &
      & iLevel  = iLevel,                               &
      & idx     = me%bc_states%pressure                 &
      &           %pntIndex%indexLvl(iLevel)            &
      &           %val(1:globBC%nElems(iLevel)),        &
      & nVals   = globBC%nElems(iLevel),                &
      & res     = rhoDef                                )

    ! convert physical pressure into LB density
    rhoDef = rhoDef * inv_rho_phy

    do iElem = 1, globBC%nElems( iLevel )
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp( (iElem-1)*nScalars+1: (iElem-1)*nScalars+QQ ) &
        &       = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                   (posInBuffer-1)*nScalars+varPos(QQ)  )
    end do

    ! Get my velocity
    call derVarPos%velFromState( state  = fTmp,                  &
      &                          iField = iField,                &
      &                          nElems = globBC%nElems(iLevel), &
      &                          varSys = varSys,                &
      &                          layout = layout,                &
      &                          res    = uxF                    )


    do iElem = 1, globBC%nElems( iLevel )

      ! my density and velocity
      rhoLodi   = me%elemLvl( iLevel )%lodi( 1,   1, iElem )
      uLodi     = me%elemLvl( iLevel )%lodi( 2:4, 1, iElem )
      ! density and velocity of my first neighbor from last iteration
      rhoN1Prev = me%elemLvl( iLevel )%lodi( 1,   2, iElem )
      uN1Prev   = me%elemLvl( iLevel )%lodi( 2:4, 2, iElem )
      ! density and velocity of my secnd neighbor from last iteration
      rhoN2Prev = me%elemLvl( iLevel )%lodi( 1,   3, iElem )
      uN2Prev   = me%elemLvl( iLevel )%lodi( 2:4, 3, iElem )
      ! Use a swapping variable to treat +x and -x directions
      ! Outgoing wave amplitudes L2, L3 and L4 in Eq (9) point outwards
      ! and incoming wave amplitude L1 points inwards.
      ! The derivatives at +x and +y direction are computed by second-order
      ! backward stencil i.e. Eq (15) and at -x and -y, the derivatives
      ! are computed by second-order forward stencil i.e. Eq (21).
      !
      ! Swap is used to switch between forward and backward stencil.
      ! It is positive for backward stencil and negative for forward stencil.
      ! NormalInd refers to incoming direction so use cxDirInv to get normal
      ! pointing away from boundary.
      normal = layout%fStencil%cxDir(:,                                   &
        &                      layout%fStencil%cxDirInv(globBC%normalInd) )
!write(dbgUnit(6),*) 'normal ', normal
      ! x-direction
      if (abs(normal(1)) == 1) then
        swap = real(normal(1),rk)

        ! LODI Treatment
        ! calculate spatial derivatives according to (15) which will be used to
        ! compute the characteristics
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(8._rk*rhoLodi  - 9._rk*rhoN1Prev  + rhoN2Prev)
        dudx   =  swap*div1_3*(8._rk*uLodi(1) - 9._rk*uN1Prev(1) + uN2Prev(1))
        dvdx   =  swap*div1_3*(8._rk*uLodi(2) - 9._rk*uN1Prev(2) + uN2Prev(2))
        dwdx   =  swap*div1_3*(8._rk*uLodi(3) - 9._rk*uN1Prev(3) + uN2Prev(3))

        ! wave amplitudes
        ! Incoming wave L1 is approached with a linear relaxation model (16)
        ! \todo 21012021 KM:check whether swap needs to be multiplied in L1
        L1 = k_i * cs2 * ( rhoLodi-rhoDef(iElem) )
        ! Outgoing wave amplitudes L2, L3, L4 and L5
        L2 = uLodi(1) * dvdx  ! equation (9) second row
        L3 = uLodi(1) * dwdx  ! equation (9) second row adapted for z-direction
        ! equation (9) third row with EOS equation (3) used
        L4 = 0.0_rk
        ! equation (9) fourth row with EOS equation (3) used
        ! MH: Here rhoLodi should be replaced with rho0
        L5 = ( uLodi(1) + swap * cs )                     &
          & * ( cs2 * drhodx + swap * cs * rho0 * dudx )

        ! LODI Equations (18)
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        ! Eq. (18b) for x-velocity
        uLodi(1) = uLodi(1) - (L5-L1)/(rho0*swap*cs)*0.5_rk
        uLodi(2) = uLodi(2) - L2  ! equation (18c) for y-direction
        uLodi(3) = uLodi(3) - L3  ! equation (18c) for z-direction
      else if( abs(normal(2)) == 1 ) then
        ! Treat y-direction
        swap = real(normal(2),rk)

        ! LODI Treatment
        ! calculate derivativies along y-direction
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(8._rk*rhoLodi  - 9._rk*rhoN1Prev  + rhoN2Prev)
        dudx   =  swap*div1_3*(8._rk*uLodi(1) - 9._rk*uN1Prev(1) + uN2Prev(1))
        dvdx   =  swap*div1_3*(8._rk*uLodi(2) - 9._rk*uN1Prev(2) + uN2Prev(2))
        dwdx   =  swap*div1_3*(8._rk*uLodi(3) - 9._rk*uN1Prev(3) + uN2Prev(3))

        ! wave amplitudes
        L1 = k_i * cs2 * ( rhoLodi-rhoDef(iElem) )
        L2 = ( uLodi(2) ) * dudx
        L3 = ( uLodi(2) ) * dwdx
        L4 = 0.0_rk
        L5 = ( uLodi(2) + swap * cs )                         &
          &* ( cs2 * drhodx + swap * cs * rho0 * dvdx )

        ! LODI Equations
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        uLodi(1) = uLodi(1) - L2
        uLodi(2) = uLodi(2) - (L5-L1)/(rho0*swap*cs)*0.5_rk
        uLodi(3) = uLodi(3) - L3
      else if( abs(normal(3)) == 1 ) then
        ! Treat z-direction
        swap = real( normal(3),rk )

        ! LODI Treatment
        ! calculate derivatives in z-direction
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(8._rk*rhoLodi  - 9._rk*rhoN1Prev  + rhoN2Prev)
        dudx   =  swap*div1_3*(8._rk*uLodi(1) - 9._rk*uN1Prev(1) + uN2Prev(1))
        dvdx   =  swap*div1_3*(8._rk*uLodi(2) - 9._rk*uN1Prev(2) + uN2Prev(2))
        dwdx   =  swap*div1_3*(8._rk*uLodi(3) - 9._rk*uN1Prev(3) + uN2Prev(3))

        ! wave amplitudes
        L1 = k_i * cs2 * ( rhoLodi-rhoDef(iElem) )
        L2 = ( uLodi(3) ) * dudx
        L3 = ( uLodi(3) ) * dvdx
        L4 = 0.0_rk
        L5 = ( uLodi(3) + swap * cs )                    &
          &* ( cs2 * drhodx + swap * cs * rho0 * dwdx )

        ! LODI Equations
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        uLodi(1) = uLodi(1) - L2
        uLodi(2) = uLodi(2) - L3
        uLodi(3) = uLodi(3) - (L5-L1)/(rho0*swap*cs)*0.5_rk
      end if

      ! Now store the previous variables for next iteration
      me%elemLvl(iLevel)%lodi(1,   1, iElem) = rhoLodi
      me%elemLvl(iLevel)%lodi(2:4, 1, iElem) = uLodi

      do iDir = 1,QQ
        fTmpPre(iDir) = me%neigh(iLevel)%neighBufferPre_nNext(1, iDir+(iElem-1)*QQ)
        fTmpPre_2(iDir) = me%neigh(iLevel)%neighBufferPre_nNext(2, iDir+(iElem-1)*QQ)
      end do

      ! 1st neighbor
      ! density
      me%elemLvl(iLevel)%lodi(1, 2, iElem) = sum( fTmpPre(1:QQ) )
      ! Calculate my neighbor velocity for next iteration
      call derVarPos%velFromState( state  = fTmpPre,             &
        &                          iField = iField,              &
        &                          nElems = 1,                   &
        &                          varSys = varSys,              &
        &                          layout = layout,              &
        &                          res    = me%elemLvl(iLevel)   &
        &                                     %lodi(2:4,2,iElem) )

      ! 2nd neighbor
      ! density
      me%elemLvl(iLevel)%lodi(1, 3, iElem) = sum( fTmpPre_2(1:QQ) )
      ! Calculate my neighbor velocity for next iteration
      call derVarPos%velFromState( state  = fTmpPre_2,           &
        &                          iField = iField,              &
        &                          nElems = 1,                   &
        &                          varSys = varSys,              &
        &                          layout = layout,              &
        &                          res    = me%elemLvl(iLevel)   &
        &                                     %lodi(2:4,3,iElem) )


      ! Fluid density
      rhoF = sum(fTmp( (iElem-1)*nScalars+1: (iElem-1)*nScalars+QQ ))

      usqB = uLodi(1)*uLodi(1) + uLodi(2)*uLodi(2) + uLodi(3)*uLodi(3)
      usqF =   uxF((iElem-1)*3+1) * uxF((iElem-1)*3+1) &
        &    + uxF((iElem-1)*3+2) * uxF((iElem-1)*3+2) &
        &    + uxF((iElem-1)*3+3) * uxF((iElem-1)*3+3)

      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! omega is STfun so get omega for current element
      omega = fieldProp%fluid%viscKine%omLvl(iLevel)%val(elemPos)

      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          ! Calc the correction term and the Dirichlet (Eq) term

          invDir = layout%fStencil%cxDirInv(iDir)

          fEqPlusFluid = layout%weight( iDir )*rhoF                  &
            &          + 4.5_rk*layout%weight( iDir )*rho0*((        &
            &            layout%fStencil%cxDirRK( 1, invDir )*uxF((iElem-1)*3+1) &
            &          + layout%fStencil%cxDirRK( 2, invDir )*uxF((iElem-1)*3+2) &
            &          + layout%fStencil%cxDirRK( 3, invDir )*uxF((iElem-1)*3+3) &
            &          )**2                                       &
            &          - div1_3*usqF)

          fEqPlus = layout%weight( iDir )*rhoLodi                &
            &     + 4.5_rk*layout%weight( iDir )*rho0*((         &
            &       layout%fStencil%cxDirRK( 1, InvDir)*uLodi(1) &
            &     + layout%fStencil%cxDirRK( 2, InvDir)*uLodi(2) &
            &     + layout%fStencil%cxDirRK( 3, InvDir)*uLodi(3) &
            &     )**2                                           &
            &     - div1_3*usqB)

          fPlusFluid = 0.5_rk*(   fTmp( (iElem-1)*nScalars + iDir )  &
            &                   + fTmp( (iElem-1)*nScalars + invDir ) )

          ! Now assign the values
          ! Actually, we should take fPlusLAst and fEqPlusLast.
          ! But this introduces errors in parallel
          state( neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) = &
          ! antibounceback term
          ! We need to get post-collision pdf in direction
          ! alpha-, which is the inverse direction of bitmask
          ! For PULL this means, get the outgoing one, as this is the one
          ! which will be bounced back
          ! For PUSH this means, get the already bounced back pdf back, so take the incoming
    &   - fTmp( (iElem-1)*nScalars + invDir ) &
    &   + 2._rk*fEqPlus     &                          ! Dirichlet pressure term
    &   + (2._rk - omega )*( fPlusFluid - fEqPlusFluid ) ! Correction term

        end if ! bitMask
      end do ! iDir
    end do ! iElem


  end subroutine outlet_nrbc_incomp
! **************************************************************************** !


! ****************************************************************************** !
  !> author: Manuel Hasert, Kannan Masilamani
  !! Characteristic-based non-reflective open boundary conditions
  !! @todo: add explaination of these steps
  !!
  !! These boundary conditions are taken from the paper:
  !!   S. Izquierdo and N. Fueyo,
  !!   "Characteristic nonreflecting boundary
  !!   conditions for open boundaries in lattice Boltzmann methods,"
  !!   Physical Review E, vol. 78, no. 46707, 2008.
  !!
  !! > Note: unstable behavior for high omega relaxation parameters.
  !! >       It is recommended to use a sponge layer with a
  !! >       \ref mus_aux_module::mus_setspatialomega "spatial omega"
  !! >       where the omega is decreased towards the boundaries to increase
  !! >       stability
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine outlet_nrbc( me, state, bcBuffer, globBC, levelDesc, &
    &                     tree, nSize, iLevel, sim_time, neigh,   &
    &                     layout, fieldProp, varPos, nScalars,    &
    &                     varSys, derVarPos, physics, iField,     &
    &                     mixture                                 )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rhoDef(globBC%nElems(iLevel)) ! Density on boundary element
    real(kind=rk) :: rhoF    ! density on fluid element
    real(kind=rk) :: uxF(globBC%nElems(iLevel)*3)  ! Velocity on fluid element
    real(kind=rk) :: usqB, usqF ! Squared velocity magnitude
    real(kind=rk) :: fPlusFluid ! f Plus for anti-bounceback correction
    ! Equilibrium Plus for anti-bounceback correction
    real(kind=rk) :: fEqPlus
    ! Equlibrium Plus from the old time step
    real(kind=rk) :: fEqPlusFluid
    real(kind=rk) :: inv_rho_phy
    real(kind=rk) :: fTmp( nScalars * globBC%nElems(iLevel) )
    real(kind=rk) :: fTmpPre( nScalars )
    real(kind=rk) :: fTmpPre_2( nScalars )
    real(kind=rk) :: omega       ! Relaxation parameter
    integer :: iDir, iElem, QQ, invDir, elemPos
    integer :: posInBuffer, bcPress_pos
    real(kind=rk) :: sigma, kappa, csMod, Ma_L, lodi_length, k_i, swap
    real(kind=rk) :: L1, L2, L3, L4, L5
    real(kind=rk) :: dudx, dvdx, dwdx, drhodx
    real(kind=rk) :: rhoN1Prev, uN1Prev(3), rhoN2Prev, uN2Prev(3)
    real(kind=rk) :: rhoLodi, uLodi(3)
    integer :: normal(3)
    ! ---------------------------------------------------------------------------

    ! write(dbgUnit(6), *) ''
    ! write(dbgUnit(6), *) ' do non-reflecting pressure_antiBounceBack for ' &
    !   &                //'bc label: '//trim(globBC%label)

    sigma = me%nrbc%sigma
    kappa = me%nrbc%kappa  ! =1 for D2Q9, =1.66666??? for D3Q19
    csMod = me%nrbc%cs_mod
    Ma_L  = me%nrbc%Ma_L
    lodi_length = me%nrbc%lodi_length
    ! Eq (17)
    k_i = sigma*(1.0_rk - Ma_L**2)*cs/lodi_length
!write(dngUnit(6),*) 'k_i ', k_i

    QQ = layout%fStencil%QQ
    inv_rho_phy = 1.0_rk / physics%fac(iLevel)%press * cs2inv

    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure variable from spacetime function
    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
      & varSys  = varSys,                               &
      & time    = sim_time,                             &
      & iLevel  = iLevel,                               &
      & idx     = me%bc_states%pressure                 &
      &           %pntIndex%indexLvl(iLevel)            &
      &           %val(1:globBC%nElems(iLevel)),        &
      & nVals   = globBC%nElems(iLevel),                &
      & res     = rhoDef                                )

    ! convert physical pressure into LB density
    rhoDef = rhoDef * inv_rho_phy

    do iElem = 1, globBC%nElems( iLevel )
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp( (iElem-1)*nScalars+1: (iElem-1)*nScalars+QQ ) &
        &       = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                   (posInBuffer-1)*nScalars+varPos(QQ)  )
    end do

    ! Get my velocity
    call derVarPos%velFromState( state  = fTmp,                  &
      &                          iField = iField,                &
      &                          nElems = globBC%nElems(iLevel), &
      &                          varSys = varSys,                &
      &                          layout = layout,                &
      &                          res    = uxF                    )


    do iElem = 1, globBC%nElems( iLevel )

      ! my density and velocity
      rhoLodi   = me%elemLvl( iLevel )%lodi( 1,   1, iElem )
      uLodi     = me%elemLvl( iLevel )%lodi( 2:4, 1, iElem )
      ! density and velocity of my first neighbor from last iteration
      rhoN1Prev = me%elemLvl( iLevel )%lodi( 1,   2, iElem )
      uN1Prev   = me%elemLvl( iLevel )%lodi( 2:4, 2, iElem )
      ! density and velocity of my secnd neighbor from last iteration
      rhoN2Prev = me%elemLvl( iLevel )%lodi( 1,   3, iElem )
      uN2Prev   = me%elemLvl( iLevel )%lodi( 2:4, 3, iElem )
      ! Use a swapping variable to treat +x and -x directions
      ! Outgoing wave amplitudes L2, L3 and L4 in Eq (9) point outwards
      ! and incoming wave amplitude L1 points inwards.
      ! The derivatives at +x and +y direction are computed by second-order
      ! backward stencil i.e. Eq (15) and at -x and -y, the derivatives
      ! are computed by second-order forward stencil i.e. Eq (21).
      !
      ! Swap is used to switch between forward and backward stencil.
      ! It is positive for backward stencil and negative for forward stencil.
      ! NormalInd refers to incoming direction so use cxDirInv to get normal
      ! pointing away from boundary.
      normal = layout%fStencil%cxDir(:,                                   &
        &                      layout%fStencil%cxDirInv(globBC%normalInd) )
!write(dbgUnit(6),*) 'normal ', normal
      ! x-direction
      if (abs(normal(1)) == 1) then
        swap = real(normal(1),rk)

        ! LODI Treatment
        ! calculate spatial derivatives according to (15) which will be used to
        ! compute the characteristics
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(8._rk*rhoLodi  - 9._rk*rhoN1Prev  + rhoN2Prev)
        dudx   =  swap*div1_3*(8._rk*uLodi(1) - 9._rk*uN1Prev(1) + uN2Prev(1))
        dvdx   =  swap*div1_3*(8._rk*uLodi(2) - 9._rk*uN1Prev(2) + uN2Prev(2))
        dwdx   =  swap*div1_3*(8._rk*uLodi(3) - 9._rk*uN1Prev(3) + uN2Prev(3))

        ! wave amplitudes
        ! Incoming wave L1 is approached with a linear relaxation model (16)
        ! \todo 21012021 KM:check whether swap needs to be multiplied in L1
        L1 = k_i * cs2 * ( rhoLodi-rhoDef(iElem) )
        ! Outgoing wave amplitudes L2, L3, L4 and L5
        L2 = uLodi(1) * dvdx  ! equation (9) second row
        L3 = uLodi(1) * dwdx  ! equation (9) second row adapted for z-direction
        ! equation (9) third row with EOS equation (3) used
        L4 = 0.0_rk
        ! equation (9) fourth row with EOS equation (3) used
        ! MH: Here rhoLodi should be replaced with rho0
        L5 = ( uLodi(1) + swap * cs )                     &
          & * ( cs2 * drhodx + swap * cs * rhoLodi * dudx )

        ! LODI Equations (18)
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        ! Eq. (18b) for x-velocity
        uLodi(1) = uLodi(1) - (L5-L1)/(rhoLodi*swap*cs)*0.5_rk
        uLodi(2) = uLodi(2) - L2  ! equation (18c) for y-direction
        uLodi(3) = uLodi(3) - L3  ! equation (18c) for z-direction
      else if( abs(normal(2)) == 1 ) then
        ! Treat y-direction
        swap = real(normal(2),rk)

        ! LODI Treatment
        ! calculate derivativies along y-direction
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(8._rk*rhoLodi  - 9._rk*rhoN1Prev  + rhoN2Prev)
        dudx   =  swap*div1_3*(8._rk*uLodi(1) - 9._rk*uN1Prev(1) + uN2Prev(1))
        dvdx   =  swap*div1_3*(8._rk*uLodi(2) - 9._rk*uN1Prev(2) + uN2Prev(2))
        dwdx   =  swap*div1_3*(8._rk*uLodi(3) - 9._rk*uN1Prev(3) + uN2Prev(3))

        ! wave amplitudes
        L1 = k_i * cs2 * ( rhoLodi-rhoDef(iElem) )
        L2 = ( uLodi(2) ) * dudx
        L3 = ( uLodi(2) ) * dwdx
        L4 = 0.0_rk
        L5 = ( uLodi(2) + swap * cs )                         &
          &* ( cs2 * drhodx + swap * cs * rhoLodi * dvdx )

        ! LODI Equations
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        uLodi(1) = uLodi(1) - L2
        uLodi(2) = uLodi(2) - (L5-L1)/(rhoLodi*swap*cs)*0.5_rk
        uLodi(3) = uLodi(3) - L3
      else if( abs(normal(3)) == 1 ) then
        ! Treat z-direction
        swap = real( normal(3),rk )

        ! LODI Treatment
        ! calculate derivatives in z-direction
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(8._rk*rhoLodi  - 9._rk*rhoN1Prev  + rhoN2Prev)
        dudx   =  swap*div1_3*(8._rk*uLodi(1) - 9._rk*uN1Prev(1) + uN2Prev(1))
        dvdx   =  swap*div1_3*(8._rk*uLodi(2) - 9._rk*uN1Prev(2) + uN2Prev(2))
        dwdx   =  swap*div1_3*(8._rk*uLodi(3) - 9._rk*uN1Prev(3) + uN2Prev(3))

        ! wave amplitudes
        L1 = k_i * cs2 * ( rhoLodi-rhoDef(iElem) )
        L2 = ( uLodi(3) ) * dudx
        L3 = ( uLodi(3) ) * dvdx
        L4 = 0.0_rk
        L5 = ( uLodi(3) + swap * cs )                    &
          &* ( cs2 * drhodx + swap * cs * rhoLodi * dwdx )

        ! LODI Equations
        rhoLodi  = rhoLodi  - ( L4 + 0.5_rk*(L5+L1) ) * cs2inv
        uLodi(1) = uLodi(1) - L2
        uLodi(2) = uLodi(2) - L3
        uLodi(3) = uLodi(3) - (L5-L1)/(rhoLodi*swap*cs)*0.5_rk
      end if

      ! Now store the previous variables for next iteration
      me%elemLvl(iLevel)%lodi(1,   1, iElem) = rhoLodi
      me%elemLvl(iLevel)%lodi(2:4, 1, iElem) = uLodi

      do iDir = 1,QQ
        fTmpPre(iDir) = me%neigh(iLevel)%neighBufferPre_nNext(1, iDir+(iElem-1)*QQ)
        fTmpPre_2(iDir) = me%neigh(iLevel)%neighBufferPre_nNext(2, iDir+(iElem-1)*QQ)
      end do

      ! 1st neighbor
      ! density
      me%elemLvl(iLevel)%lodi(1, 2, iElem) = sum( fTmpPre(1:QQ) )
      ! Calculate my neighbor velocity for next iteration
      call derVarPos%velFromState( state  = fTmpPre,             &
        &                          iField = iField,              &
        &                          nElems = 1,                   &
        &                          varSys = varSys,              &
        &                          layout = layout,              &
        &                          res    = me%elemLvl(iLevel)   &
        &                                     %lodi(2:4,2,iElem) )

      ! 2nd neighbor
      ! density
      me%elemLvl(iLevel)%lodi(1, 3, iElem) = sum( fTmpPre_2(1:QQ) )
      ! Calculate my neighbor velocity for next iteration
      call derVarPos%velFromState( state  = fTmpPre_2,           &
        &                          iField = iField,              &
        &                          nElems = 1,                   &
        &                          varSys = varSys,              &
        &                          layout = layout,              &
        &                          res    = me%elemLvl(iLevel)   &
        &                                     %lodi(2:4,3,iElem) )


      ! Fluid density
      rhoF = sum(fTmp( (iElem-1)*nScalars+1: (iElem-1)*nScalars+QQ ))

      usqB = uLodi(1)*uLodi(1) + uLodi(2)*uLodi(2) + uLodi(3)*uLodi(3)
      usqF =   uxF((iElem-1)*3+1) * uxF((iElem-1)*3+1) &
        &    + uxF((iElem-1)*3+2) * uxF((iElem-1)*3+2) &
        &    + uxF((iElem-1)*3+3) * uxF((iElem-1)*3+3)

      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! omega is STfun so get omega for current element
      omega = fieldProp%fluid%viscKine%omLvl(iLevel)%val(elemPos)

      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          ! Calc the correction term and the Dirichlet (Eq) term

          invDir = layout%fStencil%cxDirInv(iDir)

          fEqPlusFluid = layout%weight( iDir )*rhoF                  &
            &          + 4.5_rk*layout%weight( iDir )*rhoF*((        &
            &            layout%fStencil%cxDirRK( 1, invDir )*uxF((iElem-1)*3+1) &
            &          + layout%fStencil%cxDirRK( 2, invDir )*uxF((iElem-1)*3+2) &
            &          + layout%fStencil%cxDirRK( 3, invDir )*uxF((iElem-1)*3+3) &
            &          )**2                                       &
            &          - div1_3*usqF)

          fEqPlus = layout%weight( iDir )*rhoLodi                &
            &     + 4.5_rk*layout%weight( iDir )*rhoLodi*((      &
            &       layout%fStencil%cxDirRK( 1, InvDir)*uLodi(1) &
            &     + layout%fStencil%cxDirRK( 2, InvDir)*uLodi(2) &
            &     + layout%fStencil%cxDirRK( 3, InvDir)*uLodi(3) &
            &     )**2                                           &
            &     - div1_3*usqB)

          fPlusFluid = 0.5_rk*(   fTmp( (iElem-1)*nScalars + iDir )  &
            &                   + fTmp( (iElem-1)*nScalars + invDir ) )

          ! Now assign the values
          ! Actually, we should take fPlusLAst and fEqPlusLast.
          ! But this introduces errors in parallel
          state( neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) = &
          ! antibounceback term
          ! We need to get post-collision pdf in direction
          ! alpha-, which is the inverse direction of bitmask
          ! For PULL this means, get the outgoing one, as this is the one
          ! which will be bounced back
          ! For PUSH this means, get the already bounced back pdf back, so take the incoming
    &   - fTmp( (iElem-1)*nScalars + invDir ) &
    &   + 2._rk*fEqPlus     &                          ! Dirichlet pressure term
    &   + (2._rk - omega )*( fPlusFluid - fEqPlusFluid ) ! Correction term

        end if ! bitMask
      end do ! iDir
    end do ! iElem


  end subroutine outlet_nrbc
! **************************************************************************** !



! ****************************************************************************** !
  !> author: Manuel Hasert
  !! Characteristic-based non-reflective open boundary conditions
  !!
  !! These boundary conditions are taken from the paper:
  !!   S. Izquierdo and N. Fueyo,
  !!   "Characteristic nonreflecting boundary
  !!   conditions for open boundaries in lattice Boltzmann methods,"
  !!   Physical Review E, vol. 78, no. 46707, 2008.
  !!
  !! > Note: unstable behavior for high omega relaxation parameters.
  !! >       It is recommended to use a sponge layer with a
  !! >       \ref mus_aux_module::mus_setspatialomega "spatial omega"
  !! >       where the omega is decreased towards the boundaries to increase
  !! >       stability
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine outlet_nrbc_eq( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                        nSize, iLevel, sim_time, neigh, layout,         &
    &                        fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                        physics, iField, mixture                        )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rhoDef(globBC%nElems( iLevel ))
    real(kind=rk) :: Ma, swap, sigma, kappa, K_mod
    real(kind=rk) :: L1, L2, L3, L4, L5
    real(kind=rk) :: dudx, dvdx, dwdx, drhodx
    real(kind=rk) :: rhoN1Prev, uN1Prev(3), rhoN2Prev, uN2Prev(3)
    real(kind=rk) :: rhoLodi, uLodi(3), csMod
    real(kind=rk) :: fEq( layout%fStencil%QQ )
    real(kind=rk) :: fTmpPost( nScalars )
    real(kind=rk) :: fTmpPost_2( nScalars )
    real(kind=rk) :: inv_rho_phy
    integer :: iELem, iDir, QQ
    integer :: posInBuffer, bcPress_pos
    ! ---------------------------------------------------------------------------

    QQ = layout%fStencil%QQ
    inv_rho_phy = 1.0_rk / physics%fac(iLevel)%press * cs2inv

    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure variable from spacetime function
    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
      & varSys  = varSys,                               &
      & time    = sim_time,                             &
      & iLevel  = iLevel,                               &
      & idx     = me%bc_states%pressure                 &
      &           %pntIndex%indexLvl(iLevel)            &
      &           %val(1:globBC%nElems(iLevel)),        &
      & nVals   = globBC%nElems(iLevel),                &
      & res     = rhoDef                                )

    ! convert physical pressure to LB density
    rhoDef = rhoDef * inv_rho_phy

    sigma = me%nrbc%sigma
    kappa = me%nrbc%kappa  ! =1 for D2Q9, D3Q19
    csMod = me%nrbc%cs_mod

    do iElem = 1, globBC%nElems( iLevel )
      rhoLodi = me%elemLvl( iLevel )%lodi( 1,   1, iElem )
      uLodi   = me%elemLvl( iLevel )%lodi( 2:4, 1, iElem )
      ! density and velocity of my first neighbor from last iteration
      rhoN1Prev = me%elemLvl( iLevel )%lodi( 1,   2, iElem )
      uN1Prev   = me%elemLvl( iLevel )%lodi( 2:4, 2, iElem )
      ! density and velocity of my secnd neighbor from last iteration
      rhoN2Prev = me%elemLvl( iLevel )%lodi( 1,   3, iElem )
      uN2Prev   = me%elemLvl( iLevel )%lodi( 2:4, 3, iElem )

      ! Use a swapping variable to treat +x and -x directions
      swap = 1._rk
      ! x-direction
      if( globBC%elemLvl( iLevel )%normal%val( 1, iElem ) /= 0 ) then
        ! Use a swapping variable to treat +x and -x directions
        swap = 1._rk*(-real( globBC%elemLvl( iLevel )%normal%val( 1, iElem ),kind=rk))

        ! LODI Treatment
        Ma   = abs(uLodi(1))*csInv ! X direction
        ! Equation (17) for the linear relaxation model constant
        ! lodi_length= distance outlet-inlet
        K_mod = swap*sigma*(1._rk - Ma*Ma)*cs/(me%nrbc%lodi_length)

        ! calculate spatial derivatives according to (15) which will be used to
        ! compute the characteristics
        ! This is a quadratic extrapolation
        drhodx =  swap*div1_3*(8._rk*rhoLodi  - 9._rk*rhoN1Prev  + rhoN2Prev)
        dudx   =  swap*div1_3*(8._rk*uLodi(1) - 9._rk*uN1Prev(1) + uN2Prev(1))
        dvdx   =  swap*div1_3*(8._rk*uLodi(2) - 9._rk*uN1Prev(2) + uN2Prev(2))
        dwdx   =  swap*div1_3*(8._rk*uLodi(3) - 9._rk*uN1Prev(3) + uN2Prev(3))

        ! wave amplitudes
        ! Incoming wave L1 is approached with a linear relaxation model (16)
        L1 = K_mod*csMod*csMod*(rhoLodi-rhoDef(iElem))
        L2 = (uLodi(1))*dvdx  ! equation (9) second row
        L3 = (uLodi(1))*dwdx  ! equation (9) second row adapted for z-direction
        ! equation (9) third row with EOS equation (3) used
        L4 = (uLodi(1))*csMod*csMod*drhodx*(1._rk-1._rk/kappa)
        ! equation (9) fourth row with EOS equation (3) used
        ! MH: Here rhoLodi should be replaced with rho0
        L5 = ( uLodi(1) + swap * csMod )                                       &
          &* ( csMod * csMod * drhodx / kappa                                  &
          &+   swap * csMod * dudx )

        ! LODI Equations (18)
        rhoLodi  = rhoLodi  - (L4+0.5_rk*(L5+L1))/(csMod*csMod)
        ! Eq. (18b) for x-velocity
        uLodi(1) = uLodi(1) - ((L5-L1))/(swap*csMod)*0.5_rk
        uLodi(2) = uLodi(2) - L2  ! equation (18c) for y-direction
        uLodi(3) = uLodi(3) - L3  ! equation (18c) for z-direction
      else
        ! Treat y-direction
        if( globBC%elemLvl( iLevel )%normal%val( 2, iElem ) /= 0 ) then
          swap = 1._rk*(-real( globBC%elemLvl( iLevel )%normal%val(2,iElem),kind=rk ))

          ! LODI Treatment
          Ma   = abs(uLodi(2))*csInv ! X direction
          ! lodi_length= distance outlet-inlet
          K_mod = swap*sigma*(1._rk - Ma*Ma)*cs/(me%nrbc%lodi_length)

          ! calculate derivatives
          ! This is a quadratic extrapolation
          drhodx =  swap*div1_3*(8._rk*rhoLodi  - 9._rk*rhoN1Prev  + rhoN2Prev)
          dudx   =  swap*div1_3*(8._rk*uLodi(1) - 9._rk*uN1Prev(1) + uN2Prev(1))
          dvdx   =  swap*div1_3*(8._rk*uLodi(2) - 9._rk*uN1Prev(2) + uN2Prev(2))
          dwdx   =  swap*div1_3*(8._rk*uLodi(3) - 9._rk*uN1Prev(3) + uN2Prev(3))

          ! wave amplitudes
          L1 = K_mod*csMod*csMod*(rhoLodi-rhoDef(iElem))
          L2 = (uLodi(2))*dudx
          L3 = (uLodi(2))*dwdx
          L4 = (uLodi(2))*csMod*csMod*drhodx*(1._rk-1._rk/kappa)
          L5 = ( uLodi(2) + swap * csMod )                                     &
            &* ( csMod * csMod * drhodx / kappa                                &
            &+   swap * csMod * dvdx )

          ! LODI Equations
          rhoLodi  = rhoLodi  - (L4+0.5_rk*(L5+L1))/(csMod*csMod)
          uLodi(1) = uLodi(1) - L2
          uLodi(2) = uLodi(2) - ((L5-L1))/(swap*csMod)*0.5_rk
          uLodi(3) = uLodi(3) - L3
        else
          ! Treat z-direction
          if( globBC%elemLvl( iLevel )%normal%val( 3, iElem ) /= 0 ) then
            swap = 1._rk*(-real(globBC%elemLvl( iLevel )%normal%val(3,iElem),kind=rk))

            ! LODI Treatment
            Ma   = abs(uLodi(3))*csInv ! X direction
            !lodi_length= distance outlet-inlet
            K_mod = swap*sigma*(1._rk - Ma*Ma)*cs/(me%nrbc%lodi_length)

            ! calculate derivatives
           ! This is a quadratic extrapolation
            drhodx =  swap*div1_3*(8._rk*rhoLodi  - 9._rk*rhoN1Prev  + rhoN2Prev)
            dudx   =  swap*div1_3*(8._rk*uLodi(1) - 9._rk*uN1Prev(1) + uN2Prev(1))
            dvdx   =  swap*div1_3*(8._rk*uLodi(2) - 9._rk*uN1Prev(2) + uN2Prev(2))
            dwdx   =  swap*div1_3*(8._rk*uLodi(3) - 9._rk*uN1Prev(3) + uN2Prev(3))

            ! wave amplitudes
            L1 = K_mod*csMod*csMod*(rhoLodi-rhoDef(iElem))
            L2 = (uLodi(3))*dudx
            L3 = (uLodi(3))*dvdx
            L4 = (uLodi(3))*csMod*csMod*drhodx*(1._rk-1._rk/kappa)
            L5 = ( uLodi(3) + swap * csMod )                                   &
              &* ( csMod * csMod * drhodx / kappa                              &
              &+   swap * csMod * dwdx )

            ! LODI Equations
            rhoLodi  = rhoLodi  - (L4+0.5_rk*(L5+L1))/(csMod*csMod)
            uLodi(1) = uLodi(1) - L2
            uLodi(2) = uLodi(2) - L3
            uLodi(3) = uLodi(3) - ((L5-L1))/(swap*csMod)*0.5_rk
          end if
        end if
      end if

      ! Now store the previous variables for next iteration
      me%elemLvl( iLevel )%lodi( 1,   1, iElem ) = rhoLodi
      me%elemLvl( iLevel )%lodi( 2:4, 1, iElem ) = uLodi
      ! density and velocity of my first neighbor
      ! MH: Velocity and Density should be calculated with access as SAVE
      ! For pull, this is fine however, as IDX and SAVE coincide
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmpPost(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                        (posInBuffer-1)*nScalars+varPos(QQ)  )
      do iDir = 1, QQ
        fTmpPost_2(iDir) = me%neigh( iLevel )%neighBufferPost( 1,iDir+(iElem-1)*QQ)
      end do
      !density
      me%elemLvl( iLevel )%lodi( 1, 2, iElem ) = sum( fTmpPost(1:QQ) )
      !velocity
      call derVarPos%velFromState( state  = fTmpPost,   &
        &                          iField = iField,     &
        &                          nElems = 1,      &
        &                          varSys = varSys, &
        &                          layout = layout, &
        &                          res    = me%elemLvl(iLevel)%lodi(2:4,2,iElem) )

      ! density and velocity of my secnd neighbor
      me%elemLvl( iLevel )%lodi( 1, 3, iElem )  = sum( fTmpPost_2(1:QQ) )
      call derVarPos%velFromState( state  = fTmpPost_2, &
        &                          iField = iField,     &
        &                          nElems = 1,          &
        &                          varSys = varSys,     &
        &                          layout = layout,     &
        &                          res    = me%elemLvl(iLevel)%lodi(2:4,3,iElem) )

      ! Calculate the equilibrium distribution
      call derVarPos%equilFromState( state  = fTmpPost_2, &
        &                            iField = iField,     &
        &                            nElems = 1,          &
        &                            varSys = varSys,     &
        &                            layout = layout,     &
        &                            res     = fEq )

      ! Now do the anti-bounce back
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          ! Now assign the values
          state(                                                                &
& neigh (( idir-1)* nsize+ globbc%elemlvl( ilevel )%elem%val( ielem ))+( ifield-1)* qq+ nscalars*0 &
&         ) = fEq( iDir )
        end if
      end do ! iDir
    end do

  end subroutine outlet_nrbc_eq
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Moment based open boundary condition from Sam Bennent PhD thesis
  !! "A Lattice Boltzmann Model for Diffusion of Binary Gas Mixtures"
  !! for weakly compressible LBM model
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'inlet',
  !!    kind = 'pressure_momentsbased',
  !!    pressure = 'p0' }
  !!}
  !! variable = {
  !!   name = 'p0',
  !!   ncomponents = 1,
  !!   vartype = 'st_fun',
  !!   st_fun = 1.0
  !! }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type::fnct]] function pointer.
  subroutine pressure_momentsbased( me, state, bcBuffer, globBC, levelDesc, &
    &                               tree, nSize, iLevel, sim_time, neigh,   &
    &                               layout, fieldProp, varPos, nScalars,    &
    &                               varSys, derVarPos, physics, iField,     &
    &                               mixture                                 )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: pressB( globBC%nElems(iLevel) ), rho, vel(3), inv_press
    integer :: iELem, iDir
    integer :: nLinks, iLink, elemPos, QQ
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp(:)
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: first_moments(3), second_moments(6), third_moments(6), &
      &              fourth_moments(3)
    integer :: nMom(4)
    integer :: posInBuffer, bcPress_pos
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    inv_press = 1.0_rk / physics%fac( iLevel )%press

!    write(*,*) 'Boundary label ', trim(me%label)
    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure variable from spacetime function
    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
      & varSys  = varSys,                               &
      & time    = sim_time,                             &
      & iLevel  = iLevel,                               &
      & idx     = me%bc_states%pressure                 &
      &           %pntIndex%indexLvl(iLevel)            &
      &           %val(1:globBC%nElems(iLevel)),        &
      & nVals   = globBC%nElems(iLevel),                &
      & res     = pressB                                )

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
    pressB = pressB * inv_press

    nMom(1) = 1 + size(layout%moment%first_moments)
    nMom(2) = nMom(1) + size(layout%moment%second_moments)
    nMom(3) = nMom(2) + size(layout%moment%third_moments)
    nMom(4) = nMom(3) + size(layout%moment%fourth_moments)

    do iElem = 1, globBC%nElems( iLevel )
      fTmp = 0.0_rk
      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(unKnown_fTmp(nLinks))
      allocate(rhs(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )

      iLink = 0
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink + 1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = &
            &state(  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = &
&state(  neigh (( layout%fstencil%restposition-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )

      ! density from defined pressure
      rho = cs2inv*pressB(iElem)

      !local velocity
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      vel(1) = sum ( bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                      (posInBuffer-1)*nScalars+varPos(QQ) ) &
        &            * layout%fStencil%cxDir(1,:) )
      vel(2) = sum ( bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                      (posInBuffer-1)*nScalars+varPos(QQ) ) &
        &            * layout%fStencil%cxDir(2,:) )
      vel(3) = sum ( bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                      (posInBuffer-1)*nScalars+varPos(QQ) ) &
        &            * layout%fStencil%cxDir(3,:) )
      vel = vel/rho

      ! moments for in-compressible model
      first_moments = rho*vel !momX, momY, momZ
      second_moments = (/ pressB(iElem) + rho*vel(1)*vel(1), & !momXX
        &                 pressB(iElem) + rho*vel(2)*vel(2), & !momYY
        &                 pressB(iElem) + rho*vel(3)*vel(3), & !momZZ
        &                 rho*vel(1)*vel(2),                 & !momXY
        &                 rho*vel(2)*vel(3),                 & !momYZ
        &                 rho*vel(3)*vel(1) /)                 !momXZ
      third_moments = (/ cs2*rho*vel(2), & !momXXY
        &                cs2*rho*vel(3), & !momXXZ
        &                cs2*rho*vel(1), & !momYYX
        &                cs2*rho*vel(3), & !momYYZ
        &                cs2*rho*vel(1), & !momZZX
        &                cs2*rho*vel(2) /) !momZZY
      !fourth_moments = cs2*pressB(iElem) !momXXYY, momYYZZ, momZZXX
      fourth_moments = (/ &
        & cs2*(pressB(iElem)+rho*( vel(1)**2 + vel(2)**2 - 0.5_rk*vel(3)**2 )), &!momXXYY
        & cs2*(pressB(iElem)+rho*(-0.5_rk*vel(1)**2 + vel(2)**2 + vel(3)**2 )), &!momYYZZ
        & cs2*(pressB(iElem)+rho*( vel(1)**2 - 0.5_rk*vel(2)**2 + vel(3)**2 ))  &!momZZXX
        & /)

      ! moments for pressure-wall corner
      ! KM: for pressure-velocity, velocity at boundary node must be defined
      moments = 0.0_rk
      moments(1) = rho
      moments(2:nMom(1)) = &
        & first_moments(layout%moment%first_moments(:))
      moments(nMom(1)+1:nMom(2)) = &
        & second_moments(layout%moment%second_moments(:))
      moments(nMom(2)+1:nMom(3)) = &
        & third_moments(layout%moment%third_moments(:))
      moments(nMom(3)+1:nMom(4)) = &
        & fourth_moments(layout%moment%fourth_moments(:))

      ! compute rhs of lse to compute unknown pdfs
      do iLink = 1, nLinks
        rhs(iLink) = &
          & moments( me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink) ) &
          ! convert pdf to moments
          & - dot_product(layout%moment%toMoments%A(                        &
          &    me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink), :), fTmp)
      end do

      ! unknown pdfs
      unKnown_fTmp = matmul( &
        & me%elemLvl(iLevel)%moments(iElem)%unKnownPdfs_MatInv, rhs )

      do iLink = 1, nLinks
        state(                                                                 &
          & neigh (( missing_links(ilink)-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
          & )  = unKnown_fTmp(iLink)
      end do

     deallocate(missing_links)
     deallocate(unKnown_fTmp)
     deallocate(rhs)
   end do

  end subroutine pressure_momentsbased
! ****************************************************************************** !



! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Moment based open boundary condition from Sam Bennent PhD thesis
  !! "A Lattice Boltzmann Model for Diffusion of Binary Gas Mixtures"
  !! for incompressible LBM model
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'inlet',
  !!    kind = 'pressure_momentsbased',
  !!    pressure = 'p0'}
  !!}
  !! variable = {
  !!   name = 'p0',
  !!   ncomponents = 1,
  !!   vartype = 'st_fun',
  !!   st_fun = 1.0
  !! }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine pressure_momentsbased_incomp( me, state, bcBuffer, globBC,        &
    &                                      levelDesc, tree, nSize, iLevel,     &
    &                                      sim_time, neigh, layout, fieldProp, &
    &                                      varPos, nScalars, varSys,           &
    &                                      derVarPos, physics, iField, mixture )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: pressB( globBC%nElems(iLevel) ), rho, vel(3)
    real(kind=rk) :: inv_press
    integer :: iELem, iDir
    integer :: nLinks, iLink, QQ, elemPos
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp(:)
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: first_moments(3), second_moments(6), third_moments(6), &
      &              fourth_moments(3)
    integer :: nMom(4)
    integer :: posInBuffer, bcPress_pos
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    inv_press = 1.0_rk / physics%fac( iLevel )%press

    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure variable from spacetime function
    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
      & varSys  = varSys,                               &
      & time    = sim_time,                             &
      & iLevel  = iLevel,                               &
      & idx     = me%bc_states%pressure                 &
      &           %pntIndex%indexLvl(iLevel)            &
      &           %val(1:globBC%nElems(iLevel)),        &
      & nVals   = globBC%nElems(iLevel),                &
      & res     = pressB                                )

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
    pressB = pressB * inv_press

    nMom(1) = 1 + size(layout%moment%first_moments)
    nMom(2) = nMom(1) + size(layout%moment%second_moments)
    nMom(3) = nMom(2) + size(layout%moment%third_moments)
    nMom(4) = nMom(3) + size(layout%moment%fourth_moments)

    do iElem = 1, globBC%nElems( iLevel )
      fTmp = 0.0_rk
      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(unKnown_fTmp(nLinks))
      allocate(rhs(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )

      iLink = 0
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink + 1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = &
            & state( neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = &
& state( neigh (( layout%fstencil%restposition-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )

      !local velocity
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      vel(1) = sum ( bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                      (posInBuffer-1)*nScalars+varPos(QQ) ) &
        &            * layout%fStencil%cxDir(1,:) )
      vel(2) = sum ( bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                      (posInBuffer-1)*nScalars+varPos(QQ) ) &
        &            * layout%fStencil%cxDir(2,:) )
      vel(3) = sum ( bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                      (posInBuffer-1)*nScalars+varPos(QQ) ) &
        &            * layout%fStencil%cxDir(3,:) )
      vel = vel*rho0Inv

      ! density from defined pressure
      rho = cs2inv*pressB(iElem)

      ! moments for in-compressible model
      !first_moments = (/ 0.0_rk, 0.0_rk, 0.0_rk /) !momX, momY, momZ
      first_moments = rho0*vel !momX, momY, momZ
      second_moments = (/ pressB(iElem) + rho0*vel(1)*vel(1), & !momXX
        &                 pressB(iElem) + rho0*vel(2)*vel(2), & !momYY
        &                 pressB(iElem) + rho0*vel(3)*vel(3), & !momZZ
        &                 rho0*vel(1)*vel(2),                 & !momXY
        &                 rho0*vel(2)*vel(3),                 & !momYZ
        &                 rho0*vel(3)*vel(1) /)                 !momXZ
      third_moments = (/ cs2*rho0*vel(2), & !momXXY
        &                cs2*rho0*vel(3), & !momXXZ
        &                cs2*rho0*vel(1), & !momYYX
        &                cs2*rho0*vel(3), & !momYYZ
        &                cs2*rho0*vel(1), & !momZZX
        &                cs2*rho0*vel(2) /) !momZZY
      !fourth_moments = cs2*pressB(iElem) !momXXYY, momYYZZ, momZZXX
      fourth_moments = (/ &
        & cs2*(pressB(iElem)+rho0*( vel(1)**2 + vel(2)**2 - 0.5_rk*vel(3)**2 )), &!momXXYY
        & cs2*(pressB(iElem)+rho0*(-0.5_rk*vel(1)**2 + vel(2)**2 + vel(3)**2 )), &!momYYZZ
        & cs2*(pressB(iElem)+rho0*( vel(1)**2 - 0.5_rk*vel(2)**2 + vel(3)**2 ))  &!momZZXX
        & /)

      ! moments for pressure-wall corner
      ! KM: for pressure-velocity, velocity at boundary node must be defined
      moments = 0.0_rk
      moments(1) = rho
      moments(2:nMom(1)) = &
        & first_moments(layout%moment%first_moments(:))
      moments(nMom(1)+1:nMom(2)) = &
        & second_moments(layout%moment%second_moments(:))
      moments(nMom(2)+1:nMom(3)) = &
        & third_moments(layout%moment%third_moments(:))
      moments(nMom(3)+1:nMom(4)) = &
        & fourth_moments(layout%moment%fourth_moments(:))

      ! compute rhs of lse to compute unknown pdfs
      do iLink = 1, nLinks
        rhs(iLink) = &
          & moments( me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink) ) &
          ! convert pdf to moments
          & - dot_product(layout%moment%toMoments%A(                        &
          &    me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink), :), fTmp)
      end do

      ! unknown pdfs
      unKnown_fTmp = matmul( &
        & me%elemLvl(iLevel)%moments(iElem)%unKnownPdfs_MatInv, rhs )

      do iLink = 1, nLinks
        state(                                                                 &
          & neigh (( missing_links(ilink)-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
          & )  = unKnown_fTmp(iLink)
      end do

     deallocate(missing_links)
     deallocate(unKnown_fTmp)
     deallocate(rhs)
   end do

  end subroutine pressure_momentsbased_incomp
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Moment based wall boundary condition from Sam Bennent PhD thesis
  !! "A Lattice Boltzmann Model for Diffusion of Binary Gas Mixtures"
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'inlet',
  !!    kind = 'moments_wall'}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine moments_wall( me, state, bcBuffer, globBC, levelDesc, tree,      &
    &                      nSize, iLevel, sim_time, neigh, layout, fieldProp, &
    &                      varPos, nScalars, varSys, derVarPos, physics,      &
    &                      iField, mixture                                    )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rhoB, press
    integer :: iELem, iDir, tmpDir
    integer :: nLinks, iLink, QQ, elemPos
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp( : )
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: first_moments(3), second_moments(6), third_moments(6), &
      &              fourth_moments(3)
    integer :: nMom(4)
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ

    nMom(1) = 1 + size(layout%moment%first_moments)
    nMom(2) = nMom(1) + size(layout%moment%second_moments)
    nMom(3) = nMom(2) + size(layout%moment%third_moments)
    nMom(4) = nMom(3) + size(layout%moment%fourth_moments)

    do iElem = 1, globBC%nElems( iLevel )
      fTmp = 0.0_rk
      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(rhs(nLinks))
      allocate(unKnown_fTmp(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )

      iLink = 0
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink + 1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = &
            &state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = &
& state( neigh (( layout%fstencil%restposition-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )

      rhoB = fTmp(layout%fStencil%restPosition)
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          tmpDir = layout%fStencil%cxDirInv(iDir)
          if( globBC%elemLvl(iLevel)%bitmask%val( tmpDir, iElem )) then
            tmpDir = layout%fStencil%cxDirInv(globBC%elemLvl(iLevel)%          &
                &                                         normalInd%val(iElem))
          end if
        else
          tmpDir = iDir
        end if
        rhoB  = rhoB + fTmp(tmpDir)
      end do
      press = cs2 * rhoB

      ! moments for in-compressible model
      first_moments = (/ 0.0_rk, 0.0_rk, 0.0_rk /) !momX, momY, momZ
      second_moments = (/ press, & !momXX
        &                 press, & !momYY
        &                 press, & !momZZ
        &                 0.0_rk,& !momXY
        &                 0.0_rk,& !momYZ
        &                 0.0_rk /)!momXZ
      third_moments = 0.0_rk
      fourth_moments = cs2*press !momXXYY, momYYZZ, momZZXX

      ! moments for pressure-wall corner
      moments = 0.0_rk
      moments(1) = rhoB
      moments(2:nMom(1)) = &
        & first_moments(layout%moment%first_moments(:))
      moments(nMom(1)+1:nMom(2)) = &
        & second_moments(layout%moment%second_moments(:))
      moments(nMom(2)+1:nMom(3)) = &
        & third_moments(layout%moment%third_moments(:))
      moments(nMom(3)+1:nMom(4)) = &
        & fourth_moments(layout%moment%fourth_moments(:))

!     write(*,*) 'Moments ', moments

      do iLink = 1, nLinks
        rhs(iLink) = &
          & moments( me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink) ) &
          ! convert pdf to moments
          & - dot_product(layout%moment%toMoments%A(                        &
          &    me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink), :), fTmp)
      end do
     !write(*,*) 'rhs ', rhs

     unKnown_fTmp = matmul( &
       & me%elemLvl(iLevel)%moments(iElem)%unKnownPdfs_MatInv, rhs )

     !write(*,*) 'unknown fTmp ', unKnown_fTmp

      do iLink = 1, nLinks
        state(                                                         &
&neigh (( missing_links(ilink)-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
         & )  = unKnown_fTmp(iLink)
      end do

      deallocate(missing_links)
      deallocate(rhs)
      deallocate(unKnown_fTmp)
    end do

  end subroutine moments_wall
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Moment based velocity boundary condition for weakly compressible LBM model
  !! Based on Sam Bennent PhD thesis
  !! "A Lattice Boltzmann Model for Diffusion of Binary Gas Mixtures"
  !!  Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'inlet',
  !!    kind = 'velocity_momentsbased'
  !!    velocityX = 0.06,
  !!    velocityY = 0.0,
  !!    velocityZ = 0.0 }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine velocity_momentsbased( me, state, bcBuffer, globBC, levelDesc, &
    &                               tree, nSize, iLevel, sim_time, neigh,   &
    &                               layout, fieldProp, varPos, nScalars,    &
    &                               varSys, derVarPos, physics, iField,     &
    &                               mixture                                 )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rho, press, dens_correction, vel(3)
    real(kind=rk) :: vel_b(globBC%nElems(iLevel)*3), inv_vel
    integer :: iELem, iDir, tmpDir
    integer :: nLinks, iLink, QQ, elemPos
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp( : )
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: first_moments(3), second_moments(6), third_moments(6), &
      &              fourth_moments(3)
    integer :: nMom(4)
    integer :: bcVel_pos
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel
    !write(*,*) 'Boundary label ', trim(me%label)

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_b                               )

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
    !write(*,*) 'phy Vel', ux, uy, uz
    vel_b = vel_b * inv_vel

    nMom(1) = 1 + size(layout%moment%first_moments)
    nMom(2) = nMom(1) + size(layout%moment%second_moments)
    nMom(3) = nMom(2) + size(layout%moment%third_moments)
    nMom(4) = nMom(3) + size(layout%moment%fourth_moments)

    do iElem = 1, globBC%nElems( iLevel )

      vel(:) = vel_b( (iElem-1)*3+1 : iElem*3 )

      fTmp = 0.0_rk
      dens_correction = 0.0_rk
      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(rhs(nLinks))
      allocate(unKnown_fTmp(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )

      iLink = 0
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink +1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = &
            &state( neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = &
& state( neigh (( layout%fstencil%restposition-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )

      rho = fTmp(layout%fStencil%restPosition)
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          tmpDir = layout%fStencil%cxDirInv(iDir)
          if( globBC%elemLvl(iLevel)%bitmask%val( tmpDir, iElem )) then
            tmpDir = layout%fStencil%cxDirInv(globBC%elemLvl(iLevel)%          &
                &                                         normalInd%val(iElem))
          end if
        else
          tmpDir = iDir
        end if
        rho  = rho + fTmp(tmpDir)
      end do

      ! add small correction term to density which comes from momentum in the
      ! flow direction while substituing unknowns pdfs terms obtained
      ! from moments BC derivation to compute density.
      ! With this correction term density is recovered correctly at this node
      dens_correction = dot_product(vel, globBC%elemLvl(iLevel)%normal%val(:, iElem))
      ! for compressible model, rho_new = rho_old + rho_new*dens_correction
      ! thus, rho_new = rho_old/(1-dens_correction)
      rho = rho/(1.0_rk-dens_correction)

      press = cs2 * rho

      ! moments for in-compressible model
      first_moments = vel !momX, momY, momZ
      second_moments = (/ press + rho*vel(1)*vel(1), & !momXX
        &                 press + rho*vel(2)*vel(2), & !momYY
        &                 press + rho*vel(3)*vel(3), & !momZZ
        &                 rho*vel(1)*vel(2),                 & !momXY
        &                 rho*vel(2)*vel(3),                 & !momYZ
        &                 rho*vel(3)*vel(1) /)                 !momXZ
      third_moments = (/ cs2*rho*vel(2), & !momXXY
        &                cs2*rho*vel(3), & !momXXZ
        &                cs2*rho*vel(1), & !momYYX
        &                cs2*rho*vel(3), & !momYYZ
        &                cs2*rho*vel(1), & !momZZX
        &                cs2*rho*vel(2) /) !momZZY
      !fourth_moments = cs2*pressB(iElem) !momXXYY, momYYZZ, momZZXX
      fourth_moments = (/ &
        & cs2*(press+rho*( vel(1)**2 + vel(2)**2 - 0.5_rk*vel(3)**2 )), &!momXXYY
        & cs2*(press+rho*(-0.5_rk*vel(1)**2 + vel(2)**2 + vel(3)**2 )), &!momYYZZ
        & cs2*(press+rho*( vel(1)**2 - 0.5_rk*vel(2)**2 + vel(3)**2 ))  &!momZZXX
        & /)

      ! moments for pressure-wall corner
      ! KM: for pressure-velocity, velocity at boundary node must be defined
      moments = 0.0_rk
      moments(1) = rho
      moments(2:nMom(1)) = &
        & first_moments(layout%moment%first_moments(:))
      moments(nMom(1)+1:nMom(2)) = &
        & second_moments(layout%moment%second_moments(:))
      moments(nMom(2)+1:nMom(3)) = &
        & third_moments(layout%moment%third_moments(:))
      moments(nMom(3)+1:nMom(4)) = &
        & fourth_moments(layout%moment%fourth_moments(:))

      do iLink = 1, nLinks
        rhs(iLink) = &
          & moments( me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink) ) &
          ! convert pdf to moments
          & - dot_product(layout%moment%toMoments%A(                        &
          &    me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink), :), fTmp)
      end do

      unKnown_fTmp = matmul( &
        & me%elemLvl(iLevel)%moments(iElem)%unKnownPdfs_MatInv, rhs )

      do iLink = 1, nLinks
        state(                                                                 &
&neigh (( missing_links(ilink)-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
         & )  = unKnown_fTmp(iLink)
      end do

      deallocate(missing_links)
      deallocate(rhs)
      deallocate(unKnown_fTmp)

    end do

  end subroutine velocity_momentsbased
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Moment based velocity boundary condition for incompressible LBM model
  !! Based on Sam Bennent PhD thesis
  !! "A Lattice Boltzmann Model for Diffusion of Binary Gas Mixtures"
  !!  Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'inlet',
  !!    kind = 'velocity_momentsbased'
  !!    velocityX = 0.06,
  !!    velocityY = 0.0,
  !!    velocityZ = 0.0 }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine velocity_momentsbased_incomp( me, state, bcBuffer, globBC,        &
    &                                      levelDesc, tree, nSize, iLevel,     &
    &                                      sim_time, neigh, layout, fieldProp, &
    &                                      varPos, nScalars, varSys,           &
    &                                      derVarPos, physics, iField, mixture )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: rho, press, dens_correction, vel(3)
    real(kind=rk) :: vel_b(globBC%nElems(iLevel)*3), inv_vel
    integer :: iELem, iDir, tmpDir
    integer :: nLinks, iLink, QQ, elemPos
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp( : )
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: first_moments(3), second_moments(6), third_moments(6), &
      &              fourth_moments(3)
    integer :: nMom(4)
    integer :: bcVel_pos
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel
    !write(*,*) 'Boundary label ', trim(me%label)

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_b                               )

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
    vel_b = vel_b * inv_vel

    nMom(1) = 1 + size(layout%moment%first_moments)
    nMom(2) = nMom(1) + size(layout%moment%second_moments)
    nMom(3) = nMom(2) + size(layout%moment%third_moments)
    nMom(4) = nMom(3) + size(layout%moment%fourth_moments)

    do iElem = 1, globBC%nElems( iLevel )

      vel(:) = vel_b( (iElem-1)*3+1 : iElem*3 )

      fTmp = 0.0_rk
      dens_correction = 0.0_rk
      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(rhs(nLinks))
      allocate(unKnown_fTmp(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )

      iLink = 0
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink +1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = &
            &  state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = &
& state( neigh (( layout%fstencil%restposition-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )

      rho = fTmp( layout%fStencil%restPosition )
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          tmpDir = layout%fStencil%cxDirInv(iDir)
          if( globBC%elemLvl(iLevel)%bitmask%val( tmpDir, iElem )) then
            tmpDir = layout%fStencil%cxDirInv(globBC%elemLvl(iLevel)%          &
                &                                         normalInd%val(iElem))
          end if
        else
          tmpDir = iDir
        end if
        rho  = rho + fTmp(tmpDir)
      end do

      ! add small correction term to density which comes from momentum in the
      ! flow direction while substituing unknowns pdfs terms obtained
      ! from moments BC derivation to compute density.
      ! With this correction term density is recovered correctly at this node
      dens_correction = dot_product(vel, globBC%elemLvl(iLevel)%normal%val(:, iElem))
      ! for compressible model, rho_new = rho_old + rho_new*dens_correction
      ! thus, rho_new = rho_old/(1-dens_correction)
      ! correction term for in-compressible mode
      rho = rho + rho0*dens_correction

      press = cs2 * rho

      ! moments for in-compressible model
      first_moments = vel !momX, momY, momZ
      second_moments = (/ press + rho0*vel(1)*vel(1), & !momXX
        &                 press + rho0*vel(2)*vel(2), & !momYY
        &                 press + rho0*vel(3)*vel(3), & !momZZ
        &                 rho0*vel(1)*vel(2),                 & !momXY
        &                 rho0*vel(2)*vel(3),                 & !momYZ
        &                 rho0*vel(3)*vel(1) /)                 !momXZ
      third_moments = (/ cs2*rho0*vel(2), & !momXXY
        &                cs2*rho0*vel(3), & !momXXZ
        &                cs2*rho0*vel(1), & !momYYX
        &                cs2*rho0*vel(3), & !momYYZ
        &                cs2*rho0*vel(1), & !momZZX
        &                cs2*rho0*vel(2) /) !momZZY
      !fourth_moments = cs2*pressB(iElem) !momXXYY, momYYZZ, momZZXX
      fourth_moments = (/ &
        & cs2*(press+rho0*( vel(1)**2 + vel(2)**2 - 0.5_rk*vel(3)**2 )), &!momXXYY
        & cs2*(press+rho0*(-0.5_rk*vel(1)**2 + vel(2)**2 + vel(3)**2 )), &!momYYZZ
        & cs2*(press+rho0*( vel(1)**2 - 0.5_rk*vel(2)**2 + vel(3)**2 ))  &!momZZXX
        & /)

      ! moments for pressure-wall corner
      ! KM: for pressure-velocity, velocity at boundary node must be defined
      moments = 0.0_rk
      moments(1) = rho
      moments(2:nMom(1)) = &
        & first_moments(layout%moment%first_moments(:))
      moments(nMom(1)+1:nMom(2)) = &
        & second_moments(layout%moment%second_moments(:))
      moments(nMom(2)+1:nMom(3)) = &
        & third_moments(layout%moment%third_moments(:))
      moments(nMom(3)+1:nMom(4)) = &
        & fourth_moments(layout%moment%fourth_moments(:))

      do iLink = 1, nLinks
        rhs(iLink) = &
          & moments( me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink) ) &
          ! convert pdf to moments
          & - dot_product(layout%moment%toMoments%A(                        &
          &    me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink), :), fTmp)
      end do

      unKnown_fTmp = matmul( &
        & me%elemLvl(iLevel)%moments(iElem)%unKnownPdfs_MatInv, rhs )

      do iLink = 1, nLinks
        state(                                                         &
&neigh (( missing_links(ilink)-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
         & )  = unKnown_fTmp(iLink)
      end do

      deallocate(missing_links)
      deallocate(rhs)
      deallocate(unKnown_fTmp)
    end do

  end subroutine velocity_momentsbased_incomp
! ****************************************************************************** !


! ****************************************************************************** !
  !> Boundary condition for APESMATE that supports pdfs as input variable
  !! Makes it possible to use pdf from, for example, the left to the right domain
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'some_label',
  !!   kind  = 'bc_pdf',
  !!   pdf   = 'some_value'
  !! },
  !! ...
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine bc_pdf( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
    &                iLevel, sim_time, neigh, layout, fieldProp, varPos,   &
    &                nScalars, varSys, derVarPos, physics, iField, mixture )
    ! -------------------------------------------------------------------- !
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
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: in_pdf( layout%fStencil%QQ * globBC%nElems(iLevel) )
    !> get number of directions from stencil
    integer :: QQ
    !> position of pdf in variable system
    integer :: bcPdf_pos
    integer :: iELem, iDir, elemPos
    ! ---------------------------------------------------------------------------

    QQ = layout%fStencil%QQ

    ! position of boundary pdf in varSys
    bcPdf_pos = me%bc_states%pdf%varPos
    ! Get pdf
    call varSys%method%val(bcPdf_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%pdf                    &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = in_pdf                              )

    do iElem = 1, globBC%nElems( iLevel )
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          state(  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
            & )  = in_pdf((iElem-1)*QQ+iDir)
        end if
      end do
    end do

  end subroutine bc_pdf
! ****************************************************************************** !


end module mus_bc_fluid_module
! ****************************************************************************** !
