! Copyright (c) 2012-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> Boundary condition treatment routines for multispecies simulation
!!
module mus_bc_species_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,               only: rk, labelLen
  use tem_param_module,         only: div1_3, div1_36, div1_8, div3_4h,        &
    &                                 div1_4, div3_8, div3_2, div9_16, div3_16,&
    &                                 cs2inv, div1_6, cs4inv,                  &
    &                                 t2cs2inv, t2cs4inv, div1_18, cs2
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_math_module,          only: invert_matrix
  use tem_spatial_module,       only: tem_spatial_for
  use tem_spacetime_fun_module, only: tem_spacetime_for
  use tem_isNaN_module,         only: tem_isNaN
  use tem_debug_module,         only: dbgUnit
  use tem_property_module,      only: prp_solid
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_stencil_module,       only: tem_stencilHeader_type

  ! include musubi modules
  use mus_bc_header_module,       only: boundary_type, glob_boundary_type
  use mus_scheme_layout_module,   only: mus_scheme_layout_type
  use mus_field_prop_module,      only: mus_field_prop_type
  use mus_derVarPos_module,       only: mus_derVarPos_type
  use mus_param_module,           only: mus_param_type
  use mus_physics_module,         only: mus_physics_type
  use mus_mixture_module,         only: mus_mixture_type
  use mus_species_module,         only: mus_species_type
  use mus_eNRTL_module,           only: mus_calc_thermFactor,                  &
    &                                   mus_calc_MS_DiffMatrix
  use mus_varSys_module,          only: mus_varSys_data_type
  use mus_derQuanMSLiquid_module, only: momentumFromMacroLSE


  implicit none

  private

  public :: spc_inlet_eq, spc_inlet, spc_vel_bb

  public :: spc_outlet_zero_prsgrd
  public :: spc_outlet_eq, spc_outlet_vel, spc_outlet_expol
  public :: spc_moleFrac, spc_moleDiff_Flux
  public :: spc_moleFlux, spc_moleFlux_eq
  public :: spc_moleDens_eq
  public :: spc_blackbox_mem_ion, spc_blackbox_mem_solvent
  public :: spc_moleFrac_wtdf, spc_moleFrac_eq
  public :: spc_outflow, spc_inflow
  public :: spc_solvent_outflow, spc_solvent_inflow
  public :: spc_velocity_noneq_expol, spc_mole_fraction_noneq_expol
  ! moments BC
  public :: spc_moments_moleFrac
  public :: spc_moments_moleFlux
  public :: spc_moments_wall
  public :: spc_moments_vel

contains
! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Outlet boundary conditions with zero pressure gradient.
  !!
  !! These boundary conditions use the neighbor cells density and velocity in
  !! the equilibrium function. It is not necessary to specify density at
  !! boundary in the lua configuration file
  !! \( f = f^{eq}(\rho_{b-1},u_{b-1}) \)
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_outlet_zero_prsgrd( me, state, bcBuffer, globBC, levelDesc, &
    &                                tree, nSize, iLevel, sim_time, neigh,   &
    &                                layout, fieldProp, varPos, nScalars,    &
    &                                varSys, derVarPos, physics, iField,     &
    &                                mixture                                 )
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
    real(kind=rk) :: rho, ux(3), usq, m_ratio
    integer :: iElem, iDir, iField_2, pos, QQ
    real(kind=rk) :: fTmp( layout%fStencil%QQ * varSys%nStateVars )
    integer :: elemPos, posInBuffer
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ

    m_ratio = fieldProp%species%molWeigRatio

    do iElem = 1, globBC%nElems( iLevel )
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
      do iField_2 = 1, varSys%nStateVars
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iField_2)%state_varPos(iDir)
          fTmp(pos) = bcBuffer( &
          & ( posinbuffer-1)* nscalars+ pos )
        end do
      end do

      !Calculate local density and velocity moments
      rho = 0._rk
      do iDir = 1,layout%fStencil%QQ
        rho = rho + fTmp(varPos(iDir))
      end do

      call derVarPos%velFromState( state  = fTmp,   &
        &                          iField = iField, &
        &                          nElems = 1,      &
        &                          varSys = varSys, &
        &                          layout = layout, &
        &                          res    = ux      )

      usq = ux(1)*ux(1) + ux(2)*ux(2) + ux(3)*ux(3)

      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
        ! Write the values
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          ! Only changed from save to fetch and added currT
          state( neigh((idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) =  &
            &       layout%weight( iDir )*rho*(2.0_rk*m_ratio+9._rk*( (        &
            &       layout%fStencil%cxDir( 1, layout%fStencil                  &
            &                                       %cxDirInv( iDir ))*ux(1)   &
            &    +  layout%fStencil%cxDir( 2, layout%fStencil                  &
            &                                       %cxDirInv( iDir ))*ux(2)   &
            &    +  layout%fStencil%cxDir( 3, layout%fStencil                  &
            &                                       %cxDirInv( iDir ))*ux(3)   &
            &    )**2 -  div1_3*usq))                                          &
            &    - state(                                                      &
            &  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)
        end if
      end do
    end do

  end subroutine spc_outlet_zero_prsgrd
! ****************************************************************************** !

! ****************************************************************************** !
  !> Inlet boundary condition for defined species velocity and mole fraction
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_inlet',
  !!    velocity = 'inlet_vel',
  !!    mole_fraction = 'inlet_mole'
  !! }
  !!}
  !!variable = {
  !! {
  !!   label = 'inlet_vel',
  !!   ncomponents = 3,
  !!   vartype = 'st_fun',
  !!   st_fun = {0.1,0.0,0.0}
  !! },
  !! {
  !!   label = 'inlet_mole',
  !!   ncomponents = 1,
  !!   vartype = 'st_fun',
  !!   st_fun = 0.1
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_inlet( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
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
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: mass_dens
    real(kind=rk) :: vel_b(globBC%nElems(iLevel)*3), inv_vel
    real(kind=rk) :: moleFrac(globBC%nElems(iLevel))
    real(kind=rk) :: massFlux(3)
    integer :: iElem, iDir, QQ, posInBuffer, bcVel_pos, offset, elemPos
    integer :: bcMoleFrac_pos
    ! ------------------------------------------------------------------------
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

    bcMoleFrac_pos = me%bc_states%moleFrac%varPos
    ! mole fraction
    call varSys%method%val(bcMoleFrac_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFrac                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFrac(:)                              )

    do iElem = 1, globBC%nElems(iLevel)
      offset = (iElem-1)*3

      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                    (posInBuffer-1)*nScalars+varPos(QQ)  )

      ! local density
      !mass_dens = sum(fTmp)
      mass_dens = mixture%moleDens0LB *  moleFrac(iElem) &
        &       * fieldProp%species%molWeight

      massFlux = mass_dens * vel_b(offset+1:offset+3)

      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state( neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
              & = fTmp(layout%fStencil%cxDirInv( iDir ))                      &
              & + layout%weight( iDir )*2._rk*cs2inv                          &
              &    * ( layout%fStencil%cxDirRK( 1, iDir )*massFlux(1)         &
              &    +   layout%fStencil%cxDirRK( 2, iDir )*massFlux(2)         &
              &    +   layout%fStencil%cxDirRK( 3, iDir )*massFlux(3)         )
        end if
      end do
    end do
  end subroutine spc_inlet
! ****************************************************************************** !


! ****************************************************************************** !
  !> Inlet species velocity equilibrium boundary with specified
  !! mixture averaged mass velocity and its molefraction
  !! mixture kinematic pressure is extrapolated here.
  !! Density and velocity of all fields are used to compute equilibrium
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_inlet_eq',
  !!    velocity = 'inlet_vel',
  !!    mole_fraction = 'inlet_mole'
  !! }
  !!}
  !!variable = {
  !! {
  !!   label = 'inlet_vel',
  !!   ncomponents = 3,
  !!   vartype = 'st_fun',
  !!   st_fun = {0.1,0.0,0.0}
  !! },
  !! {
  !!   label = 'inlet_mole',
  !!   ncomponents = 1,
  !!   vartype = 'st_fun',
  !!   st_fun = 0.1
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_inlet_eq( me, state, bcBuffer, globBC, levelDesc, tree,      &
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
    real(kind=rk) :: fTmp( layout%fStencil%QQ*globBC%nElems(iLevel) &
      &              * varSys%nStateVars )
    ! fEq uses AOS layout
    real(kind=rk) :: fEq( layout%fStencil%QQ*globBC%nElems(iLevel) )
    real(kind=rk) :: mass_dens( globBC%nElems(iLevel)*varSys%nStateVars )
    !real(kind=rk) :: uxB(3*globBC%nElems(iLevel))   !< Velocity on boundary element
    !> Velocity on boundary element
    real(kind=rk) :: uxB(globBC%nElems(iLevel)*varSys%nStateVars*3)
    real(kind=rk) :: velocity(3,globBC%nElems(iLevel)*varSys%nStateVars)
    real(kind=rk) :: spc_vel(globBC%nElems(iLevel)*3), inv_vel
    real(kind=rk) :: moleFrac(globBC%nElems(iLevel))
    integer :: iElem, iDir, iFieldLoc, nFields, pos, QQ
    integer :: offset, bcVel_pos, bcMoleFrac_pos, elemPos, posInBuffer
    ! ------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars
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
      & res     = spc_vel                             )

    ! convert physical velocity into LB velocity
    spc_vel = spc_vel * inv_vel

    ! position of molefraction spacetime variable in varSys
    bcMoleFrac_pos = me%bc_states%moleFrac%varPos
    ! mole fraction
    call varSys%method%val(bcMoleFrac_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFrac                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFrac                                 )

    ! copy state
    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
      do iFieldLoc = 1, nFields
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp( ( ielem-1)* nscalars+ pos ) = &
            & bcBuffer(                                                 &
            & ( posinbuffer-1)* nscalars+ pos)
        end do
      end do !iField
    end do !iElem

    ! Derive all species velocities at once since
    ! it is inefficient to derive each species velocity
    call derVarPos%velocitiesFromState( state  = fTmp,                  &
      &                                 iField = iField,                &
      &                                 nElems = globBC%nElems(iLevel), &
      &                                 varSys = varSys,                &
      &                                 layout = layout,                &
      &                                 res    = uxB                    )

    ! get density and velocity.
    ! if current field, use velocity and density defined in lua file.
    ! for other fields, derive density and velocity from state
    mass_dens = 0.0_rk
    do iFieldLoc = 1, nFields
      if (iFieldLoc == iField) then
        ! store velocity in input_loc array
        do iElem = 1, globBC%nElems(iLevel)
          offset = (iElem-1)*nFields + iFieldLoc
          velocity(:,offset) = spc_vel((iElem-1)*3+1 : iElem*3)

          ! compute current species mass density from specified molefraction at
          ! boundary
          ! rho = n_t * chi_i * m_i + massFrac_i * rho0 * KinePress/(cs2*phi_i)
          ! massFrac_i = rho_i/rho
          ! rho = n_t * chi_i * m_i + rho_i * KinePress/(cs2*phi_i)
          ! 1st term is zero order density
          ! 2nd term is second order density with kinematic mixture pressure,
          ! p = cs2*(sum(phi_k*rho_k) - min_a(m_a)*n0)/rho0
          !KM: using inital mixture number density rho0/mixtureMOlWeight
          !Using local tot_NuMdens increases density over time
          mass_dens( offset ) = mixture%moleDens0LB * moleFrac(iElem) &
            &                 * fieldProp%species%molWeight
        end do
      else
        do iElem = 1, globBC%nElems(iLevel)
          offset = (iElem-1)*nFields + iFieldLoc
          velocity(:,offset) =                             &
            & uxB((iElem-1)*nFields*3 + (iFieldLoc-1)*3 + 1 : &
            &     (iElem-1)*nFields*3 + iFieldLoc*3 )
          ! species density
          do iDir = 1, layout%fStencil%QQ
            pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
            mass_dens(offset) = mass_dens(offset) &
              & + fTmp(( ielem-1)* nscalars+ pos )
          end do !iDir
        end do !iElem
      end if
    end do !iField

    ! compute equilibrium
    call derVarPos%equilFromMacro( density  = mass_dens,             &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq                    )

    do iElem = 1, globBC%nElems(iLevel)
      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq( iDir+(iElem-1)*QQ )
        end if
      end do
    end do

  end subroutine spc_inlet_eq
! ****************************************************************************** !

! ****************************************************************************** !
  !> Inlet species velocity bounce back boundary with specified
  !! mixture averaged mass velocity and its molefraction
  !! mixture kinematic pressure is extrapolated here.
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_vel_bb',
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
  subroutine spc_vel_bb( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
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
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: mass_dens
    real(kind=rk) :: vel_b(globBC%nElems(iLevel)*3), inv_vel
    integer :: iElem, iDir, QQ, posInBuffer, bcVel_pos, offset, elemPos
    ! ------------------------------------------------------------------------
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

    do iElem = 1, globBC%nElems(iLevel)
      offset = (iElem-1)*3

      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                    (posInBuffer-1)*nScalars+varPos(QQ)  )
      ! local density
      mass_dens = sum(fTmp)

      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state( neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
              & = fTmp(layout%fStencil%cxDirInv( iDir ))                      &
              & + layout%weight( iDir )*2._rk*cs2inv*mass_dens                &
              &    * ( layout%fStencil%cxDirRK( 1, iDir )*vel_b(offset+1)     &
              &    +   layout%fStencil%cxDirRK( 2, iDir )*vel_b(offset+2)     &
              &    +   layout%fStencil%cxDirRK( 3, iDir )*vel_b(offset+3)     )
        end if
      end do
    end do
  end subroutine spc_vel_bb
! ****************************************************************************** !

  ! ************************************************************************ !
  !> species Outlet Pressure extrapolation boundary. NOT VERIFIED
  !!
  !! This is taken from the paper:
  !! M. Junk and Z. Yang,
  !! Asymptotic Analysis of Lattice Boltzmann Outflow Treatments
  !! Commun. Comput. Phys., pp. 1–11, 2011.
  !!
  !! * Pressure as defined in the configuration file
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_outlet_expol',
  !!    pressure = 'p1'
  !!  }
  !!}
  !!variable = {
  !!  name = 'p1',
  !!  ncomponents = 1,
  !!  vartype = 'st_fun',
  !!  st_fun = 1.0
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_outlet_expol( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                          nSize, iLevel, sim_time, neigh, layout,         &
    &                          fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                          physics, iField, mixture                        )
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
    real(kind=rk) :: fTmpAll( layout%fStencil%QQ*globBC%nElems(iLevel) &
      &              * varSys%nStateVars )
    real(kind=rk) :: fTmp_1
    real(kind=rk) :: fTmp_2
    ! fEq uses AOS layout
    real(kind=rk) :: fEq( layout%fStencil%QQ*globBC%nElems(iLevel) )
    integer :: iDir              !< Direction counter
    integer :: iElem             !< Element counter
    integer :: pos, ifieldloc, nFields, neighPos, QQ, elemPos
    ! ------------------------------------------------------------------------

    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    ! Run over all elements in list for this boundary
    do iElem = 1, globBC%nElems(iLevel)
      neighPos = me%neigh(iLevel)%posInState(1,iElem)
      do iFieldLoc = 1, nFields
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          !neighbor element
          fTmpAll( ( ielem-1)* nscalars+ pos )     &
            & = state((neighpos-1)*nscalars+ idir+(ifieldloc-1)*qq)
          !local element
          !fTmpAll( ( ielem-1)* nscalars+ pos ) = bcBuffer(                    &
          !  & ( globbc%elemlvl( ilevel )%posinbcelembuf%val( ielem )-1)* nscalars+ pos )
        end do
      end do !iField
    end do !iElem

    ! Calculate the equilibrium distribution
    call derVarPos%equilFromState( state  = fTmpAll,                 &
      &                            iField = iField,                  &
      &                            nElems = globBC%nElems( iLevel ), &
      &                            varSys = varSys,                  &
      &                            layout = layout,                  &
      &                            res    = fEq                      )


    ! Run over all elements in list for this boundary
    do iElem = 1, globBC%nElems(iLevel)
      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      ! Treat all directions, but actually we only want to treat
      ! non-normal directions
      ! Equation (3.2b)
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
        ! KM: pre-collision of current time step must be used for extrapolation
          fTmp_1 = me%neigh( iLevel )%neighBufferPre_nNext( 1,  iDir+(iElem-1)*QQ)
          fTmp_2 = me%neigh( iLevel )%neighBufferPre_nNext( 2,  iDir+(iElem-1)*QQ)

          ! update the incoming velocity direction
           state(neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = 1.5_rk*fTmp_1 - 0.5_rk*fTmp_2
            !& = feq( iDir+(iElem-1)*QQ )
        end if
      end do

      ! then overwrite the normal direction with special treatment
      ! Equation (3.2a)
      iDir = globBC%elemLvl( iLevel )%normalInd%val( iElem )
state( neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)  = &
        &  feq( ( ielem-1)* qq+idir ) !&
        !& +(fTmpAll(( ielem-1)* nscalars+varpos(layout%fstencil%cxdirinv(idir)))&
        !& -feq( ( ielem-1)* layout%fstencil%qq+layout%fstencil%cxdirinv(idir) ) )
        !& fTmpAll(( ielem-1)* nscalars+varpos(idir))
    end do !iElem

  end subroutine spc_outlet_expol
! ****************************************************************************** !


! ****************************************************************************** !
  !> Outlet species velocity equilibrium boundary with specified
  !! mixture averaged mass velocity.
  !! molefraction is extrapolated here.
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_outlet_vel',
  !!    velocity = 'inlet_vel',
  !! }
  !!}
  !!variable = {
  !! {
  !!   label = 'inlet_vel',
  !!   ncomponents = 3,
  !!   vartype = 'st_fun',
  !!   st_fun = {0.1,0.0,0.0}
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_outlet_vel( me, state, bcBuffer, globBC, levelDesc, tree,   &
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
    real(kind=rk) :: fTmp( layout%fStencil%QQ*globBC%nElems(iLevel) &
      &              * varSys%nStateVars )
    ! fEq uses AOS layout
    real(kind=rk) :: fEq( layout%fStencil%QQ*globBC%nElems(iLevel) )
    real(kind=rk) :: mass_dens( globBC%nElems(iLevel)*varSys%nStateVars )
    !> Velocity on boundary element
    real(kind=rk) :: uxB(globBC%nElems(iLevel)*varSys%nStateVars*3)
    real(kind=rk) :: velocity(3,globBC%nElems(iLevel)*varSys%nStateVars)
    real(kind=rk) :: spc_vel(globBC%nElems(iLevel)*3), inv_vel
    integer :: iElem, iDir, iFieldLoc, nFields, pos, QQ
    integer :: offset, bcVel_pos, elemPos, posInBuffer
    ! ------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars
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
      & res     = spc_vel                             )

    ! convert physical velocity into LB velocity
    spc_vel = spc_vel * inv_vel

    ! copy state
    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      do iFieldLoc = 1, nFields
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp( ( ielem-1)* nscalars+ pos )   &
            & = bcBuffer(                                               &
            & ( posinbuffer-1)* nscalars+ pos )
        end do
      end do !iField
    end do !iElem

    ! Derive all species velocities at once since
    ! it is inefficient to derive each species velocity
    call derVarPos%velocitiesFromState( state  = fTmp,                  &
      &                                 iField = iField,                &
      &                                 nElems = globBC%nElems(iLevel), &
      &                                 varSys = varSys,                &
      &                                 layout = layout,                &
      &                                 res    = uxB                    )

    ! use same velocity for all species and extrapolate density
    do iFieldLoc = 1, nFields
        do iElem = 1, globBC%nElems(iLevel)
          offset = (iElem-1)*nFields + iFieldLoc
          velocity(:,offset) = spc_vel((iElem-1)*3+1 : iElem*3)
        end do
    end do !iField

    ! species density
    mass_dens = 0.0_rk
    do iElem = 1, globBC%nElems(iLevel)
      do iFieldLoc = 1, nFields
        offset = (iElem-1)*nFields + iFieldLoc
        do iDir = 1, layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          mass_dens(offset) = mass_dens(offset) &
            & + fTmp(( ielem-1)* nscalars+ pos )
        end do !iDir
      end do !iField
    end do !iElem

    ! compute equilibrium
    call derVarPos%equilFromMacro( density  = mass_dens,             &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq                    )

    do iElem = 1, globBC%nElems(iLevel)
      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq( iDir+(iElem-1)*QQ )
        end if
      end do
    end do

  end subroutine spc_outlet_vel
! ****************************************************************************** !


! ****************************************************************************** !
  !> Outlet mixture pressure species equilibrium boundary
  !! kinematic pressure is computed from pressure
  !! species density and velocity are extrapolated
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_outlet_eq',
  !!    pressure = 'press0'
  !! }
  !!}
  !!variable = {
  !! {
  !!   label = 'press0',
  !!   ncomponents = 1,
  !!   vartype = 'st_fun',
  !!   st_fun = 1.0
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_outlet_eq( me, state, bcBuffer, globBC, levelDesc, tree,      &
    &                       nSize, iLevel, sim_time, neigh, layout, fieldProp, &
    &                       varPos, nScalars, varSys, derVarPos, physics,      &
    &                       iField, mixture                                    )
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
    ! fEq uses AOS layout
    real(kind=rk) :: fEq( layout%fStencil%QQ*globBC%nElems(iLevel) )
    !store pdf of all fields for derive function pointer input
    real(kind=rk) :: fTmp( layout%fStencil%QQ*globBC%nElems(iLevel) &
      &              * varSys%nStateVars )
    integer :: iDir              ! Direction counter
    integer :: iElem             ! Element counter
    integer :: nFields, iFieldLoc, pos, offset, QQ!, bcPress_pos
!    real(kind=rk) :: press(globBC%nElems(iLevel))
    integer :: elemPos, posInBuffer, neighPos
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    ! @todo KM: use pressure term to set species density at boundary
    ! get density value from the definition in the lua file
    ! position of boundary pressure in varSys
!KM!    bcPress_pos = me%bc_states%pressure%varPos
!KM!    ! get pressure variable from spacetime function
!KM!    call varSys%method%val(bcPress_pos)%get_valOfIndex( &
!KM!      & varSys  = varSys,                               &
!KM!      & time    = sim_time,                             &
!KM!      & iLevel  = iLevel,                               &
!KM!      & idx     = me%bc_states%pressure                 &
!KM!      &           %pntIndex%indexLvl(iLevel)            &
!KM!      &           %val(1:globBC%nElems(iLevel)),        &
!KM!      & nVals   = globBC%nElems(iLevel),                &
!KM!      & res     = press                                 )
!KM!
!KM!    ! If physical quantities are given, transform to lattice units by division
!KM!    ! with the conversion factor
!KM!    ! kinematic pressure = pressure/rho0
!KM!    press = (press/mixture%rho0)/(physics%fac( iLevel )%press/physics%rho0)

    ! Calculate the density of current element
    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      neighPos = me%neigh(iLevel)%posInState(1,iElem)
      do iFieldLoc = 1, nFields
        offset = (iElem-1)*nFields + iFieldLoc
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          !fTmp( ( ielem-1)* nscalars+ pos )   &
          !  & = bcBuffer(                                               &
          !  & ( posinbuffer-1)* nscalars+ pos )
          fTmp( ( ielem-1)* nscalars+ pos )     &
            & = state((neighpos-1)*nscalars+ idir+(ifieldloc-1)*qq)
        end do
      end do
    end do !iElem

    ! Extrapolate equilibrium from local
    ! Calculate the equilibrium distribution
    call derVarPos%equilFromState( state  = fTmp,                    &
      &                            iField = iField,                  &
      &                            nElems = globBC%nElems( iLevel ), &
      &                            varSys = varSys,                  &
      &                            layout = layout,                  &
      &                            res    = fEq                      )

    ! Run over all elements in list for this boundary
    do iElem = 1, globBC%nElems( iLevel )
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Now assign the values
          state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 ) &
            & = fEq( iDir+(iElem-1)*QQ )
        end if
      end do
    end do

  end subroutine spc_outlet_eq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Mole fraction boundary condition
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moleFrac',
  !!    moleFraction = 0.0
  !!     }
  !!```
  !! Post collision pdf of incoming link is updated with equilibrium functions
  !! which is similar to initial condition.
  !! \f$ \bar{f^c_k} = f^{eq}_k(\rho_k, u_k) + \frac{\lambda}{2}
  !!             (f^{eq}_k(\rho_k, u_k) + f^{eq}_k(\rho_k, u^{eq}_k})) \f$
  !! Here,
  !! $u_k$ - velocity of species k from original pdf computed from LSE
  !! $u^{eq}_k$ - equilibrium velocity of species k
  !!
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type::fnct]] function pointer.
  subroutine spc_moleFrac( me, state, bcBuffer, globBC, levelDesc, tree,      &
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
    real(kind=rk) :: moleFrac(globBC%nElems(iLevel))
    real(kind=rk) :: rho, press
    integer :: iElem, iDir, iFieldLoc, nFields, pos, QQ
    integer :: bcMoleFrac_pos, elemPos, posInBuffer
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*varSys%nStateVars )
    real(kind=rk) :: fEq, fEqStar, ucx, ucxStar, ucxQuad, ucxQuadStar
    real(kind=rk) :: usq, usqStar
    real(kind=rk) :: velAvg(3), velQuad(3), velQuadStar(3), eqVel(3)
    real(kind=rk) :: vel(3,varSys%nStateVars)
    real(kind=rk) :: mom(3,varSys%nStateVars)
    real(kind=rk) :: mass_dens(varSys%nStateVars), num_dens(varSys%nStateVars)
    real(kind=rk) :: totMassDens
    real(kind=rk) :: moleFrac_loc(varSys%nStateVars)
    real(kind=rk) :: molWeightInv(varSys%nStateVars), phi(varSys%nStateVars)
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: weight0_inv, paramBInv
    type(mus_varSys_data_type), pointer :: fPtr
    !---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
!    write(*,*) 'boundary label ', trim(me%label)
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    do iFieldLoc = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeightInv(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                           %fieldProp%species%molWeightInv
      ! molecular weight ratio
      phi(iFieldLoc) = fPtr%solverData%scheme &
        &                  %field(iFieldLoc)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFieldLoc, :) = fPtr%solverData%scheme          &
        &                            %field(iFieldLoc)%fieldProp &
        &                            %species%resi_coeff(:)
    end do
    weight0_inv = 1.0_rk/layout%weight(layout%fStencil%restPosition)

    paramBInv = 1.0_rk / mixture%paramB

    ! position of molefraction spacetime variable in varSys
    bcMoleFrac_pos = me%bc_states%moleFrac%varPos
    ! mole fraction
    call varSys%method%val(bcMoleFrac_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFrac                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFrac                                 )

    ! Calculate the density of current element
    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
      ! local state vector to compute density and velocity of other species
      do iFieldLoc = 1, varSys%nStateVars
        do iDir = 1, QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all( pos ) = bcBuffer(                                  &
            & ( posinbuffer-1)* nscalars+ pos)
        end do
      end do

      !KM: using inital mixture number density rho0/mixtureMOlWeight
      rho = mixture%moleDens0LB *  moleFrac(iElem) &
        & * fieldProp%species%molWeight
      press = cs2 * rho * fieldProp%species%molWeigRatio

      !local density and momentum
      mass_dens = 0.0_rk
      vel = 0.0_rk
      do iFieldLoc = 1, nFields
        mass_dens(iFieldLoc) = &
          & sum( fTmp_all( varSys%method%val(iFieldLoc)%state_varPos(:) ) )
        num_dens(iFieldLoc) = mass_dens(iFieldLoc) * molWeightInv(iFieldLoc)
        !velocity
        call derVarPos%momFromState( state  = fTmp_all,              &
          &                          iField = iFieldLoc,             &
          &                          nElems = 1,                     &
          &                          varSys = varSys,                &
          &                          layout = layout,                &
          &                          res    = mom(:, iFieldLoc)      )
      end do

      mass_dens(iField) = rho
      num_dens(iField) = rho / fieldProp%species%molWeight

      !mass flux
      do iFieldLoc = 1, nFields
        vel(:, iFieldLoc) = mom(:, iFieldLoc) / mass_dens(iFieldLoc)
      end do

      ! total mass density
      totMassDens = sum(mass_dens)

      ! mole fraction
      moleFrac_loc = num_dens/sum(num_dens)

      ! equilibrium velocity
      eqVel = vel(:, iField)
      do iFieldLoc = 1, nFields
        eqVel(:) = eqVel(:) + resi_coeff(iField, iFieldLoc)               &
          &      * phi(iField) * moleFrac_loc(iFieldLoc)                  &
          &      * ( vel(:, iFieldLoc) - vel( :, iField ) )               &
          &      / mixture%paramB
      end do

      ! mass averaged mixture velocity
      velAvg(1) = sum(mom(1,:)) / totMassDens
      velAvg(2) = sum(mom(2,:)) / totMassDens
      velAvg(3) = sum(mom(3,:)) / totMassDens

      velQuadStar(:) = mixture%theta_eq*velAvg(:)                       &
        &            + (1.0_rk-mixture%theta_eq) * eqVel(:)
      velQuad(:) = mixture%theta_eq*velAvg(:)                           &
        &            + (1.0_rk-mixture%theta_eq)*vel(:, iField)

      usqStar = dot_product(velQuadStar, velQuadStar)*t2cs2inv
      usq = dot_product(velQuad, velQuad)*t2cs2inv

      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir =1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ucx = dot_product(layout%fStencil%cxDir(:, iDir),             &
            &               vel(:,iField))
          ucxQuad = dot_product(layout%fStencil%cxDir(:, iDir),         &
            &               velQuad)

          ucxStar = dot_product(layout%fStencil%cxDir(:, iDir),         &
            &                   eqVel)
          ucxQuadStar = dot_product(layout%fStencil%cxDir(:, iDir),     &
            &                   velQuadStar)

          ! eqVel is actually is rho_i*eqVel so ucxStar is not multiplied
          ! with rho in below equation
          fEqStar = layout%weight(iDir) * rho                  &
            & * ( phi(iField) + ucxStar * cs2inv               &
            & + ucxQuadStar * ucxQuadStar * t2cs4inv - usqStar )

          fEq = layout%weight(iDir) * rho          &
            & * ( phi(iField) + ucx * cs2inv       &
            & + ucxQuad * ucxQuad * t2cs4inv - usq )

          if ( iDir == layout%fStencil%restPosition ) then
            ! equilibrium at rest
            fEq = layout%weight( iDir ) * rho * (                          &
              & ( weight0_inv + (1.0_rk-weight0_inv) * phi(iField) ) - usq )
            fEqStar = layout%weight( iDir ) * rho * (                        &
              & ( weight0_inv + (1.0_rk-weight0_inv) * phi(iField) ) - usqStar )
          end if

          ! set equilibrium
          ! setting transformed pdf in similar to initial condition i.e
          ! \bar{f^c} = f^{eq}(\rho_spc, u_spc) + \frac{\lambda}{2}
          !             (f^{eq}(\rho_spc, u_spc) + f^{eq}(\rho_spc, u^{eq}_spc}))
          state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq + mixture%omega_diff*0.5_rk*( fEq - fEqStar )
        end if
      end do! iDir
    end do

  end subroutine spc_moleFrac
! ****************************************************************************** !

! ****************************************************************************** !
  !> Mole fraction boundary condition
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_molefrac_eq',
  !!    moleFraction = 0.0
  !!     }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moleFrac_eq( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                         nSize, iLevel, sim_time, neigh, layout,         &
    &                         fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                         physics, iField, mixture                        )
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
    ! fEq uses AOS layout
    real(kind=rk) :: fEq( layout%fStencil%QQ*globBC%nElems(iLevel) )
    !store pdf of all fields for derive function pointer input
    real(kind=rk) :: fTmp( layout%fStencil%QQ*globBC%nElems(iLevel) &
      &              * varSys%nStateVars )
    real(kind=rk) :: mass_dens( globBC%nElems(iLevel)*varSys%nStateVars )
    real(kind=rk) :: moleFrac(globBC%nElems(iLevel))
    integer :: iElem, iDir, iFieldLoc, nFields, pos
    !> Velocity on boundary element
    real(kind=rk) :: uxB(globBC%nElems(iLevel)*varSys%nStateVars*3)
    real(kind=rk) :: velocity(3,globBC%nElems(iLevel)*varSys%nStateVars)
    integer :: offset, QQ, bcMoleFrac_pos, elemPos, posInBuffer
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    ! position of molefraction spacetime variable in varSys
    bcMoleFrac_pos = me%bc_states%moleFrac%varPos
    ! mole fraction
    call varSys%method%val(bcMoleFrac_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFrac                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFrac                                 )

    ! copy state
    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      do iFieldLoc = 1, nFields
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp( ( ielem-1)* nscalars+ pos )   &
            & = bcBuffer(                                               &
            & ( posinbuffer-1)* nscalars+ pos )
        end do
      end do !iField
    end do !iElem

    ! get density and velocity.
    ! if current field, use velocity and density defined in lua file.
    ! for other fields, derive density and velocity from state
    mass_dens = 0.0_rk
    do iFieldLoc = 1, nFields
      if (iFieldLoc == iField) then
        ! store velocity in input_loc array
        do iElem = 1, globBC%nElems(iLevel)
          offset = (iElem-1)*nFields + iFieldLoc
          ! compute current species mass density from specified molefraction at
          ! boundary
          ! rho = n_t * chi_i * m_i + massFrac_i * rho0 * KinePress/(cs2*phi_i)
          ! massFrac_i = rho_i/rho
          ! rho = n_t * chi_i * m_i + rho_i * KinePress/(cs2*phi_i)
          ! 1st term is zero order density
          ! 2nd term is second order density with kinematic mixture pressure,
          ! p = cs2*(sum(phi_k*rho_k) - min_a(m_a)*n0)/rho0
          !KM: using inital mixture number density rho0/mixtureMOlWeight
          !Using local tot_NuMdens increases density over time
          !kinePress = cs2 * ( dot_product(phi, mass_dens)                     &
          !  &       - minval(molWeight)*mixture%moleDens0LB )/mixture%rho0LB
          mass_dens( offset ) = mixture%moleDens0LB * moleFrac(iElem) &
            &                 * fieldProp%species%molWeight
        end do
      else
        do iElem = 1, globBC%nElems(iLevel)
          offset = (iElem-1)*nFields + iFieldLoc
          ! species density
          do iDir = 1, layout%fStencil%QQ
            pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
            mass_dens(offset) = mass_dens(offset) &
              & + fTmp(( ielem-1)* nscalars+ pos )
          end do !iDir
        end do !iElem
      end if
    end do !iField

    ! Derive all species velocities at once since
    ! it is inefficient to derive each species velocity
    call derVarPos%velocitiesFromState( state  = fTmp,                  &
      &                                 iField = iField,                &
      &                                 nElems = globBC%nElems(iLevel), &
      &                                 varSys = varSys,                &
      &                                 layout = layout,                &
      &                                 res    = uxB                    )

    do iElem = 1, globBC%nElems(iLevel)
      offset = (iElem-1)*nFields
      do iFieldLoc = 1, nFields
        velocity(:, offset+iFieldLoc) =      &
          & uxB(offset*3 + (iFieldLoc-1)*3 + 1 : &
          &     offset*3 + iFieldLoc*3 )
      end do !iField
    end do !iElem

    fEq = 0.0_rk
    ! Calculate the equilibrium distribution
    call derVarPos%equilFromMacro( density  = mass_dens,             &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq                    )

    do iElem = 1, globBC%nElems(iLevel)
      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq( iDir+(iElem-1)*QQ )
        end if
      end do
    end do

  end subroutine spc_moleFrac_eq
! ****************************************************************************** !

! ****************************************************************************** !
  !> Mole density boundary condition
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moledens_eq',
  !!    moleFraction = 0.0
  !!     }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moleDens_eq( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                         nSize, iLevel, sim_time, neigh, layout,         &
    &                         fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                         physics, iField, mixture                        )
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
    ! ---------------------------------------------------------------------------
    ! fEq uses AOS layout
    real(kind=rk) :: fEq( layout%fStencil%QQ*globBC%nElems(iLevel) )
    !store pdf of all fields for derive function pointer input
    real(kind=rk) :: fTmp( layout%fStencil%QQ*globBC%nElems(iLevel) &
      &              * varSys%nStateVars )
    real(kind=rk) :: mass_dens( globBC%nElems(iLevel)*varSys%nStateVars )
    real(kind=rk) :: moleDens(globBC%nElems(iLevel))
    integer :: iElem, iDir, iFieldLoc, nFields, pos
    !> Velocity on boundary element
    real(kind=rk) :: uxB(globBC%nElems(iLevel)*varSys%nStateVars*3)
    real(kind=rk) :: velocity(3,globBC%nElems(iLevel)*varSys%nStateVars)
    integer :: offset, QQ, bcMoleDens_pos, elemPos, posInBuffer
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    ! position of molefraction spacetime variable in varSys
    bcMoleDens_pos = me%bc_states%moleDens%varPos
    ! mole fraction
    call varSys%method%val(bcMoleDens_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFrac                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleDens                                 )

    ! Convert to lattice Unit
    moleDens = moleDens / physics%moleDens0

    ! copy state
    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      do iFieldLoc = 1, nFields
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp( ( ielem-1)* nscalars+ pos )   &
            & = bcBuffer(                                               &
            & ( posinbuffer-1)* nscalars+ pos )
        end do
      end do !iField
    end do !iElem

    ! get density and velocity.
    ! if current field, use velocity and density defined in lua file.
    ! for other fields, derive density and velocity from state
    mass_dens = 0.0_rk
    do iFieldLoc = 1, nFields
      if (iFieldLoc == iField) then
        ! store velocity in input_loc array
        do iElem = 1, globBC%nElems(iLevel)
          offset = (iElem-1)*nFields + iFieldLoc
          ! compute current species mass density from specified molefraction at
          ! boundary
          ! rho = n_t * chi_i * m_i + massFrac_i * rho0 * KinePress/(cs2*phi_i)
          ! massFrac_i = rho_i/rho
          ! rho = n_t * chi_i * m_i + rho_i * KinePress/(cs2*phi_i)
          ! 1st term is zero order density
          ! 2nd term is second order density with kinematic mixture pressure,
          ! p = cs2*(sum(phi_k*rho_k) - min_a(m_a)*n0)/rho0
          !KM: using inital mixture number density rho0/mixtureMOlWeight
          !Using local tot_NuMdens increases density over time
          !kinePress = cs2 * ( dot_product(phi, mass_dens)                     &
          !  &       - minval(molWeight)*mixture%moleDens0LB )/mixture%rho0LB
          mass_dens( offset ) = moleDens(iElem) * fieldProp%species%molWeight
        end do
      else
        do iElem = 1, globBC%nElems(iLevel)
          offset = (iElem-1)*nFields + iFieldLoc
          ! species density
          do iDir = 1, layout%fStencil%QQ
            pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
            mass_dens(offset) = mass_dens(offset) &
              & + fTmp(( ielem-1)* nscalars+ pos )
          end do !iDir
        end do !iElem
      end if
    end do !iField

    ! Derive all species velocities at once since
    ! it is inefficient to derive each species velocity
    call derVarPos%velocitiesFromState( state  = fTmp,                  &
      &                                 iField = iField,                &
      &                                 nElems = globBC%nElems(iLevel), &
      &                                 varSys = varSys,                &
      &                                 layout = layout,                &
      &                                 res    = uxB                    )

    do iElem = 1, globBC%nElems(iLevel)
      offset = (iElem-1)*nFields
      do iFieldLoc = 1, nFields
        velocity(:, offset+iFieldLoc) =      &
          & uxB(offset*3 + (iFieldLoc-1)*3 + 1 : &
          &     offset*3 + iFieldLoc*3 )
      end do !iField
    end do !iElem

    fEq = 0.0_rk
    ! Calculate the equilibrium distribution
    call derVarPos%equilFromMacro( density  = mass_dens,             &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq                    )

    do iElem = 1, globBC%nElems(iLevel)
      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq( iDir+(iElem-1)*QQ )
        end if
      end do
    end do

  end subroutine spc_moleDens_eq
! ****************************************************************************** !



! ****************************************************************************** !
  !> Mole fraction boundary condition with thermodynamic factor
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moleFrac_wtdf',
  !!    moleFraction = 0.0
  !!     }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moleFrac_wtdf( me, state, bcBuffer, globBC, levelDesc, tree, &
    &                           nSize, iLevel, sim_time, neigh, layout,       &
    &                           fieldProp, varPos, nScalars, varSys,          &
    &                           derVarPos, physics, iField, mixture           )
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
    real(kind=rk) :: fEq
    !store pdf of all fields for derive function pointer input
    real(kind=rk) :: fTmp( varSys%nScalars )
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: num_dens( varSys%nStateVars )
    real(kind=rk) :: tot_numDens, tot_massDens
    real(kind=rk) :: moleFrac(globBC%nElems(iLevel))
    integer :: iElem, iDir, iFieldLoc, iField_2, iField_3, nFields, pos, QQ
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk), dimension(varSys%nStateVars) :: moleFrac_loc
    real(kind=rk), dimension(3, varSys%nStateVars) :: first_moments, velocity
    real(kind=rk), &
      & dimension(varSys%nStateVars, varSys%nStateVars) :: matA, invA, &
      & resi_coeff, thermodynamic_fac, inv_thermodyn_fac, diff_coeff
    real(kind=rk) :: temp, press, phy_moleDens_fac
    real(kind=rk) :: omega, paramB, theta_eq
    real(kind=rk) :: eqVel(3), velAvg(3), velQuad(3), usqr, ucx, ucxQuadTerm
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: bcMoleFrac_pos, posInBuffer, posInState
    ! ------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    do iFieldLoc = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc)   &
        &                                   %fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc)      &
        &                             %fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFieldLoc, :) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                            %fieldProp%species%resi_coeff(:)
    end do

    !KM @todo check moleDens for multilevel
    phy_moleDens_fac = physics%moleDens0
    !fixed parameter
    paramB = mixture%paramB
    !equilibrium theta
    theta_eq = mixture%theta_eq

    omega = mixture%relaxLvl(iLevel)%omega_diff
    ! temperature
    temp = mixture%temp0
    ! atmospheric pressure
    press = mixture%atm_press
!write(*,*) 'iField ', iField
!write(*,*) 'resi_coeff ', resi_coeff
!write(*,*) 'molWeight moleFrac BC', molWeight
!write(*,*) 'phi ', phi
!stop

    ! position of molefraction spacetime variable in varSys
    bcMoleFrac_pos = me%bc_states%moleFrac%varPos
    ! mole fraction
    call varSys%method%val(bcMoleFrac_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFrac                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFrac                                 )

    ! Calculate the density of current element
    do iElem = 1, globBC%nElems(iLevel)
!write(*,*) 'IElem ', iElem
      mass_dens = 0.0_rk
      first_moments = 0.0_rk
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      do iFieldLoc = 1, varSys%nStateVars
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp( pos ) = bcBuffer(                                   &
            & ( posinbuffer-1)* nscalars+ pos )
          mass_dens(iFieldLoc) = mass_dens(iFieldLoc) + fTmp( pos )

          !field momentum (rho*u)
          first_moments( 1, iFieldLoc ) = first_moments( 1, iFieldLoc ) &
            & + fTmp( pos ) * layout%fStencil%cxDir( 1, iDir )

          first_moments( 2, iFieldLoc ) = first_moments( 2, iFieldLoc ) &
            & + fTmp( pos ) * layout%fStencil%cxDir( 2, iDir )

          first_moments( 3, iFieldLoc ) = first_moments( 3, iFieldLoc ) &
            & + fTmp( pos ) * layout%fStencil%cxDir( 3, iDir )

        end do
        num_dens(iFieldLoc) = mass_dens(iFieldLoc)/molWeight(iFieldLoc)
      end do

      ! update num_dens, massDens and moleFrac of current species with specified
      ! molefraction
      num_dens(iField) = moleFrac(iElem) * mixture%moleDens0LB
      mass_dens(iField) = num_dens(iField) * molWeight(iField)

      !total number density
      tot_NumDens = sum(num_dens)

      !total mass density
      tot_massDens = sum(mass_dens)

      !mole fraction
      moleFrac_loc(:) =  num_dens(:)/tot_NumDens
!write(*,*) 'num_dens ', num_dens
!write(*,*) 'num_dens*phy_moleDens_fac ', num_dens*phy_moleDens_fac

      ! MS-Diff coeff matrix from C++ code
      call mus_calc_MS_DiffMatrix( nFields, temp, press,                       &
        &                          num_dens*phy_moleDens_fac, diff_coeff )
!write(*,*) 'diff_coeff ', diff_coeff
      ! Convert to lattice unit
      resi_coeff = physics%fac(iLevel)%diffusivity/diff_coeff
!write(*,*) 'resi_coeff ', resi_coeff

      call mus_calc_thermFactor( nFields, temp, press, moleFrac_loc,           &
        &                        thermodynamic_fac )
!write(*,*) 'thermodyn_fac ', thermodynamic_fac

      inv_thermodyn_fac = invert_matrix( thermodynamic_fac )
!write(*,*) 'inv_thermodyn_fac ', inv_thermodyn_fac

      matA = 0.0_rk
      !build up matrix to solver LSE for actual velocity
      do iFieldLoc = 1, nFields
        !set diagonal part
        matA(iFieldLoc, iFieldLoc) = 1.0_rk
        do iField_2 = 1, nFields
          do iField_3 = 1, nFields
            matA(iFieldLoc, iField_2) = matA(iFieldLoc, iField_2)              &
              &                       + omega * 0.5_rk                         &
              &                       * inv_thermodyn_fac(iFieldLoc, iField_2) &
              &                       * resi_coeff(iField_2, iField_3)         &
              &                       * phi(iField_2) * moleFrac_loc(iField_3) &
              &                       / paramB
          end do
        end do
        !set non-diagonal part
        do iField_2 = 1, nFields
          do iField_3 = 1, nFields
            matA(iFieldLoc, iField_3) = matA(iFieldLoc, iField_3)              &
              &                       - omega * 0.5_rk                         &
              &                       * inv_thermodyn_fac(iFieldLoc, iField_2) &
              &                       * resi_coeff(iField_2, iField_3)         &
              &                       * phi(iField_3) * moleFrac_loc(iField_2) &
              &                       / paramB
          end do
        end do
      end do
!write(*,*) 'matA ', matA

      ! invert matrix
      invA = invert_matrix( matA )
!write(*,*) 'invA ', invA

      !actual velocity of all species
      velocity(1, :) = matmul( invA, first_moments(1,:) ) / mass_dens(:)
      velocity(2, :) = matmul( invA, first_moments(2,:) ) / mass_dens(:)
      velocity(3, :) = matmul( invA, first_moments(3,:) ) / mass_dens(:)
!write(*,*) 'velocity ', velocity

      ! equilibrium velocity with thermodynamic factor
      eqVel( : ) = mass_dens(iField)*velocity( :, iField )
      do iField_2 = 1, nFields
        do iField_3 = 1, nFields
          eqVel( : ) = eqVel( : )                                      &
            &        + inv_thermodyn_fac(iField, iField_2)             &
            &        * mass_dens(iField_2)                             &
            &        * resi_coeff( iField_2, iField_3 ) * phi(iField_2)&
            &        * moleFrac_loc(iField_3)                          &
            &        * (velocity(:, iField_3) - velocity(:,iField_2))  &
            &        / paramB
        end do
      end do
!write(*,*) 'eqVel ', eqVel
      !compute mass averaged mixture velocity
      velAvg(1) = dot_product( mass_dens, velocity(1,:) )/tot_massDens
      velAvg(2) = dot_product( mass_dens, velocity(2,:) )/tot_massDens
      velAvg(3) = dot_product( mass_dens, velocity(3,:) )/tot_massDens

      velQuad = theta_eq*velAvg                                              &
        &     + (1.0_rk - theta_eq) * eqVel(:) / mass_dens(iField)

      usqr = dot_product( velQuad, velQuad ) * t2cs2inv

      posInState = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
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

          ucx = dot_product(                                                   &
            & dble(layout%fStencil%cxDir(:, iDir)), eqVel(:) )

          ucxQuadTerm = dot_product(                                           &
            & dble(layout%fStencil%cxDir(:, iDir)), velQuad )

          ! Equilibrium
          feq = layout%weight(iDir) * ( mass_dens(iField) * ( phi(iField)      &
            & + ucxQuadTerm * ucxQuadTerm * t2cs4inv - usqr ) + ucx * cs2inv )

          if ( iDir == layout%fStencil%restPosition ) then
          ! equilibrium at rest
            select case( trim(layout%fStencil%label) )
            case('d2q9')
              feq = layout%weight( iDir ) * mass_dens(iField) * ( &
                    & ( 9._rk - 5._rk * phi(iField) )/4._rk - usqr )
            case('d3q19')
              feq = layout%weight( iDir ) * mass_dens(iField) * ( &
                    & ( 3._rk - 2._rk * phi(iField) ) - usqr )
            end select
          end if

          ! set equilibrium
          state(                                                               &
            &  neigh (( idir-1)* nsize+ posinstate)+( ifield-1)* qq+ nscalars*0)&
            & = fEq
!write(*,*) 'fEq ', fEq
        end if
      end do
    end do

  end subroutine spc_moleFrac_wtdf
! ****************************************************************************** !


! ****************************************************************************** !
  !> molar flux equilibrium boundary condition
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moleflux',
  !!    moleflux = 0.0
  !!     }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moleFlux_eq( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                         nSize, iLevel, sim_time, neigh, layout,         &
    &                         fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                         physics, iField, mixture                        )
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
    !store pdf of all fields for derive function pointer input
    real(kind=rk) :: fTmp( layout%fStencil%QQ*globBC%nElems(iLevel) &
      &              * varSys%nStateVars )
    real(kind=rk) :: moleFlux(globBC%nElems(iLevel)*3), inv_flux
    !> Velocity on boundary element
    real(kind=rk) :: uxB(globBC%nElems(iLevel)*varSys%nStateVars*3)
    real(kind=rk) :: mass_dens( globBC%nElems(iLevel)*varSys%nStateVars )
    real(kind=rk) :: velocity(3,globBC%nElems(iLevel)*varSys%nStateVars)
    real(kind=rk) :: num_dens
    integer :: iElem, iDir, iFieldLoc, nFields, pos, QQ
    integer :: offset, posInBuffer, bcMoleFlux_pos, elemPos
    ! ------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars
    inv_flux = 1.0_rk / physics%fac( iLevel )%moleFlux

    ! position of boundary moleflux in varSys
    bcMoleFlux_pos = me%bc_states%moleFlux%varPos
    ! Get moleflux
    call varSys%method%val(bcMoleFlux_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFlux                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFlux                                 )

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
    moleFlux = moleFlux * inv_flux

    ! copy state
    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      do iFieldLoc = 1, nFields
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp( ( ielem-1)* nscalars+ pos ) = &
            & bcBuffer(                                                 &
            & ( posinbuffer-1)* nscalars+ pos )
        end do
      end do !iField
    end do !iElem

    ! get density and velocity.
    ! if current field, use velocity and density defined in lua file.
    ! for other fields, derive density and velocity from state
    ! species density
    mass_dens = 0.0_rk
    do iElem = 1, globBC%nElems(iLevel)
      do iFieldLoc = 1, nFields
        offset = (iElem-1)*nFields + iFieldLoc
        do iDir = 1, layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          mass_dens(offset) = mass_dens(offset) &
            & + fTmp(( ielem-1)* nscalars+ pos )
        end do !iDir
      end do !iField
    end do !iElem

    ! Derive all species velocities at once since
    ! it is inefficient to derive each species velocity
    call derVarPos%velocitiesFromState( state  = fTmp,                  &
      &                                 iField = iField,                &
      &                                 nElems = globBC%nElems(iLevel), &
      &                                 varSys = varSys,                &
      &                                 layout = layout,                &
      &                                 res    = uxB                    )

    do iFieldLoc = 1, nFields
      if (iFieldLoc == iField) then
        ! store velocity in input_loc array
        do iElem = 1, globBC%nElems(iLevel)
          offset = (iElem-1)*nFields + iFieldLoc
          num_dens = mass_dens(offset) * fieldProp%species%molWeightInv
          velocity(:,offset) = moleFlux((iElem-1)*3+1:iElem*3) / num_dens
        end do
      else
        do iElem = 1, globBC%nElems(iLevel)
          offset = (iElem-1)*nFields + iFieldLoc
          velocity(:,offset) =                             &
            & uxB((iElem-1)*nFields*3 + (iFieldLoc-1)*3 + 1 : &
            &     (iElem-1)*nFields*3 + iFieldLoc*3 )
        end do !iElem
      end if
    end do !iField

    ! Calculate the equilibrium distribution
    call derVarPos%equilFromMacro( density  = mass_dens,             &
      &                            velocity = velocity,              &
      &                            iField   = iField,                &
      &                            nElems   = globBC%nElems(iLevel), &
      &                            varSys   = varSys,                &
      &                            layout   = layout,                &
      &                            res      = fEq                    )

    do iElem = 1, globBC%nElems(iLevel)
      if( .not. btest( levelDesc%property(                        &
        &              globBC%elemLvl(iLevel)%elem%val(iElem)), prp_solid))then
      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
            & = fEq( iDir+(iElem-1)*QQ )
        end if
      end do
      end if
    end do

  end subroutine spc_moleFlux_eq
! ****************************************************************************** !


! ****************************************************************************** !
  !> molar flux boundary condition like velocity bounce back bc type
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moleflux',
  !!    mole_flux = 'mole_flux'
  !! }
  !!}
  !!variable = {
  !!  name = 'mole_flux',
  !!  ncomponents = 3,
  !!  vartype = 'st_fun',
  !!  st_fun = {0.06, 0.0, 0.0}
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moleFlux( me, state, bcBuffer, globBC, levelDesc, tree,      &
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
    ! ------------------------------------------------------------------------
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: molWeight
    real(kind=rk) :: massFlux(3)
    real(kind=rk) :: moleFlux(globBC%nElems(iLevel)*3), inv_flux
    integer :: iElem, iDir, QQ, posInBuffer, offset, bcMoleFlux_pos
    integer :: elemPos
    ! ------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    molWeight = fieldProp%species%molWeight
    inv_flux = 1.0_rk / physics%fac( iLevel )%moleFlux

    ! position of boundary velocity in varSys
    bcMoleFlux_pos = me%bc_states%moleFlux%varPos
    ! Get moleflux
    call varSys%method%val(bcMoleFlux_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFlux                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFlux                                 )

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
    moleFlux = moleFlux * inv_flux
    !moleFlux = moleFlux / physics%fac(iLevel)%flux

    ! Calculate the density of current element
    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      do iDir = 1,layout%fStencil%QQ
        fTmp( iDir ) = bcBuffer(                    &
       & ( posinbuffer-1)* nscalars+ varpos(idir) )
      end do
      !caulate mass flux from moleFlux
      ! massFlux = moleflux * molWeight
      offset = (iElem-1)*3
      massFlux(1) = moleFlux(offset+1) * molWeight
      massFlux(2) = moleFlux(offset+2) * molWeight
      massFlux(3) = moleFlux(offset+3) * molWeight
      !write(dbgUnit(1),*) 'massFlux ', massFlux
      !write(*,*) iElem, 'mass_flux ', massFlux

      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
              & = fTmp(layout%fStencil%cxDirInv( iDir ))            &
              & + layout%weight( iDir )*2._rk*cs2inv                &
              &    * ( layout%fStencil%cxDir( 1, iDir )*massFlux(1) &
              &    +   layout%fStencil%cxDir( 2, iDir )*massFlux(2) &
              &    +   layout%fStencil%cxDir( 3, iDir )*massFlux(3) )
        end if
      end do
    end do

  end subroutine spc_moleFlux
! ****************************************************************************** !

! **************************************************************************** !
  !> Inflow boundary condition for solvent based on non-Equilbrium
  !! extrapolation method. Similar to spc_velocity_noneq_expol except
  !! the mass density for solvent is enforced such that total moleDens0 is
  !! maintained.
  !! Default qVal=1.0.
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'inlet',
  !!    kind = 'spc_solvent_inflow',
  !!    velocity = {0.1,0.0,0.0},
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_solvent_inflow( me, state, bcBuffer, globBC, levelDesc, tree, &
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
    real(kind=rk) :: fTmpAll_1f(nScalars), fTmpAll_2f(nScalars)
    real(kind=rk) :: fEq_b, fEq_1f
    real(kind=rk) :: velAvg_1f(3), vel_b(3), eqVel_1f(3)
    real(kind=rk) :: usq_1f, ucx_1f, usq_b, ucx_b, ucxQuad_1f
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens_1f, mass_dens_2f
    real(kind=rk), dimension(varSys%nStateVars) :: num_dens_1f, num_dens_2f
    real(kind=rk) :: resi_coeff(varSys%nStateVars, varSys%nStateVars)
    real(kind=rk) :: mass_dens_b, num_dens_b, tot_NumDens_b
    real(kind=rk) :: omega_fac, paramBInv
    integer :: iElem, iDir, QQ, iFld, vPos
    integer :: bcVel_pos
    integer :: elemPos, elemOff, posInBuffer, neighPos
    real(kind=rk) :: vel_w(3*globBC%nElems(iLevel)), inv_vel
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: nFields
    ! -------------------------------------------------------------------------
    nFields = varSys%nStateVars
    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! get velocity_phy on boundary (surface)
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_w                               )

    ! convert physical velocity into LB velocity
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel
    vel_w = vel_w * inv_vel

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    ! resitivity coefficient of all species
    do iFld = 1, nFields
      !resi_coeff(iFld, :) = species(iFld)%resi_coeff(:)
      resi_coeff(iFld, :) = fPtr%solverData%scheme                          &
        &                                  %field(iFld)%fieldProp           &
        &                                              %species%resi_coeff(:)
    end do
    ! omega factor of LSE
    omega_fac = mixture%omega_diff * 0.5_rk
    paramBInv = 1.0_rk / mixture%paramB

    associate( auxField => fPtr%solverData%scheme%auxField(iLevel)%val,      &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species, &
      &        neighBufferPre_nNext => me%neigh(iLevel)%neighBufferPre_nNext,            &
      &        stencil => layout%fStencil                                    )

      QQ = layout%fStencil%QQ

      ! Calculate the density of current element
      do iElem = 1, globBC%nElems(iLevel)

        ! Position of current boundary element in bcBuffer
        posInBuffer = globBC%elemLvl(iLevel)%posInBcElemBuf%val( iElem )
        elemOff = (posInBuffer-1)*nScalars
        ! Position of 2nd fluid neighbor in state array
        neighPos = me%neigh(iLevel)%posInState(1, iElem)
        mass_dens_2f = 0.0_rk
        do iFld = 1, nFields
          do iDir = 1, QQ
            ! Position of current field variable current dir in state array
            vPos = varSys%method%val(iFld)%state_varPos(iDir)
            ! State array of 1st fluid element
            fTmpAll_1f(vPos) = bcBuffer( elemOff + vPos )
            ! Pre-Collision values of second fluid element
            fTmpAll_2f(vPos) = state(                                       &
              &  neigh((idir-1)* nsize+ neighpos)+( ifld-1)* qq+ nscalars*0 )
            ! mass density of second neighbor
            mass_dens_2f(iFld) = mass_dens_2f(iFld) + fTmpAll_2f(vPos)
          end do
          num_dens_2f(iFld) = mass_dens_2f(iFld) * species(iFld)%molWeightInv
        end do

        ! calculate spc density and velAvg from PDF of first fluid
        call calcDensAndVelsFromPDF( nFields    = nFields,         &
          &                          iField     = iField,          &
          &                          mass_dens  = mass_dens_1f,    &
          &                          num_dens   = num_dens_1f,     &
          &                          velAvg     = velAvg_1f,       &
          &                          eqVel      = eqVel_1f,        &
          &                          varSys     = varSys,          &
          &                          pdf        = fTmpAll_1f,      &
          &                          stencil    = layout%fStencil, &
          &                          species    = species,         &
          &                          resi_coeff = resi_coeff,      &
          &                          omega_fac  = omega_fac,       &
          &                          paramBInv  = paramBInv        )

        ! fluid node velocity square term of equilibrium function
        usq_1f = dot_product( velAvg_1f(1:3), velAvg_1f(1:3) ) * t2cs2inv

        ! Second-order extrapolation of number density of current species
        num_dens_b = ( 4.0_rk * num_dens_1f(iField) - num_dens_2f(iField) ) &
          &         / 3.0_rk
        ! Second-order extrapolation of numerical total number density
        tot_NumDens_b = ( 4.0_rk * sum(num_dens_1f) - sum(num_dens_2f) ) &
          &           / 3.0_rk
        ! Compute number density of current species at boundary element
        num_dens_b = mixture%moleDens0LB - ( tot_NumDens_b - num_dens_b )
        ! Mass density of current species at boundary element
        mass_dens_b = num_dens_b * species(iField)%molWeight
        !mass_dens_b = mixture%rho0LB - ( sum(mass_dens_1f) - mass_dens_1f(iField) )

        ! velocity at boundary node for qVal = 1.0
        vel_b = vel_w( (iElem-1)*3 + 1: iElem*3 )

        ! boundary node velocity square term of equilibrium function
        usq_b = dot_product( vel_b(1:3), vel_b(1:3) ) * t2cs2inv

        ! Position of this boundary element in state array
        elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)

        do iDir = 1, layout%fStencil%QQN
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            ! Compute equilibrium on fluid node
            ucx_1f = dot_product( stencil%cxDirRK(1:3, iDir), velAvg_1f )
            ucxQuad_1f = dot_product( stencil%cxDirRK(1:3, iDir), velAvg_1f )

            fEq_1f = layout%weight(iDir) * mass_dens_1f(iField)      &
              &   * ( species(iField)%molWeigRatio + ucx_1f * cs2inv &
              &   + ucxQuad_1f * ucxQuad_1f * t2cs4inv - usq_1f      )

            ! Compute equilibrium on boundary node
            ucx_b = dot_product( stencil%cxDirRK(1:3, iDir), vel_b )

            fEq_b = layout%weight(iDir) * mass_dens_b               &
              &   * ( species(iField)%molWeigRatio + ucx_b * cs2inv &
              &   + ucx_b * ucx_b * t2cs4inv - usq_b                )

            state( neigh((idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)  &
              & = fEq_b + ( fTmpAll_1f(varPos(iDir)) - fEq_1f )

          end if
        end do
      end do
    end associate

  end subroutine spc_solvent_inflow
! **************************************************************************** !


! **************************************************************************** !
  !> Velocity boundary condition based on non-Equilbrium extrapolation method.
  !! Default qVal=1.0.
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'inlet',
  !!    kind = 'spc_velocity_noneq_expol',
  !!    velocity = {0.1,0.0,0.0},
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_velocity_noneq_expol( me, state, bcBuffer, globBC, levelDesc, &
    &                                  tree, nSize, iLevel, sim_time, neigh,   &
    &                                  layout, fieldProp, varPos, nScalars,    &
    &                                  varSys, derVarPos, physics, iField,     &
    &                                  mixture                                 )
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
    real(kind=rk) :: fTmpAll_1f(nScalars), fTmpAll_2f(nScalars)
    real(kind=rk) :: fEq_b, fEq_1f
    real(kind=rk) :: velAvg_1f(3), vel_b(3), eqVel_1f(3)
    real(kind=rk) :: usq_1f, ucx_1f, usq_b, ucx_b
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens_1f, mass_dens_2f, &
      &                                            num_dens_1f
    real(kind=rk) :: resi_coeff(varSys%nStateVars, varSys%nStateVars)
    real(kind=rk) :: mass_dens_b
    real(kind=rk) :: omega_fac, paramBInv
    integer :: iElem, iDir, QQ, iFld, vPos
    integer :: bcVel_pos
    integer :: elemPos, elemOff, posInBuffer, neighPos
    real(kind=rk) :: vel_w(3*globBC%nElems(iLevel)), inv_vel
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: nFields
    ! -------------------------------------------------------------------------
    nFields = varSys%nStateVars
    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! get velocity_phy on boundary (surface)
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_w                               )

    ! convert physical velocity into LB velocity
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel
    vel_w = vel_w * inv_vel

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    ! resitivity coefficient of all species
    do iFld = 1, nFields
      !resi_coeff(iFld, :) = species(iFld)%resi_coeff(:)
      resi_coeff(iFld, :) = fPtr%solverData%scheme                          &
        &                                  %field(iFld)%fieldProp           &
        &                                              %species%resi_coeff(:)
    end do
    ! omega factor of LSE
    omega_fac = mixture%omega_diff * 0.5_rk
    paramBInv = 1.0_rk / mixture%paramB

    associate( auxField => fPtr%solverData%scheme%auxField(iLevel)%val,      &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species, &
      &        neighBufferPre_nNext => me%neigh(iLevel)%neighBufferPre_nNext,            &
      &        stencil => layout%fStencil                                    )

      QQ = layout%fStencil%QQ

      ! Calculate the density of current element
      do iElem = 1, globBC%nElems(iLevel)

        ! Position of current boundary element in bcBuffer
        posInBuffer = globBC%elemLvl(iLevel)%posInBcElemBuf%val( iElem )
        elemOff = (posInBuffer-1)*nScalars
        ! Position of 2nd fluid neighbor in state array
        neighPos = me%neigh(iLevel)%posInState(1, iElem)
        mass_dens_2f = 0.0_rk
        do iFld = 1, nFields
          do iDir = 1, QQ
            ! Position of current field variable current dir in state array
            vPos = varSys%method%val(iFld)%state_varPos(iDir)
            ! State array of 1st fluid element
            fTmpAll_1f(vPos) = bcBuffer( elemOff + vPos )
            ! Pre-Collision values of second fluid element
            fTmpAll_2f(vPos) = state(                                       &
              &  neigh((idir-1)* nsize+ neighpos)+( ifld-1)* qq+ nscalars*0 )
            ! mass density of second neighbor
            mass_dens_2f(iFld) = mass_dens_2f(iFld) + fTmpAll_2f(vPos)
          end do
        end do

        ! calculate spc density and velAvg from PDF of first fluid
        call calcDensAndVelsFromPDF( nFields    = nFields,         &
          &                          iField     = iField,          &
          &                          mass_dens  = mass_dens_1f,    &
          &                          num_dens   = num_dens_1f,     &
          &                          velAvg     = velAvg_1f,       &
          &                          eqVel      = eqVel_1f,        &
          &                          varSys     = varSys,          &
          &                          pdf        = fTmpAll_1f,      &
          &                          stencil    = layout%fStencil, &
          &                          species    = species,         &
          &                          resi_coeff = resi_coeff,      &
          &                          omega_fac  = omega_fac,       &
          &                          paramBInv  = paramBInv        )

        ! fluid node velocity square term of equilibrium function
        usq_1f = dot_product( velAvg_1f(1:3), velAvg_1f(1:3) ) * t2cs2inv

        ! density of current species at boundary node extrapolated from
        ! fluid elements
        mass_dens_b = ( 4.0_rk * mass_dens_1f(iField) - mass_dens_2f(iField) ) &
          &         / 3.0_rk

        ! velocity at boundary node for qVal = 1.0
        vel_b = vel_w( (iElem-1)*3 + 1: iElem*3 )

        ! boundary node velocity square term of equilibrium function
        usq_b = dot_product( vel_b(1:3), vel_b(1:3) ) * t2cs2inv

        ! Position of this boundary element in state array
        elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)

        do iDir = 1, layout%fStencil%QQN
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            ! Compute equilibrium on fluid node
            ucx_1f = dot_product( stencil%cxDirRK(1:3, iDir), velAvg_1f )

            fEq_1f = layout%weight(iDir) * mass_dens_1f(iField)      &
              &   * ( species(iField)%molWeigRatio + ucx_1f * cs2inv &
              &   + ucx_1f * ucx_1f * t2cs4inv - usq_1f              )

            ! Compute equilibrium on boundary node
            ucx_b = dot_product( stencil%cxDirRK(1:3, iDir), vel_b )

            fEq_b = layout%weight(iDir) * mass_dens_b               &
              &   * ( species(iField)%molWeigRatio + ucx_b * cs2inv &
              &   + ucx_b * ucx_b * t2cs4inv - usq_b                )

            state( neigh((idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)  &
              & = fEq_b + ( fTmpAll_1f(varPos(iDir)) - fEq_1f )

          end if
        end do
      end do
    end associate

  end subroutine spc_velocity_noneq_expol
! **************************************************************************** !


! **************************************************************************** !
  !> Mole fraction boundary condition for nonequilibrium extrapolation based.
  !! Default qVal=0.0.
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'spc_mole_fraction_noneq_expol',
  !!    kind = 'mole_fraction',
  !!    mole_fraction = 0.1
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_mole_fraction_noneq_expol( me, state, bcBuffer, globBC, &
    &          levelDesc, tree, nSize, iLevel, sim_time, neigh, layout,  &
    &          fieldProp, varPos, nScalars, varSys, derVarPos, physics,  &
    &          iField, mixture                                           )
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
    real(kind=rk) :: fTmpAll_2f(nScalars)
    real(kind=rk) :: fEq_b, fEq_2f
    real(kind=rk) :: velAvg_2f(3), vel_b(3), eqVel_2f(3)
    real(kind=rk) :: usq_b, ucx_b, usq_2f, ucx_2f
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens_2f
    real(kind=rk), dimension(varSys%nStateVars) :: num_dens_2f
    real(kind=rk) :: resi_coeff(varSys%nStateVars, varSys%nStateVars)
    real(kind=rk) :: omega_fac, paramBInv
    integer :: iElem, iDir, QQ, iFld, vPos
    integer :: bcMoleFrac_pos
    integer :: elemPos, elemOff, neighPos, posInBuffer
    real(kind=rk) :: moleFrac_w(globBC%nElems(iLevel)), mass_dens_b
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: nFields
    !> Inverse of lattice weight ar restPosition
    real(kind=rk) :: weight0Inv
    ! -------------------------------------------------------------------------
    weight0Inv = 1.0_rk / layout%weight(layout%fStencil%restPosition)
    nFields = varSys%nStateVars

    ! position of molefraction spacetime variable in varSys
    bcMoleFrac_pos = me%bc_states%moleFrac%varPos
    ! mole fraction
    call varSys%method%val(bcMoleFrac_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFrac                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFrac_w                               )

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )

    ! resitivity coefficient of all species
    do iFld = 1, nFields
      !resi_coeff(iFld, :) = species(iFld)%resi_coeff(:)
      resi_coeff(iFld, :) = fPtr%solverData%scheme                          &
        &                                  %field(iFld)%fieldProp           &
        &                                              %species%resi_coeff(:)
    end do
    ! omega factor of LSE
    omega_fac = mixture%omega_diff * 0.5_rk
    paramBInv = 1.0_rk / mixture%paramB

    associate( auxField => fPtr%solverData%scheme%auxField(iLevel)%val,      &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species, &
      &        stencil => layout%fStencil                                    )

      QQ = layout%fStencil%QQ

      ! Calculate the density of current element
      do iElem = 1, globBC%nElems(iLevel)

        ! Position of current boundary element in bcBuffer
        posInBuffer = globBC%elemLvl(iLevel)%posInBcElemBuf%val( iElem )
        elemOff = (posInBuffer-1)*nScalars
        ! Position of neighbor element in state array
        neighPos = me%neigh(iLevel)%posInState(1, iElem)

        do iFld = 1, nFields
          do iDir = 1, QQ
            ! Position of current field variable current dir in state array
            vPos = varSys%method%val(iFld)%state_varPos(iDir)
            ! Pre-Collision values of neighbor fluid element
            fTmpAll_2f(vPos) = state(                                    &
              &  neigh((idir-1)* nsize+ neighpos)+( ifld-1)* qq+ nscalars*0 )
            ! State array of 1st fluid element
            !fTmpAll_1f(vPos) = bcBuffer( elemOff + vPos )
          end do
        end do

        ! calculate spc density and velAvg from PDF of second fluid element
        call calcDensAndVelsFromPDF( nFields    = nFields,         &
          &                          iField     = iField,          &
          &                          mass_dens  = mass_dens_2f,    &
          &                          num_dens   = num_dens_2f,     &
          &                          velAvg     = velAvg_2f,       &
          &                          eqVel      = eqVel_2f,        &
          &                          varSys     = varSys,          &
          &                          pdf        = fTmpAll_2f,      &
          &                          stencil    = layout%fStencil, &
          &                          species    = species,         &
          &                          resi_coeff = resi_coeff,      &
          &                          omega_fac  = omega_fac,       &
          &                          paramBInv  = paramBInv        )

        ! fluid node velocity square term of equilibrium function
        usq_2f = dot_product( velAvg_2f(1:3), velAvg_2f(1:3) ) * t2cs2inv

        ! density of current species at boundary node
        mass_dens_b = mixture%moleDens0LB * moleFrac_w(iElem) &
          &         * species(iField)%molWeight

        ! velocity at boundary node extrapolated from 2nd fluid element
        vel_b = velAvg_2f

        ! boundary node velocity square term of equilibrium function
        usq_b = dot_product( vel_b(1:3), vel_b(1:3) ) * t2cs2inv

        ! Position of this boundary element in state array
        elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)

        do iDir = 1, layout%fStencil%QQN
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            ! Compute equilibrium on fluid node
            ucx_2f = dot_product( stencil%cxDirRK(1:3, iDir), velAvg_2f )

            fEq_2f = layout%weight(iDir) * mass_dens_2f(iField)      &
              &   * ( species(iField)%molWeigRatio + ucx_2f * cs2inv &
              &   + ucx_2f * ucx_2f * t2cs4inv - usq_2f              )

            ! Compute equilibrium on boundary node
            ucx_b = dot_product( stencil%cxDirRK(1:3, iDir), vel_b )

            fEq_b = layout%weight(iDir) * mass_dens_b               &
              &   * ( species(iField)%molWeigRatio + ucx_b * cs2inv &
              &   + ucx_b * ucx_b * t2cs4inv - usq_b                )

            ! update incomping link
            state( neigh((idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)  &
              & = fEq_b + ( fTmpAll_2f(varPos(iDir)) - fEq_2f )
          end if
        end do
      end do
    end associate

  end subroutine spc_mole_fraction_noneq_expol
! **************************************************************************** !

! **************************************************************************** !
  !> Open outflow boundary condition based on nonequilibrium extrapolation
  !! method. Default qVal = 0.0
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_outflow',
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_outflow( me, state, bcBuffer, globBC, levelDesc, tree, nSize, &
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
    real(kind=rk) :: fTmpAll_2f(nScalars), fTmpAll_3f(nScalars)
    real(kind=rk) :: fEq_b, fEq_2f
    real(kind=rk) :: velAvg_2f(3), vel_b(3), eqVel_2f(3)
    real(kind=rk) :: usq_2f, ucx_2f, usq_b, ucx_b
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens_2f, mass_dens_3f
    real(kind=rk), dimension(varSys%nStateVars) :: num_dens_2f
    real(kind=rk) :: resi_coeff(varSys%nStateVars, varSys%nStateVars)
    real(kind=rk) :: omega_fac, paramBInv, mass_dens_b
    integer :: iElem, iDir, QQ, iFld, vPos
    integer :: elemPos, neighPos_2f, neighPos_3f
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: nFields
    ! -------------------------------------------------------------------------
    nFields = varSys%nStateVars

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    ! resitivity coefficient of all species
    do iFld = 1, nFields
      !resi_coeff(iFld, :) = species(iFld)%resi_coeff(:)
      resi_coeff(iFld, :) = fPtr%solverData%scheme                          &
        &                                  %field(iFld)%fieldProp           &
        &                                              %species%resi_coeff(:)
    end do
    ! omega factor of LSE
    omega_fac = mixture%omega_diff * 0.5_rk
    paramBInv = 1.0_rk / mixture%paramB

    associate( auxField => fPtr%solverData%scheme%auxField(iLevel)%val,      &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species, &
      &        neighBufferPre_nNext => me%neigh(iLevel)%neighBufferPre_nNext,            &
      &        stencil => layout%fStencil                                    )

      QQ = layout%fStencil%QQ

      ! Calculate the density of current element
      do iElem = 1, globBC%nElems(iLevel)

        ! Position of second neighbor element in state array
        neighPos_2f = me%neigh(iLevel)%posInState(1, iElem)
        ! Position of third neighbor element in state array
        neighPos_3f = me%neigh(iLevel)%posInState(2, iElem)
        ! Compute only mass density from third neighbor
        mass_dens_3f = 0.0_rk
        do iFld = 1, nFields
          do iDir = 1, QQ
            ! Position of current field variable current dir in state array
            vPos = varSys%method%val(iFld)%state_varPos(iDir)
            ! Pre-Collision values of second fluid element
            fTmpAll_2f(vPos) = state(                                        &
              &  neigh((idir-1)* nsize+ neighpos_2f)+( ifld-1)* qq+ nscalars*0 )
            ! Pre-Collision values of third fluid element
            fTmpAll_3f(vPos) = state(                                        &
              &  neigh((idir-1)* nsize+ neighpos_3f)+( ifld-1)* qq+ nscalars*0 )
            ! mass density of third neighbor
            mass_dens_3f(iFld) = mass_dens_3f(iFld) + fTmpAll_3f(vPos)
          end do
        end do

        ! calculate spc density and velAvg from PDF of second fluid
        call calcDensAndVelsFromPDF( nFields    = nFields,         &
          &                          iField     = iField,          &
          &                          mass_dens  = mass_dens_2f,    &
          &                          num_dens   = num_dens_2f,     &
          &                          velAvg     = velAvg_2f,       &
          &                          eqVel      = eqVel_2f,        &
          &                          varSys     = varSys,          &
          &                          pdf        = fTmpAll_2f,      &
          &                          stencil    = layout%fStencil, &
          &                          species    = species,         &
          &                          resi_coeff = resi_coeff,      &
          &                          omega_fac  = omega_fac,       &
          &                          paramBInv  = paramBInv        )

        ! fluid node velocity square term of equilibrium function
        usq_2f = dot_product( velAvg_2f(1:3), velAvg_2f(1:3) ) * t2cs2inv

        ! second-order extrapolation of density of current species at boundary
        ! node from second and third fluid element
        mass_dens_b = ( 4.0_rk * mass_dens_2f(iField) - mass_dens_3f(iField) ) &
          &         / 3.0_rk

        ! first-order extrapolation of velocity at boundary node from second
        ! fluid element
        vel_b = velAvg_2f
        ! second-order extrapolation of velocity at boundary node
        ! vel_b =  ( 4.0_rk * velAvg_2f - velAvg_3f ) / 3.0_rk

        ! boundary node velocity square term of equilibrium function
        usq_b = dot_product( vel_b(1:3), vel_b(1:3) ) * t2cs2inv

        ! Position of this boundary element in state array
        elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)

        do iDir = 1, layout%fStencil%QQN
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            ! Compute equilibrium on second fluid node
            ucx_2f = dot_product( stencil%cxDirRK(1:3, iDir), velAvg_2f )

            fEq_2f = layout%weight(iDir) * mass_dens_2f(iField)      &
              &   * ( species(iField)%molWeigRatio + ucx_2f * cs2inv &
              &   + ucx_2f * ucx_2f * t2cs4inv - usq_2f              )

            ! Compute equilibrium on boundary node
            ucx_b = dot_product( stencil%cxDirRK(1:3, iDir), vel_b )

            fEq_b = layout%weight(iDir) * mass_dens_b               &
              &   * ( species(iField)%molWeigRatio + ucx_b * cs2inv &
              &   + ucx_b * ucx_b * t2cs4inv - usq_b                )

            ! update incomping direction
            state( neigh((idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)  &
              & = fEq_b + ( fTmpAll_2f(varPos(iDir)) - fEq_2f )
          end if
        end do
      end do
    end associate

  end subroutine spc_outflow
! **************************************************************************** !

! **************************************************************************** !
  !> Open outflow boundary condition for solvent based on nonequilibrium
  !! extrapolation. total moledens at boundary is enforced.
  !! method. Default qVal = 0.0
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_solvent_outflow',
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_solvent_outflow( me, state, bcBuffer, globBC, levelDesc,      &
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
    real(kind=rk) :: fTmpAll_2f(nScalars), fTmpAll_3f(nScalars)
    real(kind=rk) :: fEq_b, fEq_2f
    real(kind=rk) :: velAvg_2f(3), vel_b(3), eqVel_2f(3)
    real(kind=rk) :: usq_2f, ucx_2f, usq_b, ucx_b
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens_2f, mass_dens_3f
    real(kind=rk), dimension(varSys%nStateVars) :: num_dens_2f, num_dens_3f
    real(kind=rk) :: resi_coeff(varSys%nStateVars, varSys%nStateVars)
    real(kind=rk) :: omega_fac, paramBInv, mass_dens_b, num_dens_b, &
      &              tot_NumDens_b
    integer :: iElem, iDir, QQ, iFld, vPos
    integer :: elemPos, neighPos_2f, neighPos_3f
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: nFields
    ! -------------------------------------------------------------------------
    nFields = varSys%nStateVars

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    ! resitivity coefficient of all species
    do iFld = 1, nFields
      !resi_coeff(iFld, :) = species(iFld)%resi_coeff(:)
      resi_coeff(iFld, :) = fPtr%solverData%scheme                          &
        &                                  %field(iFld)%fieldProp           &
        &                                              %species%resi_coeff(:)
    end do
    ! omega factor of LSE
    omega_fac = mixture%omega_diff * 0.5_rk
    paramBInv = 1.0_rk / mixture%paramB

    associate( auxField => fPtr%solverData%scheme%auxField(iLevel)%val,      &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species, &
      &        neighBufferPre_nNext => me%neigh(iLevel)%neighBufferPre_nNext,            &
      &        stencil => layout%fStencil                                    )

      QQ = layout%fStencil%QQ

      ! Calculate the density of current element
      do iElem = 1, globBC%nElems(iLevel)

        ! Position of second neighbor element in state array
        neighPos_2f = me%neigh(iLevel)%posInState(1, iElem)
        ! Position of third neighbor element in state array
        neighPos_3f = me%neigh(iLevel)%posInState(2, iElem)
        ! Compute only mass density from third neighbor
        mass_dens_3f = 0.0_rk
        do iFld = 1, nFields
          do iDir = 1, QQ
            ! Position of current field variable current dir in state array
            vPos = varSys%method%val(iFld)%state_varPos(iDir)
            ! Pre-Collision values of second fluid element
            fTmpAll_2f(vPos) = state(                                       &
              &  neigh((idir-1)* nsize+ neighpos_2f)+( ifld-1)* qq+ nscalars*0 )
            ! Pre-Collision values of third fluid element
            fTmpAll_3f(vPos) = state(                                       &
              &  neigh((idir-1)* nsize+ neighpos_3f)+( ifld-1)* qq+ nscalars*0 )
            ! mass density of third neighbor
            mass_dens_3f(iFld) = mass_dens_3f(iFld) + fTmpAll_3f(vPos)
          end do
        end do

        ! calculate spc density and velAvg from PDF of second fluid
        call calcDensAndVelsFromPDF( nFields    = nFields,         &
          &                          iField     = iField,          &
          &                          mass_dens  = mass_dens_2f,    &
          &                          num_dens   = num_dens_2f,     &
          &                          velAvg     = velAvg_2f,       &
          &                          eqVel      = eqVel_2f,        &
          &                          varSys     = varSys,          &
          &                          pdf        = fTmpAll_2f,      &
          &                          stencil    = layout%fStencil, &
          &                          species    = species,         &
          &                          resi_coeff = resi_coeff,      &
          &                          omega_fac  = omega_fac,       &
          &                          paramBInv  = paramBInv        )

        ! fluid node velocity square term of equilibrium function
        usq_2f = dot_product( velAvg_2f(1:3), velAvg_2f(1:3) ) * t2cs2inv

        ! Second-order extrapolation of number density of current species
        num_dens_b = ( 4.0_rk * num_dens_2f(iField) - num_dens_3f(iField) ) &
          &         / 3.0_rk
        ! Second-order extrapolation of numerical total number density
        tot_NumDens_b = ( 4.0_rk * sum(num_dens_2f) - sum(num_dens_3f) ) &
          &           / 3.0_rk
        ! Compute number density of current species at boundary element
        num_dens_b = mixture%moleDens0LB - ( tot_NumDens_b - num_dens_b )
        ! Mass density of current species at boundary element
        mass_dens_b = num_dens_b * species(iField)%molWeight

        ! first-order extrapolation of velocity at boundary node from second
        ! fluid element
        vel_b = velAvg_2f
        ! second-order extrapolation of velocity at boundary node
        ! vel_b =  ( 4.0_rk * velAvg_2f - velAvg_3f ) / 3.0_rk

        ! boundary node velocity square term of equilibrium function
        usq_b = dot_product( vel_b(1:3), vel_b(1:3) ) * t2cs2inv

        ! Position of this boundary element in state array
        elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)

        do iDir = 1, layout%fStencil%QQN
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            ! Compute equilibrium on second fluid node
            ucx_2f = dot_product( stencil%cxDirRK(1:3, iDir), velAvg_2f )

            fEq_2f = layout%weight(iDir) * mass_dens_2f(iField)      &
              &   * ( species(iField)%molWeigRatio + ucx_2f * cs2inv &
              &   + ucx_2f * ucx_2f * t2cs4inv - usq_2f              )

            ! Compute equilibrium on boundary node
            ucx_b = dot_product( stencil%cxDirRK(1:3, iDir), vel_b )

            fEq_b = layout%weight(iDir) * mass_dens_b               &
              &   * ( species(iField)%molWeigRatio + ucx_b * cs2inv &
              &   + ucx_b * ucx_b * t2cs4inv - usq_b                )

            ! update incomping direction
            state( neigh((idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)  &
              & = fEq_b + ( fTmpAll_2f(varPos(iDir)) - fEq_2f )
          end if
        end do
      end do
    end associate

  end subroutine spc_solvent_outflow
! **************************************************************************** !

! ****************************************************************************** !
  !> Inflow boundary condition based on non-Equilbrium extrapolation method.
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_inflow',
  !!    mole_fraction = 0.01,
  !!    velocity = {0.1,0.0,0.0},
  !! }
  !!}
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_inflow( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
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
    real(kind=rk) :: fTmpAll_f(nScalars), fTmpAll_ff(nScalars)
    real(kind=rk) :: fEq_b, fEq_f, fEq_ff, fnEq_f, fnEq_ff, fnEq_b
    real(kind=rk) :: eqVel_f(3), velAvg_f(3), eqVel_ff(3), velAvg_ff(3), &
      &              vel_b(3)
    real(kind=rk) :: usq_f, ucx_f, ucxQuad_f, usq_b, ucx_b, usq_ff, ucx_ff, &
      &              ucxQuad_ff
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens_f, mass_dens_ff
    real(kind=rk), dimension(3, varSys%nStateVars) :: velocity_f, velocity_ff
    real(kind=rk) :: resi_coeff(varSys%nStateVars, varSys%nStateVars)
    real(kind=rk) :: omega_fac, paramBInv
    integer :: iElem, iDir, QQ, iFld, vPos, posInBuffer
    integer :: bcMoleFrac_pos, bcVel_pos
    integer :: elemPos, neighPos_1, neighPos_2
    real(kind=rk) :: vel_w(3*globBC%nElems(iLevel)), inv_vel
    real(kind=rk) :: moleFrac_w(globBC%nElems(iLevel)), mass_dens_b
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: coeff_w, coeff_f, coeff_ff
    integer :: nFields
    !> Inverse of lattice weight ar restPosition
    real(kind=rk) :: weight0Inv
    ! ------------------------------------------------------------------------
    weight0Inv = 1.0_rk / layout%weight(layout%fStencil%restPosition)
    nFields = varSys%nStateVars
    ! coefficient for extrapolation of macroscopic quantities for qVal=0.5
    coeff_w = 5.0_rk/3.0_rk
    coeff_f = -0.5_rk
    coeff_ff = -1.0_rk/6.0_rk
    ! position of molefraction spacetime variable in varSys
    bcMoleFrac_pos = me%bc_states%moleFrac%varPos
    ! mole fraction
    call varSys%method%val(bcMoleFrac_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFrac                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFrac_w                               )

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! get velocity_phy on boundary (surface)
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_w                               )

    ! convert physical velocity into LB velocity
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel
    vel_w = vel_w * inv_vel

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    ! resitivity coefficient of all species
    do iFld = 1, nFields
      !resi_coeff(iFld, :) = species(iFld)%resi_coeff(:)
      resi_coeff(iFld, :) = fPtr%solverData%scheme                          &
        &                                  %field(iFld)%fieldProp           &
        &                                              %species%resi_coeff(:)
    end do
    ! omega factor of LSE
    omega_fac = mixture%omega_diff * 0.5_rk
    paramBInv = 1.0_rk / mixture%paramB

    associate( auxField => fPtr%solverData%scheme%auxField(iLevel)%val,      &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species, &
      &        stencil => layout%fStencil                                    )

      QQ = layout%fStencil%QQ

!write(dbgUnit(1),*) 'globBC%nElems ', globBC%nElems(ilevel)
!write(dbgUnit(1),*) 'iField ', iField
      ! Calculate the density of current element
      do iElem = 1, globBC%nElems(iLevel)

        ! Position of current boundary element in bcBuffer
        posInBuffer = globBC%elemLvl(iLevel)%posInBcElemBuf%val( iElem )
        ! Position of neighbor element in state array
        neighPos_1 = me%neigh(iLevel)%posInState(1, iElem)
        neighPos_2 = me%neigh(iLevel)%posInState(2, iElem)

        do iFld = 1, nFields
          do iDir = 1, QQ
            ! Position of current field variable current dir in state array
            vPos = varSys%method%val(iFld)%state_varPos(iDir)
            ! State array of current element
            fTmpAll_f(vPos) = bcBuffer( (posInBuffer-1)*nScalars + vPos )
            !fTmpAll_f(vPos) = state(                                    &
            !  &  neigh((idir-1)* nsize+ neighpos_1)+( ifld-1)* qq+ nscalars*0 )
            ! Post-Collision values of neighbor element
            fTmpAll_ff(vPos) = state(                                    &
              &  neigh((idir-1)* nsize+ neighpos_1)+( ifld-1)* qq+ nscalars*0 )
          end do
        end do

        ! calculate spc density, velocity, and velAvg from PDF of first fluid
        call calcMacrosFromPDF( varSys     = varSys,      &
          &                     pdf        = fTmpAll_f,   &
          &                     mass_dens  = mass_dens_f, &
          &                     velocity   = velocity_f,  &
          &                     velAvg     = velAvg_f,    &
          &                     eqVel      = eqVel_f,     &
          &                     species    = species      )


        ! calculate spc density, velocity, and velAvg from PDF of second fluid
        call calcMacrosFromPDF( varSys     = varSys,       &
          &                     pdf        = fTmpAll_ff,   &
          &                     mass_dens  = mass_dens_ff, &
          &                     velocity   = velocity_ff,  &
          &                     velAvg     = velAvg_ff,    &
          &                     eqVel      = eqVel_ff,     &
          &                     species    = species       )

        ! fluid node velocity square term of equilibrium function
        usq_f = dot_product( velAvg_f(1:3), velAvg_f(1:3) ) * t2cs2inv
        usq_ff = dot_product( velAvg_ff(1:3), velAvg_ff(1:3) ) * t2cs2inv

        ! density of current species at boundary node
        !mass_dens_b = (4.0_rk * mass_dens_f(iField) - mass_dens_ff(iField)) / 3.0_rk
        !if (mass_dens_b < 0.0_rk) then
        !  mass_dens_b = mixture%moleDens0LB * 0.001 &
        !    &         * species(iField)%molWeight
        !end if
        mass_dens_b = mixture%moleDens0LB * moleFrac_w(iElem) &
          &         * species(iField)%molWeight
        !mass_dens_b = ( 4.0_rk * mass_dens_f(iField) - mass_dens_ff(iField) ) &
        !  &         / 3.0_rk
        !mass_dens_b = mixture%rho0LB - (sum(mass_dens_f) - mass_dens_b)
        !mass_dens_b = 3.0_rk * mixture%moleDens0LB * moleFrac_w(iElem) &
        !  &         * species(iField)%molWeight                         &
        !  &         - mass_dens_f(iField) - mass_dens_ff(iField)
        !mass_dens_b = coeff_w * mixture%moleDens0LB * moleFrac_w(iElem) &
        !  &         * species(iField)%molWeight                         &
        !  &         + coeff_f * mass_dens_f(iField)                     &
        !  &         + coeff_ff * mass_dens_ff(iField)

        ! velocity at boundary node
        vel_b = vel_w( (iElem-1)*3 + 1: iElem*3 )
        !vel_b = velocity_ff(1:3, iField)
        !vel_b = coeff_w * vel_w( (iElem-1)*3 + 1: iElem*3 ) &
        !  &   + coeff_f * velocity_f(1:3, iField)           &
        !  &   + coeff_ff * velocity_ff(1:3, iField)
        !vel_b = coeff_w * vel_w( (iElem-1)*3 + 1: iElem*3 ) &
        !  &   + coeff_f * velAvg_f(1:3) + coeff_ff * velAvg_ff(1:3)

!write(dbgUnit(1),*) 'vel_w ', vel_w( (iElem-1)*3 + 1: iElem*3 )
!write(dbgUnit(1),*) 'mass_dens_b ', mass_dens_b, ' vel_b ', vel_b

        ! boundary node velocity square term of equilibrium function
        usq_b = dot_product( vel_b(1:3), vel_b(1:3) ) * t2cs2inv

        ! Position of this boundary element in state array
        elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)

        do iDir = 1, layout%fStencil%QQN
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            ! Compute equilibrium on fluid node
            ucx_f = dot_product( stencil%cxDirRK(1:3, iDir), velAvg_f )
            ucxQuad_f = dot_product( stencil%cxDirRK(1:3, iDir), velAvg_f(1:3) )

            fEq_f = layout%weight(iDir) * mass_dens_f(iField)       &
              &   * ( species(iField)%molWeigRatio + ucx_f * cs2inv &
              &   + ucxQuad_f * ucxQuad_f * t2cs4inv - usq_f        )

            ! Compute equilibrium on fluid node
            ucx_ff = dot_product( stencil%cxDirRK(1:3, iDir), velAvg_ff )
            ucxQuad_ff = dot_product( stencil%cxDirRK(1:3, iDir), &
              &                       velAvg_ff(1:3) )

            fEq_ff = layout%weight(iDir) * mass_dens_ff(iField)      &
              &   * ( species(iField)%molWeigRatio + ucx_ff * cs2inv &
              &   + ucxQuad_ff * ucxQuad_ff * t2cs4inv - usq_ff      )

            ! Compute equilibrium on boundary node
            ucx_b = dot_product( stencil%cxDirRK(1:3, iDir), vel_b )

            fEq_b = layout%weight(iDir) * mass_dens_b               &
              &   * ( species(iField)%molWeigRatio + ucx_b * cs2inv &
              &   + ucx_b * ucx_b * t2cs4inv - usq_b                )

            ! Position of current field variable current dir in state array
            !pdf_f = bcBuffer( (posInBuffer-1)*nScalars + varPos(iDir) )

            fnEq_f = fTmpAll_f(varPos(iDir)) - fEq_f
            fnEq_ff = fTmpAll_ff(varPos(iDir)) - fEq_ff
            fnEq_b = 0.5_rk * fnEq_f + 0.5_rk * fnEq_ff

            state( neigh((idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)  &
              !& = fEq_b + ( pdf_f - fEq_f )
              !& = fEq_b + fnEq_b
              & = fEq_b + fnEq_f
          end if
        end do

        !iDir = layout%fStencil%restPosition
        !fEq_b = layout%weight(iDir) * mass_dens_b * ( weight0Inv             &
        ! &    + (1.0_rk - weight0Inv) * species(iField)%molWeigRatio - usq_b )

        !state( neigh((idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
        !  & = fEq_b
      end do
    end associate
!    stop

  contains
    !> This routine computes macroscopic quantities like density, velocity
    !! equilibrium velocity, velocity average for given iField from PDF
    subroutine calcMacrosFromPDF( varSys, pdf, mass_dens, velocity, velAvg, &
      &                           eqVel, species                            )
      ! -----------------------------------------------------------------------
      type(tem_varSys_type), intent(in) :: varSys
      real(kind=rk), intent(in) :: pdf(varSys%nScalars)
      real(kind=rk), intent(out) :: mass_dens(nFields)
      real(kind=rk), intent(out) :: velocity(nFields, nFields)
      real(kind=rk), intent(out) :: velAvg(3)
      real(kind=rk), intent(out) :: eqVel(3)
      type(mus_species_type), intent(in) :: species(nFields)
      ! -----------------------------------------------------------------------
      integer :: iFld, iDir, vPos
      real(kind=rk) :: fTmp(layout%fStencil%QQ)
      real(kind=rk) :: tot_numDens, tot_massDens
      real(kind=rk) :: momentum(3, nFields), first_moments(3, nFields)
      real(kind=rk) :: num_dens(nFields), moleFrac(nFields)
      ! -----------------------------------------------------------------------
      do iFld = 1, nFields
        do iDir = 1, layout%fStencil%QQ
          ! Position of current field variable current dir in state array
          vPos = varSys%method%val(iFld)%state_varPos(iDir)
          ! State array of current element
          fTmp(iDir) = pdf(vPos)
        end do
        ! species density
        mass_dens(iFld) = sum(fTmp)
        ! number density
        num_dens(iFld) = mass_dens(iFld) * species(iFld)%molWeightInv
        ! first moments: momentum of transformed PDF
        first_moments(1, iFld) = sum( fTmp * layout%fStencil%cxDirRK(1, :) )
        first_moments(2, iFld) = sum( fTmp * layout%fStencil%cxDirRK(2, :) )
        first_moments(3, iFld) = sum( fTmp * layout%fStencil%cxDirRK(3, :) )
      end do

      !total number density
      tot_numDens = sum(num_dens)
      !total mass density
      tot_massDens = sum(mass_dens)
      ! mole fraction
      moleFrac(:) = num_dens(:) / tot_numDens

      ! actual momentum of all species obtained from solving
      ! Linear Equation system
      momentum = momentumFromMacroLSE(                      &
        &          moleFraction  = moleFrac,                &
        &          first_moments = first_moments,           &
        &          nFields       = nFields,                 &
        &          phi           = species(:)%molWeigRatio, &
        &          resi_coeff    = resi_coeff,              &
        &          omega_fac     = omega_fac,               &
        &          paramBInv     = paramBInv                )

      ! velocity of all species on current fluid node
      do iFld = 1, nFields
        velocity(:, iFld) = momentum(:, iFld) / mass_dens(iFld)
      end do

      ! equilibrium velocity for current species
      eqVel(:) = velocity(:, iField)
      do iFld = 1, nFields
        eqVel(:) = eqVel(:) + resi_coeff(iField, iFld)           &
          &      * species(iField)%molWeigRatio * moleFrac(iFld) &
          &      * ( velocity(:, iFld) - velocity(:, iField) )   &
          &      * paramBInv
      end do

      ! Mass averaged mixture velocity
      velAvg(1) = dot_product( mass_dens, velocity(1, :) ) / tot_massDens
      velAvg(2) = dot_product( mass_dens, velocity(2, :) ) / tot_massDens
      velAvg(3) = dot_product( mass_dens, velocity(3, :) ) / tot_massDens

    end subroutine calcMacrosFromPDF

  end subroutine spc_inflow
! ****************************************************************************** !



! ****************************************************************************** !
  !> molar diffusion flux boundary condition
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moleDiff_flux',
  !!    mole_diff_flux = 'mole_diff_flux'
  !! }
  !!}
  !!variable = {
  !!  name = 'mole_diff_flux',
  !!  ncomponents = 3,
  !!  vartype = 'st_fun',
  !!  st_fun = {0.06, 0.0, 0.0}
  !!}
  !!```
  !! @todo KM Adapt moleDiff flux boundary input_loc for equilibrium function
  !! pointer
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moleDiff_Flux( me, state, bcBuffer, globBC, levelDesc, tree, &
    &                           nSize, iLevel, sim_time, neigh, layout,       &
    &                           fieldProp, varPos, nScalars, varSys,          &
    &                           derVarPos, physics, iField, mixture           )
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
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    real(kind=rk) :: mixMolvel( 3 )
    real(kind=rk) :: mass_dens( varSys%nStateVars ), &
      &              num_dens( varSys%nStateVars )
    real(kind=rk) :: tot_NumDens, tot_MassDens
    real(kind=rk) :: moleDiff_Flux(globBC%nElems(iLevel)*3), inv_flux
    real(kind=rk) :: moleFrac(varSys%nStateVars)
    integer :: iElem, iDir, iFieldLoc, nFields, iComp, pos, QQ
    integer :: bcMoleDiffFlux_pos, posInBuffer, posInState
    real(kind=rk) :: molWeightInv(varSys%nStateVars), phi(varSys%nStateVars)
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: paramBInv, omega_fac
    type(mus_varSys_data_type), pointer :: fPtr
    !------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )

    do iFieldLoc = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeightInv(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                           %fieldProp%species%molWeightInv
      ! molecular weight ratio
      phi(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                  %fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFieldLoc, :) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                            %fieldProp%species%resi_coeff(:)
    end do

    omega_fac = mixture%omega_diff * 0.5_rk
    paramBInv = 1.0_rk/mixture%paramB

    inv_flux = 1.0_rk / physics%fac( iLevel )%moleFlux
    ! position of boundary moleflux in varSys
    bcMoleDiffFlux_pos = me%bc_states%moleDiff_flux%varPos
    ! Get moleflux
    call varSys%method%val(bcMoleDiffFlux_pos)%get_valOfIndex( &
      & varSys  = varSys,                                      &
      & time    = sim_time,                                    &
      & iLevel  = iLevel,                                      &
      & idx     = me%bc_states%moleDiff_flux                   &
      &           %pntIndex%indexLvl(iLevel)                   &
      &           %val(1:globBC%nElems(iLevel)),               &
      & nVals   = globBC%nElems(iLevel),                       &
      & res     = moleDiff_Flux                                )

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
    moleDiff_Flux = moleDiff_Flux * inv_flux

    ! Calculate the mole averaged mixture velocity w = sum(n_i v_i)/n
    do iElem = 1, globBC%nElems(iLevel)
      ! initialize the local velocity for every element
      first_moments = 0.0_rk
      momentum = 0.0_rk
      mass_dens = 0.0_rk
      num_dens = 0.0_rk
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      do iFieldLoc = 1, nFields
        ! all fields have some nComponents
        do iDir = 1,layout%fStencil%QQ
          ! get the position of the current direction of the depending variable
          ! in the global system (= position in the input)
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          ! density
          mass_dens(iFieldLoc) = mass_dens(iFieldLoc) + bcBuffer(       &
            & ( posinbuffer-1)* nscalars+ pos )

          first_moments( :, iFieldLoc ) = first_moments( :, iFieldLoc )   &
& + bcBuffer( ( posinbuffer-1)* nscalars+ pos ) &
              & * layout%fStencil%cxDirRK( :, iDir )
        end do !iDir
        ! number density
        num_dens(iFieldLoc) = mass_dens(iFieldLoc) * molWeightInv(iFieldLoc)
      end do !iField

      !total number density
      tot_NumDens = sum(num_dens)
      !total mass density
      tot_massDens = sum(mass_dens)
      ! mole fraction
      moleFrac(:) = num_dens(:) / tot_NumDens

      ! momentum of all species
      momentum = momentumFromMacroLSE( moleFraction  = moleFrac,       &
        &                              first_moments = first_moments,  &
        &                              nFields       = nFields,        &
        &                              phi           = phi,            &
        &                              resi_coeff    = resi_coeff,     &
        &                              omega_fac     = omega_fac,      &
        &                              paramBInv     = paramBInv       )

      !velocity of all species
      do iFieldLoc = 1,nFields
        vel( :, iFieldLoc) = momentum( :, iFieldLoc) / mass_dens(iFieldLoc)
      end do

      !mole averaged mixture velocity, w = sum(n_i*v_i)/n_T
      do iComp = 1, 3
        mixMolVel(iComp) = dot_product( num_dens, vel(iComp,:) )/tot_numDens
      end do

      ! update current species velocity with defined mole diffusion flux
      ! species velocity at boundary element
      ! moleDiff_flux, J_i = n_i(v_i - w )
      ! v_i = J_i/n_i + w
      vel(:, iField) = moleDiff_Flux((iElem-1)*3+1:iElem*3) &
        &            / num_dens(iField)                     &
        &            + mixMolVel

      ! compute equilibrium
      call derVarPos%equilFromMacro( density  = mass_dens, &
        &                            velocity = vel,       &
        &                            iField   = iField,    &
        &                            nElems   = 1,         &
        &                            varSys   = varSys,    &
        &                            layout   = layout,    &
        &                            res      = fEq        )

      posInState = globBC%elemLvl( iLevel )%elem%val( iElem )
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          state(neigh (( idir-1)* nsize+ posinstate)+( ifield-1)* qq+ nscalars*0) &
            & = fEq( iDir )
        end if
      end do
    end do

  end subroutine spc_moleDiff_Flux
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine computes mole diffusion flux of the ionic species at the
  !! membrance using black box model and then mass density at the membrane
  !! boundary from mole diffusion flux. Then equilibrium is set at the boundary
  !! which is computed from mass density and velocity
  !!
  !! Black box model: \n
  !!  - for ionic species mole diffusion flux:
  !! $J^m_k (x^m,t) = \frac{t^m_k(x^m,t)}{z_k F} i^m(x^m,t)$
  !! $J^m_k$ - mole diffusion flux of ion in the membrane
  !! $t^m_k$ - transference number of ion of the membrane
  !!           (selective membrane property)
  !! $z_k$ - charge number of ion
  !! $F$ - Faraday constant
  !! $i^m(x^m,t) = n^l_x_1 \cdot i^l(x,t) = - n^r_x_1 \cdot i^r(x,t)$
  !!    - orthogonal projection of the current density at the left $i^l$ and
  !!      right $i^r$ side of the membrance.
  !! The normal vectors $n^l_1$ and $n^r_1$ point into the membrane surface
  !! $i(x,t) = F \sum_{j=1}^{n} z_j J^e_j(x,t)$
  !!    - local current density [A/m^3]
  !! $J^e_k = c^e_k(v^e_k - v) $ - mole diffusion flux of ion in the electrolyte
  !! In black box model: mixture averaged velocity at the membrane is assumed to
  !! be zero
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_blackbox_mem_ion( me, state, bcBuffer, globBC, levelDesc,     &
    &                              tree, nSize, iLevel, sim_time, neigh,       &
    &                              layout, fieldProp, varPos, nScalars,        &
    &                              varSys, derVarPos, physics, iField, mixture )
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
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*varSys%nStateVars )
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: num_dens( varSys%nStateVars )
    integer :: iElem, iDir, iFieldLoc, nFields, pos, QQ
    integer :: posInBuffer, posInState
    real(kind=rk) :: molWeight(varSys%nStateVars), chargeNr(varSys%nStateVars)
    real(kind=rk), dimension(3, varSys%nStateVars):: molar_flux, vel, momentum
    real(kind=rk) :: loc_curr_dens(3)
    real(kind=rk) :: mem_curr_dens, mem_molar_flux
    real(kind=rk) :: paramBInv, omega_fac
    real(kind=rk) :: molefrac(3), mem_mass_flux, elec_mass_flux(3)
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    real(kind=rk) :: tot_NumDens, tot_MassDens
    real(kind=rk) :: molWeightInv(varSys%nStateVars), phi(varSys%nStateVars)
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    type(mus_varSys_data_type), pointer :: fPtr
    ! --------------------------------------------------------------------------
!write(dbgUnit(1),*) 'Boundary label ', trim(me%label)
!write(dbgUnit(1),*) 'iField ', iField
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )

    do iFieldLoc = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFieldLoc) =  fPtr%solverData%scheme%field(iFieldLoc)   &
        &                                    %fieldProp%species%molWeight
      molWeightInv(iFieldLoc) =  fPtr%solverData%scheme%field(iFieldLoc)   &
        &                                    %fieldProp%species%molWeightInv
      ! charge number
      chargeNr(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                                  %fieldProp%species%chargeNr
      ! molecular weight ratio
      phi(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                  %fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFieldLoc, :) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                            %fieldProp%species%resi_coeff(:)
    end do

    omega_fac = mixture%omega_diff * 0.5_rk
    paramBInv = 1.0_rk/mixture%paramB

    ! Calculate the mole averaged mixture velocity w = sum(n_i v_i)/n
    do iElem = 1, globBC%nElems(iLevel)
!write(dbgUnit(1),*) 'iElem ', iElem
      ! initialize the local velocity for every element
      first_moments = 0.0_rk
      mass_dens = 0.0_rk
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
      do iFieldLoc = 1, nFields
        ! all fields have some nComponents
        do iDir = 1,layout%fStencil%QQ
          ! get the position of the current direction of the depending variable
          ! in the global system (= position in the input)
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all( pos ) = bcBuffer(                                   &
            & ( posinbuffer-1)* nscalars+ pos )
          ! density
          mass_dens(iFieldLoc) = mass_dens(iFieldLoc) + fTmp_all(pos)
          first_moments( :, iFieldLoc ) = first_moments( :, iFieldLoc ) &
            & + fTmp_all(pos) * layout%fStencil%cxDirRK( :, iDir )
        end do !iDir
        ! number density
        num_dens(iFieldLoc) = mass_dens(iFieldLoc) * molWeightInv(iFieldLoc)
      end do !iField

      !total number density
      tot_NumDens = sum(num_dens)
      !total mass density
      tot_massDens = sum(mass_dens)
      ! mole fraction
      moleFrac(:) = num_dens(:) / tot_NumDens

      ! momentum of all species
      momentum = momentumFromMacroLSE( moleFraction  = moleFrac,      &
        &                              first_moments = first_moments, &
        &                              nFields       = nFields,       &
        &                              phi           = phi,           &
        &                              resi_coeff    = resi_coeff,    &
        &                              omega_fac     = omega_fac,     &
        &                              paramBInv     = paramBInv      )

      !velocity of all species
      do iFieldLoc = 1,nFields
        vel( :, iFieldLoc) = momentum( :, iFieldLoc) / mass_dens(iFieldLoc)
      end do

      ! molar diffusive flux = molar flux: Assuming mixture velocity is zero
      ! at membrane surface
      do iFieldLoc = 1, nFields
        molar_flux(:, iFieldLoc) = num_dens(iFieldLoc) * vel(:, iFieldLoc)
      end do

      ! KM: \todo 2016_08_05 Use current density variable
      ! local current density
      loc_curr_dens = 0.0_rk
      do iFieldLoc = 1, nFields
        loc_curr_dens(:) = loc_curr_dens(:) + chargeNr(iFieldLoc) &
          &              * molar_flux(:, iFieldLoc)
      end do

      loc_curr_dens(:) = loc_curr_dens * mixture%faradayLB

      ! project local current density on the membrane
      mem_curr_dens = dot_product(abs(globBC%elemLvl(iLevel)%normal%val(:, iElem)), &
        &                             loc_curr_dens)

      ! membrance molar flux of ionic species
      mem_molar_flux = me%blackbox_mem%transNr * mem_curr_dens                 &
        &            / ( chargeNr(iField) * mixture%faradayLB )

      ! membrance mass flux
      mem_mass_flux = mem_molar_flux*molWeight(iField)
      ! set normal direction of membrance flux
      elec_mass_flux = 0.0_rk
      if( abs(globBC%elemLvl(iLevel)%normal%val(1, iElem)) == 1 )                &
        & elec_mass_flux(1) = mem_mass_flux
      if( abs(globBC%elemLvl(iLevel)%normal%val(2, iElem)) == 1 )                &
        & elec_mass_flux(2) = mem_mass_flux
      if( abs(globBC%elemLvl(iLevel)%normal%val(3, iElem)) == 1 )                &
        & elec_mass_flux(3) = mem_mass_flux

!write(dbgUnit(1),*) 'transNr ', me%blackbox_mem%transNr
!write(dbgUnit(1),*) 'moleFlux ', molar_flux(:, iField)
!write(dbgUnit(1),*) 'moleFluxAll ', molar_flux(:, :)
!write(dbgUnit(1),*) 'loc_curr_dens ', loc_curr_dens
!write(dbgUnit(1),*) 'mem_curr_dens ', mem_curr_dens
!write(dbgUnit(1),*) 'mem_mass_flux ', mem_mass_flux
!write(dbgUnit(1),*) 'elec_mass_flux ', elec_mass_flux
      posInState = globBC%elemLvl(iLevel)%elem%val(iElem)
      do iDir = 1,layout%fStencil%QQN
        ! Write the values
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          state(                                                               &
            & neigh((idir-1)* nsize+ posinstate)+( ifield-1)* qq+ nscalars*0) = &
          ! We need to get post-collision pdf in direction
          ! alpha- outgoing direction, which is the inverse direction of bitmask
          ! For PULL this means, get the outgoing one, as this is the one which
          ! will be bounced back
          ! For PUSH this means, get the already bounced back pdf back, so take
          ! the incoming
            & fTmp_all(varPos(layout%fStencil%cxDirInv(iDir)))               &
            &       + layout%weight( iDir )*2._rk*cs2inv                     &
            &       * ( layout%fStencil%cxDir( 1, iDir )*elec_mass_flux(1)   &
            &       +   layout%fStencil%cxDir( 2, iDir )*elec_mass_flux(2)   &
            &       +   layout%fStencil%cxDir( 3, iDir )*elec_mass_flux(3) )
        end if
      end do
    end do !iElem

  end subroutine spc_blackbox_mem_ion
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine computes mole diffusion flux of the solvent at the membrance
  !! using black box model and then mass density at the membrane boundary
  !! from mole diffusion flux. Then equilibrium is set at the boundary which is
  !! computed from mass density and velocity
  !! @todo KM: Not completeed. Requires ionic strength at left and right side of
  !! membrane which is not available at the moment
  !!
  !! Black box model: \n
  !!  - for solvent mole diffusion flux:
  !! $J^m_s (x^m,t) = \frac{t^m_s(x^m,t)}{F} i^m(x^m,t)
  !!                + LW(IS^l(x^m,t)-IS^r(x^m,t))$
  !!
  !! $J^m_s$ - mole diffusion flux of solvent through the membrane
  !! $t^m_k$ - transference number of solvent on the membrane
  !!           (selective membrane property)
  !! $F$ - Faraday constant
  !! $LW$ - osmotic solvent transport coefficient
  !! $IS^l$ and $IS^r$ are ionic strength of electrolyte solution at the left
  !! and right side of membrance
  !! $IS^l(x^m,t)=\frac{1}{2}\sum^n_{k=1} z^2_k c^m_k(x^m,t)$
  !! $c^m_k = \rho^m_k/m_k$ - molar concentration
  !! $i^m(x^m,t) = n^l_x_1 \cdot i^l(x,t) = - n^r_x_1 \cdot i^r(x,t)$
  !!    - orthogonal projection of the current density at the left $i^l$ and
  !!      right $i^r$ side of the membrance.
  !! The normal vectors $n^l_1$ and $n^r_1$ point into the membrane surface
  !! $i(x,t) = F \sum_{j=1}^{n} z_j J^e_j(x,t)$
  !!    - local current density [A/m^3]
  !! $J^e_k = c^e_k(v^e_k - v) $ - mole diffusion flux of ion in the electrolyte
  !! In black box model: mixture averaged velocity at the membrane is assumed to
  !! be zero
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_blackbox_mem_solvent( me, state, bcBuffer, globBC, levelDesc, &
    &                                  tree, nSize, iLevel, sim_time, neigh,   &
    &                                  layout, fieldProp, varPos, nScalars,    &
    &                                  varSys, derVarPos, physics, iField,     &
    &                                  mixture                                 )
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
    !store pdf of all fields for derive function pointer input
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*varSys%nStateVars )
    real(kind=rk) :: mass_dens( varSys%nStateVars ), &
      &              num_dens( varSys%nStateVars )
    integer :: iElem, iDir, iFieldLoc, nFields, pos, bndNormalDir, QQ
    real(kind=rk) :: molWeight(varSys%nStateVars), chargeNr(varSys%nStateVars)
    real(kind=rk), dimension(3, varSys%nStateVars):: molar_flux, vel, momentum
    real(kind=rk) :: loc_curr_dens(3)
    real(kind=rk) :: mem_curr_dens, mem_molar_flux
    real(kind=rk) :: molefrac(3), mem_mass_flux, elec_mass_flux(3)
    real(kind=rk) :: tot_NumDens, tot_MassDens
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    real(kind=rk) :: molWeightInv(varSys%nStateVars), phi(varSys%nStateVars)
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: paramBInv, omega_fac
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: posInBuffer, posInState
    ! ------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )

    do iFieldLoc = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFieldLoc) =  fPtr%solverData%scheme%field(iFieldLoc) &
        &                         %fieldProp%species%molWeight
      molWeightInv(iFieldLoc) =  fPtr%solverData%scheme%field(iFieldLoc)   &
        &                                    %fieldProp%species%molWeightInv
      ! charge number
      chargeNr(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                       %fieldProp%species%chargeNr
      ! molecular weight ratio
      phi(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                  %fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFieldLoc, :) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                            %fieldProp%species%resi_coeff(:)
    end do

    omega_fac = mixture%omega_diff * 0.5_rk
    paramBInv = 1.0_rk/mixture%paramB

    ! Calculate the mole averaged mixture velocity w = sum(n_i v_i)/n
    do iElem = 1, globBC%nElems(iLevel)
      ! initialize the local velocity for every element
      first_moments = 0.0_rk
      mass_dens = 0.0_rk
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
      do iFieldLoc = 1, nFields
        ! all fields have some nComponents
        do iDir = 1,layout%fStencil%QQ
          ! get the position of the current direction of the depending variable
          ! in the global system (= position in the input)
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all( pos ) = bcBuffer(                                  &
            & ( posinbuffer-1)* nscalars+ pos)
          ! density
          mass_dens(iFieldLoc) = mass_dens(iFieldLoc) + fTmp_all(pos)
          first_moments( :, iFieldLoc ) = first_moments( :, iFieldLoc ) &
            & + fTmp_all(pos) * layout%fStencil%cxDirRK( :, iDir )
        end do !iDir
        ! number density
        num_dens(iFieldLoc) = mass_dens(iFieldLoc) * molWeightInv(iFieldLoc)
      end do !iField

      !total number density
      tot_NumDens = sum(num_dens)
      !total mass density
      tot_massDens = sum(mass_dens)
      ! mole fraction
      moleFrac(:) = num_dens(:) / tot_NumDens

      ! momentum of all species
      momentum = momentumFromMacroLSE( moleFraction  = moleFrac,      &
        &                              first_moments = first_moments, &
        &                              nFields       = nFields,       &
        &                              phi           = phi,           &
        &                              resi_coeff    = resi_coeff,    &
        &                              omega_fac     = omega_fac,     &
        &                              paramBInv     = paramBInv      )

      !velocity of all species
      do iFieldLoc = 1,nFields
        vel( :, iFieldLoc) = momentum( :, iFieldLoc) / mass_dens(iFieldLoc)
      end do

      ! molar diffusive flux = molar flux
      do iFieldLoc = 1, nFields
        molar_flux(:, iFieldLoc) = num_dens(iFieldLoc) &
          & * (vel(:, iFieldLoc))
      end do

      ! local current density
      loc_curr_dens = 0.0_rk
      do iFieldLoc = 1, nFields
        loc_curr_dens(:) = loc_curr_dens(:) + chargeNr(iFieldLoc) &
          &              * molar_flux(:, iFieldLoc)
      end do

      loc_curr_dens = loc_curr_dens * mixture%faradayLB

      ! normal direction point to boundary
      bndNormalDir = layout%fStencil%cxDirInv( globBC%elemLvl( iLevel )%       &
        &                                               normalInd%val( iElem ))
      ! project local current density on the membrane
      mem_curr_dens = dot_product(loc_curr_dens,                          &
        &                         layout%fStencil%cxDirRK(:, bndNormalDir))

      ! membrance molar flux of ionic species
      mem_molar_flux = me%blackbox_mem%transNr * mem_curr_dens                 &
        &            / ( mixture%faradayLB )

      ! membrance mass flux
      mem_mass_flux = mem_molar_flux*molWeight(iField)
      elec_mass_flux = 0.0_rk
      ! set normal direction of membrance flux
      if( abs(globBC%elemLvl(iLevel)%normal%val(1, iElem)) == 1 )                &
        & elec_mass_flux(1) = mem_mass_flux
      if( abs(globBC%elemLvl(iLevel)%normal%val(2, iElem)) == 1 )                &
        & elec_mass_flux(2) = mem_mass_flux
      if( abs(globBC%elemLvl(iLevel)%normal%val(3, iElem)) == 1 )                &
        & elec_mass_flux(3) = mem_mass_flux

      posInState = globBC%elemLvl(iLevel)%elem%val(iElem)
      do iDir = 1,layout%fStencil%QQN
        ! Write the values
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state(                                                               &
            & neigh((idir-1)* nsize+ posinstate)+( ifield-1)* qq+ nscalars*0) = &
          ! We need to get post-collision pdf in direction
          ! alpha- outgoing direction, which is the inverse direction of bitmask
          ! For PULL this means, get the outgoing one, as this is the one which
          ! will be bounced back
          ! For PUSH this means, get the already bounced back pdf back, so take
          ! the incoming
            & fTmp_all(varPos(layout%fStencil%cxDirInv(iDir)))             &
            &       + layout%weight( iDir )*2._rk*cs2inv                   &
            &       * ( layout%fStencil%cxDir( 1, iDir )*elec_mass_flux(1) &
            &       +   layout%fStencil%cxDir( 2, iDir )*elec_mass_flux(2) &
            &       +   layout%fStencil%cxDir( 3, iDir )*elec_mass_flux(3) )
        end if
      end do
    end do !iElem

  end subroutine spc_blackbox_mem_solvent
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Moment based wall boundary condition from Sam Bennent PhD thesis
  !! "A Lattice Boltzmann Model for Diffusion of Binary Gas Mixtures"
  !!
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moments_wall( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                          nSize, iLevel, sim_time, neigh, layout,         &
    &                          fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                          physics, iField, mixture                        )
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
    integer :: iElem, iDir, tmpDir, QQ
    integer :: nLinks, iLink, elemPos
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

!    write(*,*) 'Boundary label ', trim(me%label)
    do iElem = 1, globBC%nElems( iLevel )
      if( .not. btest( levelDesc%property(                        &
        &              globBC%elemLvl(iLevel)%elem%val(iElem)), prp_solid))then

      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(rhs(nLinks))
      allocate(unKnown_fTmp(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
!      write(*,*) 'elemPos ', elemPos
!      write(*,*) 'treeID ', levelDesc%total(elemPos)

      iLink = 0
      fTmp = 0.0_rk
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink +1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = &
            &state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = &
& state(neigh (( layout%fstencil%restposition-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)
!        write(*,*) 'fTmp0 ',  fTmp

      rhoB = fTmp(layout%fStencil%restPosition)
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          tmpDir = layout%fStencil%cxDirInv(iDir)
          if( globBC%elemLvl(iLevel)%bitmask%val( tmpDir, iElem )) then
            tmpDir = layout%fStencil%cxDirInv( &
              &            globBC%elemLvl(iLevel)%normalInd%val(iElem))
          end if
        else
          tmpDir = iDir
        end if
        rhoB  = rhoB + fTmp(tmpDir)
      end do

      press = cs2*fieldProp%species%molWeigRatio*rhoB
!write(*,*) 'press ', press

      ! moments for multispecies LBM model of transformed variable
      first_moments = (/ 0.0_rk, 0.0_rk, 0.0_rk /) !momX, momY, momZ
      second_moments = (/ press, & !momXX
        &                 press, & !momYY
        &                 press, & !momZZ
        &                 0.0_rk, & !momXY
        &                 0.0_rk, & !momYZ
        &                 0.0_rk /) !momXZ
      third_moments = 0.0_rk
      fourth_moments = cs2*press !momXXYY, momYYZZ, momZZXX

      ! moments for liquid mixture
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
       state(                                                                 &
&neigh (( missing_links(ilink)-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
        & )  = unKnown_fTmp(iLink)
     end do

      deallocate(missing_links)
      deallocate(rhs)
      deallocate(unKnown_fTmp)
      end if
    end do

  end subroutine spc_moments_wall
! ****************************************************************************** !


! ****************************************************************************** !
  !> Mole fraction boundary condition
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moments_moleFrac',
  !!    moleFraction = 0.0
  !!     }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moments_moleFrac( me, state, bcBuffer, globBC, levelDesc,     &
    &                              tree, nSize, iLevel, sim_time, neigh,       &
    &                              layout, fieldProp, varPos, nScalars,        &
    &                              varSys, derVarPos, physics, iField, mixture )
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
    real(kind=rk) :: rho, press
    real(kind=rk) :: moleFrac(globBC%nElems(iLevel))
    integer :: iElem, iDir, iFieldLoc, nFields, pos
    real(kind=rk) :: molWeight(varSys%nStateVars), phi(varSys%nStateVars)
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    integer :: nLinks, iLink, elemPos, QQ
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp(:)
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*varSys%nStateVars )
    real(kind=rk) :: vel_inG(3), velAvg(3), velQuad(3)
    real(kind=rk) :: mom(3,varSys%nStateVars), mom_inG(3)
    real(kind=rk) :: mass_dens(varSys%nStateVars), num_dens(varSys%nStateVars)
    real(kind=rk) :: totMassDens
    real(kind=rk) :: moleFrac_loc(varSys%nStateVars)
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: first_moments(3), second_moments(6), third_moments(6), &
      &              fourth_moments(3)
    integer :: nMom(4)
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: bcMoleFrac_pos, posInBuffer
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )

    do iFieldLoc = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                        %fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                  %fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFieldLoc, :) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                            %fieldProp%species%resi_coeff(:)
    end do

    ! position of molefraction spacetime variable in varSys
    bcMoleFrac_pos = me%bc_states%moleFrac%varPos
    ! mole fraction
    call varSys%method%val(bcMoleFrac_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFrac                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFrac                                 )

    nMom(1) = 1 + size(layout%moment%first_moments)
    nMom(2) = nMom(1) + size(layout%moment%second_moments)
    nMom(3) = nMom(2) + size(layout%moment%third_moments)
    nMom(4) = nMom(3) + size(layout%moment%fourth_moments)

    ! Calculate the density of current element
    do iElem = 1, globBC%nElems(iLevel)
      fTmp = 0.0_rk
      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(unKnown_fTmp(nLinks))
      allocate(rhs(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )

      iLink = 0
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink + 1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = &
            &state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = &
&state(neigh (( layout%fstencil%restposition-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)

      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
      ! local state vector to compute density and velocity of other species
      do iFieldLoc = 1, varSys%nStateVars
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all( pos ) = bcBuffer(                                  &
            & ( posinbuffer-1)* nscalars+ pos)
        end do
      end do

      !KM: using inital mixture number density rho0/mixtureMOlWeight
      rho = mixture%moleDens0LB *  moleFrac(iElem) &
        & * fieldProp%species%molWeight
      press = cs2 * rho * fieldProp%species%molWeigRatio

      !local density and momentum
      mass_dens = 0.0_rk
      do iFieldLoc = 1, nFields
        mass_dens(iFieldLoc) = &
          & sum(fTmp_all(varSys%method%val(iFieldLoc)%state_varPos(:)))
        num_dens(iFieldLoc) = mass_dens(iFieldLoc) &
          &                 / molWeight(iFieldLoc)
        !velocity
        call derVarPos%momFromState( state  = fTmp_all,         &
          &                          iField = iFieldLoc,        &
          &                          nElems = 1,                &
          &                          varSys = varSys,           &
          &                          layout = layout,           &
          &                          res    = mom(:, iFieldLoc) )
      end do

      mass_dens(iField) = rho
      num_dens(iField) = rho / fieldProp%species%molWeight

      ! total mass density
      totMassDens = sum(mass_dens)

      ! mole fraction
      moleFrac_loc = num_dens/sum(num_dens)

!do iFieldLoc = 1, nFields
!  write(dbgUnit(1),*) iFieldLoc, 'vel ', mom(:, iFieldLoc)/mass_dens(iFieldLoc), &
!    & 'dens ', mass_dens(iFieldLoc), 'moleFrac ', moleFrac(iFieldLoc)
!end do

      ! compute velocity of transformed variable g. Eq. which is normally
      ! used to compute actual velocity by solving LSE
      mom_inG = mom(:, iField)
      do iFieldLoc = 1, nFields
        mom_inG = mom_inG + ( mixture%omega_diff * 0.5_rk                    &
            &               * resi_coeff(iField, iFieldLoc) * phi(iField)    &
            &               * moleFrac_loc(iFieldLoc) / mixture%paramB )     &
            &               * mom(:, iField)
      end do
      !set non-diagonal part
      do iFieldLoc = 1, nFields
        mom_inG = mom_inG - ( mixture%omega_diff * 0.5_rk                    &
          &                 * resi_coeff(iField, iFieldLoc) * phi(iFieldLoc) &
          &                 * moleFrac_loc(iField) / mixture%paramB )        &
          &                 * mom(:,iFieldLoc)
      end do

      vel_inG = mom_inG/rho
!write(dbgUnit(1),*) 'transformed momentum ', mom_inG
!write(dbgUnit(1),*) 'transformed velocity ', vel_inG

      ! mass averaged mixture velocity
      velAvg(1) = sum(mom(1,:)) / totMassDens
      velAvg(2) = sum(mom(2,:)) / totMassDens
      velAvg(3) = sum(mom(3,:)) / totMassDens

!      write(*,*) 'velAvg ', velAvg
      velQuad(:) = mixture%theta_eq*velAvg(:) &
        &        + (1.0_rk - mixture%theta_eq)*vel_inG

      ! moments for multispecies LBM model of transformed variable
      first_moments = mom_inG !momX, momY, momZ
      second_moments = (/ press + rho*velQuad(1)*velQuad(1), & !momXX
        &                 press + rho*velQuad(2)*velQuad(2), & !momYY
        &                 press + rho*velQuad(3)*velQuad(3), & !momZZ
        &                 rho*velQuad(1)*velQuad(2),         & !momXY
        &                 rho*velQuad(2)*velQuad(3),         & !momYZ
        &                 rho*velQuad(3)*velQuad(1) /)         !momXZ
      third_moments = (/ cs2*mom_inG(2), & !momXXY
        &                cs2*mom_inG(3), & !momXXZ
        &                cs2*mom_inG(1), & !momYYX
        &                cs2*mom_inG(3), & !momYYZ
        &                cs2*mom_inG(1), & !momZZX
        &                cs2*mom_inG(2) /) !momZZY
      fourth_moments = (/ &
        & cs2*(press+rho*(        velQuad(1)**2 +        velQuad(2)**2 &
        &                - 0.5_rk*velQuad(3)**2 )), &!momXXYY
        & cs2*(press+rho*(-0.5_rk*velQuad(1)**2 +        velQuad(2)**2 &
        &                + velQuad(3)**2 )), &!momYYZZ
        & cs2*(press+rho*(        velQuad(1)**2 - 0.5_rk*velQuad(2)**2 &
        &                + velQuad(3)**2 ))  &!momZZXX
        & /)

      ! moments for liquid mixture
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
        state(                                                         &
&neigh (( missing_links(ilink)-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
        & )  = unKnown_fTmp(iLink)
      end do

     deallocate(missing_links)
     deallocate(unKnown_fTmp)
     deallocate(rhs)
   end do

  end subroutine spc_moments_moleFrac
! ****************************************************************************** !


! ****************************************************************************** !
  !> molar flux boundary condition like moments velocity bc type
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moments_moleflux',
  !!    moleflux = 0.0
  !!     }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moments_moleFlux( me, state, bcBuffer, globBC, levelDesc,     &
    &                              tree, nSize, iLevel, sim_time, neigh,       &
    &                              layout, fieldProp, varPos, nScalars,        &
    &                              varSys, derVarPos, physics, iField, mixture )
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
    real(kind=rk) :: moleFlux(globBC%nElems(iLevel)*3), inv_flux
    integer :: iElem, iDir, iFieldLoc, nFields, pos, tmpDir
    real(kind=rk) :: rho, press, dens_correction
    integer :: nLinks, iLink, elemPos, QQ
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp(:)
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*varSys%nStateVars )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: massFlux(3), vel_inG(3)
    real(kind=rk) :: mom(3,varSys%nStateVars), mom_inG(3)
    real(kind=rk) :: mass_dens(varSys%nStateVars), num_dens(varSys%nStateVars)
    real(kind=rk) :: totMassDens, velAvg(3), velQuad(3)
    real(kind=rk) :: moleFrac(varSys%nStateVars)
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: molWeight(varSys%nStateVars), phi(varSys%nStateVars)
    real(kind=rk) :: molWeightInv
    real(kind=rk) :: first_moments(3), second_moments(6), third_moments(6), &
      &              fourth_moments(3)
    integer :: nMom(4)
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: offset, bcMoleFlux_pos, posInBuffer
    ! ------------------------------------------------------------------------
!write(dbgUnit(1),*) 'Boundary label ', trim(me%label)
!write(dbgUnit(1),*) 'iField ', iField

    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars
    molWeightInv = fieldProp%species%molWeightInv

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )

    do iFieldLoc = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                        %fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                  %fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFieldLoc, :) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                            %fieldProp%species%resi_coeff(:)
    end do

    inv_flux = 1.0_rk / physics%fac( iLevel )%moleFlux
    ! position of boundary moleflux in varSys
    bcMoleFlux_pos = me%bc_states%moleFlux%varPos
    ! Get moleflux
    call varSys%method%val(bcMoleFlux_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleFlux                    &
      &           %pntIndex%indexLvl(iLevel)               &
      &           %val(1:globBC%nElems(iLevel)),           &
      & nVals   = globBC%nElems(iLevel),                   &
      & res     = moleFlux                                 )

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
    moleFlux = moleFlux * inv_flux
!write(*,*) 'physical moleFlux ', physics%fac( iLevel )%moleFlux
!write(*,*) 'moleFlux ', moleFlux
!stop

    nMom(1) = 1 + size(layout%moment%first_moments)
    nMom(2) = nMom(1) + size(layout%moment%second_moments)
    nMom(3) = nMom(2) + size(layout%moment%third_moments)
    nMom(4) = nMom(3) + size(layout%moment%fourth_moments)

!    write(*,*) 'iField ', iField
    ! Calculate the density of current element
    do iElem = 1, globBC%nElems(iLevel)
!write(dbgUnit(1),*) 'iElem ', iElem

      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(rhs(nLinks))
      allocate(unKnown_fTmp(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )

      iLink = 0
      fTmp = 0.0_rk
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemlvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink +1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = &
            &state(neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = &
&state(neigh (( layout%fstencil%restposition-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)
      !write(*,*) 'fTmp ', fTmp

      rho = fTmp(layout%fStencil%restPosition)
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          tmpDir = layout%fStencil%cxDirInv(iDir)
          if( globBC%elemLvl(iLevel)%bitmask%val( tmpDir, iElem )) then
            tmpDir = layout%fStencil                                      &
              &          %cxDirInv(globBC%elemLvl(iLevel)%normalInd%val(iElem))
          end if
        else
          tmpDir = iDir
        end if
        rho  = rho + fTmp(tmpDir)
      end do
!write(dbgUnit(1),*) 'density ', rho

      !caulate mass flux from moleFlux
      ! massFlux = moleflux * molWeight
      offset = (iElem-1)*3
      massFlux(1) = moleFlux(offset+1) * molWeightInv
      massFlux(2) = moleFlux(offset+2) * molWeightInv
      massFlux(3) = moleFlux(offset+3) * molWeightInv
!write(dbgUnit(1),*) 'massFlux ', massFlux

      ! add small correction term to density which comes from momentum in the
      ! flow direction while substituing unknowns pdfs terms obtained
      ! from moments BC derivation to compute density.
      ! With this correction term density is recovered correctly at this node
      dens_correction = dot_product(massFlux,                            &
        &                           globBC%elemLvl(iLevel)%normal%val(:, iElem))
      ! correction term for Multispecies when flux is defined
      ! rho_new = rho_old + dens_correction
      rho = rho + dens_correction
      !write(*,*) 'density after correction', rho

      !pressure
      press = cs2*rho*fieldProp%species%molWeigRatio
!write(dbgUnit(1),*) 'rho ', rho, 'press', press

      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
      ! local state vector to compute density and velocity of other species
      do iFieldLoc = 1, varSys%nStateVars
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all( pos ) = bcBuffer(                                  &
            & ( posinbuffer-1)* nscalars+ pos)
        end do
      end do

      !local density and momentum
      mass_dens = 0.0_rk
      do iFieldLoc = 1, nFields
        if(iField == iFieldLoc) then
          ! massflux = moleflux * molecular_weight
          mom(:, iField) = massFlux
          mass_dens(iField) = rho
          num_dens(iField) = rho * molWeightInv
        else
          ! mass density
          mass_dens(iFieldLoc) = &
            & sum(fTmp_all(varSys%method%val(iFieldLoc)%state_varPos(:)))
          num_dens(iFieldLoc) = mass_dens(iFieldLoc) / molWeight(iFieldLoc)
          !mass flux
          call derVarPos%momFromState( state  = fTmp_all,         &
            &                          iField = iFieldLoc,        &
            &                          nElems = 1,                &
            &                          varSys = varSys,           &
            &                          layout = layout,           &
            &                          res    = mom(:, iFieldLoc) )
        end if
      end do

      ! mass density
      mass_dens = num_dens * molWeight
      ! total mass density
      totMassDens = sum(mass_dens)

      ! mole fraction
      moleFrac = num_dens/sum(num_dens)

      ! compute velocity of transformed variable g. Eq. which is normally
      ! used to compute actual velocity by solving LSE
      mom_inG = mom(:, iField)
      do iFieldLoc = 1, nFields
        mom_inG = mom_inG + ( mixture%omega_diff * 0.5_rk                    &
            &               * resi_coeff(iField, iFieldLoc) * phi(iField)    &
            &               * moleFrac(iFieldLoc) / mixture%paramB )         &
            &               * mom(:, iField)
      end do
      !set non-diagonal part
      do iFieldLoc = 1, nFields
        mom_inG = mom_inG - ( mixture%omega_diff * 0.5_rk                    &
          &                 * resi_coeff(iField, iFieldLoc) * phi(iFieldLoc) &
          &                 * moleFrac(iField) / mixture%paramB )            &
          &                 * mom(:,iFieldLoc)
      end do

      vel_inG = mom_inG/rho
!write(dbgUnit(1),*) 'defined momentum ', massFlux
!write(dbgUnit(1),*) 'transformed momentum ', mom_inG
!write(dbgUnit(1),*) 'transformed velocity ', vel_inG

      ! mass averaged mixture velocity
      velAvg(1) = sum(mom(1,:)) / totMassDens
      velAvg(2) = sum(mom(2,:)) / totMassDens
      velAvg(3) = sum(mom(3,:)) / totMassDens
!      write(*,*) 'velAvg ', velAvg
      velQuad(:) = mixture%theta_eq*velAvg(:) &
        &        + (1.0_rk - mixture%theta_eq)*vel_inG

      ! moments for multispecies LBM model of transformed variable
      first_moments = mom_inG !momX, momY, momZ
      second_moments = (/ press + rho*velQuad(1)*velQuad(1), & !momXX
        &                 press + rho*velQuad(2)*velQuad(2), & !momYY
        &                 press + rho*velQuad(3)*velQuad(3), & !momZZ
        &                 rho*velQuad(1)*velQuad(2),         & !momXY
        &                 rho*velQuad(2)*velQuad(3),         & !momYZ
        &                 rho*velQuad(3)*velQuad(1) /)         !momXZ
      third_moments = (/ cs2*mom_inG(2), & !momXXY
        &                cs2*mom_inG(3), & !momXXZ
        &                cs2*mom_inG(1), & !momYYX
        &                cs2*mom_inG(3), & !momYYZ
        &                cs2*mom_inG(1), & !momZZX
        &                cs2*mom_inG(2) /) !momZZY
      fourth_moments = (/ &
        & cs2*(press+rho*(        velQuad(1)**2 +        velQuad(2)**2 &
        &                - 0.5_rk*velQuad(3)**2 )), &!momXXYY
        & cs2*(press+rho*(-0.5_rk*velQuad(1)**2 +        velQuad(2)**2 &
        &                + velQuad(3)**2 )), &!momYYZZ
        & cs2*(press+rho*(        velQuad(1)**2 - 0.5_rk*velQuad(2)**2 &
        &                + velQuad(3)**2 ))  &!momZZXX
        & /)

      ! moments for liquid mixture
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

      !write(*,*) 'Moments ', moments
      !write(*,*) 'knownMomPos ', me%elemLvl(iLevel)%moments(iElem)%knownMom_pos

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
        state(                                                                 &
&neigh (( missing_links(ilink)-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
         & )  = unKnown_fTmp(iLink)
      end do

      deallocate(missing_links)
      deallocate(rhs)
      deallocate(unKnown_fTmp)
    end do

  end subroutine spc_moments_moleFlux
! ****************************************************************************** !


! ****************************************************************************** !
  !> velocity boundary condition like moments velocity bc type
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moments_vel',
  !!    moleflux = 0.0
  !!     }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moments_vel( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                         nSize, iLevel, sim_time, neigh, layout,         &
    &                         fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                         physics, iField, mixture                        )
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
    real(kind=rk) :: vel_b(globBC%nElems(iLevel)*3), inv_vel
    integer :: iElem, iDir, iFieldLoc, nFields, pos, tmpDir
    real(kind=rk) :: rho, press, dens_correction
    integer :: nLinks, iLink, elemPos, QQ
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp(:)
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*varSys%nStateVars )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: vel(3), vel_inG(3)
    real(kind=rk) :: mom(3,varSys%nStateVars), mom_inG(3)
    real(kind=rk) :: mass_dens(varSys%nStateVars), num_dens(varSys%nStateVars)
    real(kind=rk) :: velAvg(3), velQuad(3)
    real(kind=rk) :: totMassDens
    real(kind=rk) :: moleFrac(varSys%nStateVars)
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: molWeight(varSys%nStateVars), phi(varSys%nStateVars)
    real(kind=rk) :: first_moments(3), second_moments(6), third_moments(6), &
      &              fourth_moments(3)
    integer :: nMom(4)
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: bcVel_pos, posInBuffer
    ! ------------------------------------------------------------------------
    !write(*,*) 'Boundary label ', trim(me%label)

    !write(*,*) 'iField ', iField

    QQ = layout%fStencil%QQ
    nFields = varSys%nStateVars
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )

    do iFieldLoc = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                        %fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFieldLoc) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                  %fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFieldLoc, :) = fPtr%solverData%scheme%field(iFieldLoc) &
        &                            %fieldProp%species%resi_coeff(:)
    end do

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

    nMom(1) = 1 + size(layout%moment%first_moments)
    nMom(2) = nMom(1) + size(layout%moment%second_moments)
    nMom(3) = nMom(2) + size(layout%moment%third_moments)
    nMom(4) = nMom(3) + size(layout%moment%fourth_moments)

!KM!write(dbgUnit(1),*) 'iField ', iField
    do iElem = 1, globBC%nElems(iLevel)
!KM!write(dbgUnit(1),*) 'iElem ', iElem
!KM!write(dbgUnit(1),*) 'vel ', ux(iElem), uy(iElem), uz(iElem)
      vel = vel_b( (iElem-1)*3+1 : iElem*3 )

      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(rhs(nLinks))
      allocate(unKnown_fTmp(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )

      iLink = 0
      fTmp = 0.0_rk
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink +1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = state(                                           &
            & neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = state(                            &
& neigh (( layout%fstencil%restposition-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)

      ! Calculate the density of current element
      rho = fTmp(layout%fStencil%restPosition)
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          tmpDir = layout%fStencil%cxDirInv(iDir)
          if( globBC%elemLvl(iLevel)%bitmask%val( tmpDir, iElem )) then
            tmpDir = layout%fStencil                                      &
              &          %cxDirInv(globBC%elemLvl(iLevel)%normalInd%val(iElem))
          end if
        else
          tmpDir = iDir
        end if
        rho  = rho + fTmp(tmpDir)
      end do
!KM!write(dbgUnit(1),*) 'density ', rho

      ! add small correction term to density which comes from momentum in the
      ! flow direction while substituing unknowns pdfs terms obtained
      ! from moments BC derivation to compute density.
      ! With this correction term density is recovered correctly at this node
      dens_correction = dot_product(vel, globBC%elemLvl(iLevel)%normal%val(:, iElem))
      ! correction term for Multispecies when vel is defined
      ! rho_new = rho_old/(1.0 - dens_correction)
      rho = rho/(1.0_rk-dens_correction)

!KM!write(dbgUnit(1),*) 'density after correction', rho
      !pressure
      press = cs2*rho*fieldProp%species%molWeigRatio
!KM!write(dbgUnit(1),*) 'press ', press

      ! local state vector to compute density and velocity of other species
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
      do iFieldLoc = 1, varSys%nStateVars
        do iDir = 1,layout%fStencil%QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all( pos ) = bcBuffer(                                  &
            & ( posinbuffer-1)* nscalars+ pos)
        end do
      end do

      !local density and momentum
      mass_dens = 0.0_rk
      do iFieldLoc = 1, nFields
        if(iField == iFieldLoc) then
          ! massflux = moleflux * molecular_weight
          mom(:, iField) = rho * vel
          mass_dens(iField) = rho
          num_dens(iField) = rho / fieldProp%species%molWeight
        else
          mass_dens(iFieldLoc) =                                         &
            & sum(fTmp_all(varSys%method%val(iFieldLoc)%state_varPos(:)))
          num_dens(iFieldLoc) = mass_dens(iFieldLoc) &
            &                 / molWeight(iFieldLoc)

          ! mass flux
          call derVarPos%momFromState( state  = fTmp_all,         &
            &                          iField = iFieldLoc,        &
            &                          nElems = 1,                &
            &                          varSys = varSys,           &
            &                          layout = layout,           &
            &                          res    = mom(:, iFieldLoc) )
        end if
!KM!write(dbgUnit(1),*) iFieldLoc, ' mom ', mom(:, iFieldLoc)
      end do

      ! total mass density
      totMassDens = sum(mass_dens)

      ! mole fraction
      moleFrac = num_dens/sum(num_dens)

!KM!do iFieldLoc = 1, nFields
!KM!  write(dbgUnit(1),*) iFieldLoc, 'vel ', mom(:, iFieldLoc)/mass_dens(iFieldLoc), &
!KM!    & 'dens ', mass_dens(iFieldLoc), 'moleFrac ', moleFrac(iFieldLoc)
!KM!end do

      ! compute velocity of transformed variable g. Eq. which is normally
      ! used to compute actual velocity by solving LSE
      mom_inG = mom(:, iField)
      do iFieldLoc = 1, nFields
        mom_inG = mom_inG + ( mixture%omega_diff * 0.5_rk                    &
            &               * resi_coeff(iField, iFieldLoc) * phi(iField)    &
            &               * moleFrac(iFieldLoc) / mixture%paramB )         &
            &               * mom(:, iField)
      end do
      !set non-diagonal part
      do iFieldLoc = 1, nFields
        mom_inG = mom_inG - ( mixture%omega_diff * 0.5_rk                    &
          &                 * resi_coeff(iField, iFieldLoc) * phi(iFieldLoc) &
          &                 * moleFrac(iField) / mixture%paramB )            &
          &                 * mom(:,iFieldLoc)
      end do

      vel_inG = mom_inG/rho
!write(dbgUnit(1),*) 'transformed momentum ', mom_inG
!write(dbgUnit(1),*) 'transformed velocity ', vel_inG

      ! mass averaged mixture velocity
      velAvg(1) = sum(mom(1,:)) / totMassDens
      velAvg(2) = sum(mom(2,:)) / totMassDens
      velAvg(3) = sum(mom(3,:)) / totMassDens

      velQuad(:) = mixture%theta_eq*velAvg(:) &
        &        + (1.0_rk - mixture%theta_eq)*vel_inG
!      write(*,*) 'velAvg ', velAvg

      ! moments for multispecies LBM model of transformed variable
      first_moments = mom_inG !momX, momY, momZ
      second_moments = (/ press + rho*velQuad(1)*velQuad(1), & !momXX
        &                 press + rho*velQuad(2)*velQuad(2), & !momYY
        &                 press + rho*velQuad(3)*velQuad(3), & !momZZ
        &                 rho*velQuad(1)*velQuad(2),         & !momXY
        &                 rho*velQuad(2)*velQuad(3),         & !momYZ
        &                 rho*velQuad(3)*velQuad(1) /)         !momXZ
      third_moments = (/ cs2*mom_inG(2), & !momXXY
        &                cs2*mom_inG(3), & !momXXZ
        &                cs2*mom_inG(1), & !momYYX
        &                cs2*mom_inG(3), & !momYYZ
        &                cs2*mom_inG(1), & !momZZX
        &                cs2*mom_inG(2) /) !momZZY
      fourth_moments = (/ &
        & cs2*(press+rho*(        velQuad(1)**2 +        velQuad(2)**2 &
        &                - 0.5_rk*velQuad(3)**2 )), &!momXXYY
        & cs2*(press+rho*(-0.5_rk*velQuad(1)**2 +        velQuad(2)**2 &
        &                + velQuad(3)**2 )), &!momYYZZ
        & cs2*(press+rho*(        velQuad(1)**2 - 0.5_rk*velQuad(2)**2 &
        &                + velQuad(3)**2 ))  &!momZZXX
        & /)

      ! moments for liquid mixture
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

      !write(*,*) 'Moments ', moments
      !write(*,*) 'knownMomPos ', me%elemLvl(iLevel)%moments(iElem)%knownMom_pos

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
        state(                                                                &
&neigh ((missing_links(ilink)-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
         & = unKnown_fTmp(iLink)
      end do

      deallocate(missing_links)
      deallocate(rhs)
      deallocate(unKnown_fTmp)
    end do

  end subroutine spc_moments_vel
! ****************************************************************************** !

! **************************************************************************** !
  !> This routine returns mass density of all species and mass averaged mixture
  !! velocity from given pdf of all species for single element.
  !! It is used in Nonequilbrium extrapolation based boundary conditions.
  subroutine calcDensAndVelsFromPDF( nFields, iField, mass_dens, num_dens,     &
    &                                velAvg, eqVel, varSys, pdf, stencil,      &
    &                                species, resi_coeff, omega_fac, paramBInv )
    ! --------------------------------------------------------------------------
    !> Number of fields
    integer, intent(in) :: nFields
    !> current field to compute eqVel
    integer, intent(in) :: iField
    !> Mass density of all species
    real(kind=rk), intent(out) :: mass_dens(nFields)
    !> Number density of all species
    real(kind=rk), intent(out) :: num_dens(nFields)
    !> Mass averaged velocity
    real(kind=rk), intent(out) :: velAvg(3)
    !> Variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> State array of single element for all species
    real(kind=rk), intent(in) :: pdf(varSys%nScalars)
    !> Stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    ! Equilibrium velocity of current species
    real(kind=rk), intent(out) :: eqVel(3)
    !> Species information
    type(mus_species_type), intent(in) :: species(nFields)
    !> Resisitivity coefficient matrix
    real(kind=rk), intent(in) :: resi_coeff(nFields, nFields)
    !> relaxation parameter, omega_diff*0.5_rk
    real(kind=rk), intent(in) :: omega_fac
    !> inverse of free parameter B
    real(kind=rk), intent(in) :: paramBInv
    ! --------------------------------------------------------------------------
    integer :: iFld, iDir
    real(kind=rk) :: fTmp(stencil%QQ)
    real(kind=rk) :: tot_numDensInv, tot_massDensInv
    real(kind=rk) :: momentum(3, nFields), first_moments(3, nFields)
    real(kind=rk) :: moleFrac(nFields), massFrac(nFields)
    real(kind=rk) :: velocity(nFields, nFields)
    ! --------------------------------------------------------------------------
    do iFld = 1, nFields
      do iDir = 1, stencil%QQ
        ! State array of iFld
        fTmp(iDir) = pdf( varSys%method%val(iFld)%state_varPos(iDir) )
      end do
      ! species density
      mass_dens(iFld) = sum(fTmp)
      ! number density
      num_dens(iFld) = mass_dens(iFld) * species(iFld)%molWeightInv
      ! first moments: momentum of transformed PDF
      first_moments(1, iFld) = sum( fTmp * stencil%cxDirRK(1, :) )
      first_moments(2, iFld) = sum( fTmp * stencil%cxDirRK(2, :) )
      first_moments(3, iFld) = sum( fTmp * stencil%cxDirRK(3, :) )
    end do

    !total mass density
    tot_massDensInv = 1.0_rk / sum(mass_dens)
    ! mass fraction
    massFrac(:) = mass_dens(:) * tot_massDensInv
    !total number density
    tot_numDensInv = 1.0_rk / sum(num_dens)
    ! mole fraction
    moleFrac(:) = num_dens(:) * tot_numDensInv

    ! actual momentum of all species obtained from solving
    ! Linear Equation system
    momentum = momentumFromMacroLSE(                      &
      &          moleFraction  = moleFrac,                &
      &          first_moments = first_moments,           &
      &          nFields       = nFields,                 &
      &          phi           = species(:)%molWeigRatio, &
      &          resi_coeff    = resi_coeff,              &
      &          omega_fac     = omega_fac,               &
      &          paramBInv     = paramBInv                )

    ! velocity of all species on current fluid node
    do iFld = 1, nFields
      velocity(:, iFld) = momentum(:, iFld) / mass_dens(iFld)
    end do

    ! equilibrium velocity for current species
    eqVel(:) = velocity(:, iField)
    do iFld = 1, nFields
      eqVel(:) = eqVel(:) + resi_coeff(iField, iFld)           &
        &      * species(iField)%molWeigRatio * moleFrac(iFld) &
        &      * ( velocity(:, iFld) - velocity(:, iField) )   &
        &      * paramBInv
    end do

    ! Mass averaged mixture velocity
    velAvg(1) = dot_product( massFrac, velocity(1, :) )
    velAvg(2) = dot_product( massFrac, velocity(2, :) )
    velAvg(3) = dot_product( massFrac, velocity(3, :) )

  end subroutine calcDensAndVelsFromPDF
! **************************************************************************** !

end module mus_bc_species_module
! ****************************************************************************** !
