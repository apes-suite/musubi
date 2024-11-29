! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017-2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!> Boundary condition treatment routines for fluid simulation
!!
!! A detailed description on the implementation details are given in
!! \ref boundary_implementation
!!
module mus_bc_fluid_experimental_module
  use iso_c_binding, only: c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,               only: rk
  use tem_param_module,         only: cs, cs2, cs2inv, rho0
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_debug_module,         only: dbgUnit
  use tem_logging_module,       only: logUnit
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_property_module,      only: prp_hasQVal
  use tem_aux_module,           only: tem_abort

  ! include musubi modules
  use mus_bc_header_module,      only: boundary_type, glob_boundary_type
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_field_prop_module,     only: mus_field_prop_type
  use mus_derVarPos_module,      only: mus_derVarPos_type
  use mus_param_module,          only: mus_param_type
  use mus_physics_module,        only: mus_physics_type
  use mus_mixture_module,        only: mus_mixture_type
  use mus_varSys_module,         only: mus_varSys_data_type


  implicit none

  private

  public :: pressure_expol_slow

  public :: moments_inflow, moments_outflow
  public :: spc_moments_outflow
  public :: spc_bb_wall
  public :: spc_bb_vel_test


contains

! ****************************************************************************** !
  !> species bounce back velocity boundary
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!   kind = 'spc_bb_vel_test',
  !!   velocity = 'inlet_vel',
  !! }
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
  subroutine spc_bb_vel_test( me, state, bcBuffer, globBC, levelDesc, tree,   &
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
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: fEq( layout%fStencil%QQ, varSys%nStateVars )
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*varSys%nStateVars )
    integer :: iELem, iDir, iFieldLoc, nFields, pos, invPos
    integer :: iNeigh, elemPos, fluidNeighbor, QQ
    real(kind=rk) :: mass_dens(varSys%nStateVars), vel(3, varSys%nStateVars)
    real(kind=rk) :: fTmp_allNeigh( layout%fStencil%QQ*varSys%nStateVars )
    real(kind=rk) :: fTmp_Next( layout%fStencil%QQ*varSys%nStateVars )
    !real(kind=rk) :: fEqNeigh( layout%fStencil%QQ )
    real(kind=rk) :: fEqNeigh( layout%fStencil%QQ, varSys%nStateVars )
    real(kind=rk) :: fEqNext( layout%fStencil%QQ )
    real(kind=rk) :: ucx
    real(kind=rk) :: vel_b(globBC%nElems(iLevel)*3), inv_vel
    integer :: posInBuffer, bcVel_pos, offset
    ! ------------------------------------------------------------------------
    write(logUnit(1),*) 'WARNING: Experimental BC spc_bb_vel_test'
!!write(dbgUnit(1),*) 'Boundary label ', trim(me%label)
!!write(dbgUnit(1),*) 'iField ', iField

    nFields = varSys%nStateVars
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

    ! If physical quantities are given, transform to lattice units by division
    ! with the conversion factor
    !write(*,*) 'phy Vel', ux, uy, uz
    vel_b = vel_b * inv_vel

    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
!!write(dbgUnit(1),*) 'iElem ', iElem
      fTmp_all = 0.0_rk
      ! Calculate the density of current element
      do iFieldLoc = 1, nFields
        do iDir = 1, QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all(pos) = bcBuffer( pos + ( posInBuffer-1)*nScalars )
        end do
      end do !iField

      fTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                    (posInBuffer-1)*nScalars+varPos(1)+QQ-1 )

!!write(dbgUnit(1),*) 'fTmp_all ', fTmp_all
      fTmp_Next = 0.0_rk
      do iFieldLoc = 1, nFields
        call derVarPos%equilFromState( state  = fTmp_All,         &
          &                            iField = iFieldLoc,        &
          &                            nElems = 1,                &
          &                            varSys = varSys,           &
          &                            layout = layout,           &
          &                            res    = fEq(:, iFieldLoc) )

        !write(*,*) 'iField ', iFieldLoc, 'local fEq ', fEq
        call tem_abort('Spc_vel_bb_test: bitmast is allcoated with QQN')
        do iDir = 1,layout%fStencil%QQ
          if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem ) .or. &
            & iDir == layout%fStencil%restPosition ) then
            pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
            invPos = varSys%method%val(iFieldLoc)%state_varPos( &
              & layout%fStencil%cxDirInv(iDir) )
            !write(*,*) 'iDir ', iDir, pos, invPos
            fTmp_Next(pos) = ( fTmp_all( invPos )                      &
              &            + mixture%omega_diff * 0.5_rk *             &
              &              fEq( layout%fStencil%cxDirInv(iDir), iFieldLoc ) ) &
              &            / (1.0_rk+mixture%omega_diff*0.5_rk)
          end if
        end do
      end do

!!write(dbgUnit(1),*) 'local fTmp_Next ', fTmp_Next

      ! current element position in state array
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
      offset = (iElem-1)*3
!!write(dbgUnit(1),*) 'elemPos', elemPos, 'treeID ', levelDesc%total(elemPos)

      ! Transform g to f in each neighbor node
      do iNeigh = 1,layout%fStencil%QQN
        ! fluid neighbor in inverse direction of iNeigh (x-c_i)
        ! @todo: fluid neighbor should be found before
        fluidNeighbor = &
 & int((neigh( ( ineigh-1)* nsize + elempos)-1)/ nscalars)+1

        !if (elemPos /= fluidNeighbor) then
        if( .not. globBC%elemLvl(iLevel)%bitmask%val( iNeigh, iElem ) ) then
          if (elemPos == fluidNeighbor) then
            fTmp_allNeigh =  fTmp_all
          else
            ! get pdf of neighbor in iNeigh
            do iFieldLoc = 1, nFields
              do iDir = 1, QQ
                pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
                fTmp_allNeigh(pos) = &
& state( ( fluidneighbor-1)* nscalars+ idir+( ifieldloc-1)* qq )
              end do
            end do !iFieldLoc
          end if

!!write(dbgUnit(1),*) 'iDir ', iNeigh, 'fluidNeighbor ', fluidNeighbor
!!write(dbgUnit(1),*) 'neighID ', levelDesc%total(fluidNeighbor)
!!write(dbgUnit(1),*) 'fTmp_allNeigh ', fTmp_allNeigh

          ! compute equilibrium of neighbor
          fEqNeigh = 0.0_rk
          do iFieldLoc = 1, nFields
            call derVarPos%equilFromState( state  = fTmp_AllNeigh,         &
              &                            iField = iFieldLoc,             &
              &                            nElems = 1,                     &
              &                            varSys = varSys,                &
              &                            layout = layout,                &
              &                            res    = fEqNeigh(:, iFieldLoc) )
          !write(*,*) 'iFieldLoc ', iFieldLoc, 'fEqNeigh ', fEqNeigh(:, iFieldLoc)
          !write(*,*) 'iFieldLoc ', iFieldLoc, 'sum(fEqNeigh) ', sum(fEqNeigh(:, iFieldLoc))
          end do

          do iFieldLoc = 1, nFields
            pos = varSys%method%val(iFieldLoc)%state_varPos(iNeigh)
            fTmp_Next(pos) = ( fTmp_allNeigh(pos) &
              &               + mixture%omega_diff*0.5_rk*fEqNeigh(iNeigh, iFieldLoc) ) &
              &               / (1.0_rk+mixture%omega_diff*0.5_rk)
          !fTmp_Next(varSys%method%val(iFieldLoc)%state_varPos(tmpDir)) = &
          !  & ( fTmp_allNeigh(varSys%method%val(iFieldLoc)%state_varPos(iNeigh)) &
          !  &               + mixture%omega_diff*0.5_rk*fEqNeigh(iNeigh, iFieldLoc) ) &
          !  &               / (1.0_rk+mixture%omega_diff*0.5_rk)
          end do
        end if
      end do !iNeigh
!write(dbgUnit(1),*) 'fTmp_Next ', fTmp_Next

      mass_dens = 0.0_rk
      do iFieldLoc = 1, nFields
        do iDir = 1,QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          mass_dens(iFieldLoc) = mass_dens(iFieldLoc) + fTmp_Next(pos)
        end do
      end do
!write(dbgUnit(1),*) 'mass_dens ', mass_dens
!write(dbgUnit(1),*) 'vel ', ux(iElem), uy(iElem), uz(iElem)

      ! update untransformed f with velocity bounce back
      do iDir = 1,QQ
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem ) ) then
!write(dbgUnit(1),*) 'iDir', iDir
          ucx = layout%fStencil%cxDir(1,iDir) * vel_b(offset+1)   &
            & + layout%fStencil%cxDir(2,iDir) * vel_b(offset+2)   &
            & + layout%fStencil%cxDir(3,iDir) * vel_b(offset+3)
          do iFieldLoc = 1, nFields
            pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
            fTmp_Next(pos) = fTmp_Next(pos) &
              & + layout%weight(iDir) * 2._rk * cs2inv * mass_Dens(iFieldLoc) &
              & * ucx
          end do
        end if
      end do
!write(dbgUnit(1),*) 'after ubb fTmp_Next ', fTmp_Next

      mass_dens = 0.0_rk
      vel = 0.0_rk
      do iFieldLoc = 1, nFields
        do iDir = 1,QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          mass_dens(iFieldLoc) = mass_dens(iFieldLoc) + fTmp_Next(pos)
          vel(:, iFieldLoc)  = vel(:, iFieldLoc) &
            & + real(layout%fStencil%cxDir( :, iDir), kind=rk) * fTmp_Next(pos)
        end do
        !vel(:, iFieldLoc)  = (/ ux(iElem), uy(iElem), uz(iElem) /)
      end do
      do iFieldLoc = 1, nFields
        vel(:, iFieldLoc) = vel(:,iFieldLoc)/mass_dens(iFieldLoc)
      end do

!write(dbgUnit(1),*) 'after ubb mass_dens ', mass_dens
!write(dbgUnit(1),*) 'after ubb vel ', vel

!write(dbgUnit(1),*) 'fEqNext ', fEqNext
      ! Calculate the equilibrium distribution
      call derVarPos%equilFromMacro( density  = mass_dens, &
        &                            velocity = vel,       &
        &                            iField   = iField,    &
        &                            nElems   = 1,         &
        &                            varSys   = varSys,    &
        &                            layout   = layout,    &
        &                            res      = fEqNext    )

!KM!      ! density
!KM!      rho_spc(1) = sum(fTmp_all(varSys%method%val(iField)%state_varPos(:)))
!KM!
!KM!      ! equilibrium velocity
!KM!      call derVarPos%velocitiesFromState( state  = fTmp_Next, &
!KM!        &                                 iField = iField,    &
!KM!        &                                 nElems = 1,         &
!KM!        &                                 varSys = varSys,    &
!KM!        &                                 layout = layout,    &
!KM!        &                                 res    = velocities )
!KM!
!KM!      eqVel = equilVelFromMacro( iField, moleFraction, velocities, nFields, &
!KM!        &                        paramBInv, phi, resi_coeff )
!KM!
      ! set boundary
      do iDir = 1, QQ
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state(                                                               &
& neigh (( idir-1)*nsize+ globbc%elemlvl(ilevel)%elem%val(ielem))+( ifield-1)* qq+ nscalars*0 ) &
            !! simple Bounce Back
            !& = fTmp(layout%fStencil%cxDirInv( iDir ))

            ! using f in next time step
            & = fTmp_Next(varPos(iDir)) &
            & + mixture%omega_diff*0.5_rk*( fTmp_Next(varPos(iDir)) &
            & - fEqNext(iDir) )

            ! Bounce back with correction term for g with fEq in next time step
            !& = fTmp(layout%fStencil%cxDirInv( iDir ))  &
            !& + mixture%omega_diff*0.5_rk*(fEq(layout%fStencil%cxDirInv( iDir ), iField)&
            !& - fEqNext(iDir) )

            !& = fTmp(layout%fStencil%cxDirInv( iDir ))  &
            !& - layout%weight( iDir )&
            !& * mixture%omega_diff*cs2inv*rho_spc(1)                        &
            !&    * ( layout%fStencil%cxDir( 1, iDir )*eqVel(1)            &
            !&    +   layout%fStencil%cxDir( 2, iDir )*eqVel(2)            &
            !&    +   layout%fStencil%cxDir( 3, iDir )*eqVel(3))

        end if
      end do
    end do

!KM!    ! debug output
!KM!    if (iField == nFields) then
!KM!    do iElem = 1, globBC%nElems( iLevel )
!KM!      !ux = 0.0_rk
!KM!      do iFieldLoc = 1, nFields
!KM!        do iDir = 1,layout%fStencil%QQ
!KM!          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
!KM!          fTmp_all(pos) = state(                                                  &
!KM!& &
!KM!        & )
!KM!        end do
!KM!      end do
!KM!      bound(1)=iElem
!KM!      bound(2)=iElem
!KM!      write(dbgUnit(1),*) 'iElem ', iElem
!KM!      write(dbgUnit(1),*) 'fTmp_all ', fTmp_all
!KM!      ! density
!KM!      do iFieldLoc = 1, nFields
!KM!        rhoTmp = sum(fTmp_all(varSys%method%val(iFieldLoc)%state_varPos(:)))
!KM!
!KM!        ! mass flux
!KM!        call derVarPos%momFromState( state  = fTmp_all,  &
!KM!          &                          iField = iFieldLoc, &
!KM!          &                          nElems = 1,         &
!KM!          &                          varSys = varSys,    &
!KM!          &                          layout = layout,    &
!KM!          &                          res    = uxTmp      )
!KM!        write(dbgUnit(1),*) 'iField ', iFieldLoc, 'density ', rhoTmp
!KM!        write(dbgUnit(1),*) 'iField ', iFieldLoc, 'velocity ', uxTmp
!KM!      end do
!KM!    end do
!KM!    end if
  end subroutine spc_bb_vel_test
! ****************************************************************************** !


! ****************************************************************************** !
  !> species bounce back wall boundary
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_bb_wall',
  !!     }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_bb_wall( me, state, bcBuffer, globBC, levelDesc, tree, nSize, &
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
!    real(kind=rk) :: fEq( layout%fStencil%QQ*globBC%nElems(iLevel) )
!    real(kind=rk) :: fEq
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: fEq( layout%fStencil%QQ, varSys%nStateVars )
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*varSys%nStateVars )
!    real(kind=rk) :: eqVel(3)
    integer :: iELem, iDir, iFieldLoc, nFields, pos, invPos, QQ
    integer :: iNeigh, elemPos, fluidNeighbor
    real(kind=rk) :: mass_dens(varSys%nStateVars), vel(3, varSys%nStateVars)
    real(kind=rk) :: fTmp_allNeigh( layout%fStencil%QQ*varSys%nStateVars )
    real(kind=rk) :: fTmp_Next( layout%fStencil%QQ*varSys%nStateVars )
    !real(kind=rk) :: fEqNeigh( layout%fStencil%QQ )
    real(kind=rk) :: fEqNeigh( layout%fStencil%QQ, varSys%nStateVars )
    real(kind=rk) :: fEqNext( layout%fStencil%QQ )
    integer :: posInBuffer
    ! ------------------------------------------------------------------------
    write(logUnit(1),*) 'WARNING: Experimental BC spc_bb_wall'
!write(dbgUnit(1),*) 'Boundary label ', trim(me%label)
!write(dbgUnit(1),*) 'iField ', iField
    nFields = varSys%nStateVars
    QQ = layout%fStencil%QQ

    do iElem = 1, globBC%nElems(iLevel)

      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
!write(dbgUnit(1),*) 'iElem ', iElem
      fTmp_all = 0.0_rk
      ! Calculate the density of current element
      do iFieldLoc = 1, nFields
        do iDir = 1,QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all(pos) = bcBuffer( pos+(posInBuffer-1)*nScalars )
        end do
      end do !iField

      fTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                    (posInBuffer-1)*nScalars+varPos(1)+QQ-1 )

!write(dbgUnit(1),*) 'fTmp_all ', fTmp_all
      fTmp_Next = 0.0_rk
      do iFieldLoc = 1, nFields
        call derVarPos%equilFromState( state  = fTmp_All,         &
          &                            iField = iFieldLoc,        &
          &                            nElems = 1,                &
          &                            varSys = varSys,           &
          &                            layout = layout,           &
          &                            res    = fEq(:, iFieldLoc) )

!write(dbgUnit(1),*) 'iFieldLoc ', iFieldLoc, 'fEqLocal ', fEq(:, iFieldLoc)

        !write(*,*) 'iField ', iFieldLoc, 'local fEq ', fEq
        call tem_abort('Spc_vel_wall: bitmast is allcoated with QQN')
        do iDir = 1,QQ
          if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem ) .or. &
            & iDir == layout%fStencil%restPosition ) then
            pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
            invPos = varSys%method%val(iFieldLoc)%state_varPos( &
              & layout%fStencil%cxDirInv(iDir) )
            !write(*,*) 'iDir ', iDir, pos, invPos
            fTmp_Next(pos) = ( fTmp_all( invPos )                              &
              &            + mixture%omega_diff * 0.5_rk *                     &
              &              fEq( layout%fStencil%cxDirInv(iDir), iFieldLoc ) )&
              &            / (1.0_rk+mixture%omega_diff*0.5_rk)
          end if
        end do
      end do

!write(dbgUnit(1),*) 'local fTmp_Next ', fTmp_Next

      ! current element position in state array
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
!write(dbgUnit(1),*) 'elemPos', elemPos, 'treeID ', levelDesc%total(elemPos)

      ! Transform g to f in each neighbor node
      do iNeigh = 1,layout%fStencil%QQN
        ! fluid neighbor in inverse direction of iNeigh (x-c_i)
        fluidNeighbor = &
 & int((neigh( ( ineigh-1)* nsize + elempos)-1)/ nscalars)+1

        !if (elemPos /= fluidNeighbor) then
        if( .not. globBC%elemLvl(iLevel)%bitmask%val( iNeigh, iElem ) ) then
          if (elemPos == fluidNeighbor) then
            fTmp_allNeigh =  fTmp_all
          else
            ! get pdf of neighbor in iNeigh
            do iFieldLoc = 1, nFields
              do iDir = 1,QQ
                pos = varSys%method%val(iFieldLoc) &
                  &         %state_varPos(iDir)
                fTmp_allNeigh(pos) =  &
& state( ( fluidneighbor-1)* nscalars+ idir+( ifieldloc-1)* qq )
              end do
            end do !iFieldLoc
          end if

!write(dbgUnit(1),*) 'iDir ', iNeigh, 'fluidNeighbor ', fluidNeighbor
!write(dbgUnit(1),*) 'neighID ', levelDesc%total(fluidNeighbor)
!write(dbgUnit(1),*) 'fTmp_allNeigh ', fTmp_allNeigh

          ! compute equilibrium of neighbor
          fEqNeigh = 0.0_rk
          do iFieldLoc = 1, nFields
            call derVarPos%equilFromState( state  = fTmp_AllNeigh,         &
              &                            iField = iFieldLoc,             &
              &                            nElems = 1,                     &
              &                            varSys = varSys,                &
              &                            layout = layout,                &
              &                            res    = fEqNeigh(:, iFieldLoc) )

!write(dbgUnit(1),*) 'iFieldLoc ', iFieldLoc, 'fEqNeigh ', fEqNeigh(:, iFieldLoc)
!write(dbgUnit(1),*) 'iFieldLoc ', iFieldLoc, 'sum(fEqNeigh) ', sum(fEqNeigh(:, iFieldLoc))
          end do

          do iFieldLoc = 1, nFields
            pos = varSys%method%val(iFieldLoc)%state_varPos(iNeigh)
            !fTmp_Next(pos) = fTmp_allNeigh(pos)
            fTmp_Next(pos) = ( fTmp_allNeigh(pos) &
              &               + mixture%omega_diff*0.5_rk*fEqNeigh(iNeigh, iFieldLoc) ) &
              &               / (1.0_rk+mixture%omega_diff*0.5_rk)
          !fTmp_Next(varSys%method%val(iFieldLoc)%state_varPos(tmpDir)) = &
          !  & ( fTmp_allNeigh(varSys%method%val(iFieldLoc)%state_varPos(iNeigh)) &
          !  &               + mixture%omega_diff*0.5_rk*fEqNeigh(iNeigh, iFieldLoc) ) &
          !  &               / (1.0_rk+mixture%omega_diff*0.5_rk)
          end do
        end if
      end do !iNeigh
!write(dbgUnit(1),*) 'fTmp_Next ', fTmp_Next

      mass_dens = 0.0_rk
      vel = 0.0_rk
      do iFieldLoc = 1, nFields
        do iDir = 1,QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          mass_dens(iFieldLoc) = mass_dens(iFieldLoc) + fTmp_Next(pos)
          vel(:, iFieldLoc)  = vel(:, iFieldLoc) &
            & + real(layout%fStencil%cxDir( :, iDir), kind=rk) * fTmp_Next(pos)
        end do
      end do
      do iFieldLoc = 1, nFields
        vel(:, iFieldLoc) = vel(:,iFieldLoc)/mass_dens(iFieldLoc)
      end do

!write(dbgUnit(1),*) 'mass_dens ', mass_dens
!write(dbgUnit(1),*) 'vel ', vel

      ! Calculate the equilibrium distribution
      call derVarPos%equilFromMacro( density  = mass_dens, &
        &                            velocity = vel,       &
        &                            iField   = iField,    &
        &                            nElems   = 1,         &
        &                            varSys   = varSys,    &
        &                            layout   = layout,    &
        &                            res      = fEqNext    )

!write(dbgUnit(1),*) 'fEqNext ', fEqNext

!KM!      ! density
!KM!      rho_spc(1) = sum(fTmp_all(varSys%method%val(iField)%state_varPos(:)))
!KM!
!KM!      ! equilibrium velocity
!KM!      call derVarPos%velocitiesFromState( state  = fTmp_Next, &
!KM!        &                                 iField = iField,    &
!KM!        &                                 nElems = 1,         &
!KM!        &                                 varSys = varSys,    &
!KM!        &                                 layout = layout,    &
!KM!        &                                 res    = velocities )
!KM!
!KM!      eqVel = equilVelFromMacro( iField, moleFraction, velocities, nFields, &
!KM!        &                        paramBInv, phi, resi_coeff )
!KM!
      ! set boundary
      do iDir = 1,QQ
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          state(                                                               &
& neigh (( idir-1)*nsize+ globbc%elemlvl( ilevel )%elem%val( ielem ))+( ifield-1)* qq+ nscalars*0 ) &
              !! simple Bounce Back
              !& = fTmp(layout%fStencil%cxDirInv( iDir ))

              ! using f in next time step
              & = fTmp_Next(varPos(iDir)) &
              & + mixture%omega_diff*0.5_rk*( fTmp_Next(varPos(iDir)) &
              & - fEqNext(iDir) )

              ! Bounce back with correction term for g with fEq in next time step
              !& = fTmp(layout%fStencil%cxDirInv( iDir ))  &
              !& + mixture%omega_diff*0.5_rk*(fEq(layout%fStencil%cxDirInv( iDir ), iField)&
              !& - fEqNext(iDir) )

              !& = fTmp(layout%fStencil%cxDirInv( iDir ))  &
              !& + mixture%omega_diff*0.5_rk*(fEq(layout%fStencil%cxDirInv( iDir ), iField)&
              !& - fEq(iDir, iField) )

              !& = fTmp(layout%fStencil%cxDirInv( iDir ))  &
              !& - layout%weight( iDir )&
              !& * mixture%omega_diff*cs2inv*rho_spc(1)                        &
              !&    * ( layout%fStencil%cxDir( 1, iDir )*eqVel(1)            &
              !&    +   layout%fStencil%cxDir( 2, iDir )*eqVel(2)            &
              !&    +   layout%fStencil%cxDir( 3, iDir )*eqVel(3))

        end if
      end do
    end do

    ! debug output
!KM!    if (iField == nFields) then
!KM!    do iElem = 1, globBC%nElems( iLevel )
!KM!      !ux = 0.0_rk
!KM!      do iFieldLoc = 1, nFields
!KM!        do iDir = 1,layout%fStencil%QQ
!KM!          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
!KM!          fTmp_all(pos) = state(                                                  &
!KM!& &
!KM!        & )
!KM!        end do
!KM!      end do
!KM!      bound(1)=iElem
!KM!      bound(2)=iElem
!KM!!write(dbgUnit(1),*) 'iElem ', iElem
!KM!!write(dbgUnit(1),*) 'fTmp_all ', fTmp_all
!KM!      ! density
!KM!      do iFieldLoc = 1, nFields
!KM!        rhoTmp = sum(fTmp_all(varSys%method%val(iFieldLoc)%state_varPos(:)))
!KM!
!KM!        ! mass flux
!KM!        call derVarPos%momFromState( state  = fTmp_all,  &
!KM!          &                          iField = iFieldLoc, &
!KM!          &                          nElems = 1,         &
!KM!          &                          varSys = varSys,    &
!KM!          &                          layout = layout,    &
!KM!          &                          res    = uxTmp      )
!KM!!write(dbgUnit(1),*) 'iField ', iFieldLoc, 'density ', rhoTmp
!KM!!write(dbgUnit(1),*) 'iField ', iFieldLoc, 'velocity ', uxTmp
!KM!      end do
!KM!    end do
!KM!    end if
  end subroutine spc_bb_wall
! ****************************************************************************** !


! ****************************************************************************** !
  !> molar flux boundary condition like moments velocity bc type
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !! { label = 'outlet',
  !!    kind = 'spc_moments_outflow',
  !!    moleflux = 0.0
  !!     }
  !!```
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_moments_outflow( me, state, bcBuffer, globBC, levelDesc,      &
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
    integer :: iELem, iDir, iFieldLoc, nFields, pos
    real(kind=rk) :: rho, press
    integer :: nLinks, iLink, elemPos, neighPos, QQ
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp(:)
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*varSys%nStateVars )
    real(kind=rk) :: toMoments( layout%fStencil%QQ )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: vel_inG(3)
    real(kind=rk) :: mom(3,varSys%nStateVars), mom_inG(3)
    real(kind=rk) :: mass_dens(varSys%nStateVars), num_dens(varSys%nStateVars)
    real(kind=rk) :: totMassDens, velAvg(3)
    real(kind=rk) :: moleFrac(varSys%nStateVars)
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: molWeight(varSys%nStateVars), phi(varSys%nStateVars)
    type(mus_varSys_data_type), pointer :: fPtr
    ! ------------------------------------------------------------------------
    write(logUnit(1),*) 'WARNING: Experimental BC spc_moments_outflow'
    !write(*,*) 'Boundary label ', trim(me%label)

    nFields = varSys%nStateVars
    QQ = layout%fStencil%QQ

    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )

    do iFieldLoc = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFieldLoc) =                                      &
        &  fPtr%solverData%scheme%field(iFieldLoc)%fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFieldLoc) = &
        & fPtr%solverData%scheme%field(iFieldLoc)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFieldLoc, :) =                                     &
        & fPtr%solverData%scheme%field(iFieldLoc)%fieldProp%species%resi_coeff(:)
    end do

!    write(*,*) 'iField ', iField
    ! Calculate the density of current element
    do iElem = 1, globBC%nElems(iLevel)
!      write(*,*) 'iElem ', iElem

      fTmp = 0.0_rk
      rho = 0
      nLinks = me%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs
      allocate(missing_links(nLinks))
      allocate(rhs(nLinks))
      allocate(unKnown_fTmp(nLinks))
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      neighPos = me%neigh(iLevel)%posInState(1,iElem)

      iLink = 0
      do iDir = 1,layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          iLink = iLink +1
          missing_links(iLink) = iDir
        else
          fTmp( iDir ) = &
            & state( neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )
        end if
      end do
      fTmp( layout%fStencil%restPosition ) = &
& state( neigh (( layout%fstencil%restposition-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )

      ! local state vector to compute density and velocity of other species
      do iFieldLoc = 1, varSys%nStateVars
        do iDir = 1,QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all( pos ) = &
!& state( )
& state(( neighpos-1)* nscalars+ idir+( ifieldloc-1)* qq )
        end do
      end do

      do iFieldLoc = 1, nFields
        mass_dens(iFieldLoc) =                                          &
          & sum(fTmp_all( varSys%method%val(iFieldLoc)%state_varPos(:) ))
        num_dens(iFieldLoc) = mass_dens(iFieldLoc) &
          &                 / molWeight(iFieldLoc)

        ! mass flux
        call derVarPos%momFromState( state  = fTmp_all,         &
          &                          iField = iFieldLoc,        &
          &                          nElems = 1,                &
          &                          varSys = varSys,           &
          &                          layout = layout,           &
          &                          res    = mom(:, iFieldLoc) )
      end do

      ! mass density
      mass_dens = num_dens * molWeight
      totMassDens = sum(mass_dens)

      rho = mass_dens(iField)
!      write(*,*) 'density ', rho

      ! mole fraction
      moleFrac = num_dens/sum(num_dens)

!      do iFieldLoc = 1, nFields
!        write(*,*) iFieldLoc, 'vel ', vel(:, iFieldLoc), &
!          & 'dens ', mass_dens(iFieldLoc), 'moleFrac ', moleFrac(iFieldLoc)
!      end do

      ! compute velocity of transformed variable g. Eq. which is normally
      ! used to compute actual velocity by solving LSE
      mom_inG = mom(:, iField)
      do iFieldLoc = 1, nFields
        mom_inG = mom_inG + ( mixture%omega_diff * 0.5_rk                    &
            &               * resi_coeff(iField, iFieldLoc) * phi(iField)    &
            &               * moleFrac(iFieldLoc) / mixture%paramB )         &
            &             * mom(:, iField)
      end do
      !set non-diagonal part
      do iFieldLoc = 1, nFields
        mom_inG = mom_inG - ( mixture%omega_diff * 0.5_rk                    &
          &                 * resi_coeff(iField, iFieldLoc) * phi(iFieldLoc) &
          &                 * moleFrac(iField) / mixture%paramB )            &
          &               * mom(:,iFieldLoc)
      end do

      vel_inG = mom_inG/rho
!      write(*,*) 'transformed momentum ', mom_inG

      ! mass averaged mixture velocity
      velAvg(1) = sum(mom(1,:)) / totMassDens
      velAvg(2) = sum(mom(2,:)) / totMassDens
      velAvg(3) = sum(mom(3,:)) / totMassDens

!      write(*,*) 'velAvg ', velAvg

      press = cs2*rho*fieldProp%species%molWeigRatio
      !write(*,*) 'rho ', rho, 'press', press

      toMoments = matmul(layout%moment%toMoments%A, fTmp)
      !write(*,*) 'fTmp ', fTmp
      !write(*,*) 'toMoments ', toMoments

      ! moments for liquid mixture
      moments = (/ rho, mom_inG(1), mom_inG(2), press+rho*velAvg(1)**2,    &
        & press+rho*velAvg(2)**2, rho*velAvg(1)*velAvg(2), cs2*mom_inG(2), &
        & cs2*mom_inG(1), cs2*(press+rho*velAvg(1)**2+rho*velAvg(2)**2) /)

      !write(*,*) 'Moments ', moments
      !write(*,*) 'knownMomPos ', me%elemLvl(iLevel)%moments(iElem)%knownMom_pos

      do iLink = 1, nLinks
        rhs(iLink) = &
          & moments( me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink) ) &
          & - toMoments( me%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink) )
      end do
      !write(*,*) 'rhs ', rhs

      unKnown_fTmp = matmul( &
        & me%elemLvl(iLevel)%moments(iElem)%unKnownPdfs_MatInv, rhs )
      !write(*,*) 'unknown fTmp ', unKnown_fTmp

      do iLink = 1, nLinks
        state(                                                                 &
& neigh (( missing_links(ilink)-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )  = &
&  unKnown_fTmp(iLink)
      end do

      deallocate(missing_links)
      deallocate(rhs)
      deallocate(unKnown_fTmp)
    end do

  end subroutine spc_moments_outflow
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Moment based velocity boundary condition from Sam Bennent PhD thesis
  !! "A Lattice Boltzmann Model for Diffusion of Binary Gas Mixtures"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine moments_inflow( me, state, bcBuffer, globBC, levelDesc, tree,   &
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
    real(kind=rk) :: rho, vel(3), dens_correction
    real(kind=rk) :: vel_b(globBC%nElems(iLevel)*3), &
      &              pressB(globBC%nElems(iLevel))
    integer :: iELem, iDir, tmpDir
    integer :: nLinks, iLink, elemPos, QQ
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp( : )
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    real(kind=rk) :: fTmp_2( layout%fStencil%QQ )
    real(kind=rk) :: rhoTmp, uxTmp(3), inv_vel, inv_press
    integer :: bcVel_pos, bcPress_pos
    ! ---------------------------------------------------------------------------
    write(logUnit(1),*) 'WARNING: Experimental BC moments_inflow'
    !write(*,*) 'Boundary label ', trim(me%label)
    QQ = layout%fStencil%QQ
    inv_vel = 1.0_rk / physics%fac( iLevel )%vel
    inv_press = 1.0_rk / physics%fac( iLevel )%press

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

    ! position of boundary pressure in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! Get pressure
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

    do iElem = 1, globBC%nElems( iLevel )

      write(*,*) 'iElem ', iElem
      vel(:) = vel_b((iElem-1)*3+1:iElem*3)
      write(*,*) 'velocity ', vel

      fTmp = 0.0_rk
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
!           & bcBuffer( ( globbc%elemlvl( ilevel )%posinbcelembuf%val(ielem)-1)* nscalars + varpos(idir))
            &state( neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0)
        end if
!        write(*,*) 'fTmp0 ',  fTmp(iDir)
      end do
      fTmp( layout%fStencil%restPosition ) = &
& state( neigh (( layout%fstencil%restposition-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )

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

write(*,*) 'dens before correction ', rho

      ! add small correction term to density which comes from momentum in the
      ! flow direction while substituing unknowns pdfs terms obtained
      ! from moments BC derivation to compute density.
      ! With this correction term density is recovered correctly at this node
      dens_correction = dot_product(vel, globBC%elemLvl(iLevel)%normal%val(:, iElem))
      ! for compressible model, rho_new = rho_old + rho_new*dens_correction
      ! thus, rho_new = rho_old/(1-dens_correction)
      ! correction term for in-compressible mode
      rho = rho + rho0*dens_correction
write(*,*) 'dens after correction ', rho, rho*cs2
write(*,*) 'dens from press', pressB(iElem)*cs2inv, pressB(iElem)

      rho = pressB(iElem) * cs2inv

      !!! add small correction term to density which comes from momentum in the
      !!! flow direction while substituing unknowns pdfs terms obtained
      !!! from moments BC derivation to compute density.
      !!! With this correction term density is recovered correctly at this node
      !!dens_correction = dot_product(vel, globBC%elemLvl(iLevel)%normal%val(:, iElem))
      !!! correction term for Multispecies
      !!! rho_new = rho_old/(1.0 - dens_correction)
      !!rho = rho + rho0*dens_correction

      !write(*,*) 'fTmp ', fTmp

      ! moments for compressible model
      !moments = (/ rho, rho*ux(iElem), rho*uy(iElem), press+rho*ux(iElem)**2, &
      !  & press+rho*uy(iElem)**2, rho*ux(iElem)*uy(iElem), cs2*rho*uy(iElem), &
      !  & cs2*rho*ux(iElem), cs2*(press+rho*ux(iElem)**2+rho*uy(iElem)**2) /)
      ! moments for in-compressible model
      moments = (/ rho, rho0*vel(1), rho0*vel(2), pressB(iElem)+rho0*vel(1)**2,&
        & pressB(iElem)+rho0*vel(2)**2, rho0*vel(1)*vel(2), cs2*rho0*vel(2),   &
        & cs2*rho0*vel(1), cs2*(pressB(iElem)+rho0*vel(1)**2+rho0*vel(2)**2) /)
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
& neigh(( missing_links(ilink)-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
         & )  = unKnown_fTmp(iLink)
      end do

      deallocate(missing_links)
      deallocate(rhs)
      deallocate(unKnown_fTmp)

    end do

    do iElem = 1, globBC%nElems( iLevel )
      !ux = 0.0_rk
      do iDir = 1,QQ
        fTmp_2(iDir) = state(                                                  &
& neigh (( idir-1)*nsize+ globbc%elemlvl(ilevel)%elem%val( ielem ))+( ifield-1)* qq+ nscalars*0 &
        & )
      end do
      rhoTmp = sum(fTmp_2)
      call derVarPos%velFromState( state  = fTmp_2, &
        &                          iField = iField, &
        &                          nElems = 1,      &
        &                          varSys = varSys, &
        &                          layout = layout, &
        &                          res    = uxTmp   )
      write(*,*) 'iElem ', iElem
      write(*,*) 'density ', rhoTmp, 'press ', cs2*rhoTmp, cs2*rhoTmp*physics%fac( iLevel )%press
      write(*,*) 'velocity ', uxTmp
!      write(*,*) 'velocity ', uxTmp*rhoTmp
    end do
!
!  stop

  end subroutine moments_inflow
! ****************************************************************************** !


! ****************************************************************************** !
  !> author: Kannan Masilamani
  !! Moment based velocity boundary condition from Sam Bennent PhD thesis
  !! "A Lattice Boltzmann Model for Diffusion of Binary Gas Mixtures"
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine moments_outflow( me, state, bcBuffer, globBC, levelDesc, tree,   &
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
    real(kind=rk) :: rho, press, vel(3)
    integer :: iELem, iDir, QQ
    integer :: nLinks, iLink, elemPos
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: rhs(:)
    real(kind=rk), allocatable :: unKnown_fTmp( : )
    real(kind=rk) :: fTmp( layout%fStencil%QQ )
    real(kind=rk) :: moments( layout%fStencil%QQ )
    ! Vel on frst neighbor element
    real(kind=rk) :: uxB_1(3*globBC%nElems(iLevel))
    ! density on first neighbor element
    real(kind=rk) :: rho_1(globBC%nElems(iLevel))
    ! ---------------------------------------------------------------------------
    write(logUnit(1),*) 'WARNING: Experimental BC moments_outflow'
    !write(*,*) 'Boundary label ', trim(me%label)

    QQ = layout%fStencil%QQ

    ! velocity is taken from precollision of current
    ! @todo: what layout is used inside this routine
    call derVarPos%velFromState(                            &
      &     state  = me%neigh(iLevel)%neighBufferPre_nNext(1, :), &
      &     iField = iField,                                &
      &     nElems = globBC%nElems(iLevel),                 &
      &     varSys = varSys,                                &
      &     layout = layout,                                &
      &     res    = uxB_1                                  )

    rho_1 = 0.0_rk
    ! density is taken from precollision of current
    do iElem = 1, globBC%nElems( iLevel )
      do iDir = 1,QQ
        rho_1 = rho_1 + me%neigh(iLevel)%neighBufferPre_nNext( 1, (iElem-1)*QQ+iDir )
      end do
    end do

    do iElem = 1, globBC%nElems( iLevel )

      !write(*,*) 'iElem ', iElem
      rho = rho_1(iElem)
      vel = (/ uxB_1((iElem-1)*3+1), uxB_1((iElem-1)*3+2), uxB_1((iElem-1)*3+3) /)
      !write(*,*) 'dens', rho, 'velocity ', vel

      fTmp = 0.0_rk
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
          fTmp(iDir) = state( neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0)
        end if
!        write(*,*) 'fTmp0 ',  fTmp(iDir)
      end do
      fTmp( layout%fStencil%restPosition ) = &
& state( neigh (( layout%fstencil%restposition-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0 )

      press = rho * cs2

      !write(*,*) 'fTmp ', fTmp

      ! moments for compressible model
      !moments = (/ rho, rho*ux(iElem), rho*uy(iElem), press+rho*ux(iElem)**2, &
      !  & press+rho*uy(iElem)**2, rho*ux(iElem)*uy(iElem), cs2*rho*uy(iElem), &
      !  & cs2*rho*ux(iElem), cs2*(press+rho*ux(iElem)**2+rho*uy(iElem)**2) /)
      ! moments for in-compressible model
      moments = (/ rho, rho0*vel(1), rho0*vel(2), press+rho0*vel(1)**2, &
        & press+rho0*vel(2)**2, rho0*vel(1)*vel(2), cs2*rho0*vel(2), &
        & cs2*rho0*vel(1), cs2*(press+rho0*vel(1)**2+rho0*vel(2)**2) /)
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
& neigh(( missing_links(ilink)-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
         & )  = unKnown_fTmp(iLink)
      end do

      deallocate(missing_links)
      deallocate(rhs)
      deallocate(unKnown_fTmp)

    end do

!!    do iElem = 1, globBC%nElems( iLevel )
!!      !ux = 0.0_rk
!!      do iDir = 1,QQ
!!        fTmp_2(iDir) = state(                                                  &
!!&nelems (( idir-1)* neigh+ globbc%elemlvl(ilevel)%elem%val( ielem ))+( ifield-1)* qq+ nscalars*0 &
!!        & )
!!      end do
!!      rhoTmp = sum(fTmp_2)
!!      call derVarPos%velFromState( state  = fTmp_2, &
!!        &                          iField = iField, &
!!        &                          nElems = 1,      &
!!        &                          varSys = varSys, &
!!        &                          layout = layout, &
!!        &                          res    = uxTmp   )
!!      write(*,*) 'iElem ', iElem
!!      write(*,*) 'density ', rhoTmp, cs2*rhoTmp
!!!      write(*,*) 'velocity ', uxTmp
!!!      write(*,*) 'velocity ', uxTmp*rhoTmp
!!    end do
!
!  stop

  end subroutine moments_outflow
! ****************************************************************************** !

! ****************************************************************************** !
  !> No comment yet!
  !!
  !! @TODO add comment
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine pressure_expol_slow( me, state, bcBuffer, globBC, levelDesc,      &
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
    real(kind=rk) :: fPreCol_prev
    real(kind=rk) :: fTmp_1 ! first  neighbor pdf
    real(kind=rk) :: fTmp_2 ! second neighbor pdf
    real(kind=rk) :: fEq(  layout%fStencil%QQ * globBC%nElems(iLevel) )
    real(kind=rk) :: fEq0( layout%fStencil%QQ * globBC%nElems(iLevel) )
    real(kind=rk) :: rho( globBC%nElems(iLevel) ), inv_rho_phy
    real(kind=rk) :: velocity( 3, globBC%nElems(iLevel) )
    integer :: iDir, iElem, QQ, elemPos, invDir
    integer :: iLink, bcPress_pos, ii, ff
    logical :: axisNormal
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos, vel_pos(3), elemOff
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    dens_pos = varSys%method%val(derVarPos%density)%auxField_varPos(1)
    vel_Pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)

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
        &                            iField   = ff,                    &
        &                            nElems   = globBC%nElems(iLevel), &
        &                            varSys   = varSys,                &
        &                            layout   = layout,                &
        &                            res      = fEq0                   )
    ! Here ends the copied part of pressure_expol

    ! linkwise treatment for links that have boundary
do ii = 1, QQ * 4
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
        state( me%links(iLevel)%val(iLink) ) = fEq0( me%outletExpol(iLevel) &
          &                                            %statePos(iLink)     )
      else
        fTmp_1 = me%neigh(iLevel)%neighBufferPre_nNext(1,                  &
          &                         me%outletExpol(iLevel)%statePos(iLink) )
        fTmp_2 = me%neigh(iLevel)%neighBufferPre_nNext(2,                  &
          &                         me%outletExpol(iLevel)%statePos(iLink) )
        state( me%links(iLevel)%val(iLink) ) = 1.5_rk*fTmp_1 - 0.5_rk*fTmp_2

      end if
      state( me%links(iLevel)%val(iLink) ) = &
        & state( me%links(iLevel)%val(iLink) ) + dble(QQ*ii - nScalars*ii)
    end do ! iLink, linkwise treatment ends here
end do

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
        fPreCol_prev = fPtr%solverData%scheme%state(iLevel)%val(          &
          &  neigh (( invdir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0, &
          & fPtr%solverData%scheme%pdf(iLevel)%nNow )
        state(  neigh (( idir-1)* nsize+ elempos)+( ifield-1)* qq+ nscalars*0 &
          & )  = &
          ! KM:we use fEq0 since this equilibrium has been computed from
          ! defined density at boundary
          & fEq0((iElem-1)*QQ+iDir) + fPreCol_prev            &
          &                       -  fEq( (iElem-1)*QQ+invDir )
      end if
    end do ! iElem

  end subroutine pressure_expol_slow
! ****************************************************************************** !

end module mus_bc_fluid_experimental_module
! ****************************************************************************** !
