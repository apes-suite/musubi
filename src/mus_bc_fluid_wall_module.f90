! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2017, 2019-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2015-2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017-2018, 2020 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2020 Jana Gericke <jana.gericke@uni-siegen.de>
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
!> Boundary condition wall treatment routines
!!
!! This module contains higher order wall treatments
!! A detailed description on the implementation details are given
!! in [[tem_bc_module]].
!!
module mus_bc_fluid_wall_module

  ! include treelm modules
  use env_module,               only: rk
  use tem_param_module,         only: cs2inv, cs2, rho0, rho0Inv
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_debug_module,         only: dbgUnit
  use tem_geometry_module,      only: tem_ElemSizeLevel
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
  use mus_varSys_module,          only: mus_varSys_data_type
  use mus_relaxationParam_module, only: mus_viscosity_type

  implicit none

  private

  ! public :: wall_multiReflection
  public :: slip_wall, spc_slip_wall
  public :: wall_libb
  public :: do_nothing

contains

! ****************************************************************************** !
  !> slip-wall boundary condition. Slip defined by a slip factor
  !!
  !! \li Normal velocity,\f$ u_n = 0 \f$
  !! \li Tangential velocity, \f$ \frac{\partial u_t}{\partial n} = 0 \f$
  !! \li Pressure, \f$ \frac{\partial P}{\partial n} = 0 \f$
  !! For slip-wall boundary, the slip factor will be multiplied by the velocity
  !! if slip factor = 1, then it is full/free-slip and if slip factor = 0, then
  !! it is no-slip
  !!
  !! @todo KM: Currently, free-slip boundary works only for axis-parallel planes.
  !!           Need to extend it for arbitrary geometries
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine slip_wall( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
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
    ! defining local variables
    real(kind=rk) :: fTmp( layout%fStencil%QQ*globBC%nElems(iLevel) )
    real(kind=rk) :: vel(3*globBC%nElems(iLevel)) ! Velocity on boundary element
    real(kind=rk) :: velTmp(3), rho
    integer :: iELem, iDir, bndNormalDir, QQ, posInBuffer
    ! ---------------------------------------------------------------------------

    QQ = layout%fStencil%QQ

    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp( (iElem-1)*QQ+1: (iElem-1)*QQ+QQ ) &
        &       = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                   (posInBuffer-1)*nScalars+varPos(1)+QQ-1 )
    end do

    ! Get local velocity
    call derVarPos%velFromState( state  = fTmp ,                 &
      &                          iField = iField,                &
      &                          nElems = globBC%nElems(iLevel), &
      &                          varSys = varSys,                &
      &                          layout = layout,                &
      &                          res    = vel                    )

    do iElem = 1, globBC%nElems(iLevel)
      velTmp = vel((iElem-1)*3+1 : iELem*3) * me%slip_fac
      rho = sum(fTmp( (iElem-1)*QQ+1: (iElem-1)*QQ+QQ ))
      bndNormalDir = layout%fStencil%cxDirInv( globBC%elemLvl( iLevel )%       &
        &                                                normalInd%val( iElem ))
      !write(dbgUnit(1),*) 'bndNormalDir ',  bndNormalDir
      if( abs(layout%fStencil%cxDir( 1, bndNormalDir )) == 1) velTmp(1) = 0.0_rk
      if( abs(layout%fStencil%cxDir( 2, bndNormalDir )) == 1) velTmp(2) = 0.0_rk
      if( abs(layout%fStencil%cxDir( 3, bndNormalDir )) == 1) velTmp(3) = 0.0_rk
      !write(dbgUnit(1),*) 'velTmp ', velTmp

      do iDir = 1, layout%fStencil%QQN
        ! Write the values
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state(                                                               &
& neigh((idir-1)* nsize+ globbc%elemlvl(ilevel)%elem%val(ielem))+( ifield-1)* qq+ nscalars*0)=&
          ! We need to get post-collision pdf in direction
          ! alpha- outgoing direction, which is the inverse direction of bitmask
          ! For PULL this means, get the outgoing one, as this is the one which
          ! will be bounced back
          ! For PUSH this means, get the already bounced back pdf back, so take
          ! the incoming
            & fTmp((ielem-1)*qq+layout%fstencil%cxdirinv(idir))&
            &       - layout%weight( iDir )*6._rk*rho                          &
            &       * ( layout%fStencil%cxDir( 1, layout%fStencil%             &
            &                                       cxDirInv( iDir ))*velTmp(1)&
            &       +   layout%fStencil%cxDir( 2, layout%fStencil%             &
            &                                       cxDirInv( iDir ))*velTmp(2)&
            &       +   layout%fStencil%cxDir( 3, layout%fStencil%             &
            &                                       cxDirInv( iDir ))*velTmp(3))
        end if
      end do
    end do !iElem

  end subroutine slip_wall
! ****************************************************************************** !


! ****************************************************************************** !
  !> slip-wall boundary condition. Slip defined by a slip factor
  !!
  !! * Normal velocity,\( u_n = 0 \)
  !! * Tangential velocity, \( \frac{\partial u_t}{\partial n} = 0 \)
  !! * Pressure, \( \frac{\partial P}{\partial n} = 0 \)
  !!
  !! For slip-wall boundary, the slip factor will be multiplied by the velocity
  !! if slip factor = 1, then it is full/free-slip and if slip factor = 0, then
  !! it is no-slip
  !!
  !! @todo KM: Currently, free-slip boundary works only for axis-parallel planes.
  !!           Need to extend it for arbitrary geometries
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine spc_slip_wall( me, state, bcBuffer, globBC, levelDesc, tree,      &
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
    real(kind=rk) :: fTmp_all( layout%fStencil%QQ*globBC%nElems(iLevel) &
      &              * varSys%nStateVars )
    real(kind=rk) :: fTmp(layout%fStencil%QQ)
    real(kind=rk) :: mom(3*globBC%nElems(iLevel)), momTmp(3)
    integer :: iELem, iDir, bndNormalDir, pos, iFieldLoc, QQ, posInBuffer
    ! ------------------------------------------------------------------------
    QQ = layout%fStencil%QQ

    do iElem = 1, globBC%nElems(iLevel)
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      do iFieldLoc = 1, varSys%nStateVars
        do iDir = 1, QQ
          pos = varSys%method%val(iFieldLoc)%state_varPos(iDir)
          fTmp_all( pos+(iElem-1)*nScalars ) &
            &  = bcBuffer( pos+(posInBuffer-1)*nScalars )
        end do
      end do
    end do

    call derVarPos%momFromState( state  = fTmp_all,                &
      &                          iField = iField,                  &
      &                          nElems = globBC%nElems( iLevel ), &
      &                          varSys = varSys,                  &
      &                          layout = layout,                  &
      &                          res    = mom                      )


    do iElem = 1, globBC%nElems(iLevel)
      if( .not. btest( levelDesc%property(                        &
        &              globBC%elemLvl(iLevel)%elem%val(iElem)), prp_solid))then

      momTmp(1) = mom((iElem-1)*3+1) * me%slip_fac
      momTmp(2) = mom((iElem-1)*3+2) * me%slip_fac
      momTmp(3) = mom((iElem-1)*3+3) * me%slip_fac

      bndNormalDir = layout%fStencil%cxDirInv( globBC%elemLvl( iLevel )%       &
        &                                                normalInd%val( iElem ))
      if( abs(layout%fStencil%cxDir( 1, bndNormalDir )) == 1) momTmp(1) = 0.0_rk
      if( abs(layout%fStencil%cxDir( 2, bndNormalDir )) == 1) momTmp(2) = 0.0_rk
      if( abs(layout%fStencil%cxDir( 3, bndNormalDir )) == 1) momTmp(3) = 0.0_rk

      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                    (posInBuffer-1)*nScalars+varPos(1)+QQ-1 )

      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! Depending on PUSH or pull, use + or - for cxDir, because directions
          ! are inverted
          state(                                                               &
&neigh(( idir-1)* nsize+ globbc%elemlvl(ilevel)%elem%val(ielem))+( ifield-1)* qq+ nscalars*0)&
              & = fTmp(layout%fStencil%cxDirInv( iDir ))  &
              &    + layout%weight( iDir )*2._rk*cs2inv &
              &    * ( layout%fStencil%cxDir( 1, iDir )*momTmp(1) &
              &    +   layout%fStencil%cxDir( 2, iDir )*momTmp(2) &
              &    +   layout%fStencil%cxDir( 3, iDir )*momTmp(3) )
        end if
      end do
      end if
    end do

  end subroutine spc_slip_wall
! ****************************************************************************** !

! ****************************************************************************** !
  !> No comment yet!
  !!
  !! @TODO add comment
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine wall_libb( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
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
    real(kind=rk) :: fIn, fOut, fNgh
    real(kind=rk) :: cIn, cOut, cNgh
    integer :: iLink
    ! ---------------------------------------------------------------------------

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

  end subroutine wall_libb
! ****************************************************************************** !

! ****************************************************************************** !
  !> No comment yet!
  !!
  !! @TODO add comment
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine do_nothing( me, state, bcBuffer, globBC, levelDesc, tree, nSize,  &
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
  end subroutine do_nothing
! ****************************************************************************** !

end module mus_bc_fluid_wall_module
! ****************************************************************************** !
