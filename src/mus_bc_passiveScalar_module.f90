! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2015, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!! boundary_implementation
!!
module mus_bc_passiveScalar_module

  ! include treelm modules
  use env_module,               only: rk
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_param_module,         only: div1_3, div4_3, csInv, cs, cs2, cs2inv

  ! include musubi modules
  use mus_bc_header_module,      only: boundary_type, glob_boundary_type
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_field_prop_module,     only: mus_field_prop_type
  use mus_derVarPos_module,      only: mus_derVarPos_type
  use mus_param_module,          only: mus_param_type
  use mus_physics_module,        only: mus_physics_type
  use mus_mixture_module,        only: mus_mixture_type

  implicit none

  private

  public :: inlet_pasScal
  public :: outlet_pasScal
  public :: pressure_antiBounceBack_pasScal

contains

! ****************************************************************************** !
  !> Inlet boundary conditions for passive scalar transport (Flekkoy).
  !!
  !! All links of the inlet boundary elements are set to 0.0
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine inlet_pasScal( me, state, bcBuffer, globBC, levelDesc, tree,      &
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
    integer :: iElem, iDir, QQ
    integer :: nElems
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nElems = size( neigh ) / QQ

    do iElem = 1, globBC%nElems( iLevel )
      do iDir = 1, layout%fStencil%QQ
        state(                                                                 &
& ( globbc%elemlvl(ilevel)%elem%val(ielem)-1)* nscalars+idir+( ifield-1)* qq)&
          & = 0._rk
      end do
    end do

  end subroutine inlet_pasScal
! ****************************************************************************** !


! ****************************************************************************** !
  !> Outlet boundary conditions for passive scalar transport (Flekkoy).
  !!
  !!  currently no obstacles allowed lx-2 upstream fluid outlet nodes
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine outlet_pasScal( me, state, bcBuffer, globBC, levelDesc, tree,   &
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
    integer :: iElem, iDir, QQ
    integer :: nElems
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    nElems = size( neigh ) / QQ

    do iElem = 1, globBC%nElems( iLevel )
      do iDir = 1, layout%fStencil%QQ
        state(                                                                 &
& ( globbc%elemlvl( ilevel )%elem%val( ielem )-1)* nscalars+ idir+( ifield-1)* qq)  = &
          & ( state(                                                           &
& ( me%neigh( ilevel )%posinstate( 1, ielem )-1)* nscalars+ idir+( ifield-1)* qq)  + &
          &   state(                                                           &
& ( me%neigh( ilevel )%posinstate( 2, ielem )-1)* nscalars+ idir+( ifield-1)* qq)) / &
          &   2.0_rk
      end do
    end do

  end subroutine outlet_pasScal
! ****************************************************************************** !

  ! ****************************************************************************** !
  !> Dirichlet stationary boundary conditions for passive scalar transport in
  !!
  !! Irina Ginzburg (2005), "Generic boundary conditions for lattice Boltzmann
  !! models and their application to advection and anisotropic dispersion
  !! equations.", Advances in Water Resources, Volume 28, Issue 11
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine pressure_antiBounceBack_pasScal( me, state, bcBuffer, globBC, levelDesc, tree,   &
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
    integer :: iElem, iDir, QQ, invDir, elemPos
    real(kind=rk) :: rhoDef(globBC%nElems(iLevel)) ! Density on boundary element
    real(kind=rk) :: inv_rho_phy
    real(kind=rk) :: fTmp( nScalars * globBC%nElems(iLevel) )
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

    ! convert physical pressure into LB density
    rhoDef = rhoDef * inv_rho_phy

    do iElem = 1, globBC%nElems( iLevel )
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
      fTmp( (iElem-1)*nScalars+1: (iElem-1)*nScalars+QQ ) &
        &       = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) : &
        &                   (posInBuffer-1)*nScalars+varPos(QQ)  )
    end do

    do iElem = 1, globBC%nElems( iLevel )
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      do iDir = 1, layout%fStencil%QQN
        if( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem )) then
          invDir = layout%fStencil%cxDirInv(iDir)

          state( neigh (( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) = &
            &   - fTmp( (iElem-1)*nScalars + invDir ) &
            &   + 2._rk*layout%weight( invDir ) * rhoDef(iElem)

        end if ! bitMask
      end do ! iDir
    end do ! iElem

  end subroutine pressure_antiBounceBack_pasScal
! ****************************************************************************** !


end module mus_bc_passiveScalar_module
! ****************************************************************************** !
