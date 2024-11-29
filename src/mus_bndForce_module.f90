! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
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
!> This module contains routines to bnd_force type and routines to initialize
!! bndForce array and compute bndForce on all boundary elements
!!
!! author: Kannan Masilamani
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
module mus_bndForce_module
  ! include treelm modules
  use env_module,               only: rk
  use tem_debug_module,         only: dbgUnit
  use tem_logging_module,       only: logUnit
  use tem_aux_module,           only: tem_abort
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_math_module,          only: cross_product3D

  ! include musubi modules
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_scheme_header_module,  only: mus_scheme_header_type
  use mus_scheme_type_module,    only: array2D_type
  use mus_pdf_module,            only: pdf_data_type
  use mus_field_module,          only: mus_field_type
  use mus_bc_header_module,      only: boundary_type, glob_boundary_type
  use mus_physics_module,        only: mus_physics_type, mus_convertFac_type

  implicit none
  private

  public :: mus_init_bndForce
  public :: mus_calcBndForce

contains

  ! ************************************************************************** !
  !> This routine initialize bndForce array assign function pointer to
  !! calculate bndForce
  subroutine mus_init_bndForce(bndForce, bndMoment, bc_prop, schemeHeader, bc)
    !---------------------------------------------------------------------------
    !> Bnd force to allocate
    real(kind=rk), allocatable, intent(out) :: bndForce(:,:)
    !> Bnd moment to allocate
    real(kind=rk), allocatable, intent(out) :: bndMoment(:,:)
    !> Boundary property
    type(tem_bc_prop_type), intent(in) :: bc_prop
    !> scheme header info
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    !> global array boundary type
    type(boundary_type), intent(inout) :: bc(:)
    !---------------------------------------------------------------------------
    integer :: iBnd
    !---------------------------------------------------------------------------
    write(logUnit(1),'(A)') 'Initialize BndForce'

    select case( trim(schemeHeader%kind) )
    case('fluid', 'fluid_incompressible', 'isotherm_acEq')
      ! deallocated for dynamic load balancing
      if (allocated(bndForce)) deallocate(bndForce)
      allocate(bndForce(bc_prop%property%nElems, 3))
      bndForce = 0.0_rk

      if (allocated(bndMoment)) deallocate(bndMoment)
      allocate(bndMoment(bc_prop%property%nElems, 3))
      bndMoment = 0.0_rk

      do iBnd = 1, bc_prop%nBCtypes
        select case (trim(bc(iBnd)%BC_kind))
        case('wall')
          bc(iBnd)%calcBndForce => mus_calcBndForce_wall
        case('wall_libb')
          bc(iBnd)%calcBndForce => mus_calcBndForce_wall_libb
        case('turbulent_wall', 'turbulent_wall_noneq_expol', &
          & 'turbulent_wall_eq')
          bc(iBnd)%calcBndForce => mus_calcBndForce_turbWall
        case default
          bc(iBnd)%calcBndForce => mus_calcBndForce_dummy
        end select
      end do
    case default
      write(logUnit(1),'(A)') 'WARNING: BndForce calculation is not supported' &
        &                    //' for '//trim(schemeHeader%kind)
      do iBnd = 1, bc_prop%nBCtypes
        bc(iBnd)%calcBndForce => mus_calcBndForce_dummy
      end do
    end select

  end subroutine mus_init_bndForce
  ! ************************************************************************** !


  ! -------------------------------------------------------------------------- !
  !> This routine computes force on boundary elements which are used to
  !! compute lift and drag coefficient on boundary
  subroutine mus_calcBndForce( bndForce, bndMoment, posInBndID, nBCs, field, &
    &                          globBC, minLevel, maxLevel, state, pdf,    &
    &                          levelDesc, layout, varSys, physics )
    ! -------------------------------------------------------------------- !
    !> Boundary force on wall boundary
    real(kind=rk), intent(inout) :: bndForce(:,:)
    !> Boundary moment on wall boundary
    real(kind=rk), intent(inout) :: bndMoment(:,:)
    !> Mapping from global tree%treeid to global boundary%boundaryID
    integer, intent(in) :: posInBndID(:)
    !> number of BC
    integer, intent(in) :: nBCs
    !> fluid parameters and properties
    type(mus_field_type), intent(in) :: field(:)
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC(:)
    !> minlevel and maxLevel
    integer, intent(in) :: minLevel, maxLevel
    !> contains global state vector
    type(pdf_data_type), intent(in) :: pdf(minLevel: maxLevel)
    !> state arrays fo current iLevel both now and next
    type(array2D_type), intent(in) :: state(minLevel: maxLevel)
    !> Level Descriptor
    type(tem_leveldesc_type), intent(in) :: levelDesc(minLevel: maxLevel)
    !> scheme layout type
    type(mus_scheme_layout_type), intent(in) ::layout
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> contains physics conversion factors
    type(mus_physics_type), intent(in) :: physics
    ! -------------------------------------------------------------------- !
    integer :: iLevel, iBnd, nFields
    ! -------------------------------------------------------------------- !
    nFields = size( field )
    ! Boundary force calculation is valid only for single field schemes
    ! like fluid and fluid_incompressible
    if ( nBCs > 0 .and. nFields == 1) then
      do iLevel = minLevel, maxLevel
        ! Treat all boundary conditions
        do iBnd = 1, nBCs
          ! Calculate force on wall and wall_libb boundaries
          call field( nFields )%bc( iBnd )%calcBndForce(                  &
            & bndForce   = bndForce,                                      &
            & bndMoment  = bndMoment,                                     &
            & posInBndID = posInBndID,                                    &
            & globBC     = globBC( iBnd ),                                &
            & currState  = state( iLevel )%val( :, pdf( iLevel )%nNext ), &
            & levelDesc  = levelDesc( iLevel ),                           &
            & nSize      = pdf( iLevel )%nSize,                           &
            & iLevel     = iLevel,                                        &
            & neigh      = pdf( iLevel )%neigh,                           &
            & layout     = layout,                                        &
            & nScalars   = varSys%nScalars,                               &
            & phyConvFac = physics%fac(iLevel)                            )
        end do ! iBnd
      end do ! iLevel
    end if ! nBCs>0


  end subroutine mus_calcBndForce
  ! -------------------------------------------------------------------------- !


  ! ************************************************************************** !
  !> Dummy routine for calcBndForce
  subroutine mus_calcBndForce_dummy( me, bndForce, bndMoment, posInBndID, &
    &                                globBC, currState, levelDesc, nSize, &
    &                                iLevel, neigh, layout, nScalars,     &
    &                                phyConvFac)
    ! --------------------------------------------------------------------- !
    !> field boundary type
    class( boundary_type ), intent(in) :: me
    !> bndForce to fill
    real(kind=rk), intent(inout) :: bndForce(:,:)
    !> Boundary moment on wall boundary
    real(kind=rk), intent(inout) :: bndMoment(:,:)
    ! position of boundary element in boundary%bcID
    integer, intent(in) :: posInBndID(:)
    !> scheme global boundary type
    type( glob_boundary_type ), intent(in) :: globBC
    !> current state array to access post-collision values
    real(kind=rk), intent(in) :: currState(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> level which invokes boundary
    integer,intent(in) :: iLevel
    !> global parameters
    ! type(mus_param_type),intent(in) :: params
    integer,intent(in)  :: neigh(:)  !< connectivity array
    !> scheme layout
    type( mus_scheme_layout_type ),intent(in) :: layout
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> physics conversion factor
    type(mus_convertFac_type), intent(in) :: phyConvFac
    ! --------------------------------------------------------------------- !
    !call tem_abort('Dummy routine for calcBndForce')
  end subroutine mus_calcBndForce_dummy
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine computes bndForce on wall boundary elements
  subroutine mus_calcBndForce_wall( me, bndForce, bndMoment, posInBndID, &
    &                               globBC, currState, levelDesc, nSize, &
    &                               iLevel, neigh, layout, nScalars,     &
    &                               phyConvFac)
    ! --------------------------------------------------------------------- !
    !> field boundary type
    class( boundary_type ), intent(in) :: me
    !> bndForce to fill
    real(kind=rk), intent(inout) :: bndForce(:,:)
    !> Boundary moment on wall boundary
    real(kind=rk), intent(inout) :: bndMoment(:,:)
    ! position of boundary element in boundary%bcID
    integer, intent(in) :: posInBndID(:)
    !> scheme global boundary type
    type( glob_boundary_type ), intent(in) :: globBC
    !> current state array to access post-collision values
    real(kind=rk), intent(in) :: currState(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> level which invokes boundary
    integer,intent(in) :: iLevel
    !> global parameters
    ! type(mus_param_type),intent(in) :: params
    integer,intent(in)  :: neigh(:)  !< connectivity array
    !> scheme layout
    type( mus_scheme_layout_type ),intent(in) :: layout
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> physics conversion factor
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !---------------------------------------------------------------------------
    integer :: elemPos, iElem, iDir, QQN, QQ
    integer :: invDir
    real(kind=rk) :: force(3), fOut
    real(kind=rk) :: delta_x(3), conv2LatLen
    !---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    QQN    = layout%fStencil%QQN
    conv2LatLen = 1.0_rk / phyConvFac%length

    do iElem = 1, globBC%nElems_Fluid(iLevel)
      force = 0.0_rk

      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! arm of the momentum
      delta_x = (levelDesc%baryOfTotal(elemPos, :) - me%bndMomRefPnt(:)) &
        &     * conv2LatLen

      do iDir = 1, QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          ! direction towards boundary
          invDir = layout%fStencil%cxDirInv( iDir )

          ! Post-collision PDF
          fOut = currState(                                          &
            & ( elempos-1)* nscalars+invdir+( 1-1)* qq )
          ! For wall, qVal=0.5 so fIn = fOut
          force = force + layout%fStencil%cxDirRK(:,invDir) * 2.0_rk * fOut
        end if
      end do !iComp

      ! store computed force on bndForce val
      bndForce(posInBndID( levelDesc%pntTID(elemPos) ), :) = force
      ! moment = arm x force
      bndMoment(posInBndID( levelDesc%pntTID(elemPos) ), :) &
        & = cross_product3D(delta_x, force)
    end do

  end subroutine mus_calcBndForce_wall
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine computes bndForce on wall_libb boundary elements
  subroutine mus_calcBndForce_wall_libb( me, bndForce, bndMoment, posInBndID, &
    &                                    globBC, currState, levelDesc, nSize, &
    &                                    iLevel, neigh, layout, nScalars,     &
    &                                    phyConvFac)
    ! --------------------------------------------------------------------- !
    !> field boundary type
    class( boundary_type ), intent(in) :: me
    !> bndForce to fill
    real(kind=rk), intent(inout) :: bndForce(:,:)
    !> Boundary moment on wall boundary
    real(kind=rk), intent(inout) :: bndMoment(:,:)
    ! position of boundary element in boundary%boundaryID
    integer, intent(in) :: posInBndID(:)
    !> scheme global boundary type
    type( glob_boundary_type ), intent(in) :: globBC
    !> current state array to access post-collision values
    real(kind=rk), intent(in) :: currState(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> level which invokes boundary
    integer,intent(in) :: iLevel
    !> global parameters
    ! type(mus_param_type),intent(in) :: params
    integer,intent(in)  :: neigh(:)  !< connectivity array
    !> scheme layout
    type( mus_scheme_layout_type ),intent(in) :: layout
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> physics conversion factor
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !---------------------------------------------------------------------------
    integer :: elemPos, posInBuffer, iElem, iDir, QQN, QQ
    integer :: invDir, iLink
    real(kind=rk) :: force(3), fIn, fOut, fOutNeigh, fIntp
    real(kind=rk) :: delta_x(3), conv2LatLen
    real(kind=rk) :: qVal, cIn, cOut, cNgh
    !---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    QQN    = layout%fStencil%QQN
    iLink = 0

    conv2LatLen = 1.0_rk / phyConvFac%length
    do iElem = 1, globBC%nElems_Fluid(iLevel)
      force = 0.0_rk

      elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )
      posInBuffer = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )

      ! arm of the momentum
      delta_x = (levelDesc%baryOfTotal(elemPos, :) - me%bndMomRefPnt(:)) &
        &     * conv2LatLen

      do iDir = 1, QQN
        if( globBC%elemLvl( iLevel )%bitmask%val( iDir, iElem )) then
          iLink = iLink + 1
          ! direction towards boundary
          invDir = layout%fStencil%cxDirInv( iDir )
          ! qValues
          qVal = globBC%elemLvl( iLevel )%qVal%val( invDir, iElem )

          cIn  = me%bouzidi(iLevel)% cIn( iLink )
          cOut = me%bouzidi(iLevel)%cOut( iLink )
          cNgh = me%bouzidi(iLevel)%cNgh( iLink )

          ! Incoming state after collision
          fIn = currState(                                           &
            & ( elempos-1)* nscalars+idir+( 1-1)* qq )
          ! Outgoing state after collision
          fOut = currState(                                          &
            & ( elempos-1)* nscalars+invdir+( 1-1)* qq )

          ! Outgoing direction of neighbor after collision can be obtained
          ! using FETCH in invDir
          fOutNeigh = currState(                                          &
            &  neigh((invdir-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0 )

          ! Incoming interpolated state according to bouzidi rule
          fIntp = cIn*fIn + cOut*fOut + cNgh*fOutNeigh

!          ! For outgoing direction, post-collision is extracted from
!          ! currState using SAVE macro in invDir since this value is
!          ! not modified by any BC. So no need to use bcBuffer.
!          fOut = currState(                                          &
!            & ( elempos-1)* nscalars+invdir+( 1-1)* qq )
!
          force = force + layout%fStencil%cxDirRK(:, invDir) * (fOut + fIntp)
        end if
      end do !iComp

      ! store computed force on bndForce val
      bndForce(posInBndID( levelDesc%pntTID(elemPos) ), :) = force
      ! moment = arm x force
      bndMoment(posInBndID( levelDesc%pntTID(elemPos) ), :) &
        & = cross_product3D(delta_x, force)
    end do

  end subroutine mus_calcBndForce_wall_libb
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine access bndForce from turbulent wall function boundary elements
  subroutine mus_calcBndForce_turbWall( me, bndForce, bndMoment, posInBndID, &
    &                                   globBC, currState, levelDesc, nSize, &
    &                                   iLevel, neigh, layout, nScalars,     &
    &                                   phyConvFac)
    ! --------------------------------------------------------------------- !
    !> field boundary type
    class( boundary_type ), intent(in) :: me
    !> bndForce to fill
    real(kind=rk), intent(inout) :: bndForce(:,:)
    !> Boundary moment on wall boundary
    real(kind=rk), intent(inout) :: bndMoment(:,:)
    ! position of boundary element in boundary%bcID
    integer, intent(in) :: posInBndID(:)
    !> scheme global boundary type
    type( glob_boundary_type ), intent(in) :: globBC
    !> current state array to access post-collision values
    real(kind=rk), intent(in) :: currState(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> level which invokes boundary
    integer,intent(in) :: iLevel
    !> global parameters
    ! type(mus_param_type),intent(in) :: params
    integer,intent(in)  :: neigh(:)  !< connectivity array
    !> scheme layout
    type( mus_scheme_layout_type ),intent(in) :: layout
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> physics conversion factor
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !---------------------------------------------------------------------------
    integer :: elemPos, iElem
    !---------------------------------------------------------------------------
    do iElem = 1, globBC%nElems_Fluid(iLevel)
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! access computed force on turbWallFunc bndForce
      bndForce(posInBndID( levelDesc%pntTID(elemPos) ), :) &
        & = me%turbWallFunc%dataOnLvl(iLevel)%bndForce(:, iElem)
      bndMoment(posInBndID( levelDesc%pntTID(elemPos) ), :) &
        & = me%turbWallFunc%dataOnLvl(iLevel)%bndMoment(:, iElem)
    end do

  end subroutine mus_calcBndForce_turbWall
  ! ************************************************************************** !


end module mus_bndForce_module
