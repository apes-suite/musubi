! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012, 2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2016, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Monika Harlacher <monika.harlacher@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> This module contains the definition of geometry type and routines to
!! geometry information like mesh, boundary, immersed_boundary and restart
module mus_geom_module

  ! include treelm modules
  use env_module,              only: rk, long_k, globalMaxLevels
  use tem_global_module,       only: tem_global_mesh_read
  use tem_property_module,     only: prp_hasBnd, prp_hasQVal
  use treelmesh_module,        only: treelmesh_type, load_tem
  use tem_bc_prop_module,      only: tem_bc_prop_type, init_tem_bc_prop, &
    &                                load_tem_bc_qVal
  use tem_geometry_module,     only: tem_build_treeToProp_pointer
  use tem_solveHead_module,    only: tem_solveHead_type
  use tem_timeControl_module,  only: tem_timeControl_start_at_sim
  use tem_restart_module,      only: tem_restart_type, tem_load_restart
  use tem_comm_env_module,     only: tem_comm_env_type
  use tem_simControl_module,   only: tem_simControl_type
  use tem_tools_module,        only: tem_horizontalSpacer
  use tem_logging_module,      only: logUnit

  ! include musubi modules
  use mus_geomIncrHead_module, only: mus_geomIncrHead_type, &
    &                                mus_geomIncrHead_load
  use mus_IBM_module,          only: mus_IBM_globType, mus_load_IBM

  implicit none

  private

  public :: mus_geom_type
  public :: mus_build_posInProp
  public :: mus_load_geom
  public :: mus_load_bc_data

  !> Geometric information and definitions
  type mus_geom_type

    !> tree data type
    type( treelmesh_type ) :: tree

    !> boundary information as stored on disk
    type( tem_bc_prop_type ) :: boundary

    !> The header type containing all the geometry increase information
    type( mus_geomIncrHead_type ), allocatable :: geomIncr(:)

    !> Logical to define whether geometry increase is active or not
    logical :: dynamicGeom = .false.

    !> Tree element position in the boundary_ID( nDir, nElems) in bc_prop_type
    !! it has a size of tree%nElems
    !! How to use:
    !! do iElem = 1, tree%nElems
    !!   posInBndID = posInBndID( iElem )
    !!   ! current element has boundary only if posInBndID>0
    !!   ! else posInBndID = -1
    !!   if (posInBnd > 0 )
    !!     bnd_ID(1:nDir) = bc_prop%boundary_ID( 1:nDir, posInBndID )
    !!   end if
    !! end do
    integer, allocatable :: posInBndID(:)

    !> Tree element position in the qVal( nDir, nElems) in bc_prop_type
    !! it has a size of tree%nElems
    !! How to use:
    !! do iElem = 1, tree%nElems
    !!   posInQVal = posInQVal( iElem )
    !!   ! current element has qVal if posInQVal>0 else posInQVal = -1
    !!   if (posInQVal > 0 )
    !!     qVal(1:nDir) = bc_prop%qVal( 1:nDir, posInQVal )
    !!   end if
    !! end do
    integer, allocatable :: posInQVal(:)

    !> tree element position in level descriptor total list
    !! it has a size of tree%nElems
    !! How to use:
    !! do iElem = 1, tree%nElems
    !!   treeID = tree%treeID( iElem )
    !!   level = tem_levelOf( treeID )
    !!   posInTotal = levelPointer( iElem )
    !!   treeID = LevelDesc( iLevel )%total( posInTotal )
    !! end do
    integer, allocatable :: levelPointer(:)

    !> Boundary element poisition in the levelwise globBC%elemLvl(:)%elem%val
    !! It has a size of geometry%boundary%property%nElems.
    !! It is used in tracking to extract value stored in boundary types.
    !! Hot to use this access normal direction of boundary element:
    !! do iElem = 1, tree%nElems
    !!   level = tem_levelOf( treeID )
    !!   posInBndID = posInBndID(iElem)
    !!   if (posInBndID > 0) then
    !!     BCIDs = bc_prop%boundary_ID(:, posInBndID)
    !!     minBCID = minval(BCIDs, BCIDs > 0)
    !!     posInBcElem = bcLevelPointer(posInBndID)
    !!     normal = globBC%elemLvl(iLevel)%normal(posInBcElem)
    !!   end if
    !! end do
    integer, allocatable :: bcLevelPointer(:)

    !> Minimum bcID for each boundary element.
    !! if a element has more than one boundary then use minBcID which depends
    !! on boundary order in seeder configuration.
    integer, allocatable :: minBcID(:)

    !> immersed boundary data
    type(mus_IBM_globType) :: globIBM

    !> Contains Forces on boundary elements computed using momentum exchange
    !! method. This will be used to derive_bndForce routine to compute force
    !! of certain boundaries.
    !! Forces are stored in level-independent fashion as geometry%boundaryID
    !! loaded from mesh files.
    !! Dim1: geometry%boundaryi%property%nElems
    !! Dim2: 3
    real(kind=rk), allocatable :: bndForce(:,:)

    !> Contains Moments on boundary elements computed using momentum exchange
    !! method. This will be used to derive_bndMoment routine to compute moment
    !! of certain boundaries.
    !! Forces are stored in level-independent fashion as geometry%boundaryID
    !! loaded from mesh files.
    !! Dim1: geometry%boundaryi%property%nElems
    !! Dim2: 3
    real(kind=rk), allocatable :: bndMoment(:,:)
  end type mus_geom_type

contains

  ! ************************************************************************** !
  !> This routine load all geometry related datas like mesh, boundary
  !! and immersed_boundary. Restart is also loaded here because mesh is loaded
  !! in tem_load_restart if restart read is defined.
  subroutine mus_load_geom(me, restart, solverHead, simControl, proc, &
    &                      scaleFactor, initial_balance)
    ! --------------------------------------------------------------------------
    !< contains geometry information which are loaded in this routine
    type(mus_geom_type), intent(out) :: me
    !> contains restart information
    type(tem_restart_type), intent(out) :: restart
    !> contains general description of the solver including flu_state
    type(tem_solveHead_type), intent(inout) :: solverHead
    !> contains simulation time control information
    type(tem_simControl_type), intent(inout) :: simControl
    !> contains MPI communication environment
    type(tem_comm_env_type), intent(in) :: proc
    !> Temporal scaling factor for multilevel mesh
    integer, intent(in) :: scaleFactor
    !> If true, do initial balancing using level_weights
    logical, intent(in) :: initial_balance
    ! --------------------------------------------------------------------------
    integer :: iLevel, minLevel, maxLevel
    real(kind=rk) :: level_weights(globalMaxLevels)
    ! --------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    !  Load mesh                                                             !
    ! -------------------------------------------------------------------------
    ! First check, if we are starting from a restart
    ! KJ: Now the restart is read from initial conditions, which are loaded
    ! much later than the config in the process flow
    call tem_load_restart( me       = restart,            &
      &                    conf     = solverHead%conf(1), &
      &                    tree     = me%tree,            &
      &                    timing   = simControl%now,     &
      &                    globProc = proc                )

    if( restart%controller%readRestart ) then
      ! If there is a restart, the timings in the params type have to be
      ! updated to those read from the restart
      call tem_timeControl_start_at_sim( me  = simControl%timeControl, &
        &                                now = simControl%now          )
    else

      write(logUnit(1),*) 'No read restart given. Loading mesh file'

      ! First load global info
      call tem_global_mesh_read( me     = me%tree%global,     &
        &                        conf   = solverHead%conf(1), &
        &                        myPart = proc%rank,          &
        &                        nParts = proc%comm_size,     &
        &                        comm   = proc%comm           )

      minLevel = me%tree%global%minLevel
      maxLevel = me%tree%global%maxLevel

      if( (minLevel /= maxLevel) .and. initial_balance ) then

        call tem_horizontalSpacer(fUnit = logUnit(1))
        write(logUnit(1),"(A)") 'Do initial level-wise load balancing.'

        ! Calculate level wise weights to scaleFactor^( level - minlevel )
        !                 acoustic      diffusive
        !  ------------------------------------------------
        !  fac =         |   2              4
        !  ------------------------------------------------
        !  minLevel      |   1              1
        !  minLevel + 1  |   2              4
        !  minLevel + 2  |   4              16
        !  ...
        !  ------------------------------------------------
        do iLevel = minLevel, maxLevel
          level_weights( iLevel ) = dble( scaleFactor ** (iLevel-minLevel) )
        end do

        call load_tem( me          = me%tree,              &
          &            conf        = solverHead%conf(1),   &
          &            myPart      = proc%rank,            &
          &            nParts      = proc%comm_size,       &
          &            comm        = proc%comm,            &
          &            levelWeight = level_weights,        &
          &            meshDir     = solverHead%meshFolder )
      else
        ! load the tree from the mesh = '' definition in case no restartRead is
        ! given in this case the mesh from the restart header is read
        ! (in tem_load_restart)
        call load_tem( me      = me%tree,              &
          &            conf    = solverHead%conf(1),   &
          &            myPart  = proc%rank,            &
          &            nParts  = proc%comm_size,       &
          &            comm    = proc%comm,            &
          &            meshDir = solverHead%meshFolder )
      end if ! minLevel /= maxLevel .and. initial_balance
    end if ! readRestart
    ! ---------------------------------------------------------------------------
    ! Done loading mesh.                                                       !
    ! ---------------------------------------------------------------------------

    ! Load boundary and qval
    call mus_load_bc_data( geometry = me,        &
      &                    rank     = proc%rank, &
      &                    comm     = proc%comm  )

    ! Load solidification/fluidification settings
    call mus_geomIncrHead_load( me          = me%geomIncr,        &
      &                         conf        = solverHead%conf(1), &
      &                         dynamicGeom = me%dynamicGeom      )

    ! load IBM data
    call mus_load_IBM( me     = me%globIBM,         &
      &                conf   = solverHead%conf(1), &
      &                rank   = proc%rank           )

  end subroutine mus_load_geom
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine invokes the treelm routines to load the boundary conditions
  !!
  subroutine mus_load_bc_data( geometry, rank, comm )
    ! --------------------------------------------------------------------------
    type( mus_geom_type ),  intent(inout)  :: geometry !< Treelmesh data
    integer, intent(in) :: rank, comm
    ! --------------------------------------------------------------------------
    integer :: iProp
    ! --------------------------------------------------------------------------

    ! ----------------------Load the boundary conditions------------------------
    ! Boundary conditions are loaded even in the case of dynamic load balancing
    ! when the restart is done at requisite load balancing intervals
    call init_tem_bc_prop( geometry%tree, rank, &
      &                    comm, geometry%boundary )

    do iProp = 1, geometry%tree%global%nProperties
      if( geometry%tree%global%property( iProp )%nElems > 0 ) then
        select case( geometry%tree%global%property( iProp )%bitpos )
        case( prp_hasBnd )
          ! Already loaded
        case( prp_hasQVal )
          ! Load qVal from disk
          ! prp_hasQVal is the 2nd property in mesh
          write(logUnit(1),"(A)") 'Loading qVal data from directory: '&
            &                     //trim(geometry%tree%global%dirname)
          call load_tem_BC_qVal(                                               &
            &           me       = geometry%boundary,                          &
            &           offset   = geometry%tree%Property(iProp)%Offset,       &
            &           nElems   = geometry%tree%Property(iProp)%nElems,       &
            &           basename = trim(geometry%tree%global%dirname)//'qval', &
            &           mypart   = rank,                   &
            &           comm     = comm                    )
          write(logUnit(1),*)'Done, reading the qVal!'

        end select ! property( iProp )%bitpos
      endif ! property( iProp )%nElems > 0
    enddo ! iProp

    call mus_build_posInProp( geometry )

    ! when qVal exist, it is allocated inside load_tem_BC_qVal
    ! otherwise allocate it with size 0
    if ( .not. allocated( geometry%boundary%qVal ) ) then
      allocate( geometry%boundary%qVAl(0,0) )
    end if

  end subroutine mus_load_bc_data
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine builds mapping between elements in tree to to propery list
  subroutine mus_build_posInProp( me )
    ! --------------------------------------------------------------------------
    type( mus_geom_type ), intent(inout) :: me !< Treelmesh data
    ! --------------------------------------------------------------------------
    integer :: iProp
    ! --------------------------------------------------------------------------

    do iProp = 1, me%tree%global%nProperties
      if( me%tree%global%property( iProp )%nElems > 0 ) then
        select case( me%tree%global%property( iProp )%bitpos )
        case( prp_hasBnd )
          if (  allocated( me%posInBndID )) deallocate( me%posInBndID )
          allocate( me%posInBndID(me%tree%nElems) )
          ! build me%posInBnd which maps tree to boundary%boundary_ID
          call tem_build_treeToProp_pointer( &
            &        treeToProp = me%posInBndID, &
            &        nElems     = me%tree%nElems,&
            &        ElemPropertyBits = me%tree%ElemPropertyBits,&
            &        prp_bit    = prp_hasBnd           )

          call mus_build_minBcID( minBcID    = me%minBcID,   &
            &                     bc_prop    = me%boundary,  &
            &                     posInBndID = me%posInBndID )

        case( prp_hasQVal )
          if ( allocated( me%posInQval )) deallocate( me%posInQval )
          allocate( me%posInQval(me%tree%nElems) )
          ! build me%posInQVal which maps tree to boundary%qVal
          call tem_build_treeToProp_pointer( &
            &        treeToProp = me%posInQVal, &
            &        nElems     = me%tree%nElems,&
            &        ElemPropertyBits = me%tree%ElemPropertyBits,&
            &        prp_bit    = prp_hasQVal         )

        end select ! property( iProp )%bitpos
      endif ! property( iProp )%nElems > 0
    enddo ! iProp

  end subroutine mus_build_posInProp
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine build minBcID for boundary elements, it is required if a
  !! element has more than one boundary in its directions.
  !! if a element has more than one boundary then use minBcID which depends
  !! on boundary order in seeder configuration
  subroutine mus_build_minBcID(minBcID, bc_prop, posInBndID)
    ! --------------------------------------------------------------------------
    integer, allocatable, intent(out) :: minBcID(:)
    !> boundary information from mesh
    type(tem_bc_prop_type), intent(in) :: bc_prop
    !> tree element position in boundaryID
    integer, intent(in) :: posInBndID(:)
    ! --------------------------------------------------------------------------
    integer :: iElem
    integer(kind=long_k) :: bcIDs(bc_prop%nSides)
    ! --------------------------------------------------------------------------
    allocate(minBcID(bc_prop%property%nElems))

    do iElem = 1, bc_prop%property%nElems
      bcIDs = bc_prop%boundary_ID(:, iElem)
      minBcID(iElem) = int(minval(bcIDs, bcIDs > 0_long_k))
    end do
  end subroutine mus_build_minBcID
  ! ************************************************************************** !

end module mus_geom_module
! ****************************************************************************** !
