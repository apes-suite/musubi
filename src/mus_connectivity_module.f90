! Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
! **************************************************************************** !
!> This module contains routines related to neighbor connectivity
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
module mus_connectivity_module
  use env_module,              only: rk, long_k
  use treelmesh_module,        only: treelmesh_type
  use tem_debug_module,        only: main_debug, dbgUnit
  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit, tem_toStr
  use tem_varSys_module,       only: tem_varSys_type
  use tem_varMap_module,       only: tem_varMap_type
  use tem_stencil_module,      only: tem_stencilHeader_type,  &
    &                                tem_stencil_findIndexOfDir
  use tem_geometry_module,     only: tem_determine_discreteVector,          &
    &                                tem_eligibleChildren, tem_CoordOfReal, &
    &                                tem_PosofId, tem_BaryOfId
  use tem_topology_module,     only: tem_directChildren, tem_levelOf, &
    &                                tem_parentOf, tem_IdOfCoord,     &
    &                                tem_FirstIdAtLevel
  use tem_construction_module, only: tem_levelDesc_type
  use tem_property_module,     only: prp_solid
  use tem_float_module,        only: operator(.fne.)
  use tem_element_module,      only: eT_fluid, eT_halo, eT_ghostFromFiner, &
    &                                eT_ghostFromCoarser

  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_bc_header_module,     only: glob_boundary_type

  implicit none
  private

  public :: mus_construct_connectivity
  public :: mus_updateConnectivity_forSymmetricBC
  public :: mus_intp_getSrcElemPosInLevelDesc
  public :: mus_intp_getSrcElemPosInTree


contains


! **************************************************************************** !
  !> Construct the propagation list for each element of 1st field.
  !!
  !! the bounce back rules have to be applied here
  !!
  !! Create connectivity array for one field such that it points to poisition
  !! of 1st Field in state array and use macro NGOFFSET to access other field
  !! states.
  !!
  subroutine mus_construct_connectivity(neigh, nSize, nElems, levelDesc, &
    &                                   stencil, varSys, stateVarMap     )
    ! --------------------------------------------------------------------------
    !> connectivity array
    integer, intent(out) :: neigh(:)
    !> number of elements in state array
    integer, intent(in) :: nSize
    !> number of elements in local partition
    integer, intent(in) :: nElems
    !> current level description
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> fluid stencil
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> global variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> state varMap
    type(tem_varMap_type), intent(in) :: stateVarMap
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, zeroPos, nghDir
    ! result: from which direction to take the pdf from (inverse for bounceback)
    integer :: GetFromDir
    ! result: from which element to take the pdf from (local for bounceback)
    integer :: GetFromPos
    integer :: neighPos      ! position of current neighElem in treeID list
    integer :: sourceDir
    integer :: stateVarPos(stencil%QQ)
    integer(kind=long_k) :: neighProp, elemProp
    logical :: missing_neigh_for_nonghost
    logical :: solidified
    ! --------------------------------------------------------------------------

    write(logUnit(1),*) 'Building state array connectivity ...'
!KM!    write(dbgUnit(1),*) 'Building state array connectivity ...'

    ! neigh array points to actual position of 1st field in the state
    ! array with size of nFields*QQ*nElems and
    ! position of current field in the state array can be accessed
    ! by using offset (iField-1)*QQ

    ! state varpos for 1st field since neigh is created only for 1st field
    stateVarPos(:stencil%QQ) = varSys%method%val(stateVarMap%varPos%val(1))  &
      &                                             %state_varPos(:stencil%QQ)

    ! set zero position direction to itself
    ! (e.g. the d3q6 stencil does not have a rest direction at all)
    if ( stencil%restPosition > 0 ) then
      zeroPos = stencil%restPosition
      do iElem = 1, nElems
        neigh( ( zeropos-1)* nsize + ielem) &
          & = ( ielem-1)* varsys%nscalars+ zeropos
      end do ! iElem
    end if

    ! set connectivity for non-rest positions
    elemLoop: do iElem = 1, nElems
      elemProp = levelDesc%property( iElem )
!KM!      write(dbgUnit(1),*) 'iElem: ', iElem, 'treeID: ', &
!KM!        & levelDesc%total(iElem), 'prop ', elemProp

      ! Loop over every link
      neighLoop: do iDir  = 1, stencil%QQN
        ! For push we need to find different neighbors than for pull
        ! -> ngDir
        nghDir =  stencil%cxdirinv(idir)
        neighPos = levelDesc%neigh(1)%nghElems(nghDir, iElem)

        ! JQ: get neighbor's property
        if (neighPos > 0) then
          neighProp = levelDesc%property(neighPos)
        else
          neighProp = 0_long_k
        end if
!KM!        write(dbgUnit(1),*) iDir, nghDir, 'neighID', &
!KM!          & levelDesc%total(neighPos), 'neighProp ', neighProp


        missing_neigh_for_nonghost                               &
          &  = (neighPos <= 0)                                   &
          &  .and. ( (iElem <= levelDesc%offset(2, eT_fluid) )   &
          &          .or. (iElem > levelDesc%offset(1, eT_halo)) )

        solidified = ( btest(neighProp, prp_solid)     &
          &            .or. btest(elemProp, prp_solid) )

        sourceDir = iDir
        if (missing_neigh_for_nonghost .or. solidified) then
          ! Treat missing or solid neighbor element for non-ghosts as solid
          ! with Bounce-Back rule (invDir).
          ! Get inverse direction of LOCAL element:
          sourceDir = stencil%cxDirInv( iDir )
        end if

        ! link to get from neighbor
        GetFromDir = stateVarPos(sourceDir)

        GetFromPos = neighPos
        if ((neighPos <= 0) .or. solidified) then
          GetFromPos = iElem
        end if

        neigh( (idir-1)* nsize+ ielem )                       &
          &  = ( getfrompos-1)* varsys%nscalars+getfromdir

      end do neighLoop
    end do elemLoop
    write(logUnit(1),*) 'Done with state array connectivity.'
  end subroutine mus_construct_connectivity
! **************************************************************************** !


! **************************************************************************** !
  !> Update the connectivity for elements with symmetric boundary condition
  !! such that during the propagation they are applied implicitly.
  !!
  !! update connectivity only for neighbors along edge and corner because
  !! face neighbor boundary are treated as bounce back
  subroutine mus_updateConnectivity_forSymmetricBC(neigh, nSize, iLevel,  &
    & levelDesc, layout, varSys, stateVarMap, nBCs, globBC, nSymBCs,      &
    & symmetricBCs)
    ! --------------------------------------------------------------------------
    !> connectivity array
    integer, intent(inout) :: neigh(:)
    !> number of elements in state array
    integer, intent(in) :: nSize
    !> current level
    integer, intent(in) :: iLevel
    !> current level description
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> scheme layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> global variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> state varMap
    type(tem_varMap_type), intent(in) :: stateVarMap
    integer, intent(in) :: nBCs !< number of boundary conditiosn
    !> global boundary information
    type(glob_boundary_type), intent(in) :: globBC(nBCs)
    !> number of symmetric boundaries
    integer, intent(in) :: nSymBCs
    !> symmetric boundary ids
    integer, intent(in) :: symmetricBCs(nSymBCs)
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, nScalars
    ! result: from which direction to take the pdf from (inverse for bounceback)
    integer :: GetFromDir
    integer :: neighPos      ! position of current neighElem in treeID list
    integer :: stateVarPos(layout%fStencil%QQ)
    integer :: nElems
    integer :: iBC, elemPos, neighDirInd, neighDirInd_inv, symDirInd
    integer :: iSymBC
    integer :: normal(3), symDir(3), neighDir(3)
    ! --------------------------------------------------------------------------
    write(logUnit(1),*) 'Updating state array connectivity for symmetric BC...'
!KM!    write(dbgUnit(1),*) 'Updating state array connectivity for symmetric BC...'
!KM!    write(dbgUnit(1),*) 'Number of symmetric boundaries to update: ', nSymBCs
!KM!    flush(dbgUnit(1))

    write(logUnit(1),*) 'Number of symmetric boundaries to update: ', nSymBCs

    ! neigh array points to actual position of 1st field in the state
    ! array with size of nFields*QQ*nElems and
    ! position of current field in the state array can be accessed
    ! by using offset (iField-1)*QQ

    QQ       = layout%fStencil%QQ
    nScalars = varSys%nScalars
    ! state varpos for 1st field since neigh is created only for 1st field
    stateVarPos(:QQ) = varSys%method%val(stateVarMap%varPos%val(1))  &
      &                                             %state_varPos(:QQ)

    do iSymBC = 1, nSymBCs
      iBC = symmetricBCs(iSymBC)
!KM!      write(dbgUnit(1),*) 'Symmetric BC: ', trim(scheme%globBC(iBC)%label)
      nElems = globBC(iBC)%nElems(iLevel)
!KM!        write(dbgUnit(1),*) 'iLevel ', iLevel, ' nElems:', nElems
      elemLoop: do iElem = 1, nElems
        elemPos = globBC(iBC)%elemLvl(iLevel)%elem%val(iElem)
!KM!        write(dbgUnit(1),*) 'iElem: ', iElem, 'treeID: ', &
!KM!          & scheme%levelDesc(iLevel)%total(elemPos)
          ! Loop over every link
        neighLoop: do iDir  = 1, layout%fStencil%QQN
          if ( globBC(iBC)%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
!KM!            write(dbgUnit(1), *) 'iDir ', iDir
            ! Do bounce back for axis-aligned links and links along the normal
            ! direction of boundary

            ! this sum is 1 for axis-aligned and they are already treated
            ! with bounce back rule in construct_connectivity routine.
            ! so goto next dir
            if (sum( abs(layout%fStencil%cxDir( :, iDir ))) == 1) cycle

            ! also do nothing if a link is same as normal direction of
            ! boundary
            if (globBC(iBC)%elemLvl(iLevel)%normalInd%val(iElem) == iDir) cycle

            ! for rest of the directions compute the symDir of iDir
            ! with respect to the boundary normal and also find from
            ! which neighbor the symdir must be taken.
            ! symDir are outgoing directions pointing towards the boundary
            normal = globBC(iBC)%elemLvl(iLevel)%normal%val(:,iElem)
            symDir = layout%fStencil%cxDir(:, iDir) - 2*normal
            call tem_determine_discreteVector(symDir, layout%prevailDir)
            symDirInd = tem_stencil_findIndexOfDir( symDir,               &
              &                                     layout%fStencil%cxDir )

            ! neighbor direction
            neighDir = normal - layout%fStencil%cxDir(:, iDir)
            call tem_determine_discreteVector( neighDir,         &
              &                                layout%prevailDir )

            neighDirInd = tem_stencil_findIndexOfDir( neighDir,             &
              &                                       layout%fStencil%cxDir )

            ! neighbor to get symDir
            ! NgDir uses inverse direction so use inverse of neighDirInd
            neighDirInd_inv = layout%fStencil%cxDirInv(neighDirInd)
            neighPos = levelDesc%neigh(1)%nghElems(                  &
              &  layout%fstencil %cxdirinv( neighdirind_inv), elemPos )

!KM!write(dbgUnit(1),*) 'cxDir ', scheme%layout%fStencil%cxDir(:, iDir)
!KM!write(dbgUnit(1),*) 'symDir: ', symDir
!KM!write(dbgUnit(1),*) 'neighDir', neighDir, ' neighPos: ', neighPos
              ! update connectivity only if the neighbor element exist
            if (neighPos>0) then
!KM!write(dbgUnit(1),*) 'neighID: ', scheme%levelDesc(iLevel)%total(neighPos)
              getFromDir = stateVarPos(symDirInd)
              neigh( ( idir-1)* nsize + elempos)  &
                & = ( neighpos-1)* nscalars+ getfromdir
            end if

          end if !bitmask
        end do neighLoop
!KM!        write(dbgUnit(1),*)
      end do elemLoop
!KM!      write(dbgUnit(1),*)
    end do ! iBC

    write(logUnit(1),*) 'Done with state array connectivity for symmetric BC.'
!KM!    write(dbgUnit(1),*) 'Done with state array connectivity for symmetric BC.'
!KM!    flush(dbgUnit(1))
!KM!    stop

  end subroutine mus_updateConnectivity_forSymmetricBC
! **************************************************************************** !

! **************************************************************************** !
  ! This routine provides the source element position in level descriptor
  ! for interpolation. Here, Ghost elements are considered as ghost.
  ! This routine is used in setup_indices for state and derive variable where
  ! point value is computed using elements in the same level.
  subroutine mus_intp_getSrcElemPosInLevelDesc(srcElemPos, weights, nSrcElems, &
    & point, statePos, neigh, baryOfTotal, nElems, nSolve, stencil, nScalars,  &
    & excludeHalo)
    ! --------------------------------------------------------------------------
    !> position of source element in the levelwise state array
    integer, intent(out) :: srcElemPos(:)
    !> weights for interpolation
    real(kind=rk), intent(out) :: weights(:)
    !> number of source elements found
    integer, intent(out) :: nSrcElems
    !> target point
    real(kind=rk), intent(in) :: point(3)
    !> position of element which contains target point on the levelwise state
    !! array
    integer, intent(in) :: statePos
    !> neighbor connectivity on the level in which the element of the point
    !! is found
    integer, intent(in) :: neigh(:)
    !> bary of elements in the total list
    real(kind=rk), intent(in) :: baryOfTotal(:,:)
    !> nElems in the connectivity array
    integer, intent(in) :: nElems
    !> number of elements excluding halos
    integer, intent(in) :: nSolve
    !> stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> number of scalars in the varSys
    integer, intent(in) :: nScalars
    !> excludeHalo element for interpolation. True of auxFieldVar.
    !! True for state var if reducedComm is False else otherwise.
    !! For reducedComm only part of state array is communicated so exclude
    !! halo elements for state vars
    logical, intent(in) :: excludeHalo
    ! --------------------------------------------------------------------------
    integer :: iNeigh, neighPos
    integer :: iSrc
    real(kind=rk) :: bary(3), dist(3), dist_len
    ! --------------------------------------------------------------------------
    nSrcElems = 0
    srcElemPos = 0
    ! use statePos as 1st source element for interpolation
    ! if restPosition is part of Stencil
    if (stencil%QQ /= stencil%QQN) then
      nSrcElems = 1
      ! use position in levelwise for interpolation
      srcElemPos(nSrcElems)  = statePos
    end if

    bary = baryOfTotal(statePos, :)
    dist = abs(bary - point(:))
    dist_len = sqrt(dot_product(dist, dist))
    ! if point is exact bary center of current element then
    ! no need to do interpolation
    if ( dist_len .fne. 0.0_rk ) then

      ! get existing neighbor elements in actual tree
      do iNeigh = 1, stencil%QQN

        ! Find neighor using neigh array because for boundary, neigh array
        ! points to current element so no need to check for existence of the
        ! neighbor element
        neighPos =                                                         &
          & int((neigh(( stencil%cxdirinv( ineigh)-1)* nelems+ statepos)-1)/ nscalars)+1

        ! skip this neighbor and goto next neighbor
        if (excludeHalo .and. neighPos > nSolve) cycle

        ! append to list of source element only if neighPos is not statePos
        ! if neighPos is statePos then neighbor is boundary
        if (neighPos /= statePos) then
          nSrcElems = nSrcElems + 1
          srcElemPos(nSrcElems)  = neighPos
        end if
      end do !iNeigh

      ! calculate weight
      do iSrc = 1, nSrcElems
        bary = baryOfTotal(srcElemPos(iSrc),:)
        dist = abs(bary - point(:))
        weights(iSrc) = 1.0_rk/sqrt(dot_product( dist, dist ))
      end do
    else
      ! if point is exact bary center of current element then
      ! no need to do interpolation
      weights(nSrcElems) = 1.0_rk
    end if

    ! normalize weights
    weights = weights/sum(weights)

  end subroutine  mus_intp_getSrcElemPosInLevelDesc
! **************************************************************************** !


! **************************************************************************** !
  ! This routine provides source element position in the global tree. It is
  ! used in state and derive variable get_point routines to interpolate to a
  ! point using elements in global fluid tree.
  ! NOTE: Ghost and halo elements are not part of source elements.
  subroutine mus_intp_getSrcElemPosinTree(srcElemPos, weights, nSrcElems, &
    & point, stencil, tree, levelPointer, levelDesc)
    ! --------------------------------------------------------------------------
    !> position of source element in the level independent list
    integer, intent(out) :: srcElemPos(:)
    !> weights for interpolation
    real(kind=rk), intent(out) :: weights(:)
    !> number of source elements found
    integer, intent(out) :: nSrcElems
    !> target point
    real(kind=rk), intent(in) :: point(3)
    !> stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree
    !> Pointer from treeIDlist entry to level-wise fluid part of total list
    integer, intent(in)         :: levelPointer(:)
    !> level description of all levels
    type(tem_levelDesc_type), intent(in) :: levelDesc(tree%global%minLevel:)
    ! --------------------------------------------------------------------------
    integer :: elemPos, neighPos, statePos, childPos
    integer :: iNeigh, iSrc, iLevel, iChild
    real(kind=rk) :: bary(3), dist(3), dist_len
    integer :: coord(4)
    integer(kind=long_k) :: neighID, childID, parentID
    integer, allocatable :: eligible_childs(:)
    logical :: skip_neigh
    real(kind=rk) :: mindist
    integer :: closest_neigh
    ! --------------------------------------------------------------------------
    nSrcElems = 0
    srcElemPos = 0
    weights = 0.0_rk
    mindist = huge(dist)
    closest_neigh = 0

    ! Coordinate of requested point in maxlevel
    coord =  tem_CoordOfReal(tree, point(:), tree%global%maxLevel )
    ! returns position of existing element in tree which contains point
    ! 0 if no corresponding node is found,
    ! or the negative of the found ID, if it is a virtual node.
    elemPos = abs(tem_PosofId( tem_IdOfCoord(coord), tree%treeID ))

    ! Extrapolate if point is outside fluid domain.
    ! Use closest neighbor as elemPos and uses its neighbors for interpolation
    ! Using algorithm from tem_cano_checkNeigh
    if (elemPos == 0) then
      ! Assuming the point is near boundary at finest level
      iLevel = tree%global%maxLevel
      directionLoop: do iNeigh = 1, stencil%QQN
        neighID = tem_IdOfCoord(                              &
          &         [ coord(1) + stencil%cxDir( 1, iNeigh ),  &
          &           coord(2) + stencil%cxDir( 2, iNeigh ),  &
          &           coord(3) + stencil%cxDir( 3, iNeigh ),  &
          &           coord(4) ], tem_FirstIdAtLevel(iLevel) )
        neighPos = abs( tem_PosOfId( neighID, tree%treeID) )
        if (neighPos > 0) then
          bary = tem_BaryOfId(tree, neighID)
          dist = bary - point
          dist_len = sqrt(dot_product(dist, dist))
          if (dist_len < mindist) then
            mindist = dist_len
            closest_neigh = neighPos
          else if (dist_len <= mindist) then
            ! Make sure to pick the first element in the SFC ordering
            ! as the nearest neighbor, if there are multiple with the
            ! same distance (<= comparison to avoid == comparison for
            ! reals)
            closest_neigh = min(neighPos, closest_neigh)
          end if
        end if
      end do directionLoop
      elemPos = closest_neigh
    end if

    ! level of existing element
    iLevel = tem_levelOf( tree%treeID( elemPos ) )
    ! position of element in levelDesc total list
    statePos = levelPointer( elemPos )

    ! use elemPos as 1st source element for interpolation
    ! if restPosition is part of Stencil
    if (stencil%QQ /= stencil%QQN) then
      nSrcElems = 1
      ! use position in levelwise for interpolation
      srcElemPos(nSrcElems)  = elemPos
    end if

    bary = tem_BaryOfID(tree, tree%treeID(elemPos) )
    dist = abs(bary - point(:))
    dist_len = sqrt(dot_product(dist, dist))
    ! if point is exact bary center of current element then
    ! no need to do interpolation
    if ( dist_len .fne. 0.0_rk ) then

      ! get existing neighbor elements in actual tree
      do iNeigh = 1, stencil%QQN

        ! neighbor position in levelwise list
        neighPos = levelDesc(iLevel)%neigh(1)%nghElems(iNeigh, statePos)
        ! skip this neighbor if it is halo or boundary and goto next neighbor
        skip_neigh                                                 &
          & = (neighPos <= 0)                                      &
          & .or. ( neighPos > levelDesc(iLevel)%offset(1, eT_halo) )

        if ( skip_neigh ) cycle

        neighID = levelDesc(iLevel)%total(neighPos)
        ! KM: DO NOT CHANGE THIS ORDER OF IF because it relys on order in total
        ! list i.e. ghostFromFiner are added after ghostFromCoarser.
        if (neighPos > levelDesc(iLevel)%offset(1, eT_ghostFromFiner) .and. &
          & neighPos < levelDesc(iLevel)%offset(2, eT_ghostFromFiner) ) then
          ! neighbor is ghostFromFiner, use eligible children in direction
          ! iNeigh as source elements
          call tem_eligibleChildren( eligible_childs,    &
            &                        stencil%map(iNeigh) )

          do iChild = 1, size(eligible_childs)
            childID = neighID*8_long_k + eligible_childs(iChild)
            childPos = tem_PosOfId(childID, tree%treeID)
            ! children could be distributed so addd to source list only if
            ! it is local
            if (childPos > 0) then
              nSrcElems = nSrcElems + 1
              srcElemPos(nSrcElems)  = childPos
            end if
          end do
        else if ( neighPos > levelDesc(iLevel)                          &
          &                    %offset(1, eT_ghostFromCoarser) .and.    &
          & neighPos < levelDesc(iLevel)%offset(2, eT_ghostFromCoarser) ) then
          ! neighbor is ghostFromCoarser, use parent as source element
          nSrcElems = nSrcElems + 1
          parentID = tem_parentOf( neighID )
          srcElemPos(nSrcElems) = tem_PosOfId(parentID, tree%treeID)
        else if (neighPos /= statePos) then
          ! append to list of source element only if neighPos is not elemPos
          ! if neighPos is elemPos then neighbor is boundary
          nSrcElems = nSrcElems + 1
          srcElemPos(nSrcElems)  = levelDesc(iLevel)%pntTID(neighPos)
        end if
      end do !iNeigh
    end if ! dist_len /= 0

    ! calculate weight
    if (nSrcElems == 1) then
      ! if point is exact bary center of current element then
      ! no need to do interpolation
      weights(nSrcElems) = 1.0_rk
    else
      do iSrc = 1, nSrcElems
        bary = tem_BaryOfID(tree, tree%treeID( srcElemPos(iSrc) ))
        dist = abs(bary - point(:))
        weights(iSrc) = 1.0_rk/sqrt(dot_product( dist, dist ))
      end do
      ! normalize weightss
      weights = weights/sum(weights)
    end if

  end subroutine  mus_intp_getSrcElemPosInTree
! **************************************************************************** !

end module mus_connectivity_module
