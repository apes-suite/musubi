! Copyright (c) 2017, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
! **************************************************************************** !
!> author: Kannan Masilamani, Verena Krupp
!! This module contains routine to retrieve state variables for getElement,
!! getPoint, setupIndices and getValOfIndex
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
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015, 2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
!--------------------------------------------
!    A O S - Array of structures layout new
!-------------------------------------------
! Access to get_point value output
! Access to get_element value output
module mus_stateVar_module

  use, intrinsic :: iso_c_binding,  only: c_ptr, c_f_pointer, c_loc
  use env_module,                   only: rk, long_k

  use tem_logging_module,           only: logUnit
  use tem_debug_module,             only: dbgUnit
  use tem_aux_module,               only: tem_abort
  use tem_varSys_module,            only: tem_varSys_type,                    &
    &                                     tem_varSys_op_type
  use tem_grow_array_module,        only: grw_realArray_type, &
    &                                     grw_intArray_type, &
    &                                     append, init, truncate
  use tem_pointData_module,         only: tem_pointData_list_type, &
    &                                     init, append, truncate
  use treelmesh_module,             only: treelmesh_type
  use tem_geometry_module,          only: tem_CoordOfReal, &
    &                                     tem_PosofId, tem_BaryOfId
  use tem_time_module,              only: tem_time_type
  use tem_topology_module,          only: tem_IdOfCoord
  use tem_topology_module,          only: tem_levelOf
  use tem_element_module,           only: eT_minRelevant, eT_maxRelevant, &
    &                                     eT_fluid

  ! include musubi modules
  use mus_scheme_type_module,  only: mus_scheme_type
  use mus_varSys_module,       only: mus_varSys_data_type,            &
    &                                mus_createSrcElemInTreeForGetPoint
  use mus_connectivity_module, only: mus_intp_getSrcElemPosInLevelDesc, &
    &                                mus_intp_getSrcElemPosInTree


  implicit none

  private

  public :: mus_access_state_forElement
  public :: mus_access_stateFetch_forElement
  public :: mus_access_stateFetch_now_forElement
  public :: mus_stateVar_forPoint
  public :: mus_accessVar_setupIndices
  public :: mus_stateVar_fromIndex
  public :: mus_stateVar_Fetch_fromIndex
  public :: mus_stateVar_Fetch_now_fromIndex

contains

  ! ************************************************************************* !
  !> Return the solver state variable for a given set of elements
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine mus_access_state_forElement(fun, varsys, elempos, time, &
    &                                              tree, nElems, nDofs, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: statePos, iElem, iComp, iLevel
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: nSize, iFld, QQ, nScalars
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    ! myPos can be only used to access state variable since state variables are
    ! added to varSys first.
    iFld = fun%myPos
    QQ = fPtr%solverData%scheme%layout%fStencil%QQ
    nScalars = varSys%nSCalars

    ! res is always AOS layout
    res = 0.0_rk
    associate( state => fPtr%solverData%scheme%state,                &
      &        pdf => fPtr%solverData%scheme%pdf,                    &
      &        levelPointer => fPtr%solverData%geometry%levelPointer )
      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
        nSize = pdf( iLevel )%nSize
        do iComp = 1, fun%nComponents
          res( (iElem-1)*fun%nComponents+iComp ) =                       &
            & state( iLevel )%val(                                       &
            ! position of this state variable in the state array
& ( statepos-1)* nscalars+icomp+( ifld-1)* qq, &
            & pdf( iLevel )%nNext )
        end do !iComp
      end do !iElem
    end associate

  end subroutine mus_access_state_forElement
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> Return the solver state variable for a given set of elements by using
  !! FETCH macro for nNow
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine mus_access_stateFetch_now_forElement( &
    & fun, varsys, elempos, time, tree, nElems, nDofs, res   )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: statePos, iElem, iComp, iLevel
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: nSize, iFld, QQ, nScalars
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    QQ = fPtr%solverData%scheme%layout%fStencil%QQ
    nScalars = varSys%nSCalars
    iFld = fun%input_varPos(1)

    ! res is always AOS layout
    res = 0.0_rk
    associate( state => fPtr%solverData%scheme%state,                &
      &        pdf => fPtr%solverData%scheme%pdf,                    &
      &        levelPointer => fPtr%solverData%geometry%levelPointer )
      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
        nSize = pdf( iLevel )%nSize
        do iComp = 1, fun%nComponents
          res( (iElem-1)*fun%nComponents+iComp ) =                        &
            & state( iLevel )%val(                                        &
            ! position of this state variable in the state array
&  pdf(ilevel)%neigh((icomp-1)* nsize+ statepos)+( ifld-1)* qq+ nscalars*0, &
            & pdf( iLevel )%nNow )
        end do !iComp
      end do !iElem
    end associate

  end subroutine mus_access_stateFetch_now_forElement
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> Return the solver state variable for a given set of elements by using
  !! FETCH macro for nNext
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine mus_access_stateFetch_forElement(fun, varsys, elempos, &
    &                                                   time, tree, nElems,   &
    &                                                   nDofs, res            )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: statePos, iElem, iComp, iLevel
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: nSize, iFld, QQ, nScalars
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    QQ = fPtr%solverData%scheme%layout%fStencil%QQ
    nScalars = varSys%nSCalars
    iFld = fun%input_varPos(1)

    ! res is always AOS layout
    res = 0.0_rk
    associate( state => fPtr%solverData%scheme%state,                &
      &        pdf => fPtr%solverData%scheme%pdf,                    &
      &        levelPointer => fPtr%solverData%geometry%levelPointer )
      do iElem = 1, nElems
        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
        nSize = pdf( iLevel )%nSize
        do iComp = 1, fun%nComponents
          res( (iElem-1)*fun%nComponents+iComp ) =                        &
            & state( iLevel )%val(                                        &
            ! position of this state variable in the state array
&  pdf(ilevel)%neigh((icomp-1)* nsize+ statepos)+( ifld-1)* qq+ nscalars*0, &
            & pdf( iLevel )%nNext )
        end do !iComp
      end do !iElem
    end associate

  end subroutine mus_access_stateFetch_forElement
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> State variable for a given set of points using linear interpolation.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_point]].
  !!
  recursive subroutine mus_stateVar_forPoint(fun, varsys, point, time, tree, &
    &                                        nPnts, res                      )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: elemPos, statePos, iPnt, iLevel
    real(kind=rk), allocatable :: srcRes(:), pntVal(:), weights(:)
    integer :: iSrc, iComp, nSize, iFld, QQ, nScalars, nSrcElems
    integer, allocatable :: srcElemPos(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    ! myPos can be only used to access state variable since state variables are
    ! added to varSys first.
    iFld = fun%myPos
    nScalars = varSys%nScalars
    QQ = scheme%layout%fStencil%QQ
    allocate(srcRes(QQ*fun%nComponents))
    allocate(pntVal(fun%nComponents))
    allocate(weights(QQ))
    allocate(srcElemPos(QQ))

    res = 0.0_rk
    do iPnt = 1, nPnts
      srcRes = 0.0_rk

      ! get position of source element position in global tree for
      ! interpolation.
      ! Also calculate weights for interpolation using distance between the
      ! point and source element barycenter
      call mus_intp_getSrcElemPosInTree(          &
        & srcElemPos   = srcElemPos,              &
        & weights      = weights,                 &
        & nSrcElems    = nSrcElems,               &
        & point        = point(iPnt,:),           &
        & stencil      = scheme%layout%fStencil,  &
        & tree         = tree,                    &
        & levelPointer = fPtr%solverData%geometry &
        &                    %levelPointer,       &
        & levelDesc    = scheme%levelDesc         )

      ! get source element values
      do iSrc = 1, nSrcElems
        ! position in global tree
        elemPos = srcElemPos(iSrc)
        ! position of element in levelDesc total list
        statePos = fPtr%solverData%geometry%levelPointer( elemPos )
        iLevel = tem_levelOf( tree%treeID( elemPos ) )
        nSize = scheme%pdf( iLevel )%nSize
        do iComp = 1, fun%nComponents
          srcRes( (iSrc-1)*fun%nComponents + iComp )                           &
            & = scheme%state( iLevel )%val(                                    &
& ( statepos-1)* nscalars+icomp+( ifld-1)* qq,&
            & scheme%pdf( iLevel )%nNext )
        end do !iComp
      end do !iSrc

      ! Linear interpolation res = sum(weight_i*phi_i)
      pntVal = 0.0_rk
      do iSrc = 1, nSrcElems
        pntVal(:) = pntVal(:) + weights(iSrc)                           &
          & * srcRes( (iSrc-1)*fun%nComponents+1 : iSrc*fun%nComponents )
      end do
      res( (iPnt-1)*fun%nComponents+1 : iPnt*fun%nComponents ) = pntVal

    end do !iPnt
  end subroutine mus_stateVar_forPoint
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> This routine takes points coordinates, stores them in the method_data and
  !! return indices where points are located in the growing array of points or
  !! values ( sometimes we do not need to store the points )
  !! It is need to setup points for every variable. Points will be provided by
  !! boundaries or sources depends on what uses the variable. This points do not
  !! change with time . This indices will be stored in corresponding boundary
  !! or source to evaluate a value on that point later using
  !! tem_varSys_proc_getValOfIndex.
  subroutine mus_accessVar_setupIndices( fun, varSys, point, offset_bit, &
    &                                   iLevel, tree, nPnts, idx        )
    ! -------------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, for this routine
    !! we need the location where to store the points.
    class(tem_varSys_op_type), intent(in) :: fun
    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys
    !>  arrays of points for which the indices are returned
    real(kind=rk), intent(in) :: point(:,:)
    !> Offset bit encoded as character for every point.
    !! If not present default is to center i.e offset_bit = achar(1+4+16)
    character, optional, intent(in) :: offset_bit(:)
    !> the point data need to be loaded levelwise, we need the current level
    integer, intent(in) :: iLevel
    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree
    !> number of points
    integer, intent(in) :: nPnts
    integer, intent(out) :: idx(:)
    ! -------------------------------------------------------------------------- !
    integer :: iPoint
    integer :: coord(4)
    integer :: elemPos, loc_level, statePos
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk), allocatable :: weights(:)
    integer, allocatable :: srcElemPos(:)
    integer :: nSrcElems
    integer(kind=long_k) :: treeID
    ! -------------------------------------------------------------------------- !

!    write(dbgUnit(4),*) 'setup indices for the points of variable ', &
!      &                 trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    allocate(srcElemPos(scheme%layout%fStencil%QQ))
    allocate(weights(scheme%layout%fStencil%QQ))

    idx = 0
    ! store elemPos of each point in global treeID list
    do iPoint = 1, nPnts
      weights = 0.0_rk

      ! Get the treeID of the element and then get the element position
      ! in the levelwise array
      coord =  tem_CoordOfReal( mesh  = tree,                &
        &                       point = point(iPoint,:),     &
        &                       level = tree%global%maxLevel )

      treeID = tem_IdOfCoord(coord)
      ! returns position of existing element in tree which contains point
      ! 0 if no corresponding node is found,
      ! or the negative of the found ID, if it is a virtual node.
      elemPos = abs(tem_PosofId( treeID, tree%treeID ))

!      write(dbgUnit(1),*) 'iPnt ', iPoint, ' elemPos ', elemPos, &
!        & 'treeID ', treeID, 'point ', point(iPoint,:)

      ! skip this point if no corresponding element is found
      if (elemPos == 0) cycle
      ! store point data only if element position is found in tree

      ! level of existing element
      loc_level = tem_levelOf( tree%treeID( elemPos ) )

      ! append grw array for this point's element position in tree
      call append(me  = fPtr%pointData%pntLvl(iLevel)%elemPos, &
        &         val = elemPos                                )
      ! Actual level on which this point was found
      call append(me  = fPtr%pointData%pntLvl(iLevel)%pntLevel, &
        &         val = loc_level                               )
      ! Append points to growing array of points
      call append(me  = fPtr%pointData%pntLvl(iLevel)%grwPnt, &
        &         val = point(iPoint,:)                       )

      ! position of element in levelDesc total list
      statePos = fPtr%solverData%geometry%levelPointer( elemPos )

      ! get position of source element position in state array for interpolation
      ! using neighbor connectivity array.
      ! Also calculate weights for interpolation using distance between the
      ! point and source element barycenter
      ! Excluded halo elements as they require communication of all links
      call mus_intp_getSrcElemPosInLevelDesc(                        &
        & srcElemPos  = srcElemPos,                                  &
        & weights     = weights,                                     &
        & nSrcElems   = nSrcElems,                                   &
        & point       = point(iPoint,:),                             &
        & statePos    = statePos,                                    &
        & neigh       = fPtr%solverData%scheme%pdf(loc_level)%neigh, &
        & baryOfTotal = scheme%levelDesc(loc_level)%baryOfTotal,     &
        & nElems      = scheme%pdf(loc_level)%nSize,                 &
        & nSolve      = scheme%pdf(loc_level)%nElems_computed,       &
        & stencil     = scheme%layout%fStencil,                      &
        & excludeHalo = .true.,                                      &
        & nScalars    = varSys%nScalars                              )

      ! First position of source elements
      call append(me  = fPtr%pointData%pntLvl(iLevel)%srcElem%first,  &
        &         val = fPtr%pointData%pntLvl(iLevel)%srcElem%elemPos &
        &               %nVals + 1                                    )

      ! Append all src elemPos
      call append(me  = fPtr%pointData%pntLvl(iLevel)%srcElem%elemPos, &
        &         val = srcElemPos(1:nSrcElems)                        )

      ! last position of source elements
      call append(me  = fPtr%pointData%pntLvl(iLevel)%srcElem%last,   &
        &         val = fPtr%pointData%pntLvl(iLevel)%srcElem%elemPos &
        &               %nVals                                        )

      ! weights for srcElements
      call append(me  = fPtr%pointData%pntLvl(iLevel)%srcElem%weight, &
        &         val = weights(1:nSrcElems)                          )

      ! store the index,
      idx(iPoint)= fPtr%pointData%pntLvl(iLevel)%elemPos%nVals
    end do
!    write(dbgUnit(1),*) 'varName ', trim(varSys%varName%val(fun%myPos))
!    write(dbgUnit(1),*) 'idx ', idx

    fPtr%pointData%pntLvl(iLevel)%nPnts = fPtr%pointData%pntLvl(iLevel)%grwPnt &
      &                                       %coordX%nVals
    call truncate(fPtr%pointData%pntLvl(iLevel)%grwPnt)
    call truncate(fPtr%pointData%pntLvl(iLevel)%elemPos)
    call truncate(fPtr%pointData%pntLvl(iLevel)%srcElem%first)
    call truncate(fPtr%pointData%pntLvl(iLevel)%srcElem%last)
    call truncate(fPtr%pointData%pntLvl(iLevel)%srcElem%elemPos)
    call truncate(fPtr%pointData%pntLvl(iLevel)%srcElem%weight)

  end subroutine mus_accessVar_setupIndices
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> Routine to get the actual value for a given array of indices.
  !! The indices belong to the grwarray of points storing levelwise in
  !! Pointdata%pntLvl(iLevel).
  !! Hence this routines takes the indeices as input, can refer to the pointData
  !! and evaluate the variable and returns the values
  subroutine mus_stateVar_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                idxLen,  nVals,  res            )
    ! -------------------------------------------------------------------------- !
    !> Description of the method to obtain the variables,
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: n
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: dependes on number of first index for contiguous array,
    !! but the sum of all idxLen is equal to nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: iComp, iSrc, iVal, elemPos, first, last, nSize, loc_level
    integer :: iSrcElems
    real(kind=rk), allocatable :: srcRes(:), pntVal(:)
    real(kind=rk) :: weight
    integer :: iFld, QQ, nScalars
    ! -------------------------------------------------------------------------- !

!    write(dbgUnit(4),*) 'get the values of indices for variable ',  &
!      &                  trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    ! myPos can be only used to access state variable since state variables are
    ! added to varSys first.
    iFld = fun%myPos
    nScalars = varSys%nScalars
    QQ = scheme%layout%fStencil%QQ
    allocate(srcRes(QQ*fun%nComponents))
    allocate(pntVal(fun%nComponents))

    res = 0.0_rk

    ! distinguish if we have an array of index or we have contingous memory
    ! access where index are always first entries!
    if (present(idxLen)) then
      call tem_abort('Error: idxLen is not supported in get_valOfIndex for ' &
        &          //'state variable')
    else
      ! Get the state value at the specific point
      do iVal = 1, nVals
        if (idx(iVal)>0) then

          ! elemPos in tree
          ! elemPos = fPtr%pointData%pntLvl(iLevel)%elemPos%val(idx(iVal))

          first = fPtr%pointData%pntLvl(iLevel)%srcElem%first%val(idx(iVal))
          last = fPtr%pointData%pntLvl(iLevel)%srcElem%last%val(idx(iVal))
          loc_level = fPtr%pointData%pntLvl(iLevel)%pntLevel%val(idx(iVal))
          nSize = scheme%pdf( loc_level )%nSize

          ! get pdf's of source elements
          srcRes = 0.0_rk
          iSrcElems = 0
          do iSrc = first, last
            iSrcElems = iSrcElems + 1

            ! position of element in levelDesc total list
            elemPos = fPtr%pointData%pntLvl(iLevel)%srcElem%elemPos%val(iSrc)

            do iComp = 1, fun%nComponents
              srcRes( (iSrcElems-1)*fun%nComponents + iComp )                  &
                & = scheme%state( loc_level )%val(                             &
& ( elempos-1)* nscalars+icomp+( ifld-1)* qq, &
                & scheme%pdf( loc_level )%nNext )
            end do !iComp
          end do !iSrc

          ! Linear interpolation res = sum(weight_i*phi_i)
          pntVal = 0.0_rk
          iSrcElems = 0
          do iSrc = first, last
            weight = fPtr%pointData%pntLvl(iLevel)%srcElem%weight%val(iSrc)
            iSrcElems = iSrcElems + 1
            pntVal(:) = pntVal(:) + weight &
              & * srcRes( (iSrcElems-1)*fun%nComponents+1 &
              &         : iSrcElems*fun%nComponents )
          end do

          ! get the state value for each component of this
          res( (iVal-1)*fun%nComponents+1: iVal*fun%nComponents ) = pntVal
        end if !idx>0
      end do !iVal
    end if ! idx_len

  end subroutine mus_stateVar_fromIndex
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> Routine to get the actual value for a given array of indices.
  !! The indices belong to the grwarray of points storing levelwise in
  !! Pointdata%pntLvl(iLevel).
  !! Hence this routines takes the indeices as input, can refer to the pointData
  !! and evaluate the variable and returns the values
  subroutine mus_stateVar_Fetch_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                      idxLen,  nVals,  res            )
    ! -------------------------------------------------------------------------- !
    !> Description of the method to obtain the variables,
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: n
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: dependes on number of first index for contiguous array,
    !! but the sum of all idxLen is equal to nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: iComp, iSrc, iVal, elemPos, first, last, nSize, loc_level
    integer :: nSrcElems, iFld, QQ, nScalars
    real(kind=rk), allocatable :: srcRes(:), pntVal(:)
    real(kind=rk) :: weight
    ! -------------------------------------------------------------------------- !

!    write(dbgUnit(4),*) 'get the values of indices for variable ',  &
!      &                  trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nSCalars
    allocate(srcRes(QQ*fun%nComponents))
    allocate(pntVal(fun%nComponents))
    iFld = fun%input_varPos(1)

    res = 0.0_rk

    ! distinguish if we have an array of index or we have contingous memory
    ! access where index are always first entries!
    if (present(idxLen)) then
      call tem_abort('Error: idxLen is not supported in get_valOfIndex for ' &
        &          //'state variable')
    else
      ! Get the state value at the specific point
      do iVal = 1, nVals
        if (idx(iVal)>0) then

          ! elemPos in tree
          ! elemPos = fPtr%pointData%pntLvl(iLevel)%elemPos%val(idx(iVal))
          first = fPtr%pointData%pntLvl(iLevel)%srcElem%first%val(idx(iVal))
          last = fPtr%pointData%pntLvl(iLevel)%srcElem%last%val(idx(iVal))
          loc_level = fPtr%pointData%pntLvl(iLevel)%pntLevel%val(idx(iVal))
          nSize = scheme%pdf( loc_level )%nSize

          ! get pdf's of source elements
          srcRes = 0.0_rk
          nSrcElems = 0
          do iSrc = first, last
            nSrcElems = nSrcElems + 1

            ! position of element in levelDesc total list
            elemPos = fPtr%pointData%pntLvl(iLevel)%srcElem%elemPos%val(iSrc)

            do iComp = 1, fun%nComponents
              srcRes( (nSrcElems-1)*fun%nComponents + iComp )                &
                & = scheme%state( loc_level )%val(                           &
& scheme%pdf(loc_level)%neigh((icomp-1)*nsize+elempos)+(ifld-1)*qq+nscalars*0, &
& scheme%pdf( loc_level )%nNext )
            end do !iComp
          end do !iSrc

          ! Linear interpolation res = sum(weight_i*phi_i)
          pntVal = 0.0_rk
          nSrcElems = 0
          do iSrc = first, last
            weight = fPtr%pointData%pntLvl(iLevel)%srcElem%weight%val(iSrc)
            nSrcElems = nSrcElems + 1
            pntVal(:) = pntVal(:) + weight &
              & * srcRes( (nSrcElems-1)*fun%nComponents+1 &
              &         : nSrcElems*fun%nComponents )
          end do

          ! get the state value for each component of this
          ! get the state value for each component of this
          res( (iVal-1)*fun%nComponents+1: iVal*fun%nComponents ) = pntVal

        end if !idx>0
      end do !iVal
    end if ! idx_len

  end subroutine mus_stateVar_Fetch_fromIndex
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> Routine to get the actual value for a given array of indices.
  !! The indices belong to the grwarray of points storing levelwise in
  !! Pointdata%pntLvl(iLevel).
  !! Hence this routines takes the indeices as input, can refer to the pointData
  !! and evaluate the variable and returns the values
  subroutine mus_stateVar_Fetch_now_fromIndex( fun, varSys, time, iLevel,  &
    &                                          idx, idxLen,  nVals,  res   )
    ! -------------------------------------------------------------------------- !
    !> Description of the method to obtain the variables,
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: n
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: dependes on number of first index for contiguous array,
    !! but the sum of all idxLen is equal to nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: iComp, iSrc, iVal, elemPos, first, last, nSize, loc_level
    integer :: nSrcElems, iFld, QQ, nScalars
    real(kind=rk), allocatable :: srcRes(:), pntVal(:)
    real(kind=rk) :: weight
    ! -------------------------------------------------------------------------- !

!    write(dbgUnit(4),*) 'get the values of indices for variable ',  &
!      &                  trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nSCalars
    allocate(srcRes(QQ*fun%nComponents))
    allocate(pntVal(fun%nComponents))
    iFld = fun%input_varPos(1)

    res = 0.0_rk

    ! distinguish if we have an array of index or we have contingous memory
    ! access where index are always first entries!
    if (present(idxLen)) then
      call tem_abort('Error: idxLen is not supported in get_valOfIndex for ' &
        &          //'state variable')
    else
      ! Get the state value at the specific point
      do iVal = 1, nVals
        if (idx(iVal)>0) then

          ! elemPos in tree
          ! elemPos = fPtr%pointData%pntLvl(iLevel)%elemPos%val(idx(iVal))
          first = fPtr%pointData%pntLvl(iLevel)%srcElem%first%val(idx(iVal))
          last = fPtr%pointData%pntLvl(iLevel)%srcElem%last%val(idx(iVal))
          loc_level = fPtr%pointData%pntLvl(iLevel)%pntLevel%val(idx(iVal))
          nSize = scheme%pdf( loc_level )%nSize

          ! get pdf's of source elements
          srcRes = 0.0_rk
          nSrcElems = 0
          do iSrc = first, last
            nSrcElems = nSrcElems + 1

            ! position of element in levelDesc total list
            elemPos = fPtr%pointData%pntLvl(iLevel)%srcElem%elemPos%val(iSrc)

            do iComp = 1, fun%nComponents
              srcRes( (nSrcElems-1)*fun%nComponents + iComp )          &
                & = scheme%state( loc_level )%val(                        &
& scheme%pdf(loc_level)%neigh((icomp-1)*nsize+elempos)+(ifld-1)*qq+nscalars*0, &
& scheme%pdf( loc_level )%nNow )
            end do !iComp
          end do !iSrc

          ! Linear interpolation res = sum(weight_i*phi_i)
          pntVal = 0.0_rk
          nSrcElems = 0
          do iSrc = first, last
            weight = fPtr%pointData%pntLvl(iLevel)%srcElem%weight%val(iSrc)
            nSrcElems = nSrcElems + 1
            pntVal(:) = pntVal(:) + weight &
              & * srcRes( (nSrcElems-1)*fun%nComponents+1 &
              &         : nSrcElems*fun%nComponents )
          end do

          ! get the state value for each component of this
          ! get the state value for each component of this
          res( (iVal-1)*fun%nComponents+1: iVal*fun%nComponents ) = pntVal

        end if !idx>0
      end do !iVal
    end if ! idx_len

  end subroutine mus_stateVar_Fetch_now_fromIndex
  ! ************************************************************************* !

end module mus_stateVar_module
! **************************************************************************** !

