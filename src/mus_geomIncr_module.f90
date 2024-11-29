! Copyright (c) 2012, 2014-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2018 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
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
!> author: Kartik Jain
!! author: Jiaxing Qi
!! This module provides functionality to perform geometry change
!! depending on the flow properties. Fluidify/Solidify, proximity
!! condition variable are supplied in the lua
!! file. Checks are performed against those parameters and geometry is
!! accordingly changed (an element solidified or fluidified)
!!
!> Geometry Increment module
!!
module mus_geomIncr_module

  ! include treelm modules
  use env_module,              only: rk, long_k, io_buffer_size, LabelLen
  use tem_topology_module,     only: tem_LevelOf
  use tem_property_module,     only: prp_fluidify, prp_solid, prp_fluid, &
    &                                prp_noSolidification
  use tem_comm_module,         only: tem_commPattern_type
  use tem_timeControl_module,  only: tem_timeControl_check
  use tem_logging_module,      only: logUnit
  use tem_general_module,      only: tem_general_type
  use tem_dyn_array_module,    only: PositionOfVal
  use tem_varSys_module,       only: tem_varSys_type
  use tem_construction_module, only: tem_levelDesc_type
  use treelmesh_module,        only: treelmesh_type
  use tem_aux_module,          only: tem_abort

  ! include musubi modules
  use mus_geomIncrHead_module, only: mus_geomIncrHead_type
  use mus_geom_module,         only: mus_geom_type
  use mus_scheme_type_module,  only: mus_scheme_type
  use mus_time_module,         only: mus_timeControl_homogenize



  implicit none

  private

  public :: mus_geomIncr, proximity
  public :: mus_init_geomIncr
  public :: prepare_target_elem_list
  ! public :: update_connectivity

contains

! ****************************************************************************** !
  !> This subroutine initializes the geometry increment.
  !!
  subroutine mus_init_geomIncr(me, varSys, dt, reqInt)
    ! ---------------------------------------------------------------------------
    !> geometry increment type
    type( mus_geomIncrHead_type ), intent(inout) :: me(:)
    !> Global variable system
    type( tem_varSys_type ), intent(in)          :: varSys
    !> dt of maxlevel or smallest dt
    real(kind=rk), intent(in) :: dt
    !> Required interval, in which the update MUST occur.
    !! This is required for the musubi multilevel, where the time step should
    !! only be determined active, when the end of the largest cycle is reached.
    integer, intent(in)                          :: reqInt
    ! ---------------------------------------------------------------------------
    integer :: iGInc
    ! ---------------------------------------------------------------------------
    do iGInc = 1, size(me)
      ! Initialize the timer
      call mus_timeControl_homogenize( me     = me( iGInc )%timeControl, &
        &                              dt     = dt,                      &
        &                              reqInt = reqInt                   )

      me(iGInc)%cond_varPos = PositionOfVal(me  = varSys%varname,             &
        &                                   val = trim(me(iGInc)%cond_varName))
      if (me(iGInc)%cond_varPos < 1) then
        write(logUnit(1),*) 'Error: variable '//trim(me(iGInc)%cond_varName) &
          &              // ' is not added to variable system'
        call tem_abort()
      end if

    end do

  end subroutine mus_init_geomIncr
! ******************************************************************************!


! ******************************************************************************!
  !> This subroutine checks for various conditions defined in the geomIncr table
  !! within the lua file, calculates the requested macroscpoic variables and
  !! then compares them against the specified threshold.
  !! Accordingly then solidification or fluidification of elements is performed.
  !!
  subroutine mus_geomIncr( geometry, scheme, commPattern, general)
    ! ---------------------------------------------------------------------------
    !>
    type( mus_geom_type )                     :: geometry
    !>
    type( mus_scheme_type ), intent(inout)    :: scheme
    !>
    type( tem_commPattern_type ), intent(in)  :: commPattern
    !> global parameters
    type( tem_general_type ),intent(in)       :: general
    ! ---------------------------------------------------------------------------
    integer :: pot_nElems ! potential elements for gi temporarily
    integer :: nElems_gi  ! number of elements which participate in geomincr
    integer :: nElems ! Total number of elements in the tree
    integer :: varPos
    ! levelPointer list of element to be solidify
    integer, allocatable :: LPlist(:), pntTreeID(:)
    ! List of target iElems
    integer, allocatable :: target_ielem_list(:)
    real(kind=rk), allocatable :: chunk(:)
    integer :: iGInc, iChunk
    integer :: nComp, chunkSize, nChunks
    integer :: nChunkElems, iLevel
    integer :: elemoff ! starting index of iElem for different chunks
    real(kind=rk), allocatable :: chunk_res(:)
    integer :: buf_start, buf_end
    integer :: res_size
    logical :: triggered
    ! ---------------------------------------------------------------------------

    ! The number of elements in the tree
    nElems = geometry%tree%nElems

    chunkSize = io_buffer_size

    allocate( LPlist( size(geometry%tree%treeID)) )
    allocate( pntTreeID( size(geometry%tree%treeID)) )
    allocate( target_ielem_list( size(geometry%tree%treeID)) )

    allocate(chunk(io_buffer_size))

    ! Loop over the number of geomIncrs defined in the Lua file
    do iGInc = 1, size( geometry%geomIncr )
      call tem_timeControl_check(                                         &
        &               me        = geometry%geomIncr(iGInc)%timeControl, &
        &               now       = general%simControl%now,               &
        &               comm      = general%proc%comm,                    &
        &               triggered = triggered                             )

      ! time_modulo is to do only when one complete cyle is done with multilevel
      ! when read from restart trigger is active at 1st time step
      if( triggered ) then

        ! Initialize the total number of target elements to zero
        nElems_gi = 0
        ! Call the proximity routine to generate
        ! a list of treeIDs which pass the proximity check
        call proximity( scheme       = scheme,      &
          &             geometry     = geometry,    &
          &             iGInc        = iGInc,       &
          &             pot_nElems   = pot_nElems,  &
          &             pot_lp       = LPlist(:),   &
          &             pot_lp_tree  = pntTreeID(:) )

        nChunks = ceiling( real(pot_nElems, kind=rk) &
          &              / real(chunkSize, kind=rk) )

        do iChunk = 1, nChunks
          elemOff = ( (iChunk-1)*chunkSize )

          ! number of elements written to THIS chunk
          nChunkElems = min(chunkSize, pot_nElems-elemOff)

          ! Compute the element lower and upper bound for the current chunk
          buf_start = (iChunk-1)*chunkSize+1
          buf_end   = (iChunk-1)*chunkSize+nChunkElems

          ! position of condition variable in varSys
          varPos = geometry%geomIncr(iGInc)%cond_varPos
          ! get the no. of components for each variable system
          nComp = scheme%varSys%method%val(varPos)%nComponents

          ! allocate local array where the results are placed
          res_size = nComp*nChunkElems
          allocate( chunk_res(res_size) )

          ! Use the function pointers to calculate variables identified
          ! from Lua results are placed in array chunk_res.
          ! Variable returns 1 if condition satisfies else 0
          call scheme%varSys%method%val(varpos)%get_element(         &
            &               varSys  = scheme%varSys,                 &
            &               elemPos = pntTreeID(buf_start:buf_end),  &
            &               time    = general%simControl%now, &
            &               tree    = geometry%tree,                 &
            &               nElems  = nChunkElems,                   &
            &               nDofs   = 1,                             &
            &               res     = chunk_res(1:res_size)          )

          ! Prepare the list of target elements and set the property bits
          ! of these elements
          call prepare_target_elem_list(                    &
            & nComp             = nComp,                    &
            & elemoff           = elemoff,                  &
            & nChunkElems       = nChunkElems,              &
            & res               = chunk_res,                &
            & LPList            = LPList,                   &
            & nElems_gi         = nElems_gi,                &
            & pntTreeID         = pntTreeID,                &
            & target_ielem_list = target_ielem_list,        &
            & scheme            = scheme,                   &
            & tree              = geometry%tree,            &
            & geomIncr          = geometry%geomIncr(iGInc)  )

          ! Deallocate the array containing computed results
          deallocate( chunk_res )
        end do ! iChunk

        if(nElems_gi > 0) then
          if (geometry%geomIncr(iGInc)%solidify) then
            write(logUnit(1),*)'Elements found for solidification:', nElems_gi
          elseif (geometry%geomIncr(iGInc)%fluidify) then
            write(logUnit(1),*)'Elements found for fluidification:', nElems_gi
          end if
          nElems_gi = 0
          ! and update the flag meshChange because something changed
          geometry%tree%global%meshChange = .true.
        end if

        do iLevel = geometry%tree%global%minLevel, geometry%tree%global%maxLevel
          ! Communicate the property bits of halo elements
          call commPattern%exchange_long(                            &
            &  send         = scheme%levelDesc( iLevel )%sendbuffer, &
            &  recv         = scheme%levelDesc( iLevel )%recvbuffer, &
            &  state        = scheme%levelDesc( iLevel )%property,   &
            &  message_flag = iLevel,                                &
            &  comm         = general%proc%comm                    )

          ! Invoke the update connectivity routine which updates the
          ! neighbor array of target elements and the neighbor array of
          ! their neighbors.
          call update_connectivity( scheme, iLevel)
        end do
      end if ! Timing check
    end do ! iGInc
    ! Deallocate the chunk array
    deallocate( chunk )
    deallocate( LPlist )
    deallocate( target_ielem_list )
  end subroutine mus_geomIncr
! ****************************************************************************** !


! ****************************************************************************** !
  !> A subroutine which checks the proximity condition and generates a new
  !! treeID list of potential target elements which can then be checked against
  !! various thresholds as defined in the lua configuration file.
  !!
  !! Proximity condition means that an element can be solidified only if it has
  !! at least one wall as a neighboring element.
  !! In biological clotting models it ensures that no isolated clot is produced
  !! in the fluid domain.
  !! Furthermore the noSolidification property is checked for every element.
  !! It is used to define spaces where solidification is not allowed.
  !!
  subroutine proximity( scheme, geometry, iGInc, pot_nELems, pot_lp,      &
    &                   pot_lp_tree)
    ! ---------------------------------------------------------------------------
    !> scheme information
    type(mus_scheme_type), intent(inout)  :: scheme
    !> scheme information
    type(mus_geom_type), intent(in)       :: geometry
    !> Geomtry icrement
    integer, intent(in)                   :: iGInc
    !> number of potential elements found
    integer, intent(inout)                :: pot_nElems
    !> level pointer of potential elements
    integer, intent(inout)                :: pot_lp(:)
    !> position of potential elements in the treeID list
    integer, intent(inout)                :: pot_lp_tree(:)
    ! ---------------------------------------------------------------------------
    ! local variables
    integer :: iElem, iLevel, iField, iDir
    integer :: neighPos, elemPos
    integer :: QQN, QQ, neighbor, nScalars
    integer(kind=long_k) :: elemProp
    ! ---------------------------------------------------------------------------

    pot_nElems = 0
    pot_lp = 0
    pot_lp_tree = 0
    QQN = scheme%layout%fStencil%QQN
    QQ  = scheme%layout%fStencil%QQ
    nScalars = scheme%varSys%nScalars

    if ( geometry%geomIncr(iGInc)%proximity ) then
      do iElem = 1, size( geometry%tree%treeID )
        elemProp = geometry%tree%ElemPropertyBits( iElem )
        ! Only consider element if solidification is allowed
        if ( .NOT. btest ( elemProp, prp_noSolidification ) ) then
          do iField = 1, scheme%nFields
            iLevel = tem_levelOf( geometry%tree%treeID( iElem ))
            elemPos = geometry%levelPointer( iElem )
            do iDir = 1, QQN
              ! Get the neighbor
              neighbor = scheme%pdf( iLevel )%neigh(               &
                & (idir-1)* scheme%pdf(ilevel)%nsize+ elempos )
              ! position in total list
              neighPos = int(( neighbor-1)/ nscalars)+1

              !> Acc the proximity condition, as soon as its found that any of the
              !! neighbors of current element is solid, that element is stored in
              !! the new list of potential elements. Loop is terminated for this
              !! element and check is performed for other element.
              if (neighPos == elemPos) then
                pot_nElems              = pot_nElems + 1
                pot_lp(pot_nElems)      = elemPos
                pot_lp_tree(pot_nElems) = iElem
                exit
              end if
            end do ! iDir
          end do ! iField
        end if
      end do ! iElem
    else
      do iElem = 1, size( geometry%tree%treeID )
        elemProp = geometry%tree%ElemPropertyBits( iElem )
        if ( .NOT. btest ( elemProp, prp_noSolidification ) ) then
          elemPos = geometry%levelPointer( iElem )
          pot_nElems              = pot_nElems + 1
          pot_lp(pot_nElems)      = elemPos
          pot_lp_tree(pot_nElems) = iElem
        end if
      end do
    end if

  end subroutine proximity
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine compares the macroscopic quantity obtained for an element
  !! against specified threshold and populates the list of target elements
  !! which participate in geometry change. In addition it sets the property
  !! bit(s) of target elements to fluidify in case of solidification and clears
  !! this bit in case of fluidification
  !!
  subroutine prepare_target_elem_list( nComp, elemoff, nChunkElems, res, &
    &                                  LPList, nElems_gi, pntTreeID,     &
    &                                  target_ielem_list, scheme, tree,  &
    &                                  geomIncr                          )
    ! ---------------------------------------------------------------------------
    !> total number of components
    integer, intent(in)    :: nComp
    !> offset of elements in chunk
    integer, intent(in)    :: elemoff
    !> number of elements in chunk
    integer, intent(in)    :: nChunkElems
    !> result from condition variable
    real(kind=rk), intent(in) :: res(:)
    !> list of level pointers
    integer, intent(in)    :: LPList(:)
    !> number of final elements
    integer, intent(inout) :: nElems_gi
    !> element position in global treeID list tree%treeID
    integer, intent(in)    :: pntTreeID(:)
    !> target element list for second iGInc to reduce computation effort
    integer, intent(inout) :: target_ielem_list(:)
    !> scheme information
    type(mus_scheme_type), intent(inout)  :: scheme
    !> fluid tree from mesh, element property bits are set to solidify
    !! if condition variable is 1
    type( treelmesh_type ), intent(inout)  :: tree
    !> current geometry increment
    type( mus_geomIncrHead_type ), intent(inout) :: geomIncr
    ! ---------------------------------------------------------------------------
    ! local variables
    integer                :: local_elem
    integer                :: iElem, iLevel, elemPos, neighPos, iField, iDir
    integer(kind=long_k)   :: neighProp
    logical                :: isSolid, cond
    ! ---------------------------------------------------------------------------

    ! Loop over all elements in chunks for the condition variable in var system
    ! For solidify (solidify=true), fluid elements will be solidified if the
    ! condition is satisfied and for fluidify (fluidify=true) solid elements
    ! will be fluidified.
      local_elem = 0
      do iElem = elemoff+1, elemoff+nChunkElems
        local_elem = local_elem + 1
        elemPos = LPlist( iElem ) ! elem position in total list
        iLevel = tem_LevelOf( tree%treeID( elemPos ))
        if( geomIncr%solidify ) then
          ! if condition variable is a vector then all index must be satisfied
          cond = all(res( (local_elem-1)*nComp + 1 : local_elem*nComp) > 0.0_rk)
          if ( cond ) then
            ! Update the property bit
            call mus_setProp( property = scheme%levelDesc(iLevel) &
              &                          %property(elemPos),      &
              &               geomIncr = geomIncr                 )
            ! and the ones in the global tree
            call mus_setProp( property = tree%ElemPropertyBits &
              &                          (elemPos),            &
              &               geomIncr = geomIncr              )
            nElems_gi = nElems_gi + 1
          end if ! condition satisfies
        elseif ( geomIncr%fluidify ) then
          ! if condition variable is a vector then all index must be satisfied
          cond = all(res( (local_elem-1)*nComp + 1 : local_elem*nComp) > 0.0_rk)
          if ( cond ) then
            ! Get the direct neighbor elements for fluidification
            do iField = 1, scheme%nFields
              do iDir = 1, 6
                neighPos = scheme%levelDesc( iLevel )%neigh(1)%               &
                  &        nghElems(  scheme%layout%fstencil %cxdirinv( idir), &
                  &                  elemPos )
                ! neighPos < 0 means boundary in that direction
                if ( neighPos > 0 ) then
                  neighProp = scheme%levelDesc(iLevel)%property( neighPos )

                  isSolid = btest( neighProp, prp_Solid )

                  ! Only need to fluidify if neighbor is solid
                  if( isSolid ) then
                    ! Update the property bit
                    call mus_setProp( property = scheme%levelDesc(iLevel) &
                      &                          %property(neighPos),     &
                      &               geomIncr = geomIncr                 )
                    ! and the ones in the global tree
                    call mus_setProp( property = tree%ElemPropertyBits  &
                      &                          (neighPos),            &
                      &               geomIncr = geomIncr               )
                    nElems_gi = nElems_gi + 1
                    ! As soon as one solid neighbor is fluidifed
                    ! dir-loop is terminated. This will keep a balance
                    ! between soldification and fluidification.
          !!          EXIT
                  end if
                end if
              end do ! iDir
            end do ! iField
          end if ! condition satisfies
        end if ! is my target
      end do ! iElem

! DF: Following code is only used to improve performance.
!     Still needs a couple of changes to run with
!     solidification AND fluidification.

!!    else ! iVar != 1 or iDepend != 1
!!      ! loop over the tree ID list which was found for the first geomIncr loop.
!!      ! Doing this we ensure the AND condition among
!!      ! different variables within a depend table and as the number of target
!!      ! elements is hardly 5% of the total geometry, it improves computational
!!      ! efficiency of the code.
!!      local_elem = 0
!!      do iElem = 1, nElems_gi
!!        elemPos = LPlist( target_ielem_list( iElem ))
!!        iLevel = tem_LevelOf( tree%treeID( &
!!          &                   pntTreeID( target_ielem_list(iElem) ) ) )
!!        elemProp = levelDesc( iLevel )%property( elemPos )
!!        isSolid = btest( elemProp, prp_fluidify )
!!        if( isSolid .neqv. geomIncr%solidify ) then
!!          ! if condition variable is a vector then all index must be satisfied
!!          cond = all(res( (target_ielem_list(iElem)-1)*nComp + 1 &
!!            &            : target_ielem_list(iElem)*nComp) > 0.0_rk)
!!
!!          ! Update the local tree ID list
!!          if ( cond ) then
!!            ! Update the property bit
!!            call mus_setSolidProp(                                    &
!!              &       property = levelDesc(iLevel)%property(elemPos), &
!!              &       solidify = geomIncr%solidify                    )
!!            ! and the ones in the global tree
!!            call mus_setSolidProp(                             &
!!              &   property = tree%ElemPropertyBits(pntTreeID(iElem)), &
!!              &   solidify = geomIncr%solidify                        )
!!            ! and update the flag meshChange because something changed
!!            ! geometry%tree%global%meshChange = .true.
!!            local_elem = local_elem + 1
!!            target_ielem_list(local_elem) = target_ielem_list(iElem)
!!          end if ! condition satisfies
!!        end if ! is my target
!!      end do ! iElem
!!      ! number of target elements is updated
!!      nElems_gi = local_elem
!!    end if ! nElems_gi == 0

  end subroutine prepare_target_elem_list
! ****************************************************************************** !


! ****************************************************************************** !
  !> Construct the propagation stencil for each element
  !!
  !! the bounce back rules have to be applied here
  !!
  subroutine update_connectivity( scheme, iLevel)
    ! --------------------------------------------------------------------------
    !> scheme information
    type( mus_scheme_type ), intent(inout) :: scheme
    !> current level
    integer, intent(in)     :: iLevel
    ! --------------------------------------------------------------------------
    ! defining local variables
    integer :: iElem  ! element counter for current element in treeID list
    integer :: iField ! field counter
    integer :: idx_idir ! current direction in global varsys
    integer :: QQ
    integer :: iDir ! counter for current direction in iElem
    integer :: invDir ! inverse direction for iDir
    ! result: from which direction to take pdf from (inverse for bounceback)
    integer :: GetFromDir
    ! result: from which element to take the pdf from (local for bounceback)
    integer :: GetFromPos
    integer :: neighPos      ! position of current neighElem in treeID list
    integer(kind=long_k)  :: neighProp, elemProp
    integer :: nElems, nSize
    ! ---------------------------------------------------------------------------
    QQ = scheme%layout%fStencil%QQ

    nElems = scheme%pdf( iLevel )%nElems_local
    nSize  = scheme%pdf( iLevel )%nSize
    do iElem = 1, nElems
      elemProp = scheme%levelDesc( iLevel )%property( iElem )

      ! Loop over every link
      do iDir  = 1, QQ
        if ( iDir == QQ  ) then
          ! Rest density takes from itself, iElem and iDir
          GetFromPos = iElem
          GetFromDir = iDir
        else ! iDir not rest density
          ! For push we need to find different neighbors than for pull
          ! -> ngDir
          neighPos = scheme%levelDesc( iLevel )%neigh(1)%               &
            &        nghElems(  scheme%layout%fstencil %cxdirinv( idir),    &
            &                  iElem )

          neighProp = 0_long_k
          ! get neighbor's property
          if( neighPos > 0 )                                                 &
            & neighProp = scheme%levelDesc( iLevel )%property( neighPos )

          if( (neighPos <= 0) .or. (btest(neighProp, prp_fluidify))        &
              &  .or. (btest(elemProp, prp_fluidify))) then
            ! 1. Fluid element has boundary in current link direction
            ! 2. neighbor element has fluidify property (still solid now)
            ! Get inverse direction of LOCAL element
            invDir = scheme%layout%fStencil%cxDirInv( iDir )
            GetFromDir = invDir
            GetFromPos = iElem
            ! I set the source to the target link for halo cells
            ! with no neighbor in that direction
          else
            ! neighElem = scheme%levelDesc( iLevel )%total( neighPos )
            GetFromDir = iDir
            GetFromPos = neighPos
          endif ! has no boundary -> regular treatment
        endif ! not the 0 velocity

        scheme%pdf( iLevel )%neigh( ( idir-1)* nelems + ielem) =    &
          &  ( getfrompos-1)* scheme%varsys%nscalars+ getfromdir

        do iField = 1, scheme%nFields
          idx_idir = scheme%varSys%method%val(                      &
            &               scheme%stateVarMap%varPos%val(iField) ) &
            &                     %state_varPos(iDir)
!!          ! Set all the pdf links of target element to an infinitesimally
!!          ! small value. This will prevent density from becoming NaN if pdfs
!!          ! are set to zero and besides changed elements will have a very
!!          ! small density which will help in visualization and differentiate
!!          ! them from other fluid elements
!!          if (btest(elemProp, prp_fluidify) ) then
!!             scheme%state(iLevel)%val(                                         &
!!               &   ( ielem-1)* scheme%varsys%nscalars+idx_idir,:) =  &
!!               &   0.000000001_rk
!!          end if
        end do

      enddo ! neighLoop
    enddo ! elemLoop

  end subroutine update_connectivity
! ****************************************************************************** !

! ****************************************************************************** !
  !> This routine updates the propertybits of an element
  !! If solidify == true, set prp_fluidify and set prp_solid
  !! If fluidify == treu, clean prp_fluidify and set prp_fluid
  subroutine mus_setProp( property, geomIncr )
    ! ---------------------------------------------------------------------------
    integer( kind=long_k ) :: property
    type( mus_geomIncrHead_type ), intent(in) :: geomIncr
    ! ---------------------------------------------------------------------------

    if( geomIncr%solidify ) then
      ! set solid property and clean fluid property
      property = ibset( property, prp_fluidify )
      property = ibset( property, prp_solid )
      property = ibclr( property, prp_fluid )
    elseif( geomIncr%fluidify ) then
      ! set fluid property and clean solid property
      property = ibclr( property, prp_fluidify )
      property = ibset( property, prp_fluid )
      property = ibclr( property, prp_solid )
    endif

  end subroutine mus_setProp
! ****************************************************************************** !


end module mus_geomIncr_module
! ****************************************************************************** !
