! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012, 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
! **************************************************************************** !
!> Collection of debugging functions
!!
!!
module mus_debug_module

  ! include treelm subroutines
  use env_module,              only: long_k
  use tem_param_module,        only: qDirName
  use tem_construction_module, only: tem_levelDesc_type
  use treelmesh_module,        only: treelmesh_type
  use tem_bc_prop_module,      only: tem_bc_prop_type
  use tem_element_module,      only: eT_ghostFromCoarser, &
    &                                eT_ghostFromFiner
  use tem_property_module,     only: prp_hasBnd
  use tem_topology_module,     only: tem_CoordOfId, tem_IdOfCoord, &
    &                                tem_firstIdAtLevel
  use tem_geometry_module,     only: tem_baryOfId
  use tem_debug_module,        only: dbgUnit, main_debug
  use tem_tools_module,        only: tem_horizontalSpacer
  use tem_logging_module,      only: logUnit, tem_toStr
  use tem_time_module,         only: tem_time_advance

  ! include musubi subroutines
  use mus_param_module,              only: mus_param_type
  use mus_geom_module,               only: mus_geom_type
  ! use mus_comm_module,               only: mus_exchange
  use mus_bc_general_module,         only: set_boundary
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_debug_tools_module,        only: dumpPdfAll

  implicit none
  private

  public :: debug_dependenciesFromFiner
  public :: debug_dependenciesFromCoarser
  public :: debug_connectivity
  public :: debug_normals
  public :: dump_debug_info

  contains

! **************************************************************************** !
  !> Dump pdf values into debug files
  !!
  subroutine dump_debug_info( buffer, scheme, params, iLevel, iTime )
    ! --------------------------------------------------------------------------
    character(len=*), intent(in)  :: buffer
    !> scheme type contains pdf
    type( mus_scheme_type ), intent(in) :: scheme
    !> global parameters
    type( mus_param_type ), intent(in) :: params
    !> Level counter variable
    integer, intent(in)           :: iLevel
    !> integer time counter variable
    integer, intent(in), optional :: iTime
    ! --------------------------------------------------------------------------
    if ( main_debug%dumpLevelwiseState ) then
      ! dump current time and level
      write(dbgUnit(1),*) trim( buffer ), params%general%simControl%now%sim
      write(dbgUnit(1),*)'            level      ', iLevel

      ! dump pdf values
      call dumpPDFAll( pdf       = scheme%pdf(iLevel),              &
        &              levelDesc = scheme%leveldesc(iLevel),        &
        &              QQ        = scheme%layout%fStencil%QQ,       &
        &              text      = buffer,                          &
        &              iTime     = iTime,                           &
        &              state     = scheme%state(iLevel)%val,        &
        &              nFields   = scheme%nFields,                  &
        &              nScalars  = scheme%varSys%nScalars,          &
        &              dumpHalos = main_debug%dumpHaloState )
    end if
  end subroutine dump_debug_info
! **************************************************************************** !


! **************************************************************************** !
  !> Debug the dependencies for the ghostFromFiner elements
  !!
  subroutine debug_dependenciesFromFiner( levelDesc, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    !> Level range
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    type( tem_levelDesc_type ),  intent(in) :: levelDesc(minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: sourceLevel  ! level of source elements
    integer :: targetElem   ! treeId of current source element
    integer :: iElem, iLevel
    integer :: nSources ! number of source elements for the current target
    integer :: mySources(8)
    ! --------------------------------------------------------------------------

    call tem_horizontalSpacer( dbgUnit(5) )
    write(dbgUnit(5),'(A)') 'Inside routine: debug_dependenciesFromFiner'

    do iLevel = minLevel, maxLevel-1

write(dbgUnit(5),"(A,I0)") 'Ghost on level: ', iLevel
write(dbgUnit(5),"(A,I0)") 'Number of Ghost from Finer: ', &
  &                        levelDesc( iLevel )%elem%nElems( eT_ghostFromFiner )
write(dbgUnit(5),"(4A12)") 'iElem', 'ghostID', 'nSources', 'sourceID'

      ! Treat every coarse target element:
      do iElem = 1, levelDesc( iLevel )%elem%nElems( eT_ghostFromFiner )

        ! Check source level
        sourceLevel = levelDesc( iLevel )%depFromFiner( iElem )%dependencyLevel
        if( sourceLevel - iLevel /= 1 ) then
          write(dbgUnit(5),*)'error in level difference on level ', iLevel, &
            &                ' elem ', iElem, ' sourceLevel ', sourceLevel
          stop
        endif

        targetElem = iElem + levelDesc( iLevel )%offset( 1, eT_ghostFromFiner)
        nSources = levelDesc( iLevel )%depFromFiner( iElem )%elem%nVals
        mySources(1:nSources) = &
          & levelDesc( iLevel )%depFromFiner( iElem )%elem%val(1:nSources)

! write(dbgUnit(1),"(A,3F6.3)")  'coord: ', &
!   &   levelDesc( iLevel )%depFromFiner( iElem )%coord(1:3)
write(dbgUnit(5),'(3I12,A)') &
  &  iElem, levelDesc( iLevel )%total( targetElem ), nSources, &
  &  trim(tem_toStr( levelDesc(sourceLevel)%total(mySources(1:nSources)),','))

      enddo ! iElem
    enddo ! iLevel

    write(dbgUnit(5),"(A)") 'Out of routine: debug_dependenciesFromFiner'
    call tem_horizontalSpacer( dbgUnit(5) )

  end subroutine debug_dependenciesFromFiner
! **************************************************************************** !


! **************************************************************************** !
  !> Interpolation routine that is based on a simple weighted average of source
  !! nodes. This is the interpolation coarse-> fine. Weights are needed here,
  !! as the distances source <-> target are different for the source nodes.
  !!
  subroutine debug_dependenciesFromCoarser( levelDesc, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    !> Level range
    integer,                    intent(in) :: minLevel, maxLevel
    !> Level descriptor
    type( tem_levelDesc_type ), intent(in) :: levelDesc( minlevel:maxLevel )
    ! --------------------------------------------------------------------------
    integer :: sourceLevel   ! level of source elements
    integer :: targetElem    ! treeId of current source element
    integer :: iElem         ! current target element (for outer loop)
    integer :: nSources  ! number of source elements for the current target
    integer :: iLevel
    integer :: mySources(27)
    ! --------------------------------------------------------------------------

    call tem_horizontalSpacer( dbgUnit(5) )
    write(dbgUnit(5),"(A)") 'Inside routine: debug_dependenciesFromCoarser'

    do iLevel = minLevel+1, maxLevel

      sourceLevel = iLevel - 1

write(dbgUnit(5),"(A,I0)") 'Ghost on level: ', iLevel
write(dbgUnit(5),"(A,I0)") 'Number of Ghost from Coarser: ', &
  &                         levelDesc( iLevel )%elem         &
  &                           %nElems( eT_ghostFromCoarser )
write(dbgUnit(5),"(5A12)") 'iElem', 'ghostID', 'nSources', 'childNum', 'coord'

      do iElem = 1, levelDesc( iLevel )%elem%nElems( eT_ghostFromCoarser )

        ! find out the treeId of the current target element
        targetElem = iElem + levelDesc( iLevel )%offset( 1, eT_ghostFromCoarser)

        nSources = levelDesc( iLevel )%depFromCoarser( iElem )%elem%nVals
        mySources(1:nSources) = &
          & levelDesc( iLevel )%depFromCoarser( iElem )%elem%val(1:nSources)

write(dbgUnit(5),'(4I12,3F6.3)') &
  &  iElem, levelDesc( iLevel )%total( targetElem ), nSources, &
  &  levelDesc( iLevel )%depFromCoarser( iElem )%childNum, &
  &  levelDesc( iLevel )%depFromCoarser( iElem )%coord(1:3)
write(dbgUnit(5),'(2A)')  'sourceID: ', &
  &  trim(tem_toStr( levelDesc(sourceLevel)%total(mySources(1:nSources)),','))

if ( allocated( levelDesc( iLevel )%depFromCoarser( iElem )%weight) ) then
write(dbgUnit(5),'(2A)')  'weights: ', &
  &  trim(tem_toStr(levelDesc(iLevel)%depFromCoarser(iElem) &
  &                   %weight(1:nSources),',',main_debug%logger))
end if
    write(dbgUnit(5),"(A)") ''

      end do ! iElem
    end do ! iLevel

    write(dbgUnit(5),"(A)") 'Out of routine: debug_dependenciesFromCoarser'
    call tem_horizontalSpacer( dbgUnit(5) )
  end subroutine debug_dependenciesFromCoarser
! **************************************************************************** !


! **************************************************************************** !
  !> summary: Dump the neighbors found for the boundaries with neighbor elements
  !!
  subroutine debug_normals( tree, scheme, nBCs )
    ! --------------------------------------------------------------------------
    type( treelmesh_type ),  intent(in) :: tree
    integer,                 intent(in) :: nBCs
    type( mus_scheme_type ), intent(in) :: scheme
    ! -------------------------------------------------------------------------
    integer :: iElem, iLevel, iBnd, iNeigh
    integer :: minLevel, maxLevel
    integer(kind=long_k) :: nTreeID
    ! -------------------------------------------------------------------------

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel

    write(dbgUnit(1),*) 'Debugging Boundary elements normals...'

    do iLevel = minLevel, maxLevel
      do iBnd = 1, nBCs

        write(dbgUnit(1),*) 'BC ID: ', iBnd, ', at Level: ', iLevel, &
          &                 ', nElems: ', scheme%globBC( iBnd )%nElems( iLevel )
        if ( .not. scheme%globBC(iBnd)%isWall ) then

        do iElem = 1, scheme%globBC( iBnd )%nElems( iLevel )

          nTreeID = scheme%levelDesc( iLevel )%total(              &
            &              scheme%globBC( iBnd )%elemLvl( iLevel ) &
            &                %elem%val( iElem ))

          write(dbgUnit(1),*) 'tID ', nTreeID, ' pos ',                        &
            &                  tem_BaryOfId( tree, nTreeID )
          write(dbgUnit(1),*) 'normal ', scheme%globBC( iBnd )%                &
            &                        elemLvl( iLevel )%normal%val( :, iElem ), &
            &                ' nrmInd ', scheme%globBC( iBnd )%                &
            &                        elemLvl( iLevel )%normalInd%val( iElem ), &
            &                ' nNeighs ', scheme%field(1)%bc( iBnd )%nNeighs

          do iNeigh = 1, scheme%field(1)%bc( iBnd )%nNeighs
            write(dbgUnit(1),*) '   ', iNeigh, ' neigh ',                     &
              &   scheme%levelDesc(iLevel)                                    &
              &     %total( scheme%field(1)                                   &
              &               %bc(iBnd)%neigh(iLevel)                         &
              &               %posInState(iNeigh, iElem) ),                   &
              &   'bardofID', tem_baryOfId( tree, scheme%levelDesc( iLevel )  &
              &                                     %total( scheme%           &
              &                   field(1)%bc( iBnd )%neigh( iLevel )%        &
              &                             posInState( iNeigh, iElem )) )
          enddo ! iNeigh
        enddo ! iElem = 1, scheme%globBC( iBnd )%nElems( iLevel )
        else
          write(dbgUnit(5),*) 'this is wall'
        end if ! isWall
      enddo ! iBnd = 1, nBCs
    enddo ! iLevel = minLevel, maxLevel
  end subroutine debug_normals
! **************************************************************************** !


! **************************************************************************** !
  !> Detailed Debug output to the PDF neighbor (connectivity) list
  !!
  !! currently the dubug connectivity is done only for 1st field
  subroutine debug_connectivity( scheme, minLevel, maxLevel)
    ! --------------------------------------------------------------------------
    type(mus_scheme_type), intent(in) :: scheme
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    ! --------------------------------------------------------------------------
    integer :: iElem, iLevel, iDir
    integer :: posInNeigh, posInState, posInTotal
    integer(kind=long_k) :: treeID
    integer :: QQ, varPos(scheme%layout%fStencil%QQ), nScalars
    integer :: nElems, nSize
    ! --------------------------------------------------------------------------

    QQ = scheme%layout%fStencil%QQ
    ! take varpos of 1st field
    nScalars = scheme%varSys%nScalars
    varPos = scheme%varSys%method%val(               &
      &            scheme%stateVarMap%varPos%val(1)) &
      &                  %state_varPos

    write(dbgUnit(1),*)
    write(dbgUnit(1),*) ' Debugging connectivity for the First Field'
    write(dbgUnit(1),*) ' level from: ', minLevel, ' to ', maxLevel

    do iLevel = minLevel, maxLevel

      nSize  = scheme%pdf( iLevel )%nSize
      nElems = scheme%pdf( iLevel )%nElems_local

      write(dbgUnit(1),*)
      write(dbgUnit(1),"(A,I0)") 'LEVEL: ', iLevel
      write(dbgUnit(1),"(A,I0)") 'nElems fluid ', scheme%pdf( iLevel ) &
        &                                               %nElems_fluid
      write(dbgUnit(1),"(A,I0)") 'nElems solve ', scheme%pdf( iLevel ) &
        &                                               %nElems_solve
      write(dbgUnit(1),"(A,I0)") 'nElems halo  ', scheme%pdf( iLevel ) &
        &                                               %nElems_halo
      write(dbgUnit(1),"(A,I0)") 'nElems ghost ', scheme%pdf( iLevel ) &
        &                                               %nElems_ghost
      write(dbgUnit(1),"(A,I0)") 'nElems local ', scheme%pdf( iLevel ) &
        &                                               %nElems_local
      write(dbgUnit(1),"(A,I0)") 'nElems size  ', nSize
      write(dbgUnit(1),*) 'size of total ', &
        &                 size(scheme%levelDesc( iLevel )%total)
      write(dbgUnit(1),*) 'size of state ', &
        &                 size(scheme%state( iLevel )%val(:,1))
      write(dbgUnit(1),*) 'size of neigh ', size(scheme%pdf( iLevel )%neigh )

      write(dbgUnit(1),*) 'Loop over all local elements:'
      write(dbgUnit(1),"(A10,A20,4A8)") 'iElem', 'treeID', 'X', 'Y', 'Z', &
        &                               'level'
      do iElem = 1, nElems

        ! write element info
        write(dbgUnit(1),"(I10,I20,4I8)") iElem, &
          &  scheme%levelDesc( iLevel )%total( iElem ), &
          &  tem_CoordOfId( scheme%levelDesc( iLevel )%total( iElem ))

        write(dbgUnit(1),"(3A10,4A8)") 'varPos', 'qDirName', 'treeID', 'X', &
          &                            'Y', 'Z', 'level'
        do iDir = 1, QQ
          posInNeigh = ( idir-1)* nsize + ielem
          ! pos in State vector for the 1st field
          posInState = scheme%pdf( iLevel )%neigh( posInNeigh )

          if( posInState <= 0 ) then
            write(dbgUnit(1),*) 'ERROR in qDirName: ',qDirName( iDir )
          else
            posInTotal = int((posinstate-1)/ nscalars)+1
            treeID     = scheme%levelDesc( iLevel )%total( posInTotal )
            write(dbgUnit(1),"(I10,A10,I10,4I8)") varPos( iDir ), &
              &        qDirName( iDir ), treeID, tem_CoordOfId( treeID )
          endif
        end do ! iDir

        write(dbgUnit(1),*) ''

      end do ! iElem

    end do ! iLevel

    write(dbgUnit(1),*) ' End Debugging connectivity '

    call tem_horizontalSpacer( fUnit = dbgUnit(1))
  end subroutine debug_connectivity
! **************************************************************************** !


end module mus_debug_module
! **************************************************************************** !
