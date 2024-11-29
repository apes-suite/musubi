! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2015, 2018-2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2021-2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> author: Manuel Hasert
!! Interpolation scheme tools
!!
!! For an overview over implemented interpolation methods, see
!! [Interpolation methods](../page/features/intp_methods.html)
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
module mus_interpolate_tools_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,              only: rk, long_k, newUnit, pathLen
  use treelmesh_module,        only: treelmesh_type
  use tem_grow_array_module,   only: grw_intArray_type
  use tem_construction_module, only: tem_levelDesc_type
  use tem_element_module,      only: eT_GhostFromCoarser, &
    &                                eT_ghostFromFiner
  use tem_geometry_module,     only: tem_baryOfId, &
    &                                tem_elemSize
  use tem_param_module,        only: cs2inv, cs2, PI, div1_2, div1_9, div4_9,  &
    &                                div1_36, div2_3, div1_18, div2_9, div1_18,&
    &                                rho0, rho0Inv
  use tem_matrix_module,       only: tem_matrix_type
  use tem_varSys_module,       only: tem_varSys_type
  use tem_property_module,      only: prp_fluid
  !use tem_construction_module,  only: tem_levelDesc_type
  use tem_debug_module, only: dbgUnit

  ! include musubi modules
  use mus_interpolate_header_module, only: mus_interpolation_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_pdf_module,                only: pdf_data_type
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_varSys_module,        only: mus_varSys_data_type
  use mus_gradData_module,      only: mus_gradData_type
  use mus_scheme_type_module,        only: mus_scheme_type

  implicit none

  private

  public :: dump_intpLists
  public :: debug_dependencies

  contains
! ****************************************************************************** !
  !> check the dependencies from Finer
  !!
  subroutine debug_dependencies( intp, leveldesc, tree, rank )
    ! ---------------------------------------------------------------------------
    !> state properties
    type( tem_leveldesc_type ), intent(in) :: levelDesc(:)
    !> interpolation method info
    type( mus_interpolation_type ), intent(in) :: intp
    !> global tree information
    type( treelmesh_type ), intent(in) :: tree
    !> musubi mpi communicator environment
    integer, intent(in) :: rank
    ! ---------------------------------------------------------------------------
    integer :: nUnit
    character(len=17)  :: buffer
    ! ---------------------------------------------------------------------------
    nUnit = newunit()
    write(buffer,'(a,i7.7)') 'debugDep. rank: ', rank
    open(file = trim(buffer), unit= nUnit, recl=1024)

    write(nUnit,*)
    write(nUnit,*) 'GhostFromFiner'
    write(nUnit,*)
    write(nUnit,*) '-----------------------------------------------------------'
    call dump_MyGhostsFromFiner( intp, levelDesc, nUnit, tree )
    write(nUnit,*) '-----------------------------------------------------------'
    write(nUnit,*) '-----------------------------------------------------------'
    write(nUnit,*)
    write(nUnit,*) 'GhostFromCoarser'
    write(nUnit,*)
    write(nUnit,*) '-----------------------------------------------------------'
    call dump_FinerGhostsFromMe( intp, levelDesc, nUnit, tree )
    write(nUnit,*) '-----------------------------------------------------------'
    write(nUnit,*) '-----------------------------------------------------------'
    write(nUnit,*)
    write(nUnit,*) 'GhostFromCoarserBuffer'
    write(nUnit,*)
    write(nUnit,*) '-----------------------------------------------------------'
    call dump_FinerGhostsFromMeBuffer( intp, leveldesc, nUnit, tree )
    write(nUnit,*) '-----------------------------------------------------------'

    close( nUnit )

  end subroutine debug_dependencies
! ****************************************************************************** !

! ****************************************************************************** !
  !> check the dependencies from Finer
  !!
  subroutine dump_intpLists( minLevel, maxLevel, order, levelDesc, rank )
    ! ---------------------------------------------------------------------------
    !> global pdf information
    integer, intent(in) :: minLevel, maxLevel
    !> state properties
    type( tem_leveldesc_type ), intent(in) :: levelDesc(minLevel:maxLevel)
    integer, intent(in) :: order
    !> musubi mpi communicator environment
    integer, intent(in) :: rank
    ! ---------------------------------------------------------------------------
    integer :: nUnit, iLevel, iOrder
    character(len=17)  :: buffer
    ! ---------------------------------------------------------------------------
    nUnit = newunit()
    write(buffer,'(a4,i7.7)') 'dumpIntpLists.', rank
    open(file = trim(buffer), unit= nUnit, recl=1024)

    do iLevel = minLevel, maxLevel
      write(nUnit,*)
      write(nUnit,*) 'GhostFromFiner      '
      write(nUnit,*)
      write(nUnit,*) '---------------------------------------------------------'
      call dump_intpList( eType = eT_ghostFromFiner, &
        &    levelDesc = levelDesc(iLevel), &
        &    ind = levelDesc( iLevel)%intpFromFiner, nUnit = nUnit )
      write(nUnit,*) '---------------------------------------------------------'
      write(nUnit,*) '---------------------------------------------------------'
      write(nUnit,*)
      write(nUnit,*) 'GhostFromCoarser    '
      write(nUnit,*)
      do iOrder = 0, order
        write(nUnit,*) '--------------------------------------------------------'
        write(nUnit,*) 'order: ',iOrder, '--------------------------------------'
        write(nUnit,*) '--------------------------------------------------------'
        call dump_intpList( eType = eT_ghostFromCoarser, &
          &    levelDesc = levelDesc(iLevel), &
          &    ind = levelDesc(iLevel)%intpFromCoarser(iOrder), nUnit = nUnit )
        write(nUnit,*) '--------------------------------------------------------'
      end do
    end do
    close( nUnit )

  end subroutine dump_intpLists
! ****************************************************************************** !

! ****************************************************************************** !
  !> check the dependencies from Finer and write them out so we can compare
  !!
  subroutine dump_intpList( eType, levelDesc, ind, nUnit )
    ! ---------------------------------------------------------------------------
    !> state properties
    type( tem_levelDesc_type ), intent(in) :: levelDesc
    !> indirectio list
    type( grw_intArray_type ), intent(in)  :: ind
    !>
    integer, intent(in) :: nUnit
    !>
    integer, intent(in) :: eType
    ! ---------------------------------------------------------------------------
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: targetElem     ! element counter for indirection list
    integer(kind=long_k) :: tID
    character(len=pathLen) :: buffer
    ! ---------------------------------------------------------------------------

    buffer = ''
    do indElem = 1, ind%nVals
      iElem = ind%val( indElem )
      targetElem = iElem + levelDesc%offset( 1, eType )
      tID = levelDesc%total( targetElem )
      write(buffer, '(2a, i9)') trim(buffer), ' ', tID
      if( mod( indElem, 32 ) == 0 .or. indElem == ind%nVals ) then
        write(nUnit, *) trim(buffer)
        buffer = ''
      end if
    enddo

  end subroutine dump_intpList
! ****************************************************************************** !

! ****************************************************************************** !
  !> check the dependencies from Finer and write them out so we can compare
  !!
  subroutine dump_MyGhostsFromFiner( intp, levelDesc, nUnit, tree)
    ! ---------------------------------------------------------------------------
    !> state properties
    type( tem_leveldesc_type ), intent(in) :: leveldesc(:)
    !> interpolation method info
    type( mus_interpolation_type ), intent(in) :: intp
    !> global tree information
    type( treelmesh_type ), intent(in) :: tree
    !> unit to write to
    integer, intent(in) :: nUnit
    ! --------------------------------------------------------------------------
    integer :: ilevel         ! grid refinement level
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iElem          ! current target element (for outer loop)
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    integer(kind=long_k) :: tID(0:intp%fillMineFromFiner%nMaxSources )
    ! --------------------------------------------------------------------------

    do iLevel = tree%global%minLevel, tree%global%maxLevel - 1

      sourceLevel = iLevel + 1
      targetLevel = iLevel

      ! Treat all coarse target elements
      do iElem = 1, levelDesc( targetLevel )%elem%nElems( eT_ghostFromFiner )
        ! Read the target element treeId
        targetElem = iElem + levelDesc( targetlevel )%offset( 1, eT_ghostFromFiner )
        ! Find out how many fine source elements we have for interpolation.
        ! Usually 8, but can differ at corners, obstacles, boundaries...
        nSourceELems = levelDesc( ilevel )%depFromFiner( iElem )%elem%nVals
        tID = 0
        tID( 0 ) = levelDesc( targetLevel )%total( targetElem )

        write(nUnit, '(a,i9,3f10.5)') ' targetGhost  ', tID( 0 ),              &
          &                           tem_baryOfId( tree, tID(0) )

        ! Now loop over all fine source elements for this target:
        do iSourceElem = 1, nSourceElems

          ! Get the source element's treeId
          sourceElem = levelDesc( targetLevel )%depFromFiner( iElem )     &
            &                                       %elem%val( iSourceElem )
          tID( iSourceElem ) = levelDesc( sourceLevel )%total( sourceElem )

          call dump_elemDep( targetElem = tID(0),                              &
            &                sourceElem = tID( iSourceElem ),                  &
            &                nUnit = nUnit,                                    &
            &                tree = tree)
        end do  ! iSourceElem
      enddo
    enddo

  end subroutine dump_MyGhostsFromFiner
! ****************************************************************************** !

! ****************************************************************************** !
  !> check the dependencies from Coarser
  !!
  subroutine dump_FinerGhostsFromMe( intp, leveldesc, nUnit, tree)
    ! ---------------------------------------------------------------------------
    !> state properties
    type( tem_levelDesc_type ), intent(in) :: levelDesc(:)
    !> interpolation method info
    type( mus_interpolation_type ), intent(in) :: intp
    !> global tree information
    type( treelmesh_type ), intent(in) :: tree
    !> unit to write to
    integer, intent(in) :: nUnit
    ! ---------------------------------------------------------------------------
    integer :: ilevel         ! grid refinement level
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iElem          ! current target element (for outer loop)
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    integer(kind=long_k), allocatable :: tID(:)
    logical :: weights
    integer :: nMaxSources
    ! ---------------------------------------------------------------------------
    nMaxSources = maxval(intp%fillFinerFromMe(:)%nMaxSources)
    allocate(tID(0:nMaxSources))

    do iLevel = tree%global%minLevel, tree%global%maxLevel - 1

      sourceLevel = ilevel
      targetLevel = ilevel + 1

      ! Treat all coarse target elements
      do iElem = 1, levelDesc( targetLevel )%elem%nElems( eT_ghostFromCoarser )
        ! Read the target element treeId
        targetElem = iElem + levelDesc( targetlevel )%offset( 1, eT_ghostFromCoarser)
        ! Find out how many fine source elements we have for interpolation.
        ! Usually 8, but can differ at corners, obstacles, boundaries...
        nSourceELems = levelDesc( targetlevel )%depFromCoarser( iElem )%  &
          &                                                          elem%nVals
        tID = 0
        tID( 0 ) = levelDesc( targetLevel )%total( targetElem )
        write(nUnit, '(a,i9,3f10.5)') ' targetGhost  ', tID( 0 ),              &
          &                                        tem_baryOfId( tree, tID(0) )

        ! Now loop over all fine source elements for this target:
        do iSourceElem = 1, nSourceElems

          ! Get the source element's treeId
          sourceElem = levelDesc( targetLevel )%depFromCoarser( iElem )   &
            &                                       %elem%val( iSourceElem )
          tID( iSourceElem ) = levelDesc( sourceLevel )%total( sourceElem )

          if( allocated( levelDesc( targetLevel )                         &
              &       %depFromCoarser( iElem )%weight)) then
            weights = .true.
          else
            weights = .false.
          end if
          if( weights ) then
            call dump_elemDep( targetElem = tID(0),                            &
              &                sourceElem = tID( iSourceElem ),                &
              &                nUnit = nUnit,                                  &
              &                tree = tree,                                    &
              &                weight = levelDesc( targetLevel )          &
              &                   %depFromCoarser( iElem )%weight( iSourceElem))
          else
            call dump_elemDep( targetElem = tID(0),                            &
              &                sourceElem = tID( iSourceElem ),                &
              &                nUnit = nUnit,                                  &
              &                tree = tree )
          end if

        end do ! iSourceElem
      end do ! iElem
    end do ! iLevel

  end subroutine dump_FinerGhostsFromMe
! ****************************************************************************** !

! ****************************************************************************** !
  !> check the dependencies from Coarser
  !!
  subroutine dump_FinerGhostsFromMeBuffer( intp, levelDesc, nUnit, tree )
    ! ---------------------------------------------------------------------------
    !> state properties
    type( tem_leveldesc_type ), intent(in) :: leveldesc(:)
    !> interpolation method info
    type( mus_interpolation_type ), intent(in) :: intp
    !> global tree information
    type( treelmesh_type ), intent(in) :: tree
    !> unit to write to
    integer, intent(in) :: nUnit
    ! ---------------------------------------------------------------------------
    integer :: ilevel         ! grid refinement level
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iElem          ! current target element (for outer loop)
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    integer :: totalPos ! Position of the element in the total list
    integer(kind=long_k), allocatable :: tID(:)
    logical :: weights
    integer :: nMaxSources
    ! ---------------------------------------------------------------------------
    nMaxSources = maxval(intp%fillFinerFromMe(:)%nMaxSources)
    allocate(tID(0:nMaxSources))

    do iLevel = tree%global%minLevel, tree%global%maxLevel - 1
      sourceLevel = ilevel
      targetLevel = ilevel + 1

      ! Treat all coarse target elements
      do iElem = 1, levelDesc( targetLevel )%elem%nElems( eT_ghostFromCoarser )
        ! Read the target element treeId
        targetElem = iElem + levelDesc( targetlevel )                     &
          &                      %offset( 1, eT_ghostFromCoarser)
        ! Find out how many fine source elements we have for interpolation.
        ! Usually 8, but can differ at corners, obstacles, boundaries...
        nSourceELems = levelDesc( targetlevel )%depFromCoarser( iElem )%  &
          &                                                          elem%nVals
        tID = 0
        tID( 0 ) = levelDesc( targetLevel )%total( targetElem )
        write(nUnit, '(a,i9,3f10.5)') ' targetGhost  ', tID( 0 ),              &
          &                           tem_baryOfId( tree, tID(0) )

        ! Now loop over all fine source elements for this target:
        do iSourceElem = 1, nSourceElems

          ! Get the source element's treeId
          sourceElem = levelDesc( targetLevel )%depFromCoarser( iElem )   &
            &                                    %elemBuffer%val( iSourceElem )
          totalPos = levelDesc( targetLevel )%sourceFromCoarser%          &
            &                                                 val( sourceElem )
          tID( iSourceElem ) = levelDesc( sourceLevel )%total( totalPos )

          if( allocated( levelDesc( targetLevel )                         &
              &       %depFromCoarser( iElem )%weight)) then
              weights = .true.
          else
            weights = .false.
          end if
          if( weights ) then
            call dump_elemDep( targetElem = tID(0),                            &
              &                sourceElem = tID( iSourceElem ),                &
              &                nUnit = nUnit,                                  &
              &                tree = tree,                                    &
              &                weight = levelDesc( targetLevel )          &
              &                   %depFromCoarser( iElem )%weight( iSourceElem))
          else
            call dump_elemDep( targetElem = tID(0),                            &
              &                sourceElem = tID( iSourceElem ),                &
              &                nUnit = nUnit,                                  &
              &                tree = tree )
          end if

        end do  ! iSourceElem
      enddo
    enddo

  end subroutine dump_FinerGhostsFromMeBuffer
! ****************************************************************************** !

! ****************************************************************************** !
  !> dump dependencies for one element
  !!
  subroutine dump_elemDep( targetElem, sourceElem, nUnit, tree, weight )
    ! ---------------------------------------------------------------------------
    !>
    integer(kind=long_k), intent(in) :: targetElem
    !>
    integer(kind=long_k), intent(in) :: sourceElem
    !>
    integer, intent(in) :: nUnit
    !>
    type( treelmesh_type ), intent(in) :: tree
    !>
    real(kind=rk), optional :: weight
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: xTarget(3), xSource(3)
    real(kind=rk) :: elemSize
    character(len=pathLen) :: buffer
    ! ---------------------------------------------------------------------------
    buffer = ''
    elemSize = tem_elemSize( tree, targetElem )
    xTarget = tem_baryOfId( tree, targetElem )
    xSource = tem_baryOfId( tree, sourceElem )
    write(buffer,'(a,i9, 6f10.5)') '  ', sourceElem,                           &
      &          ((xSource(1:3) - xTarget(1:3))/elemSize), xSource(1:3)
    if( present( weight )) then
      write(buffer, '(a,a,f8.5)') trim( buffer), '    weight: ', weight
    end if
    write(nUnit, '(a)') trim(buffer)

  end subroutine dump_elemDep
! ****************************************************************************** !

end module mus_interpolate_tools_module
! ****************************************************************************** !
