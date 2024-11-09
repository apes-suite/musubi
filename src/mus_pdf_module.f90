! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2016, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2017 Jana Gericke <jana.gericke@uni-siegen.de>
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
!> Definitions for the main state and neighbor arrays
!!
module mus_pdf_module

  ! include treelm modules
  use env_module,              only: rk, globalMaxLevels, long_k
  use tem_construction_module, only: tem_levelDesc_type
  use tem_debug_module,        only: dbgUnit
  use tem_logging_module,      only: logUnit

  implicit none
  private

  public :: pdf_data_type
  public :: mus_calc_nElems
  public :: mus_pdf_allocate
  public :: mus_swap_Now_Next
  public :: allocate_momBuf

  !> This level-wise data type contains the PDF arrays which are fed into the
  !! kernel.
  !! The solver updates the state vector and finds the position of the neighbor
  !! vectors by looking it up in the neigh array
  type pdf_data_type
    ! Per level element numbers
    !> number of local fluid elements
    integer :: nElems_fluid
    !> number of halo elements (from remote processes)
    integer :: nElems_halo
    !> number of ghost elements (from other levels)
    integer :: nElems_ghostFromCoarser
    !> number of ghost elements (from other levels)
    integer :: nElems_ghostFromFiner

    !> number of ghost elements (from other levels)
    integer :: nElems_ghost
    !> number of local elements (fluid+halos+ghost)
    integer :: nElems_local
    !> fluid elements + ghostFromCoarser elements (elements for solver).
    !! PDF and auxFied are interpolated for ghostFromFiner and there is no
    !! need to do collision on these elements.
    integer :: nElems_solve
    !> fluid elements + ghostFromCoarse + ghostFromFiner
    integer :: nElems_computed

    !> number of elements padded to 4
    integer :: nSize

    !> which buffer to use for current time step
    integer :: nNow  = 1
    !> which buffer to use for next time step
    integer :: nNext = 2

    !> Connectivity array
    !! Points to where to send respective pdfs
    !! Access in a neigh way
    !! Size: QQ * nSize
    !! allocated in routine: mus_pdf_allocate
    integer,allocatable,dimension(:) :: neigh
    !dir$ attributes align : 32 :: neigh

    !> containing state vector values of elements which have a boundary
    !! It always uses AOS data layout
    !! allocated in routine: mus_pdf_allocate
    !! filled for each iteration in routine: fill_bcBuffer
    real(kind=rk),allocatable, dimension(:) :: bcBuffer

    !> Buffer storing the moments of all source from Coarser
    real(kind=rk), allocatable :: momBuf(:,:)

  end type pdf_data_type


contains


! **************************************************************************** !
  !> Compute nElems for different types
  subroutine mus_calc_nElems(me, nFluids, nGhostFromCoarser, nGhostFromFiner, &
    &                        nHalos)
    ! --------------------------------------------------------------------------
    type( pdf_data_type ) :: me
    integer, intent(in) :: nFluids, nGhostFromCoarser, nGhostFromFiner, nHalos
    ! --------------------------------------------------------------------------
    integer :: remainder

    me%nElems_fluid            = nFluids
    me%nElems_ghostFromCoarser = nGhostFromCoarser
    me%nElems_ghostFromFiner   = nGhostFromFiner
    me%nElems_halo             = nHalos

    me%nElems_local = nFluids + nGhostFromCoarser + nGhostFromFiner + nHalos
    me%nElems_ghost = nGhostFromCoarser + nGhostFromFiner
    me%nElems_solve = nFluids + nGhostFromCoarser
    me%nElems_computed = nFluids + nGhostFromCoarser + nGhostFromFiner

    remainder = mod( me%nElems_local, 4 )
    me%nSize  = me%nElems_local + mod(4 - remainder, 4)

    ! print out number of elements of all kinds
    write(dbgUnit(1),"(A   )") ''
    write(dbgUnit(1),"(A   )") ' Number of local elements: '
    write(dbgUnit(1),"(A,I0)") '   fluid: ', me%nElems_fluid
    write(dbgUnit(1),"(A,I0)") '   ghostFromCoarser: ', &
      &                        me%nElems_ghostFromCoarser
    write(dbgUnit(1),"(A,I0)") '   ghostFromFiner:   ', me%nElems_ghostFromFiner
    write(dbgUnit(1),"(A,I0)") '   halo: ', me%nElems_halo
    write(dbgUnit(1),"(A,I0)") '   size: ', me%nSize
    write(dbgUnit(1),"(A   )") ''
    write(logUnit(1),"(A,I0)") ' nSize of state array: ', me%nSize

  end subroutine mus_calc_nElems
! **************************************************************************** !


! **************************************************************************** !
  subroutine mus_pdf_allocate(me, nScalars, QQ, nElems_bcBuffer, isPDF )
    ! --------------------------------------------------------------------------
    type( pdf_data_type ), intent(inout) :: me
    integer, intent(in) :: nScalars, QQ, nElems_bcBuffer
    logical, intent(in) :: isPDF
    ! --------------------------------------------------------------------------

    if ( isPDF ) then
      allocate( me%neigh( me%nSize * QQ ) )
      allocate( me%bcBuffer( nElems_bcBuffer * nScalars ) )

      if ( nElems_bcBuffer /= 0 ) me%bcBuffer = -1.0_rk
    end if

  end subroutine mus_pdf_allocate
! **************************************************************************** !


! **************************************************************************** !
  subroutine mus_swap_Now_Next( me )
    ! --------------------------------------------------------------------------
    type( pdf_data_type ), intent(inout) :: me
    ! --------------------------------------------------------------------------

    me%nNow  = mod( me%nNow ,2 )+1
    me%nNext = mod( me%nNext,2 )+1

  end subroutine mus_swap_Now_Next
! **************************************************************************** !

! **************************************************************************** !
  subroutine allocate_momBuf( me, nVals )
    type( pdf_data_type ), intent(inout) :: me
    integer, intent(in) :: nVals

    if ( nVals > 0 ) then
      allocate( me%momBuf(10,nVals) )
    else
      allocate( me%momBuf(0,0) )
    end if

  end subroutine allocate_momBuf
! **************************************************************************** !

end module mus_pdf_module
! **************************************************************************** !


! **************************************************************************** !
!>\page datastructures Data structures
!! The Octree data structure is mapped to a one-dimensional array in order to
!! have an efficient data structure on which the solver can act on in a
!! performant way. An efficient representation of the elements and their
!! neighbor relations is chosen.
!! The fluid elements are mapped to a one-dimensional array and the neighboring
!! relations are introduced by an additional connectivity array. The access of a
!! neighbor element is performed by looking up the correct position of an
!! element's link neighbor in the connectivity array, thus constituting an
!! indirect access.
!!
!! \image html  statevector_neighborlist.png
!!
!! The different dependencies of each link in an element require a thorough
!! treatment, when data is exchanged at domain boundaries. Only the links, which
!! point outside the domain have to be sent to neighbor partitions, and the
!! links pointing inwards have to be filled with valid values from these.
