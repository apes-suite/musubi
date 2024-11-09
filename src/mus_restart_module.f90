! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011-2012, 2016, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
! **************************************************************************** !
!> author: Simon Zimny
!! This module provides the MUSUBI subroutines needed for the restart
!! functionality.
!!
!! This module makes the restart IO functionality available for MUSUBI. When
!! using the write functionality the workflow is as follows
!! - open the file and write the 'normal header'
!! - run over the chunks and
!!   + serialize the data
!!   + write them to file
!! - close the file and write the 'last header'
!! For reading the restart files the workflow is as follows:
!! - open the file
!! - read the solver specific character
!! - loop over the chunks and
!!   + read the data from disc
!!   + unserialize the data and store it in the state vector
!! - close the file
!!
!! Further information on the treelm functions and the usage can be found in
!! the [[tem_restart_module]] and inside the Restart Usage section itself.
!!
module mus_restart_module

  ! include treelm modules
  use mpi
  use env_module,            only: rk, io_buffer_size
  use treelmesh_module,      only: treelmesh_type
  use tem_time_module,       only: tem_time_type
  use tem_timer_module,      only: tem_timer_type, tem_startTimer, tem_stopTimer
  use tem_restart_module,    only: tem_restart_type, tem_restart_openRead,     &
    &                              tem_restart_closeRead, tem_restart_readData,&
    &                              tem_restart_openWrite,                      &
    &                              tem_restart_closeWrite,                     &
    &                              tem_restart_writeData
  use tem_debug_module,      only: dbgUnit, main_debug
  use tem_logging_module,    only: logUnit

  ! include musubi modules
  use mus_scheme_type_module, only: mus_scheme_type
  use mus_buffer_module,      only: mus_pdf_unserialize, mus_pdf_serialize

  implicit none
  private

  public :: mus_writeRestart
  public :: mus_readRestart


contains


! **************************************************************************** !
  !> Write the serialized buffer assembled in mus_serializeData to disk
  !!
  subroutine mus_writeRestart( levelPointer, restart, scheme, tree, timing, &
    &                          timerHandle, suffix)
    ! --------------------------------------------------------------------------
    !> global pdf info
    integer, intent(in) :: levelPointer(:)
    !> restart information
    type(tem_restart_type), intent(inout)  :: restart
    !> array of schemes including the data to be serialized and dumped
    type( mus_scheme_type ), intent(inout) :: scheme
    !> mesh, provided in treelm format
    type(treelmesh_type), intent(in)       :: tree
    !> current simulation time information
    type(tem_time_type), intent(inout)     :: timing
    !> Timer handle
    integer, intent(in) :: timerHandle
    !> optional suffix (if present NO timestamp will be added!!!!)
    character(len=*), optional, intent(in) :: suffix
    ! --------------------------------------------------------------------------
    ! array for chunkwise memory bounded output
    real(kind=rk), allocatable :: buffer(:)
    ! counter variables
    integer :: elemOff
    ! --------------------------------------------------------------------------
    allocate( buffer( io_buffer_size ))

    ! open the file, write header and prepare buffering
    call tem_restart_openWrite( me     = restart,       &
      &                         tree   = tree,          &
      &                         timing = timing,        &
      &                         varSys = scheme%varSys, &
      &                         suffix = suffix         )

    ! debug output
    if (main_debug%debugRestart) then
      write(dbgUnit(1),*) 'DEBUG: WRITE RESTART'
      write(dbgUnit(1),*) 'minLevel ', tree%global%minLevel, &
        & ' PDF state array size ',                          &
        & size( scheme%state( tree%global%minLevel )%val, 1)
      write(dbgUnit(1),*) 'maxLevel ', tree%global%maxLevel, &
        & ' PDF state array size ', &
        & size( scheme%state( tree%global%maxLevel )%val, 1)
    end if

    elemOff = 0
    ! loop over the chunks
    do while ( elemOff < tree%nElems )

      ! set the number of elements that are on the stack
      restart%nChunkElems = min( restart%write_file%chunkSize, &
        &                        tree%nElems-elemOff           )

      ! prepare the data
      call mus_pdf_serialize(                                     &
        &    scheme,                                              &
        &    tree%treeID(elemOff+1:elemOff+restart%nChunkElems),  &
        &    levelPointer(elemOff+1:elemOff+restart%nChunkElems), &
        &    restart%varMap,                                      &
        &    restart%nChunkElems,                                 &
        &    buffer,                                              &
        &    tree%global%minLevel,                                &
        &    tree%global%maxLevel                                 )

      call tem_startTimer( timerHandle = timerHandle )

      ! write the data to file
      call tem_restart_writeData( restart, buffer )

      call tem_stopTimer( timerHandle = timerHandle )

      ! set elements offset for next chunk
      elemOff = elemOff + restart%nChunkElems
    end do

    ! close the file and write last header
    call tem_restart_closeWrite( me     = restart,      &
      &                          timing = timing,       &
      &                          tree   = tree,         &
      &                          varSys = scheme%varSys )

    deallocate( buffer )

    ! debug output
    if( main_debug%debugRestart )then
      write(dbgUnit(1),*) 'END DEBUG: WRITE RESTART'
    end if

  end subroutine mus_writeRestart
! **************************************************************************** !


! **************************************************************************** !
  !> Read the serialized restart file into the state vectors
  !!
  subroutine mus_readRestart( levelPointer, restart, scheme, tree )
    ! -------------------------------------------------------------------------
    !> restart information
    type(tem_restart_type), intent(inout)  :: restart
    !> mesh, provided in treelm format
    type(treelmesh_type), intent(in) :: tree
    !> Level pointer, from tree mesh to level descriptor
    integer, intent(in) :: levelPointer(tree%nElems)
    !> array of schemes including the data to be serialized and dumped
    type( mus_scheme_type ), intent(inout) :: scheme
    ! -------------------------------------------------------------------------
    ! local variables
    real(kind=rk), allocatable :: buffer(:)
    ! integer :: iChunk
    integer :: elemOff ! Offset in the overall number of elements
    ! -------------------------------------------------------------------------
    write(logUnit(1),*)'Reading restart...'

    allocate( buffer( io_buffer_size ))

    if( main_debug%debugRestart ) then
      write(dbgUnit(1),*) 'DEBUG: READ RESTART'
    end if

    ! open the file, read header and prepare buffering
    call tem_restart_openRead( me = restart )

    elemOff = 0
    ! loop over the chunks
    do while ( elemOff < tree%nElems )
      ! set the number of elements that are on the stack
      restart%nChunkElems = min( restart%read_file%chunkSize, &
        &                        tree%nElems-elemOff          )

      ! fill buffer from restart file
      call tem_restart_readData( restart, buffer )

      ! transfer buffer to level-wise PDF array
      call mus_pdf_unserialize(                                   &
        &    scheme,                                              &
        &    tree%treeID(elemOff+1:elemOff+restart%nChunkElems),  &
        &    levelPointer(elemOff+1:elemOff+restart%nChunkElems), &
        &    restart%varMap,                                      &
        &    restart%nChunkElems,                                 &
        &    buffer,                                              &
        &    tree%global%minLevel,                                &
        &    tree%global%maxLevel                                 )
      ! Set treeID offset for next buffer
      elemOff = elemOff + restart%nChunkElems
    end do

    ! close the file
    call tem_restart_closeRead( me = restart )

    deallocate( buffer )

    ! debug output
    if( main_debug%debugRestart )then
      write(dbgUnit(1),*) 'END DEBUG: READ RESTART'
    end if
    write(logUnit(1),*) 'Done reading restart.'

  end subroutine mus_readRestart
! **************************************************************************** !

end module mus_restart_module
! **************************************************************************** !
