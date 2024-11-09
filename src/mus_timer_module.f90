! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2021 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> This module contains data types and routines used to measure time spend on
!! musubi routines.
!! author: Kannan Masilamani
!!
module mus_timer_module
  use env_module,         only: rk
  use tem_timer_module,   only: tem_resetTimer, tem_addTimer,&
    &                           tem_getTimerVal, tem_getNTimers

  implicit none

  private

  public :: mus_timerHandles
  public :: mus_timer_handle_type
  public :: mus_init_mainTimer
  public :: mus_init_levelTimer
  public :: mus_init_bcTimer
  public :: mus_reset_levelTimer
  public :: mus_reset_mainTimer
  public :: mus_reset_bcTimer
  public :: mus_get_timerHandles
  public :: mus_set_timerHandles

  public :: get_mainLoopTime

  public :: get_computeTime
  public :: get_computeRatio

  public :: get_auxTime
  public :: get_relaxTime

  public :: get_intpFromCoarserTime
  public :: get_intpFromCoarserRatio
  public :: get_intpFromFinerTime
  public :: get_intpFromFinerRatio
  public :: get_intpRatio

  public :: get_communicateTime
  public :: get_communicateRatio

  public :: get_bcBufferTime
  public :: get_bcBufferRatio

  public :: get_boundaryTime
  public :: get_boundaryRatio

  public :: get_stageTime
  public :: get_stageRatio

  integer, parameter, public :: nStages = 12

  interface get_boundaryTime
      module procedure get_boundaryTime_total
      module procedure get_boundaryTime_byID
  end interface get_boundaryTime

  interface get_computeTime
      module procedure get_computeTime_total
      module procedure get_computeTime_atLevel
  end interface get_computeTime

  interface get_computeRatio
      module procedure get_computeRatio_total
      module procedure get_computeRatio_atLevel
  end interface get_computeRatio

  interface get_auxTime
      module procedure get_auxTime_total
      module procedure get_auxTime_atLevel
  end interface get_auxTime

  interface get_relaxTime
      module procedure get_relaxTime_total
      module procedure get_relaxTime_atLevel
  end interface get_relaxTime

  interface get_intpFromCoarserTime
      module procedure get_intpFromCoarserTime_total
      module procedure get_intpFromCoarserTime_atLevel
  end interface get_intpFromCoarserTime

  interface get_intpFromCoarserRatio
      module procedure get_intpFromCoarserRatio_total
      module procedure get_intpFromCoarserRatio_atLevel
  end interface get_intpFromCoarserRatio

  interface get_intpFromFinerTime
      module procedure get_intpFromFinerTime_total
      module procedure get_intpFromFinerTime_atLevel
  end interface get_intpFromFinerTime

  interface get_intpFromFinerRatio
      module procedure get_intpFromFinerRatio_total
      module procedure get_intpFromFinerRatio_atLevel
  end interface get_intpFromFinerRatio

  interface get_bcBufferTime
      module procedure get_bcBufferTime_total
      module procedure get_bcBufferTime_atLevel
  end interface get_bcBufferTime

  interface get_bcBufferRatio
      module procedure get_bcBufferRatio_total
      module procedure get_bcBufferRatio_atLevel
  end interface get_bcBufferRatio

  !> Musubi timer type --------------------------------------------------
  type mus_timer_handle_type

    !> handle for the complete mainloop
    integer :: mainloop

    !> handle for loading / creating the mesh and config
    integer :: loadMesh

    !> handle for initialising the levelDescriptor
    integer :: initLvlD

    !> handle for writing restart
    integer :: wRestart

    !> handle for the dyn_loadBal routine
    integer :: balance

    !> handle for source terms
    integer :: source

    ! GGS: erase the following
    !integer :: kappa_abg, inv_kappa_abg, inner_loop

    !> First main handle position in treelm timer object
    integer :: first = 0

    !> Last main handle position in treelm timer object
    integer :: last = -1

    ! level-wise timers:

    !> handle for advection relaxation
    integer, allocatable :: compute(:)

    !> handle for auxfield calculation
    integer, allocatable :: aux(:)

    !> handle for relax parameter update
    integer, allocatable :: relax(:)

    !> handle for communicate
    integer, allocatable :: comm(:)

    !> handle for interpolation and communicate
    integer, allocatable :: intpFromCoarser(:)
    integer, allocatable :: intpFromFiner(:)
    integer, allocatable :: commFromCoarser(:)
    integer, allocatable :: commFromFiner(:)

    !> handle for setboundary
    integer, allocatable :: setBnd(:)
    integer, allocatable :: bcBuffer(:)

    !> handle for immersed boundary method
    integer, allocatable :: doIBM(:)

    !> Stage timers for multi level recursive algorithm
    integer :: stage(nStages)

    !> min. level in mesh
    integer :: minLevel
    !> max. level in mesh
    integer :: maxLevel

    !> number of BCs
    integer :: nBCs

  end type mus_timer_handle_type

  !> Musubi timer type --------------------------------------------------
  type( mus_timer_handle_type ), save :: mus_timerHandles


contains


  ! ------------------------------------------------------------------------ !
  !> Timers initialization routine for whatever
  !!
  subroutine mus_init_mainTimer()
    ! -------------------------------------------------------------------- !
    character(len=2) :: buffer
    integer :: ii
    ! -------------------------------------------------------------------- !

    ! add timer handles to measure wall clock time of different routines
    call tem_addTimer( timerHandle = mus_timerHandles%mainLoop, &
      &                timerName   = 'MainLoop'         )

    ! Set position of the fisrt handle
    mus_timerHandles%first = tem_getNTimers()

    call tem_addTimer( timerHandle = mus_timerHandles%loadMesh,&
      &                timerName   = 'LoadMesh'        )

    call tem_addTimer( timerHandle = mus_timerHandles%initLvlD,&
      &                timerName   = 'InitLvlD' )

    call tem_addTimer( timerHandle = mus_timerHandles%wRestart, &
      &                timerName   = 'wRestart' )

    call tem_addTimer( timerHandle = mus_timerHandles%balance, &
      &                timerName   = 'Balance' )

    call tem_addTimer( timerHandle = mus_timerHandles%source, &
      &                timerName   = 'Source' )

    ! GGS: erase the following 3
    !call tem_addTimer( timerHandle = mus_timerHandles%kappa_abg, &
    !  &                timerName   = 'Chimera_ops' )
    !
    !call tem_addTimer( timerHandle = mus_timerHandles%inv_kappa_abg, &
    !  &                timerName   = 'Backward_Chimera_ops' )
    !
    !call tem_addTimer( timerHandle = mus_timerHandles%inner_loop, &
    !  &                timerName   = 'Inner_Loop' )

    ! Set position of the last handle
    mus_timerHandles%last = tem_getNTimers()

    do ii = 1, nStages
      write(buffer, "(I2.2)") ii
      call tem_addTimer( timerHandle = mus_timerHandles%stage(ii), &
        &                timerName   = 'Stage'//(buffer)           )
    end do

  end subroutine mus_init_mainTimer
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_init_levelTimer(minLevel, maxLevel)
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: minLevel, maxLevel
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    character(len=3) :: level_string
    ! -------------------------------------------------------------------- !

    mus_timerHandles%minLevel = minLevel
    mus_timerHandles%maxLevel = maxLevel
    ! Allocate level-wise handles
    allocate( mus_timerHandles%compute( minLevel:maxLevel ) )
    allocate( mus_timerHandles%aux( minLevel:maxLevel ) )
    allocate( mus_timerHandles%relax( minLevel:maxLevel ) )
    allocate( mus_timerHandles%comm( minLevel:maxLevel ) )
    allocate( mus_timerHandles%intpFromCoarser( minLevel:maxLevel ) )
    allocate( mus_timerHandles%intpFromFiner( minLevel:maxLevel ) )
    allocate( mus_timerHandles%commFromCoarser( minLevel:maxLevel ) )
    allocate( mus_timerHandles%commFromFiner( minLevel:maxLevel ) )
    allocate( mus_timerHandles%doIBM( minLevel:maxLevel ) )
    allocate( mus_timerHandles%bcBuffer( minLevel:maxLevel ) )

    do iLevel = mus_timerHandles%minLevel, mus_timerHandles%maxLevel

      ! level prependix
      write( level_string, "(A, I2.2)") "L", iLevel

      call tem_addTimer( timerHandle = mus_timerHandles%compute(iLevel), &
        &                timerName   = level_string//'_compute' )

      call tem_addTimer( timerHandle = mus_timerHandles%aux(iLevel), &
        &                timerName   = level_string//'_aux' )

      call tem_addTimer( timerHandle = mus_timerHandles%relax(iLevel), &
        &                timerName   = level_string//'_relax' )

      call tem_addTimer( timerHandle = mus_timerHandles%comm(iLevel), &
        &                timerName   = level_string//'_comm' )

      call tem_addTimer( timerHandle = mus_timerHandles                 &
        &                              %intpFromCoarser(iLevel),        &
        &                timerName   = level_string//'_intpFromCoarser' )

      call tem_addTimer( timerHandle = mus_timerHandles%intpFromFiner(iLevel), &
        &                timerName   = level_string//'_intpFromFiner' )

      call tem_addTimer( timerHandle = mus_timerHandles                 &
        &                              %commFromCoarser(iLevel),        &
        &                timerName   = level_string//'_commFromCoarser' )

      call tem_addTimer( timerHandle = mus_timerHandles%commFromFiner(iLevel), &
        &                timerName   = level_string//'_commFromFiner' )

      call tem_addTimer( timerHandle = mus_timerHandles%doIBM(iLevel), &
        &                timerName   = level_string//'_doIBM' )

      call tem_addTimer( timerHandle = mus_timerHandles%bcBuffer(iLevel), &
        &                timerName   = level_string//'_bcBuffer' )

    end do

  end subroutine mus_init_levelTimer
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_reset_levelTimer()
    ! -------------------------------------------------------------------- !
    integer :: iLevel, ii
    ! -------------------------------------------------------------------- !

    do iLevel = mus_timerHandles%minLevel, mus_timerHandles%maxLevel
      call tem_resetTimer( timerhandle = mus_timerHandles%compute(iLevel) )
      call tem_resetTimer( timerhandle = mus_timerHandles%aux(iLevel) )
      call tem_resetTimer( timerhandle = mus_timerHandles%relax(iLevel) )
      call tem_resetTimer( timerhandle = mus_timerHandles%comm(iLevel) )
      call tem_resetTimer(                                            &
        &      timerhandle = mus_timerHandles%intpFromCoarser(iLevel) )
      call tem_resetTimer(                                          &
        &      timerhandle = mus_timerHandles%intpFromFiner(iLevel) )
      call tem_resetTimer(                                            &
        &      timerhandle = mus_timerHandles%commFromCoarser(iLevel) )
      call tem_resetTimer(                                          &
        &      timerhandle = mus_timerHandles%commFromFiner(iLevel) )
      ! call tem_resetTimer( timerhandle = mus_timerHandles%setBnd(iLevel) )
      call tem_resetTimer( timerhandle = mus_timerHandles%doIBM(iLevel)  )
      call tem_resetTimer( timerhandle = mus_timerHandles%bcBuffer(iLevel)  )
    end do

    do ii = 1, nStages
      call tem_resetTimer( timerhandle = mus_timerHandles%stage(ii) )
    end do

  end subroutine mus_reset_levelTimer
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_reset_mainTimer()
    ! -------------------------------------------------------------------- !

    call tem_resetTimer( timerhandle = mus_timerHandles%mainLoop )
    call tem_resetTimer( timerhandle = mus_timerHandles%loadMesh )
    call tem_resetTimer( timerhandle = mus_timerHandles%initLvlD )
    call tem_resetTimer( timerhandle = mus_timerHandles%wRestart )
    call tem_resetTimer( timerhandle = mus_timerHandles%balance  )
    call tem_resetTimer( timerhandle = mus_timerHandles%source   )

  end subroutine mus_reset_mainTimer
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function returns musubi modular variable mus_timerHandles to apesmate
  !! and deallocate mus_timerHandles level timers.
  function mus_get_timerHandles() result(timerHandles)
    ! -------------------------------------------------------------------- !
    type(mus_timer_handle_type) :: timerHandles
    ! -------------------------------------------------------------------- !
    timerHandles = mus_timerHandles

    deallocate( mus_timerHandles%compute )
    deallocate( mus_timerHandles%aux )
    deallocate( mus_timerHandles%relax )
    deallocate( mus_timerHandles%comm )
    deallocate( mus_timerHandles%intpFromCoarser )
    deallocate( mus_timerHandles%intpFromFiner )
    deallocate( mus_timerHandles%commFromCoarser )
    deallocate( mus_timerHandles%commFromFiner )
    deallocate( mus_timerHandles%setBnd )
    deallocate( mus_timerHandles%doIBM )
    deallocate( mus_timerHandles%bcBuffer )

  end function mus_get_timerHandles
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine sets mus_timerHandles passed by apesmate
  subroutine mus_set_timerHandles(timerHandles)
    ! -------------------------------------------------------------------- !
    type(mus_timer_handle_type), intent(in) :: timerHandles
    ! -------------------------------------------------------------------- !
    mus_timerHandles = timerHandles
  end subroutine mus_set_timerHandles
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_computeTime_total() result ( time )
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = 0.0_rk
    do iLevel = mus_timerHandles%minLevel, mus_timerHandles%maxLevel
      time = time + get_computeTime( iLevel )
    end do

  end function get_computeTime_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_computeTime_atLevel( level ) result ( time )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: level
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = tem_getTimerVal( timerhandle = mus_timerHandles%compute(level) )

  end function get_computeTime_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_computeRatio_total() result ( ratio )
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_computeTime() / get_mainLoopTime() * 100._rk

  end function get_computeRatio_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_computeRatio_atLevel( level ) result ( ratio )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: level
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_computeTime( level ) / get_mainLoopTime() * 100._rk

  end function get_computeRatio_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_auxTime_total() result ( time )
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = 0.0_rk
    do iLevel = mus_timerHandles%minLevel, mus_timerHandles%maxLevel
      time = time + get_auxTime( iLevel )
    end do

  end function get_auxTime_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_auxTime_atLevel( level ) result ( time )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: level
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = tem_getTimerVal( timerhandle = mus_timerHandles%aux(level) )

  end function get_auxTime_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_relaxTime_total() result ( time )
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = 0.0_rk
    do iLevel = mus_timerHandles%minLevel, mus_timerHandles%maxLevel
      time = time + get_relaxTime( iLevel )
    end do

  end function get_relaxTime_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_relaxTime_atLevel( level ) result ( time )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: level
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = tem_getTimerVal( timerhandle = mus_timerHandles%relax(level) )

  end function get_relaxTime_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpFromCoarserTime_total() result ( time )
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = 0.0_rk
    do iLevel = mus_timerHandles%minLevel, mus_timerHandles%maxLevel
      time = time + get_intpFromCoarserTime( iLevel )
    end do

  end function get_intpFromCoarserTime_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpFromCoarserTime_atLevel( level ) result ( time )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: level
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = tem_getTimerVal(                                          &
      &        timerhandle = mus_timerHandles%intpFromCoarser(level) )

  end function get_intpFromCoarserTime_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpFromCoarserRatio_total() result ( ratio )
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_intpFromCoarserTime() / get_mainLoopTime() * 100._rk

  end function get_intpFromCoarserratio_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpFromCoarserRatio_atLevel( level ) result ( ratio )
    ! -------------------------------------------------------------------- !
    integer :: level
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_intpFromCoarserTime( level ) / get_mainLoopTime() * 100._rk

  end function get_intpFromCoarserratio_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpFromFinerTime_total() result ( time )
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = 0.0_rk
    do iLevel = mus_timerHandles%minLevel, mus_timerHandles%maxLevel
      time = time + tem_getTimerVal(                                         &
        &               timerhandle = mus_timerHandles%intpFromFiner(iLevel) )
    end do

  end function get_intpFromFinerTime_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpFromFinerTime_atLevel( level ) result ( time )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: level
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = tem_getTimerVal(timerhandle = mus_timerHandles%intpFromFiner(level))

  end function get_intpFromFinerTime_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpFromFinerRatio_total() result ( ratio )
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_intpFromFinerTime() / get_mainLoopTime() * 100._rk

  end function get_intpFromFinerRatio_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpFromFinerRatio_atLevel( level ) result ( ratio )
    ! -------------------------------------------------------------------- !
    integer :: level
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_intpFromFinerTime( level ) / get_mainLoopTime() * 100._rk

  end function get_intpFromFinerRatio_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_communicateTime() result ( time )
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = 0.0_rk
    do iLevel = mus_timerHandles%minLevel, mus_timerHandles%maxLevel
      time = time + tem_getTimerVal(                                &
        &               timerhandle = mus_timerHandles%comm(iLevel) )
      time = time + tem_getTimerVal(                                         &
        &               timerhandle = mus_timerHandles%commFromFiner(iLevel) )
      time = time + tem_getTimerVal(                                           &
        &               timerhandle = mus_timerHandles%commFromCoarser(iLevel) )
    end do

  end function get_communicateTime
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_communicateRatio() result ( ratio )
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_communicateTime() / get_mainLoopTime() * 100._rk

  end function get_communicateRatio
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_boundaryTime_byID( bcID ) result ( time )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: bcID
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    if ( bcID > mus_timerHandles%nBCs ) then
      time = 0._rk
    else
      time = tem_getTimerVal( timerhandle = mus_timerHandles%setBnd(bcID) )
    end if

  end function get_boundaryTime_byID
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_boundaryTime_total() result ( time )
    ! -------------------------------------------------------------------- !
    integer :: ii
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = 0.0_rk
    do ii = 1, mus_timerHandles%nBCs
      time = time + get_boundaryTime( ii )
    end do

  end function get_boundaryTime_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_boundaryRatio() result ( ratio )
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_boundaryTime() / get_mainLoopTime() * 100._rk

  end function get_boundaryRatio
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpTime() result ( time )
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = get_intpFromCoarserTime() + get_intpFromFinerTime()

  end function get_intpTime
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_intpRatio() result ( ratio )
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_intpTime() / get_mainLoopTime() * 100.0_rk

  end function get_intpRatio
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_mainLoopTime() result ( time )
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = tem_getTimerVal( timerhandle = mus_timerHandles%mainLoop )

  end function get_mainLoopTime
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_stageTime(ii) result ( time )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: ii
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    if ( ii <= nStages ) then
      time = tem_getTimerVal( timerhandle = mus_timerHandles%stage(ii) )
    else
      time = - get_mainLoopTime()
    end if

  end function get_stageTime
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_stageRatio( ii ) result ( ratio )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: ii
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_stageTime(ii) / get_mainLoopTime() * 100._rk

  end function get_stageRatio
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_init_bcTimer( nBCs )
    integer, intent(in) :: nBCs
    character(len=4) :: bc_string
    integer :: ii

    mus_timerHandles%nBCs = nBCs
    allocate( mus_timerHandles%setBnd( nBCs ) )

    do ii = 1, nBCs
      write( bc_string, "(A, I2.2)") "BC", ii
      call tem_addTimer( timerHandle = mus_timerHandles%setBnd(ii), &
        &                timerName   = bc_string//'_setBnd' )
    end do

  end subroutine mus_init_bcTimer
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_reset_bcTimer( )
    integer :: ii

    do ii = 1, mus_timerHandles%nBCs
      call tem_resetTimer( timerhandle = mus_timerHandles%setBnd(ii) )
    end do
  end subroutine mus_reset_bcTimer
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_bcBufferTime_atLevel( level ) result ( time )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: level
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = tem_getTimerVal( timerhandle = mus_timerHandles%bcBuffer(level) )

  end function get_bcBufferTime_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

  ! ------------------------------------------------------------------------ !
  function get_bcBufferTime_total() result ( time )
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    real(kind=rk) :: time
    ! -------------------------------------------------------------------- !

    time = 0.0_rk
    do iLevel = mus_timerHandles%minLevel, mus_timerHandles%maxLevel
      time = time + get_bcBufferTime( iLevel )
    end do

  end function get_bcBufferTime_total
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_bcBufferRatio_atLevel( level ) result ( ratio )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: level
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_bcBufferTime( level ) / get_mainLoopTime() * 100._rk

  end function get_bcBufferRatio_atLevel
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function get_bcBufferRatio_total() result ( ratio )
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: ratio
    ! -------------------------------------------------------------------- !

    ratio = get_bcBufferTime() / get_mainLoopTime() * 100._rk

  end function get_bcBufferRatio_total
  ! ------------------------------------------------------------------------ !

end module mus_timer_module
! *************************************************************************** !
