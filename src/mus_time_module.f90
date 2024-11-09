! Copyright (c) 2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2020 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> This module contains time definition needed for musubi
module mus_time_module

  ! include treelm modules
  use env_module,             only: rk
  use tem_time_module,        only: tem_time_type, tem_time_dump
  use tem_timeControl_module, only: tem_timeControl_type, tem_timeControl_dump,&
    &                               tem_timeControl_update,                    &
    &                               tem_timeControl_globalTriggered
  use tem_logging_module,     only: logUnit

  implicit none
  private

  public :: mus_time_homogenize
  public :: mus_timeControl_homogenize
  public :: mus_time_modulo


contains


! **************************************************************************** !
  !> Convert itime from restart to real time
  !!
  subroutine mus_time_homogenize( me, dt, readRestart )
    ! -------------------------------------------------------------------------
    type(tem_time_type), intent(inout) :: me
    real(kind=rk), intent(in) :: dt
    logical, intent(in) :: readRestart
    ! -------------------------------------------------------------------------

    ! if readRestart is active, set sim real time from iTime
    if ( readRestart ) then
      me%sim = real(me%iter, kind=rk)*dt
    end if

    ! debug time output
    write(logUnit(1),*) 'Simulation time homogenize:'
    call tem_time_dump( me, logUnit(1) )

  end subroutine mus_time_homogenize
! **************************************************************************** !

! **************************************************************************** !
  !> Converts sim time to iter and vice versa depends on which one is defined
  !! in the configuration file
  subroutine mus_timeControl_homogenize( me, dt, reqInt )
    ! -------------------------------------------------------------------------
    !> simulation time control
    type(tem_timeControl_type), intent(inout) :: me
    !> dt of maxlevel or smallest dt
    real(kind=rk), intent(in) :: dt
    !> Required interval, in which the update MUST occur.
    !! This is required for the musubi multilevel, where the time step should
    !! only be determined active, when the end of the largest cycle is reached.
    integer, intent(in) :: reqInt
    ! -------------------------------------------------------------------------

    ! if sim time is defined set itime
    if ( me%min%sim < huge(me%min%sim) .and.         &
      &  me%min%iter == huge(me%min%iter)) then
      me%min%iter = ceiling(me%min%sim/dt)
      ! if ( mod(me%min%iter, reqInt) /= 0 ) &
      !   & me%min%iter = me%min%iter + (reqInt - mod(me%min%iter, reqInt)
      me%min%iter = me%min%iter + mod((reqInt - mod(me%min%iter, reqInt)), &
        &                             reqInt)
    end if
    if(me%max%sim < huge(me%max%sim) .and.           &
      & me%max%iter == huge(me%max%iter)) then
      me%max%iter = ceiling(me%max%sim/dt)
      ! if ( mod(me%max%iter, reqInt) /= 0 ) &
      !   & me%max%iter = me%max%iter + reqInt - mod(me%max%iter, reqInt)
      me%max%iter = me%max%iter + mod((reqInt - mod(me%max%iter, reqInt)), &
        &                             reqInt)
    end if
    if(me%interval%sim < huge(me%interval%sim) .and. &
      & me%interval%iter == huge(me%interval%iter)) then
      me%interval%iter = int(me%interval%sim/dt)

      ! update trigger if suppose min is defined by simulation time
      me%trigger%iter = ceiling(me%trigger%sim/dt)
    end if

    ! set simulation time to never thus timeControl trigger is only
    ! based on iterations and not on simulation time
    me%min%sim = huge(me%min%sim)
    me%max%sim = huge(me%max%sim)
    me%interval%sim = huge(me%interval%sim)

    ! if iTime is defined set sim time
    !if(me%min%iter < huge(me%min%iter)) then
    !  me%min%sim = real(me%min%iter, kind=rk)*dt
    !endif
    !if(me%max%iter < huge(me%max%iter)) then
    !  me%max%sim = real(me%max%iter, kind=rk)*dt
    !endif

    write(logUnit(1),*) 'Time control homogenize:'
    call tem_timeControl_dump(me, logUnit(1))

  end subroutine mus_timeControl_homogenize
! **************************************************************************** !

! **************************************************************************** !
  !> Check for multilevel cycle complete by modulo of nIters by scaleFactor
  !! depends on acoustic or diffusive scaling.
  !! Acoustic scaling: scale factor = 2
  !! Diffusive scaling: scale factor = 4
  pure function mus_time_modulo(now, reqInt) result(triggered)
    ! -------------------------------------------------------------------------
    !> current simulation time
    type(tem_time_type), intent(in) :: now
    !> Required interval, in which the update MUST occur.
    !! This is required for the musubi multilevel, where the time step should
    !! only be determined active, when the end of the largest cycle is reached.
    integer, intent(in) :: reqInt
    ! -------------------------------------------------------------------------
    logical :: triggered
    ! -------------------------------------------------------------------------
    triggered = (mod(now%iter, reqInt) == 0)
  end function
! **************************************************************************** !


end module mus_time_module
! **************************************************************************** !
