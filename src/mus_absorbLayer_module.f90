! Copyright (c) 2021 Kannan Masilamani <kannan.masilamani@dlr.de>
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
!> author: Kannan Masilamani
!! Module containing subroutines for building MUSUBI specific absorbLayer source
!! variables
!!
module mus_absorbLayer_module
  ! include treelm modules
  use env_module,               only: rk, labelLen
  use tem_tools_module,         only: upper_to_lower
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit
  use tem_tools_module,         only: tem_horizontalSpacer

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_Fatal, &
    &                         aoterr_NonExistent, aoterr_WrongType
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val

  implicit none
  private

  public :: mus_absorbLayer_type
  public :: mus_absorbLayer_dynAvg_type
  public :: mus_load_absorbLayer
  public :: mus_init_absorbLayer


  ! ************************************************************************** !
  !! Make target_pressure and target_velocity as stFun.
  !> Contains additional information for absorblayer source
  type absorbLayer_config_type
    !> target pressure
    real(kind=rk) :: target_pressure
    !> target velocityX, velocityY and velocityZ
    real(kind=rk) :: target_velocity(3)
    !> Use time average for pressure. Default: false.
    logical :: isPressDyn = .false.
    !> Use time average for Velocity. Default: false.
    logical :: isVelDyn = .false.
    !> Number of iterations to record for time-averaging
    integer :: nRecord = 100
  end type absorbLayer_config_type
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Stores time average values of density and velocity for dynamic absorbLayer
  type mus_absorbLayer_dynAvg_type
    !> density time average in lattice unit
    real(kind=rk), allocatable :: dens(:)
    !> velocity time average in lattice unit
    real(kind=rk), allocatable :: velX(:)
    real(kind=rk), allocatable :: velY(:)
    real(kind=rk), allocatable :: velZ(:)
    !> It is used to initialiye dynamic average density with initial condition
    logical :: isInitDens = .true.

    !> It is used to initialiye dynamic average velocity with initial condition
    logical :: isInitVel = .true.
  end type mus_absorbLayer_dynAvg_type
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Contains information for absorblayer
  type mus_absorbLayer_type
    !> Information loaded from configuration file
    type(absorbLayer_config_type) :: config

    !> Smoothing factor for expoential moving average
    !! = 2 / (nRecord+1)
    real(kind=rk) :: smoothFac
  end type mus_absorbLayer_type
  ! ************************************************************************** !


contains


  ! ************************************************************************** !
  !> This routine load additional information for absorblayer
  subroutine mus_load_absorbLayer(me, conf, key, parent, loadPres, loadVel)
    ! -------------------------------------------------------------------------!
    !> Absorb layer
    type(absorbLayer_config_type), intent(out) :: me
    !> flu state
    type( flu_State ) :: conf
    !> Table name to load target states
    character(len=*), intent(in) :: key
    !> parent source handle
    integer, intent(in) :: parent
    !> Load pressure if true else set to dynamic
    logical, intent(in) :: loadPres
    !> Load velocity if true else set to dynamic
    logical, intent(in) :: loadVel
    ! -------------------------------------------------------------------------!
    integer :: target_handle, iError, vError(3), errFatal(3)
    character(len=labelLen) :: tarStateAsStr
    ! -------------------------------------------------------------------------!
    errfatal = aotErr_Fatal

    call aot_table_open( L       = conf,                 &
      &                  parent  = parent,               &
      &                  thandle = target_handle,        &
      &                  key     = trim(key)             )

    ! Load target pressure and velocity only for static absorblayer
    write(logUnit(1),*) ' * Target state:'
    ! target state: pressure
    if (loadPres) then
      call aot_get_val( L       = conf,               &
        &               thandle = target_handle,      &
        &               key     = 'pressure',         &
        &               val     = me%target_pressure, &
        &               ErrCode = iError              )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1), *) 'Error loading target pressure as value'
        write(logUnit(5), *) 'Trying to load as string'
        call aot_get_val( L       = conf,          &
          &               thandle = target_handle, &
          &               key     = 'pressure',    &
          &               val     = tarStateAsStr, &
          &               ErrCode = iError         )
        if (btest(iError, aoterr_Fatal)) then
          write(logUnit(1),*) 'Error loading target pressure as string'
          call tem_abort()
        end if

        if (upper_to_lower(trim(tarStateAsStr)) == 'dynamic') then
          me%isPressDyn = .true.
          me%target_pressure = 0.0_rk
          write(logUnit(1),*) "    pressure = dynamic'"
        else
          call tem_abort('Target pressure is neither value or string "dynamic"')
        end if
      else
        write(logUnit(1),*) '    pressure =', me%target_pressure
      end if
    else
      me%isPressDyn = .true.
      me%target_pressure = 0.0_rk
      write(logUnit(1),*) "    pressure = 'dynamic'"
    end if

    ! target state: velocity
    if (loadVel) then
      call aot_get_val( L       = conf,               &
        &               thandle = target_handle,      &
        &               key     = 'velocity',         &
        &               val     = me%target_velocity, &
        &               ErrCode = vError              )
      if (any(btest( vError, errFatal )) ) then
        write(logUnit(1), *) 'Error loading target velocity as table'
        write(logUnit(5), *) 'Trying to load as string'
        call aot_get_val( L       = conf,          &
          &               thandle = target_handle, &
          &               key     = 'velocity',    &
          &               val     = tarStateAsStr, &
          &               ErrCode = iError         )
        if (btest(iError, aoterr_Fatal)) then
          write(logUnit(1),*) 'Error loading target velocity as string'
          call tem_abort()
        end if

        if (upper_to_lower(trim(tarStateAsStr)) == 'dynamic') then
          me%isVelDyn = .true.
          me%target_velocity = 0.0_rk
          write(logUnit(1),*) "    velocity = dynamic'"
        else
          call tem_abort('Target velocity is neither value or string "dynamic"')
        end if

      else
        write(logUnit(1),*) '    velocity =', me%target_velocity(1:3)
      end if
    else
      me%isVelDyn = .true.
      me%target_velocity = 0.0_rk
      write(logUnit(1),*) "    velocity = 'dynamic'"
    end if

    ! Load nRecord for time averaging
    if (me%isPressDyn .or. me%isVelDyn) then
      call aot_get_val( L       = conf,          &
        &               thandle = target_handle, &
        &               key     = 'nrecord',     &
        &               val     = me%nRecord,    &
        &               ErrCode = iError         )
      if (btest(iError, aoterr_Fatal)) then
        write(*,*) 'FATAL Error occured, when loading nrecord for absorb layer'
        call tem_abort()
      end if
      write(logUnit(1),*) '    nRecord =', me%nRecord
    end if

    call aot_table_close(conf, target_handle)

  end subroutine mus_load_absorbLayer
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Initialize arrays to store time average density and velocity for
  !! dynamic absorbing layer.
  !! \todo KM: 20210301 Allocate only pressure or velocity depending on
  !! absorb_layer_inlet or absorb_layer_outlet
  subroutine mus_init_absorbLayer(absLayer, dynAvg, nElems)
    ! --------------------------------------------------------------------------
    !> Absorblayer type
    type(mus_absorbLayer_type), intent(inout) :: absLayer
    !> Contains dynamic average density and velocity for absorblayer
    type(mus_absorbLayer_dynAvg_type), intent(inout) :: dynAvg
    !> Number of source elements
    integer, intent(in) :: nElems
    ! --------------------------------------------------------------------------
    if (absLayer%config%isPressDyn .or. absLayer%config%isVelDyn) then
      allocate(dynAvg%dens(nElems))
      allocate(dynAvg%velX(nElems))
      allocate(dynAvg%velY(nElems))
      allocate(dynAvg%velZ(nElems))
      dynAvg%dens(:) = 0.0_rk
      dynAvg%velX(:) = 0.0_rk
      dynAvg%velY(:) = 0.0_rk
      dynAvg%velZ(:) = 0.0_rk

      dynAvg%isInitDens = .true.
      dynAvg%isInitVel = .true.
      ! smooth factor for expontential average
      absLayer%smoothFac = 2.0_rk/(absLayer%config%nRecord+1.0_rk)
    end if


  end subroutine mus_init_absorbLayer
  ! ************************************************************************** !

end module mus_absorbLayer_module
