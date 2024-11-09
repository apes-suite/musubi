! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!! This module contains the datatypes for the geometry increase table within the
!! lua configuration file. It also contains the routines to load various
!! variables, parameters from the lua file which are then used in performing
!! geometry changes
!!
module mus_geomIncrHead_module

  ! include treelm modules
  use env_module,             only: labelLen
  use tem_timeControl_module, only: tem_timeControl_type, tem_timeControl_load
  use tem_tools_module,       only: tem_horizontalSpacer
  use tem_depend_module,      only: tem_load_depend, tem_depend_type
  use tem_logging_module,     only: logUnit

  ! include aotus modules
  use aotus_module,      only: flu_State, aot_get_val,           &
    &                          aoterr_NonExistent, aoterr_Fatal, &
    &                          aoterr_WrongType
  use aot_table_module,  only: aot_get_val, aot_table_open, aot_table_close,   &
    &                          aot_table_length
  use aot_vector_module, only: aot_get_val

  ! include musubi modules
  use mus_param_module,       only: mus_param_type

  implicit none
  private

  public :: mus_geomIncrHead_type
  public :: mus_geomIncrHead_load

  type mus_geomIncrHead_type
    logical :: active     = .false.  !< geometry increment active?
    logical :: solidify   = .false.  !< if true solidify, if false fluidify
    logical :: fluidify   = .false.
    logical :: proximity  = .false.  !< Check proximity condition if true
    type(tem_timeControl_type) :: timeControl

    !> Name of the variable defined for condition varname in config file.
    !! Variable refered should return 0 for if condition is false and
    !  1 for true.
    !! If there are more than one condition variable required then they must
    !! be combined via variable definitions in config file.
    character(len=labelLen) :: cond_varName
    !> Position of variable defined for the condition varname in the varSys
    integer :: cond_varPos
  end type mus_geomIncrHead_type

contains

! ****************************************************************************** !
  !> Read all the necessary information for the geometry increase from the lua
  !! config file. This routine basically provides as a wrapper to the routine
  !! which reads single values
  !!
  !! Example to geomIncr table:
  !!```lua
  !! variable = {
  !!   { name = 'vel_threshold',
  !!     ncomponents = 1,
  !!     var_type = st_fun,
  !!     st_fun = 0.01
  !!   },
  !!   { name = 'incr_condition',
  !!     ncomponents = 1,
  !!     var_type = 'operation' ,
  !!     operation = {
  !!       kind = '<',
  !!       input_varname = {'vel_mag','vel_threshold'}
  !!     }
  !!   },
  !! }
  !! geomIncr = {
  !!  condition = 'incr_condition'
  !! }
  !!```
  subroutine mus_geomIncrHead_load( me, conf, parent, dynamicGeom)
    ! ---------------------------------------------------------------------------
    type( mus_geomIncrHead_type), allocatable, intent(inout) :: me(:)
    type( flu_state ) :: conf
    integer, optional, intent(in) :: parent
    logical, intent(inout)        :: dynamicGeom
    ! ---------------------------------------------------------------------------
    integer :: tc_handle, sub_handle
    integer :: nGeomIncrs
    integer :: iGInc
    ! ---------------------------------------------------------------------------

    ! Attempt to open the geomIncr table (within another table, if a parent is
    ! given )
    call aot_table_open( L       = conf,      &
      &                  parent  = parent,    &
      &                  thandle = tc_handle, &
      &                  key     = 'geomIncr' )

    ! Check if geomIncr is actually defined
    if (tc_handle /= 0) then
      ! Set the dynamicGeom flag to true
      dynamicGeom = .true.
      ! check whether there are other members inside geomIncr
      call aot_table_open( L       = conf,       &
        &                  parent  = tc_handle,  &
        &                  thandle = sub_handle, &
        &                  pos     = 1           )
      ! If there is only one member in geomIncr, call the load routine once
      if(sub_handle == 0) then
        allocate(me(1))
        call aot_table_close( L = conf, thandle = sub_handle)
        call mus_geomIncrHead_load_single( me      = me(1),    &
          &                                conf    = conf,     &
          &                                thandle = tc_handle )
      else
        ! multiple definitions inside geomIncr
        call aot_table_close( L = conf, thandle = sub_handle)
        nGeomIncrs = aot_table_length( L = conf, thandle = tc_handle)
        allocate(me(nGeomIncrs))

        ! Read the sub tables individually
        do iGInc = 1, nGeomIncrs
          call aot_table_open( L       = conf,       &
            &                  parent  = tc_handle,  &
            &                  thandle = sub_handle, &
            &                  pos     = iGInc       )
          call mus_geomIncrHead_load_single( me      = me( iGInc ), &
            &                                conf    = conf,        &
            &                                thandle = sub_handle   )
        end do
      end if
    else
      dynamicGeom = .false.
    end if

    call aot_table_close(L=conf, thandle=tc_handle)

  end subroutine mus_geomIncrHead_load
! ****************************************************************************** !


! ****************************************************************************** !
  !> Reads various parameters from the lua file defined for geometry increase
  !! This routine reads single values and is wrapped around in another function
  !! where it is called multiple times as required
  !!
  subroutine mus_geomIncrHead_load_single( me, conf, thandle )
    ! ---------------------------------------------------------------------------
    type( mus_geomIncrHead_type),intent(inout)  :: me
    type( flu_state ), intent(in) :: conf
    integer, intent(in) :: thandle
    ! ---------------------------------------------------------------------------
    integer :: iError
    ! ---------------------------------------------------------------------------

    me%active = .true.
    ! Read the solidify, fluidify and proximity flags
    call aot_get_val( L       = conf,        &
      &               thandle = thandle,     &
      &               val     = me%solidify, &
      &               ErrCode = iError,      &
      &               key     = 'solidify',  &
      &               default = .false.      )

    call aot_get_val( L       = conf,        &
      &               thandle = thandle,     &
      &               val     = me%fluidify, &
      &               ErrCode = iError,      &
      &               key     = 'fluidify',  &
      &               default = .false.      )

    call aot_get_val( L       = conf,         &
      &               thandle = thandle,      &
      &               val     = me%proximity, &
      &               ErrCode = iError,       &
      &               key     = 'proximity',  &
      &               default = .false.       )

    ! Load variable name for condition.
    ! This variable should be logical operation variable which
    ! returns 0 for false and 1 for true.
    call aot_get_val( L       = conf,            &
      &               thandle = thandle,         &
      &               val     = me%cond_varName, &
      &               ErrCode = iError,          &
      &               key     = 'condition'      )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving "condition" ' &
        &                //'variable name from operation table :'
      if ( btest( iError, aotErr_NonExistent ))                                &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
    end if

    ! load time control to perform geomIncr
    call tem_timeControl_load( conf           = conf,                          &
      &                        parent         = thandle,                       &
      &                        me             = me%timeControl )

    ! If geomIncr is activited, output its configuration to console
    if( me%active ) then
      write(logUnit(1),*)'Geometry Increment ACTIVATED! '
      if( me%solidify ) then
        write(logUnit(1),*)'  Solidification is activited!'
      else
        write(logUnit(1),*)'  Fluidification is activited!'
      endif
      if( me%proximity) then
        write(logUnit(1),*)'  Proximity is activited.'
      else
        write(logUnit(1),*)'  Proximity is activited.'
      endif
      call tem_horizontalSpacer(fUnit = logUnit(1))
    endif

  end subroutine mus_geomIncrHead_load_single
! ****************************************************************************** !


end module mus_geomIncrHead_module
! ****************************************************************************** !
