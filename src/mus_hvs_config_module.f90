! Copyright (c) 2015-2016, 2018, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! ***************************************************************************** !
!> In this module, all parameter files are read in
!! as lua script or a sample configuration is being loaded
!!
!! Possible Parameter configuration
!! - General Parameters
!! - [[mus_scheme_layout_module]] "Scheme Definitions"
!! - [[tem_debug_module]] "Debug Parameters"
!!
module mus_hvs_config_module

  ! include musubi modules
  use mpi
  use mus_param_module,              only: mus_param_type, mus_load_param
  use mus_geom_module,               only: mus_geom_type, mus_load_bc_data
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_scheme_module,             only: mus_load_scheme
  use mus_geomIncrHead_module,       only: mus_geomIncrHead_load
  use mus_physics_module,            only: mus_load_physics
  use mus_varSys_module,             only: mus_varSys_solverData_type
  use mus_config_module,             only: mus_open_config
  use mus_timer_module,              only: mus_timerHandles

  ! include treelm modules
  use env_module,             only: pathLen, rk, long_k, rk_mpi
  use tem_property_module,    only: prp_hasBnd, prp_hasQVal
  use tem_bc_prop_module,     only: init_tem_bc_prop, load_tem_BC_prop, &
    &                               load_tem_bc_qVal
  use tem_restart_module,     only: tem_load_restart
  use tem_timer_module,       only: tem_startTimer, &
    &                               tem_stopTimer
  use tem_tools_module,       only: tem_horizontalSpacer
  use tem_geometry_module,    only: tem_setEffBoundingBox
  use tem_general_module,     only: tem_load_general
  use tem_logging_module,     only: logUnit, tem_logging_load_primary
  use tem_timeControl_module, only: tem_timeControl_start_at_sim
  use tem_debug_module,       only: dbgUnit, tem_debug_load_main
  use tem_aux_module,         only: tem_abort, tem_open_distconf
  use tem_tracking_module,    only: tem_load_tracking

  use hvs_output_module, only: hvs_output_config_type, hvs_output_load

  ! include aotus modules
  use aotus_module,    only: flu_State, open_config_chunk, close_config,       &
    &                        aot_get_val

  implicit none

  private

  public :: mus_hvs_config_load
  public :: mus_hvs_config_type

  !> This datatype describes the various settings to load from the configuration
  !! file.
  type mus_hvs_config_type
    !> Location on disk to load the mesh data from.
    !!
    !! This prefix will be put before the various filenames of the individual
    !! mesh data files.
    character(len=pathLen) :: prefix

    !> Description of how the visualization output should be done.
    type(hvs_output_config_type) :: output

  end type mus_hvs_config_type

contains

! ****************************************************************************** !
  !> Read in LUA parameter file
  !! See http://www.lua.org for a reference on how to use
  !! Lua is a scripting language in itself which allows
  !! more complex parameter files including comments
  !! And load / create the mesh depending on the configuration
  subroutine mus_hvs_config_load( me, scheme, solverData, geometry, params )
    ! ---------------------------------------------------------------------------
    !> Musubi harvesting configuration to load when no tracking table is defined
    type( mus_hvs_config_type ), intent(out) :: me
    !> scheme type
    type( mus_scheme_type ),    target :: scheme
    !> contains pointer to scheme, physics types
    type( mus_varSys_solverData_type ), target :: solverData
    !> Treelmesh data
    type( mus_geom_type ), intent(inout), target    :: geometry
    !> Global parameters
    type( mus_param_type ), target, intent(inout)   :: params
    ! ---------------------------------------------------------------------------
    character(len=PathLen) :: filename
    integer :: minLevel, maxLevel
    integer :: iError
    ! ---------------------------------------------------------------------------

    call tem_startTimer( timerHandle =  mus_timerHandles%loadMesh )

    ! Load configuration data according to command line arguments.
    filename = ''

    ! Get filename from command line argument
    call get_command_argument( 1,filename )

    if( trim(filename) == '' ) then
      ! Default to harvester.lua, if no filename is provided on the command
      ! line.
      filename = 'harvester.lua'
    end if

    if (params%general%proc%rank == 0) then
      write(logUnit(1),*) "Loading configuration file: "//trim( filename )
      call tem_horizontalSpacer(fUnit = logUnit(1))
    end if

    ! open musubi config file and solver specific lua functions as chunk
    call mus_open_config( conf     = params%general%solver%conf,             &
      &                   filename = filename,                               &
      &                   proc     = params%general%proc )

    ! load and initialize logUnit
    call tem_logging_load_primary(conf = params%general%solver%conf(1), &
      &                           rank = params%general%proc%rank       )
    ! load and initialize debug unit
    call tem_debug_load_main(conf = params%general%solver%conf(1), &
      &                      rank = params%general%proc%rank       )

    ! load general information
    call tem_load_general( me   = params%general,                            &
      &                    conf = params%general%solver%conf(1))
    ! -------------------------------------------------------------------------
    ! First check, if we are starting from a restart
    call tem_load_restart( me             = params%general%restart,          &
      &                    conf           = params%general%solver%conf(1),   &
      &                    tree           = geometry%tree,                   &
      &                    timing         = params%general%simControl%now,   &
      &                    globProc       = params%general%proc )

    if ( .not. params%general%restart%controller%readRestart ) then
      write(logUnit(1),*) 'ERROR: No read restart given.'
      write(logUnit(1),*) 'Solution: Provide restart file in: '
      write(logUnit(1),*) '          restart = { read = <filename> }'
      call tem_abort()
    end if

    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel

    ! If there is a restart, the timings in the params type have to be
    ! updated to those read from the restart
    call tem_timeControl_start_at_sim(              &
      & me = params%general%simControl%timeControl, &
      & now = params%general%simControl%now         )

    ! Load boundary and qval
    call mus_load_bc_data( geometry, params%general%proc%rank, &
      &                              params%general%proc%comm  )

    ! load params, physics and solver specific info
    call mus_hvs_load_solverData( scheme     = scheme,     &
      &                           solverData = solverData, &
      &                           geometry   = geometry,   &
      &                           params     = params      )

    ! Initialize requiredInterval for multilevel to determine one complete cycle
    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel
    params%reqInterval =  params%scaleFactor**(maxLevel-minLevel)

    ! If tracking table is not defined, load output_folder key
    ! to dump restart input or mesh to disk
    if ( .not. scheme%track%control%active ) then
      ! load output format and other config for output from output table
      call hvs_output_load( me   = me%output,                    &
        &                   conf = params%general%solver%conf(1),&
        &                   isReduce = .false. )

      ! Load output folder
      call aot_get_val( L       = params%general%solver%conf(1), &
        &               key     = 'output_folder',               &
        &               val     = me%prefix,                     &
        &               default = './',                          &
        &               ErrCode = iError                         )
      write(logUnit(1),*) 'Output folder: '//trim(me%prefix)
    end if

    call tem_horizontalSpacer(fUnit = logUnit(1))

    call tem_stopTimer( timerHandle =  mus_timerHandles%loadMesh )

  end subroutine mus_hvs_config_load
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routines load solver data from config file except tracking
  subroutine mus_hvs_load_solverData( scheme, solverData, geometry, params )
    ! ---------------------------------------------------------------------------
    !> scheme type
    type( mus_scheme_type ),    target :: scheme
    !> contains pointer to scheme, physics types
    type( mus_varSys_solverData_type ), target :: solverData
    !> Treelmesh data
    type( mus_geom_type ), intent(inout), target    :: geometry
    !> Global parameters
    type( mus_param_type ), target, intent(inout)   :: params
    ! ---------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = logUnit(1))

    ! load global musubi params
    call mus_load_param( params = params, &
      &                  conf   = params%general%solver%conf(1) )

    ! load physics table for unit converstion
    call mus_load_physics( me          = params%physics,                &
      &                    conf        = params%general%solver%conf(1), &
      &                    tree        = geometry%tree,                 &
      &                    scaleFactor = params%scaleFactor             )

    ! Load basic scheme information from the restart file
    call mus_load_scheme( me         = scheme,                        &
      &                   solverData = solverData,                    &
      &                   geometry   = geometry,                      &
      &                   conf       = params%general%solver%conf(1), &
      &                   params     = params,                        &
      &                   isMusHvs   = .true.                         )

    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_hvs_load_solverData

! ****************************************************************************** !

end module mus_hvs_config_module
! ****************************************************************************** !

