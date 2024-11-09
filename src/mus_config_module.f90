! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2012,2020-2021,2024 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011-2016, 2019-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
!> In this module, all parameter files are read in
!! as lua script or a sample configuration is being loaded
!!
!! Possible Parameter configuration
!!
!! - General Parameters
!! - [[mus_scheme_layout_module]] for Scheme Definitions
!! - [[tem_debug_module]] for Debug Parameters
!!
module mus_config_module

  ! include musubi modules
  use mpi
  use mus_param_module,         only: mus_param_type, mus_load_param, &
    &                                 mus_init_latticeUnit
  use mus_geom_module,          only: mus_geom_type, mus_load_geom
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_scheme_module,        only: mus_load_scheme
  use mus_physics_module,       only: mus_create_funcStr, mus_load_physics
  use mus_varSys_module,        only: mus_varSys_solverData_type
  use mus_timer_module,         only: mus_timerHandles
  use mus_tools_module,         only: dump_linear_partition

  ! include treelm modules
  use env_module,               only: pathLen, rk, long_k, globalMaxLevels
  use tem_aux_module,           only: tem_open_distconf_array, tem_abort
  use tem_timer_module,         only: tem_startTimer, &
    &                                 tem_stopTimer
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_general_module,       only: tem_load_general
  use tem_logging_module,       only: logUnit, tem_logging_load_primary
  use tem_comm_env_module,      only: tem_comm_env_type
  use tem_debug_module,         only: dbgUnit, tem_debug_load_main

  use tem_adaptation_config_module, only: tem_adapt_type, tem_load_adapt

  ! include aotus modules
  use aotus_module,    only: flu_State, open_config_chunk

  implicit none

  private

  public :: mus_load_config
  public :: mus_open_config


contains


  ! ************************************************************************** !
  !> Read in LUA parameter file
  !! See http://www.lua.org for a reference on how to use
  !! Lua is a scripting language in itself which allows
  !! more complex parameter files including comments
  !! And load / create the mesh depending on the configuration
  subroutine mus_load_config( scheme, solverData, geometry, params, adapt )
    ! ---------------------------------------------------------------------- !
    !> scheme type
    type( mus_scheme_type ), target :: scheme
    !> contains pointer to scheme, physics types
    type( mus_varSys_solverData_type ), target :: solverData
    !> Treelmesh data
    type( mus_geom_type ), intent(out), target :: geometry
    !> Global parameters
    type( mus_param_type ), target, intent(inout) :: params
    !> mesh adaptation
    type(tem_adapt_type), intent(inout) :: adapt
    ! ---------------------------------------------------------------------- !
    character(len=PathLen) :: filename
    integer :: minLevel, maxLevel
    ! ---------------------------------------------------------------------- !

    call tem_startTimer( timerHandle = mus_timerHandles%loadMesh )

    ! check whether params%general%solver%configFile is defined. When loading
    ! musubi from apes, musubi config file is provided via params
    if ( trim(params%general%solver%configFile) == '' ) then
      ! Get filename from command line argument
      call get_command_argument( 1,filename )
      if( trim(filename) == '' ) then
        filename = 'musubi.lua'
      end if

      params%general%solver%configFile = filename
    else
      filename = params%general%solver%configFile
    endif

    if (params%general%proc%rank == 0) then
      write(logUnit(1),*) "Loading configuration file: "//trim( filename )
      call tem_horizontalSpacer(fUnit = logUnit(1))
    end if

    ! open musubi config file and solver specific lua functions as chunk
    call mus_open_config( conf     = params%general%solver%conf, &
      &                   filename = filename,                   &
      &                   proc     = params%general%proc         )

    ! load and initialize logUnit
    call tem_logging_load_primary(conf = params%general%solver%conf(1), &
      &                           rank = params%general%proc%rank       )

    ! load and initialize debug unit
    call tem_debug_load_main(conf = params%general%solver%conf(1), &
      &                      rank = params%general%proc%rank       )

    ! load general information
    call tem_load_general(                               &
      &    me           = params%general,                &
      &    conf         = params%general%solver%conf(1), &
      &    solverAborts = params%mus_Aborts              )

    ! load global musubi params
    call mus_load_param( params = params,                       &
      &                  conf   = params%general%solver%conf(1) )

    ! load geometry information like mesh, boundary, immersed_boundary and
    ! restart. If restart read is defined then simControl is updated
    call mus_load_geom( me              = geometry,                  &
      &                 restart         = params%general%restart,    &
      &                 solverHead      = params%general%solver,     &
      &                 simControl      = params%general%simControl, &
      &                 proc            = params%general%proc,       &
      &                 scaleFactor     = params%scaleFactor,        &
      &                 initial_balance = params%initial_balance     )


    if ( params%dump_linear_partition ) then
      call dump_linear_partition( treeID = geometry%tree%treeID,     &
        &                         nElems = geometry%tree%nElems,     &
        &                         offset = geometry%tree%elemOffset, &
        &                         myRank = params%general%proc%rank, &
        &                         iter   = 0                         )
    end if

    ! minlevel and maxlevel of the mesh
    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel

    ! Initialize requiredInterval for multilevel to determine one complete cycle
    params%reqInterval = params%scaleFactor**(maxLevel-minLevel)

    ! Set lattice dx and dt according scaling type
    call mus_init_latticeUnit(params%lattice, minLevel, maxLevel, &
      &                       params%scaleFactor                  )

    ! load physics table for unit converstion
    call mus_load_physics( me          = params%physics,                &
      &                    conf        = params%general%solver%conf(1), &
      &                    tree        = geometry%tree,                 &
      &                    scaleFactor = params%scaleFactor             )

    ! The scheme information is loaded from the configuration file by invoking
    ! the function. THIS IS REQUIRED EVEN IN THE CASE OF DYNAMIC LOAD BALANCING
    ! DUE TO PRESENCE OF DERIVED VARIABLE LISTS BUILDING INSIDE THE LOAD SCHEME
    ! SINGLE ROUTINE
    call mus_load_scheme( me         = scheme,                        &
      &                   solverData = solverData,                    &
      &                   geometry   = geometry,                      &
      &                   conf       = params%general%solver%conf(1), &
      &                   params     = params                         )

    call tem_horizontalSpacer(fUnit = logUnit(1))

    ! load adaptation table
    call tem_load_adapt( me   = adapt,                        &
      &                  conf = params%general%solver%conf(1) )


    call tem_stopTimer( timerHandle = mus_timerHandles%loadMesh )

  end subroutine mus_load_config
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine loads musubi specific lua function from string and musubi
  !! input configuration file
  subroutine mus_open_config( conf, filename, proc )
    ! ---------------------------------------------------------------------- !
    !> lua state to be stored
    type(flu_State), allocatable :: conf(:)
    !> name of the config file to be opened
    character(len=*), intent(in) :: filename
    !> process description to use
    type(tem_comm_env_type), intent(in) :: proc
    ! ---------------------------------------------------------------------- !
    character(len=2048) :: fun_str !contains function to be used in musubi
    integer :: iThread
    ! ---------------------------------------------------------------------- !

    ! allocate the array of lua states
    allocate( conf( proc%nThreads ))

    call mus_create_funcStr( fun_str = fun_str )

    !pre-loading musubi specific lua functions using same conf used to load
    !main musubi config file
    do iThread=1, proc%nThreads
      call open_config_chunk(L = conf(iThread), chunk = fun_str)
    end do

    ! Open the lua config
    call tem_open_distconf_array( L        = conf,           &
      &                           filename = trim(filename), &
      &                           proc     = proc            )

  end subroutine mus_open_config
  ! ************************************************************************** !

end module mus_config_module
! **************************************************************************** !

