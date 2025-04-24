! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011, 2013, 2015, 2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011-2016, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Nikhil Anand <n.anand@grs-sim.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
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
!> This module provides the definition and methods for
!! boundarys.
module mus_param_module

  ! include musubi modules
  use mus_physics_module,     only: mus_physics_type, set_values_by_levels
  use mus_abortCriteria_module, only: mus_abortCriteria_type

  ! include treelm modules
  use env_module,             only: labelLen, rk, globalMaxLevels
  use tem_logging_module,     only: logUnit
  use tem_general_module,     only: tem_general_type
  use tem_tools_module,       only: tem_horizontalSpacer

  ! include aotus modules
  use aotus_module,   only: flu_state, aot_get_val
  use aot_out_module, only: aot_out_type, aot_out_val

  implicit none

  private

  public :: mus_param_type
  public :: mus_load_param
  public :: mus_param_out
  public :: mus_init_latticeUnit
  public :: mus_latticeUnit_type

  integer, parameter :: block_size = 576

  !> lattice dx and dt on each level
  type mus_latticeUnit_type
    !> dt in LB unit, dtLvl(minLevel) = 1.0
    real(kind=rk), allocatable :: dtLvl(:)
    !> dx in LB unit, dxLvl(minLevel) = 1.0
    real(kind=rk), allocatable :: dxLvl(:)
  end type mus_latticeUnit_type

  !> Global parameter type definition, filled with
  type mus_param_type

    !> Treelm param parameter type
    type( tem_general_type ) :: general

    !> Musubi specific abort criteria
    type( mus_abortCriteria_type) :: mus_Aborts

    !> Lattice dx and dt on each level
    type(mus_latticeUnit_type) :: lattice

    !> contains basic SI units to convert from lattice to physical and
    !! vice versa
    type(mus_physics_type) :: physics

    !> Type of particle modeling to use:
    !!
    !! Available kinds:
    !! * Momentum-exchange method (kind = 'MEM')
    !! * Discrete Particle Simulations (kind = 'DPS')
    !! * One-way coupled Discrete Particle Simulations (kind = 'DPS_oneway')
    !! Default: 'none'
    character(len=labelLen) :: particle_kind = 'none'

    !> type of the control routine
    character(len=labelLen) :: controlRoutine
    !> initialize all elements with valid entries?
    !! This should only be activated for debugging,
    !! as it needs to be ensured that all helper elements
    !! are filled by communication and interpolation instead
    !! of filling initial values (consider restart!)
    logical :: init_allElems = .false.
    !> type of the scaling:
    !! * acoustic
    !! * diffusive
    character(len=labelLen) :: scaling
    !> Temporal scaling factor for the scaling. Acoustic = 2, Diffusive = 4
    integer :: scaleFactor
    !> Nesting is 2: acoustic, 4: diffusive
    !! To calculate turbulent viscosity, velocity on buffer ghost elements
    !! should be valid to nesting is set to same as scaling Factor
    integer :: nNesting = 2
    !> Required interval, in which the update MUST occur.
    !! This is required for the musubi multilevel, where the time step should
    !! only be determined active, when the end of the largest cycle is reached.
    integer :: reqInterval
    logical :: comm_reduced = .true. !< Communicate all links?
    !> need to set solver version in  general%solver%version
    character(len=labelLen) :: version = 'v2.0'
    !> active when restart is triggered by restart timeControl
    !! dump restart when simulation reached end only when
    !! restart is not triggered by its timeControl before
    logical :: restart_triggered = .false.

    !> remove solid from BC list
    logical :: remove_solid = .true.

    !> Block size for compute kernel
    integer :: block = block_size

    !> Initial balance
    logical :: initial_balance = .false.

    !> scratch file unit contains solver specific info in dump in restart header
    !! This file should contain the information in form of a Lua script.
    !! KM: Not required anymore. Load config file name from restart header
    !! integer :: solSpec_unit = -1

    !> Dump level timing
    logical :: dump_level_timing = .false.

    !> Dump linear partition
    logical :: dump_linear_partition = .false.

    !> Dump computation and bc timing information for all ranks
    logical :: dump_bc_timing = .false.

  end type mus_param_type


contains


! **************************************************************************** !
  !> load global parameter from conf
  subroutine mus_load_param( params, conf )
    ! --------------------------------------------------------------------------
    !> global parameter info
    type(mus_param_type), intent(inout)  :: params
    !> lua state
    type(flu_state) :: conf
    ! --------------------------------------------------------------------------
    integer :: iError
    ! --------------------------------------------------------------------------
    write(logUnit(1),*) 'Loading general solver params: '
    ! Load the control routine
    call aot_get_val( L       = conf,                                          &
      &               key     = 'control_routine',                             &
      &               val     = params%controlRoutine,                         &
      &               ErrCode = iError,                                        &
      &               default = 'standard' )
    write(logUnit(1),*)'Choosing Control Routine: ', trim(params%controlRoutine)

    ! load scaling type
    ! KM: @todo move scaling inside physics table because for lbm
    ! default scaling is acoustic and for multispecies default is diffusive
    call aot_get_val( L       = conf,           &
      &               key     = 'scaling',      &
      &               val     = params%scaling, &
      &               ErrCode = iError,         &
      &               default = 'acoustic'      )
    write(logUnit(1),*) 'Setting SCALING to:',trim( params%scaling )

    ! Set acoustic or diffusive scaling parameters
    !! KM: To calculate turbulent viscosity, velocity on buffer ghost elements
    !! should be valid to nesting is set to same as scaling Factor
    select case( trim( params%scaling ))
    case('acoustic')
      params%nNesting = 2
      params%scaleFactor = 2
    case('diffusive')
      params%nNesting = 4
      params%scaleFactor = 4
    end select

    ! Communicate reduced set of links only?
    call aot_get_val( L       = conf,                                          &
      &               key     = 'comm_reduced',                                &
      &               val     = params%comm_reduced,                           &
      &               ErrCode = iError,                                        &
      &               default = .true.)

    ! initialize all elements? this includes not only fluid but also ghost and
    ! halo
    call aot_get_val( L       = conf,                                          &
      &               key     = 'init_allElems',                               &
      &               val     = params%init_allElems,                          &
      &               ErrCode = iError,                                        &
      &               default = .false.)
    if( params%init_allElems ) then
      write(logUnit(1),*) 'WARNING: Initializing all (including helper '//     &
        &                           'elements) with valid data'
    end if

    call aot_get_val( L       = conf,                                          &
      &               key     = 'remove_solid',                                &
      &               val     = params%remove_solid,                           &
      &               ErrCode = iError,                                        &
      &               default = .true.)
    if( params%remove_solid ) then
      write(logUnit(1),"(A)") 'Solid elements will be removed from BC list'
    end if

    call aot_get_val( L       = conf,                     &
      &               key     = 'initial_balance',        &
      &               val     = params%initial_balance,   &
      &               ErrCode = iError,                   &
      &               default = .false.)

    call aot_get_val( L       = conf,                     &
      &               key     = 'dump_level_timing',      &
      &               val     = params%dump_level_timing, &
      &               ErrCode = iError,                   &
      &               default = .false.                   )

    call aot_get_val( L       = conf,                  &
      &               key     = 'dump_bc_timing',      &
      &               val     = params%dump_bc_timing, &
      &               ErrCode = iError,                &
      &               default = .false.                )

    call aot_get_val( L       = conf,                     &
      &               key     = 'dump_linear_partition',      &
      &               val     = params%dump_linear_partition, &
      &               ErrCode = iError,                   &
      &               default = .false.                   )

    call aot_get_val( L       = conf,         &
      &               key     = 'block',      &
      &               val     = params%block, &
      &               ErrCode = iError,       &
      &               default = block_size    )
    write(logUnit(1),"(A,I0)") 'Block size (elements): ', params%block

    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_load_param
! **************************************************************************** !


! **************************************************************************** !
  !> This routine writes global parameter into solver specific string in lua
  !! format
  subroutine mus_param_out( me, conf )
    ! --------------------------------------------------------------------------
    type( mus_param_type ), intent(in) :: me !< params
    type( aot_out_type ) :: conf
    ! --------------------------------------------------------------------------
    call aot_out_val( put_conf = conf, vname = 'scaling', &
      &               val = trim(me%scaling) )

  end subroutine mus_param_out
! **************************************************************************** !

! **************************************************************************** !
  !> This routine initialize lattice dx and dt
  subroutine mus_init_latticeUnit(lattice, minLevel, maxLevel, scaleFactor)
    ! -------------------------------------------------------------------------
    !> Lattice unit
    type(mus_latticeUnit_type), intent(out) :: lattice
    !> minlevel and maxlevel
    integer, intent(in) :: minLevel, maxLevel
    !> scaleFactor depending on acoustic or diffusive scaling
    integer, intent(in) :: scaleFactor
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------
    allocate(lattice%dxLvl(minLevel:maxLevel))
    allocate(lattice%dtLvl(minLevel:maxLevel))
    lattice%dxLvl( minLevel:maxLevel ) =                    &
      & set_values_by_levels( 1.0_rk, minLevel, maxLevel, 2 )

    lattice%dtLvl( minLevel:maxLevel ) =                              &
      & set_values_by_levels( 1.0_rk, minLevel, maxLevel, scaleFactor )

    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1), '(A)') 'Lattice dx and dt on each level'
    do iLevel = minLevel, maxLevel
      write(logUnit(1), '(A,I0,A,F10.5,A,F10.5)' ) 'level=', iLevel,     &
        & ', dxL=', lattice%dxLvl(iLevel), ', dtL=', lattice%dtLvl(iLevel)
    end do
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_init_latticeUnit
! **************************************************************************** !

end module mus_param_module
! **************************************************************************** !

