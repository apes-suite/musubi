! Copyright (c) 2015-2016, 2018-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> author: Kannan Masilamani
!! Auxiliary functionality for musubi harvesting
!!
module mus_hvs_aux_module

  ! include musubi modules
  use tem_aux_module,                only: tem_abort
  use mus_param_module,              only: mus_param_type
  use mus_geom_module,               only: mus_geom_type
  use mus_tools_module,              only: check_streaming_layout
  use mus_interpolate_module,        only: mus_init_interpolate
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_interpolate_tools_module,  only: debug_dependencies, dump_intpLists
  use mus_time_module,               only: mus_time_homogenize
  use mus_fluid_module,              only: mus_init_fluid
  use mus_gradData_module,           only: mus_init_gradData
  use mus_tracking_module,           only: mus_init_tracker
  use mus_bndForce_module,           only: mus_init_BndForce

  ! include treelm modules
  use tem_debug_module,         only: dbgUnit
  use tem_restart_module,       only: tem_init_restart
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_tracking_module,      only: tem_tracker
  use tem_debug_module,         only: main_debug
  use tem_logging_module,       only: logUnit
  use tem_spacetime_fun_module, only: tem_create_subTree_of_st_funList

  implicit none
  private

  public :: mus_hvs_init_aux

contains

! ****************************************************************************** !
  !> Init auxiliary features such as interpolation boundaries, restart and
  !! the tracker
  subroutine mus_hvs_init_aux( scheme, geometry, params )
    ! ---------------------------------------------------------------------------
    !> container for the scheme
    type(mus_scheme_type), target, intent(inout)     :: scheme
    !> geometry infomation
    type(mus_geom_type), intent(inout)       :: geometry
    !> global parameters
    type(mus_param_type),intent(inout)       :: params
    ! ---------------------------------------------------------------------------
    integer :: iLevel, minLevel, maxLevel
    ! ---------------------------------------------------------------------------
    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel

    ! create subTree for all spacetime function in the linked list of
    ! spacetime function
    call tem_create_subTree_of_st_funList(     &
      &       me      = scheme%st_funList,     &
      &       tree    = geometry%tree,         &
      &       bc_prop = geometry%boundary,     &
      &       stencil = scheme%layout%fStencil )

    call tem_horizontalSpacer(fUnit = logUnit(1))

    ! verify correct settings for the streaming layout
    call check_streaming_layout( minLevel, maxLevel )

    ! When the restart is read from separate restart table, we start the
    ! simulation from the time step given in restart file. In the case,
    ! when restart is read from initial condition table, the simulation start
    ! time step is taken as the one defined in configuration file
    call mus_time_homogenize( me = params%general%simControl%now,              &
      &                       dt = params%physics%dtLvl( maxLevel ),           &
      &              readRestart = params%general%restart%controller           &
      &                             %readRestart )

    ! initialize the restart
    write(logUnit(1),*) 'Initializing restart...'
    call tem_init_restart( me           = params%general%restart, &
      &                    solver       = params%general%solver,  &
      ! &                    varSys       = scheme%varSys,          &
      &                    varMap       = scheme%stateVarMap,     &
      &                    tree         = geometry%tree           )

    ! scheme loop to initialize tracking, boundary, interpolation and output

    !> initialize fluid type which contains relaxation parameter
    !! and function pointers to get mrt paramter and nonEqScaling factor
    !! for interpolation
    select case( trim(scheme%header%kind) )
    case('fluid', 'fluid_incompressible', 'isotherm_acEq')
      if (scheme%nFields > 1) then
        call tem_abort('chosen scheme kind supports only one field')
      end if
      ! initialize fluid viscosity relaxation paramters
      call mus_init_fluid(                                &
        & me           = scheme%field(1)%fieldProp%fluid, &
        & physics      = params%physics,                  &
        & schemeHeader = scheme%header,                   &
        & minLevel     = minLevel,                        &
        & maxLevel     = maxLevel,                        &
        & levelDesc    = scheme%levelDesc(:),             &
        & pdf          = scheme%pdf(:),                   &
        & stencil      = scheme%layout%fStencil,          &
        & general      = params%general,                  &
        & tNow         = params%general%simControl%now    )
    end select

    ! Initialize gradient data. Required for LES tuburbulent and evaluating
    ! gradient of a variable
    if (scheme%readVarIsPdf) then
      allocate(scheme%gradData(minLevel:maxLevel))
      do iLevel = minLevel, maxLevel
        call mus_init_gradData( me        = scheme%gradData(iLevel),         &
          &                     neigh     = scheme%pdf(iLevel)%neigh(:),     &
          !&                     levelDesc = scheme%levelDesc(iLevel),        &
          &                     stencil   = scheme%layout%fStencil,          &
          &                     nSize     = scheme%pdf(iLevel)%nSize,        &
          &                     nSolve    = scheme%pdf(iLevel)%nElems_solve, &
          &                     nScalars  = scheme%varSys%nScalars           )
      end do
    end if

    ! initialize tracking objects.
    call mus_init_tracker( scheme    = scheme,   &
      &                    geometry  = geometry, &
      &                    params    = params    )

    if( minLevel /= maxlevel ) then
      write(logUnit(1),*) 'Initializing interpolation...'
      ! initialize the interpolation
      call mus_init_interpolate(                               &
        &             intp         = scheme%intp,              &
        &             leveldesc    = scheme%leveldesc,         &
        &             schemeHeader = scheme%header,            &
        &             stencil      = scheme%layout%fStencil,   &
        &             minLevel     = minLevel,                 &
        &             maxLevel     = maxLevel,                 &
        &             fieldProp    = scheme%field(:)%fieldProp )
    end if
    if( main_debug%debugDependencies) then
      call debug_dependencies( scheme%intp, scheme%levelDesc, &
        &                      geometry%tree, params%general%proc%rank )
      call dump_intpLists( minLevel, maxLevel, scheme%intp%config%order,   &
        &                  scheme%levelDesc, params%general%proc%rank )
    end if

    ! Boundary force calculation is valid only for single field schemes
    ! like fluid and fluid_incompressible so initialize only if nFields = 1
    if (geometry%boundary%nBCtypes > 0 .and. scheme%nFields==1) then
      call mus_init_bndForce(bndForce     = geometry%bndForce,  &
        &                    bndMoment    = geometry%bndMoment, &
        &                    bc_prop      = geometry%boundary,  &
        &                    schemeHeader = scheme%header,      &
        &                    bc           = scheme%field(1)%bc  )
    end if

    call tem_horizontalSpacer( after = 1, fUnit = logUnit(1) )

  end subroutine mus_hvs_init_aux
! ****************************************************************************** !


end module mus_hvs_aux_module
! ****************************************************************************** !
