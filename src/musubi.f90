! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2016, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2013, 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012-2013, 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2025 Tristan Vlogman <t.g.vlogman@utwente.nl>
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
!> M U S U B I
!! The lattice Boltzmann solver within the APES suite
!! (c) 2012 German Research School for Simulation Sciences GmbH
!! (c) 2013 University of Siegen.
!!
!! For a documentation, run ./waf gendoxy and find the documentation at
!! ./Documentation/html/index.html
program musubi
  ! include treelm modules
  use mpi
  use tem_general_module,           only: tem_start, tem_finalize
  use tem_adaptation_config_module, only: tem_adapt_type

  ! include musubi modules
  ! RESPECT THE ORDER !!!!!!!!
  use mus_aux_module,                only: mus_banner
  use mus_config_module,             only: mus_load_config
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_param_module,              only: mus_param_type
  use mus_timer_module,              only: mus_init_mainTimer,  &
    &                                      mus_init_levelTimer, &
    &                                      mus_init_bcTimer
  use mus_geom_module,               only: mus_geom_type
  ! ESPECIALLY OF THE FOLLOWING MODULES
  use mus_control_module,            only: mus_control_type
  use mus_program_module,            only: mus_initialize, mus_solve, &
    &                                      mus_finalize
  use mus_varSys_module,             only: mus_varSys_solverData_type

  ! include modules for coupled LBM-DEM simulations of solid particles
  use mus_particle_type_module,      only: mus_particle_group_type
  use mus_particle_timer_module,     only: mus_init_particleTimer

  implicit none
  ! -------------------------------------------------------------------------- !
  !> scheme types
  type(mus_scheme_type),            target :: scheme
  type(mus_geom_type),              target :: geometry
  type(mus_param_type),             target :: params
  type(mus_particle_group_type),    target :: particleGroup
  type(mus_varSys_solverData_type), target :: solverData
  type(tem_adapt_type)                     :: adapt
  type(mus_control_type)                   :: control
  integer :: ierr
  ! -------------------------------------------------------------------------- !

  control = mus_control_type( scheme        = scheme,       &
    &                         geometry      = geometry,     &
    &                         params        = params,       &
    &                         particleGroup = particleGroup )

  ! Initialize environment
  call tem_start(codeName   = 'Musubi',                 &
    &            general    = params%general,           &
    &            simControl = params%general%simControl )

  if (params%general%proc%rank == 0) then
    call mus_banner(solver = params%general%solver)
  end if

  call mus_init_mainTimer()
  call mus_init_particleTimer

  ! load configuration file
  call mus_load_config( scheme        = scheme,       &
    &                   solverData    = solverData,   &
    &                   geometry      = geometry,     &
    &                   params        = params,       &
    &                   adapt         = adapt,        &
    &                   particleGroup = particleGroup )

  ! KM: Do not move this init_levelTimer and init_bcTimer from here,
  ! Need to be here for apesmate
  call mus_init_levelTimer( geometry%tree%global%minLevel, &
    &                       geometry%tree%global%maxLevel  )
  call mus_init_bcTimer( geometry%boundary%nBCtypes )

  ! initialize musubi
  call mus_initialize( scheme        = scheme,        &
    &                  solverData    = solverData,    &
    &                  geometry      = geometry,      &
    &                  params        = params,        &
    &                  particleGroup = particleGroup, &
    &                  control       = control        )

  call mpi_barrier( MPI_COMM_WORLD, ierr )

  ! do main loop
  call mus_solve( scheme     = scheme,     &
    &             solverData = solverData, &
    &             geometry   = geometry,   &
    &             params     = params,     &
    &             control    = control,    &
    &             adapt      = adapt       )

  ! finalize musubi
  call mus_finalize( scheme        = scheme,                     &
    &                params        = params,                     &
    &                particleGroup = particleGroup,              &
    &                tree          = geometry%tree,              &
    &                nBCs          = geometry%boundary%nBCtypes, &
    &                levelPointer  = geometry%levelPointer,      &
    &                globIBM       = geometry%globIBM            )

  ! finalize treelm function like print run time info and mpi
  call tem_finalize(params%general)

end program musubi
