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
!> mus_particle_config_module contains routines for loading the data for LBM-DEM
!! simulations of particles in a flow from the lua script.
!! All particle data is read from the "particles" table in the lua script.
!! This table should look as follows:
!! particles = {
!!  nParticles = 1,        -- Total number of particles in the simulation

!!  -- particle kind, can be:
!!  -- * 'MEM' for fully-resolved particles using the Momentum-Exchange Method
!!  -- * 'DPS' for four-way coupled unresolved particles based on the
!!  --   Volume-Averaged Navier-Stokes equations
!!  -- * 'DPS_twoway' for two-way coupled unresolved particles, neglecting the
!!  --   effect of volume fraction
!!  -- * 'DPS_oneway' for one-way coupled unresolved particles. That is, particles
!!  --   do not affect the flow
!!  kind = 'DPS',
!!
!!  -- Describe how particles should interact with boundaries. They can either bounce
!!  -- off them ('wall') or have periodic boundaries ('periodic'). Currently only
!!  -- interactions with prismatic boundaries are supported. If not specified,
!!  -- particles will be removed from the domain once they hit a boundary of the
!!  -- fluid domain (also works for arbitrary geometries).
!!  boundaries = {
!!    domain_bnd = { 0.0, L, 0.0, W, 0.0, H },
!!    bnd_kind = {'wall', 'wall', 'wall', 'wall', 'wall', 'wall'}
!!  },
!!
!!  -- Number of DEM subcycles per LBM time step to use
!!  nDEMsubcycles = 50,
!!
!!  -- How often to log the particle data (1 = log every LBM time step)
!!  particleLogInterval = 1,
!!
!!  -- Size of the buffers in the particle communication routines. Should be
!!  -- sufficient to hold all particles anticipated to be on a domain boundary
!!  -- at any time during the simulation. Safest option is to make this equal
!!  -- to or greater than nParticles.
!!  particleBufferSize = 100,
!!
!!  -- Collision time (in physical units) of particle-particle collisions
!!  particle_collision_time = 0.1*dt,
!!
!!  -- Collision tolerance, i.e. what should the gap between particles be
!!  -- before we call it a collision
!!  particle_collision_tol = 0.0,
!!
!!  -- Size of the gap before capping the calculation of lubrication forces
!!  critical_gap_Flub = 0.01,
!!
!!  -- Particle positions. Should contain nParticles tables of the form
!!  -- {x,y,z,rx,ry,rz} with x, y, z the translational positions and
!!  -- rx, ry, rz the rotation angles about the corresponding axes.
!!  position = { { x_p_phy, y_p_phy, z_p_phy, 0.0, 0.0, 0.0} },
!!
!!  -- Particle initial velocities. Should contain nParticles tables of the form
!!  -- {ux, uy, uz, urx, ury, urz}
!!  velocity = { {0.0, 0.0, 0.0, 0.0, 0.0, 0.0} },
!!
!!  -- External forces on particles. Should contain nParticles tables of the form
!!  -- {Fx, Fy, Fz, Frx, Fry, Frz}
!!  force    = { {0.0, 0.0, -F_gravity + F_buoyancy, 0.0, 0.0, 0.0} },
!!
!!  -- Particle radii, should have nParticles entries
!!  radius = { 0.5*Dia_p_phy },
!!
!!  -- Particle masses, should have nParticles entries
!!  mass = { m_p_phy }
!!}















! --- Macros for loading particle data from lua script --- !



! Module to load particle parameters from lua file
module mus_particle_config_module
  use env_module, only: rk, long_k, labelLen

  use tem_grow_array_module, only: append, destroy, empty
  use tem_aux_module, only: tem_abort, check_aot_error
  use tem_timer_module, only: tem_startTimer, tem_stopTimer
  use tem_logging_module, only: logUnit

  use flu_kinds_module, only: double_k
  use aotus_module, only: flu_State, open_config_file,   &
    &                     close_config, aot_get_val,     &
    &                     aot_top_get_val, aoterr_Fatal, &
    &                     aoterr_WrongType, aoterr_NonExistent
  use aot_table_module, only: aot_table_open, aot_table_close, &
    &                         aot_table_length
  use aot_fun_module, only: aot_fun_type, aot_fun_open, &
    &                       aot_fun_put,                &
    &                       aot_fun_do, aot_fun_close
  use aot_out_module, only: aot_out_type, aot_out_open,        &
    &                       aot_out_close, aot_out_open_table, &
    &                       aot_out_close_table, aot_out_val

  use mus_scheme_type_module, only: mus_scheme_type
  use mus_geom_module,        only: mus_geom_type
  use mus_param_module,       only: mus_param_type

  use mus_particle_logging_type_module, only: &
    &   mus_particle_debugtracking_type,      &
    &   debugTracker,                         &
    &   printDebugTrackerData
  use mus_particle_type_module, only: &
    &   mus_particle_MEM_type,        &
    &   mus_particle_DPS_type,        &
    &   mus_particle_group_type,      &
    &   dyn_particle_MEM_array_type,  &
    &   dyn_particle_DPS_array_type,  &
    &   init_da_particle_MEM,         &
    &   append_da_particle_MEM,       &
    &   destroy_da_particle_MEM,      &
    &   init_da_particle_DPS,         &
    &   append_da_particle_DPS,       &
    &   destroy_da_particle_DPS,      &
    &   maxContainerSize
  use mus_particle_boundary_module, only: &
    &   pgBndData, &
    &   mus_particle_boundarydata_type
  use mus_particle_interpolation_module, only: &
    &   mus_particle_interpolator_type, &
    &   intp_1D_delta, &
    &   intp_1D_peskin, &
    &   getwght1d_linear, &
    &   one, &
    &   printParticleInterpolator
  use mus_particle_timer_module, only: mus_particle_timerHandles
  use mus_particle_aux_module, only: getProcessBoundingBox, &
    &   positionLocalOnMyRank
  use mus_particle_blob_module, only:    &
    &   mus_particle_blob_prob_type,     &
    &   mus_particle_blob_cylinder_type, &
    &   mus_particle_blob_prism_type,    &
    &   particleblob,                    &
    &   particleblob_prism,              &
    &   print_particleblob,              &
    &   print_particleblob_prism
  use mus_particle_creator_module, only:       &
    &   mus_particle_creator_type,             &
    &   print_particle_creator,                &
    &   print_particle_creator_positions,      &
    &   init_particle_creator,                 &
    &   init_particle_creator_from_blob,       &
    &   init_particle_creator_from_blob_prism, &
    &   particle_creator

  implicit none


contains


  subroutine mus_load_particle_creator_timing(    &
    &          conf, p_thandle, particle_creator, &
    &          Nparticles, chunkSize, scheme,     &
    &          geometry, myRank                   )
    !> configuration
    type(flu_State) :: conf
    !> Handle to particle table
    integer, intent(in) :: p_thandle
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: particle_creator
    !> Number of particles in lua table
    integer, intent(in) :: Nparticles
    !> Size of the chunks of particle data to read from the lua script
    integer, intent(in) :: chunkSize
    !> Scheme
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry
    type(mus_geom_type), intent(in) :: geometry
    !> This MPI process rank
    integer, intent(in) :: myRank
    !> logical
    ! --------------------------------------------- !
    integer :: thandle ! handle to timing subtable
    integer :: ierror
    ! --------------------------------------------- !
    ! Open "timing" subtable within particles table
    call aot_table_open( l       = conf,      &
                       & thandle = thandle,   &
                       & parent  = p_thandle, &
                       & key     = "timing"   )

    ! Load the particle creator timing data
    if (thandle /= 0) then
      call aot_get_val( l       = conf,                        &
                      & thandle = thandle,                     &
                      & val     = particle_creator%iter_start, &
                      & errcode = ierror,                      &
                      & pos     = 1,                           &
                      & default = 0                            )

      call aot_get_val( l       = conf,                      &
                      & thandle = thandle,                   &
                      & val     = particle_creator%iter_end, &
                      & errcode = ierror,                    &
                      & pos     = 2,                         &
                      & default = 0                          )

      ! Interval: value of 0 means only initialize particles once
      call aot_get_val( l       = conf,                           &
                      & thandle = thandle,                        &
                      & val     = particle_creator%iter_interval, &
                      & errcode = ierror,                         &
                      & pos     = 3,                              &
                      & default = 1                               )
    else
      ! If no timing subtable is present, load the default values
      ! This will only create particles once, at the start of the simulation
      particle_creator%iter_start = 1
      particle_creator%iter_end = 1
      particle_creator%iter_interval = 1
    end if ! thandle /= 0

    ! Allocate space in particle creator arrays to hold particle positions, velocity data etc.
    call init_particle_creator( me = particle_creator )

  end subroutine mus_load_particle_creator_timing


  subroutine mus_load_particle_interpolation(conf, p_thandle, layout, interpolator)
      !> configuration
      type(flu_State) :: conf
      !> Handle to parent table
      integer, intent(in) :: p_thandle
      !> Layout label, e.g. d3q19, d2q9
      character(len=labelLen) :: layout
      !> tracker to initialize values of
      type(mus_particle_interpolator_type), intent(inout) :: interpolator
      ! -----------------------------------------!
      character(len=labelLen) :: interpolation_kind
      integer :: iError
      ! -----------------------------------------!
      select case( trim(layout) )
        case( 'd3q19', 'd3q27' )
          interpolator%bnd_x(1) = -1
          interpolator%bnd_x(2) = 1
          interpolator%bnd_y(1) = -1
          interpolator%bnd_y(2) = 1
          interpolator%bnd_z(1) = -1
          interpolator%bnd_z(2) = 1
        case( 'd2q9' )
          interpolator%bnd_x(1) = -1
          interpolator%bnd_x(2) = 1
          interpolator%bnd_y(1) = -1
          interpolator%bnd_y(2) = 1
          interpolator%bnd_z(1) = 0
          interpolator%bnd_z(2) = 0
        case default
          write(logUnit(1),*) 'ERROR mus_load_particle_interpolation layout kind unknown'
          write(logUnit(1),*) 'ABORTING'
          call tem_abort()
        end select

      call aot_get_val( L       = conf,                 &
        &               thandle = p_thandle,            &
        &               key     = 'interpolation_kind', &
        &               val     = interpolation_kind,   &
        &               default = 'delta',              &
        &               ErrCode = iError                )

      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*) 'FATAL Error occured, while retrieving particle interpolation kind:'
        if (btest(iError, aoterr_NonExistent)) write(*,*) 'Variable nonexistent!'
        if (btest(iError, aoterr_WrongType)) write(*,*) 'Variable of wrong type!'
      end if

      interpolator%interpolation_kind = interpolation_kind

      select case( trim(interpolator%interpolation_kind) )
        case( 'delta' )
          interpolator%getWght_x => intp_1D_delta
          interpolator%getWght_y => intp_1D_delta
          if( trim(layout) == 'd2q9' ) then
            interpolator%getWght_z => one
          else
            interpolator%getWght_z => intp_1D_delta
          end if

        case( 'linear' )
          interpolator%getWght_x => getwght1d_linear
          interpolator%getWght_y => getwght1d_linear
          if( trim(layout) == 'd2q9' ) then
            interpolator%getWght_z => one
          else
            interpolator%getWght_z => getwght1d_linear
          end if
        case( 'peskin' )
          ! For peskin stencil also interpolation boundaries change
          ! depending on layout kind
          if( trim(layout) == 'd2q9' ) then
            interpolator%bnd_x(1) = -2
            interpolator%bnd_x(2) = 2
            interpolator%bnd_y(1) = -2
            interpolator%bnd_y(2) = 2
            interpolator%bnd_z(1) = 0
            interpolator%bnd_z(2) = 0
            interpolator%getWght_x => intp_1D_peskin
            interpolator%getWght_y => intp_1D_peskin
            interpolator%getWght_z => one
          else
            interpolator%bnd_x(1) = -2
            interpolator%bnd_x(2) = 2
            interpolator%bnd_y(1) = -2
            interpolator%bnd_y(2) = 2
            interpolator%bnd_z(1) = -2
            interpolator%bnd_z(2) = 2
            interpolator%getWght_x => intp_1D_peskin
            interpolator%getWght_y => intp_1D_peskin
            interpolator%getWght_z => intp_1D_peskin
          end if
        case default
          write(logUnit(1),*) 'ERROR mus_load_particle_interpolation layout kind unknown'
          write(logUnit(1),*) 'ABORTING'
          call tem_abort()
      end select


  end subroutine mus_load_particle_interpolation

  subroutine mus_load_debugtracker(conf, p_thandle, tracker)
    !> configuration
    type(flu_State) :: conf
    !> Handle to parent table
    integer, intent(in) :: p_thandle
    !> tracker to initialize values of
    type(mus_particle_debugtracking_type), intent(inout) :: tracker
    ! -----------------------------------------!
    integer :: thandle, tthandle, iError
    integer :: nVals, iVal
    real(kind=rk) :: realBuffer
    integer :: intBuffer
    character(len=labelLen) :: stringBuffer

    ! -----------------------------------------!
    call aot_table_open(L       = conf,      &
      &                 thandle = thandle,   &
      &                 parent  = p_thandle, &
      &                 key     = 'tracker'  )

    ! If tracker table is available, load its data
    if (thandle /= 0) then
      ! ---------------- LOAD FILE NAME PREFIX ----------------- !
      ! Actual filename will be this with _tXXXX.dat appended to it
      call aot_get_val(L = conf, thandle = thandle, key = 'name', &
        &              val = stringBuffer, ErrCode = iError       )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*) 'FATAL Error occured, while retrieving tracker{name}:'
        if (btest(iError, aoterr_NonExistent)) write(*,*) 'Variable nonexistent!'
        if (btest(iError, aoterr_WrongType)) write(*,*) 'Variable of wrong type!'
      end if

      tracker%lfile = stringBuffer
      ! ---------------- LOAD TRACKING TYPE ----------------- !
      ! Can be 'line' or 'plane'
      call aot_get_val(L = conf, thandle = thandle, key = 'trackertype', &
        &              val = stringBuffer, ErrCode = iError           )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*) 'FATAL Error occured, while retrieving tracker{trackertype}:'
        if (btest(iError, aoterr_NonExistent)) write(*,*) 'Variable nonexistent!'
        if (btest(iError, aoterr_WrongType)) write(*,*) 'Variable of wrong type!'
      end if

      tracker%trackertype = stringBuffer

      ! ----------- OPEN DIR1 TABLE ---------- !
      call aot_table_open(L       = conf,     &
        &                 thandle = tthandle, &
        &                 parent  = thandle,  &
        &                 key     = 'dir1'    )
      if (tthandle /= 0) then
        ! get the number of values in dir1 table, should be 3
        nVals = aot_table_length(L=conf, thandle=tthandle)
        if(nVals /= 3) then
          write(*,*) 'FATAL Error occured, dir1 table is not length 3'
          write(*,*) 'use e.g. dir1 = {1, 0, 0} for x-direction'
          call tem_abort()
        end if ! nVals /= 3

        do iVal = 1, nVals
          call aot_get_val(L = conf, thandle = tthandle, &
            &              val = intBuffer, ErrCode = iError, &
            &              pos = iVal)
          if (btest(iError, aoterr_Fatal)) then
            write(*,*) 'FATAL Error occured, while retrieving dir1'
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not existent!'
            if (btest(iError, aoterr_WrongType)) write(*,*) &
              &  'Variable has wrong type!'
          else
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not set in' &
              &  // ' config, Using default value!'

          end if ! btest(iError, aoterror_fatal)

          ! Store table value in tracker
          tracker%dir1(iVal) = intBuffer

        end do ! Loop over dir1 table

      end if ! tthandle /= 0
      call aot_table_close(L = conf, thandle = tthandle)

      ! ----------- OPEN DIR2 TABLE ---------- !
      call aot_table_open(L       = conf,     &
        &                 thandle = tthandle, &
        &                 parent  = thandle,  &
        &                 key     = 'dir2'    )
      if (tthandle /= 0) then
        ! get the number of values in dir2 table, should be 3
        nVals = aot_table_length(L=conf, thandle=tthandle)
        if(nVals /= 3) then
          write(*,*) 'FATAL Error occured, dir2 table is not length 3'
          write(*,*) 'use e.g. dir2 = {1, 0, 0} for x-direction'
          call tem_abort()
        end if ! nVals /= 3

        do iVal = 1, nVals
          call aot_get_val(L = conf, thandle = tthandle,      &
            &              val = intBuffer, ErrCode = iError, &
            &              pos = iVal                         )
          if (btest(iError, aoterr_Fatal)) then
            write(*,*) 'FATAL Error occured, while retrieving dir2'
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not existent!'
            if (btest(iError, aoterr_WrongType)) write(*,*) &
              &  'Variable has wrong type!'
          else
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not set in' &
              &  // ' config, Using default value!'

          end if ! btest(iError, aoterror_fatal)

          ! Store table value in tracker
          tracker%dir2(iVal) = intBuffer

        end do ! Loop over dir2 table

      end if ! tthandle /= 0
      call aot_table_close(L = conf, thandle = tthandle)

      ! ---------------- LOAD XSTART ----------------- !
      call aot_table_open(L       = conf,     &
        &                 thandle = tthandle, &
        &                 parent  = thandle,  &
        &                 key     = 'xstart'  )
      if (tthandle /= 0) then
        ! get the number of values in xstart table, should be 3
        nVals = aot_table_length(L=conf, thandle=tthandle)
        if(nVals /= 3) then
          write(*,*) 'FATAL Error occured, xstart table is not length 3'
          write(*,*) 'use e.g. xstart = {0.5, 0.5, 0.5}'
          call tem_abort()
        end if ! nVals /= 3

        do iVal = 1, nVals
          call aot_get_val(L = conf, thandle = tthandle, &
            &              val = realBuffer, ErrCode = iError, &
            &              pos = iVal)
          if (btest(iError, aoterr_Fatal)) then
            write(*,*) 'FATAL Error occured, while retrieving xstart'
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not existent!'
            if (btest(iError, aoterr_WrongType)) write(*,*) &
              &  'Variable has wrong type!'
          else
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not set in' &
              &  // ' config, Using default value!'

          end if ! btest(iError, aoterror_fatal)

          ! Store table value in tracker
          tracker%xstart(iVal) = realBuffer

        end do ! Loop over xstart table

      end if ! tthandle /= 0
      call aot_table_close(L = conf, thandle = tthandle)


      ! ---------------- LOAD LENGTH1 ----------------- !
      call aot_get_val(L = conf, thandle = thandle, key = 'length1', &
        &              val = realBuffer, ErrCode = iError            )

      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*) 'FATAL Error occured, while retrieving tracker{length1}:'
        if (btest(iError, aoterr_NonExistent)) write(*,*) 'Variable nonexistent!'
        if (btest(iError, aoterr_WrongType)) write(*,*) 'Variable of wrong type!'
      end if

      tracker%length1 = realBuffer
      ! ---------------- LOAD LENGTH2 ----------------- !
      call aot_get_val(L = conf, thandle = thandle, key = 'length2', &
        &              val = realBuffer, ErrCode = iError            )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*) 'FATAL Error occured, while retrieving tracker{length2}:'
        if (btest(iError, aoterr_NonExistent)) write(*,*) 'Variable nonexistent!'
        if (btest(iError, aoterr_WrongType)) write(*,*) 'Variable of wrong type!'
      end if

      tracker%length2 = realBuffer

      ! Set debug tracker%active to ON
      tracker%active = .TRUE.

      ! Print tracker data to make sure it was loaded correctly
      call printDebugTrackerData( debugTracker = tracker, logUnit = logUnit(6) )

    else
      ! Cannot find debug tracking part of particles table.
      ! Set debug tracking to OFF
      tracker%active = .FALSE.
    end if ! thandle /= 0

  end subroutine mus_load_debugtracker

  subroutine mus_load_particle_boundaries(conf, p_thandle, bndData)
    !> configuration
    type(flu_State) :: conf
    !> Handle to parent table
    integer, intent(in) :: p_thandle
    !> Particle boundary data
    type(mus_particle_boundarydata_type), intent(inout) :: bndData
    ! -----------------------------------------------!
    integer :: boundaries_handle, domain_bnd_handle
    integer :: bnd_kind_handle
    integer :: iError
    integer :: iVal, nVals
    real(kind=rk) :: realBuffer
    character(len=labelLen) :: bnd_kind
    ! -----------------------------------------------!

    call aot_table_open(L       = conf,              &
      &                 thandle = boundaries_handle, &
      &                 parent  = p_thandle,         &
      &                 key     = 'boundaries'       )

    ! If boundaries table is available, load particle domain boundary data
    if (boundaries_handle /= 0) then

      ! ---------------- Start reading domain_bnd table --------------- !
      ! Get cartesian coordinates of prismatic domain boundaries
      call aot_table_open(L       = conf,              &
        &                 thandle = domain_bnd_handle, &
        &                 parent  = boundaries_handle, &
        &                 key     = 'domain_bnd'       )
      if (domain_bnd_handle /= 0) then
        ! get the number of values in domain_bnd, should be 6
        nVals = aot_table_length(L=conf, thandle=domain_bnd_handle)
        if(nVals /= 6) then
          write(*,*) 'FATAL Error occured, domain_bnd table is not length 6'
          write(*,*) 'use domain_bnd = [xmin, xmax, ymin, ymax, zmin, zmax]'
          call tem_abort()
        end if ! nVals /= 6

        do iVal = 1, nVals
          call aot_get_val(L = conf, thandle = domain_bnd_handle, &
            &              val = realBuffer, ErrCode = iError, &
            &              pos = iVal)
          if (btest(iError, aoterr_Fatal)) then
            write(*,*) 'FATAL Error occured, while retrieving domain_bnd'
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not existent!'
            if (btest(iError, aoterr_WrongType)) write(*,*) &
              &  'Variable has wrong type!'
          else
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not set in' &
              &  // ' config, Using default value!'

          end if ! btest(iError, aoterror_fatal)

          ! Use the values read from domain_bnd table to set pgBndData%bnd
          pgBndData%bnd(iVal) = realBuffer

        end do ! loop over values in domain_bnd table


      else
        write(*,*) 'FATAL Error occured, domain_bnd table is not defined'
        call tem_abort()
      end if ! domain_bnd_handle =/ 0

      call aot_table_close(L = conf, thandle = domain_bnd_handle)
      ! -------- Finished reading domain_bnd table --------- !

      ! Get value for boundary kinds
      call aot_table_open(L       = conf,              &
        &                 thandle = bnd_kind_handle,   &
        &                 parent  = boundaries_handle, &
        &                 key     = 'bnd_kind'         )
      if( bnd_kind_handle /= 0 ) then
        ! get the number of values in bnd_kind, should be 6
        nVals = aot_table_length(L=conf, thandle=bnd_kind_handle)
        if(nVals /= 6) then
          write(*,*) 'FATAL Error occured, domain_bnd table is not length 6'
          write(*,*) 'use e.g. bnd_kind = [wall open periodic periodic wall wall]'
          call tem_abort()
        end if ! nVals /= 6

        do iVal = 1, nVals
          call aot_get_val( L = conf, thandle = bnd_kind_handle,              &
            &               val = bnd_kind, default='wall', ErrCode = iError, &
            &               pos = iVal                                        )

          if (btest(iError, aoterr_Fatal)) then
            write(*,*) 'FATAL Error occured, while retrieving bnd_kind'
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not existent!'
            if (btest(iError, aoterr_WrongType)) write(*,*) &
              &  'Variable has wrong type!'
          else
            if (btest(iError, aoterr_NonExistent)) write(*,*) &
              &  'Variable not set in' &
              &  // ' config, Using default value!'
          end if ! btest(iError, aoterror_fatal)

          select case(bnd_kind)
            case('wall')
              pgBndData%periodicBnd(iVal) = .FALSE.
              pgBndData%wallBnd(iVal) = .TRUE.
            case('periodic')
              pgBndData%periodicBnd(iVal) = .TRUE.
              pgBndData%wallBnd(iVal) = .FALSE.
            case default
              ! Default case is no boundary treatment: particles are removed from
              ! the simulation if they exit the domain
              pgBndData%periodicBnd(iVal) = .FALSE.
              pgBndData%wallBnd(iVal) = .FALSE.
          end select

        end do

      end if ! bnd_kind_handle /= 0

      ! Indicate that we are using boundary interactions in this simulation
      pgBndData%useBnd = .TRUE.

    else
        ! If boundaries_handle == 0 then turn off useBnd to indicate we are not
        ! considering any particle-boundary interactions in this simulation.
        write(logUnit(1),*) "Particle boundary data NOT loaded, particles will be removed ", &
        & "from simulation upon exiting domain"
        pgBndData%useBnd = .FALSE.
    end if ! boundaries_handle /= 0
    call aot_table_close(L = conf, thandle = boundaries_handle)

  end subroutine mus_load_particle_boundaries


  !> Get the particle kind from the configuration
  !!
  !! Resulting string is one of:
  !! * Momentum-exchange method (kind = 'MEM')
  !! * Discrete Particle Simulations (kind = 'DPS')
  !! * One-way coupled Discrete Particle Simulations (kind = 'DPS_oneway')
  !! or 'none', if there are no particle to be modeled
  subroutine mus_load_particlekind(particle_kind, conf)
    !> particleGroup to add particles loaded from the lua script to
    character(len=*), intent(inout) :: particle_kind
    !> configuration
    type(flu_State) :: conf
    ! -----------------------------------------!
    integer :: p_thandle
    integer :: iError
    ! -----------------------------------------!
    !-- Open particle table if it exists --!
    call aot_table_open( L=conf, thandle=p_thandle, key='particles' )

    particle_kind = 'none'

    if( p_thandle > 0 ) then

      ! Set default particle kind to Momentum Exchange Method (MEM)
      particle_kind = 'DPS'

      ! -- GET PARTICLE KIND TO USE --!
      ! Available kinds:
      ! * Momentum-exchange method (kind = 'MEM')
      ! * Discrete Particle Simulations (kind = 'DPS')
      ! * One-way coupled Discrete Particle Simulations (kind = 'DPS_oneway')
      call aot_get_val( L       = conf,          &
        &               thandle = p_thandle,     &
        &               key     = 'kind',        &
        &               val     = particle_kind, &
        &               default = 'DPS',         &
        &               ErrCode = iError         )

      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*) 'FATAL Error occured, while retrieving particle kind:'
        if (btest(iError, aoterr_NonExistent)) write(*,*) 'Variable nonexistent!'
        if (btest(iError, aoterr_WrongType)) write(*,*) 'Variable of wrong type!'
      end if

    end if

    select case(trim(particle_kind))
    case('MEM', 'MEM_unittest')
      write(logUnit(1),*) 'MEM particle modelling:', &
        &                 trim(particle_kind)
    case('DPS', 'DPS_twoway', 'DPS_oneway', 'DPS_unittest')
      write(logUnit(1),*) 'DPS particle modelling:', &
        &                 trim(particle_kind)
    case('none')
      write(logUnit(1),*) 'NO particle modelling'

    case default
      write(logUnit(1),*) 'ERROR: unknown particle model: ', &
        &                 trim(particle_kind)
      write(logUnit(1),*) 'particles.kind needs to be one of:'
      write(logUnit(1),*) ' * MEM, MEM_unittest'
      write(logUnit(1),*) ' * DPS, DPS_oneway, DPS_twoway, DPS_unittest'
      write(logUnit(1),*) ' * none'
      write(logUnit(1),*) 'Aborting...'
      call tem_abort()
    end select

  end subroutine mus_load_particlekind

  function check_particle_scheme_kind_compatibility(particle_kind, scheme_kind) result(is_compatible)
    character(len=*), intent(in) :: particle_kind
    character(len=*), intent(in) :: scheme_kind
    logical :: is_compatible
    ! ----------------------------------------------- !

    is_compatible = .FALSE.
    select case( trim(particle_kind) )
    case('DPS')
      select case(trim(scheme_kind))
      case('fluid_GNS', 'fluid_incompressible_GNS')
        is_compatible = .TRUE.
      case default
        is_compatible = .FALSE.
      end select
    case('DPS_oneway','DPS_twoway')
      select case(trim(scheme_kind))
      case('fluid', 'fluid_incompressible')
        is_compatible = .TRUE.
      case default
        is_compatible = .FALSE.
      end select
    case('MEM')
      select case(trim(scheme_kind))
      case('fluid', 'fluid_incompressible')
        is_compatible = .TRUE.
      case default
        is_compatible = .FALSE.
      end select
    end select

  end function check_particle_scheme_kind_compatibility

  subroutine mus_load_particle_collisions(particleGroup, conf, p_thandle)
    !> particleGroup to add particles loaded from the lua script to
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> configuration
    type(flu_State) :: conf
    !> Handle to particle table (which should be opened before calling this routine)
    integer, intent(in) :: p_thandle
    ! ----------------------------------------------------------- !
    real(kind=rk) :: realBuffer
    integer :: iError
    ! ----------------------------------------------------------- !
    ! Get number of DEM subcycles to use per LBM time step
    call aot_get_val( L = conf, thandle = p_thandle,  &
      &               key = 'nDEMsubcycles',          &
      &               val = particleGroup%Nsubcycles, &
      &               default = 50,                   &
      &               ErrCode = iError                )
    call check_aot_error(iError, key = "nDEMsubcycles")

    ! Get particle collision time.
    ! If this is not specified or negative particles will NOT
    ! collide at all during the simulation. We specify a default value so
    ! the simulation does not abort if we don't manage to read it. In this
    ! case we will just disable collisions
    call aot_get_val( L = conf, thandle = p_thandle,   &
      &               key = 'particle_collision_time', &
      &               val = realBuffer,                &
      &               default = -1.0_rk,               &
      &               ErrCode = iError                 )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1), *) 'ERROR in reading particle_collision_time!'
      write(logUnit(1), *) 'Aborting!'
      call tem_abort()
    end if

    particleGroup%enableCollisions = (realBuffer > 0.0_rk)

    if (particleGroup%enableCollisions) then
      ! Successfully read particle collision time.
      ! Enable collisions with the specified time
      particleGroup%collision_time = realBuffer

      ! If particle collisions are enabled, get particle threshold collision gap
      ! This is the distance between particle surfaces at which we consider two
      ! particles colliding
      call aot_get_val( L = conf, thandle = p_thandle,     &
        &               key = 'particle_collision_tol',    &
        &               val = particleGroup%collision_tol, &
        &               ErrCode = iError                   )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1), *) 'ERROR in reading particle_collision_tol!'
        if (btest(iError, aoterr_NonExistent)) then
          write(logUnit(1), *) '    setting not found!'
        else
          write(logUnit(1), *) '    needs to be a real number!'
        end if
        write(logUnit(1), *) 'Aborting!'
        call tem_abort()
      end if
    end if
  end subroutine mus_load_particle_collisions

  !> Routine to load the particle data from the musubi.lua file
  !! meant to be called from within mus_load_config after the lua
  !! file has already been opened
  subroutine mus_load_particles( particleGroup, particle_kind, conf, &
    &                            chunkSize, scheme, geometry, myRank )
    !> particleGroup to add particles loaded from the lua script to
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Kind of particle modelling to use
    character(len=*), intent(in) :: particle_kind
    !> configuration
    type(flu_State) :: conf
    !> Size of the number of particles to be read in one chunk from the lua table
    integer, intent(in) :: chunkSize
    !> Scheme to determine if particles belong on this process
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry to determine if particles belong on this process
    type(mus_geom_type), intent(in) :: geometry
    !> This process rank
    integer, intent(in) :: myRank
    ! -----------------------------------------!
    !> Scheme kind to check compatibility with particle kind
    character(len=labelLen) :: schemeKind
    !> Layout kind (e.g. d2q9 or d3q19) to load correct interpolation routines
    character(len=labelLen) :: layout

    !!integer :: particleLogInterval, particleBufferSize
    !!integer :: intBuffer
    integer :: iError
    !> Handle to particle table
    integer :: p_thandle
    integer :: Nparticles

    integer :: lev
    real(kind=rk) :: dx

    logical :: predefined
    logical :: flag
    character(len=labelLen) :: blob_type
    ! -----------------------------------------!
    schemeKind = scheme%header%kind
    layout = scheme%layout%fStencil%label
    lev = geometry%tree%global%maxLevel
    dx = geometry%tree%global%BoundingCubeLength / 2**lev


    !-- Open particle table if it exists --!
    call aot_table_open( L=conf, thandle=p_thandle, key='particles' )

    if ((p_thandle > 0) .and. (trim(particle_kind) /= 'none')) then

      write(logUnit(1),*) '---- PARTICLE MUSUBI SETTINGS ----'

      ! Check to make sure particle kind is compatible with fluid kind
      flag = check_particle_scheme_kind_compatibility( &
        &      particle_kind = particle_kind,          &
        &      scheme_kind   = schemeKind              )
      if (.NOT. flag) then
        write(logUnit(1),*) 'ERROR: particle kind ', trim(particle_kind), &
          &                 ' is not compatible with schemeKind ', &
          &                 trim(schemeKind), '!'
        write(logUnit(1),*) 'ABORTING!'
        call tem_abort()
      end if

      ! ---- LOAD NUMBER OF PARTICLES ---- !
      ! This should be the total amount of particles defined in the
      ! 'particles' table in musubi.lua
      ! This is required if particle table is present.
      call aot_get_val( L = conf, thandle = p_thandle, &
        &               key = 'nParticles',            &
        &               val = nParticles,              &
        &               ErrCode = iError               )
      call check_aot_error(iError, key = 'nParticles', &
        &    event_string = 'loading particles'        )
      particleGroup%nParticles = nParticles

      !-- INITIALIZE MAXIMUM SIZE OF PARTICLE DYN ARRAYS --!
      ! This should be the greater than the total amount of particles defined
      ! in the 'particles' table in musubi.lua
      ! This is required if particle table is present.
      call aot_get_val( L = conf, thandle = p_thandle,      &
        &               key     = 'maxDynArraySize',        &
        &               val     = maxContainerSize,         &
        &               default = particleGroup%nParticles, &
        &               ErrCode = iError                    )
      call check_aot_error(iError, key = 'maxDynArraySize', &
        &    event_string = 'loading particles'             )

      ! ---- LOAD PARTICLE LOG INTERVAL ---- !
      ! Optional, default value is logging every LBM time step
      call aot_get_val( L = conf, thandle = p_thandle,           &
        &               key = 'particleLogInterval',             &
        &               val = particleGroup%particleLogInterval, &
        &               default = 1,                             &
        &               ErrCode = iError                         )
      call check_aot_error(iError, key = 'particleLogInterval', &
        &    event_string = 'loading particles'             )

      select case(particle_kind)
        case( 'MEM', 'MEM_unittest' )
          ! For MEM particles there is no default halo_distance. It must be
          ! specified by the user, one typically takes it equal to the particle
          ! diameter.
          call aot_get_val( L = conf, thandle = p_thandle,     &
            &               key = 'halo_distance',             &
            &               val = particleGroup%halo_distance, &
            &               ErrCode = iError                   )

        case('DPS', 'DPS_twoway', 'DPS_oneway', 'DPS_unittest')
          call aot_get_val( L = conf, thandle = p_thandle,     &
            &               key = 'halo_distance',             &
            &               val = particleGroup%halo_distance, &
            &               default = dx,                      &
            &               ErrCode = iError                   )
      end select
      call check_aot_error(                                          &
        &    iError, key = 'halo_distance',                          &
        &    event_string = 'reading '//trim(particle_kind)//' data' )

      call mus_load_particle_collisions( particleGroup = particleGroup, &
        &                                conf          = conf,          &
        &                                p_thandle     = p_thandle      )

      !-- INITIALIZE DEBUG TRACKERS --!
      call mus_load_debugtracker( conf = conf, p_thandle = p_thandle, tracker = debugTracker )

      !-- LOAD PARTICLE DOMAIN BOUNDARIES --!
      ! (FOR SIMPLE PRISMATIC DOMAIN CASES ONLY)
      ! By default (if domain_bnd table is not given) particles will be deleted upon leaving
      ! the domain. Prescribing domain_bnd allows us to handle interactions with periodic
      ! or wall BC's provided that the geometry is simple prismatic.

      call mus_load_particle_boundaries( conf      = conf,      &
        &                                p_thandle = p_thandle, &
        &                                bndData   = pgBndData  )

      !-- INITIALIZE PARTICLE BUFFER SIZE --!
      call aot_get_val( L       = conf,                             &
        &               thandle = p_thandle,                        &
        &               key     = 'particleBufferSize',             &
        &               val     = particleGroup%particleBufferSize, &
        &               default = 100,                              &
        &               ErrCode = iError                            )
      call check_aot_error(iError, key = 'particleBufferSize')

      !-- LOAD ACTUAL PARTICLE DATA --!
      ! Load the timing settings for the particle creator. This tells us at what
      ! iterations particles need to be created during the simulation.
      call mus_load_particle_creator_timing( conf             = conf,             &
        &                                    p_thandle        = p_thandle,        &
        &                                    particle_creator = particle_creator, &
        &                                    Nparticles       = Nparticles,       &
        &                                    chunkSize        = chunkSize,        &
        &                                    scheme           = scheme,           &
        &                                    geometry         = geometry,         &
        &                                    myRank           = myRank            )

      ! See if we are using a predefined shape to load a "blob" of particles
      ! or if we are loading each particle position individually.
      call aot_get_val( L       = conf,         &
        &               thandle = p_thandle,    &
        &               key     = 'predefined', &
        &               val     = predefined,   &
        &               default = .FALSE.,      &
        &               ErrCode = iError        )

      call aot_get_val( L       = conf,                   &
        &               thandle = p_thandle,              &
        &               key     = 'rho0_lat',             &
        &               val     = particleGroup%rho0_lat, &
        &               default = 1.0_rk,                 &
        &               ErrCode = iError                  )

      ! Initialize the particle dynamic arrays
      select case(particle_kind)
        case( 'MEM', 'MEM_unittest' )
            call init_da_particle_MEM(particleGroup%particles_MEM, 1)

        case('DPS', 'DPS_twoway', 'DPS_oneway', 'DPS_unittest')
          ! ------ LOAD INTERPOLATION STENCIL -------!
          call mus_load_particle_interpolation(            &
            &    conf         = conf,                      &
            &    p_thandle    = p_thandle,                 &
            &    layout       = layout,                    &
            &    interpolator = particleGroup%interpolator )
          write(logUnit(1),*) '-- Settings for interpolation of fluid props --'
          call printParticleInterpolator(                   &
            &    interpolator = particleGroup%interpolator, &
            &    logUnit      = logUnit(1)                  )

          call init_da_particle_DPS(particleGroup%particles_DPS, 1)
        end select


      ! Now load the particle positions into the particle creator. We can do this
      ! using either predefined shape (called a particleblob) or by loading
      ! individual positions from the lua file.
      if (predefined) then
        ! Load parameters of the cylindrical particleblob (length, radius, etc.)
        call mus_load_predefined_particleblob( conf         = conf,      &
          &                                    parent       = p_thandle, &
          &                                    blob_type    = blob_type, &
          &                                    flag         = flag       )

        select case(blob_type)
        case('cylinder')
          call print_particleblob( particleblob = particleblob, &
            &                      logUnit      = logUnit(1)    )
        case('prism')
          call print_particleblob_prism( particleblob = particleblob_prism, &
            &                            logUnit      = logUnit(1)          )
        end select

        ! Now initialize the particle creator using the data inside particleBlob
        select case(blob_type)
        case('cylinder')
          ! Generate the particle initial positions inside this blob
          call init_particle_creator_from_blob(       &
            &    particle_creator = particle_creator, &
            &    particleblob     = particleblob,     &
            &    Nparticles       = nParticles,       &
            &    scheme           = scheme,           &
            &    geometry         = geometry,         &
            &    myRank           = myRank            )
        case('prism')
          ! Generate the particle initial positions inside this blob
          call init_particle_creator_from_blob_prism(   &
            &    particle_creator = particle_creator,   &
            &    particleblob     = particleblob_prism, &
            &    Nparticles       = nParticles,         &
            &    scheme           = scheme,             &
            &    geometry         = geometry,           &
            &    myRank           = myRank              )
        end select

      else
        ! If we do NOT use a predefined shape, load the particle positions
        ! individually
        select case(particle_kind)
          case( 'MEM', 'MEM_unittest' )
            write(logUnit(1),*) 'Loading particle kind: MEM'
            ! Load the individual particle initial positions from lua script
            call load_particle_mem_creator_data(        &
              &    conf             = conf,             &
              &    parent           = p_thandle,        &
              &    particle_creator = particle_creator, &
              &    Nparticles       = Nparticles,       &
              &    chunksize        = chunkSize,        &
              &    scheme           = scheme,           &
              &    geometry         = geometry,         &
              &    myrank           = myrank            )

          case('DPS', 'DPS_twoway', 'DPS_oneway', 'DPS_unittest')
            ! Load the individual particle initial positions from lua script
            call load_particle_dps_creator_data(        &
              &    conf             = conf,             &
              &    parent           = p_thandle,        &
              &    particle_creator = particle_creator, &
              &    Nparticles       = Nparticles,       &
              &    chunksize        = chunkSize,        &
              &    scheme           = scheme,           &
              &    geometry         = geometry,         &
              &    myrank           = myrank            )

          case default
            write(logUnit(1),*) 'FATAL Error occurred, particle kind unknown:'
            call tem_abort()
        end select ! particle kind

      end if ! use predefined shape for particle 'blob' initialization

      call aot_table_close( L=conf, thandle=p_thandle )

    else
      particleGroup%nParticles = 0
      write(logUnit(1),*) 'mus_load_particles: no particle table found'
    end if ! particle table exists

    call mus_particles_print_config( particleGroup = particleGroup, &
      &                              logUnit = logUnit(1)           )

  end subroutine mus_load_particles


  subroutine load_particle_mem_data( conf, parent, particles,     &
                                   &  nparticles, chunksize  )
    type(flu_state), intent(in) :: conf

    !> handle to parent table if position, velocity tables etc are inside
    !! another table
    integer, intent(in), optional :: parent

    !> dynamic particle array to append particles read from file to
    type(dyn_particle_mem_array_type), intent(inout) :: particles

    !> total number of particles to read
    integer, intent(in) :: nparticles

    ! number of particles to read as one "chunk"
    integer, intent(in) :: chunksize
    !--------------------------------------------!
    type(mus_particle_mem_type), allocatable :: particlebuffer(:)
    integer :: istart, ichunk, iparticle
    integer :: nchunks, nchunkvals
    logical :: wasadded
    !--------------------------------------------!
    nchunks = nparticles/chunksize + 1
    istart = 1
    allocate( particlebuffer(chunksize) )

    chunkloop: do ichunk = 1, nchunks
      ! load one chunk into the particlebuffer
      ! write(logunit(1),*) 'loading chunk ', ichunk
      ! note: loading position values also sets the particle id's!
      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'position',     &
                                       & istart         = istart          )

      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'velocity',     &
                                       & istart         = istart          )

      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'force',        &
                                       & istart         = istart          )

      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'radius',       &
                                       & istart         = istart          )

      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'mass',         &
                                       & istart         = istart          )


      ! nchunkvals is the number of particles actually read from the lua file
      ! can be less than chunksize if there are fewer than chunksize particles
      ! left in lua table
      ! write(logunit(1),*) 'nchunkvals = ', nchunkvals

      !-- set particle rotational inertia --!
      do iparticle = 1, nchunkvals
        particlebuffer( iparticle )%rotinertia &
          & = 0.4 * particlebuffer( iparticle )%mass &
          &       * particlebuffer( iparticle )%radius**2

      end do

      ! append particles in particlebuffer to particlegroup
      ! write(logunit(1),*) 'loading particles: nchunkvals = ', nchunkvals
      do iparticle = 1, nchunkvals
        wasadded = .false.
        ! write(logunit(1),*) 'appending particle'

        call append_da_particle_mem( me        = particles,                       &
                              & particle  = particlebuffer(iparticle),   &
                              & length    = 1,                               &
                              & wasadded  = wasadded                         )

      end do

      istart = istart + chunksize
    end do chunkloop

    deallocate( particlebuffer )

  end subroutine load_particle_mem_data



  subroutine load_particle_mem_creator_data( conf, parent, particle_creator,     &
                                   &  nparticles, chunksize, scheme,        &
                                   & geometry, myrank  )
    type(flu_state), intent(in) :: conf

    !> handle to parent table if position, velocity tables etc are inside
    !! another table
    integer, intent(in), optional :: parent

    !> dynamic particle array to append particles read from file to
    type(mus_particle_creator_type), intent(inout) :: particle_creator

    !> total number of particles to read
    integer, intent(in) :: nparticles

    ! number of particles to read as one "chunk"
    integer, intent(in) :: chunksize

    !> scheme
    type(mus_scheme_type), intent(in) :: scheme

    !> geometry
    type(mus_geom_type), intent(in) :: geometry

    !> this mpi process rank
    integer, intent(in) :: myrank
    !--------------------------------------------!
    type(mus_particle_mem_type), allocatable :: particlebuffer(:)
    integer :: istart, ichunk, iparticle, kparticle
    integer :: nchunks, nchunkvals
    logical :: wasadded, islocal
    !--------------------------------------------!
    nchunks = nparticles/chunksize + 1
    istart = 1
    kparticle = 1
    allocate( particlebuffer(chunksize) )
    particle_creator%global_nparticles = 0
    particle_creator%n_times_called = 0

    chunkloop: do ichunk = 1, nchunks
      ! load one chunk into the particlebuffer
      write(logunit(1),*) 'loading chunk ', ichunk
      ! note: loading position values also sets the particle id's!
      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'position',     &
                                       & istart         = istart          )

      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'velocity',     &
                                       & istart         = istart          )

      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'force',        &
                                       & istart         = istart          )

      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'radius',       &
                                       & istart         = istart          )

      call mus_load_particle_mem_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'mass',         &
                                       & istart         = istart          )


      ! nchunkvals is the number of particles actually read from the lua file
      ! can be less than chunksize if there are fewer than chunksize particles
      ! left in lua table
      ! write(logunit(1),*) 'nchunkvals = ', nchunkvals

      !-- set particle rotational inertia --!
      do iparticle = 1, nchunkvals
        particlebuffer( iparticle )%rotinertia &
          & = 0.4 * particlebuffer( iparticle )%mass &
          &       * particlebuffer( iparticle )%radius**2

      end do

      ! append particles in particlebuffer to particlegroup
      write(logunit(1),*) 'loading particles: nchunkvals = ', nchunkvals
      do iparticle = 1, nchunkvals
        wasadded = .false.
        ! write(logunit(1),*) 'appending particle'

        ! check if this particle position is local to my process
        islocal = positionlocalonmyrank( pos = particlebuffer(iparticle)%pos(1:3), &
                                       & geometry = geometry,                      &
                                       & scheme = scheme,                          &
                                       & myrank = myrank                           )

        if(islocal) then
          call append( me = particle_creator%position, val = particlebuffer(iparticle)%pos(1:6) )
          call append( me = particle_creator%velocity, val = particlebuffer(iparticle)%vel(1:6) )
          call append( me = particle_creator%force, val = particlebuffer(iparticle)%fext(1:6) )
          call append( me = particle_creator%radius, val = particlebuffer(iparticle)%radius )
          call append( me = particle_creator%mass, val = particlebuffer(iparticle)%mass )
          call append( me = particle_creator%idoffset, val = kparticle )
        end if

        ! increment global number of particles
        particle_creator%global_nparticles = particle_creator%global_nparticles + 1

        ! kparticle keeps track of the number of particles that have been added to the
        ! particle creators on all processes. so after we've loaded all the particles
        ! it will be equal to the number of particles specified in the lua script.
        kparticle = kparticle + 1

      end do

      istart = istart + chunksize
    end do chunkloop

    particle_creator%nparticles = particle_creator%radius%nvals
    deallocate( particlebuffer )

  end subroutine load_particle_mem_creator_data


  subroutine mus_load_particle_mem_data_chunk( conf, parent, particlebuffer, &
                                         & nchunkvals, key, istart       )

    type(flu_state), intent(in) :: conf

    !> handle to parent table if position, velocity tables etc are inside
    !! another table
    integer, intent(in), optional :: parent

    ! buffer to hold the data read from chunk
    type(mus_particle_mem_type), allocatable, intent(inout) :: particlebuffer(:)

    ! number of values actually read from table
    ! if entire buffer is filled, nchunkvals = size(particlebuffer)
    integer, intent(out) :: nchunkvals

    ! key = 'position', 'velocity', 'radius' or 'mass'
    character(len=*), intent(in) :: key
    ! index to start reading the lua array at
    integer, intent(in) :: istart

    !--------------------------------------------!
    integer :: ierror
    integer :: verr(6)
    real(kind=rk) :: vecbuffer(6)
    real(kind=rk) :: realbuffer

    integer :: nchunk

    integer :: nvals
    integer :: thandle
    integer :: iparticle, ichunk
    !--------------------------------------------!
    ierror = 0
    nchunkvals = 0

    if( .not.( allocated( particlebuffer ) ) ) then
      write(logunit(1),*) 'error mus_load_particle_data_chunk:'
      call tem_abort( 'particlebuffer not allocated')
    end if

    ! number of elements to read = size of particlebuffer
    nchunk = size(particlebuffer)

    ! determine whether we're loading a vector or scalar quantity
    if( key == 'position' .or. key == 'velocity' .or. key == 'force' ) then
      ! initialize vector particle quantity
      ! call aot_table_open(l = conf, thandle = thandle, key = key)
      call aot_table_open(l       = conf,    &
        &                 thandle = thandle, &
        &                 parent  = parent,  &
        &                 key     = key      )
      if (thandle /= 0) then
        ! get the number of position vectors, should be equal to nparticles
        nvals = aot_table_length(l=conf, thandle=thandle)

        ! loop over all the position vals
        ichunk = 1
        do iparticle=istart, istart + nchunk - 1

          ! check to make sure we haven't reached the end of the lua table
          if( iparticle > nvals ) exit

          ! get the entire position vector
          call aot_get_val( l = conf, thandle = thandle, &
                          & pos = iparticle, val = vecbuffer, &
                          & errcode = verr, &
                          & default = [0.0_double_k, 0.0_double_k, &
                          &            0.0_double_k, 0.0_double_k, &
                          &            0.0_double_k, 0.0_double_k ])

          if (btest(ierror, aoterr_fatal)) then
            write(*,*) 'fatal error occured, while retrieving particle pos'
            if (btest(ierror, aoterr_nonexistent)) write(*,*) &
              &  'variable not existent!'
            if (btest(ierror, aoterr_wrongtype)) write(*,*) &
              &  'variable has wrong type!'
          else
            if (btest(ierror, aoterr_nonexistent)) write(*,*) &
              &  'variable not set in' &
              &  // ' config, using default value!'
          end if

          ! copy property to the particlebuffer
          if( key == 'position' ) then
            particlebuffer(ichunk)%pos(1:6) = vecbuffer(1:6)
            particlebuffer(ichunk)%particleid = iparticle
          else if( key == 'velocity' ) then
            particlebuffer(ichunk)%vel(1:6) = vecbuffer(1:6)
          else if( key == 'force' ) then
            particlebuffer(ichunk)%fext(1:6) = vecbuffer(1:6)
          end if

          ichunk = ichunk + 1

          ! increment nchunkvals to signify we've read another chunk
          nchunkvals = nchunkvals + 1
        end do

      end if
      call aot_table_close(l = conf, thandle = thandle)


    else if ( key == 'radius' .or. key == 'mass' ) then
      ! initialize scalar particle quantity
      ! call aot_table_open(l = conf, thandle = thandle, key = key)
      call aot_table_open(l       = conf,    &
        &                 thandle = thandle, &
        &                 parent  = parent,  &
        &                 key     = key      )
      if (thandle /= 0) then
        ! get the number of position vectors, should be equal to nparticles
        nvals = aot_table_length(l=conf, thandle=thandle)

        ichunk = 1
        do iparticle = istart, istart + nchunk - 1

          ! check to make sure we haven't reached the end of the lua table
          if( iparticle > nvals ) exit

          call aot_get_val(l = conf, thandle = thandle, &
            &              val = realbuffer, errcode = ierror, &
            &              pos = iparticle)
          if (btest(ierror, aoterr_fatal)) then
            write(*,*) 'fatal error occured, while retrieving radius'
            if (btest(ierror, aoterr_nonexistent)) write(*,*) &
              &  'variable not existent!'
            if (btest(ierror, aoterr_wrongtype)) write(*,*) &
              &  'variable has wrong type!'
          else
            if (btest(ierror, aoterr_nonexistent)) write(*,*) &
              &  'variable not set in' &
              &  // ' config, using default value!'
          end if

          ! copy property to the particlebuffer
          if( key == 'radius' ) then
            particlebuffer(ichunk)%radius = realbuffer
          else if ( key == 'mass' ) then
            particlebuffer(ichunk)%mass = realbuffer
          end if

          ichunk = ichunk + 1
          nchunkvals = nchunkvals + 1

        end do
      end if
      call aot_table_close(l = conf, thandle = thandle)

    else
      write(logunit(1),'(a)') &
        & 'error mus_load_particle_data: unrecognized property'
    end if ! particle property


  end subroutine mus_load_particle_mem_data_chunk




  subroutine load_particle_dps_data( conf, parent, particles,     &
                                   &  nparticles, chunksize  )
    type(flu_state), intent(in) :: conf

    !> handle to parent table if position, velocity tables etc are inside
    !! another table
    integer, intent(in), optional :: parent

    !> dynamic particle array to append particles read from file to
    type(dyn_particle_dps_array_type), intent(inout) :: particles

    !> total number of particles to read
    integer, intent(in) :: nparticles

    ! number of particles to read as one "chunk"
    integer, intent(in) :: chunksize
    !--------------------------------------------!
    type(mus_particle_dps_type), allocatable :: particlebuffer(:)
    integer :: istart, ichunk, iparticle
    integer :: nchunks, nchunkvals
    logical :: wasadded
    !--------------------------------------------!
    nchunks = nparticles/chunksize + 1
    istart = 1
    allocate( particlebuffer(chunksize) )

    chunkloop: do ichunk = 1, nchunks
      ! load one chunk into the particlebuffer
      ! write(logunit(1),*) 'loading chunk ', ichunk
      ! note: loading position values also sets the particle id's!
      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'position',     &
                                       & istart         = istart          )

      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'velocity',     &
                                       & istart         = istart          )

      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'force',        &
                                       & istart         = istart          )

      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'radius',       &
                                       & istart         = istart          )

      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'mass',         &
                                       & istart         = istart          )


      ! nchunkvals is the number of particles actually read from the lua file
      ! can be less than chunksize if there are fewer than chunksize particles
      ! left in lua table
      ! write(logunit(1),*) 'nchunkvals = ', nchunkvals

      !-- set particle rotational inertia --!
      do iparticle = 1, nchunkvals
        particlebuffer( iparticle )%rotinertia &
          & = 0.4 * particlebuffer( iparticle )%mass &
          &       * particlebuffer( iparticle )%radius**2

      end do

      ! append particles in particlebuffer to particlegroup
      ! write(logunit(1),*) 'loading particles: nchunkvals = ', nchunkvals
      do iparticle = 1, nchunkvals
        wasadded = .false.
        ! write(logunit(1),*) 'appending particle'

        call append_da_particle_dps( me        = particles,                       &
                              & particle  = particlebuffer(iparticle),   &
                              & length    = 1,                               &
                              & wasadded  = wasadded                         )

      end do

      istart = istart + chunksize
    end do chunkloop

    deallocate( particlebuffer )

  end subroutine load_particle_dps_data



  subroutine load_particle_dps_creator_data( conf, parent, particle_creator,     &
                                   &  nparticles, chunksize, scheme,        &
                                   & geometry, myrank  )
    type(flu_state), intent(in) :: conf

    !> handle to parent table if position, velocity tables etc are inside
    !! another table
    integer, intent(in), optional :: parent

    !> dynamic particle array to append particles read from file to
    type(mus_particle_creator_type), intent(inout) :: particle_creator

    !> total number of particles to read
    integer, intent(in) :: nparticles

    ! number of particles to read as one "chunk"
    integer, intent(in) :: chunksize

    !> scheme
    type(mus_scheme_type), intent(in) :: scheme

    !> geometry
    type(mus_geom_type), intent(in) :: geometry

    !> this mpi process rank
    integer, intent(in) :: myrank
    !--------------------------------------------!
    type(mus_particle_dps_type), allocatable :: particlebuffer(:)
    integer :: istart, ichunk, iparticle, kparticle
    integer :: nchunks, nchunkvals
    logical :: wasadded, islocal
    !--------------------------------------------!
    nchunks = nparticles/chunksize + 1
    istart = 1
    kparticle = 1
    allocate( particlebuffer(chunksize) )
    particle_creator%global_nparticles = 0
    particle_creator%n_times_called = 0

    chunkloop: do ichunk = 1, nchunks
      ! load one chunk into the particlebuffer
      write(logunit(1),*) 'loading chunk ', ichunk
      ! note: loading position values also sets the particle id's!
      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'position',     &
                                       & istart         = istart          )

      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'velocity',     &
                                       & istart         = istart          )

      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'force',        &
                                       & istart         = istart          )

      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'radius',       &
                                       & istart         = istart          )

      call mus_load_particle_dps_data_chunk( conf           = conf,           &
                                       & parent         = parent,      &
                                       & particlebuffer = particlebuffer, &
                                       & nchunkvals     = nchunkvals,     &
                                       & key            = 'mass',         &
                                       & istart         = istart          )


      ! nchunkvals is the number of particles actually read from the lua file
      ! can be less than chunksize if there are fewer than chunksize particles
      ! left in lua table
      ! write(logunit(1),*) 'nchunkvals = ', nchunkvals

      !-- set particle rotational inertia --!
      do iparticle = 1, nchunkvals
        particlebuffer( iparticle )%rotinertia &
          & = 0.4 * particlebuffer( iparticle )%mass &
          &       * particlebuffer( iparticle )%radius**2

      end do

      ! append particles in particlebuffer to particlegroup
      write(logunit(1),*) 'loading particles: nchunkvals = ', nchunkvals
      do iparticle = 1, nchunkvals
        wasadded = .false.
        ! write(logunit(1),*) 'appending particle'

        ! check if this particle position is local to my process
        islocal = positionlocalonmyrank( pos = particlebuffer(iparticle)%pos(1:3), &
                                       & geometry = geometry,                      &
                                       & scheme = scheme,                          &
                                       & myrank = myrank                           )

        if(islocal) then
          call append( me = particle_creator%position, val = particlebuffer(iparticle)%pos(1:6) )
          call append( me = particle_creator%velocity, val = particlebuffer(iparticle)%vel(1:6) )
          call append( me = particle_creator%force, val = particlebuffer(iparticle)%fext(1:6) )
          call append( me = particle_creator%radius, val = particlebuffer(iparticle)%radius )
          call append( me = particle_creator%mass, val = particlebuffer(iparticle)%mass )
          call append( me = particle_creator%idoffset, val = kparticle )
        end if

        ! increment global number of particles
        particle_creator%global_nparticles = particle_creator%global_nparticles + 1

        ! kparticle keeps track of the number of particles that have been added to the
        ! particle creators on all processes. so after we've loaded all the particles
        ! it will be equal to the number of particles specified in the lua script.
        kparticle = kparticle + 1

      end do

      istart = istart + chunksize
    end do chunkloop

    particle_creator%nparticles = particle_creator%radius%nvals
    deallocate( particlebuffer )

  end subroutine load_particle_dps_creator_data


  subroutine mus_load_particle_dps_data_chunk( conf, parent, particlebuffer, &
                                         & nchunkvals, key, istart       )

    type(flu_state), intent(in) :: conf

    !> handle to parent table if position, velocity tables etc are inside
    !! another table
    integer, intent(in), optional :: parent

    ! buffer to hold the data read from chunk
    type(mus_particle_dps_type), allocatable, intent(inout) :: particlebuffer(:)

    ! number of values actually read from table
    ! if entire buffer is filled, nchunkvals = size(particlebuffer)
    integer, intent(out) :: nchunkvals

    ! key = 'position', 'velocity', 'radius' or 'mass'
    character(len=*), intent(in) :: key
    ! index to start reading the lua array at
    integer, intent(in) :: istart

    !--------------------------------------------!
    integer :: ierror
    integer :: verr(6)
    real(kind=rk) :: vecbuffer(6)
    real(kind=rk) :: realbuffer

    integer :: nchunk

    integer :: nvals
    integer :: thandle
    integer :: iparticle, ichunk
    !--------------------------------------------!
    ierror = 0
    nchunkvals = 0

    if( .not.( allocated( particlebuffer ) ) ) then
      write(logunit(1),*) 'error mus_load_particle_data_chunk:'
      call tem_abort( 'particlebuffer not allocated')
    end if

    ! number of elements to read = size of particlebuffer
    nchunk = size(particlebuffer)

    ! determine whether we're loading a vector or scalar quantity
    if( key == 'position' .or. key == 'velocity' .or. key == 'force' ) then
      ! initialize vector particle quantity
      ! call aot_table_open(l = conf, thandle = thandle, key = key)
      call aot_table_open(l       = conf,    &
        &                 thandle = thandle, &
        &                 parent  = parent,  &
        &                 key     = key      )
      if (thandle /= 0) then
        ! get the number of position vectors, should be equal to nparticles
        nvals = aot_table_length(l=conf, thandle=thandle)

        ! loop over all the position vals
        ichunk = 1
        do iparticle=istart, istart + nchunk - 1

          ! check to make sure we haven't reached the end of the lua table
          if( iparticle > nvals ) exit

          ! get the entire position vector
          call aot_get_val( l = conf, thandle = thandle, &
                          & pos = iparticle, val = vecbuffer, &
                          & errcode = verr, &
                          & default = [0.0_double_k, 0.0_double_k, &
                          &            0.0_double_k, 0.0_double_k, &
                          &            0.0_double_k, 0.0_double_k ])

          if (btest(ierror, aoterr_fatal)) then
            write(*,*) 'fatal error occured, while retrieving particle pos'
            if (btest(ierror, aoterr_nonexistent)) write(*,*) &
              &  'variable not existent!'
            if (btest(ierror, aoterr_wrongtype)) write(*,*) &
              &  'variable has wrong type!'
          else
            if (btest(ierror, aoterr_nonexistent)) write(*,*) &
              &  'variable not set in' &
              &  // ' config, using default value!'
          end if

          ! copy property to the particlebuffer
          if( key == 'position' ) then
            particlebuffer(ichunk)%pos(1:6) = vecbuffer(1:6)
            particlebuffer(ichunk)%particleid = iparticle
          else if( key == 'velocity' ) then
            particlebuffer(ichunk)%vel(1:6) = vecbuffer(1:6)
          else if( key == 'force' ) then
            particlebuffer(ichunk)%fext(1:6) = vecbuffer(1:6)
          end if

          ichunk = ichunk + 1

          ! increment nchunkvals to signify we've read another chunk
          nchunkvals = nchunkvals + 1
        end do

      end if
      call aot_table_close(l = conf, thandle = thandle)


    else if ( key == 'radius' .or. key == 'mass' ) then
      ! initialize scalar particle quantity
      ! call aot_table_open(l = conf, thandle = thandle, key = key)
      call aot_table_open(l       = conf,    &
        &                 thandle = thandle, &
        &                 parent  = parent,  &
        &                 key     = key      )
      if (thandle /= 0) then
        ! get the number of position vectors, should be equal to nparticles
        nvals = aot_table_length(l=conf, thandle=thandle)

        ichunk = 1
        do iparticle = istart, istart + nchunk - 1

          ! check to make sure we haven't reached the end of the lua table
          if( iparticle > nvals ) exit

          call aot_get_val(l = conf, thandle = thandle, &
            &              val = realbuffer, errcode = ierror, &
            &              pos = iparticle)
          if (btest(ierror, aoterr_fatal)) then
            write(*,*) 'fatal error occured, while retrieving radius'
            if (btest(ierror, aoterr_nonexistent)) write(*,*) &
              &  'variable not existent!'
            if (btest(ierror, aoterr_wrongtype)) write(*,*) &
              &  'variable has wrong type!'
          else
            if (btest(ierror, aoterr_nonexistent)) write(*,*) &
              &  'variable not set in' &
              &  // ' config, using default value!'
          end if

          ! copy property to the particlebuffer
          if( key == 'radius' ) then
            particlebuffer(ichunk)%radius = realbuffer
          else if ( key == 'mass' ) then
            particlebuffer(ichunk)%mass = realbuffer
          end if

          ichunk = ichunk + 1
          nchunkvals = nchunkvals + 1

        end do
      end if
      call aot_table_close(l = conf, thandle = thandle)

    else
      write(logunit(1),'(a)') &
        & 'error mus_load_particle_data: unrecognized property'
    end if ! particle property


  end subroutine mus_load_particle_dps_data_chunk



  subroutine mus_finalize_particleGroup( particleGroup )

    type(mus_particle_group_type), intent(inout) :: particleGroup
    !--------------------------------------------!

    ! Destroy particleGroup's dynamic array of particles
    call destroy_da_particle_MEM(particleGroup%particles_MEM)

  end subroutine mus_finalize_particleGroup

  subroutine mus_particles_print_config(particleGroup, logUnit)
    type(mus_particle_group_type), intent(in) :: particleGroup
    integer, intent(in) :: logUnit
    ! ------------------------------------ !
    if (particleGroup%nParticles > 0) then
      write(logUnit,'(A)') '---- PARTICLEGROUP SETTINGS ----'
      write(logUnit,'(A,L2)')    'enableCollisions = ', particleGroup%enableCollisions
      write(logUnit,'(A,E17.9)') 'collision_tol    = ', particleGroup%collision_tol
      write(logUnit,'(A,E17.9)') 'collision_time   = ', particleGroup%collision_time
      write(logUnit,'(A,I12)') 'Number of DEM subcycles = ', particleGroup%Nsubcycles
      write(logUnit,'(A,I12)') 'Particle buffer size    = ', particleGroup%particleBufferSize
      write(logUnit,'(A,I12)') 'Particle log interval   = ', particleGroup%particleLogInterval
    end if
    ! ------------------------------------ !
  end subroutine mus_particles_print_config

  subroutine mus_load_predefined_particleblob( conf, parent, blob_type, flag)
    !> configuration
    type(flu_State) :: conf
    !> Handle to parent table in which the "shape" table exists
    integer, intent(in) :: parent
    !> Output string containing the blob type
    character(len=LabelLen), intent(inout) :: blob_type
    !> Output flag set to TRUE when reading was successful, FALSE if not
    logical, intent(out) :: flag
    ! ---------------------------------------------------- !
    integer :: shape_thandle, iError
    integer :: vError(6)
    real(kind=rk) :: default_force(6)
    real(kind=rk) :: default_vel(6)
    real(kind=rk) :: force(6)
    real(kind=rk) :: velocity(6)
    real(kind=rk) :: radius
    real(kind=rk) :: mass
    logical :: has_init_velocity
    ! ---------------------------------------------------- !
    default_force = 0.0_rk
    default_vel = 0.0_rk

    call aot_table_open(L       = conf,          &
      &                 thandle = shape_thandle, &
      &                 parent  = parent,        &
      &                 key     = 'shape'        )

    if (shape_thandle /= 0) then
      ! Load shape kind
      call aot_get_val( L       = conf,          &
        &               thandle = shape_thandle, &
        &               key     = 'kind',        &
        &               val     = blob_type,     &
        &               ErrCode = iError         )
      call check_aot_error( iError, key = 'kind',                   &
        &                   event_string = 'getting particle shape' )

      ! Load the particle properties from the particle table (handle: parent).
      ! These are the same for each particle in the "blob"
      call aot_get_val( L       = conf,                      &
        &               thandle = parent,                    &
        &               key     = 'velocity',                &
        &               val     = particleblob%particle_vel, &
        &               default = default_vel,               &
        &               ErrCode = vError                     )
      if (any(btest(vError, aoterr_Fatal))) then
        write(logUnit(0), *) "ERROR reading the velocity for " &
          &                  // "particleblob cylinder"
        call tem_abort()
      end if
      has_init_velocity = (sum(abs(velocity)) > 8*tiny(1.0_rk))

      call aot_get_val( L       = conf,                        &
        &               thandle = parent,                      &
        &               key     = 'force',                     &
        &               val     = particleblob%particle_force, &
        &               default = default_vel,                 &
        &               ErrCode = vError                       )
      if (any(btest(vError, aoterr_Fatal))) then
        write(logUnit(0), *) "ERROR reading the force for " &
          &                  // "particleblob cylinder"
        call tem_abort()
      end if

      call aot_get_val( L       = conf,                         &
        &               thandle = parent,                       &
        &               key     = 'radius',                     &
        &               val     = particleblob%particle_radius, &
        &               ErrCode = iError                        )
      call check_aot_error(iError, key = 'radius', &
        &                  event_string = 'reading particleblob')

      call aot_get_val( L       = conf,                       &
        &               thandle = parent,                     &
        &               key     = 'mass',                     &
        &               val     = particleblob%particle_mass, &
        &               ErrCode = iError                      )
      call check_aot_error(iError, key = 'mass', &
        &                  event_string = 'reading particleblob')

      select case(trim(blob_type))
      case('cylinder')
        call mus_load_particleblob_cylinder( particleblob = particleblob, &
          &                                  conf         = conf,         &
          &                                  parent       = shape_thandle )
        particleblob%particle_vel = velocity
        particleblob%particle_force = force
        particleblob%init_particles_to_fluid_vel = has_init_velocity
        particleblob%particle_radius = radius
        particleblob%particle_mass = mass

      case('prism')
        call mus_load_particleblob_prism( particleblob = particleblob_prism, &
          &                               conf         = conf,               &
          &                               parent       = shape_thandle       )
        particleblob_prism%particle_vel = velocity
        particleblob_prism%particle_force = force
        particleblob_prism%init_particles_to_fluid_vel = has_init_velocity
        particleblob_prism%particle_radius = radius
        particleblob_prism%particle_mass = mass

      case default
        write(logUnit(1),*) "ERROR: Predefined particle shape kind not recognized, aborting!"
        write(logUnit(1),*) "       Has to be one of:"
        write(logUnit(1),*) "       * cylinder"
        write(logUnit(1),*) "       * prism"
        call tem_abort()
      end select

      flag = .TRUE.
    else
      flag = .FALSE.
    end if
    call aot_table_close(L = conf, thandle = shape_thandle)

  end subroutine mus_load_predefined_particleblob

  subroutine mus_load_particle_distribution(distribution, conf, parent)
    type(mus_particle_blob_prob_type), intent(inout) :: distribution
    !> configuration
    type(flu_State) :: conf
    !> Handle to parent table in which the "shape" table exists
    integer, intent(in) :: parent
    ! --------------------------------------- !
    integer :: iError, pError(2)
    ! --------------------------------------- !
    call aot_get_val( L       = conf,              &
      &               thandle = parent,            &
      &               key     = 'distribution',    &
      &               val     = distribution%kind, &
      &               default = 'uniform',         &
      &               ErrCode = iError             )
    call check_aot_error(iError, key = "distribution",   &
      &    event_string = "reading particle distribution")

    if ( trim(particleblob%distribution%kind) == 'gaussian' ) then
      ! Load the random seed to be used for generating random positions
      call aot_get_val( L       = conf,              &
        &               thandle = parent,            &
        &               key     = 'seed',            &
        &               val     = distribution%seed, &
        &               default = 10,                &
        &               ErrCode = iError             )
      call check_aot_error(iError, key = "seed",   &
        &    event_string = "reading gaussian distribution")

      ! Load the mean and standard deviation of the distribution
      call aot_get_val( L       = conf,            &
        &               thandle = parent,          &
        &               key     = 'mu',            &
        &               val     = distribution%mu, &
        &               ErrCode = pError           )
      if (any(btest(pError, aoterr_Fatal))) then
        write(logUnit(0),*) "ERROR reading mu for gaussian distribution"
      end if

      call aot_get_val( L       = conf,               &
        &               thandle = parent,             &
        &               key     = 'sigma',            &
        &               val     = distribution%sigma, &
        &               ErrCode = pError              )
      if (any(btest(pError, aoterr_Fatal))) then
        write(logUnit(0),*) "ERROR reading sigma for gaussian distribution"
      end if

    end if

  end subroutine mus_load_particle_distribution

  subroutine mus_load_particleblob_cylinder(particleblob, conf, parent)
    type(mus_particle_blob_cylinder_type), intent(inout) :: particleblob
    !> configuration
    type(flu_State) :: conf
    !> Handle to parent table in which the "shape" table exists
    integer, intent(in) :: parent
    ! --------------------------------------- !
    integer :: iError, vError(3)
    logical :: failure
    ! --------------------------------------- !
    failure = .false.
    ! load origin
    call aot_get_val( L       = conf,                &
      &               thandle = parent,              &
      &               key     = 'origin',            &
      &               val     = particleblob%origin, &
      &               ErrCode = vError               )
    failure = failure .or. any(btest(vError, aoterr_Fatal))

    ! load cylinder axis
    call aot_get_val( L       = conf,             &
      &               thandle = parent,           &
      &               key     = 'vec',            &
      &               val     = particleblob%vec, &
      &               ErrCode = vError            )
    failure = failure .or. any(btest(vError, aoterr_Fatal))

    ! load cylinder radius
    call aot_get_val( L       = conf,                &
      &               thandle = parent,              &
      &               key     = 'cylinder_radius',   &
      &               val     = particleblob%radius, &
      &               ErrCode = iError               )
    failure = failure .or. btest(iError, aoterr_Fatal)

    call mus_load_particle_distribution(particleblob%distribution, conf, parent)

    if (failure) then
      write(logUnit(1),*) "ERROR mus_load_particleblob_cylinder: could not read cylinder params"
      call tem_abort()
    end if

    ! --------------------------------------- !
  end subroutine mus_load_particleblob_cylinder

  subroutine mus_load_particleblob_prism(particleblob, conf, parent)
    type(mus_particle_blob_prism_type), intent(inout) :: particleblob
    !> configuration
    type(flu_State) :: conf
    !> Handle to parent table in which the "shape" table exists
    integer, intent(in) :: parent
    ! --------------------------------------- !
    integer :: iError, vError(3)
    logical :: failure
    ! --------------------------------------- !
    failure = .false.
    ! load origin
    call aot_get_val( L       = conf,                &
      &               thandle = parent,              &
      &               key     = 'origin',            &
      &               val     = particleblob%origin, &
      &               ErrCode = vError               )
    failure = failure .or. any(btest(vError, aoterr_Fatal))

    ! load x-direction axis
    call aot_get_val( L       = conf,               &
      &               thandle = parent,             &
      &               key     = 'vec_x',            &
      &               val     = particleblob%vec_x, &
      &               ErrCode = vError              )
    failure = failure .or. any(btest(vError, aoterr_Fatal))

    ! load y-direction axis
    call aot_get_val( L       = conf,               &
      &               thandle = parent,             &
      &               key     = 'vec_y',            &
      &               val     = particleblob%vec_y, &
      &               ErrCode = vError              )
    failure = failure .or. any(btest(vError, aoterr_Fatal))

    ! load z-direction axis
    call aot_get_val( L       = conf,               &
      &               thandle = parent,             &
      &               key     = 'vec_z',            &
      &               val     = particleblob%vec_z, &
      &               ErrCode = vError              )
    failure = failure .or. any(btest(vError, aoterr_Fatal))

    call aot_get_val( L       = conf,            &
      &               thandle = parent,          &
      &               key     = 'nx',            &
      &               val     = particleblob%nx, &
      &               default = -1,              &
      &               ErrCode = iError           )
    failure = failure .or. btest(iError, aoterr_Fatal)

    call aot_get_val( L       = conf,            &
      &               thandle = parent,          &
      &               key     = 'ny',            &
      &               val     = particleblob%ny, &
      &               default = -1,              &
      &               ErrCode = iError           )
    failure = failure .or. btest(iError, aoterr_Fatal)

    call aot_get_val( L       = conf,            &
      &               thandle = parent,          &
      &               key     = 'nz',            &
      &               val     = particleblob%nz, &
      &               default = -1,              &
      &               ErrCode = iError           )
    failure = failure .or. btest(iError, aoterr_Fatal)

    ! Optional: load the probability distribution of particles inside the blob
    call mus_load_particle_distribution(particleblob%distribution, conf, parent)

    if (failure) then
      write(logUnit(1),*) "ERROR mus_load_particleblob_prism: could not read prism params"
      call tem_abort()
    end if

    ! --------------------------------------- !
  end subroutine mus_load_particleblob_prism

end module mus_particle_config_module
