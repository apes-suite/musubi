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
!> mus_particle_creator_module contains data types and routines for the particle 
!! creator object which is used to initialize new particles at specified locations 
!! at specified time steps.

module mus_particle_creator_module
  use mpi
  use env_module,                       only: rk, newunit
  use tem_param_module,                 only: PI
  use tem_aux_module,                   only: tem_abort
  use tem_logging_module,               only: logUnit
  use tem_geometry_module,              only: tem_CoordOfReal
  use tem_grow_array_module,            only: init, append, destroy, empty, &
    &                                         grw_realarray_type,           &
    &                                         grw_intarray_type,            &
    &                                         grw_real2darray_type
  use mus_geom_module,                  only: mus_geom_type
  use mus_scheme_type_module,           only: mus_scheme_type
  use mus_param_module,                 only: mus_param_type
  use mus_particle_logging_type_module, only: pgDebugLog
  use mus_particle_aux_module,          only: positionLocalOnMyRank
  use mus_particle_comm_module,         only: exchangeNewParticles_DPS, &
    &                                         exchangeNewParticles_MEM
  use mus_particle_type_module, only:         &
    &   mus_particle_group_type,              &
    &   mus_particle_DPS_type,                &
    &   mus_particle_MEM_type,                &
    &   allocateProcessMasks,                 &
    &   printParticleGroup2_DPS,              &
    &   printParticleGroup2_MEM,              &
    &   remove_particle_from_da_particle_DPS, &
    &   remove_particle_from_da_particle_MEM, &
    &   append_da_particle_DPS,               &
    &   append_da_particle_MEM
  use mus_particle_DPS_module,  only: initParticle_DPS, &
    &                                 interpolateFluidProps
  use mus_particle_MEM_module,  only: updateCoordOfOrigin, &
    &                                 initParticle_MEM
  use mus_particle_blob_module, only: mus_particle_blob_cylinder_type, &
    &                                 mus_particle_blob_prism_type, &
    &                                 fill_cylinder, &
    &                                 fill_prism, &
    &                                 set_random_seed, &
    &                                 init_particle_blob_prob, &
    &                                 pick_random_position, &
    &                                 print_particleblob_prob, &
    &                                 print_positions, &
    &                                 translate_positions, &
    &                                 rotate_positions
  
  implicit none
  
  !> Data type used to create particles with certain properties at specified time 
  !! steps
  type mus_particle_creator_type
    ! ---------- PARTICLE DATA -----------!
    !> Number of particles in this particle creator object.
    !! This is the amount of particles that will be created at each interval.
    integer :: Nparticles
    !> Dynamic array containing positions of particles to be created at each iter
    !! First index corresponds to a certain particle, second to x, y, z, rx, ry, rz coordinate
    type(grw_real2darray_type) :: position
    !> Dynamic array containing velocities of particles to be created at each iter
    type(grw_real2darray_type) :: velocity
    !> Dynamic array containing forces of particles to be created at each iter
    type(grw_real2darray_type) :: force
    !> Dynamic array containing radii of particles to be created at each iter
    type(grw_realarray_type) :: radius
    !> Dynamic array containing mass of particles to be created at each iter
    type(grw_realarray_type) :: mass
    !> Dynamic array containing the particleID offsets for each particle
    !! These are used to create unique (across all processes) ID's for each particle
    type(grw_intarray_type) :: IDoffset
    
    !> Logical to indicate whether to initialize particles to the local fluid velocity
    !! If this is set to TRUE, particles will be initialized with the local fluid velocity 
    !! when they are created. In this case the particle_creator%velocities array is not used.
    logical :: init_particle_to_fluid_vel
  
  
    ! ---------- GLOBAL DATA ------------!
    ! We need this to create unique particleID's for each particle we create
  
    !> Number of particles in the particle creators on all processes
    integer :: global_Nparticles
  
    !> Times that the particleCreator has been called
    !! This is the same on all processes
    integer :: N_times_called = 0
  
    ! ---------- TEMPORAL DATA start-----------!
    !> Iteration number to start creating particles
    integer :: iter_start
    !> Iteration number to stop creating particles
    integer :: iter_end
    !> Create particles every this many time steps
    integer :: iter_interval
  end type mus_particle_creator_type
  
  type(mus_particle_creator_type), save :: particle_creator
  

contains
  

  subroutine init_particle_creator(me, Nparticles )
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: me
    !> Number of particles to create at each time step
    integer, optional, intent(in) :: Nparticles
    ! ---------------------------------------------------------!
    ! Initialize the 2D growing arrays, giving them a first dimension of 6
    ! e.g. me%position(:,iParticle) = (x, y, z, rx, ry, rz)
    if (present(Nparticles)) then
      call init( me = me%position, width = 6, length = Nparticles )
      call init( me = me%velocity, width = 6, length = Nparticles )
      call init( me = me%force, width = 6, length = Nparticles )
      call init( me = me%radius, length = Nparticles )
      call init( me = me%mass, length = Nparticles )
      call init( me = me%IDoffset, length = Nparticles )
    else
      call init( me = me%position, width = 6 )
      call init( me = me%velocity, width = 6 )
      call init( me = me%force, width = 6 )
    end if
  
  end subroutine init_particle_creator
  
  !> Routine to create a new particle using the information in the particleCreator object
  subroutine createNewParticle_MEM( pos, vel, F, radius, mass, particleID, &
                                  & particleGroup, geometry, scheme, myRank )
    !> Create particle at this position
    real(kind=rk), intent(in) :: pos(6)
    !> Initial velocity of particle
    real(kind=rk), intent(in) :: vel(6)
    !> External forces on particle
    real(kind=rk), intent(in) :: F(6)
    !> Particle radius
    real(kind=rk), intent(in) :: radius 
    !> Particle mass
    real(kind=rk), intent(in) :: mass 
    !> ID that should be assigned to this particle
    integer, intent(in) :: particleID
    !> ParticleGroup to add this particle to
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme to initialize particles on the lattice
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry to initialize particles on the lattice 
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------ !
    type(mus_particle_MEM_type) :: dummyParticle
    real(kind=rk) :: dx
    logical :: wasAdded, rmflag
    integer :: iParticle, lev
    ! ------------------------------------------ !
    lev = geometry%tree%global%maxLevel
    dx = geometry%tree%global%BoundingCubeLength / 2**lev
    ! Set the dummy particle's properties
    dummyParticle%pos(1:6) = pos(1:6)
    dummyParticle%vel(1:6) = vel(1:6)
    dummyParticle%Fext(1:6) = F(1:6)
    dummyParticle%radius = radius
    dummyParticle%rotInertia = 0.4*mass*radius**2
    dummyParticle%Rn = ceiling( radius / dx )
    dummyParticle%mass = mass
    dummyParticle%particleID = particleID
  
    ! Append it to the particleGroup
    call append_DA_particle_MEM( me        = particleGroup%particles_MEM, &
                               & particle  = dummyParticle,               &
                               & length    = 1,                           &
                               & wasAdded  = wasAdded                     )
  
    iParticle = particleGroup%particles_MEM%nvals
    ! Initialize the communication masks
    call allocateProcessMasks(                                    &
      &    particle = particleGroup%particles_MEM%val(iParticle), &
      &    nProcs   = particleGroup%send%nProcs                   )
  
    ! Set initial coordOfOrigin
    particleGroup%particles_MEM%val(iParticle)%coordOfOrigin &
      &  = tem_CoordOfReal( &
      &      mesh  = geometry%tree,                                  &
      &      point = particleGroup%particles_MEM%val(iParticle)%pos, &
      &      level = lev                                             ) 
  
    ! Update particle coordOfOrigin: this ensures that both
    ! particle%coordOfOrigin and particle%oldCoordOfOrigin are set
    call updateCoordOfOrigin(                                  &
      &     this = particleGroup%particles_MEM%val(iParticle), &
      &     geometry = geometry                                )
  
    ! Set initial existsOnProc to FALSE. Will be set to true in initParticle_MEM
    ! If particle actually exists on proc
    particleGroup%particles_MEM%val(iParticle)%existsOnProc = .FALSE.
    
    ! Initialize particle representation on the grid. Also sets owner.
    call initParticle_MEM(                                            &
      &     particle    = particleGroup%particles_MEM%val(iParticle), &
      &     particleID  = iParticle,                                  &
      &     geometry    = geometry,                                   & 
      &     scheme      = scheme,                                     &
      &     myRank      = myRank,                                     &
      &     comm        = particleGroup%send,                         &
      &     rmflag      = rmflag                                      )
      
  
    if(rmflag) then
      call remove_particle_from_da_particle_MEM( particleGroup%particles_MEM, &
        &                                        iParticle                    )
    end if
    
    ! Set particle force to 0 initially
    particleGroup%particles_MEM%val(iParticle)%F(1:6) = 0.0_rk
  
  end subroutine createNewParticle_MEM
  
  !> Routine to create a new particle using the information in the
  !! particleCreator object
  subroutine createNewParticle_DPS( pos, vel, F, radius, mass, particleID,  &
    &                               particleGroup, geometry, scheme, myRank )
    !> Create particle at this position
    real(kind=rk), intent(in) :: pos(6)
    !> Initial velocity of particle
    real(kind=rk), intent(in) :: vel(6)
    !> External forces on particle
    real(kind=rk), intent(in) :: F(6)
    !> Particle radius
    real(kind=rk), intent(in) :: radius 
    !> Particle mass
    real(kind=rk), intent(in) :: mass 
    !> ID that should be assigned to this particle
    integer, intent(in) :: particleID
    !> ParticleGroup to add this particle to
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme to initialize particles on the lattice
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry to initialize particles on the lattice 
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------ !
    type(mus_particle_DPS_type) :: dummyParticle
    logical :: wasAdded
    integer :: iParticle
    ! ------------------------------------------ !
    ! Set the dummy particle's properties
    dummyParticle%pos(1:6) = pos(1:6)
    dummyParticle%vel(1:6) = vel(1:6)
    dummyParticle%Fext(1:6) = F(1:6)
    dummyParticle%radius = radius
    dummyParticle%mass = mass
    dummyParticle%rotInertia = 0.4*mass*radius**2
    dummyParticle%particleID = particleID
  
    ! Append it to the particleGroup
    call append_DA_particle_DPS( me        = particleGroup%particles_DPS, &
                               & particle  = dummyParticle,               &
                               & length    = 1,                           &
                               & wasAdded  = wasAdded                     )
  
    iParticle = particleGroup%particles_DPS%nvals
  
    ! Initialize the communication masks for this particle
    call allocateProcessMasks(                                    &
      &    particle = particleGroup%particles_DPS%val(iParticle), &
      &    nProcs   = particleGroup%send%nProcs                   )
  
    ! Initialize the new particle on the lattice
    call initParticle_DPS(                                             &
      &    particle     = particleGroup%particles_DPS%val(iParticle),  &
      &    interpolator = particleGroup%interpolator,                  &
      &    particleID   = particleID,                                  &
      &    geometry     = geometry,                                    & 
      &    scheme       = scheme,                                      &
      &    myRank       = myRank,                                      &
      &    comm         = particleGroup%send                           )
  
    if( particleGroup%particles_DPS%val(iParticle)%removeParticle_local ) then
      open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
        write(pgDebugLog%lu,*) "CREATE_NEW_PARTICLE DPS: Removing particle with ID ", &
        & particleGroup%particles_DPS%val(iParticle)%particleID
      call remove_particle_from_da_particle_DPS( particleGroup%particles_DPS, iParticle)
      close(pgDebugLog%lu)
    end if
  
  end subroutine createNewParticle_DPS
  

  !> Routine to generate a unique particleID for a particle
  function getNewParticleID(particle_creator, iParticle)
    !> Particle creator object
    type(mus_particle_creator_type), intent(in) :: particle_creator
    !> Index of the particle data we want to use to create this particle 
    !! in the particle creator object 
    integer :: iParticle
    ! ------------------------------ !
    ! Output particleID
    integer :: getNewParticleID
    ! ------------------------------ !
    getNewParticleID = particle_creator%IDoffset%val(iParticle) &
    & + particle_creator%N_times_called*particle_creator%global_Nparticles
  
  end function getNewParticleID
  

  !> Routine to create new fully resolved MEM particles
  subroutine create_particles_MEM( particle_creator, particleGroup, scheme, &
                                 & geometry, params, myRank                 )
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: particle_creator
    !> ParticleGroup to add this particle to
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme to initialize particles on the lattice
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry to initialize particles on the lattice 
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------- !
    integer :: k, particleID
    real(kind=rk) :: particle_vel(6)
    real(kind=rk) :: dx
    integer :: lev
    ! ------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    dx = params%physics%dxLvl(lev) 
  
    ! Re-initialize the particle creator positions with new normally distributed random values.
  
    open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
    do k = 1, particle_creator%Nparticles
      ! Generate a unique ID for this particle
      particleID = getNewParticleID( particle_creator = particle_creator, &
                                   & iParticle        = k                 )
      write(pgDebugLog%lu,*) "Creating new particle with ID = ", particleID
  
      if(particle_creator%init_particle_to_fluid_vel) then
        write(logUnit(1),* ) "ERROR: initializing particles to fluid velocity ", &
                              & "is not supported for MEM particles."
        call tem_abort()
      else
        ! Use the velocity stored in the particle_creator%velocity array
        particle_vel = particle_creator%velocity%val(1:6,k)
      end if
  
      ! Create the particle
      call createNewParticle_MEM( pos = particle_creator%position%val(1:6,k), &
                                & vel = particle_vel, & 
                                & F = particle_creator%force%val(1:6,k),      &
                                & radius = particle_creator%radius%val(k),    &
                                & mass = particle_creator%mass%val(k),        &
                                & particleID = particleID,                    &
                                & particleGroup = particleGroup,              &
                                & geometry = geometry,                        &
                                & scheme = scheme,                            &
                                & myRank = myRank                             )
  
    end do
    close(pgDebugLog%lu)
  
    particle_creator%N_times_called = particle_creator%N_times_called + 1
  
  end subroutine create_particles_MEM
  
  !> Routine to create new unresolved DPS particles
  subroutine create_particles_DPS( particle_creator, particleGroup, scheme, &
                                 & geometry, params, myRank                 )
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: particle_creator
    !> ParticleGroup to add this particle to
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme to initialize particles on the lattice
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry to initialize particles on the lattice 
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------- !
    integer :: k, particleID
    real(kind=rk) :: particle_vel(6)
    integer :: coord(4)
    real(kind=rk) :: vel_tmp(3), rho_tmp, eps_f_tmp
    real(kind=rk) :: dx
    integer :: lev
    ! ------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    dx = params%physics%dxLvl(lev) 
  
    ! Re-initialize the particle creator positions with new normally distributed random values.
  
    open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
    do k = 1, particle_creator%Nparticles
      ! Generate a unique ID for this particle
      particleID = getNewParticleID( particle_creator = particle_creator, &
                                   & iParticle        = k                 )
      write(pgDebugLog%lu,*) "Creating new particle with ID = ", particleID
  
      if(particle_creator%init_particle_to_fluid_vel) then
        ! Interpolate the local fluid velocity to the particle position
        coord = tem_CoordOfReal( mesh  = geometry%tree,                        &
                               & point = particle_creator%position%val(1:3,k), &
                               & level = geometry%tree%global%maxLevel         ) 
  
        call particleGroup%intp(                                   &
            & xp           = particle_creator%position%val(1:6,k), &
            & coord_xp     = coord,                                &
            & scheme       = scheme,                               &
            & geom_origin  = geometry%tree%global%origin,          &
            & dx           = dx,                                   &
            & interpolator = particleGroup%interpolator,           &
            & vel_xp       = vel_tmp,                              &
            & rho_xp       = rho_tmp,                              &
            & eps_f_xp     = eps_f_tmp                             )
  
        particle_vel(1:3) = vel_tmp
        particle_vel(4:6) = 0.0_rk
  
      else
        ! Use the velocity stored in the particle_creator%velocity array
        particle_vel = particle_creator%velocity%val(1:6,k)
      end if
  
      ! Create the particle
      call createNewParticle_DPS( pos = particle_creator%position%val(1:6,k), &
                                & vel = particle_vel,                         & 
                                & F = particle_creator%force%val(1:6,k),      &
                                & radius = particle_creator%radius%val(k),    &
                                & mass = particle_creator%mass%val(k),        &
                                & particleID = particleID,                    &
                                & particleGroup = particleGroup,              &
                                & geometry = geometry,                        &
                                & scheme = scheme,                            &
                                & myRank = myRank                             )
  
    end do
    close(pgDebugLog%lu)
  
    particle_creator%N_times_called = particle_creator%N_times_called + 1
  
  end subroutine create_particles_DPS
  

  !> Routine that checks if new particles should be created and creates them if so.
  subroutine check_and_create_new_particles_MEM( particle_creator, iter, particleGroup, &
    &                                        scheme, geometry, params, myRank           )
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: particle_creator
    !> Current LBM iteration
    integer, intent(in) :: iter
    !> particleGroup to add particles to
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme to initialize particles on the lattice
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry to initialize particles on the lattice 
    type(mus_geom_type), intent(in) :: geometry
    !> Params
    type(mus_param_type), intent(in) ::params
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------ !
    integer :: iParticle, lev
    logical :: createParticles, rmflag
    ! ------------------------------------------ !
    lev = geometry%tree%global%maxLevel
    rmflag = .FALSE.
  
    createParticles = must_create_new_particles(    &
        & iter     = iter,                          &
        & i_start  = particle_creator%iter_start,   &
        & i_end    = particle_creator%iter_end,     &
        & interval = particle_creator%iter_interval ) 
  
    if( createParticles ) then                    
      call create_particles_MEM( particle_creator = particle_creator, &
                              & particleGroup    = particleGroup, &
                              & scheme           = scheme, &
                              & geometry         = geometry, &
                              & params           = params, &
                              & myRank           = myRank )
  
      ! At this point particles have only been created on the process to which 
      ! they are local. So we need to exchange them to all other processes which are 
      ! within one lattice site of the particle position.
      call exchangeNewParticles_MEM( this         = particleGroup,       &
                                   & send         = particleGroup%send,  &
                                   & recv         = particleGroup%recv,  &
                                   & comm         = MPI_COMM_WORLD,      &
                                   & myRank       = myRank,              &
                                   & message_flag = 1                    )
  
      ! Initialize the new particles I have received
      do iParticle = 1, particleGroup%particles_MEM%nvals
        if( particleGroup%particles_MEM%val(iParticle)%newForMe ) then
  
          call allocateProcessMasks(                                    &
            &    particle = particleGroup%particles_MEM%val(iParticle), &
            &    nProcs   = particleGroup%send%nProcs                   )
  
          ! Set initial coordOfOrigin
          particleGroup%particles_MEM%val(iParticle)%coordOfOrigin &
            &  = tem_CoordOfReal( mesh  = geometry%tree,                               &
                                & point = particleGroup%particles_MEM%val(iParticle)%pos,  &
                                & level = lev                                          ) 
  
          ! Update particle coordOfOrigin: this ensures that both particle%coordOfOrigin and 
          ! particle%oldCoordOfOrigin are set
          call updateCoordOfOrigin( this = particleGroup%particles_MEM%val(iParticle), &
                                  & geometry = geometry                            )
  
          ! Set initial existsOnProc to FALSE. Will be set to true in initParticle_MEM
          ! If particle actually exists on proc
          particleGroup%particles_MEM%val(iParticle)%existsOnProc = .FALSE.
          
          ! Initialize particle representation on the grid. Also sets owner.
          call initParticle_MEM( particle    = particleGroup%particles_MEM%val(iParticle),  &
                                & particleID  = iParticle,                               &
                                & geometry    = geometry,                                & 
                                & scheme      = scheme,                                  &
                                & myRank      = myRank,                                  &
                                & comm        = particleGroup%send,                      &
                                & rmflag      = rmflag                                   )
            
  
          if(rmflag) then
            call remove_particle_from_da_particle_MEM( particleGroup%particles_MEM, iParticle )
          end if
          
          ! Set particle force to 0 initially
          particleGroup%particles_MEM%val(iParticle)%F(1:6) = 0.0_rk
  
          particleGroup%particles_MEM%val(iParticle)%newForMe = .FALSE. 
  
          open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
            write(pgDebugLog%lu,*) 'Initializing new particle with ID ', &
              & particleGroup%particles_MEM%val(iParticle)%particleID
            write(pgDebugLog%lu,*) 'pos = ', particleGroup%particles_MEM%val(iParticle)%pos(1:6)
            write(pgDebugLog%lu,*) 'vel = ', particleGroup%particles_MEM%val(iParticle)%vel(1:6)
            write(pgDebugLog%lu,*) 'F = ', particleGroup%particles_MEM%val(iParticle)%F(1:6)
            write(pgDebugLog%lu,*) 'F_DEM(1,:) = ', particleGroup%particles_MEM%val(iParticle)%Fbuff(1,1:6)
            write(pgDebugLog%lu,*) 'F_DEM(2,:) = ', particleGroup%particles_MEM%val(iParticle)%Fbuff(2,1:6)
            write(pgDebugLog%lu,*) 'coordOfOrigin = ', particleGroup%particles_MEM%val(iParticle)%coordOfOrigin(1:4)
          close( pgDebugLog%lu )
  
        end if
      end do ! initialization of new particles received from other processes.
  
  
      open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
      write(pgDebugLog%lu,*) "CREATED NEW PARTICLES iter = ", iter
      call printParticleGroup2_MEM(particleGroup, pgDebugLog%lu, myRank, iter) 
      close(pgDebugLog%lu)
    end if
  end subroutine check_and_create_new_particles_MEM
  
  !> Routine that checks if new particles should be created and creates them if so.
  subroutine check_and_create_new_particles_DPS( particle_creator, iter, particleGroup, &
    &                                        scheme, geometry, params, myRank           )
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: particle_creator
    !> Current LBM iteration
    integer, intent(in) :: iter
    !> particleGroup to add particles to
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme to initialize particles on the lattice
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry to initialize particles on the lattice 
    type(mus_geom_type), intent(in) :: geometry
    !> Params
    type(mus_param_type), intent(in) ::params
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------ !
    integer :: iParticle
    logical :: createParticles
    ! ------------------------------------------ !
    createParticles = must_create_new_particles(    &
        & iter     = iter,                          &
        & i_start  = particle_creator%iter_start,   &
        & i_end    = particle_creator%iter_end,     &
        & interval = particle_creator%iter_interval ) 
  
    if( createParticles ) then                    
      call create_particles_DPS( particle_creator = particle_creator, &
                              & particleGroup    = particleGroup, &
                              & scheme           = scheme, &
                              & geometry         = geometry, &
                              & params           = params, &
                              & myRank           = myRank )
  
      ! At this point particles have only been created on the process to which 
      ! they are local. So we need to exchange them to all other processes which are 
      ! within one lattice site of the particle position.
      call exchangeNewParticles_DPS( this         = particleGroup,       &
                                   & send         = particleGroup%send,  &
                                   & recv         = particleGroup%recv,  &
                                   & comm         = MPI_COMM_WORLD,      &
                                   & myRank       = myRank,              &
                                   & message_flag = 1                    )
  
      ! Initialize the new particles I have received
      do iParticle = 1, particleGroup%particles_DPS%nvals
        if( particleGroup%particles_DPS%val(iParticle)%newForMe ) then
  
  
          call allocateProcessMasks(                                    &
            &    particle = particleGroup%particles_DPS%val(iParticle), &
            &    nProcs   = particleGroup%send%nProcs                   )
  
          call initParticle_DPS( &
            & particle    = particleGroup%particles_DPS%val(iParticle),            &
            & interpolator = particleGroup%interpolator,                           &
            & particleID  = particleGroup%particles_DPS%val(iParticle)%particleID, &
            & geometry    = geometry,                                              &
            & scheme      = scheme,                                                &
            & myRank      = myRank,                                                &
            & comm        = particleGroup%send                                     )
  
          particleGroup%particles_DPS%val(iParticle)%newForMe = .FALSE. 
  
          open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
            write(pgDebugLog%lu,*) 'Initializing new particle with ID ', &
              & particleGroup%particles_DPS%val(iParticle)%particleID
            write(pgDebugLog%lu,*) 'pos = ', particleGroup%particles_DPS%val(iParticle)%pos(1:6)
            write(pgDebugLog%lu,*) 'vel = ', particleGroup%particles_DPS%val(iParticle)%vel(1:6)
            write(pgDebugLog%lu,*) 'F = ', particleGroup%particles_DPS%val(iParticle)%F(1:6)
            write(pgDebugLog%lu,*) 'F_DEM(1,:) = ', particleGroup%particles_DPS%val(iParticle)%F_DEM(1,1:6)
            write(pgDebugLog%lu,*) 'F_DEM(2,:) = ', particleGroup%particles_DPS%val(iParticle)%F_DEM(2,1:6)
            write(pgDebugLog%lu,*) 'coordOfOrigin = ', particleGroup%particles_DPS%val(iParticle)%coordOfOrigin(1:4)
          close( pgDebugLog%lu )
  
        end if
      end do ! initialization of new particles received from other processes.
  
  
      open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
      write(pgDebugLog%lu,*) "CREATED NEW PARTICLES iter = ", iter
      call printParticleGroup2_DPS(particleGroup, pgDebugLog%lu, myRank, iter) 
      close(pgDebugLog%lu)
    end if
  end subroutine check_and_create_new_particles_DPS
  

  !> Function to check whether particles should be created using the particleCreator 
  !! object at the current iteration iter.
  function must_create_new_particles(iter, i_start, i_end, interval)
    !> Current simulation time (in iterations)
    integer, intent(in) :: iter
    !> Start creating new particles at this time (in iterations)
    integer, intent(in) :: i_start
    !> Max time (in iterations) to create new particles
    integer, intent(in) :: i_end 
    !> Interval (in iterations) at which to create new particles
    integer, intent(in) :: interval
    !> Output logical: TRUE if particles should be created at this instant
    logical :: must_create_new_particles
    ! ----------------------------------------- !
    integer :: k
    ! ----------------------------------------- !
    if( iter < i_start .OR. iter > i_end ) then
      must_create_new_particles = .FALSE.
    else
      k = iter - i_start
      if( mod(k, interval) == 0 ) then
        must_create_new_particles = .TRUE. 
      else
        must_create_new_particles = .FALSE.
      end if
    end if
  end function
  
  !> Routine to initialize the particle creator object from a particle 
  !! blob cylinder object
  subroutine init_particle_creator_from_blob(           &
          & particle_creator, particleblob, Nparticles, &
          & scheme, geometry, myRank                    )
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: particle_creator
    !> Particle blob object describing the shape of the initial "blob" of particles
    type(mus_particle_blob_cylinder_type), intent(inout) :: particleblob
    !> Desired number of particles to fill the ENTIRE blob with (so not just the 
    !! part of it on this rank, but on all ranks)
    integer, intent(in) :: Nparticles
    !> Scheme to initialize particles on the lattice
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry to initialize particles on the lattice 
    type(mus_geom_type), intent(in) :: geometry
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------- !
    real(kind=rk) :: blob_length  ! length of particle blob cylinder
    real(kind=rk) :: d ! distance between particles
    real(kind=rk) :: R ! particle blob radius
    real(kind=rk) :: n_cylinder(3) ! unit vector in particleblob cylinder axis direction
    integer :: iPos, kParticle
    logical :: isLocal
  
    ! Array holding coordinates of all positions generated by the particleblob
    type(grw_real2darray_type) :: positions
    ! ------------------------------------------- !
    ! Initialize counter for global number of particles in domain (across all procs)
    kParticle = 1
    ! Get cylinder length L
    blob_length = dot_product( particleblob%vec, particleblob%vec ) 
    blob_length = sqrt(blob_length)
    n_cylinder = particleblob%vec/blob_length
    R = particleblob%radius
  
    ! Estimate the distance between particles to get Nparticles in total
    d = ( (blob_length*PI*R**2)/(Nparticles) )**(1.0_rk/3.0_rk)
  
    ! Check to make sure we are not placing particles at overlapping positions
    if( d < 2*particleblob%particle_radius ) then
      call tem_abort("ERROR init_particle_creator_from_blob: particles placed at overlapping positions!")
    end if
  
    call init( me = positions, width = 6 )
    ! Fill the cylinder with particles according to the chosen distribution kind
    select case( trim(particleblob%distribution%kind) )
    case ('uniform')
      call fill_cylinder( R         = particleblob%radius, &
                        & L         = blob_length,         &
                        & d         = d,                   &
                        & positions = positions            ) 
  
      ! Rotate the cylinder to the desired orientation 
      call rotate_positions( positions = positions,                 &
                          & n1        = [ 0.0_rk, 0.0_rk, 1.0_rk ], &
                          & n2        = n_cylinder                  )
  
      ! Translate the cylinder to the desired position
      call translate_positions( positions = positions, &
                            & translation_vec = particleblob%origin )
    case ('gaussian')
      call fill_blob_positions_gauss( particleblob = particleblob, &
                                    & positions    = positions,    &
                                    & Nparticles   = Nparticles    )
  
    end select
  
    !> Check if particles should be initialized to the local fluid velocity.
    if(particleblob%init_particles_to_fluid_vel) then
      particle_creator%init_particle_to_fluid_vel = .TRUE.
    else
      particle_creator%init_particle_to_fluid_vel = .FALSE.
    end if
  
  
    ! By now the positions of all particles are determined. These positions should 
    ! only be added to the particle_creator if they are on a local fluid cell of 
    ! this MPI rank. 
    
    ! Initialize the particle creators position array
    call init( me=particle_creator%position, width=6 )
  
    ! Set velocity, force, radius and mass of particles
    do iPos = 1, positions%nvals
      islocal = positionLocalOnMyRank(                       &
                      & pos      = positions%val(1:3, iPos), &
                      & geometry = geometry,                 &
                      & scheme   = scheme,                   &
                      & myRank   = myRank                    )
      if(isLocal) then
        call append( me  = particle_creator%position, &
                   & val = positions%val(1:6,iPos)    )
  
        call append( me  = particle_creator%velocity, &
                   & val = particleblob%particle_vel  )
  
        call append( me  = particle_creator%force,     &
                   & val = particleblob%particle_force )
  
        call append( me = particle_creator%radius,      &
                   & val = particleblob%particle_radius )
  
        call append( me  = particle_creator%mass,     &
                   & val = particleblob%particle_mass )
  
        call append( me  = particle_creator%IDoffset, &
                   & val = kparticle                  )
      end if
      kParticle = kParticle + 1
    end do
    
    ! Set number of particles in the creator = equal to the number of entries 
    ! in the growing arrays.
    particle_creator%Nparticles = particle_creator%radius%nvals
  
    ! Destroy the positions array to free memory
    call destroy(me = positions )
  
  end subroutine init_particle_creator_from_blob
  
  !> Routine to initialize the particle creator object from a particle 
  !! blob prism object
  subroutine init_particle_creator_from_blob_prism(     &
          & particle_creator, particleblob, Nparticles, &
          & scheme, geometry, myRank                    )
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: particle_creator
    !> Particle blob object describing the shape of the initial "blob" of particles
    type(mus_particle_blob_prism_type), intent(inout) :: particleblob
    !> Desired number of particles to fill the ENTIRE blob with (so not just the 
    !! part of it on this rank, but on all ranks)
    integer, intent(in) :: Nparticles
    !> Scheme to initialize particles on the lattice
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry to initialize particles on the lattice 
    type(mus_geom_type), intent(in) :: geometry
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------- !
    real(kind=rk) :: d ! distance between particles
    real(kind=rk) :: lx, ly, lz  ! lengths of the blob prism
    real(kind=rk) :: n_prism(3) ! unit vector in particleblob cylinder axis direction
    integer :: iPos, kParticle
    logical :: isLocal
  
    ! Array holding coordinates of all positions generated by the particleblob
    type(grw_real2darray_type) :: positions
    ! ------------------------------------------- !
    ! Initialize counter for global number of particles in domain (across all procs)
    kParticle = 1
  
    call init( me = positions, width = 6 )
  
    ! Do particleblob prism stuff here
    lx = sqrt( dot_product(particleblob%vec_x, particleblob%vec_x) )
    ly = sqrt( dot_product(particleblob%vec_y, particleblob%vec_y) )
    lz = sqrt( dot_product(particleblob%vec_z, particleblob%vec_z) )
    n_prism = particleblob%vec_x / lx
  
    d = ( lx*ly*lz / Nparticles )**(1.0_rk/3.0_rk)
    ! Check to make sure we are not placing particles at overlapping positions
    if( d < 2*particleblob%particle_radius ) then
      call tem_abort("ERROR init_particle_creator_from_blob: particles placed at overlapping positions!")
    end if
  
    select case( trim(particleblob%distribution%kind) )
    case ('uniform')
      if(particleblob%nx < 0 .OR. particleblob%ny < 0 .OR. particleblob%nz < 0) then
        call fill_prism( lx       = lx,       &
                      & ly        = ly,       &
                      & lz        = lz,       &
                      & d         = d,        &
                      & positions = positions )
      else 
        call fill_prism( lx       = lx,              &
                      & ly        = ly,              &
                      & lz        = lz,              &
                      & nx        = particleblob%nx, &
                      & ny        = particleblob%ny, &
                      & nz        = particleblob%nz, &
                      & positions = positions        )
      end if
    case('gaussian')
      call fill_blob_positions_gauss_prism( particleblob = particleblob, &
                                          & positions    = positions,    &
                                          & Nparticles   = Nparticles    )
    end select
  
    !> Check if particles should be initialized to the local fluid velocity.
    if(particleblob%init_particles_to_fluid_vel) then
      particle_creator%init_particle_to_fluid_vel = .TRUE.
    else
      particle_creator%init_particle_to_fluid_vel = .FALSE.
    end if
  
    ! By now the positions of all particles are determined. These positions should 
    ! only be added to the particle_creator if they are on a local fluid cell of 
    ! this MPI rank. 
  
    ! Initialize the particle creators position array
    call init( me=particle_creator%position, width=6 )
    
    ! Set velocity, force, radius and mass of particles
    do iPos = 1, positions%nvals
      islocal = positionLocalOnMyRank(                       &
                      & pos      = positions%val(1:3, iPos), &
                      & geometry = geometry,                 &
                      & scheme   = scheme,                   &
                      & myRank   = myRank                    )
      if(isLocal) then
        call append( me  = particle_creator%position, &
                   & val = positions%val(1:6,iPos)    )
  
        call append( me  = particle_creator%velocity, &
                   & val = particleblob%particle_vel  )
  
        call append( me  = particle_creator%force,     &
                   & val = particleblob%particle_force )
  
        call append( me = particle_creator%radius,      &
                   & val = particleblob%particle_radius )
  
        call append( me  = particle_creator%mass,     &
                   & val = particleblob%particle_mass )
  
        call append( me  = particle_creator%IDoffset, &
                   & val = kparticle                  )
      end if
      kParticle = kParticle + 1
    end do
  
    ! Set number of particles in the creator = equal to the number of entries 
    ! in the growing arrays.
    particle_creator%Nparticles = particle_creator%radius%nvals
    write(logUnit(1),*) "Particle creator Nparticles = ", particle_creator%Nparticles
  
    ! Destroy the positions array to free memory
    call destroy(me = positions )
  
  end subroutine init_particle_creator_from_blob_prism
  

  !> Routine to initialize particle creator object with random positions 
  !! inside a cylinder described by blob_cylinder type
  subroutine fill_blob_positions_gauss( particleblob, positions, Nparticles )
    !> Particleblob type
    type(mus_particle_blob_cylinder_type), intent(inout) :: particleblob
    !> Positions array to append randomly created values to
    type(grw_real2darray_type), intent(inout) :: positions
    !> Desired number of particles
    integer, intent(in) :: Nparticles
    ! --------------------------------------------------- !
    integer :: iSlice, iParticle, Nslices, Nparticles_per_slice
    real(kind=rk) :: blob_length, pos(6), n_cylinder(3)
    real(kind=rk) :: d, R
    ! --------------------------------------------------- !
    R = particleblob%radius
    ! Initialize the random positiong generator
    blob_length = dot_product( particleblob%vec, particleblob%vec ) 
    blob_length = sqrt(blob_length)
    n_cylinder = particleblob%vec/blob_length
  
    ! Estimate the distance between particles to get Nparticles in total
    d = ( (blob_length*PI*R**2)/(Nparticles) )**(1.0_rk/3.0_rk)
  
    Nslices = floor(blob_length/d) 
    Nparticles_per_slice = Nparticles/Nslices
  
    ! Set the random number generator
    call set_random_seed( particleblob%distribution%seed )
  
    write( logUnit(1),* ) "Random seed set to ", particleblob%distribution%seed
    ! Construct the 3D positions by stacking several of these slices together 
    do iSlice = 1, Nslices
      call init_particle_blob_prob( blob    = particleblob%distribution,       &
                                  & d       = 2*particleblob%particle_radius,  &
                                  & R       = particleblob%radius,             &
                                  & mu      = particleblob%distribution%mu,    &
                                  & sigma   = particleblob%distribution%sigma, &
                                  & Nchosen = Nparticles_per_slice             )
                                  
      ! Pick (almost) random non-overlapping positions in the x, y plane
      do iParticle = 1, Nparticles_per_slice
        call pick_random_position( blob = particleblob%distribution )
      end do
  
      ! Now add these random positions to the big positions array
      do iParticle = 1, particleblob%distribution%chosen_positions%nvals
        pos(1) = particleblob%distribution%chosen_positions%val(1,iParticle)
        pos(2) = particleblob%distribution%chosen_positions%val(2,iParticle)
        pos(3) = d*(iSlice - 1)
        pos(4:6) = 0.0_rk
        call append( me = positions, val = pos)
      end do ! iParticle
  
    end do ! iSlice
  
    ! Rotate the positions inside the cylinder to the desired orientation 
    call rotate_positions( positions = positions,                  &
                          & n1        = [ 0.0_rk, 0.0_rk, 1.0_rk ], &
                          & n2        = n_cylinder                  )
  
    ! Translate the cylinder to the desired position
    call translate_positions( positions       = positions,          &
                            & translation_vec = particleblob%origin )
  
  end subroutine fill_blob_positions_gauss
  

  !> Routine to initialize particle creator object with random positions 
  !! inside a prism described by blob_prism type
  subroutine fill_blob_positions_gauss_prism( particleblob, positions, Nparticles )
    !> Particleblob type
    type(mus_particle_blob_prism_type), intent(inout) :: particleblob
    !> Positions array to append randomly created values to
    type(grw_real2darray_type), intent(inout) :: positions
    !> Desired number of particles
    integer, intent(in) :: Nparticles
    ! --------------------------------------------------- !
    integer :: iSlice, iParticle, Nslices, Nparticles_per_slice
    real(kind=rk) :: Lx, Ly, Lz, pos(6)
    real(kind=rk) :: n_z(3)
    real(kind=rk) :: d, dslice
    ! --------------------------------------------------- !
    ! Initialize the random positiong generator
    Lx = dot_product( particleblob%vec_x, particleblob%vec_x ) 
    Lx = sqrt(Lx)
    Ly = dot_product( particleblob%vec_y, particleblob%vec_y ) 
    Ly = sqrt(Ly)
    Lz = dot_product( particleblob%vec_z, particleblob%vec_z ) 
    Lz = sqrt(Lz)
  
    write(logUnit(1),*) "Lx = ", Lx
    write(logUnit(1),*) "Ly = ", Ly
    write(logUnit(1),*) "Lz = ", Lz
  
    ! Calculate the normal vector pointing in the z-direction of 
    ! the prism
    n_z = particleblob%vec_z/Lz
  
    ! Estimate the distance between particles to get Nparticles in total
    d = ( Lx*Ly*Lz/Nparticles )**(1.0_rk/3.0_rk)
  
    Nslices = floor(Lz/d) 
    Nparticles_per_slice = Nparticles/Nslices
    dslice = Lz/Nslices
  
    ! Set the random number generator
    call set_random_seed( particleblob%distribution%seed )
  
    write( logUnit(1),* ) "Random seed set to ", particleblob%distribution%seed
    ! Construct the 3D positions by stacking several of these slices together 
    do iSlice = 1, Nslices
      call init_particle_blob_prob( blob    = particleblob%distribution,       &
                                  & d       = 2*particleblob%particle_radius,  &
                                  & xlim    = [-0.5*Lx, 0.5*Lx],               &
                                  & ylim    = [-0.5*Ly, 0.5*Ly],               &
                                  & mu      = particleblob%distribution%mu,    &
                                  & sigma   = particleblob%distribution%sigma, &
                                  & Nchosen = Nparticles_per_slice             )
                                  
      ! Pick (almost) random non-overlapping positions in the x, y plane
      do iParticle = 1, Nparticles_per_slice
        call pick_random_position( blob = particleblob%distribution )
      end do
  
      ! Now add these random positions to the big positions array
      do iParticle = 1, particleblob%distribution%chosen_positions%nvals
        pos(1) = particleblob%distribution%chosen_positions%val(1,iParticle)
        pos(2) = particleblob%distribution%chosen_positions%val(2,iParticle)
        pos(3) = dslice*(iSlice - 1)
        pos(4:6) = 0.0_rk
        call append( me = positions, val = pos)
      end do ! iParticle
  
    end do ! iSlice
  
    ! Rotate the positions inside the cylinder to the desired orientation 
    call rotate_positions( positions = positions,                  &
                          & n1       = [ 0.0_rk, 0.0_rk, 1.0_rk ], &
                          & n2       = n_z                         )
  
    ! Translate the cylinder to the desired position
    call translate_positions( positions       = positions,          &
                            & translation_vec = particleblob%origin )
  
  end subroutine fill_blob_positions_gauss_prism
  
  !> Print the data in particle creator object, used for debugging
  subroutine print_particle_creator(particle_creator, logUnit)
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: particle_creator
    !> Unit to write to
    integer, intent(in) :: logUnit
    ! ------------------------------------------- !
    integer :: iParticle, Nparticles, i
    ! ------------------------------------------- !
    Nparticles = particle_creator%radius%nVals
  
    write(logUnit,'(A)') "--------- Particle creator data ---------"
    write(logUnit,'(A)') "--- TIMING ---"
    write(logUnit,'(A,I8)') "iter_start = ", particle_creator%iter_start
    write(logUnit,'(A,I8)') "iter_end = ", particle_creator%iter_end
    write(logUnit,'(A,I8)') "iter_interval = ", particle_creator%iter_interval
    write(logUnit,'(A,I8)') "N_times_called = ", particle_creator%N_times_called
    write(logUnit,'(A,I8)') "global_Nparticles = ", particle_creator%global_Nparticles
  
    write(logUnit,'(A,L3)') "init_particles_to_fluid_vel = ", particle_creator%init_particle_to_fluid_vel
  
    write(logUnit,*) "--- PARTICLES ---"
    do iParticle = 1, Nparticles
      write(logUnit,'(A,I4,A)') "-----------", iParticle, "-----------"
      write(logUnit,'(A)') "Position"
      do i = 1, 6
        write(logUnit,'(E17.9)', advance='no') particle_creator%position%val(i,iParticle)
        write(logUnit,'(A)', advance='no') " "
      end do
      write(logUnit,'(A)') 
      write(logUnit,'(A)') "velocity"
      do i = 1, 6
        write(logUnit,'(E17.9)', advance='no') particle_creator%velocity%val(i,iParticle)
        write(logUnit,'(A)', advance='no') " "
      end do
      write(logUnit,'(A)') 
      write(logUnit,'(A)') "force"
      do i = 1, 6
        write(logUnit,'(E17.9)', advance='no') particle_creator%force%val(i,iParticle)
        write(logUnit,'(A)', advance='no') " "
      end do
      write(logUnit,'(A)') 
      write(logUnit,'(A, E17.9)') "radius ", particle_creator%radius%val(iParticle)
      write(logUnit,'(A, E17.9)') "mass ", particle_creator%mass%val(iParticle)
      write(logUnit,'(A, I0.9)') "IDoffset ", particle_creator%IDoffset%val(iParticle)
    end do
    write(logUnit,'(A)') "---------------------------------------"
  end subroutine print_particle_creator
  
  !> Print the positions in particle creator object, used for debugging
  subroutine print_particle_creator_positions(particle_creator, logUnit)
    !> Particle creator object
    type(mus_particle_creator_type), intent(inout) :: particle_creator
    !> Unit to write to
    integer, intent(in) :: logUnit
    ! ------------------------------------------- !
    integer :: iParticle, Nparticles, i
    ! ------------------------------------------- !
    Nparticles = particle_creator%radius%nVals
  
    write(logUnit,'(A)') "--------- Particle creator data ---------"
    write(logUnit,'(A)') "--- TIMING ---"
    write(logUnit,'(A,I8)') "iter_start = ", particle_creator%iter_start
    write(logUnit,'(A,I8)') "iter_end = ", particle_creator%iter_end
    write(logUnit,'(A,I8)') "iter_interval = ", particle_creator%iter_interval
    write(logUnit,'(A,I8)') "N_times_called = ", particle_creator%N_times_called
    write(logUnit,'(A,I8)') "global_Nparticles = ", particle_creator%global_Nparticles
  
    write(logUnit,'(A,L3)') "init_particles_to_fluid_vel = ", particle_creator%init_particle_to_fluid_vel
  
    write(logUnit,*) "--- PARTICLE POSITIONS ---"
    do iParticle = 1, Nparticles
      do i = 1, 6
        write(logUnit,'(E17.9)', advance='no') particle_creator%position%val(i,iParticle)
        write(logUnit,'(A)', advance='no') " "
      end do
      write(logUnit,'(A)') 
    end do
  end subroutine print_particle_creator_positions

end module mus_particle_creator_module
