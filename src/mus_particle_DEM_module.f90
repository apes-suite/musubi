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
!> In this module the core routines of the Discrete Element Method (DEM) for 
!! simulating particle-laden flows are implemented.
!!
!! Two types of particle simulations are supported: unresolved Discrete Particle
!! Simulations (DPS) and fully resolved Momentum-Exchange Method simulations.

module mus_particle_DEM_module

  use env_module,           only: rk, long_k, newUnit
  use tem_param_module,     only: PI
  use tem_logging_module,   only: logUnit
  use tem_aux_module,       only: tem_abort
  use tem_dyn_array_module, only: init, append, destroy, &
    &                             empty, dyn_intArray_type
  use tem_timer_module,     only: tem_startTimer, tem_stopTimer

  use mus_geom_module,          only: mus_geom_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_param_module,         only: mus_param_type

  use mus_particle_type_module, only: mus_particle_MEM_type, &
    &                                 mus_particle_DPS_type, &
    &                                 mus_particle_group_type
  use mus_particle_comm_type_module, only: mus_particles_communication_type 
  use mus_particle_comm_module, only: exchangeParticlesToRemove, &
    &                                 DEM_exchangeForces,        &
    &                                 DEM_exchangeForces_DPS,    &
    &                                 exchangePositions,         &
    &                                 exchangePositions_DPS,     &
    &                                 exchangeMomInc_DPS,        &
    &                                 exchangeVelocities,        &
    &                                 exchangeVelocities_DPS,    &
    &                                 exchangeHydroforces_DPS,   &
    &                                 addCollisionForceToBuffer, &
    &                                 resetForceBuffers,         &
    &                                 DEM_exchangeWallPositions, &
    &                                 DEM_exchangeWallPositions_DPS
  use mus_particle_logging_module, only: openLogFile, &
    &                                    dumpdata,    &
    &                                    mus_particles_log_total_momentum
  use mus_particle_logging_type_module, only: mus_particle_logging_type, &
    &                                         pgDebugLog
  use mus_particle_MEM_module, only: updateCoordOfOrigin       
  use mus_particle_DPS_module, only: updateCoordOfOrigin_DPS,        &
    &                                applyDragForce_DPS,             &
    &                                applyDragForce_DPS_noeps,       &
    &                                applyLiftForce_DPS,             &
    &                                mapParticlesToLattice_DPS,      &
    &                                interpolateFluidProperties_DPS, &
    &                                incrementAuxField_DPS,          &
    &                                recalculate_auxField_DPS,       &
    &                                applyPressureForce_DPS,         &
    &                                mus_particles_DPS_interpolateFluidProperties
  use mus_particle_interactions_module, only: isLocalCollision,      &
    &                                         computeWallPosSum,     &
    &                                         DEM_isRemoteCollision, &
    &                                         checkAndCollideDEM,    &
    &                                         DEM_collideWithWall,   &
    &                                         computeWallForce_1D,   &
    &                                         DEM_computeWallForces_MEM
  
  use mus_particle_boundary_module, only: pgBndData,                      & 
    &                                     mus_particle_boundaryData_type, &
    &                                     wrapPeriodicPos,                &
    &                                     computeDisplacement
  use mus_particle_timer_module, only: mus_particle_timerHandles

  implicit none
  
  interface updatePositionVerlet
    !> Interface for routines to update particle positions
    !! using Verlet integration
    module procedure updatePositionVerlet_MEM
    module procedure updatePositionVerlet_DPS
  end interface
  
  interface updatePositionEuler
    !> Interface for routines to update particle positions
    !! using Euler integration
    module procedure updatePositionEuler_MEM
    module procedure updatePositionEuler_DPS
  end interface
  
  interface updateVelocityEuler
    !> Interface for routines to update particle velocities
    !! using Euler integration
    module procedure updateVelocityEuler_MEM
    module procedure updateVelocityEuler_DPS
  end interface
  
  interface updateVelocityVerlet
    !> Interface for routines to update particle velocities
    !! using Verlet integration
    module procedure updateVelocityVerlet_MEM
    module procedure updateVelocityVerlet_DPS
  end interface
  
  interface DEM_swapFnowFnext
    !> Interface for routines to swap the indices pointing 
    !! to the force in the current (Fnow) and next time step 
    !! (Fnext). Used for averaging force over two time steps.
    module procedure DEM_swapFnowFnext_MEM
    module procedure DEM_swapFnowFnext_DPS
  end interface
  
  ! For debugging: we temporarily store the values of the computed lift and 
  ! drag force in these variables so we can log them.
  real(kind=rk) :: logger_Flift(3) = 0.0_rk
  real(kind=rk) :: logger_Fdrag(3) = 0.0_rk


contains


  !> DEMsubcycles runs Nsubcycles time integration steps per LBM time step 
  !! In this routine simple Forward Euler integration of particle position and
  !! velocity is used.
  subroutine DEMSubcycles_MEM( particleGroup, scheme, geometry, params, Nsubcycles )
    !> Array of particles
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Number of subcycles
    integer, intent(in) :: Nsubcycles
    ! ------------------------------------------!
    real(kind=rk) :: dt, dx       ! LBM time step
    real(kind=rk) :: dt_sub       ! subcycle time step
    real(kind=rk) :: eps, Tc      ! collision tolerance and collision time
    real(kind=rk) :: mu, nu_lat   ! fluid dynamic viscosity (physical units)
    integer :: ksub               ! iterator over subcycle loops
    integer :: iParticle
    integer :: myRank             ! This process rank
    integer :: lev
    ! ------------------------------------------!
    lev = geometry%tree%global%maxLevel
    dt = params%physics%dtLvl(lev)
    dx = params%physics%dxLvl(lev)
    myRank = params%general%proc%rank
    eps = particleGroup%collision_tol
    Tc = particleGroup%collision_time
    dt_sub = dt/Nsubcycles
  
    ! Get viscosity values
    nu_lat = scheme%field(1)%fieldProp%fluid%viscKine%dataOnLvl(lev)%val(1)
    ! Note dynamic viscosity in lattice units is same as kinematic visc in lat units
    mu = nu_lat * params%physics%fac(lev)%viscDyna
  
    ! ---- INITIALIZE F_DEM ---- !
    call DEM_storeWallPositions_MEM( particleGroup, scheme, geometry, &
                                     & params, Tc, 0.0_rk         ) 
  
    call DEM_resetFnext_MEM( particleGroup = particleGroup )
  
    call DEM_computeLocalCollisionForces( particleGroup = particleGroup, &
                                        & myRank        = myRank,        &
                                        & eps           = eps,        &
                                        & Tc            = Tc,            &
                                        & mu            = mu             )
  
  
    ! Compute collision forces where one of the particles is remote
    call DEM_computeRemoteCollisionForces( particleGroup = particleGroup, &
                                         & myRank        = myRank,        &
                                         & eps           = eps,           &
                                         & Tc            = Tc,            &
                                         & mu            = mu             )
  
    ! Compute collision forces between particles and walls
    call DEM_computeWallForces_MEM( particleGroup = particleGroup, &
                                  & scheme        = scheme,        &
                                  & geometry      = geometry,      &
                                  & params        = params,        &
                                  & Tc            = Tc,            &
                                  & eps           = eps            )
  
  
    ! Add external forces (hydrodynamic and body forces) to F_DEM(next,:)
    call DEM_computeExternalForces( particleGroup = particleGroup )
  
    ! Communicate forces: owner of the particle gets all DEM force contributions
    ! from other processes
    call DEM_exchangeForces( this         = particleGroup,            &
                           & send         = particleGroup%send,       &
                           & recv         = particleGroup%recv,       &
                           & comm         = params%general%proc%comm, &
                           & myRank       = myRank,                   &
                           & message_flag = 1                         )
  
    ! Swap force buffer for subcycle (collision) forces
    do iParticle = 1, particleGroup%particles_MEM%nvals
      call DEM_swapFnowFnext( this = particleGroup%particles_MEM%val(iParticle) )
    end do
  
    do ksub = 1, Nsubcycles
  
      ! --- POSITION UPDATE --- !
      do iParticle = 1, particleGroup%particles_MEM%nvals
        if(particleGroup%particles_MEM%val(iParticle)%owner &
          &                         == params%general%proc%rank) then
          ! Owner updates the position and integer coordOfOrigin of the particle
          call updatePositionVerlet( particleGroup%particles_MEM%val(iParticle), dt_sub)
  
          call wrapPeriodicPos( pos  =  particleGroup%particles_MEM%val(iParticle)%pos(1:3), &
                              & boundaryData = pgBndData                        )
  
          call updateCoordOfOrigin( this = particleGroup%particles_MEM%val(iParticle), &
                                      & geometry = geometry                            )
  
        end if
      end do
  
      call exchangePositions( this = particleGroup,              &
                            & send = particleGroup%send,         &
                            & recv = particleGroup%recv,         &
                            & comm = params%general%proc%comm,   &
                            & myRank = params%general%proc%rank, &
                            & message_flag = 1                   )
  
      ! --- DEM FORCE COMPUTATION --- !
      call DEM_resetFnext_MEM( particleGroup = particleGroup )
  
      ! Reset the force buffers which will be filled in 
      ! the DEM collision routines for particles and walls 
      call resetForceBuffers( particleGroup )
  
      call DEM_computeLocalCollisionForces( particleGroup = particleGroup, &
                                          & myRank        = myRank,        &
                                          & eps           = eps,           &
                                          & Tc            = Tc,            &
                                          & mu            = mu             )
  
  
      ! Compute collision forces where one of the particles is remote
      call DEM_computeRemoteCollisionForces( particleGroup = particleGroup, &
                                          & myRank        = myRank,        &
                                          & eps           = eps,           &
                                          & Tc            = Tc,            &
                                          & mu            = mu             )
  
      ! Compute collision forces between particles and walls
      call DEM_computeWallForces_MEM( particleGroup = particleGroup, &
                                    & scheme        = scheme,        &
                                    & geometry      = geometry,      &
                                    & params        = params,        &
                                    & Tc            = Tc,            &
                                    & eps           = eps            )
  
  
      ! Add external forces (hydrodynamic and body forces) to F_DEM(next,:)
      ! NOTE: These are not added to force buffer
      call DEM_computeExternalForces( particleGroup = particleGroup )
  
      ! Communicate forces: owner of the particle gets all DEM force contributions
      ! from other processes
      call DEM_exchangeForces( this         = particleGroup,            &
                             & send         = particleGroup%send,       &
                             & recv         = particleGroup%recv,       &
                             & comm         = params%general%proc%comm, &
                             & myRank       = myRank,                   &
                             & message_flag = 1                         )
  
  
      ! --- VELOCITY UPDATE --- !
      do iParticle = 1, particleGroup%particles_MEM%nvals
        if(particleGroup%particles_MEM%val(iParticle)%owner &
          &                         == params%general%proc%rank) then
  
            ! Only update velocities of particle for which I am owner
            call updateVelocityVerlet(particleGroup%particles_MEM%val(iParticle), dt_sub )
  
        end if
      end do
  
      ! --- SYNCHRONIZE VELOCITIES --- !
      ! * Owner sends updated velocities to all other processes on which
      !   the particle exists
      call  exchangeVelocities( this = particleGroup,              &
                              & send = particleGroup%send,         &
                              & recv = particleGroup%recv,         &
                              & comm = params%general%proc%comm,   &
                              & myRank = params%general%proc%rank, &
                              & message_flag = 1                   )
  
  
      ! --- SWAP SUBCYCLE FORCE BUFFER --- !
      do iParticle = 1, particleGroup%particles_MEM%nvals
        call DEM_swapFnowFnext( this = particleGroup%particles_MEM%val(iParticle) )
      end do
    end do ! ksub
  
  end subroutine DEMSubcycles_MEM
  
  !> DEMsubcycles runs Nsubcycles time integration steps per LBM time step 
  !! In this routine velocity verlet integration of particle position and
  !! velocity is used. 
  subroutine DEMSubcycles_DPS( particleGroup, scheme, geometry, &
                             & params, Nsubcycles               )
    !> Array of particles
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Number of subcycles
    integer, intent(in) :: Nsubcycles
    ! ------------------------------------------!
    real(kind=rk) :: dt, dx       ! LBM time step
    real(kind=rk) :: dt_sub       ! subcycle time step
    real(kind=rk) :: nu, nu_lat   ! dynamic viscosity for computing hydro forces
    integer :: ksub               ! iterator over subcycle loops
    integer :: iParticle
    integer :: myRank             ! This process rank
    integer :: lev
  
    real(kind=rk) :: t
    real(kind=rk) :: dt_sub_lat
    
    ! ------------------------------------------!
    lev = geometry%tree%global%maxLevel
    dt = params%physics%dtLvl(lev)
    dx = params%physics%dxLvl(lev)
    myRank = params%general%proc%rank
    dt_sub = dt/Nsubcycles
    dt_sub_lat = 1.0_rk/Nsubcycles
  
    ! Get viscosity values
    nu_lat = scheme%field(1)%fieldProp%fluid%viscKine%dataOnLvl(lev)%val(1)
    nu = nu_lat * params%physics%fac(lev)%visc
  
    ! ---- INITIALIZE SUBCYCLING LOOP ---- !
  
    ! Fill neighbor list for efficient collision detection during subcycles
    if(particleGroup%enableCollisions) then
      call DEM_fillNeighborList( particleGroup = particleGroup, &
                               & d0            = 3*dx           )
    end if
  
    ! Calculate initial forces on particles
    call updateParticleForces( particleGroup = particleGroup, &
                             & scheme        = scheme,        &
                             & geometry      = geometry,      &
                             & params        = params,        &
                             & boundaryData  = pgBndData      )
  
  
    ! Swap force buffer for subcycle (collision) forces
    do iParticle = 1, particleGroup%particles_DPS%nvals
      call DEM_swapFnowFnext( this = particleGroup%particles_DPS%val(iParticle) )
    end do
  
    ! --- START OF ACTUAL SUBCYCLING LOOP --- !
    ! Start subcycle timer
    call tem_startTimer(timerHandle = mus_particle_timerHandles%subcycleTimer)
    do ksub = 1, Nsubcycles
      t = params%general%simControl%now%sim + dt_sub*ksub
  
      
      call tem_startTimer(timerHandle = mus_particle_timerHandles%positionTimer)
      ! --- POSITION UPDATE --- !
      do iParticle = 1, particleGroup%particles_DPS%nvals
        if(particleGroup%particles_DPS%val(iParticle)%owner &
          &                         == params%general%proc%rank) then
          ! Owner updates the position of the particle
          call updatePositionVerlet( particleGroup%particles_DPS%val(iParticle), dt_sub)
  
          call wrapPeriodicPos( pos  =  particleGroup%particles_DPS%val(iParticle)%pos(1:3), &
                              & boundaryData = pgBndData                        )
  
  
        end if
      end do
  
      call tem_startTimer(timerHandle = mus_particle_timerHandles%exchangePositionsTimer)
      call exchangePositions_DPS( this = particleGroup,              &
                                & send = particleGroup%send,         &
                                & recv = particleGroup%recv,         &
                                & comm = params%general%proc%comm,   &
                                & myRank = params%general%proc%rank, &
                                & message_flag = 1                   )
      call tem_stopTimer(timerHandle = mus_particle_timerHandles%exchangePositionsTimer)
  
      ! --- MAP PARTICLES TO LATTICE --- !
      ! And add/remove particles to/from processes as needed
      call mapParticlesToLattice_DPS( particleGroup = particleGroup, &
                                    & scheme        = scheme,        &
                                    & geometry      = geometry,      &
                                    & params        = params         )
  
      call tem_stopTimer(timerHandle = mus_particle_timerHandles%positionTimer)
  
      call tem_startTimer(timerHandle = mus_particle_timerHandles%exchangeMomIncTimer)
      call exchangeMomInc_DPS( this = particleGroup,              &
                             & send = particleGroup%send,         &
                             & recv = particleGroup%recv,         &
                             & comm = params%general%proc%comm,   &
                             & myRank = params%general%proc%rank, &
                             & message_flag = 1                   )
      call tem_stopTimer(timerHandle = mus_particle_timerHandles%exchangeMomIncTimer)
  
  
      call tem_startTimer(timerHandle = mus_particle_timerHandles%forceTimer)
      ! Update DEM forces on all particles and 
      ! synchronize across processes
  
      call updateParticleForces( particleGroup = particleGroup, &
                               & scheme        = scheme,        &
                               & geometry      = geometry,      &
                               & params        = params,        &
                               & boundaryData  = pgBndData      )
  
  
      ! Increment the fluid velocity using the forces by particles on fluid
      call tem_startTimer(timerHandle = mus_particle_timerHandles%incrementAuxFieldTimer)
  
      !$OMP PARALLEL
      !$OMP DO
      do iParticle = 1, particleGroup%particles_DPS%nvals
        call incrementAuxField_DPS(                                            &
            & particle         = particleGroup%particles_DPS%val( iParticle ), &
            & interpolator     = particleGroup%interpolator,                   & 
            & scheme           = scheme,                                       &
            & geometry         = geometry,                                     &
            & params           = params,                                       &
            & dt_DEM_lat       = dt_sub_lat                                    )
      end do ! iParticle
      !$OMP END DO
      !$OMP END PARALLEL
      call tem_stopTimer(timerHandle = mus_particle_timerHandles%incrementAuxFieldTimer)
  
      call tem_stopTimer(timerHandle = mus_particle_timerHandles%forceTimer)
  
      call tem_startTimer(timerHandle = mus_particle_timerHandles%velocityTimer)
      do iParticle = 1, particleGroup%particles_DPS%nvals
        if(particleGroup%particles_DPS%val(iParticle)%owner &
          &                         == params%general%proc%rank) then
  
            ! Only update velocities of particle for which I am owner
            call updateVelocityVerlet(particleGroup%particles_DPS%val(iParticle), dt_sub)
  
        end if
      end do
  
      ! --- SYNCHRONIZE VELOCITIES --- !
      ! * Owner sends updated velocities to all other processes on which
      !   the particle exists
      call tem_startTimer(timerHandle = mus_particle_timerHandles%exchangeVelocitiesTimer)
      call exchangeVelocities_DPS( this = particleGroup,              &
                                 & send = particleGroup%send,         &
                                 & recv = particleGroup%recv,         &
                                 & comm = params%general%proc%comm,   &
                                 & myRank = params%general%proc%rank, &
                                 & message_flag = 1                   )
      call tem_stopTimer(timerHandle = mus_particle_timerHandles%exchangeVelocitiesTimer)
      ! --- SWAP SUBCYCLE FORCE BUFFER --- !
      do iParticle = 1, particleGroup%particles_DPS%nvals
        call DEM_swapFnowFnext( this = particleGroup%particles_DPS%val(iParticle) )
      end do
  
      call tem_stopTimer(timerHandle = mus_particle_timerHandles%velocityTimer)
    end do ! ksub
  
    !$OMP PARALLEL
    !$OMP DO
    do iParticle = 1, particleGroup%particles_DPS%nvals
  
      call recalculate_auxField_DPS(                                             &
              & particle         = particleGroup%particles_DPS%val( iParticle ), &
              & interpolator     = particleGroup%interpolator,                   & 
              & scheme           = scheme,                                       &
              & geometry         = geometry,                                     &
              & params           = params                                        )
    end do ! iParticle
    !$OMP END DO
    !$OMP END PARALLEL
    
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%subcycleTimer)
  
  end subroutine DEMSubcycles_DPS
  
  !> DEMsubcycles runs Nsubcycles time integration steps per LBM time step 
  !! In this routine velocity verlet integration of particle position and
  !! velocity is used. This routine is for one-way coupled DPS particles. 
  subroutine DEMSubcycles_DPS_onewaycoupled( particleGroup, scheme, geometry, params, Nsubcycles )
    !> Array of particles
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Number of subcycles
    integer, intent(in) :: Nsubcycles
    ! ------------------------------------------!
    real(kind=rk) :: dt, dx       ! LBM time step
    real(kind=rk) :: dt_sub       ! subcycle time step
    real(kind=rk) :: nu, nu_lat   ! dynamic viscosity for computing hydro forces
    integer :: ksub               ! iterator over subcycle loops
    integer :: iParticle
    integer :: myRank             ! This process rank
    integer :: lev
  
    real(kind=rk) :: dt_sub_lat
    
    ! ------------------------------------------!
    lev = geometry%tree%global%maxLevel
    dt = params%physics%dtLvl(lev)
    dx = params%physics%dxLvl(lev)
    myRank = params%general%proc%rank
    dt_sub = dt/Nsubcycles
    dt_sub_lat = 1.0_rk/Nsubcycles
  
    ! Get viscosity values
    nu_lat = scheme%field(1)%fieldProp%fluid%viscKine%dataOnLvl(lev)%val(1)
    nu = nu_lat * params%physics%fac(lev)%visc
  
  
    ! ---- INITIALIZE SUBCYCLING LOOP ---- !
  
    ! Fill neighbor list for efficient collision detection during subcycles
    if(particleGroup%enableCollisions) then
      call DEM_fillNeighborList( particleGroup = particleGroup, &
                              & d0            = 3*dx           )
    end if
  
    ! Calculate initial forces on particles
    call updateParticleForces( particleGroup = particleGroup, &
                             & scheme        = scheme,        &
                             & geometry      = geometry,      &
                             & params        = params,        &
                             & boundaryData  = pgBndData      )
  
  
    ! Swap force buffer for subcycle (collision) forces
    do iParticle = 1, particleGroup%particles_DPS%nvals
      call DEM_swapFnowFnext( this = particleGroup%particles_DPS%val(iParticle) )
    end do
  
    ! --- START OF ACTUAL SUBCYCLING LOOP --- !
    do ksub = 1, Nsubcycles
  
      ! --- POSITION UPDATE --- !
      do iParticle = 1, particleGroup%particles_DPS%nvals
        if(particleGroup%particles_DPS%val(iParticle)%owner &
          &                         == params%general%proc%rank) then
          ! Owner updates the position of the particle
          call updatePositionVerlet( & 
                        & this = particleGroup%particles_DPS%val(iParticle), &
                        & dt   = dt_sub                                      )
  
          ! Apply the effect of periodic boundaries (if used) on the position
          call wrapPeriodicPos( &
            & pos =  particleGroup%particles_DPS%val(iParticle)%pos(1:3), &
            & boundaryData = pgBndData                                    )
        end if
      end do
  
      call exchangePositions_DPS( this = particleGroup,              &
                                & send = particleGroup%send,         &
                                & recv = particleGroup%recv,         &
                                & comm = params%general%proc%comm,   &
                                & myRank = params%general%proc%rank, &
                                & message_flag = 1                   )
  
      ! --- MAP PARTICLES TO LATTICE --- !
      ! And add/remove particles to/from processes as needed
      call mapParticlesToLattice_DPS( particleGroup = particleGroup, &
                                    & scheme        = scheme,        &
                                    & geometry      = geometry,      &
                                    & params        = params         )
  
      ! Update forces (DEM and hydrodynamic) on all particles and 
      ! synchronize across processes
      call updateParticleForces( particleGroup = particleGroup, &
                               & scheme        = scheme,        &
                               & geometry      = geometry,      &
                               & params        = params,        &
                               & boundaryData  = pgBndData      )
  
  
      do iParticle = 1, particleGroup%particles_DPS%nvals
        if(particleGroup%particles_DPS%val(iParticle)%owner &
          &                         == params%general%proc%rank) then
  
            ! Only update velocities of particle for which I am owner
            call updateVelocityVerlet(                             & 
              & this = particleGroup%particles_DPS%val(iParticle), &
              & dt   = dt_sub                                      )
  
        end if
      end do
  
      ! --- SYNCHRONIZE VELOCITIES --- !
      ! * Owner sends updated velocities to all other processes on which
      !   the particle exists
      call exchangeVelocities_DPS( this = particleGroup,              &
                                 & send = particleGroup%send,         &
                                 & recv = particleGroup%recv,         &
                                 & comm = params%general%proc%comm,   &
                                 & myRank = params%general%proc%rank, &
                                 & message_flag = 1                   )
  
      ! --- SWAP SUBCYCLE FORCE BUFFER --- !
      do iParticle = 1, particleGroup%particles_DPS%nvals
        call DEM_swapFnowFnext(this = particleGroup%particles_DPS%val(iParticle))
      end do
    end do ! ksub
  end subroutine DEMSubcycles_DPS_onewaycoupled
  
  ! ---------------- TIME INTEGRATION ROUTINES ----------------- !
  !> Update particle continuous position using Euler integration 
  subroutine updatePositionEuler_MEM(this, dt)
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Time step
    real(kind=rk), intent(in) :: dt
    ! ------------------------------------------- !
    ! Store old position for handling possible collisions later
    this%oldPos = this%pos
  
    ! Update continuous position
    this%pos = this%pos + dt*this%vel
  
  end subroutine updatePositionEuler_MEM
  
  !> Update particle continuous position using Euler integration 
  subroutine updatePositionEuler_DPS(this, dt)
    type(mus_particle_DPS_type), intent(inout) :: this
    !> Time step
    real(kind=rk), intent(in) :: dt
    ! ------------------------------------------- !
    ! Store old position for handling possible collisions later
    this%oldPos = this%pos
  
    ! Update continuous position
    this%pos = this%pos + dt*this%vel
  
  end subroutine updatePositionEuler_DPS
  
  !> Update particle velocity according to current forces on particle
  !! using Euler integration
  subroutine updateVelocityEuler_MEM(this, dt)
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Time step
    real(kind=rk), intent(in) :: dt
    ! ------------------------------------------!
    integer :: now, next
  
    ! ------------------------------------------!
    now = this%F_DEM_now
    next = this%F_DEM_next
  
  
    this%vel(1:3) = this%vel(1:3) &
      &           + dt * (this%F_DEM(now,1:3) + this%F_DEM(next,1:3) ) / this%mass
  
    this%vel(4:6) = this%vel(4:6) &
    &             + dt * ( this%F_DEM(now,4:6) + this%F_DEM(next,4:6) ) / this%rotInertia
  
  end subroutine updateVelocityEuler_MEM
  
  !> Update particle velocity according to current forces on particle
  !! using Euler integration
  subroutine updateVelocityEuler_DPS(this, dt)
    type(mus_particle_DPS_type), intent(inout) :: this
    !> Time step
    real(kind=rk), intent(in) :: dt
    ! ------------------------------------------!
    integer :: now, next
  
    ! ------------------------------------------!
    now = this%F_DEM_now
    next = this%F_DEM_next
  
  
    this%vel(1:3) = this%vel(1:3) &
      &           + dt * 0.5*(this%F_DEM(now,1:3) + this%F_DEM(next,1:3) ) / this%mass
  
    this%vel(4:6) = this%vel(4:6) &
    &             + dt * 0.5*( this%F_DEM(now,4:6) + this%F_DEM(next,4:6) ) / this%rotInertia
  
  end subroutine updateVelocityEuler_DPS
  
  !> Update particle position using Verlet integration
  subroutine updatePositionVerlet_MEM( this, dt )
    type(mus_particle_MEM_type), intent(inout) :: this
    real(kind=rk),intent(in) :: dt
    ! -----------------------------------------!
    ! Store old position
    this%oldPos = this%pos
  
    this%pos(1:3) = this%pos(1:3) + dt*this%vel(1:3) &
      & + 0.5*dt**2 * ( this%F_DEM(this%F_DEM_now,1:3) )/this%mass
  
    this%pos(4:6) = this%pos(4:6) + dt*this%vel(4:6) &
      & + 0.5*dt**2 * this%F_DEM(this%F_DEM_now,4:6)/this%rotInertia
  
  end subroutine updatePositionVerlet_MEM
  
  !> Update particle position using Verlet integration
  subroutine updatePositionVerlet_DPS( this, dt )
    type(mus_particle_DPS_type), intent(inout) :: this
    real(kind=rk),intent(in) :: dt
    ! -----------------------------------------!
    real(kind=rk) :: Ftot(6)
    ! -----------------------------------------!
    ! Store old position
    this%oldPos = this%pos
  
    Ftot(1:6) = 0.5*(this%F_DEM(this%F_DEM_next,1:6) &
              & + this%F_DEM(this%F_DEM_now,1:6) )   &
              & + this%Fext(1:6)                     &
              & + this%F(1:6)
  
    this%pos(1:3) = this%pos(1:3) + dt*this%vel(1:3) &
      & + 0.5*dt**2 * Ftot(1:3)/this%mass
  
    this%pos(4:6) = this%pos(4:6) + dt*this%vel(4:6) &
      & + 0.5*dt**2 * Ftot(4:6)/this%rotInertia
  
  end subroutine updatePositionVerlet_DPS
  
  !> Update particle velocity using Verlet integration
  subroutine updateVelocityVerlet_MEM( this, dt )
    type(mus_particle_MEM_type), intent(inout) :: this
    real(kind=rk),intent(in) :: dt
    ! -----------------------------------------!
  
    this%vel(1:3) = this%vel(1:3) &
      & + 0.5*dt*( this%F_DEM(this%F_DEM_now,1:3) &
      &          + this%F_DEM(this%F_DEM_next,1:3) )/this%mass
  
    this%vel(4:6) = this%vel(4:6) &
      & + 0.5*dt*( this%F_DEM(this%F_DEM_now,4:6) &
      &          + this%F_DEM(this%F_DEM_next,4:6) )/this%rotInertia
  
  end subroutine updateVelocityVerlet_MEM
  
  !> Update particle velocity using Verlet integration
  subroutine updateVelocityVerlet_DPS( this, dt )
    type(mus_particle_DPS_type), intent(inout) :: this
    real(kind=rk),intent(in) :: dt
    ! -----------------------------------------!
    real(kind=rk) :: infty
    real(kind=rk) :: Ftot(6)
    integer :: now, next
    ! -----------------------------------------!
    now = this%F_DEM_now
    next = this%F_DEM_next
  
    infty = huge(infty)
  
    ! Check for infty forces
    if( any( this%F_DEM(this%F_DEM_now,1:6) > infty )&
      & .OR. any( this%F_DEM(this%F_DEM_next,1:6) > infty ) ) then
      write(logUnit(1),*) "ERROR updateVelocityVerlet: infty force"
      call tem_abort()
    end if
  
    Ftot(1:6) = 0.5*(this%F_DEM(next,1:6) + this%F_DEM(now,1:6) ) &
              & + this%Fext(1:6) &
              & + this%F(1:6)
  
    this%vel(1:3) = this%vel(1:3) + dt*Ftot(1:3) / this%mass
    this%vel(4:6) = this%vel(4:6) + dt*Ftot(4:6) / this%rotInertia
  
  end subroutine updateVelocityVerlet_DPS
  
  ! updateParticleForces updates forces (DEM and hydrodynamic) on all particles 
  ! and synchronizes these across processes
  subroutine updateParticleForces( particleGroup, scheme, geometry, &
                                 & params, boundaryData             )
    !> Array of particles
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Domain boundary information for particle-wall collisions
    type(mus_particle_boundarydata_type) :: boundaryData
    ! ------------------------------------------ !
    real(kind=rk) :: collision_tol, collision_time
    real(kind=rk) :: nu_lat, nu, mu
    integer :: myRank, lev
    ! ------------------------------------------ !
    lev = geometry%tree%global%maxlevel
    collision_tol = particleGroup%collision_tol
    collision_time = particleGroup%collision_time
    myRank =  params%general%proc%rank
  
    ! Get viscosity values
    nu_lat = scheme%field(1)%fieldProp%fluid%viscKine%dataOnLvl(lev)%val(1)
    nu = nu_lat * params%physics%fac(lev)%visc
  
    ! Note dynamic viscosity in lattice units is same as kinematic visc 
    ! in lattice units
    mu = nu_lat * params%physics%fac(lev)%viscDyna
  
    ! --- DEM FORCE COMPUTATION --- !
    ! Reset particle%F_DEM(next,:) = 0.0 for all particles
    call DEM_resetFnext_DPS( particleGroup = particleGroup )
  
    if(particleGroup%enableCollisions) then
      ! Reset the force buffers which will be filled in 
      ! the DEM collision routines for particles and walls 
      call resetForceBuffers( particleGroup )
  
    call tem_startTimer(timerHandle = mus_particle_timerHandles%collisionHandlingTimer)
    call DEM_computeLocalCollisionForces_DPS(                          &
                            & particleGroup = particleGroup,           &
                            & myRank        = myRank,                  &
                            & eps           = collision_tol,           &
                            & Tc            = collision_time,          &
                            & mu            = mu                       )
  
      ! Compute collision forces where one of the particles is remote
      call DEM_computeRemoteCollisionForces_DPS( &
                                          & particleGroup = particleGroup,  &
                                          & myRank        = myRank,         &
                                          & eps           = collision_tol,  &
                                          & Tc            = collision_time, &
                                          & mu            = mu              )
    end if
  
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%collisionHandlingTimer)
    ! Add external forces (hydrodynamic and body forces) to F_DEM(next,:)
    ! NOTE: These are not added to force buffer
  
    ! --- INTERPOLATE FLUID PROPERTIES --- !
    call tem_startTimer(timerHandle = mus_particle_timerHandles%interpolateFluidTimer)
    call interpolateFluidProperties_DPS( particleGroup = particleGroup, &
                                       & scheme        = scheme,        &
                                       & geometry      = geometry,      &
                                       & params        = params         )
    call tem_stopTimer(timerHandle = mus_particle_timerHandles%interpolateFluidTimer)
  
  
    call computeHydroForces_DPS( particleGroup = particleGroup, &
                               & nu            = nu,            &
                               & myRank        = myRank         )
  
    if( boundaryData%useBnd ) then
      call DEM_computeWallForces_DPS( particleGroup = particleGroup,  &
                                    & boundaryData  = boundaryData,   & 
                                    & eps           = collision_tol,  &
                                    & Tc            = collision_time, &
                                    & myRank        = myRank          )
    end if
  
  
    ! Communicate forces: owner of the particle gets all DEM force contributions
    ! from other processes
    if(particleGroup%enableCollisions) then
    call tem_startTimer( timerHandle = mus_particle_timerHandles%exchangeDEMForcesTimer )
      call DEM_exchangeForces_DPS( this         = particleGroup,           &
                                & send         = particleGroup%send,       &
                                & recv         = particleGroup%recv,       &
                                & comm         = params%general%proc%comm, &
                                & myRank       = myRank,                   &
                                & message_flag = 1                         )
    call tem_stopTimer( timerHandle = mus_particle_timerHandles%exchangeDEMForcesTimer )
    end if
  
    ! Synchronize hydrodynamic forces across processes
    call tem_startTimer( timerHandle = mus_particle_timerHandles%exchangeHydroForcesTimer )
    call exchangeHydroforces_DPS( this = particleGroup,              &
                                & send = particleGroup%send,         &
                                & recv = particleGroup%recv,         &
                                & comm = params%general%proc%comm,   &
                                & myRank = params%general%proc%rank, &
                                & message_flag = 1                   )
    call tem_stopTimer( timerHandle = mus_particle_timerHandles%exchangeHydroForcesTimer )
    
  end subroutine updateParticleForces
  
  ! ------ ROUTINES FOR COMPUTING COLLISION FORCES BETWEEN PARTICLES ---------- !
  !> Compute collision forces using Discrete Element Method between particles 
  !! local to this process i.e. both particles are owned by this process.
  subroutine DEM_computeLocalCollisionForces(particleGroup, myRank, eps, Tc, mu)
    !> particleGroup to search for collisions in
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> This proc's rank
    integer, intent(in) :: myRank
    !> Threshold gap at which to call it a collision
    real(kind=rk), intent(in) :: eps
    !> DEM collision time, used for calculating spring and damper coefficients
    real(kind=rk), intent(in) :: Tc
    !> Fluid dynamic viscosity (physical units)
    real(kind=rk), intent(in) :: mu
    ! ------------------------------------------------------ !
    integer :: iParticle, jParticle
    logical :: collision
    real(kind=rk) :: Fcoll(3)
    ! ------------------------------------------------------ !
    Fcoll = 0.0_rk
  
    ! reset hasCollided
    do iParticle = 1, particleGroup%particles_MEM%nvals
      particleGroup%particles_MEM%val(iParticle)%hasCollided = .FALSE.
    end do
  
    ! First handle LOCAL collisions, where I own both of the particles
    do iParticle = 1, particleGroup%particles_MEM%nvals
      do jParticle = iParticle + 1, particleGroup%particles_MEM%nvals
  
        ! Check if the collision should be resolved on this process
        if( isLocalCollision( particleGroup%particles_MEM%val(iParticle), &
                            & particleGroup%particles_MEM%val(jParticle), &
                            & myRank                               )  ) then
          collision = .FALSE.
  
          call checkAndCollideDEM(                                          & 
                  & particleA = particleGroup%particles_MEM%val(iParticle), &
                  & particleB = particleGroup%particles_MEM%val(jParticle), &
                  & hasCollided = collision,                                &
                  & eps       = eps,                                        &
                  & Tc        = Tc,                                         &
                  & mu        = mu                                          )
  
  
          ! Set particle%hasCollided to true if collision has taken place
          if(collision) then
            ! write(logUnit(1),*) "Local collision"
            particleGroup%particles_MEM%val(iParticle)%hasCollided = .TRUE.
            particleGroup%particles_MEM%val(jParticle)%hasCollided = .TRUE.
          end if ! hasCollided
        end if ! collideOnThisProc
      end do ! jParticle
    end do ! iParticle
  end subroutine DEM_computeLocalCollisionForces
  
  !> Compute collision forces using Discrete Element Method between particles 
  !! local to this process i.e. both particles are owned by this process.
  subroutine DEM_computeLocalCollisionForces_DPS( particleGroup, myRank, &
                                                & eps, Tc, mu            )
    !> particleGroup to search for collisions in
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> This proc's rank
    integer, intent(in) :: myRank
    !> Threshold gap at which to call it a collision
    real(kind=rk), intent(in) :: eps
    !> DEM collision time, used for calculating spring and damper coefficients
    real(kind=rk), intent(in) :: Tc
    !> Fluid dynamic viscosity, in physical units
    real(kind=rk), intent(in) :: mu
    ! ------------------------------------------------------ !
    integer :: iParticle, jParticle
    integer :: iNgh
    logical :: collision
    real(kind=rk) :: Fcoll(3)
    ! ------------------------------------------------------ !
    Fcoll = 0.0_rk
  
    ! reset hasCollided
    do iParticle = 1, particleGroup%particles_DPS%nvals
      particleGroup%particles_DPS%val(iParticle)%hasCollided = .FALSE.
    end do
  
    ! First handle LOCAL collisions, where I own both of the particles
    do iParticle = 1, particleGroup%particles_DPS%nvals
       do iNgh = 1, particleGroup%particles_DPS%val(iParticle)%DEM_neighborList%nvals
        jParticle = particleGroup%particles_DPS%val(iParticle)%DEM_neighborList%val(iNgh)
  
  
        ! Check if the collision should be resolved on this process
        if( isLocalCollision( particleGroup%particles_DPS%val(iParticle), &
                            & particleGroup%particles_DPS%val(jParticle), &
                            & myRank                               )  ) then
          collision = .FALSE.
  
          call checkAndCollideDEM(                                              &
                      & particleA = particleGroup%particles_DPS%val(iParticle), &
                      & particleB = particleGroup%particles_DPS%val(jParticle), &
                      & hasCollided = collision,                                &
                      & eps       = eps,                                        &
                      & Tc        = Tc,                                         &
                      & mu        = mu                                          )
  
  
          ! Set particle%hasCollided to true if collision has taken place
          if(collision) then
            ! write(logUnit(1),*) "Local collision"
            particleGroup%particles_DPS%val(iParticle)%hasCollided = .TRUE.
            particleGroup%particles_DPS%val(jParticle)%hasCollided = .TRUE.
          end if ! hasCollided
        end if ! collideOnThisProc
      end do ! jParticle
    end do ! iParticle
  end subroutine DEM_computeLocalCollisionForces_DPS
  
  
  !> Compute collision forces using Discrete Element Method between particle
  !! pairs where particles belong to different processes.
  subroutine DEM_computeRemoteCollisionForces(particleGroup, myRank, eps, Tc, mu)
    !> particleGroup to search for collisions in
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> This proc's rank
    integer, intent(in) :: myRank
    !> Threshold gap at which to call it a collision
    real(kind=rk), intent(in) :: eps
    !> DEM collision time, used for calculating spring and damper coefficients
    real(kind=rk), intent(in) :: Tc
    !> Fluid dynamic viscosity, physical units
    real(kind=rk), intent(in) :: mu
    ! ------------------------------------------------------ !
    integer :: iParticle, jParticle
    integer :: otherRank             ! rank of the particle owned by other proc
    integer :: otherRankIndex        ! index of otherRank in send%proc
    integer :: iproc
    logical :: isRemoteCollision
    logical :: collision
    real(kind=rk) :: Fcoll(3)        ! Collision force
    ! ------------------------------------------------------ !
    Fcoll = 0.0_rk
    otherRank = -1
    isRemoteCollision = .FALSE.
  
    ! reset hasCollided
    do iParticle = 1, particleGroup%particles_MEM%nvals
      particleGroup%particles_MEM%val(iParticle)%hasCollided = .FALSE.
    end do
  
    ! Clear the send buffer for DEM forces to be sent to other processes
    do iproc = 1, particleGroup%send%nProcs
      particleGroup%send%buf_force(iproc)%nParticles = 0
    end do ! iproc
  
    ! First handle LOCAL collisions, where I own both of the particles
    do iParticle = 1, particleGroup%particles_MEM%nvals
      do jParticle = iParticle + 1, particleGroup%particles_MEM%nvals
  
        ! Check if the collision should be resolved on this process
        call DEM_isRemoteCollision( particleGroup%particles_MEM%val(iParticle), &
                                  & particleGroup%particles_MEM%val(jParticle), &
                                  & myRank,                                 &
                                  & particleGroup%send,                     &
                                  & isRemoteCollision,                      &
                                  & otherRank,                              &
                                  & otherRankIndex                          )
  
        if(isRemoteCollision) then
          collision = .FALSE.
  
          ! Check to see if iParticle and jParticle collide. If so, return the 
          ! collision force on iParticle in Fcoll.
          call checkAndCollideDEM(                                          &
                  & particleA = particleGroup%particles_MEM%val(iParticle), &
                  & particleB = particleGroup%particles_MEM%val(jParticle), &
                  & hasCollided = collision,                                &
                  & eps       = eps,                                        &
                  & Tc        = Tc,                                         &
                  & mu        = mu                                          )
  
          ! Set particle%hasCollided to true if collision has taken place
          if(collision) then
            ! Send over the force on the remotely-owned particle.
  
            ! If jParticle is the remotely-owned particle, calculate the force 
            ! on it (force on jParticle = - force on iParticle)
            if( particleGroup%particles_MEM%val(jParticle)%owner == otherRank ) then
              call addCollisionForceToBuffer(                                           &
                & particleGroup = particleGroup,                                        &
                & recvRankIndex = otherRankIndex,                                       &
                & Fcoll         = -Fcoll,                                               &
                & particleID    = particleGroup%particles_MEM%val(jParticle)%particleID )
  
            else if( particleGroup%particles_MEM%val(iParticle)%owner == otherRank ) then
              ! iParticle is the remotely-owned particle
              call addCollisionForceToBuffer(                                           &
                & particleGroup = particleGroup,                                        &
                & recvRankIndex = otherRankIndex,                                       &
                & Fcoll         = Fcoll,                                                &
                & particleID    = particleGroup%particles_MEM%val(iParticle)%particleID )
            else
              write(logUnit(1),*) "ERROR DEM_computeRemoteCollisionForces_MEM: could not identify remote particle"
              call tem_abort()
            end if
  
            particleGroup%particles_MEM%val(iParticle)%hasCollided = .TRUE.
            particleGroup%particles_MEM%val(jParticle)%hasCollided = .TRUE.
          end if ! hasCollided
        end if ! collideOnThisProc
      end do ! jParticle
    end do ! iParticle
  
  end subroutine DEM_computeRemoteCollisionForces
  
  !> Compute collision forces using Discrete Element Method between particle
  !! pairs where particles belong to different processes.
  subroutine DEM_computeRemoteCollisionForces_DPS( particleGroup, myRank, &
                                                 & eps, Tc, mu            )
    !> particleGroup to search for collisions in
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> This proc's rank
    integer, intent(in) :: myRank
    !> Threshold gap at which to call it a collision
    real(kind=rk), intent(in) :: eps
    !> DEM collision time, used for calculating spring and damper coefficients
    real(kind=rk), intent(in) :: Tc
    !> Fluid dynamic viscosity, physical units
    real(kind=rk), intent(in) :: mu
    ! ------------------------------------------------------ !
    integer :: iParticle, jParticle
    integer :: iNgh
    integer :: otherRank             ! rank of the particle owned by other proc
    integer :: otherRankIndex        ! index of otherRank in send%proc
    integer :: iproc
    logical :: isRemoteCollision
    logical :: collision
    real(kind=rk) :: Fcoll(3)        ! Collision force
    ! ------------------------------------------------------ !
    Fcoll = 0.0_rk
    otherRank = -1
    isRemoteCollision = .FALSE.
  
    ! reset hasCollided
    do iParticle = 1, particleGroup%particles_DPS%nvals
      particleGroup%particles_DPS%val(iParticle)%hasCollided = .FALSE.
    end do
  
    ! Clear the send buffer for DEM forces to be sent to other processes
    do iproc = 1, particleGroup%send%nProcs
      particleGroup%send%buf_force(iproc)%nParticles = 0
    end do ! iproc
  
    ! Loop over particle pairs to detect collisions
    do iParticle = 1, particleGroup%particles_DPS%nvals
      do iNgh = 1, particleGroup%particles_DPS%val(iParticle)%DEM_neighborList%nvals
        jParticle = particleGroup%particles_DPS%val(iParticle)%DEM_neighborList%val(iNgh)
  
        ! It is possible for the DEM_neighborlist to include particles which were 
        ! removed from the simulation during the DEM subcycles within the current 
        ! LBM time step. If this is the case the ID of such a particle will have 
        ! been set to the negative of what it was originally. If we encounter such 
        ! a particle we should skip it.
        if( particleGroup%particles_DPS%val(jParticle)%particleID <= 0 ) cycle
      
  
        ! Check if the collision should be resolved on this process
        call DEM_isRemoteCollision( particleGroup%particles_DPS%val(iParticle), &
                                  & particleGroup%particles_DPS%val(jParticle), &
                                  & myRank,                                 &
                                  & particleGroup%send,                     &
                                  & isRemoteCollision,                      &
                                  & otherRank,                              &
                                  & otherRankIndex                          )
  
        if(isRemoteCollision) then
          collision = .FALSE.
  
          call checkAndCollideDEM(                                            &
                  & particleA   = particleGroup%particles_DPS%val(iParticle), &
                  & particleB   = particleGroup%particles_DPS%val(jParticle), &
                  & hasCollided = collision,                                  &
                  & eps         = eps,                                        &
                  & Tc          = Tc,                                         &
                  & mu          = mu                                          )
  
  
          ! Set particle%hasCollided to true if collision has taken place
          if(collision) then
            ! Send over the force on the remotely-owned particle.
  
            ! If jParticle is the remotely-owned particle, calculate the force on it
            ! Force on jParticle = - force on iParticle
            if( particleGroup%particles_DPS%val(jParticle)%owner == otherRank ) then
              call addCollisionForceToBuffer(                                      &
                  & particleGroup = particleGroup,                                 &
                  & recvRankIndex = otherRankIndex,                                &
                  & Fcoll         = -Fcoll,                                        &
                  & particleID    = particleGroup%particles_DPS%val(jParticle)%particleID )
  
            else if( particleGroup%particles_DPS%val(iParticle)%owner == otherRank ) then
              ! iParticle is the remotely-owned particle
              call addCollisionForceToBuffer(                                               &
                    & particleGroup = particleGroup,                                        &
                    & recvRankIndex = otherRankIndex,                                       &
                    & Fcoll         = Fcoll,                                                &
                    & particleID    = particleGroup%particles_DPS%val(iParticle)%particleID )
            else
              write(logUnit(1),*) "ERROR DEM_computeRemoteCollisionForces_DPS: could not identify remote particle"
              call tem_abort()
            end if
  
            particleGroup%particles_DPS%val(iParticle)%hasCollided = .TRUE.
            particleGroup%particles_DPS%val(jParticle)%hasCollided = .TRUE.
          end if ! hasCollided
        end if ! collideOnThisProc
      end do ! jParticle
    end do ! iParticle
  end subroutine DEM_computeRemoteCollisionForces_DPS
  
  ! ----------- ROUTINES FOR ADDING EXTERNAL AND HYDRODYNAMIC FORCES ---------- !
  !> DEM_computeExternalForces adds the non-DEM forces (including
  !! hydrodynamic force from the LBM) to F_DEM which is used to
  !! update velocity within the DEM subcycle.
  subroutine DEM_computeExternalForces( particleGroup )
    !> particleGroup to search for collisions in
    type(mus_particle_group_type), intent(inout) :: particleGroup
    ! ------------------------------------------!
    integer :: iParticle
    integer :: F_DEM_next
    ! ------------------------------------------!
    do iParticle = 1, particleGroup%particles_MEM%nvals
      F_DEM_next = particleGroup%particles_MEM%val(iParticle)%F_DEM_next
  
      particleGroup%particles_MEM%val(iParticle)%F_DEM(F_DEM_next,1:6) &
        & = particleGroup%particles_MEM%val(iParticle)%F_DEM(F_DEM_next,1:6) &
        & + particleGroup%particles_MEM%val(iParticle)%F(1:6) &
        & + particleGroup%particles_MEM%val(iParticle)%Fext(1:6)
  
    end do
  end subroutine DEM_computeExternalForces
  
  !> DEM_computeExternalForces adds the non-DEM forces (including
  !! hydrodynamic force from the LBM) to F_DEM which is used to
  !! update velocity within the DEM subcycle.
  subroutine computeHydroForces_DPS( particleGroup, nu, myRank )
    !> particleGroup to search for collisions in
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Fluid dynamic viscosity (physical units) for computing hydro forces
    real(kind=rk), intent(in) :: nu
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------!
    integer :: iParticle
    real(kind=rk) :: Fd(3), Fp(3), Flift(3)
    real(kind=rk) :: eps_p
    ! ------------------------------------------!
    do iParticle = 1, particleGroup%particles_DPS%nvals
      if(  particleGroup%particles_DPS%val(iParticle)%owner == myRank ) then 
  
        ! Compute solid volume fraction
        eps_p = 1.0_rk - particleGroup%particles_DPS%val(iParticle)%eps_f_fluid
  
        call particleGroup%calcDragForce(                            &
            & particle = particleGroup%particles_DPS%val(iParticle), &
            & eps_p    = eps_p,                                      &
            & nu       = nu,                                         &
            & Fd       = Fd                                          )
  
        call particleGroup%calcPressureForce(                            &
            & particle = particleGroup%particles_DPS%val(iParticle),     &
            & Fp           = Fp                                          )
  
        call particleGroup%calcLiftForce(                            &
            & particle = particleGroup%particles_DPS%val(iParticle), &
            & nu       = nu,                                         &
            & Flift    = Flift                                       )
  
        ! particleGroup%particles_DPS%val(iParticle)%F(1:3) = Fd + Fp + Flift 
        ! particleGroup%particles_DPS%val(iParticle)%F(1:3) = Fd + Flift
        particleGroup%particles_DPS%val(iParticle)%F(1:3) = Fd 
  
      end if
    end do
  end subroutine computeHydroForces_DPS
  
  !> DEM_computeExternalForces adds the non-DEM forces (including
  !! hydrodynamic force from the LBM) to F_DEM which is used to
  !! update velocity within the DEM subcycle.
  subroutine DEM_computeExternalForces_DPS_oneway( particleGroup, nu, myRank )
    !> particleGroup to search for collisions in
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Fluid dynamic viscosity (physical units) for computing hydro forces
    real(kind=rk), intent(in) :: nu
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------!
    integer :: iParticle
    integer :: F_DEM_next
    real(kind=rk) :: Fd(3), Fp(3), Flift(3)
    ! ------------------------------------------!
    do iParticle = 1, particleGroup%particles_DPS%nvals
      if(  particleGroup%particles_DPS%val(iParticle)%owner == myRank ) then 
        F_DEM_next = particleGroup%particles_DPS%val(iParticle)%F_DEM_next
  
        ! Compute hydrodynamic forces
        call applyDragForce_DPS(                                     &
            & particle = particleGroup%particles_DPS%val(iParticle), &
            & eps_p    = 0.0_rk,                                     &
            & nu       = nu,                                         &
            & Fd       = Fd                                          )
  
        call applyPressureForce_DPS( & 
              & particle     = particleGroup%particles_DPS%val(iParticle), &
              & Fp           = Fp                                          )
  
        call applyLiftForce_DPS(                                     &
            & particle = particleGroup%particles_DPS%val(iParticle), &
            & nu       = nu,                                         &
            & Flift    = Flift                                       )
  
        ! particleGroup%particles_DPS%val(iParticle)%F(1:3) = Fd + Fp + Flift 
        particleGroup%particles_DPS%val(iParticle)%F(1:3) = Fd
  
        ! Compute total force (hydro + DEM) on particle and store in F_DEM
        particleGroup%particles_DPS%val(iParticle)%F_DEM(F_DEM_next,1:6) &
          & = particleGroup%particles_DPS%val(iParticle)%F_DEM(F_DEM_next,1:6) &
          & + particleGroup%particles_DPS%val(iParticle)%F(1:6)                &
          & + particleGroup%particles_DPS%val(iParticle)%Fext(1:6)
  
      end if
    end do
  end subroutine DEM_computeExternalForces_DPS_oneway
  
  ! -------- ROUTINES FOR INTERACTIONS BETWEEN PARTICLES AND WALLS ----------- !
  
  !> DEM_computeWallForces_DPS computes the forces on particles as a result of 
  !! collisions with walls described in the boundaryData object. These forces 
  !! are only computed for particles owned by this process.
  subroutine DEM_computeWallForces_DPS( particleGroup, boundaryData, eps, &
                                      & Tc, myRank                        )
    !> particleGroup to search for collisions in
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Data about particle domain boundaries
    type(mus_particle_boundaryData_type), intent(in) :: boundaryData
    !> threshold at which gap to call it a collision
    real(kind=rk), intent(in) :: eps
    !> Collision time
    real(kind=rk), intent(in) :: Tc
    !> This process rank
    integer, intent(in) :: myRank
    ! ------------------------------------------!
    integer :: iParticle
    integer :: F_DEM_next
    real(kind=rk) :: Fwall(3)
    !> particle mass, spring constant, damping coefficient
    real(kind=rk) :: mp, kn, dn
    !> dry resitution coefficient
    real(kind=rk) :: e_dry
    ! ------------------------------------------!
    e_dry = 1.0_rk
  
    do iParticle = 1, particleGroup%particles_DPS%nvals
      if(  particleGroup%particles_DPS%val(iParticle)%owner == myRank ) then 
        F_DEM_next = particleGroup%particles_DPS%val(iParticle)%F_DEM_next
  
        ! Compute spring and damper coefficient
        mp = particleGroup%particles_DPS%val(iParticle)%mass
        kn = mp * ( PI**2 + (log(e_dry) )**2 )/Tc**2
        dn = (-2.0*mp*log(e_dry))/Tc
  
        Fwall = computeWallForce_DPS(                                  &
                  & this = particleGroup%particles_DPS%val(iParticle), &
                  & boundaryData = boundaryData,                       &
                  & eps          = eps,                                &
                  & kn           = kn,                                 &
                  & dn           = dn                                  )
  
  
        ! Add wall interaction force to total F_DEM
        particleGroup%particles_DPS%val(iParticle)%F_DEM(F_DEM_next,1:3) &
          & = particleGroup%particles_DPS%val(iParticle)%F_DEM(F_DEM_next,1:3) &
          & + Fwall(1:3)
  
      end if
    end do
  end subroutine DEM_computeWallForces_DPS
  
  
  
  !> Routine for computing the collision force between particles and walls for the 
  !! case of a simple prismatic domain with walls aligned with the Cartesian axes.
  function computeWallForce_DPS( this, boundaryData, eps, kn, dn) result ( Fwall )
    !> Particle to collide
    type(mus_particle_DPS_type), intent(inout) :: this
    !> Boundary data containing information on wall location
    type(mus_particle_boundaryData_type) :: boundaryData
    !> Threshold gap at which to call it a collision
    real(kind=rk), intent(in) :: eps
    !> Spring coefficient
    real(kind=rk), intent(in) :: kn
    !> Damping coefficient
    real(kind=rk), intent(in) :: dn
    !> Output: collision force on particleA = negative of collision force on particleB
    real(kind=rk) :: Fwall(3)
    ! ------------------------------------------!
    real(kind=rk) :: Fwall_1D
    real(kind=rk) :: xwall
    integer :: iDir, iBnd, i
    ! ------------------------------------------!
  
    ! For each wall in boundaryData, compute the force from that wall 
    ! force = 0 if no collision
    Fwall = 0.0_rk
  
    do iDir = 1,3 
      do i = 1,2
        iBnd = i + (iDir-1)*2
  
          if( boundaryData%wallBnd(iBnd) ) then
            xwall = boundaryData%bnd(iBnd)
  
            Fwall_1D = computeWallForce_1D( xp    = this%pos(iDir), &
                                          & up    = this%vel(iDir), &
                                          & Rp    = this%radius,    &
                                          & xwall = xwall,          &
                                          & kn    = kn,             &
                                          & dn    = dn,             &
                                          & eps   = eps             )
          else
            Fwall_1D = 0.0_rk
          end if
  
          Fwall(iDir) = Fwall(iDir) + Fwall_1D
        end do
    end do
     
    
  end function computeWallForce_DPS
  
  ! -------- ROUTINES FOR FAST COLLISION DETECTION BETWEEN PARTICLES ----------- !
  subroutine DEM_fillNeighborList( particleGroup, d0 )
    !> particleGroup
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Radius of search area: particles closer than this distance will be added
    !! to the neighbor list
    real(kind=rk), intent(in) :: d0
    ! ------------------------------------------!
    !> Positions of two particles we are checking
    real(kind=rk) :: xpi(3), xpj(3), rij(3)
    !> Square of distance between particles
    real(kind=rk) :: d2
    integer :: iParticle, jParticle
    ! ------------------------------------------!
    do iParticle = 1, particleGroup%particles_DPS%nvals
      !  First clear iParticle's neighborList
      call empty( me = particleGroup%particles_DPS%val(iParticle)%DEM_neighborList )
  
      ! Now fill iParticle's DEM_neighborlist with new neighbors
      do jParticle = iParticle + 1, particleGroup%particles_DPS%nvals
        xpi = particleGroup%particles_DPS%val(iParticle)%pos(1:3)
        xpj = particleGroup%particles_DPS%val(jParticle)%pos(1:3)
        ! Compute distance between these two particles, taking periodic boundaries 
        ! into account if present.
        rij = computeDisplacement( x1           = xpi,      &
                                 & x2           = xpj,      &
                                 & boundaryData = pgBndData )
  
        d2 = dot_product(rij,rij) 
        if( d2 < d0**2 ) then
          ! Add jParticle to iParticle's neighbor list
          call append( me  = particleGroup%particles_DPS%val(iParticle)%DEM_neighborList, &
                     & val = jParticle                                                    )
        end if
        
      end do
    end do
    
  end subroutine DEM_fillNeighborList
  
  
  
  ! Code below is old stuff that is not currently used
  ! These routines compute wall interactions based on "average" wall 
  ! positions calculated from the treeID's. Works OK for MEM 
  ! if dx << particle diameter. Works horrible for DPS. 
  subroutine DEM_storeWallPositions_MEM( particleGroup, scheme, geometry, &
                                       & params, Tc, eps                  )
    !> particleGroup
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Collision time
    real(kind=rk) :: Tc
    !> Threshold gap for DEM collisions
    real(kind=rk) :: eps
    ! ------------------------------------------!
    integer :: iParticle, iproc
    real(kind=rk) :: wallPosSum(3) ! vector from center of particle to wall
  
    integer :: lev, myRank
    integer :: owner
    logical :: foundWall, rmflag
    integer :: nWallPos
  
    ! ------------------------------------------!
    myRank = params%general%proc%rank
  
    ! First clear wall position buffers
    do iproc = 1, particleGroup%send%nProcs
      particleGroup%send%buf_wall(iproc)%nParticles = 0
    end do
  
    do iParticle = 1, particleGroup%particles_MEM%nvals
      lev = particleGroup%particles_MEM%val(iParticle)%coordOfOrigin(4)
      rmflag = .FALSE.
      foundWall = .FALSE.
      nWallPos = 0
      wallPosSum = 0.0_rk
  
      ! For each particle, compute the LOCAL wallPosSum = the sum of 
      ! the barycenter coordinates of all lattice sites in the neighborhood
      ! of the particle
      call computeWallPosSum( this = particleGroup%particles_MEM%val(iParticle), &
                            & BCinteraction = particleGroup%BC_interaction,  &
                            & scheme = scheme,                               &
                            & stencil = scheme%layout%fStencil,              &
                            & geometry = geometry,                           &
                            & params = params,                               &
                            & rmflag = rmflag,                               &
                            & foundWall = foundWall,                         &
                            & wallPosSum = wallPosSum,                       &
                            & nWallPos   = nWallPos                          )
      if(foundWall) then
        particleGroup%particles_MEM%val(iParticle)%interactWithWall = .TRUE.
  
        ! Add wallPosSum to buffer for particle owner
        owner = particleGroup%particles_MEM%val(iParticle)%owner
  
        if( owner == params%general%proc%rank ) then
          particleGroup%particles_MEM%val(iParticle)%avgWallPos(1:3) = wallPosSum
          particleGroup%particles_MEM%val(iParticle)%nWallPos = nWallPos
        else
          do iproc = 1, particleGroup%send%nProcs
            if (particleGroup%send%proc(iproc) == owner ) then
              ! Check if wall buffer is large enough
              if( particleGroup%send%buf_wall(iproc)%nParticles+1 &
                &     >= particleGroup%send%buf_wall(iproc)%maxParticles  ) then
                call tem_abort('Error DEM_storeWallPositions: wall buffer too small.')
              end if
  
              ! If wall buffer is large enough, add wallPosSum
              particleGroup%send%buf_wall(iproc)%nParticles &
                & = particleGroup%send%buf_wall(iproc)%nParticles+1
  
              particleGroup%send%buf_wall(iproc)%val( & 
                & particleGroup%send%buf_wall(iproc)%nParticles )%V(1:3) = wallPosSum(1:3)
  
              particleGroup%send%buf_wall(iproc)%val( &
                & particleGroup%send%buf_wall(iproc)%nParticles )%I(1) &
                  & = nWallPos 
              
              particleGroup%send%buf_wall(iproc)%val( & 
                & particleGroup%send%buf_wall(iproc)%nParticles )%I(2) &
                  & = particleGroup%particles_MEM%val(iParticle)%particleID
            end if
          end do
        end if ! owner /= myRank
      else
        particleGroup%particles_MEM%val(iParticle)%interactWithWall = .FALSE.
      end if !foundWall
  
      if( rmflag ) then
        ! particleGroup%particles%val(iParticle)%removeParticle_global = .TRUE.
        particleGroup%particles_MEM%val(iParticle)%removeParticle_global = .FALSE.
      else
        particleGroup%particles_MEM%val(iParticle)%removeParticle_global = .FALSE.
      end if
    end do
  
  
    call DEM_exchangeWallPositions( this = particleGroup,            &
                                  & send = particleGroup%send,       &
                                  & recv = particleGroup%recv,       &
                                  & comm = params%general%proc%comm, &
                                  & myRank = myRank,                 &
                                  & message_flag = 1                 )
  
  end subroutine DEM_storeWallPositions_MEM
  
  ! ------------- MISCELLANEOUS ROUTINES ------------------!
  
  !> Swap the indices of Fnow and Fnext pointing to the 
  !! force at the current and next time step respectively
  subroutine DEM_swapFnowFnext_MEM(this)
    type(mus_particle_MEM_type), intent(inout) :: this
    ! ------------------------------------------!
    ! Set Fnow for next time step
    if(this%F_DEM_now == 1) then
      this%F_DEM_now = 2
      this%F_DEM_next = 1
    else
      this%F_DEM_now = 1
      this%F_DEM_next = 2
    end if
  end subroutine DEM_swapFnowFnext_MEM
  
  !> Swap the indices of Fnow and Fnext pointing to the 
  !! force at the current and next time step respectively
  subroutine DEM_swapFnowFnext_DPS(this)
    type(mus_particle_DPS_type), intent(inout) :: this
    ! ------------------------------------------!
    ! Set Fnow for next time step
    if(this%F_DEM_now == 1) then
      this%F_DEM_now = 2
      this%F_DEM_next = 1
    else
      this%F_DEM_now = 1
      this%F_DEM_next = 2
    end if
  end subroutine DEM_swapFnowFnext_DPS
  
  !> Reset the force at the next time step to 0
  subroutine DEM_resetFnext_MEM( particleGroup )
    !> particleGroup to reset F_DEM(F_DEM_next,:) for
    type(mus_particle_group_type), intent(inout) :: particleGroup
    ! ------------------------------------------------------ !
    integer :: iParticle
    integer :: F_DEM_next
    ! ------------------------------------------------------ !
    ! reset DEM forces
    do iParticle = 1, particleGroup%particles_MEM%nvals
      F_DEM_next = particleGroup%particles_MEM%val(iParticle)%F_DEM_next
      particleGroup%particles_MEM%val(iParticle)%F_DEM(F_DEM_next,:) = 0.0_rk
    end do
  
  end subroutine DEM_resetFnext_MEM
  
  !> Reset the force at the next time step to 0
  subroutine DEM_resetFnext_DPS( particleGroup )
    !> particleGroup to reset F_DEM(F_DEM_next,:) for
    type(mus_particle_group_type), intent(inout) :: particleGroup
    ! ------------------------------------------------------ !
    integer :: iParticle
    integer :: F_DEM_next
    ! ------------------------------------------------------ !
    ! reset DEM forces
    do iParticle = 1, particleGroup%particles_DPS%nvals
      F_DEM_next = particleGroup%particles_DPS%val(iParticle)%F_DEM_next
      particleGroup%particles_DPS%val(iParticle)%F_DEM(F_DEM_next,:) = 0.0_rk
    end do
  
  end subroutine DEM_resetFnext_DPS

end module mus_particle_DEM_module
