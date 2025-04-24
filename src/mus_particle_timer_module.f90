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
!> mus_particle_timer_module contains routines and data types to time various 
!! parts of the particle code to analyze its performance

module mus_particle_timer_module
  use env_module,         only: rk
  use tem_timer_module,   only: tem_resetTimer, tem_addTimer,&
   &                            tem_getTimerVal, tem_getNTimers
  implicit none

  type mus_particles_timerhandle_type
    !> Main timer that measures the execution time of the entire particle code 
    integer :: mainParticleTimer
    !> Timer for loading of particles from the lua script
    integer :: loadParticleTimer
    !> Timer for full DEM subcycles
    integer :: subcycleTimer
    !> Timer for update and exchange of particle positions
    integer :: positionTimer
    !> Timer for momInc, force computation and communication
    integer :: forceTimer
    !> Timer for velocity update and exchange
    integer :: velocityTimer
    !> Timer measuring the idle time in the communication routines
    !! i.e. the time spent in MPI_Wait() across all particle communication
    integer :: idleTimer
    !> Time spent on exchangePositions routine
    integer :: exchangePositionsTimer
    !> Time spent on exchangeVelocities routine
    integer :: exchangeVelocitiesTimer
    !> Time spent on exchangeMomInc routine
    integer :: exchangeMomIncTimer
    !> Time spent on exchanging DEM forces
    integer :: exchangeDEMForcesTimer
    !> Time spent on exchanging hydrodynamic forces
    integer :: exchangeHydroForcesTimer
    !> Time spent interpolating fluid properties
    integer :: interpolateFluidTimer
    !> Time spent incrementing the auxfield
    integer :: incrementAuxFieldTimer
    !> Time spent in collision detection and adding particles to buffer
    integer :: collisionHandlingTimer
    !> First main handle position in treelm timer object
    integer :: first = 0
    !> Last main handle position in treelm timer object
    integer :: last = -1
  end type mus_particles_timerhandle_type

  type(mus_particles_timerhandle_type), save :: mus_particle_timerHandles


contains


  subroutine mus_init_particleTimer()
    ! add timer handles to measure wall clock time of different routines
    call tem_addTimer( timerHandle = mus_particle_TimerHandles%mainParticleTimer, &
      &                timerName   = 'mainParticleTimer'         )

    ! Set position of the first handle
    mus_particle_TimerHandles%first = tem_getNTimers()

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%loadParticleTimer, &
      &                timerName   = 'loadParticleTimer'         )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%subcycleTimer, &
      &                timerName   = 'subcycleTimer'         )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%positionTimer, &
      &                timerName   = 'positionTimer'         )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%velocityTimer, &
      &                timerName   = 'velocityTimer'         )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%forceTimer, &
      &                timerName   = 'forceTimer'         )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%idleTimer, &
      &                timerName   = 'idleTime'        )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%exchangePositionsTimer, &
      &                timerName   = 'exchangePositions' )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%exchangeVelocitiesTimer, &
      &                timerName   = 'exchangeVelocities' )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%exchangeMomIncTimer, &
      &                timerName   = 'exchangeMomInc' )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%exchangeDEMForcesTimer, &
      &                timerName   = 'exchangeDEMForces' )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%interpolateFluidTimer, &
      &                timerName   = 'interpolateFluid' )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%incrementAuxFieldTimer, &
      &                timerName   = 'incrementAuxField' )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%exchangeHydroForcesTimer, &
      &                timerName   = 'exchangeHydroForces' )

    call tem_addTimer( timerHandle = mus_particle_TimerHandles%collisionHandlingTimer, &
      &                timerName   = 'collisionHandling' )

    ! Set position of the last handle
    mus_particle_TimerHandles%last = tem_getNTimers()

  end subroutine mus_init_particleTimer
end module mus_particle_timer_module
