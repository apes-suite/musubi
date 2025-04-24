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
!> mus_particle_checks_module contains a number of debugging routines to 
!! calculate the fluid and particle momentum so that the conservation of 
!! total momentum can be verified. 

module mus_particle_checks_module
  ! use statements go here
  use mpi
  use env_module,               only: rk, rk_mpi, long_k
  use tem_property_module,      only: prp_solid, prp_particle, &
    &                                 prp_hasBnd, prp_sendHalo
  use mus_geom_module,          only: mus_geom_type
  use mus_param_module,         only: mus_param_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_particle_type_module, only: mus_particle_group_type
  implicit none

  public :: iMomNow
  public :: iMomLast

  integer, save :: iMomNow = 2
  integer, save :: iMomLast = 1


contains


  !> compute_fluid_momentum calculates the total momentum of the fluid for all 
  !! fluid elements in the global domain. The result "totalMomentum" is returned 
  !! in physical units.
  subroutine compute_fluid_momentum(scheme, lev, params, totalMomentum, & 
    &                               comm, myRank, nProcs                )
    !> scheme for access to auxfield
    type(mus_scheme_type), intent(in) :: scheme
    !> level in the octree to compute total momentum on
    integer, intent(in) :: lev
    !> params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> output: total momentum on ENTIRE fluid domain (across all procs)
    real(kind=rk), intent(out) :: totalMomentum(3)
    !> MPI communicator
    integer :: comm
    !> My MPI rank
    integer :: myRank
    !> Number of MPI processes in communicator
    integer :: nProcs
    ! ---------------------------------------- !
    integer :: iElem, elemOff 
    integer :: dens_pos, vel_pos(3)
    integer :: iError
    integer(kind=long_k) :: elemProp
    real(kind=rk) :: rho, u(3), cellMomentum(3), dx
    real(kind=rk) :: processMomentum(3)
    real(kind=rk) :: convertFac_mom
    ! ---------------------------------------- !
    dx = params%physics%dxLvl(lev)
    dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density) &
      &                  %auxField_varPos(1)
    vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity) &
      &                  %auxField_varPos(1:3)
    convertFac_mom = params%physics%rho0 * dx**3 * params%physics%fac(lev)%vel

    totalMomentum = 0.0_rk
    processMomentum = 0.0_rk

    ! Loop over all local fluid cells 
    do iElem = 1, scheme%pdf( lev )%nElems_fluid
      elemOff = (iElem-1)*scheme%varSys%nAuxScalars
      elemProp = scheme%levelDesc(lev)%property(iElem)
      ! Calculate fluid momentum of this cell and add to total
      if( btest(elemProp, prp_solid) ) then
        cycle
      end if
      rho  = scheme%auxField(lev)%val(elemOff + dens_pos)
      u(1) = scheme%auxField(lev)%val(elemOff + vel_pos(1))
      u(2) = scheme%auxField(lev)%val(elemOff + vel_pos(2))
      u(3) = scheme%auxField(lev)%val(elemOff + vel_pos(3))
      cellMomentum = rho*u
      processMomentum = processMomentum + cellMomentum
      
    end do

    ! If we have multiple processes, sum the momentum on each of them
    if(nProcs > 1) then
      call MPI_REDUCE( sendbuf   = processMomentum, &
                     & recvbuf   = totalMomentum,   &
                     & count     = 3,               &
                     & datatype  = rk_mpi,          &
                     & op        = MPI_SUM,         &
                     & root      = 0,               &
                     & comm      = comm,            &
                     & iError    = iError           )
    else
      totalMomentum = processMomentum
    end if

    ! Finally convert the total momentum from lattice units to physical units
    totalMomentum = totalMomentum * convertFac_mom
      
  end subroutine compute_fluid_momentum

  subroutine compute_particle_momentum(particleGroup, lev, params, &
                                      & totalMomentum, comm, myRank, nProcs)
    type(mus_particle_group_type), intent(in) :: particleGroup
    !> Level in the octree that the particles are represented on
    integer, intent(in) :: lev
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> output: total momentum on ENTIRE fluid domain (across all procs)
    real(kind=rk), intent(out) :: totalMomentum(3)
    !> MPI communicator
    integer, intent(in) :: comm
    !> My MPI rank
    integer, intent(in) :: myRank
    !> Number of MPI processes in communicator
    integer, intent(in) :: nProcs
    ! ---------------------------------------- !
    integer :: iParticle, iError
    real(kind=rk) :: mp, up(3), processMomentum(3)
    ! ---------------------------------------- !
    totalMomentum = 0.0_rk
    processMomentum = 0.0_rk

    ! Loop over particles and sum the momentum of particles owned by this process
    if( params%particle_kind == 'DPS'             &
      & .OR. params%particle_kind == 'DPS_twoway' ) then
      do iParticle = 1, particleGroup%particles_DPS%nvals
        if( particleGroup%particles_DPS%val(iParticle)%owner == myRank) then
          mp =  particleGroup%particles_DPS%val(iParticle)%mass
          up =  particleGroup%particles_DPS%val(iParticle)%vel(1:3)
          processMomentum = processMomentum + mp*up
        end if
      end do
    else ! If not DPS particles we must be using MEM particles
      do iParticle = 1, particleGroup%particles_MEM%nvals
        if( particleGroup%particles_MEM%val(iParticle)%owner == myRank) then
          mp =  particleGroup%particles_MEM%val(iParticle)%mass
          up =  particleGroup%particles_MEM%val(iParticle)%vel(1:3)
          processMomentum = processMomentum + mp*up
        end if
      end do
    end if

    ! If we have multiple processes, sum the momentum on each of them
    if(nProcs > 1) then
      call MPI_REDUCE( sendbuf   = processMomentum, &
                     & recvbuf   = totalMomentum,   &
                     & count     = 3,               &
                     & datatype  = rk_mpi,          &
                     & op        = MPI_SUM,         &
                     & root      = 0,               &
                     & comm      = comm,            &
                     & iError    = iError           )
    else
      totalMomentum = processMomentum
    end if
  end subroutine compute_particle_momentum

  subroutine swapMomNowMomLast(iMomNow, iMomLast)
    integer, intent(inout) :: iMomNow
    integer, intent(inout) :: iMomLast
    ! --------------------------------- !
    iMomNow = mod(iMomNow,2) + 1
    iMomLast = mod(iMomLast,2) + 1
  end subroutine swapMomNowMomLast

end module mus_particle_checks_module
