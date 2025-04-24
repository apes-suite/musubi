! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011-2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013-2015, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
!> This module contains the routines related to the LBM-DEM coupling for
!! particle-laden flows using unresolved Discrete Particle Simulations (DPS)

module mus_particle_DPS_module

  use env_module,          only: rk, long_k
  use tem_param_module,    only: div1_3, div1_6, PI, cs2inv, &
    &                            cs4inv, rho0Inv
  use tem_aux_module,      only: tem_abort
  use tem_logging_module,  only: logUnit
  use tem_geometry_module, only: tem_CoordOfReal, tem_posOfId
  use tem_topology_module, only: tem_IdOfCoord, tem_FirstIdAtLevel
  use tem_varSys_module,   only: tem_varSys_type

  use mus_geom_module,        only: mus_geom_type
  use mus_scheme_type_module, only: mus_scheme_type
  use mus_param_module,       only: mus_param_type
  use mus_auxField_module,    only: mus_auxFieldVar_type
  use mus_derVarPos_module,   only: mus_derVarPos_type

  use mus_particle_type_module, only: mus_particle_DPS_type,   &
    &                                 mus_particle_group_type, &
    &                                 allocateProcessMasks
  use mus_particle_comm_type_module, only: mus_particles_communication_type
  use mus_particle_comm_module, only: exchangeNewParticles_DPS, &
    &                                 exchangeParticlesToRemove_DPS
  use mus_particle_type_module, only: remove_particle_from_da_particle_DPS, &
    &                                 interpolateFluidPropFunc, &
    &                                 calcVelAndPGradFunc
  use mus_particle_aux_module, only: getBaryOfCoord,  &
    &                                cross_product,   &
    &                                followNeighPath, &
    &                                getPosOfCoord,   &
    &                                findPartitionOfTreeID
  use mus_particle_logging_module, only: getParticleLogUnit, &
    &                                    logParticleData,    &
    &                                    closeParticleLog
  use mus_particle_logging_type_module, only: pgDebugLog
  use mus_particle_boundary_module, only: pgBndData,        &
    &                                     getNeighborCoord, &
    &                                     wrapPeriodicCoord
  use mus_particle_interpolation_module, only: intp_1D_delta,                  &
    &                                          mus_particle_interpolator_type, &
    &                                          intp_1D_peskin,                 &
    &                                          intp1D_linear,                  &
    &                                          getwght1d_linear,               &
    &                                          get_xd,                         &
    &                                          getInterpolationBnds

  implicit none

  !> Interface for routine to get the value of the velocity and
  !! pressure field at a specified integer coordinate
  interface grabValueAtCoord
    module procedure grabValueAtCoord_twowaycoupled
    module procedure grabValueAtCoord_onewaycoupled
  end interface


contains


  !> initParticle_DPS performs the initialization of a DPS particle once
  !! it has been created. This consists of determining the coordinates of
  !! particle origin and determining which processes this particle exists
  !! on.
  subroutine initParticle_DPS( particle, interpolator, particleID, geometry, &
                             & scheme, myRank, comm                          )
    !> Particle to initialize
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Interpolation data type. We don't interpolate stuff here but we use it
    !! to determine the directions we need to search for local elements in
    !! (only x, y for d2q9, x, y, z for d3q19)
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Number to uniquely identify particle particle
    integer, intent(in) :: particleID
    !> Geometry to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> MPI rank of this process
    integer, intent(in) :: myRank
    !> Communication data type for particles
    type(mus_particles_communication_type) :: comm
    ! -----------------------------------------------!
    real(kind=rk) :: dx
    integer :: lev                  ! level to map particle on
    ! -----------------------------------------------!

    lev = geometry%tree%global%maxLevel
    dx = geometry%tree%global%BoundingCubeLength / 2**lev

    ! particle radius in lattice units, rounded up
    particle%Rn = ceiling(particle%radius / dx)

    ! Set initial coordOfOrigin
    particle%coordOfOrigin = tem_CoordOfReal( mesh  = geometry%tree,     &
                                            & point = particle%pos(1:3), &
                                            & level = lev                )
    particle%posOfOrigin = getPosOfCoord( coord = particle%coordOfOrigin, &
                                        & scheme = scheme )

    ! Call updateCoordOfOrigin so that also oldCoordOfOrigin is set
    call updateCoordOfOrigin_DPS( this     = particle, &
                                & scheme   = scheme,   &
                                & geometry = geometry, &
                                & myRank   = myRank    )


    call updateParticleOwner( this     = particle,     &
                            & scheme   = scheme,       &
                            & geometry = geometry,     &
                            & myrank   = myRank,       &
                            & procs    = comm%proc,    &
                            & nProcs   = comm%nProcs   )

    call updateExistsOnProc_DPS( this     = particle,     &
                               & interpolator = interpolator, &
                               & scheme = scheme,     &
                               & geometry = geometry,     &
                               & procs    = comm%proc,    &
                               & nProcs   = comm%nProcs,  &
                               & myRank   = myRank        )

  end subroutine initParticle_DPS

  !> updateCoordOfOrigin updates the integer coordinate of
  !! the origin of a particle
  subroutine updateCoordOfOrigin_DPS( this, scheme, geometry, myRank )
    !> Particle to update coordOfOrigin of
    type(mus_particle_DPS_type), intent(inout) :: this
    !> Scheme to look up posOfOrigin
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry information to determine new coord of origin
    type(mus_geom_type), intent(in) :: geometry
    !> This process's rank
    integer, intent(in) :: myRank
    ! --------------------------------!
    integer :: moveDir(3)
    integer :: lev
    integer :: newCoordOfOrigin(4)
    ! --------------------------------!
    lev = this%coordOfOrigin(4)
    ! Check if particle has moved more than one lattice site
    newCoordOfOrigin = tem_CoordOfReal( mesh  = geometry%tree, &
                                      & point = this%pos(1:3), &
                                      & level = lev            )

    moveDir = newCoordOfOrigin(1:3) - this%coordOfOrigin(1:3)


    if( .NOT. all(moveDir == 0) ) then
      ! If particle has moved to a different coordinate, update
      ! the posOfOrigin
      this%posOfOrigin =  getPosOfCoord( coord  = newCoordOfOrigin, &
                                       & scheme = scheme            )
    end if
    ! Update coordinate of origin
    this%oldCoordOfOrigin(:) = this%coordOfOrigin(:)
    this%coordOfOrigin(:) = newCoordOfOrigin(:)

  end subroutine updateCoordOfOrigin_DPS

  !> UpdateParticleOwner updates the "owner process" of each particle,
  !! which is the process responsible for performing operations
  !! (e.g. modifying velocity and position) on the particle and sending
  !! updates to other processes.
  subroutine updateParticleOwner( this, scheme, geometry, myRank, &
                                & procs, nProcs                   )
    !> Particle to update owner of
    type(mus_particle_DPS_type), intent(inout) :: this
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> This process's rank
    integer, intent(in) :: myRank
    !> Array of neighbor procs on which particle could possibly exist
    integer, intent(in) :: procs(:)
    !> Number of procs in array procs
    integer, intent(in) :: nProcs
    ! ------------------------------------- !
    integer(kind=long_k) :: TreeID
    integer :: ldPos
    integer :: iElemProc
    integer :: lev
    ! ------------------------------------- !
    lev = this%coordOfOrigin(4)

    TreeID = tem_IdOfCoord(coord = this%coordOfOrigin)
    ! First check whether this process is the owner
    ! In that case the TreeID should be within this process's
    ! part_first and part_last
    if( TreeID >= geometry%tree%Part_First( myRank + 1) &
      & .AND. TreeID <= geometry%tree%Part_Last( myRank + 1) ) then
      ! I am the owner of this particle. Now check whether the particle
      ! is within the actual fluid domain by checking if we can find the
      ! TreeID in this process's total list
      ldPos = tem_PosOfId( sTreeID    = TreeID,                      &
                         & treeIDlist = scheme%levelDesc(lev)%total, &
                         & lower      = 1,                           &
                         & upper      = scheme%pdf(lev)%nElems_fluid )
      if(ldPos > 0) then
        ! Found the TreeID in total list. Set particle owner to myRank and
        ! exit subroutine
        this%owner = myRank
        this%removeParticle_global = .FALSE.
        return
      else
        this%owner = myRank
        open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
        write(pgDebugLog%lu,*) "updateParticleOwner: Particle ID ", &
          & this%particleID, " moved outside domain"
          write(pgDebugLog%lu,*) "particle pos = ", this%pos(1:3)
          this%removeParticle_global = .TRUE.
        return
        close(pgDebugLog%lu)
      end if ! ldPos > 0
    else
      ! TreeID does not belong to our process, so find out to which other
      ! process it belongs

      ! First try to look only in my neighbor procs. Most likely we'll find
      ! the particle owner there. If owner is NOT found on neighbor procs,
      ! this%owner will be set to -1 to indicate this
      call findPartitionOfTreeID( sTreeID   = TreeID,     &
                                & geometry  = geometry,   &
                                & procs     = procs,      &
                                & nProcs    = nProcs,     &
                                & outProc   = this%owner, &
                                & procIndex = ielemProc   )

      ! If we still cannot find the owner in neighbor procs, look in all processes
      if(this%owner < 0) then
        call findPartitionOfTreeID( sTreeID   = TreeID,     &
                                  & geometry  = geometry,   &
                                  & nProcs    = geometry%tree%global%nParts,    &
                                  & outProc   = this%owner, &
                                  & procIndex = ielemProc   )

        ! If we cannot find the TreeID on ANY process, remove particle from the
        ! global domain
        if(this%owner < 0 ) then
          this%owner = myRank
          open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
          write(pgDebugLog%lu,*) "updateParticleOwner: Could not find particle ID ", &
            & "on any process, removing particle ", &
            & this%particleID, " from global domain"
          write(pgDebugLog%lu,*) "particle pos = ", this%pos(1:3)
          close(pgDebugLog%lu)
          this%removeParticle_global = .TRUE.
        end if
      end if
    end if ! TreeID between myRank's part_first and part_last

  end subroutine updateParticleOwner

  !> updateExistsOnProc updates the boolean values of the array of neighbor procs
  !! that tell us whether the particle exists on that process or not.
  subroutine updateExistsOnProc_DPS( this, interpolator, scheme, geometry, &
                                   & procs, nProcs, myRank                 )
    !> Particle to update owner of
    type(mus_particle_DPS_type), intent(inout) :: this
    !> Interpolation data type. We don't interpolate stuff here but we use it to
    !! determine the directions we need to search for local elements in
    !! (only x, y for d2q9, x, y, z for d3q19)
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
      !> Array of neighbor procs on which particle could possibly exist
    integer, intent(in) :: procs(:)
    !> Number of procs in array procs
    integer, intent(in) :: nProcs
    !> This process's rank
    integer, intent(in) :: myRank
    ! ------------------------------------- !
    integer :: iDir
    integer :: coord(4)
    integer :: lev, upperBound, iproc, iElemProc, elemProc
    integer(kind=long_k) :: TreeID
    integer :: elemPos, nghPos
    logical :: oldExistsOnProc( 1:nProcs )
    ! ------------------------------------- !
    lev = geometry%tree%global%maxLevel
    upperBound = 2**lev

    this%addToProc = .FALSE.
    this%removeFromProc = .FALSE.

    ! Store the existsOnProc mask from the previous time step
    do iproc = 1, nProcs
      oldExistsOnProc(iproc) = this%existsOnProc(iproc)
    end do

    this%existsOnProc = .FALSE.

    ! Set removeParticle_local to true initially. Will be set to false once we encounter a
    ! particle that is local to our proc within a bounding box of 1 lattice
    ! site around this particle's coordOfOrigin
    this%removeParticle_local = .TRUE.

    ! Get position of this element in the total list
    elemPos = this%posOfOrigin

    do iDir = 1, interpolator%Nelems
      ! Get coordinate of this neighbor
      coord(1:3) = this%coordOfOrigin(1:3) + interpolator%neighDirs(1:3,iDir)
      coord(4) = this%coordOfOrigin(4)

      call wrapPeriodicCoord( coord        = coord,    &
                            & boundaryData = pgBndData )

      ! Check if coordinate is within actual simulation domain
      if( any(coord(1:3) < 0) .OR. any(coord(1:3) > upperBound) ) then
        cycle
      else
        ! Update removeParticle_local if it is not yet set to FALSE for this iteration
        nghPos = followNeighPath(                                        &
                  & neighpath    = interpolator%neighPaths(1:3,idir),    &
                  & startPos     = elemPos,                              &
                  & neigh        = scheme%levelDesc(lev)%neigh(1),       &
                  & nElems_fluid = scheme%pdf(lev)%nElems_fluid,         &
                  & zeroDir      = scheme%layout%stencil(1)%restPosition )
        if(nghPos < 1) then
          ! If we cannot reach nghPos by followNeighPath (for example because )
          ! elemPos points to a halo element, use the (expensive) Treelm routines
          TreeID = tem_IdOfCoord(coord = coord)
          nghPos = tem_PosOfId( sTreeID    = TreeID,                      &
                              & treeIDlist = scheme%levelDesc(lev)%total, &
                              & lower      = 1,                           &
                              & upper      = scheme%pdf(lev)%nElems_fluid )
        else
          TreeID = scheme%levelDesc(lev)%total(nghPos)
        end if

        if(nghPos > 0 .AND. nghPos <= scheme%pdf(lev)%nElems_fluid) then
          this%removeParticle_local = .FALSE.
        end if

        ! Check which process this element belongs to to update
        ! existsOnProc masks
        if( TreeID < geometry%tree%Part_First(myRank + 1) &
          & .OR. TreeID > geometry%tree%Part_Last(myRank + 1) ) then
          ! Element is not local so find the neighbor process that
          ! this element is on
          call findPartitionOfTreeID(                                         &
                          & sTreeID = TreeID,                                 &
                          & geometry = geometry,                              &
                          & procs = procs,    &
                          & nProcs = nProcs, &
                          & outProc = elemProc,                               &
                          & procIndex = iElemProc                             )
          if(elemProc >= 0 ) then
            this%existsOnProc( iElemProc ) = .TRUE.
          end if ! elemProc >= 0
        end if ! element is NOT local

      end if ! coordinate within simulation domain

    end do ! iDir

    ! Set addToProc and removeFromProc based on old and new values
    ! of existsOnProc
    do iproc = 1, nProcs
      if( (.NOT. oldExistsOnProc(iproc)) .AND. this%existsOnProc(iproc)  ) then
        this%addToProc(iproc) = .TRUE.
      else if( oldExistsOnProc(iproc) .AND. ( .NOT. this%existsOnProc(iproc) ) ) then
        this%removeFromProc(iproc) = .TRUE.
      end if
    end do ! iproc

  end subroutine updateExistsOnProc_DPS

  !> Routine that applies forces from particles to the fluid for unresolved
  !! DPS particles based on the VANS equations
  subroutine transferMomentumToFluid_DPS( particle, interpolator, scheme, &
                                               & geometry, params, Ftot                )
    !> Particle to interpolate fluid properties to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> interpolator object containing stencil and weight function info
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Output, total force exerted on fluid by particle
    real(kind=rk), intent(out) :: Ftot(3)
    ! ----------------------------------------------- !
    integer :: lev, nghCoord(4), iNghDir, nghPos
    real(kind=rk) :: dx, bary(3), rbary_lat(3), r_lat(3)
    real(kind=rk) :: wght_x, wght_y, wght_z, G_lat(3)
    ! ----------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    dx = params%physics%dxLvl(lev)

    bary = getBaryOfCoord( coord  = particle%coordOfOrigin,      &
                         & origin = geometry%tree%global%Origin, &
                         & dx     = dx                           )

    rbary_lat = (particle%pos(1:3) - bary)/dx

    ! Loop over neighboring cells
    do iNghDir = 1, interpolator%Nelems
      ! Get coordinate of this neighbor element
      nghCoord(1:3) = particle%coordOfOrigin(1:3) + interpolator%neighDirs(1:3,iNghDir)
      nghCoord(4) = particle%coordOfOrigin(4)

      ! Get position of neighbor element in total list
      nghPos = followNeighPath( neighPath    = interpolator%neighPaths(1:3, iNghDir), &
                              & startPos     = particle%posOfOrigin,                  &
                              & neigh        = scheme%levelDesc(lev)%neigh(1), &
                              & nElems_fluid = scheme%pdf(lev)%nElems_fluid, &
                              & zeroDir      = scheme%layout%stencil(1)%restPosition )
      if(nghPos < 1) then
        ! If we cannot follow neighPath to get this element's position, get its position
        ! using (expensive) call to getPosOfCoord which uses the tem_topology routines
        nghPos = getPosOfCoord( coord = nghCoord, scheme = scheme )
      end if

      if( nghPos > 0 .AND. nghPos < scheme%pdf(lev)%nElems_fluid ) then
        ! If this is a local fluid element, modify the PDF values

        ! Compute vector from particle to baryCenter of sample point
        r_lat = interpolator%neighDirs(1:3,iNghDir) - rbary_lat
        r_lat = abs(r_lat)
        ! Get weights in each Cartesian direction for this coord
        wght_x = interpolator%getWght_x( r_lat(1) )
        wght_y = interpolator%getWght_y( r_lat(2) )
        wght_z = interpolator%getWght_z( r_lat(3) )
        G_lat = -particle%Favg(1:3)*wght_x*wght_y*wght_z
        G_lat = G_lat / dx**3
        G_lat = G_lat / params%physics%fac(lev)%body_force

        call applySrc_toElem_eps( G = G_lat,           &
                                & posInTotal = nghPos, &
                                & scheme = scheme,     &
                                & iLevel = lev         )
      end if
    end do

    particle%Favg = 0.0_rk
  end subroutine transferMomentumToFluid_DPS

  !> Routine that applies forces from particles to the fluid for unresolved
  !! DPS two-way coupled particles
  subroutine transferMomentumToFluid_DPS_twoway( particle, interpolator, scheme, &
                                               & geometry, params, Ftot                )
    !> Particle to interpolate fluid properties to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> interpolator object containing stencil and weight function info
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Output, total force exerted on fluid by particle
    real(kind=rk), intent(out) :: Ftot(3)
    ! ----------------------------------------------- !
    integer :: lev, nghCoord(4), iNghDir, nghPos
    real(kind=rk) :: dx, bary(3), rbary_lat(3), r_lat(3)
    real(kind=rk) :: wght_x, wght_y, wght_z, G_lat(3)
    ! ----------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    dx = params%physics%dxLvl(lev)

    bary = getBaryOfCoord( coord  = particle%coordOfOrigin,      &
                         & origin = geometry%tree%global%Origin, &
                         & dx     = dx                           )

    rbary_lat = (particle%pos(1:3) - bary)/dx

    ! Loop over neighboring cells
    do iNghDir = 1, interpolator%Nelems
      ! Get coordinate of this neighbor element
      nghCoord(1:3) = particle%coordOfOrigin(1:3) + interpolator%neighDirs(1:3,iNghDir)
      nghCoord(4) = particle%coordOfOrigin(4)

      ! Get position of neighbor element in total list
      nghPos = followNeighPath( neighPath    = interpolator%neighPaths(1:3, iNghDir), &
                              & startPos     = particle%posOfOrigin,                  &
                              & neigh        = scheme%levelDesc(lev)%neigh(1), &
                              & nElems_fluid = scheme%pdf(lev)%nElems_fluid, &
                              & zeroDir      = scheme%layout%stencil(1)%restPosition )
      if(nghPos < 1) then
        ! If we cannot follow neighPath to get this element's position, get its position
        ! using (expensive) call to getPosOfCoord which uses the tem_topology routines
        nghPos = getPosOfCoord( coord = nghCoord, scheme = scheme )
      end if

      if( nghPos > 0 .AND. nghPos < scheme%pdf(lev)%nElems_fluid ) then
        ! If this is a local fluid element, modify the PDF values

        ! Compute vector from particle to baryCenter of sample point
        r_lat = interpolator%neighDirs(1:3,iNghDir) - rbary_lat
        r_lat = abs(r_lat)
        ! Get weights in each Cartesian direction for this coord
        wght_x = interpolator%getWght_x( r_lat(1) )
        wght_y = interpolator%getWght_y( r_lat(2) )
        wght_z = interpolator%getWght_z( r_lat(3) )
        G_lat = -particle%Favg(1:3)*wght_x*wght_y*wght_z
        G_lat = G_lat / dx**3
        G_lat = G_lat / params%physics%fac(lev)%body_force

        call applySrc_toElem( G = G_lat,           &
                            & posInTotal = nghPos, &
                            & scheme = scheme,     &
                            & iLevel = lev         )
      end if
    end do
    particle%Favg = 0.0_rk
  end subroutine transferMomentumToFluid_DPS_twoway

  !> Routine to modify the auxField (velocity) with the forces exerted
  !! by particles on the fluid
  subroutine addParticleSourceToAuxfield_DPS( particle, interpolator, scheme, &
                                               & geometry, params, Ftot                )
    !> Particle to interpolate fluid properties to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> interpolator object containing stencil and weight function info
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Output, total force exerted on fluid by particle
    real(kind=rk), intent(out) :: Ftot(3)
    ! ----------------------------------------------- !
    integer :: lev, nghCoord(4), iNghDir, nghPos
    real(kind=rk) :: dx, bary(3), rbary_lat(3), r_lat(3)
    real(kind=rk) :: wght_x, wght_y, wght_z, G_lat(3)
    ! ----------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    dx = params%physics%dxLvl(lev)

    bary = getBaryOfCoord( coord  = particle%coordOfOrigin,      &
                         & origin = geometry%tree%global%Origin, &
                         & dx     = dx                           )

    rbary_lat = (particle%pos(1:3) - bary)/dx

    ! Loop over neighboring cells
    do iNghDir = 1, interpolator%Nelems
      ! Get coordinate of this neighbor element
      nghCoord(1:3) = particle%coordOfOrigin(1:3) + interpolator%neighDirs(1:3,iNghDir)
      nghCoord(4) = particle%coordOfOrigin(4)

      ! Get position of neighbor element in total list
      nghPos = followNeighPath( neighPath    = interpolator%neighPaths(1:3, iNghDir), &
                              & startPos     = particle%posOfOrigin,                  &
                              & neigh        = scheme%levelDesc(lev)%neigh(1), &
                              & nElems_fluid = scheme%pdf(lev)%nElems_fluid, &
                              & zeroDir      = scheme%layout%stencil(1)%restPosition )
      if(nghPos < 1) then
        ! If we cannot follow neighPath to get this element's position, get its position
        ! using (expensive) call to getPosOfCoord which uses the tem_topology routines
        nghPos = getPosOfCoord( coord = nghCoord, scheme = scheme )
      end if

      if( nghPos > 0 .AND. nghPos < scheme%pdf(lev)%nElems_fluid ) then
        ! If this is a local fluid element, modify the PDF values

        ! Compute vector from particle to baryCenter of sample point
        r_lat = interpolator%neighDirs(1:3,iNghDir) - rbary_lat
        r_lat = abs(r_lat)
        ! Get weights in each Cartesian direction for this coord
        wght_x = interpolator%getWght_x( r_lat(1) )
        wght_y = interpolator%getWght_y( r_lat(2) )
        wght_z = interpolator%getWght_z( r_lat(3) )

        G_lat = -particle%F(1:3)*wght_x*wght_y*wght_z
        G_lat = G_lat / dx**3
        G_lat = G_lat / params%physics%fac(lev)%body_force

        ! Call modifyAuxFieldCell()
        call modify_AuxField_of_Elem( posInTotal = nghPos,                   &
                                    & scheme     = scheme,                   &
                                    & auxField   = scheme%auxField(lev)%val, &
                                    & iLevel     = lev,                      &
                                    & varSys     = scheme%varSys,            &
                                    & derVarPos  = scheme%derVarPos,         &
                                    & G_lat      = G_lat                     )
      end if
    end do
  end subroutine addParticleSourceToAuxfield_DPS

  !> Remove the momentum increments added during subcycles to prevent double-adding
  !! the momentum of particles to the fluid.
  subroutine recalculate_auxField_DPS( particle, interpolator, scheme, &
                                  & geometry, params                   )
    !> Particle to interpolate fluid properties to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> interpolator object containing stencil and weight function info
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    ! ----------------------------------------------- !
    integer :: lev, nghCoord(4), iNghDir, nghPos
    real(kind=rk) :: dx, bary(3), rbary_lat(3), r_lat(3)
    real(kind=rk) :: wght_x, wght_y, wght_z, momentumIncrement(3)
    ! ----------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    dx = params%physics%dxLvl(lev)


    bary = getBaryOfCoord( coord  = particle%coordOfOrigin,      &
                         & origin = geometry%tree%global%Origin, &
                         & dx     = dx                           )

    rbary_lat = (particle%pos(1:3) - bary)/dx

    ! Loop over neighboring cells
    do iNghDir = 1, interpolator%Nelems
      ! Get coordinate of this neighbor element
      nghCoord(1:3) = particle%coordOfOrigin(1:3) + interpolator%neighDirs(1:3,iNghDir)
      nghCoord(4) = particle%coordOfOrigin(4)

      ! Get position of neighbor element in total list
      nghPos = followNeighPath( neighPath    = interpolator%neighPaths(1:3, iNghDir), &
                              & startPos     = particle%posOfOrigin,                  &
                              & neigh        = scheme%levelDesc(lev)%neigh(1), &
                              & nElems_fluid = scheme%pdf(lev)%nElems_fluid, &
                              & zeroDir      = scheme%layout%stencil(1)%restPosition )
      if(nghPos < 1) then
        ! If we cannot follow neighPath to get this element's position, get its position
        ! using (expensive) call to getPosOfCoord which uses the tem_topology routines
        nghPos = getPosOfCoord( coord = nghCoord, scheme = scheme )
      end if

      if( nghPos > 0 .AND. nghPos < scheme%pdf(lev)%nElems_fluid ) then
        ! If this is a local fluid element, modify the PDF values

        ! Compute vector from particle to baryCenter of sample point
        r_lat = interpolator%neighDirs(1:3,iNghDir) - rbary_lat
        r_lat = abs(r_lat)
        ! Get weights in each Cartesian direction for this coord
        wght_x = interpolator%getWght_x( r_lat(1) )
        wght_y = interpolator%getWght_y( r_lat(2) )
        wght_z = interpolator%getWght_z( r_lat(3) )

        ! Calculate avg particle force during last LBM time step on this cell in lattice units
        ! Force on fluid by particle is -particle%Favg and we SUBTRACT this, so total term is positive
        ! Favg.
        momentumIncrement = 0.5 * particle%Favg(1:3)*wght_x*wght_y*wght_z / params%physics%fac(lev)%force

        call add_Momentum_Increment_to_Auxfield(                                 &
                                & posInTotal         = nghPos,                   &
                                & scheme             = scheme,                   &
                                & auxField           = scheme%auxField(lev)%val, &
                                & iLevel             = lev,                      &
                                & varSys             = scheme%varSys,            &
                                & derVarPos          = scheme%derVarPos,         &
                                & momentumIncrement  = momentumIncrement         )

      end if
    end do
  end subroutine recalculate_auxField_DPS

  !> Add momentum increments due to particles within DEM subcycles to the
  !! auxField
  subroutine incrementAuxField_DPS( particle, interpolator, scheme, &
                                  & geometry, params, dt_DEM_lat                )
    !> Particle to interpolate fluid properties to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> interpolator object containing stencil and weight function info
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Time step of DEM subcycles, in lattice units
    real(kind=rk), intent(in) :: dt_DEM_lat
    ! ----------------------------------------------- !
    integer :: lev, nghCoord(4), iNghDir, nghPos
    real(kind=rk) :: dx, bary(3), rbary_lat(3), r_lat(3)
    real(kind=rk) :: wght_x, wght_y, wght_z, j_lat(3)
    ! ----------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    dx = params%physics%dxLvl(lev)


    bary = getBaryOfCoord( coord  = particle%coordOfOrigin,      &
                         & origin = geometry%tree%global%Origin, &
                         & dx     = dx                           )

    rbary_lat = (particle%pos(1:3) - bary)/dx

    particle%Favg = particle%Favg + particle%F(1:3)*dt_DEM_lat

    ! Loop over neighboring cells
    do iNghDir = 1, interpolator%Nelems
      ! Get coordinate of this neighbor element
      nghCoord(1:3) = particle%coordOfOrigin(1:3) + interpolator%neighDirs(1:3,iNghDir)
      nghCoord(4) = particle%coordOfOrigin(4)

      ! Get position of neighbor element in total list
      nghPos = followNeighPath( neighPath    = interpolator%neighPaths(1:3, iNghDir), &
                              & startPos     = particle%posOfOrigin,                  &
                              & neigh        = scheme%levelDesc(lev)%neigh(1), &
                              & nElems_fluid = scheme%pdf(lev)%nElems_fluid, &
                              & zeroDir      = scheme%layout%stencil(1)%restPosition )
      if(nghPos < 1) then
        ! If we cannot follow neighPath to get this element's position, get its position
        ! using (expensive) call to getPosOfCoord which uses the tem_topology routines
        nghPos = getPosOfCoord( coord = nghCoord, scheme = scheme )
      end if

      if( nghPos > 0 .AND. nghPos < scheme%pdf(lev)%nElems_fluid ) then
        ! If this is a local fluid element, modify the PDF values

        ! Compute vector from particle to baryCenter of sample point
        r_lat = interpolator%neighDirs(1:3,iNghDir) - rbary_lat
        r_lat = abs(r_lat)
        ! Get weights in each Cartesian direction for this coord
        wght_x = interpolator%getWght_x( r_lat(1) )
        wght_y = interpolator%getWght_y( r_lat(2) )
        wght_z = interpolator%getWght_z( r_lat(3) )

        ! Calculate momentum increment on this cell in lattice units
        j_lat = dt_DEM_lat * (-particle%F(1:3)*wght_x*wght_y*wght_z) / params%physics%fac(lev)%force

        call add_Momentum_Increment_to_Auxfield(                                 &
                                & posInTotal         = nghPos,                   &
                                & scheme             = scheme,                   &
                                & auxField           = scheme%auxField(lev)%val, &
                                & iLevel             = lev,                      &
                                & varSys             = scheme%varSys,            &
                                & derVarPos          = scheme%derVarPos,         &
                                & momentumIncrement  = j_lat                     )

      end if
    end do
  end subroutine incrementAuxField_DPS

  !> ModifyAuxFieldCell is a helper routine that modifies the auxfield in the cell with
  !! position elemPos in the total list. It takes the total force F_particle and the vector
  !! r_lat (in lattice units) pointing from the particle location to the barycenter of the
  !! fluid element, then computes the weight and modifies the auxField.
  !! Performing addForceToAuxFieldCell over every neighboring cell around a particle
  !! (including the one the particle is on) distributes the entire force.
  subroutine addForceToAuxFieldCell( auxField, vel_pos, dens_pos, nAuxScalars, &
    & elemPos, F, r_lat, convFac_force, dt_DEM_lat, interpolator  )
    type(mus_auxFieldVar_type), intent(inout) :: auxField
    !> Position of velocity in auxField
    integer, intent(in) :: vel_pos(3)
    !> Position of density in auxField
    integer, intent(in) :: dens_pos
    !> Number of scalars in auxField
    integer, intent(in) :: nAuxScalars
    !> Position of element to add force contribution to in auxField
    !! inside total list
    integer, intent(in) :: elemPos
    !> Total force to apply to auxField (only a fraction of this will
    !! be applied to the element at elemPos)
    real(kind=rk), intent(in) :: F(3)
    !> Vector (in lattice units) pointing from the location where F is
    !! applied to the barycenter of the cell at elemPos
    real(kind=rk), intent(in) :: r_lat(3)
    !> Conversion factor between physical and lattice units
    real(kind=rk) :: convFac_force
    !> DEM subcycle time step, in lattice units
    real(kind=rk) :: dt_DEM_lat
    !> Interpolation object containing stencil and weight info
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    ! --------------------------------------------------------- !
    real(kind=rk) :: wght_x, wght_y, wght_z
    real(kind=rk) :: F_cell(3)
    integer :: elemOff
    real(kind=rk) :: inv_rho
    real(kind=rk) :: velocityIncrement(3)
    ! --------------------------------------------------------- !
    ! Get weights
    wght_x = interpolator%getWght_x( r_lat(1) )
    wght_y = interpolator%getWght_y( r_lat(2) )
    wght_z = interpolator%getWght_z( r_lat(3) )

    ! Compute portion of particle hydrodynamic force to distribute to this cell
    F_cell = wght_x * wght_y * wght_z * F(1:3)

    ! Convert to lattice units
    F_cell = F_cell / convFac_force

    ! Transfer momentum from particle to fluid by modifying the velocity in auxField
    elemOff = (elemPos-1)*nAuxScalars
    inv_rho = 1.0_rk / auxField%val(elemOff + dens_pos)

    velocityIncrement = F_cell*dt_DEM_lat*inv_rho

    !$OMP ATOMIC UPDATE
    ! Add this term to the velocity field
    auxField%val(elemOff+vel_pos(1)) &
      & = auxField%val(elemOff+vel_pos(1)) + velocityIncrement(1)
    !$OMP END ATOMIC

    !$OMP ATOMIC UPDATE
    auxField%val(elemOff+vel_pos(2)) &
      & = auxField%val(elemOff+vel_pos(2)) + velocityIncrement(2)
    !$OMP END ATOMIC

    !$OMP ATOMIC UPDATE
    auxField%val(elemOff+vel_pos(3)) &
      & = auxField%val(elemOff+vel_pos(3)) + velocityIncrement(3)
    !$OMP END ATOMIC
  end subroutine addForceToAuxFieldCell

  !> Distribute the volume of a particle to update the fluid volume fraction field
  subroutine distributeParticleVolume( particle, interpolator, scheme, geometry, params )
    !> Particle to interpolate fluid properties to
    type(mus_particle_DPS_type), intent(in) :: particle
    !> Interpolation object containing stencil and weight info
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    ! ------------------------------------------!
    integer :: iDir, lev
    integer :: coord(4)
    integer :: elemOff, elemPos, nghPos, vol_frac_pos
    real(kind=rk) :: baryOfOrigin(3), rbary_lat(3), r_lat(3)
    real(kind=rk) :: geom_origin(3), dx
    real(kind=rk) :: wght_x, wght_y, wght_z
    real(kind=rk) :: Vparticle

    ! ------------------------------------------!
    ! Find neighboring fluid cells
    lev = particle%coordOfOrigin(4)
    dx = params%physics%dxLvl(lev)
    geom_origin = geometry%tree%global%origin
    vol_frac_pos = scheme%varSys%method%val(scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)

    elemPos = particle%posOfOrigin

    ! Find out what the nearest 8 barycenters are to which we will
    ! distribute the particle volume
    baryOfOrigin = getBaryOfCoord( coord  = particle%coordOfOrigin(1:3), &
                                 & origin = geom_origin,             &
                                 & dx     = dx                       )

    ! rbary_lat = vector FROM particle CoordOfOrigin barycenter to particle
    rbary_lat = (particle%pos(1:3) - baryOfOrigin)/dx

    ! Calculate volume of the particle
    Vparticle = (4.0/3.0)*PI*particle%radius**3

    ! Loop over lattice cells neighboring the particle and distribute
    ! particle volume over them

    do iDir = 1, interpolator%Nelems
      coord(1:3) = particle%coordOfOrigin(1:3) + interpolator%neighDirs(1:3,idir)
      coord(4) = particle%coordOfOrigin(4)
      ! Get the position of this neighbor element in the total list
      ! First try to do this by following the interpolator's neighPaths
      nghPos = followNeighPath(                                        &
                & neighpath    = interpolator%neighPaths(1:3,idir),    &
                & startPos     = elemPos,                              &
                & neigh        = scheme%levelDesc(lev)%neigh(1),       &
                & nElems_fluid = scheme%pdf(lev)%nElems_fluid,         &
                & zeroDir      = scheme%layout%stencil(1)%restPosition )

      if(nghPos < 1) then
        ! If we cannot follow neighPath to get this element's position, get its position
        ! using (expensive) call to getPosOfCoord which uses the tem_topology routines
        nghPos = getPosOfCoord( coord = coord, scheme = scheme )
      end if

      if(nghPos < 1) then
        ! If neighPos is still less than 0, skip it.
        ! This can be the case for example near walls.
        cycle
      end if ! neighPos < 0
      elemOff = (nghPos-1)*scheme%varSys%nAuxScalars

      ! Compute vector from particle to baryCenter of sample point
      r_lat = interpolator%neighDirs(1:3,iDir) - rbary_lat
      r_lat = abs(r_lat)

      ! Get weights in each Cartesian direction for this coord
      wght_x = interpolator%getWght_x( r_lat(1) )
      wght_y = interpolator%getWght_y( r_lat(2) )
      wght_z = interpolator%getWght_z( r_lat(3) )

      !$OMP ATOMIC UPDATE
      scheme%auxField(lev)%val(elemOff + vol_frac_pos) =                 &
                      & scheme%auxField(lev)%val(elemOff + vol_frac_pos) &
                      &       - (wght_x*wght_y*wght_z*Vparticle) / dx**3
      !$OMP END ATOMIC
    end do ! iDir

  end subroutine distributeParticleVolume

  !> Get the index in the total list of a specified integer coordinate
  subroutine getIndexOfCoordInTotal(scheme, coord, ldPos)
    integer, intent(in) :: coord(4)
    type(mus_scheme_type), intent(in) :: scheme
    integer, intent(out) :: ldPos
    ! ------------------------------------------!
    integer(kind=long_k) :: TreeID
    integer :: lev
    ! ------------------------------------------!
    lev = coord(4)
    ldPos = -1

    ! get TreeID and position in complete tree
    TreeID = tem_IdOfCoord(coord = coord)

    ldPos = tem_PosOfId( sTreeID = TreeID,                        &
                      & treeIDlist = scheme%levelDesc(lev)%total, &
                      & lower      = 1,                           &
                      & upper      = scheme%pdf(lev)%nElems_fluid )
    if( ldPos <= 0) then
      ! If we could not find ldPos in local elems, look in halos
      ldPos = tem_PosOfId( sTreeID = TreeID,                               &
                        & treeIDlist = scheme%levelDesc(lev)%total,        &
                        & lower      = scheme%pdf(lev)%nElems_fluid + 1,  &
                        & upper      = scheme%pdf(lev)%nElems_local        )

    end if ! could not find element in local elems

  end subroutine getIndexOfCoordInTotal

  !> Main routine to update the fluid volume fraction in the auxField
  subroutine mus_particles_updateFluidVolumeFraction( particleGroup, scheme, geometry, params, nElems )
    !> Array of particles
    type(mus_particle_group_type), intent(in) :: particleGroup
    !> scheme for access to varSys
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Number of elements including halos (multilevel not supported yet so no ghosts)
    integer, intent(in) :: nElems
    ! ------------------------------------- !
    integer :: iElem, vol_frac_pos, elemOff
    integer :: iParticle, lev, myRank
    ! ------------------------------------- !
    myRank = params%general%proc%rank
    lev = geometry%tree%global%maxLevel
    vol_frac_pos = scheme%varSys%method%val(scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)

    ! First reset the fluid volume fraction field to all 1.0
    do iElem = 1, nElems
      elemOff = (iElem-1)*scheme%varSys%nAuxScalars
      scheme%auxField(lev)%val(elemOff+vol_frac_pos) = 1.0_rk
    end do

    ! Loop over particles and distribute particle volume over nearby
    ! lattice sites.
    !$OMP PARALLEL
    !$OMP DO
    do iParticle = 1, particleGroup%particles_DPS%nvals
      call distributeParticleVolume( particle     = particleGroup               &
                                   &                %particles_DPS              &
                                   &                %val(iParticle),            &
                                   & interpolator = particleGroup%interpolator, &
                                   & scheme       = scheme,                     &
                                   & geometry     = geometry,                   &
                                   & params       = params                      )
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine mus_particles_updateFluidVolumeFraction

  !> Routine to initialize the fluid volume fraction field
  subroutine mus_particles_initFluidVolumeFraction( scheme, geometry, nElems )
    !> scheme for access to varSys
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Number of elements including halos (multilevel not supported yet so no ghosts)
    integer, intent(in) :: nElems
    ! ------------------------------------- !
    integer :: iElem, vol_frac_pos, elemOff
    integer :: lev
    ! ------------------------------------- !
    lev = geometry%tree%global%maxLevel
    vol_frac_pos = scheme%varSys%method%val(scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)

    ! First reset the fluid volume fraction field to all 1.0
    do iElem = 1, nElems
      elemOff = (iElem-1)*scheme%varSys%nAuxScalars
      scheme%auxField(lev)%val(elemOff+vol_frac_pos) = 1.0_rk
    end do

  end subroutine mus_particles_initFluidVolumeFraction

  !> Routine to map particles to the lattice (update coordinate of the origin etc.)
  subroutine mapToLattice_DPS(particle, interpolator, scheme, geometry, params, &
                             & comm, particlelogInterval  )
    !> Array of particles
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Interpolation data type. We don't interpolate stuff here but we use it to
    !! determine the directions we need to search for local elements in
    !! (only x, y for d2q9, x, y, z for d3q19)
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Communication data type for particles
    type(mus_particles_communication_type), intent(inout) :: comm
    !> Log the particle data to file every this many iterations
    integer, intent(in) :: particleLogInterval
    ! ---------------------------------------------------------------- !

    ! Update coordOfOrigin and particle owner
    call updateCoordOfOrigin_DPS( this =  particle,                   &
                                & scheme = scheme,                    &
                                & geometry = geometry,                &
                                & myRank = params%general%proc%rank   )

    ! Store previous owner of particle, then update the owner
    particle%previousOwner = particle%owner

    call updateParticleOwner( this     = particle,                       &
                            & scheme   = scheme,                         &
                            & geometry = geometry,                       &
                            & myrank   = params%general%proc%rank,       &
                            & procs    = comm%proc,                      &
                            & nProcs   = comm%nProcs                     )

    ! Reset the existsOnProc masks, then update them
    call updateExistsOnProc_DPS( this     = particle,                    &
                               & interpolator = interpolator,            &
                               & scheme = scheme,                        &
                               & geometry = geometry,                    &
                               & procs    = comm%proc,                   &
                               & nProcs   = comm%nProcs,                 &
                               & myRank   = params%general%proc%rank     )

  end subroutine mapToLattice_DPS

  !> Main control routine which maps each DPS particle to the lattice
  subroutine mapParticlesToLattice_DPS(particleGroup, scheme, geometry, params)
    !> Array of particles
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    ! -----------------------------------!
    integer :: iParticle
    ! -----------------------------------!

    !$OMP PARALLEL
    !$OMP DO
    do iParticle = 1, particleGroup%particles_DPS%nvals

      call mapToLattice_DPS( particle = particleGroup%particles_DPS%val(iParticle),   &
                           & interpolator = particleGroup%interpolator,               &
                           & scheme   = scheme,                                       &
                           & geometry = geometry,                                     &
                           & params   = params,                                       &
                           & comm     = particleGroup%send,                           &
                           & particleLogInterval = particleGroup%particleLogInterval  )


    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! -------------- Exchange new particles with neighboring processes --------------!
    call exchangeNewParticles_DPS( this = particleGroup,              &
                                 & send = particleGroup%send,         &
                                 & recv = particleGroup%recv,         &
                                 & comm = params%general%proc%comm,   &
                                 & myRank = params%general%proc%rank, &
                                 & message_flag = 1                   )

    ! Initialize the new particles I have received
    do iParticle = 1, particleGroup%particles_DPS%nvals
      if( particleGroup%particles_DPS%val(iParticle)%newForMe ) then


        call allocateProcessMasks(                                    &
          &    particle = particleGroup%particles_DPS%val(iParticle), &
          &    nProcs   = particleGroup%send%nProcs                   )

        call initParticle_DPS( &
          &    particle     = particleGroup%particles_DPS%val(iParticle), &
          &    interpolator = particleGroup%interpolator,                 &
          &    particleID   = particleGroup%particles_DPS%val(iParticle)  &
          &                                %particleID,                   &
          &    geometry     = geometry,                                   &
          &    scheme       = scheme,                                     &
          &    myRank       = params%general%proc%rank,                   &
          &    comm         = particleGroup%send                          )

        particleGroup%particles_DPS%val(iParticle)%newForMe = .FALSE.

        ! write(logUnit(1),*) 'mus_particles_mapping_DPS: initialized new particle on proc', &
        !   & params%general%proc%rank

         open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', &
           &   position='append' )
           write(pgDebugLog%lu,*) 'Initializing new particle with ID ', &
             & particleGroup%particles_DPS%val(iParticle)%particleID, &
             & ' iter = ', params%general%simcontrol%now%iter
           write(pgDebugLog%lu,*) 'pos = ',                                  &
             &                    particleGroup%particles_DPS%val(iParticle) &
             &                                 %pos(1:6)
           write(pgDebugLog%lu,*) 'vel = ',                                  &
             &                    particleGroup%particles_DPS%val(iParticle) &
             &                                 %vel(1:6)
           write(pgDebugLog%lu,*) 'F = ',                                    &
             &                    particleGroup%particles_DPS%val(iParticle) &
             &                                 %F(1:6)
           write(pgDebugLog%lu,*) 'F_DEM(1,:) = ',                           &
             &                    particleGroup%particles_DPS%val(iParticle) &
             &                                 %F_DEM(1,1:6)
           write(pgDebugLog%lu,*) 'F_DEM(2,:) = ',                           &
             &                    particleGroup%particles_DPS%val(iParticle) &
             &                                 %F_DEM(2,1:6)
           write(pgDebugLog%lu,*) 'coordOfOrigin = ',                        &
             &                    particleGroup%particles_DPS%val(iParticle) &
             &                                 %coordOfOrigin(1:4)
         close( pgDebugLog%lu )

      end if
    end do

    ! -------------- Exchange particles that have left the global domain --------------!
    call exchangeParticlesToRemove_DPS( this         = particleGroup,            &
                                      & send         = particleGroup%send,       &
                                      & recv         = particleGroup%recv,       &
                                      & comm         = params%general%proc%comm, &
                                      & myRank       = params%general%proc%rank, &
                                      & message_flag = 1                         )

    ! Remove particles that have left the global domain or which should be removed
    ! locally from this process
    iParticle = 1
    do while ( iParticle <= particleGroup%particles_DPS%nvals )
      if( particleGroup%particles_DPS%val(iParticle)%removeParticle_global ) then
        open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
          write(pgDebugLog%lu,*) 'Destroying particle with ID ', &
            & particleGroup%particles_DPS%val(iParticle)%particleID, &
            & ' iter = ', params%general%simcontrol%now%iter, 'from GLOBAL domain'
        close( pgDebugLog%lu )

        call remove_particle_from_da_particle_DPS( particleGroup%particles_DPS, iParticle)
        cycle
      else if ( particleGroup%particles_DPS%val(iParticle)%removeParticle_local ) then
        open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
          write(pgDebugLog%lu,*) 'Destroying particle with ID ', &
            & particleGroup%particles_DPS%val(iParticle)%particleID, &
            & ' iter = ', params%general%simcontrol%now%iter, 'from LOCAL domain'
        close( pgDebugLog%lu )

        call remove_particle_from_da_particle_DPS( particleGroup%particles_DPS, iParticle)
        cycle
      end if
      iParticle = iParticle + 1
    end do

  end subroutine mapParticlesToLattice_DPS

  !> Main control routine which interpolates and stores the fluid properties at
  !! the particle locations
  subroutine interpolateFluidProperties_DPS(particleGroup, scheme, geometry, params)
    !> Array of particles
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    ! --------------------------------------------- !
    integer :: iParticle, myRank
    ! --------------------------------------------- !
    myRank = params%general%proc%rank
    !$OMP PARALLEL
    !$OMP DO
    do iParticle = 1, particleGroup%particles_DPS%nvals
      if(particleGroup%particles_DPS%val(iParticle)%owner == myRank) then
      ! Owner of particle interpolates fluid properties from auxField to particle position
      call mus_particles_DPS_interpolateFluidProperties(                    &
                  & particle     = particleGroup                            &
                  &                %particles_DPS                           &
                  &                %val(iParticle),                         &
                  & interpolator = particleGroup%interpolator,              &
                  & scheme       = scheme,                                  &
                  & geometry     = geometry,                                &
                  & params       = params,                                  &
                  & intp         = particleGroup%intp,                      &
                  & calc_vel_and_p_grad = particleGroup%calc_vel_and_p_grad )
       end if
     end do
     !$OMP END DO
     !$OMP END PARALLEL
  end subroutine interpolateFluidProperties_DPS

  !> Interpolates fluid properties from neighboring lattice sites to determine the
  !! fluid density and velocity at the location of the particle
  subroutine mus_particles_DPS_interpolateFluidProperties( particle, interpolator, &
                                                         & scheme, geometry, params, &
                                                         & intp,  &
                                                         & calc_vel_and_p_grad )
    !> Particle to interpolate fluid properties to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Interpolation object containing interpolation bounds and weigth function info
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> Procedure to interpolate fluid props (different for Generalized Navier-Stokes)
    !! and regular Navier-Stokes
    procedure(interpolateFluidPropFunc) :: intp
    !> Procedure to interpolate fluid props (different for Generalized Navier-Stokes)
    !! and regular Navier-Stokes
    procedure(calcVelAndPGradFunc) :: calc_vel_and_p_grad
    ! ------------------------------------------!
    integer :: lev
    real(kind=rk) :: vel_tmp(3), rho_tmp, eps_f_tmp
    real(kind=rk) :: grad_p_tmp(3), curl_u_tmp(3)
    real(kind=rk) :: dx
    logical :: failedToGrabValue
    ! ------------------------------------------!
    failedToGrabValue = .FALSE.
    lev = particle%coordOfOrigin(4)
    dx = params%physics%dxLvl(lev)

    call intp( xp           = particle%pos(1:3),           &
             & coord_xp     = particle%coordOfOrigin,      &
             & scheme       = scheme,                      &
             & geom_origin  = geometry%tree%global%origin, &
             & dx           = dx,                          &
             & interpolator = interpolator,                &
             & vel_xp       = vel_tmp,                     &
             & rho_xp       = rho_tmp,                     &
             & eps_f_xp     = eps_f_tmp,                   &
             & posOfCoord   = particle%posOfOrigin         )

    ! Due to rounding errors, for small particles eps_f may turn out greater than one
    ! In that case, set it equal to one
    if( eps_f_tmp > 1.0_rk ) then
      eps_f_tmp = 1.0_rk
    end if

    ! Convert to physical units
    particle%rho_fluid   = params%physics%rho0
    particle%u_fluid     = vel_tmp * params%physics%fac(lev)%vel
    particle%eps_f_fluid = eps_f_tmp

    ! Calculate gradient of pressure and curl of velocity
    ! These are used
    ! Right now we assume grad(p) and curl(u) at the particle location are
    ! equal to the value at the barycenter of its coordOfOrigin.
    ! This can be changed to interpolation later but for that we would need
    ! non-local values of rho and u.

    call calc_vel_and_p_grad( coord    = particle%coordOfOrigin, &
                            & scheme   = scheme,                 &
                            & grad_p   = grad_p_tmp,             &
                            & curl_u   = curl_u_tmp,             &
                            & err      = failedToGrabValue,      &
                            & posOfCoord = particle%posOfOrigin  )
    if(failedToGrabValue) then
      write(logUnit(1),*) "ERROR interpolateFluidProperties: could not grab value of gradients!"
      write(logUnit(1),*) "CoordOfOrigin = ", particle%coordOfOrigin
      write(logUnit(1),*) "Particle owner = ", particle%owner
      write(logUnit(1),*) "myRank = ", params%general%proc%rank

      call tem_abort()
    end if

    ! Convert to physical units
    particle%grad_p_fluid(1:3) = grad_p_tmp(1:3) * params%physics%fac(lev)%press / dx
    particle%curl_u_fluid(1:3) = curl_u_tmp(1:3) * params%physics%fac(lev)%vel / dx

  end subroutine mus_particles_DPS_interpolateFluidProperties

  !> Interpolate the fluid properties to a point xp
  subroutine interpolateFluidProps(xp, coord_xp, scheme, geom_origin, dx, &
                                        & interpolator, vel_xp, rho_xp, eps_f_xp, posOfCoord )
    !> Query point (x,y,z) to interpolate fluid properties to
    real(kind=rk), intent(in) :: xp(3)
    !> Coordinate in the tree (ix,iy,iz,level) of xp
    integer, intent(in) :: coord_xp(4)
    !> Scheme for access to fluid data
    type(mus_scheme_type), intent(in) :: scheme
    !> Origin of bounding cube (geometry%tree%global%origin)
    real(kind=rk), intent(in) :: geom_origin(3)
    !> Mesh size
    real(kind=rk), intent(in) :: dx
    !> Interpolator object containing stencil and weight function information
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Fluid velocity at xp
    real(kind=rk), intent(out) :: vel_xp(3)
    !> Fluid density at xp
    real(kind=rk), intent(out) :: rho_xp
    !> Fluid volume fraction at xp
    real(kind=rk), intent(out) :: eps_f_xp
    !> Position of coord_xp in the total list, if available.
    !! Supplying this makes the routine run faster.
    integer, intent(in), optional :: posOfCoord
    ! ------------------------------------------!
    integer :: idir, lev
    integer :: coord(4)
    integer :: nghPos, elemPos, elemOff
    integer :: vel_pos(3), dens_pos, vol_frac_pos
    real(kind=rk) :: del_x, del_y, del_z, wght
    real(kind=rk) :: r_lat(3), bary(3)
    real(kind=rk) :: u_tmp(3), rho_tmp, eps_f_tmp
    ! ------------------------------------------!
    lev = coord_xp(4)
    vol_frac_pos = scheme%varSys%method%val(scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)
    dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
    vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    ! Get position of particle coordinate in the total list
    if( present(posOfCoord) ) then
      elemPos = posOfCoord
    else
      elemPos = getPosOfCoord( coord = coord_xp, scheme = scheme )
    end if

    ! Initialize output variables to 0
    vel_xp = 0.0_rk
    rho_xp = 0.0_rk
    eps_f_xp = 0.0_rk

    ! Loop over sample points (= neigboring lattice cells)
    do idir = 1, interpolator%Nelems
      coord(1:3) = coord_xp(1:3) + interpolator%neighDirs(1:3,idir)
      coord(4) = coord_xp(4)
      ! Get the position of this neighbor element in the total list
      ! First try to do this by following the interpolator's neighPaths
      nghPos = followNeighPath(                                        &
                & neighpath    = interpolator%neighPaths(1:3,idir),    &
                & startPos     = elemPos,                              &
                & neigh        = scheme%levelDesc(lev)%neigh(1),       &
                & nElems_fluid = scheme%pdf(lev)%nElems_fluid,         &
                & zeroDir      = scheme%layout%stencil(1)%restPosition )

      if(nghPos < 1) then
        ! If we cannot follow neighPath to get this element's position, get its position
        ! using (expensive) call to getPosOfCoord which uses the tem_topology routines
        nghPos = getPosOfCoord( coord = coord, scheme = scheme )
      end if

      if(nghPos < 1) then
        ! If neighPos is still less than 0, use the value at the current element instead of
        ! its neighbor in the interpolation. This can be the case for example near walls.
        elemOff = (elemPos-1)*scheme%varSys%nAuxScalars
      else
        ! Grab the flow variables from the neighboring lattice site.
        elemOff = (nghPos-1)*scheme%varSys%nAuxScalars

      end if ! neighPos < 0

      rho_tmp = scheme%auxField(lev)%val( elemOff + dens_pos )
      u_tmp(1) = scheme%auxField(lev)%val( elemOff + vel_pos(1) )
      u_tmp(2) = scheme%auxField(lev)%val( elemOff + vel_pos(2) )
      u_tmp(3) = scheme%auxField(lev)%val( elemOff + vel_pos(3) )
      eps_f_tmp = scheme%auxField(lev)%val( elemOff + vol_frac_pos )

      ! Compute weights
      bary = getBaryOfCoord( coord  = coord(1:3),  &
                           & origin = geom_origin, &
                           & dx     = dx           )

      ! Get normalized distance from query point to current sample point
      r_lat = abs( (xp - bary)/dx )

      del_x = interpolator%getWght_x( r_lat(1) )
      del_y = interpolator%getWght_y( r_lat(2) )
      del_z = interpolator%getWght_z( r_lat(3) )
      wght  = del_x*del_y*del_z

      ! Add contribution of this value to vel_xp, rho_xp, eps_xp
      vel_xp = vel_xp + wght*u_tmp
      rho_xp = rho_xp + wght*rho_tmp
      eps_f_xp = eps_f_xp + wght*eps_f_tmp

    end do ! idir
  end subroutine interpolateFluidProps

  !> Interpolate the fluid properties to a point xp
  subroutine interpolateFluidProps_onewaycoupled(xp, coord_xp, scheme, geom_origin, dx, &
                                        & interpolator, vel_xp, rho_xp, eps_f_xp, posOfCoord )
    !> Query point (x,y,z) to interpolate fluid properties to
    real(kind=rk), intent(in) :: xp(3)
    !> Coordinate in the tree (ix,iy,iz,level) of xp
    integer, intent(in) :: coord_xp(4)
    !> Scheme for access to fluid data
    type(mus_scheme_type), intent(in) :: scheme
    !> Origin of bounding cube (geometry%tree%global%origin)
    real(kind=rk), intent(in) :: geom_origin(3)
    !> Mesh size
    real(kind=rk), intent(in) :: dx
    !> Interpolator object containing stencil and weight function information
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Fluid velocity at xp
    real(kind=rk), intent(out) :: vel_xp(3)
    !> Fluid density at xp
    real(kind=rk), intent(out) :: rho_xp
    !> Fluid volume fraction at xp
    real(kind=rk), intent(out) :: eps_f_xp
    !> Position of coord_xp in the total list, if available.
    !! Supplying this makes the routine run faster.
    integer, intent(in), optional :: posOfCoord
    ! ------------------------------------------!
    integer :: idir, lev
    integer :: coord(4)
    integer :: nghPos, elemPos, elemOff
    integer :: vel_pos(3), dens_pos
    real(kind=rk) :: del_x, del_y, del_z, wght
    real(kind=rk) :: r_lat(3), bary(3)
    real(kind=rk) :: u_tmp(3), rho_tmp
    ! ------------------------------------------!
    lev = coord_xp(4)
    dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
    vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    ! Get position of particle coordinate in the total list
    if( present(posOfCoord) ) then
      elemPos = posOfCoord
    else
      elemPos = getPosOfCoord( coord = coord_xp, scheme = scheme )
    end if

    ! Initialize output variables to 0
    vel_xp = 0.0_rk
    rho_xp = 0.0_rk
    eps_f_xp = 0.0_rk

    ! Loop over sample points (= neigboring lattice cells)
    do idir = 1, interpolator%Nelems
      coord(1:3) = coord_xp(1:3) + interpolator%neighDirs(1:3,idir)
      coord(4) = coord_xp(4)
      ! Get the position of this neighbor element in the total list
      ! First try to do this by following the interpolator's neighPaths
      nghPos = followNeighPath(                                        &
                & neighpath    = interpolator%neighPaths(1:3,idir),    &
                & startPos     = elemPos,                              &
                & neigh        = scheme%levelDesc(lev)%neigh(1),       &
                & nElems_fluid = scheme%pdf(lev)%nElems_fluid,         &
                & zeroDir      = scheme%layout%stencil(1)%restPosition )

      if(nghPos < 1) then
        ! If we cannot follow neighPath to get this element's position, get its position
        ! using (expensive) call to getPosOfCoord which uses the tem_topology routines
        nghPos = getPosOfCoord( coord = coord, scheme = scheme )
      end if

      if(nghPos < 1) then
        ! If neighPos is still less than 0, use the value at the current element instead of
        ! its neighbor in the interpolation. This can be the case for example near walls.
        elemOff = (elemPos-1)*scheme%varSys%nAuxScalars
      else
        ! Grab the flow variables from the neighboring lattice site.
        elemOff = (nghPos-1)*scheme%varSys%nAuxScalars

      end if ! neighPos < 0

      rho_tmp = scheme%auxField(lev)%val( elemOff + dens_pos )
      u_tmp(1) = scheme%auxField(lev)%val( elemOff + vel_pos(1) )
      u_tmp(2) = scheme%auxField(lev)%val( elemOff + vel_pos(2) )
      u_tmp(3) = scheme%auxField(lev)%val( elemOff + vel_pos(3) )

      ! Compute weights
      bary = getBaryOfCoord( coord  = coord(1:3),  &
                           & origin = geom_origin, &
                           & dx     = dx           )

      ! Get normalized distance from query point to current sample point
      r_lat = abs( (xp - bary)/dx )

      del_x = interpolator%getWght_x( r_lat(1) )
      del_y = interpolator%getWght_y( r_lat(2) )
      del_z = interpolator%getWght_z( r_lat(3) )
      wght  = del_x*del_y*del_z

      ! Add contribution of this value to vel_xp, rho_xp, eps_xp
      vel_xp = vel_xp + wght*u_tmp
      rho_xp = rho_xp + wght*rho_tmp

    end do ! idir

    ! Set fluid volume fraction to 1 for one way coupled simulations
    eps_f_xp = 1.0_rk
  end subroutine interpolateFluidProps_onewaycoupled

  subroutine interpolateFluidProps_onewaycoupled_old(xp, coord_xp, scheme, geom_origin, dx, &
    & interpolator, vel_xp, rho_xp )
  !> Query point (x,y,z) to interpolate fluid properties to
  real(kind=rk), intent(in) :: xp(3)
  !> Coordinate in the tree (ix,iy,iz,level) of xp
  integer, intent(in) :: coord_xp(4)
  !> Scheme for access to fluid data
  type(mus_scheme_type), intent(inout) :: scheme
  !> Origin of bounding cube (geometry%tree%global%origin)
  real(kind=rk), intent(in) :: geom_origin(3)
  !> Mesh size
  real(kind=rk), intent(in) :: dx
  !> Interpolator object containing stencil and weight function information
  type(mus_particle_interpolator_type), intent(in) :: interpolator
  !> Fluid velocity at xp
  real(kind=rk), intent(out) :: vel_xp(3)
  !> Fluid density at xp
  real(kind=rk), intent(out) :: rho_xp
  ! ------------------------------------------!
  integer :: nx, ny, nz
  integer :: coord(4)
  integer :: vel_pos(3), dens_pos
  real(kind=rk) :: del_x, del_y, del_z, wght
  real(kind=rk) :: r_lat(3), bary(3)
  real(kind=rk) :: u_tmp(3), rho_tmp
  logical :: failedToGrabValue
  ! ------------------------------------------!
  dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
  vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

  ! Initialize output variables to 0
  vel_xp = 0.0_rk
  rho_xp = 0.0_rk

  ! Loop over sample points (= neigboring lattice cells)
  do nx = interpolator%bnd_x(1), interpolator%bnd_x(2)
    do ny = interpolator%bnd_y(1), interpolator%bnd_y(2)
      do nz = interpolator%bnd_z(1), interpolator%bnd_z(2)
      coord = getNeighborCoord( coord        = coord_xp, &
      & nx           = nx,               &
      & ny           = ny,               &
      & nz           = nz,               &
      & boundaryData = pgBndData               )

      bary = getBaryOfCoord( coord  = coord(1:3), &
      & origin = geom_origin,   &
      & dx     = dx             )

      ! Get normalized distance from query point to current sample point
      r_lat = abs( (xp - bary)/dx )

      ! Compute weights
      del_x = interpolator%getWght_x( r_lat(1) )
      del_y = interpolator%getWght_y( r_lat(2) )
      del_z = interpolator%getWght_z( r_lat(3) )
      wght  = del_x*del_y*del_z

      call grabValueAtCoord( coord        = coord,            &
      & scheme       = scheme,           &
      & vel_pos      = vel_pos,          &
      & dens_pos     = dens_pos,         &
      & rho          = rho_tmp,          &
      & u            = u_tmp,            &
      & err          = failedToGrabValue )

      if(failedToGrabValue) then
        call grabValueAtCoord( coord        = coord_xp,         &
        & scheme       = scheme,           &
        & vel_pos      = vel_pos,          &
        & dens_pos     = dens_pos,         &
        & rho          = rho_tmp,          &
        & u            = u_tmp,            &
        & err          = failedToGrabValue )
        if(failedToGrabValue) then
          write(logUnit(1),*) 'ERROR interpolateFluidProps_delta: ', &
          & 'could not grab fluid prop values at coord', coord_xp
          call tem_abort()
        end if

      end if

      ! Add contribution of this value to vel_xp, rho_xp, eps_xp
      vel_xp = vel_xp + wght*u_tmp
      rho_xp = rho_xp + wght*rho_tmp

      end do
    end do
  end do
  end subroutine interpolateFluidProps_onewaycoupled_old

  !> Interpolates fluid properties from neighboring lattice sites to determine the
  !! fluid density and velocity at the location of the particle
  subroutine mus_particles_DPS_interpolateFluidProperties_onewaycoupled( particle, interpolator, &
                                                         & scheme, geometry, params)
    !> Particle to interpolate fluid properties to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Interpolation object containing interpolation bounds and weigth function info
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    !> Scheme
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    ! ------------------------------------------!
    integer :: lev
    real(kind=rk) :: vel_tmp(3), rho_tmp
    real(kind=rk) :: grad_p_tmp(3), curl_u_tmp(3)
    real(kind=rk) :: dx
    logical :: failedToGrabValue
    ! ------------------------------------------!
    failedToGrabValue = .FALSE.
    lev = particle%coordOfOrigin(4)
    dx = params%physics%dxLvl(lev)

    call interpolateFluidProps_onewaycoupled_old(                       &
                      & xp           = particle%pos(1:3),           &
                      & coord_xp     = particle%coordOfOrigin,      &
                      & scheme       = scheme,                      &
                      & geom_origin  = geometry%tree%global%origin, &
                      & dx           = dx,                          &
                      & interpolator = interpolator,                &
                      & vel_xp       = vel_tmp,                     &
                      & rho_xp       = rho_tmp                      )

    ! Convert to physical units
    particle%rho_fluid   = rho_tmp * params%physics%rho0
    ! Note the division by eps_f because we want to interpolate the ACTUAL fluid velocity,
    ! whereas in auxField the VOLUME-AVERAGED fluid velocity is stored.
    particle%u_fluid     = (vel_tmp * params%physics%fac(lev)%vel)

    ! Calculate gradient of pressure and curl of velocity
    ! These are used
    ! Right now we assume grad(p) and curl(u) at the particle location are
    ! equal to the value at the barycenter of its coordOfOrigin.
    ! This can be changed to interpolation later but for that we would need
    ! non-local values of rho and u.

    call calcVelocityAndPressureGradient_onewaycoupled_old(        &
                          & coord    = particle%coordOfOrigin, &
                          & scheme   = scheme,                 &
                          & grad_p   = grad_p_tmp,             &
                          & curl_u   = curl_u_tmp,             &
                          & err      = failedToGrabValue       )

    if(failedToGrabValue) then
      write(logUnit(1),*) "ERROR interpolateFluidProperties: could not grab value of gradients!"
      write(logUnit(1),*) "CoordOfOrigin = ", particle%coordOfOrigin
      write(logUnit(1),*) "Particle owner = ", particle%owner
      write(logUnit(1),*) "myRank = ", params%general%proc%rank

      call tem_abort()
    end if

    ! Convert to physical units
    particle%grad_p_fluid(1:3) = grad_p_tmp(1:3) * params%physics%fac(lev)%press / dx
    particle%curl_u_fluid(1:3) = curl_u_tmp(1:3) * params%physics%fac(lev)%vel / dx

  end subroutine mus_particles_DPS_interpolateFluidProperties_onewaycoupled


  !> calcVelocityAndPressureGradient calculates the gradient of pressure and
  !! curl of the velocity field (in lattice units) at coordinate coord
  !! Coord must be a local element on this process!
  subroutine calcVelocityAndPressureGradient(coord, scheme, grad_p, curl_u, err, posOfCoord )
    !> Integer coordinate to calculate grad(p) and curl(u) at
    integer, intent(in) :: coord(4)
    !> Scheme for access to total list
    type(mus_scheme_type), intent(in) :: scheme
    !> Output pressure gradient at coord
    real(kind=rk), intent(out) :: grad_p(3)
    !> Output velocity curl at coord
    real(kind=rk), intent(out) :: curl_u(3)
    !> Error code: set to TRUE if we could not find ldPos of coord
    !! we are trying to grab value at (for example because it is
    !! outside domain and not a halo.
    logical, intent(out) :: err
    !> Optional: position of coord in the total list. This can speed up execution
    !! of this routine if this is known.
    integer, intent(in), optional :: posOfCoord
    ! ------------------------------------------!
    integer :: dens_pos, vel_pos(3), vol_frac_pos
    integer :: elemPos, nghPos, elemOff, lev, Nelems_fluid
    integer :: iDir, QQN
    integer :: invDir
    integer :: cq(3)
    real(kind=rk) :: wq
    real(kind=rk) :: rho_tmp, u_tmp(3), curl_tmp(3)
    real(kind=rk) :: rho_coord, u_coord(3)
    real(kind=rk) :: eps_f_inv         ! 1 / fluid volume fraction
    ! ------------------------------------------!
    lev = coord(4)
    QQN = scheme%layout%fStencil%QQN
    Nelems_fluid = scheme%pdf(lev)%nElems_fluid

    ! Position of volume frac in varSys
    vol_frac_pos = scheme%varSys%method%val(scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)
    dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
    vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    ! Get fluid volume fraction at coord
    if( present(posOfCoord) ) then
      elemPos = posOfCoord
    else
      call getIndexOfCoordInTotal( coord  = coord,  &
                                & scheme = scheme, &
                                & ldPos  = elemPos   )
    end if

    ! If this element is not a local fluid on our process domain, skip it
    if( elemPos <= 0 .OR. elemPos > scheme%pdf(lev)%nElems_fluid) then
      err = .TRUE.
      return
    else
      err = .FALSE.
    end if

    elemOff = (elemPos-1)*scheme%varSys%nAuxScalars

    ! Get flow variables at the location we are calculating grad p and curl u at.
    eps_f_inv = 1.0_rk / scheme%auxField(lev)%val(elemOff + vol_frac_pos)
    rho_coord = scheme%auxField(lev)%val( elemOff + dens_pos )
    u_coord(1) = scheme%auxField(lev)%val( elemOff + vel_pos(1) )
    u_coord(2) = scheme%auxField(lev)%val( elemOff + vel_pos(2) )
    u_coord(3) = scheme%auxField(lev)%val( elemOff + vel_pos(3) )

    ! Loop over points around coord and grab the values of velocity
    ! and pressure at those points
    grad_p = 0.0_rk
    curl_u = 0.0_rk

    do iDir = 1, QQN
      ! Get lattice direction and weight
      cq = scheme%layout%fStencil%cxDir(:,iDir)
      wq = scheme%layout%weight(iDir)

      ! Get position of neighboring lattice site in total list
      nghPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, elemPos)

      ! If nghPos is positive, we have a neighboring fluid element in this direction
      if(nghPos > 0) then
        ! Grab the values of pressure and density
        elemOff = (nghPos-1)*scheme%varSys%nAuxScalars
        rho_tmp = scheme%auxField(lev)%val( elemOff + dens_pos )
        u_tmp(1) = scheme%auxField(lev)%val( elemOff + vel_pos(1) )
        u_tmp(2) = scheme%auxField(lev)%val( elemOff + vel_pos(2) )
        u_tmp(3) = scheme%auxField(lev)%val( elemOff + vel_pos(3) )
      else
        ! If nghPos is negative, we have a boundary in this direction. In that case,
        ! we use the neighbor in the OPPOSITE direction to extrapolate the values we need.
        invDir = scheme%layout%fStencil%cxDirInv(iDir)
        nghPos = scheme%levelDesc(lev)%neigh(1)%nghElems(invDir, elemPos)

        if(nghPos < 1) then
          ! If the element in opposite direction is also not available,
          ! use the value at the current element
          elemOff = (elemPos-1)*scheme%varSys%nAuxScalars
        else
          elemOff = (nghPos-1)*scheme%varSys%nAuxScalars
        end if

        ! Grab the values of pressure and density at the
        rho_tmp = scheme%auxField(lev)%val( elemOff + dens_pos )
        u_tmp(1) = scheme%auxField(lev)%val( elemOff + vel_pos(1) )
        u_tmp(2) = scheme%auxField(lev)%val( elemOff + vel_pos(2) )
        u_tmp(3) = scheme%auxField(lev)%val( elemOff + vel_pos(3) )

        rho_tmp = rho_coord - (rho_tmp - rho_coord)
        u_tmp   = u_coord - ( u_tmp - u_coord )
      end if

      ! Add this coordinate's contribution to the gradients
      ! according to the lattice differential operators
      ! note that rho_tmp = cs^2 * p. This cancels out with the
      ! division by cs^2 in the definition of the lattice differential
      ! operator, which is why cs^2 is not in this term
      grad_p = grad_p + wq * cq * rho_tmp

      ! Compute the cross-product term, then add it to curl_u
      ! Convert the volume-averaged fluid velocity (= what grabValueAtCoord gives)
      ! to the fluid phase velocity for Generalized Navier-Stokes
      call cross_product( wq*cq, u_tmp, curl_tmp )
      curl_u = curl_u + curl_tmp
    end do ! iDir

    curl_u = curl_u * cs2inv * eps_f_inv

  end subroutine calcVelocityAndPressureGradient


  !> calcVelocityAndPressureGradient calculates the gradient of pressure and
  !! curl of the velocity field (in lattice units) at coordinate coord
  !! Coord must be a local element on this process!
  subroutine calcVelocityAndPressureGradient_onewaycoupled(coord, scheme, grad_p, curl_u, err, posOfCoord )
    !> Integer coordinate to calculate grad(p) and curl(u) at
    integer, intent(in) :: coord(4)
    !> Scheme for access to total list
    type(mus_scheme_type), intent(in) :: scheme
    !> Output pressure gradient at coord
    real(kind=rk), intent(out) :: grad_p(3)
    !> Output velocity curl at coord
    real(kind=rk), intent(out) :: curl_u(3)
    !> Error code: set to TRUE if we could not find ldPos of coord
    !! we are trying to grab value at (for example because it is
    !! outside domain and not a halo.
    logical, intent(out) :: err
    !> Optional: position of coord in the total list. This can speed up execution
    !! of this routine if this is known.
    integer, intent(in), optional :: posOfCoord
    ! ------------------------------------------!
    integer :: dens_pos, vel_pos(3)
    integer :: elemPos, nghPos, elemOff, lev, Nelems_fluid
    integer :: iDir, QQN
    integer :: invDir
    integer :: cq(3)
    real(kind=rk) :: wq
    real(kind=rk) :: rho_tmp, u_tmp(3), curl_tmp(3)
    real(kind=rk) :: rho_coord, u_coord(3)
    ! ------------------------------------------!
    lev = coord(4)
    QQN = scheme%layout%fStencil%QQN
    Nelems_fluid = scheme%pdf(lev)%nElems_fluid

    ! Position of volume frac in varSys
    dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
    vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    ! Get fluid volume fraction at coord
    if( present(posOfCoord) ) then
      elemPos = posOfCoord
    else
      call getIndexOfCoordInTotal( coord  = coord,  &
                                & scheme = scheme, &
                                & ldPos  = elemPos   )
    end if

    ! If this element is not a local fluid on our process domain, skip it
    if( elemPos <= 0 .OR. elemPos > scheme%pdf(lev)%nElems_fluid) then
      err = .TRUE.
      return
    else
      err = .FALSE.
    end if

    elemOff = (elemPos-1)*scheme%varSys%nAuxScalars

    ! Get flow variables at the location we are calculating grad p and curl u at.
    rho_coord = scheme%auxField(lev)%val( elemOff + dens_pos )
    u_coord(1) = scheme%auxField(lev)%val( elemOff + vel_pos(1) )
    u_coord(2) = scheme%auxField(lev)%val( elemOff + vel_pos(2) )
    u_coord(3) = scheme%auxField(lev)%val( elemOff + vel_pos(3) )

    ! Loop over points around coord and grab the values of velocity
    ! and pressure at those points
    grad_p = 0.0_rk
    curl_u = 0.0_rk

    do iDir = 1, QQN
      ! Get lattice direction and weight
      cq = scheme%layout%fStencil%cxDir(:,iDir)
      wq = scheme%layout%weight(iDir)

      ! Get position of neighboring lattice site in total list
      nghPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, elemPos)

      ! If nghPos is positive, we have a neighboring fluid element in this direction
      if(nghPos > 0) then
        ! Grab the values of pressure and density
        elemOff = (nghPos-1)*scheme%varSys%nAuxScalars
        rho_tmp = scheme%auxField(lev)%val( elemOff + dens_pos )
        u_tmp(1) = scheme%auxField(lev)%val( elemOff + vel_pos(1) )
        u_tmp(2) = scheme%auxField(lev)%val( elemOff + vel_pos(2) )
        u_tmp(3) = scheme%auxField(lev)%val( elemOff + vel_pos(3) )
      else
        ! If nghPos is negative, we have a boundary in this direction. In that case,
        ! we use the neighbor in the OPPOSITE direction to extrapolate the values we need.
        invDir = scheme%layout%fStencil%cxDirInv(iDir)
        nghPos = scheme%levelDesc(lev)%neigh(1)%nghElems(invDir, elemPos)

        if(nghPos < 1) then
          ! If the element in opposite direction is also not available,
          ! use the value at the current element
          elemOff = (elemPos-1)*scheme%varSys%nAuxScalars
        else
          elemOff = (nghPos-1)*scheme%varSys%nAuxScalars
        end if

        ! Grab the values of pressure and density at the
        rho_tmp = scheme%auxField(lev)%val( elemOff + dens_pos )
        u_tmp(1) = scheme%auxField(lev)%val( elemOff + vel_pos(1) )
        u_tmp(2) = scheme%auxField(lev)%val( elemOff + vel_pos(2) )
        u_tmp(3) = scheme%auxField(lev)%val( elemOff + vel_pos(3) )

        rho_tmp = rho_coord - (rho_tmp - rho_coord)
        u_tmp   = u_coord - ( u_tmp - u_coord )
      end if

      ! Add this coordinate's contribution to the gradients
      ! according to the lattice differential operators
      ! note that rho_tmp = cs^2 * p. This cancels out with the
      ! division by cs^2 in the definition of the lattice differential
      ! operator, which is why cs^2 is not in this term
      grad_p = grad_p + wq * cq * rho_tmp

      ! Compute the cross-product term, then add it to curl_u
      ! Convert the volume-averaged fluid velocity (= what grabValueAtCoord gives)
      ! to the fluid phase velocity for Generalized Navier-Stokes
      call cross_product( wq*cq, u_tmp, curl_tmp )
      curl_u = curl_u + curl_tmp
    end do ! iDir

    curl_u = curl_u * cs2inv

  end subroutine calcVelocityAndPressureGradient_onewaycoupled

  !> calcVelocityAndPressureGradient calculates the gradient of pressure and
  !! curl of the velocity field (in lattice units) at coordinate coord
  !! Coord must be a local element on this process!
  subroutine calcVelocityAndPressureGradient_onewaycoupled_old(coord, scheme, grad_p, curl_u, err )
    !> Integer coordinate to calculate grad(p) and curl(u) at
    integer, intent(in) :: coord(4)
    !> Scheme for access to total list
    type(mus_scheme_type), intent(inout) :: scheme
    !> Output pressure gradient at coord
    real(kind=rk), intent(out) :: grad_p(3)
    !> Output velocity curl at coord
    real(kind=rk), intent(out) :: curl_u(3)
    !> Error code: set to TRUE if we could not find ldPos of coord
    !! we are trying to grab value at (for example because it is
    !! outside domain and not a halo.
    logical, intent(out) :: err
    ! ------------------------------------------!
    integer :: dens_pos, vel_pos(3)
    integer :: ldPos, elemOff, lev
    integer :: iDir, QQ
    integer :: invDir
    integer :: cq(3), cqInv(3), neighborCoord(4)
    real(kind=rk) :: wq
    real(kind=rk) :: rho_tmp, u_tmp(3), curl_tmp(3)
    real(kind=rk) :: rho_coord, u_coord(3)
    logical :: failedToGrabValue
    ! ------------------------------------------!
    lev = coord(4)
    QQ = scheme%layout%fStencil%QQ

    ! Position of volume frac in varSys
    dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
    vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    ! Get fluid volume fraction at coord
    call getIndexOfCoordInTotal( coord  = coord,  &
                             & scheme = scheme, &
                             & ldPos  = ldPos   )

    ! If this element is not a local fluid on our process domain, skip it
    if( ldPos <= 0 .OR. ldPos > scheme%pdf(lev)%nElems_fluid) then
      err = .TRUE.
      return
    else
      err = .FALSE.
    end if

    elemOff = (ldPos-1)*scheme%varSys%nAuxScalars

    ! Loop over points around coord and grab the values of velocity
    ! and pressure at those points
    grad_p = 0.0_rk
    curl_u = 0.0_rk

    do iDir = 1, QQ
      ! Get lattice direction and weight
      cq = scheme%layout%fStencil%cxDir(:,iDir)
      wq = scheme%layout%weight(iDir)

      ! Calculate neighbor coord
      neighborCoord = getNeighborCoord( coord       = coord,    &
                                     & nx           = cq(1),    &
                                     & ny           = cq(2),    &
                                     & nz           = cq(3),    &
                                     & boundaryData = pgBndData )

      ! Grab value of density and velocity at neighbor
      call grabValueAtCoord( coord    = neighborCoord,    &
                           & scheme   = scheme,           &
                           & vel_pos  = vel_pos,          &
                           & dens_pos = dens_pos,         &
                           & rho      = rho_tmp,          &
                           & u        = u_tmp,            &
                           & err      = failedToGrabValue )
      if( failedToGrabValue ) then
        ! If we could not get density and velocity values at neighbor
        ! i.e. if neighbor is outside domain, use the neighbor in the
        ! OTHER direction to extrapolate the values we need.
        ! write(logUnit(1),*) 'calcVelocityAndPressureGradient: extrapolating ', &
        !   & 'from neighbor in opposite direction'

        ! First grab value at coord itself
        call grabValueAtCoord( coord    = coord,    &
                             & scheme   = scheme,           &
                             & vel_pos  = vel_pos,          &
                             & dens_pos = dens_pos,         &
                             & rho      = rho_coord,         &
                             & u        = u_coord,            &
                             & err      = failedToGrabValue )
        if( failedToGrabValue ) then
          write(logUnit(1),*) 'ERROR calcVelocityAndPressureGradient: ', &
            & 'could not grab pressure value at coord', coord
          call tem_abort()
        end if

        ! Then get the value at neighbor in OPPOSITE direction
        invDir = scheme%layout%fStencil%cxDirInv(iDir)
        cqInv = scheme%layout%fStencil%cxDir(:,invDir)


        neighborCoord = getNeighborCoord( coord       = coord,       &
                                       & nx           = cqInv(1),    &
                                       & ny           = cqInv(2),    &
                                       & nz           = cqInv(3),    &
                                       & boundaryData = pgBndData )

        ! Grab value of density and velocity at neighbor
        ! in INVERSE direction
        call grabValueAtCoord( coord    = neighborCoord,    &
                             & scheme   = scheme,           &
                             & vel_pos  = vel_pos,          &
                             & dens_pos = dens_pos,         &
                             & rho      = rho_tmp,          &
                             & u        = u_tmp,            &
                             & err      = failedToGrabValue )
        if( failedToGrabValue ) then
          ! write(logUnit(1),*) 'WARNING calcVelocityAndPressureGradient: ',  &
          !   & 'could not grab value at neighbor coord. Assuming value at ', &
          !   & 'coordOfOrigin'
          rho_tmp = rho_coord
          u_tmp   = u_coord
        end if

        ! Extrapolate the values of rho and u at the "missing" lattice
        ! site using those in of the neighbor in opposite direction.
        rho_tmp = rho_coord - (rho_tmp - rho_coord)
        u_tmp   = u_coord - ( u_tmp - u_coord )
      end if

      ! Add this coordinate's contribution to the gradients
      ! according to the lattice differential operators
      ! note that rho_tmp = cs^2 * p. This cancels out with the
      ! division by cs^2 in the definition of the lattice differential
      ! operator, which is why cs^2 is not in this term
      grad_p = grad_p + wq * cq * rho_tmp

      ! Compute the cross-product term, then add it to curl_u
      ! Convert the volume-averaged fluid velocity (= what grabValueAtCoord gives)
      ! to the fluid phase velocity for Generalized Navier-Stokes
      u_tmp = u_tmp

      call cross_product( wq*cq, u_tmp, curl_tmp )
      curl_u = curl_u + curl_tmp
    end do ! iDir

    curl_u = curl_u * cs2inv

  end subroutine calcVelocityAndPressureGradient_onewaycoupled_old

  !> Get the value of the density and velocity at a certain integer coordinate
  subroutine grabValueAtCoord_twowaycoupled(coord, scheme, vel_pos, dens_pos, vol_frac_pos, rho, u, eps_f, err )
    integer, intent(in) :: coord(4)
    type(mus_scheme_type), intent(inout) :: scheme
    ! Position of density in varSys
    integer, intent(in) :: dens_pos
    ! Position of velocity in varSys
    integer, intent(in) :: vel_pos(3)
    ! Position of volume fraction in varSys
    integer, intent(in) :: vol_frac_pos
    ! Output: density at coord
    real(kind=rk), intent(out) :: rho
    ! Output: velocity at coord
    real(kind=rk), intent(out) :: u(3)
    ! Output: fluid volume fraction at coord
    real(kind=rk), intent(out) :: eps_f
    !> Error code: set to TRUE if we could not find ldPos of coord
    !! we are trying to grab value at (for example because it is
    !! outside domain and not a halo.
    logical, intent(out) :: err
    ! ------------------------------------------!
    integer(kind=long_k) :: TreeID
    integer :: ldPos, elemOff, lev
    ! ------------------------------------------!
    lev = coord(4)
    err = .FALSE.

    ! get TreeID and position in complete tree
    TreeID = tem_IdOfCoord(coord = coord)

    ldPos = tem_PosOfId( sTreeID = TreeID,                        &
                      & treeIDlist = scheme%levelDesc(lev)%total, &
                      & lower      = 1,                           &
                      & upper      = scheme%pdf(lev)%nElems_fluid )
    if( ldPos <= 0) then
      ! If we could not find ldPos in local elems, look in halos
      ldPos = tem_PosOfId( sTreeID = TreeID,                               &
                         & treeIDlist = scheme%levelDesc(lev)%total,        &
                         & lower      = scheme%pdf(lev)%nElems_fluid + 1,  &
                         & upper      = scheme%pdf(lev)%nElems_local        )

      if (ldPos <= 0) then
        ! If we couldn't find the element on this proc, set err to TRUE
        ! write(logUnit(1),*) 'WARNING rank mus_particles_DPS_grabValueAtCoord: ldPos <= 0'
        ! write(logUnit(1),*) 'ldPos =',  ldPos
        ! write(logUnit(1),*) 'TreeID =', TreeID
        ! write(logUnit(1),*) 'coord = ', coord
        err = .TRUE.
        rho = -1.0_rk
        u = -1.0_rk
        eps_f = -1.0_rk
        return
      end if ! could not find element in halos
    end if

    ! If we were able to find the element on this proc, get
    ! density and velocity values from auxField
    elemOff = (ldPos-1)*scheme%varSys%nAuxScalars

    rho  = scheme%auxField(lev)%val(elemOff + dens_pos)
    u(1) = scheme%auxField(lev)%val(elemOff + vel_pos(1))
    u(2) = scheme%auxField(lev)%val(elemOff + vel_pos(2))
    u(3) = scheme%auxField(lev)%val(elemOff + vel_pos(3))
    eps_f  = scheme%auxField(lev)%val(elemOff + vol_frac_pos)

  end subroutine grabValueAtCoord_twowaycoupled

  !> Get the value of the density and velocity at a certain integer coordinate
  subroutine grabValueAtCoord_onewaycoupled(coord, scheme, vel_pos, dens_pos, rho, u, err )
    integer, intent(in) :: coord(4)
    type(mus_scheme_type), intent(inout) :: scheme
    ! Position of density in varSys
    integer, intent(in) :: dens_pos
    ! Position of velocity in varSys
    integer, intent(in) :: vel_pos(3)
    ! Output: density at coord
    real(kind=rk), intent(out) :: rho
    ! Output: velocity at coord
    real(kind=rk), intent(out) :: u(3)
    !> Error code: set to TRUE if we could not find ldPos of coord
    !! we are trying to grab value at (for example because it is
    !! outside domain and not a halo.
    logical, intent(out) :: err
    ! ------------------------------------------!
    integer(kind=long_k) :: TreeID
    integer :: ldPos, elemOff, lev
    ! ------------------------------------------!
    lev = coord(4)
    err = .FALSE.

    ! get TreeID and position in complete tree
    TreeID = tem_IdOfCoord(coord = coord)

    ldPos = tem_PosOfId( sTreeID = TreeID,                        &
                      & treeIDlist = scheme%levelDesc(lev)%total, &
                      & lower      = 1,                           &
                      & upper      = scheme%pdf(lev)%nElems_fluid )
    if( ldPos <= 0) then
      ! If we could not find ldPos in local elems, look in halos
      ldPos = tem_PosOfId( sTreeID = TreeID,                               &
                         & treeIDlist = scheme%levelDesc(lev)%total,        &
                         & lower      = scheme%pdf(lev)%nElems_fluid + 1,  &
                         & upper      = scheme%pdf(lev)%nElems_local        )

      if (ldPos <= 0) then
        ! If we couldn't find the element on this proc, set err to TRUE
        ! write(logUnit(1),*) 'WARNING rank mus_particles_DPS_grabValueAtCoord: ldPos <= 0'
        ! write(logUnit(1),*) 'ldPos =',  ldPos
        ! write(logUnit(1),*) 'TreeID =', TreeID
        ! write(logUnit(1),*) 'coord = ', coord
        err = .TRUE.
        rho = -1.0_rk
        u = -1.0_rk
        return
      end if ! could not find element in halos
    end if

    ! If we were able to find the element on this proc, get
    ! density and velocity values from auxField
    elemOff = (ldPos-1)*scheme%varSys%nAuxScalars

    rho  = scheme%auxField(lev)%val(elemOff + dens_pos)
    u(1) = scheme%auxField(lev)%val(elemOff + vel_pos(1))
    u(2) = scheme%auxField(lev)%val(elemOff + vel_pos(2))
    u(3) = scheme%auxField(lev)%val(elemOff + vel_pos(3))

  end subroutine grabValueAtCoord_onewaycoupled

  ! ***** Routines for hydrodynamic force computation ***** !

  !> Routine to calculate the drag force according to
  !! [1] S. Tenneti, R. Garg, and S. Subramaniam, “Drag law for monodisperse
  !! gas–solid systems using particle-resolved direct numerical simulation of flow
  !! past fixed assemblies of spheres,” International Journal of Multiphase Flow,
  !! vol. 37, no. 9, pp. 1072–1092, Nov. 2011, doi:
  !! 10.1016/j.ijmultiphaseflow.2011.05.010.

  subroutine applyDragForce_DPS( particle, eps_p, nu, Fd )
    !> Particle to apply force to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Solid volume fraction interpolated to location of the particle
    real(kind=rk), intent(in) :: eps_p
    !> Fluid kinematic viscosity (phy)
    real(kind=rk), intent(in) :: nu
    !> Output: drag force on particle
    real(kind=rk), intent(out) :: Fd(3)
    ! ------------------------------------------!
    real(kind=rk) :: F_stokes(3)
    real(kind=rk) :: umag, Re, Cd0, A, B
    real(kind=rk) :: c_eps   ! c_eps = 1-eps_p
    real(kind=rk) :: u_rel(3)
    ! ------------------------------------------!
    ! 1) Compute Reynolds number
    c_eps = abs(1.0_rk - eps_p)
    u_rel = particle%u_fluid(1:3) - particle%vel(1:3)
    umag = dot_product(u_rel,u_rel)
    umag = sqrt(umag)
    Re = (umag*c_eps*2.0*particle%radius)/nu

    ! 2) Compute drag coefficient

    Cd0 = (1.0_rk + 0.15*Re**0.687)/c_eps**3

    A   = 5.81*eps_p / c_eps**3 + 0.48*( eps_p**(1.0/3.0) / c_eps**4 )
    B   = eps_p**3 * Re * ( 0.95 + 0.61*eps_p**3 / c_eps**2 )

    ! 3) Compute drag force
    F_stokes = 6.0*PI*particle%radius*nu*particle%rho_fluid*u_rel*c_eps
    Fd = (Cd0 + A + B)*F_stokes

  end subroutine applyDragForce_DPS


  subroutine applyDragForce_DPS_noeps( particle, eps_p, nu, Fd )
    !> Particle to apply force to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Solid volume fraction interpolated to location of the particle
    real(kind=rk), intent(in) :: eps_p
    !> Fluid kinematic viscosity (phy)
    real(kind=rk), intent(in) :: nu
    !> Output: drag force on particle
    real(kind=rk), intent(out) :: Fd(3)
    ! ------------------------------------------!
    real(kind=rk) :: umag, Re, Cd0
    real(kind=rk) :: u_rel(3)
    ! ------------------------------------------!
    ! 1) Compute Reynolds number
    u_rel = particle%u_fluid(1:3) - particle%vel(1:3)
    umag = dot_product(u_rel,u_rel)
    umag = sqrt(umag)
    Re = (umag*2.0*particle%radius)/nu

    Cd0 = 1.0_rk + 0.15*Re**0.687

    ! 3) Compute drag force
    Fd = 6.0*PI*particle%radius*nu*particle%rho_fluid*Cd0*u_rel

  end subroutine applyDragForce_DPS_noeps

  !> applyLiftForce_DPS computes the Saffman lift force on a particle according to
  !! the paper :
  !!   Saffman PG. The lift on a small sphere in a slow shear flow.
  !!   J Fluid Mech 1965;22(2):385–400. doi:10.1017/S0 022112065000824.
  subroutine applyLiftForce_DPS( particle, nu, Flift )
    !> Particle to apply force to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Fluid kinematic viscosity (phy)
    real(kind=rk), intent(in) :: nu
    !> Output: drag force on particle
    real(kind=rk), intent(out) :: Flift(3)
    ! ------------------------------------------!
    real(kind=rk) :: u_rel(3)
    real(kind=rk) :: curl_u_mag
    real(kind=rk) :: fac1, fac2(3)
    real(kind=rk) :: tol ! tolerance when curl_u = 0
    ! ------------------------------------------!
    tol = 1.0e-6
    ! Compute relative velocity between particle and fluid
    u_rel = particle%u_fluid(1:3) - particle%vel(1:3)

    ! Compute magnitude of curl of fluid velocity
    curl_u_mag = dot_product( particle%curl_u_fluid, particle%curl_u_fluid )
    curl_u_mag = sqrt(curl_u_mag)
    if( curl_u_mag < tol ) then
      Flift = (/ 0.0_rk, 0.0_rk, 0.0_rk /)
    else
      ! Compute square root term
      fac1 = sqrt( nu * particle%rho_fluid / curl_u_mag )

      ! Compute the cross-product term and store result in fac2
      call cross_product( u_rel, particle%curl_u_fluid, fac2 )

      Flift = 1.61 * (2.0*particle%radius)**2 * fac1 * fac2
    end if

  end subroutine applyLiftForce_DPS

  !> Calculate the pressure gradient force on a particle
  subroutine applyPressureForce_DPS( particle, Fp )
    !> Particle to apply force to
    type(mus_particle_DPS_type), intent(inout) :: particle
    !> Output: pressure gradient force
    real(kind=rk), intent(out) :: Fp(3)
    ! ------------------------------------------!
    real(kind=rk) :: Vp      ! particle volume
    ! ------------------------------------------!
    Vp = (4.0_rk/3.0_rk) * PI * particle%radius**3
    Fp = -Vp * particle%grad_p_fluid(1:3)
  end subroutine applyPressureForce_DPS

  ! ************************************************************************** !
  !> Modify the auxField of an element with posInTotal in the total list
  subroutine modify_AuxField_of_Elem(posInTotal, scheme, auxField, iLevel, &
    & varSys, derVarPos, G_lat)
    ! ------------------------------------------------------------------------ !
    !> Position of element to modify in total list
    integer, intent(in) :: posInTotal
    !> Scheme
    type(mus_scheme_type), intent(inout) :: scheme
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Force to apply, in lattice units
    real(kind=rk), intent(in) :: G_lat(3)
    ! ------------------------------------------------------------------------ !
    integer :: vel_pos(3), elemOff
    real(kind=rk) :: forceTerm(3)
    ! ------------------------------------------------------------------------ !
    ! position of velocity field in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    elemoff = (posInTotal-1)*varSys%nAuxScalars
    ! forceterm to add to velocity: F/(2*rho)
    forceTerm = G_lat(1:3)*0.5_rk*rho0Inv

    ! add force to velocity
    auxField(elemOff+vel_pos(1)) = auxField(elemOff+vel_pos(1)) + forceTerm(1)
    auxField(elemOff+vel_pos(2)) = auxField(elemOff+vel_pos(2)) + forceTerm(2)
    auxField(elemOff+vel_pos(3)) = auxField(elemOff+vel_pos(3)) + forceTerm(3)

  end subroutine modify_AuxField_of_Elem

  !> Add a momentum increment to element with PosInTotal in total list
  subroutine add_Momentum_Increment_to_Auxfield(posInTotal, scheme, auxField,  &
    & iLevel, varSys, derVarPos, momentumIncrement)
    ! ------------------------------------------------------------------------ !
    !> Position of element to modify in total list
    integer, intent(in) :: posInTotal
    !> Scheme
    type(mus_scheme_type), intent(inout) :: scheme
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Momentum increment to apply, in lattice units
    real(kind=rk), intent(in) :: momentumIncrement(3)
    ! ------------------------------------------------------------------------ !
    integer :: vel_pos(3), elemOff
    real(kind=rk) :: velocityIncrement(3)
    ! ------------------------------------------------------------------------ !
    ! position of velocity field in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    elemoff = (posInTotal-1)*varSys%nAuxScalars

    velocityIncrement = momentumIncrement(1:3)*rho0Inv

    ! add force to velocity
    auxField(elemOff+vel_pos(1)) = auxField(elemOff+vel_pos(1)) + velocityIncrement(1)
    auxField(elemOff+vel_pos(2)) = auxField(elemOff+vel_pos(2)) + velocityIncrement(2)
    auxField(elemOff+vel_pos(3)) = auxField(elemOff+vel_pos(3)) + velocityIncrement(3)

  end subroutine add_Momentum_Increment_to_Auxfield

  subroutine applySrc_toElem( G, posInTotal, scheme, iLevel )
    ! -------------------------------------------------------------------- !
    !> Body force to apply to cell, in lattice units
    real(kind=rk), intent(in) :: G(3)

    !> Position of this element in the total list
    integer, intent(in) :: posInTotal

    !> scheme type
    type(mus_scheme_type), intent(inout) :: scheme

    !> current level
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: velocity(3), ucx, uMinusCX(3), forceTerm
    integer :: nPdfSize, iDir, QQ, nScalars, statePos, elemOff, next
    integer :: vel_pos(3)
    real(kind=rk) :: omega, omega_fac
    ! ---------------------------------------------------------------------------
    QQ = scheme%layout%fStencil%QQ
    nScalars = scheme%varSys%nScalars
    nPdfSize = scheme%pdf(iLevel)%nSize
    next = scheme%pdf(iLevel)%nNext

    ! Position of velocity variable in auxField
    vel_pos = scheme%varSys%method%val( &
    & scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    ! element offset
    elemoff = (posInTotal-1)*scheme%varSys%nAuxScalars

    ! obtain velocity from auxField
    velocity(1) = scheme%auxField(iLevel)%val(elemOff + vel_pos(1))
    velocity(2) = scheme%auxField(iLevel)%val(elemOff + vel_pos(2))
    velocity(3) = scheme%auxField(iLevel)%val(elemOff + vel_pos(3))

    ! get the correct omega value
    omega = scheme%field(1)%fieldProp%fluid%viscKine              &
      &                              %omLvl(iLevel)%val(posInTotal)
    omega_fac = 1.0_rk - omega * 0.5_rk

    ! force term:
    ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
    !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}

    ! write(logUnit(1),*) "posInTotal = ", posInTotal
    do iDir = 1, QQ
      ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
        &                velocity )
      uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - velocity

      forceTerm = dot_product( uMinusCx * cs2inv               &
        &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
        &       * cs4inv, G )

      ! Get position in state array
      statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq

      ! Update outstate
      scheme%state(iLevel)%val(statePos, next) = scheme%state(iLevel)%val(statePos, next)                         &
        & + omega_fac * scheme%layout%weight( iDir ) * forceTerm
    end do
  end subroutine applySrc_toElem

  subroutine applySrc_toElem_eps( G, posInTotal, scheme, iLevel )
    ! -------------------------------------------------------------------- !
    !> Body force to apply to cell, in lattice units
    real(kind=rk), intent(in) :: G(3)

    !> Position of this element in the total list
    integer, intent(in) :: posInTotal

    !> scheme type
    type(mus_scheme_type), intent(inout) :: scheme

    !> current level
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: velocity(3), inv_vol_frac, ucx, uMinusCX(3), forceTerm
    integer :: nPdfSize, iDir, QQ, nScalars, statePos, elemOff, next
    integer :: vel_pos(3), vol_frac_pos
    real(kind=rk) :: omega, omega_fac
    ! ---------------------------------------------------------------------------
    QQ = scheme%layout%fStencil%QQ
    nScalars = scheme%varSys%nScalars
    nPdfSize = scheme%pdf(iLevel)%nSize
    next = scheme%pdf(iLevel)%nNext

    ! Position of velocity variable in auxField
    vel_pos = scheme%varSys%method%val( &
    & scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)

    vol_frac_pos = scheme%varSys%method%val( &
    & scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)

    ! element offset
    elemoff = (posInTotal-1)*scheme%varSys%nAuxScalars

    ! obtain velocity from auxField
    velocity(1) = scheme%auxField(iLevel)%val(elemOff + vel_pos(1))
    velocity(2) = scheme%auxField(iLevel)%val(elemOff + vel_pos(2))
    velocity(3) = scheme%auxField(iLevel)%val(elemOff + vel_pos(3))
    inv_vol_frac = 1.0_rk / scheme%auxField(iLevel)%val(elemOff + vol_frac_pos)

    ! get the correct omega value
    omega = scheme%field(1)%fieldProp%fluid%viscKine              &
      &                              %omLvl(iLevel)%val(posInTotal)
    omega_fac = 1.0_rk - omega * 0.5_rk

    ! force term:
    ! F_i = w_i( (\vec{e}_i-\vec{u}*)/cs2 +
    !       (\vec{e}_i \cdot \vec{u}*)\vec{e}_i/cs4) \cdot \vec{F}

    ! write(logUnit(1),*) "posInTotal = ", posInTotal
    do iDir = 1, QQ
      ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), &
        &                velocity )
      uMinusCx = scheme%layout%fStencil%cxDirRK(:, iDir) - velocity*inv_vol_frac

      forceTerm = dot_product( uMinusCx * cs2inv               &
        &       + ucx * scheme%layout%fStencil%cxDirRK(:,iDir) &
        &       * cs4inv * inv_vol_frac, G )

      ! Get position in state array
      statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq

      ! Update outstate
      scheme%state(iLevel)%val(statePos, next) &
        &    = scheme%state(iLevel)%val(statePos, next) &
        &    + omega_fac * scheme%layout%weight( iDir ) * forceTerm
    end do
  end subroutine applySrc_toElem_eps

end module mus_particle_DPS_module
