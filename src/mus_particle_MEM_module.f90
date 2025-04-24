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
!> mus_particle_MEM module contains routines needed for the simulation of fully
!! resolved Momentum Exchange Method particles

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

module mus_particle_MEM_module

  use env_module,              only: rk, long_k
  use tem_param_module,        only: rho0_lat => rho0, cs2, cs2inv
  use tem_logging_module,      only: logUnit
  use tem_aux_module,          only: tem_abort
  use tem_geometry_module,     only: tem_CoordOfReal, tem_PosOfId
  use tem_topology_module,     only: tem_IdOfCoord, tem_FirstIdAtLevel
  use tem_construction_module, only: tem_levelDesc_type
  use tem_varSys_module,       only: tem_varSys_type
  use tem_varMap_module,       only: tem_varMap_type
  use tem_stencil_module,      only: tem_stencilHeader_type,  &
    &                                tem_stencil_findIndexOfDir
  use tem_property_module,     only: prp_solid, prp_particle, &
    &                                prp_hasBnd, prp_sendHalo
  use tem_dyn_array_module,    only: init, append, destroy,             &
    &                                empty, dyn_intArray_type,          &
    &                                dyn_longArray_type, PosOfVal_long, &
    &                                SortPosOfVal_long
  use tem_grow_array_module,   only: init, append, destroy, empty, &
    &                                grw_int2darray_type,          &
    &                                grw_logical2darray_type,      &
    &                                grw_intarray_type,            &
    &                                grw_longarray_type,           &
    &                                grw_real2darray_type

  use mus_geom_module,        only: mus_geom_type
  use mus_scheme_type_module, only: mus_scheme_type
  use mus_param_module,       only: mus_param_type

  use mus_particle_type_module, only: mus_particle_MEM_type,   &
    &                                 mus_particle_group_type, &
    &                                 printParticleGroup,      &
    &                                 printParticleGroup2_MEM, &
    &                                 printpIDlist,            &
    &                                 remove_particle_from_da_particle_MEM
  use mus_particle_logging_module, only: logParticleData, getParticleLogUnit, &
    &                                    closeParticleLog
  use mus_particle_logging_type_module, only: mus_particle_logging_type, &
    &                                         pgDebugLog
  use mus_particle_comm_type_module, only: mus_particles_communication_type
  use mus_particle_aux_module, only: findPartitionOfTreeID, &
    &                                getBaryOfCoord,        &
    &                                cross_product
  use mus_particle_boundary_module,  only: pgBndData, wrapPeriodicPos, &
    &                                      getNeighborCoord,           &
    &                                      calcPeriodicRsurface,       &
    &                                      computeDisplacement

  implicit none
  private

  public :: initParticle_MEM
  public :: applyHydrodynamicForces
  public :: mapToLattice
  public :: updateCoordOfOrigin
  public :: updateExclusionList
  public :: updateSolidNodes
  public :: updateFluidNeighbors
  public :: updateNewFluidNodes
  public :: destroyParticle_MEM
  public :: applyVelocityBounceback

  public :: checkForParticleOverlap
  public :: setToEquilibrium
  public :: make_pdf_tiny


contains


  !> Routine to initialize a newly created MEM particle
  subroutine initParticle_MEM( particle, particleID, geometry, scheme, myRank, comm, rmflag  )
    !> Particle to initialize
    type(mus_particle_MEM_type), intent(inout) :: particle
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
    !> flag to indicate whether particle exists on this process or should be removed
    logical, intent(out) :: rmflag
    ! -----------------------------------------------!
    real(kind=rk) :: dx
    integer :: lev                  ! level to map particle on
    ! -----------------------------------------------!

    lev = geometry%tree%global%maxLevel
    dx = geometry%tree%global%BoundingCubeLength / 2**lev

    ! particle radius in lattice units, rounded up
    particle%Rn = ceiling(particle%radius / dx)

    ! ---- INITIALIZE PARTICLE REPRESENTATION ON GRID ---- !
    call updateExclusionList( this     = particle,    &
                            & scheme   = scheme,      &
                            & geometry = geometry,    &
                            & myRank   = myRank   ,   &
                            & procs    = comm%proc,   &
                            & nProcs   = comm%nProcs, &
                            & dx       = dx,          &
                            & rmflag   = rmflag       )

    ! Set old particle owner to same value as current owner initially
    if( particle%previousOwner < 0 ) then
      particle%previousOwner = particle%owner
    end if

    ! Make PDF values of particle tiny for visualization
    call make_pdf_tiny(scheme, particle, lev)

    ! intialize exclusionListBuffer to same size as exclusionList
    call init( me     = particle%exclusionListBuffer,        &
      &        length = particle%exclusionList%containersize )

    ! initialize the makeFluidList and newSolidList
    call init( me     = particle%makeFluidList,              &
      &        length = particle%exclusionList%containersize )

    ! set newForMe to false to indicate particle has been initialized
    particle%newForMe = .FALSE.


    ! Modify the connectivity (scheme%pdf%neigh) list for elements belonging
    ! to particle (exclusionList)
    ! * sets the streaming directions for all elements in exclusionList
    ! * sets the prp_solid bits

    call updateSolidNodes( this = particle,                  &
                         & scheme = scheme,                  &
                         & stencil = scheme%layout%fStencil  )

    ! Modify the connectivity (scheme%pdf%neigh) list for fluid elements
    ! adjacent to particle
    ! * loops over all particle elements in exclusionList
    ! * for each element it checks each lattice direction for neighboring
    !   fluid elems (w/o prp_solid)
    ! * then sets bounceback for those fluid elems

    call updateFluidNeighbors( this = particle,                   &
                             & scheme = scheme,                   &
                             & stencil = scheme%layout%fStencil   )

  end subroutine initParticle_MEM

  !> ApplyHydrodynamicForces performs the momentum transfer
  !! from fluid TO particle
  subroutine applyHydrodynamicForces( this, scheme, stencil, params, rho_p_lat )
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> fluid stencil
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> Parameters for access to conversion factors between units
    type(mus_param_type), intent(in) :: params
    !> Density value to compute equilibrium PDFs with in case local fluid density
    !! is not available because i.e. the adjacent element is a wall or
    !! another particle.
    real(kind=rk), intent(in), optional :: rho_p_lat
    ! ------------------------------------------!

    integer :: lev
    integer :: elemPos, neighPos, fluidPos
    ! integer :: elemPos, neighPos, fluidPos
    integer :: iDir, fluidDir, iField, idx_idir
    integer(kind=long_k) :: neighProp
    integer :: iElem, nLocalElems, nSize
    integer :: stateVarPos(stencil%QQ)
    integer :: nNow, nNext
    real(kind=rk) :: dt, dx, inv_dx, inv_dt, inv_vel, pdf_val
    real(kind=rk) :: wq, MoBoFactor, movingBoundaryTerm
    real(kind=rk) :: u_surf_lat(3)
    real(kind=rk) :: baryOfOrigin(3), baryOfSurface(3), cq(3)

    ! position vector particle origin to surface elem, angular velocity,
    ! force and torque per link. All in lattice units
    real(kind=rk) :: r(3), r_lat(3), om_lat(3), F_q_lat(3),  dM_lat(3)

    ! For computing the equilibrium distribution in case our
    ! particle borders a solid element
    logical :: linkBordersSolid
    real(kind=rk) :: feq( scheme%layout%fStencil%QQ )
    real(kind=rk) :: eq_rho( 1 )
    real(kind=rk) :: eq_vel( 3, 1 )

    ! ------------------------------------------!
    lev = this%coordOfOrigin(4)
    iField = 1

    dt = params%physics%dtLvl(lev)
    inv_dt = 1.0 / dt
    dx = params%physics%dxLvl(lev)
    inv_dx = 1.0_rk / dx
    inv_vel = 1.0_rk / params%physics%fac(lev)%vel
    om_lat = this%vel(4:6) * dt
    MoBoFactor = -2.0 * cs2inv * rho0_lat ! prefactor for moving boundary term

    baryOfOrigin = this%pos(1:3)

    nLocalElems = scheme%pdf( lev )%nElems_fluid
    nSize  = scheme%pdf( lev )%nSize
    nNow   = scheme%pdf( lev )%nNow
    nNext  = scheme%pdf( lev )%nNext

    stateVarPos(:stencil%QQ) = scheme%varSys                    &
      &                              %method                    &
      &                              %val(scheme%stateVarMap    &
      &                              %varPos                    &
      &                              %val(1))                   &
      &                              %state_varPos( :stencil%QQ )

    ! Initialize current force and torque in buffer to 0
    this%Fbuff(this%Fnow,1:6) = (/ 0.0_rk, 0.0_rk, 0.0_rk, &
       &                          0.0_rk, 0.0_rk, 0.0_rk /)

    ! Loop over the elements belonging to this particle (those in exclusionList)
    ! Then identify fluid neighbors from there.
    particleElemLoop: do iElem = 1, this%exclusionList%nVals
      ! position of this node in levelDesc state vector
      elemPos = this%exclusionList%val(iElem)

      ! Only compute force contribution from links connecting local fluid elems
      ! not halos!
      if(elemPos > nLocalElems ) then
        cycle
      end if

      baryOfSurface = scheme%levelDesc(lev)%baryOfTotal(elemPos,1:3)



      linkLoop: do iDir = 1, stencil%QQN
        ! Get the index of the neighbor element in this direction.
        neighPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, elemPos)

        ! direction FROM fluid TO particle = cq in Bartuschat eqn 5.2
        fluidDir = scheme%layout%fStencil%cxDirInv(stateVarPos(iDir))

        ! get vector for stencil direction fluidDir
        ! NB: this is in lattice units!
        cq = stencil%cxDirRK(:,fluidDir)

        ! compute vector from particle origin to current point on the surface
        ! subtract 0.5*cq as surface is halfway between solid node and fluid node
        r = baryOfSurface - baryOfOrigin
        call calcPeriodicRsurface( r            = r,           &
                                 & R_particle   = this%radius, &
                                 & boundaryData = pgBndData    )
        r_lat = r * inv_dx - 0.5*cq

        ! compute rotational component of surface velocity using cross product
        call cross_product( om_lat, r_lat, u_surf_lat )

        ! add translational component
        u_surf_lat = u_surf_lat + this%vel(1:3) * inv_vel

        ! Stencil weight in this direction
        wq = scheme%layout%weight(fluidDir)

        ! Initialize neighbor's property
        neighProp = 0_long_k

        ! If neighbor element exists, get its property
        if( neighPos > 0 ) then
          neighProp = scheme%levelDesc(lev)%property(neighPos)
            if ( btest(neighProp, prp_solid) ) then
              ! Check if this solid element is part of our particle or a different one
              if( ( any( neighPos &
              &      == this%exclusionList%val(1:this%exclusionList%nVals)  ) )  ) then
                ! Solid element is part of this particle, skip to next link.
                cycle
              else
                linkBordersSolid = .TRUE.
              end if ! element is part of this particle
            else
              ! No prp_solid means link points to a fluid neighbor
              linkBordersSolid = .FALSE.
            end if ! prp_solid
        else
          ! neighPos <= 0 means we have a boundary in this direction and neighPos = -BCid
          linkBordersSolid = .TRUE.
        end if ! neighPos > 0

        if( linkBordersSolid ) then
          ! fill in the missing links for the force computation using the equilibrium distribution
          ! calculate eq dist (note: rho and velocity must be in lattice units)
          if(present(rho_p_lat)) then
            eq_rho(1) = rho_p_lat
          else
            eq_rho(1) = rho0_lat
          end if
          eq_vel(1:3,1) = u_surf_lat(1:3)

          call scheme%derVarPos(1)%equilFromMacro( density  = eq_rho,         &
                                                 & velocity = eq_vel,         &
                                                 & iField   = 1,              &
                                                 & nElems   = 1,              &
                                                 & varSys   = scheme%varSys,  &
                                                 & layout   = scheme%layout,  &
                                                 & res      = fEq             )

          ! Use the equilibrium distribution as pdf value for momentum exchange
          pdf_val = fEq(fluidDir)

        else
          ! this link points to a fluid neighbor
          fluidPos = neighPos

          idx_idir = scheme%varSys%method%val(                      &
            &               scheme%stateVarMap%varPos%val(iField) ) &
            &                     %state_varPos(fluidDir)

          ! Get the pdf value in the direction from fluid to particle
          pdf_val =  scheme%state(lev)%val(                                    &
            &   ( fluidpos-1)* scheme%varsys%nscalars+idx_idir, nNext)
        end if ! link borders solid

        ! compute momentum exchange term due to moving boundary
        movingBoundaryTerm = MoBoFactor * wq  * dot_product(cq, u_surf_lat)

        ! Now transfer momentum FROM fluid TO particle:
        ! Compute contribution to total force by this link and this particle
        F_q_lat = (2.0 * pdf_val + movingBoundaryTerm) * cq

        ! Add this link's contribution to the particle force and torque
        ! in PHYSICAL units
        this%Fbuff( this%Fnow, 1:3 ) = this%Fbuff( this%Fnow, 1:3 ) &
          & + F_q_lat * params%physics%fac(lev)%force

        call cross_product(r_lat, F_q_lat, dM_lat)

        this%Fbuff( this%Fnow, 4:6 ) = this%Fbuff( this%Fnow, 4:6) &
          & + dM_lat * params%physics%fac(lev)%force * dx

      end do linkLoop
    end do particleElemLoop

  end subroutine applyHydrodynamicForces


  ! -------- ROUTINES FOR MAPPING PARTICLE TO LATTICE -------!
  !> mapToLattice is the high-level routine performs a full re-mapping each time step
  subroutine mapToLattice( this, particleGroup, scheme, stencil, &
      &                     geometry, params, rmflag              )
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Particle group in which to search for collisions
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> fluid stencil
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> Geometry to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> Parameters for access to e.g. time step dt
    type(mus_param_type), intent(in) :: params
    !> Flag that tells whether to remove this particle from this process
    logical, intent(out) :: rmflag

    ! ------------------------------------------!
    ! property bits for neighbor element
    integer :: elemPos
    integer :: iElem
    integer :: lev

    logical :: overlapsOtherParticle

    ! macroscopic quantities for intializing new fluid elems to eq PDF
    real(kind=rk), allocatable :: rho_lat(:)            ! size (Nelems)
    real(kind=rk), allocatable :: vel_surf_lat(:,:)     ! size (3,Nelems)

    real(kind=rk) :: dt, dx
    real(kind=rk) :: inv_dx, inv_dt, inv_vel
    real(kind=rk) :: baryOfOrigin(3)
    real(kind=rk) :: baryOfSurface(3)

    ! vector from particle origin to surface element, lattice units
    real(kind=rk) :: r_lat(3)

    ! particle angular velocity vector in lattice units
    real(kind=rk) :: om_lat(3)

    ! ------------------------------------------!
    ! 0. Initialize parameters from scheme, stencil, etc.
    ! write(logUnit(1),'(A)') 'Enter mapToLattice'

    ! Set rmflag to false by default
    rmflag = .FALSE.

    lev = this%coordOfOrigin(4)
    dt = params%physics%dtLvl(lev)
    dx = params%physics%dxLvl(lev)
    inv_dx = 1.0_rk / dx
    inv_dt = 1.0_rk / dt
    inv_vel = 1.0_rk / params%physics%fac(lev)%vel

    ! --- Update exclusionList using neighbor information --- !
    ! copy old exclusionList into exclusionListBuffer for creating makeFluidList later
    do iElem = 1, this%exclusionList%nVals
      call append( me  = this%exclusionListBuffer,     &
                 & val = this%exclusionList%val(iElem) )
    end do

    ! Then empty the exclusionList and re-fill with new values
    ! This also sets rmflag to TRUE if the particle "bounding box"
    ! No longer overlaps any elements of this process.
    call updateExclusionList( this     = this,                      &
                            & scheme   = scheme,                    &
                            & geometry = geometry,                  &
                            & myRank   = params%general%proc%rank,  &
                            & procs    = particleGroup%send%proc,   &
                            & nProcs   = particleGroup%send%nProcs, &
                            & dx       = dx,                        &
                            & rmflag   = rmflag                     )

    ! ----- CHECK IF RE-MAPPING LEADS TO OVERLAP ----- !
    ! If so:
    ! * Reset all element positions to the ones prior to moving particle
    ! * Then exit mapToLattice function

    ! Check for overlap with discrete representations of other particles
    ! If overlap is detected, overlapsOtherParticle will be set to .TRUE.
    overlapsOtherParticle = .FALSE.
    call checkForParticleOverlap(this, particleGroup, scheme, overlapsOtherParticle)

    if(overlapsOtherParticle) then

      open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
        write(pgDebugLog%lu,'(A)') 'OverlapsOtherParticle = true at t = '
        write(pgDebugLog%lu, *) params%general%simcontrol%now%sim
      close( pgDebugLog%lu )


      call empty( me = this%exclusionList )

      do iElem = 1, this%exclusionListBuffer%nVals
        call append( me  = this%exclusionList,                 &
                   & val = this%exclusionListBuffer%val(iElem) )
      end do

      ! Empty the exclusionListBuffer for next time step and exit mapToLattice
      call empty( me = this%exclusionListBuffer )
      return

    end if ! overlapsOtherParticle

    ! --- NO OVERLAP --- !
    ! Update the  solid node connectivity in kernel list and set prp_solid bit
    call updateSolidNodes(this, scheme, scheme%layout%fStencil)

    ! Set pdf values of solid nodes to tiny for visualization
    call make_pdf_tiny(scheme, this, lev)

    ! ------- INITIALIZE NEW FLUID ELEMENTS -----------!
    ! * initialize the new fluid elements (former particle elements) to fEq
    !   * find out what the new fluid elements are
    !   * clear their prp_solid bits
    !   * set their pdf to equilibrium distribution fEq(rho, us)
    !     where us is the particle surface velocity
    !   * restore their connectivity to neighboring fluid elements

    call empty( me = this%makeFluidList )

    ! * MakeFluidList will contain all elements which were in EL before
    !   but are no longer after updateExclusionList
    do iElem = 1, this%exclusionListBuffer%nVals
      if( .not.( any( this%exclusionListBuffer%val(iElem) &
        &      == this%exclusionList%val(1:this%exclusionList%nVals)  ) )  ) then

        call append( me   = this%makeFluidList,                 &
                   &  val = this%exclusionListBuffer%val(iElem) )
      end if
    end do

    baryOfOrigin = this%pos(1:3)

    ! Rotational velocity of particle
    om_lat = this%vel(4:6) * dt

    allocate( rho_lat( this%makeFluidList%nVals ) )
    allocate( vel_surf_lat( 3, this%makeFluidList%nVals ) )

    ! Compute surface velocity values at each new fluid node
    do iElem = 1,this%makeFluidList%nVals
      ! get element position in total list
      elemPos = this%makeFluidList%val(iElem)


      ! --- STEP 1: compute density and surface velocity to use --- !
      ! ---         in feq calculation for new fluid nodes      --- !
      baryOfSurface = scheme%levelDesc(lev)%baryOfTotal(elemPos,1:3)
      r_lat = computeDisplacement( baryOfOrigin, baryOfSurface, pgBndData )
      r_lat = r_lat * inv_dx

      ! Compute surface velocity and store in vel_surf_lat
      call cross_product(om_lat, r_lat, vel_surf_lat(:,iElem))

      vel_surf_lat(:,iElem) = vel_surf_lat(:,iElem) + this%vel(1:3) * inv_vel

      ! Compute average density from neighboring fluid nodes
      rho_lat(iElem) = particleGroup%rho0_lat
    end do

    ! Clear the prp_solid bit for all elements
    do iElem = 1,this%makeFluidList%nVals
      ! get element position in total list
        elemPos = this%makeFluidList%val(iElem)

      scheme%levelDesc(lev)%property( elemPos ) = &
          &  ibclr( scheme%levelDesc(lev)%property( elemPos ), prp_solid )
    end do

    ! set PDF of new fluid elems to fEq value using surface velocity values
    call setToEquilibrium( elemList = this%makeFluidList                      &
                         &                %val(1:this%makeFluidList%nVals),   &
                         &  N       = this%makeFluidList%nVals,               &
                         &  lev     = lev,                                    &
                         &  scheme  = scheme,                                 &
                         &  rho     = rho_lat,                                &
                         &  vel     = vel_surf_lat                            )

    deallocate(rho_lat)
    deallocate(vel_surf_lat)

    ! restore connectivity of new fluid elements to neighbor elems and vice-versa
    call updateNewFluidNodes(this, scheme, scheme%layout%fStencil)

    ! -------------------------------------------------!
    ! * modify connectivity of elements adjacent to the updated particle elems
    !   so that they perform bounceback at the particle walls

    call updateFluidNeighbors(this, scheme, scheme%layout%fStencil)
    ! -------------------------------------------------!
    ! clean up: empty the exclusionListBuffer
    call empty( me = this%exclusionListBuffer )


  end subroutine mapToLattice

  !> updateCoordOfOrigin updates the integer coordinate of
  !! the origin of a particle
  subroutine updateCoordOfOrigin( this, geometry )
    !> Particle to update coordOfOrigin of
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Geometry information to determine new coord of origin
    type(mus_geom_type), intent(in) :: geometry
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

    if( any( abs(moveDir) > 1 ) ) then
      write(logUnit(1), *) 'ERROR: particle moved more than one lattice site'
      write(logUnit(1), '(A)',advance='no') 'Particle ID '
      write(logUnit(1), *) this%particleID
      write(logUnit(1), '(A)',advance='no') 'pos'
      write(logUnit(1), *) this%pos(1:3)
      write(logUnit(1), '(A)',advance='no') 'vel'
      write(logUnit(1), *) this%vel(1:3)
      write(logUnit(1), '(A)',advance='no') 'F'
      write(logUnit(1), *) this%F(1:3)
      write(logUnit(1), '(A)',advance='no') 'MoveDir'
      write(logUnit(1), *) moveDir
      write(logUnit(1), '(A)',advance='no') 'newCoordOfOrigin'
      write(logUnit(1), *) newCoordOfOrigin
      write(logUnit(1), '(A)',advance='no') 'oldCoordOfOrigin'
      write(logUnit(1), *) this%coordOfOrigin
    !  call tem_abort()
    end if

    ! Update coordinate of origin
    this%oldCoordOfOrigin(:) = this%coordOfOrigin(:)
    this%coordOfOrigin(:) = newCoordOfOrigin(:)

  end subroutine updateCoordOfOrigin

  ! UpdateParticleOwner
  ! * updates the MPI process that owns this particle by looking to
  !   find the process which contains the particle center of mass in its local fluid cells.
  ! * Stores the previous owner of the particle
  ! In case the owner of the particle is myRank, then the routine checks if the
  ! TreeID of the particle coordOfOrigin exists as a local fluid element on this process.
  ! If not, it throws an error.
  subroutine updateParticleOwner(this, scheme, geometry, myRank, procs, nProcs)
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry information to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> This process's rank
    integer, intent(in) :: myRank
    !> Array of neighbor procs on which particle could possibly exist
    integer, intent(in) :: procs(:)
    !> Number of procs in array procs
    integer, intent(in) :: nProcs
    ! ----------------------------------------- !
    integer(kind=long_k) :: TreeID
    integer :: ielemproc, lev, ldPos
    ! ----------------------------------------- !
    lev = geometry%tree%global%maxLevel
    TreeID = tem_IdOfCoord(coord = this%coordOfOrigin)

    ! Store the previous owner of the particle
    this%previousOwner = this%owner

    ! First try to look only in my neighbor procs. Most likely we'll find
    ! the particle owner there. If owner is NOT found on neighbor procs,
    ! this%owner will be set to -1 to indicate this
    call findPartitionOfTreeID( sTreeID   = TreeID,     &
                              & geometry  = geometry,   &
                              & procs     = procs,      &
                              & nProcs    = nProcs,     &
                              & outProc   = this%owner, &
                              & procIndex = ielemProc   )

    ! If we cannot find the owner in our neighbor procs, look in
    ! list of all procs
    if( this%owner < 0 ) then
      call findPartitionOfTreeID( sTreeID   = TreeID,     &
                                & geometry  = geometry,   &
                                & nProcs    = geometry%tree%global%nParts,    &
                                & outProc   = this%owner, &
                                & procIndex = ielemProc   )
    end if

    if( this%owner < 0 ) then
      write(logUnit(1),*) "Error: proc ", myRank, &
      & "could not find owner of particle with ID", this%particleID
      call tem_abort()
    end if

    ! At this point we should have found the owner of the particle.
    ! If I am the owner, check to make sure that the TreeID of
    ! the particle origin is actually a fluid element

    if( this%owner == myRank ) then
      ldPos = tem_PosOfId( sTreeID   = TreeID,                      &
                        & treeIDlist = scheme%levelDesc(lev)%total, &
                        & lower      = 1,                           &
                        & upper      = scheme%pdf(lev)%nElems_fluid )
      if(ldPos <= 0) then
        write(logUnit(1),*) "Error: proc ", myRank, &
        & ": origin of particle with ID", this%particleID, "is not on a lattice site!"
        write(logUnit(1),*) "Origin = ", this%pos(1:3)

        call tem_abort()
      end if
    end if

  end subroutine updateParticleOwner

  !> UpdateExclusionList uses the continuous particle position and particle
  !! radius to determine which lattice elements belong to the particle.
  !! These are stored in the 'ExclusionList' of elements which do not participate
  !! in stream-and-collide. It also checks if particles exist or should exist
  !! on some other process and stores this information in the logical masks of
  !! the particle.
  subroutine updateExclusionList(this, scheme, geometry, myRank, procs, nProcs, dx, rmflag)
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry information to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> This process's rank
    integer, intent(in) :: myRank
    !> Array of neighbor procs on which particle could possibly exist
    integer, intent(in) :: procs(:)
    !> Number of procs in array procs
    integer, intent(in) :: nProcs
    !> Mesh size
    real(kind=rk), intent(in) :: dx
    !> Logical to indicate whether particle should be removed from this process
    logical, intent(out) :: rmflag

    ! ------------------------------------------- !
    real(kind=rk) :: x, y, z
    real(kind=rk) :: bary(3)
    real(kind=rk) :: rbary(3)          ! Vector from barycenter of coordOfOrigin to particle position
    integer :: currentCoord(4)
    integer(kind=long_k) :: TreeID
    integer :: ldPos
    integer :: nx, ny, nz, lev

    ! These variables are used to check elements added to dynamic arrays
    integer :: newPos
    logical :: wasAdded
    ! ------------------------------------------ !
    lev = geometry%tree%global%maxLevel

    ! Set previous owner before updating owner in this time step
    this%previousOwner = this%owner

    call updateParticleOwner( this = this, &
                            & scheme = scheme, &
                            & geometry = geometry, &
                            & myRank = myRank, &
                            & procs = procs, &
                            & nProcs = nProcs )

    ! -- UPDATE THE EXCLUSIONLIST -- !
    ! Empty the current exclusionList
    call empty( me = this%exclusionList )

    bary = getBaryOfCoord( coord  = this%coordOfOrigin,          &
                         & origin = geometry%tree%global%origin, &
                         & dx     = dx                           )

    rbary = this%pos(1:3) - bary(1:3)

    ! Identify nodes in neighborhood of origin belonging to particle
    do nx = -this%Rn, this%Rn
      do ny = -this%Rn, this%Rn
        do nz = -this%Rn, this%Rn

          ! integer coordinate of this element
          currentCoord = getNeighborCoord( this%coordOfOrigin, nx, ny, nz, pgBndData )

          ! get TreeID and position in complete tree
          TreeID = tem_IdOfCoord(coord = currentCoord)

          ! Look for the position of this TreeID in levelDescriptor total list
          ! This has to be done in two stages because tem_PosOfId only works for
          ! sorted lists. The total lists is a concatenation of two sorted lists,
          ! one for local fluids and one for halos but the complete list is NOT sorted.
          ldPos = tem_PosOfId( sTreeID = TreeID,                         &
                             & treeIDlist = scheme%levelDesc(lev)%total, &
                             & lower      = 1,                           &
                             & upper      = scheme%pdf(lev)%nElems_fluid )
          if( ldPos <= 0) then
            ! Can't find in local fluids so look in halos
            ldPos = tem_PosOfId( sTreeID    = TreeID,                           &
                               & treeIDlist = scheme%levelDesc(lev)%total,      &
                               & lower      = scheme%pdf(lev)%nElems_fluid + 1, &
                               & upper      = scheme%pdf(lev)%nElems_local      )
          end if ! ldPos <= 0

          if( ldPos > 0 ) then
            ! Compute cartesian coordinates relative to particle origin
            x = nx*dx - rbary(1)
            y = ny*dx - rbary(2)
            z = nz*dx - rbary(3)

            if( x**2 + y**2 + z**2 < this%radius**2 ) then
              ! If element belongs to this particle's local or halo elems,
              ! add it to the exclusionList
              call append( me       = this%exclusionList, &
                        & val      = ldPos,              &
                        & pos      = newPos,             &
                        & wasAdded = wasAdded            )
              if( .NOT. wasAdded ) then
                write(logUnit(1),'(A)') &
                & 'ERROR updateExclusionList: could not add element'
              end if
            end if ! belongs to particle
          end if
        end do
      end do
    end do

    ! Update the existsOnProc mask and check if this particle should be removed from
    ! this group (in that case rmflag is set to .TRUE.)
    call updateExistsOnProc( this, scheme, geometry, myRank, procs, nProcs, rmflag )

  end subroutine updateExclusionList


  !> updateExistsOnProc updates the mask that tells us whether a particle exists
  !! on our neighboring procs. It does this by checking which processes are intersected
  !! by the bounding box of the particle + a "bumper" of one lattice cell to account
  !! for the fact that a particle can also exist on a proc as purely halos
  subroutine updateExistsOnProc( this, scheme, geometry, myRank, procs, nProcs, rmflag)
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry information to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> This process's rank
    integer, intent(in) :: myRank
    !> Array of neighbor procs on which particle could possibly exist
    integer, intent(in) :: procs(:)
    !> Number of procs in array procs
    integer, intent(in) :: nProcs
    !> logical to set whether particle should be removed from this proc
    logical, intent(out) :: rmflag
    ! ------------------------------------------- !
    integer :: nx, ny, nz
    integer(kind=long_k) :: TreeID
    integer :: lev, upperBound
    integer :: coord(4)
    integer :: elemProc, iElemProc
    integer :: iproc

    logical :: oldExistsOnProc( 1:nProcs)
    ! ------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    upperBound = 2**lev
    ! Set rmflag to true initially. Will be set to false once we see the
    ! particle bounding box intersect this rank's elements
    rmflag = .TRUE.

    ! addToProc will be set to true if existsOnProc was FALSE in the last time step
    ! (so before updating it) and is TRUE during this time step
    this%addToProc = .FALSE.
    this%removeFromProc = .FALSE.

    ! removeFromProc should be set to true if this particle existed on
    ! proc in the last time step but no longer does in this time step
    ! So set it to true initially if existsOnProc = TRUE in the last time step
    do iproc = 1, nProcs
      oldExistsOnProc(iproc) = this%existsOnProc(iproc)
    end do ! iproc

    ! Now initialize existsOnProc for THIS time step to false
    this%existsOnProc = .FALSE.

    do nx = -this%Rn-1, this%Rn+1
      do ny = -this%Rn-1, this%Rn+1
        do nz = -this%Rn-1, this%Rn+1
          elemProc = -1
          iElemProc = -1

          ! integer coordinate of this element
          ! coord = this%coordOfOrigin(:) + (/ nx,ny,nz,0 /)
          coord = getNeighborCoord( this%coordOfOrigin, nx, ny, nz, pgBndData  )

          ! Check if coordinate is within actual simulation domain
          if( any(coord(1:3) < 0) .OR. any(coord(1:3) > upperBound) ) then
            cycle
          else
            ! get TreeID
            TreeID = tem_IdOfCoord(coord = coord)

            ! Check if element is local
            if( TreeID >= geometry%tree%Part_First(myRank + 1) &
              & .AND. TreeID <= geometry%tree%Part_Last(myRank + 1) ) then
              rmflag = .FALSE.
            else ! Element is not local, either halo or not on this proc
              ! Find the neighbor process that this element is on
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
            end if ! element is local
          end if ! coord within simulation domain
        end do
      end do
    end do

    ! Set addToProc and removeFromProc based on old and new values
    ! of existsOnProc
    do iproc = 1, nProcs
      if( (.NOT. oldExistsOnProc(iproc)) .AND. this%existsOnProc(iproc)  ) then
        this%addToProc(iproc) = .TRUE.
      else if( oldExistsOnProc(iproc) .AND. ( .NOT. this%existsOnProc(iproc) ) ) then
        this%removeFromProc(iproc) = .TRUE.
      end if
    end do ! iproc

  end subroutine updateExistsOnProc

  !> updateSolidNodes  modifies connectivity of scheme%pdf%neigh such
  !! that the particle cells do not participate in the streaming process.
  !! Specifically it:
  !!  * runs over the elements in this%exclusionList
  !!  * sets kernel neighbor list values to the particle elements themselves.
  !!  * This way they are effectively excluded from the streaming process.
  !!  * Sets the property bit prp_solid to 1 to indicate elements
  !!    belonging to a particle
  subroutine updateSolidNodes(this, scheme, stencil)
    ! ---------------------------------------------------------------------------
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Scheme for access to pdf%neigh
    type(mus_scheme_type), intent(inout) :: scheme
    !> fluid stencil
    type(tem_stencilHeader_type), intent(in) :: stencil

    ! ---------------------------------------------------------------------------
    integer :: iElem, iDir, zeroPos
    integer :: lev
    integer :: GetFromDir
    integer :: elemPos
    integer :: stateVarPos(stencil%QQ)
    integer :: nParticleElems
    integer :: nElems, nSize
    ! ---------------------------------------------------------------------------
    !write(logUnit(1),*) 'Updating state array connectivity for particles ...'
    ! FOR PARTICLES:
    ! Map PDF values of elements with particle property to themselves
    ! This way the compute kernel leaves them unchanged.

    ! write(logUnit(1), '(A)') 'updateSolidNodes'
    lev = this%coordOfOrigin(4)
    nElems = scheme%pdf( lev )%nElems_local
    nSize  = scheme%pdf( lev )%nSize
    ! write(logUnit(1), '(A,I12)') 'nSize = ', nSize

    ! state varpos for 1st field since neigh is created only for 1st field
    stateVarPos(:stencil%QQ) = scheme%varSys                    &
      &                              %method                    &
      &                              %val(scheme%stateVarMap    &
      &                              %varPos                    &
      &                              %val(1))                   &
      &                              %state_varPos( :stencil%QQ )

    nParticleElems = this%exclusionList%nvals

    ! First set the rest position for all elements belonging to this particle
    if ( stencil%restPosition > 0 ) then
      zeroPos = stencil%restPosition
      do iElem = 1, nParticleElems
        elemPos = this%exclusionList%val(iElem)

        scheme%pdf(lev)%neigh( ( zeropos-1)* nsize + elempos) &
          & = ( elempos-1)* scheme%varsys%nscalars+ zeropos
      end do
    end if

    ! Then loop over the other links for each element belonging to this particle
    elemLoop: do iElem = 1, nParticleElems

      ! Position of this element in the state array
      elemPos = this%exclusionList%val(iElem)


      ! Set property bit prp_solid to 1 to indicate element belongs to particle
      scheme%levelDesc(lev)%property( elemPos ) &
        & = ibset( scheme%levelDesc(lev)%property( elemPos ), prp_solid )

      neighLoop: do iDir = 1, stencil%QQN
        ! set direction to get the current direction itself
        GetFromDir = stateVarPos(iDir)


        scheme%pdf(lev)%neigh( ( idir-1)* nsize + elempos) &
          & = ( elempos-1)* scheme%varsys%nscalars+ getfromdir

      end do neighloop
    end do elemLoop
  end subroutine updateSolidNodes

  !> updateFluidNeighbors Runs over all particle elements
  !! (those in exclusionList) and then over each stencil direction.
  !!  When it finds a fluid neighbor in a direction it sets
  !!  the bounceback condition for that neighbor
  subroutine updateFluidNeighbors(this, scheme, stencil)
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> fluid stencil
    type(tem_stencilHeader_type), intent(in) :: stencil

    integer :: lev
    integer :: elemPos, neighPos, fluidPos
    integer :: iDir, fluidDir, fluidInvDir
    integer(kind=long_k) :: neighProp
    integer :: iElem, nElems, nSize
    integer :: stateVarPos(stencil%QQ)

    ! write(logUnit(1), '(A)') 'updateFluidNeighbors'

    lev = this%coordOfOrigin(4)
    nElems = scheme%pdf( lev )%nElems_local
    nSize  = scheme%pdf( lev )%nSize

    stateVarPos(:stencil%QQ) = scheme%varSys%method%val(           &
      &  scheme%stateVarMap%varPos%val(1))%state_varPos(:stencil%QQ)

    ! Loop over the elements belonging to this particle (those in exclusionList)
    ! Then identify fluid neighbors from there.
    particleElemLoop: do iElem = 1, this%exclusionList%nVals
      ! position of this node in levelDesc state vector
      elemPos = this%exclusionList%val(iElem)

      linkLoop: do iDir = 1, stencil%QQN
        ! Get the index of the neighbor element in this direction.
        ! -> ngDir because we have different neighbors for PUSH and PULL
        neighPos = scheme%levelDesc( lev )%neigh(1)%               &
          &        nghElems(  scheme%layout%fstencil %cxdirinv( idir),    &
          &                  elemPos )

        ! Initialize neighbor's property
        neighProp = 0_long_k

        ! If neighbor element exists, get its property
        if( neighPos > 0 ) neighProp = scheme%levelDesc(lev)%property(neighPos)

        if ( neighPos <= 0 .or. (btest(neighProp, prp_solid)) ) then
          ! this link either points to a BC or another solid element so skip it
          cycle
        else
          ! This link points to a fluid neighbor
          ! We now modify the entry in neigh(:) for the fluid neighbor
          fluidPos = neighPos

          ! direction to set bounceback FROM fluid TO solid = inverse of iDir
          fluidDir = scheme%layout%fStencil%cxDirInv(stateVarPos(iDir))

          ! set bounceback: stream from link in inverse direction of same element
          fluidInvDir = scheme%layout%fStencil%cxDirInv(fluidDir)
          ! scheme%pdf( lev )%neigh( ( fluiddir-1)* nelems + fluidpos) =    &
          !   &  ( fluidpos-1)* scheme%varsys%nscalars+ fluidinvdir
          scheme%pdf( lev )%neigh( ( fluiddir-1)* nsize + fluidpos) =    &
            &  ( fluidpos-1)* scheme%varsys%nscalars+ fluidinvdir

        end if
      end do linkLoop
    end do particleElemLoop
  end subroutine updateFluidNeighbors

  !> updateNewFluidNodes restores the correct connectivity for
  !! new fluid elements which were particle elements in the previous time step.
  subroutine updateNewFluidNodes(this, scheme, stencil)
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> fluid stencil
    type(tem_stencilHeader_type), intent(in) :: stencil

    integer :: lev
    integer :: elemPos, neighPos
    integer :: iDir, inv_iDir, fluidDir
    integer(kind=long_k) :: neighProp
    integer :: iElem, nElems, nSize
    integer :: stateVarPos(stencil%QQ)

    ! write(logUnit(1), '(A)') 'updateNewFluidNodes'
    lev = this%coordOfOrigin(4)
    nElems = scheme%pdf( lev )%nElems_local
    nSize  = scheme%pdf( lev )%nSize

    stateVarPos(:stencil%QQ) = scheme%varSys                    &
      &                              %method                    &
      &                              %val(scheme%stateVarMap    &
      &                              %varPos                    &
      &                              %val(1))                   &
      &                              %state_varPos( :stencil%QQ )

    particleElemLoop: do iElem = 1, this%makeFluidList%nVals
      ! position of this node in levelDesc state vector
      elemPos = this%makeFluidList%val(iElem)

      linkLoop: do iDir = 1, stencil%QQN
        ! Get the index of the neighbor element in this direction.
        ! -> ngDir because we have different neighbors for PUSH and PULL
        neighPos = scheme%levelDesc( lev )%neigh(1)%               &
          &        nghElems(  scheme%layout%fstencil %cxdirinv( idir),    &
          &                  elemPos )

        neighProp = 0_long_k

        if( neighPos > 0 ) then
          neighProp = scheme%levelDesc(lev)%property(neighPos)
        end if
        ! We also need to set the connectivity of the NEIGHBORS
        ! of the new fluid element. This allows the surrounding elements
        ! to stream from the new fluid element

        if ( neighPos <= 0 .or. (btest(neighProp, prp_solid)) ) then
          ! set bounceback for the new fluid element and then go to next
          ! iDir iteration
          inv_iDir = scheme%layout%fStencil%cxDirInv(stateVarPos(iDir))

          scheme%pdf( lev )%neigh( ( idir-1)* nsize + elempos) =    &
                  &  ( elempos-1)* scheme%varsys%nscalars+ inv_idir
          cycle
        else
          ! Set connectivity for the NEW fluid element
          ! This allows the new fluid element to stream/pull from its
          ! surrounding elements

          scheme%pdf( lev )%neigh( ( idir-1)* nsize + elempos) =    &
                  &  ( neighpos-1)* scheme%varsys%nscalars+ idir

          ! Now modify the entry in neigh(:) for the neighbor of our new
          ! fluid element. The direction to set FROM neighbor TO new
          ! fluid element = inverse of iDir

          fluidDir = scheme%layout%fStencil%cxDirInv(stateVarPos(iDir))

          scheme%pdf( lev )%neigh( ( fluiddir-1)* nsize + neighpos) =    &
            &          ( elempos-1)* scheme%varsys%nscalars+ fluiddir

        end if
      end do linkLoop
    end do particleElemLoop
  end subroutine updateNewFluidNodes


  !> destroyParticle_MEM removes a particle from the lattice and the particleGroup
  !! It will
  !! * empty the exclusionList
  !! * initialize all the former particle elements to new fluid elements
  !! * restore connectivity for all these new fluid elements
  !! This routine should be followed up by a call to remove_particle_from_da_particle_MEM!
  subroutine destroyParticle_MEM(this, particleGroup, scheme, params )
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Particle group in which to search for collisions
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Parameters for access to dx and dt
    type(mus_param_type), intent(in) :: params
    ! ---------------------------------------------------!
    integer :: lev
    integer :: elemPos
    integer :: iElem

    ! macroscopic quantities for intializing new fluid elems to eq PDF
    real(kind=rk), allocatable :: rho(:)                ! size (Nelems)
    real(kind=rk), allocatable :: vel_surf_lat(:,:)     ! size (3,Nelems)

    real(kind=rk) :: dt, dx
    real(kind=rk) :: inv_dx, inv_vel
    real(kind=rk) :: baryOfOrigin(3)
    real(kind=rk) :: baryOfSurface(3)

    ! vector from particle origin to surface element, lattice units
    real(kind=rk) :: r_lat(3)
    ! particle angular velocity vector in lattice units
    real(kind=rk) :: om_lat(3)
    ! ---------------------------------------------------!
    lev = this%coordOfOrigin(4)
    dt = params%physics%dtLvl(lev)
    dx = params%physics%dxLvl(lev)
    inv_dx = 1.0_rk / dx
    inv_vel = 1.0_rk / params%physics%fac(lev)%vel

    call empty( me = this%makeFluidList )

    ! To destroy the particle, all elements should be initialized to
    ! new fluid elements so add them to the makeFluidList
    do iElem = 1, this%exclusionList%nVals
        call append( me   = this%makeFluidList,           &
                   &  val = this%exclusionList%val(iElem) )
    end do

    baryOfOrigin = this%pos(1:3)

    ! Compute rotational velocity of particle
    om_lat = this%vel(4:6) * dt

    allocate( rho( this%makeFluidList%nVals ) )
    allocate( vel_surf_lat( 3, this%makeFluidList%nVals ) )

    ! Use reference density for feq calculation as in Bartuschat
    rho(1:this%makeFluidList%nVals) = rho0_lat

    ! Compute surface velocity values at each new fluid node
    do iElem = 1,this%makeFluidList%nVals
      elemPos = this%makeFluidList%val(iElem)
      scheme%levelDesc(lev)%property( elemPos ) = &
          &  ibclr( scheme%levelDesc(lev)%property( elemPos ), prp_solid )

      baryOfSurface = scheme%levelDesc(lev)%baryOfTotal(elemPos,1:3)
      r_lat = computeDisplacement(baryOfOrigin, baryOfSurface, pgBndData)
      r_lat = r_lat * inv_dx

      ! Compute surface velocity and store in vel_surf_lat
      call cross_product(om_lat, r_lat, vel_surf_lat(:,iElem))

      vel_surf_lat(:,iElem) = vel_surf_lat(:,iElem) + this%vel(1:3) * inv_vel
    end do

    ! set PDF of new fluid elems to fEq value using surface velocity values
    call setToEquilibrium( elemList = this%makeFluidList                      &
                         &                %val(1:this%makeFluidList%nVals),   &
                         &  N       = this%makeFluidList%nVals,               &
                         &  lev     = lev,                                    &
                         &  scheme  = scheme,                                 &
                         &  rho     = rho,                                    &
                         &  vel     = vel_surf_lat                            )

    deallocate(rho)
    deallocate(vel_surf_lat)

    ! restore connectivity of new fluid elements to neighbor elems and vice-versa
    call updateNewFluidNodes(this, scheme, scheme%layout%fStencil)
    ! -------------------------------------------------!

  end subroutine destroyParticle_MEM


  !> ApplyVelocityBounceback adds the momentum exchange term to the
  !! distribution function at the particle  neighboring fluid nodes
  !! this term accounts for the moving boundary of the particle.
  subroutine applyVelocityBounceback(this, scheme, stencil, params)
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> fluid stencil
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> Parameters for access to dt, dx, conversion factors between units
    type(mus_param_type), intent(in) :: params
    ! ------------------------------------------!
    ! All variables with suffix _lat in lattice units
    integer :: lev
    integer :: elemPos, neighPos, fluidPos
    integer :: iDir, fluidDir, iField, idx_idir
    integer(kind=long_k) :: neighProp
    integer :: iElem, nElems, nSize
    integer :: stateVarPos(stencil%QQ)
    integer :: nNow, nNext
    real(kind=rk) :: dt, dx, inv_dx, inv_vel
    real(kind=rk) :: wq, MoBoFactor, movingBoundaryTerm
    real(kind=rk) :: u_surf_lat(3)
    real(kind=rk) :: baryOfOrigin(3), baryOfSurface(3), cq(3)

    ! vector particle origin to surface elem and angular velocity
    real(kind=rk) :: r(3), r_lat(3), om_lat(3)
    ! ------------------------------------------!
    ! write(logUnit(1),'(A)') 'Inside applyParticleVelocityBounceback'
    lev = this%coordOfOrigin(4)
    iField = 1

    dt = params%physics%dtLvl(lev)
    dx = params%physics%dxLvl(lev)
    inv_dx = 1.0_rk / dx
    inv_vel = 1.0_rk / params%physics%fac(lev)%vel
    om_lat = this%vel(4:6) * dt
    MoBoFactor = -2.0 * cs2inv * rho0_lat ! prefactor for moving boundary term

    baryOfOrigin = this%pos(1:3)

    nElems = scheme%pdf( lev )%nElems_local
    nSize  = scheme%pdf( lev )%nSize
    nNow   = scheme%pdf( lev )%nNow
    nNext  = scheme%pdf( lev )%nNext

    stateVarPos(:stencil%QQ) = scheme%varSys                    &
      &                              %method                    &
      &                              %val(scheme%stateVarMap    &
      &                              %varPos                    &
      &                              %val(1))                   &
      &                              %state_varPos( :stencil%QQ )

    ! Loop over the elements belonging to this particle (those in exclusionList)
    ! Then identify fluid neighbors from there.
    particleElemLoop: do iElem = 1, this%exclusionList%nVals

    ! position of this node in levelDesc state vector
    elemPos = this%exclusionList%val(iElem)

    baryOfSurface = scheme%levelDesc(lev)%baryOfTotal(elemPos,1:3)
      linkLoop: do iDir = 1, stencil%QQN
      ! Get the index of the neighbor element in this direction.
      neighPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, elemPos)
      ! Initialize neighbor's property
      neighProp = 0_long_k

      ! If neighbor element exists, get its property
      if( neighPos > 0 ) neighProp = scheme%levelDesc(lev)%property(neighPos)

      if ( neighPos <= 0 .or. (btest(neighProp, prp_solid)) ) then
      ! this link either points to a BC or another solid element so skip it
        cycle
      else
        ! this link points to a fluid neighbor
        fluidPos = neighPos

        ! fluidDir = direction FROM fluid TO particle = cq in Bartuschat eqn 5.2
        fluidDir = scheme%layout%fStencil%cxDirInv(stateVarPos(iDir))

        ! get vector for stencil direction fluidDir
        ! NB: this is in lattice units!
        cq = stencil%cxDirRK(:,fluidDir)

        ! vector from particle origin to current point on the surface
        r = baryOfSurface - baryOfOrigin
        call calcPeriodicRsurface( r            = r,           &
                                 & R_particle   = this%radius, &
                                 & boundaryData = pgBndData    )

        ! Convert to lattice units
        r_lat = r * inv_dx - 0.5*cq

        ! Correct for possible periodic boundary conditions

        ! Stencil weight in this direction
        wq = scheme%layout%weight(fluidDir)

        idx_idir = scheme%varSys%method%val(                      &
          &               scheme%stateVarMap%varPos%val(iField) ) &
          &                     %state_varPos(fluidDir)


        ! compute rotational component of surface velocity using cross product
        call cross_product( om_lat, r_lat, u_surf_lat )

        ! add translational component
        u_surf_lat = u_surf_lat + this%vel(1:3) * inv_vel

        ! compute momentum exchange term due to moving boundary
        movingBoundaryTerm = MoBoFactor * wq  * dot_product(cq, u_surf_lat)

        ! Add the extra term to the pdf from fluid to particle
        ! to transfer momentum FROM particle TO fluid

        scheme%state(lev)%val(                                                  &
          & ( fluidpos-1)* scheme%varsys%nscalars+idx_idir,nNext)     &
          &   = scheme%state(lev)%val(                                          &
          &     ( fluidpos-1)* scheme%varsys%nscalars+idx_idir,nNext) &
          &     + movingBoundaryTerm


      end if
      end do linkLoop
    end do particleElemLoop

  end subroutine applyVelocityBounceback

  !> handleParticleOverlap checks for overlap with discrete representations
  !! of other particles. Returns overlapsOtherParticle = TRUE if
  !! there is overlap. This does not require communication with other
  !! processes.
  subroutine checkForParticleOverlap(this, particleGroup, scheme, &
      &                            overlapsOtherParticle        )
    !> Particle we are currently checking overlap for
    type(mus_particle_MEM_type), intent(inout) :: this
    !> particleGroup to search for overlapping particles in
    type(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> logical indicating if a overlap has been detected (1) or not (0)
    logical, intent(out) :: overlapsOtherParticle
    ! ------------------------------------------------------ !

    integer :: lev, iElem
    integer :: collPartIndex, collPartID
    integer(kind=long_k) :: elemProp
    integer :: elemPos

    lev = this%coordOfOrigin(4)
    elemProp = 0_long_k
    overlapsOtherParticle = .FALSE.
    ! Check the new solid elements (those in exclusionList but NOT in
    ! exclusionListBuffer) for prp_solid
    ! If any of the new solid elements already has prp_solid this indicates
    ! overlap with another particle
    checkOverlapLoop: do iElem = 1, this%exclusionList%nVals
      if( .not.( any( this%exclusionList%val(iElem) &
        &             == this%exclusionListBuffer                    &
        &                    %val(1:this%exclusionListBuffer%nVals)))) then
        ! For elements that were newly added to this particle,
        ! check if they already belong to some other particle

        elemProp = scheme%levelDesc(lev)%property( this%exclusionList%val(iElem) )
        elemPos = this%exclusionList%val(iElem)

        if( btest(elemProp, prp_solid) ) then
          ! Determine origin and radius of neighbor particle
          call findParticleFromElem( particleGroup     = particleGroup,   &
                                   & elemPos           = elemPos,         &
                                   & thisParticleIndex = this%particleID, &
                                   & findParticleIndex = collPartIndex,   &
                                   & particleID        = collPartID       )

          open( pgDebugLog%lu, file=pgDebugLog%lfile, &
              & status='old', position='append' )
            write(pgDebugLog%lu,'(A)') 'Overlap with particle at index:'
            write(pgDebugLog%lu, *) collPartIndex
          close( pgDebugLog%lu )

          if (collPartIndex == 0) then
            call tem_abort( 'Error: could not find overlapping particle')
          end if

          overlapsOtherParticle = .TRUE.
          exit checkOverlapLoop
        end if
      end if
    end do checkOverlapLoop
  end subroutine checkForParticleOverlap


  !! -------------------- AUXILIARY ROUTINES ------------------------!!
  !> setToEquilibrium sets elements to equilibrium distribution
  subroutine setToEquilibrium(elemList, N, lev, scheme, rho, vel)
    !> list containing positions in state vector of elements
    !! to set to equilibrium distribution
    integer, intent(in) :: elemList(:)
    !> number of values in elemList
    integer, intent(in) :: N
    !> level
    integer, intent(in) :: lev
    !> scheme for access to levelDesc and pdf
    type(mus_scheme_type), intent(inout) :: scheme
    !> density value for computing fEq
    real(kind=rk), intent(in) :: rho(:)
    !> velocity vector
    real(kind=rk), intent(in) :: vel(:,:)

    real(kind=rk), allocatable :: fEq(:)
    integer :: nSize
    integer :: QQ,QQN ! number of stencil directions
    integer :: iElem, iDir, iField, zeroPos, idx_idir, offset
    integer :: elemPos
    integer :: nNext

    ! write(logUnit(1), '(A)') 'setToEquilibrium'
    ! --- Compute equilibrium distribution from prescribed rho, vel_arr
    zeroPos = scheme%layout%fStencil%restPosition
    QQN = scheme%layout%fStencil%QQN
    QQ = scheme%layout%fStencil%QQ
    iField = 1
    nNext = scheme%pdf( lev )%nNext

    allocate( fEq( N * QQ ) )

    ! Use rho(:) and vel_arr(:,:) to compute fEq(:)
    call scheme%derVarPos(iField)%equilFromMacro( density  = rho(1:N),      &
                                                & velocity = vel(1:3,1:N),  &
                                                & iField   = iField,        &
                                                & nElems   = N,             &
                                                & varSys   = scheme%varSys, &
                                                & layout   = scheme%layout, &
                                                & res      = fEq            )

    ! Now set elems in elemList to computed fEq
    nSize  = scheme%pdf( lev )%nSize

    do iElem = 1, N
      elemPos = elemList(iElem)   ! position of element in state array
      offset = (iElem-1)*QQ

      ! Set equilibrium for zero (rest) direction
      idx_idir = scheme%varSys%method%val(                      &
          &             scheme%stateVarMap%varPos%val(iField) ) &
          &                   %state_varPos(zeroPos)

      scheme%state(lev)%val( &
        &  ( elempos-1)* scheme%varsys%nscalars+idx_idir, nNext)  &
        &     =  fEq(offset + zeroPos)


      ! Then for the other directions
      do iDir = 1,QQN
        idx_idir = scheme%varSys%method%val(                        &
            &               scheme%stateVarMap%varPos%val(iField) ) &
            &                     %state_varPos(iDir)


        scheme%state(lev)%val(                                              &
          & ( elempos-1)* scheme%varsys%nscalars+idx_idir, nNext) &
          &    =  fEq(offset + iDir)

      end do
    end do

    deallocate(fEq)

  end subroutine setToEquilibrium


  !> make_pdf_tiny sets the state vector (particle distribution function) values
  !! to some tiny value. This is done for visualization purposes. It has no effect
  !! on the flow as the particle elements do not participate in the
  !! stream-and-collide process.
  subroutine make_pdf_tiny(scheme,particle,lev)
    type(mus_scheme_type) :: scheme
    type(mus_particle_MEM_type) :: particle
    integer :: lev


    integer :: nElems, nSize, nScalars
    integer :: iDir, iField, iElem, idx_idir, QQ
    nElems = particle%exclusionList%nvals

    nSize  = scheme%pdf( lev )%nSize
    nScalars  = scheme%varSys%nScalars
    QQ = scheme%layout%fStencil%QQ

    do iElem = 1,nElems
      do iDir = 1,QQ
        do iField = 1, scheme%nFields
          idx_idir = scheme%varSys%method%val(                      &
            &               scheme%stateVarMap%varPos%val(iField) ) &
            &                     %state_varPos(iDir)
         ! Set all the pdf links of target element to an infinitesimally
         ! small value. This will prevent density from becoming NaN if pdfs
         ! are set to zero and besides changed elements will have a very
         ! small density which will help in visualization and differentiate
         ! them from other fluid elements

         scheme%state(lev)%val(                                          &
         & ( particle%exclusionlist%val(ielem)-1)* nscalars+idx_idir,:) &
         & = 0.000000001_rk
        end do
      end do
    end do
  end subroutine make_pdf_tiny

  !> find particle index and ID in particleGroup array
  !  for an element with position elemPos in state array
  subroutine findParticleFromElem( particleGroup, elemPos, thisParticleIndex, &
                                 & findParticleIndex, particleID              )
    !> Array of particles
    type(mus_particle_group_type), intent(in) :: particleGroup
    !> Index of element for which to find particle in state array
    integer, intent(in) :: elemPos

    !> Index of particle that was found in particleGroup%particles
    !  Returns 0 if not found
    integer :: thisParticleIndex
    integer :: findParticleIndex
    !> ID of particle that was found
    integer :: particleID
    !---------------------------------------------------!

    integer iParticle, iElem

    findParticleIndex = 0
    particleID = 0

    do iParticle = 1, particleGroup%particles_MEM%nvals
      if (iParticle == thisParticleIndex) then
        cycle
      end if

      do iElem = 1, particleGroup%particles_MEM%val(iParticle)%exclusionList%nVals
        if( elemPos == particleGroup%particles_MEM%val(iParticle)%exclusionList &
          &                                              %val(iElem)    ) then
          ! If particle to which this element belongs is found, get its
          ! index and ID, then exit function
          findParticleIndex = iParticle
          particleID = particleGroup%particles_MEM%val(iParticle)%particleID
          return
        end if
      end do
    end do

  end subroutine findParticleFromElem

  ! !!!! Unused Code ???
  !HK! !> computeCellMomentum computes the momentum (first-order moment) of
  !HK! !! a fluid cell with position elemPos in the total list
  !HK! subroutine computeCellMomentum( elemPos, scheme, stencil, params, lev, j_phy )
  !HK!   !> Position of this element in the levelDesc total list
  !HK!   integer, intent(in) :: elemPos
  !HK!   !> Scheme for access to level descriptor
  !HK!   type(mus_scheme_type), intent(inout) :: scheme
  !HK!   !> fluid stencil
  !HK!   type(tem_stencilHeader_type), intent(in) :: stencil
  !HK!   !> Parameters for access to conversion factors
  !HK!   type(mus_param_type), intent(in) :: params
  !HK!   !> Level of cell
  !HK!   integer :: lev
  !HK!   !> Output momentum vector j_phy = (jx, jy, jz). Physical units.
  !HK!   real(kind=rk) :: j_phy(3)
  !HK!   ! --------------------------------------------- !
  !HK!   integer :: iDir, idx_idir
  !HK!   integer :: nSize, nNext
  !HK!   real(kind=rk) :: pdf_val
  !HK!   real(kind=rk) :: cq_lat(3)    ! stencil direction in lat units
  !HK!   real(kind=rk) :: j_lat(3)     ! momentum vector in lat units
  !HK!   real(kind=rk) :: dx           ! mesh size
  !HK!   ! --------------------------------------------- !
  !HK!   nSize  = scheme%pdf( lev )%nSize
  !HK!   nNext  = scheme%pdf( lev )%nNext
  !HK!   dx = params%physics%dxLvl(lev)
  !HK!
  !HK!   j_lat = 0.0_rk
  !HK!   ! Loop over directions
  !HK!   do iDir = 1, stencil%QQ
  !HK!     cq_lat = stencil%cxDirRK(:,iDir)
  !HK!
  !HK!     idx_idir = scheme%varSys%method%val(                      &
  !HK!     &               scheme%stateVarMap%varPos%val(1) ) &
  !HK!     &                     %state_varPos(iDir)
  !HK!
  !HK!     ! Get the pdf value in this direction
  !HK!     pdf_val =  scheme%state(lev)%val(                                    &
  !HK!       &   ( elempos-1)* scheme%varsys%nscalars+idx_idir, nNext)
  !HK!
  !HK!     j_lat = j_lat + cq_lat(1:3) * pdf_val
  !HK!   end do
  !HK!
  !HK!   j_phy = j_lat * dx**3 * params%physics%rho0 * params%physics%fac(lev)%vel
  !HK!
  !HK! end subroutine computeCellMomentum
  !HK!
  !HK! subroutine getAverageNeighborDensity(elemPos, scheme, stencil, lev, rho_lat )
  !HK!   !> Position of this element in the levelDesc total list
  !HK!   integer, intent(in) :: elemPos
  !HK!   !> Scheme for access to level descriptor
  !HK!   type(mus_scheme_type), intent(inout) :: scheme
  !HK!   !> fluid stencil
  !HK!   type(tem_stencilHeader_type), intent(in) :: stencil
  !HK!   !> Level that particle is on
  !HK!   integer, intent(in) :: lev
  !HK!   !> Output: average density, in lattice units
  !HK!   real(kind=rk), intent(out) :: rho_lat
  !HK!   ! --------------------------------------------!
  !HK!   integer :: iDir
  !HK!   integer :: neighPos
  !HK!   integer(kind=long_k) :: neighProp
  !HK!   integer :: nLocalElems
  !HK!   integer :: Navg          ! counter how many elems in our average
  !HK!   integer :: dens_pos, elemOff
  !HK!   real(kind=rk) :: rho_el_lat
  !HK!   ! --------------------------------------------!
  !HK!   dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
  !HK!   nLocalElems = scheme%pdf( lev )%nElems_fluid
  !HK!
  !HK!   Navg = 0
  !HK!   rho_lat = 0.0_rk
  !HK!
  !HK!   do iDir = 1, stencil%QQN
  !HK!       ! Get the index of the neighbor element in this direction.
  !HK!       neighPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, elemPos)
  !HK!             ! Initialize neighbor's property
  !HK!       neighProp = 0_long_k
  !HK!
  !HK!       ! If neighbor element exists, get its property
  !HK!       if( neighPos > 0 ) neighProp = scheme%levelDesc(lev)%property(neighPos)
  !HK!
  !HK!       if ( neighPos <= 0 .or. neighPos > nLocalElems  &
  !HK!         & .or. (btest(neighProp, prp_solid)) ) then
  !HK!         ! this link either points to a BC or a halo or another solid element
  !HK!         ! so skip it
  !HK!         cycle
  !HK!
  !HK!       else
  !HK!         ! Get density of this fluid node
  !HK!         elemOff = (neighPos-1)*scheme%varSys%nAuxScalars
  !HK!         rho_el_lat = scheme%auxField(lev)%val( elemOff + dens_pos )
  !HK!         rho_lat = rho_lat + rho_el_lat
  !HK!         Navg = Navg + 1
  !HK!       end if
  !HK!   end do ! iDir
  !HK!
  !HK!   ! Take the average
  !HK!   if( Navg > 0 ) then
  !HK!     rho_lat = rho_lat / Navg
  !HK!   else
  !HK!     ! If there are no fluid elems to take the average from,
  !HK!     ! use the reference density.
  !HK!     rho_lat = rho0_lat
  !HK!   end if
  !HK!
  !HK! end subroutine getAverageNeighborDensity
end module mus_particle_MEM_module
