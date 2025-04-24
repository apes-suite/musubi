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
!> mus_particle_module contains the main control routines for LBM-DEM simulations
!! of particles in a flow.

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

module mus_particle_module

  use mpi

  use env_module, only : rk, double_k, long_k, newunit

  use tem_param_module,        only: rho0_lat => rho0, cs2, cs2inv
  use tem_logging_module,      only: logUnit
  use tem_aux_module,          only: tem_abort
  use tem_geometry_module,     only: tem_CoordOfReal, tem_PosOfId
  use tem_topology_module,     only: tem_IdOfCoord, tem_coordOfId, &
    &                                tem_FirstIdAtLevel
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

  use mus_particle_type_module, only:         &
    &   mus_particle_MEM_type,                &
    &   mus_particle_DPS_type,                &
    &   init_da_particle_DPS,                 &
    &   mus_particle_group_type,              &
    &   allocateProcessMasks,                 &
    &   printParticleGroup,                   &
    &   printParticleGroup2_MEM,              &
    &   printParticleGroup2_DPS,              &
    &   printpIDlist,                         &
    &   remove_particle_from_da_particle_MEM, &
    &   remove_particle_from_da_particle_DPS

  use mus_particle_comm_type_module, only:          &
    &   init_mus_particles_comm_type,               &
    &   mus_particles_comm_init_vectorbuffer,       &
    &   mus_particles_comm_init_wallbuffer,         &
    &   mus_particles_comm_init_statebuffer,        &
    &   mus_particles_comm_init_posbuffer,          &
    &   mus_particles_comm_init_IDbuffer,           &
    &   mus_particles_comm_init_particlebuffer,     &
    &   mus_particles_initForceContributionMPItype, &
    &   mus_particles_initParticleStateMPItype,     &
    &   mus_particles_initPositionUpdateMPItype,    &
    &   mus_particles_initWallPosMPItype,           &
    &   mus_particles_initParticleInfoMPItype,      &
    &   print_particles_comm,                       &
    &   print_particles_pIDvectorbuffer,            &
    &   find_particle_comm_procs4,                  &
    &   mus_particles_communication_type,           &
    &   mus_pIDvector_type,                         &
    &   mus_particleState_type,                     &
    &   mus_particleInfo_type

  use mus_particle_comm_module, only:  &
    &   exchangeForces,                &
    &   exchangeParticlesToRemove,     &
    &   exchangeParticlesToRemove_DPS, &
    &   exchangeVelocities,            &
    &   exchangeParticleStates,        &
    &   exchangeNewParticles_MEM,      &
    &   exchangeNewParticles_DPS,      &
    &   exchangeHydroForces_DPS,       &
    &   mus_particles_initialize_communication

  use mus_particle_MEM_module, only: &
    &   mapToLattice,                &
    &   initParticle_MEM,            &
    &   applyHydrodynamicForces,     &
    &   updateCoordOfOrigin,         &
    &   updateExclusionList,         &
    &   updateSolidNodes,            &
    &   updateFluidNeighbors,        &
    &   updateNewFluidNodes,         &
    &   destroyParticle_MEM,         &
    &   applyVelocityBounceback,     &
    &   checkForParticleOverlap,     &
    &   setToEquilibrium,            &
    &   make_pdf_tiny

  use mus_particle_DPS_module, only:                   &
    &   initParticle_DPS,                              &
    &   mapToLattice_DPS,                              &
    &   applyDragForce_DPS,                            &
    &   applyDragForce_DPS_noeps,                      &
    &   applyLiftForce_DPS,                            &
    &   applyPressureForce_DPS,                        &
    &   transferMomentumToFluid_DPS,                   &
    &   mus_particles_updateFluidVolumeFraction,       &
    &   transferMomentumToFluid_DPS_twoway,            &
    &   interpolateFluidProps,                         &
    &   interpolateFluidProps_onewaycoupled,           &
    &   calcVelocityAndPressureGradient,               &
    &   calcVelocityAndPressureGradient_onewaycoupled, &
    &   addParticleSourceToAuxField_DPS,               &
    &   mus_particles_DPS_interpolateFluidProperties,  &
    &   mus_particles_DPS_interpolateFluidProperties_onewaycoupled

  use mus_particle_interpolation_module, only: init_particle_interpolator

  use mus_particle_aux_module, only: &
    &   findPartitionOfTreeID,       &
    &   cross_product,               &
    &   getProcessBoundingBox

  use mus_particle_logging_module, only: &
    &   closeParticleLog,                &
    &   getParticleLogUnit,              &
    &   generateElemListLine

  use mus_particle_logging_type_module, only: &
    &   mus_particle_logging_type,            &
    &   pgDebugLog,                           &
    &   init_particle_logger

  use mus_particle_DEM_module, only: &
    &   DEMSubcycles_MEM,            &
    &   DEMSubcycles_DPS,            &
    &   DEMSubcycles_DPS_onewaycoupled

  use mus_particle_boundary_module, only: pgBndData
  use mus_particle_creator_module, only: create_particles_DPS, particle_creator

  implicit none


contains


!> Swap index of the particle force buffer
subroutine swapFBuff(this)
  type(mus_particle_MEM_type), intent(inout) :: this

  ! ------------------------------------------!

  ! Set Fnow for next time step
  if(this%Fnow == 1) then
    this%Fnow = 2
    this%Flast = 1
  else
    this%Fnow = 1
    this%Flast = 2
  end if
end subroutine swapFBuff

! *************************************************************************** !
!> Initialization for particleGroup and all the particles in it
!! Includes:
!! * Assigning the required procedure pointers for particleGroup and particles
!! * Initializing loggers
!! * Initializing communication routines
!! * Building the representation of the particles on the grid
subroutine mus_particles_initialize( particleGroup, scheme, &
                                   & geometry, params       )
  !> Array of particles
  type(mus_particle_group_type), target :: particleGroup
  !> Scheme for access to leveldescriptor
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  !---------------------------------------------------!
  integer :: lev, nBCs, BCid
  integer :: iParticle, iBC
  real(kind=rk) :: dt

  ! For cubic periodic validation cases only
  integer :: bnd_coord(4)
  real(kind=rk) :: bnd(3), dx
  ! --------------------------------------------------!
  lev = geometry%tree%global%maxLevel
  dt = params%physics%dtLvl(lev)
  dx = params%physics%dxLvl(lev)

  ! ------ SET PROCEDURES TO USE FOR PARTICLEGROUP ------ !
  select case( trim(params%particle_kind) )
    case('MEM')
      particleGroup%applyHydrodynamicForces => mus_particles_applyHydrodynamicForces_MEM
      particleGroup%moveParticles => mus_particles_move
      particleGroup%mapParticles => mus_particles_mapping_MEM
      particleGroup%transferMomentumToFluid => mus_particles_transferMomentumToFluid_MEM
    case('DPS', 'DPS_unittest')
      particleGroup%applyHydrodynamicForces => mus_particles_applyHydrodynamicForces_DPS
      particleGroup%transferMomentumToFluid => mus_particles_transferMomentumToFluid_DPS
      particleGroup%addParticleSourcesToAuxField => mus_particles_addSourceTermsToAuxField_DPS
      particleGroup%transfer_momentum => transferMomentumToFluid_DPS
      particleGroup%modify_auxfield => addParticleSourceToAuxfield_DPS
      particleGroup%intp => interpolateFluidProps
      particleGroup%calc_vel_and_p_grad => calcVelocityAndPressureGradient
      particleGroup%moveParticles => mus_particles_move_DPS
      particleGroup%mapParticles => mus_particles_mapping_DPS
      particleGroup%calcDragForce => applyDragForce_DPS
      particleGroup%calcLiftForce => applyLiftForce_DPS
      particleGroup%calcPressureForce => applyPressureForce_DPS
    case('DPS_twoway')
      particleGroup%applyHydrodynamicForces => null()
      particleGroup%moveParticles => mus_particles_move_DPS
      particleGroup%mapParticles => mus_particles_mapping_DPS
      particleGroup%transferMomentumToFluid => mus_particles_transferMomentumToFluid_DPS
      particleGroup%addParticleSourcesToAuxField => mus_particles_addSourceTermsToAuxField_DPS
      particleGroup%transfer_momentum => transferMomentumToFluid_DPS_twoway
      particleGroup%modify_auxfield => addParticleSourceToAuxfield_DPS
      particleGroup%intp => interpolateFluidProps_onewaycoupled
      particleGroup%calc_vel_and_p_grad => calcVelocityAndPressureGradient_onewaycoupled
      particleGroup%calcDragForce => applyDragForce_DPS_noeps
      particleGroup%calcLiftForce => applyLiftForce_DPS
      particleGroup%calcPressureForce => applyPressureForce_DPS
    case('DPS_oneway')
      particleGroup%applyHydrodynamicForces => mus_particles_applyHydrodynamicForces_DPS_onewaycoupled
      particleGroup%moveParticles => mus_particles_move_DPS_onewaycoupled
      particleGroup%mapParticles => mus_particles_mapping_DPS
      particleGroup%transferMomentumToFluid => null()
      particleGroup%transfer_momentum => null()
      particleGroup%modify_auxfield => null()
      particleGroup%intp => interpolateFluidProps_onewaycoupled
      particleGroup%calc_vel_and_p_grad => calcVelocityAndPressureGradient_onewaycoupled
      particleGroup%calcDragForce => applyDragForce_DPS_noeps
      particleGroup%calcLiftForce => applyLiftForce_DPS
      particleGroup%calcPressureForce => applyPressureForce_DPS
    case default
      write(logUnit(1),*) 'ERROR mus_particles_initialize: unknown particle kind!'
      call tem_abort()
      particleGroup%applyHydrodynamicForces => mus_particles_applyHydrodynamicForces_DPS
      particleGroup%moveParticles => mus_particles_move_DPS
      particleGroup%mapParticles => mus_particles_mapping_DPS
      particleGroup%transferMomentumToFluid => mus_particles_transferMomentumToFluid_DPS
      particleGroup%intp => interpolateFluidProps
      particleGroup%calc_vel_and_p_grad => calcVelocityAndPressureGradient
      particleGroup%calcDragForce => applyDragForce_DPS_noeps
      particleGroup%calcLiftForce => applyLiftForce_DPS
      particleGroup%calcPressureForce => applyPressureForce_DPS
    end select ! select particle kind

  particleGroup%exchangeForces => exchangeForces
  particleGroup%exchangeParticleStates => exchangeParticleStates
  particleGroup%exchangeNewParticles => exchangeNewParticles_MEM


  ! ------ INITIALIZE PARTICLE DEBUG LOGGER  -------- !
  call init_particle_logger( pgDebugLog, params%general%proc%rank )

  ! ------ INITIALIZE BOUNDARY INFORMATION FOR PARTICLES ------ !
  ! Set integer coordinates of the particle domain boundaries
  if( pgBndData%useBnd ) then
    bnd = (/ pgBndData%bnd(1), pgBndData%bnd(3), pgBndData%bnd(5) /) + 0.5*dx
    bnd_coord = tem_coordOfReal( mesh  = geometry%tree,   &
                               & point = bnd,             &
                               & level = lev              )
    pgBndData%bnd_coord(1) = bnd_coord(1)
    pgBndData%bnd_coord(3) = bnd_coord(2)
    pgBndData%bnd_coord(5) = bnd_coord(3)

    bnd = (/ pgBndData%bnd(2), pgBndData%bnd(4), pgBndData%bnd(6) /) - 0.5*dx
    bnd_coord = tem_coordOfReal( mesh  = geometry%tree, &
                               & point = bnd,           &
                               & level = lev            )
    pgBndData%bnd_coord(2) = bnd_coord(1)
    pgBndData%bnd_coord(4) = bnd_coord(2)
    pgBndData%bnd_coord(6) = bnd_coord(3)

    pgBndData%domain_size(1) = pgBndData%bnd(2) - pgBndData%bnd(1)
    pgBndData%domain_size(2) = pgBndData%bnd(4) - pgBndData%bnd(3)
    pgBndData%domain_size(3) = pgBndData%bnd(6) - pgBndData%bnd(5)
  end if ! pgBndData%useBnd

  ! NB: I think boundary interaction code below is old stuff and can be removed
  nBCs = geometry%boundary%nBCtypes
  if ( allocated(particleGroup%BC_interaction) ) deallocate( particleGroup%BC_interaction )
  allocate( particleGroup%BC_interaction(1:nBCs) )


  ! write(logUnit(1),*) '---- INIT PARTICLE BOUNDARY INTERACTIONS ------'
  open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
  do iBC = 1, nBCs
    BCid = scheme%field(1)%bc(iBC)%bc_id
    ! write(pgDebugLog%lu,*) 'BCid = ', BCid, ' kind ', trim(scheme%field(1)%bc(iBC)%BC_kind)
    ! Determine whether this BC is a wall or not
    if(trim(scheme%field(1)%bc(iBC)%BC_kind) == 'wall') then
      particleGroup%BC_interaction( BCid ) = 0
    else
      particleGroup%BC_interaction( BCid ) = 1
      write(pgDebugLog%lu,*) '--- Particle BC treatment BC with ID: ', BCid, ' ---'
      write(pgDebugLog%lu,*) 'kind: ', scheme%field(1)%bc(iBC)%BC_kind
      write(pgDebugLog%lu,*) 'interaction set to: ', particleGroup%BC_interaction( BCid )
    end if
  end do ! iBC
  close( pgDebugLog%lu )

  ! ------ INITIALIZE COMMUNICATION BUFFERS AND DATATYPES------ !
  call mus_particles_initialize_communication( particleGroup = particleGroup, &
                                             & scheme = scheme,               &
                                             & geometry = geometry,           &
                                             & params = params                )

  ! ------ INITIALIZE PARTICLE INTERPOLATOR ------- !
  call init_particle_interpolator( interpolator = particleGroup%interpolator, &
                                 & stencil      = scheme%layout%stencil(1)    )

  ! ------ INITIALIZE PARTICLE REPRESENTATION ON GRID ------ !
  ! call printTotalElemList(scheme, geometry, lev, params%general%proc%rank, 90+params%general%proc%rank)
  ! call printNeighList(scheme, geometry, lev, params%general%proc%rank, 95+params%general%proc%rank)

  !-- Initialize the individual particles in this group --!
  iParticle = 1

  select case ( trim(params%particle_kind) )
    case ('MEM')
      open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
      write( pgDebugLog%lu, *) "Initial particles in group: ", &
        & particleGroup%particles_MEM%nvals
      close( pgDebugLog%lu )

    case('DPS','DPS_oneway','DPS_unittest')
      ! If we read particle domain boundary data from lua file, print it here.
      if( pgBndData%useBnd ) then
        open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
          write( pgDebugLog%lu, *) "Particle domain bounds defined as: "
          write( pgDebugLog%lu, *) "xmin = ", pgBndData%bnd(1)
          write( pgDebugLog%lu, *) "xmax = ", pgBndData%bnd(2)
          write( pgDebugLog%lu, *) "ymin = ", pgBndData%bnd(3)
          write( pgDebugLog%lu, *) "ymax = ", pgBndData%bnd(4)
          write( pgDebugLog%lu, *) "zmin = ", pgBndData%bnd(5)
          write( pgDebugLog%lu, *) "zmax = ", pgBndData%bnd(6)
          write( pgDebugLog%lu, *) "xmin_coord = ", pgBndData%bnd_coord(1)
          write( pgDebugLog%lu, *) "xmax_coord = ", pgBndData%bnd_coord(2)
          write( pgDebugLog%lu, *) "ymin_coord = ", pgBndData%bnd_coord(3)
          write( pgDebugLog%lu, *) "ymax_coord = ", pgBndData%bnd_coord(4)
          write( pgDebugLog%lu, *) "zmin_coord = ", pgBndData%bnd_coord(5)
          write( pgDebugLog%lu, *) "zmax_coord = ", pgBndData%bnd_coord(6)
          write( pgDebugLog%lu, *) "periodicBnd = ", pgBndData%periodicBnd
          write( pgDebugLog%lu, *) "wallBnd = ", pgBndData%wallBnd
        close( pgDebugLog%lu )
      end if

      open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
      write( pgDebugLog%lu, *) "Initial particles in group: ", &
        & particleGroup%particles_DPS%nvals
      close( pgDebugLog%lu )

  end select

  ! If doing unresolved particle simulations using the Generalized Navier Stokes (GNS)
  ! equations, update the initial fluid volume fraction field.
  if( trim(scheme%header%kind) == 'fluid_GNS' &
  & .OR. trim(scheme%header%kind) == 'fluid_incompressible_GNS') then
    call mus_particles_updateFluidVolumeFraction(                                &
      & particleGroup = particleGroup,                                           &
      & scheme        = scheme,                                                  &
      & geometry      = geometry,                                                &
      & params        = params,                                                  &
      & nElems        = scheme%pdf( geometry%tree%global%maxLevel )%nElems_local )
  end if ! scheme kind == fluid_GNS

end subroutine mus_particles_initialize

!> This routine moves all particles in particleArray.
!  Meant to be called once per time step.
subroutine mus_particles_move(particleGroup, scheme, geometry, params)
    !> Array of particles
    class(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !---------------------------------------------------!


    call DEMsubcycles_MEM( particleGroup = particleGroup,           &
                         & scheme        = scheme,                  &
                         & geometry      = geometry,                &
                         & params        = params,                  &
                         & Nsubcycles    = particleGroup%Nsubcycles )


end subroutine mus_particles_move

!> This routine moves all particles in particleArray.
!  Meant to be called once per time step.
subroutine mus_particles_move_DPS(particleGroup, scheme, geometry, params)
    !> Array of particles
    class(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !---------------------------------------------------!

    call DEMSubcycles_DPS( particleGroup = particleGroup,             &
                           & scheme        = scheme,                  &
                           & geometry      = geometry,                &
                           & params        = params,                  &
                           & Nsubcycles    = particleGroup%nSubcycles )

end subroutine mus_particles_move_DPS

!> This routine moves all particles in particleArray.
!  Meant to be called once per time step.
subroutine mus_particles_move_DPS_onewaycoupled(particleGroup, scheme, geometry, params)
    !> Array of particles
    class(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !---------------------------------------------------!

    call DEMSubcycles_DPS_onewaycoupled( particleGroup = particleGroup,           &
                                       & scheme        = scheme,                  &
                                       & geometry      = geometry,                &
                                       & params        = params,                  &
                                       & Nsubcycles    = particleGroup%Nsubcycles )

end subroutine mus_particles_move_DPS_onewaycoupled

!> mus_particles_mapping maps the current particle positions to the lattice
!! This means the exclusionLists are updated, connectivity of solid and
!! fluid neighbor particles is modified, new fluid particles are initialized
!! and have their connectivity restored. Also particles that have no more elements
!! on this process (either local or halo) get removed from this process.
!! Should be called once per LBM time step.
subroutine mus_particles_mapping_MEM(particleGroup, scheme, geometry, params)
    !> Array of particles
    class(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params

    !---------------------------------------------------!
    integer :: iParticle, lev
    real(kind=rk) :: dt, dx
    logical :: removeFromMyRank

    !---------------------------------------------------!
    lev = geometry%tree%global%maxLevel
    dt = params%physics%dtLvl(lev)
    dx = params%physics%dxLvl(lev)



    ! --- UPDATE PARTICLE REPRESENTATION ON GRID --- !
    iParticle = 1
    do
      if(iParticle > particleGroup%particles_MEM%nvals) exit
      removeFromMyRank = .FALSE.

      call mapToLattice( this = particleGroup%particles_MEM%val(iParticle), &
                       & particleGroup = particleGroup,                 &
                       & scheme = scheme,                               &
                       & stencil = scheme%layout%fstencil,              &
                       & geometry = geometry,                           &
                       & params = params,                               &
                       & rmflag = removeFromMyRank                      )


      if(removeFromMyRank) then
        ! write(logUnit(1),*) 'Removing particle pID', &
        !   & particleGroup%particles%val(iParticle)%particleID, ' from process ', &
        !   & params%general%proc%rank
        ! write(logUnit(1),*) 'iter = ', params%general%simcontrol%now%iter
        ! write(logUnit(1), '(A)', advance='no') 'coordOfOrigin = [ '
        ! write(logUnit(1), '(4I3)', advance = 'no') ( particleGroup%particles &
        !   &                                          %val(iParticle)      &
        !   &                                          %coordOfOrigin(i), i = 1,4)
        ! write(logUnit(1), '(A)') ' ]'

        call remove_particle_from_da_particle_MEM( particleGroup%particles_MEM, iParticle)
        ! Note: iParticle is NOT incremented after removing a particle from
        ! group because after removal the index iParticle corresponds to a
        ! different particle

      else
        iParticle = iParticle + 1
      end if
    end do

end subroutine mus_particles_mapping_MEM

subroutine mus_particles_mapping_DPS(particleGroup, scheme, geometry, params)
  !> Array of particles
  class(mus_particle_group_type), intent(inout) :: particleGroup
  !> Scheme for access to leveldescriptor
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  ! -----------------------------------!
  integer :: iParticle
  ! -----------------------------------!

  do iParticle = 1, particleGroup%particles_DPS%nvals
    call mapToLattice_DPS( particle = particleGroup%particles_DPS%val(iParticle),   &
                         & interpolator = particleGroup%interpolator, &
                         & scheme   = scheme,                                       &
                         & geometry = geometry,                                     &
                         & params   = params,                                       &
                         & comm     = particleGroup%send,                           &
                         & particleLogInterval = particleGroup%particleLogInterval  )
  end do

  ! -------------- Exchange new particles with neighboring processes --------------!
  call exchangeNewParticles_DPS( this = particleGroup,              &
                               & send = particleGroup%send,         &
                               & recv = particleGroup%recv,         &
                               & comm = MPI_COMM_WORLD,             &
                               & myRank = params%general%proc%rank, &
                               & message_flag = 1                   )

  ! Initialize the new particles I have received
  do iParticle = 1, particleGroup%particles_DPS%nvals
    if( particleGroup%particles_DPS%val(iParticle)%newForMe ) then


      call allocateProcessMasks(                                    &
        &    particle = particleGroup%particles_DPS%val(iParticle), &
        &    nProcs   = particleGroup%send%nProcs                   )

      call initParticle_DPS( &
                 & particle    = particleGroup%particles_DPS%val(iParticle),            &
                 & interpolator = particleGroup%interpolator,            &
                 & particleID  = particleGroup%particles_DPS%val(iParticle)%particleID, &
                 & geometry    = geometry,                                              &
                 & scheme      = scheme,                                                &
                 & myRank      = params%general%proc%rank,                              &
                 & comm        = particleGroup%send                                     )

      particleGroup%particles_DPS%val(iParticle)%newForMe = .FALSE.

      write(logUnit(1),*) 'mus_particles_mapping_DPS: initialized new particle on proc', &
        & params%general%proc%rank

      open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
        write(pgDebugLog%lu,*) 'Initializing new particle with ID ', &
          & particleGroup%particles_DPS%val(iParticle)%particleID, &
          & ' iter = ', params%general%simcontrol%now%iter
        write(pgDebugLog%lu,*) 'pos = ', particleGroup%particles_DPS%val(iParticle)%pos(1:6)
        write(pgDebugLog%lu,*) 'vel = ', particleGroup%particles_DPS%val(iParticle)%vel(1:6)
        write(pgDebugLog%lu,*) 'F = ', particleGroup%particles_DPS%val(iParticle)%F(1:6)
        write(pgDebugLog%lu,*) 'F_DEM(1,:) = ', particleGroup%particles_DPS%val(iParticle)%F_DEM(1,1:6)
        write(pgDebugLog%lu,*) 'F_DEM(2,:) = ', particleGroup%particles_DPS%val(iParticle)%F_DEM(2,1:6)
        write(pgDebugLog%lu,*) 'coordOfOrigin = ', particleGroup%particles_DPS%val(iParticle)%coordOfOrigin(1:4)
      close( pgDebugLog%lu )

    end if
  end do

  ! -------------- Exchange particles that have left the global domain --------------!
  call exchangeParticlesToRemove_DPS( this         = particleGroup,            &
                                    & send         = particleGroup%send,       &
                                    & recv         = particleGroup%recv,       &
                                    & comm         = MPI_COMM_WORLD,           &
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

end subroutine mus_particles_mapping_DPS



subroutine mus_particles_applyHydrodynamicForces_MEM( &
  &                                 particleGroup, scheme, geometry, params)
  !> Array of particles
  class(mus_particle_group_type), intent(inout) :: particleGroup
  !> Scheme for access to leveldescriptor
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params


  integer :: iParticle, lev
  ! integer :: particleLogUnit
  logical :: rmflag       ! .TRUE. if particle should be removed from group
  real(kind=rk) :: dt
  !---------------------------------------------------!
  lev = geometry%tree%global%maxLevel
  dt = params%physics%dtLvl(lev)

  iParticle = 1
  do
    if(iParticle > particleGroup%particles_MEM%nvals) exit
    ! Compute hydrodynamic force of fluid on particle and store the result
    ! in particle%Fbuff(particle%Fnow,:)
    call applyHydrodynamicForces(                                 &
            & this = particleGroup%particles_MEM%val(iParticle),  &
            & scheme = scheme,                                    &
            & stencil = scheme%layout%fStencil,                   &
            & params = params,                                    &
            & rho_p_lat = particleGroup%rho0_lat                  )

    ! call particleGroup%particles%val(iParticle)%applyParticleWallForces(     &
    ! &            BCinteraction = particleGroup%BC_interaction, &
    ! &            scheme        = scheme,                       &
    ! &            stencil       = scheme%layout%fStencil,       &
    ! &            geometry      = geometry,                     &
    ! &            params        = params,                       &
    ! &            rmflag        = rmflag                        )

    ! if(rmflag) then
    !   particleGroup%particles%val(iParticle)%removeParticle_global = .TRUE.
    ! else
    !   particleGroup%particles%val(iParticle)%removeParticle_global = .FALSE.
    ! end if

    iParticle = iParticle + 1

  end do

  ! --- EXCHANGE FORCE CONTRIBUTIONS WITH NEIGHBOR PROCESSES --- !
  call  particleGroup%exchangeForces( send = particleGroup%send,         &
                                    & recv = particleGroup%recv,         &
                                    & comm = MPI_COMM_WORLD,   &
                                    & myRank = params%general%proc%rank, &
                                    & message_flag = 1                   )

  ! NOTE 25-04-2022: I think a call to exchangeParticlesToRemove is missing here!

  ! --- REMOVE PARTICLES THAT HAVE LEFT THE GLOBAL DOMAIN --- !
  iParticle = 1
  do
    if(iParticle > particleGroup%particles_MEM%nvals) exit

    if( particleGroup%particles_MEM%val(iParticle)%removeParticle_global ) then
      open( pgDebugLog%lu, file=pgDebugLog%lfile, status='old', position='append' )
        write(pgDebugLog%lu,*) '--- Rank ', params%general%proc%rank, ' particle hit wall! ---'
        write(pgDebugLog%lu,*) 'Destroying particle with ID ', &
          & particleGroup%particles_MEM%val(iParticle)%particleID
      close( pgDebugLog%lu )
      call destroyParticle_MEM( this = particleGroup%particles_MEM%val(iParticle), &
                              & particleGroup = particleGroup,                 &
                              & scheme        = scheme,                        &
                              & params        = params                         )

      call remove_particle_from_da_particle_MEM( particleGroup%particles_MEM, iParticle )
      ! Note: iParticle is NOT incremented after removing a particle from
      ! group because after removal the index iParticle corresponds to a
      ! different particle
    else
      iParticle = iParticle + 1
    end if
  end do

  ! --- UPDATE VELOCITIES OF PARTICLES --- !
  do iParticle = 1, particleGroup%particles_MEM%nvals
    if(iParticle > particleGroup%particles_MEM%nvals) exit
    ! Average hydrodynamic forces over the past two time steps
    particleGroup%particles_MEM%val(iParticle)%F = &
      & 0.5 * ( particleGroup%particles_MEM%val(iParticle)%Fbuff(1,:) &
      & + particleGroup%particles_MEM%val(iParticle)%Fbuff(2,:) )

    call swapFBuff( particleGroup%particles_MEM%val(iParticle) )
  end do


  ! --- EXCHANGE INCOMING/OUTGOING PARTICLES WITH NEIGHBOR PROCESSES
  call  particleGroup%exchangeNewParticles( send = particleGroup%send,         &
                                          & recv = particleGroup%recv,         &
                                          & comm = MPI_COMM_WORLD,             &
                                          & myRank = params%general%proc%rank, &
                                          & message_flag = 1                   )

  ! Initialize all the new particles we might have received
  do iParticle = 1, particleGroup%particles_MEM%nvals
    if( particleGroup%particles_MEM%val(iParticle)%newForMe ) then

      ! write(logUnit(1),*) '--- INIT NEW PARTICLE ON PROC ', params%general%proc%rank,  '---'
      ! write(logUnit(1),*) 'Particle ID =', particleGroup%particles%val(iParticle)%particleID
      ! write(logUnit(1),*) 'iter = ', params%general%simcontrol%now%iter
      ! write(logUnit(1),*) '---------------------------------------------------------'

      ! Set the particle procedure pointers and initialize the communication masks
      call allocateProcessMasks(                                    &
        &    particle = particleGroup%particles_MEM%val(iParticle), &
        &    nProcs   = particleGroup%send%nProcs                   )

      ! Map particle to lattice using existing coordOfOrigin
      call initParticle_MEM( particle   = particleGroup%particles_MEM%val(iParticle), &
                           & particleID = iParticle,                              &
                           & geometry   = geometry,                               &
                           & scheme     = scheme,                                 &
                           & myRank     = params%general%proc%rank,               &
                           & comm       = particleGroup%send,                     &
                           & rmflag     = rmflag                                  )

    end if
  end do

end subroutine mus_particles_applyHydrodynamicForces_MEM

subroutine mus_particles_applyHydrodynamicForces_DPS( particleGroup, scheme, &
                                                    & geometry, params )
  !> Array of particles
  class(mus_particle_group_type), intent(inout) :: particleGroup
  !> Scheme for access to leveldescriptor
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  ! ------------------------------------------!
  integer :: iParticle
  ! ------------------------------------------!

  do iParticle = 1, particleGroup%particles_DPS%nvals
    if( particleGroup%particles_DPS%val(iParticle)%owner == params%general%proc%rank ) then
      ! First interpolate the fluid properties at the position of the particle
      call mus_particles_DPS_interpolateFluidProperties( particle = particleGroup      &
                                                       &            %particles_DPS     &
                                                       &            %val(iParticle),   &
                                                       & interpolator = particleGroup  &
                                                       &                %interpolator, &
                                                       & scheme   = scheme,            &
                                                       & geometry = geometry,          &
                                                       & params   = params, &
                                                       & intp = interpolateFluidProps, &
                                                       & calc_vel_and_p_grad = calcVelocityAndPressureGradient )

      ! Actual hydrodynamic forces are computed (using the values interpolated here)
      ! inside the DEM subcycling routine.

    end if ! particle%owner = myRank

  end do

end subroutine mus_particles_applyHydrodynamicForces_DPS

subroutine mus_particles_applyHydrodynamicForces_DPS_onewaycoupled( particleGroup, scheme, &
                                                    & geometry, params )
  !> Array of particles
  class(mus_particle_group_type), intent(inout) :: particleGroup
  !> Scheme for access to leveldescriptor
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  ! ------------------------------------------!
  integer :: iParticle
  ! ------------------------------------------!
  do iParticle = 1, particleGroup%particles_DPS%nvals
    if( particleGroup%particles_DPS%val(iParticle)%owner == params%general%proc%rank ) then
      ! First interpolate the fluid properties at the position of the particle
      call mus_particles_DPS_interpolateFluidProperties_onewaycoupled(       &
                                & particle     = particleGroup               &
                                &                %particles_DPS              &
                                &                %val(iParticle),            &
                                & interpolator = particleGroup%interpolator, &
                                & scheme       = scheme,                     &
                                & geometry     = geometry,                   &
                                & params       = params                      )

      ! Actual hydrodynamic forces are computed (using the values interpolated here)
      ! inside the DEM subcycling routine.

    end if ! particle%owner = myRank

  end do

end subroutine mus_particles_applyHydrodynamicForces_DPS_onewaycoupled


subroutine mus_particles_transferMomentumToFluid_MEM( &
    &                                 particleGroup, scheme, geometry, params)
    !> Array of particles
    class(mus_particle_group_type), intent(inout) :: particleGroup
    !> Scheme for access to leveldescriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !---------------------------------------------------!
    integer :: iParticle
    !---------------------------------------------------!

    do iParticle = 1, particleGroup%particles_MEM%nvals
        call applyVelocityBounceback( this = particleGroup%particles_MEM%val(iParticle), &
                                    & scheme = scheme,                               &
                                    & stencil = scheme%layout%fStencil,              &
                                    & params = params                                )
    end do

end subroutine mus_particles_transferMomentumToFluid_MEM

subroutine mus_particles_transferMomentumToFluid_DPS( particleGroup, scheme, &
                                                    & geometry, params       )
  !> Array of particles
  class(mus_particle_group_type), intent(inout) :: particleGroup
  !> Scheme for access to leveldescriptor
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  ! ------------------------------------------!
  integer :: iParticle
  real(kind=rk) :: F_particle(3), Fparticle_fluid_total(3)
  ! ------------------------------------------!
  Fparticle_fluid_total = 0.0_rk
  ! Before transferring momentum to fluid, we need to ensure the hydrodynamic force
  ! on the particles is up-to-date on each process:
  call exchangeHydroForces_DPS( this         = particleGroup,            &
                              & send         = particleGroup%send,       &
                              & recv         = particleGroup%recv,       &
                              & comm         = MPI_COMM_WORLD,           &
                              & myRank       = params%general%proc%rank, &
                              & message_flag = 1                         )

  do iParticle = 1, particleGroup%particles_DPS%nvals
    F_particle = 0.0_rk
    call particleGroup%transfer_momentum( &
                  & particle     = particleGroup%particles_DPS%val(iParticle), &
                  & interpolator = particleGroup%interpolator,                 &
                  & scheme       = scheme,                                     &
                  & geometry     = geometry,                                   &
                  & params       = params,                                     &
                  & Ftot         = F_particle                                  )
  end do



end subroutine mus_particles_transferMomentumToFluid_DPS

subroutine mus_particles_addSourceTermsToAuxField_DPS( particleGroup, scheme, &
                                                    & geometry, params       )
  !> Array of particles
  class(mus_particle_group_type), intent(inout) :: particleGroup
  !> Scheme for access to leveldescriptor
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  ! ------------------------------------------!
  integer :: iParticle
  real(kind=rk) :: auxFieldForceSum(3), tmp(3)
  ! ------------------------------------------!
  tmp = 0.0_rk
  auxFieldForceSum = 0.0_rk

  do iParticle = 1, particleGroup%particles_DPS%nvals
    call particleGroup%modify_auxfield(                                  &
            & particle     = particleGroup%particles_DPS%val(iParticle), &
            & interpolator = particleGroup%interpolator,                 &
            & scheme       = scheme,                                     &
            & geometry     = geometry,                                   &
            & params       = params,                                     &
            & Ftot         = tmp                                         )

    auxFieldForceSum = auxFieldForceSum + tmp
  end do

end subroutine mus_particles_addSourceTermsToAuxField_DPS

! ---------------- ROUTINES FOR DEBUGGING/PRINTING STUFF ------------------ !
!> Routine for manually checking connectivity of particle elements,
!! to be removed later
subroutine testParticleConnectivity(scheme, particle, lev)
  type(mus_scheme_type) :: scheme
  type(mus_particle_MEM_type) :: particle
  integer :: lev


  integer :: nElems, nStateElems, nSize, nScalars
  integer :: iDir, iField, iElem, idx_idir, QQ
  integer :: nghElem, nghDir
  integer :: stateVectorPos, neighPos, ldNeighPos
  integer :: elemPos

  nElems = particle%exclusionList%nvals
  nStateElems = scheme%pdf( lev )%nElems_local
  nScalars = scheme%varSys%nScalars

  nSize  = scheme%pdf( lev )%nSize
  QQ = scheme%layout%fStencil%QQ

  do iElem = 1, 1! nElems
  elemPos = particle%exclusionList%Val(iElem)
    do iDir = 1,QQ
      do iField = 1, scheme%nFields
        idx_idir = scheme%varSys%method%val(                      &
          &               scheme%stateVarMap%varPos%val(iField) ) &
          &                     %state_varPos(iDir)



        ! Position in state vector of this direction and this element
        stateVectorPos = ( elempos-1)* nscalars+idx_idir

        ! Position of neighbor PDF to stream from in state vector.
        ! Should be identical to stateVectorPos for elems belonging to particle

        ! NB: Compare this to the other calls to pdf%neigh!
        !     nStateElems should be nSize!!!
        neighPos = scheme%pdf(lev)%neigh( ( idir-1)* nsize + elempos)
        !neighPos = scheme%pdf(lev)%neigh( ( idir-1)* nstateelems + elempos)

        write(logUnit(1),'(A)') ' '
        write(logUnit(1), '(A,I12)') 'idir = ', iDir
        write(logUnit(1), '(A,I12)') 'idx_idir = ', idx_idir
        write(logUnit(1), '(A,I12)') 'elemPos = ', elemPos
        write(logUnit(1), '(A,I12)') 'stateVectorPos = ', stateVectorPos
        write(logUnit(1), '(A,I12)') 'PDF neigh = ', neighPos

        ! ldNeighElem = scheme%levelDesc( lev )                                   &
        !   &              %neigh(1)                                           &
        !   &              %nghElems(  scheme%layout%fstencil %cxdirinv( idir), &
        !   &                          elemPos                                 )

        ! Now compare to neighbor index found from level descriptor.
        ! These should be identical before invoking particle connectivity
        ! routines

        ! Direction to stream from
        nghDir =  scheme%layout%fstencil%cxdirinv(idx_idir)
        ! Neighboring element in that direction
        nghElem = scheme%levelDesc(lev)%neigh(1)%nghElems(nghDir, elemPos)

        ! Position of link to stream from
        ldNeighPos = ( nghelem-1)* nscalars+ idx_idir

        write(logUnit(1), '(A,I12)') 'levelDesc neigh = ', ldNeighPos

      end do
    end do
  end do
end subroutine testParticleConnectivity


subroutine printTotalElemList(scheme, geometry, lev, proc, plogUnit)
  type(mus_scheme_type), intent(in) :: scheme
  type(mus_geom_type), intent(in) :: geometry
  integer, intent(in) :: lev
  integer, intent(in) :: proc
  integer, intent(in) :: plogUnit
  ! --------------------------------------------!
  integer :: iElem, posInBnd
  integer(kind=long_k) :: bcID
  integer :: iDir
  integer(kind=long_k) :: elemTreeID
  integer(kind=long_k) :: elemProp
  ! --------------------------------------------!

  character(len=1024) :: filename

  write (filename, "(A9,I0.4)") "totallist", proc

  filename = trim(filename)//'.dat'

  open(plogUnit, file = filename, status = 'new')

  write(plogUnit, *) 'Local fluid TreeIDs'
  write(plogUnit,'(A12)', advance='no') 'TreeID'
  write(plogUnit,'(A10)', advance='no') 'hasBnd?'
  write(plogUnit,'(A)') 'Boundary IDs in each of the stencil directions'

  do iElem = 1, scheme%pdf(lev)%nElems_fluid
    elemTreeID = scheme%levelDesc(lev)%total(iElem)
    elemProp = scheme%levelDesc(lev)%property(iElem)
    if( btest(elemProp, prp_hasBnd) ) then
      write(plogUnit, '(I12)', advance='no') elemTreeID
      write(plogUnit, '(A10)', advance='no') ' hasBnd'
      ! Find out what the BCid is
      posInBnd = geometry%posInBndID(iElem)
      do iDir = 1,scheme%layout%fstencil%QQN
        bcID = geometry%boundary%boundary_ID(iDir,posInBnd)
        write(plogUnit, '(I12)', advance='no') bcID
      end do
      write(plogUnit,'(A)') ' '
    else
      write(plogUnit, '(I12)', advance='no') elemTreeID
      write(plogUnit, '(A10)') ' '
    end if
  end do

  write(plogUnit, *) 'Halo TreeIDs'
  do iElem = scheme%pdf(lev)%nElems_fluid + 1, scheme%pdf(lev)%nElems_fluid + scheme%pdf(lev)%nElems_halo
    elemTreeID = scheme%levelDesc(lev)%total(iElem)
    elemProp = scheme%levelDesc(lev)%property(iElem)
    if( btest(elemProp, prp_hasBnd) ) then
      write(plogUnit, '(I12)', advance='no') elemTreeID
      write(plogUnit, '(A)') 'hasBnd'
    else
      write(plogUnit, '(I12)') elemTreeID
    end if
  end do
  ! do iElem = 1, ubound( scheme%levelDesc(lev)%total, 1 )
  !   elemTreeID = scheme%levelDesc(lev)%total(iElem)
  !   write(logUnit, '(I12)') elemTreeID
  ! end do

  close(unit = plogUnit, status='KEEP')
end subroutine printTotalElemList

subroutine printNeighList(scheme, geometry, lev, proc, plogUnit)
  type(mus_scheme_type), intent(in) :: scheme
  type(mus_geom_type), intent(in) :: geometry
  integer, intent(in) :: lev
  integer, intent(in) :: proc
  integer, intent(in) :: plogUnit
  ! --------------------------------------------!
  integer :: iElem, iDir, iBnd
  integer :: QQN
  integer :: neighPos
  integer(kind=long_k) :: elemTreeID

  ! --------------------------------------------!

  character(len=1024) :: filename
  QQN = scheme%layout%fStencil%QQN

  write (filename, "(A9,I0.4)") "neighlist", proc

  filename = trim(filename)//'.dat'

  open(plogUnit, file = filename, status = 'new')

  write(plogUnit, *) 'Local fluid TreeIDs'
  do iElem = 1, scheme%pdf(lev)%nElems_fluid
    elemTreeID = scheme%levelDesc(lev)%total(iElem)
    write(plogUnit, '(I8)', advance='no') elemTreeID
    do iDir = 1, QQN
      neighPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, iElem)
      if(neighPos > 0) then
        elemTreeID = scheme%levelDesc(lev)%total(neighPos)
        write(plogUnit, '(I8)', advance='no') elemTreeID
      else
        write(plogUnit, '(I8)', advance='no') neighPos
      end if
    end do
    write(plogUnit, '(A)') ' '
  end do

  write(plogUnit, *) 'Halo TreeIDs'
  do iElem = scheme%pdf(lev)%nElems_fluid + 1, scheme%pdf(lev)%nElems_fluid + scheme%pdf(lev)%nElems_halo
    elemTreeID = scheme%levelDesc(lev)%total(iElem)
    write(plogUnit, '(I8)', advance='no') elemTreeID
    do iDir = 1, QQN
      neighPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, iElem)
      if(neighPos > 0) then
        elemTreeID = scheme%levelDesc(lev)%total(neighPos)
        write(plogUnit, '(I8)', advance='no') elemTreeID
      else
        write(plogUnit, '(I8)', advance='no') neighPos
      end if

    end do
    write(plogUnit, '(A)') ' '
  end do

  write(plogUnit,*) 'This mesh has ', geometry%boundary%nBCtypes,  ' boundary labels'
  do iBnd = 1, geometry%boundary%nBCtypes
    write(plogUnit,*) iBnd, ': ', geometry%boundary%BC_label(iBnd)
  end do

  close(unit = plogUnit, status='KEEP')
end subroutine printNeighList

! This is to test whether applyHydrodynamicForces indeed only loops over local fluid links
subroutine test_loopOverLocalLinks( this, scheme, stencil, params )
  type(mus_particle_MEM_type), intent(inout) :: this
  !> Scheme for access to level descriptor
  type(mus_scheme_type), intent(inout) :: scheme
  !> fluid stencil
  type(tem_stencilHeader_type), intent(in) :: stencil
  !> Parameters for access to conversion factors between units
  type(mus_param_type), intent(in) :: params
  ! ------------------------------------------!

  integer :: lev
  integer :: elemPos, neighPos
  integer :: iDir, iField
  integer(kind=long_k) :: neighProp
  integer :: iElem, nLocalElems
  integer :: nLinks


  ! ------------------------------------------!
  lev = this%coordOfOrigin(4)
  iField = 1

  nLocalElems = scheme%pdf( lev )%nElems_fluid


  ! Loop over the elements belonging to this particle (those in exclusionList)
  ! Then identify fluid neighbors from there.
  nLinks = 0
  particleElemLoop: do iElem = 1, this%exclusionList%nVals
    ! position of this node in levelDesc state vector
    elemPos = this%exclusionList%val(iElem)

    ! Only compute force contribution from links connecting local fluid elems
    ! not halos!
    if(elemPos > nLocalElems ) then
      cycle
    end if

    linkLoop: do iDir = 1, stencil%QQN
      ! Get the index of the neighbor element in this direction.
      neighPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, elemPos)
      ! Initialize neighbor's property
      neighProp = 0_long_k

      ! If neighbor element exists, get its property
      if( neighPos > 0 ) neighProp = scheme%levelDesc(lev)%property(neighPos)

      if ( neighPos <= 0 .or. neighPos > nLocalElems .or. (btest(neighProp, prp_solid)) ) then
        ! this link either points to a BC, a halo or another solid element so skip it
        cycle
      else

      nLinks = nLinks + 1

      end if
    end do linkLoop
  end do particleElemLoop

  write(logUnit(1), *) 'Process ', params%general%proc%rank, ' force comp loops over ', nLinks, ' links'

end subroutine test_loopOverLocalLinks

end module mus_particle_module
