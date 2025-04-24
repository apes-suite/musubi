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
!> mus_particle_tests_module contains a number of tests for routines of the
!! particle Musubi code. Right now this is a work-in-progress.

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

module mus_particle_tests_module

use mpi

use env_module,                        only : rk, long_k, newUnit
use tem_param_module,                  only : div1_3, div1_6, PI, cs2inv, cs4inv
use tem_aux_module,                    only : tem_abort
use tem_logging_module,                only : logUnit
use tem_dyn_array_module,              only : init, append, destroy,             &
                                            & empty, dyn_intArray_type,          &
                                            & dyn_longArray_type
use tem_geometry_module,               only : tem_CoordOfReal, tem_posOfId, tem_BaryOfId
use tem_topology_module,               only : tem_IdOfCoord, tem_FirstIdAtLevel, &
                                            & tem_coordOfId
use tem_varSys_module,                 only : tem_varSys_type
use mus_geom_module,                   only : mus_geom_type
use mus_scheme_type_module,            only : mus_scheme_type
use mus_param_module,                  only : mus_param_type
use mus_auxField_module,               only : mus_auxFieldVar_type
use mus_derVarPos_module,              only : mus_derVarPos_type
use mus_particle_type_module,          only : mus_particle_DPS_type,        &
                                            & mus_particle_group_type,      &
                                            & allocateProcessMasks,         &
                                            & append_da_particle_DPS, &
                                            & init_da_particle_DPS
use mus_particle_comm_type_module,     only : mus_particles_communication_type
use mus_particle_comm_module,          only : exchangeNewParticles_DPS
use mus_particle_aux_module,           only : getBaryOfCoord,  &
                                            & cross_product,   &
                                            & mus_particles_global_errorcheck, &
                                            & findPartitionOfTreeID
use mus_particle_logging_module,       only : getParticleLogUnit, &
                                            & logParticleData, &
                                            & closeParticleLog, &
                                            & generateElemListLine, &
                                            & openLogFile
use mus_particle_logging_type_module,  only : pgDebugLog
use mus_particle_DPS_module,           only : grabValueAtCoord, &
                                            & getInterpolationBnds,  &
                                            & getIndexOfCoordInTotal, &
                                            & interpolateFluidProps, &
                                            & calcVelocityAndPressureGradient
use mus_particle_boundary_module,      only : mus_particle_boundaryData_type, &
                                            & computeDisplacement
use mus_particle_DEM_module,           only : computeWallForce_DPS, &
                                            & DEM_fillNeighborList
use mus_particle_interpolation_module, only : mus_particle_interpolator_type, &
                                            & intp1d_linear, &
                                            & intp_1D_peskin, &
                                            & intp_1D_delta, &
                                            & get_xd, &
                                            & getWght1D_linear
use mus_particle_checks_module,        only : compute_fluid_momentum, &
                                            & compute_particle_momentum


implicit none

! -- Procedures and subroutines -- !
contains


subroutine mus_particles_runtests( particleGroup, scheme, geometry, params )
  type(mus_particle_group_type), intent(inout) :: particleGroup
  type(mus_scheme_type), intent(inout) :: scheme
  type(mus_geom_type), intent(in) :: geometry
  type(mus_param_type), intent(in) :: params
  ! ------------------------------------------!
  integer :: plogunit
  integer :: myRank, nProcs
  character(len=1024) :: fileName
  ! real(kind=rk) :: xp(3)
  ! integer :: iParticle
  ! logical :: flag_computeDisplacement, flag_computeWallForces
  ! logical :: flag_fillNeighborlist, flag_peskin, flag_compute_fluid_momentum
  ! ------------------------------------------!
  myRank = params%general%proc%rank
  nProcs = params%general%proc%comm_size
  write(fileName,'(A11,I1,A4)') 'testresults',myRank,'.txt'

  call openLogFile(fileName, plogunit)

  call test_global_errorcheck(myRank, plogunit)

  ! Only root process opens and writes to the log file
  ! if(myRank == 0) then
  !   call openLogFile( fileName, plogunit )
  ! end if

  ! xp = (/ 0.51_rk, 0.47_rk, 0.000001_rk /)

  ! call test_compute_fluid_momentum( scheme = scheme, &
  !                                 & geometry = geometry, &
  !                                 & params = params,&
  !                                 & comm = MPI_COMM_WORLD, &
  !                                 & myRank = myRank, &
  !                                 & nProcs = nProcs, &
  !                                 & plogunit = plogunit, &
  !                                 & flag = flag_compute_fluid_momentum)

  ! call test_compute_particle_momentum( lev = geometry%tree%global%maxLevel, &
  !                                    & params = params,&
  !                                    & comm = MPI_COMM_WORLD, &
  !                                    & myRank = myRank, &
  !                                    & nProcs = nProcs, &
  !                                    & plogunit = plogunit, &
  !                                    & flag = flag_compute_fluid_momentum)

  ! do iParticle = 1, particleGroup%particles_DPS%nvals, 1

    ! call test_distribution_DPS( particle = particleGroup%particles_DPS%val(iParticle), &
    !                            & scheme   = scheme,                                     &
    !                            & geometry = geometry,                                   &
    !                            & params   = params,                                     &
    !                            & plogunit  = plogunit                                     )

    ! call test_interpolation_delta_DPS_d2q9( xp           = xp,                         &
    !                                  & interpolator = particleGroup%interpolator, &
    !                                  & scheme       = scheme,                     &
    !                                  & geometry     = geometry,                   &
    !                                  & params       = params,                     &
    !                                  & plogunit      = plogunit                     )

  ! end do
  ! call test_computeDisplacement(flag_computeDisplacement)
  ! write( plogunit, * ) "test_computeDisplacement: ", flag_computeDisplacement

  ! call test_computeWallForces(flag_computeWallForces)
  ! write( plogunit, * ) "test_computeWallForces: ", flag_computeWallForces

  ! call test_DEM_fillNeighborList( flag_fillNeighborList )

  ! call test_calcVelocityAndPressuregradient_DPS( xp = xp, &
  !                                              & scheme = scheme, &
  !                                              & geometry = geometry, &
  !                                              & params = params, &
  !                                              & plogunit = plogunit )

  ! call test_generateElemListLine( scheme = scheme, &
  !                               & geometry = geometry, &
  !                               & params = params )

  ! call test_calcVelocityAndPressuregradient_DPS( xp = xp, &
  !                                              & scheme = scheme, &
  !                                              & geometry = geometry, &
  !                                              & params = params, &
  !                                              & plogunit = plogunit )

  close(plogunit)

end subroutine mus_particles_runtests


subroutine test_interpolation_delta_DPS_d3q19(xp, interpolator, scheme, geometry, params, plogunit)
  !> Query point (x,y,z) to interpolate fluid properties to
  real(kind=rk), intent(in) :: xp(3)
  !> Object containing interpolation information
  type(mus_particle_interpolator_type) :: interpolator
  !> Scheme for access to fluid data
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  !> Unit to write output to
  integer, intent(in) :: plogunit
  ! ----------------------------------------!
  integer :: lev
  integer :: vel_pos(3), dens_pos, vol_frac_pos
  integer :: elemOff
  integer :: iElem
  integer :: coord_xp(4)
  integer(kind=long_k) :: TreeID
  real(kind=rk) :: baryOfCoord(3), dx
  real(kind=rk) :: vel_tmp(3), rho_tmp, eps_f_tmp
  real(kind=rk) :: x, y, z
  ! ----------------------------------------!

  lev = geometry%tree%global%maxLevel
  dx = params%physics%dxLvl(lev)
  coord_xp = tem_coordOfReal( mesh  = geometry%tree, &
                            & point = xp,            &
                            & level = lev            )

  ! Position of velocity variable in auxField
  vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)
  dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
  vol_frac_pos = scheme%varSys%method%val(scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)


  ! Set the auxField quantities to some known value
  ! NOTE: known value should be prescribed in lattice units!
  do iElem = 1, scheme%pdf(lev)%nElems_local
    elemOff = (iElem-1)*scheme%varSys%nAuxScalars
    TreeID  = scheme%levelDesc(lev)%total(iElem)

    baryOfCoord = tem_baryOfId( tree   = geometry%tree, &
                              & TreeID = TreeID         )
    x = baryOfCoord(1)
    y = baryOfCoord(2)
    z = baryOfCoord(3)

    ! Add known (e.g. sine) function here
    scheme%auxField(lev)%val(elemOff+vol_frac_pos) = x + 2*y + 3*z
    scheme%auxField(lev)%val(elemOff+dens_pos)     = 1.0_rk
    scheme%auxField(lev)%val(elemOff+vel_pos(1))   = z
    scheme%auxField(lev)%val(elemOff+vel_pos(2))   = x
    scheme%auxField(lev)%val(elemOff+vel_pos(3))   = y

  end do

  ! Now interpolate the fluid props
  call interpolateFluidProps( xp           = xp, &
                            & interpolator = interpolator, &
                            & coord_xp     = coord_xp, &
                            & scheme       = scheme, &
                            & geom_origin  = geometry%tree%global%origin, &
                            & dx           = dx, &
                            & vel_xp       = vel_tmp,    &
                            & rho_xp       = rho_tmp,    &
                            & eps_f_xp     = eps_f_tmp )


  write(plogunit,'(A)') 'Interpolated values: '
  write(plogunit,'(A,E17.9)') 'eps_intp = ', eps_f_tmp
  write(plogunit,'(A,E17.9)') 'rho_intp = ', rho_tmp
  write(plogunit,'(A,E17.9)') 'ux_intp   = ', vel_tmp(1)
  write(plogunit,'(A,E17.9)') 'uy_intp   = ', vel_tmp(2)
  write(plogunit,'(A,E17.9)') 'uz_intp   = ', vel_tmp(3)

  write(plogunit,'(A)') 'Actual values: '
  write(plogunit,'(A,E17.9)') 'eps_intp = ', xp(1) + 2*xp(2) + 3*xp(3)
  write(plogunit,'(A,E17.9)') 'rho_intp = ', 1.0_rk
  write(plogunit,'(A,E17.9)') 'ux_intp   = ', xp(3)
  write(plogunit,'(A,E17.9)') 'uy_intp   = ', xp(1)
  write(plogunit,'(A,E17.9)') 'uz_intp   = ', xp(2)

end subroutine test_interpolation_delta_DPS_d3q19

subroutine test_interpolation_delta_DPS_d2q9(xp, interpolator, scheme, geometry, params, plogunit)
  !> Query point (x,y,z) to interpolate fluid properties to
  real(kind=rk), intent(in) :: xp(3)
  !> Object containing interpolation information
  type(mus_particle_interpolator_type) :: interpolator
  !> Scheme for access to fluid data
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  !> Unit to write output to
  integer, intent(in) :: plogunit
  ! ----------------------------------------!
  integer :: lev
  integer :: vel_pos(3), dens_pos, vol_frac_pos
  integer :: elemOff
  integer :: iElem
  integer :: coord_xp(4)
  integer(kind=long_k) :: TreeID
  real(kind=rk) :: baryOfCoord(3), dx
  real(kind=rk) :: vel_tmp(3), rho_tmp, eps_f_tmp
  real(kind=rk) :: x, y, z
  ! ----------------------------------------!

  lev = geometry%tree%global%maxLevel
  dx = params%physics%dxLvl(lev)
  coord_xp = tem_coordOfReal( mesh  = geometry%tree, &
                            & point = xp,            &
                            & level = lev            )

  ! Position of velocity variable in auxField
  vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)
  dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
  vol_frac_pos = scheme%varSys%method%val(scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)


  ! Set the auxField quantities to some known value
  ! NOTE: known value should be prescribed in lattice units!
  do iElem = 1, scheme%pdf(lev)%nElems_local
    elemOff = (iElem-1)*scheme%varSys%nAuxScalars
    TreeID  = scheme%levelDesc(lev)%total(iElem)

    baryOfCoord = tem_baryOfId( tree   = geometry%tree, &
                              & TreeID = TreeID         )
    x = baryOfCoord(1)
    y = baryOfCoord(2)
    z = baryOfCoord(3)

    ! Add known (e.g. sine) function here
    scheme%auxField(lev)%val(elemOff+vol_frac_pos) = x + 2*y
    scheme%auxField(lev)%val(elemOff+dens_pos)     = 1.0_rk
    scheme%auxField(lev)%val(elemOff+vel_pos(1))   = x + y
    scheme%auxField(lev)%val(elemOff+vel_pos(2))   = x
    scheme%auxField(lev)%val(elemOff+vel_pos(3))   = y

  end do

  ! Now interpolate the fluid props
  call interpolateFluidProps( xp           = xp, &
                            & interpolator = interpolator, &
                            & coord_xp     = coord_xp, &
                            & scheme       = scheme, &
                            & geom_origin  = geometry%tree%global%origin, &
                            & dx           = dx, &
                            & vel_xp       = vel_tmp,    &
                            & rho_xp       = rho_tmp,    &
                            & eps_f_xp     = eps_f_tmp )


  write(plogunit,'(A)') 'Interpolated values: '
  write(plogunit,'(A,E17.9)') 'eps_intp = ', eps_f_tmp
  write(plogunit,'(A,E17.9)') 'rho_intp = ', rho_tmp
  write(plogunit,'(A,E17.9)') 'ux_intp   = ', vel_tmp(1)
  write(plogunit,'(A,E17.9)') 'uy_intp   = ', vel_tmp(2)
  write(plogunit,'(A,E17.9)') 'uz_intp   = ', vel_tmp(3)

  write(plogunit,'(A)') 'Actual values: '
  write(plogunit,'(A,E17.9)') 'eps_intp = ', xp(1) + 2*xp(2)
  write(plogunit,'(A,E17.9)') 'rho_intp = ', 1.0_rk
  write(plogunit,'(A,E17.9)') 'ux_intp   = ', xp(1) + xp(2)
  write(plogunit,'(A,E17.9)') 'uy_intp   = ', xp(1)
  write(plogunit,'(A,E17.9)') 'uz_intp   = ', xp(2)

end subroutine test_interpolation_delta_DPS_d2q9


subroutine test_calcVelocityAndPressureGradient_DPS(xp, scheme, geometry, params, plogunit)
  !> Query point (x,y,z) to interpolate fluid properties to
  real(kind=rk), intent(in) :: xp(3)
  !> Scheme for access to fluid data
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  !> Unit to write output to
  integer, intent(in) :: plogunit
  ! ----------------------------------------!
  integer :: lev
  integer :: vel_pos(3), dens_pos, vol_frac_pos
  integer :: elemOff
  integer :: iElem
  integer :: coord_xp(4)
  integer(kind=long_k) :: TreeID
  real(kind=rk) :: baryOfCoord(3), dx
  real(kind=rk) :: curl_u_tmp(3), grad_p_tmp(3)
  real(kind=rk) :: x, y, z
  real(kind=rk) :: dt
  logical :: err
  ! ----------------------------------------!
  err = .FALSE.

  lev = geometry%tree%global%maxLevel
  dx  = params%physics%dxLvl(lev)
  dt  = params%physics%dtLvl(lev)
  coord_xp = tem_coordOfReal( mesh  = geometry%tree, &
                            & point = xp,            &
                            & level = lev            )

  ! Position of variables in auxField
  vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)
  dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
  vol_frac_pos = scheme%varSys%method%val(scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)


  ! Set the auxField quantities to some known value
  ! NOTE: known value should be prescribed in lattice units!
  do iElem = 1, scheme%pdf(lev)%nElems_local
    elemOff = (iElem-1)*scheme%varSys%nAuxScalars
    TreeID  = scheme%levelDesc(lev)%total(iElem)

    baryOfCoord = tem_baryOfId( tree   = geometry%tree, &
                              & TreeID = TreeID         )

    ! coordinates of our cell (in physical units!)
    x = baryOfCoord(1)
    y = baryOfCoord(2)
    z = baryOfCoord(3)

    ! Add known (e.g. sine) function here
    scheme%auxField(lev)%val(elemOff+vol_frac_pos) = 1.0_rk
    ! scheme%auxField(lev)%val(elemOff+dens_pos)     = x*params%physics%rho0
    ! scheme%auxField(lev)%val(elemOff+vel_pos(1))   = -y*(dt/dx)
    ! scheme%auxField(lev)%val(elemOff+vel_pos(2))   = x*(dt/dx)
    ! scheme%auxField(lev)%val(elemOff+vel_pos(3))   = 0.0_rk

    ! Base our prescribed density and velocity fields off of position in LATTICE units
    ! Note p = cs2*rho -> to prescribe p = x, rho = x/cs2inv
    scheme%auxField(lev)%val(elemOff+dens_pos)     = (cs2inv*x)/dx
    scheme%auxField(lev)%val(elemOff+vel_pos(1))   = -y/dx
    scheme%auxField(lev)%val(elemOff+vel_pos(2))   = x/dx
    scheme%auxField(lev)%val(elemOff+vel_pos(3))   = 0.0_rk

  end do

  ! Calculate the curl of velocity and pressure gradient
  call calcVelocityAndPressureGradient( coord  = coord_xp,   &
                                      & scheme = scheme,     &
                                      & grad_p = grad_p_tmp, &
                                      & curl_u = curl_u_tmp, &
                                      & err    = err         )

  if(err) then
    write(logUnit(1),*) "Test calcVelocityAndPressureGradient FAILED, aborting"
    call tem_abort()
  end if

  ! Compare the calculate curl to the analytical result
  ! grad_p_tmp = grad_p_tmp*params%physics%fac(lev)%press*(1.0/dx)
  ! curl_u_tmp = curl_u_tmp*dt

  write(plogunit,*) "test_calcVelocityAndPressureGradient_DPS"
  write(plogunit,*) "grad p = "
  write(plogunit,*) grad_p_tmp
  write(plogunit,*) "curl u = "
  write(plogunit,*) curl_u_tmp

  write(plogunit,*) "analytic grad p = "
  write(plogunit,*) 1.0_rk
  write(plogunit,*) "analytic curl u = "
  write(plogunit,*) 2.0_rk

end subroutine test_calcVelocityAndPressureGradient_DPS
! ------------- TESTS FOR mus_particle_boundary_module ---------------- !
subroutine test_computeDisplacement(flag)
  logical, intent(out) :: flag
  ! ---------------------------- !
  type( mus_particle_boundarydata_type)  :: boundaryData
  real(kind=rk) :: x1(3)
  real(kind=rk) :: x2(3)
  real(kind=rk) :: r12(3)
  real(kind=rk) :: r12_exact(3)
  !> tolerance
  real(kind=rk) :: tol  = 1.0e-5
  ! ---------------------------- !
  ! Initialize the boundary data object
  boundaryData%bnd = (/ 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk /)
  boundaryData%domain_size(1) = boundaryData%bnd(2)-boundaryData%bnd(1)
  boundaryData%domain_size(2) = boundaryData%bnd(4)-boundaryData%bnd(3)
  boundaryData%domain_size(3) = boundaryData%bnd(6)-boundaryData%bnd(5)
  boundaryData%periodicBnd(1:6) = .TRUE.

  ! Set sample positions x1 and x2
  x1 = (/ 0.2_rk, 0.2_rk, 0.8_rk /)
  x2 = (/ 0.3_rk, 0.1_rk, 0.1_rk /)

  ! Set exact answer
  r12_exact(1) = 0.1_rk
  r12_exact(2) = -0.1_rk
  r12_exact(3) = 0.3_rk

  ! Compute the displacement
  r12 = computeDisplacement( x1, x2, boundaryData )

  ! See if it matches exact answer
  if( all( abs(r12 - r12_exact) < tol ) ) then
    flag = .TRUE.
  else
    flag = .FALSE.
    write(logUnit(1),*) "test_computeDisplacement FAILED"
    write(logUnit(1),*) "x1 = ", x1
    write(logUnit(1),*) "x2 = ", x2
    write(logUnit(1),*) "calculated displacement = ", r12
    write(logUnit(1),*) "Error = ", abs(r12 - r12_exact)

  end if


end subroutine test_computeDisplacement

subroutine test_computeWallForces(flag)
  logical, intent(out) :: flag
  ! ---------------------------- !
  type(mus_particle_DPS_type) :: particle
  !> Data for boundaries
  type( mus_particle_boundarydata_type)  :: boundaryData
  !> Particle position
  real(kind=rk) :: eps, tol
  real(kind=rk) :: Fcoll(3), Fcoll_exact(3)
  real(kind=rk) :: kn, dn
  ! ---------------------------- !
  Fcoll = 0.0_rk
  tol   = 1.0e-5

  ! Initialize the boundary data object
  boundaryData%bnd = (/ 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk /)
  boundaryData%domain_size(1) = boundaryData%bnd(2)-boundaryData%bnd(1)
  boundaryData%domain_size(2) = boundaryData%bnd(4)-boundaryData%bnd(3)
  boundaryData%domain_size(3) = boundaryData%bnd(6)-boundaryData%bnd(5)
  boundaryData%periodicBnd(1:6) = .FALSE.

  ! --- TEST #1 --- !
  ! Initialize mock particle object
  particle%pos(1:3) = (/ 0.95_rk, 0.05_rk, 0.5_rk /)
  particle%vel(1:3) = (/ -0.1_rk, 0.1_rk, 0.0_rk /)
  particle%radius = 0.05_rk
  particle%mass = 1.0_rk

  ! Collision parameters
  eps = 0.1_rk
  ! Tc  = PI
  ! e_dry = 1.0_rk
  ! kn = this%mass*( PI**2 + (log(e_dry) )**2 )/Tc**2
  ! dn = (-2.0*this%mass*log(e_dry))/Tc

  ! ---- TEST #1: no damping
  ! kn = 1.0
  ! dn = 0.0

  ! Fcoll = computeWallForce_DPS( this         = particle,     &
      !                           & boundaryData = boundaryData, &
      !                           & eps          = eps,          &
      !                           & kn           = kn,           &
      !                           & dn           = dn            )

  ! ! Overlap should be 0.1 in x-direction. So Fcoll should be -0.1
  ! ! See if it matches exact answer
  ! Fcoll_exact = (/ -0.1_rk, 0.1_rk, 0.0_rk /)

  ! if( all( abs( Fcoll - Fcoll_exact) < tol )  ) then
  !   flag = .TRUE.
  !   write(logUnit(1),*) "test_computeWallForces SUCCEEDED"
  !   write(logUnit(1),*) "calculated force = ", Fcoll
  ! else
  !   flag = .FALSE.
  !   write(logUnit(1),*) "test_computeWallForces FAILED"
  !   write(logUnit(1),*) "calculated force = ", Fcoll
  !   write(logUnit(1),*) "Error = ", abs(Fcoll - Fcoll_exact)
  ! end if

  ! ---- TEST #2 damping forces ---- !
  kn = 0.0
  dn = 1.0

  Fcoll = computeWallForce_DPS( this         = particle,     &
                           & boundaryData = boundaryData, &
                           & eps          = eps,          &
                           & kn           = kn,           &
                           & dn           = dn            )

  ! Overlap should be 0.1 in x-direction. So Fcoll should be -0.1
  ! See if it matches exact answer
  Fcoll_exact = (/ 0.1_rk, -0.1_rk, 0.0_rk /)

  if( all( abs( Fcoll - Fcoll_exact) < tol )  ) then
    flag = .TRUE.
    write(logUnit(1),*) "test_computeWallForces SUCCEEDED"
    write(logUnit(1),*) "calculated force = ", Fcoll
  else
    flag = .FALSE.
    write(logUnit(1),*) "test_computeWallForces FAILED"
    write(logUnit(1),*) "calculated force = ", Fcoll
    write(logUnit(1),*) "Error = ", abs(Fcoll - Fcoll_exact)
  end if




end subroutine test_computeWallForces

subroutine test_DEM_fillNeighborList(flag)
  logical, intent(out) :: flag
  ! ------------------------------------------!
  type(mus_particle_group_type) :: particleGroup
  type(mus_particle_DPS_type) :: particle
  real(kind=rk) :: xp0(3), xp1(3), xp2(3), xp3(3)
  real(kind=rk) :: d0, Rp
  integer :: iParticle, iNgh
  logical :: wasAdded

  call init_da_particle_DPS( particleGroup%particles_DPS, 4 )

  ! Position of particle whose neighbor list we will fill
  xp0 = (/ 0.5_rk, 0.5_rk, 0.5_rk /)
  xp1 = (/ 0.3_rk, 0.5_rk, 0.5_rk /)
  xp2 = (/ 0.5_rk, 0.3_rk, 0.7_rk /)
  xp3 = (/ 0.2_rk, 0.2_rk, 0.2_rk /)

  ! Particle radius (same for all particles)
  Rp = 0.05_rk

  ! Include particles in neighborlist if they're less than this distance away
  d0 = 0.25_rk

  ! Append the mock particles to particleGroup
  particle%pos(1:3) = xp0
  particle%pos(4:6) = 0.0_rk
  particle%particleID = 1
  particle%radius   = Rp

  call append_da_particle_DPS( me     = particleGroup%particles_DPS, &
                           & particle = particle,           &
                           & length   = 1,                  &
                           & wasAdded = wasAdded            )

  ! Append the mock particles to particleGroup
  particle%pos(1:3) = xp1
  particle%pos(4:6) = 0.0_rk
  particle%particleID = 2
  particle%radius   = Rp

  call append_da_particle_DPS( me     = particleGroup%particles_DPS, &
                           & particle = particle,           &
                           & length   = 1,                  &
                           & wasAdded = wasAdded            )

  ! Append the mock particles to particleGroup
  particle%pos(1:3) = xp2
  particle%pos(4:6) = 0.0_rk
  particle%particleID = 3
  particle%radius   = Rp

  call append_da_particle_DPS( me     = particleGroup%particles_DPS, &
                           & particle = particle,           &
                           & length   = 1,                  &
                           & wasAdded = wasAdded            )

  ! Append the mock particles to particleGroup
  particle%pos(1:3) = xp3
  particle%pos(4:6) = 0.0_rk
  particle%particleID = 4
  particle%radius   = Rp

  call append_da_particle_DPS( me     = particleGroup%particles_DPS, &
                           & particle = particle,           &
                           & length   = 1,                  &
                           & wasAdded = wasAdded            )

  ! Check if adding the particles went OK
  do iParticle = 1, 4
    write(logUnit(1),*) "Particle ", iParticle
    write(logUnit(1),*) "pos ", particleGroup%particles_DPS%val(iParticle)%pos(1:3)
    write(logUnit(1),*) "radius ", particleGroup%particles_DPS%val(iParticle)%radius
  end do

  call DEM_fillNeighborList( particleGroup = particleGroup, &
                           & d0            = d0             )

  ! Print out the neighbor lists
  do iParticle = 1, 4
    write(logUnit(1),*) "Particle ", iParticle
    write(logUnit(1),*) "Neighbor list: "
    do iNgh = 1, particleGroup%particles_DPS%val(iParticle)%DEM_neighborList%nvals
      write(logUnit(1),*) particleGroup%particles_DPS%val(iParticle) &
        &                                            %DEM_neighborList%val(iNgh)
    end do
  end do


end subroutine test_DEM_fillNeighborList

subroutine test_generateElemListLine(scheme, geometry, params)
  !> Scheme for access to fluid data
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to tree
  type(mus_geom_type), intent(in) :: geometry
  !> Params
  type(mus_param_type), intent(in) :: params
  ! ------------------------------ !
  type(dyn_intArray_type) :: elemList
  integer :: dir(3)
  real(kind=rk) :: xstart(3), length
  integer :: iElem, lev, coord(4)
  integer(kind=long_k) :: TIDoffset
  integer :: ldPos
  real(kind=rk) :: x(3), dx
  ! ------------------------------ !
  lev = geometry%tree%global%maxlevel
  dx  = params%physics%dxLvl(lev)
  TIDoffset = tem_firstIdAtLevel(lev)

  dir = (/ 0, 0, 1 /)
  xstart = (/ 0.5_rk, 0.5_rk, 0.1_rk /)
  length = 0.7_rk

  ! Create the elemList
  call generateElemListLine( dir = dir, &
                           & xstart = xstart,     &
                           & length = length,     &
                           & scheme = scheme,     &
                           & geometry = geometry, &
                           & elemList = elemList  )

  ! To test, print out the coordinates of elems in the list
  do iElem = 1, elemList%nvals
    ldPos =  elemList%val(iElem)
    coord = tem_coordOfID( TreeID = scheme%levelDesc(lev)%total(ldPos), &
                         & offset = TIDoffset)

    x = getBaryOfCoord( coord  = coord,                       &
                      & origin = geometry%tree%global%origin, &
                      & dx     = dx                           )
    write(logUnit(1), *) "x = ", x

  end do

end subroutine test_generateElemListLine

subroutine test_intp1D_peskin(flag, plogunit)
  ! logical to indicate whether test passed
  ! set to TRUE if test is successful
  logical :: flag
  integer :: plogunit
  ! -------------------- !
  ! query points r_lat
  real(kind=rk) :: a, b, c
  ! weights:
  real(kind=rk) :: wa, wb, wc
  ! tolerance when comparing to known result
  real(kind=rk) :: tol
  ! -------------------- !

  ! Set query values and tolerance
  a = 0.6_rk
  b = 1.43_rk
  c = 2.2_rk
  tol = 1.0e-3_rk

  ! Compute weights and compare to known results
  wa = intp_1D_peskin( r_lat = a )
  wb = intp_1D_peskin( r_lat = b )
  wc = intp_1D_peskin( r_lat = c )

  write(plogunit,*) "test_intp1D_peskin"
  write(plogunit,*) "wa = ", wa
  write(plogunit,*) "wb = ", wb
  write(plogunit,*) "wc = ", wc

  if( abs(wa - 0.4_rk) < tol      .AND. &
    & abs(wb - 0.091592_rk) < tol .AND. &
    & abs(wc - 0.0_rk) < tol            ) then
    flag = .TRUE.
  else
    flag = .FALSE.
  end if

end subroutine test_intp1D_peskin

subroutine test_compute_fluid_momentum( scheme, geometry, params, &
                                      & comm, myRank, nProcs, plogunit, flag )
  !> Scheme for access to auxField
  type(mus_scheme_type), intent(inout) :: scheme
  !> Geometry for access to global number of elems
  type(mus_geom_type), intent(in) :: geometry
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  !> MPI communicator
  integer, intent(in) :: comm
  !> My MPI rank
  integer, intent(in) :: myRank
  !> Number of MPI processes in communicator
  integer, intent(in) :: nProcs
  !> Unit to log test result to
  integer, intent(in) :: plogunit
  !> flag gets set to TRUE if test succeeded, FALSE if it failed
  logical, intent(out) :: flag
  ! ---------------------------------------- !
  real(kind=rk) :: totalMomentum(3)
  real(kind=rk) :: correctMomentum(3)
  integer(kind=long_k) :: Nelems_global
  integer :: iElem
  integer :: elemOff
  integer :: lev
  integer :: dens_pos, vel_pos(3)
  real(kind=rk) :: tol, dx, convertFac_mom
  ! ---------------------------------------- !
  Nelems_global = geometry%tree%global%nElems
  tol = 1.0e-5_rk
  lev = geometry%tree%global%maxLevel
  dx = params%physics%dxLvl(lev)
  dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
  vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)
  convertFac_mom = params%physics%rho0 * dx**3 * params%physics%fac(lev)%vel

  ! First initialize the auxField rho and u to all 1.0
  ! Loop over all local fluid cells
  do iElem = 1, scheme%pdf( lev )%nElems_fluid
    elemOff = (iElem-1)*scheme%varSys%nAuxScalars
    ! Calculate fluid momentum of this cell and add to total
    scheme%auxField(lev)%val(elemOff + dens_pos) = 1.0_rk
    scheme%auxField(lev)%val(elemOff + vel_pos(1)) = 1.0_rk
    scheme%auxField(lev)%val(elemOff + vel_pos(2)) = 0.0_rk
    scheme%auxField(lev)%val(elemOff + vel_pos(3)) = 0.0_rk
  end do

  correctMomentum(1) = real(Nelems_global) * convertFac_mom
  correctMomentum(2) = 0.0_rk
  correctMomentum(3) = 0.0_rk

  ! Call compute_fluid_momentum to compute the total fluid momentum across all procs
  call compute_fluid_momentum( scheme = scheme, &
                              & lev = lev, &
                              & params = params, &
                              & totalMomentum = totalMomentum, &
                              & comm = comm, &
                              & myRank = myRank, &
                              & nProcs = nProcs )

  ! The result should be equal to the global number of elements in the mesh.
  if(myRank == 0) then
    if( any(abs(totalMomentum - correctMomentum) > tol) ) then
      flag = .FALSE.
      write(plogunit, *) "test_compute_fluid_momentum FAILED"
      write(plogunit, *) "totalMomentum = ", totalMomentum
      write(plogunit, *) "correctMomentum = ", correctMomentum
    else
      flag = .TRUE.
      write(plogunit, *) "test_compute_fluid_momentum SUCCEEDED"
      write(plogunit, *) "totalMomentum = ", totalMomentum
      write(plogunit, *) "correctMomentum = ", correctMomentum
    end if
  end if

end subroutine test_compute_fluid_momentum

subroutine test_compute_particle_momentum(lev, params, comm, myRank, nProcs, &
                                         &  plogunit, flag                    )
  integer, intent(in) :: lev
  !> Params for access to dt, dx, etc.
  type(mus_param_type), intent(in) :: params
  !> MPI communicator
  integer, intent(in) :: comm
  !> My MPI rank
  integer, intent(in) :: myRank
  !> Number of MPI processes in communicator
  integer, intent(in) :: nProcs
  !> Unit to log test result to
  integer, intent(in) :: plogunit
  !> flag gets set to TRUE if test succeeded, FALSE if it failed
  logical, intent(out) :: flag
  ! ------------------------------------------ !
  type(mus_particle_group_type) :: particleGroup
  type(mus_particle_DPS_type) :: particle
  integer :: Nparticles, iParticle
  logical :: wasAdded
  real(kind=rk) :: totalMomentum(3), correctMomentum(3), tol
  ! ------------------------------------------ !
  ! Tolerance to compare computed momentum to the correct result
  tol = 1.0e-5

  ! Initialize the particleGroup on each process with 5 mock particle objects
  Nparticles = 5

  ! The total momentum of the mock particles on each process should be equal to
  ! the total number of mock particles over all processes.
  correctMomentum(1) = nProcs*Nparticles
  correctMomentum(2) = 0.0_rk
  correctMomentum(3) = 0.0_rk

  call init_da_particle_DPS(particleGroup%particles_DPS, 1)

  do iParticle = 1, Nparticles
    ! Initialize properties of some mock particles
    particle%vel(1:3) = (/ 1.0_rk, 0.0_rk, 0.0_rk /)
    particle%mass = 1.0_rk
    particle%owner = myRank
    particle%particleID = iParticle

    ! Add them to particleGroup
    call append_da_particle_dps( me        = particleGroup%particles_DPS,     &
                               & particle  = particle,                        &
                               & length    = 1,                               &
                               & wasadded  = wasadded                         )
  end do
  write(logUnit(1),*) "Mock particle nvals = ", particleGroup%particles_DPS%nvals

  call compute_particle_momentum( particleGroup = particleGroup, &
                                & lev = lev, &
                                & params = params, &
                                & totalMomentum = totalMomentum, &
                                & comm = comm, &
                                & myRank = myRank, &
                                & nProcs = nProcs )

  ! Check if total momentum is correct
  ! The result should be equal to the global number of elements in the mesh.
  if(myRank == 0) then
    if( any(abs(totalMomentum - correctMomentum) > tol) ) then
      flag = .FALSE.
      write(plogunit, *) "test_compute_particle_momentum FAILED"
      write(plogunit, *) "totalMomentum = ", totalMomentum
      write(plogunit, *) "correctMomentum = ", correctMomentum
    else
      flag = .TRUE.
      write(plogunit, *) "test_compute_particle_momentum SUCCEEDED"
      write(plogunit, *) "totalMomentum = ", totalMomentum
      write(plogunit, *) "correctMomentum = ", correctMomentum
    end if
  end if

end subroutine test_compute_particle_momentum

subroutine test_global_errorcheck(myRank, plogunit)
  integer, intent(in) :: myRank
  integer, intent(in) :: plogunit
  ! --------------------------------- !
  logical :: flag1, flag2
  ! --------------------------------- !
  ! First test: set local error flag to FALSE for all processes
  ! in this case the flag should remain false after global errorcheck
  flag1 = .FALSE.
  write(plogunit,*) "Rank", myRank, "flag 1 before global errorcheck: ", flag1
  call mus_particles_global_errorcheck(flag1, MPI_COMM_WORLD)
  write(plogunit,*) "Rank", myRank, "flag 1 after global errorcheck: ", flag1

  ! Second test: set local error flag to TRUE on one process (the root)
  ! in this case the flag should be true on all procs after global errorcheck
  if(myRank == 0) then
    flag2 = .TRUE.
  else
    flag2 = .FALSE.
  end if

  write(plogunit,*) "Rank", myRank, "flag 2 before global errorcheck: ", flag2
  call mus_particles_global_errorcheck(flag2, MPI_COMM_WORLD)
  write(plogunit,*) "Rank", myRank, "flag 2 after global errorcheck: ", flag2

  ! Test if both flags are what they should be after the global errorcheck
  if( (.NOT. flag1) .AND. flag2 ) then
    write(plogunit,*) "test_global_errorcheck on rank", myRank, " success!"
  end if

end subroutine test_global_errorcheck

end module mus_particle_tests_module
