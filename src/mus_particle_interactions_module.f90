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
!> In this module the routines for interactions between particles in the DEM
!! solver are implemented.

module mus_particle_interactions_module

  use env_module,                        only: rk, long_k
  use tem_param_module,                  only: PI, div1_6, div1_12
  use tem_logging_module,                only: logUnit
  use tem_aux_module,                    only: tem_abort
  use tem_logging_module,                only: logUnit
  use tem_property_module,               only: prp_solid, prp_particle, &
    &                                          prp_hasBnd, prp_sendHalo
  use tem_geometry_module,               only: tem_CoordOfReal, tem_PosOfId
  use tem_topology_module,               only: tem_IdOfCoord, tem_FirstIdAtLevel
  use tem_stencil_module,                only: tem_stencilHeader_type
  use mus_geom_module,                   only: mus_geom_type
  use mus_scheme_type_module,            only: mus_scheme_type
  use mus_param_module,                  only: mus_param_type
  use mus_particle_type_module,          only: mus_particle_MEM_type, &
    &                                          mus_particle_DPS_type, &
    &                                          mus_particle_group_type
  use mus_particle_comm_type_module,     only: mus_particles_communication_type
  use mus_particle_logging_type_module,  only: mus_particle_logging_type, &
    &                                          pgDebugLog
  use mus_particle_aux_module,           only: cross_product
  use mus_particle_boundary_module,      only: computeDisplacement, pgBndData

  implicit none

  ! --------- INTERFACE BLOCKS FOR GENERIC PROCEDURES -------- !
  ! This is used so that names for similar procedures for different kinds of
  ! particles can be referenced using the same name

  interface isLocalCollision
    module procedure isLocalCollision_MEM
    module procedure isLocalCollision_DPS
  end interface

  interface DEM_isRemoteCollision
    module procedure DEM_isRemoteCollision_MEM
    module procedure DEM_isRemoteCollision_DPS
  end interface

  interface checkAndCollideDEM
    module procedure checkAndCollideDEM_MEM
    module procedure checkAndCollideDEM_DPS
  end interface

  interface computeWallPosSum
    module procedure computeWallPosSum_MEM
    module procedure computeWallPosSum_DPS
  end interface


  interface DEM_collideWithWall
    module procedure DEM_collideWithWall_MEM
    module procedure DEM_collideWithWall_DPS
  end interface


contains


  !> isLocalCollision checks if collision between two particles is local
  !! This is the case when they are both owned by this process
  pure logical function isLocalCollision_MEM(particleA, particleB, myRank)
    !> First particle to collide
    type(mus_particle_MEM_type), intent(in) :: particleA
    !> Second particle to collide
    type(mus_particle_MEM_type), intent(in) :: particleB
    !> This process's rank
    integer, intent(in) :: myRank
    ! ------------------------------------------------------ !

    if(particleA%owner == myRank .AND. particleB%owner == myRank ) then
      isLocalCollision_MEM = .TRUE.
    else
      isLocalCollision_MEM = .FALSE.
    end if
  end function isLocalCollision_MEM

  !> isLocalCollision checks if collision between two particles is local
  !! This is the case when they are both owned by this process
  pure logical function isLocalCollision_DPS(particleA, particleB, myRank)
    !> First particle to collide
    type(mus_particle_DPS_type), intent(in) :: particleA
    !> Second particle to collide
    type(mus_particle_DPS_type), intent(in) :: particleB
    !> This process's rank
    integer, intent(in) :: myRank
    ! ------------------------------------------------------ !

    if(particleA%owner == myRank .AND. particleB%owner == myRank ) then
      isLocalCollision_DPS = .TRUE.
    else
      isLocalCollision_DPS = .FALSE.
    end if
  end function isLocalCollision_DPS

  subroutine DEM_isRemoteCollision_MEM(particleA, particleB, myRank, send, &
                                  & isRemoteCollision, otherRank, otherRankIndex)
    !> First particle to collide
    type(mus_particle_MEM_type), intent(in) :: particleA
    !> Second particle to collide
    type(mus_particle_MEM_type), intent(in) :: particleB
    !> This process's rank
    integer, intent(in) :: myRank
    !> Communication type for sending
    type(mus_particles_communication_type), intent(in) :: send
    !> logical indicating if this is a remote collision that needs to be
    !! resolved on this process
    logical, intent(out) :: isRemoteCollision
    !> Other particle's rank as an output
    integer, intent(out) :: otherRank
    !> Position of other particle's rank in send%proc
    integer, intent(out) :: otherRankIndex
    ! ------------------------------------------------------ !
    integer :: iproc
    integer :: i
    ! ------------------------------------------------------ !
    otherRank = -1
    otherRankIndex = -1

    if(particleA%owner == myRank .OR. particleB%owner == myRank ) then
      ! I should execute a remote collision only if I own at least one of
      ! the participating particles

      if( particleA%owner == particleB%owner) then
        ! If I own both particles it is not a remote collision
        isRemoteCollision = .FALSE.
        return
      end if ! both particles have the same owner

      ! If I own at least one of the particles, find out who owns the other
      if( particleA%owner == myRank ) then
        otherRank = particleB%owner
      else
        otherRank = particleA%owner
      end if

      ! Find the position of the other rank in send%proc
      do iproc = 1, send%nProcs
        if( otherRank == send%proc(iproc) ) then
          otherRankIndex = iproc
          ! Found the position of the other particle's owner in send%proc
          exit
        end if
      end do

      ! Check to make sure we found the proc
      if( otherRank /= send%proc(iproc) ) then
        write(logUnit(1), '(A)') '------------------------------------'
        write(logUnit(1),*) 'ERROR isRemoteCollision myRank = ', myRank, ' : could not find particle rank'
        write(logUnit(1),*) 'particle A with ID = ', particleA%particleID, ' owner ', particleA%owner
        write(logUnit(1),*) 'particle B with ID = ', particleB%particleID, ' owner ', particleB%owner
        write(logUnit(1),*) 'otherRank = ', otherRank
        write(logUnit(1), '(A)', advance = 'no') 'send%proc = ['
        write(logUnit(1), '(4I6)', advance = 'no') ( send%proc(i), i = 1, send%nProcs)
        write(logUnit(1), '(A)') ']'
        write(logUnit(1), '(A)') '------------------------------------'

        call tem_abort()
      end if

      ! Check if the collision can be detected on the other process
      if( particleA%existsOnProc(iproc) .AND. particleB%existsOnProc(iproc) ) then
        ! Collision can also be detected on other proc. In that case the
        ! proc with highest rank wll execute the collision
        if( otherRank < myRank) then
          isRemoteCollision = .TRUE.
        else
          isRemoteCollision = .FALSE.
        end if
      else
        ! Collision cannot be detected on other proc
        ! In this case I definitely have to handle the collision
        isRemoteCollision = .TRUE.
      end if
    else
      ! If neither of the two particles are owned by me, another rank
      ! will handle the collision
      isRemoteCollision = .FALSE.
    end if ! I own at least one particle
  end subroutine DEM_isRemoteCollision_MEM


  subroutine DEM_isRemoteCollision_DPS(particleA, particleB, myRank, send, &
                                  & isRemoteCollision, otherRank, otherRankIndex)
    !> First particle to collide
    type(mus_particle_DPS_type), intent(in) :: particleA
    !> Second particle to collide
    type(mus_particle_DPS_type), intent(in) :: particleB
    !> This process's rank
    integer, intent(in) :: myRank
    !> Communication type for sending
    type(mus_particles_communication_type), intent(in) :: send
    !> logical indicating if this is a remote collision that needs to be
    !! resolved on this process
    logical, intent(out) :: isRemoteCollision
    !> Other particle's rank as an output
    integer, intent(out) :: otherRank
    !> Position of other particle's rank in send%proc
    integer, intent(out) :: otherRankIndex
    ! ------------------------------------------------------ !
    integer :: iproc
    integer :: i
    ! ------------------------------------------------------ !
    otherRank = -1
    otherRankIndex = -1

    if(particleA%owner == myRank .OR. particleB%owner == myRank ) then
      ! I should execute a remote collision only if I own at least one of
      ! the participating particles

      if( particleA%owner == particleB%owner) then
        ! If I own both particles it is not a remote collision
        isRemoteCollision = .FALSE.
        return
      end if ! both particles have the same owner

      ! If I own at least one of the particles, find out who owns the other
      if( particleA%owner == myRank ) then
        otherRank = particleB%owner
      else
        otherRank = particleA%owner
      end if

      ! Find the position of the other rank in send%proc
      do iproc = 1, send%nProcs
        if( otherRank == send%proc(iproc) ) then
          otherRankIndex = iproc
          ! Found the position of the other particle's owner in send%proc
          exit
        end if
      end do

      ! Check to make sure we found the proc
      if( otherRank /= send%proc(iproc) ) then
        write(logUnit(1), '(A)') '------------------------------------'
        write(logUnit(1),*) 'ERROR isRemoteCollision myRank = ', myRank, ' : could not find particle rank'
        write(logUnit(1),*) 'particle A with ID = ', particleA%particleID, ' owner ', particleA%owner
        write(logUnit(1),*) 'particle B with ID = ', particleB%particleID, ' owner ', particleB%owner
        write(logUnit(1),*) 'otherRank = ', otherRank
        write(logUnit(1), '(A)', advance = 'no') 'send%proc = ['
        write(logUnit(1), '(4I6)', advance = 'no') ( send%proc(i), i = 1, send%nProcs)
        write(logUnit(1), '(A)') ']'
        write(logUnit(1), '(A)') '------------------------------------'

        call tem_abort()
      end if

      ! Check if the collision can be detected on the other process
      if( particleA%existsOnProc(iproc) .AND. particleB%existsOnProc(iproc) ) then
        ! Collision can also be detected on other proc. In that case the
        ! proc with highest rank wll execute the collision
        if( otherRank < myRank) then
          isRemoteCollision = .TRUE.
        else
          isRemoteCollision = .FALSE.
        end if
      else
        ! Collision cannot be detected on other proc
        ! In this case I definitely have to handle the collision
        isRemoteCollision = .TRUE.
      end if
    else
      ! If neither of the two particles are owned by me, another rank
      ! will handle the collision
      isRemoteCollision = .FALSE.
    end if ! I own at least one particle
  end subroutine DEM_isRemoteCollision_DPS

  ! ---------------- COLLISION DETECTION ROUTINES ----------------- !
  !> checkAndCollideDEM checks if two particles A and B collide
  !! This is the case if the continuous representations overlap
  !! If there is a collision, a collision force is applied to each particle
  !> checkAndCollideDEM checks if two particles A and B collide
  !! This is the case if the continuous representations overlap
  !! If there is a collision, a collision force is applied to each particle
  subroutine checkAndCollideDEM_MEM( particleA, particleB, hasCollided, eps, Tc, mu )
    !> First particle to collide
    type(mus_particle_MEM_type), intent(inout) :: particleA
    !> Second particle to collide
    type(mus_particle_MEM_type), intent(inout) :: particleB
    !> Logical which indicates whether particles have indeed collided
    logical, intent(out) :: hasCollided
    !> Threshold gap at which to call it a collision
    real(kind=rk), intent(in) :: eps
    !> Collision time, used to compute DEM spring and damper constant
    real(kind=rk), intent(in) :: Tc
    !> Dynamic viscosity, in physical units
    real(kind=rk), intent(in) :: mu
    ! ------------------------------------------------------ !
    ! vector joining the origins of two colliding particles
    real(kind=rk) :: rab(3), length_rab, norm_rab(3)
    ! Relative velocity of particle B w.r.t. A
    real(kind=rk) :: uab(3)
    ! normal velocities of a and b, un = relative normal velocity
    ! ut = relative tangential velocity
    real(kind=rk) :: un, ut(3)
    ! Surface velocities due to particle rotational motion
    ! and relative velocity usurf_ab of B w.r.t. A
    real(kind=rk) :: usurf_a(3), usurf_b(3), usurf_ab(3)
    ! gap between particles
    real(kind=rk) :: dab
    ! Average radius of the two particles
    real(kind=rk) :: Ravg
    ! Dry coefficient of restitution
    real(kind=rk) :: e_dry
    ! mass and effective mass of the colliding particles
    real(kind=rk) :: m_a, m_b, m_eff
    ! Collision force, spring and damper constant
    real(kind=rk) :: kn, dn
    ! Normalized gap h = dab/R and critical gap hc
    ! (below this we don't consider lubrication)
    ! as the lubrication force diverges for h -> 0.
    ! hmax is the minimum gap at which we compute lubrication
    ! forces
    real(kind=rk) :: h, hmax, hc
    ! Lubrication force in normal direction
    real(kind=rk) :: Flub_n(3)
    ! Lubrication force in tangential direction
    real(kind=rk) :: Flub_t(3)
    ! Collision force
    real(kind=rk) :: Fcoll(3)
    ! ------------------------------------------------------ !
    hc = 0.01_rk
    hmax = 1.0_rk
    hasCollided = .FALSE.
    Ravg = 0.5*(particleA%radius + particleB%radius)
    rab =  computeDisplacement( x1           = particleA%pos(1:3), &
                              & x2           = particleB%pos(1:3), &
                              & boundaryData = pgBndData           )
    e_dry = 1.0_rk

    length_rab = dot_product(rab, rab)
    length_rab = sqrt( length_rab )

    ! compute separation distance
    dab = length_rab - (particleA%radius + particleB%radius + eps)
    h = 2.0*dab/(particleA%radius + particleB%radius)

    if( h < hmax ) then
      ! Normal vector to plane of collision with length 1
      norm_rab = rab / length_rab

      ! get normal and tangential velocity magnitudes
      uab = particleB%vel(1:3) - particleA%vel(1:3)
      un = dot_product( uab, norm_rab )
      ut = uab - un*rab

      usurf_a = calc_surface_vel( particleA%vel(4:6), particleA%radius, norm_rab )
      usurf_b = calc_surface_vel( particleB%vel(4:6), particleB%radius, -norm_rab )
      usurf_ab = usurf_b - usurf_a

      Flub_n = computeLubForce_normal( h        = h,       &
                                     & hc       = hc,      &
                                     & mu       = mu,      &
                                     & r        = Ravg,    &
                                     & un       = un,      &
                                     & norm_rab = norm_rab )

      Flub_t = computeLubForce_tangential( h   = h,       &
                                         & hc  = hc,      &
                                         & mu  = mu,      &
                                         & r   = Ravg,    &
                                         & ut  = ut,      &
                                         & utr = usurf_ab )

      ! Apply forces to the particle
      particleA%F_DEM( particleA%F_DEM_next,1:3 ) &
        & = particleA%F_DEM( particleA%F_DEM_next,1:3 ) - Flub_n - Flub_t

      particleB%F_DEM( particleB%F_DEM_next,1:3 ) &
        & = particleB%F_DEM( particleB%F_DEM_next,1:3 ) + Flub_n + Flub_t

      if( dab < 0.0_rk ) then
        ! Calculate collision forces
        m_a = particleA%mass
        m_b = particleB%mass
        m_eff = m_a*m_b/(m_a+m_b)

        ! Compute spring and damper constants
        ! both kn and dn are always > 0
        kn = m_eff*( PI**2 + (log(e_dry) )**2 )/Tc**2
        dn = (-2.0*m_eff*log(e_dry))/Tc

        Fcoll = kn*dab*norm_rab

        ! Apply forces to the particle
        particleA%F_DEM( particleA%F_DEM_next,1:3 ) &
          & = particleA%F_DEM( particleA%F_DEM_next,1:3 ) + Fcoll

        particleB%F_DEM( particleB%F_DEM_next,1:3 ) &
          & = particleB%F_DEM( particleB%F_DEM_next,1:3 ) - Fcoll

      end if ! dab < 0.0

      hasCollided = .TRUE.

    end if ! h < hc

  end subroutine checkAndCollideDEM_MEM


  !> checkAndCollideDEM checks if two particles A and B collide
  !! This is the case if the continuous representations overlap
  !! If there is a collision, a collision force is applied to each particle
  subroutine checkAndCollideDEM_DPS( particleA, particleB, hasCollided, eps, Tc, mu )
    !> First particle to collide
    type(mus_particle_DPS_type), intent(inout) :: particleA
    !> Second particle to collide
    type(mus_particle_DPS_type), intent(inout) :: particleB
    !> Logical which indicates whether particles have indeed collided
    logical, intent(out) :: hasCollided
    !> Threshold gap at which to call it a collision
    real(kind=rk), intent(in) :: eps
    !> Collision time, used to compute DEM spring and damper constant
    real(kind=rk), intent(in) :: Tc
    !> Dynamic viscosity, in physical units
    real(kind=rk), intent(in) :: mu
    ! ------------------------------------------------------ !
    ! vector joining the origins of two colliding particles
    real(kind=rk) :: rab(3), length_rab, norm_rab(3)
    ! Relative velocity of particle B w.r.t. A
    real(kind=rk) :: uab(3)
    ! normal velocities of a and b, un = relative normal velocity
    ! ut = relative tangential velocity
    real(kind=rk) :: un, ut(3)
    ! Surface velocities due to particle rotational motion
    ! and relative velocity usurf_ab of B w.r.t. A
    real(kind=rk) :: usurf_a(3), usurf_b(3), usurf_ab(3)
    ! gap between particles
    real(kind=rk) :: dab
    ! Average radius of the two particles
    real(kind=rk) :: Ravg
    ! Dry coefficient of restitution
    real(kind=rk) :: e_dry
    ! mass and effective mass of the colliding particles
    real(kind=rk) :: m_a, m_b, m_eff
    ! Collision force, spring and damper constant
    real(kind=rk) :: kn, dn
    ! Normalized gap h = dab/R and critical gap hc
    ! (below this we don't consider lubrication)
    ! as the lubrication force diverges for h -> 0.
    ! hmax is the minimum gap at which we compute lubrication
    ! forces
    real(kind=rk) :: h, hmax, hc
    ! Lubrication force in normal direction
    real(kind=rk) :: Flub_n(3)
    ! Lubrication force in tangential direction
    real(kind=rk) :: Flub_t(3)
    ! Collision force
    real(kind=rk) :: Fcoll(3)
    ! ------------------------------------------------------ !
    hc = 0.01_rk
    hmax = 1.0_rk
    hasCollided = .FALSE.
    Ravg = 0.5*(particleA%radius + particleB%radius)
    rab =  computeDisplacement( x1           = particleA%pos(1:3), &
                              & x2           = particleB%pos(1:3), &
                              & boundaryData = pgBndData           )
    e_dry = 1.0_rk

    length_rab = dot_product(rab, rab)
    length_rab = sqrt( length_rab )

    ! compute separation distance
    dab = length_rab - (particleA%radius + particleB%radius + eps)
    h = 2.0*dab/(particleA%radius + particleB%radius)

    if( h < hmax ) then
      ! Normal vector to plane of collision with length 1
      norm_rab = rab / length_rab

      ! get normal and tangential velocity magnitudes
      uab = particleB%vel(1:3) - particleA%vel(1:3)
      un = dot_product( uab, norm_rab )
      ut = uab - un*rab

      usurf_a = calc_surface_vel( particleA%vel(4:6), particleA%radius, norm_rab )
      usurf_b = calc_surface_vel( particleB%vel(4:6), particleB%radius, -norm_rab )
      usurf_ab = usurf_b - usurf_a

      Flub_n = computeLubForce_normal( h        = h,       &
                                     & hc       = hc,      &
                                     & mu       = mu,      &
                                     & r        = Ravg,    &
                                     & un       = un,      &
                                     & norm_rab = norm_rab )

      Flub_t = computeLubForce_tangential( h   = h,       &
                                         & hc  = hc,      &
                                         & mu  = mu,      &
                                         & r   = Ravg,    &
                                         & ut  = ut,      &
                                         & utr = usurf_ab )

      ! Apply forces to the particle
      particleA%F_DEM( particleA%F_DEM_next,1:3 ) &
        & = particleA%F_DEM( particleA%F_DEM_next,1:3 ) - Flub_n - Flub_t

      particleB%F_DEM( particleB%F_DEM_next,1:3 ) &
        & = particleB%F_DEM( particleB%F_DEM_next,1:3 ) + Flub_n + Flub_t

      if( dab < 0.0_rk ) then
        ! Calculate collision forces
        m_a = particleA%mass
        m_b = particleB%mass
        m_eff = m_a*m_b/(m_a+m_b)

        ! Compute spring and damper constants
        ! both kn and dn are always > 0
        kn = m_eff*( PI**2 + (log(e_dry) )**2 )/Tc**2
        dn = (-2.0*m_eff*log(e_dry))/Tc

        Fcoll = kn*dab*norm_rab

        ! Apply forces to the particle
        particleA%F_DEM( particleA%F_DEM_next,1:3 ) &
          & = particleA%F_DEM( particleA%F_DEM_next,1:3 ) + Fcoll

        particleB%F_DEM( particleB%F_DEM_next,1:3 ) &
          & = particleB%F_DEM( particleB%F_DEM_next,1:3 ) - Fcoll

      end if ! dab < 0.0

      hasCollided = .TRUE.

    end if ! h < hc

  end subroutine checkAndCollideDEM_DPS

  function calc_surface_vel( om, particle_radius, n ) result(usurf)
    real(kind=rk), intent(in) :: om(3)
    real(kind=rk), intent(in) :: particle_radius
    real(kind=rk), intent(in) :: n(3)
    real(kind=rk) :: usurf(3)
    ! ----------------------------------------- !
    real(kind=rk) :: r(3)
    ! ----------------------------------------- !
    r = particle_radius*n
    call cross_product( om, r, usurf)
  end function calc_surface_vel

  function computeLubForce_normal(h, hc, mu, R, un, norm_rab) result(Flub)
    !> Gap between particles (physical units)
    real(kind=rk), intent(in) :: h
    !> Critical gap (physical units)
    !! If h < hc we use hc to compute the lubrication forces
    !! to prevent them from diverging as h -> 0.
    real(kind=rk), intent(in) :: hc
    !> Dynamic viscosity (physical units)
    real(kind=rk), intent(in) :: mu
    !> Radius of the particles
    real(kind=rk), intent(in) :: R
    !> Relative velocity in normal direction between the two particles
    real(kind=rk), intent(in) :: un
    !> Unit vector pointing in the direction from particle A to particle B
    real(kind=rk), intent(in) :: norm_rab(3)
    !> Output: lubrication force in normal direction
    real(kind=rk) :: Flub(3)
    ! ------------------------------------------------ !
    if( h > hc ) then
      Flub = -6*PI*R*mu*un*(0.25/h &
          &   - (9.0/40.0)*log10(h) - (3.0/112.0)*h*log10(h))*norm_rab
    else
      Flub = -6*PI*R*mu*un*(0.25/hc &
          &   - (9.0/40.0)*log10(hc) - (3.0/112.0)*hc*log10(hc))*norm_rab
    end if
  end function computeLubForce_normal

  function computeLubForce_tangential(h, hc, mu, R, ut, utr) result(Flub)
    !> Gap between particles (physical units)
    real(kind=rk), intent(in) :: h
    !> Critical gap (physical units)
    !! If h < hc we use hc to compute the lubrication forces
    !! to prevent them from diverging as h -> 0.
    real(kind=rk), intent(in) :: hc
    !> Dynamic viscosity (physical units)
    real(kind=rk), intent(in) :: mu
    !> Radius of the particles
    real(kind=rk), intent(in) :: R
    !> Relative velocity in tangential direction between the two particles
    real(kind=rk), intent(in) :: ut(3)
    !> Relative velocity in tangential direction caused by rotational motion of particles
    real(kind=rk), intent(in) :: utr(3)
    !> Output: lubrication force in tangential direction
    real(kind=rk) :: Flub(3)
    ! ------------------------------------------------ !
    if( h > hc ) then
      Flub = -6*PI*R*mu*( ut*(-div1_6*log10(h)) &
        &                + utr*( -div1_6*log10(h) - div1_12*h*log10(h) ) )
    else
      Flub = -6*PI*R*mu*( ut*(-div1_6*log10(hc)) &
        &                + utr*( -div1_6*log10(hc) - div1_12*hc*log10(hc) ) )
    end if
  end function computeLubForce_tangential

  !> computeWallForce_1D computes the 1D DEM force between a particle at position
  !! xp with velocity un to a wall at position xwall.
  function computeWallForce_1D(xp, up, Rp, xwall, kn, dn, eps) result(Fwall)
    !> Particle position
    real(kind=rk), intent(in) :: xp
    !> Particle velocity
    real(kind=rk), intent(in) :: up
    !> Particle radius
    real(kind=rk), intent(in) :: Rp
    !> Wall position
    real(kind=rk), intent(in) :: xwall
    !> Spring constant
    real(kind=rk), intent(in) :: kn
    !> Damping coefficient
    real(kind=rk), intent(in) :: dn
    !> Threshold distance at which to execute a collision
    real(kind=rk), intent(in) :: eps
    ! Output: 1D wall force
    real(kind=rk) :: Fwall
    ! ------------------------------------------!
    real(kind=rk) :: rwall, nwall, del, d, un

    rwall = xwall - xp
    nwall = rwall/abs(rwall)

    del = rwall - (Rp + eps)*nwall
    d   = del*nwall
    if( d < 0.0_rk ) then
      un = up*nwall
      Fwall = (kn*d - dn*un)*nwall
    else
      Fwall = 0.0_rk
    end if
  end function computeWallForce_1D



  !> computeWallPosSum computes the sum of the locations of wall boundaries in the
  !! area surrounding the particle. It also sets rmflag to true if an open boundary
  !! is detected, indicating that this particle should be removed from the global
  !! domain (in a different routine)
  subroutine computeWallPosSum_MEM( this, BCinteraction, scheme, stencil, geometry,    &
                              & params, rmflag, foundWall, wallPosSum, nWallPos )
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Array of length nBCs which tells us how particles should interact with
    !! each boundary. Indexed using boundary ID!
    integer, intent(in) :: BCinteraction(:)
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Stencil for access to cq
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> Geometry information to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> Parameters for access to conversion factors between units
    type(mus_param_type), intent(in) :: params
    !> Flag to indicate whether particle should be removed because it just
    !! hit a non-wall BC
    logical, intent(out) :: rmflag
    !> Flag to indicate whether we found a wall or not
    logical, intent(out) :: foundWall
    !> Sum of locations (x,y,z) of the detected walls
    real(kind=rk), intent(out) :: wallPosSum(3)
    !> Number of elements in wallPosSum
    integer, intent(out) :: nWallPos
    ! ------------------------------------------- !
    real(kind=rk) :: x, y, z
    integer :: currentCoord(4)
    integer(kind=long_k) :: TreeID
    integer :: lev
    integer(kind=long_k) :: BCid
    integer :: ldPos, neighPos
    integer :: posInBnd
    integer(kind=long_k) :: treePos
    integer :: nx, ny, nz, iDir
    integer(kind=long_k) :: elemProp


    ! Number of lattice cells to check in each direction starting
    ! from the particle origin
    integer :: searchRadius

    ! vector of current lattice direction
    real(kind=rk) :: cq_phy(3)
    ! position and average position of all the wall points we find
    real(kind=rk) :: wallPos(3)

    real(kind=rk) :: dx
    ! ------------------------------------------ !
    rmflag = .FALSE.

    lev = this%coordOfOrigin(4)
    dx = params%physics%dxLvl(lev)
    nWallPos = 0
    WallPos = 0.0_rk
    wallPosSum = 0.0_rk


    searchRadius = this%Rn + 3
    ! Identify other nodes in neighborhood of origin belonging to particle
    do nx = -searchRadius, searchRadius
      do ny = -searchRadius, searchRadius
        do nz = -searchRadius, searchRadius

          if( nx**2 + ny**2 + nz**2 <= (this%Rn - 1)**2 ) then
            ! These elements are definitely part of the particle itself.
            ! So we don't need to check them for walls
            cycle
          end if
          elemProp = 0.0_rk
          ! integer coordinate of this element
          currentCoord = this%coordOfOrigin(:) + (/ nx,ny,nz,0 /)

          ! get TreeID and position in complete tree
          TreeID = tem_IdOfCoord(coord = currentCoord)
          treePos = tem_PosOfId( sTreeID    = TreeID,             &
                               & treeIDlist = geometry%tree%treeID)

          if( treePos == 0 ) then
            ! If element is not found in treeID total list
            ! it is probably outside the domain so skip it.
            cycle
          end if

          ! Get position of this elemenet in level descriptor list
          ldPos = geometry%levelPointer(treePos)

          ! Compute cartesian coordinates relative to particle origin
          x = scheme%levelDesc(lev)%baryOfTotal(ldPos, 1) - this%pos(1)
          y = scheme%levelDesc(lev)%baryOfTotal(ldPos, 2) - this%pos(2)
          z = scheme%levelDesc(lev)%baryOfTotal(ldPos, 3) - this%pos(3)

          if( x**2 + y**2 + z**2 >= this%radius**2 ) then
            ! if node is NOT within particle radius (so does not belong
            ! to particle) then check it for a boundary property

            elemProp = scheme%levelDesc(lev)%property(ldPos)
            if ( btest( elemProp, prp_hasBnd ) ) then

              ! --- Check if this is a wall or other boundary ---
              ! Get position of element in boundary ID array
              posInBnd = geometry%posInBndID(ldPos)

              ! For all the elements with prp_hasBnd, loop over the neighboring
              ! elements to find in which direction the boundary is
              linkLoop: do iDir = 1, stencil%QQN
                bcID = geometry%boundary%boundary_ID(iDir,posInBnd)

                ! bcID = 0 if no boundary in this direction
                ! bcID > 0 if boundary (wall or open)
                ! bcID < 0 if periodic boundary
                if(bcID > 0_long_k) then
                ! Check how particles should interact with this boundary
                if( BCinteraction(int(bcID)) == 1 ) then
                  ! If we hit an open boundary, exit the routine and remove this
                  ! particle from the group in a subsequent call to destroyParticle
                  ! and remove_particle_from_da_particle_MEM
                  rmflag = .TRUE.
                  return
                else
                  neighPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, ldPos)

                  ! NeighPos = 0 indicates a regular wall boundary condition
                  ! NeighPos < 0 indicates a different boundary condition
                  if(neighPos <= 0) then
                    ! There is a wall in this direction.
                    nWallPos = nWallPos + 1


                    ! Stencil direction vector in physical coordinates
                    cq_phy = stencil%cxDirRK(:,iDir) * dx

                    ! Estimate its location as half a lattice spacing away
                    wallPos = scheme%levelDesc(lev)%baryOfTotal(ldPos,1:3) &
                      & + 0.5 * cq_phy

                    wallPosSum = wallPosSum + wallPos

                  end if ! neighPos <= 0
                end if ! BCinteraction(bcID) == 1
                end if ! bcID > 0
              end do linkLoop
            end if ! hasBnd
          end if ! x^2 + y^2 + z^2 > radius^2
        end do ! nz
      end do ! ny
    end do ! nx

    ! If we found a wall, set foundWall to true
    if( nWallPos > 0 ) then
      foundWall = .TRUE.
    else
      ! We didn't find any walls so set foundWall to false
      foundWall = .FALSE.
      return
    end if

  end subroutine computeWallPosSum_MEM



  !> computeWallPosSum computes the sum of the locations of wall boundaries in the
  !! area surrounding the particle. It also sets rmflag to true if an open boundary
  !! is detected, indicating that this particle should be removed from the global
  !! domain (in a different routine)
  subroutine computeWallPosSum_DPS( this, BCinteraction, scheme, stencil, geometry,    &
                              & params, rmflag, foundWall, wallPosSum, nWallPos )
    type(mus_particle_DPS_type), intent(inout) :: this
    !> Array of length nBCs which tells us how particles should interact with
    !! each boundary. Indexed using boundary ID!
    integer, intent(in) :: BCinteraction(:)
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Stencil for access to cq
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> Geometry information to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> Parameters for access to conversion factors between units
    type(mus_param_type), intent(in) :: params
    !> Flag to indicate whether particle should be removed because it just
    !! hit a non-wall BC
    logical, intent(out) :: rmflag
    !> Flag to indicate whether we found a wall or not
    logical, intent(out) :: foundWall
    !> Sum of locations (x,y,z) of the detected walls
    real(kind=rk), intent(out) :: wallPosSum(3)
    !> Number of elements in wallPosSum
    integer, intent(out) :: nWallPos
    ! ------------------------------------------- !
    real(kind=rk) :: x, y, z
    integer :: currentCoord(4)
    integer(kind=long_k) :: TreeID
    integer :: lev
    integer(kind=long_k) :: BCid
    integer :: ldPos, neighPos
    integer :: posInBnd
    integer(kind=long_k) :: treePos
    integer :: nx, ny, nz, iDir
    integer(kind=long_k) :: elemProp


    ! Number of lattice cells to check in each direction starting
    ! from the particle origin
    integer :: searchRadius

    ! vector of current lattice direction
    real(kind=rk) :: cq_phy(3)
    ! position and average position of all the wall points we find
    real(kind=rk) :: wallPos(3)

    real(kind=rk) :: dx
    ! ------------------------------------------ !
    rmflag = .FALSE.

    lev = this%coordOfOrigin(4)
    dx = params%physics%dxLvl(lev)
    nWallPos = 0
    WallPos = 0.0_rk
    wallPosSum = 0.0_rk

    this%nWallPos = 0

    searchRadius = this%Rn + 3
    ! Identify other nodes in neighborhood of origin belonging to particle
    do nx = -searchRadius, searchRadius
      do ny = -searchRadius, searchRadius
        do nz = -searchRadius, searchRadius

          elemProp = 0.0_rk
          ! integer coordinate of this element
          currentCoord = this%coordOfOrigin(:) + (/ nx,ny,nz,0 /)

          ! get TreeID and position in complete tree
          TreeID = tem_IdOfCoord(coord = currentCoord)
          treePos = tem_PosOfId( sTreeID    = TreeID,             &
                               & treeIDlist = geometry%tree%treeID)

          if( treePos == 0 ) then
            ! If element is not found in treeID total list
            ! it is probably outside the domain so skip it.
            cycle
          end if

          ! Get position of this elemenet in level descriptor list
          ldPos = geometry%levelPointer(treePos)

          ! Compute cartesian coordinates relative to particle origin
          x = scheme%levelDesc(lev)%baryOfTotal(ldPos, 1) - this%pos(1)
          y = scheme%levelDesc(lev)%baryOfTotal(ldPos, 2) - this%pos(2)
          z = scheme%levelDesc(lev)%baryOfTotal(ldPos, 3) - this%pos(3)

          if( x**2 + y**2 + z**2 >= this%radius**2 ) then
            ! if node is NOT within particle radius (so does not belong
            ! to particle) then check it for a boundary property

            elemProp = scheme%levelDesc(lev)%property(ldPos)
            if ( btest( elemProp, prp_hasBnd ) ) then

              ! --- Check if this is a wall or other boundary ---
              ! Get position of element in boundary ID array
              posInBnd = geometry%posInBndID(ldPos)

              ! For all the elements with prp_hasBnd, loop over the neighboring
              ! elements to find in which direction the boundary is
              linkLoop: do iDir = 1, stencil%QQN
                bcID = geometry%boundary%boundary_ID(iDir,posInBnd)

                ! bcID = 0 if no boundary in this direction
                ! bcID > 0 if boundary (wall or open)
                ! bcID < 0 if periodic boundary
                if(bcID > 0_long_k) then
                ! Check how particles should interact with this boundary
                if( BCinteraction(int(bcID)) == 1 ) then
                  ! If we hit an open boundary, exit the routine and remove this
                  ! particle from the group in a subsequent call to destroyParticle
                  ! and remove_particle_from_da_particle_MEM
                  rmflag = .TRUE.
                  return
                else
                  neighPos = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, ldPos)

                  ! NeighPos = 0 indicates a regular wall boundary condition
                  ! NeighPos < 0 indicates a different boundary condition
                  if(neighPos <= 0) then
                    ! There is a wall in this direction.
                    nWallPos = nWallPos + 1

                    ! Stencil direction vector in physical coordinates
                    cq_phy = stencil%cxDirRK(:,iDir) * dx

                    ! Estimate its location as half a lattice spacing away
                    wallPos = scheme%levelDesc(lev)%baryOfTotal(ldPos,1:3) &
                      & + 0.5 * cq_phy

                    wallPosSum = wallPosSum + wallPos

                  end if ! neighPos <= 0
                end if ! BCinteraction(bcID) == 1
                end if ! bcID > 0
              end do linkLoop
            end if ! hasBnd
          end if ! x^2 + y^2 + z^2 > radius^2
        end do ! nz
      end do ! ny
    end do ! nx

    ! If we found a wall, set foundWall to true
    if( nWallPos > 0 ) then
      foundWall = .TRUE.
    else
      ! We didn't find any walls so set foundWall to false
      foundWall = .FALSE.
      return
    end if

  end subroutine computeWallPosSum_DPS

  subroutine DEM_collideWithWall_MEM( this, BCinteraction, scheme,          &
                                & geometry, params, Tc, eps )
    type(mus_particle_MEM_type), intent(inout) :: this
    !> Array of length nBCs which tells us how particles should interact with
    !! each boundary. Indexed using boundary ID!
    integer, intent(in) :: BCinteraction(:)
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry information to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> Parameters for access to conversion factors between units
    type(mus_param_type), intent(in) :: params
    !> Collision time
    real(kind=rk), intent(in) :: Tc
    !> threshold gap for collision
    real(kind=rk), intent(in) :: eps

    ! ------------------------------------------- !
    real(kind=rk) :: norm_rwall(3) ! vector from center of particle to wall

    real(kind=rk) :: un                ! particle velocity normal to wall
    real(kind=rk) :: length_rwall, dx  ! distance to wall, mesh size
    real(kind=rk) :: del(3)            ! displacement of DEM spring w.r.t. eq position
    real(kind=rk) :: dab            ! displacement of DEM spring w.r.t. eq position
    real(kind=rk) :: Fwall(3);         ! Wall force
    real(kind=rk) :: kn                ! spring constant
    real(kind=rk) :: dn                ! damping constant
    real(kind=rk) :: e_dry             ! dry restitution coefficient for collision

    integer :: lev
    logical :: foundWall

    ! ------------------------------------------ !
    foundWall = .FALSE.

    lev = this%coordOfOrigin(4)
    dx = params%physics%dxLvl(lev)
    e_dry = 1.0_rk

    ! Update dwall and norm_rwall using the particles displacement since last DEM time step
    this%rwall = this%rwall - ( this%pos(1:3) - this%oldPos(1:3) )

    length_rwall = dot_product(this%rwall,this%rwall)
    length_rwall = sqrt(length_rwall)
    norm_rwall = this%rwall/length_rwall

    ! If we did find a wall, calculate DEM force
    kn = this%mass*( PI**2 + (log(e_dry) )**2 )/Tc**2
    dn = (-2.0*this%mass*log(e_dry))/Tc


    del = this%rwall - ( this%radius + eps)*norm_rwall
    dab = dot_product( del, norm_rwall )

    if( dab < 0.0_rk ) then
      un = dot_product(this%vel(1:3), norm_rwall)

      Fwall = (kn*dab - dn*un)*norm_rwall
      this%F_DEM(this%F_DEM_next,1:3) = this%F_DEM(this%F_DEM_next,1:3) + Fwall

    end if

  end subroutine DEM_collideWithWall_MEM

  subroutine DEM_collideWithWall_DPS( this, BCinteraction, scheme,          &
                                & geometry, params, Tc, eps )
    type(mus_particle_DPS_type), intent(inout) :: this
    !> Array of length nBCs which tells us how particles should interact with
    !! each boundary. Indexed using boundary ID!
    integer, intent(in) :: BCinteraction(:)
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry information to determine TreeIDs of elements 'covered' by particle
    type(mus_geom_type), intent(in) :: geometry
    !> Parameters for access to conversion factors between units
    type(mus_param_type), intent(in) :: params
    !> Collision time
    real(kind=rk), intent(in) :: Tc
    !> threshold gap for collision
    real(kind=rk), intent(in) :: eps

    ! ------------------------------------------- !
    real(kind=rk) :: norm_rwall(3) ! vector from center of particle to wall

    real(kind=rk) :: un                ! particle velocity normal to wall
    real(kind=rk) :: length_rwall, dx  ! distance to wall, mesh size
    real(kind=rk) :: del(3)            ! displacement of DEM spring w.r.t. eq position
    real(kind=rk) :: dab            ! displacement of DEM spring w.r.t. eq position
    real(kind=rk) :: Fwall(3);         ! Wall force
    real(kind=rk) :: kn                ! spring constant
    real(kind=rk) :: dn                ! damping constant
    real(kind=rk) :: e_dry             ! dry restitution coefficient for collision

    integer :: lev
    logical :: foundWall

    ! ------------------------------------------ !
    foundWall = .FALSE.

    lev = this%coordOfOrigin(4)
    dx = params%physics%dxLvl(lev)
    e_dry = 1.0_rk

    ! Update dwall and norm_rwall using the particles displacement since last DEM time step
    this%rwall = this%rwall - ( this%pos(1:3) - this%oldPos(1:3) )

    length_rwall = dot_product(this%rwall,this%rwall)
    length_rwall = sqrt(length_rwall)
    norm_rwall = this%rwall/length_rwall

    ! If we did find a wall, calculate DEM force
    kn = this%mass*( PI**2 + (log(e_dry) )**2 )/Tc**2
    dn = (-2.0*this%mass*log(e_dry))/Tc


    del = this%rwall - ( this%radius + eps)*norm_rwall
    dab = dot_product( del, norm_rwall )

    if( dab < 0.0_rk ) then
      un = dot_product(this%vel(1:3), norm_rwall)

      Fwall = (kn*dab - dn*un)*norm_rwall
      this%F_DEM(this%F_DEM_next,1:3) = this%F_DEM(this%F_DEM_next,1:3) + Fwall
    end if

  end subroutine DEM_collideWithWall_DPS

  subroutine DEM_collideWithPlanarWall_DPS( this, dir, xwall, Tc, eps )
    !> Particle to collide with wall
    type(mus_particle_DPS_type), intent(inout) :: this
    !> Cartesian direction to collide with (x,y,z) = (1,2,3)
    integer, intent(in) :: dir
    !> Wall position (in one of three Cartesian coordinates
    real(kind=rk), intent(in) :: xwall
    !> Collision time
    real(kind=rk), intent(in) :: Tc
    !> threshold gap for collision
    real(kind=rk), intent(in) :: eps
    ! ------------------------------------------- !
    real(kind=rk) :: xp, rwall, norm_rwall
    real(kind=rk) :: kn, dn, Fwall
    real(kind=rk) :: del, un
    real(kind=rk) :: e_dry   ! dry coefficient of restitution
    ! ------------------------------------------- !
    e_dry = 1.0_rk

    xp = this%pos(dir)
    rwall = xwall - xp
    norm_rwall = rwall / abs(rwall)

    del = rwall - (this%radius + eps)*norm_rwall

    if( del*norm_rwall < 0.0_rk ) then
      kn = this%mass*( PI**2 + (log(e_dry) )**2 )/Tc**2
      dn = (-2.0*this%mass*log(e_dry))/Tc

      un = this%vel(dir)
      Fwall = kn*del - dn*un

      this%F_DEM(this%F_DEM_next, dir) =  this%F_DEM(this%F_DEM_next, dir) + Fwall
    end if
  end subroutine DEM_collideWithPlanarWall_DPS



  subroutine DEM_computeWallForces_MEM( particleGroup, scheme, geometry, &
                                  & params, Tc, eps              )
    !> particleGroup to search for collisions in
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
    integer :: iParticle
    integer :: myRank
    logical :: rmflag       ! Flag to indicate whether particle should be removed
    ! ------------------------------------------!
    rmflag = .FALSE.
    myRank = params%general%proc%rank

    do iParticle = 1, particleGroup%particles_MEM%nvals

      ! Only particle owner computes wall interaction force using this%rwall
      ! this%rwall itself is set at the beginning of the subcycling loop using
      ! information from ALL processes.
      if( particleGroup%particles_MEM%val(iParticle)%interactWithWall &
        & .AND. particleGroup%particles_MEM%val(iParticle)%owner == myRank) then
        call DEM_collideWithWall( this = particleGroup%particles_MEM%val(iParticle), &
                                & BCinteraction = particleGroup%BC_interaction,  &
                                & scheme        = scheme,                        &
                                & geometry      = geometry,                      &
                                & params        = params,                        &
                                & Tc            = Tc,                            &
                                & eps           = eps                            )
        if( rmflag ) then
          particleGroup%particles_MEM%val(iParticle)%removeParticle_global = .TRUE.
        else
          particleGroup%particles_MEM%val(iParticle)%removeParticle_global = .FALSE.
        end if
      end if
    end do
  end subroutine DEM_computeWallForces_MEM

end module mus_particle_interactions_module
