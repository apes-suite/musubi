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
!> mus_particle_boundary_module contains routines and data types to handle the 
!! interaction of particles with (periodic) boundaries

module mus_particle_boundary_module
  use env_module, only: rk, long_k, newunit
  
  implicit none
  private

  public :: mus_particle_boundarydata_type
  public :: pgBndData
  public :: getNeighborCoord
  public :: wrapPeriodicCoord
  public :: wrapPeriodicPos
  public :: computeDisplacement
  public :: calcPeriodicRsurface

  type mus_particle_boundarydata_type
    !> boundary locations [xmin, xmax, ymin, ymax, zmin, zma]
    real(kind=rk) :: bnd(6)
    !> length of domain in x, y, z directions
    real(kind=rk) :: domain_size(3)
    !> integer coordinates of bnd in x, y, z directions
    integer :: bnd_coord(6)
    !> logical set to TRUE if particle domain boundaries are active
    logical :: useBnd
    ! logical set to TRUE if periodic boundary 
    ! periodicBnd(1) = periodicBnd(2) = TRUE if periodic in x
    ! periodicBnd(3) = periodicBnd(4) = TRUE if periodic in y 
    ! periodicBnd(5) = periodicBnd(6) = TRUE if periodic in z
    logical :: periodicBnd(6)
    ! logical array set to true if boundary should be treated as a wall in 
    ! the corresponding direction
    ! wallBnd(i) = TRUE if bnd(i) should be treated as a wall
    logical :: wallBnd(6)
  end type mus_particle_boundarydata_type
  
  type(mus_particle_boundarydata_type), save :: pgBndData


contains


  !> getNeighborCoord gets the coordinate of the element
  !! offset from input coord by (nx,ny,nz) while taking into
  !! account periodicity
  function getNeighborCoord( coord, nx, ny, nz, boundaryData ) &
    &        result(neighborCoord)
    integer, intent(in) :: coord(4)
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    type(mus_particle_boundarydata_type), intent(in) :: boundaryData
    integer :: neighborCoord(4)
    ! ----------------------------- ! 
    integer :: xsize, ysize, zsize
    integer :: xmin, xmax, ymin, ymax, zmin, zmax
    ! ----------------------------- ! 
    neighborCoord = coord + (/ nx, ny, nz, 0 /)
  
    if (boundaryData%useBnd) then
      xmin = boundaryData%bnd_coord(1)
      xmax = boundaryData%bnd_coord(2)
      ymin = boundaryData%bnd_coord(3)
      ymax = boundaryData%bnd_coord(4)
      zmin = boundaryData%bnd_coord(5)
      zmax = boundaryData%bnd_coord(6)
  
      xsize = xmax-xmin + 1
      ysize = ymax-ymin + 1
      zsize = zmax-zmin + 1
  
      ! If we have periodic boundaries, wrap the coordinate
      ! x - direction
      if( boundaryData%periodicBnd(1) ) then
        if( neighborCoord(1) < xmin) then
          neighborCoord(1) = neighborCoord(1) + xsize
        else if( neighborCoord(1) > xmax) then
          neighborCoord(1) = neighborCoord(1) - xsize
        end if
      end if
  
      ! y - direction
      if( boundaryData%periodicBnd(3) ) then
        if( neighborCoord(2) < ymin) then
          neighborCoord(2) = neighborCoord(2) + ysize
        else if( neighborCoord(2) > ymax) then
          neighborCoord(2) = neighborCoord(2) - ysize
        end if
      end if
  
      ! z - direction
      if( boundaryData%periodicBnd(5) ) then
        if( neighborCoord(3) < zmin) then
          neighborCoord(3) = neighborCoord(3) + zsize
        else if( neighborCoord(3) > zmax) then
          neighborCoord(3) = neighborCoord(3) - zsize
        end if
      end if
    end if ! boundaryData%useBnd
  end function getNeighborCoord
  
  !> wrapPeriodicCoord modifies the input coord to take into
  !! account periodicity
  subroutine wrapPeriodicCoord(coord, boundaryData)
    integer, intent(inout) :: coord(4)
    type(mus_particle_boundarydata_type), intent(in) :: boundaryData
    ! ----------------------------- ! 
    integer :: xsize, ysize, zsize
    integer :: xmin, xmax, ymin, ymax, zmin, zmax
    ! ----------------------------- ! 
    if (boundaryData%useBnd) then
      xmin = boundaryData%bnd_coord(1)
      xmax = boundaryData%bnd_coord(2)
      ymin = boundaryData%bnd_coord(3)
      ymax = boundaryData%bnd_coord(4)
      zmin = boundaryData%bnd_coord(5)
      zmax = boundaryData%bnd_coord(6)
  
      xsize = xmax-xmin + 1
      ysize = ymax-ymin + 1
      zsize = zmax-zmin + 1
  
      ! If we have periodic boundaries, wrap the coordinate
      ! x - direction
      if (boundaryData%periodicBnd(1)) then
        if (coord(1) < xmin) then
          coord(1) = coord(1) + xsize
        else if( coord(1) > xmax) then
          coord(1) = coord(1) - xsize
        end if
      end if
  
      ! y - direction
      if (boundaryData%periodicBnd(3)) then
        if (coord(2) < ymin) then
          coord(2) = coord(2) + ysize
        else if( coord(2) > ymax) then
          coord(2) = coord(2) - ysize
        end if
      end if
  
      ! z - direction
      if (boundaryData%periodicBnd(5)) then
        if (coord(3) < zmin) then
          coord(3) = coord(3) + zsize
        else if( coord(3) > zmax) then
          coord(3) = coord(3) - zsize
        end if
      end if
    end if ! boundaryData%useBnd
  end subroutine wrapPeriodicCoord
  
  ! -------------- FUNCTIONS FOR PERIODIC BOUNDARIES -------------- !
  
  subroutine wrapPeriodicPos(pos, boundaryData)
    real(kind=rk), intent(inout) :: pos(3)
    type(mus_particle_boundarydata_type) :: boundaryData
    ! ------------------------------- !
    real(kind=rk) :: xmin 
    real(kind=rk) :: xmax 
    real(kind=rk) :: ymin 
    real(kind=rk) :: ymax 
    real(kind=rk) :: zmin 
    real(kind=rk) :: zmax 
    ! ------------------------------- !
    ! If particle boundary data is active, use it to incorporate the effect of 
    ! periodic boundaries
    if (boundaryData%useBnd) then
      xmin = boundaryData%bnd(1)
      xmax = boundaryData%bnd(2)
      ymin = boundaryData%bnd(3)
      ymax = boundaryData%bnd(4)
      zmin = boundaryData%bnd(5)
      zmax = boundaryData%bnd(6)
  
      if (pgBndData%periodicBnd(1)) then
        ! x - direction
        if (pos(1) < xmin) then
          pos(1) = xmax - ( xmin - pos(1) )
        else if( pos(1) > xmax) then
          pos(1) = xmin + ( pos(1) - xmax )
        end if
      end if
  
      ! y - direction
      if (pgBndData%periodicBnd(3)) then
        if (pos(2) < ymin) then
          pos(2) = ymax - ( ymin - pos(2) )
        else if( pos(2) > ymax) then
          pos(2) = ymin + ( pos(2) - ymax )
        end if
      end if
  
      ! z - direction
      if (pgBndData%periodicBnd(5)) then
        if (pos(3) < zmin) then
          pos(3) = zmax - ( zmin - pos(3) )
        else if( pos(3) > zmax) then
          pos(3) = zmin + ( pos(3) - zmax )
        end if
      end if
    end if ! boundaryData%useBnd
  end subroutine wrapPeriodicPos
  
  !> computeDistance computes the shortest distance between points x1 and x2
  !! In doing so it takes possible periodic boundaries into account.
  function computeDisplacement(x1, x2, boundaryData) result( r12 )
    !> xyz coordinates of point 1
    real(kind=rk) :: x1(3)
    !> xyz coordinates of point 2
    real(kind=rk) :: x2(3)
    !> Boundary data tells us the domain bounds and if we have periodic bounds
    type(mus_particle_boundarydata_type), intent(in) :: boundaryData
    !> Output: shortest vector pointing from x1 to x2
    real(kind=rk) :: r12(3)
    ! ------------------------- !
    !> magnitude of r12
    real(kind=rk) :: r12_mag(3)  
    integer :: i, iBnd
    ! ------------------------- !
  
    r12 = x2 - x1
    r12_mag = abs(r12)
  
    ! If particle boundary data is active, use it to incorporate the effect of 
    ! periodic boundaries
    if (boundaryData%useBnd) then
      ! If domain is periodic in x, y, z direction adjust distance computation
      ! in each dimension accordingly
      do i = 1,3
        iBnd = 1 + 2*(i-1)
        if (boundaryData%periodicBnd(iBnd)) then
          if (r12_mag(i) > 0.5*boundaryData%domain_size(i)) then
            r12(i) = -1 * ( r12(i) / r12_mag(i) ) &
              &         * ( boundaryData%domain_size(i) - r12_mag(i) )
          end if
        end if
      end do
    end if
    
  end function computeDisplacement
  
  !>  wrap_periodic checks whether a distance r from x_particle to 
  !   bary_of_surface r is less than the radius of the particle R. 
  !   If not, this indicates that the barycenter should wrap around 
  !   the periodic boundary i.e. we should add or subtract the domain 
  !   length L from it.
  subroutine calcPeriodicDistanceToSurface( ri, R, L )
    !> Distance (in one of the Cartesian directions xi) from the particle 
    !! origin to the barycenter of a point on the surface.
    real(kind=rk), intent(inout) :: ri
    !> Radius of the particle
    real(kind=rk), intent(in) :: R
    !> Length of the periodic domain (in direction xi )
    real(kind=rk), intent(in) :: L
    ! --------------------------------------------------------- !
    if( abs(ri) > R ) then
      if( ri < 0.0_rk ) then
        ri = ri + L
      else
        ri = ri - L
      end if
    end if
  
  end subroutine calcPeriodicDistanceToSurface
  
  !> calcPeriodicRsurface is used to calculate the vector from the 
  !! particle origin to a point on the surface in presence of periodic 
  !! boundaries.
  !! Usage: first calculate the distance from the particle origin to a 
  !! surface element using r = baryOfSurface - x_origin. If the particle is 
  !! close to a periodic boundary this vector may not be correct.
  !! In that case a call to calcPeriodicRsurface modifies the vector r to 
  !! take into account the periodicity. This is done by checking if the 
  !! magnitude of r is less than the particle radius R_particle (which can 
  !! be modified with some tolerance if needed).
  !! This routine is used for fully resolved (MEM) particles only.
  subroutine calcPeriodicRsurface( r, R_particle, boundaryData )
    !> Original vector r
    real(kind=rk), intent(inout) :: r(3)
    !> Particle radius + a tolerance if desired
    real(kind=rk), intent(inout) :: R_particle
    !> Datatype containing periodic boundary data
    type(mus_particle_boundarydata_type), intent(in) :: boundaryData
    ! --------------------------------------------------------- !
    integer :: k, i
    ! --------------------------------------------------------- !
  
    i = 1
    do k = 1, 5, 2
      if (boundaryData%periodicBnd(k)) then
        call calcPeriodicDistanceToSurface( ri = r(i),                       &
          &                                 R  = R_particle,                 &
          &                                 L  = boundaryData%domain_size(i) )
      end if
      i = i + 1
    end do
  
  end subroutine calcPeriodicRsurface

end module mus_particle_boundary_module
