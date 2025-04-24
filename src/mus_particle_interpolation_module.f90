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
!> In this module the routines and data types for interpolation of fluid properties
!! at particle locations are implemented.

module mus_particle_interpolation_module
  use env_module,              only: rk, labelLen
  use tem_logging_module,      only: logUnit
  use tem_param_module,        only: div1_3, div1_6, PI, cs2inv, cs4inv
  use tem_stencil_module,      only: tem_stencilHeader_type, &
    &                                d3q27_cxDir, &
    &                                d2q9_cxDir
  use mus_particle_aux_module, only: getPathToVec

  implicit none

  abstract interface
    function wghtFunc(r_lat) result(wght)
      import :: rk
      real(kind=rk), intent(in) :: r_lat
      real(kind=rk) :: wght
    end function wghtFunc
  end interface

  !> Data type containing information required by routines to interpolate
  !! fluid properties to particle locations
  type mus_particle_interpolator_type
    !> Boundaries of interpolation stencil
    !! This indicates the offset of cells to each side
    !! of the cell containing the point to be interpolated
    !! intpBnds = (/ x_lo, x_hi, y_lo, y_hi, z_lo, z_hi /)
    !! For example: intpBnds = (/ -1, 1, -1, 1, -1, 1 /) will
    !! loop over all neighbor cells, intpBnds = (/ -1 1, -1, 1, 0, 0 /)
    !! will only loop over neighbors in x and y direction etc.
    integer :: bnd_x(2)
    integer :: bnd_y(2)
    integer :: bnd_z(2)

    !> Directions of neighboring elements to interpolate from
    integer, allocatable :: neighDirs(:,:)

    !> "Paths" in terms of iDirs with which to reach neighboring stencil elements
    !! This depends on the stencil being used. We need this because not every element
    !! to interpolate from is accessible using just one neighDir e.g. when using the
    !! D3Q19 stencil the "corner elements" cannot be reached using just one neighDir.
    integer, allocatable :: neighPaths(:,:)

    !> Number of elements (neighbors + the center element) to interpolate from.
    !! If all direct neighbors are used this is 27 in 3D and 9 in 2D.
    integer :: Nelems

    !> Functions to calculate interpolation weights
    procedure(wghtFunc), pointer, nopass :: getWght_x
    procedure(wghtFunc), pointer, nopass :: getWght_y
    procedure(wghtFunc), pointer, nopass :: getWght_z

    !> Label of interpolation scheme used
    !! can be e.g. 'delta', 'linear', or 'peskin'
    !! Peskin stencil currently only works in serial.
    character(len=labelLen) :: interpolation_kind
  end type mus_particle_interpolator_type


contains


  subroutine init_particle_interpolator( interpolator, stencil )
    type(mus_particle_interpolator_type), intent(inout) :: interpolator
    type(tem_stencilHeader_type), intent(in) :: stencil
    ! ----------------------------------------- !
    integer :: iDir, path(3), vec(3)
    ! ----------------------------------------- !

    select case( trim(stencil%label) )
    case('d3q27')
      interpolator%Nelems = 27
      ! Allocate space for neighPaths array
      allocate(interpolator%neighPaths( 3,interpolator%Nelems ))
      allocate(interpolator%neighDirs( 3,interpolator%Nelems ))
      interpolator%neighDirs = d3q27_cxDir

      do iDir = 1, interpolator%Nelems
        vec = d3q27_cxDir(1:3,iDir)
        interpolator%neighPaths(1,iDir) = iDir
        interpolator%neighPaths(2:3,iDir) = stencil%restPosition
      end do

    case('d3q19')
      interpolator%Nelems = 27
      ! Allocate space for neighPaths array
      allocate(interpolator%neighPaths( 3,interpolator%Nelems ))
      allocate(interpolator%neighDirs( 3,interpolator%Nelems ))
      interpolator%neighDirs = d3q27_cxDir

      do iDir = 1, interpolator%Nelems
        vec = d3q27_cxDir(1:3,iDir)
        call getPathToVec( vec, stencil%cxDir, path )
        interpolator%neighPaths(1:3,iDir) = path
      end do
    case('d2q9')
      interpolator%Nelems = 9
      ! Allocate space for neighPaths array
      allocate(interpolator%neighPaths( 3,interpolator%Nelems ))
      allocate(interpolator%neighDirs( 3,interpolator%Nelems ))
      interpolator%neighDirs = d2q9_cxDir

      do iDir = 1, interpolator%Nelems
        vec = d2q9_cxDir(1:3,iDir)
        call getPathToVec( vec, stencil%cxDir, path )
        interpolator%neighPaths(1:3,iDir) = path
      end do

    case default
      write(logUnit(1),*) "Error init_particle_interpolator: stencil kind not supported"
    end select
  end subroutine init_particle_interpolator

  subroutine finalize_particle_interpolator(interpolator)
    type(mus_particle_interpolator_type), intent(inout) :: interpolator
    ! ----------------------------------------- !
    deallocate(interpolator%neighPaths)
    deallocate(interpolator%neighDirs)
  end subroutine finalize_particle_interpolator

  !> 1D discrete delta interpolation function. Used to interpolate the
  !! fluid property of a lattice cell with barycenter xf to the position
  !! of a particle xp.
  pure function intp_1D_delta(r_lat) result(del)
    !> Distance from particle to barycenter of fluid cell
    real(kind=rk), intent(in) :: r_lat
    !> Output: weighting factor
    real(kind=rk) :: del
    ! ------------------------------------------!
    if( r_lat < 0.5_rk ) then
      del = div1_3*( 1.0 + sqrt(-3.0*r_lat**2 + 1.0_rk) )
    else if( r_lat < 1.5_rk) then
      del = div1_6*( 5.0 - 3.0*r_lat - sqrt( -3.0*(1-r_lat)**2 + 1.0) )
    else
      del = 0.0_rk
    end if

  end function intp_1D_delta

  !> 1D discrete delta interpolation function according to Peskin.
  !! Used to interpolate the
  !! fluid property of a lattice cell with barycenter xf to the position
  !! of a particle xp.
  !! Note that this function has a support of 2*dx so requires both neighbor
  !! and next-neighbor for interpolation
  pure function intp_1D_peskin(r_lat) result(del)
    !> Distance from particle to barycenter of fluid cell (must be positive)
    real(kind=rk), intent(in) :: r_lat
    !> Output: weighting factor
    real(kind=rk) :: del
    ! ------------------------------------------!
    if( r_lat <= 1.0_rk ) then
      del = 0.125*( 3.0_rk - 2*r_lat + sqrt(1.0 + 4*r_lat - 4*r_lat**2 ) );
    else if( r_lat <= 2.0_rk) then
      del = 0.125*( 5.0_rk - 2*r_lat - sqrt(-7.0 + 12*r_lat - 4*r_lat**2 ) );
    else
      del = 0.0_rk
    end if

  end function intp_1D_peskin

  !> Weight function for 1d linear interpolation or distribution.
  !! Given the distance from particle to fluid cell barycenter,
  !! it returns weight, where weight is a linear function
  !! of the distance between the position of the particle and the
  !! barycenter of the lattice site we are distributing to.
  pure function getwght1d_linear(r) result(weight)
    !> Distance r = xp - x/dx from particle to lattice barycenter
    real(kind=rk), intent(in) :: r
    !> Output: part of f distributed to this lattice site
    real(kind=rk) :: weight
    ! ------------------------------------------!
    real(kind=rk) :: abs_r
    ! ------------------------------------------!
    abs_r = abs(r)

    if(abs_r > 1.0) then
      weight = 0.0_rk
    else
      weight = (1.0 - abs_r)
    end if
  end function

  !> Function to return interpolation weight 1.0 regardless of input.
  !! We need this as the weight for the z-direction interpolation in
  !! case of d2q9 stencil.
  pure function one(r_lat) result(del)
    !> Distance from particle to barycenter of fluid cell
    real(kind=rk), intent(in) :: r_lat
    !> Output: weighting factor
    real(kind=rk) :: del
    ! ------------------------------------------!
    del = 1.0_rk
  end function one

  pure function intp1d_linear(xd, f0, f1) result(f)
    ! Nondimensionalized coordinate xd = (x-x0)/(x1-x0)
    real(kind=rk), intent(in) :: xd
    ! Value of interpolant at f0
    real(kind=rk), intent(in) :: f0(0:)
    ! Value of interpolant at f0
    real(kind=rk), intent(in) :: f1(0:)
    ! Output: interpolated value
    real(kind=rk) :: f( lbound(f0,1):ubound(f0,1))
    ! ------------------------------------------!
    f = f0*(1-xd) + f1*xd
  end function

  pure function get_xd(x_lat) result(xd)
    !> Input: distance from particle coordOfOrigin barycenter
    !! to particle position. Choose x, y or z component depending
    !! on which direction we are interpolating
    real(kind=rk), intent(in) :: x_lat
    !> Output: normalized interpolation coordinate
    !! xd = (x-x0)/(x1-x0) = (x-x0)/dx
    real(kind=rk) :: xd
    ! ------------------------------------------!
    if(x_lat >= 0.0_rk) then
      xd = x_lat
    else
      xd = 1.0_rk - abs(x_lat)
    end if
  end function

  !> getInterpolationBnds is used to determine which cells to interpolate
  !! fluid quantities from for a particle positioned at its coordOfOrigin + r_lat.
  !! The 8 cells obtained are the cells whose barycenter forms the tightest
  !! bounding cube around the particle.
  subroutine getInterpolationBnds(r_lat, lo_x, up_x, lo_y, up_y, lo_z, up_z)
    real(kind=rk), intent(in) :: r_lat(3)
    integer, intent(out) :: lo_x
    integer, intent(out) :: up_x
    integer, intent(out) :: lo_y
    integer, intent(out) :: up_y
    integer, intent(out) :: lo_z
    integer, intent(out) :: up_z
    ! ------------------------------------------!
    if( r_lat(1) >= 0 ) then
      lo_x = 0
      up_x = 1
    else
      lo_x = -1
      up_x = 0
    end if

    if( r_lat(2) >= 0 ) then
      lo_y = 0
      up_y = 1
    else
      lo_y = -1
      up_y = 0
    end if

    if( r_lat(3) >= 0 ) then
      lo_z = 0
      up_z = 1
    else
      lo_z = -1
      up_z = 0
    end if
  end subroutine getInterpolationBnds

  subroutine printParticleInterpolator(interpolator, logUnit)
    type(mus_particle_interpolator_type), intent(in) :: interpolator
    integer, intent(in) :: logUnit
    ! ------------------------------------------!

    write(logUnit,*) '---- Particle interpolator settings ----'
    write(logUnit,*) 'bnd_x = ', interpolator%bnd_x
    write(logUnit,*) 'bnd_y = ', interpolator%bnd_y
    write(logUnit,*) 'bnd_z = ', interpolator%bnd_z
    write(logUnit,*) 'interpolation_kind = ', &
                      & trim(interpolator%interpolation_kind)

  end subroutine printParticleInterpolator

end module mus_particle_interpolation_module
