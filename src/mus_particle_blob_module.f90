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
!> mus_particle_blob_module contains types for particle "blobs"
!! representing data that can be used to quickly create groups
!! (blobs) of particles with certain shapes like prisms, cylinders
!! etc.

module mus_particle_blob_module
  use env_module,               only: rk, LabelLen
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit
  use tem_grow_array_module,    only: init, append, destroy, empty, &
    &                                 grw_realarray_type,           &
    &                                 grw_real2darray_type
  use mus_particle_aux_module,  only: cross_product
  use mus_particle_prob_module, only: normcdf

  implicit none

  interface fill_prism
    module procedure fill_prism_by_distance, fill_prism_by_number
  end interface

  type mus_particle_blob_prob_type
    character(len=labelLen) :: kind

    integer :: seed

    real(kind=rk) :: mu(2)
    real(kind=rk) :: sigma(2)
    !> Positions that can be picked at random. Stored in a 2D array as
    !! [ [x1, x2, ... xn], [y1, y2, ...yn] ]
    !! So the i-th position is given by [x,y] = available_positions%val(1:2,i)
    type(grw_real2darray_type) :: available_positions
    type(grw_real2darray_type) :: chosen_positions
    type(grw_realarray_type) :: probabilities
    type(grw_realarray_type) :: bins
  end type mus_particle_blob_prob_type

  type mus_particle_blob_cylinder_type
    real(kind=rk) :: origin(3)
    real(kind=rk) :: vec(3)
    real(kind=rk) :: radius
    ! Properties for particles inside blob
    real(kind=rk) :: particle_radius
    real(kind=rk) :: particle_mass
    real(kind=rk) :: particle_vel(6)
    real(kind=rk) :: particle_force(6)

    !> If this is TRUE, particles will be initialized to the
    !! local fluid velocity instead of particle_vel
    logical :: init_particles_to_fluid_vel

    !> Probability distribution to use to place particles inside
    !! the cylindrical blob
    type(mus_particle_blob_prob_type) :: distribution
  end type mus_particle_blob_cylinder_type

  type mus_particle_blob_prism_type
    real(kind=rk) :: origin(3)
    real(kind=rk) :: vec_x(3)
    real(kind=rk) :: vec_y(3)
    real(kind=rk) :: vec_z(3)

    integer :: nx
    integer :: ny
    integer :: nz

    ! Properties for particles inside blob
    real(kind=rk) :: particle_radius
    real(kind=rk) :: particle_mass
    real(kind=rk) :: particle_vel(6)
    real(kind=rk) :: particle_force(6)

    !> If this is TRUE, particles will be initialized to the
    !! local fluid velocity instead of particle_vel
    logical :: init_particles_to_fluid_vel

    !> Probability distribution to use to place particles inside
    !! the cylindrical blob
    type(mus_particle_blob_prob_type) :: distribution
  end type mus_particle_blob_prism_type


  type(mus_particle_blob_cylinder_type), save :: particleblob
  type(mus_particle_blob_prism_type), save :: particleblob_prism

  interface init_particle_blob_prob
    module procedure init_particle_blob_prob_circle
    module procedure init_particle_blob_prob_square
  end interface


contains


  subroutine set_random_seed(s)
    integer, intent(in) :: s
    ! --------------------------------------------- !
    integer :: seed_length
    integer, allocatable :: seed(:)
    ! --------------------------------------------- !
    ! Initialize the random seed to a hard-coded value.
    ! This ensures that each MPI process picks the same sequence of
    ! random numbers. That is necessary for the procedure to create
    ! unique particle IDs for each particle to work.
    call random_seed(size=seed_length)
    allocate(seed(seed_length))
    seed(1:seed_length) = s ! set seed to some hard-coded integer value
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine set_random_seed

  subroutine init_particle_blob_prob_circle(blob, d, R, mu, sigma, Nchosen)
    type(mus_particle_blob_prob_type), intent(inout) :: blob
    !> Distance between particle positions. This should be chosen larger
    !! than one particle diameter to prevent overlap.
    real(kind=rk), intent(in) :: d
    !> Radius within which all positions must be located
    real(kind=rk) :: R
    !> Mean value of the probability distributions in x and y directions (mu_x, mu_y)
    real(kind=rk) :: mu(2)
    !> Standard deviation of the probability distributions in x and y directions (sig_x, sig_y)
    real(kind=rk) :: sigma(2)
    !> Number of particles we plan to choose. This determines the size
    !! to which the chosen_positions array will be allocated
    integer, intent(in) :: Nchosen
    ! --------------------------------------------- !
    integer :: Nx, Ny, Npositions
    integer :: ix, iy
    real(kind=rk) :: pos(2)
    real(kind=rk) :: xlim(2), ylim(2)
    real(kind=rk) :: x2, y2, R2
    real(kind=rk) :: prob_x, prob_y, prob
    real(kind=rk) :: hi, lo, p_hi, p_lo, ignore
    ! --------------------------------------------- !

    ! Allocate space for the arrays
    Nx = floor( 2*R / d )
    Ny = floor( 2*R / d )
    xlim = [-R, R]
    ylim = [-R, R]
    Npositions = Nx*Ny

    R2 = R**2

    call init( me = blob%available_positions, width = 2, length = Npositions )
    call init( me = blob%chosen_positions, width = 2, length = Nchosen )
    call init( me = blob%probabilities, length = Nchosen )
    call init( me = blob%bins, length = Nchosen+1 )

    ! Fill the array of available positions and the array of
    ! probabilities for these positions
    do ix = 0, Nx - 1
      pos(1) = xlim(1) + (ix+0.5)*d
      x2 = pos(1)**2
      hi = pos(1) + 0.5*d
      lo = pos(1) - 0.5*d
      call normcdf( x     = hi ,      &
                  & mu    = mu(1),    &
                  & sigma = sigma(1), &
                  & cum   = p_hi,     &
                  & ccum  = ignore    )

      call normcdf( x     = lo ,      &
                  & mu    = mu(1),    &
                  & sigma = sigma(1), &
                  & cum   = p_lo,     &
                  & ccum  = ignore    )

      prob_x = p_hi-p_lo
      do iy = 0, Ny - 1
        pos(2) = ylim(1) + (iy+0.5)*d
        y2 = pos(2)**2
        if( x2+y2 >= R2 ) then
          cycle
        end if

        hi = pos(2) + 0.5*d
        lo = pos(2) - 0.5*d
        call normcdf( x     = hi ,      &
                    & mu    = mu(2),    &
                    & sigma = sigma(2), &
                    & cum   = p_hi,     &
                    & ccum  = ignore    )

        call normcdf( x     = lo ,      &
                    & mu    = mu(2),    &
                    & sigma = sigma(2), &
                    & cum   = p_lo,     &
                    & ccum  = ignore    )

        prob_y = p_hi - p_lo

        ! Calculate probability of this position being chosen
        prob = prob_x*prob_y

        call append( me = blob%available_positions, val = pos )
        call append( me = blob%probabilities, val = prob )
      end do
    end do

    ! Right now the probabilities do not sum to 1 exactly if the variance is large
    ! (because part of the Gaussian curve extends beyond xlim and ylim).
    ! So we normalize the probabilities such that their sum = 1
    call normalize_probabilities(blob = blob)

    ! Now fill the bins according to the probabilities we just computed
    call fill_bins( blob = blob )

  end subroutine init_particle_blob_prob_circle

  subroutine init_particle_blob_prob_square(blob, d, xlim, ylim, mu, sigma, Nchosen)
    type(mus_particle_blob_prob_type), intent(inout) :: blob
    !> Distance between particle positions. This should be chosen larger
    !! than one particle diameter to prevent overlap.
    real(kind=rk), intent(in) :: d
    !> Range of x-positions particles can have [xmin, xmax]
    real(kind=rk) :: xlim(2)
    !> Range of y-positions particles can have [ymin, ymax]
    real(kind=rk) :: ylim(2)
    !> Mean value of the probability distributions in x and y directions (mu_x, mu_y)
    real(kind=rk) :: mu(2)
    !> Standard deviation of the probability distributions in x and y directions (sig_x, sig_y)
    real(kind=rk) :: sigma(2)
    !> Number of particles we plan to choose. This determines the size
    !! to which the chosen_positions array will be allocated
    integer, intent(in) :: Nchosen
    ! --------------------------------------------- !
    integer :: Nx, Ny, Npositions
    integer :: ix, iy
    real(kind=rk) :: pos(2)
    real(kind=rk) :: prob_x, prob_y, prob
    real(kind=rk) :: hi, lo, p_hi, p_lo, ignore
    ! --------------------------------------------- !

    ! Allocate space for the arrays
    Nx = floor( (xlim(2) - xlim(1)) / d )
    Ny = floor( (ylim(2) - ylim(1)) / d )
    Npositions = Nx*Ny

    call init( me = blob%available_positions, width = 2, length = Npositions )
    call init( me = blob%chosen_positions, width = 2, length = Nchosen )
    call init( me = blob%probabilities, length = Nchosen )
    call init( me = blob%bins, length = Nchosen+1 )

    ! Fill the array of available positions and the array of
    ! probabilities for these positions
    do ix = 0, Nx - 1
      pos(1) = xlim(1) + (ix+0.5)*d
      hi = pos(1) + 0.5*d
      lo = pos(1) - 0.5*d
      call normcdf( x     = hi ,      &
                  & mu    = mu(1),    &
                  & sigma = sigma(1), &
                  & cum   = p_hi,     &
                  & ccum  = ignore    )

      call normcdf( x     = lo ,      &
                  & mu    = mu(1),    &
                  & sigma = sigma(1), &
                  & cum   = p_lo,     &
                  & ccum  = ignore    )

      prob_x = p_hi-p_lo
      do iy = 0, Ny - 1
        pos(2) = ylim(1) + (iy+0.5)*d
        hi = pos(2) + 0.5*d
        lo = pos(2) - 0.5*d
        call normcdf( x     = hi ,      &
                    & mu    = mu(2),    &
                    & sigma = sigma(2), &
                    & cum   = p_hi,     &
                    & ccum  = ignore    )

        call normcdf( x     = lo ,      &
                    & mu    = mu(2),    &
                    & sigma = sigma(2), &
                    & cum   = p_lo,     &
                    & ccum  = ignore    )

        prob_y = p_hi - p_lo

        ! Calculate probability of this position being chosen
        prob = prob_x*prob_y

        call append( me = blob%available_positions, val = pos )
        call append( me = blob%probabilities, val = prob )
      end do
    end do

    ! Right now the probabilities do not sum to 1 exactly if the variance is large
    ! (because part of the Gaussian curve extends beyond xlim and ylim).
    ! So we normalize the probabilities such that their sum = 1
    call normalize_probabilities(blob = blob)

    ! Now fill the bins according to the probabilities we just computed
    call fill_bins( blob = blob )

  end subroutine init_particle_blob_prob_square

  !> Routine which fills the bins in mus_particle_blob_prob_type
  subroutine fill_bins( blob )
    type(mus_particle_blob_prob_type), intent(inout) :: blob
    ! --------------------------------------- !
    integer :: ipos
    real(kind=rk) :: bin_bound
    ! --------------------------------------- !
    ! Set the first bin boundary to 0
    call append( me=blob%bins, val = 0.0_rk )
    do ipos = 1, blob%probabilities%nvals
      ! Then append the bin boundaries for all the other positions
      bin_bound = blob%bins%val(ipos) + blob%probabilities%val(ipos)
      call append( me=blob%bins, val = bin_bound )
    end do
  end subroutine fill_bins

  subroutine normalize_probabilities(blob)
    type(mus_particle_blob_prob_type), intent(inout) :: blob
    ! --------------------------------- !
    real(kind=rk) :: sum_of_probs
    integer :: nvals
    ! --------------------------------- !
    nvals = blob%probabilities%nvals
    sum_of_probs = sum( blob%probabilities%val(1:nvals) )

    blob%probabilities%val(1:nvals) = blob%probabilities%val(1:nvals)/sum_of_probs

  end subroutine normalize_probabilities

  subroutine pick_random_position(blob)
    type(mus_particle_blob_prob_type), intent(inout) :: blob
    ! -------------------------------- !
    real(kind=rk) :: r
    integer :: nbins, ibin
    ! -------------------------------- !
    nbins = blob%bins%nvals

    ! Check that there are still available positions left to pick
    if( blob%available_positions%nvals < 1 )  then
      write(*,*) "ERROR pick_random_position: no more available positions left!"
      call tem_abort()
    end if

    ! Generate a random number between 0 and 1 (uniform distribution)
    ! until we get one that's inside the bin range.
    do
      call random_number(r)

      ! Determine what bin it falls into
      call map_to_bin( x    = r, &
                    & bins = blob%bins%val(1:nbins), &
                    & N    = nbins, &
                    & ibin = ibin )
      if(ibin > 0) exit
    end do

    ! ibin is now the index of the position we will pick
    ! Add it to the list of chosen positions
    call append( me  = blob%chosen_positions,                 &
               & val = blob%available_positions%val(1:2,ibin) )

    ! Now remove this position from the list of available positions
    call remove_elem_grw_real2darray( me          = blob%available_positions, &
                                    & i           = ibin,                     &
                                    & array_width = 2                         )

    ! Also remove the probability associated with this position
    ! from the probabilities array
    call remove_elem_grw_realarray( me = blob%probabilities, &
                                  & i  = ibin                )

    ! Normalize the probabilities of the remaining elements so their sum equals 1
    call normalize_probabilities(blob)

    ! Recalculate the bins
    blob%bins%nvals = 0
    call fill_bins( blob = blob )

  end subroutine pick_random_position

  subroutine map_to_bin( x, bins, N, ibin)
    !> Number to map to a bin
    real(kind=rk), intent(in) :: x
    !> Array containing the boundaries of each bin
    !! Should be monotonic increasing i.e. bin(0) < bin(1) < ... bin(Nbins)
    real(kind=rk), intent(in) :: bins(:)
    !> Number of elements in bins array
    integer, intent(in) :: N
    !> Output: index of the LOWER bound of the bin containing x
    !! so that bins(ibin) < x <= bins(ibin+1)
    integer, intent(out) :: ibin

    !! EXAMPLE
    !! if bins = [0.0, 0.1, 0.3, 0.7, 1.0]
    !! and x = 0.32
    !! then this routine returns ibin = 3
    !! since bins(ibin) = 0.3 < x <= bins(ibin+1)
    ! ------------------------------- !
    integer :: low, high, mid

    ! First check that x is actually within the bin bounds
    if( bins(1) >= x .OR. x > bins(N) ) then
      write(logUnit(1),*) "WARNING mapToBins: value outside range of bins"
      write(logUnit(1),*) "x = ", x
      write(logUnit(1),*) "bins(1) = ", bins(1)
      write(logUnit(1),*) "bins(N) = ", bins(N)

      ibin = -1
      return
    end if

    low = 1
    high = N
    do
      mid = (low+high)/2
      if( bins(mid) < x .AND. x <= bins(mid+1) ) then
        exit
      else if( x <= bins(mid) ) then
        ! x is to the left of bins(mid)
        high = mid - 1
      else ! x > bins(mid+1)
        low = mid + 1
      end if

    end do

    ibin = mid
  end subroutine map_to_bin

  subroutine swap_grw_real2darray(me,i,j,array_width)
    type( grw_real2darray_type ), intent(inout) :: me
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: array_width
    ! -------------------------------- !
    real(kind=rk) :: tmp(array_width)

    tmp = me%val(1:array_width, i)
    me%val(1:array_width, i) = me%val(1:array_width, j)
    me%val(1:array_width, j) = tmp

  end subroutine swap_grw_real2darray

  subroutine swap_grw_realarray(me,i,j)
    type( grw_realarray_type ), intent(inout) :: me
    integer, intent(in) :: i
    integer, intent(in) :: j
    ! -------------------------------- !
    real(kind=rk) :: tmp

    tmp = me%val(i)
    me%val(i) = me%val(j)
    me%val(j) = tmp

  end subroutine swap_grw_realarray


  subroutine remove_elem_grw_real2darray(me, i, array_width)
    type( grw_real2darray_type ), intent(inout) :: me
    integer, intent(in) :: i
    integer, intent(in) :: array_width
    ! -------------------------------- !
    integer :: Nvals
    ! -------------------------------- !
    Nvals = me%nvals

    ! First put the value to be "deleted"
    ! (it is not actually deleted, but looks that way)
    ! at the end.
    call swap_grw_real2darray( me          = me,         &
                             & i           = i,          &
                             & j           = Nvals,      &
                             & array_width = array_width )

    ! Now decrease Nvals so it seems as if the element has disappeared
    me%nvals = me%nvals - 1

  end subroutine remove_elem_grw_real2darray


  subroutine remove_elem_grw_realarray(me, i)
    type( grw_realarray_type ), intent(inout) :: me
    integer, intent(in) :: i
    ! -------------------------------- !
    integer :: Nvals
    ! -------------------------------- !
    Nvals = me%nvals

    ! First put the value to be "deleted"
    ! (it is not actually deleted, but looks that way)
    ! at the end.
    call swap_grw_realarray( me          = me,         &
                           & i           = i,          &
                           & j           = Nvals       )

    ! Now decrease Nvals so it seems as if the element has disappeared
    me%nvals = me%nvals - 1

  end subroutine remove_elem_grw_realarray

  !> Fill_cylinder fills the dynamic array positions with coordinates
  !! of particles spaced a distance d apart inside a cylinder with
  !! faces z = 0 and z = L and radius R.
  !! For cylinders of different position and orientation, these coordinates
  !! can be transformed.
  subroutine fill_cylinder(R, L, d, positions)
    !> Radius of the cylinder
    real(kind=rk), intent(in) :: R
    !> Length of the cylinder
    real(kind=rk), intent(in) :: L
    !> Distance between particles
    real(kind=rk) :: d
    !> Output: coordinates of the particles
    type(grw_real2darray_type), intent(inout) :: positions
    ! -------------------------------------------- !
    real(kind=rk) :: x, y, z, r2
    integer :: Nr, Nz  ! number of particles per length 2R
    integer :: ix, iy, iz
    real(kind=rk) :: xs, ys, zs
    real(kind=rk) :: pos(6) ! x, y, z, rx, ry, rz
    ! -------------------------------------------- !
    call init( me = positions, width = 6 )

    ! Determine the number of particles per side
    Nr = floor(2*R/d)
    Nz = floor(L/d)

    ! Coordinates of the bottom-left corner of the circle's bounding square
    xs = -R + 0.5*d
    ys = -R + 0.5*d
    zs = 0.5*d

    ! Create a meshgrid of x and y coordinates
    do iz = 1, Nz
      z = zs + (iz-1)*d
      do ix = 1, Nr
        x = xs + (ix-1)*d
        do iy = 1, Nr
          y = ys + (iy-1)*d

          ! Compute distance of this point from circle origin
          r2 = x**2 + y**2

          if( r2 < R**2 ) then
            ! Append this point to the array of particle coordinates
            pos = 0.0_rk
            pos(1) = x
            pos(2) = y
            pos(3) = z
            call append( me = positions, val =  pos(1:6) )

          end if
        end do ! iy
      end do ! ix
    end do ! iz
  end subroutine fill_cylinder

  subroutine fill_prism_by_distance(lx, ly, lz, d, positions)
    !> Vector spanning the prism in x-direction
    real(kind=rk), intent(in) :: lx
    !> Vector spanning the prism in y-direction
    real(kind=rk), intent(in) :: ly
    !> Vector spanning the prism in z-direction
    real(kind=rk), intent(in) :: lz
    !> Distance between particles
    real(kind=rk), intent(in) :: d
    !> Output: coordinates of the particles
    type(grw_real2darray_type), intent(inout) :: positions
    ! ------------------------------------ !
    integer :: Nx, Ny, Nz
    integer :: ix, iy, iz
    real(kind=rk) :: pos(6)
    real(kind=rk) :: x, y, z
    real(kind=rk) :: xs, ys, zs
    ! ------------------------------------ !
    call init( me = positions, width = 6 )
    Nx = floor(lx/d)
    Ny = floor(ly/d)
    Nz = floor(lz/d)

    ! Initialize pos in a corner of the prism
    xs = 0.5*d
    ys = 0.5*d
    zs = 0.5*d

    do ix = 1, Nx
      x = xs + (ix-1)*d
      do iy = 1, Ny
        y = ys + (iy-1)*d
        do iz = 1, Nz
          z = zs + (iz-1)*d
          pos = 0.0_rk
          pos(1) = x
          pos(2) = y
          pos(3) = z
          call append( me = positions, val = pos(1:6) )
        end do
      end do
    end do

  end subroutine fill_prism_by_distance

  subroutine fill_prism_by_number(lx, ly, lz, Nx, Ny, Nz, positions)
    !> Vector spanning the prism in x-direction
    real(kind=rk), intent(in) :: lx
    !> Vector spanning the prism in y-direction
    real(kind=rk), intent(in) :: ly
    !> Vector spanning the prism in z-direction
    real(kind=rk), intent(in) :: lz
    !> Number of particles in x-direction
    integer, intent(in) :: Nx
    !> Number of particles in y-direction
    integer, intent(in) :: Ny
    !> Number of particles in z-direction
    integer, intent(in) :: Nz
    !> Output: coordinates of the particles
    type(grw_real2darray_type), intent(inout) :: positions
    ! ------------------------------------ !
    integer :: ix, iy, iz
    real(kind=rk) :: pos(6)
    real(kind=rk) :: x, y, z
    real(kind=rk) :: dx, dy, dz
    real(kind=rk) :: xs, ys, zs
    ! ------------------------------------ !
    call init( me = positions, width = 6 )
    dx = lx/Nx
    dy = ly/Ny
    dz = lz/Nz

    ! Initialize pos in a corner of the prism
    xs = 0.5*dx
    ys = 0.5*dy
    zs = 0.5*dz

    do ix = 1, Nx
      x = xs + (ix-1)*dx
      do iy = 1, Ny
        y = ys + (iy-1)*dy
        do iz = 1, Nz
          z = zs + (iz-1)*dz
          pos = 0.0_rk
          pos(1) = x
          pos(2) = y
          pos(3) = z
          call append( me = positions, val = pos(1:6) )
        end do
      end do
    end do

  end subroutine fill_prism_by_number

  function rodriguez_rotation( n1, n2, v ) result(v_rot)
    !> Unit vector describing the orientation to rotate from
    real(kind=rk), intent(in) :: n1(3)
    !> Unit vector describing orientation to rotate to
    real(kind=rk), intent(in) :: n2(3)
    !> Vector to rotate
    real(kind=rk), intent(in) :: v(3)
    !> Output: rotated vector
    real(kind=rk) :: v_rot(3)
    ! ------------------------------------- !
    real(kind=rk) :: n2_x_n1(3) ! cross product n2 x n1
    real(kind=rk) :: k_x_v(3) ! cross product n2 x n1
    real(kind=rk) :: k(3) ! unit vector describing axis of rotation
    real(kind=rk) :: c, s  ! cosine and sine of angle between n1 and n2
    ! ----------------------------------------- !
    c = dot_product(n1,n2)

    ! Compute the cross product which defines the axis of rotation
    call cross_product(n1, n2, n2_x_n1)

    ! Use cross product to determine the sine of the angle between n1 and n2
    s = dot_product(n2_x_n1, n2_x_n1)
    s = sqrt(s)

    ! Calculate axis of rotation
    k = n2_x_n1/s
    call cross_product(k, v, k_x_v)

    ! Compute rotated vector
    v_rot = v*c + (k_x_v)*s + k*( dot_product(k, v) )*(1 - c)

  end function rodriguez_rotation

  subroutine rotate_positions( positions, n1, n2 )
    !> Output: coordinates of the particles
    type(grw_real2darray_type), intent(inout) :: positions
    !> Unit vector describing the input cylinder axis
    real(kind=rk), intent(in) :: n1(3)
    !> Unit vector describing the output cylinder axis.
    real(kind=rk), intent(in) :: n2(3)
    ! ----------------------------------------- !
    integer :: iParticle
    real(kind=rk) :: xp_rot(3)
    ! ----------------------------------------- !
    ! Loop over particles in positions array
    do iParticle = 1, positions%nvals
      xp_rot = rodriguez_rotation( n1 = n1,                           &
                                 & n2 = n2,                           &
                                 & v  = positions%val(1:3, iParticle) )
      positions%val(1:3, iParticle) = xp_rot
    end do ! iParticle

  end subroutine rotate_positions


  subroutine translate_positions( positions, translation_vec )
    !> Output: coordinates of the particles
    type(grw_real2darray_type), intent(inout) :: positions
    !> Vector with which to translate positions by
    real(kind=rk), intent(in) :: translation_vec(3)
    ! ----------------------------------------- !
    integer :: iParticle
    ! ----------------------------------------- !
    ! Loop over particles in positions array
    do iParticle = 1, positions%nvals
      positions%val(1:3, iParticle) = positions%val(1:3, iParticle) &
       &                              + translation_vec(1:3)
    end do ! iParticle

  end subroutine translate_positions

  subroutine print_positions(positions, logUnit)
    !> Input: coordinates of the particles
    type(grw_real2darray_type), intent(in) :: positions
    !> logUnit: if this is a file, it must be opened prior to calling this routine
    integer, intent(in) :: logUnit
    ! --------------------------------------- !
    integer :: iParticle, k
    ! --------------------------------------- !
    do iParticle = 1, positions%nvals
      do k = 1, 6
        write(logUnit, '(E17.9)', advance='no') positions%val(k, iParticle)
      end do
      write(logUnit,'(A)')
    end do

  end subroutine print_positions

  subroutine print_particleblob(particleblob, logUnit)
    type(mus_particle_blob_cylinder_type), intent(in) :: particleblob
    integer, intent(in) :: logUnit
    ! --------------------------------------------- !
    write(logUnit,'(A)') '----- PARTICLE BLOB PROPERTIES -----'
    write(logUnit,'(A,3E17.9)') 'origin = ', particleblob%origin(1:3)
    write(logUnit,'(A,3E17.9)') 'vec    = ', particleblob%vec(1:3)
    write(logUnit,'(A,E17.9)') 'radius = ', particleblob%radius
    write(logUnit,'(A,E17.9)') 'particle radius = ', particleblob%particle_radius
    write(logUnit,'(A,E17.9)') 'particle mass   = ', particleblob%particle_mass
    write(logUnit,'(A,6E17.9)') 'particle_vel   = ', particleblob%particle_vel(1:6)
    write(logUnit,'(A,6E17.9)') 'particle_force = ', particleblob%particle_force(1:6)


  end subroutine print_particleblob

  subroutine print_particleblob_prism(particleblob, logUnit)
    type(mus_particle_blob_prism_type), intent(in) :: particleblob
    integer, intent(in) :: logUnit
    ! --------------------------------------------- !
    write(logUnit,'(A)') '----- PARTICLE BLOB PROPERTIES -----'
    write(logUnit,'(A,3E17.9)') 'origin = ', particleblob%origin(1:3)
    write(logUnit,'(A,3E17.9)') 'vec_x  = ', particleblob%vec_x(1:3)
    write(logUnit,'(A,3E17.9)') 'vec_y  = ', particleblob%vec_y(1:3)
    write(logUnit,'(A,3E17.9)') 'vec_z  = ', particleblob%vec_z(1:3)
    write(logUnit,'(A,E17.9)') 'particle radius = ', particleblob%particle_radius
    write(logUnit,'(A,E17.9)') 'particle mass   = ', particleblob%particle_mass
    write(logUnit,'(A,6E17.9)') 'particle_vel   = ', particleblob%particle_vel(1:6)
    write(logUnit,'(A,6E17.9)') 'particle_force = ', particleblob%particle_force(1:6)


  end subroutine print_particleblob_prism

  subroutine print_particleblob_prob(particleblob, logUnit)
    type(mus_particle_blob_prob_type), intent(in) :: particleblob
    integer, intent(in) :: logUnit
    ! --------------------------------------------- !
    integer :: i
    ! --------------------------------------------- !
    write(logUnit,'(A)') '----- PARTICLE BLOB AVAILABLE POSITIONS -----'
    do i = 1, particleblob%available_positions%nvals
      write(logUnit,'(2E17.9)') particleblob%available_positions%val(1:2,i)
    end do

    write(logUnit,'(A)') '----- PARTICLE BLOB POSITION PROBABILITIES -----'
    do i = 1, particleblob%probabilities%nvals
      write(logUnit,'(E17.9)') particleblob%probabilities%val(i)
    end do

    write(logUnit,'(A)') '----- PARTICLE BLOB BINS -----'
    do i = 1, particleblob%bins%nvals
      write(logUnit,'(E17.9)') particleblob%bins%val(i)
    end do

    write(logUnit,'(A)') '----- PARTICLE BLOB CHOSEN POSITIONS -----'
    do i = 1, particleblob%chosen_positions%nvals
      write(logUnit,'(2E17.9)') particleblob%chosen_positions%val(1:2,i)
    end do

    write(logUnit,'(A)') '----- PARTICLE BLOB SUM OF PROBABILITIES -----'
    write(logUnit,'(E17.9)') sum(particleblob%probabilities%val(1:particleblob%probabilities%nvals ))

  end subroutine print_particleblob_prob

end module mus_particle_blob_module
