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
!> Module containing some auxiliary routines for LBM-DEM 
!! simulations for particles in a flow.

module mus_particle_aux_module

  use mpi

  use env_module,              only: rk, long_k
  use tem_geometry_module,     only: tem_PosOfId
  use tem_global_module,       only: tem_global_type
  use tem_propHead_module,     only: tem_propHead_type 
  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit
  use tem_stencil_module,      only: tem_stencil_findIndexOfDir
  use tem_property_module,     only: prp_hasBnd, prp_hasRemoteNgh, &
    &                                prp_sendHalo
  use tem_topology_module,     only: tem_firstIdAtLevel,        &
    &                                tem_coordOfId, tem_IdOfCoord
  use tem_geometry_module,     only: tem_CoordOfReal
  use tem_construction_module, only: tem_levelNeighbor_type
  use mus_geom_module,         only: mus_geom_type
  use mus_scheme_type_module,  only: mus_scheme_type

  implicit none
  private

  public :: findPartitionOfTreeID
  public :: coordExistsOnMyRank
  public :: coordLocalOnMyRank
  public :: cross_product
  public :: getBaryOfCoord
  public :: getCartesianStencilIndices
  public :: getPathToVec
  public :: getPosOfCoord
  public :: followNeighPath
  public :: findPropIndex
  public :: getProcessBoundingBox
  public :: mus_particles_global_errorcheck
  public :: positionLocalOnMyRank


contains

  ! ************************************************************************** !
  !> Routine to find the process a certain treeID is on.
  !! Looks in the ranks in the array of argument procs
  !! Returns the INDEX of the found proc in this array
  !! Example: procs = [4 8 7] and we find our treeID on rank 8
  !!          will return procIndex = 2
  subroutine findPartitionOfTreeID( sTreeID, geometry, procs, nProcs, &
                                  & outProc, procIndex                )
    ! ---------------------------------------------------------------------- !
    !> TreeID to look for
    integer(kind=long_k), intent(in) :: sTreeID
    !> Geometry to look in
    type(mus_geom_type), intent(in) :: geometry
    !> Only look on these ranks (look in all ranks if omitted)
    integer, optional, intent(in) :: procs(:)
    !> Number of elements in procs
    integer, intent(in) :: nProcs
    !> Process which the treeID is on
    integer, intent(out) :: outProc
    !> Index of that process
    integer, intent(out) :: procIndex
    ! ---------------------------------------------------------------------- !
    integer :: iproc, proc
    ! ---------------------------------------------------------------------- !
    procIndex = -1
    outProc = -1

    do iproc = 1, nProcs
      if( present(procs) ) then
        proc = procs(iproc)
      else
        proc = iproc-1
      end if

      ! NOTE: the index of proc in Part_First and Part_Last
      !       will be proc + 1 since the indices start at 1
      !       whereas process numbering starts at 0
      if( sTreeID > geometry%tree%Part_Last( proc+1 ) ) then
        cycle
      else if( sTreeID >= geometry%tree%Part_First( proc + 1) ) then
        ! TreeID is between first and last TreeID on partition iproc - 1
        ! So we found the partition we're looking for
        outProc = proc
        procIndex = iproc
        return
      end if
    end do ! iproc

  end subroutine findPartitionOfTreeID
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Routine which checks whether integer coordinate is local 
  !! on the current rank
  function coordLocalOnMyRank( coord, scheme, myRank )
    ! ---------------------------------------------------------------------- !
    !> integer coordinate (x,y,z,level)
    integer, intent(in) :: coord(4)
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> This process rank
    integer, intent(in) :: myRank
    ! ---------------------------------------------------------------------- !
    logical :: coordLocalOnMyRank
    integer :: lev
    integer(kind=long_k) :: TreeID
    integer :: ldPos
    ! ---------------------------------------------------------------------- !
    lev = coord(4)
    TreeID = tem_IdOfCoord(coord = coord)

    ldPos = tem_PosOfId( sTreeID    = TreeID,                      &
      &                  treeIDlist = scheme%levelDesc(lev)%total, &
      &                  lower      = 1,                           &
      &                  upper      = scheme%pdf(lev)%nElems_fluid )
    if (ldPos > 0) then
      coordLocalOnMyRank = .TRUE.
    else
      coordLocalOnMyRank = .FALSE.
    end if

  end function coordLocalOnMyRank
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Routine which checks whether a spatial position is local 
  !! on the current rank
  function positionLocalOnMyRank( pos, geometry, scheme, myRank )
    ! ---------------------------------------------------------------------- !
    !> Cartesian position (x, y, z)
    real(kind=rk), intent(in) :: pos(3)
    !> Geometry to look in
    type(mus_geom_type), intent(in) :: geometry
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> This process rank
    integer, intent(in) :: myRank
    ! ---------------------------------------------------------------------- !
    integer :: coord(4), lev
    logical :: positionLocalOnMyRank
    ! ---------------------------------------------------------------------- !
    lev = geometry%tree%global%maxLevel

    ! Get integer coordinate of this position
    coord = tem_coordOfReal( mesh  = geometry%tree,   &
      &                      point = pos,             &
      &                      level = lev              )

    positionLocalOnMyRank = coordLocalOnMyRank( coord = coord,   &
      &                                         scheme = scheme, &
      &                                         myRank = myRank  )
  end function positionLocalOnMyRank
  ! ************************************************************************** !

  ! ************************************************************************** !
  function treeIDlocalOnMyRank( TreeID, scheme, lev )
    ! ---------------------------------------------------------------------- !
    !> TreeID to check
    integer(kind=long_k), intent(in) :: TreeID
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> Level
    integer, intent(in) :: lev
    !> Output: TRUE if this TreeID is a local fluid element on my rank
    logical :: treeIDlocalOnMyRank
    ! ---------------------------------------------------------------------- !
    integer :: ldPos
    ! ---------------------------------------------------------------------- !
    ldPos = tem_PosOfId( sTreeID    = TreeID,                      &
      &                  treeIDlist = scheme%levelDesc(lev)%total, &
      &                  lower      = 1,                           &
      &                  upper      = scheme%pdf(lev)%nElems_fluid )
    if(ldPos > 0) then
      TreeIDlocalOnMyRank = .TRUE.
    else
      TreeIDlocalOnMyRank = .FALSE.
    end if
  end function treeIDlocalOnMyRank
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Function which checks if element has remote neighbors
  logical function has_remote_neighbors( iElem, scheme, lev )
    ! ---------------------------------------------------------------------- !
    !> Index of this element in the total list
    integer, intent(in) :: iElem
    !> Scheme for access to levelDescriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> Level
    integer, intent(in) :: lev
    ! ---------------------------------------------------------------------- !
    integer :: iDir
    integer :: iNeigh
    integer(kind=long_k) :: TreeID
    logical :: elementIsLocal
    ! ---------------------------------------------------------------------- !
    ! Initially set has_remote_neighbors to .FALSE.
    has_remote_neighbors = .FALSE.

    do iDir = 1, scheme%layout%fStencil%QQN
      ! Get the position of neighbor element in this direction
      ! in total list
      iNeigh = scheme%levelDesc(lev)%neigh(1)%nghElems(iDir, iElem)
      if(iNeigh <= 0) then
        ! iNeigh <= 0 indicates the neighbor of this element is a boundary. 
        ! In that case we also do not want to place a particle here.
        has_remote_neighbors = .TRUE.
        return
      end if

      TreeID = scheme%levelDesc(lev)%total(iNeigh)

      ! Check whether this neighbor element is a local fluid
      elementIsLocal = treeIDlocalOnMyRank( TreeID = TreeID, &
        &                                   scheme = scheme, &
        &                                   lev    = lev     )

      ! As soon as we find one non-local neighbor we 
      ! can set has_remote_neighbors to true
      if( .NOT. elementIsLocal ) then
        has_remote_neighbors = .TRUE.
        return
      end if
    end do ! iDir
  end function has_remote_neighbors
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Checks if integer coordinate is in the total list of this rank
  function coordExistsOnMyRank( coord, scheme, myRank )
    ! ---------------------------------------------------------------------- !
    !> integer coordinate (x,y,z,level)
    integer, intent(in) :: coord(4)
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> This process rank
    integer, intent(in) :: myRank
    ! ---------------------------------------------------------------------- !
    logical :: coordExistsOnMyRank
    integer :: lev
    integer(kind=long_k) :: TreeID
    integer :: ldPos
    ! ---------------------------------------------------------------------- !
    lev = coord(4)
    TreeID = tem_IdOfCoord(coord = coord)

    ldPos = tem_PosOfId( sTreeID    = TreeID,                      &
      &                  treeIDlist = scheme%levelDesc(lev)%total, &
      &                  lower      = 1,                           &
      &                  upper      = scheme%pdf(lev)%nElems_local )
    if(ldPos > 0) then
      coordExistsOnMyRank = .TRUE.
    else
      coordExistsOnMyRank = .FALSE.
    end if

  end function coordExistsOnMyRank
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> computes cross-product a x b = res
  subroutine cross_product(a,b,res)
    ! ---------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: a(3)
    real(kind=rk), intent(in) :: b(3)
    real(kind=rk), intent(out) :: res(3)
    ! ---------------------------------------------------------------------- !

    res(1) = a(2) * b(3) - a(3) * b(2)
    res(2) = -a(1) * b(3) + a(3) * b(1)
    res(3) = a(1) * b(2) - a(2) * b(1)

  end subroutine cross_product
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Convenience function to get barycenter of integer coordinate
  pure function getBaryOfCoord( coord, origin, dx ) result(bary)
    ! ---------------------------------------------------------------------- !
    integer, intent(in) :: coord(3)
    real(kind=rk), intent(in) :: origin(3)
    real(kind=rk), intent(in) :: dx
    real(kind=rk) :: bary(3)
    ! ---------------------------------------------------------------------- !
    real(kind=rk) :: c(3)
    ! ---------------------------------------------------------------------- !
    c(:) = origin(:) + 0.5_rk * dx
    bary(1) = c(1) + real(coord(1), kind=rk) * dx
    bary(2) = c(2) + real(coord(2), kind=rk) * dx
    bary(3) = c(3) + real(coord(3), kind=rk) * dx

  end function getBaryOfCoord
  ! ************************************************************************** !

  ! ************************************************************************** !
  subroutine getCartesianStencilIndices( cxDir, cartIndices )
    ! ---------------------------------------------------------------------- !
    !> Integer stencil direction
    integer, intent(in) :: cxDir(:,:)
    !> Output indices of the Cartesian directions
    !! [-x +x -y +y -z +z] in cxDir
    integer, intent(out) :: cartIndices(6)
    ! ---------------------------------------------------------------------- !
    cartIndices(1) = tem_stencil_findIndexOfDir( findDir = (/ -1, 0, 0 /), &
      &                                          cxDir   = cxDir           )

    cartIndices(2) = tem_stencil_findIndexOfDir( findDir =  (/ 1, 0, 0 /), &
      &                                          cxDir   = cxDir           )

    cartIndices(3) = tem_stencil_findIndexOfDir( findDir =  (/ 0, -1, 0 /), &
      &                                          cxDir   = cxDir            )

    cartIndices(4) = tem_stencil_findIndexOfDir( findDir =  (/ 0, 1, 0 /), &
      &                                          cxDir   = cxDir           )

    cartIndices(5) = tem_stencil_findIndexOfDir( findDir =  (/ 0, 0, -1 /), &
      &                                          cxDir   = cxDir            )

    cartIndices(6) = tem_stencil_findIndexOfDir( findDir =  (/ 0, 0, 1 /), &
      &                                          cxDir   = cxDir           )

    ! if( any( cartIndices < 0 ) ) then
    !   write(logUnit(1),*) "ERROR getCartesianStencilIndices: could not find dirs"
    !   call tem_abort()
    ! end if

  end subroutine getCartesianStencilIndices
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> getPathToVec builds a series of 3 indices which can be followed using the 
  !! levelDesc(lev)%neigh array to get the position of a neighboring element 
  !! in the total list.
  !! Usage: given an element at coordinate [0,0,0] with position startPos in the
  !! total list, to reach the element at [1,1,1] we do:
  !! call getPathToVec([1,1,1], cxDir, path)
  !! pos = startPos 
  !! do k = 1,3
  !!   pos = levelDesc(lev)%neigh%nghElems(path(k),pos)
  !! end do
  !! At the end of the loop pos will point to the requested element
  subroutine getPathToVec( vec, cxDir, path )
    ! ---------------------------------------------------------------------- !
    !> Vector we'd like to build the path to
    integer, intent(in) :: vec(3)
    !> Integer stencil direction
    integer, intent(in) :: cxDir(:,:)
    !> Output path of iDirs
    integer, intent(out) :: path(3)
    ! ---------------------------------------------------------------------- !
    integer :: uvec(3)
    integer :: k
    ! ---------------------------------------------------------------------- !
    ! Break vec down into unit vectors and construct path of iDir's 
    ! that need to be followed to reach it
    do k = 1,3 
      uvec = 0
      uvec(k) = vec(k)
      path(k) = tem_stencil_findIndexOfDir( findDir = uvec,  &
        &                                   cxDir   = cxDir  )
    end do
  end subroutine getPathToVec
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Get the position in the levelDesc%total list of an element with 
  !! integer coordinate coord. This routine returns ldPos <= 0 if 
  !! the element could not be found on this process.
  function getPosOfCoord(coord, scheme) result(ldPos)
    ! ---------------------------------------------------------------------- !
    integer, intent(in) :: coord(4)
    type(mus_scheme_type), intent(in) :: scheme
    integer :: ldPos
    ! ---------------------------------------------------------------------- !
    integer(kind=long_k) :: TreeID
    integer :: lev
    ! ---------------------------------------------------------------------- !
    lev = coord(4)

    ! get TreeID and position in complete tree
    TreeID = tem_IdOfCoord(coord = coord)

    ldPos = tem_PosOfId( sTreeID    = TreeID,                      &
      &                  treeIDlist = scheme%levelDesc(lev)%total, &
      &                  lower      = 1,                           &
      &                  upper      = scheme%pdf(lev)%nElems_fluid )
    if( ldPos <= 0) then
      ! If we could not find ldPos in local elems, look in halos
      ldPos = tem_PosOfId( sTreeID    = TreeID,                           &
        &                  treeIDlist = scheme%levelDesc(lev)%total,      &
        &                  lower      = scheme%pdf(lev)%nElems_fluid + 1, &
        &                  upper      = scheme%pdf(lev)%nElems_local      )
    end if

  end function getPosOfCoord
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> followNeighPath can be used to determine the position of a neighboring
  !! lattice site in the total list. For this we follow one or more entries of
  !! the levelDesc%neigh array. Since this array is only defined for local fluid
  !! elements, startPos must point to a local fluid element.
  function followNeighPath( neighPath, startPos, neigh, nElems_fluid, &
    &                       zeroDir ) result( endPos )
    ! ---------------------------------------------------------------------- !
    !> Array containing the sequence of indices of directions that must be
    !! followed to get to the direction described by neighDir. This depends on
    !! the stencil being used for the simulation. For example in the d3q19
    !! stencil, to reach the neighbor located with displacement (1,1,1) from the
    !! element at startPos, we follow first (1,0,0) and then (0,1,1). For d3q27
    !! neighPath will only be one entry (with value 26) corresponding to the
    !! index (1,1,1).
    integer, intent(in) :: neighPath(:)
    !> Position of the starting element in the total list
    !! This element must be a local fluid!
    integer, intent(in) :: startPos
    !> Neighbor list
    type( tem_levelNeighbor_type ), intent(in) :: neigh
    !> Number of local fluid elements on this process
    !! We need this to check if intermediate values of endPos 
    !! are still local fluids
    integer, intent(in) :: nElems_fluid
    !> Index of the "zero direction" of the stencil i.e. 
    !! stencil%cxdir(1:3,zeroDir) = [0,0,0]
    integer, intent(in) :: zeroDir
    !> Output: position of the neighbor element located by following 
    !! consecutively neighPath(1), neighPath(2) etc.
    integer :: endPos
    ! ---------------------------------------------------------------------- !
    integer :: i, iDir, N
    ! ---------------------------------------------------------------------- !
    endPos = startPos
    N = ubound(neighPath, 1)

    do i = 1, N 
      iDir = neighPath(i)

      ! If the iDir in this step of neighPath corresponds to [0,0,0],
      ! we do not need to do anything so skip to next iteration.
      if(iDir == zeroDir) cycle

      ! NOTE: this will break if at any point endPos points at a halo
      ! element. In this case abort the routine and set endPos = -1 to 
      ! indicate that execution was unsuccesful.
      if(endPos > nElems_fluid .OR. endPos < 1 ) then
        endPos = -1
        return
      end if

      endPos = neigh%nghElems(iDir, endPos) 
    end do

  end function followNeighPath
  ! ************************************************************************** !

  ! ************************************************************************** !
  ! Function to find the index of a property with prp_bitpos
  function findPropIndex(prp_bitpos, glob) result(propIndex)
    ! ---------------------------------------------------------------------- !
    !> Property bit of property to find index for
    integer, intent(in) :: prp_bitpos
    !> tree%global mesh information
    type(tem_global_type), intent(in) :: glob
    ! Output index of property with prp_bitpos in propHead array
    integer :: propIndex
    ! ---------------------------------------------------------------------- !
    integer iProp
    ! ---------------------------------------------------------------------- !
    propIndex = -1
    write(logUnit(1),*) 'findPropIndex glob%nproperties = ', glob%nProperties
    do iProp = 1, glob%nProperties 
      write(logUnit(1),*) 'glob%property(iProp) = ', glob%property(iProp)%label
      if( glob%property( iProp )%bitpos == prp_bitpos ) then
        propIndex = iProp
      end if
    end do

  end function findPropIndex
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Routine to find the bounding box of the part of the spatial 
  !! domain on this process
  subroutine getProcessBoundingBox( scheme, geometry, boundingBox )
    ! ---------------------------------------------------------------------- !
    !> Scheme for access to level descriptor
    type(mus_scheme_type), intent(in) :: scheme
    !> Geometry information to determine TreeIDs of elements 'covered' by 
    !! particle
    type(mus_geom_type), intent(in) :: geometry
    !> Output: bounding box (xmin, xmax, ymin, ymax, zmin, zmax)
    real(kind=rk), intent(out) :: boundingBox(6)
    ! ---------------------------------------------------------------------- !
    integer :: intBoundingBox(6)
    integer :: lev, upperBound
    integer :: iElem, nFluids
    integer :: elemCoord(4), coord(4)
    integer(kind=long_k) :: TIDoffset
    integer :: nProcs
    real(kind=rk) :: dx
    real(kind=rk) :: coordX(3)
    ! ---------------------------------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    upperBound = 2**lev
    TIDoffset = tem_firstIdAtLevel(lev)
    dx = geometry%tree%global%BoundingCubeLength / 2**lev
    nFluids = scheme%pdf(lev)%nElems_fluid

    ! Check if we are running the serial case
    nProcs = geometry%tree%global%nParts
    if(nProcs < 2) then
      boundingBox(1) = geometry%tree%global%Origin(1)
      boundingBox(2) = geometry%tree%global%Origin(1) &
        &            + geometry%tree%global%BoundingCubeLength
      boundingBox(3) = geometry%tree%global%Origin(2)
      boundingBox(4) = geometry%tree%global%Origin(2) &
        &            + geometry%tree%global%BoundingCubeLength
      boundingBox(5) = geometry%tree%global%Origin(3)
      boundingBox(6) = geometry%tree%global%Origin(3) &
        &            + geometry%tree%global%BoundingCubeLength
    else
      elemCoord = tem_coordOfID( TreeID = scheme%levelDesc(lev)%total(1), &
        &                        offset = TIDoffset                       )

      intBoundingBox(1:2) = elemCoord(1)
      intBoundingBox(3:4) = elemCoord(2)
      intBoundingBox(5:6) = elemCoord(3)

      ! Loop over all local fluid elements in this process
      do iElem = 1, nFluids
        ! Check if elem has prp_sendhalo
        if( btest( scheme%levelDesc(lev)%property(iElem), prp_sendHalo) .OR. &
          & btest( scheme%levelDesc(lev)%property(iElem), prp_hasBnd )  ) then
          ! For elems with prp_sendhalo, we will look to the surrounding
          ! elements to find which remote processes they belong to

          ! Get the searchBoxLimits
          elemCoord = tem_coordOfID(                                 &
            &           TreeID = scheme%levelDesc(lev)%total(iElem), &
            &           offset = TIDoffset                           )
          if( any(elemCoord(1:3) < 0)                    &
            & .OR. any(elemCoord(1:3) > upperBound) ) then
            cycle
          else ! elemCoordinate is within domain
            if( elemCoord(1) < intBoundingBox(1) ) then
              intBoundingBox(1) = elemCoord(1)
            end if
            if( elemCoord(2) < intBoundingBox(3) ) then
              intBoundingBox(3) = elemCoord(2)
            end if
            if( elemCoord(3) < intBoundingBox(5) ) then
              intBoundingBox(5) = elemCoord(3)
            end if
            if( elemCoord(1) > intBoundingBox(2) ) then
              intBoundingBox(2) = elemCoord(1)
            end if
            if( elemCoord(2) > intBoundingBox(4) ) then
              intBoundingBox(4) = elemCoord(2)
            end if
            if( elemCoord(3) > intBoundingBox(6) ) then
              intBoundingBox(6) = elemCoord(3)
            end if
          end if
        end if ! elem has prp_sendHalo

      end do ! iElem

      ! Convert intBoundingBox integer coordinates to cartesian and write to log
      coord(1) = intBoundingBox(1)
      coord(2) = intBoundingBox(3)
      coord(3) = intBoundingBox(5)
      coordX = getBaryOfCoord( coord  = coord(1:3),                  &
        &                      origin = geometry%tree%global%origin, &
        &                      dx     = dx                           )

      boundingBox(1) = coordX(1)
      boundingBox(3) = coordX(2)
      boundingBox(5) = coordX(3)

      coord(1) = intBoundingBox(2)
      coord(2) = intBoundingBox(4)
      coord(3) = intBoundingBox(6)
      coordX = getBaryOfCoord( coord  = coord(1:3),                  &
        &                      origin = geometry%tree%global%origin, &
        &                      dx     = dx                           )

      boundingBox(2) = coordX(1)
      boundingBox(4) = coordX(2)
      boundingBox(6) = coordX(3)

    end if

  end subroutine getProcessBoundingBox
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Routine to aggregate local error flags and use MPI allreduce 
  !! to set a global flag if any of the local flags is true
  subroutine mus_particles_global_errorcheck( flag, comm )
    ! ---------------------------------------------------------------------- !
    !> Error flag: at input this will be the local error code of 
    !! this rank only which is set to TRUE if there is an error and 
    !! FALSE if not. Upon output this will be set to flag = TRUE if 
    !! ANY of the local error codes (on all processes in comm) were TRUE.
    logical, intent(inout) :: flag
    !> MPI communicator
    integer, intent(in) :: comm
    ! ---------------------------------------------------------------------- !
    integer :: iErr
    logical :: flag_local, flag_global
    ! ---------------------------------------------------------------------- !
    flag_local = flag

    ! Use MPI reduce to check if error code on ANY process was TRUE
    call MPI_allreduce( &
      &    flag_local,  & ! send buffer
      &    flag_global, & ! receive buffer
      &    1,           & ! count
      &    MPI_LOGICAL, & ! datatype
      &    MPI_LOR,     & ! reduction operation
      &    comm,        & ! MPI communicator
      &    iErr         ) ! iError

    ! Set flag equal to the global error code.
    flag = flag_global

  end subroutine mus_particles_global_errorcheck
  ! ************************************************************************** !

end module mus_particle_aux_module
