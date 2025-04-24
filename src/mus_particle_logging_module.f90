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
!> Routines containing routines for logging particle data.

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

module mus_particle_logging_module
  use env_module,                       only: rk, long_k, newunit, labelLen
  use tem_aux_module,                   only: tem_abort
  use tem_logging_module,               only: logUnit
  use tem_dyn_array_module,             only: init, append, destroy,    &
    &                                         empty, dyn_intArray_type, &
    &                                         dyn_longArray_type
  use tem_geometry_module,              only: tem_CoordOfReal, tem_PosOfId
  use tem_topology_module,              only: tem_IdOfCoord, tem_FirstIdAtLevel, &
    &                                         tem_coordOfId
  use tem_timer_module,                 only: tem_getMaxTimerVal, &
    &                                         tem_getTimerVal


  use mus_geom_module,                  only: mus_geom_type
  use mus_scheme_type_module,           only: mus_scheme_type
  use mus_param_module,                 only: mus_param_type

  use mus_particle_logging_type_module, only: mus_particle_debugtracking_type
  use mus_particle_type_module,         only: mus_particle_MEM_type, &
    &                                         mus_particle_DPS_type, &
    &                                         mus_particle_group_type
  use mus_particle_aux_module,          only: getBaryOfCoord
  use mus_particle_checks_module,       only: compute_fluid_momentum, compute_particle_momentum
  use mus_particle_timer_module,        only: mus_particle_timerHandles

  implicit none

  interface logParticleData
    module procedure logParticleData_MEM
    module procedure logParticleData_DPS
  end interface

  public :: initParticleLog_DPS
  public :: initParticleLog_MEM
  public :: closeParticleLog


contains


  !> Main routine to log particle data for MEM particles
  subroutine mus_particles_logdata_MEM( particleGroup, params, t )
    !> Particle group to log
    type(mus_particle_group_type), intent(in) :: particleGroup
    !> Params for access to time control
    type(mus_param_type), intent(in) :: params
    !> Time to log. If not provided this will be taken from params
    real(kind=rk), optional :: t
    ! --------------------------------------------!
    integer :: particleLogInterval, iParticle, particleLogUnit
    ! --------------------------------------------!
    particleLogInterval = particleGroup%particleLogInterval

    do iParticle = 1, particleGroup%particles_MEM%nvals
      if( mod(params%general%simcontrol%now%iter, particleLogInterval ) == 0  ) then
        if( particleGroup%particles_MEM%val(iParticle)%owner == params%general%proc%rank ) then

          particleLogUnit = getParticleLogUnit(                                        &
                              & particleGroup%particles_MEM%val(iParticle)%particleID, &
                              & params%general%proc%rank )

          if( present(t) ) then
            call logParticleData( particle = particleGroup%particles_MEM%val(iParticle), &
                                & plogUnit  = particleLogUnit,                           &
                                & myRank   = params%general%proc%rank,                   &
                                & t        = t                                           )
          else
            call logParticleData( particle = particleGroup%particles_MEM%val(iParticle), &
                                & plogUnit  = particleLogUnit,                           &
                                & myRank   = params%general%proc%rank,                   &
                                & t        = params%general%simcontrol%now%sim           )
          end if
          call closeParticleLog(particleLogUnit)
        end if
      end if
    end do
  end subroutine mus_particles_logdata_MEM

  !> Main routine to log particle data for DPS particles
  subroutine mus_particles_logdata_DPS( particleGroup, params )
    !> Particle group to log
    type(mus_particle_group_type), intent(in) :: particleGroup
    !> Params for access to time control
    type(mus_param_type), intent(in) :: params
    ! --------------------------------------------!
    integer :: particleLogInterval, iParticle, particleLogUnit
    ! --------------------------------------------!
    particleLogInterval = particleGroup%particleLogInterval

    do iParticle = 1, particleGroup%particles_DPS%nvals
      if( mod(params%general%simcontrol%now%iter, particleLogInterval ) == 0  ) then
        if( particleGroup%particles_DPS%val(iParticle)%owner == params%general%proc%rank ) then

          particleLogUnit = getParticleLogUnit(                                        &
                              & particleGroup%particles_DPS%val(iParticle)%particleID, &
                              & params%general%proc%rank )

          call logParticleData( particle = particleGroup%particles_DPS%val(iParticle), &
                              & plogUnit = particleLogUnit,                            &
                              & myRank   = params%general%proc%rank,                   &
                              & t        = params%general%simcontrol%now%sim           )
          call closeParticleLog(particleLogUnit)
        end if
      end if
    end do
  end subroutine mus_particles_logdata_DPS

  !> Routine to create a particle logunit based on particle ID
  pure function getParticleLogUnit( particleID, myRank ) result(plogUnit)
      integer, intent(in) :: particleID
      integer, intent(in) :: myRank
      integer :: plogUnit
      ! --------------------------------------------!
      plogUnit = 150 + particleID

  end function getParticleLogUnit

  !> Routine to log MEM particle data
  subroutine logParticleData_MEM( particle, plogUnit, myRank, t )
      type(mus_particle_MEM_type), intent(in) :: particle
      integer, intent(in) :: plogUnit
      integer, intent(in) :: myRank
      real(kind=rk), intent(in)  :: t
      ! --------------------------------------------!
      logical :: fileExists, fileIsOpen
      logical :: logUnitTaken
      integer :: existingLogUnit
      character(len=1024) :: filename
      integer :: i
      ! --------------------------------------------!

      ! The desired file name determined by particleID and rank
      ! write (filename, "(A8,I0.4,A1,I0.4)") "particle", particle%particleID, &
      !   & "p", myRank
      write (filename, "(A8,I0.4)") "particle", particle%particleID

      filename = trim(filename)//'.dat'

      inquire(file=filename, exist=fileExists, opened=fileIsOpen, number=existingLogUnit)

      if(fileIsOpen) then
        if(existingLogUnit /= plogUnit) then
          write(logUnit(1),*) 'ERROR logParticleData: wrong log unit connected'
          call tem_abort()
        end if ! plogUnit does not match
      else ! File is not open
        ! Check if the log unit is available or used for some other file
        inquire(unit=plogUnit, opened=logUnitTaken)
        if(logUnitTaken) then
          ! write(logUnit(1),*) 'ERROR logParticleData: log unit taken'
          call tem_abort()
        else
          ! write(logUnit(1),*) 'logParticleData: initializing particle log'
          call initParticleLog_MEM( particle%particleID, filename, plogUnit, fileExists )
        end if ! logUnitTaken
      end if ! fileIsOpen

      ! After all the checks above, the file should be opened and the correct
      ! logUnit attached
      write(plogUnit, '(E17.9)', advance='no') t
      write(plogUnit, '(A)', advance='no') ' '

      ! Positions
      do i = 1, 3
          write(plogUnit, '(E17.9)', advance='no') particle%pos(i)
          write(plogUnit, '(A)', advance='no') ' '
      end do

      ! Velocities
      do i = 1, 6
          write(plogUnit, '(E17.9)', advance='no') particle%vel(i)
          write(plogUnit, '(A)', advance='no') ' '
      end do

      ! Hydrodynamic forces
      do i = 1, 6
          write(plogUnit, '(E17.9)', advance='no') &
            & ( particle%F(i) + particle%Fext(i) )
          write(plogUnit, '(A)', advance='no') ' '
      end do

      ! Collision forces
      do i = 1, 3
          write(plogUnit, '(E17.9)', advance='no') &
            & ( 0.5*(particle%F_DEM(particle%F_DEM_now,i) + &
            &   particle%F_DEM(particle%F_DEM_next,i) ) &
            & - particle%F(i) - particle%Fext(i) )
          write(plogUnit, '(A)', advance='no') ' '
      end do

      write(plogUnit, '(A)')
  end subroutine logParticleData_MEM

  !> Routine to log DPS particle data
  subroutine logParticleData_DPS( particle, plogUnit, myRank, t )
      type(mus_particle_DPS_type), intent(in) :: particle
      integer, intent(in) :: plogUnit
      integer, intent(in) :: myRank
      real(kind=rk), intent(in)  :: t
      ! --------------------------------------------!
      logical :: fileExists, fileIsOpen
      logical :: logUnitTaken
      integer :: existingLogUnit
      character(len=1024) :: filename
      integer :: i
      ! --------------------------------------------!

      ! The desired file name determined by particleID and rank
      ! write (filename, "(A8,I0.4,A1,I0.4)") "particle", particle%particleID, &
      !   & "p", myRank
      write (filename, "(A8,I0.7)") "particle", particle%particleID

      filename = trim(filename)//'.dat'

      inquire(file=filename, exist=fileExists, opened=fileIsOpen, number=existingLogUnit)

      if(fileIsOpen) then
        if(existingLogUnit /= plogUnit) then
          write(logUnit(1),*) 'ERROR logParticleData: wrong log unit connected'
          call tem_abort()
        end if ! plogUnit does not match
      else ! File is not open
        ! Check if the log unit is available or used for some other file
        inquire(unit=plogUnit, opened=logUnitTaken)
        if(logUnitTaken) then
          ! write(logUnit(1),*) 'ERROR logParticleData: log unit taken'
          call tem_abort()
        else
          ! write(logUnit(1),*) 'logParticleData: initializing particle log'
          call initParticleLog_DPS( particle%particleID, filename, plogUnit, fileExists )
        end if ! logUnitTaken
      end if ! fileIsOpen

      ! After all the checks above, the file should be opened and the correct
      ! plogUnit attached
      write(plogUnit, '(E17.9)', advance='no') t
      write(plogUnit, '(A)', advance='no') ' '

      ! Positions
      do i = 1, 3
          write(plogUnit, '(E17.9)', advance='no') particle%pos(i)
          write(plogUnit, '(A)', advance='no') ' '
      end do

      ! Velocities
      do i = 1, 6
          write(plogUnit, '(E17.9)', advance='no') particle%vel(i)
          write(plogUnit, '(A)', advance='no') ' '
      end do

      ! Forces
      do i = 1, 6
          write(plogUnit, '(E17.9)', advance='no') &
            & ( particle%F_DEM( particle%F_DEM_now, i ) + particle%F(i) )
          write(plogUnit, '(A)', advance='no') ' '
      end do
      ! Average forces
      do i = 1, 3
          write(plogUnit, '(E17.9)', advance='no') &
            & particle%Favg(i)
          write(plogUnit, '(A)', advance='no') ' '
      end do

      write(plogUnit, '(A)')

      close(plogUnit)

  end subroutine logParticleData_DPS

  !> Routine to initialize MEM particle log file (e.g. print header)
  subroutine initParticleLog_MEM( particleID, fileName, plogUnit, fileExists )
    integer, intent(in) :: particleID
    character(*), intent(in) :: fileName
    integer, intent(in) :: plogUnit
    logical, intent(in) :: fileExists
    ! --------------------------------------------!

    ! write(logUnit(1),*) 'Calling initParticleLog!'
    if(fileExists) then
      open(plogUnit, file = filename, status = 'old', position='append')
    else
      ! write(logUnit(1),*) 'INITPARTICLELOG logUnit', logUnit
      open(plogUnit, file = filename, status = 'new')
      write(plogUnit, '(A,I5)') 'Particle ID = ', particleID
      write(plogUnit, '(A18)', advance='no') 't'
      write(plogUnit, '(A18)', advance='no') 'x'
      write(plogUnit, '(A18)', advance='no') 'y'
      write(plogUnit, '(A18)', advance='no') 'z'
      write(plogUnit, '(A18)', advance='no') 'ux'
      write(plogUnit, '(A18)', advance='no') 'uy'
      write(plogUnit, '(A18)', advance='no') 'uz'
      write(plogUnit, '(A18)', advance='no') 'urx'
      write(plogUnit, '(A18)', advance='no') 'ury'
      write(plogUnit, '(A18)', advance='no') 'urz'
      write(plogUnit, '(A18)', advance='no') 'Fx_hydro'
      write(plogUnit, '(A18)', advance='no') 'Fy_hydro'
      write(plogUnit, '(A18)', advance='no') 'Fz_hydro'
      write(plogUnit, '(A18)', advance='no') 'RFx_hydro'
      write(plogUnit, '(A18)', advance='no') 'RFy_hydro'
      write(plogUnit, '(A18)', advance='no') 'RFz_hydro'
      write(plogUnit, '(A18)', advance='no') 'Fx_coll'
      write(plogUnit, '(A18)', advance='no') 'Fy_coll'
      write(plogUnit, '(A18)') 'Fz_coll'
    end if

  end subroutine initParticleLog_MEM

  !> Routine to initialize DPS particle log file (e.g. print header)
  subroutine initParticleLog_DPS( particleID, fileName, plogUnit, fileExists )
    integer, intent(in) :: particleID
    character(*), intent(in) :: fileName
    integer, intent(in) :: plogUnit
    logical, intent(in) :: fileExists
    ! --------------------------------------------!

    ! write(logUnit(1),*) 'Calling initParticleLog!'
    if(fileExists) then
      open(plogUnit, file = filename, status = 'old', position='append')
    else
      ! write(logUnit(1),*) 'INITPARTICLELOG logUnit', plogUnit
      open(plogUnit, file = filename, status = 'new')
      write(plogUnit, '(A,I5)') 'Particle ID = ', particleID
      write(plogUnit, '(A18)', advance='no') 't'
      write(plogUnit, '(A18)', advance='no') 'x'
      write(plogUnit, '(A18)', advance='no') 'y'
      write(plogUnit, '(A18)', advance='no') 'z'
      write(plogUnit, '(A18)', advance='no') 'ux'
      write(plogUnit, '(A18)', advance='no') 'uy'
      write(plogUnit, '(A18)', advance='no') 'uz'
      write(plogUnit, '(A18)', advance='no') 'urx'
      write(plogUnit, '(A18)', advance='no') 'ury'
      write(plogUnit, '(A18)', advance='no') 'urz'
      write(plogUnit, '(A18)', advance='no') 'Fx'
      write(plogUnit, '(A18)', advance='no') 'Fy'
      write(plogUnit, '(A18)', advance='no') 'Fz'
      write(plogUnit, '(A18)', advance='no') 'RFx'
      write(plogUnit, '(A18)', advance='no') 'RFy'
      write(plogUnit, '(A18)', advance='no') 'RFz'
      write(plogUnit, '(A18)', advance='no') 'Fxa'
      write(plogUnit, '(A18)', advance='no') 'Fya'
      write(plogUnit, '(A18)') 'Fza'
    end if

  end subroutine initParticleLog_DPS


  subroutine closeParticleLog(plogUnit)
      integer, intent(in) :: plogUnit
      ! --------------------------------------------!
      close(plogUnit)
  end subroutine closeParticleLog

  !> openLogFile opens a file with name fileName and returns the unit
  !! attached to it. It checks whether the file exists and if not
  subroutine openLogFile( fileName, plogUnit, isNewFile )
      character(len=*), intent(inout) :: filename
      integer, intent(out) :: plogUnit
      logical, intent(out), optional :: isNewFile
      ! --------------------------------------------!
      logical :: fileExists, fileIsOpen
      ! --------------------------------------------!
      ! Check if file with fileName already exists and whether it is already open
      inquire(file=trim(filename), exist=fileExists, opened=fileIsOpen, number=plogUnit)

      if(fileIsOpen) then
        ! File is open and we can write to it on unit logUnit
        if(present(isNewFile)) isNewFile = .FALSE.
        return
      else ! file is not open
        plogUnit = newunit()
        if( fileExists ) then
          open(plogUnit, file = trim(fileName), status = 'old', position='append')
          if(present(isNewFile)) isNewFile = .FALSE.
        else
          open(plogUnit, file = trim(fileName), status = 'new')
          if(present(isNewFile)) isNewFile = .TRUE.
        end if
      end if
  end subroutine openLogFile

  !> Debugging routine to create a list of elements along a line so
  !! that the properties of these elements can be printed
  subroutine generateElemListLine( dir, xstart, length, scheme, geometry, elemList )
    ! Cartesian direction vector, i.e. dir = (/ 1,0,0 /) for line in x-direction
    integer, intent(in) :: dir(3)
    ! Cartesian starting coordinate
    real(kind=rk), intent(in) :: xstart(3)
    ! Line will run from xstart to xstart + dir*length
    real(kind=rk), intent(in) :: length
    !> scheme
    type(mus_scheme_type), intent(in) :: scheme
    !> geometry
    type(mus_geom_type), intent(in) :: geometry
    !> output: element list
    type(dyn_intArray_type), intent(inout) :: elemList
    ! -------------------------------------------- !
    integer :: coord(4), lev
    integer(kind=long_k) :: TreeID
    integer :: ldPos
    real(kind=rk) :: x(3), s, dx
    ! -------------------------------------------- !
    lev = geometry%tree%global%maxLevel
    dx = geometry%tree%global%BoundingCubeLength / 2**lev

    ! Empty the array
    call empty(me = elemList)

    x = xstart
    s = 0.0_rk

    do while( s < length )
      ! Get coordinate
      coord = tem_CoordOfReal( mesh  = geometry%tree,     &
                             & point = x,                 &
                             & level = lev                )

      ! Get TreeID and position in total list
      TreeID = tem_IdOfCoord(coord = coord)
      ldPos = tem_PosOfId( sTreeID    = TreeID,                      &
                         & treeIDlist = scheme%levelDesc(lev)%total, &
                         & lower      = 1,                           &
                         & upper      = scheme%pdf(lev)%nElems_fluid )

      ! If we found the element on this proc, add it to elemList
      if(ldPos > 0) then
        call append( me       = elemList, &
                   & val      = ldPos     )
      end if

      x = x + dir*dx
      s = s + dx

    end do

  end subroutine generateElemListLine

  !> Debugging routine to dump debug tracking data
  subroutine dumpdata(tracker, t, scheme, geometry, params)
    !> Tracker containing file name and elemList
    type(mus_particle_debugtracking_type), intent(inout) :: tracker
    !> Current simulation time
    real(kind=rk), intent(in) :: t
    !> Scheme for access to fluid data
    type(mus_scheme_type), intent(inout) :: scheme
    !> Geometry for access to tree
    type(mus_geom_type), intent(in) :: geometry
    !> Params
    type(mus_param_type), intent(in) :: params
    ! -------------------------------------------- !
    integer :: plogUnit
    integer :: lev
    integer :: vel_pos(3), dens_pos, vol_frac_pos
    integer :: elemOff
    integer :: iElem, i
    integer :: coord(4)
    integer :: ldPos
    integer(kind=long_k) :: TIDoffset
    real(kind=rk) :: x(3), dx
    real(kind=rk) :: u_fluid(3)
    real(kind=rk) :: eps_f
    character(len=1024) :: fileName
    ! -------------------------------------------- !
    ! If tracker is inactive, do nothing and return from the routine
    if( .NOT. tracker%active ) then
      return
    end if

    ! If tracker is active, log data
    lev = geometry%tree%global%maxlevel
    dx  = params%physics%dxLvl(lev)
    TIDoffset = tem_firstIdAtLevel(lev)

    write(fileName,'(A,A,E11.5,A)') trim(tracker%lfile), '_t', t, '.dat'


    ! Position of variables in auxField
    vel_pos      = scheme%varSys%method%val(scheme%derVarPos(1)%velocity)%auxField_varPos(1:3)
    dens_pos     = scheme%varSys%method%val(scheme%derVarPos(1)%density)%auxField_varPos(1)
    vol_frac_pos = scheme%varSys%method%val(scheme%derVarPos(1)%vol_frac)%auxField_varPos(1)

    ! Open log file
    call openLogFile( fileName, plogUnit )

    do iElem = 1, tracker%elemList%nvals
      ldPos = int(tracker%elemList%val(iElem))
      elemOff = (ldPos-1)*scheme%varSys%nAuxScalars

      ! Get cartesian coordinates of current element
      coord = tem_coordOfID( TreeID = scheme%levelDesc(lev)%total(ldPos), &
                           & offset = TIDoffset)

      x = getBaryOfCoord( coord  = coord,                       &
                        & origin = geometry%tree%global%origin, &
                        & dx     = dx                           )

      ! Get auxField data to write to logUnit
      u_fluid(1) = scheme%auxField(lev)%val(elemOff+vel_pos(1))
      u_fluid(2) = scheme%auxField(lev)%val(elemOff+vel_pos(2))
      u_fluid(3) = scheme%auxField(lev)%val(elemOff+vel_pos(3))
      eps_f      = scheme%auxField(lev)%val(elemOff+vol_frac_pos)

      u_fluid = u_fluid*params%physics%fac(lev)%vel / eps_f

      ! Write data to logUnit

      ! Time
      write(plogUnit, '(E17.9)', advance='no') t
      write(plogUnit, '(A)', advance='no') ' '

      ! Element position
      do i = 1, 3
          write(plogUnit, '(E17.9)', advance='no') x(i)
          write(plogUnit, '(A)', advance='no') ' '
      end do

      ! Fluid velocity
      do i = 1, 3
          write(plogUnit, '(E17.9)', advance='no') u_fluid(i)
          write(plogUnit, '(A)', advance='no') ' '
      end do

      ! Newline
      write(plogUnit, '(A)')
    end do

    close(plogUnit)

  end subroutine dumpdata

  !> perform_particle_checks computes the total momentum of the fluid and particle
  !! phases and logs this to a file
  subroutine mus_log_fluid_momentum(scheme, lev, params, comm)
    !> scheme for access to auxfield
    type(mus_scheme_type), intent(inout) :: scheme
    !> level in the octree to compute total momentum on
    integer, intent(in) :: lev
    !> params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> MPI communicator
    integer :: comm
    ! ------------------------------------------------------- !
    integer :: plogUnit, k, myRank, nProcs
    character(len=1024) :: fileName
    real(kind=rk) :: fluidMomentum(3)
    real(kind=rk) :: t
    logical :: isNewFile

    ! ------------------------------------------------------- !
    myRank = params%general%proc%rank
    nProcs = params%general%proc%comm_size

    t = params%general%simcontrol%now%sim
    write(fileName,'(A)') 'momentumcheck.txt'

    call compute_fluid_momentum( scheme        = scheme,         &
                                & lev           = lev,           &
                                & params        = params,        &
                                & totalMomentum = fluidMomentum, &
                                & comm          = comm,          &
                                & myRank        = myRank,        &
                                & nProcs        = nProcs         )

    ! Only root process opens and writes to the log file
    if(myRank == 0) then
      call openLogFile( fileName, plogUnit, isNewFile )
      if(isNewFile) then
      ! Write header if the file is new
        write(plogUnit, '(A18)', advance='no') 't'
        write(plogUnit, '(A18)', advance='no') 'fluidMomX'
        write(plogUnit, '(A18)', advance='no') 'fluidMomY'
        write(plogUnit, '(A18)', advance='no') 'fluidMomZ'
        write(plogUnit, '(A)')
      end if ! isNewFile

      ! Write time and particle and fluid momentum to log file
      write(plogUnit, '(E17.9)', advance='no') t

      do k = 1,3
        write(plogUnit, '(E17.9)', advance='no') fluidMomentum(k)
      end do

      ! Write a newline
      write(plogUnit, '(A)', advance='no')

      close(plogUnit)
    end if ! myRank == 0

  end subroutine mus_log_fluid_momentum

  !> perform_particle_checks computes the total momentum of the fluid and particle
  !! phases and logs this to a file
  subroutine mus_particles_log_total_momentum(particleGroup, scheme, lev, params, comm, t)
    !> Particle group on this process
    type(mus_particle_group_type), intent(in) :: particleGroup
    !> scheme for access to auxfield
    type(mus_scheme_type), intent(inout) :: scheme
    !> level in the octree to compute total momentum on
    integer, intent(in) :: lev
    !> params for access to dt, dx, etc.
    type(mus_param_type), intent(in) :: params
    !> MPI communicator
    integer :: comm
    !> Current (physical) time
    real(kind=rk) :: t
    ! ------------------------------------------------------- !
    integer :: plogUnit, k, myRank, nProcs
    character(len=1024) :: fileName
    real(kind=rk) :: fluidMomentum(3)
    real(kind=rk) :: particleMomentum(3)
    logical :: isNewFile

    ! ------------------------------------------------------- !
    myRank = params%general%proc%rank
    nProcs = params%general%proc%comm_size

    write(fileName,'(A)') 'momentumcheck.txt'

    call compute_fluid_momentum( scheme        = scheme,         &
                                & lev           = lev,           &
                                & params        = params,        &
                                & totalMomentum = fluidMomentum, &
                                & comm          = comm,          &
                                & myRank        = myRank,        &
                                & nProcs        = nProcs         )

    call compute_particle_momentum( particleGroup = particleGroup,    &
                                  & lev           = lev,              &
                                  & params        = params,           &
                                  & totalMomentum = particleMomentum, &
                                  & comm          = comm,             &
                                  & myRank        = myRank,           &
                                  & nProcs        = nProcs            )

    ! Only root process opens and writes to the log file
    if(myRank == 0) then
      call openLogFile( fileName, plogUnit, isNewFile )
      if(isNewFile) then
      ! Write header if the file is new
        write(plogUnit, '(A18)', advance='no') 't'
        write(plogUnit, '(A18)', advance='no') 'fluidMomX'
        write(plogUnit, '(A18)', advance='no') 'fluidMomY'
        write(plogUnit, '(A18)', advance='no') 'fluidMomZ'
        write(plogUnit, '(A18)', advance='no') 'particleMomX'
        write(plogUnit, '(A18)', advance='no') 'particleMomY'
        write(plogUnit, '(A18)', advance='no') 'particleMomZ'
        write(plogUnit, '(A18)', advance='no') 'totalMomX'
        write(plogUnit, '(A18)', advance='no') 'totalMomY'
        write(plogUnit, '(A18)', advance='no') 'totalMomZ'
        write(plogUnit, '(A)')
      end if ! isNewFile

      ! Write time and particle and fluid momentum to log file
      write(plogUnit, '(E17.9)', advance='no') t

      do k = 1,3
        write(plogUnit, '(E17.9)', advance='no') fluidMomentum(k)
      end do

      do k = 1,3
        write(plogUnit, '(E17.9)', advance='no') particleMomentum(k)
      end do

      do k = 1,3
        write(plogUnit, '(E17.9)', advance='no') fluidMomentum(k) + particleMomentum(k)
      end do

      ! Write a newline
      write(plogUnit, '(A)', advance='no')

      close(plogUnit)
    end if ! myRank == 0

  end subroutine mus_particles_log_total_momentum

  !> Routine to dump particle timing data
  subroutine dump_particle_timing(proc_logUnit)
    !> Unit to log local timing data to
    integer, intent(in) :: proc_logUnit
    ! ------------------------------------- !
    real(kind=rk) :: t_upPos, t_upVel, t_Fhydro, t_FDEM
    real(kind=rk) :: t_idle, t_MomInc, t_coll, t_particle, t_intp, t_aux
    real(kind=rk) :: t_pos, t_vel, t_force, t_sub, t_loadParticles
    ! real(kind=rk) :: tParticles
    ! character(len=100) :: fileName
    ! integer :: plogUnit, iError
    ! integer :: Nprocs, myRank
    ! ------------------------------------- !

    ! Get local timer values
    t_upPos = tem_getTimerVal( timerHandle = mus_particle_timerHandles%exchangePositionsTimer )
    t_upVel = tem_getTimerVal( timerHandle = mus_particle_timerHandles%exchangeVelocitiesTimer )
    t_Fhydro = tem_getTimerVal( timerHandle = mus_particle_timerHandles%exchangeHydroForcesTimer )
    t_FDEM = tem_getTimerVal( timerHandle = mus_particle_timerHandles%exchangeDEMForcesTimer )
    t_idle = tem_getTimerVal( timerHandle = mus_particle_timerHandles%idleTimer )
    t_MomInc = tem_getTimerVal( timerHandle = mus_particle_timerHandles%exchangeMomIncTimer )
    t_coll = tem_getTimerVal( timerHandle = mus_particle_timerHandles%collisionHandlingTimer )
    t_intp = tem_getTimerVal( timerHandle = mus_particle_timerHandles%interpolateFluidTimer )
    t_aux = tem_getTimerVal( timerHandle = mus_particle_timerHandles%incrementAuxFieldTimer )
    t_particle = tem_getTimerVal( timerHandle = mus_particle_timerHandles%mainParticleTimer )

    t_pos = tem_getTimerVal( timerHandle = mus_particle_timerHandles%positionTimer )
    t_vel = tem_getTimerVal( timerHandle = mus_particle_timerHandles%velocityTimer )
    t_force = tem_getTimerVal( timerHandle = mus_particle_timerHandles%forceTimer )
    t_sub = tem_getTimerVal( timerHandle = mus_particle_timerHandles%subcycleTimer )
    t_loadParticles = tem_getTimerVal( timerHandle = mus_particle_timerHandles%loadParticleTimer )


    ! Write to plogUnit (should be opened before calling this routine)
    write(proc_logUnit,'(A)') "------- PARTICLE CODE TIMING -------"
    write(proc_logUnit,'(A25,E17.9)') "loadParticles = ", t_loadParticles
    write(proc_logUnit,'(A25,E17.9)') "exchangePositions = ", t_upPos
    write(proc_logUnit,'(A25,E17.9)') "exchangeVelocities = ", t_upVel
    write(proc_logUnit,'(A25,E17.9)') "exchangeHydroForces = ", t_Fhydro
    write(proc_logUnit,'(A25,E17.9)') "exchangeDEMForces = ", t_FDEM
    write(proc_logUnit,'(A25,E17.9)') "idle time comm = ", t_idle
    write(proc_logUnit,'(A25,E17.9)') "exchangeMomInc = ", t_MomInc
    write(proc_logUnit,'(A25,E17.9)') "collision detection = ", t_coll
    write(proc_logUnit,'(A25,E17.9)') "Interpolate fluid props = ", t_intp
    write(proc_logUnit,'(A25,E17.9)') "Increment auxField = ", t_aux
    write(proc_logUnit,'(A)') "------------------------------"
    write(proc_logUnit,'(A25,E17.9)') "Total positions = ", t_pos
    write(proc_logUnit,'(A25,E17.9)') "Total velocities = ", t_vel
    write(proc_logUnit,'(A25,E17.9)') "Total forces = ", t_force
    write(proc_logUnit,'(A25,E17.9)') "Total DEM subcycles = ", t_sub
    write(proc_logUnit,'(A)') "------------------------------"
    write(proc_logUnit,'(A25,E17.9)') "Total particle time = ", t_particle

    ! Calculate the maximum value of the timer.
    ! tParticles = tem_getMaxTimerVal( timerHandle = timerHandle, &
    !                                & comm        = comm         )

    ! call MPI_COMM_SIZE( comm, Nprocs, iError )
    ! call MPI_COMM_RANK( comm, myRank, iError )

    ! write(fileName,'(A)') 'mus_particle_timing.res'
    ! ! Open file to output timing to
    ! if(myRank == 0) then
    !   call openLogFile( fileName, plogUnit)

    !   ! Write timing data to file
    !   write(plogUnit,'(A18)', advance='no') 'tParticles'
    !   write(plogUnit,'(A18)', advance='no') 'Nprocs'
    !   write(plogUnit,'(A)')
    !   write(plogUnit,'(E17.9)', advance='no') tParticles
    !   write(plogUnit,'(I18)', advance='no') Nprocs

    !   close(plogUnit)
    ! end if

  end subroutine dump_particle_timing

end module mus_particle_logging_module
