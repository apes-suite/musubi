! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2014,2020-2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2017, 2019, 2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Nikhil Anand <n.anand@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!See copyright notice in the COPYRIGHT file.
!! *************************************************************************** !
!> Some generic handy check routines to check the properties of the flow field
!! and the current run
!!
!! Performance estimation, check of the total density, check for NaNs
!!
module mus_tools_module

  ! include treelm modules
  use mpi
!$ use omp_lib
  use env_module,             only: rk, long_k, newunit, PathLen, labelLen, &
    &                               rk_mpi, outLen
  use tem_time_module,        only: tem_time_dump
  use tem_timeControl_module, only: tem_timeControl_globalTriggered
  use tem_timer_module,       only: tem_getTimerVal, tem_getMaxTimerVal, &
    &                               tem_getTimerName
  use tem_status_module,      only: tem_stat_nan_detected, tem_stat_nonPhysical
  use tem_isNaN_module,       only: tem_isNaN
  use tem_debug_module,       only: dbgUnit
  use tem_logging_module,     only: logUnit
  use tem_topology_module,    only: tem_LevelOf
  use tem_comm_env_module,    only: tem_comm_env_type
  use tem_comm_module,        only: tem_communication_type
  use tem_general_module,     only: tem_general_type
  use tem_balance_module,     only: tem_calc_imbalance
  use tem_aux_module,         only: check_mpi_error, tem_abort
  use tem_tools_module,       only: tem_horizontalSpacer
  use tem_param_module,       only: cs

  ! include musubi modules
  use mus_abortCriteria_module, only: mus_abortCriteria_type
  use mus_param_module,       only: mus_param_type, mus_param_out
  use mus_scheme_type_module, only: mus_scheme_type
  use mus_IBM_module,         only: mus_finishIBM
  use mus_scheme_module,      only: mus_scheme_out
  use mus_physics_module,     only: mus_physics_out
  use mus_timer_module,       only: get_computeRatio, get_communicateRatio, &
    &                               get_intpRatio, get_boundaryRatio,       &
    &                               get_bcBufferTime, get_bcBufferRatio,    &
    &                               get_computeTime, get_communicateTime,   &
    &                               get_boundaryTime,                       &
    &                               get_auxTime,                            &
    &                               get_relaxTime,                          &
    &                               get_mainLoopTime, mus_timerHandles
  use mus_relaxationParam_module, only: mus_check_omegaKine

  use aotus_module, only: flu_state, aot_get_val
  use aot_out_module, only: aot_out_type, aot_out_open, aot_out_close, &
    &                       aot_out_val

  implicit none

  private

  public :: perform_checks
  public :: check_streaming_layout
  public :: check_density
  public :: check_potential
  public :: mus_writeSolverSpecInfo
  public :: dump_linear_partition
  public :: mus_perf_measure
  public :: mus_BC_timing
  public :: dump_bc_timing


contains


  ! ------------------------------------------------------------------------ !
  !> Perform run-time checks if interval is active
  !!
  subroutine perform_checks( scheme, minLevel, maxLevel, general, &
    &                        mus_aborts, initCheck                )
    ! -------------------------------------------------------------------- !
    type(mus_scheme_type), intent(in) :: scheme
    integer, intent(in) :: minLevel, maxLevel
    type(tem_general_type), intent(inout) :: general
    type(mus_abortCriteria_type), intent(in) :: mus_aborts
    !> True for initial check before main time toop
    logical, intent(in) :: initCheck
    ! -------------------------------------------------------------------- !
    !  Check and output in the main intervals
    if ( tem_timeControl_globalTriggered(       &
      &    me = general%simControl%timeControl, &
      &    now = general%simControl%now,        &
      &    comm = general%proc%comm )           &
      & .or. initCheck                          ) then

      select case (trim(scheme%header%kind))
      case ('fluid', 'fluid_incompressible',       &
        &   'fluid_GNS', 'fluid_incompressible_GNS')
        ! check total density
        call check_density( scheme, minLevel, maxLevel, general )
        ! check whether lattice velocity above stability limit
        call check_velocityFluid( scheme, minLevel, maxLevel, &
          &                       general, mus_aborts         )
        ! Check omega range only if viscKine STfun is not constant or
        ! it is initial check.
        associate(fluid => scheme%field(1)%fieldProp%fluid)
          if ( (trim(fluid%viscKine%STfun%fun_kind) /= 'const') &
            & .or. initCheck                                    ) then
            call mus_check_omegaKine(                            &
              &    omLvlKine   = fluid%viscKine%omLvl,           &
              &    nSolve      = scheme%pdf(:)%nElems_solve,     &
              &    schemeRelax = trim(scheme%header%relaxation), &
              &    minLevel    = minLevel,                       &
              &    maxLevel    = maxLevel,                       &
              &    general     = general                         )
          end if
        end associate

      case ('multispecies_gas', 'multispecies_liquid')
        call check_density( scheme, minLevel, maxLevel, general )
        call check_velocityMS( scheme, minLevel, maxLevel, general, mus_aborts )
      case ('passive_scalar', 'nernst_planck', 'isotherm_acEq')
        call check_density( scheme, minLevel, maxLevel, general )
      case ('poisson', 'poisson_boltzmann_linear', &
        &  'poisson_boltzmann_nonlinear')
        call check_potential( scheme, minlevel, maxLevel, general )
      case default
        write(logUnit(1),*) 'Unknown scheme kind: '//trim(scheme%header%kind)
        call tem_abort()
      end select
    end if

  end subroutine perform_checks
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Check the total potential for poisson scheme
  !!
  !! The output might be delayed by using arrays which are then dumped
  !! to keep disc access more restricted
  !!
  subroutine check_potential( scheme, minLevel, maxLevel, general, &
    &                         total_potential )
    ! -------------------------------------------------------------------- !
    !> scheme type
    type(mus_scheme_type), intent(in) :: scheme
    !> global scheme independent information
    integer, intent(in) :: minLevel, maxLevel
    type(tem_general_type), intent(inout) :: general
    real(kind=rk), intent(out), optional :: total_potential
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: pot, total_pot
    integer iLevel, ierror, iElem
    ! -------------------------------------------------------------------- !
    total_pot = 0._rk
    pot = 0._rk
    ! Sum up all values in the complete vector to get the total density
    lvlLoop: do iLevel = minLevel, maxLevel
      do iElem = 1, scheme%pdf( iLevel )%nElems_fluid
        ! Potential is stored in auxField. Use it to check total potential
        pot = pot + scheme%auxField(iLevel)%val(iElem)
      end do ! iElem
    end do lvlLoop !iLevel

    ! All ranks check, if there are NaN entries
    if ( tem_isnan(pot) ) then
      write(logUnit(1),*) ' Error: Total potential is NaN on rank ', &
        &                 general%proc%rank
      general%simControl%status%bits(tem_stat_nan_detected) = .true.
    endif

    ! compute total_pot
    call mpi_reduce( pot, total_pot, 1, mpi_double_precision, &
                     mpi_sum, 0, general%proc%comm, ierror    )

    call tem_time_dump(me = general%simControl%now, outUnit = logUnit(1))

    ! dump total potential for each field
    write(logUnit(1),'(a11,a,a20)')        'field',' |', 'total potential'
    write(logUnit(1), '(i11,a,en20.10)' ) scheme%nFields,' |', total_pot
    write(logUnit(1), '(a)' ) ''
    if (present(total_potential)) then
      total_potential = total_pot
    end if

  end subroutine check_potential
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Check the total density for a selected scheme and write to unit
  !!
  !! The output might be delayed by using arrays which are then dumped
  !! to keep disc access more restricted
  !!
  subroutine check_density( scheme, minLevel, maxLevel, general, total_density )
    ! -------------------------------------------------------------------- !
    !> scheme type
    type(mus_scheme_type), intent(in) :: scheme
    !> global scheme independent information
    integer, intent(in) :: minLevel, maxLevel
    type(tem_general_type), intent(inout) :: general
    real(kind=rk), intent(out), optional :: total_density
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: dens(:) , total_dens(:)
    real(kind=rk) :: dens_field
    integer :: iLevel, ierror, iElem, iField
    integer :: elemOff, dens_pos(scheme%nFields)
    ! -------------------------------------------------------------------- !
    ! density is stored in auxField. Use it to check density
    select case (trim(scheme%header%kind))
    case ('nernst_planck')
      do iField = 1, scheme%nFields
        dens_pos(iField) = scheme%varSys                                    &
          &                      %method                                    &
          &                      %val(scheme%derVarpos(iField)%moleDensity) &
          &                      %auxField_varPos(1)
      end do
    case default
      do iField = 1, scheme%nFields
        dens_pos(iField) = scheme%varSys%method                        &
          &                     %val(scheme%derVarpos(iField)%density) &
          &                     %auxField_varPos(1)
      end do
    end select

    allocate( dens( scheme%nFields ) )
    allocate( total_dens( scheme%nFields ) )
    total_dens = 0._rk
    dens = 0._rk
    ! Sum up all values in the complete vector to get the total density
    lvlLoop: do iLevel = minLevel, maxLevel
      do iElem = 1, scheme%pdf( iLevel )%nElems_fluid
        ! element offset for auxField
        elemOff = (iElem-1)*scheme%varSys%nAuxScalars
        ! compute density for each field
        do iField = 1, scheme%nFields
          ! field density
          dens_field = scheme%auxField(iLevel)%val(elemOff + dens_pos(iField))
          ! set nonPhysical if local field density is negative
          ! density can never be negative
          if (dens_field <= 0.0_rk) then
            write(logUnit(1),'(a,i3)') 'ERROR: Negative density for field', &
              &                        iField
            write(logUnit(1),*) 'Terminating check total density calculation'
            write(logUnit(1),*) 'non-physical density=',dens_field

            write(dbgUnit(1),'(a,i3)') 'ERROR: Negative density for field', &
              &                        iField
            write(dbgUnit(1),*) 'Terminating check total density calculation'
            write(dbgUnit(1),*) 'non-physical density=',dens_field
            write(dbgUnit(1),*) 'TreeID ', scheme%levelDesc(iLevel)%total(iElem)
            general%simControl%status%bits(tem_stat_nonPhysical) = .true.
            exit lvlLoop
          end if
          ! total density
          dens( iField ) = dens( iField ) + dens_field
        end do ! iField
      end do ! iElem
    end do lvlLoop ! iLevel

    ! All ranks check, if there are NaN entries
    if ( any( tem_isnan(dens) )  ) then
      write(logUnit(1),*) ' Error: Total density is NaN on rank ', &
        &                 general%proc%rank
      general%simControl%status%bits(tem_stat_nan_detected) = .true.
    end if

    ! compute total_dens
    call mpi_reduce( dens, total_dens,scheme%nFields, mpi_double_precision, &
                     mpi_sum, 0, general%proc%comm, ierror                  )

    call tem_time_dump(me = general%simControl%now, outUnit = logUnit(1))

    ! dump total density for each field
    write(logUnit(1),'(a11,a,a20)')        'field',' |', 'total density'
    do iField = 1, scheme%nFields
      write(logUnit(1), '(i11,a,en20.10)' ) iField,' |', total_dens(iField)
    end do
    write(logUnit(1), '(a)' ) ''
    if (present(total_density)) then
      total_density = sum(total_dens)
    end if

  end subroutine check_density
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Check the maximum velocity whether it is above Ma>0.1
  !!
  subroutine check_velocityFluid( scheme, minLevel, maxLevel, general, &
    &                             mus_aborts )
    ! -------------------------------------------------------------------- !
    !> scheme type
    type(mus_scheme_type), intent(in) :: scheme
    !> global scheme independent information
    integer, intent(in) :: minLevel, maxLevel
    type(tem_general_type), intent(inout) :: general
    type(mus_abortCriteria_type), intent(in) :: mus_aborts
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: vel(3) , velMag, maxVel, glob_maxVel
    integer :: iLevel, ierror, iElem
    integer :: elemOff, vel_pos(3)
    ! -------------------------------------------------------------------- !
    ! Sum up all values in the complete vector to get the total density
    maxVel = 0.0_rk
    lvlLoop: do iLevel = minLevel, maxLevel
      do iElem = 1, scheme%pdf(iLevel)%nElems_fluid
        ! element offset for auxField
        elemOff = (iElem-1)*scheme%varSys%nAuxScalars
        ! velocity is stored in auxField. Use it to check velocity
        vel_pos = scheme%varSys%method                            &
          &                    %val(scheme%derVarpos(1)%velocity) &
          &                    %auxField_varPos(1:3)
        vel(1) = scheme%auxField(iLevel)%val(elemOff + vel_pos(1))
        vel(2) = scheme%auxField(iLevel)%val(elemOff + vel_pos(2))
        vel(3) = scheme%auxField(iLevel)%val(elemOff + vel_pos(3))
        velMag = sqrt(dot_product(vel, vel))
        maxVel = max(maxVel, velMag)
     end do
   end do lvlLoop

   ! All ranks check, if there are NaN entries
   if (tem_isnan(maxVel)) then
     write(logUnit(1),*) ' Error: Lattice velocity is NaN on rank ', &
       &                 general%proc%rank
     general%simControl%status%bits(tem_stat_nan_detected) = .true.
   endif

   ! compute maximum of all processes
   call mpi_reduce( maxVel, glob_maxVel, 1, mpi_double_precision, &
                    mpi_max, 0, general%proc%comm, ierror )

   if (general%proc%isRoot) then
     if (glob_maxVel > mus_aborts%velLat_max) then
       write(logUnit(1),'(2(A,F6.3))')                      &
         &   'Error: Maximum lattice velocity in domain: ', &
         &   glob_maxVel , ' > ', mus_aborts%velLat_max
       general%simControl%status%bits(tem_stat_nonPhysical) = .true.
     else if (glob_maxVel > 0.1_rk) then
       write(logUnit(1),'(A,F6.3)') 'WARNING: Maximum lattice velocity in ' &
         &                          //'domain:', glob_maxVel
     end if
   end if

  end subroutine check_velocityFluid
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Check the maximum velocity whether it is above Ma>0.1
  !!
  subroutine check_velocityMS( scheme, minLevel, maxLevel, general, &
    &                          mus_aborts )
    ! -------------------------------------------------------------------- !
    !> scheme type
    type(mus_scheme_type), intent(in) :: scheme
    !> global scheme independent information
    integer, intent(in) :: minLevel, maxLevel
    type(tem_general_type), intent(inout) :: general
    type(mus_abortCriteria_type), intent(in) :: mus_aborts
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: velMag, maxVel, glob_maxVel
    real(kind=rk) :: dens(scheme%nFields), mom(3, scheme%nFields)
    real(kind=rk) :: mixVel(3), inv_tRho
    integer :: iLevel, ierror, iElem, iField
    integer :: elemOff, dens_pos, mom_pos(3)
    ! -------------------------------------------------------------------- !
    ! Sum up all values in the complete vector to get the total density
    maxVel = 0.0_rk
    lvlLoop: do iLevel = minLevel, maxLevel
      do iElem = 1, scheme%pdf(iLevel)%nElems_fluid
        ! element offset for auxField
        elemOff = (iElem-1)*scheme%varSys%nAuxScalars
        do iField = 1, scheme%nFields
          ! velocity is stored in auxField. Use it to check velocity
          dens_pos = scheme%varSys%method                                &
            &                     %val(scheme%derVarpos(iField)%density) &
            &                     %auxField_varPos(1)
          mom_pos = scheme%varSys%method                                 &
            &                    %val(scheme%derVarpos(iField)%momentum) &
            &                    %auxField_varPos(1:3)
          dens(iField) = scheme%auxField(iLevel)%val(elemOff + dens_pos)
          mom(1, iField) = scheme%auxField(iLevel)%val(elemOff + mom_pos(1))
          mom(2, iField) = scheme%auxField(iLevel)%val(elemOff + mom_pos(2))
          mom(3, iField) = scheme%auxField(iLevel)%val(elemOff + mom_pos(3))
        end do
        inv_tRho = 1.0_rk / sum(dens)
        ! mixture averaged mass velocity
        mixVel(1) = sum(mom(1,:)) * inv_tRho
        mixVel(2) = sum(mom(2,:)) * inv_tRho
        mixVel(3) = sum(mom(3,:)) * inv_tRho
        velMag = sqrt(dot_product(mixVel, mixVel))
        maxVel = max(maxVel, velMag)
     end do
   end do lvlLoop

   ! All ranks check, if there are NaN entries
   if (tem_isnan(maxVel)) then
     write(logUnit(1),*) ' Error: Lattice velocity is NaN on rank ', &
       &                 general%proc%rank
     general%simControl%status%bits(tem_stat_nan_detected) = .true.
   endif

   ! compute maximum of all processes
   call mpi_reduce( maxVel, glob_maxVel, 1, mpi_double_precision, &
                    mpi_max, 0, general%proc%comm, ierror )

   if (general%proc%isRoot) then
     if (glob_maxVel > mus_aborts%velLat_max) then
       write(logUnit(1),'(2(A,F6.3))')                      &
         &   'Error: Maximum lattice velocity in domain: ', &
         &   glob_maxVel , ' > ', mus_aborts%velLat_max
       general%simControl%status%bits(tem_stat_nonPhysical) = .true.
     else if (glob_maxVel > 0.1_rk) then
       write(logUnit(1),'(A,F6.3)') 'WARNING: Maximum lattice velocity in ' &
         &                          //'domain:', glob_maxVel
     end if
   end if

  end subroutine check_velocityMS
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine measures performance imbalance, MLUPS and dumps timings
  !! to disk
  subroutine mus_perf_measure( totalDens, DomSize, minLevel,          &
    &                          maxLevel, nElems, scaleFactor, general )
    ! -------------------------------------------------------------------- !
    !> Total density from check_density
    real(kind=rk), intent(in) :: totalDens
    !> Total number of elements in tree
    integer(kind=long_k), intent(in) :: DomSize
    !> level range
    integer, intent(in) :: minLevel, maxLevel
    !> array of nElems levelwise
    integer, intent(in) :: nElems(minLevel:maxLevel)
    !> global parameter
    integer, intent(in) :: scaleFactor
    !> Contains proc, simControl, solveHead
    type(tem_general_type), intent(in) :: general
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: tMainLoop, tCompute, mlups, mlups_kernel
    real(kind=rk) :: cpuCost, imbalance
    integer       :: iter, iTimer, nTimers, counter, iErr, nLevels
    real(kind=rk), allocatable :: timerVal(:)
    integer(kind=long_k) :: nTotals(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    cpuCost = get_mainLoopTime() - get_communicateTime()
    imbalance = tem_calc_imbalance( cpuCost,                &
      &                             general%proc%comm,      &
      &                             general%proc%comm_size, &
      &                             general%proc%isRoot     )

    ! first and last handle convers all musubi handles which are added
    ! contigously
    nTimers = mus_timerHandles%last - mus_timerHandles%first + 1
    allocate( timerVal( nTimers ) )
    ! Get maxTimer of all process
    counter = 0
    do iTimer = mus_timerHandles%first, mus_timerHandles%last
      counter = counter+1
      timerVal(counter) = tem_getMaxTimerVal( timerHandle = iTimer,           &
        &                                     comm        = general%proc%comm )
    end do

    nLevels = maxLevel - minLevel + 1
    call mpi_allreduce( int(nElems(minLevel:maxLevel), long_k), &
      &                 nTotals(minLevel:maxLevel),             &
      &                 nLevels, mpi_integer8, mpi_sum,         &
      &                 general%proc%comm, iErr                 )

    if ( general%proc%isRoot ) then

      iter = ( general%simControl%now%iter               &
        &    - general%simControl%timeControl%min%iter ) &
        &    / ( scaleFactor ** ( maxLevel - minLevel )  )

      ! calculate MLUPS and MLUPS kernel
      tMainLoop = get_mainLoopTime()
      mlups = calc_MLUPS( nElems      = nTotals,     &
        &                 minLevel    = minLevel,    &
        &                 maxLevel    = maxLevel,    &
        &                 scaleFactor = scaleFactor, &
        &                 iter        = iter,        &
        &                 time        = tMainLoop    )

!write(logUnit(7),"(A)") "Calculate MLUPS f r this stage."
!write(logUnit(7),"(A,I0)") " Now iter: ", general%simControl%now%iter
!write(logUnit(7),"(A,I0)") " Ran iter: ", iter
!write(logUnit(7),"(A,F10.2)") " Main Loop time: ", tMainLoop
!write(logUnit(7),"(A,F10.2)") " MLUPS: ", mlups

      ! Using compute kernel timer
      tCompute = get_computeTime()
      mlups_kernel = calc_MLUPS( nElems       = nTotals,     &
        &                        minLevel     = minLevel,    &
        &                        maxLevel     = maxLevel,    &
        &                        scaleFactor  = scaleFactor, &
        &                        iter         = iter,        &
        &                        time         = tCompute     )

      ! dump timing
      call dump_timing( totalDens    = totalDens,    &
        &               DomSize      = DomSize,      &
        &               MLUPS        = mlups,        &
        &               MLUPS_kernel = mlups_kernel, &
        &               imbalance    = imbalance,    &
        &               timerVal     = timerVal,     &
        &               general      = general       )
    end if

  end subroutine mus_perf_measure
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Performance results (MLUPs) are written to a file for statistical review
  !! The file-format is simple can be evaluated with gnuplot
  !!
  subroutine dump_timing( totalDens, DomSize, MLUPS, MLUPS_kernel, &
    &                     timerVal, imbalance, general )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: totalDens
    integer(kind=long_k), intent(in) :: DomSize
    real(kind=rk), intent(in) :: MLUPS, MLUPS_kernel, imbalance
    !> Max. timers of all process
    real(kind=rk), intent(in) :: timerVal(:)
    type(tem_general_type), intent(in) :: general
    ! -------------------------------------------------------------------- !
    real(kind=rk)          :: tMusubi
    logical                :: file_exists
    character(len=PathLen) :: filename
    integer                :: fileunit, iTimer, counter
    character(len=PathLen) :: header, output
    ! -------------------------------------------------------------------- !

    ! Header and output should match each other!
    header = ''
    output = ''

    write(header,'(a,a1,a15)' ) trim(header), '#', 'Revision'
    write(output,'(a,a16)')     trim(output), trim(general%solver%revision)

    write(header,'(a,a21)') trim(header), 'SimName'
    write(output,'(a,a21)') trim(output), ' '//trim(general%solver%simName)

    write(header,'(a,a15)') trim(header), 'DomSize'
    write(output,'(a,i15)') trim(output), DomSize

    write(header,'(a,a10)') trim(header), 'nProcs'
    write(output,'(a,i10)') trim(output), general%proc%comm_size
!$  write(header,'(a,a10)') trim(header), 'nThreads'
!$  write(output,'(a,i10)') trim(output), general%proc%nThreads

    write(header,'(a,a14)')    trim(header), 'MLUPs'
    write(output,'(a,en14.2)') trim(output), MLUPS

    write(header,'(a,a14)')    trim(header), 'MLUPs_kernel'
    write(output,'(a,en14.2)') trim(output), mlups_kernel

    write(header,'(a,a14)')   trim(header), 'imbalance(%)'
    write(output,'(a,F14.2)') trim(output), imbalance

    tMusubi = tem_getTimerVal( timerHandle = general%solver%timerHandle )
    write(header,'(a,a14)')    trim(header), 'timeMusubi'
    write(output,'(a,en14.4)') trim(output), tMusubi

    write(header,'(a,a10)') trim(header), 'maxIter'
    write(output,'(a,I10)') trim(output), general%simControl%now%iter

    write(header,'(a,a19)')    trim(header), 'totalDens'
    write(output,'(a,en19.9)') trim(output), totalDens

    ! Append all main timer values
    counter = 0
    do iTimer = mus_timerHandles%first, mus_timerHandles%last
      write(header,'(a,a16)') trim(header), &
        &           'time'//trim(tem_getTimerName(timerHandle = iTimer))
      counter = counter + 1
      write(output,'(a,f16.4)') &
        &           trim(output), timerVal(counter)
    end do

                                          !  12345678901234567890
    write(header,'(a,a12)')   trim(header), 'timeAux'
    write(output,'(a,f12.2)') trim(output), get_auxTime()

    write(header,'(a,a12)')   trim(header), 'timeRelax'
    write(output,'(a,f12.2)') trim(output), get_relaxTime()

    write(header,'(a,a12)')   trim(header), 'Comp(%)'
    write(output,'(a,f12.2)') trim(output), get_computeRatio()

    write(header,'(a,a12)')   trim(header), 'Comm(%)'
    write(output,'(a,f12.2)') trim(output), get_communicateRatio()

    write(header,'(a,a15)')   trim(header), 'BCbuffer(%)'
    write(output,'(a,f15.2)') trim(output), get_bcBufferRatio()

    write(header,'(a,a12)')   trim(header), 'BC(%)'
    write(output,'(a,f12.2)') trim(output), get_boundaryRatio()

    write(header,'(a,a12)')   trim(header), 'Intp(%)'
    write(output,'(a,f12.2)') trim(output), get_intpRatio()

    ! open and write to file
    filename=trim(general%timingFile)

    inquire(file=trim(filename),exist=file_exists)
    fileunit = newunit()
    open(unit=fileunit,file=trim(filename),position='append')

    if( .not. file_exists ) then
       write(fileunit,'(a)') trim(header)
    end if
    write(fileunit,'(a)') trim(output)
    close(fileunit)
    ! Write timing.res --------------------------------------------------------

  end subroutine dump_timing
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Calculate the MLUPS or MFLUPS for the simulation
  !!
  pure function calc_MLUPS( minLevel, maxLevel, scaleFactor, nElems, iter, &
    &                       time ) result( mlups )
    ! -------------------------------------------------------------------- !
    !> level range
    integer, intent(in) :: minLevel, maxLevel
    !> global parameter
    integer, intent(in) :: scaleFactor
    !> array of nElems levelwise
    integer(kind=long_k), intent(in) :: nElems(minLevel:maxLevel)
    !> number of iterations on maxLevel
    !! number of iteration on iLevel = iter / scaleFactor**(maxLevel-iLevel)
    integer, intent(in) :: iter
    !> time consumed for running iter iterations
    real(kind=rk), intent(in) :: time
    !> resulting mlups
    real(kind=rk) :: mlups
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    integer(kind=long_k) :: cellUpdates
    ! -------------------------------------------------------------------- !

    cellUpdates = 0_long_k

    ! MLUPS = sum( nElems(iLevel) * iter(iLevel) ) / mainLoopTime
    !       = sum( nElems(iLevel) * iter/scaleFactor(iLevel) ) / mainLoopTime
    !       = sum( nElems(iLevel) / scaleFactor(iLevel) ) * iter / mainLoopTime
    do iLevel = minLevel, maxLevel
      cellUpdates = cellUpdates + nElems(iLevel)                       &
        &         / int( scaleFactor**(maxLevel - iLevel), kind=long_k )
    enddo

    mlups = dble(cellUpdates*iter) / (time*1000000._rk)

  end function calc_MLUPS
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Check for the streaming layout.
  !!
  subroutine check_streaming_layout( minLevel, maxLevel )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: minLevel, maxLevel
    ! -------------------------------------------------------------------- !

    ! Check, if advection layout is PUSH. if so, abort, because the current code
    ! can only handle PUSH with single level
    if( "?StreamName?" == 'push') then
      write(logUnit(1),*)' Pushing pdfs to neighbor elements after collision '
      if( minLevel /= maxLevel ) then
        write(logUnit(1),*) ' ----------------------------------------------' &
          &                 //'------------ '
        write(logUnit(1),*) '           WARNING !!'
        write(logUnit(1),*) '    Push is not possible with multi-level ' &
          &                 //'execution.'
        write(logUnit(1),*) ' ----------------------------------------------' &
          &                 //'------------ '
        write(logUnit(1),*)

      endif
    else
      write(logUnit(1),*)' Pulling pdfs from neighbor elements before collision'
    endif

  end subroutine check_streaming_layout
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Write solver specific info to scratch file
  subroutine mus_writeSolverSpecInfo( scheme, params, rank, outUnit )
    ! -------------------------------------------------------------------- !
    !> scheme type
    type(mus_scheme_type), intent(in) :: scheme
    !> Contains physical convertion info and scaling type
    type(mus_param_type), intent(in) :: params
    !> Rank of the process either from
    !! global communicator or tracking output communicator
    integer, intent(in) :: rank
    !> unit to output solver info in lua format
    integer, intent(inout) :: outUnit
    ! -------------------------------------------------------------------- !
    type(aot_out_type) :: conf
    ! -------------------------------------------------------------------- !
    ! Write scratch file only from root process and
    ! outUnit is not initialized before
    if (rank == 0 .and. outUnit == -1) then

      outUnit = newUnit()
      open(unit=outUnit, status='scratch')
      call aot_out_open( put_conf = conf, outUnit = outUnit )

      ! write global parameter info
      call mus_param_out( me = params, conf = conf )

      ! write physics
      call mus_physics_out( me   = params%physics, &
        &                   conf = conf            )

      ! write schemes
      call mus_scheme_out( me   = scheme, &
        &                  conf = conf    )

      ! close conf
      call aot_out_close(put_conf = conf)
    end if

  end subroutine mus_writeSolverSpecInfo
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine dump_linear_partition( treeID, nElems, offset, myRank, iter )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: nElems
    integer(kind=long_k), intent(in) :: offset
    integer(kind=long_k), intent(in) :: treeID(1:nElems)
    integer, intent(in) :: myRank
    integer, intent(in) :: iter
    ! -------------------------------------------------------------------- !
    integer :: ii, level, previous, fileunit, minElem, maxElem, L
    integer :: block = 1024
    character(len=PathLen) :: filename
    logical                :: file_exists, toWrite
    character(len=65536) :: buffer
    ! -------------------------------------------------------------------- !

    write(logUnit(1), "(A)") 'Dump linear partition file.'
    write( filename, "(A,I5.5,A,I6.6,A)") &
      &              'elemlist_partition_p', myRank,'_t', iter, '.res'

    inquire(file=trim(filename),exist=file_exists)
    fileunit = newunit()
    open(unit=fileunit,file=trim(filename),position='append')

    write(fileunit,"(A10, A8)" ) '# iElem', 'level'

    buffer = ''

    ! always write the first element
    level = tem_LevelOf( treeID(1) )
    previous = level
    write(fileunit,"(A,I10,I8)") trim(buffer), 1+offset, level

    ! write elements by blocks
    do minElem = 2, nElems, block

      maxElem = min( minElem + block - 1, nElems )
      buffer = ''
      toWrite = .false.
      do ii = minElem, maxElem
        level = tem_LevelOf( treeID(ii) )
        ! only dump when my level is different from the level of my previous
        ! neighbor
        if ( level /= previous ) then
          write(buffer,"(A,I10,I8,A)") trim(buffer), ii+offset-1, previous, &
            &                          new_line("A")
          write(buffer,"(A,I10,I8,A)") trim(buffer), ii+offset, level, &
            &                          new_line("A")
          previous = level
          toWrite = .true.
        end if
      end do ! ii = minElem, maxElem
      ! the last character is a new_line, do not need to write it
      L = len(trim(buffer))
      if ( toWrite ) write(fileunit,"(A)") buffer(1:L-1)

    end do ! minElem = 2, nElems

    ! always write the last element
    level = tem_LevelOf( treeID(nElems) )
    write(fileunit,"(A,I10,I8)") trim(buffer), nElems+offset, level

    close(fileunit)

  end subroutine dump_linear_partition
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Output the min and max time across all ranks,
  !! which are spent on each boundary condition.
  subroutine mus_BC_timing( nBCs, bc_labels, comm )
    integer, intent(in) :: nBCs
    character(len=labelLen), intent(in) :: bc_labels(nBCs)
    integer, intent(in) :: comm

    real(kind=rk) :: bc_t(nBCs), min_t(nBCs), max_t(nBCs)
    integer :: ii, iError

    if ( nBCs > 0 ) then
      do ii = 1, nBCs
        bc_t(ii) = get_boundaryTime(ii)
      end do

      call mpi_reduce( bc_t, min_t, nBCs, rk_mpi, mpi_min, 0, comm, ierror )
      call mpi_reduce( bc_t, max_t, nBCs, rk_mpi, mpi_max, 0, comm, ierror )

      call tem_horizontalSpacer(fUnit=logUnit(5))
      write(logUnit(5), "(A)") 'Boundary timing information:'
      do ii = 1, nBCs
        write(logUnit(5), "(A12,2(A,F8.2))") trim(bc_labels(ii)), &
          &       ', min time: ', min_t(ii), ', max time: ', max_t(ii)
      end do
      call tem_horizontalSpacer(fUnit=logUnit(5))
    end if

  end subroutine mus_BC_timing
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine dump compute and BC timing for all ranks
  !! rank    nFluids     tCompute     nBCElems     tBC    tCPU    tMainLoop
  subroutine dump_bc_timing( proc, nFluids, nBCElems, DomSize, iter, send )
    type( tem_comm_env_type ), intent(in) :: proc
    integer, intent(in) :: nFluids
    integer, intent(in) :: nBCElems
    integer(kind=long_k), intent(in) :: DomSize
    integer, intent(in) :: iter
    type( tem_communication_type ), intent(in) :: send
    ! --------------------------------------------------------------------------
    character(len=PathLen) :: filename
    real(kind=rk) :: tMainLoop, tBC, tCompute, tBCBuffer, tComm, tCPU
    real(kind=rk) :: totalCPU, maxCPU
    integer :: fileunit, iError, msgsize, nLinks, totalProcs, totalLinks, ii
    character(len=OutLen) :: output
    integer :: ioStatus( mpi_status_size )
    ! --------------------------------------------------------------------------

    tMainLoop = get_mainLoopTime()
    tCompute  = get_computeTime()
    tBCBuffer = get_bcBufferTime()
    tBC       = get_boundaryTime()
    tComm     = get_communicateTime()
    tCPU      = tBC + tCompute + tBCBuffer

    nLinks = 0
    do ii = 1, send%nProcs
      nLinks = nLinks + send%buf_real(ii)%nVals
    end do

    call mpi_reduce( tCPU, totalCPU, 1, mpi_double_precision, &
                     mpi_sum, 0, proc%comm, ierror            )
    call mpi_reduce( tCPU,   maxCPU, 1, mpi_double_precision, &
                     mpi_max, 0, proc%comm, ierror            )
    totalCPU = totalCPU / dble( proc%comm_size )

    call mpi_reduce( send%nProcs, totalProcs, 1, mpi_integer, &
                     mpi_sum, 0, proc%comm, ierror            )
    call mpi_reduce( nLinks, totalLinks, 1, mpi_integer, &
                     mpi_sum, 0, proc%comm, ierror       )
    totalProcs = totalProcs / proc%comm_size
    totalLinks = totalLinks / proc%comm_size

    write(filename, "(A)") 'bc_timing.res'
    if ( proc%isRoot ) then

      fileUnit = newunit()
      open(unit=fileunit,file=trim(filename),position='append')

      ! Write Header ---------------------------------------------------
      !                       12345678901234567890
      write(fileUnit, '(A)') ''
      write(fileUnit, '(A,I0)') '# DomSize: ', DomSize
      write(fileUnit, '(A,I0)') '# Iteration: ', iter
      write(fileUnit, '(A,ES10.2)') '# tMainLoop: ', tMainLoop
      write(fileUnit, '(2(A,ES10.2),A,F5.1)') '# tCPU average/max = ', &
        &                              totalCPU, ' / ',maxCPU, &
        &                              ' = ', (totalCPU/maxCPU)*100._rk
      write(fileUnit, '(A,I0,A,I0)') '# Send mean nProcs: ', totalProcs, &
        &                           ', mean nLinks: ', totalLinks
      write(fileUnit, '(A6,7A10)') &
        &                    '# rank', &
        &                    'nFluids', &
        &                    'tCompute', &
        &                    'nBCElems', &
        &                    'tBCBuffer', &
        &                    'tBC', &
        &                    'tCPU', &
        &                    'tComm'
      close(fileunit)
    end if

    ! Make sure root finishes writing header
    call MPI_BARRIER( proc%comm, iError )

    ! Write Data ---------------------------------------------------
    output = ''
    ! Write data into string
    write(output,"(A,I6)")  trim(output), proc%rank
    write(output,"(A,I10)") trim(output), nFluids
    write(output,"(A,ES10.2)") trim(output), tCompute
    write(output,"(A,I10)") trim(output), nBCElems
    write(output,"(A,ES10.2)") trim(output), tBCBuffer
    write(output,"(A,ES10.2)") trim(output), tBC
    write(output,"(A,ES10.2)") trim(output), (tBC + tCompute + tBCBuffer)
    write(output,"(A,ES10.2)") trim(output), tComm
    write(output,"(A,A)")  trim(output), new_line('A')

    call MPI_File_open( proc%comm, trim(filename),                           &
      &                 ior(MPI_MODE_APPEND,MPI_MODE_WRONLY), MPI_INFO_NULL, &
      &                 fileunit, iError                                     )
    call check_MPI_error( iError, 'Open BC timing file')

    msgsize = len_trim( output )
    call MPI_File_write_ordered(fileunit, trim(output), msgsize, &
      &                         MPI_CHARACTER, iostatus, iError  )
    call MPI_File_close(fileunit, iError)
    call check_MPI_error( iError, 'Close BC timing file')

  end subroutine dump_bc_timing
  ! ------------------------------------------------------------------------ !

end module mus_tools_module
! **************************************************************************** !
