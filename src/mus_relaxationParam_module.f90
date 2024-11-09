! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2021 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2020-2021 Harald Klimach <harald.klimach@dlr.de>
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
!> This module contains the data type for MRT.
!!
!! Also provides functions and routines to set relaxation parameters for
!! single-component and multispecies.
!!
module mus_relaxationParam_module
  use mpi
  ! include treelm modules
  use env_module,               only: rk, io_buffer_size, eps_single_k
  use tem_param_module,         only: cs2inv, div2_3
  use tem_time_module,          only: tem_time_type
  use tem_temporal_module,      only: tem_temporal_type, tem_temporal_for
  use tem_varSys_module,        only: tem_varSys_type
  use tem_spacetime_fun_module, only: tem_spacetime_fun_type, tem_spacetime_for
  use tem_logging_module,       only: logUnit
  use tem_aux_module,           only: tem_abort
  use tem_grow_array_module,    only: grw_realArray_type, init
  use tem_comm_env_module,      only: tem_comm_env_type
  use tem_general_module,       only: tem_general_type
  use tem_status_module,        only: tem_stat_nonPhysical

  ! include musubi modules
  use mus_physics_module,       only: mus_convertFac_type
  use mus_scheme_header_module, only: mus_scheme_header_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_moments_type_module,  only: mus_moment_type
  use mus_derVarPos_module,     only: mus_derVarPos_type
  use mus_turbulence_module,    only: mus_turbulence_type, mus_turb_calcVisc
  use mus_gradData_module,      only: mus_gradData_type
  use mus_nonNewtonian_module,  only: mus_nNwtn_type
  use mus_gradData_module,      only: mus_Grad_type

  implicit none

  private

  public :: mus_viscosity_type
  public :: mus_update_viscKine
  public :: mus_update_relaxParamKine
  public :: mus_update_relaxParamFromViscSTfun
  public :: mus_init_relaxParam
  public :: mus_calcOmegaFromVisc
  public :: mus_check_omegaKine

  !> Contains relaxation parameter for a level
  type mus_relaxationParam_type
    !> Relaxation parameter computed from viscosity
    !! For kinematic viscosity, if turbulence is active, this omega refers to
    !! effective omega which is omega_bg + omega_turb
    !! size: nElems_solve
    real(kind=rk), allocatable :: val(:)
  end type mus_relaxationParam_type

  !> Contains STfun of viscosity variable and relaxation parameter for each
  !! level
  type mus_viscosity_type
    !> space-time function
    type(tem_spacetime_fun_type) :: STfun
    !> viscosity value evaluated from STfun
    type(grw_realArray_type), allocatable :: dataOnLvl(:)
    !> relaxation paramter omega for each level
    type(mus_relaxationParam_type), allocatable :: omLvl(:)
  end type mus_viscosity_type


contains


  ! ------------------------------------------------------------------------ !
  !> This routine initialize relaxation parameter
  subroutine mus_init_relaxParam( omLvl, minLevel, maxLevel, nElems )
    ! -------------------------------------------------------------------- !
    !> relaxation paramter
    type(mus_relaxationParam_type), allocatable, intent(out) :: omLvl(:)
    !> minlevel and maxLevel
    integer, intent(in) :: minLevel, maxLevel
    !> number of local elements per level
    integer, intent(in) :: nElems(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    ! -------------------------------------------------------------------- !
    allocate(omLvl(minLevel:maxLevel))
    do iLevel = minLevel, maxLevel
      allocate( omLvl(iLevel)%val( nElems(iLevel) ) )
      omLvl(iLevel)%val = -100.0_rk
    end do

  end subroutine mus_init_relaxParam
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Update kinematic viscosity from STfun and calculate turbulent viscosity
  !! from velocity gradient or nonEqPDF
  !! Viscosity obtained from this routine are normalized to the level
  subroutine mus_update_viscKine( viscKine, state, neigh, auxField, gradData, &
    &                             layout, baryOfTotal, tNow, nSize, nFluids,  &
    &                             nGhostFromCoarser, nGhostFromFiner, nHalo,  &
    &                             varSys, iLevel, convFac, dxL, dtL,          &
    &                             derVarPos, turb, nNwtn, Grad                )
    ! ------------------------------------------------------------------------ !
    !> Kinematic viscosity
    type(mus_viscosity_type), intent(inout) :: viscKine
    !> bary of treeID in total list
    real(kind=rk), intent(in) :: baryOfTotal(:,:)
    !> number of elements in state array
    integer, intent(in) :: nSize
    !> number of fluid elements in state array
    integer, intent(in) :: nFluids
    !> Number of ghostFromCoarser element in state array
    integer, intent(in) :: nGhostFromCoarser
    !> Number of ghostFromFiner element in state array
    integer, intent(in) :: nGhostFromFiner
    !> Number of halo element in state array
    integer, intent(in) :: nHalo
    !> state array
    real(kind=rk), intent(in) :: state(:)
    !> neighbor connectivity array
    integer, intent(in) :: neigh(:)
    !> Auxiliary field variable array
    real(kind=rk), intent(in) :: auxField(:)
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> stencil layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> current level
    integer, intent(in) :: iLevel
    !> current simulation time
    type(tem_time_type), intent(in) :: tNow
    !> reference physical conversion factors for current level
    type(mus_convertFac_type), intent(in) :: convFac
    !> lattice element size in current level
    real(kind=rk), intent(in) :: dxL
    !> lattice time step size in current level
    real(kind=rk), intent(in) :: dtL
    !> variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> contains position of in varSys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> turbulence type
    type(mus_turbulence_type), intent(inout) :: turb
    !> non-Newtonian type
    type(mus_nNwtn_type), intent(in) :: nNwtn
    !> Object that contains pointers to calculate gradients
    type(mus_Grad_type), intent(in) :: Grad
    ! -------------------------------------------------------------------- !
    integer :: nScalars, nAuxScalars, nSolve, densPos, velPos(3)
    ! -------------------------------------------------------------------- !
    nSolve = nFluids + nGhostFromCoarser + nGhostFromFiner + nHalo
    nScalars = varSys%nScalars
    nAuxScalars = varSys%nAuxScalars
    densPos = varSys%method%val(derVarPos%density)%auxField_varPos(1)
    velPos = varSys%method%val(derVarPos%velocity)%auxField_varPos(:)

    ! Update physical kinematic viscosity acccording to non-Newtonian model
    if (nNwtn%active) then
      ! calculate lattice kinemaic viscosity from local shear rate according to
      ! non-Newtonian model.
      ! Output of this routine is vL_l/dtL_l
      call nNwtn%calcVisc(viscKine     = viscKine%dataOnLvl(iLevel) &
        &                                        %val(1:nSolve),    &
        &                 omega        = viscKine%omLvl(iLevel)     &
        &                                        %val(1:nSolve),    &
        &                 state        = state,                     &
        &                 neigh        = neigh,                     &
        &                 auxField     = auxField,                  &
        &                 densPos      = densPos,                   &
        &                 velPos       = velPos,                    &
        &                 nSize        = nSize,                     &
        &                 nSolve       = nSolve,                    &
        &                 nScalars     = nScalars,                  &
        &                 nAuxScalars  = nAuxScalars,               &
        &                 layout       = layout,                    &
        &                 convFac      = convFac                    )

    else
      ! background fluid kinematic viscosity for nFluid+nGhost
      viscKine%dataOnLvl(iLevel)%val(1:nSolve)                    &
        & = tem_spacetime_for( me      = viscKine%STfun,          &
        &                      coord   = baryOfTotal(1:nSolve,:), &
        &                      time    = tNow,                    &
        &                      n       = nSolve                   )

      ! convert physical viscosity to lattice unit
      ! convFac%visc is the reference physical viscosity on current level
      ! i.e. (dxP_l)^2/dtP_l
      ! Dividing physical viscosity with the viscRef gives vL_l/dtL_l
      viscKine%dataOnLvl(iLevel)%val(:) = viscKine%dataOnLvl(iLevel)%val(:) &
        &                               / convFac%visc
    end if


    if (turb%active) then
      !> calculate turbulence viscosity for nFluids and ghostFromCoarser,
      !! for ghostFromFiner elements, turb visc is interpolated
      nSolve = nFluids + nGhostFromCoarser
      call mus_turb_calcVisc(turbData     = turb%dataOnLvl(iLevel),    &
        &                    turbConfig   = turb%config,               &
        &                    calcTurbVisc = turb%calcVisc,             &
        &                    state        = state,                     &
        &                    neigh        = neigh,                     &
        &                    auxField     = auxField,                  &
        &                    gradData     = gradData,                  &
        &                    densPos      = densPos,                   &
        &                    velPos       = velPos,                    &
        &                    nSize        = nSize,                     &
        &                    nSolve       = nSolve,                    &
        &                    nScalars     = nScalars,                  &
        &                    nAuxScalars  = nAuxScalars,               &
        &                    layout       = layout,                    &
        &                    viscKine     = viscKine%dataOnLvl(iLevel) &
        &                                           %val(:),           &
        &                    dxL          = dxL,                       &
        &                    dtL          = dtL,                       &
        &                    Grad         = Grad                       )
    end if

  end subroutine mus_update_viscKine
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Update kinematic relaxation parameter from viscosity and check omega
  subroutine mus_update_relaxParamKine(viscKine, turb, nSolve, iLevel)
    ! -------------------------------------------------------------------- !
    !> Kinematic viscosity
    type(mus_viscosity_type), intent(inout) :: viscKine
    !> Number of elements to solve in compute kernel
    integer, intent(in) :: nSolve
    !> turbulence type
    type(mus_turbulence_type), intent(in) :: turb
    !> current level
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: iElem
    real(kind=rk) :: tot_visc
    ! -------------------------------------------------------------------- !

    if (turb%active) then
      do iElem = 1, nSolve
        ! total normalized viscosity
        tot_visc = viscKine%dataOnLvl(iLevel)%val(iElem) &
          &      + turb%dataOnLvl(iLevel)%visc(iElem)
        ! compute omega
        viscKine%omLvl(iLevel)%val(iElem) = 1.0_rk / ( cs2inv * tot_visc &
          &                                           + 0.5_rk           )
      end do
    else

      ! compute omega from kinematic viscosity
      do iElem = 1, nSolve
        viscKine%omLvl(iLevel)%val(iElem)                               &
          & = 1.0_rk / ( cs2inv * viscKine%dataOnLvl(iLevel)%val(iElem) &
          &             + 0.5_rk                                        )
      end do
    end if

  end subroutine mus_update_relaxParamKine
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine is used to initialize relaxation paramter and update
  !! bulk viscosity at every time step
  !! Bulk visocisty is defined as space-time function to apply
  !! ramping and spatial sponge in bulk viscosity
  subroutine mus_update_relaxParamFromViscSTfun( omega, visc, viscSTfun,    &
    &                                            nElems, baryOfTotal, tNow, &
    &                                            viscRef                    )
    ! -------------------------------------------------------------------- !
    !> relaxation parameter
    real(kind=rk), intent(inout) :: omega(:)
    !> Kinematic viscosity
    real(kind=rk), intent(inout) :: visc(:)
    !> viscosity space-time function
    type(tem_spacetime_fun_type), intent(in) :: viscSTfun
    !> Number of local elements including halos
    integer, intent(in) :: nElems
    !> baryID of total list
    real(kind=rk), intent(in) :: baryOfTotal(:,:)
    !> current simulation time
    type(tem_time_type), intent(in) :: tNow
    !> reference physical viscosity on current level i.e. (dxP_l)^2/dtP_l
    !! Dividing physical viscosity with the viscRef gives vL_l/dtL_l
    real(kind=rk), intent(in) :: viscRef
    !> lattice time step size in current level
    ! real(kind=rk), intent(in) :: dtL
    ! -------------------------------------------------------------------- !
    ! viscosity values from space-time function
    integer :: iChunk, nChunks, nChunkElems, elemoff
    integer :: minBuf, maxBuf
    ! -------------------------------------------------------------------- !
    !\todo KM: 'Optimize this routine for constant viscosity'
    ! find chunksize and number of chunks required for initialzation
    nChunks = ceiling( real(nElems, rk) / real(io_buffer_size, rk) )

    do iChunk = 1, nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ( (iChunk-1)*io_buffer_size )
      nChunkElems = min(io_buffer_size, nElems - elemOff)
      minBuf = elemOff+1
      maxBuf = elemOff+nChunkElems

      ! background viscosity
      visc(minBuf:maxBuf) = tem_spacetime_for(                           &
        &                       me      = viscSTfun,                     &
        &                       coord   = baryOfTotal(minBuf:maxBuf, :), &
        &                       time    = tNow,                          &
        &                       n       = nChunkElems                    )

      ! convert viscosity to lattice unit
      ! viscRef is scaled with current level dtL
      visc(minBuf:maxBuf) = visc(minBuf:maxBuf) / viscRef

      ! compute omega
      !do iElem = 1, nChunkElems
      !  omega(elemOff+iElem) = mus_calcOmegaFromVisc(visc(iElem))
      !end do
      omega(minBuf:maxBuf) = mus_calcOmegaFromVisc(visc(minBuf:maxBuf))

    end do

  end subroutine mus_update_relaxParamFromViscSTfun
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This function compute relaxation paramter omega from viscosity
  elemental function mus_calcOmegaFromVisc(visc) result(omega)
    ! -------------------------------------------------------------------- !
    !> scaled lattice viscosity i.e vL_c/dtL
    real(kind=rk), intent(in) :: visc
    !> lattice time step size in current level
    !!real(kind=rk), intent(in) :: dtL
    !> output: relaxation parameter omega
    real(kind=rk) :: omega
    ! -------------------------------------------------------------------- !
    !omega = 1.0_rk / ( cs2inv * visc / dtL + 0.5_rk )
    omega = 1.0_rk / ( cs2inv * visc + 0.5_rk )
  end function mus_calcOmegaFromVisc
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine checks whether omega is within the stability limit.
  !! If not it will terminate the simulation with error message.
  !! Using limits given in
  !! Tölke, J., Freudiger, S., & Krafczyk, M. (2006). An adaptive scheme using
  !! hierarchical grids for lattice Boltzmann multi-phase flow simulations.
  !! Computers & Fluids, 35(8–9), 820–830.
  !! For BGK: 2/3 < omega < 1.976
  !! For MRT: 2/3 < omega < 1.999
  subroutine mus_check_omegaKine( schemeRelax, omLvlKine, nSolve, minLevel, &
    &                             maxLevel, general                         )
    ! -------------------------------------------------------------------- !
    !> minlevel and maxlevel
    integer, intent(in) :: minLevel, maxLevel
    !> scheme relaxation type
    character(len=*), intent(in) :: schemeRelax
    !> array of kinematic relaxation parameter on all levels
    type(mus_relaxationParam_type), intent(in) :: omLvlKine(minLevel:maxLevel)
    !> Number of elements to solve in compute kernel
    integer, intent(in) :: nSolve(minLevel:maxLevel)
    !> Contains proc, simControl, solveHead
    type(tem_general_type), intent(inout) :: general
    ! -------------------------------------------------------------------- !
    integer :: iLevel, iError
    !> minimum omega
    real(kind=rk) :: om_min, om_max, glob_om_min, glob_om_max
    ! -------------------------------------------------------------------- !
    do iLevel = minLevel, maxLevel
      ! minimum omega
      om_min = minval( omLvlKine(iLevel)%val( 1 : nSolve(iLevel) ) )

      ! maximum omega
      om_max = maxval( omLvlKine(iLevel)%val( 1 : nSolve(iLevel) ) )

      ! get global min and max omega
      call mpi_reduce( om_min, glob_om_min, 1, mpi_double_precision, &
                       mpi_min, 0, general%proc%comm, ierror         )
      call mpi_reduce( om_max, glob_om_max, 1, mpi_double_precision, &
                       mpi_max, 0, general%proc%comm, ierror         )

      if (general%proc%isRoot) then
        if (glob_om_min < div2_3) then
          write(logUnit(1),'(A,F10.5)') 'Error: Kinematic omega < 2/3:', &
            &                           glob_om_min
          write(logUnit(1),'(A,I0)') 'On level:', iLevel
          write(logUnit(1),*) 'Solution: Decrease kinematic viscosity'
          general%simControl%status%bits(tem_stat_nonPhysical) = .true.
        end if

        select case (trim(schemeRelax))
        case ('bgk')
          if (glob_om_max > 1.976_rk) then
            write(logUnit(1),'(A,F10.5)') 'WARNING: Kinematic omega > 1.976:', &
              &                           glob_om_max
            write(logUnit(1),'(A,I0)') 'On level:', iLevel
            write(logUnit(1),*) 'Solution: Increase kinematic viscosity,'
            write(logUnit(1),*) 'especially near BC and level jumps.'
          end if
        case ( 'mrt','trt','cumulant','cumulant_extended','hrr_bgk','rr_bgk', &
          &    'prr_bgk','r_bgk', 'drt_bgk', 'rr_bgk_corrected',              &
          &    'cumulant_extended_generic', 'hrr_bgk_corrected',              &
          &    'prr_bgk_corrected'                                            )
          if (glob_om_max > 1.999_rk) then
            write(logUnit(1),'(A,F10.5)') 'WARNING: Kinematic omega > 1.999:', &
              &                           glob_om_max
            write(logUnit(1),'(A,I0)') 'On level:', iLevel
            write(logUnit(1),*) 'Solution: Increase kinematic viscosity,'
            write(logUnit(1),*) 'especially near BC and level jumps.'
          end if
        case default
          call tem_abort('Error: Unknown scheme relaxation in check omega')
        end select

        write(logUnit(10), '(A,I0)') 'On level:', iLevel
        write(logUnit(10),'(A,F10.5,A,F10.5,A)') &
          &  'Kinematic omega (min,max):(', glob_om_min, ',', glob_om_max, ')'
      end if

      if ( any(abs( omLvlKine(iLevel)%val(1:nSolve(iLevel)) - 1.0_rk ) &
        &      < 1000 * eps_single_k ) ) then
        write(*,"(A,I0)")' Attention: Omega on level: ', iLevel
        write(*,"(A)")   ' is chosen close to 1.0'
        write(*,"(A)")   ' This can cause program to crash '  &
          &            //'when it calculates strain rate or ' &
          &            //'shear stress.'
      end if
    end do

  end subroutine mus_check_omegaKine
  ! ------------------------------------------------------------------------ !

end module mus_relaxationParam_module
