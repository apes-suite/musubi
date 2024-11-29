! Copyright (c) 2018-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Jana Gericke <jana.gericke@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2021-2023 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
! ************************************************************************** !
!> author: Kannan Masilamani
!! Average Interpolation of flow quantities between different grid levels
!!
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
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! Make sure loglvl is defined to an integer value.
! Usually this should be defined on the command line with -Dloglvl=
! to some value.
module mus_interpolate_average_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,              only: rk
  use tem_aux_module,          only: tem_abort
  use tem_element_module,      only: eT_GhostFromCoarser, &
    &                                eT_ghostFromFiner
  use tem_param_module,        only: cs2inv, cs2, PI, div1_2, div1_9, div4_9,&
    &                                div1_36, div2_3, div2_9, div1_18, rho0, rho0Inv
  use tem_comm_env_module,     only: tem_comm_env_type
  use tem_debug_module,        only: dbgUnit
  use tem_construction_module, only: tem_levelDesc_type
  use tem_matrix_module,       only: tem_matrix_type
  use tem_stencil_module,      only: tem_stencilHeader_type
  use tem_logging_module,      only: logUnit
  use tem_varSys_module,       only: tem_varSys_type
  use tem_time_module,         only: tem_time_type


  ! include musubi modules
  use mus_interpolate_debug_module,  only: TGV_2D
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_physics_module,            only: mus_physics_type
  use mus_interpolate_header_module, only: mus_interpolation_method_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_fluid_module,              only: mus_fluid_type
  use mus_relaxationParam_module,    only: mus_calcOmegaFromVisc
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_derivedQuantities_module2,  only: getNonEqFac_intp_fine_to_coarse, &
    &                                       getNonEqFac_intp_coarse_to_fine
  use mus_varSys_module,        only: mus_varSys_data_type

  implicit none

  private

  public :: fillArbiMyGhostsFromFiner_avg
  public :: fillMyGhostsFromFiner_avg_feq_fneq
  public :: fillMyGhostsFromFiner_avgLES_feq_fneq
  public :: fillMyGhostsFromFiner_avg2D_feq_fneq
  public :: fillArbiFinerGhostsFromMe_weighAvg
  public :: fillFinerGhostsFromMe_weighAvg_feq_fneq
  public :: fillFinerGhostsFromMe_weighAvgLES_feq_fneq
  public :: fillFinerGhostsFromMe_weighAvg2D_feq_fneq

  contains


! **************************************************************************** !
  !> Interpolate auxiliary field from fine source to coarse target
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine_arbitraryVal]] in intp/[[mus_interpolate_header_module]].f90
  !! in order to be callable via
  !! [[mus_interpolation_method_type:do_intpArbiVal]] function pointer.
  subroutine fillArbiMyGhostsFromFiner_avg( method, tLevelDesc, level,     &
    &                                       stencil, sVal, tVal, nTargets, &
    &                                       targetList, nScalars           )
    ! ------------------------------------------------------------------ !
    class(mus_interpolation_method_type), intent(inout) :: method

    !> my refinement level
    integer, intent(in) :: level

    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil

    !> State vector of SOURCE FLUID elements
    real(kind=rk), intent(in) :: sVal(:)

    !> State vector of TARGET GHOST elements
    real(kind=rk), intent(inout) :: tVal(:)

    !> level descriptor on target level
    type( tem_levelDesc_type ), intent(in) :: tLevelDesc

    !> List of target elements ( their position in depSource list )
    integer, intent(in) :: nTargets
    integer, intent(in) :: targetList(nTargets)

    !> number of scalars to interpolate
    integer, intent(in) :: nScalars
    ! ------------------------------------------------------------------ !
    ! --------------------------------------------------------------------------
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    real(kind=rk) :: inv_nSourceElems
    real(kind=rk) :: tArbi(nScalars)  ! target auxField
    real(kind=rk) :: sArbi(nScalars)  ! temp source ArbiField
    ! --------------------------------------------------------------------------
    sourceLevel = level + 1
    targetLevel = level

    ! Treat all coarse target elements
    do indElem = 1, nTargets

      iElem = targetList( indElem )

      ! Read the target element
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromFiner)

      ! Get how many fine source elements we have for interpolation.
      nSourceElems = tLevelDesc%depFromFiner( iElem )%elem%nVals
      inv_nSourceElems = 1.0_rk / dble(nSourceElems)

      tArbi = 0.0_rk
      ! Now loop over all fine source elements for this target:
      do iSourceElem = 1, nSourceElems

        ! Get the source element
        sourceElem = tLevelDesc%depFromFiner( iElem ) &
          &                    %elem%val( iSourceElem )

        ! Get souce auxilary variables
        sArbi(1:nScalars) = sVal((sourceElem-1)*nScalars+1:sourceElem*nScalars)

        ! target aux
        tArbi = sArbi + tArbi

      end do  ! iSourceElem

      ! interpolate all auxiliary variables by average
      tVal((targetElem-1)*nScalars+1: targetElem*nScalars) &
        & = tArbi(:) * inv_nSourceElems

    enddo

  end subroutine fillArbiMyGhostsFromFiner_avg
! **************************************************************************** !

! **************************************************************************** !
  !> [Average interpolation](../page/features/intp_methods.html) of ghostFromFiner
  !! The interpolation procedure used in this routine is:\n
  !! 1. Calculate Equilibrium and nonEquilibrium
  !! 2. Compute scaling factor
  !! 3. calculate target: Eq + Scale * nonEquilibrium
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine fillMyGhostsFromFiner_avg_feq_fneq( method, fieldProp,           &
    &      tLevelDesc, level, sState, sNeigh, snSize, sAuxField, tState,      &
    &      tNeigh, tnSize, layout, nTargets, targetList, physics, time,       &
    &      varSys, derVarPos                                                  )
    ! -------------------------------------------------------------------- !
    class(mus_interpolation_method_type), intent(inout) :: method

    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), target, intent(in) :: fieldProp(:)

    !> level descriptor on target level
    type( tem_levelDesc_type ), intent(in) :: tLevelDesc

    !> my refinement level
    integer, intent(in) :: level

    !> State vector of SOURCE FLUID elements
    real(kind=rk), intent(in) :: sState(:)
    integer, intent(in) :: sNeigh(:)
    integer, intent(in) :: snSize

    !> State vector of TARGET GHOST elements
    real(kind=rk), intent(inout) :: tState(:)
    integer, intent(in) :: tNeigh(:)
    integer, intent(in) :: tnSize

    !> AuxField variable to read rho and vel from source elements
    real(kind=rk), intent(inout) :: sAuxField(:)

    !> the layout used
    type( mus_scheme_layout_type ), intent(in) :: layout

    !> List of target elements ( their position in depSource list )
    integer, intent(in) :: nTargets
    integer, intent(in) :: targetList(nTargets)

    !> physics type to convert lattice to physics SI unit and vice versa
    !! @todo: This can be replaced by scale factor by level
    type( mus_physics_type ), intent(in) :: physics

    !> time required to compute viscosity on target element barycenter
    type(tem_time_type), intent(in) :: time

    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys

    !> position of all derive variable in varSys for all fields
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iDir           ! current direction (discrete velocity) for loop
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    real(kind=rk) :: f( layout%fStencil%QQ ), f_eq( layout%fStencil%QQ )
    integer :: QQ
    real(kind=rk) :: inv_nSourceElems
    real(kind=rk) :: t_f_eq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: t_f_neq( layout%fStencil%QQ )  ! temp pdf calculation
    type(mus_fluid_type), pointer :: fluid
    real(kind=rk) :: cVisc, cOmegaKine, fOmegaKine, rho, vel(3), fac
    integer :: elemOff, dens_pos, vel_pos(3), nScalars
    ! --------------------------------------------------------------------------


    fluid => fieldProp(1)%fluid
    nScalars = varSys%nScalars
    QQ = layout%fStencil%QQ

    sourceLevel = level + 1
    targetLevel = level

    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)


    ! Treat all coarse target elements
    do indElem = 1, nTargets

      iElem = targetList( indElem )

      ! Read the target element
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromFiner)


      ! Get how many fine source elements we have for interpolation.
      nSourceElems = tLevelDesc%depFromFiner( iElem )%elem%nVals
      inv_nSourceElems = 1.0_rk / dble(nSourceElems)

      t_f_eq = 0._rk
      t_f_neq = 0._rk
      ! Now loop over all fine source elements for this target:
      do iSourceElem = 1, nSourceElems

        ! Get the source element
        sourceElem = tLevelDesc%depFromFiner( iElem ) &
          &                    %elem%val( iSourceElem )

        do iDir = 1, QQ
        ! this is post collision
          f(iDir) = sState( ( sourceelem-1)* nscalars+ idir+( 1-1)* qq)
        enddo

        ! element offset for auxField array
        elemOff = (sourceElem-1)*varSys%nAuxScalars

        ! local density
        rho = sAuxField(elemOff + dens_pos)
        ! local x-, y-, z-velocity
        vel(1) = sAuxField(elemOff + vel_pos(1))
        vel(2) = sAuxField(elemOff + vel_pos(2))
        vel(3) = sAuxField(elemOff + vel_pos(3))


        f_eq(:) = layout%quantities%pdfEq_ptr( rho = rho,               &
          &                                    vel = vel,               &
          &                                    QQ = QQ                  )

        t_f_eq(:) = t_f_eq(:) + f_eq(:)
        t_f_neq(:) = t_f_neq(:) + (f(:) - f_eq(:))

      end do  ! iSourceElem

      ! get normalized kinematic viscosity on target (coarse) element
      cVisc = fluid%viscKine%dataOnLvl(targetLevel)%val(targetElem)

      ! relation between coarse and fine grid kinematic viscosity:
      ! v^s_f = 2 v^s_c
      ! calculate omega on source and target level
      fOmegaKine = mus_calcOmegaFromVisc(2_rk*cVisc)
      cOmegaKine = mus_calcOmegaFromVisc(cVisc)

      !evaluate scaling factor
      fac = getNonEqFac_intp_fine_to_coarse( cOmegaKine, fOmegaKine)

      ! interpolate all pdfs by average and scale accordingly
      do iDir = 1, QQ
        t_f_eq(iDir) = t_f_eq(iDir) * inv_nSourceElems
        t_f_neq(iDir) = t_f_neq(iDir) * inv_nSourceElems * fac
        f(iDir) = t_f_eq(iDir) + t_f_neq(iDir)
        ! Now write the resulting pdf in the current direction to the target
        tState( ( targetelem-1)* nscalars+ idir+( 1-1)* qq) = &
          & f(iDir)
      enddo


    enddo

  end subroutine fillMyGhostsFromFiner_avg_feq_fneq
! **************************************************************************** !

! **************************************************************************** !
  !> [Average interpolation](../page/features/intp_methods.html) of ghostFromFiner
  !! The interpolation procedure used in this routine is:\n
  !! 1. Calculate Equilibrium and nonEquilibrium
  !! 2. Compute scaling factor
  !! 3. calculate target: Eq + Scale * nonEquilibrium
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine fillMyGhostsFromFiner_avgLES_feq_fneq( method, fieldProp,     &
    &      tLevelDesc, level, sState, sNeigh, snSize, sAuxField, tState,   &
    &      tNeigh, tnSize, layout, nTargets, targetList, physics, time,    &
    &      varSys, derVarPos                                               )
    ! -------------------------------------------------------------------- !
    class(mus_interpolation_method_type), intent(inout) :: method

    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), target, intent(in) :: fieldProp(:)

    !> level descriptor on target level
    type( tem_levelDesc_type ), intent(in) :: tLevelDesc

    !> my refinement level
    integer, intent(in) :: level

    !> State vector of SOURCE FLUID elements
    real(kind=rk), intent(in) :: sState(:)
    integer, intent(in) :: sNeigh(:)
    integer, intent(in) :: snSize

    !> State vector of TARGET GHOST elements
    real(kind=rk), intent(inout) :: tState(:)
    integer, intent(in) :: tNeigh(:)
    integer, intent(in) :: tnSize

    !> AuxField variable to read rho and vel from source elements
    real(kind=rk), intent(inout) :: sAuxField(:)

    !> the layout used
    type( mus_scheme_layout_type ), intent(in) :: layout

    !> List of target elements ( their position in depSource list )
    integer, intent(in) :: nTargets
    integer, intent(in) :: targetList(nTargets)

    !> physics type to convert lattice to physics SI unit and vice versa
    !! @todo: This can be replaced by scale factor by level
    type( mus_physics_type ), intent(in) :: physics

    !> time required to compute viscosity on target element barycenter
    type(tem_time_type), intent(in) :: time

    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys

    !> position of all derive variable in varSys for all fields
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iDir           ! current direction (discrete velocity) for loop
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    real(kind=rk) :: f( layout%fStencil%QQ ), f_eq( layout%fStencil%QQ )
    integer :: QQ
    real(kind=rk) :: inv_nSourceElems
    real(kind=rk) :: t_f_eq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: t_f_neq( layout%fStencil%QQ )  ! temp pdf calculation
    type(mus_fluid_type), pointer :: fluid
    real(kind=rk) :: cVisc, cOmegaKine, fOmegaKine, rho, vel(3), fac
    real(kind=rk) :: tTurbVisc
    integer :: elemOff, dens_pos, vel_pos(3), nScalars
    ! --------------------------------------------------------------------------


    fluid => fieldProp(1)%fluid
    nScalars = varSys%nScalars
    QQ = layout%fStencil%QQ

    sourceLevel = level + 1
    targetLevel = level

    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)


    ! Treat all coarse target elements
    do indElem = 1, nTargets

      iElem = targetList( indElem )

      ! Read the target element
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromFiner)


      ! Get how many fine source elements we have for interpolation.
      nSourceElems = tLevelDesc%depFromFiner( iElem )%elem%nVals
      inv_nSourceElems = 1.0_rk / dble(nSourceElems)

      t_f_eq = 0._rk
      t_f_neq = 0._rk
      tTurbVisc = 0.0_rk
      ! Now loop over all fine source elements for this target:
      do iSourceElem = 1, nSourceElems

        ! Get the source element
        sourceElem = tLevelDesc%depFromFiner( iElem ) &
          &                    %elem%val( iSourceElem )

        do iDir = 1, QQ
        ! this is post collision
          f(iDir) = sState( ( sourceelem-1)* nscalars+ idir+( 1-1)* qq)
        enddo

        ! element offset for auxField array
        elemOff = (sourceElem-1)*varSys%nAuxScalars

        ! local density
        rho = sAuxField(elemOff + dens_pos)
        ! local x-, y-, z-velocity
        vel(1) = sAuxField(elemOff + vel_pos(1))
        vel(2) = sAuxField(elemOff + vel_pos(2))
        vel(3) = sAuxField(elemOff + vel_pos(3))


        f_eq(:) = layout%quantities%pdfEq_ptr( rho = rho,               &
          &                                    vel = vel,               &
          &                                    QQ = QQ                  )

        t_f_eq(:) = t_f_eq(:) + f_eq(:)
        t_f_neq(:) = t_f_neq(:) + (f(:) - f_eq(:))

        tTurbVisc = tTurbVisc + fluid%turbulence%dataOnLvl(sourceLevel) &
          &                                      %visc(sourceElem)

      end do  ! iSourceElem

      ! interpolate turbulent viscosity on target element
      tTurbVisc = tTurbVisc * inv_nSourceElems
      ! scale interpolated turbulent viscosity to target element
      fluid%turbulence%dataOnLvl(targetLevel)%visc(targetElem) &
        & = fluid%turbulence%fac_f2c*tTurbVisc

      ! get normalized kinematic viscosity on target (coarse) element
      cVisc = fluid%viscKine%dataOnLvl(targetLevel)%val(targetElem)

      ! relation between coarse and fine grid kinematic viscosity:
      ! v^s_f = 2 v^s_c
      ! calculate omega on source and target level
      fOmegaKine = mus_calcOmegaFromVisc(2_rk*cVisc + tTurbVisc)
      cOmegaKine = mus_calcOmegaFromVisc(cVisc + fluid%turbulence%fac_f2c*tTurbVisc)

      !evaluate scaling factor
      fac = getNonEqFac_intp_fine_to_coarse(cOmegaKine, fOmegaKine)

      ! interpolate all pdfs by average and scale accordingly
      do iDir = 1, QQ
        t_f_eq(iDir) = t_f_eq(iDir) * inv_nSourceElems
        t_f_neq(iDir) = t_f_neq(iDir) * inv_nSourceElems * fac
        f(iDir) = t_f_eq(iDir) + t_f_neq(iDir)
        ! Now write the resulting pdf in the current direction to the target
        tState( ( targetelem-1)* nscalars+ idir+( 1-1)* qq) = f(iDir)
      enddo


    enddo

  end subroutine fillMyGhostsFromFiner_avgLES_feq_fneq
! **************************************************************************** !


! **************************************************************************** !
  !> Fill coarse target ghost from fine source fluid by average interpolation.
  !! 1. Calculate Equilibrium and nonEquilibrium
  !! 2. Compute scaling factor
  !! 3. calculate target: Eq + Scale * nonEquilibrium
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine fillMyGhostsFromFiner_avg2D_feq_fneq( method, fieldProp,         &
    &      tLevelDesc, level, sState, sNeigh, snSize, sAuxField, tState,      &
    &      tNeigh, tnSize, layout, nTargets, targetList, physics, time,       &
    &      varSys, derVarPos                                                  )
    ! -------------------------------------------------------------------- !
    class(mus_interpolation_method_type), intent(inout) :: method

    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), target, intent(in) :: fieldProp(:)

    !> level descriptor on target level
    type( tem_levelDesc_type ), intent(in) :: tLevelDesc

    !> my refinement level
    integer, intent(in) :: level

    !> State vector of SOURCE FLUID elements
    real(kind=rk), intent(in) :: sState(:)
    integer, intent(in) :: sNeigh(:)
    integer, intent(in) :: snSize

    !> State vector of TARGET GHOST elements
    real(kind=rk), intent(inout) :: tState(:)
    integer, intent(in) :: tNeigh(:)
    integer, intent(in) :: tnSize

    !> AuxField variable to read rho and vel from source elements
    real(kind=rk), intent(inout) :: sAuxField(:)

    !> the layout used
    type( mus_scheme_layout_type ), intent(in) :: layout

    !> List of target elements ( their position in depSource list )
    integer, intent(in) :: nTargets
    integer, intent(in) :: targetList(nTargets)

    !> physics type to convert lattice to physics SI unit and vice versa
    !! @todo: This can be replaced by scale factor by level
    type( mus_physics_type ), intent(in) :: physics

    !> time required to compute viscosity on target element barycenter
    type(tem_time_type), intent(in) :: time

    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys

    !> position of all derive variable in varSys for all fields
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iDir           ! current direction (discrete velocity) for loop
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    real(kind=rk) :: f( layout%fStencil%QQ ), f_eq( layout%fStencil%QQ )
    integer :: QQ
    real(kind=rk) :: inv_nSourceElems
    ! source elements' pdf
    real(kind=rk) :: t_f_eq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: t_f_neq( layout%fStencil%QQ )  ! temp pdf calculation
    type(mus_fluid_type), pointer :: fluid
    real(kind=rk) :: cVisc, cOmegaKine, fOmegaKine, rho, vel(3), fac
    integer :: elemOff, dens_pos, vel_pos(3), nScalars
    ! --------------------------------------------------------------------------
    ! for initialization do something else


    vel = 0._rk
    fluid => fieldProp(1)%fluid
    nScalars = varSys%nScalars
    QQ = layout%fStencil%QQ

    sourceLevel = level + 1
    targetLevel = level

    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)


    ! Treat all coarse target elements
    do indElem = 1, nTargets

      iElem = targetList( indElem )

      ! Read the target element treeId
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromFiner)

      ! Find out how many fine source elements we have for interpolation.
      nSourceELems = tLevelDesc%depFromFiner( iElem )%elem%nVals
      inv_nSourceElems = 1.0_rk / dble(nSourceElems)

      ! Now loop over all fine source elements for this target:
      t_f_eq = 0._rk
      t_f_neq = 0._rk
      do iSourceElem = 1, nSourceElems

        ! Get the source element position
        sourceElem = tLevelDesc%depFromFiner( iElem ) &
          &                    %elem%val( iSourceElem )

        do iDir = 1, QQ
        ! this is post collision
          f(iDir) = sState( ( sourceelem-1)* nscalars+ idir+( 1-1)* qq)
        enddo

        ! element offset for auxField array
        elemOff = (sourceElem-1)*varSys%nAuxScalars

        ! local density
        rho = sAuxField(elemOff + dens_pos)
        ! local x-, y-velocity
        vel(1) = sAuxField(elemOff + vel_pos(1))
        vel(2) = sAuxField(elemOff + vel_pos(2))


      f_eq(:) = layout%quantities%pdfEq_ptr( rho = rho,               &
        &                                    vel = vel,               &
        &                                    QQ = QQ                  )

      t_f_eq(:) = t_f_eq(:) + f_eq(:)
      t_f_neq(:) = t_f_neq(:) + (f(:) - f_eq(:))

      end do ! iSourceElem

      ! get normalized kinematic viscosity on target (coarse) element
      cVisc = fluid%viscKine%dataOnLvl(targetLevel)%val(targetElem)

      ! relation between coarse and fine grid kinematic viscosity:
      ! v^s_f = 2 v^s_c
      ! calculate omega on source and target level
      fOmegaKine = mus_calcOmegaFromVisc(2_rk*cVisc)
      cOmegaKine = mus_calcOmegaFromVisc(cVisc)

      !evaluate scaling factor
      fac = getNonEqFac_intp_fine_to_coarse( cOmegaKine, fOmegaKine)

      ! interpolate all pdfs by average and scale accordingly
      do iDir = 1, QQ
        t_f_eq(iDir) = t_f_eq(iDir) * inv_nSourceElems
        t_f_neq(iDir) = t_f_neq(iDir) * inv_nSourceElems * fac
        f(iDir) = t_f_eq(iDir) + t_f_neq(iDir)
        ! Now write the resulting pdf in the current direction to the target
        tState( ( targetelem-1)* nscalars+ idir+( 1-1)* qq) = f(iDir)
      enddo


    end do

  end subroutine fillMyGhostsFromFiner_avg2D_feq_fneq
! **************************************************************************** !


! **************************************************************************** !
  !> Interpolate auxiliary field from coarse source to fine target
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine_arbitraryVal]] in intp/[[mus_interpolate_header_module]].f90
  !! in order to be callable via
  !! [[mus_interpolation_method_type:do_intpArbiVal]] function pointer.
  subroutine fillArbiFinerGhostsFromMe_weighAvg( method, tLevelDesc, level,    &
    &                                            stencil, sVal, tVal, nTargets,&
    &                                            targetList, nScalars          )
    ! ------------------------------------------------------------------ !
    class(mus_interpolation_method_type), intent(inout) :: method

    !> my refinement level
    integer, intent(in) :: level

    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil

    !> State vector of SOURCE FLUID elements
    real(kind=rk), intent(in) :: sVal(:)

    !> State vector of TARGET GHOST elements
    real(kind=rk), intent(inout) :: tVal(:)

    !> level descriptor on target level
    type( tem_levelDesc_type ), intent(in) :: tLevelDesc

    !> List of target elements ( their position in depSource list )
    integer, intent(in) :: nTargets
    integer, intent(in) :: targetList(nTargets)

    !> number of scalars to interpolate
    integer, intent(in) :: nScalars
    ! ------------------------------------------------------------------ !
    ! --------------------------------------------------------------------------
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    real(kind=rk) :: weight( stencil%QQ )
    real(kind=rk) :: sArbi(nScalars, stencil%QQ)  ! temp source ArbiField
    integer :: iVar
    ! --------------------------------------------------------------------------
    sourceLevel = level
    targetLevel = level + 1

    ! Treat all coarse target elements
    do indElem = 1, nTargets

      iElem = targetList( indElem )

      ! Read the target element
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromCoarser)

      ! Get how many fine source elements we have for interpolation.
      nSourceElems = tLevelDesc%depFromCoarser( iElem )%elem%nVals

      ! Read the pre-calculated weights (like in Average intp.)
      weight(1:nSourceElems) = tLevelDesc%depFromCoarser( iElem ) &
        &                                %weight(1:nSourceElems)

      ! Now loop over all fine source elements for this target:
      do iSourceElem = 1, nSourceElems

        ! Get the source element
        sourceElem = tLevelDesc%depFromCoarser( iElem ) &
          &                    %elem%val( iSourceElem )

        ! Get souce auxilary variables
        sArbi(:, iSourceElem) = sVal( (sourceElem-1)*nScalars+1 &
          &                           : sourceElem*nScalars       )

      end do  ! iSourceElem

      ! interpolate all auxiliary variables by average
      do iVar = 1, nScalars
        tVal((targetElem-1)*nScalars+iVar)                    &
          & = sum( weight(1:nSourceElems)*sArbi(iVar,1:nSourceElems) )
      end do

    enddo

  end subroutine fillArbiFinerGhostsFromMe_weighAvg
! **************************************************************************** !

! **************************************************************************** !
  !> [Linear interpolation](../page/features/intp_methods.html) of ghostFromFiner
  !! 1. Calculate Equilibrium and nonEquilibrium
  !! 2. Compute scaling factor
  !! 3. calculate target: Eq + Scale * nonEquilibrium
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine fillFinerGhostsFromMe_weighAvg_feq_fneq( method, fieldProp,      &
    &      tLevelDesc, level, sState, sNeigh, snSize, sAuxField, tState,      &
    &      tNeigh, tnSize, layout, nTargets, targetList, physics, time,       &
    &      varSys, derVarPos                                                  )
    ! -------------------------------------------------------------------- !
    class(mus_interpolation_method_type), intent(inout) :: method

    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), target, intent(in) :: fieldProp(:)

    !> level descriptor on target level
    type( tem_levelDesc_type ), intent(in) :: tLevelDesc

    !> my refinement level
    integer, intent(in) :: level

    !> State vector of SOURCE FLUID elements
    real(kind=rk), intent(in) :: sState(:)
    integer, intent(in) :: sNeigh(:)
    integer, intent(in) :: snSize

    !> State vector of TARGET GHOST elements
    real(kind=rk), intent(inout) :: tState(:)
    integer, intent(in) :: tNeigh(:)
    integer, intent(in) :: tnSize

    !> AuxField variable to read rho and vel from source elements
    real(kind=rk), intent(inout) :: sAuxField(:)

    !> the layout used
    type( mus_scheme_layout_type ), intent(in) :: layout

    !> List of target elements ( their position in depSource list )
    integer, intent(in) :: nTargets
    integer, intent(in) :: targetList(nTargets)

    !> physics type to convert lattice to physics SI unit and vice versa
    !! @todo: This can be replaced by scale factor by level
    type( mus_physics_type ), intent(in) :: physics

    !> time required to compute viscosity on target element barycenter
    type(tem_time_type), intent(in) :: time

    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys

    !> position of all derive variable in varSys for all fields
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iDir           ! current direction (discrete velocity) for loop
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    real(kind=rk) :: f_neq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: f_eq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: weight(layout%fStencil%QQ)
    real(kind=rk) :: t_f_eq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: t_f_neq( layout%fStencil%QQ )  ! temp pdf calculation
    ! pdf to reconstruct from
    real(kind=rk) :: f( layout%fStencil%QQ )
    integer :: QQ
    type(mus_fluid_type), pointer :: fluid
    real(kind=rk) :: fVisc, cOmegaKine, fOmegaKine, rho, vel(3), fac
    integer :: elemOff, dens_pos, vel_pos(3), nScalars
    ! --------------------------------------------------------------------------


    fluid => fieldProp(1)%fluid
    nScalars = varSys%nScalars
    QQ = layout%fStencil%QQ

    sourceLevel = level
    targetLevel = level + 1

    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)


    ! Treat all fine target elements:
    do indElem = 1, nTargets

      iElem = targetList( indElem )

      ! Read the target element treeId
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromCoarser)

      ! Find out how many fine source elements we have for interpolation.
      nSourceElems = tLevelDesc%depFromCoarser( iElem )%elem%nVals

      ! Read the pre-calculated weights (like in Average intp.)
      weight(1:nSourceElems) = tLevelDesc%depFromCoarser( iElem ) &
        &                                %weight(1:nSourceElems)

      ! Now loop over all fine source elements for this target:
      do iSourceElem = 1, nSourceElems

        ! Get the source element position
        sourceElem = tLevelDesc%depFromCoarser( iElem ) &
          &                    %elem%val( iSourceElem )

        do iDir = 1, QQ
        ! this is post collision
          f(iDir) = sState( ( sourceelem-1)* nscalars+ idir+( 1-1)* qq)
        enddo

        ! element offset for auxField array
        elemOff = (sourceElem-1)*varSys%nAuxScalars

        ! local density
        rho = sAuxField(elemOff + dens_pos)
        ! local x-, y-, z-velocity
        vel(1) = sAuxField(elemOff + vel_pos(1))
        vel(2) = sAuxField(elemOff + vel_pos(2))
        vel(3) = sAuxField(elemOff + vel_pos(3))


        f_eq(:, iSourceElem) = layout%quantities%pdfEq_ptr( rho = rho,  &
          &                                    vel = vel,               &
          &                                    QQ = QQ                  )

        f_neq(:, iSourceElem) = f(:) - f_eq(:, iSourceElem)

      end do ! iSourceElem

      ! interpolate all pdfs by weighted average
      do iDir = 1,QQ
        t_f_eq(iDir) = sum(weight(1:nSourceElems)*f_eq(iDir,1:nSourceElems))
        t_f_neq(iDir) = sum(weight(1:nSourceElems)*f_neq(iDir,1:nSourceElems))
      end do

      ! get normalized kinematic viscosity on target (fine) element
      fVisc = fluid%viscKine%dataOnLvl(targetLevel)%val(targetElem)

      ! relation between coarse and fine grid kinematic viscosity:
      ! v^s_c = 0.5 v^s_f
      ! calculate omega on source and target level
      fOmegaKine = mus_calcOmegaFromVisc(fVisc)
      cOmegaKine = mus_calcOmegaFromVisc(0.5_rk*fVisc)

      !evaluate scaling factor
      fac = getNonEqFac_intp_coarse_to_fine( cOmegaKine, fOmegaKine )

      ! Rescale the non eq pdfs
      do iDir=1, QQ
        t_f_neq(iDir) = t_f_neq(iDir) * fac
        f(iDir) = t_f_neq(iDir) + t_f_eq(iDir)
        ! Now write the resulting pdf in the current direction to the target
        tState( ( targetelem-1)* nscalars+ idir+( 1-1)* qq) = f(iDir)
      enddo


    end do ! indElem

  end subroutine fillFinerGhostsFromMe_weighAvg_feq_fneq
! **************************************************************************** !

! **************************************************************************** !
  !> [Linear interpolation](../page/features/intp_methods.html) of ghostFromFiner
  !! 1. Calculate Equilibrium and nonEquilibrium
  !! 2. Compute scaling factor
  !! 3. calculate target: Eq + Scale * nonEquilibrium
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine fillFinerGhostsFromMe_weighAvgLES_feq_fneq( method, fieldProp,   &
    &      tLevelDesc, level, sState, sNeigh, snSize, sAuxField, tState,      &
    &      tNeigh, tnSize, layout, nTargets, targetList, physics, time,       &
    &      varSys, derVarPos                                                  )
    ! -------------------------------------------------------------------- !
    class(mus_interpolation_method_type), intent(inout) :: method

    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), target, intent(in) :: fieldProp(:)

    !> level descriptor on target level
    type( tem_levelDesc_type ), intent(in) :: tLevelDesc

    !> my refinement level
    integer, intent(in) :: level

    !> State vector of SOURCE FLUID elements
    real(kind=rk), intent(in) :: sState(:)
    integer, intent(in) :: sNeigh(:)
    integer, intent(in) :: snSize

    !> State vector of TARGET GHOST elements
    real(kind=rk), intent(inout) :: tState(:)
    integer, intent(in) :: tNeigh(:)
    integer, intent(in) :: tnSize

    !> AuxField variable to read rho and vel from source elements
    real(kind=rk), intent(inout) :: sAuxField(:)

    !> the layout used
    type( mus_scheme_layout_type ), intent(in) :: layout

    !> List of target elements ( their position in depSource list )
    integer, intent(in) :: nTargets
    integer, intent(in) :: targetList(nTargets)

    !> physics type to convert lattice to physics SI unit and vice versa
    !! @todo: This can be replaced by scale factor by level
    type( mus_physics_type ), intent(in) :: physics

    !> time required to compute viscosity on target element barycenter
    type(tem_time_type), intent(in) :: time

    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys

    !> position of all derive variable in varSys for all fields
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iDir           ! current direction (discrete velocity) for loop
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    real(kind=rk) :: f_neq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: f_eq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: weight(layout%fStencil%QQ), sTurbVisc(layout%fStencil%QQ)
    real(kind=rk) :: t_f_eq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: t_f_neq( layout%fStencil%QQ )  ! temp pdf calculation
    ! pdf to reconstruct from
    real(kind=rk) :: f( layout%fStencil%QQ ), tTurbVisc
    integer :: QQ
    type(mus_fluid_type), pointer :: fluid
    real(kind=rk) :: fVisc, cOmegaKine, fOmegaKine, rho, vel(3), fac
    integer :: elemOff, dens_pos, vel_pos(3), nScalars
    ! --------------------------------------------------------------------------


    fluid => fieldProp(1)%fluid
    nScalars = varSys%nScalars
    QQ = layout%fStencil%QQ

    sourceLevel = level
    targetLevel = level + 1

    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)


    ! Treat all fine target elements:
    do indElem = 1, nTargets

      iElem = targetList( indElem )

      ! Read the target element treeId
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromCoarser)

      ! Find out how many fine source elements we have for interpolation.
      nSourceElems = tLevelDesc%depFromCoarser( iElem )%elem%nVals

      ! Read the pre-calculated weights (like in Average intp.)
      weight(1:nSourceElems) = tLevelDesc%depFromCoarser( iElem ) &
        &                                %weight(1:nSourceElems)

      ! Now loop over all fine source elements for this target:
      do iSourceElem = 1, nSourceElems

        ! Get the source element position
        sourceElem = tLevelDesc%depFromCoarser( iElem ) &
          &                    %elem%val( iSourceElem )

        do iDir = 1, QQ
        ! this is post collision
          f(iDir) = sState( ( sourceelem-1)* nscalars+ idir+( 1-1)* qq)
        enddo

        ! element offset for auxField array
        elemOff = (sourceElem-1)*varSys%nAuxScalars

        ! local density
        rho = sAuxField(elemOff + dens_pos)
        ! local x-, y-, z-velocity
        vel(1) = sAuxField(elemOff + vel_pos(1))
        vel(2) = sAuxField(elemOff + vel_pos(2))
        vel(3) = sAuxField(elemOff + vel_pos(3))


        f_eq(:, iSourceElem) = layout%quantities%pdfEq_ptr( rho = rho,  &
          &                                    vel = vel,               &
          &                                    QQ = QQ                  )

        f_neq(:, iSourceElem) = f(:) - f_eq(:, iSourceElem)

        ! get turbulent viscosity
        sTurbVisc(iSourceElem) = fluid%turbulence%dataOnLvl(sourceLevel) &
          &                                       %visc(sourceElem)

      end do ! iSourceElem

      ! interpolate all pdfs by weighted average
      tTurbVisc = sum(weight(1:nSourceElems)*sTurbVisc(1:nSourceElems))
      do iDir = 1,QQ
        t_f_eq(iDir) = sum(weight(1:nSourceElems)*f_eq(iDir,1:nSourceElems))
        t_f_neq(iDir) = sum(weight(1:nSourceElems)*f_neq(iDir,1:nSourceElems))
      end do

      ! scale interpolated turbulent viscosity to target element
      fluid%turbulence%dataOnLvl(targetLevel)%visc(targetElem) &
        & = fluid%turbulence%fac_c2f*tTurbVisc

      ! get normalized kinematic viscosity on target (fine) element
      fVisc = fluid%viscKine%dataOnLvl(targetLevel)%val(targetElem)

      ! relation between coarse and fine grid kinematic viscosity:
      ! v^s_c = 0.5 v^s_f
      ! calculate omega on source and target level
      fOmegaKine = mus_calcOmegaFromVisc(fVisc + fluid%turbulence%fac_c2f*tTurbVisc)
      cOmegaKine = mus_calcOmegaFromVisc(0.5_rk*fVisc + tTurbVisc)

      !evaluate scaling factor
      fac = getNonEqFac_intp_coarse_to_fine( cOmegaKine, fOmegaKine )

      ! Rescale the non eq pdfs
      do iDir=1, QQ
        t_f_neq(iDir) = t_f_neq(iDir) * fac
        f(iDir) = t_f_neq(iDir) + t_f_eq(iDir)
        ! Now write the resulting pdf in the current direction to the target
        tState( ( targetelem-1)* nscalars+ idir+( 1-1)* qq) = f(iDir)
      enddo


    end do ! indElem

  end subroutine fillFinerGhostsFromMe_weighAvgLES_feq_fneq
! **************************************************************************** !

! **************************************************************************** !
  !> [Linear interpolation](../page/features/intp_methods.html) of ghostFromFiner
  !! 1. Calculate Equilibrium and nonEquilibrium
  !! 2. Compute scaling factor
  !! 3. calculate target: Eq + Scale * nonEquilibrium
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine fillFinerGhostsFromMe_weighAvg2D_feq_fneq( method, fieldProp,    &
    &      tLevelDesc, level, sState, sNeigh, snSize, sAuxField, tState,      &
    &      tNeigh, tnSize, layout, nTargets, targetList, physics, time,       &
    &      varSys, derVarPos                                                  )
    ! -------------------------------------------------------------------- !
    class(mus_interpolation_method_type), intent(inout) :: method

    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), target, intent(in) :: fieldProp(:)

    !> level descriptor on target level
    type( tem_levelDesc_type ), intent(in) :: tLevelDesc

    !> my refinement level
    integer, intent(in) :: level

    !> State vector of SOURCE FLUID elements
    real(kind=rk), intent(in) :: sState(:)
    integer, intent(in) :: sNeigh(:)
    integer, intent(in) :: snSize

    !> State vector of TARGET GHOST elements
    real(kind=rk), intent(inout) :: tState(:)
    integer, intent(in) :: tNeigh(:)
    integer, intent(in) :: tnSize

    !> AuxField variable to read rho and vel from source elements
    real(kind=rk), intent(inout) :: sAuxField(:)

    !> the layout used
    type( mus_scheme_layout_type ), intent(in) :: layout

    !> List of target elements ( their position in depSource list )
    integer, intent(in) :: nTargets
    integer, intent(in) :: targetList(nTargets)

    !> physics type to convert lattice to physics SI unit and vice versa
    !! @todo: This can be replaced by scale factor by level
    type( mus_physics_type ), intent(in) :: physics

    !> time required to compute viscosity on target element barycenter
    type(tem_time_type), intent(in) :: time

    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys

    !> position of all derive variable in varSys for all fields
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iDir           ! current direction (discrete velocity) for loop
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    real(kind=rk) :: f_neq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: f_eq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: weight(layout%fStencil%QQ)
    real(kind=rk) :: t_f_eq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: t_f_neq( layout%fStencil%QQ )  ! temp pdf calculation
    ! pdf to reconstruct from
    real(kind=rk) :: f( layout%fStencil%QQ )
    integer :: QQ
    type(mus_fluid_type), pointer :: fluid
    real(kind=rk) :: fVisc, cOmegaKine, fOmegaKine, rho, vel(3), fac
    integer :: elemOff, dens_pos, vel_pos(3), nScalars
    ! --------------------------------------------------------------------------


    vel = 0._rk
    fluid => fieldProp(1)%fluid
    nScalars = varSys%nScalars
    QQ = layout%fStencil%QQ

    sourceLevel = level
    targetLevel = level + 1

    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)


    ! Treat all fine target elements:
    do indElem = 1, nTargets

      iElem = targetList( indElem )

      ! Read the target element treeId
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromCoarser)

      ! Find out how many fine source elements we have for interpolation.
      nSourceElems = tLevelDesc%depFromCoarser( iElem )%elem%nVals

      ! Read the pre-calculated weights (like in Average intp.)
      weight(1:nSourceElems) = tLevelDesc%depFromCoarser( iElem ) &
        &                                %weight(1:nSourceElems)

      ! Now loop over all fine source elements for this target:
      do iSourceElem = 1, nSourceElems

        ! Get the source element position
        sourceElem = tLevelDesc%depFromCoarser( iElem ) &
          &                    %elem%val( iSourceElem )

        do iDir = 1, QQ
        ! this is post collision
          f(iDir) = sState( ( sourceelem-1)* nscalars+ idir+( 1-1)* qq)
        enddo

        ! element offset for auxField array
        elemOff = (sourceElem-1)*varSys%nAuxScalars

        ! local density
        rho = sAuxField(elemOff + dens_pos)
        ! local x-, y-, z-velocity
        vel(1) = sAuxField(elemOff + vel_pos(1))
        vel(2) = sAuxField(elemOff + vel_pos(2))


        f_eq(:, iSourceElem) = layout%quantities%pdfEq_ptr( rho = rho,  &
          &                                    vel = vel,               &
          &                                    QQ = QQ                  )

        f_neq(:, iSourceElem) = f(:) - f_eq(:, iSourceElem)

      end do ! iSourceElem

      ! interpolate all pdfs by weighted average
      do iDir = 1,QQ
        t_f_eq(iDir) = sum(weight(1:nSourceElems)*f_eq(iDir,1:nSourceElems))
        t_f_neq(iDir) = sum(weight(1:nSourceElems)*f_neq(iDir,1:nSourceElems))
      end do

      ! get normalized kinematic viscosity on target (fine) element
      fVisc = fluid%viscKine%dataOnLvl(targetLevel)%val(targetElem)

      ! relation between coarse and fine grid kinematic viscosity:
      ! v^s_c = 0.5 v^s_f
      ! calculate omega on source and target level
      fOmegaKine = mus_calcOmegaFromVisc(fVisc)
      cOmegaKine = mus_calcOmegaFromVisc(0.5_rk*fVisc)

      !evaluate scaling factor
      fac = getNonEqFac_intp_coarse_to_fine( cOmegaKine, fOmegaKine )

      ! Rescale the non eq pdfs
      do iDir=1, QQ
        t_f_neq(iDir) = t_f_neq(iDir) * fac
        f(iDir) = t_f_neq(iDir) + t_f_eq(iDir)
        ! Now write the resulting pdf in the current direction to the target
        tState( ( targetelem-1)* nscalars+ idir+( 1-1)* qq) = f(iDir)
      enddo


    end do ! indElem

  end subroutine fillFinerGhostsFromMe_weighAvg2D_feq_fneq
! **************************************************************************** !

end module mus_interpolate_average_module
! **************************************************************************** !
