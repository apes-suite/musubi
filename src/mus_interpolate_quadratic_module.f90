! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2011-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2016, 2018-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016, 2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019 Jana Gericke <jana.gericke@uni-siegen.de>
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
! ****************************************************************************** !
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
!> summary: Interpolation of flow quantities between different grid levels
!!
!! Ghost elements are employed at grid level interfaces to provide valid
!! pdf values to the neighboring fluid elements. This way, the solvers can
!! act on elements of the same size only, treating the levels successively.
!! Target elements are the ghost elements, which have to be filled with
!! valid values.
!! Source elements are the fluid elements from other levels, from where to
!! take the input values for the interpolation.
!! The target ghost elements on the target level have corresponding source
!! fluid elements on the source level.
!!
!! [[tem_topology_module]] For a detailed description of the grid
!!
module mus_interpolate_quadratic_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,              only: rk
  use tem_aux_module,          only: tem_abort
  use tem_param_module,        only: cs2inv, cs2, div1_3, div1_6, div1_9, &
    &                                div2_9, div5_9, div2_3, div1_4, &
    &                                div1_21, div4_21, div5_21, div1_7, &
    &                                div3_7, div1_42, div5_42, div1_2, &
    &                                div1_36, div3_4h, c_x, c_y, c_z, &
    &                                rho0, rho0Inv
  use tem_element_module,      only: eT_GhostFromCoarser
  use tem_construction_module, only: depSource_type, tem_levelDesc_type
  use tem_debug_module,        only: dbgUnit
  use tem_logging_module,      only: logUnit
  use tem_stencil_module,      only: tem_stencilHeader_type
  use tem_matrix_module,       only: tem_matrix_type, tem_matrix_dump
  use tem_varSys_module,       only: tem_varSys_type
  use tem_time_module,         only: tem_time_type

  ! include musubi modules
  use mus_pdf_module,                only: pdf_data_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_interpolate_header_module, only: mus_interpolation_method_type
  use mus_physics_module,            only: mus_physics_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_fluid_module,              only: mus_fluid_type
  use mus_relaxationParam_module,    only: mus_calcOmegaFromVisc
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_derivedQuantities_module2,  only: getNonEqFac_intp_fine_to_coarse, &
    &                                       getNonEqFac_intp_coarse_to_fine

  implicit none

  private

  public :: fillArbiFinerGhostsFromMe_quad
  public :: fillArbiFinerGhostsFromMe_quad2D
  public :: fillFinerGhostsFromMe_quad_feq_fneq
  public :: fillFinerGhostsFromMe_quadLES_feq_fneq
  public :: fillFinerGhostsFromMe_quad2D_feq_fneq

contains

! **************************************************************************** !
  !> Interpolate auxiliary field from coarse source to fine target
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine_arbitraryVal]] in intp/[[mus_interpolate_header_module]].f90
  !! in order to be callable via
  !! [[mus_interpolation_method_type:do_intpArbiVal]] function pointer.
  subroutine fillArbiFinerGhostsFromMe_quad( method, tLevelDesc, level,     &
    &                                        stencil, sVal, tVal, nTargets, &
    &                                        targetList, nScalars           )
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
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    integer :: posInIntpMatLSF
    real(kind=rk) :: tArbi(nScalars)
    real(kind=rk) :: sArbi(nScalars, stencil%QQ)  ! temp source ArbiField
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
      posInIntpMatLSF = tLevelDesc%depFromCoarser( iElem )%posInIntpMatLSF

      ! Now loop over all fine source elements for this target:
      do iSourceElem = 1, nSourceElems

        ! Get the source element
        sourceElem = tLevelDesc%depFromCoarser( iElem ) &
          &                    %elem%val( iSourceElem )

        ! Get souce auxilary variables
        sArbi(:, iSourceElem) = sVal( (sourceElem-1)*nScalars+1 &
          &                           : sourceElem*nScalars       )

      end do  ! iSourceElem

      ! interpolate all auxiliary variables by quadratic interpolation
      tArbi(1:nScalars) = mus_interpolate_quad3D_leastSq(         &
        &   srcMom      = sArbi(1:nScalars, 1:nSourceElems),      &
        &   targetCoord = tLevelDesc%depFromCoarser( iElem )%coord, &
        &   LSFmat      = method%intpMat_forLSF%matArray            &
        &                       %val(posInIntpMatLSF),              &
        &   nSources    = nSourceElems,                             &
        &   nVals       = nScalars                               )

      ! write interpolated value
      tVal((targetElem-1)*nScalars+1 : targetElem*nScalars) &
          & = tArbi
    enddo

  end subroutine fillArbiFinerGhostsFromMe_quad
! **************************************************************************** !

! **************************************************************************** !
  !> Interpolate auxiliary field from coarse source to fine target
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine_arbitraryVal]] in intp/[[mus_interpolate_header_module]].f90
  !! in order to be callable via
  !! [[mus_interpolation_method_type:do_intpArbiVal]] function pointer.
  subroutine fillArbiFinerGhostsFromMe_quad2D( method, tLevelDesc, level,     &
    &                                          stencil, sVal, tVal, nTargets, &
    &                                          targetList, nScalars           )
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
    integer :: sourceLevel    ! level of source elements
    integer :: sourceElem     ! treeId of current source element
    integer :: targetLevel    ! level of target elements
    integer :: targetElem     ! treeId of current source element
    integer :: iElem          ! current target element (for outer loop)
    integer :: indElem        ! element counter for indirection list
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    integer :: posInIntpMatLSF
    real(kind=rk) :: tArbi(nScalars)
    real(kind=rk) :: sArbi(nScalars, stencil%QQ)  ! temp source ArbiField
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

      posInIntpMatLSF = tLevelDesc%depFromCoarser( iElem )%posInIntpMatLSF
      ! Now loop over all fine source elements for this target:
      do iSourceElem = 1, nSourceElems

        ! Get the source element
        sourceElem = tLevelDesc%depFromCoarser( iElem ) &
          &                    %elem%val( iSourceElem )

        ! Get souce auxilary variables
        sArbi(:, iSourceElem) = sVal( (sourceElem-1)*nScalars+1 &
          &                           : sourceElem*nScalars       )

      end do  ! iSourceElem

      ! interpolate all auxiliary variables by quadratic interpolation
      tArbi(1:nScalars) = mus_interpolate_quad2D_leastSq(         &
        &   srcMom      = sArbi(1:nScalars, 1:nSourceElems),      &
        &   targetCoord = tLevelDesc%depFromCoarser( iElem )%coord, &
        &   LSFmat      = method%intpMat_forLSF%matArray            &
        &                       %val(posInIntpMatLSF),              &
        &   nSources    = nSourceElems,                             &
        &   nVals       = nScalars                               )

      ! write interpolated value
      tVal((targetElem-1)*nScalars+1 : targetElem*nScalars) &
          & = tArbi
    enddo

  end subroutine fillArbiFinerGhostsFromMe_quad2D
! **************************************************************************** !

! ************************************************************************** !
  !> Fill fine ghost from coarse fluid by quadratic interpolation for D2Q9
  !! stencil.
  !! 1. Compute moments for all source elements, save in momBuf
  !! 2. For each target, interpolate moments (den, vel, tau)
  !!    (10 moments for 3D and 6 moments for 2D)
  !! 3. calculate fEq and use it to calculate high order moments
  !! 4. convert moments to PDF
  !! This routine is used by acoustic quadratic interpolation.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine fillFinerGhostsFromMe_quad_feq_fneq( method, fieldProp,     &
    &      tLevelDesc, level, sState, sNeigh, snSize, sAuxField, tState, &
    &      tNeigh, tnSize, layout, nTargets, targetList, physics, time,  &
    &      varSys, derVarPos                                             )
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
    integer :: iElem, indElem, iDir
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    ! pdf to reconstruct from
    real(kind=rk) :: f( layout%fStencil%QQ )
    ! pdf to interpolate
    real(kind=rk) :: f_neq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: f_eq( layout%fStencil%QQ, layout%fStencil%QQ )
    ! source elements' pdf
    real(kind=rk) :: t_f_eq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: t_f_neq( layout%fStencil%QQ )  ! temp pdf calculation
    ! shear stress scaling factor
    real(kind=rk) :: fac, rho, vel(3)
    integer :: QQ
    integer :: posInIntpMatLSF
    type(mus_fluid_type), pointer :: fluid
    integer :: elemOff, dens_pos, vel_pos(3), nScalars
    real(kind=rk) :: fVisc, fOmegaKine, cOmegaKine
    ! ---------------------------------------------------------------------------


    vel = 0._rk
    fluid => fieldProp(1)%fluid
    QQ = layout%fStencil%QQ
    nScalars = varSys%nScalars

    sourceLevel = level
    targetLevel = level + 1

    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)


    ! Treat all fine target elements:
    do indElem = 1, nTargets

      iElem = targetList( indElem )
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromCoarser)
      nSourceElems = tLevelDesc%depFromCoarser( iElem )%elem%nVals
      posInIntpMatLSF = tLevelDesc%depFromCoarser( iElem )%posInIntpMatLSF

      ! First calculate all the required moments for all the source elements
      do iSourceElem = 1, nSourceElems

        ! Get the source element's position in the state vector
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


        f_eq(:, iSourceElem) = layout%quantities%pdfEq_ptr( rho = rho, &
          &                                   vel = vel,               &
          &                                   QQ = QQ                  )

        f_neq(:, iSourceElem) = f(:) - f_eq(:, iSourceElem)

      enddo

      t_f_eq = mus_interpolate_quad3D_leastSq(                      &
        &   srcMom      = f_eq(1:QQ, 1:nSourceElems),               &
        &   targetCoord = tLevelDesc%depFromCoarser( iElem )%coord, &
        &   LSFmat      = method%intpMat_forLSF%matArray            &
        &                       %val(posInIntpMatLSF),              &
        &   nSources    = nSourceElems,                             &
        &   nVals       = QQ )

      t_f_neq = mus_interpolate_quad3D_leastSq(                     &
        &   srcMom      = f_neq(1:QQ, 1:nSourceElems),              &
        &   targetCoord = tLevelDesc%depFromCoarser( iElem )%coord, &
        &   LSFmat      = method%intpMat_forLSF%matArray            &
        &                       %val(posInIntpMatLSF),              &
        &   nSources    = nSourceElems,                             &
        &   nVals       = QQ )

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

  end subroutine fillFinerGhostsFromMe_quad_feq_fneq
! ****************************************************************************** !

! ************************************************************************** !
  !> Fill fine ghost from coarse fluid by quadratic interpolation for D2Q9
  !! stencil.
  !! 1. Compute moments for all source elements, save in momBuf
  !! 2. For each target, interpolate moments (den, vel, tau)
  !!    (10 moments for 3D and 6 moments for 2D)
  !! 3. calculate fEq and use it to calculate high order moments
  !! 4. convert moments to PDF
  !! This routine is used by acoustic quadratic interpolation.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine fillFinerGhostsFromMe_quadLES_feq_fneq( method, fieldProp,     &
    &      tLevelDesc, level, sState, sNeigh, snSize, sAuxField, tState,    &
    &      tNeigh, tnSize, layout, nTargets, targetList, physics, time,     &
    &      varSys, derVarPos                                                )
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
    integer :: iElem, indElem, iDir
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    ! pdf to reconstruct from
    real(kind=rk) :: f( layout%fStencil%QQ )
    ! pdf to interpolate
    real(kind=rk) :: f_neq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: f_eq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: sTurbVisc(layout%fStencil%QQ)
    ! source elements' pdf
    real(kind=rk) :: t_f_eq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: t_f_neq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: tTurbVisc(1)
    ! shear stress scaling factor
    real(kind=rk) :: fac, rho, vel(3)
    integer :: QQ
    integer :: posInIntpMatLSF
    type(mus_fluid_type), pointer :: fluid
    integer :: elemOff, dens_pos, vel_pos(3), nScalars
    real(kind=rk) :: fVisc, fOmegaKine, cOmegaKine
    ! ---------------------------------------------------------------------------


    vel = 0._rk
    fluid => fieldProp(1)%fluid
    QQ = layout%fStencil%QQ
    nScalars = varSys%nScalars

    sourceLevel = level
    targetLevel = level + 1

    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)


    ! Treat all fine target elements:
    do indElem = 1, nTargets

      iElem = targetList( indElem )
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromCoarser)
      nSourceElems = tLevelDesc%depFromCoarser( iElem )%elem%nVals
      posInIntpMatLSF = tLevelDesc%depFromCoarser( iElem )%posInIntpMatLSF

      ! First calculate all the required moments for all the source elements
      do iSourceElem = 1, nSourceElems

        ! Get the source element's position in the state vector
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


        f_eq(:, iSourceElem) = layout%quantities%pdfEq_ptr( rho = rho, &
          &                                   vel = vel,               &
          &                                   QQ = QQ                  )

        f_neq(:, iSourceElem) = f(:) - f_eq(:, iSourceElem)

        ! get turbulent viscosity
        sTurbVisc(iSourceElem) = fluid%turbulence%dataOnLvl(sourceLevel) &
          &                                       %visc(sourceElem)
      enddo

      t_f_eq = mus_interpolate_quad3D_leastSq(                      &
        &   srcMom      = f_eq(1:QQ, 1:nSourceElems),               &
        &   targetCoord = tLevelDesc%depFromCoarser( iElem )%coord, &
        &   LSFmat      = method%intpMat_forLSF%matArray            &
        &                       %val(posInIntpMatLSF),              &
        &   nSources    = nSourceElems,                             &
        &   nVals       = QQ )

      t_f_neq = mus_interpolate_quad3D_leastSq(                     &
        &   srcMom      = f_neq(1:QQ, 1:nSourceElems),              &
        &   targetCoord = tLevelDesc%depFromCoarser( iElem )%coord, &
        &   LSFmat      = method%intpMat_forLSF%matArray            &
        &                       %val(posInIntpMatLSF),              &
        &   nSources    = nSourceElems,                             &
        &   nVals       = QQ )

      !interpolate turbulent viscosity to target element
      tTurbVisc = mus_interpolate_quad3D_leastSq(                   &
        &   srcMom      = sTurbVisc(1:nSourceElems),                &
        &   targetCoord = tLevelDesc%depFromCoarser( iElem )%coord, &
        &   LSFmat      = method%intpMat_forLSF%matArray            &
        &                       %val(posInIntpMatLSF),              &
        &   nSources    = nSourceElems,                             &
        &   nVals       = 1                                         )

      ! scale interpolated turbulent viscosity to target element
      fluid%turbulence%dataOnLvl(targetLevel)%visc(targetElem) &
        & = fluid%turbulence%fac_c2f*tTurbVisc(1)

      ! get normalized kinematic viscosity on target (fine) element
      fVisc = fluid%viscKine%dataOnLvl(targetLevel)%val(targetElem)

      ! relation between coarse and fine grid kinematic viscosity:
      ! v^s_c = 0.5 v^s_f
      ! calculate omega on source and target level
      fOmegaKine = mus_calcOmegaFromVisc(fVisc + fluid%turbulence%fac_c2f*tTurbVisc(1))
      cOmegaKine = mus_calcOmegaFromVisc(0.5_rk*fVisc + tTurbVisc(1))

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

  end subroutine fillFinerGhostsFromMe_quadLES_feq_fneq
! ****************************************************************************** !

! ************************************************************************** !
  !> Fill fine ghost from coarse fluid by quadratic interpolation for D2Q9
  !! stencil.
  !! 1. Compute moments for all source elements, save in momBuf
  !! 2. For each target, interpolate moments (den, vel, tau)
  !!    (10 moments for 3D and 6 moments for 2D)
  !! 3. calculate fEq and use it to calculate high order moments
  !! 4. convert moments to PDF
  !! This routine is used by acoustic quadratic interpolation.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine fillFinerGhostsFromMe_quad2D_feq_fneq( method, fieldProp,     &
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
    integer :: iElem, indElem, iDir
    integer :: iSourceElem    ! current source element (for inner loop)
    integer :: nSourceElems   ! number of source elements for the current target
    ! pdf to reconstruct from
    real(kind=rk) :: f( layout%fStencil%QQ )
    ! pdf to interpolate
    real(kind=rk) :: f_neq( layout%fStencil%QQ, layout%fStencil%QQ )
    real(kind=rk) :: f_eq( layout%fStencil%QQ, layout%fStencil%QQ )
    ! source elements' pdf
    real(kind=rk) :: t_f_eq( layout%fStencil%QQ )  ! temp pdf calculation
    real(kind=rk) :: t_f_neq( layout%fStencil%QQ )  ! temp pdf calculation
    ! shear stress scaling factor
    real(kind=rk) :: fac, rho, vel(3)
    integer :: QQ
    integer :: posInIntpMatLSF
    type(mus_fluid_type), pointer :: fluid
    integer :: elemOff, dens_pos, vel_pos(3), nScalars
    real(kind=rk) :: fVisc, fOmegaKine, cOmegaKine
    ! ---------------------------------------------------------------------------


    vel = 0._rk
    fluid => fieldProp(1)%fluid
    QQ = layout%fStencil%QQ
    nScalars = varSys%nScalars

    sourceLevel = level
    targetLevel = level + 1

    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)


    ! Treat all fine target elements:
    do indElem = 1, nTargets

      iElem = targetList( indElem )
      targetElem = iElem + tLevelDesc%offset( 1, eT_ghostFromCoarser)
      nSourceElems = tLevelDesc%depFromCoarser( iElem )%elem%nVals
      posInIntpMatLSF = tLevelDesc%depFromCoarser( iElem )%posInIntpMatLSF

      ! First calculate all the required moments for all the source elements
      do iSourceElem = 1, nSourceElems

        ! Get the source element's position in the state vector
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
        ! local x-, y-velocity
        vel(1) = sAuxField(elemOff + vel_pos(1))
        vel(2) = sAuxField(elemOff + vel_pos(2))


        f_eq(:, iSourceElem) = layout%quantities%pdfEq_ptr( rho = rho, &
          &                                   vel = vel,               &
          &                                   QQ = QQ                  )

        f_neq(:, iSourceElem) = f(:) - f_eq(:, iSourceElem)

      enddo

      t_f_eq = mus_interpolate_quad2D_leastSq(                      &
        &   srcMom      = f_eq(1:QQ, 1:nSourceElems),               &
        &   targetCoord = tLevelDesc%depFromCoarser( iElem )%coord, &
        &   LSFmat      = method%intpMat_forLSF%matArray            &
        &                       %val(posInIntpMatLSF),              &
        &   nSources    = nSourceElems,                             &
        &   nVals       = QQ )

      t_f_neq = mus_interpolate_quad2D_leastSq(                     &
        &   srcMom      = f_neq(1:QQ, 1:nSourceElems),              &
        &   targetCoord = tLevelDesc%depFromCoarser( iElem )%coord, &
        &   LSFmat      = method%intpMat_forLSF%matArray            &
        &                       %val(posInIntpMatLSF),              &
        &   nSources    = nSourceElems,                             &
        &   nVals       = QQ )

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

  end subroutine fillFinerGhostsFromMe_quad2D_feq_fneq
! ****************************************************************************** !

! ****************************************************************************** !
  !> Biquadratic interpolation for a vector quantity phi
  !!
  ! real(kind=rk), dimension(6,9),parameter  :: intpMaxtrix =  &
  !           reshape((/ &
  !   &  div2_9,  -div1_6,    0.0_rk,    div1_6,   -div1_3,   0.0_rk, &
  !   &  div2_9,    0.0_rk,  -div1_6,   -div1_3,    div1_6,   0.0_rk, &
  !   &  div2_9,   div1_6,    0.0_rk,    div1_6,   -div1_3,   0.0_rk, &
  !   &  div2_9,    0.0_rk,   div1_6,   -div1_3,    div1_6,   0.0_rk, &
  !   & -div1_9,  -div1_6,  -div1_6,    div1_6,    div1_6,  0.25_rk, &
  !   & -div1_9,  -div1_6,   div1_6,    div1_6,    div1_6, -0.25_rk, &
  !   & -div1_9,   div1_6,  -div1_6,    div1_6,    div1_6, -0.25_rk, &
  !   & -div1_9,   div1_6,   div1_6,    div1_6,    div1_6,  0.25_rk, &
  !   &  div5_9,    0.0_rk,    0.0_rk,   -div1_3,   -div1_3,   0.0_rk  &
  !  /),(/6,9/))
  pure function mus_interpolate_quad2D_leastSq( srcMom, targetCoord, LSFmat, &
    &                                           nSources, nVals ) result(phi)
    ! ---------------------------------------------------------------------------
    !> number of quantities to interpolation
    integer, intent(in) :: nVals
    !> Number of source elements
    integer, intent(in) :: nSources
    !> source values of the momentum on the square corners
    real(kind=rk), intent(in) :: srcMom(nVals, nSources)
    !> interpolation location within the square
    real(kind=rk), intent(in) :: targetCoord(3)
    !> matrix for least square fit
    type(tem_matrix_type), intent(in) :: LSFmat
    !> interpolated value
    real(kind=rk) :: phi( nVals )
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: a(6) ! Coefficients
    integer :: iVal
    ! ---------------------------------------------------------------------------
    ! We extract momentum information completely on the view of the source
    ! coordinate system
    ! Set the right hand side of the equation system
    ! Solve the problem, where b = rhs, x = coefficients
    ! A*x = b
    ! overdetermined, solve the least Square fit problem
    ! (A^T)A*x = (A^T)b
    ! x = ((A^T)A)^-1*(A^T)b
    ! Solve linear system of equation with inverted matrix
    do iVal = 1, nVals
      a(:) = matmul( LSFmat%A, srcMom(iVal, :) )

      ! Evaluate the bubble function with the above calculated coefficients
      ! m_pi(x) = a1+a2*xcoord+a3*ycoord+a4*xcoord*xcoord+a5*ycoord*ycoord+a6*xcoord*ycoord;
      phi( iVal ) =   a( 1) &
        &           + a( 2)*targetCoord(c_x)                  &
        &           + a( 3)*targetCoord(c_y)                  &
        &           + a( 4)*targetCoord(c_x)*targetCoord(c_x) &
        &           + a( 5)*targetCoord(c_y)*targetCoord(c_y) &
        &           + a( 6)*targetCoord(c_x)*targetCoord(c_y)
    enddo

  end function mus_interpolate_quad2D_leastSq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Triquadratic interpolation for a vector quantity phi
  !! Each phi corresponds to each moment
  ! matrixInv%A(1:10,1:19) =
  !      1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19
  ! [ 4/21  4/21  4/21  4/21  4/21  4/21 -1/21 -1/21 -1/21 -1/21 -1/21 -1/21 -1/21 -1/21 -1/21 -1/21 -1/21 -1/21   3/7]
  ! [-1/10     0     0  1/10     0     0     0     0     0     0 -1/10  1/10 -1/10  1/10 -1/10 -1/10  1/10  1/10     0]
  ! [    0 -1/10     0     0  1/10     0 -1/10 -1/10  1/10  1/10     0     0     0     0 -1/10  1/10 -1/10  1/10     0]
  ! [    0     0 -1/10     0     0  1/10 -1/10  1/10 -1/10  1/10 -1/10 -1/10  1/10  1/10     0     0     0     0     0]
  ! [ 1/42  -1/7  -1/7  1/42  -1/7  -1/7 -1/21 -1/21 -1/21 -1/21  5/42  5/42  5/42  5/42  5/42  5/42  5/42  5/42 -5/21]
  ! [ -1/7  1/42  -1/7  -1/7  1/42  -1/7  5/42  5/42  5/42  5/42 -1/21 -1/21 -1/21 -1/21  5/42  5/42  5/42  5/42 -5/21]
  ! [ -1/7  -1/7  1/42  -1/7  -1/7  1/42  5/42  5/42  5/42  5/42  5/42  5/42  5/42  5/42 -1/21 -1/21 -1/21 -1/21 -5/21]
  ! [    0     0     0     0     0     0     0     0     0     0     0     0     0     0   1/4  -1/4  -1/4   1/4     0]
  ! [    0     0     0     0     0     0   1/4  -1/4  -1/4   1/4     0     0     0     0     0     0     0     0     0]
  ! [    0     0     0     0     0     0     0     0     0     0   1/4  -1/4  -1/4   1/4     0     0     0     0     0]
  !
! ****************************************************************************** !
  function mus_interpolate_quad3D_leastSq( srcMom, targetCoord, LSFmat, &
    &                                           nSources, nVals ) result(phi)
    ! ---------------------------------------------------------------------------
    !> number of quantities to interpolation
    integer, intent(in) :: nVals
    !> Number of source elements
    integer, intent(in) :: nSources
    !> source values of the momentum on the square corners
    real(kind=rk), intent(in) :: srcMom(nVals, nSources)
    !> interpolation location within the square
    real(kind=rk), intent(in) :: targetCoord(3)
    !> matrix for least square fit
    type(tem_matrix_type), intent(in) :: LSFmat
    !> interpolated value
    real(kind=rk) :: phi( nVals )
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: a(10) ! Coefficients
    integer :: iVal
    ! ---------------------------------------------------------------------------
    ! We extract momentum information completely on the view of the source
    ! coordinate system
    ! Set the right hand side of the equation system
    ! Solve the problem, where b = rhs, x = coefficients
    ! A*x = b
    ! overdetermined, solve the least Square fit problem
    ! (A^T)A*x = (A^T)b
    ! x = ((A^T)A)^-1*(A^T)b
    ! Solve linear system of equation with inverted matrix
    do iVal = 1, nVals
      a(:) = matmul( LSFmat%A, srcMom(iVal, :) )

      ! Evaluate the bubble function with the above calculated coefficients
      ! m_pi(x) = a1+a2*xcoord+a3*ycoord+a4*xcoord*xcoord+a5*ycoord*ycoord+a6*xcoord*ycoord;
      phi( iVal ) =   a( 1) &
        &           + a( 2)*targetCoord(c_x)                  &
        &           + a( 3)*targetCoord(c_y)                  &
        &           + a( 4)*targetCoord(c_z)                  &
        &           + a( 5)*targetCoord(c_x)*targetCoord(c_x) &
        &           + a( 6)*targetCoord(c_y)*targetCoord(c_y) &
        &           + a( 7)*targetCoord(c_z)*targetCoord(c_z) &
        &           + a( 8)*targetCoord(c_x)*targetCoord(c_y) &
        &           + a( 9)*targetCoord(c_y)*targetCoord(c_z) &
        &           + a(10)*targetCoord(c_z)*targetCoord(c_x)
    enddo

  end function mus_interpolate_quad3D_leastSq
! ****************************************************************************** !

end module mus_interpolate_quadratic_module
! ****************************************************************************** !
