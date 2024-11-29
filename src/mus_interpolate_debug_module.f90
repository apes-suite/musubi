! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2015, 2018-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> author: Manuel Hasert
!! Interpolation of flow quantities between different grid levels
!!
!! # Interpolation
!! The routines defined here, fill up the ghost elements with valid data.
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
!! # Workflow
!!
!! Each interpolation routine acts on a list of ghost elements.
!! This list contains pointers to the position in the total list.
!! For each of these ghost elements, the source elements are identified.
!! Before that, the sourceLevel is identified. However, the code is restricted
!! to work with a level jump of only one level, so the sourceLevel is
!! for sourceLevel = targetLevel+1
!! sourceLevel = targetLevel-1
!!
!! For an overview over implemented interpolation methods, see
!! [Interpolation methods](../page/features/intp_methods.html)
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
module mus_interpolate_debug_module
  use iso_c_binding, only: c_f_pointer

  ! include treelm modules
  use env_module,            only: rk, long_k
  use tem_param_module,      only: cs2inv, PI
  use tem_topology_module,   only: tem_coordOfID
  use tem_element_module,    only: eT_GhostFromCoarser,eT_ghostFromFiner
  use tem_logging_module,    only: logUnit
  use tem_debug_module,      only: dbgUnit
  use tem_construction_module, only: tem_levelDesc_type
  use tem_varSys_module,       only: tem_varSys_type
  use tem_time_module,       only: tem_time_type
  use tem_stencil_module,    only: tem_stencilHeader_type
  use tem_construction_module, only: tem_levelDesc_type

  ! include musubi modules
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_physics_module,            only: mus_physics_type
  use mus_interpolate_header_module, only: mus_interpolation_method_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_fluid_module,              only: mus_fluid_type

  implicit none

  private

  public :: do_nothing
  public :: do_nothing_arbi


  public :: TGV_2D

  real(kind=rk), parameter :: debugValue = 1.0_rk

 contains

! ****************************************************************************** !
  !> Fill GhostFromFiner elements on my level with debug value
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine]] in intp/[[mus_interpolate_header_module]].f90 in order to
  !! be callable via [[mus_interpolation_method_type:do_intp]] function pointer.
  subroutine do_nothing( method, fieldProp, tLevelDesc, level, sState, sNeigh,  &
    &  snSize, sAuxField, tState, tNeigh, tnSize, layout, nTargets, targetList, &
    &  physics, time, varSys, derVarPos                                         )
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

    !> AuxField variable to read rho and vel from source elements
    real(kind=rk), intent(inout) :: sAuxField(:)

    !> State vector of TARGET GHOST elements
    real(kind=rk), intent(inout) :: tState(:)
    integer, intent(in) :: tNeigh(:)
    integer, intent(in) :: tnSize

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

  end subroutine do_nothing
! ****************************************************************************** !

! ****************************************************************************** !
  !> Fill GhostFromFiner elements on my level with debug value
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[intpRoutine_arbitraryVal]] in intp/[[mus_interpolate_header_module]].f90
  !! in order to be callable via
  !! [[mus_interpolation_method_type:do_intpArbiVal]] function pointer.
  subroutine do_nothing_arbi( method, tLevelDesc, level, stencil, sVal, &
    &                         tVal, nTargets, targetList, nScalars      )
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

  end subroutine do_nothing_arbi
! ****************************************************************************** !



! ****************************************************************************** !
  !> This routine returns the analytical solution of TGV 2D testcase for a given
  !! position and time (coord, t)
  !!
  pure function TGV_2D( coord, t ) result ( res )
    ! ---------------------------------------------------------------------------
    !> position and time
    real(kind=rk), intent(in) :: coord(3), t
    !> pressure, velX, velY, Sxx, Syy, Sxy
    real(kind=rk)             :: res(6)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: nuPhy, tD, x, y
    ! ---------------------------------------------------------------------------

    x = coord(1)
    y = coord(2)
    ! viscosity = 1 / Re
    nuPhy = 1.0 / 25.0_rk
    tD = 1.0 / nuPhy / 2

    ! pressure
    res(1) = -(cos(2.0*x)+cos(2.0*y)) * exp(-2.0*t/tD) / 4.0
    ! velocity
    res(2) = -cos(x) * sin(y) * exp(-t/tD)
    res(3) =  sin(x) * cos(y) * exp(-t/tD)
    ! strain rate S
    res(4) =  sin(x) * sin(y) * exp(-t/tD)
    res(5) = -sin(x) * sin(y) * exp(-t/tD)
    res(6) = 0.0_rk

  end function TGV_2D
! ****************************************************************************** !
end module mus_interpolate_debug_module
! ****************************************************************************** !
