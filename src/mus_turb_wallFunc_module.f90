! Copyright (c) 2022 Kannan Masilamani <kannan.masilaman@dlr.de>
! Copyright (c) 2023 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
! Copyright (c) 2023 Harald Klimach <harald.klimach@dlr.de>
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
!> This module contains turbulent wall function type and routines to calculate
!! friction velocity and stream-wise velocity component.
module mus_turb_wallFunc_module

  ! include treelm modules
  use env_module,               only: rk, long_k, labelLen
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit
  use tem_float_module,         only: operator(.fle.)

  ! include aotus modules
  use aotus_module,     only: flu_State, aoterr_Fatal, aoterr_NonExistent, &
    &                         aoterr_WrongType, aot_get_val

  use mus_wall_function_abstract_module, only: mus_wall_function_type
  use mus_wall_function_musker_module, only: mus_wall_function_musker_type
  use mus_wall_function_schmitt_module, only: mus_wall_function_schmitt_type, &
    &                                         get_uTau_logLayer,              &
    &                                         get_uTau_subVisousLayer,        &
    &                                         sc_lLmt, sc_uLmt
  use mus_wall_function_reichardt_module, only: mus_wall_function_reichardt_type

  implicit none

  private

  public :: mus_turb_wallFunc_type
  public :: mus_load_turb_wallFunc

  !> Contains friction velocity and turbulent viscosity on boundary elements on
  !! each level
  type mus_turb_wallFunc_data_type
    !> Turbulent viscosity on boundary element computed using mixing length
    !! formulation in lattice unit
    !! nu_t = (vonKarman*distToBnd)**2 * |du/dy|
    real(kind=rk), allocatable :: tVisc(:)

    !> Friction velocity on first neighbor in normal direction in lattice unit
    !! computed from wall model
    real(kind=rk), allocatable :: velTau(:)

    !> Distance to boundary from first fluid in normal direction
    !! in lattice unit.
    real(kind=rk), allocatable :: distToBnd(:)

    !> Distance to boundary from first fluid neighbor in normal direction
    !! in lattice unit.
    real(kind=rk), allocatable :: neighDistToBnd(:)

    !> Unit normal for each boundary element.
    !! Size: (3, nElems)
    real(kind=rk), allocatable :: unitNormal(:,:)

    !> Force on each boundary element.
    !! Size: (3, nElems)
    real(kind=rk), allocatable :: bndForce(:,:)

    !> Moment on each boundary element.
    !! Size: (3, nElems)
    real(kind=rk), allocatable :: bndMoment(:,:)
  end type mus_turb_wallFunc_data_type


  !> Contains function pointers to compute friction velocity and stream-wise
  !! velocity component
  type mus_turb_wallFunc_type
    !> is true if wall function is active
    logical :: isActive = .false.

    !> Wall model function
    character(len=labelLen) :: wall_func

    !> Nonlinear solver type
    character(len=labelLen) :: nonlinear_solver

    !> Von-Karman constant. Default = 0.4_rk
    real(kind=rk) :: vonKarman = 0.4_rk

    !> Use vanDriest damping function to damp turbulent viscosity
    logical :: useVanDriest = .true.

    !> Contains data computed in turbulent wall bc routine on each level
    type(mus_turb_wallFunc_data_type), allocatable :: dataOnLvl(:)

    !> Function pointer to compute friction velocity
    procedure(mus_proc_calcFricVel), pointer, pass(this) :: calcFricVel => null()

    !> Function pointer to compute strean-wise velocity component
    procedure(mus_proc_calcStreamWiseVel), pointer, nopass :: &
      & calcStreamWiseVel => null()

    !> Allocate wall function object
    class(mus_wall_function_type), allocatable :: wall_function

    !> Function pointer to the iterative method
    procedure(mus_iterative_method_interface), pointer, nopass :: &
      & iterativeMethod => null()
  end type mus_turb_wallFunc_type


  !> Interface definition for the turbulent wall bc routines
  abstract interface
    !> This abstract interface defines the interface to calculate turbulent wall
    !! friction velocity from given velocity and distance to boundary.
    !! All inputs and output are in lattice units.
    pure subroutine mus_proc_calcFricVel(this, velTau, velSW, distToBnd, &
      &             viscKine, nElems)
      import :: rk, mus_wall_function_type, mus_turb_wallFunc_type
      !> Turbulent wall model to use for the computation
      class(mus_turb_wallFunc_type), intent(in) :: this
      !> Friction velocity computed from wall model.
      !! it is inout to provide velTau from previous timestep as initial velTau
      !! for fixed-point or Newton iteration solver
      real(kind=rk), intent(inout) :: velTau(:)
      !> Stream-wise velocity component from which friction velocity is computed
      real(kind=rk), intent(in) :: velSW(:)
      !> Distance to the boundary in the discrete normal direction
      real(kind=rk), intent(in) :: distToBnd(:)
      !> Kinematic viscosity
      real(kind=rk), intent(in) :: viscKine(:)
      !> Number of elements in input and output arrays
      integer, intent(in) :: nElems
    end subroutine mus_proc_calcFricVel

    !> This abstract interface defines the interface to calculate stream-wise
    !! velocity component from friction velocity and distance to boundary.
    !! All inputs and output are in lattice units.
    pure subroutine mus_proc_calcStreamWiseVel(velSW, velTau, distToBnd,      &
      &                                   viscKine, nElems, wall_function)
      import :: rk, mus_wall_function_type
      !> Stream-wise velocity component from wall model
      real(kind=rk), intent(out) :: velSW(:)
      !> Friction velocity computd from wall model
      real(kind=rk), intent(in) :: velTau(:)
      !> Distance to the boundary in the discrete normai direction
      real(kind=rk), intent(in) :: distToBnd(:)
      !> Kinematic viscosity
      real(kind=rk), intent(in) :: viscKine(:)
      !> Number of elements in input and output arrays
      integer, intent(in) :: nElems
      !> Allocate wall function object
      class(mus_wall_function_type), intent(in) :: wall_function
    end subroutine mus_proc_calcStreamWiseVel

    !> This routine computes friction velocity from wall model profile
    !! using Newton iteration method
    pure function mus_iterative_method_interface(velTau_initialGuess, velSW, &
      &                              y, nu, wall_function) result (velTau_new)
      ! -------------------------------------------------------------------- !
      import :: rk, mus_wall_function_type
      !> Friction velocity computed from previsous time step
      real(kind=rk), intent(in) :: velTau_initialGuess
      !> Stream-wise velocity component from which friction velocity is computed
      real(kind=rk), intent(in) :: velSW
      !> vertical distance from the wall
      real(kind=rk), intent(in) :: y
      !> dynamic viscosity
      real(kind=rk), intent(in) :: nu
      !> Number of elements in input and output arrays
      class(mus_wall_function_type), intent(in) :: wall_function
      !> Friction velocity computed in this routine
      real(kind=rk) :: velTau_new
    end function mus_iterative_method_interface
  end interface

  !! Constant parameters for Implicit equation solver
  integer, parameter :: imEq_nIter = 1000
  real(kind=rk), parameter :: imEq_tol = 1e-10


contains


  ! -------------------------------------------------------------------------- !
  !> This routine loads wall model and nonlinear solver type for nonlinear
  !! equation
  subroutine mus_load_turb_wallFunc(me, conf, parent)
    ! -------------------------------------------------------------------- !
    !> Turbulent wall model type to fill assign wall model
    type(mus_turb_wallFunc_type), intent(inout) :: me
    !> lua flu state
    type( flu_state ) :: conf
    !> bc parent handle
    integer, intent(in) :: parent
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !
    write(logUnit(2), "(A)") '  Loading turbulent wall_function >>>'
    me%isActive = .true.

    call load_wall_function(me, conf, parent)
    call load_iterativeMethod(me, conf, parent)

    ! Von-Karman constant
    call aot_get_val(L       = conf,         &
      &              thandle = parent,       &
      &              key     = 'von_karman', &
      &              val     = me%vonKarman, &
      &              default = 0.4_rk,       &
      &              ErrCode = iError        )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving ' &
        &              //'von_karman:'
      if (btest(iError, aoterr_WrongType)) then
         write(logUnit(1),*)'Variable has wrong type (should be number)!'
         write(logUnit(1),*)'STOPPING'
         call tem_abort()
      end if
    end if
    write(logUnit(2),*) '  von Karman: ', me%vonKarman

    ! Use van Driest damping function
    call aot_get_val(L       = conf,             &
      &              thandle = parent,           &
      &              key     = 'use_van_driest', &
      &              val     = me%useVanDriest,  &
      &              default = .true.,           &
      &              ErrCode = iError            )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving ' &
        &              //'use_van_driest:'
      if (btest(iError, aoterr_WrongType)) then
         write(logUnit(1),*)'Variable has wrong type (should be bool)!'
         write(logUnit(1),*)'STOPPING'
         call tem_abort()
      end if
    end if
    write(logUnit(2),*) '  use van Driest: ', me%useVanDriest

    me%calcStreamWiseVel => computeStreamWiseVel
    write(logUnit(2), "(A)") '  <<< Done loading turbulent wall_function'

  end subroutine mus_load_turb_wallFunc
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> Load the iterativeMethod to use in the turbulent wall model from the user
  !! configuration.
  subroutine load_iterativeMethod(me, conf, parent)
    ! -------------------------------------------------------------------- !
    !> Turbulent wall model type to fill assign wall model
    type(mus_turb_wallFunc_type), intent(inout) :: me
    !> lua flu state
    type( flu_state ) :: conf
    !> bc parent handle
    integer, intent(in) :: parent
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    call aot_get_val(L       = conf,                &
      &              thandle = parent,              &
      &              key     = 'nonlinear_solver',  &
      &              val     = me%nonlinear_solver, &
      &              default = 'fixed_point',       &
      &              ErrCode = iError               )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving ' &
        &                 // 'nonlinear_solver:'
      if (btest(iError, aoterr_WrongType)) then
         write(logUnit(1),*)'Variable has wrong type (should be string)!'
         write(logUnit(1),*)'STOPPING'
         call tem_abort()
      end if
    end if

    select case (trim(me%nonlinear_solver))
    case ('fixed_point')
      me%iterativeMethod => fixedPoint_method

    case ('newton')
      me%iterativeMethod => newton_method

    case default
      call tem_abort('Error: Unsupported nonlinear solver '           &
        &            //trim(me%nonlinear_solver)//' for wall function')

    end select
    write(logUnit(2),*) '  nonlinear solver: ', trim(me%nonlinear_solver)

  end subroutine load_iterativeMethod
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> Load the iterativeMethod to use in the turbulent wall model from the user
  !! configuration.
  subroutine load_wall_function(me, conf, parent)
    ! -------------------------------------------------------------------- !
    !> Turbulent wall model type to fill assign wall model
    type(mus_turb_wallFunc_type), intent(inout) :: me
    !> lua flu state
    type( flu_state ) :: conf
    !> bc parent handle
    integer, intent(in) :: parent
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    call aot_get_val(L       = conf,            &
      &              thandle = parent,          &
      &              key     = 'wall_function', &
      &              val     = me%wall_func,    &
      &              default = 'musker',        &
      &              ErrCode = iError           )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving wall_function:'
      if (btest(iError, aoterr_WrongType)) then
         write(logUnit(1),*)'Variable has wrong type (should be a string)!'
         write(logUnit(1),*)'STOPPING'
         call tem_abort()
      end if
    end if

    select case (trim(me%wall_func))
    case ('musker')
      allocate(mus_wall_function_musker_type :: me%wall_function)
      me%calcFricVel => compute_fricVel

    case ('reichardt')
      allocate(mus_wall_function_reichardt_type :: me%wall_function)
      me%calcFricVel => compute_fricVel

    case ('power_law')
      allocate(mus_wall_function_schmitt_type :: me%wall_function)
      me%calcFricVel => fricVel_Schmitt

    case default
      call tem_abort("Error: Unknown turbulent wall function " &
        &            //trim(me%wall_func)                      )
    end select

    write(logUnit(2),*) '  wall function: ', trim(me%wall_func)

  end subroutine load_wall_function
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This routine computes friction velocity from Schmitt wall model.
  ! I am forced to make an extra subroutine because the wall profile is a
  ! combination of implicit and explicit functions
  pure subroutine fricVel_Schmitt(this, velTau, velSW, distToBnd, viscKine, &
    &                             nElems)
    ! -------------------------------------------------------------------- !
    !> Pass the calling object as an argument
    class(mus_turb_wallFunc_type), intent(in) :: this
    !> Friction velocity computed from wall model.
    !! it is inout to provide velTau from previous timestep as initial velTau
    !! for fixed-point or Newton iteration solver
    real(kind=rk), intent(inout) :: velTau(:)
    !> Stream-wise velocity component from which friction velocity is computed
    real(kind=rk), intent(in) :: velSW(:)
    !> Distance to the boundary in the discrete normai direction
    real(kind=rk), intent(in) :: distToBnd(:)
    !> Kinematic viscosity
    real(kind=rk), intent(in) :: viscKine(:)
    !> Number of elements in input and output arrays
    integer, intent(in) :: nElems
    ! -------------------------------------------------------------------- !
    integer :: iElem
    real(kind=rk) :: visc_div_dist, yPlus
    ! -------------------------------------------------------------------- !
    do iElem = 1, nElems
      visc_div_dist = viscKine(iElem) / distToBnd(iElem)
      yPlus = velTau(iElem) / visc_div_dist

      if (yPlus >= sc_uLmt) then
        ! Interial layer, use powerlaw profile from werner and wengle
        velTau(iElem) = get_uTau_logLayer( visc_div_dist = visc_div_dist, &
        &                                  velSW = velSW(iElem)           )
      else if ( yPlus < sc_lLmt) then
        ! viscous sublayer, use linear profile.
        velTau(iElem) = get_uTau_subVisousLayer( visc_div_dist = visc_div_dist, &
          &                                      velSW = velSW(iElem)           )
      else ! if ( yPlus >= sc_lLmt .and. yPlus < sc_uLmt)
        ! Buffer layer, use logarithmic profile
        ! the class Schmitt wall function has been modified such that the
        ! derivative of uPlus with respect to uTau is always in the buffer layer
        velTau(iElem) = this%iterativeMethod (                                &
          &                              velTau_initialGuess = velTau(iElem), &
          &                              velSW = velSW(iElem),                &
          &                              y = distToBnd(iElem),                &
          &                              nu = viscKine(iElem),                &
          &                              wall_function = this%wall_function   )
      end if
    end do

  end subroutine fricVel_Schmitt
  ! -------------------------------------------------------------------------- !


  ! -------------------------------------------------------------------------- !
  !> This routine computes friction velocity from wall model profile
  !! using Newton iteration method
  pure subroutine compute_fricVel(this, velTau, velSW, distToBnd, viscKine, &
    &                             nElems)
    ! -------------------------------------------------------------------- !
    !> Pass the calling object as an argument
    class(mus_turb_wallFunc_type), intent(in) :: this
    !> Friction velocity computed from wall model.
    !! it is inout to provide velTau from previous timestep as initial velTau
    !! for fixed-point or Newton iteration solver
    real(kind=rk), intent(inout) :: velTau(:)
    !> Stream-wise velocity component from which friction velocity is computed
    real(kind=rk), intent(in) :: velSW(:)
    !> Distance to the boundary in the discrete normai direction
    real(kind=rk), intent(in) :: distToBnd(:)
    !> Kinematic viscosity
    real(kind=rk), intent(in) :: viscKine(:)
    !> Number of elements in input and output arrays
    integer, intent(in) :: nElems
    ! -------------------------------------------------------------------- !
    integer :: iElem
    ! -------------------------------------------------------------------- !

    do iElem = 1, nElems
      velTau(iElem) = this%iterativeMethod (                                &
        &                              velTau_initialGuess = velTau(iElem), &
        &                              velSW = velSW(iElem),                &
        &                              y = distToBnd(iElem),                &
        &                              nu = viscKine(iElem),                &
        &                              wall_function = this%wall_function   )
    end do

  end subroutine compute_fricVel
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This routine computes friction velocity from wall model profile
  !! using Newton iteration method
  pure function newton_method(velTau_initialGuess, velSW, y, nu, &
    &                         wall_function) result (velTau_new)
    ! -------------------------------------------------------------------- !
    !> Friction velocity computed from previsous time step
    real(kind=rk), intent(in) :: velTau_initialGuess
    !> Stream-wise velocity component from which friction velocity is computed
    real(kind=rk), intent(in) :: velSW
    !> vertical distance from the wall
    real(kind=rk), intent(in) :: y
    !> dynamic viscosity
    real(kind=rk), intent(in) :: nu
    !> Number of elements in input and output arrays
    class(mus_wall_function_type), intent(in) :: wall_function
    !> Friction velocity computed in this routine
    real(kind=rk) :: velTau_new
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: error
    integer :: iter
    real(kind=rk) :: yPlus, yPlus_fac, velTau_old
    real(kind=rk) :: fx, dfdx
    ! -------------------------------------------------------------------- !
    iter = 0
    yPlus_fac = y / nu
    velTau_old = velTau_initialGuess

    newton: do
      yPlus = velTau_old * yPlus_fac
      ! Function of u_tau
      fx = velSW / velTau_old - wall_function%get_uPlus(yPlus)
      ! derivative of function w.r.t u_tau
      dfdx = -velSW / velTau_old**2                                  &
        &    - (wall_function%get_d_uPlus_d_uTau( y = y,             &
        &                                         uTau = velTau_old, &
        &                                         nu = nu            ) )
      velTau_new = velTau_old - (fx/dfdx)
      error = abs(velTau_new - velTau_old)
      velTau_old = velTau_new
      iter = iter+1
      if ((error .fle. imEq_tol)    &
        & .or. (iter >= imEq_nIter)) exit newton
    end do newton

  end function newton_method
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This routine computes friction velocity from wall model profile
  !! using fixed-point iterative method
  pure function fixedPoint_method(velTau_initialGuess, velSW, y, nu, &
    &                             wall_function) result (velTau_new)
    ! -------------------------------------------------------------------- !
    !> Friction velocity computed from previsous time step
    real(kind=rk), intent(in) :: velTau_initialGuess
    !> Stream-wise velocity component from which friction velocity is computed
    real(kind=rk), intent(in) :: velSW
    !> vertical distance from the wall
    real(kind=rk), intent(in) :: y
    !> dynamic viscosity
    real(kind=rk), intent(in) :: nu
    !> Number of elements in input and output arrays
    class(mus_wall_function_type), intent(in) :: wall_function
    !> Friction velocity computed in this routine
    real(kind=rk) :: velTau_new
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: error
    integer :: iter
    real(kind=rk) :: velTau_old, yPlus, yPlus_fac
    ! -------------------------------------------------------------------- !
    iter = 0
    yPlus_fac = y / nu
    velTau_old = velTau_initialGuess

    fixed_point: do
      yPlus = velTau_old * yPlus_fac
      velTau_new = velSW / wall_function%get_uPlus(yPlus)
      error = abs(velTau_new - velTau_old)
      velTau_old = velTau_new
      iter = iter+1
      if ((error .fle. imEq_tol)    &
        & .or. (iter >= imEq_nIter)) exit fixed_point
    end do fixed_point

  end function fixedPoint_method
  ! -------------------------------------------------------------------------- !

  ! -------------------------------------------------------------------------- !
  !> This routines computes streamWise velocity component from friction velocity
  !! and distance to boundary using any wall function profile.
  pure subroutine computeStreamWiseVel(velSW, velTau, distToBnd, viscKine, &
    &                                  nElems, wall_function)
    ! -------------------------------------------------------------------- !
    !> Stream-wise velocity component from wall model
    real(kind=rk), intent(out) :: velSW(:)
    !> Friction velocity computd from wall model
    real(kind=rk), intent(in) :: velTau(:)
    !> Distance to the boundary in the discrete normai direction
    real(kind=rk), intent(in) :: distToBnd(:)
    !> Kinematic viscosity
    real(kind=rk), intent(in) :: viscKine(:)
    !> Number of elements in input and output arrays
    integer, intent(in) :: nElems
    !> Allocate wall function object
    class(mus_wall_function_type), intent(in) :: wall_function
    ! -------------------------------------------------------------------- !
    integer :: iElem
    real(kind=rk) :: yPlus
    ! -------------------------------------------------------------------- !
    do iElem = 1, nElems
      yPlus = distToBnd(iElem) * velTau(iElem) / viscKine(iElem)
      velSW(iElem) = wall_function%get_uPlus(yPlus) * velTau(iElem)
    end do
  end subroutine computeStreamWiseVel
  ! -------------------------------------------------------------------------- !

end module mus_turb_wallFunc_module
