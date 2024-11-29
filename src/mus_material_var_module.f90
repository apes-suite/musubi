! Copyright (c) 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!> author: Kannan Masilamani
!! This module provides variable to extract material variable like viscosity,
!! omega, etc.
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015, 2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
!--------------------------------------------
!    A O S - Array of structures layout new
!-------------------------------------------
! Access to get_point value output
! Access to get_element value output
module mus_material_var_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,               only: rk, labelLen
  use tem_logging_module,       only: logUnit
  use tem_param_module,         only: cs2, rho0
  use tem_variable_module,      only: tem_variable_type
  use tem_topology_module,      only: tem_levelOf
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type, tem_varSys_op_type,     &
    &                                 tem_varSys_append_derVar,                &
    &                                 tem_varSys_proc_point,                   &
    &                                 tem_varSys_proc_element,                 &
    &                                 tem_varSys_proc_setParams,               &
    &                                 tem_varSys_proc_getParams,               &
    &                                 tem_varSys_proc_setupIndices,            &
    &                                 tem_varSys_proc_getValOfIndex,           &
    &                                 tem_varSys_setupIndices_dummy,           &
    &                                 tem_varSys_getValOfIndex_dummy,          &
    &                                 tem_varSys_setParams_dummy,              &
    &                                 tem_varSys_getParams_dummy
  use tem_aux_module,           only: tem_abort
  use tem_grow_array_module,    only: grw_labelarray_type, append

  ! include musubi modules
  use mus_scheme_type_module,     only: mus_scheme_type
  use mus_scheme_header_module,   only: mus_scheme_header_type
  use mus_varSys_module,          only: mus_varSys_data_type,       &
    &                                   mus_varSys_solverData_type, &
    &                                   mus_get_new_solver_ptr,     &
    &                                   mus_deriveVar_forPoint

  implicit none
  private

  public :: mus_append_materialVar

contains

  ! ************************************************************************** !
  !> subroutine to add material variable
  subroutine mus_append_materialVar(varSys, solverData, schemeHeader, derVarName)
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)  :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in)  :: schemeHeader

    !> array of derive physical variables
    type(grw_labelarray_type), intent(inout) :: derVarName
    ! --------------------------------------------------------------------------
    ! number of derive variables
    integer :: nDerVars, iVar
    integer :: nComponents, addedPos
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    character(len=labelLen), allocatable :: derVarName_loc(:)
    ! --------------------------------------------------------------------------
    select case (trim(schemeHeader%kind))
    case ('fluid', 'fluid_incompressible')
      nDerVars = 3
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'kine_viscosity    ', 'pressure_reference', &
        &                   'omega             ']
    case default
      nDerVars = 0
      write(logUnit(3),*) 'No material variables are defined for chosen '&
        &               //'scheme kind'
    end select

    do iVar = 1, nDerVars
      call append(derVarName, derVarName_loc(iVar))
      ! assign function pointers only for get_element and get_point because
      ! this variable will be only used for tracking
      get_point => mus_deriveVar_forPoint
      ! for other function pointers assign dummy routines
      setup_indices => tem_varSys_setupIndices_dummy
      get_valOfIndex => tem_varSys_getValOfIndex_dummy
      set_params => tem_varSys_setParams_dummy
      get_params => tem_varSys_getParams_dummy

      select case (trim(adjustl(derVarName_loc(iVar))))
      case ('kine_viscosity')
        get_element => access_kineVisc_forElement
        nComponents = 1
      case ('omega')
        get_element => access_kineOmega_forElement
        nComponents = 1
      case ('pressure_reference')
        get_element => access_pressRef_forElement
        nComponents = 1
      case default
          call tem_abort( 'Error: Unknown material variable')
      end select

      ! append variable to varSys
      call tem_varSys_append_derVar(                            &
        &  me             = varSys,                             &
        &  varName        = trim(derVarName_loc(iVar)),         &
        &  nComponents    = nComponents,                        &
        &  method_data    = mus_get_new_solver_ptr(solverData), &
        &  get_point      = get_point,                          &
        &  get_element    = get_element,                        &
        &  set_params     = set_params,                         &
        &  get_params     = get_params,                         &
        &  setup_indices  = setup_indices,                      &
        &  get_valOfIndex = get_valOfIndex,                     &
        &  pos            = addedPos,                           &
        &  wasAdded       = wasAdded                            )

      if (wasAdded) then
        write(logUnit(10),*) ' Appended variable: '//trim(derVarName_loc(iVar))
      else if (addedpos < 1) then
        write(logUnit(0),*) 'Error: variable '//trim(derVarName_loc(iVar)) &
          &         // ' is not added to variable system'
        call tem_abort()
      end if

    end do

  end subroutine mus_append_materialVar
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine returns the kinematic viscosity
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine access_kineVisc_forElement(fun, varsys, elempos, time, &
    &                                             tree, nElems, nDofs, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: statePos, iElem, iLevel
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    ! res is always AOS layout
    res = 0.0_rk
    do iElem = 1, nElems
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      statePos = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
      res( iElem ) = scheme%field(1)%fieldProp%fluid%viscKine &
        &                           %dataOnLvl(iLevel)%val(statePos)
    end do !iElem

  end subroutine access_kineVisc_forElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine returns the omega
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine access_kineOmega_forElement(fun, varsys, elempos, time, &
    &                                              tree, nElems, nDofs, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: statePos, iElem, iLevel
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    ! res is always AOS layout
    res = 0.0_rk
    do iElem = 1, nElems
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      statePos = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
      res( iElem ) = scheme%field(1)%fieldProp%fluid%viscKine &
        &                           %omLvl(iLevel)%val(statePos)
    end do !iElem

  end subroutine access_kineOmega_forElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine returns the reference pressure
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine access_pressRef_forElement(fun, varsys, elempos, time, &
    &                                             tree, nElems, nDofs, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Position of the TreeID of the element to get the variable for in the
    !! global treeID list.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    ! res is always AOS layout
    res = rho0 * cs2

  end subroutine access_pressRef_forElement
  ! ************************************************************************** !

end module mus_material_var_module

