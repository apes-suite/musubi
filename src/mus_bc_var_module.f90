! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
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
!> author: Kannan Masilamani
!! This module provides variable to extract from boundary condition
module mus_bc_var_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,               only: rk, labelLen
  use tem_logging_module,       only: logUnit
  use tem_variable_module,      only: tem_variable_type
  use tem_topology_module,      only: tem_levelOf
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_property_module,      only: prp_hasBnd
  use tem_varSys_module,        only: tem_varSys_type, tem_varSys_op_type,     &
    &                                 tem_varSys_append_derVar,                &
    &                                 tem_varSys_proc_point,                   &
    &                                 tem_varSys_proc_element,                 &
    &                                 tem_varSys_proc_setParams,               &
    &                                 tem_varSys_proc_getParams,               &
    &                                 tem_varSys_proc_setupIndices,            &
    &                                 tem_varSys_proc_getValOfIndex,           &
    &                                 tem_varSys_getElement_dummy,             &
    &                                 tem_varSys_setupIndices_dummy,           &
    &                                 tem_varSys_getValOfIndex_dummy,          &
    &                                 tem_varSys_setParams_dummy,              &
    &                                 tem_varSys_getParams_dummy
  use tem_aux_module,           only: tem_abort
  use tem_grow_array_module,    only: grw_labelArray_type, append
  use tem_dyn_array_module,     only: dyn_labelArray_type, init, append
  use tem_stencil_module,       only: tem_stencilHeader_type

  ! include musubi modules
  use mus_field_module,       only: mus_field_type
  use mus_varSys_module,      only: mus_varSys_data_type,             &
    &                               mus_varSys_solverData_type,       &
    &                               mus_get_new_solver_ptr,           &
    &                               mus_deriveVar_forPoint
  use mus_bc_header_module,   only: glob_boundary_type

  implicit none
  private

  public :: mus_append_bcVar

contains

  ! ************************************************************************** !
  !> This routine adds boundary variables for tracking
  subroutine mus_append_bcVar(varSys, solverData, derVarName, nFields, &
    &                         field, stencil)
    ! --------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)  :: varSys
    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData
    !> array of derive physical variables
    type(grw_labelArray_type), intent(inout) :: derVarName
    !> number of fields
    integer, intent(in)                      :: nFields
    !> Field contains sources and boundary infos
    type(mus_field_type), intent(in) :: field(nFields)
    !> Compute stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    ! --------------------------------------------------------------------------
    integer :: iField, iBC, iVar, nComponents, addedPos, iIn
    logical :: wasAdded
    character(len=labelLen), allocatable ::  input_varname(:)
    character(len=labelLen)  ::  varName
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    type(c_ptr) :: method_data
    type(dyn_labelArray_type) :: bcVarName
    integer :: nBCtypes
    ! --------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    write(logUnit(1),*) 'Appending boundary variables '
    ! all fields have same number of bc
    nBCtypes = size(field(1)%bc)

    ! First list of bc variables to append
    call init(me=bcVarName, length=0)
    do iField = 1, nFields
      do iBC = 1, nBCtypes
        select case (trim(field(iField)%bc(iBC)%BC_kind))
        case ('wall_libb', 'velocity_bounceback', 'velocity_bfl', &
          & 'velocity_noneq_expol')
          call append( me       = bcVarName,   &
            &          val      = 'bc_normal', &
            &          pos      = addedPos,    &
            &          wasAdded = wasAdded     )
          call append( me       = bcVarName,   &
            &          val      = 'bc_qval',   &
            &          pos      = addedPos,    &
            &          wasAdded = wasAdded     )

        case ('turbulent_wall', 'turbulent_wall_noneq_expol', &
          & 'turbulent_wall_eq')
          call append( me       = bcVarName,   &
            &          val      = 'bc_normal', &
            &          pos      = addedPos,    &
            &          wasAdded = wasAdded     )
          call append( me       = bcVarName,          &
            &          val      = 'bc_fric_velocity', &
            &          pos      = addedPos,           &
            &          wasAdded = wasAdded            )
          call append( me       = bcVarName,          &
            &          val      = 'bc_turb_viscosity', &
            &          pos      = addedPos,            &
            &          wasAdded = wasAdded             )
          call append( me       = bcVarName,             &
            &          val      = 'bc_norm_dist_to_bnd', &
            &          pos      = addedPos,              &
            &          wasAdded = wasAdded               )
          call append( me       = bcVarName,   &
            &          val      = 'bc_y_plus', &
            &          pos      = addedPos,    &
            &          wasAdded = wasAdded     )
          call append( me       = bcVarName,   &
            &          val      = 'bc_qval',   &
            &          pos      = addedPos,    &
            &          wasAdded = wasAdded     )
        end select
      end do
    end do

    do iField = 1, nFields
      do iVar = 1, bcVarName%nVals
        call append(me = derVarName, val = bcVarName%val(iVar))

        ! set default pointers, overwrite if neccessary
        get_element => tem_varSys_getElement_dummy
        get_point => mus_deriveVar_forPoint
        setup_indices => tem_varSys_setupIndices_dummy
        get_valOfIndex => tem_varSys_getValOfIndex_dummy
        method_data  = mus_get_new_solver_ptr(solverData)
        set_params => tem_varSys_setParams_dummy
        get_params => tem_varSys_getParams_dummy

        select case (trim(adjustl(bcVarName%val(iVar))))
        case ('bc_normal')
          get_element => access_bcNormal_forElement
          nComponents = 3
          allocate(input_varname(1))
          input_varname(1) = 'pdf'

        case ('bc_qval')
          get_element => access_qVal_forElement
          nComponents = stencil%QQN
          allocate(input_varname(1))
          input_varname(1) = 'pdf'

        case ('bc_fric_velocity')
          get_element => access_bcFricVel_forElement
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'pdf'

        case ('bc_norm_dist_to_bnd')
          get_element => access_bcNormDistToBnd_forElement
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'pdf'

        case ('bc_turb_viscosity')
          get_element => access_bcTurbVisc_forElement
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'pdf'

        case ('bc_y_plus')
          get_element => access_bcYPlus_forElement
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'pdf'

        case default
          write(logUnit(1),*) 'WARNING: Unknown variable: '// &
            &                 trim(bcVarName%val(iVar))
          cycle !go to next variable
        end select

        ! update variable names with field label
        varname = trim(field(iField)%label)//trim(adjustl(bcVarName%val(iVar)))
        do iIn = 1, size(input_varname)
          input_varname(iIn) = trim(field(iField)%label)&
            & // trim(input_varname(iIn))
        end do

        ! append variable to varSys
        call tem_varSys_append_derVar(  me             = varSys,         &
          &                             varName        = trim(varname),  &
          &                             nComponents    = nComponents,    &
          &                             input_varname  = input_varname,  &
          &                             method_data    = method_data,    &
          &                             get_point      = get_point,      &
          &                             get_element    = get_element,    &
          &                             set_params     = set_params,     &
          &                             get_params     = get_params,     &
          &                             setup_indices  = setup_indices,  &
          &                             get_valOfIndex = get_valOfIndex, &
          &                             pos            = addedPos,       &
          &                             wasAdded       = wasAdded        )

        if (wasAdded) then
          write(logUnit(10),*) ' Appended variable:'//trim(varname)
        else if (addedpos < 1) then
          write(logUnit(1),*) 'Error: variable '//trim(varname)// &
            &                 ' is not added to variable system'
        end if

        deallocate(input_varname)
      end do
    end do

  end subroutine mus_append_bcVar
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine returns the boundary normal pointing inside the domain
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine access_bcNormal_forElement(fun, varsys, elempos, time, &
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
    integer :: iElem, iLevel, posInBndID, minBcID, bcLevelPointer
    type(mus_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! res is always AOS layout
    res = 0.0_rk
    associate (globBC => fPtr%solverData%scheme%globBC)
      do iElem = 1, nElems
        ! This routine is valid only for boundary elements
        if ( btest( tree%elemPropertyBits( elemPos(iElem) ), &
          &         prp_hasBnd )  ) then
          ! position of current element in boundary_ID list
          posInBndID = fPtr%solverData%geometry%posInBndID( elempos(iElem) )
          minBcID = fPtr%solverData%geometry%minBcID(posInBndID)
          bcLevelPointer = fPtr%solverData%geometry%bcLevelPointer(posInBndID)
          iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
          if (.not. globBC(minBcID)%isWall) then
            res( (iElem-1)*fun%nComponents+1: iElem*fun%nComponents ) &
              & = globBC(minBcID)%elemLvl(iLevel)%normal%val(:, bcLevelPointer)
          end if
        end if
      end do !iElem
    end associate

  end subroutine access_bcNormal_forElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine returns the boundary qValues
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine access_qVal_forElement(fun, varsys, elempos, time, &
    &                                         tree, nElems, nDofs, res    )
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
    integer :: iElem, iLevel, posInBndID, minBcID, bcLevelPointer
    type(mus_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! res is always AOS layout
    res = 0.0_rk
    associate (globBC => fPtr%solverData%scheme%globBC)
      do iElem = 1, nElems
        ! This routine is valid only for boundary elements
        if ( btest( tree%elemPropertyBits( elemPos(iElem) ), &
          &         prp_hasBnd )  ) then
          ! position of current element in boundary_ID list
          posInBndID = fPtr%solverData%geometry%posInBndID( elempos(iElem) )
          minBcID = fPtr%solverData%geometry%minBcID(posInBndID)
          bcLevelPointer = fPtr%solverData%geometry%bcLevelPointer(posInBndID)
          iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
          if (.not. globBC(minBcID)%isWall) then
            res( (iElem-1)*fun%nComponents+1: iElem*fun%nComponents ) &
              & = globBC(minBcID)%elemLvl(iLevel)%qVal%val(:, bcLevelPointer)
          end if
        end if
      end do !iElem
    end associate

  end subroutine access_qVal_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine returns friction velocity computed in turbulent_wall bc and
  !! stored in mus_turb_wallFunc_data_type routine.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine access_bcFricVel_forElement(fun, varsys, elempos, time, &
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
    integer :: iElem, iLevel, iField, posInBndID, minBcID, bcLevelPointer
    type(mus_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    ! Position of dependent pdf is iField
    iField = fun%input_varPos(1)

    ! res is always AOS layout
    res = 0.0_rk
    associate (field => fPtr%solverData%scheme%field(iField))
      do iElem = 1, nElems
        ! This routine is valid only for turbulent wall boundary elements
        if ( btest( tree%elemPropertyBits( elemPos(iElem) ), &
          &         prp_hasBnd )  ) then
          ! position of current element in boundary_ID list
          posInBndID = fPtr%solverData%geometry%posInBndID( elempos(iElem) )
          minBcID = fPtr%solverData%geometry%minBcID(posInBndID)
          bcLevelPointer = fPtr%solverData%geometry%bcLevelPointer(posInBndID)
          iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
          if (minBcID > 0) then
            if (field%bc(minBcID)%turbWallFunc%isActive) then
              res(iElem) = field%bc(minBcID)%turbWallFunc%dataOnLvl(iLevel) &
                &                           %velTau(bcLevelPointer)
            end if
          end if
        end if
      end do !iElem
    end associate

  end subroutine access_bcFricVel_forElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine returns normal distance to boundary calculated in
  !! mus_init_turbWallFunc and stored in mus_turb_wallFunc_data_type routine.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine access_bcNormDistToBnd_forElement(fun, varsys, elempos, &
    & time, tree, nElems, nDofs, res    )
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
    integer :: iElem, iLevel, iField, posInBndID, minBcID, bcLevelPointer
    type(mus_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    ! Position of dependent pdf is iField
    iField = fun%input_varPos(1)

    ! res is always AOS layout
    res = 0.0_rk
    associate (field => fPtr%solverData%scheme%field(iField))
      do iElem = 1, nElems
        ! This routine is valid only for turbulent wall boundary elements
        if ( btest( tree%elemPropertyBits( elemPos(iElem) ), &
          &         prp_hasBnd )  ) then
          ! position of current element in boundary_ID list
          posInBndID = fPtr%solverData%geometry%posInBndID( elempos(iElem) )
          minBcID = fPtr%solverData%geometry%minBcID(posInBndID)
          bcLevelPointer = fPtr%solverData%geometry%bcLevelPointer(posInBndID)
          iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
          if (minBcID > 0) then
            if (field%bc(minBcID)%turbWallFunc%isActive) then
              res(iElem) = field%bc(minBcID)%turbWallFunc%dataOnLvl(iLevel) &
                &                           %distToBnd(bcLevelPointer)
            end if
          end if
        end if
      end do !iElem
    end associate

  end subroutine access_bcNormDistToBnd_forElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine returns turbulent viscosity computed in turbulent_wall bc
  !! according to RANS formulation and stored in mus_turb_wallFunc_data_type
  !! routine.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine access_bcTurbVisc_forElement(fun, varsys, elempos, time,&
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
    integer :: iElem, iLevel, iField, posInBndID, minBcID, bcLevelPointer
    type(mus_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    ! Position of dependent pdf is iField
    iField = fun%input_varPos(1)

    ! res is always AOS layout
    res = 0.0_rk
    associate (field => fPtr%solverData%scheme%field(iField))
      do iElem = 1, nElems
        ! This routine is valid only for turbulent wall boundary elements
        if ( btest( tree%elemPropertyBits( elemPos(iElem) ), &
          &         prp_hasBnd )  ) then
          ! position of current element in boundary_ID list
          posInBndID = fPtr%solverData%geometry%posInBndID( elempos(iElem) )
          minBcID = fPtr%solverData%geometry%minBcID(posInBndID)
          bcLevelPointer = fPtr%solverData%geometry%bcLevelPointer(posInBndID)
          iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
          if (minBcID > 0) then
            if (field%bc(minBcID)%turbWallFunc%isActive) then
              res(iElem) = field%bc(minBcID)%turbWallFunc%dataOnLvl(iLevel) &
                &                           %tVisc(bcLevelPointer)
            end if
          end if
        end if
      end do !iElem
    end associate

  end subroutine access_bcTurbVisc_forElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine returns yPlus = distToBnd * fric_vel / visc
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine access_bcYPlus_forElement(fun, varsys, elempos, time, &
    &                                            tree, nElems, nDofs, res  )
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
    integer :: iElem, iLevel, iField, posInBndID, minBcID, bcLevelPointer
    real(kind=rk) :: distToBnd, velTau, viscKine
    type(mus_varSys_data_type), pointer :: fPtr
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    ! Position of dependent pdf is iField
    iField = fun%input_varPos(1)

    ! res is always AOS layout
    res = 0.0_rk
    associate (field        => fPtr%solverData%scheme%field(iField),           &
      &        levelPointer => fPtr%solverData%geometry%levelPointer,          &
      &        fluid        => fPtr%solverData%scheme%field(1)%fieldProp%fluid )
      do iElem = 1, nElems
        ! This routine is valid only for turbulent wall boundary elements
        if ( btest( tree%elemPropertyBits( elemPos(iElem) ), &
          &         prp_hasBnd )  ) then
          ! position of current element in boundary_ID list
          posInBndID = fPtr%solverData%geometry%posInBndID( elempos(iElem) )
          minBcID = fPtr%solverData%geometry%minBcID(posInBndID)
          bcLevelPointer = fPtr%solverData%geometry%bcLevelPointer(posInBndID)
          iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
          if (minBcID > 0) then
            if (field%bc(minBcID)%turbWallFunc%isActive) then
              distToBnd = field%bc(minBcID)%turbWallFunc%dataOnLvl(iLevel) &
                &                          %distToBnd(bcLevelPointer)
              velTau = field%bc(minBcID)%turbWallFunc%dataOnLvl(iLevel) &
                &                       %velTau(bcLevelPointer)
              viscKine = fluid%viscKine%dataOnLvl(iLevel)%val( &
                &          levelPointer( elemPos(iElem) ) )
              res(iElem) = distToBnd * velTau / viscKine
            end if
          end if
        end if
      end do !iElem
    end associate

  end subroutine access_bcYPlus_forElement
  ! ************************************************************************** !

end module mus_bc_var_module
