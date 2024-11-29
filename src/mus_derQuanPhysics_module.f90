! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!! This module provides the MUSUBI specific functions for calculating
!! macroscopic quantities from the state variables in physical units.
!!
!! The depending common interface between MUSUBI and ATELES is defined in the
!! [[tem_derived_module]]. The functionality for accessing a variable from the state
!! and evaluating a lua function are also provided in the tem_derived module.
!!
!! Do not use get_Element or get_Point routines to update the state !
!!
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
module mus_derQuanPhysics_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,              only: rk, labelLen
  use tem_varSys_module,       only: tem_varSys_type, tem_varSys_op_type,      &
    &                                tem_varSys_append_derVar,                 &
    &                                tem_varSys_proc_point,                    &
    &                                tem_varSys_proc_element,                  &
    &                                tem_varSys_proc_setParams,                &
    &                                tem_varSys_proc_getParams,                &
    &                                tem_varSys_proc_setupIndices,             &
    &                                tem_varSys_proc_getValOfIndex,            &
    &                                tem_varSys_getPoint_dummy,                &
    &                                tem_varSys_getElement_dummy,              &
    &                                tem_varSys_setupIndices_dummy,            &
    &                                tem_varSys_getValOfIndex_dummy,           &
    &                                tem_varSys_setParams_dummy,               &
    &                                tem_varSys_getParams_dummy
  use tem_variable_module,     only: tem_variable_type
  use tem_topology_module,     only: tem_levelOf
  use tem_time_module,         only: tem_time_type
  use treelmesh_module,        only: treelmesh_type
  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit
  use tem_dyn_array_module,    only: PositionOfVal
  use tem_grow_array_module,   only: grw_labelarray_type

  use mus_varSys_module,        only: mus_varSys_solverData_type, &
    &                                 mus_varSys_data_type,       &
    &                                 mus_get_new_solver_ptr,     &
    &                                 mus_deriveVar_ForPoint
  use mus_operation_var_module, only: mus_opVar_setupIndices

  implicit none

  private

  public :: mus_append_derVar_physics

contains


  ! **************************************************************************** !
  !> subroutine to add derive variables for weakly compressible LBM
  !! (schemekind = 'lbm') to the varsys.
  subroutine mus_append_derVar_physics( derVarName, varSys, solverData, &
    &                                   nFields, fldLabel )
    ! --------------------------------------------------------------------------
    !> array of derive variables
    type(grw_labelarray_type),  intent(in) :: derVarName

    !> global variable system
    type(tem_varSys_type), intent(inout)      :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> number of fields
    integer, intent(in)                       :: nFields

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel(:)
    ! --------------------------------------------------------------------------
    integer :: iVar, iField, addedPos, nDerVars, nFields_loc
    integer :: input_varPos, invar_nComp
    logical :: wasAdded
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
    procedure(tem_varSys_proc_setupIndices), pointer :: &
      &                                      setup_indices => null()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
      &                                       get_valOfIndex => null()
    character(len=labelLen)  ::  phyVar_name
    character(len=labelLen) ::  input_varname(1)
    ! --------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    write(logUnit(5),*) 'Append derive physical variables to varSys'

    nDerVars = derVarName%nVals

    if (nFields > 1) then
      nFields_loc = nFields + 1 ! nSpecies + 1 mixture
    else
      nFields_loc = 1
    end if

    do iField = 1, nFields_loc
      do iVar = 1, nDerVars

        ! set default pointers, overwrite if neccessary
        get_element => tem_varSys_getElement_dummy
        get_point => mus_deriveVar_ForPoint
        setup_indices => mus_opVar_setupIndices
        get_valOfIndex => tem_varSys_getValOfIndex_dummy
        set_params => tem_varSys_setParams_dummy
        get_params => tem_varSys_getParams_dummy

        select case(trim(adjustl(derVarName%val(iVar))))
        case ('density')
          get_element => deriveDensityPhy
          get_valOfIndex => deriveDensityPhy_fromIndex
        case ('charge_density', 'charge_density_boltzmann')
          get_element => deriveChargeDensityPhy
          get_valOfIndex => deriveChargeDensityPhy_fromIndex
        case ('current_density')
          get_element => deriveCurrentDensityPhy
          get_valOfIndex => deriveCurrentDensityPhy_fromIndex
        case ('pressure', 'pressure_reference', 'pressure_deviation', &
          &   'shear_stress', 'wss', 'shear_mag')
          get_element => derivePressurePhy
          get_valOfIndex => derivePressurePhy_fromIndex
        case ('kinematic_pressure')
          get_element => deriveKinePressPhy
          get_valOfIndex => deriveKinePressPhy_fromIndex
        case ('velocity', 'vel_mag', 'bc_fric_velocity')
          get_element => deriveVelocityPhy
          get_valOfIndex => derivevelocityPhy_fromIndex
        case ('momentum')
          get_element => deriveMomentumPhy
          get_valOfIndex => deriveMomentumPhy_fromIndex
        case ('bnd_force')
          get_element => deriveForcePhy
          get_valOfIndex => deriveForcePhy_fromIndex
        case ('bnd_moment')
          get_element => deriveMomentPhy
          get_valOfIndex => deriveMomentPhy_fromIndex
        case ('strain_rate', 'shear_rate')
          get_element => deriveStrainRatePhy
          get_valOfIndex => deriveStrainRatePhy_fromIndex
        case ('kinetic_energy')
          get_element => deriveKineEnerPhy
          get_valOfIndex => deriveKineEnerPhy_fromIndex
        case ('temperature')
          get_element => deriveTempPhy
          get_valOfIndex => deriveTempPhy_fromIndex
        case ('mole_density')
          get_element => deriveMoleDensityPhy
          get_valOfIndex => deriveMoleDensityPhy_fromIndex
        case ('mole_flux')
          get_element => deriveMoleFluxPhy
          get_valOfIndex => deriveMoleFluxPhy_fromIndex
        case ('moment')
          get_element => derivePDFMomentsPhy
          get_valOfIndex => derivePDFMomentsPhy_fromIndex
        case ('potential')
          get_element => derivePotentialPhy
          get_valOfIndex => derivePotentialPhy_fromIndex
        case ('electric_field')
          get_element => deriveElectric_FieldPhy
          get_valOfIndex => deriveElectric_FieldPhy_fromIndex
        case ('kine_viscosity','turb_viscosity', 'bc_turb_viscosity')
          get_element => deriveViscosityPhy
        case ('grad_velocity') ! grad velocity has same unit as strainRate
          get_element => deriveStrainRatePhy
        case ('vorticity') ! vorticity has same unit as strainRate
          get_element => deriveStrainRatePhy
        case ('q_criterion')
          get_element => deriveQCriterionPhy
        case default
          write(logUnit(7),*) 'WARNING: Variable: '//&
            &                 trim(derVarName%val(iVar))//&
            &                 ' is not defined in physical units'
          cycle !go to next variable
        end select

        if (iField > nFields) then ! mixture
          input_varname(1) = trim(derVarName%val(iVar))
        else
          input_varname(1) = trim(fldLabel(iField))//trim(derVarName%val(iVar))
        end if
        phyvar_name = trim(input_varname(1))//'_phy'
        input_varPos = PositionOfVal( me  = varSys%varname, &
          &                        val = trim(input_varname(1)) )
        if (input_varPos > 0) then
          invar_nComp = varSys%method%val(input_varPos)%nComponents

          ! append variable to varSys
          call tem_varSys_append_derVar(                            &
            &  me             = varSys,                             &
            &  varName        = trim(phyVar_name),                  &
            &  nComponents    = invar_nComp,                        &
            &  input_varname  = input_varname,                      &
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
            write(logUnit(10),*) ' Appended variable: '//trim(phyVar_name)
          else if (addedpos < 1) then
            write(logUnit(1),*) 'Error: variable '//trim(phyVar_name)// &
              &                 ' is not added to variable system'
            call tem_abort()
          end if
        end if
      end do ! iVar
    end do ! iField

  end subroutine mus_append_derVar_physics
  ! **************************************************************************** !


! **************************************************************************** !
!       Subroutines with common interface for the function pointers            !
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the density of a given set of elements (sum up all links).
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveDensityPhy(fun, varsys, elempos, time, tree, &
    &                                   nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of density in glob system
    dens_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(dens_pos)%get_element( varSys  = varSys,  &
      &                                           elemPos = elemPos, &
      &                                           time    = time,    &
      &                                           tree    = tree,    &
      &                                           nElems  = nElems,  &
      &                                           nDofs   = nDofs,   &
      &                                           res     = res      )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%rho0

  end subroutine deriveDensityPhy
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the potential of a given set of elements
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine derivePotentialPhy(fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: pot_pos, iElem
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of potential in glob system
    pot_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(pot_pos)%get_element( varSys  = varSys,  &
      &                                          elemPos = elemPos, &
      &                                          time    = time,    &
      &                                          tree    = tree,    &
      &                                          nElems  = nElems,  &
      &                                          nDofs   = nDofs,   &
      &                                          res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical unit
    do iElem = 1,nElems
      res(iElem) = res(iElem)                                           &
        &        * fPtr%solverData%physics%fac( level(iElem) )%potential
    end do

  end subroutine derivePotentialPhy
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the electric_field of a given set of elements
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveElectric_FieldPhy(fun, varsys, elempos, time, &
    &                                          tree, nElems, nDofs, res    )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: EF_pos, iElem
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of electric_field in glob system
    EF_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(EF_pos)%get_element( varSys  = varSys,  &
      &                                         elemPos = elemPos, &
      &                                         time    = time,    &
      &                                         tree    = tree,    &
      &                                         nElems  = nElems,  &
      &                                         nDofs   = nDofs,   &
      &                                         res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical unit
    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)        &
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)  &
        &        * fPtr%solverData%physics%fac( level(iElem) )%potential &
        &        / fPtr%solverData%physics%dxLvl( level(iElem) )
    end do

  end subroutine deriveElectric_FieldPhy
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the charge density of a given set of elements
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveChargeDensityPhy(fun, varsys, elempos, time, &
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: chDens_pos, iElem
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of charge density in glob system
    chDens_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(chDens_pos)%get_element( varSys  = varSys,  &
      &                                             elemPos = elemPos, &
      &                                             time    = time,    &
      &                                             tree    = tree,    &
      &                                             nElems  = nElems,  &
      &                                             nDofs   = nDofs,   &
      &                                             res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical unit
    do iElem = 1,nElems
      res(iElem) = res(iElem)                                           &
        &        * fPtr%solverData%physics%fac( level(iElem) )%chargeDens
    end do

  end subroutine deriveChargeDensityPhy
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the charge density of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveChargeDensityPhy_fromIndex(fun, varSys, time,   &
    &                                                   iLevel, idx, idxLen, &
    &                                                   nVals, res           )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: chDens_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of charge density in glob system
    chDens_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(chDens_Pos)%get_valOfIndex( &
        & varSys  = varSys,                            &
        & time    = time,                              &
        & iLevel  = iLevel,                            &
        & idx     = fPtr%opData%input_pntIndex(1)      &
        &           %indexLvl(iLevel)%val( idx(:) ),   &
        & nVals   = nVals,                             &
        & res     = res                                )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( iLevel )%chargeDens

  end subroutine deriveChargeDensityPhy_fromIndex
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the current density of a given set of elements
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveCurrentDensityPhy(fun, varsys, elempos, time, &
    &                                          tree, nElems, nDofs, res    )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: chDens_pos, iElem
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of current density in glob system
    chDens_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(chDens_pos)%get_element( varSys  = varSys,  &
      &                                             elemPos = elemPos, &
      &                                             time    = time,    &
      &                                             tree    = tree,    &
      &                                             nElems  = nElems,  &
      &                                             nDofs   = nDofs,   &
      &                                             res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical unit
    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)      &
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        & * fPtr%solverData%physics%fac( level(iElem) )%currentDens
    end do

  end subroutine deriveCurrentDensityPhy
! **************************************************************************** !


! **************************************************************************** !
  !> Compute and convert following variable pressure/shear_stress/wss/shear_mag
  !! into physical units
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine derivePressurePhy(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: iElem, input_varPos
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of density in glob system
    input_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(input_varPos)%get_element( varSys  = varSys,  &
      &                                               elemPos = elemPos, &
      &                                               time    = time,    &
      &                                               tree    = tree,    &
      &                                               nElems  = nElems,  &
      &                                               nDofs   = nDofs,   &
      &                                               res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical unit
    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)      &
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        &   * fPtr%solverData%physics%fac( level( iElem ) )%press
    end do

    deallocate( level )

  end subroutine derivePressurePhy
! **************************************************************************** !

! **************************************************************************** !
  !> Compute and convert following variable force into physical units
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveForcePhy(fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: iElem, input_varPos
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of density in glob system
    input_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(input_varPos)%get_element( varSys  = varSys,  &
      &                                               elemPos = elemPos, &
      &                                               time    = time,    &
      &                                               tree    = tree,    &
      &                                               nElems  = nElems,  &
      &                                               nDofs   = nDofs,   &
      &                                               res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical unit
    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)      &
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        &   * fPtr%solverData%physics%fac( level( iElem ) )%force
    end do

    deallocate( level )

  end subroutine deriveForcePhy
! **************************************************************************** !


! **************************************************************************** !
  !> Compute and convert following variable moment into physical units
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveMomentPhy(fun, varsys, elempos, time, tree, &
    &                                 nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: iElem, input_varPos
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of density in glob system
    input_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(input_varPos)%get_element( varSys  = varSys,  &
      &                                               elemPos = elemPos, &
      &                                               time    = time,    &
      &                                               tree    = tree,    &
      &                                               nElems  = nElems,  &
      &                                               nDofs   = nDofs,   &
      &                                               res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical unit
    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)      &
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        &   * fPtr%solverData%physics%fac( level( iElem ) )%force &
        &   * fPtr%solverData%physics%fac( level( iElem ) )%length
    end do

    deallocate( level )

  end subroutine deriveMomentPhy
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the kinematic mixture pressure in physical
  !! of a given set of elements (sum up all links).
  !!
  !! This routine requires initial mole fraction of all species which must
  !! be provided through add_variable table in scheme and set depvar in tracking
  !! for kinematic pressure to the newly added variable.
  !! This will be used to evaluate spatial function for initial mole fraction
  !! from initial mole fraction initial mixture molar density will be computed.
  !! Formula to compute kinematic pressure
  !! \( p = c^2_s (\sum_k \rho_k \phi_k - min_l (m_l) n_0)/\rho_0 \)
  !! here, \( \rho_k \) - species density, \\
  !! \( \phi_k \) - species molecular weight ratio, \\
  !! \( n_0 \) - mixture number density,\\
  !! \( \rho_0 \) - reference density.
  !! Example:
  !! In scheme table,
  !!```lua
  !! add_variable = { name = "initialMolFrac", ncomponents=1,')
  !!   initialMolFrac = { { kind="combined", spatial=luaFunc}
  !!                      { kind="combined", spatial=luaFunc2}}
  !! }
  !!```
  !! In tracking,
  !!```lua
  !! variable = {{"kinematicpressure_phy", 1, dep = { "initialMolFrac"}}}
  !!```
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveKinePressPhy(fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: iElem, press_pos
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! point to pressure derive quantities which return \rho*cs2
    press_pos = fun%input_varPos(1)

    ! derive dependent variable
    ! mixture pressure
    call varSys%method%val(press_pos)%get_element( varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            nElems  = nElems,  &
      &                                            nDofs   = nDofs,   &
      &                                            res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    !kinematic pressure (p/rho0)[m^2/s^2]
    do iElem = 1,nElems
      res(iElem) = res(iElem)                         &
        & * ( fPtr%solverData%physics%fac( level(iElem) )%length &
        &   / fPtr%solverData%physics%fac( level(iElem) )%time )**2
    end do

  end subroutine deriveKinePressPhy
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the velocity of a given element number with the given input
  !!        vector (sum up all values)
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveVelocityPhy(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer, allocatable :: level(:)
    integer :: iElem, vel_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of vel in glob system
    vel_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(vel_pos)%get_element( varSys  = varSys,  &
      &                                          elemPos = elemPos, &
      &                                          time    = time,    &
      &                                          tree    = tree,    &
      &                                          nElems  = nElems,  &
      &                                          nDofs   = nDofs,   &
      &                                          res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        &   *fPtr%solverData%physics%fac( level( iElem ) )%vel
    end do

  end subroutine deriveVelocityPhy
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the momentum of a given element number with the given input
  !!        vector (sum up all values)
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveMomentumPhy(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer, allocatable :: level(:)
    integer :: iElem, mom_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    ! position of momentum in glob system
    mom_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(mom_pos)%get_element( varSys  = varSys,  &
      &                                          elemPos = elemPos, &
      &                                          time    = time,    &
      &                                          tree    = tree,    &
      &                                          nElems  = nElems,  &
      &                                          nDofs   = nDofs,   &
      &                                          res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        &   *fPtr%solverData%physics%fac( level( iElem ) )%flux
    end do

  end subroutine deriveMomentumPhy
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the velocity magnitude of a given element number with the given
  !! input vector (sum up all values)
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveTempPhy(fun, varsys, elempos, time, tree, nElems, &
    &                                nDofs, res                                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: temp_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position in glob system
    temp_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(temp_pos)%get_element( varSys  = varSys,  &
      &                                           elemPos = elemPos, &
      &                                           time    = time,    &
      &                                           tree    = tree,    &
      &                                           nElems  = nElems,  &
      &                                           nDofs   = nDofs,   &
      &                                           res     = res      )

    res = res*fPtr%solverData%physics%temp0

  end subroutine deriveTempPhy
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the kinetic energy in physical units
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveKineEnerPhy(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer, allocatable :: level(:)
    integer :: iElem, kinec_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position in glob system
    kinec_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(kinec_pos)%get_element( varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            nElems  = nElems,  &
      &                                            nDofs   = nDofs,   &
      &                                            res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    do iElem = 1,nElems
      res(iElem) = res(iElem)*fPtr%solverData%physics%fac( level( iElem ) )%energy
    end do

  end subroutine deriveKineEnerPhy
! **************************************************************************** !


! **************************************************************************** !
  !> author: Jiaxing Qi
  !! Calculate the strain rate (or called rate of strain)
  !!
  !! The formula is:\n
  !! \[ \tau_{\alpha \beta}=
  !! -(1-\frac{\omega}{2}) \sum_{i} f^{neq}_{i} c_{i\alpha} c_{i\beta} \]\n
  !! where \( \tau_{\alpha \beta}\) is the stress
  !! in the \(\beta\)-direction on a face normal to the \(\alpha\)-axis,\n
  !! \( f^{neq}_i = f_i - f^{eq}_i\) is the non-equilibirium density.\n
  !! For more information, please refer to:\n
  !! Krueger T, Varnik F, Raabe D. Shear stress in lattice Boltzmann
  !! simulations. Physical Review E. 2009;79(4):1-14.
  !!
  recursive subroutine deriveStrainRatePhy(fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res                )
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
    ! local variables
    type(mus_varSys_data_type), pointer :: fPtr
    integer, allocatable :: level(:) ! levels of elements
    integer :: iElem, inVar_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of input variable in glob system
    inVar_pos = fun%input_varPos(1)

    ! calculate shear strain in LB unit
    call varSys%method%val(inVar_pos)%get_element( varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            nElems  = nElems,  &
      &                                            nDofs   = nDofs,   &
      &                                            res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical units
    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)&
        &   * fPtr%solverData%physics%fac( level( iElem ) )%strainRate
    end do

  end subroutine deriveStrainRatePhy
! **************************************************************************** !


! **************************************************************************** !
!                                     MULTISPECIES                             !
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the density of a given set of elements (sum up all links).
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveMoleDensityPhy(fun, varsys, elempos, time, tree, &
    &                                       nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: moleDens_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of mole density in glob system
    moleDens_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(moleDens_pos)%get_element( varSys  = varSys,  &
      &                                               elemPos = elemPos, &
      &                                               time    = time,    &
      &                                               tree    = tree,    &
      &                                               nElems  = nElems,  &
      &                                               nDofs   = nDofs,   &
      &                                               res     = res      )

    ! convert to physical unit
    res = res*fPtr%solverData%physics%moleDens0

  end subroutine deriveMoleDensityPhy
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the density of a given set of elements (sum up all links).
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveMoleFluxPhy(fun, varsys, elempos, time, tree, &
    &                                    nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: iElem, moleFlux_pos
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of moleflux in glob system
    moleFlux_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(moleFlux_pos)%get_element( varSys  = varSys,  &
      &                                               elemPos = elemPos, &
      &                                               time    = time,    &
      &                                               tree    = tree,    &
      &                                               nElems  = nElems,  &
      &                                               nDofs   = nDofs,   &
      &                                               res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    !number density physical [mol/m^3]
    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)       &
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents) &
        & * fPtr%solverData%physics%fac( level( iElem ) )%moleFlux
    end do

  end subroutine deriveMoleFluxPhy
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the momentum of a given element number with the given input
  !!        vector (sum up all values)
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine derivePDFMomentsPhy(fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer, allocatable :: level(:)
    integer :: iElem, mom_pos, offset
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    ! position of momentum in glob system
    mom_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(mom_pos)%get_element( varSys  = varSys,  &
      &                                          elemPos = elemPos, &
      &                                          time    = time,    &
      &                                          tree    = tree,    &
      &                                          nElems  = nElems,  &
      &                                          nDofs   = nDofs,   &
      &                                          res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    do iElem = 1,nElems

      offset =  (iElem-1)*fun%nComponents
      res( offset + 1) = res( offset + 1) * fPtr%solverData%physics        &
        &                                       %fac( level( iElem ) )%press
      res( offset + 2) = res( offset + 2) * fPtr%solverData%physics        &
        &                                       %fac( level( iElem ) )%press
      res( offset + 3) = res( offset + 3) * fPtr%solverData%physics        &
        &                                       %fac( level( iElem ) )%press
      res( offset + 4) = res( offset + 4) * fPtr%solverData%physics      &
        &                                       %fac( level( iElem ) )%vel
      res( offset + 5) = res( offset + 5) * fPtr%solverData%physics      &
        &                                       %fac( level( iElem ) )%vel
      res( offset + 6) = res( offset + 6) * fPtr%solverData%physics      &
        &                                       %fac( level( iElem ) )%vel
      res( offset + 7) = res( offset + 7) * fPtr%solverData%physics      &
        &                                       %fac( level( iElem ) )%vel
      res( offset + 8) = res( offset + 8) * fPtr%solverData%physics        &
        &                                       %fac( level( iElem ) )%press
      res( offset + 9) = res( offset + 9) * fPtr%solverData%physics          &
        &                                       %fac( level( iElem ) )%press &
        & * fPtr%solverData%physics%fac( level( iElem ) )%vel
    end do

  end subroutine derivePDFMomentsPhy
! **************************************************************************** !

! **************************************************************************** !
!           Subroutines with common interface for Index routines               !
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the density of a given idx (sum up all links).
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveDensityPhy_fromIndex(fun, varSys, time, iLevel, &
    &                                             idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: dens_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of density in glob system
    dens_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( dens_pos )%get_valOfIndex( &
        & varSys  = varSys,                            &
        & time    = time,                              &
        & iLevel  = iLevel,                            &
        & idx     = fPtr%opData%input_pntIndex(1)      &
        &           %indexLvl(iLevel)%val( idx(:) ),   &
        & nVals   = nVals,                             &
        & res     = res                                )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%rho0

  end subroutine deriveDensityPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the current density of a given idx (sum up all links).
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveCurrentDensityPhy_fromIndex(fun, varSys, time,   &
    &                                                    iLevel, idx, idxLen, &
    &                                                    nVals, res           )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: chdens_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of density in glob system
    chdens_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( chdens_pos )%get_valOfIndex( &
        & varSys  = varSys,                              &
        & time    = time,                                &
        & iLevel  = iLevel,                              &
        & idx     = fPtr%opData%input_pntIndex(1)        &
        &           %indexLvl(iLevel)%val( idx(:) ),     &
        & nVals   = nVals,                               &
        & res     = res                                  )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%currentDens

  end subroutine deriveCurrentDensityPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the pressure of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine derivePressurePhy_fromIndex(fun, varSys, time, iLevel, &
    &                                              idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: input_VarPos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of density in glob system
    input_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(input_varPos)%get_valOfIndex( &
        & varSys  = varSys,                              &
        & time    = time,                                &
        & iLevel  = iLevel,                              &
        & idx     = fPtr%opData%input_pntIndex(1)        &
        &           %indexLvl(iLevel)%val( idx(:) ),     &
        & nVals   = nVals,                               &
        & res     = res                                  )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%press

  end subroutine derivePressurePhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the kinematic pressure of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveKinePressPhy_fromIndex(fun, varSys, time, iLevel, &
    &                                               idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: press_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! point to pressure derive quantities which return \rho*cs2
    press_pos = fun%input_varPos(1)

    ! derive dependent variable
    ! mixture pressure
    call varSys%method%val( press_pos )%get_valOfIndex( &
        & varSys  = varSys,                             &
        & time    = time,                               &
        & iLevel  = iLevel,                             &
        & idx     = fPtr%opData%input_pntIndex(1)       &
        &           %indexLvl(iLevel)%val( idx(:) ),    &
        & nVals   = nVals,                              &
        & res     = res                                 )

    !kinematic pressure (p/rho0)[m^2/s^2]
    res = res                                                           &
        &   * ( fPtr%solverData%physics%fac( ilevel )%length            &
        &   /   fPtr%solverData%physics%fac( ilevel )%time )**2

  end subroutine deriveKinePressPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the velocity of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveVelocityPhy_fromIndex(fun, varSys, time, iLevel, &
    &                                              idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: vel_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    vel_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( vel_pos )%get_valOfIndex( &
        & varSys  = varSys,                           &
        & time    = time,                             &
        & iLevel  = iLevel,                           &
        & idx     = fPtr%opData%input_pntIndex(1)     &
        &           %indexLvl(iLevel)%val( idx(:) ),  &
        & nVals   = nVals,                            &
        & res     = res                               )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%vel

  end subroutine deriveVelocityPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the momentum of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveMomentumPhy_fromIndex(fun, varSys, time, iLevel, &
    &                                              idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: mom_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    mom_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( mom_pos )%get_valOfIndex( &
        & varSys  = varSys,                           &
        & time    = time,                             &
        & iLevel  = iLevel,                           &
        & idx     = fPtr%opData%input_pntIndex(1)     &
        &           %indexLvl(iLevel)%val( idx(:) ),  &
        & nVals   = nVals,                            &
        & res     = res                               )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%flux

  end subroutine deriveMomentumPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the bnd_force of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveForcePhy_fromIndex(fun, varSys, time, iLevel, &
    &                                           idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: input_varpos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    input_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( input_varpos )%get_valOfIndex( &
        & varSys  = varSys,                                &
        & time    = time,                                  &
        & iLevel  = iLevel,                                &
        & idx     = fPtr%opData%input_pntIndex(1)          &
        &           %indexLvl(iLevel)%val( idx(:) ),       &
        & nVals   = nVals,                                 &
        & res     = res                                    )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%force

  end subroutine deriveForcePhy_fromIndex
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the bnd_moment of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveMomentPhy_fromIndex(fun, varSys, time, iLevel, &
    &                                            idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: input_varpos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    input_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( input_varpos )%get_valOfIndex( &
        & varSys  = varSys,                                &
        & time    = time,                                  &
        & iLevel  = iLevel,                                &
        & idx     = fPtr%opData%input_pntIndex(1)          &
        &           %indexLvl(iLevel)%val( idx(:) ),       &
        & nVals   = nVals,                                 &
        & res     = res                                    )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%force &
      & * fPtr%solverData%physics%fac( ilevel )%length

  end subroutine deriveMomentPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the strain rate of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveStrainRatePhy_fromIndex(fun, varSys, time,   &
    &                                                iLevel, idx, idxLen, &
    &                                                nVals, res           )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: strain_varpos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    strain_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( strain_varpos )%get_valOfIndex( &
        & varSys  = varSys,                                 &
        & time    = time,                                   &
        & iLevel  = iLevel,                                 &
        & idx     = fPtr%opData%input_pntIndex(1)           &
        &           %indexLvl(iLevel)%val( idx(:) ),        &
        & nVals   = nVals,                                  &
        & res     = res                                     )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%strainRate

  end subroutine deriveStrainRatePhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the kinetic energy of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveKineEnerPhy_fromIndex(fun, varSys, time, iLevel, &
    &                                              idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: kinec_varpos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    kinec_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( kinec_varpos )%get_valOfIndex( &
        & varSys  = varSys,                                &
        & time    = time,                                  &
        & iLevel  = iLevel,                                &
        & idx     = fPtr%opData%input_pntIndex(1)          &
        &           %indexLvl(iLevel)%val( idx(:) ),       &
        & nVals   = nVals,                                 &
        & res     = res                                    )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%energy

  end subroutine deriveKineEnerPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the temperature of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveTempPhy_fromIndex(fun, varSys, time, iLevel, idx, &
    &                                          idxLen, nVals, res              )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: temp_varpos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    temp_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( temp_varpos )%get_valOfIndex( &
        & varSys  = varSys,                               &
        & time    = time,                                 &
        & iLevel  = iLevel,                               &
        & idx     = fPtr%opData%input_pntIndex(1)         &
        &           %indexLvl(iLevel)%val( idx(:) ),      &
        & nVals   = nVals,                                &
        & res     = res                                   )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%temp0

  end subroutine deriveTempPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the mole density of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveMoleDensityPhy_fromIndex(fun, varSys, time,   &
    &                                                 iLevel, idx, idxLen, &
    &                                                 nVals, res           )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: moleDens_varpos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    moleDens_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( moleDens_varpos )%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fPtr%opData%input_pntIndex(1)             &
        &           %indexLvl(iLevel)%val( idx(:) ),          &
        & nVals   = nVals,                                    &
        & res     = res                                       )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%moleDens0

  end subroutine deriveMoleDensityPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the mole flux of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveMoleFluxPhy_fromIndex(fun, varSys, time, iLevel, &
    &                                              idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: moleFlux_varpos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    moleFlux_varPos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( moleFlux_varpos )%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fPtr%opData%input_pntIndex(1)             &
        &           %indexLvl(iLevel)%val( idx(:) ),          &
        & nVals   = nVals,                                    &
        & res     = res                                       )

    ! convert to physical unit [mol/m^3]
    res = res * fPtr%solverData%physics%fac( ilevel )%moleFlux

  end subroutine deriveMoleFluxPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the momentum of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine derivePDFMomentsPhy_fromIndex(fun, varSys, time, iLevel,&
    &                                                idx, idxLen, nVals, res   )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: mom_pos, iVal, offset
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    mom_Pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( mom_pos )%get_valOfIndex( &
        & varSys  = varSys,                           &
        & time    = time,                             &
        & iLevel  = iLevel,                           &
        & idx     = fPtr%opData%input_pntIndex(1)     &
        &           %indexLvl(iLevel)%val( idx(:) ),  &
        & nVals   = nVals,                            &
        & res     = res                               )

    ! convert to physical unit
    do iVal = 1, nVals
      offset =  (iVal-1)*fun%nComponents
      res( offset + 1) = res( offset + 1) * fPtr%solverData%physics          &
        &                                                 %fac( ilevel )%press
      res( offset + 2) = res( offset + 2) * fPtr%solverData%physics          &
        &                                                 %fac( ilevel )%press
      res( offset + 3) = res( offset + 3) * fPtr%solverData%physics          &
        &                                                 %fac( ilevel )%press
      res( offset + 4) = res( offset + 4) * fPtr%solverData%physics          &
        &                                                   %fac( ilevel )%vel
      res( offset + 5) = res( offset + 5) * fPtr%solverData%physics          &
        &                                                   %fac( ilevel )%vel
      res( offset + 6) = res( offset + 6) * fPtr%solverData%physics          &
        &                                                   %fac( ilevel )%vel
      res( offset + 7) = res( offset + 7) * fPtr%solverData%physics          &
        &                                                   %fac( ilevel )%vel
      res( offset + 8) = res( offset + 8) * fPtr%solverData%physics          &
        &                                                 %fac( ilevel )%press
      res( offset + 9) = res( offset + 9) * fPtr%solverData%physics          &
        &                                               %fac( ilevel )%press &
        &              * fPtr%solverData%physics%fac( ilevel )%vel
    end do

  end subroutine derivePDFMomentsPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the potential of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine derivePotentialPhy_fromIndex(fun, varSys, time, iLevel, &
    &                                               idx, idxLen, nVals, res    )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: pot_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    pot_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( pot_pos )%get_valOfIndex( &
        & varSys  = varSys,                           &
        & time    = time,                             &
        & iLevel  = iLevel,                           &
        & idx     = fPtr%opData%input_pntIndex(1)     &
        &           %indexLvl(iLevel)%val( idx(:) ),  &
        & nVals   = nVals,                            &
        & res     = res                               )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%potential

  end subroutine derivePotentialPhy_fromIndex
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the electrical field of a given idx
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  recursive subroutine deriveElectric_FieldPhy_fromIndex(fun, varSys, time,   &
    &                                                    iLevel, idx, idxLen, &
    &                                                    nVals, res           )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: nVals
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: EF_pos
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )


    ! position of vel in glob system
    EF_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val( EF_pos )%get_valOfIndex( &
        & varSys  = varSys,                          &
        & time    = time,                            &
        & iLevel  = iLevel,                          &
        & idx     = fPtr%opData%input_pntIndex(1)    &
        &           %indexLvl(iLevel)%val( idx(:) ), &
        & nVals   = nVals,                           &
        & res     = res                              )

    ! convert to physical unit
    res = res * fPtr%solverData%physics%fac( ilevel )%potential &
        &     / fPtr%solverData%physics%dxLvl( ilevel )

  end subroutine deriveElectric_FieldPhy_fromIndex
! **************************************************************************** !

! **************************************************************************** !
  !> Convert lattice viscosity to physical viscosity
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveViscosityPhy(fun, varsys, elempos, time, tree, &
    &                                     nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: visc_pos, iElem
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of viscosity var in glob system
    visc_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(visc_pos)%get_element( varSys  = varSys,  &
      &                                           elemPos = elemPos, &
      &                                           time    = time,    &
      &                                           tree    = tree,    &
      &                                           nElems  = nElems,  &
      &                                           nDofs   = nDofs,   &
      &                                           res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical unit
    do iElem = 1,nElems
      res(iElem) = res(iElem) * fPtr%solverData%physics%fac(level(iElem))%visc
    end do

  end subroutine deriveViscosityPhy
! **************************************************************************** !

! **************************************************************************** !
  !> Convert lattice q-criterion to physical unit i.e strainRate**2
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveQCriterionPhy(fun, varsys, elempos, time, tree, &
    &                                      nElems, nDofs, res                )
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
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: gradU_pos, iElem
    integer, allocatable :: level(:)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! position of grad_velocity in glob system
    gradU_pos = fun%input_varPos(1)

    ! derive dependent variable
    call varSys%method%val(gradU_pos)%get_element( varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            nElems  = nElems,  &
      &                                            nDofs   = nDofs,   &
      &                                            res     = res      )

    allocate(level(nElems))
    level(1:nElems) = tem_levelOf( tree%treeID(elemPos(1:nElems)) )

    ! convert to physical unit
    do iElem = 1,nElems
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)        &
        & = res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents)  &
        &   * fPtr%solverData%physics%fac( level( iElem ) )%strainRate**2
    end do

  end subroutine deriveQCriterionPhy
! **************************************************************************** !

end module mus_derQuanPhysics_module
! **************************************************************************** !
