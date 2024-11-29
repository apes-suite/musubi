! Copyright (c) 2013, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> author: Kannan Masilmani
!! This module provides the MUSUBI specific functions for calculating
!! macroscopic quantities from the state variables multispecies.
!!
!! The depending common interface between MUSUBI and ATELES is defined in the
!! [[tem_derived_module]]. The functionality for accessing a variable from the state
!! and evaluating a lua function are also provided in the tem_derived module.
!!
!! Do not use get_Element or get_Point routines to update the state !
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
module mus_derQuanMSGas_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,               only: rk, long_k, PathLen, labelLen
  use tem_varSys_module,        only: tem_varSys_type, tem_varSys_op_type,     &
    &                                 tem_varSys_append_derVar,                &
    &                                 tem_varSys_proc_point,                   &
    &                                 tem_varSys_proc_element,                 &
    &                                 tem_varSys_proc_setParams,               &
    &                                 tem_varSys_proc_getParams,               &
    &                                 tem_varSys_proc_setupIndices,            &
    &                                 tem_varSys_proc_getValOfIndex,           &
    &                                 tem_varSys_getPoint_dummy,               &
    &                                 tem_varSys_getElement_dummy,             &
    &                                 tem_varSys_setupIndices_dummy,           &
    &                                 tem_varSys_getValOfIndex_dummy,          &
    &                                 tem_varSys_setParams_dummy,              &
    &                                 tem_varSys_getParams_dummy
  use tem_variable_module,      only: tem_variable_type
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_param_module,         only: cs2, cs2inv, t2cs4inv, t2cs2inv
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_subTree_type_module,  only: tem_subTree_type, tem_treeIDfrom_subTree
  use tem_logging_module,       only: logUnit
  use tem_aux_module,           only: tem_abort
  use tem_math_module,          only: invert_matrix
  use tem_topology_module,      only: tem_levelOf
  use tem_operation_var_module, only: tem_evalMag_forElement,     &
    &                                 tem_evalMag_forPoint,       &
    &                                 tem_evalMag_fromIndex,      &
    &                                 tem_opVar_setupIndices,     &
    &                                 tem_get_new_varSys_data_ptr
  use tem_grow_array_module,    only: grw_labelarray_type, append

  use mus_varSys_module,          only: mus_varSys_solverData_type, &
    &                                   mus_varSys_data_type,       &
    &                                   mus_get_new_solver_ptr,     &
    &                                   mus_deriveVar_ForPoint
  use mus_scheme_header_module,   only: mus_scheme_header_type
  use mus_scheme_layout_module,   only: mus_scheme_layout_type
  use mus_scheme_type_module,     only: mus_scheme_type
  use mus_derQuanMSLiquid_module, only: deriveMoleDensityMS,           &
    &                                   derivePressureMS,              &
    &                                   deriveMoleFluxMS,              &
    &                                   deriveMoleFracMS,              &
    &                                   deriveMassFracMS,              &
    &                                   deriveVelocityMS,              &
    &                                   deriveMoleDensityMS_fromIndex, &
    &                                   deriveMoleFluxMS_fromIndex,    &
    &                                   deriveVelocityMS_fromIndex,    &
    &                                   mus_append_derMixVar_MS
  use mus_operation_var_module,   only: mus_opVar_setupIndices
  use mus_derVarPos_module,       only: mus_derVarPos_type
  use mus_physics_module,            only: mus_convertFac_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  ! include aotus modules
  use aotus_module, only: flu_State

  implicit none

  private

  ! functions for multispecies gas mixture
  public :: mus_append_derVar_MSGas
  public :: deriveAuxMSGas_fromState
  public :: deriveEquilMSGas_fromAux
  public :: deriveEquilMSGas_FromMacro
  public :: deriveVelMSGas_FromState
  public :: deriveMomMSGas_FromState
  public :: deriveEqMSGas_FromState
  public :: deriveMomentaMSGas_FromState
  public :: deriveVelocitiesMSGas_FromState

contains

  ! **************************************************************************** !
  !> subroutine to add derive variables for multispecies-liquid
  !! (schemekind = 'multispecies_gas') to the varsys.
  subroutine mus_append_derVar_MSGas( varSys, solverData, stencil,  &
    &                                 nFields, fldLabel, derVarName )
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)  :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> compute stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil

    !> number of fields
    integer, intent(in)                       :: nFields

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel(:)

    !> array of derive physical variables
    type(grw_labelarray_type), intent(inout) :: derVarName
    ! ---------------------------------------------------------------------------
    ! number of derive variables
    integer :: nDerVars, iVar, nComponents, addedPos, iIn
    integer :: iField
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
    character(len=labelLen), allocatable :: derVarName_loc(:)
    ! ---------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    nDerVars = 9
    allocate(derVarName_loc(nDerVars))
    derVarName_loc    = [ 'pressure       ', 'equilibrium    ', &
      &                   'equilibrium_vel', 'mole_density   ', &
      &                   'mole_fraction  ', 'mass_fraction  ', &
      &                   'mole_flux      ', 'velocity       ', &
      &                   'vel_mag        '  ]

    do iVar = 1, nDerVars
      call append(derVarName, derVarName_loc(iVar))
    end do

    ! add variable for each field and add mixture variable later
    do iField = 1, nFields
      do iVar = 1, nDerVars

        ! set default pointers, overwrite if neccessary
        get_element => tem_varSys_getElement_dummy
        get_point => mus_deriveVar_ForPoint
        setup_indices => mus_opVar_setupIndices
        get_valOfIndex => tem_varSys_getValOfIndex_dummy
        method_data  = mus_get_new_solver_ptr(solverData)
        set_params => tem_varSys_setParams_dummy
        get_params => tem_varSys_getParams_dummy

        select case(trim(adjustl(derVarName_loc(iVar))))
        case ('mole_density')
          get_element => deriveMoleDensityMS
          get_valOfIndex => deriveMoleDensityMS_fromIndex
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'density'

        case ('mole_fraction')
          get_element => deriveMoleFracMS
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'density'

        case ('mass_fraction')
          get_element => deriveMassFracMS
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'density'

        case ('velocity')
          get_element => deriveVelocityMS
          get_valOfIndex => deriveVelocityMS_fromIndex
          nComponents = 3
          allocate(input_varname(2))
          input_varname(1) = 'density'
          input_varname(2) = 'momentum'

        case ('mole_flux')
          get_element => deriveMoleFluxMS
          get_valOfIndex => deriveMoleFluxMS_fromIndex
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'momentum'

        case ('pressure')
          get_element => derivePressureMS
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'density'

        case ('equilibrium_vel')
          get_element => deriveEquilVelMSGas
          nComponents = 3
          allocate(input_varname(1))
          input_varname(1) = 'pdf'

        case ('equilibrium')
          get_element => deriveEquilMSGas
          nComponents = stencil%QQ
          allocate(input_varname(1))
          input_varname(1) = 'pdf'

        case ('vel_mag')
          get_element => tem_evalMag_forElement
          get_point => tem_evalMag_forPoint
          get_valOfIndex => tem_evalMag_fromIndex
          setup_indices => tem_opVar_setupIndices
          method_data = tem_get_new_varSys_data_ptr(method_data)
          nComponents = 1
          allocate(input_varname(1))
          input_varname(1) = 'velocity'

        case default
          write(logUnit(1),*) 'WARNING: Unknown variable: '//&
            &                 trim(derVarName_loc(iVar))
          cycle !go to next variable
        end select

        ! update variable names with field label
        varname = trim(fldLabel(iField))//trim(adjustl(derVarName_loc(iVar)))
        do iIn = 1, size(input_varname)
          input_varname(iIn) = trim(fldLabel(iField))//trim(input_varname(iIn))
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
          write(logUnit(10),*) ' Appended variable: '//trim(varname)
        else if (addedpos < 1) then
          write(logUnit(1),*) 'Error: variable '//trim(varname)// &
            &                 ' is not added to variable system'
        end if

        deallocate(input_varname)
      end do !iVar
    end do !iField

    ! append mixture variable for every derived species variable
    call mus_append_derMixVar_MS( varSys     = varSys,     &
      &                           solverData = solverData, &
      &                           nFields    = nFields,    &
      &                           fldLabel   = fldLabel,   &
      &                           derVarName = derVarName  )

  end subroutine mus_append_derVar_MSGas
  ! **************************************************************************** !


! ****************************************************************************** !
!       Subroutines with common interface for the function pointers            !
! ****************************************************************************** !

! ****************************************************************************** !
  !> Equilibrium velocity from state
  !! Calculate the momentum of a given element for single species or mixture
  !! from the cptr scheme state vector for gas mixture (Asinari model).
  !! Need to solve the system of equations to compute this momentum since
  !! first order moments gives only moments of transformed pdfs. Hence
  !! to obtain the first order moments of actual pdfs we need to solve
  !! system of equation of size = nSpecies for each velocity components
  !! @todo Add equation and reference
  recursive subroutine deriveEquilVelMSGas(fun, varsys, elempos, time, tree, &
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
    integer :: statePos, iElem, iComp, iLevel, iField, ifld, nFields
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: pdfPos, varPos, QQ
    real(kind=rk), allocatable :: tmpPDF(:)
    !mass fraction of nSpecies
    real(kind=rk) :: massFraction( varSys%nStateVars )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !partial pressure of  nSpecies
    real(kind=rk) :: press( varSys%nStateVars )
    real(kind=rk) :: lambda( varSys%nStateVars )
    real(kind=rk) :: pressMix, densMix, molWeightMix
    !first moments of nSpecies
    real(kind=rk) :: first_moments( fun%nComponents, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( fun%nComponents, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: chi( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    integer :: iFieldDia, iFieldNonDia
    real(kind=rk) :: matrixA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: invA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: eqVel(3)
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    integer :: nSize
    integer :: stateVarMap(varSys%nStateVars)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = scheme%layout%fStencil%QQ
    allocate(tmpPDF(QQ))

    ! equilibrium velocity can be computed only for single species
    pdfPos = fun%input_varPos(1)
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iField = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iField) = scheme%field(iField)%fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iField) = scheme%field(iField)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iField, :) =                                &
        & scheme%field(iField)%fieldProp%species%resi_coeff(:)
    end do

    do iElem = 1, nElems
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      statePos = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
      nSize = scheme%pdf( iLevel )%nSize

      do iField = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iField))%state_varPos(iComp)

          tmpPDF(iComp) = scheme%state( iLevel )%val(              &
            ! position of this state variable in the state array
            & ( statepos-1)* varsys%nscalars+varpos,    &
            & scheme%pdf( iLevel )%nNext )
        end do

        ! mass density of species
        mass_dens(iField ) = sum( tmpPDF )

        ! partial pressure of species
        press(iField) = mass_dens(iField) * phi(iField) * cs2

        ! velocity, first moments
        do iComp = 1, fun%nComponents
          first_moments(iComp, iField) = sum( tmpPDF *        &
            & scheme%layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iField

      ! total pressure
      pressMix = sum(press)

      ! total density
      densMix = sum(mass_dens)

      ! mass freaction
      massFraction(:) = mass_dens(:) / densMix

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iField = 1, nFields
        molWeightMix = molWeightMix + massFraction(iField)/molWeight(iField)
      end do
      molWeightMix = 1.0_rk/molWeightMix

      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      do iField = 1, nFields
        do iFieldDia = 1, nFields
          chi(iField,iFieldDia) = ( molWeightMix*molWeightMix           &
            & / (molWeight(iField)*molWeight(iFieldDia)) )              &
            & * (resi_coeff(iField, iFieldDia)/resi_coeff(iField,iField))
        enddo
      enddo

      !relaxation time
      do iField = 1, nFields
        lambda(iField) = pressMix*resi_coeff(iField,iField)/densMix
      enddo

      ! build up the equation system for momentum
      matrixA = 0.0_rk
      do iField = 1, nFields
        ! set diagonal part
        matrixA( iField, iField ) = 1.0_rk
        do iFieldDia = 1, nFields
          matrixA( iField, iField ) = matrixA( iField, iField )     &
            & + lambda(iField) * 0.5_rk * chi(iField, iFieldDia)    &
            & * massFraction( iFieldDia )
        end do
        ! set nonDiagonal
        do iFieldNonDia = 1,nFields
          matrixA( iField, iFieldNonDia ) = matrixA( iField, iFieldNonDia ) &
            & - lambda(iField) * 0.5_rk * chi( iField, iFieldNonDia )       &
            & * massFraction( iField )
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matrixA )

      ! momentum of all species
      momentum = 0.0_rk
      do iComp = 1, fun%nComponents
        momentum( iComp, : ) = matmul( invA, first_moments( iComp, : ) )
      end do

      !velocity of all species
      do iField = 1,nFields
        vel( :, iField) = momentum( :, iField) / mass_dens(iField)
      end do

      eqVel(:) = vel(:, pdfPos)
      do ifld = 1, nFields
        eqVel(:) = eqVel(:) + chi(ifld, pdfPos)                          &
          &        * massfraction(ifld) * ( vel(:,ifld) - vel(:, pdfPos) )
      end do

      ! copy the results to the res
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents) = eqVel

    end do !iElem

  end subroutine deriveEquilVelMSGas
! ****************************************************************************** !

! ****************************************************************************** !
  !> Calculate the equlibrium of a given element number with the given input
  !! state vector.
  !!
  !! The equilibrium distribution function is:\n
  !! \[ f^{eq}_rest = \rho w_i( (3-2 mRatio_iField ) + 3(cx . u)
  !!            + 9/2 (cx . u )^2 - 3/2 u^2 ) \]
  !! \[ f^{eq}_i = \rho w_i( mRatio_iField + 3(cx . u)
  !!            + 9/2 (cx . u )^2 - 3/2 u^2 ) \]
  !! where \(w_i\) is the weight in each direction,\n
  !! \(\rho\) is the macroscopic value of density,\n
  !! \(c_s\) is the speed of sound,\n
  !! \(\vec c_i\) is the lattice unit velocity in each direction,\n
  !! \(\vec u\) is the macroscopic value of velocity.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine deriveEquilMSGas(fun, varsys, elempos, time, tree, &
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
    integer :: statePos, iElem, iComp, iLevel, iField, iDir, nFields, iFld
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: pdfPos, varPos, QQ
    real(kind=rk), allocatable :: tmpPDF(:)
    !mass fraction of nSpecies
    real(kind=rk) :: massFraction( varSys%nStateVars )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !partial pressure of  nSpecies
    real(kind=rk) :: press( varSys%nStateVars )
    real(kind=rk) :: lambda( varSys%nStateVars )
    real(kind=rk) :: pressMix, densMix, molWeightMix
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: chi( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    integer :: iFieldDia, iFieldNonDia
    real(kind=rk) :: matrixA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: invA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: fEq(fun%nComponents)
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    real(kind=rk) :: eqVel(3)
    real(kind=rk) :: ucx, usq, weight0Inv
    integer :: nSize
    integer :: stateVarMap(varSys%nStateVars)
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = scheme%layout%fStencil%QQ

    ! equilibrium velocity can be computed only for single species
    pdfPos = fun%input_varPos(1)

    stateVarMap = scheme%stateVarMap%varPos%val(:)
    do iField = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iField) = scheme%field(iField)%fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iField) = scheme%field(iField)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iField, :) =                                &
        & scheme%field(iField)%fieldProp%species%resi_coeff(:)
    end do

    weight0Inv = 1.0_rk / scheme%layout%weight( &
      &                              scheme%layout%fStencil%restPosition)
    do iElem = 1, nElems
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      statePos = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
      nSize = scheme%pdf( iLevel )%nSize

      do iField = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iField))%state_varPos(iComp)

          tmpPDF(iComp) = scheme%state( iLevel )%val(              &
            ! position of this state variable in the state array
            & ( statepos-1)* varsys%nscalars+varpos,    &
            & scheme%pdf( iLevel )%nNext )
        end do

        ! mass density of species
        mass_dens(iField ) = sum( tmpPDF )

        ! partial pressure of species
        press(iField) = mass_dens(iField) * phi(iField) * cs2

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iField) = sum( tmpPDF * &
            & scheme%layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iField

      ! total pressure
      pressMix = sum(press)

      ! total density
      densMix = sum(mass_dens)

      ! mass freaction
      massFraction(:) = mass_dens(:) / densMix

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iField = 1, nFields
        molWeightMix = molWeightMix + massFraction(iField)/molWeight(iField)
      end do
      molWeightMix = 1.0_rk/molWeightMix

      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      do iField = 1, nFields
        do iFieldDia = 1, nFields
          chi(iField,iFieldDia) = ( molWeightMix*molWeightMix           &
            & / (molWeight(iField)*molWeight(iFieldDia)) )              &
            & * (resi_coeff(iField, iFieldDia)/resi_coeff(iField,iField))
        enddo
      enddo

      !relaxation time
      do iField = 1, nFields
        lambda(iField) = pressMix*resi_coeff(iField,iField)/densMix
      enddo

      ! build up the equation system for momentum
      matrixA = 0.0_rk
      do iField = 1, nFields
        ! set diagonal part
        matrixA( iField, iField ) = 1.0_rk
        do iFieldDia = 1, nFields
          matrixA( iField, iField ) = matrixA( iField, iField )     &
            & + lambda(iField) * 0.5_rk * chi(iField, iFieldDia)    &
            & * massFraction( iFieldDia )
        end do
        ! set nonDiagonal
        do iFieldNonDia = 1,nFields
          matrixA( iField, iFieldNonDia ) = matrixA( iField, iFieldNonDia ) &
            & - lambda(iField) * 0.5_rk * chi( iField, iFieldNonDia )       &
            & * massFraction( iField )
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matrixA )

      ! momentum of all species
      momentum = 0.0_rk
      do iComp = 1, 3
        momentum( iComp, : ) = matmul( invA, first_moments( iComp, : ) )
      end do

      !velocity of all species
      do iField = 1,nFields
        vel( :, iField) = momentum( :, iField) / mass_dens(iField)
      end do

      eqVel(:) = vel(:,pdfPos)
      do ifld = 1, nFields
        eqVel(:) = eqVel(:) + chi(ifld, pdfPos)                           &
          &          * massfraction(ifld) * ( vel(:,ifld) - vel(:,pdfPos) )
      end do

      ! Calculate the square of velocity
      usq = dot_product( eqVel,eqVel ) * t2cs2inv

      do iDir = 1, scheme%layout%fStencil%QQ
        ! Velocity times lattice unit velocity
        ucx = dot_product( scheme%layout%fStencil%cxDirRK(:, iDir), eqVel )

        ! calculate equilibrium
        fEq(iDir) = scheme%layout%weight( iDir )                          &
          &         * mass_dens(pdfPos) * ( phi(pdfPos)  + ucx * cs2inv &
          &         + ucx * ucx * t2cs4inv - usq )
      end do ! iDir

      fEq(scheme%layout%stencil%restPosition)                           &
        & = scheme%layout%weight( scheme%layout%fStencil%restPosition ) &
        &   * mass_dens(pdfPos)                                         &
        &   * ( weight0Inv + (1.0_rk - weight0Inv)*phi(pdfPos) - usq )

      ! copy the results to the res
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents) = fEq

    end do !iElem

  end subroutine deriveEquilMSGas
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine computes equilbrium from density and velocity
  !! This must comply with mus_variable_module%derive_FromMacro
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromMacro]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:equilFromMacro]] function pointer.
  subroutine deriveEquilMSGas_FromMacro( density, velocity, iField, nElems, &
    &                                    varSys, layout, res                )
    ! -------------------------------------------------------------------- !
    !> Array of density.
    !! Single species: dens_1, dens_2 .. dens_n
    !! multi-species: dens_1_sp1, dens_1_sp2, dens_2_sp1, dens_2_sp2 ...
    !!                dens_n_sp1, dens_n_sp2
    real(kind=rk), intent(in) :: density(:)

    !> Array of velocity.
    !! Size: dimension 1: n*nFields. dimension 2: 3 (nComp)
    !! 1st dimension arrangement for multi-species is same as density
    real(kind=rk), intent(in) :: velocity(:, :)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n*nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: fEq(layout%fStencil%QQ)
    integer :: QQ, iElem, iDir, iFld, nFields, offset
    type(mus_scheme_type), pointer :: scheme
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: resi_coeff(varSys%nStateVars), phi
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !mass fraction
    real(kind=rk) :: massFraction( varSys%nStateVars )
    real(kind=rk) :: totMass_densInv
    real(kind=rk) :: molWeightInv( varSys%nStateVars )
    real(kind=rk) :: vel( 3, varSys%nStateVars ), eqVel(3)
    real(kind=rk) :: ucx, usq, weight0Inv, molWeightMix
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    QQ = layout%fStencil%QQ
    nFields = scheme%nFields
    ! resistivity coefficients
    resi_coeff(:) =                                     &
      & scheme%field(iField)%fieldProp%species%resi_coeff(:)

    weight0Inv = 1.0_rk / scheme%layout%weight(                     &
      &                          scheme%layout%fStencil%restPosition)

    phi = scheme%field(iField)%fieldProp%species%molWeigRatio

    do ifld = 1, nFields
      ! molecular weight inverse
      molWeightInv(ifld) = scheme%field(ifld)%fieldProp%species%molWeightInv
    end do

    do iElem = 1, nElems
      offset = (iElem-1)*nFields
      ! get species density and velocity of iElem
      do ifld = 1, nFields
        mass_dens(ifld) = density( offset+ifld )
        vel(:, ifld) = velocity(:,offset+ifld)
      end do
      totmass_densInv = 1.0_rk/sum(mass_dens)
      ! massfraction
      massFraction = mass_dens * totmass_densInv

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iFld = 1, nFields
        molWeightMix = molWeightMix &
          &          + massFraction(iFld) * molWeightInv(iFld)
      end do
      molWeightMix = 1.0_rk/molWeightMix

      eqVel(:) = vel(:, iField)
      do ifld = 1, nFields
        eqVel(:) = eqVel(:) + molWeightMix * molWeightMix  &
        &        * molWeightInv(iField)*molWeightInv(ifld) &
        &        * resi_coeff(ifld) * massFraction(ifld)   &
        &        * ( vel(:,ifld) - vel(:,iField) )         &
        &        / resi_coeff(iField)
      end do

      ! Calculate the square of velocity
      usq = dot_product( eqVel, eqVel ) * t2cs2inv

      do iDir = 1, QQ
        ! Velocity times lattice unit velocity
        ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), eqVel )

        ! calculate equilibrium
        fEq(iDir) = layout%weight( iDir ) * mass_dens(iField)     &
          &         * ( phi + ucx * cs2inv + ucx * ucx * t2cs4inv &
          &         - usq                                         )
      end do ! iDir

      fEq(layout%stencil%restPosition) =                   &
        & layout%weight( layout%fStencil%restPosition )    &
        & * mass_dens(iField)                              &
        & * ( weight0Inv + (1.0_rk - weight0Inv)*phi - usq )

      res( (iElem-1)*QQ+1: iElem*QQ ) = fEq
    end do

  end subroutine deriveEquilMSGas_FromMacro
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine computes velocity from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveVelMSGas_FromState( state, iField, nElems, varSys, layout, &
    &                                  res                                    )
    ! -------------------------------------------------------------------- !
    !> Array of state
    !! n * layout%fStencil%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iComp, nFields
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos, QQ
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    real(kind=rk) :: vel(3)
    !mass fraction of nSpecies
    real(kind=rk) :: massFraction( varSys%nStateVars )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !partial pressure of  nSpecies
    real(kind=rk) :: press( varSys%nStateVars )
    real(kind=rk) :: lambda( varSys%nStateVars )
    real(kind=rk) :: pressMix, densMix, molWeightMix
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: chi( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    integer :: iFld, iFldDia, iFldNonDia
    real(kind=rk) :: matrixA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: invA( varSys%nStateVars, varSys%nStateVars )
    integer :: stateVarMap(varSys%nStateVars)
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = layout%fStencil%QQ
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFld) = scheme%field(iFld)%fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    do iElem = 1, nElems

      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iFld))%state_varPos(iComp)
          tmpPDF(iComp) = state(( ielem-1)* varsys%nscalars+varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! partial pressure of species
        press(iFld) = mass_dens(iFld) * phi(iFld) * cs2

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF * &
            & layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total pressure
      pressMix = sum(press)

      ! total density
      densMix = sum(mass_dens)

      ! mass freaction
      massFraction(:) = mass_dens(:) / densMix

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iFld = 1, nFields
        molWeightMix = molWeightMix + massFraction(iFld)/molWeight(iFld)
      end do
      molWeightMix = 1.0_rk / molWeightMix

      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      do iFld = 1, nFields
        do iFldDia = 1, nFields
          chi(iFld,iFldDia) = ( molWeightMix*molWeightMix       &
            & / (molWeight(iFld)*molWeight(iFldDia)) )          &
            & * (resi_coeff(iFld, iFldDia)/resi_coeff(iFld,iFld))
        enddo
      enddo

      !relaxation time
      do iFld = 1, nFields
        lambda(iFld) = pressMix*resi_coeff(iFld,iFld)/densMix
      enddo

      ! build up the equation system for momentum
      matrixA = 0.0_rk
      do iFld = 1, nFields
        ! set diagonal part
        matrixA( iFld, iFld ) = 1.0_rk
        do iFldDia = 1, nFields
          matrixA( iFld, iFld ) = matrixA( iFld, iFld )    &
            & + lambda(iFld) * 0.5_rk * chi(iFld, iFldDia) &
            & * massFraction( iFldDia )
        end do
        ! set nonDiagonal
        do iFldNonDia = 1,nFields
          matrixA( iFld, iFldNonDia ) = matrixA( iFld, iFldNonDia ) &
            & - lambda(iFld) * 0.5_rk * chi( iFld, iFldNonDia )     &
            & * massFraction( iFld )
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matrixA )

      ! momentum of all species
      momentum = 0.0_rk
      do iComp = 1, 3
        momentum( iComp, : ) = matmul( invA, first_moments( iComp, : ) )
      end do

      ! find required species velocity
      vel = momentum(:, iField)/mass_dens(iField)

      ! copy the results to the res
      res( (iElem-1)*3 + 1 : iElem*3) = vel
    end do

  end subroutine deriveVelMSGas_FromState
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine computes momentum from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveMomMSGas_FromState( state, iField, nElems, varSys, layout, &
    &                                  res                                    )
    ! -------------------------------------------------------------------- !
    !> Array of state
    !! n * layout%fStencil%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iComp, nFields
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos, QQ
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    !mass fraction of nSpecies
    real(kind=rk) :: massFraction( varSys%nStateVars )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    real(kind=rk) :: press( varSys%nStateVars )
    real(kind=rk) :: lambda( varSys%nStateVars )
    real(kind=rk) :: pressMix, densMix, molWeightMix
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: chi( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    integer :: iFld, iFldDia, iFldNonDia
    real(kind=rk) :: matrixA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: invA( varSys%nStateVars, varSys%nStateVars )
    integer :: stateVarMap(varSys%nStateVars)
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = layout%fStencil%QQ
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFld) = scheme%field(iFld)%fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    do iElem = 1, nElems
      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iFld))%state_varPos(iComp)

          tmpPDF(iComp) = state(                       &
            ! position of this state variable in the state array
            & ( ielem-1)* varsys%nscalars+varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! partial pressure of species
        press(iFld) = mass_dens(iFld) * phi(iFld) * cs2

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF * &
            & layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total pressure
      pressMix = sum(press)

      ! total density
      densMix = sum(mass_dens)

      ! mass freaction
      massFraction(:) = mass_dens(:) / densMix

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iFld = 1, nFields
        molWeightMix = molWeightMix + massFraction(iFld)/molWeight(iFld)
      end do
      molWeightMix = 1.0_rk/molWeightMix

      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      do iFld = 1, nFields
        do iFldDia = 1, nFields
          chi(iFld,iFldDia) = ( molWeightMix*molWeightMix       &
            & / (molWeight(iFld)*molWeight(iFldDia)) )          &
            & * (resi_coeff(iFld, iFldDia)/resi_coeff(iFld,iFld))
        enddo
      enddo

      !relaxation time
      do iFld = 1, nFields
        lambda(iFld) = pressMix*resi_coeff(iFld,iFld)/densMix
      enddo

      ! build up the equation system for momentum
      matrixA = 0.0_rk
      do iFld = 1, nFields
        ! set diagonal part
        matrixA( iFld, iFld ) = 1.0_rk
        do iFldDia = 1, nFields
          matrixA( iFld, iFld ) = matrixA( iFld, iFld )    &
            & + lambda(iFld) * 0.5_rk * chi(iFld, iFldDia) &
            & * massFraction( iFldDia )
        end do
        ! set nonDiagonal
        do iFldNonDia = 1,nFields
          matrixA( iFld, iFldNonDia ) = matrixA( iFld, iFldNonDia ) &
            & - lambda(iFld) * 0.5_rk * chi( iFld, iFldNonDia )     &
            & * massFraction( iFld )
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matrixA )

      ! momentum of all species
      momentum = 0.0_rk
      do iComp = 1, 3
        momentum( iComp, : ) = matmul( invA, first_moments( iComp, : ) )
      end do

      ! copy the results to the res
      res( (iElem-1)*3 + 1 : iElem*3) = momentum(:, iField)
    end do

  end subroutine deriveMomMSGas_FromState
! ***************************************************************************** !


! ****************************************************************************** !
  !> This routine computes velocity from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveVelocitiesMSGas_FromState( state, iField, nElems, varSys, &
    &                                         layout, res                    )
    ! -------------------------------------------------------------------- !
    !> Array of state
    !! n * layout%fStencil%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iComp, nFields
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos, QQ
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    real(kind=rk) :: vel(3)
    !mass fraction of nSpecies
    real(kind=rk) :: massFraction( varSys%nStateVars )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !partial pressure of  nSpecies
    !partial pressure of  nSpecies
    real(kind=rk) :: press( varSys%nStateVars )
    real(kind=rk) :: lambda( varSys%nStateVars )
    real(kind=rk) :: pressMix, densMix, molWeightMix
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: chi( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    integer :: iFld, iFldDia, iFldNonDia
    real(kind=rk) :: matrixA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: invA( varSys%nStateVars, varSys%nStateVars )
    integer :: stateVarMap(varSys%nStateVars)
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = layout%fStencil%QQ
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFld) = scheme%field(iFld)%fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    do iElem = 1, nElems
      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iFld))%state_varPos(iComp)

          tmpPDF(iComp) = state(                     &
            ! position of this state variable in the state array
            & ( ielem-1)* varsys%nscalars+varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! partial pressure of species
        press(iFld) = mass_dens(iFld) * phi(iFld) * cs2

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF * &
            & layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total pressure
      pressMix = sum(press)

      ! total density
      densMix = sum(mass_dens)

      ! mass freaction
      massFraction(:) = mass_dens(:) / densMix

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iFld = 1, nFields
        molWeightMix = molWeightMix + massFraction(iFld)/molWeight(iFld)
      end do
      molWeightMix = 1.0_rk/molWeightMix

      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      do iFld = 1, nFields
        do iFldDia = 1, nFields
          chi(iFld,iFldDia) = ( molWeightMix*molWeightMix       &
            & / (molWeight(iFld)*molWeight(iFldDia)) )          &
            & * (resi_coeff(iFld, iFldDia)/resi_coeff(iFld,iFld))
        enddo
      enddo

      !relaxation time
      do iFld = 1, nFields
        lambda(iFld) = pressMix*resi_coeff(iFld,iFld)/densMix
      enddo

      ! build up the equation system for momentum
      matrixA = 0.0_rk
      do iFld = 1, nFields
        ! set diagonal part
        matrixA( iFld, iFld ) = 1.0_rk
        do iFldDia = 1, nFields
          matrixA( iFld, iFld ) = matrixA( iFld, iFld )     &
            & + lambda(iFld) * 0.5_rk * chi(iFld, iFldDia)  &
            & * massFraction( iFldDia )
        end do
        ! set nonDiagonal
        do iFldNonDia = 1,nFields
          matrixA( iFld, iFldNonDia ) = matrixA( iFld, iFldNonDia ) &
            & - lambda(iFld) * 0.5_rk * chi( iFld, iFldNonDia )     &
            & * massFraction( iFld )
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matrixA )

      ! momentum of all species
      momentum = 0.0_rk
      do iComp = 1, 3
        momentum( iComp, : ) = matmul( invA, first_moments( iComp, : ) )
      end do


      ! copy the results to the res
      do iFld = 1, nFields
        ! find required species velocity
        vel = momentum(:, iFld)/mass_dens(iFld)
        res( (iElem-1)*nFields*3 + (iField-1)*3 + 1) = vel(1)
        res( (iElem-1)*nFields*3 + (iField-1)*3 + 2) = vel(2)
        res( (iElem-1)*nFields*3 + (iField-1)*3 + 3) = vel(3)
      end do

    end do

  end subroutine deriveVelocitiesMSGas_FromState
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine computes momentum from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveMomentaMSGas_FromState( state, iField, nElems, varSys, &
    &                                      layout, res                    )
    ! -------------------------------------------------------------------- !
    !> Array of state
    !! n * layout%fStencil%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iComp, nFields
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos, QQ
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    !mass fraction of nSpecies
    real(kind=rk) :: massFraction( varSys%nStateVars )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !partial pressure of  nSpecies
    real(kind=rk) :: press( varSys%nStateVars )
    real(kind=rk) :: lambda( varSys%nStateVars )
    real(kind=rk) :: pressMix, densMix, molWeightMix
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: chi( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    integer :: iFld, iFldDia, iFldNonDia
    real(kind=rk) :: matrixA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: invA( varSys%nStateVars, varSys%nStateVars )
    integer :: stateVarMap(varSys%nStateVars)
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = layout%fStencil%QQ
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFld) = scheme%field(iFld)%fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                     &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    do iElem = 1, nElems
      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iFld))%state_varPos(iComp)

          tmpPDF(iComp) = state(                       &
            ! position of this state variable in the state array
            & ( ielem-1)* varsys%nscalars+varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! partial pressure of species
        press(iFld) = mass_dens(iFld) * phi(iFld) * cs2

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF * &
            & layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total pressure
      pressMix = sum(press)

      ! total density
      densMix = sum(mass_dens)

      ! mass freaction
      massFraction(:) = mass_dens(:) / densMix

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iFld = 1, nFields
        molWeightMix = molWeightMix + massFraction(iFld)/molWeight(iFld)
      end do
      molWeightMix = 1.0_rk/molWeightMix

      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      do iFld = 1, nFields
        do iFldDia = 1, nFields
          chi(iFld,iFldDia) = ( molWeightMix*molWeightMix       &
            & / (molWeight(iFld)*molWeight(iFldDia)) )          &
            & * (resi_coeff(iFld, iFldDia)/resi_coeff(iFld,iFld))
        enddo
      enddo

      !relaxation time
      do iFld = 1, nFields
        lambda(iFld) = pressMix*resi_coeff(iFld,iFld)/densMix
      enddo

      ! build up the equation system for momentum
      matrixA = 0.0_rk
      do iFld = 1, nFields
        ! set diagonal part
        matrixA( iFld, iFld ) = 1.0_rk
        do iFldDia = 1, nFields
          matrixA( iFld, iFld ) = matrixA( iFld, iFld )    &
            & + lambda(iFld) * 0.5_rk * chi(iFld, iFldDia) &
            & * massFraction( iFldDia )
        end do
        ! set nonDiagonal
        do iFldNonDia = 1,nFields
          matrixA( iFld, iFldNonDia ) = matrixA( iFld, iFldNonDia ) &
            & - lambda(iFld) * 0.5_rk * chi( iFld, iFldNonDia )     &
            & * massFraction( iFld )
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matrixA )

      ! momentum of all species
      momentum = 0.0_rk
      do iComp = 1, 3
        momentum( iComp, : ) = matmul( invA, first_moments( iComp, : ) )
      end do

      ! copy the results to the res
      do iFld = 1, nFields
        res( (iElem-1)*nFields*3 + (iField-1)*3 + 1) = momentum(1, iFld)
        res( (iElem-1)*nFields*3 + (iField-1)*3 + 2) = momentum(2, iFld)
        res( (iElem-1)*nFields*3 + (iField-1)*3 + 3) = momentum(3, iFld)
      end do
    end do

  end subroutine deriveMomentaMSGas_FromState
! ***************************************************************************** !


! ****************************************************************************** !
  !> This routine computes equilibrium from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveEqMSGas_FromState( state, iField, nElems, varSys, layout, &
    &                                 res                                    )
    ! -------------------------------------------------------------------- !
    !> Array of state
    !! n * layout%fStencil%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iComp, nFields, iDir
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: varPos, QQ
    real(kind=rk) :: tmpPDF(layout%fStencil%QQ)
    !mass fraction of nSpecies
    real(kind=rk) :: massFraction( varSys%nStateVars )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !partial pressure of  nSpecies
    real(kind=rk) :: press( varSys%nStateVars )
    real(kind=rk) :: lambda( varSys%nStateVars )
    real(kind=rk) :: pressMix, densMix, molWeightMix
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum( 3, varSys%nStateVars )
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: chi( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    integer :: iFld, iFldDia, iFldNonDia
    real(kind=rk) :: matrixA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: invA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: fEq(layout%fStencil%QQ)
    real(kind=rk) :: vel( 3, varSys%nStateVars )
    real(kind=rk) :: eqVel(3)
    real(kind=rk) :: ucx, usq, weight0Inv
    integer :: stateVarMap(varSys%nStateVars)
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = layout%fStencil%QQ
    stateVarMap = scheme%stateVarMap%varPos%val(:)

    weight0Inv = 1.0_rk / layout%weight(layout%fStencil%restPosition)

    do iFld = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFld) = scheme%field(iFld)%fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                     &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    do iElem = 1, nElems
      do iFld = 1, nFields

        do iComp = 1, QQ
          varPos = varSys%method%val(stateVarMap(iFld))%state_varPos(iComp)

          tmpPDF(iComp) = state(                       &
            ! position of this state variable in the state array
            & ( ielem-1)* varsys%nscalars+varpos )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum( tmpPDF )

        ! partial pressure of species
        press(iFld) = mass_dens(iFld) * phi(iFld) * cs2

        ! velocity, first moments
        do iComp = 1, 3
          first_moments(iComp, iFld) = sum( tmpPDF * &
            & layout%fStencil%cxDirRK(iComp, :) )
        end do

      end do !iFld

      ! total pressure
      pressMix = sum(press)

      ! total density
      densMix = sum(mass_dens)

      ! mass freaction
      massFraction(:) = mass_dens(:) / densMix

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iFld = 1, nFields
        molWeightMix = molWeightMix + massFraction(iFld)/molWeight(iFld)
      end do
      molWeightMix = 1.0_rk/molWeightMix

      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      do iFld = 1, nFields
        do iFldDia = 1, nFields
          chi(iFld,iFldDia) = ( molWeightMix*molWeightMix       &
            & / (molWeight(iFld)*molWeight(iFldDia)) )          &
            & * (resi_coeff(iFld, iFldDia)/resi_coeff(iFld,iFld))
        enddo
      enddo

      !relaxation time
      do iFld = 1, nFields
        lambda(iFld) = pressMix*resi_coeff(iFld,iFld)/densMix
      enddo

      ! build up the equation system for momentum
      matrixA = 0.0_rk
      do iFld = 1, nFields
        ! set diagonal part
        matrixA( iFld, iFld ) = 1.0_rk
        do iFldDia = 1, nFields
          matrixA( iFld, iFld ) = matrixA( iFld, iFld )    &
            & + lambda(iFld) * 0.5_rk * chi(iFld, iFldDia) &
            & * massFraction( iFldDia )
        end do
        ! set nonDiagonal
        do iFldNonDia = 1,nFields
          matrixA( iFld, iFldNonDia ) = matrixA( iFld, iFldNonDia ) &
            & - lambda(iFld) * 0.5_rk * chi( iFld, iFldNonDia )     &
            & * massFraction( iFld )
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matrixA )

      ! momentum of all species
      momentum = 0.0_rk
      do iComp = 1, 3
        momentum( iComp, : ) = matmul( invA, first_moments( iComp, : ) )
      end do

      !velocity of all species
      do iFld = 1, nFields
        vel( :, iFld) = momentum( :, iFld) / mass_dens(iFld)
      end do

      eqVel(:) = vel(:,iField)
      do ifld = 1, nFields
        eqVel(:) = eqVel(:) + chi(ifld, iField)                       &
          &        * massfraction(ifld) * ( vel(:,ifld) - vel(:,iField) )
      end do

      ! Calculate the square of velocity
      usq = dot_product( eqVel,eqVel ) * t2cs2inv

      do iDir = 1, layout%fStencil%QQ
        ! Velocity times lattice unit velocity
        ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), eqVel )

        ! calculate equilibrium
        fEq(iDir) = layout%weight( iDir )            &
          &           * mass_dens(iField)            &
          &           * ( phi(iField) + ucx * cs2inv &
          &             + ucx * ucx * t2cs4inv - usq )
      end do ! iDir

      fEq(layout%stencil%restPosition) &
        & = layout%weight( layout%fStencil%restPosition ) * mass_dens(iField) &
        &   * ( weight0Inv + (1.0_rk - weight0Inv) * phi(iField) - usq )

      ! copy the results to the res
      res( (iElem-1)*QQ + 1 : iElem*QQ) = fEq

    end do !iElem

  end subroutine deriveEqMSGas_FromState
! ****************************************************************************** !


! **************************************************************************** !
  !> This routine computes auxField 'density and velocity' of given field
  !! from state array.
  !! velocity of original PDF is computed in this routine by solving LSE
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_auxFromState]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:auxFieldFromState]] function pointer.
  subroutine deriveAuxMSGas_fromState( derVarPos, state, neigh, iField,        &
    &                                  nElems, nSize, iLevel, stencil, varSys, &
    &                                  auxField, quantities                    )
    ! -------------------------------------------------------------------- !
    !> Position of derive variable in variable system
    class(mus_derVarPos_type), intent(in) :: derVarPos
    !> Array of state
    !! n * layout%stencil(1)%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> connectivity vector
    integer, intent(in) :: neigh(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> number of elements in state array
    integer, intent(in) :: nSize

    !> current level
    integer, intent(in) :: iLevel

    !> stencil header contains discrete velocity vectors
    type(tem_stencilHeader_type), intent(in) :: stencil

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities

    !> Output of this routine
    !! Size: nElems*nAuxScalars
    real(kind=rk), intent(inout) :: auxField(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iDir, nFields, dens_pos, mom_pos(3), elemOff
    real(kind=rk) :: massFraction( varSys%nStateVars )
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars ), inv_rho
    !partial pressure of  nSpecies
    real(kind=rk) :: press( varSys%nStateVars )
    real(kind=rk) :: lambda( varSys%nStateVars )
    real(kind=rk) :: pressMix, densMix, molWeightMix
    !first moments of nSpecies
    real(kind=rk) :: first_moments( 3, varSys%nStateVars )
    !momentum from linear system of equation
    real(kind=rk) :: momentum(3)
    !parameters from solver specific conf
    !field specific info from field table
    real(kind=rk), dimension(varSys%nStateVars) :: molWeight, phi
    real(kind=rk) :: resi_coeff( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: chi( varSys%nStateVars, varSys%nStateVars )
    !mixture info
    integer :: iFld, iFldDia, iFldNonDia
    real(kind=rk) :: matrixA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: invA( varSys%nStateVars, varSys%nStateVars )
    real(kind=rk) :: pdf( stencil%QQ )
    type(mus_scheme_type), pointer :: scheme
    type(mus_varSys_data_type), pointer :: fPtr
    ! ------------------------------------------------------------------------ !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields

    do iFld = 1, nFields
      ! species properties
      ! molecular weight inverse
      molWeight(iFld) = scheme%field(iFld)%fieldProp%species%molWeight
      ! molecular weight ratio
      phi(iFld) = scheme%field(iFld)%fieldProp%species%molWeigRatio
      ! resistivity coefficients
      resi_coeff(iFld, :) =                                &
        & scheme%field(iFld)%fieldProp%species%resi_coeff(:)
    end do

    !NEC$ ivdep
    do iElem = 1, nElems
      ! Compute density and momentum for all species
      !NEC$ shortloop
      do iFld = 1, nFields
        do iDir = 1, stencil%QQ
          pdf(iDir) = state(                                           &
& ( ielem-1)* varsys%nscalars+idir+( ifld-1)* stencil%qq )
        end do

        ! mass density of species
        mass_dens(iFld ) = sum(pdf)

        ! partial pressure of species
        press(iFld) = mass_dens(iFld) * phi(iFld) * cs2

        ! momentum
        first_moments(1, iFld) = sum( pdf * stencil%cxDirRK(1, :) )
        first_moments(2, iFld) = sum( pdf * stencil%cxDirRK(2, :) )
        first_moments(3, iFld) = sum( pdf * stencil%cxDirRK(3, :) )
      end do

      ! total pressure
      pressMix = sum(press)

      ! total density
      densMix = sum(mass_dens)

      ! mass freaction
      massFraction(:) = mass_dens(:) / densMix

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iFld = 1, nFields
        molWeightMix = molWeightMix + massFraction(iFld)/molWeight(iFld)
      end do
      molWeightMix = 1.0_rk / molWeightMix

      !chi = (m^2/(m_\sigma m_\varsigma))*(B_(\sigma \varsigma)
      do iFld = 1, nFields
        do iFldDia = 1, nFields
          chi(iFld,iFldDia) = ( molWeightMix*molWeightMix       &
            & / (molWeight(iFld)*molWeight(iFldDia)) )          &
            & * (resi_coeff(iFld, iFldDia)/resi_coeff(iFld,iFld))
        enddo
      enddo

      !relaxation time
      do iFld = 1, nFields
        lambda(iFld) = pressMix*resi_coeff(iFld,iFld)/densMix
      enddo

      ! build up the equation system for momentum
      matrixA = 0.0_rk
      !NEC$ shortloop
      do iFld = 1, nFields
        ! set diagonal part
        matrixA( iFld, iFld ) = 1.0_rk
        !NEC$ shortloop
        do iFldDia = 1, nFields
          matrixA( iFld, iFld ) = matrixA( iFld, iFld )    &
            & + lambda(iFld) * 0.5_rk * chi(iFld, iFldDia) &
            & * massFraction( iFldDia )
        end do
        ! set nonDiagonal
        !NEC$ shortloop
        do iFldNonDia = 1,nFields
          matrixA( iFld, iFldNonDia ) = matrixA( iFld, iFldNonDia ) &
            & - lambda(iFld) * 0.5_rk * chi( iFld, iFldNonDia )     &
            & * massFraction( iFld )
        end do
      end do

      ! invert matrix
      invA = invert_matrix( matrixA )

      ! momentum of current species
      momentum(1) = dot_product( invA(iField, :), first_moments(1, :) )
      momentum(2) = dot_product( invA(iField, :), first_moments(2, :) )
      momentum(3) = dot_product( invA(iField, :), first_moments(3, :) )

      ! position of density and velocity of current field in auxField array
      dens_pos = varSys%method%val(derVarPos%density)%auxField_varPos(1)
      mom_pos = varSys%method%val(derVarPos%momentum)%auxField_varPos(:)

      ! element offset for auxField
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! store field density
      auxField(elemOff + dens_pos) = mass_dens(iField)
      ! store field velocity
      inv_rho = 1.0_rk/mass_dens(iField)
      auxField(elemOff + mom_pos(1)) = momentum(1) * inv_rho
      auxField(elemOff + mom_pos(2)) = momentum(2) * inv_rho
      auxField(elemOff + mom_pos(3)) = momentum(3) * inv_rho
    end do

  end subroutine deriveAuxMSGas_fromState
! **************************************************************************** !

  ! ************************************************************************** !
  !> This routine computes equilbrium from auxField
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_equilFromAux]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:equilFromAux]] function pointer.
  subroutine deriveEquilMSGas_fromAux( derVarPos, auxField, iField, nElems, &
    &                                  varSys, layout, fEq                  )
    ! -------------------------------------------------------------------- !
    !> Position of derive variable in variable system
    class(mus_derVarPos_type), intent(in) :: derVarPos
    !> Array of auxField.
    !! Single species: dens_1, vel_1, dens_2, vel_2, .. dens_n, vel_n
    !! multi-species: dens_1_sp1, vel_1_spc1, dens_1_sp2, vel_1_spc2,
    !!                dens_2_sp1, vel_2_spc2, dens_2_sp2, vel_2_spc2 ...
    !!                dens_n_sp1, vel_n_sp1, dens_n_sp2, vel_n_spc2
    !! Access: (iElem-1)*nAuxScalars + auxField_varPos
    real(kind=rk), intent(in) :: auxField(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n*QQ of res
    real(kind=rk), intent(out) :: fEq(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iDir, iFld, nFields, elemOff, elemOffEq
    integer :: dens_pos, mom_pos(3)
    real(kind=rk) :: resi_coeff(varSys%nStateVars), phi
    !mass density of nSpecies
    real(kind=rk) :: mass_dens( varSys%nStateVars )
    !mass fraction
    real(kind=rk) :: massFraction( varSys%nStateVars )
    real(kind=rk) :: totMass_densInv
    real(kind=rk) :: molWeightInv( varSys%nStateVars )
    real(kind=rk) :: vel( 3, varSys%nStateVars ), eqVel(3)
    real(kind=rk) :: ucx, usq, weight0Inv, molWeightMix
    type(mus_scheme_type), pointer :: scheme
    type(mus_varSys_data_type), pointer :: fPtr
    ! ------------------------------------------------------------------------ !
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    ! resistivity coefficients
    resi_coeff(:) =                                     &
      & scheme%field(iField)%fieldProp%species%resi_coeff(:)

    weight0Inv = 1.0_rk / layout%weight(layout%fStencil%restPosition)

    phi = scheme%field(iField)%fieldProp%species%molWeigRatio

    do ifld = 1, nFields
      ! molecular weight inverse
      molWeightInv(ifld) = scheme%field(ifld)%fieldProp%species%molWeightInv
    end do

    !NEC$ ivdep
    do iElem = 1, nElems
      ! element offset
      elemoff = (iElem-1)*varSys%nAuxScalars

      ! get all field density and velocity to compute equilibrium velocity
      !NEC$ shortloop
      do iFld = 1, nFields
        ! field density
        dens_pos = varSys%method%val(scheme%derVarPos(iFld)%density) &
          &                     %auxField_varPos(1)
        mass_dens(iFld) = auxField(elemOff + dens_pos)

        ! field velocity
        mom_pos = varSys%method%val(scheme%derVarPos(iFld)%momentum) &
          &                    %auxField_varPos(:)
        vel(1, iFld) = auxField(elemOff+mom_pos(1))
        vel(2, iFld) = auxField(elemOff+mom_pos(2))
        vel(3, iFld) = auxField(elemOff+mom_pos(3))
      end do !iFld

      totmass_densInv = 1.0_rk/sum(mass_dens)
      ! massfraction
      massFraction = mass_dens * totmass_densInv

      ! mixture molecular weight
      !1/mm = \sum_\sigma massfraction_\sigma/m_\sigma
      molWeightMix = 0.0_rk
      do iFld = 1, nFields
        molWeightMix = molWeightMix &
          &          + massFraction(iFld) * molWeightInv(iFld)
      end do
      molWeightMix = 1.0_rk/molWeightMix

      eqVel(:) = vel(:, iField)
      do iFld = 1, nFields
        eqVel(:) = eqVel(:) + molWeightMix * molWeightMix  &
        &        * molWeightInv(iField)*molWeightInv(iFld) &
        &        * resi_coeff(iFld) * massFraction(iFld)   &
        &        * ( vel(:,iFld) - vel(:,iField) )         &
        &        / resi_coeff(iField)
      end do

      ! Calculate the square of velocity
      usq = ( eqVel(1)*eqVel(1) + eqVel(2)*eqVel(2)            &
        &                       + eqVel(3)*eqVel(3) ) * t2cs2inv

      ! offset of equilibrium
      elemOffEq = (iElem-1)*layout%fStencil%QQ

      !NEC$ shortloop
      do iDir = 1, layout%fStencil%QQ
        ! Velocity times lattice unit velocity
        ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), eqVel )

        ! calculate equilibrium
        fEq(elemOffEq+iDir) = layout%weight(iDir) * mass_dens(iField)     &
          &                 * ( phi + ucx * cs2inv + ucx * ucx * t2cs4inv &
          &                 - usq                                         )
      end do

      ! update rest position equilibrium
      fEq(elemOffEq + layout%stencil%restPosition) =       &
        & layout%weight( layout%fStencil%restPosition )    &
        & * mass_dens(iField)                              &
        & * ( weight0Inv + (1.0_rk - weight0Inv)*phi - usq )

    end do

  end subroutine deriveEquilMSGas_fromAux
  ! ************************************************************************** !

end module mus_derQuanMSGas_module
! ****************************************************************************** !
