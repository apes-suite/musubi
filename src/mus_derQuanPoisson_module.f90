! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2017-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!> author: sindhuja
!! This module provides the MUSUBI specific functions for calculating
!! macroscopic quantities from the state variables.
!! The depending common interface between MUSUBI and ATELES is defined in the
!! tem_derived_module. The functionality for accessing a variable from the state
!! and evaluating a lua function are also provided in the tem_derived module.
!! A Novel lattice boltzmann model for poisson equation
!! author> Zhenhua Chai , Baochang Shi
!! A Coupled Lattice Boltzmann method to solve  Nernst -Planck Model for
!! simulating Electro-Osmotic Flows
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
module mus_derQuanPoisson_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use tem_param_module,         only: cs2, cs2inv
  use env_module,               only: rk, labelLen
  use tem_variable_module,      only: tem_variable_type
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_topology_module,      only: tem_levelOf
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_logging_module,       only: logUnit
  use tem_varSys_module,        only: tem_varSys_type, tem_varSys_op_type,     &
    &                                 tem_varSys_append_derVar,                &
    &                                 tem_varSys_append_auxFieldVar,           &
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
  use tem_aux_module,           only: tem_abort
  use tem_debug_module,         only: dbgUnit
  use tem_grow_array_module,    only: grw_labelarray_type, append

  ! include musubi modules
  use mus_source_type_module,        only: mus_source_op_type
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_varSys_module,             only: mus_varSys_data_type,             &
    &                                      mus_varSys_solverData_type,       &
    &                                      mus_get_new_solver_ptr,           &
    &                                      mus_deriveVar_ForPoint,           &
    &                                      mus_generic_varFromPDF_fromIndex, &
    &                                      mus_generic_fromPDF_forElement,   &
    &                                      mus_derive_fromPDF
  use mus_stateVar_module,           only: mus_stateVar_Fetch_now_fromIndex,   &
    &                                      mus_accessVar_setupIndices,         &
    &                                      mus_access_stateFetch_now_forElement
  use mus_auxFieldVar_module,        only: mus_access_auxFieldVar_forElement, &
    &                                      mus_auxFieldVar_forPoint,          &
    &                                      mus_auxFieldVar_fromIndex
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_operation_var_module,      only: mus_opVar_setupIndices
  use mus_derivedQuantities_module2, only: convPrePost
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_physics_module,            only: mus_convertFac_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none

  private

  public :: mus_append_derVar_poisson
  public :: deriveAuxPoisson_fromState
  public :: deriveEquilPoisson_fromAux

  public :: applySrc_chargeDensity_2ndOrd
  public :: applySrc_chargeDensity_1stOrd
  public :: deriveSrc_chargeDensity

contains

  ! **************************************************************************** !
  !> subroutine to add derive variables for weakly compressible PB
  !! (schemekind = 'poisson') to the varsys.
  !! A Coupled Lattice Boltzmann Method to Solve Nernst-Planck Model
  !! for Simulating Electro-Osmotic flows
  !! author> Xuguang yang
  subroutine mus_append_derVar_poisson( varSys, solverData, fldLabel, &
    &                                   derVarName, schemeKind, stencil )
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)  :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> compute stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel

    !> array of derive physical variables
    type(grw_labelarray_type), intent(inout) :: derVarName

    !> scheme kind
    character(len=*), intent(in) :: schemeKind
    ! ---------------------------------------------------------------------------
    ! number of derive variables
    integer :: nDerVars, iVar, nComponents, addedPos, iIn
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

    if (trim(schemeKind) == 'poisson') then
      ! For poisson equation, chargedensity variable is defined as source
      ! term
      nDerVars = 2
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'fetch_pdf_now    ', &
        &                   'electric_field   ']
    else
      ! For poisson_boltzmann_linear and poisson_boltzmann_nonlinear
      nDerVars = 3
      allocate(derVarName_loc(nDerVars))
      derVarName_loc    = [ 'fetch_pdf_now    ', &
        &                   'charge_density   ', 'electric_field   ']
    end if

    do iVar = 1, nDerVars
      call append(derVarName, derVarName_loc(iVar))

      ! set default pointers, overwrite if neccessary
      get_element => tem_varSys_getElement_dummy
      get_point => mus_deriveVar_ForPoint
      setup_indices => mus_opVar_setupIndices
      get_valOfIndex => tem_varSys_getValOfIndex_dummy
      method_data  = mus_get_new_solver_ptr(solverData)
      set_params => tem_varSys_setParams_dummy
      get_params => tem_varSys_getParams_dummy

      select case(trim(adjustl(derVarName_loc(iVar))))

      case ('fetch_pdf_now')
        get_element => mus_access_stateFetch_now_forElement
        get_valOfIndex => mus_stateVar_Fetch_now_fromIndex
        setup_indices => mus_accessVar_setupIndices
        nComponents = stencil%QQ
        allocate(input_varname(1))
        input_varname(1) = 'pdf'

      case ('electric_field')
        get_element => deriveElectricfield_forElement
        get_valOfIndex => deriveElectricfield_fromIndex
        nComponents = 3
        allocate(input_varname(1))
        input_varname(1) = 'pdf'

      case ('charge_density_boltzmann')
        ! Charge density from potential using boltzmann approximation
        get_element => deriveChargeDensityBoltzAppr_forElement
        get_valOfIndex => deriveChargeDensityBoltzAppr_fromIndex
        nComponents = 1
        allocate(input_varname(1))
        input_varname(1) = 'potential'

      case default
        write(logUnit(1),*) 'WARNING: Unknown variable: '//&
          &                 trim(derVarName_loc(iVar))
        cycle !go to next variable
      end select

      ! update variable names with field label
      varname = trim(fldLabel)//trim(adjustl(derVarName_loc(iVar)))
      do iIn = 1, size(input_varname)
        input_varname(iIn) = trim(fldLabel)//trim(input_varname(iIn))
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

  end subroutine mus_append_derVar_poisson
  ! ************************************************************************** !

! ****************************************************************************** !
!       Subroutines with common interface for the function pointers            !
! ****************************************************************************** !

! **************************************************************************** !
  !> This routine computes auxField 'potential' from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_auxFromState]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:auxFieldFromState]] function pointer.
  subroutine deriveAuxPoisson_fromState( derVarPos, state, neigh, iField,  &
    &                                    nElems, nSize, iLevel, stencil,   &
    &                                    varSys, auxField, quantities )
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
    integer :: iElem, iDir
    real(kind=rk) :: pdf( stencil%QQ )
    ! ------------------------------------------------------------------------ !
    !NEC$ ivdep
    do iElem = 1, nElems
      !NEC$ shortloop
      do iDir = 1, stencil%QQ
        pdf(iDir) = state(                                                     &
          & ( ielem-1)* varsys%nscalars+idir+( 1-1)* stencil%qq)
      end do

      ! element offset is not required because poisson equation has only
      ! one aux scalar 'potential'
      auxField(iElem) = sum(pdf)
    end do

  end subroutine deriveAuxPoisson_fromState
! **************************************************************************** !

  ! ************************************************************************** !
  !> This routine computes equilbrium from auxField
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_equilFromAux]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:equilFromAux]] function pointer.
  subroutine deriveEquilPoisson_fromAux( derVarPos, auxField, iField, nElems, &
    &                                    varSys, layout, fEq                  )
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
    integer :: iElem, iDir
    ! ------------------------------------------------------------------------ !
    !NEC$ ivdep
    do iElem = 1, nElems
      !NEC$ shortloop
      do iDir = 1, layout%fStencil%QQ
        ! element offset is not required for auxField because poisson equation
        ! has only one aux scalar 'potential'
        fEq((iElem-1)*layout%fStencil%QQ+iDir) = layout%weight(iDir) &
          &                                    * auxField(iElem)
      end do
    end do

  end subroutine deriveEquilPoisson_fromAux
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate the electric_field of a given pre-collision pdfs
  !! i.e fetch_pdf_now
  !!
  !! The interface has to comply to the abstract interface
  !! [[mus_varSys_module:mus_derive_fromPDF]].
  !!
  recursive subroutine mus_deriveElectricField( fun, varsys, stencil, iLevel, &
    &                                             posInState, pdf, res, nVals )
    !> description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varsys_op_type), intent(in)     :: fun
    !> the variable system to obtain the variable from.
    type(tem_varsys_type), intent(in)         :: varsys
    !> fluid stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil
    !> current Level
    integer, intent(in)                       :: iLevel
    !> Position of element in levelwise state array
    integer, intent(in)                       :: posInState(:)
    !> pdf array contains pre-collision using FETCH and nNow
    real(kind=rk), intent(in)                 :: pdf(:)
    !> results
    real(kind=rk), intent(out)                :: res(:)
    !> nVals to get
    integer, intent(in)                       :: nVals
    ! ---------------------------------------------------------------------------
    integer :: iVal, iComp
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: nCompPDF
    real(kind=rk) :: pot, electric_prefactor, electric_field(3)
    real(kind=rk), allocatable :: pdfTmp(:)
    real(kind=rk) :: omega
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    nCompPDF = varSys%method%val(fun%input_varPos(1))%nComponents
    allocate(pdfTmp(nCompPDF))
    omega = fPtr%solverData%scheme%field(1)%fieldProp%poisson%omega
    electric_prefactor = ( omega * convPrePost(omega) ) / scheme%layout%cs**2

    ! res is always AOS layout
    res = 0.0_rk

    do iVal = 1, nVals
      ! Calculate potential
      pot = sum( pdf( (iVal-1)*nCompPDF + 1: iVal*nCompPDF) )
      ! Calculate nonequilibrium
      do iComp = 1, nCompPDF
        pdfTmp(iComp) = pdf( (iVal-1)*nCompPDF + iComp)
      end do

      ! Calculate electric_field by Multiply the difference with cxDir
      ! and additional pre-factor
      ! E = omega/cs2 * sum (cx * fPre) or
      ! E = omega/cs2/(1-omega) * sum (cx * fPost)
      electric_field(1) = sum( pdfTmp * scheme%layout%fStencil%cxDirRK(1,:) )
      electric_field(2) = sum( pdfTmp * scheme%layout%fStencil%cxDirRK(2,:) )
      electric_field(3) = sum( pdfTmp * scheme%layout%fStencil%cxDirRK(3,:) )

      electric_field = electric_prefactor * electric_field

      ! Store electric field in result output
      res((iVal-1)*3+1) = electric_field(1)
      res((iVal-1)*3+2) = electric_field(2)
      res((iVal-1)*3+3) = electric_field(3)

    end do !iVal

  end subroutine mus_deriveElectricField
! ****************************************************************************** !


  ! ************************************************************************** !
  !> Calculate the electric field of a given set of elements (sum up all links).
  !! This routine is used to compute electric field for all scheme kinds
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveElectricfield_forElement(fun, varsys, elempos, &
    &                                                 time, tree, nElems,   &
    &                                                 nDofs, res            )
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer :: fnCalcPtr
    ! -------------------------------------------------------------------- !

    fnCalcPtr => mus_deriveElectricField

    call mus_generic_fromPDF_forElement( &
      &  fun       = fun,                &
      &  varSys    = varSys,             &
      &  elempos   = elempos,            &
      &  tree      = tree,               &
      &  time      = time,               &
      &  nVals     = nElems,             &
      &  fnCalcPtr = fnCalcPtr,          &
      &  nDofs     = nDofs,              &
      &  res       = res                 )

  end subroutine deriveElectricfield_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculate the electric field of a given set of elements (sum up all links).
  !! This routine is used to compute electric field for all scheme kinds
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveElectricfield_fromIndex(fun, varSys, time,   &
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer  :: fnCalcPtr
    ! -------------------------------------------------------------------- !

    fnCalcPtr => mus_deriveElectricField

    call mus_generic_varFromPDF_fromIndex( &
      &  fun       = fun,                  &
      &  varSys    = varSys,               &
      &  time      = time,                 &
      &  iLevel    = iLevel,               &
      &  idx       = idx,                  &
      &  nVals     = nVals,                &
      &  fnCalcPtr = fnCalcPtr,            &
      &  res       = res                   )

  end subroutine deriveElectricfield_fromIndex
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculate charge density from potential field using Boltzmann approximation
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveChargeDensityBoltzAppr_forElement( &
    & fun, varsys, elempos, time, tree, nElems, nDofs, res      )
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
    integer :: iElem, iIon, potPos
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_field_prop_type), pointer :: fieldProp
    real(kind=rk) :: potential(nElems)
    real(kind=rk) :: fac, pot_fac, charge_dens
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    fieldProp => fPtr%solverData%scheme%field(1)%fieldProp

    potPos = fun%input_varPos(1)
    call varSys%method%val(potPos)%get_element( &
        & varSys  = varSys,                     &
        & elemPos = elemPos,                    &
        & time    = time,                       &
        & tree    = tree,                       &
        & nElems  = nElems,                     &
        & nDofs   = nDofs,                      &
        & res     = potential                   )

    fac = fieldProp%poisson%PB%faradayLB                                  &
      & / (fieldProp%poisson%PB%gasConst_R_LB * fieldProp%poisson%PB%temp )
    ! res is always AOS layout
    res = 0.0_rk
    do iElem = 1, nElems
      pot_fac = fac * potential(iElem)
      charge_dens = 0.0_rk
      do iIon = 1, fieldProp%poisson%PB%nIons
        charge_dens = charge_dens + fieldProp%poisson%PB%moleDens0 &
          & * fieldProp%poisson%PB%faradayLB                       &
          & * fieldProp%poisson%PB%valence(iIon)                   &
          & * exp( - fieldProp%poisson%PB%valence(iIon) * pot_fac  )
      end do
      res(iElem) = charge_dens
    end do

  end subroutine deriveChargeDensityBoltzAppr_forElement
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Calculate charge density from potential field using Boltzmann approximation
  !! from given index
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine deriveChargeDensityBoltzAppr_fromIndex( &
    & fun, varSys, time, iLevel, idx, idxLen, nVals, res       )
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
    integer :: iVal, iIon, potPos
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_field_prop_type), pointer :: fieldProp
    real(kind=rk) :: potential(nVals)
    real(kind=rk) :: fac, pot_fac, charge_dens
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
    fieldProp => fPtr%solverData%scheme%field(1)%fieldProp


    potPos = fun%input_varPos(1)

    call varSys%method%val(potPos)%get_valOfIndex(   &
        & varSys  = varSys,                          &
        & time    = time,                            &
        & iLevel  = iLevel,                          &
        & idx     = fPtr%opData%input_pntIndex(1)    &
        &           %indexLvl(iLevel)%val( idx(:) ), &
        & nVals   = nVals,                           &
        & res     = potential                        )

    fac = fieldProp%poisson%PB%faradayLB                                  &
      & / (fieldProp%poisson%PB%gasConst_R_LB * fieldProp%poisson%PB%temp )
    ! res is always AOS layout
    res = 0.0_rk
    do iVal = 1, nVals
      pot_fac = fac * potential(iVal)
      charge_dens = 0.0_rk
      do iIon = 1, fieldProp%poisson%PB%nIons
        charge_dens = charge_dens + fieldProp%poisson%PB%moleDens0 &
          & * fieldProp%poisson%PB%faradayLB                       &
          & * fieldProp%poisson%PB%valence(iIon)                   &
          & * exp( - fieldProp%poisson%PB%valence(iIon) * pot_fac  )
      end do
      res(iVal) = charge_dens
    end do

  end subroutine deriveChargeDensityBoltzAppr_fromIndex
  ! ************************************************************************** !


! ****************************************************************************** !
  !> Update state with source variable "ChargeDensity" with 2nd order
  !! integration of source Term.
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !!
  !! $$ \nabla^2 \phi = - \frac{\rho_e}{\epsilon_r \epsilon_0} $$
  !!
  !! Where \rho_e is the charge density
  !! KM: LBE solves potential equation of form
  !! $$ \partial_t \phi + \gamma \nabla^2 \phi + \gamma S = 0 $$
  !! where
  !! S = \frac{\rho_e}{\epsilon_r \epsilon_0}
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_chargeDensity_2ndOrd( fun, inState, outState, neigh, &
    &                                       auxField, nPdfSize, iLevel,    &
    &                                       varSys, time, phyConvFac,      &
    &                                       derVarPos                      )
    ! -------------------------------------------------------------------- !
    !> Description of method to apply source terms
    class(mus_source_op_type), intent(in) :: fun

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: rhs(fun%elemLvl(iLevel)%nElems)
    real(kind=rk) :: rhs_Fac
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal, statePos
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer(varSys%method%val(fun%srcTerm_varPos)%method_data, fPtr)

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems
    associate( scheme => fPtr%solverData%scheme,                             &
      &        poisson => fPtr%solverData%scheme%field(1)%fieldProp%poisson, &
      &        physics => fPtr%solverData%physics                            )
      ! factor to multiply rhs
      rhs_Fac = (1.0_rk - poisson%omega * 0.5_rk ) * poisson%pot_diff &
        &     / poisson%permittivity

      ! Get charge density which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = rhs                                       )

      ! convert physical to lattice and multiply with rhs factor
      rhs = (rhs / physics%fac(iLevel)%chargeDens) * rhs_Fac

      ! constant parameter
      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      !$omp do schedule(static)
      do iElem = 1, nElems
  !write(*,*) 'rhs ', rhs(iElem)

        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

        ! Source term: = (1-omega/2) \gamma * \rho_e / \epsilon
        do iDir = 1, QQ
          statePos = ( posintotal-1)* nscalars+idir+( 1-1)* qq

          outState(statePos) = outState(statePos) + scheme%layout%weight(iDir) &
            &                * rhs(iElem)
        end do

      end do !iElem
      !$omp end do nowait
    end associate

  end subroutine applySrc_chargeDensity_2ndOrd
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "ChargeDensity" with 1st order
  !! integration of source Term.
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  !!
  !! $$ \nabla^2 \phi = - \frac{\rho_e}{\epsilon_r \epsilon_0} $$
  !!
  !! Where \rho_e is the charge density
  !! KM: LBE solves potential equation of form
  !! $$ \partial_t \phi + \gamma \nabla^2 \phi + \gamma S = 0 $$
  !! where
  !! S = \frac{\rho_e}{\epsilon_r \epsilon_0}
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_chargeDensity_1stOrd( fun, inState, outState, neigh, &
    &                                       auxField, nPdfSize, iLevel,    &
    &                                       varSys, time, phyConvFac,      &
    &                                       derVarPos                      )
    ! -------------------------------------------------------------------- !
    !> Description of method to apply source terms
    class(mus_source_op_type), intent(in) :: fun

    !> input  pdf vector
    real(kind=rk), intent(in) :: inState(:)

    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)

    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> number of elements in state Array
    integer, intent(in) :: nPdfSize

    !> current level
    integer, intent(in) :: iLevel

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: rhs(fun%elemLvl(iLevel)%nElems)
    real(kind=rk) :: rhs_Fac
    integer :: nElems, iElem, iDir, QQ, nScalars, posInTotal
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems
    rhs_Fac = scheme%field(1)%fieldProp%poisson%pot_diff   &
      &     / scheme%field(1)%fieldProp%poisson%permittivity


    ! Get charge density which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = rhs                                       )

    ! convert physical to lattice and multiply with rhs factor
    rhs = (rhs / fPtr%solverData%physics%fac(iLevel)%chargeDens) * rhs_Fac

    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars

!$omp do schedule(static)
    do iElem = 1, nElems
!write(*,*) 'rhs ', rhs(iElem)

      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! Source term:
      ! source = \gamma * \rho_e / \epsilon
      do iDir = 1, QQ

        outState(( posintotal-1)* nscalars+idir+( 1-1)* qq ) &
& = outState(( posintotal-1)* nscalars+idir+( 1-1)* qq ) &
          & + scheme%layout%weight(iDir) * rhs(iElem)

      end do

    end do !iElem
!$omp end do nowait

  end subroutine applySrc_chargeDensity_1stOrd
! ****************************************************************************** !

  ! ************************************************************************** !
  !> Calculate charge density source variable referred in config file
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveSrc_chargeDensity(fun, varsys, elempos, time, &
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
    integer :: data_varPos
    ! -------------------------------------------------------------------- !
    ! spacetime source variable
    data_varPos = fun%input_varPos(2)
    call varSys%method%val(data_varPos)%get_element( &
        & varSys  = varSys,                          &
        & elemPos = elemPos,                         &
        & time    = time,                            &
        & tree    = tree,                            &
        & nElems  = nElems,                          &
        & nDofs   = nDofs,                           &
        & res     = res                              )

  end subroutine deriveSrc_chargeDensity
! ****************************************************************************** !


end module mus_derQuanPoisson_module
