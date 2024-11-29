! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!> author: seyfettin bilgi
!! This module provides the MUSUBI specific functions for calculating
!! macroscopic quantities from the state variables.
!! The depending common interface between MUSUBI and ATELES is defined in the
!! tem_derived_module. The functionality for accessing a variable from the state
!! and evaluating a lua function are also provided in the tem_derived module.
!! A Novel lattice boltzmann model for nernstPlanck equation
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
module mus_derQuanNernstPlanck_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
  ! include treelm modules
  use env_module,               only: rk, long_k, globalMaxLevels, labelLen
  use tem_param_module,         only: div1_2, div1_3, div1_54, div1_9, div3_4, &
    &                                 sqrt3, cs2inv, cs2, t2cs2inv, t2cs4inv,  &
    &                                 cs4inv
  use tem_aux_module,           only: tem_abort
  use tem_varSys_module,        only: tem_varSys_type, tem_varSys_op_type,     &
    &                                 tem_varSys_append_derVar,                &
    &                                 tem_varSys_proc_point,                   &
    &                                 tem_varSys_proc_element,                 &
    &                                 tem_varSys_proc_setParams,               &
    &                                 tem_varSys_proc_getParams,               &
    &                                 tem_varSys_proc_setupIndices,            &
    &                                 tem_varSys_proc_getValOfIndex
  use tem_variable_module,      only: tem_variable_type
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_logging_module,       only: logUnit
  use tem_topology_module,      only: tem_levelOf
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_debug_module,         only: dbgUnit
  use tem_operation_var_module, only: tem_opVar_setupIndices,      &
    &                                 tem_get_new_varSys_data_ptr, &
    &                                 tem_evalAdd_forElement,      &
    &                                 tem_evalAdd_fromIndex,       &
    &                                 tem_opVar_setParams,         &
    &                                 tem_opVar_getParams
  use tem_grow_array_module,    only: grw_labelarray_type, append

  ! include musubi modules
  use mus_source_type_module,        only: mus_source_op_type
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_varSys_module,             only: mus_varSys_data_type,             &
    &                                      mus_varSys_solverData_type,       &
    &                                      mus_get_new_solver_ptr,           &
    &                                      mus_deriveVar_ForPoint,           &
    &                                      mus_generic_varFromPDF_fromIndex, &
    &                                      mus_generic_fromPDF_forElement,   &
    &                                      mus_derive_fromPDF
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_operation_var_module,      only: mus_opVar_setupIndices
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_physics_module,            only: mus_convertFac_type, &
    &                                      gasConst_R, faraday
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  ! include aotus modules
  use aotus_module, only: flu_State


  implicit none

  private

!KM!  public :: mus_append_derVar_nernstPlanck
  public :: deriveAuxNP_fromState
  public :: deriveEquilNP_fromAux

  public :: mus_deriveMoleDensity
  public :: deriveMoleDensity_forElement
  public :: deriveMoleDensity_fromIndex

!!  public :: deriveSrc_electricField
  public :: applySrc_electricFieldNP

contains

  ! **************************************************************************** !
  !> subroutine to add derive variables for weakly compressible PB
  !! (schemekind = 'nernstPlanck') to the varsys.
  !! A Coupled Lattice Boltzmann Method to Solve Nernst-Planck Model
  !! for Simulating Electro-Osmotic flows
  !! author> Xuguang yang
  !KM! \todo currently this is replaced by mus_append_derMixVar_MS since
  !KM! only mixture variable is appended in this routine
!KM!  subroutine mus_append_derVar_nernstPlanck( varSys, solverData, fldLabel, &
!KM!    &                                        nFields, derVarName )
!KM!    ! ---------------------------------------------------------------------------
!KM!    !> global variable system
!KM!    type(tem_varSys_type), intent(inout)  :: varSys
!KM!
!KM!    !> Contains pointer to solver data types
!KM!    type(mus_varSys_solverData_type), target, intent(in) :: solverData
!KM!
!KM!    !> number of fields
!KM!    integer, intent(in)                      :: nFields
!KM!
!KM!    !> array of field label prefix. Size=nFields
!KM!    character(len=*), intent(in)              :: fldLabel(:)
!KM!
!KM!    !> array of derive physical variables
!KM!    type(grw_labelarray_type), intent(inout) :: derVarName
!KM!    ! ---------------------------------------------------------------------------
!KM!    ! number of derive variables
!KM!    integer :: iVar, nComponents, addedPos, iField
!KM!    logical :: wasAdded
!KM!    character(len=labelLen), allocatable ::  input_varname(:)
!KM!    character(len=labelLen)  ::  varName
!KM!    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
!KM!    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
!KM!    procedure(tem_varSys_proc_setParams), pointer :: set_params => null()
!KM!    procedure(tem_varSys_proc_getParams), pointer :: get_params => null()
!KM!    procedure(tem_varSys_proc_setupIndices), pointer :: &
!KM!      &                                      setup_indices => null()
!KM!    procedure(tem_varSys_proc_getValOfIndex), pointer :: &
!KM!      &                                       get_valOfIndex => null()
!KM!    type(c_ptr) :: method_data
!KM!    character(len=labelLen) :: derVarName_loc
!KM!    ! ---------------------------------------------------------------------------
!KM!    nullify(get_point, get_element, set_params, get_params, setup_indices, &
!KM!      &     get_valOfIndex)
!KM!
!KM!    derVarName_loc = 'mole_density'
!KM!    ! mole_density of each field is already added as auxiliary variable
!KM!    ! so only add mole_density for mixture
!KM!
!KM!    ! set pointers for mixture. mixture quantity is obtained by summing up
!KM!    ! species quantities
!KM!    get_element => tem_evalAdd_forElement
!KM!    get_point => mus_deriveVar_ForPoint
!KM!    setup_indices => tem_opVar_setupIndices
!KM!    get_valOfIndex => tem_evalAdd_fromIndex
!KM!    method_data = tem_get_new_varSys_data_ptr(method_data)
!KM!    set_params => tem_opVar_setParams
!KM!    get_params => tem_opVar_getParams
!KM!
!KM!    nComponents = 1
!KM!    varname = trim(adjustl(derVarName_loc))
!KM!    allocate(input_varname(nFields))
!KM!    do iField = 1, nFields
!KM!      input_varname(iField) = trim(fldLabel(iField))//trim(varname)
!KM!    end do
!KM!
!KM!    ! append variable to varSys
!KM!    call tem_varSys_append_derVar(  me             = varSys,         &
!KM!      &                             varName        = trim(varname),  &
!KM!      &                             nComponents    = nComponents,    &
!KM!      &                             input_varname  = input_varname,  &
!KM!      &                             method_data    = method_data,    &
!KM!      &                             get_point      = get_point,      &
!KM!      &                             get_element    = get_element,    &
!KM!      &                             set_params     = set_params,     &
!KM!      &                             get_params     = get_params,     &
!KM!      &                             setup_indices  = setup_indices,  &
!KM!      &                             get_valOfIndex = get_valOfIndex, &
!KM!      &                             pos            = addedPos,       &
!KM!      &                             wasAdded       = wasAdded        )
!KM!
!KM!    if (wasAdded) then
!KM!      write(logUnit(10),*) ' Appended variable:'//trim(varname)
!KM!    else if (addedpos < 1) then
!KM!      write(logUnit(1),*) 'Error: variable '//trim(varname)// &
!KM!        &                 ' is not added to variable system'
!KM!    end if
!KM!
!KM!    deallocate(input_varname)
!KM!
!KM!  end subroutine mus_append_derVar_nernstPlanck
  ! ************************************************************************** !


! **************************************************************************** !
  !> This routine computes auxField 'mole_density' from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_auxFromState]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:auxFieldFromState]] function pointer.
  subroutine deriveAuxNP_fromState( derVarPos, state, neigh, iField, nElems, &
    &                               nSize, iLevel, stencil, varSys, auxField, quantities )
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

      ! element offset is not required because passive scalar has only
      ! one aux scalar
      ! density
      auxField(iElem) = sum(pdf)
    end do

  end subroutine deriveAuxNP_fromState
! **************************************************************************** !

  ! ************************************************************************** !
  !> This routine computes equilbrium from auxField
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_equilFromAux]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:equilFromAux]] function pointer.
  subroutine deriveEquilNP_fromAux( derVarPos, auxField, iField, nElems, &
    &                               varSys, layout, fEq                  )
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
    !KM: \todo add transport velocity as auxField
    write(logUnit(1),*) 'ERROR: Equilibrium calculation requires transport '
    write(logUnit(1),*) 'velocity and it is not provided to this routine'
    write(logUnit(1),*) 'Solution: use deriveEquilPS_fromMacro'
    call tem_abort()
  end subroutine deriveEquilNP_fromAux
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate the potential of a given set of pdfs of elements
  !!
  !! The interface has to comply to the abstract interface
  !! [[mus_varSys_module:mus_derive_fromPDF]].
  !!
  recursive subroutine mus_deriveMoleDensity(fun, varsys, stencil, iLevel, &
    &                                          posInState, pdf, res, nVals )
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
    !> pdf array
    real(kind=rk), intent(in)                 :: pdf(:)
    !> results
    real(kind=rk), intent(out)                :: res(:)
    !> nVals to get
    integer, intent(in)                       :: nVals
    ! ---------------------------------------------------------------------------
    integer                                   :: iVal
    integer                                   :: nCompPDF
    ! ---------------------------------------------------------------------------
    nCompPDF = varSys%method%val(fun%input_varPos(1))%nComponents

    res = 0.0_rk
    do iVal = 1, nVals
      res(iVal) = sum( pdf( (iVal-1)*nCompPDF + 1: iVal*nCompPDF) )
    end do !iVal
  end subroutine mus_deriveMoleDensity
! ****************************************************************************** !


  ! ************************************************************************** !
  !> Calculate the potential of a given set of elements (sum up all links).
  !! This routine is used to compute potential for all scheme kinds
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine deriveMoleDensity_forElement(fun, varsys, elempos,      &
    &                                               time, tree, nElems, nDofs, &
    &                                               res                        )
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

    fnCalcPtr => mus_deriveMoleDensity

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

  end subroutine deriveMoleDensity_forElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate the potential of a given set of elements (sum up all links).
  !! This routine is used to compute potential for all scheme kinds
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValofIndex]].
  !!
  recursive subroutine deriveMoleDensity_fromIndex(fun, varSys, time, iLevel, &
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer  :: fnCalcPtr
    ! -------------------------------------------------------------------- !

    fnCalcPtr => mus_deriveMoleDensity

    call mus_generic_varFromPDF_fromIndex( &
      &  fun       = fun,                  &
      &  varSys    = varSys,               &
      &  time      = time,                 &
      &  iLevel    = iLevel,               &
      &  idx       = idx,                  &
      &  nVals     = nVals,                &
      &  fnCalcPtr = fnCalcPtr,            &
      &  res       = res                   )

  end subroutine deriveMoleDensity_fromIndex
  ! ************************************************************************** !


! ****************************************************************************** !
  !> Update state with source variable "electric field"
  !!
  !! $$ \S_j = \w_j*\c_j*\S $$
  !!
  !! Where \S is the source term
  !! S = (D*w/cs^2)*(z_i*F/(R*T))*n_i.E
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_electricFieldNP( fun, inState, outState, neigh,      &
    &                                  auxField, nPdfSize, iLevel, varSys, &
    &                                  time, phyConvFac, derVarPos         )
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
    real(kind=rk) :: electricField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: EF_elem(3)
    integer :: iElem, nElems, iDir, posInTotal, elemOff
    integer :: iField, nFields, depField, nScalars, QQ, nInputStates
    !number density of nSpecies
    real(kind=rk) :: num_dens( varSys%nStateVars )
    real(kind=rk) :: gasConstLB
    real(kind=rk), dimension(3, varSys%nStateVars ) :: electricTerm
    real(kind=rk), dimension(varSys%nStateVars) :: omega_fac
    real(kind=rk) :: force_fac, faradayLB
    real(kind=rk) :: forceTerm
    ! -------------------------------------------------------------------- !
!write(dbgUnit(1),*) 'source variable: ', trim(varSys%varname%val(fun%srcTerm_varPos))
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    associate( scheme => fPtr%solverData%scheme,                            &
      &        auxField => fPtr%solverData%scheme%auxField,                 &
      &        nernstPlanck => fPtr%solverData%scheme%nernstPlanck,         &
      &        physics => fPtr%solverData%physics,                          &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )


      ! Get electrical force which is refered in config file either its
      ! spacetime variable or operation variable
      call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
        & varSys  = varSys,                                   &
        & time    = time,                                     &
        & iLevel  = iLevel,                                   &
        & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
        & nVals   = nElems,                                   &
        & res     = electricField                             )

        ! convert physical to lattice
      electricField = electricField * physics%coulomb0 &
        &           / physics%fac(iLevel)%force

      ! number of pdf states this source depends on
      ! last input is spacetime function so it is neglected
      nInputStates =  varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

      ! constant parameter
      nFields = scheme%nFields
      QQ = scheme%layout%fStencil%QQ
      nScalars = varSys%nScalars

      gasConstLB = gasConst_R / physics%fac(iLevel)%gasConst
      faradayLB = faraday / physics%fac(iLevel)%faraday

      ! factor to multiply force Term
      force_fac = faradayLB / ( gasConstLB * nernstPlanck%temp )
      ! omega factor for each species
      do iField = 1, nFields
        omega_fac(iField) = (1.0_rk - species(iField)%omega*0.5_rk)
      end do

      !write(*,*) 'force_fac ', force_fac
      ! update source for each element
      do iElem = 1, nElems
        posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)
        ! element offset for auxField
        elemoff = (posInTotal-1)*varSys%nAuxScalars

        ! Get number density from auxField
        do iField = 1, nFields
          num_dens(iField) = auxField(iLevel)%val( elemOff + iField )
        end do

        ! convert physical to lattice
        EF_elem = electricField((iElem-1)*3+1 : iElem*3)

        ! S = (1-omega_i/2)*(z_i*F/(R*T))*n_i.E
        do iField = 1, nFields
          electricTerm(:, iField) = force_fac * omega_fac(iField) &
            &                     * species(iField)%chargeNr      &
            &                     * num_dens(iField) * EF_elem(:)
        end do

        !write(*,*) 'electricTerm ', electricTerm

        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)

          do iDir = 1, QQ
            ! Force on each species
            ! d^m_k = weight_m*cxDir_m ( S )
            forceTerm = scheme%layout%fStencil%cxDirRK( 1, iDir ) &
              &         * electricTerm(1, depField)               &
              &       + scheme%layout%fStencil%cxDirRK( 2, iDir ) &
              &         * electricTerm(2, depField)               &
              &       + scheme%layout%fStencil%cxDirRK( 3, iDir ) &
              &         * electricTerm(3, depField)

            outState(                                                         &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & = outState(                                                   &
              & (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
              & + scheme%layout%weight( iDir ) * forceTerm
          end do ! iDir
        end do !iField
      end do !iElem
    end associate

  end subroutine applySrc_electricFieldNP
! ****************************************************************************** !



end module mus_derQuanNernstPlanck_module
