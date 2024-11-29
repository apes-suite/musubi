! Copyright (c) 2016 Philipp Otte <otte@mathcces.rwth-aachen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016-2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2017 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> author: Philipp Otte
!! This module provides the MUSUBI specific functions for calculating
!! macroscopic quantities for the isothermal acoustic equations
!! from the state variables
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
module mus_derQuanIsothermAcEq_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use tem_param_module,         only: div1_2, div1_3, div1_54, div1_9, div3_4, &
    &                                 sqrt3, cs2inv, cs2, t2cs2inv, t2cs4inv,  &
    &                                 cs4inv, rho0, rho0Inv
  use env_module,               only: rk, long_k, labelLen
  use tem_variable_module,      only: tem_variable_type
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_topology_module,      only: tem_levelOf
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_logging_module,       only: logUnit
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
  use tem_aux_module,           only: tem_abort
  use tem_property_module,      only: prp_hasBnd, prp_hasQval
  use tem_tools_module,         only: tem_PositionInSorted
  use tem_grow_array_module,    only: grw_labelarray_type, append

  ! include musubi modules
  use mus_source_type_module,        only: mus_source_op_type
  use mus_pdf_module,                only: pdf_data_type
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_varSys_module,             only: mus_varSys_data_type,               &
    &                                      mus_varSys_solverData_type,         &
    &                                      mus_get_new_solver_ptr,             &
    &                                      mus_generic_varFromPDF_fromIndex,   &
    &                                      mus_generic_fromPDF_forElement,     &
    &                                      mus_deriveVar_ForPoint,             &
    &                                      mus_derive_fromPDF
  use mus_operation_var_module,      only: mus_opVar_setupIndices
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none

  private

  public :: mus_append_derVar_isotherm_acEq
  public :: derDensityIsothermAcEq
  public :: derVelocityIsothermAcEq
  public :: derPressureIsothermAcEq
  public :: derEquilIsothermAcEq
  public :: deriveEquil_FromMacro_IsothermAcEq
  public :: deriveEquilIsoThermAcEq_fromAux
  public :: deriveEq_FromState_IsothermAcEq
  public :: deriveVelocity_FromState_IsothermAcEq

  ! =============================================================================
  ! D3Q19 flow model
  ! =============================================================================
  !> Definition of the discrete velocity set

  ! integer,parameter :: block = 32
  integer,parameter :: QQ   = 19  !< number of pdf directions

  integer,parameter :: qN00 = 1   !< west             x-
  integer,parameter :: q0N0 = 2   !< south            y-
  integer,parameter :: q00N = 3   !< bottom           z-
  integer,parameter :: q100 = 4   !< east             x+
  integer,parameter :: q010 = 5   !< north            y+
  integer,parameter :: q001 = 6   !< top              z+
  integer,parameter :: q0NN = 7   !<                  z-,y-
  integer,parameter :: q0N1 = 8   !<                  z+,y-
  integer,parameter :: q01N = 9   !<                  z-,y+
  integer,parameter :: q011 = 10  !<                  z+,y+
  integer,parameter :: qN0N = 11  !<                  x-,z-
  integer,parameter :: q10N = 12  !<                  x+,z-
  integer,parameter :: qN01 = 13  !<                  x-,z+
  integer,parameter :: q101 = 14  !<                  x+,z+
  integer,parameter :: qNN0 = 15  !<                  y-,x-
  integer,parameter :: qN10 = 16  !<                  y+,x-
  integer,parameter :: q1N0 = 17  !<                  y-,x+
  integer,parameter :: q110 = 18  !<                  y+,x+
  integer,parameter :: q000 = 19  !< rest density is last

  real(kind=rk), parameter :: f1 = 2.0_rk / 5.0_rk
  real(kind=rk), parameter :: f2 = 1.0_rk / 30.0_rk
  real(kind=rk), parameter :: f8 = 1.0_rk / 30.0_rk

contains

  ! **************************************************************************** !
  !> subroutine to add derive variables for isothermal acoustic
  !! equations
  !! (schemekind = 'isotherm_acEq') to the varsys.
  subroutine mus_append_derVar_isotherm_acEq( varSys, solverData, schemeHeader, &
                                    stencil, fldLabel, derVarName )
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)  :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> identifier of the scheme
    type(mus_scheme_header_type), intent(in)  :: schemeHeader

    !> compute stencil defintion
    type(tem_stencilHeader_type), intent(in)  :: stencil

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel

    !> array of derive physical variables
    type(grw_labelarray_type), intent(inout) :: derVarName
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
    character(len=labelLen), allocatable :: derVarName_loc(:)
    ! ---------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    nDerVars = 2
    allocate(derVarName_loc(nDerVars))
    derVarName_loc    = [ 'pressure       ', 'equilibrium    '  ]

    do iVar = 1, nDerVars
      call append(derVarName, derVarName_loc(iVar))

      ! set default pointers, overwrite if neccessary
      get_element => tem_varSys_getElement_dummy
      get_point => mus_deriveVar_ForPoint
      setup_indices => mus_opVar_setupIndices
      get_valOfIndex => tem_varSys_getValOfIndex_dummy
      set_params => tem_varSys_setParams_dummy
      get_params => tem_varSys_getParams_dummy

      select case(trim(adjustl(derVarName_loc(iVar))))
      case ('pressure')
        get_element => derPressureIsothermAcEq
        nComponents = 1
        allocate(input_varname(1))
        input_varname(1) = 'density'

      case ('equilibrium')
        get_element => derEquilIsothermAcEq
        get_valOfIndex => derEquilIsothermAcEq_fromIndex
        nComponents = stencil%QQ
        allocate(input_varname(1))
        input_varname(1) = 'pdf'

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
      call tem_varSys_append_derVar(                            &
        &  me             = varSys,                             &
        &  varName        = trim(varname),                      &
        &  nComponents    = nComponents,                        &
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
        write(logUnit(10),*) ' Appended variable:'//trim(varname)
      else if (addedpos < 1) then
        write(logUnit(1),*) 'Error: variable '//trim(varname)// &
          &                 ' is not added to variable system'
      end if

      deallocate(input_varname)
    end do

  end subroutine mus_append_derVar_isotherm_acEq
  ! **************************************************************************** !


! ****************************************************************************** !
!       Subroutines with common interface for the function pointers            !
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the density of a given set of elements (sum up all links).
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine derDensityIsothermAcEq(fun, varsys, elempos, time, &
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
    integer :: statePos, iElem, iComp, iLevel, iDep
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: pdfPos, nCompPDF
    real(kind=rk) :: dens
    integer :: nSize
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    ! res is always AOS layout
    res = 0.0_rk

    associate( state => fPtr%solverData%scheme%state,                &
      &        pdf => fPtr%solverData%scheme%pdf,                    &
      &        QQ => fPtr%solverData%scheme%layout%fStencil%QQ,      &
      &        levelPointer => fPtr%solverData%geometry%levelPointer )
    do iElem = 1, nElems
      ! if state array is defined level wise then use levelPointer(pos)
      ! to access state array
      statePos = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )
      iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )
      nSize = fPtr%solverData%scheme%pdf( iLevel )%nSize

      dens = 0.0_rk
      ! use iDep here to use this routine to compute also mixture density
      ! in multi-species
      do iDep = 1, fun%nInputs
        pdfPos = fun%input_varPos(iDep)
        nCompPDF = varSys%method%val(pdfPos)%nComponents
        do iComp = 1, nCompPDF
          dens = dens +                                                      &
            & fPtr%solverData%scheme%state( iLevel )%val(                    &
& ( statepos-1)* varsys%nscalars+icomp+( 1-1)* qq, &
            & pdf( iLevel )%nNext )
        end do !iComp
      end do !iDep
      res( iElem ) = dens

    end do ! iElem
    end associate

  end subroutine derDensityIsothermAcEq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the pressure of a given set of elements (sum up all links).
  !!
  !! Pressure calculation according to the isentropic equation of state for
  !! the LBM \( p = \rho c_s^2 \)
  !! with the calculation of density as in deriveDensity
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine derPressureIsothermAcEq(fun, varsys, elempos, time, &
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
    integer :: dens_pos
    ! -------------------------------------------------------------------- !

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

    ! convert density to pressure
    res = res  * cs2

  end subroutine derPressureIsothermAcEq
! ****************************************************************************** !

! ****************************************************************************** !
  !> Initiates the calculation of velocity
  !! This routine sets the function Pointer for velocity calcualtion and calls
  !! the generice get Element from PDF routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine derVelocityIsothermAcEq(fun, varsys, elempos, time, &
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer :: fnCalcPtr
    ! -------------------------------------------------------------------- !
    fnCalcPtr => mus_derVelocityIsothermAcEq

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

  end subroutine derVelocityIsothermAcEq
! ****************************************************************************** !

! ****************************************************************************** !
  !> Calculate the equlibrium of given elements with the given input state
  !! array.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine derEquilIsothermAcEq(fun, varsys, elempos, time, tree, &
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
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer :: fnCalcPtr
    ! -------------------------------------------------------------------- !
    fnCalcPtr => mus_derEquilIsothermAcEq

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

  end subroutine derEquilIsothermAcEq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initiates the calculation of equilibrium.
  !! This routine sets the function Pointer for equilibrium calcualtion and calls
  !! the generice get Value of Index routine
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_getValOfIndex]].
  !!
  recursive subroutine derEquilIsothermAcEq_fromIndex( fun, varSys, time,   &
     &                                         iLevel, idx, idxLen, nVals, res )
      !> Description of the method to obtain the variables, here some preset
      !! values might be stored, like the space time function to use or the
      !! required variables.
      class(tem_varSys_op_type), intent(in) :: fun

      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in)     :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)       :: time

      !> Level on which values are requested
      integer, intent(in)                   :: iLevel

      !> Index of points in the growing array and variable val array to
      !! return.
      !! Size: most times nVals, if contiguous arrays are used it depends
      !! on the number of first indices
      integer, intent(in)                   :: idx(:)

      !> With idx as start index in contiguous memory,
      !! idxLength defines length of each contiguous memory
      !! Size: dependes on number of first index for contiguous array,
      !! but the sum of all idxLen is equal to nVals
      integer, optional, intent(in)         :: idxLen(:)

      !> Number of values to obtain for this variable (vectorized access).
      integer, intent(in)                   :: nVals

      !> Resulting values for the requested variable.
      !!
      !! Dimension: n requested entries x nComponents of this variable
      !! Access: (iElem-1)*fun%nComponents + iComp
      real(kind=rk), intent(out)            :: res(:)
    ! ---------------------------------------------------------------------------
    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer  :: fnCalcPtr
    ! ---------------------------------------------------------------------------
    fnCalcPtr => mus_derEquilIsothermAcEq

    call mus_generic_varFromPDF_fromIndex( &
      &  fun       = fun,                  &
      &  varSys    = varSys,               &
      &  time      = time,                 &
      &  iLevel    = iLevel,               &
      &  idx       = idx,                  &
      &  nVals     = nVals,                &
      &  fnCalcPtr = fnCalcPtr,            &
      &  res       = res                   )

  end subroutine derEquilIsothermAcEq_fromIndex
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
  subroutine deriveVelocity_FromState_IsothermAcEq( state, iField, nElems, &
    &                                               varSys, layout, res    )
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
    integer :: iElem, iDir
    real(kind=rk) :: pdf(layout%fStencil%QQ)
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    do iElem = 1, nElems
      do iDir = 1, layout%fStencil%QQ
        pdf(iDir) = state( iDir+(iElem-1)*varSys%nScalars )
      end do
      res( (iElem-1)*3+1 ) = sum( pdf * layout%fStencil%cxDirRK(1,:) ) * rho0Inv
      res( (iElem-1)*3+2 ) = sum( pdf * layout%fStencil%cxDirRK(2,:) ) * rho0Inv
      res( (iElem-1)*3+3 ) = sum( pdf * layout%fStencil%cxDirRK(3,:) ) * rho0Inv
    end do

  end subroutine deriveVelocity_FromState_IsothermAcEq
! ****************************************************************************** !

! ****************************************************************************** !
  !> This routine computes equil from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromState]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:velFromState]],
  !! [[mus_derVarPos_type:equilFromState]],
  !! [[mus_derVarPos_type:momFromState]],
  !! [[mus_derVarPos_type:velocitiesFromState]], and
  !! [[mus_derVarPos_type:momentaFromState]] function pointers.
  subroutine deriveEq_FromState_IsothermAcEq( state, iField, nElems, varSys, &
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
    real(kind=rk) :: pdf( layout%fStencil%QQ )
    real(kind=rk) :: fEq(layout%fStencil%QQ)
    real(kind=rk) :: rho, u_x, u_y, u_z
    integer :: QQ, iElem, iDir
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    QQ = layout%fStencil%QQ

    if(trim(layout%fStencil%label) == 'd3q19' ) then
      do iElem = 1, nElems
        do iDir = 1, QQ
          pdf(iDir) = state( iDir+(iElem-1)*varSys%nScalars )
        end do

        rho = sum( pdf )
        u_x = sum( pdf * layout%fStencil%cxDirRK(1,:) ) / rho0
        u_y = sum( pdf * layout%fStencil%cxDirRK(2,:) ) / rho0
        u_z = sum( pdf * layout%fStencil%cxDirRK(3,:) ) / rho0

        ! set equilibria
        fEq( qN00 ) = ( rho - 3 * rho0 * u_x ) * f2
        fEq( q0N0 ) = ( rho - 3 * rho0 * u_y ) * f2
        fEq( q00N ) = ( rho - 3 * rho0 * u_z ) * f2
        fEq( q100 ) = ( rho + 3 * rho0 * u_x ) * f2
        fEq( q010 ) = ( rho + 3 * rho0 * u_y ) * f2
        fEq( q001 ) = ( rho + 3 * rho0 * u_z ) * f2
        fEq( q0NN ) = ( rho - 3 * rho0 * ( u_y + u_z ) ) * f8
        fEQ( q0N1 ) = ( rho - 3 * rho0 * ( u_y - u_z ) ) * f8
        fEq( q01N ) = ( rho + 3 * rho0 * ( u_y - u_z ) ) * f8
        fEq( q011 ) = ( rho + 3 * rho0 * ( u_y + u_z ) ) * f8
        fEq( qN0N ) = ( rho - 3 * rho0 * ( u_x + u_z ) ) * f8
        fEq( qN01 ) = ( rho - 3 * rho0 * ( u_x - u_z ) ) * f8
        fEq( q10N ) = ( rho + 3 * rho0 * ( u_x - u_z ) ) * f8
        fEq( q101 ) = ( rho + 3 * rho0 * ( u_x + u_z ) ) * f8
        fEq( qNN0 ) = ( rho - 3 * rho0 * ( u_x + u_y ) ) * f8
        fEq( qN10 ) = ( rho - 3 * rho0 * ( u_x - u_y ) ) * f8
        fEq( q1N0 ) = ( rho + 3 * rho0 * ( u_x - u_y ) ) * f8
        fEq( q110 ) = ( rho + 3 * rho0 * ( u_x + u_y ) ) * f8
        fEq( q000 ) = ( rho ) * f1

        res( (iElem-1)*QQ+1: iElem*QQ ) = fEq
      end do
    else
      write(logUnit(1),*) 'stencil not supported in', &
      &                   'deriveEquil_FromState_IsothermAcEq'
    end if
  end subroutine deriveEq_FromState_IsothermAcEq
! ****************************************************************************** !

! ****************************************************************************** !
  !> This routine computes equilbrium from density and velocity
  !! This must comply with interface in mus_variable_module derive_FromMacro
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromMacro]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:equilFromMacro]] function pointer.
  subroutine deriveEquil_FromMacro_IsothermAcEq( density, velocity, iField,  &
    &                                            nElems, varSys, layout, res )
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
    real(kind=rk) :: rho, u_x, u_y, u_z
    integer :: QQ, iElem
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( varSys%method%val(iField)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    QQ = layout%fStencil%QQ

    if(trim(layout%fStencil%label) == 'd3q19' ) then
      do iElem = 1, nElems
      ! get element macroscopic quantities
        rho = density( iElem )
        u_x = velocity( 1, iElem )
        u_y = velocity( 2, iElem )
        u_z = velocity( 3, iElem )

        ! set equilibria
        fEq( qN00 ) = ( rho - 3 * rho0 * u_x ) * f2
        fEq( q0N0 ) = ( rho - 3 * rho0 * u_y ) * f2
        fEq( q00N ) = ( rho - 3 * rho0 * u_z ) * f2
        fEq( q100 ) = ( rho + 3 * rho0 * u_x ) * f2
        fEq( q010 ) = ( rho + 3 * rho0 * u_y ) * f2
        fEq( q001 ) = ( rho + 3 * rho0 * u_z ) * f2
        fEq( q0NN ) = ( rho - 3 * rho0 * ( u_y + u_z ) ) * f8
        fEQ( q0N1 ) = ( rho - 3 * rho0 * ( u_y - u_z ) ) * f8
        fEq( q01N ) = ( rho + 3 * rho0 * ( u_y - u_z ) ) * f8
        fEq( q011 ) = ( rho + 3 * rho0 * ( u_y + u_z ) ) * f8
        fEq( qN0N ) = ( rho - 3 * rho0 * ( u_x + u_z ) ) * f8
        fEq( qN01 ) = ( rho - 3 * rho0 * ( u_x - u_z ) ) * f8
        fEq( q10N ) = ( rho + 3 * rho0 * ( u_x - u_z ) ) * f8
        fEq( q101 ) = ( rho + 3 * rho0 * ( u_x + u_z ) ) * f8
        fEq( qNN0 ) = ( rho - 3 * rho0 * ( u_x + u_y ) ) * f8
        fEq( qN10 ) = ( rho - 3 * rho0 * ( u_x - u_y ) ) * f8
        fEq( q1N0 ) = ( rho + 3 * rho0 * ( u_x - u_y ) ) * f8
        fEq( q110 ) = ( rho + 3 * rho0 * ( u_x + u_y ) ) * f8
        fEq( q000 ) = ( rho ) * f1

        res( (iElem-1)*QQ+1: iElem*QQ ) = fEq
      end do
    else
      write(logUnit(1),*) 'stencil not supported in', &
      &                   'deriveEquil_FromMacro_IsothermAcEq'
    end if
  end subroutine deriveEquil_FromMacro_IsothermAcEq
! ****************************************************************************** !

! ************************************************************************** !
  !> This routine computes equilbrium from auxField
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_equilFromAux]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:equilFromAux]] function pointer.
  subroutine deriveEquilIsoThermAcEq_fromAux( derVarPos, auxField, iField, &
    &                                         nElems, varSys, layout, fEq  )
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
    real(kind=rk) :: rho, u_x, u_y, u_z, fEq_loc(layout%fStencil%QQ)
    integer :: QQ, iElem, elemOff
    integer :: dens_pos, vel_pos(3)
    ! ---------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)

    QQ = layout%fStencil%QQ
    if(trim(layout%fStencil%label) == 'd3q19' ) then
      !NEC$ ivdep
      do iElem = 1, nElems
        ! element offset
        elemoff = (iElem-1)*varSys%nAuxScalars
        ! density
        rho = auxField(elemOff + dens_pos)
        ! velocity
        u_x = auxField(elemOff + vel_pos(1))
        u_y = auxField(elemOff + vel_pos(2))
        u_z = auxField(elemOff + vel_pos(3))

        ! set equilibria
        fEq_loc( qN00 ) = ( rho - 3 * rho0 * u_x ) * f2
        fEq_loc( q0N0 ) = ( rho - 3 * rho0 * u_y ) * f2
        fEq_loc( q00N ) = ( rho - 3 * rho0 * u_z ) * f2
        fEq_loc( q100 ) = ( rho + 3 * rho0 * u_x ) * f2
        fEq_loc( q010 ) = ( rho + 3 * rho0 * u_y ) * f2
        fEq_loc( q001 ) = ( rho + 3 * rho0 * u_z ) * f2
        fEq_loc( q0NN ) = ( rho - 3 * rho0 * ( u_y + u_z ) ) * f8
        fEQ_loc( q0N1 ) = ( rho - 3 * rho0 * ( u_y - u_z ) ) * f8
        fEq_loc( q01N ) = ( rho + 3 * rho0 * ( u_y - u_z ) ) * f8
        fEq_loc( q011 ) = ( rho + 3 * rho0 * ( u_y + u_z ) ) * f8
        fEq_loc( qN0N ) = ( rho - 3 * rho0 * ( u_x + u_z ) ) * f8
        fEq_loc( qN01 ) = ( rho - 3 * rho0 * ( u_x - u_z ) ) * f8
        fEq_loc( q10N ) = ( rho + 3 * rho0 * ( u_x - u_z ) ) * f8
        fEq_loc( q101 ) = ( rho + 3 * rho0 * ( u_x + u_z ) ) * f8
        fEq_loc( qNN0 ) = ( rho - 3 * rho0 * ( u_x + u_y ) ) * f8
        fEq_loc( qN10 ) = ( rho - 3 * rho0 * ( u_x - u_y ) ) * f8
        fEq_loc( q1N0 ) = ( rho + 3 * rho0 * ( u_x - u_y ) ) * f8
        fEq_loc( q110 ) = ( rho + 3 * rho0 * ( u_x + u_y ) ) * f8
        fEq_loc( q000 ) = ( rho ) * f1

        fEq( (iElem-1)*QQ+1: iElem*QQ ) = fEq_loc
      end do
    else
      call tem_abort('stencil not supported in deriveEquilIsoThermAcEq_fromAux')
    end if
  end subroutine deriveEquilIsoThermAcEq_fromAux
! ************************************************************************** !

! ****************************************************************************** !
  !> Calculate the velocity of a given element number with the given input
  !! vector (sum up all values)
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine mus_derVelocityIsothermAcEq(fun, varsys, stencil,       &
    &                                      iLevel, posInState, pdf, res, nVals )
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
    integer :: iComp, iVal
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: pdfPos, nCompsPDF
    real(kind=rk), allocatable :: tmpPDF(:)
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    pdfPos = fun%input_varPos(1)
    nCompsPDF = varSys%method%val( pdfPos )%nComponents
    allocate( tmpPDF( nCompsPDF ) )
    res = 0.0_rk

    do iVal = 1, nVals
      tmpPDF = pdf( (iVal-1)*nCompsPDF+1 : iVal*nCompsPDF )
      do iComp = 1, fun%nComponents
        res( (iVal-1)*fun%nComponents+iComp ) =                    &
          &  sum(tmpPDF * scheme%layout%fStencil%cxDirRK(iComp,:)) &
          &  * rho0Inv
      end do
    end do ! iVal

  end subroutine mus_derVelocityIsothermAcEq
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate the equlibrium of given elements with the given input state
  !! array.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!

  recursive subroutine mus_derEquilIsothermAcEq(fun, varsys, stencil, iLevel,  &
    &                                              posInState, pdf, res, nVals )
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
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: pdfPos, nCompsPDF, iVal
    real(kind=rk) :: rho, u_x, u_y, u_z
    real(kind=rk), allocatable :: tmpPDF(:)
    real(kind=rk), allocatable :: fEq(:)
    ! ---------------------------------------------------------------------------
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    pdfPos = fun%input_varPos(1)
    nCompsPDF = varSys%method%val( pdfPos )%nComponents
    allocate( tmpPDF( nCompsPDF ) )
    allocate( fEq( fun%nComponents ) )
    res = 0.0_rk

    select case( trim(scheme%header%layout) )
    case('d3q19')
      ! for d3q19
      do iVal = 1, nVals
        tmpPDF = pdf( (iVal-1)*nCompsPDF+1 : iVal*nCompsPDF )
        ! computes density and velocity
        rho   = sum(tmpPDF)
        u_x = sum(tmpPDF * scheme%layout%fStencil%cxDirRK(1,:)) &
          &    * rho0Inv
        u_y = sum(tmpPDF * scheme%layout%fStencil%cxDirRK(2,:)) &
          &    * rho0Inv
        u_z = sum(tmpPDF * scheme%layout%fStencil%cxDirRK(3,:)) &
          &    * rho0Inv

            fEq( qN00 ) = ( rho - 3 * rho0 * u_x ) * f2
            fEq( q0N0 ) = ( rho - 3 * rho0 * u_y ) * f2
            fEq( q00N ) = ( rho - 3 * rho0 * u_z ) * f2
            fEq( q100 ) = ( rho + 3 * rho0 * u_x ) * f2
            fEq( q010 ) = ( rho + 3 * rho0 * u_y ) * f2
            fEq( q001 ) = ( rho + 3 * rho0 * u_z ) * f2
            fEq( q0NN ) = ( rho - 3 * rho0 * ( u_y + u_z ) ) * f8
            fEq( q0N1 ) = ( rho - 3 * rho0 * ( u_y - u_z ) ) * f8
            fEq( q01N ) = ( rho + 3 * rho0 * ( u_y - u_z ) ) * f8
            fEq( q011 ) = ( rho + 3 * rho0 * ( u_y + u_z ) ) * f8
            fEq( qN0N ) = ( rho - 3 * rho0 * ( u_x + u_z ) ) * f8
            fEq( qN01 ) = ( rho - 3 * rho0 * ( u_x - u_z ) ) * f8
            fEq( q10N ) = ( rho + 3 * rho0 * ( u_x - u_z ) ) * f8
            fEq( q101 ) = ( rho + 3 * rho0 * ( u_x + u_z ) ) * f8
            fEq( qNN0 ) = ( rho - 3 * rho0 * ( u_x + u_y ) ) * f8
            fEq( qN10 ) = ( rho - 3 * rho0 * ( u_x - u_y ) ) * f8
            fEq( q1N0 ) = ( rho + 3 * rho0 * ( u_x - u_y ) ) * f8
            fEq( q110 ) = ( rho + 3 * rho0 * ( u_x + u_y ) ) * f8
            fEq( q000 ) = ( rho ) * f1

        res( (iVal-1)*fun%nComponents+1: iVal*fun%nComponents ) = fEq
      end do ! iVal
    case default
      write(logUnit(1),*) 'layout not supported by',            &
      & 'derEquilIsothermAcEq'
    end select

    deallocate( tmpPDF )
    deallocate( fEq )
  end subroutine mus_derEquilIsothermAcEq
! ****************************************************************************** !
end module mus_derQuanIsothermAcEq_module
! ****************************************************************************** !
