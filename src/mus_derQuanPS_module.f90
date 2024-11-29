! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!! author: Jiaxing Qi
!! This module provides the MUSUBI specific functions for calculating
!! macroscopic quantities from the state variables.
!!
!! The depending common interface between MUSUBI and ATELES is defined in the
!! tem_derived_module. The functionality for accessing a variable from the state
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
module mus_derQuanPS_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use tem_param_module,         only: div1_2, div1_3, div1_54, div1_9, div3_4, &
    &                                 sqrt3, cs2inv, cs2, t2cs2inv, t2cs4inv,  &
    &                                 cs4inv, q000
  use env_module,               only: rk, long_k, labelLen
  use tem_float_module,         only: operator(.feq.), operator(.fge.), &
    &                                 operator(.fle.)
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
  use tem_operation_var_module, only: tem_evalMag_forElement,     &
    &                                 tem_evalMag_forPoint,       &
    &                                 tem_evalMag_fromIndex,      &
    &                                 tem_opVar_setupIndices,     &
    &                                 tem_get_new_varSys_data_ptr
  use tem_property_module,      only: prp_hasBnd, prp_hasQval
  use tem_tools_module,         only: tem_PositionInSorted
  use tem_debug_module,         only: dbgUnit
  use tem_grow_array_module,    only: grw_labelarray_type, append

  ! include musubi modules
  use mus_source_type_module,        only: mus_source_op_type
  use mus_pdf_module,                only: pdf_data_type
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_varSys_module,             only: mus_varSys_data_type,             &
    &                                      mus_varSys_solverData_type,       &
    &                                      mus_get_new_solver_ptr,           &
    &                                      mus_deriveVar_forPoint,           &
    &                                      mus_generic_varFromPDF_fromIndex, &
    &                                      mus_generic_fromPDF_forElement,   &
    &                                      mus_derive_fromPDF
  use mus_operation_var_module,      only: mus_opVar_setupIndices
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_derQuan_module,            only: derivePressure,        &
    &                                      derivePressure_fromIndex
  use mus_physics_module,            only: mus_convertFac_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none

  private

  public :: mus_append_derVar_lbmPS
  public :: deriveEquilPS_FromMacro
  public :: deriveEquilPS2ndOrder_FromMacro
  public :: deriveEquilPS_fromAux
  public :: deriveAuxPS_fromState

  ! source variable
  public :: derive_injectionPS
  public :: derive_equalInjectionPS

  ! source update
  public :: applySrc_injectionPS
  public :: applySrc_equalInjectionPS

contains

  ! **************************************************************************** !
  !> subroutine to add derive variables for weakly compressible LBM
  !! (schemekind = 'passive_scalar') to the varsys.
  !! for passive scalar contains only one derive variable:
  !! density
  subroutine mus_append_derVar_lbmPS( varSys, solverData, fldLabel, derVarName )
    ! ---------------------------------------------------------------------------
    !> global variable system
    type(tem_varSys_type), intent(inout)  :: varSys

    !> Contains pointer to solver data types
    type(mus_varSys_solverData_type), target, intent(in) :: solverData

    !> array of field label prefix. Size=nFields
    character(len=*), intent(in)              :: fldLabel

    !> array of derive physical variables
    type(grw_labelarray_type), intent(inout) :: derVarName
    ! ---------------------------------------------------------------------------
    ! number of derive variables
    integer :: addedPos
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
    character(len=labelLen) :: derVarName_loc
    ! ---------------------------------------------------------------------------
    nullify(get_point, get_element, set_params, get_params, setup_indices, &
      &     get_valOfIndex)

    derVarName_loc = 'pressure'

    call append(derVarName, derVarName_loc)

    ! set pointers for pressure variable
    get_element => derivePressure
    get_point => mus_deriveVar_forPoint
    setup_indices => mus_opVar_setupIndices
    get_valOfIndex => derivePressure_fromIndex
    set_params => tem_varSys_setParams_dummy
    get_params => tem_varSys_getParams_dummy

    ! update variable names with field label
    varname = trim(fldLabel)//trim(adjustl(derVarName_loc))
    allocate(input_varname(1))
    input_varname(1) = trim(fldLabel)//'density'

    ! append variable to varSys
    call tem_varSys_append_derVar(                            &
      &  me             = varSys,                             &
      &  varName        = trim(varname),                      &
      &  nComponents    = 1,                                  &
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
      write(logUnit(10),*) 'Appended variable: '//trim(varname)
    else if (addedpos < 1) then
      write(logUnit(1),*) 'Error: variable '//trim(varname)// &
        &                 ' is not added to variable system'
      call tem_abort()
    end if
    deallocate(input_varname)

  end subroutine mus_append_derVar_lbmPS
  ! **************************************************************************** !


! ****************************************************************************** !
!       Subroutines with common interface for the function pointers            !
! ****************************************************************************** !

! ****************************************************************************** !
  !> This routine computes equilbrium from density and velocity
  !! This must comply with mus_variable_module%derive_FromMacro
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromMacro]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:equilFromMacro]] function pointer.
  subroutine deriveEquilPS_FromMacro( density, velocity, iField, nElems, &
    &                                 varSys, layout, res                )
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
    real(kind=rk) :: fEq(layout%fStencil%QQ), vel(3), ucx
    integer :: QQ, iElem, iDir
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    do iElem = 1, nElems
      vel = velocity(:,iElem)

      do iDir = 1, QQ
        ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), vel )

        ! calculate equilibrium
        fEq( iDir ) = layout%weight( iDir ) * density(iElem) &
          & * ( 1._rk +  cs2inv * ucx )
      enddo

      res( (iElem-1)*QQ+1: iElem*QQ ) = fEq
    end do
  end subroutine deriveEquilPS_FromMacro
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine computes 2nd order equilbrium from density and velocity
  !! This must comply with mus_variable_module%derive_FromMacro
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_FromMacro]] in derived/[[mus_derVarPos_module]].f90 in order to be
  !! callable via [[mus_derVarPos_type:equilFromMacro]] function pointer.
  subroutine deriveEquilPS2ndOrder_FromMacro( density, velocity, iField, nElems, &
    &                                 varSys, layout, res                )
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
    real(kind=rk) :: fEq(layout%fStencil%QQ), vel(3), ucx, usq
    integer :: QQ, iElem, iDir
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    do iElem = 1, nElems
      vel = velocity(:,iElem)

      do iDir = 1, QQ
        ucx = dot_product( layout%fStencil%cxDirRK(:, iDir), vel )
        usq = dot_product( vel, vel )

        ! calculate equilibrium
        fEq( iDir ) = layout%weight( iDir ) * density(iElem) &
          & * ( 1._rk + cs2inv * ucx + cs2inv * cs2inv * ucx * ucx &
          & * 0.5_rk - usq * cs2inv * 0.5_rk )

      end do

      res( (iElem-1)*QQ+1: iElem*QQ ) = fEq
    end do
  end subroutine deriveEquilPS2ndOrder_FromMacro
! ****************************************************************************** !


! **************************************************************************** !
  !> This routine computes auxField from state array
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_auxFromState]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:auxFieldFromState]] function pointer.
  subroutine deriveAuxPS_fromState( derVarPos, state, neigh, iField, nElems, &
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

  end subroutine deriveAuxPS_fromState
! **************************************************************************** !

  ! ************************************************************************** !
  !> This routine computes equilbrium from auxField
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[derive_equilFromAux]] in derived/[[mus_derVarPos_module]].f90 in order to
  !! be callable via [[mus_derVarPos_type:equilFromAux]] function pointer.
  subroutine deriveEquilPS_fromAux( derVarPos, auxField, iField, nElems, &
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
    write(logUnit(1),*) 'ERROR: Equilibrium calculation requires transport '
    write(logUnit(1),*) 'velocity and it is not provided to this routine'
    write(logUnit(1),*) 'Solution: use deriveEquilPS_fromMacro'
    call tem_abort()
  end subroutine deriveEquilPS_fromAux
  ! ************************************************************************** !

! ****************************************************************************** !
   !> Derive injection variable defined as a source term for lbm passive scalar.
   !! It evaluates spacetime function defined in lua file for injection variable
   !! and convert it to state value which is to be added to the state
  recursive subroutine derive_injectionPS(fun, varsys, elempos, time, tree, &
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
    type(mus_scheme_type), pointer :: scheme
    real(kind=rk) :: luaRes(nElems)
    real(kind=rk) :: velocity(nElems*3)
    real(kind=rk) :: local_vel(3)
    real(kind=rk), allocatable :: linkwise_luaRes(:)
    integer :: trans_velPos, data_varPos
    integer :: iElem, iDir
    integer :: nCompPDF, nInputStates
    ! -------------------------------------------------------------------- !
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( fun%method_data, fPtr )
    scheme => fPtr%solverData%scheme

    ! number of pdf states this source depends on
    ! last input is spacetime function so it is neglected
    nInputStates = fun%nInputs - 1
    ! constant parameter
    nCompPDF = scheme%layout%fStencil%QQ
    allocate(linkwise_luaRes(nCompPDF))

    ! last input var is spacetime function
    data_varPos = fun%input_varPos(fun%nInputs)
    ! Get injection which might be defined to vary in space and time
    call varSys%method%val(data_varPos)%get_element( varSys  = varSys,  &
      &                                              elemPos = elemPos, &
      &                                              time    = time,    &
      &                                              tree    = tree,    &
      &                                              nElems  = nElems,  &
      &                                              ndofs   = 1,       &
      &                                              res     = luaRes   )

    ! Get transport velocity field from stFun
    trans_velPos = scheme%transVar%method(1)%data_varPos
    call varSys%method%val(trans_velPos)%get_element( varSys  = varSys,  &
      &                                               elemPos = elemPos, &
      &                                               time    = time,    &
      &                                               tree    = tree,    &
      &                                               nElems  = nElems,  &
      &                                               nDofs   = 1,       &
      &                                               res     = velocity )

    do iElem = 1, nElems

      ! velocity from dependent scheme
      local_vel = velocity((iElem-1)*3+1 : iElem*3)

      do iDir = 1, nCompPDf
        linkwise_luaRes(iDir) =                                              &
          & luaRes(iElem) /  6._rk * (1._rk + 3._rk * (                      &
          & dot_product(                                                     &
          & scheme%layout%fStencil%cxDirRK(:,iDir), local_vel(:)) ) )
      end do

      ! Update souce depends on nInputStates
      ! if nInputStates = 1, it is field source else it is global source
      res( (iElem-1)*fun%nComponents + 1 : iElem*fun%nComponents) &
        & = linkwise_luaRes

    end do !iElem

  end subroutine derive_injectionPS
! ****************************************************************************** !


! ****************************************************************************** !
   !> Derive injection variable defined as a source term for lbm passive scalar.
   !! It evaluates spacetime function defined in lua file for injection variable
   !! and convert it to state value which is to be added to the state
  recursive subroutine derive_equalInjectionPS(fun, varsys, elempos, time, &
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
    real(kind=rk) :: luaRes(nElems)
    real(kind=rk), allocatable :: linkwise_luaRes(:)
    integer :: iElem, iDir, stFun_pos
    integer :: iField, nCompPDF, nInputStates
    real(kind=rk) :: srcTerm(fun%nComponents)
    ! -------------------------------------------------------------------- !
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( fun%method_data, fPtr )

    ! number of pdf states this source depends on
    ! last input is spacetime function so it is neglected
    nInputStates = fun%nInputs - 1
    ! constant parameter
    nCompPDF = fPtr%solverData%scheme%layout%fStencil%QQ
    allocate(linkwise_luaRes(nCompPDF))

    ! last input var is spacetime function
    stFun_pos = fun%input_varPos(fun%nInputs)

    ! Get injection which might be defined to vary in space and time
    call varSys%method%val(stfun_pos)%get_element( varSys  = varSys,  &
      &                                            elemPos = elemPos, &
      &                                            time    = time,    &
      &                                            tree    = tree,    &
      &                                            ndofs   = 1,       &
      &                                            nElems  = nElems,  &
      &                                            res     = luaRes   )

    ! loop for all elements
    do iElem = 1, nElems

      ! compute the linkwise equal distribution
      do iDir = 1, nCompPDf
        linkwise_luaRes(iDir) = luaRes(iElem) /  real(nCompPDF, kind=rk)
      end do

      ! Update souce depends on nInputStates
      ! if nInputStates = 1, it is field source else it is global source
      do iField = 1, nInputStates
        ! loop over all directions
        do iDir = 1, nCompPDF
          srcTerm((iField-1)*nCompPDF + iDir) = linkwise_luaRes(iDir)
        end do
      end do
    end do !iElem

  end subroutine derive_equalInjectionPS
! ****************************************************************************** !

! ****************************************************************************** !
  !> Update state with source variable "injection"
  !! Similar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_injectionPS( fun, inState, outState, neigh, auxField,    &
    &                              nPdfSize, iLevel, varSys, time, phyConvFac, &
    &                              derVarPos                                   )
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
    real(kind=rk) :: luaRes(fun%elemLvl(iLevel)%nElems)
    real(kind=rk) :: transVel(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: local_vel(3)
    real(kind=rk), allocatable :: linkwise_luaRes(:)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, depField, nScalars, QQ, nInputStates
    integer :: trans_varPos
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )
    scheme => fPtr%solverData%scheme

    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! number of pdf states this source depends on
    ! last input is spacetime function so it is neglected
    nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1
    ! constant parameter
    QQ = scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars
    allocate(linkwise_luaRes(QQ))

    ! Get injection which might be defined to vary in space and time
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = luaRes                                    )

    ! Get transport velocity field from stFun
    ! Setup indices for transport variable is set for nSolve elements
    ! in total list so to get only velocity on active source space,
    ! use idx of source on transVar%index
    trans_varPos = scheme%transVar%method(1)%data_varPos
    call varSys%method%val(trans_varPos)%get_valOfIndex(   &
      & varSys  = varSys,                                  &
      & time    = time,                                    &
      & iLevel  = iLevel,                                  &
      & idx     = scheme%transVar%method(1)%pntIndex       &
      &           %indexLvl(iLevel)                        &
      &           %val(fun%elemLvl(iLevel)%idx(1:nElems)), &
      & nVals   = nElems,                                  &
      & res     = transVel                                 )

    ! convert physical to lattice
    transVel = transVel / fPtr%solverData%physics%fac(iLevel)%vel

    ! now update the state vector for all elements
    do iElem = 1, nElems

      ! position of source element in levelwise state array
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! velocity from dependent scheme
      local_vel = transVel((iElem-1)*3+1 : iElem*3)

      do iDir = 1, QQ
        linkwise_luaRes(iDir) =                                              &
          & luaRes(iElem) /  6._rk * (1._rk + 3._rk * (                      &
          & dot_product(                                                     &
          & fPtr%solverData%scheme%layout%fStencil%cxDirRK(:,iDir),          &
          & local_vel(:)) ) )
      end do

      ! Update souce depends on nInputStates
      ! if nInputStates = 1, it is field source else it is global source
      do iField = 1, nInputStates
        depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)
        ! loop over all directions
        do iDir = 1, QQ
          ! add the linkwise lua res for the current element
          outState( (posintotal-1)*nscalars+idir+(depfield-1)*qq )&
            & = outState(                                                      &
            & (posintotal-1)*nscalars+idir+(depfield-1)*qq )      &
            & + linkwise_luaRes(iDir)
        end do
      end do
    end do !iElem
  end subroutine applySrc_injectionPS
! ****************************************************************************** !


! ****************************************************************************** !
  !> Update state with source variable "equalInjection"
  !! Simuilar to derive routine but it updates the state whereas derive
  !! is used for tracking
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[proc_apply_source]] in derived/[[mus_source_type_module]].f90 in order to
  !! be callable via [[mus_source_op_type:applySrc]] function pointer.
  subroutine applySrc_equalInjectionPS( fun, inState, outState, neigh,      &
    &                                   auxField, nPdfSize, iLevel, varSys, &
    &                                   time, phyConvFac, derVarPos         )
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
    real(kind=rk) :: luaRes(fun%elemLvl(iLevel)%nElems)
    real(kind=rk), allocatable :: linkwise_luaRes(:)
    integer :: iElem, nElems, iDir, posInTotal
    integer :: iField, depField, nScalars, QQ, nInputStates
    ! ---------------------------------------------------------------------------
    ! convert c pointer to solver type fortran pointer
    call c_f_pointer( varSys%method%val( fun%srcTerm_varPos )%method_data, &
      &               fPtr )

    nElems = fun%elemLvl(iLevel)%nElems
    ! number of pdf states this source depends on
    ! last input is spacetime function so it is neglected
    nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1
    ! constant parameter
    QQ = fPtr%solverData%scheme%layout%fStencil%QQ
    nScalars = varSys%nScalars
    allocate(linkwise_luaRes(QQ))

    ! Get injection which might be defined to vary in space and time
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = luaRes                                    )

    ! now update the state vector for all elements
    do iElem = 1, nElems

      ! position of source element in levelwise state array
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)

      ! compute the linkwise equal distribution
      do iDir = 1, QQ
        linkwise_luaRes(iDir) = luaRes(iElem) /  real(QQ, kind=rk)
      end do

      ! Update souce depends on nInputStates
      ! if nInputStates = 1, it is field source else it is global source
      do iField = 1, nInputStates
        depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)
        ! loop over all directions
        do iDir = 1, QQ
          ! add the linkwise lua res for the current element
          outState( (posintotal-1)*nscalars+idir+(depfield-1)*qq ) &
            & = outState(                                                  &
            & (posintotal-1)*nscalars+idir+(depfield-1)*qq )       &
            & + linkwise_luaRes(iDir)
        end do
      end do
    end do !iElem
  end subroutine applySrc_equalInjectionPS
! ****************************************************************************** !


end module mus_derQuanPS_module
