! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2021-2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
! ****************************************************************************** !
!> This module provides the routine for applying operators. Currently it is
!! only implemented for 3D and needs to be extended to 2d
module mus_operation_var_module
  use, intrinsic :: iso_c_binding,   only: c_f_pointer, c_ptr
  use env_module,                    only: rk, long_k, labellen

  use tem_logging_module,            only: logUnit
  use tem_debug_module,              only: dbgUnit
  use tem_aux_module,                only: tem_abort
  use tem_float_module,              only: operator(.fne.)
  use tem_varSys_module,             only: tem_varSys_type,                    &
    &                                      tem_varSys_op_type,                 &
    &                                      tem_varSys_proc_element,            &
    &                                      tem_varSys_solverData_evalElem_type
  use tem_time_module,               only: tem_time_type
  use treelmesh_module,              only: treelmesh_type
  use tem_geometry_module,           only: tem_CoordOfReal, &
    &                                      tem_PosofId, tem_BaryOfId
  use tem_topology_module,           only: tem_IdOfCoord, tem_levelOf, &
    &                                      tem_FirstIdAtLevel, tem_CoordOfId
  use tem_grow_array_module,         only: grw_intArray_type, append, truncate
  use tem_operation_var_module,      only: tem_opVar_fill_inputIndex, &
    &                                      tem_varSys_op_Data_type
  use tem_element_module,            only: eT_fluid


  ! include musubi modules
  use mus_scheme_type_module,  only: mus_scheme_type
  use mus_varSys_module,       only: mus_varSys_data_type,            &
    &                                mus_varSys_solverData_type,      &
    &                                mus_get_new_solver_ptr,          &
    &                                mus_derVar_intpOnPoint
  use mus_connectivity_module, only: mus_intp_getSrcElemPosInTree

  implicit none

  private

  public :: mus_opVar_setupIndices
  public :: mus_set_opVar_getElement
  public :: mus_opVar_gradU_forElement
  public :: mus_opVar_vorticity_forElement
  public :: mus_opVar_QCriterion_forElement


contains


  ! ************************************************************************** !
  !> Routine to store musubi varSys Data in operation variable solver_bundle.
  !! Unline Ateles, Musubi operations does not require any special treatment so
  !! it uses to generic routines in treelm
  subroutine mus_set_opVar_getElement( solData_evalElem, fun )
    ! --------------------------------------------------------------------------
    !> Description on how to set the element retrieval function for stfuns.
    class(tem_varSys_solverData_evalElem_type), intent(in) :: solData_evalElem

    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    type(tem_varSys_op_type), intent(inout) :: fun
    ! --------------------------------------------------------------------------
    type(tem_varSys_op_Data_type), pointer :: fptr
    type(mus_varSys_solverData_type), pointer :: fSDptr
    ! --------------------------------------------------------------------------

    write(logunit(10),*) "Setting different solver_bundle and " &
      & // " get_element routine for variable at position ",    &
      & fun%myPos
    call C_F_Pointer(fun%method_data, fptr)
    call c_f_pointer(solData_evalElem%solver_bundle, fSDptr)
    fptr%solver_bundle = mus_get_new_solver_ptr( fSDptr )

    select case (trim(fun%operType))
    case ('reduction_transient')
      fun%get_point => mus_reductionTransient_forPoint
    end select

  end subroutine mus_set_opVar_getElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Same as reductionTransient_forElement except it evaluate it multiply values
  !! from points
  !!
  !! The interface has to comply to the abstract interface
  !! tem_varSys_module#tem_varSys_proc_point.
  recursive subroutine mus_reductionTransient_forPoint( fun, varsys, point, &
    &                                                time, tree, nPnts, res )
    ! ---------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    type(tem_varSys_op_Data_type), pointer :: fptr_temOpData
    type(mus_scheme_type), pointer :: scheme
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr_temOpData )
    call C_F_POINTER(fPtr_temOpData%solver_bundle, fPtr)
    scheme => fPtr%solverData%scheme
    call mus_derVar_intpOnPoint( varPos       = fun%myPos,               &
      &                          varSys       = varSys,                  &
      &                          tree         = tree,                    &
      &                          time         = time,                    &
      &                          point        = point,                   &
      &                          nPnts        = nPnts,                   &
      &                          stencil      = scheme%layout%fStencil,  &
      &                          levelPointer = fPtr%solverData%geometry &
      &                                             %levelPointer,       &
      &                          levelDesc    = scheme%levelDesc,        &
      &                          res          = res                      )

  end subroutine mus_reductionTransient_forPoint
  ! ************************************************************************** !


  ! ************************************************************************** !
  recursive subroutine mus_opVar_setupIndices( fun, varSys, point, offset_bit,&
    &                                iLevel, tree, nPnts, idx        )
    ! -------------------------`-----------------------------------------------!
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> List of space coordinate points to store as growing array in
    !! method_data
    real(kind=rk), intent(in) :: point(:,:)

    !> Offset bit encoded as character for every point.
    !!
    !! Offset integer coord(3) is converted into a character with
    !! offset_bit = achar( (coord(1)+1) + (coord(2)+1)*4 + (coord(3)+1)*16 )
    !! Backward transformation form character to 3 integer:
    !! coord(1) = mod(ichar(offset_bit),4) - 1
    !! coord(2) = mod(ichar(offset_bit),16)/4 - 1
    !! coord(3) = ichar(offset_bit)/16 - 1
    !!
    !! If not present default is to center i.e offset_bit = achar(1+4+16)
    character, optional, intent(in) :: offset_bit(:)

    !> Level to which input points belong to
    integer, intent(in) :: iLevel

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of points to add in method_data of this variable
    integer, intent(in) :: nPnts

    !> Index of points in the growing array and variable val array.
    !! Size: nPoints
    !!
    !! This must be stored in boundary or source depends on who
    !! calls this routine.
    !! This index is required to return a value using getValOfIndex.
    integer, intent(out) :: idx(:)
    ! -------------------------------------------------------------------------!
    type(mus_varSys_data_type), pointer :: fPtr
    integer :: iPnt, iDep, nVals_prev
    type(grw_intArray_type), allocatable :: inputIndex_loc(:)
    integer, allocatable :: idxPerPnt(:)
    ! --------------------------------------------------------------------------
    write(dbgUnit(4),*) 'setup indices for the points of derived variable ', &
      &                 trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr )

    ! allcoate the index array for all inpits
    if (.not. allocated(fPtr%opData%input_pntIndex)) then
      allocate( fPtr%OpData%input_pntIndex(fun%nInputs) )
    end if

    ! allocate temporary inputIndex with size of nInputs and initialize
    ! growing array with length nPnts
    allocate(inputIndex_loc(fun%nInputs))

    ! store which is the last entry in the indexLvl to contiguous fill the index
    ! array string from this position
    ! all input variables get the same points, we just take the nVals entry
    ! from the first input variable
    nVals_prev = fPtr%OpData%input_pntIndex(1)%indexLvl(iLevel)%nVals

    ! Now fill in the index arrays for the inputs
    call tem_opVar_fill_inputIndex( fun        = fun,           &
      &                             varSys     = varSys,        &
      &                             point      = point,         &
      &                             offset_bit = offset_bit,    &
      &                             iLevel     = iLevel,        &
      &                             tree       = tree,          &
      &                             nPnts      = nPnts,         &
      &                             inputIndex = inputIndex_loc )

    ! fill the index array of the derived variable, it starts with the first
    ! entry in this call = nVals_prev and is continguous until nVals_prev+nVals
    allocate(idxPerPnt(fun%nInputs))
    idx = 0
    do iPnt = 1, nPnts
      do iDep = 1, fun%nInputs
        idxPerPnt(iDep) = inputIndex_loc(iDep)%val(iPnt)
      end do
      ! set index only when any of dependent variable has valid index
      if (any(idxPerPnt > 0)) then
        do iDep = 1, fun%nInputs
          call append( me = fPtr%opData%input_pntIndex(iDep)%indexLvl(iLevel), &
            &         val = inputIndex_loc(iDep)%val(iPnt)                     )
        end do
        ! set index to last position in input_pntIndex of dep var 1 of
        ! indexLvl of iLevel
        idx(iPnt) = fPtr%opData%input_pntIndex(1)%indexLvl(iLevel)%nVals
      end if
    end do

    do iDep = 1, fun%nInputs
      call truncate (fPtr%opData%input_pntIndex(iDep)%indexLvl(iLevel) )
    end do

  end subroutine mus_opVar_setupIndices
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine returns the velocity gradient from velocity in auxField
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine mus_opVar_gradU_forElement(fun, varsys, elempos, time, &
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
    real(kind=rk) :: gradU(3,3,1)
    integer :: velPos(3), elemOff, nDims
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    associate ( levelPointer => fPtr%solverData%geometry%levelPointer,   &
      &         gradData     => fPtr%solverData%scheme%gradData,         &
      &         auxField     => fPtr%solverData%scheme%auxField,         &
      &         derVarPos    => fPtr%solverData%scheme%derVarPos(1),     &
      &         stencil      => fPtr%solverData%scheme%layout%fStencil,  &
      &         Grad         => fPtr%solverData%scheme%Grad              )


      nDims = stencil%nDims
      velPos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)

      ! res is always AOS layout
      res = 0.0_rk
      do iElem = 1, nElems
        ! to acess Grad calculation pointers
        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )

        ! compute velocity gradient from velocity field stored in auxField
        gradU = 0.0_rk

        gradU(:nDims,:nDims,:1) = Grad%U_ptr(                        &
          &                   auxField    = auxField(iLevel)%val(:), &
          &                   gradData    = gradData(iLevel),        &
          &                   velPos      = velPos,                  &
          &                   nAuxScalars = varSys%nAuxScalars,      &
          &                   nDims       = nDims,                   &
          &                   nSolve      = 1,                       &
          &                   elemOffset  = statePos - 1             )

        ! store gradU in res array
        elemOff= (iElem-1)*9
        res( elemOff + 1 ) = gradU(1,1,1) !dudx
        res( elemOff + 2 ) = gradU(1,2,1) !dudy
        res( elemOff + 3 ) = gradU(1,3,1) !dudz
        res( elemOff + 4 ) = gradU(2,1,1) !dvdx
        res( elemOff + 5 ) = gradU(2,2,1) !dvdy
        res( elemOff + 6 ) = gradU(2,3,1) !dvdz
        res( elemOff + 7 ) = gradU(3,1,1) !dwdx
        res( elemOff + 8 ) = gradU(3,2,1) !dwdy
        res( elemOff + 9 ) = gradU(3,3,1) !dwdz
      end do !iElem

    end associate

  end subroutine mus_opVar_gradU_forElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine computes vorticity from curl of velocity in auxField.
  !! $$vorticity = \nabla \times u$$
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine mus_opVar_vorticity_forElement(fun, varsys, elempos, &
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
    integer :: statePos, iElem, iLevel
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk) :: gradU(3,3,1)
    real(kind=rk) :: vorticity(3)
    integer :: velPos(3), elemOff,nDims
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    associate ( levelPointer => fPtr%solverData%geometry%levelPointer,   &
      &         gradData     => fPtr%solverData%scheme%gradData,         &
      &         auxField     => fPtr%solverData%scheme%auxField,         &
      &         derVarPos    => fPtr%solverData%scheme%derVarPos(1),     &
      &         stencil      => fPtr%solverData%scheme%layout%fStencil,  &
      &         Grad         => fPtr%solverData%scheme%Grad              )

      nDims = stencil%nDims
      velPos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)

      ! res is always AOS layout
      res = 0.0_rk
      vorticity = 0.0_rk
      do iElem = 1, nElems

        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )

        ! compute velocity gradient from velocity field stored in auxField
        gradU(:nDims,:nDims,:) = Grad%U_ptr(                     &
          &               auxField    = auxField(iLevel)%val(:), &
          &               gradData    = gradData(iLevel),        &
          &               velPos      = velPos,                  &
          &               nAuxScalars = varSys%nAuxScalars,      &
          &               nDims       = nDims,                   &
          &               nSolve      = 1,                       &
          &               elemOffset  = statePos - 1             )

        ! compute vorticity
        if (nDims > 1) then
          vorticity(3) = gradU(2,1,1) - gradU(1,2,1) !dvdx - dudy
          if (nDims > 2) then
            vorticity(1) = gradU(3,2,1) - gradU(2,3,1) !dwdy - dvdz
            vorticity(2) = gradU(1,3,1) - gradU(3,1,1) !dxdz - dwdx
          end if
        end if

        ! store gradU in res array
        elemOff= (iElem-1)*3
        res( elemOff + 1 : elemOff + 3 ) = vorticity(:)
      end do !iElem

    end associate

  end subroutine mus_opVar_vorticity_forElement
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine computes Q-criterion from velocity in auxField.
  !! $$Q = 1/2 ( (tr(\nabla u))^2 - tr(\nabla u \cdot \nabla u) )
  !! = 1/2 |\Omega|^2 - |S|^2 $$,
  !! where $$\Omega$$ and S are asymmetric (vorticity tensor) and
  !! symmetric (rate of strain) part of velocity gradient.
  !! i.e $$\Omega = 1/2(du_i/dx_j - du_j/dx_i)$$ and
  !! $$S=1/2(du_i/dx_j + du_j/dx_i)$$.
  !!
  !! Ref: http://www.ctr.maths.lu.se/media/thesis/2013/FMA820.pdf
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  recursive subroutine mus_opVar_QCriterion_forElement(fun, varsys, elempos, &
    &                                                  time, tree, nElems,   &
    &                                                  nDofs, res            )
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
    integer :: statePos, iElem, iLevel, iDir
    type(mus_varSys_data_type), pointer :: fPtr
    real(kind=rk):: gradU(3,3,1), gradU_sqr(3,3)
    real(kind=rk) :: tr_gradU, tr_gradU_sqr
    integer :: velPos(3), nDims
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )

    associate ( levelPointer => fPtr%solverData%geometry%levelPointer,  &
      &         gradData     => fPtr%solverData%scheme%gradData,        &
      &         auxField     => fPtr%solverData%scheme%auxField,        &
      &         derVarPos    => fPtr%solverData%scheme%derVarPos(1),    &
      &         stencil      => fPtr%solverData%scheme%layout%fStencil, &
      &         Grad         => fPtr%solverData%scheme%Grad             )

      nDims = stencil%nDims
      velPos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)

      ! res is always AOS layout
      res = 0.0_rk
      do iElem = 1, nElems

        iLevel = tem_levelOf( tree%treeID( elemPos(iElem) ) )

        ! if state array is defined level wise then use levelPointer(pos)
        ! to access state array
        statePos = levelPointer( elemPos(iElem) )
        ! compute velocity gradient from velocity field stored in auxField
        gradU(:nDims,:nDims,:) = Grad%U_ptr(              &
          &               auxField    = auxField(iLevel)%val(:), &
          &               gradData    = gradData(iLevel),        &
          &               velPos      = velPos,                  &
          &               nAuxScalars = varSys%nAuxScalars,      &
          &               nDims       = nDims,                   &
          &               nSolve      = 1,                       &
          &               elemOffset  = statePos - 1             )


        ! square of velocity gradient. gradU . gradU
        gradU_sqr(:nDims, :nDims) = matmul( gradU(:nDims,:nDims,1), &
          &                                 gradU(:nDims,:nDims,1)  )
        ! trace of gradU_sqr
        tr_gradU_sqr = gradU_sqr(1,1)
        ! trace of gradU
        tr_gradU = gradU(1,1,1)
        do iDir=2,nDims
          tr_gradU_sqr = tr_gradU_sqr + gradU_sqr(iDir,iDir)
          tr_gradU = tr_gradU + gradU(iDir, iDir, 1)
        end do

        ! compute and store Q-criterion
        res( iElem ) = 0.5_rk * (tr_gradU*tr_gradU - tr_gradU_sqr)

      end do !iElem

    end associate

  end subroutine mus_opVar_QCriterion_forElement
  ! ************************************************************************** !

end module mus_operation_var_module
! **************************************************************************** !
