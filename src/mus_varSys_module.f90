! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016, 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!> author: Kannan Masilamani, Verena Krupp
!! This module contains solver data type with pointers to musubi
!! [[mus_scheme_type]], [[mus_param_type]], etc which are required for
!! variable operation method data.
!! Also contains all general routines for the variable system.
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
module mus_varSys_module

  use, intrinsic :: iso_c_binding,  only: c_ptr, c_f_pointer, c_loc
  use env_module,                   only: rk, long_k
  use tem_param_module,         only: qOffset_inChar, q000
  use tem_logging_module,           only: logUnit
  use tem_debug_module,             only: dbgUnit
  use tem_aux_module,               only: tem_abort
  use tem_float_module,             only: operator(.fne.)
  use tem_varSys_module,            only: tem_varSys_type,                    &
    &                                     tem_varSys_op_type,                 &
    &                                     tem_varSys_solverData_evalElem_type
  use tem_pointData_module,         only: tem_pointData_list_type, &
    &                                     tem_pointData_type, append
  use tem_grow_array_module,        only: append, truncate
  use treelmesh_module,             only: treelmesh_type
  use tem_geometry_module,          only: tem_CoordOfReal, &
    &                                     tem_PosofId, tem_BaryOfId
  use tem_time_module,              only: tem_time_type
  use tem_topology_module,          only: tem_IdOfCoord, tem_levelOf, &
    &                                     tem_FirstIdAtLevel, tem_CoordOfId
  use tem_construction_module,      only: tem_levelDesc_type
  use tem_operation_var_module,     only: tem_varSys_op_data_type
  use tem_spacetime_fun_module,     only: tem_st_fun_listelem_type
  use tem_element_module,           only: eT_fluid
  use tem_stencil_module,           only: tem_stencilHeader_type
  use tem_dyn_array_module,         only: PositionOfVal

  ! include musubi modules
  use mus_scheme_type_module,  only: mus_scheme_type
  use mus_physics_module,      only: mus_physics_type
  use mus_geom_module,         only: mus_geom_type
  use mus_connectivity_module, only: mus_intp_getSrcElemPosInTree
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type


  implicit none

  private

  public :: mus_varSys_solverData_type
  public :: mus_varSys_data_type
  public :: mus_get_new_solver_ptr
  public :: mus_init_varSys_solverData
  public :: mus_deriveVar_forPoint
  public :: mus_generic_varFromPDF_fromIndex
  public :: mus_set_stFun_getElement
  public :: mus_derive_fromPDF
  public :: mus_generic_fromPDF_forElement
  public :: mus_createSrcElemInTreeForGetPoint
  public :: mus_derVar_intpOnPoint

  !> Method data container for every variable
  type mus_varSys_data_type
    ! type with all solver relevant data for the variable
    type(mus_varSys_solverData_type), pointer :: solverData

    ! the point_datas need to be stored levelwise
    type(tem_pointData_list_type) :: pointData

    !> data array for operation or derived varibales
    !! consists the index arrys for points stored in the
    !! poingtData of input variable
    !! size is number of input variables
    type(tem_varSys_op_data_type) :: opData
  end type mus_varSys_data_type

  !> Contains pointer to musubi data types required for variable operation
  !! method_data
  type mus_varSys_solverData_type
    !> scheme data type
    type(mus_scheme_type), pointer :: scheme => NULL()

    !> contains basic SI units to convert from lattice to physical and
    !! vice versa
    type(mus_physics_type), pointer :: physics => NULL()

    !> Contains geometry information and definitions
    type(mus_geom_type), pointer :: geometry => NULL()

  end type mus_varSys_solverData_type

  abstract interface
    !> This interface describes the arguments to be used for routines that do
    !! the derivation of variables from the variable system.
    subroutine mus_derive_fromPDF( fun, varSys, stencil, iLevel, posInState, &
      &                            pdf, res, nVals )
      import :: tem_varSys_op_type,     &
        &       tem_varSys_type,        &
        &       tem_stencilHeader_type, &
        &       rk
      !> Description of the method to obtain the variables,
      class(tem_varSys_op_type), intent(in) :: fun
      !> The variable system to obtain the variable from.
      type(tem_varSys_type), intent(in) :: varSys
      !> fluid stencil defintion
      type(tem_stencilHeader_type), intent(in)  :: stencil
      !> current level
      integer, intent(in) :: iLevel
      !> Position of element in levelwise state array
      integer, intent(in) :: posInState(:)
      !> pdf
      real(kind=rk), intent(in) :: pdf(:)
      !> results
      real(kind=rk), intent(out) :: res(:)
      !> nVals to get
      integer, intent(in) :: nVals
    end subroutine mus_derive_fromPDF
  end interface

contains


  ! ************************************************************************* !
  !> Routine to get a pointer to a new instance of mus_varSys_solverData_type
  !! to be used as method data for a variable in the variable system.
  !!
  !! A new instance is allocated and a c_ptr to this type is returned.
  !! Be aware that local pointer are not automatically deallocated when leaving
  !! the routine
  function mus_get_new_solver_ptr( solver ) result(resPtr)
    ! --------------------------------------------------------------------- !
    !> The prototype is used to initialize the new instance.
    type(mus_varSys_solverData_type), intent(in), optional, target :: solver
    !> Pointer to the newly created instance.
    type(c_ptr) :: resPtr
    ! --------------------------------------------------------------------- !
    !> Local variable to allocate a new instance.
    type(mus_varSys_data_type), pointer :: res
    ! --------------------------------------------------------------------- !

    allocate(res)
    if (present(solver)) res%solverData => solver

    !>TODO init the point arrays

    resPtr = c_loc(res)

  end function mus_get_new_solver_ptr
  ! ************************************************************************* !


  ! ************************************************************************* !
  ! This routine sets varSys_solverData pointers
  subroutine mus_init_varSys_solverData( me, scheme, physics, geometry )
    ! --------------------------------------------------------------------- !
    type(mus_varSys_solverData_type),  intent(out) :: me
    type(mus_scheme_type),  target, intent(in) :: scheme
    type(mus_physics_type), target, intent(in) :: physics
    type(mus_geom_type),    target, intent(in) :: geometry
    ! --------------------------------------------------------------------- !
    me%scheme => scheme
    me%physics => physics
    me%geometry => geometry
  end subroutine mus_init_varSys_solverData
  ! ************************************************************************* !


  ! ************************************************************************ !
  !> Routine to store musubi varSys Data in stFun variable
  subroutine mus_set_stFun_getElement(solData_evalElem, fun)
    ! -------------------------------------------------------------------- !
    !> Description on how to set the element retrieval function for stfuns.
    class(tem_varSys_solverData_evalElem_type), intent(in) :: solData_evalElem

    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    type(tem_varSys_op_type), intent(inout) :: fun
    ! -------------------------------------------------------------------- !
    type(tem_st_fun_listElem_type), pointer :: fptr
    type(mus_varSys_solverData_type), pointer :: fSDptr
    ! -------------------------------------------------------------------- !

    write(logunit(10),*) "Setting different solver_bundle and" &
      & // " get_element routine for variable at position ",   &
      & fun%myPos
    call C_F_Pointer(fun%method_data, fptr)
    call c_f_pointer(solData_evalElem%solver_bundle, fSDptr)
    fptr%solver_bundle = mus_get_new_solver_ptr( fSDptr )

  end subroutine mus_set_stFun_getElement
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Derive variable for a given set of points using linear interpolation.
  !! This is a generic routine for any variable.
  !! Limitation: If neighbor is halo element then its not considered for
  !! interpolation, only the fluid (non-ghost) elements in the local process
  !! are used for interpolation.
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_point]].
  !!
  recursive subroutine mus_deriveVar_forPoint(fun, varsys, point, time, tree, &
    &                                         nPnts, res                      )
    !--------------------------------------------------------------------- !
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
    ! -------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! -------------------------------------------------------------------- !
    call C_F_POINTER( fun%method_Data, fPtr )
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
  end subroutine mus_deriveVar_forPoint
  ! ************************************************************************* !

  ! ************************************************************************* !
  subroutine mus_derVar_intpOnPoint(varPos, varSys, tree, time, point, nPnts, &
    &                                stencil, levelPointer, levelDesc, res     )
    ! -------------------------------------------------------------------- !
    !> Position of variable in varSys
    integer, intent(in) :: varPos

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil

    !> Pointer from treeIDlist entry to level-wise fluid part of total list
    integer, intent(in)         :: levelPointer(:)

    !> level description of all levels
    type(tem_levelDesc_type), intent(in) :: levelDesc(tree%global%minLevel:)

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    integer :: iPnt, iSrc, nSrcElems, nComponents
    integer :: srcElemPos(stencil%QQ)
    real(kind=rk) :: weights(stencil%QQ)
    real(kind=rk), allocatable :: srcRes(:), pntVal(:)
    ! -------------------------------------------------------------------- !
    nComponents = varSys%method%val(varPos)%nComponents
    allocate(srcRes(stencil%QQ*nComponents))
    allocate(pntVal(nComponents))

    res = 0.0_rk
    do iPnt = 1, nPnts
      srcRes = 0.0_rk
      ! get position of source element position in global tree for
      ! interpolation.
      ! Also calculate weights for interpolation using distance between the
      ! point and source element barycenter
      call mus_intp_getSrcElemPosInTree( &
        & srcElemPos   = srcElemPos,     &
        & weights      = weights,        &
        & nSrcElems    = nSrcElems,      &
        & point        = point(iPnt,:),  &
        & stencil      = stencil,        &
        & tree         = tree,           &
        & levelPointer = levelPointer,   &
        & levelDesc    = levelDesc(:)    )

      ! Skip this element if no source element found
      if (nSrcElems == 0) then
        call tem_abort("In mus_derVar_intpOnPoint: No source element found "&
          &             //"for interpolation")
      end if

      ! get source element values
      call varSys%method%val(varPos)%get_element( &
         & varSys  = varSys,                      &
         & elemPos = srcElemPos(1:nSrcElems),     &
         & time    = time,                        &
         & tree    = tree,                        &
         & nElems  = nSrcElems,                   &
         & nDofs   = 1,                           &
         & res     = srcRes                       )

      ! Linear interpolation res = sum(weight_i*phi_i)
      pntVal = 0.0_rk
      do iSrc = 1, nSrcElems
        pntVal(:) = pntVal(:) + weights(iSrc)                         &
          & * srcRes((iSrc-1)*nComponents+1 : iSrc*nComponents)
      end do
!      write(*,*) 'Linear intp pntVal ', pntVal
!      pntVal = mus_interpolate_quadratic2d_leastSq( &
!        &          srcMom       = srcRes,           &
!        &          targetCoord  = point(iPnt,1:2),  &
!        &          nSourceElems = nSrcElems         )
!      write(*,*) 'quadratic intp pntVal ', pntVal
      res( (iPnt-1)*nComponents+1 : iPnt*nComponents ) = pntVal

    end do !iPnt
  end subroutine mus_derVar_intpOnPoint
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> Routine to get the actual value for a given array of indices for
  !! musubi derive variables
  !! The indices belong to the grwarray of points storing levelwise in
  !! Pointdata%pntLvl(iLevel).
  !! Hence this routines takes the indeices as input, can refer to the pointData
  !! and evaluate the variable and returns the values
  subroutine mus_generic_varFromPDF_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                       idxLen, nVals, fnCalcPtr, res )
    ! -------------------------------------------------------------------------- !
    !> Description of the method to obtain the variables,
    class(tem_varSys_op_type), intent(in)   :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in)       :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)         :: time

    !> Level on which values are requested
    integer, intent(in)                     :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: n
    integer, intent(in)                     :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: dependes on number of first index for contiguous array,
    !! but the sum of all idxLen is equal to nVals
    integer, optional, intent(in)           :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in)                     :: nVals

    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer  :: fnCalcPtr

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out)              :: res(:)
    ! -------------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer     :: fPtr, state_fPtr
    type(mus_scheme_type), pointer          :: scheme
    real(kind=rk), allocatable              :: pdf(:)
    integer, allocatable                    :: posInState(:)
    integer                                 :: varPos, pdfElemSize
    integer                                 :: elemPos, iVal
    ! -------------------------------------------------------------------------- !
    write(dbgUnit(4),*) 'get the derived values of indices for variable ',  &
      &                  trim(varSys%varname%val(fun%myPos))
    ! distinguish if we have an array of index or we have contingous memory
    ! access where index are always first entries!
    if (present(idxLen)) then
      call tem_abort('Error: idxLen is not supported in get_valOfIndex for ' &
        &          //'state variable')
    end if

    !convert pointer from C to Fotran
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    ! pdf entries per Element
    pdfElemSize = scheme%layout%fStencil%QQ*scheme%nFields

    allocate( pdf( pdfElemSize*nVals ) )
    allocate( posInState(nVals) )

    varPos = fun%input_varPos(1)

    ! get pdf values for IDX
    call varSys%method%val( varPos )%get_ValOfIndex(  &
      &     varSys = varSys,                          &
      &     time   = time,                            &
      &     iLevel = iLevel,                          &
      &     idx    = fPtr%opData%input_pntIndex(1)    &
      &              %indexLvl(iLevel)%val( idx(:) ), &
      &     nVals  = nVals,                           &
      &     res    = pdf                              )

    ! the pointData is stored only in state variable function pointer
    call C_F_POINTER( varSys%method%val(varPos)%method_Data, state_fPtr )
    do iVal = 1, nVals
      ! Position in original tree
      elemPos = state_fPtr%pointData%pntLvl(iLevel)%elemPos%val(idx(iVal))
      ! position in levelwise state list
      posInState(iVal) = fPtr%solverData%geometry%levelPointer(elemPos)
    end do

    ! Call the procedure that does the calculation
    call fnCalcPtr(fun        = fun,                    &
      &            varSys     = varSys,                 &
      &            iLevel     = iLevel,                 &
      &            posInState = posInState,             &
      &            stencil    = scheme%layout%fstencil, &
      &            pdf        = pdf,                    &
      &            nVals      = nVals,                  &
      &            res        = res                     )

    deallocate( pdf )
    deallocate( posInState )

  end subroutine mus_generic_varFromPDF_fromIndex
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> This routine prepares the data for variable derivation or operators. It
  !! gathers all input variables from the variable system, calls the function
  !! with the actual calculation.
  recursive subroutine mus_generic_fromPDF_forElement(fun, varSys, elempos,    &
    &  tree, time, nVals, fnCalcPtr, nDofs, res )
    ! --------------------------------------------------------------------------
    ! -------------------------------------------------------------------------- !
    !> description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varsys_op_type), intent(in)   :: fun

    !> the variable system to obtain the variable from.
    type(tem_varsys_type), intent(in)       :: varsys

    !> position of the treeid of the element to get the variable for in the
    !! global treeid list.
    integer, intent(in)                     :: elempos(:)

    !> point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)         :: time

    !> Number of values to obtain for this variable
    integer, intent(in)                     :: nVals

    !> global treelm mesh info
    type(treelmesh_type), intent(in)        :: tree

    !> Function pointer to perform specific operation.
    procedure(mus_derive_fromPDF), pointer  :: fnCalcPtr

    !> number of degrees of freedom within an element.
    integer, intent(in)                     :: ndofs

    !> resulting values for the requested variable.
    !!
    !! linearized array dimension:
    !! (n requested entries) x (ncomponents of this variable)
    !! x (ndegrees of freedom)
    !! access: (ielem-1)*fun%ncomponents*ndofs +
    !!         (idof-1)*fun%ncomponents + icomp
    real(kind=rk), intent(out)              :: res(:)
    ! --------------------------------------------------------------------------
    !> contain all pdf's for element
    real(kind=rk), allocatable              :: pdf(:)
    real(kind=rk), allocatable              :: loc_res(:)
    type(mus_varSys_data_type), pointer     :: fPtr
    type(mus_scheme_type), pointer          :: scheme
    integer                                 :: iElem, pdfpos
    integer                                 :: iLevel, posInState
    ! --------------------------------------------------------------------------
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    ! position of PDF in glob System
    pdfpos = fun%input_varPos(1)
    allocate( pdf(scheme%layout%fStencil%QQ*scheme%nFields) )
    allocate( loc_res( fun%nComponents ) )
    ! loop over all elements ( could also be done in chunks or in total )
    do iElem = 1, nVals
     ! get source element value
     call varSys%method%val(pdfpos)%get_element( &
       & varSys  = varSys,                       &
       & elemPos = (/ elempos(iElem) /),         &
       & time    = time,                         &
       & tree    = tree,                         &
       & nElems  = 1,                            &
       & nDofs   = nDofs,                        &
       & res     = pdf                           )

     ! get iLevel for element
     iLevel = tem_levelOf( tree%treeID( elemPos(iElem ) ) )
     ! get posInState
     posInState = fPtr%solverData%geometry%levelPointer( elemPos(iElem) )

     ! Call the procedure that does the calculation
     call fnCalcPtr(fun        = fun,                    &
       &            varSys     = varSys,                 &
       &            iLevel     = iLevel,                 &
       &            posInState = (/ posInState /),       &
       &            stencil    = scheme%layout%fstencil, &
       &            pdf        = pdf,                    &
       &            nVals      = 1,                      &
       &            res        = loc_res                 )
      res((iElem-1)*fun%nComponents + 1: iElem*fun%nComponents ) = loc_res
    end do !iElem
    deallocate(pdf)
    deallocate(loc_res)
  end subroutine mus_generic_fromPDF_forElement
! **************************************************************************** !


  ! ************************************************************************** !
  !> This routine creates srcElemInTree in pointData. It is called all in
  !! getPoint routine when first time the get point routine in called.
  subroutine mus_createSrcElemInTreeForGetPoint(pntDataMaptoTree, posInPntData,&
    &                                           point, nPnts, stencil, tree,   &
    &                                           levelPointer, levelDesc)
    ! ------------------------------------------------------------------------ !
    !> Contains position of source elements in Tree and weights for
    !! interpolation
    type(tem_pointData_type), intent(out) :: pntDataMapToTree
    !> Contains position of a point in pntDataMapToTree
    integer, intent(out) :: posInPntData(:)
    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)
    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts
    !> stencil definition
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree
    !> Pointer from treeIDlist entry to level-wise fluid part of total list
    integer, intent(in)         :: levelPointer(:)
    !> level description of all levels
    type(tem_levelDesc_type), intent(in) :: levelDesc(tree%global%minLevel:)
    ! ------------------------------------------------------------------------ !
    integer :: iPnt, QQ, nSrcElems
    integer, allocatable :: srcElemPos(:)
    real(kind=rk), allocatable :: weights(:)
    logical :: wasAdded
    ! ------------------------------------------------------------------------ !
    QQ = stencil%QQ
    allocate(srcElemPos(QQ))
    allocate(weights(QQ))
    do iPnt = 1, nPnts
      weights = 0.0_rk

      ! append point, offset_bit and elemPos to pointData type
      call append(me             = pntDataMapToTree,     &
        &         point          = point(iPnt,:),        &
        &         storePnt       = .false.,              &
        &         offset_bit     = qOffset_inChar(q000), &
        &         storeOffsetBit = .false.,              &
        &         elemPos        = 0,                    &
        &         tree           = tree,                 &
        &         pos            = posInPntData(iPnt),   &
        &         wasAdded       = wasAdded              )

      if (wasAdded) then
        ! get position of source element position in global tree for
        ! interpolation.
        ! Also calculate weights for interpolation using distance between the
        ! point and source element barycenter
        call mus_intp_getSrcElemPosInTree( &
          & srcElemPos   = srcElemPos,     &
          & weights      = weights,        &
          & nSrcElems    = nSrcElems,      &
          & point        = point(iPnt,:),  &
          & stencil      = stencil,        &
          & tree         = tree,           &
          & levelPointer = levelPointer,   &
          & levelDesc    = levelDesc(:)    )

        ! First position of source elements
        call append(me  = pntDataMapToTree%srcElem%first,            &
          &         val = pntDataMapToTree%srcElem%elemPos%nVals + 1 )

        ! Append all src elemPos
        call append(me  = pntDataMapToTree%srcElem%elemPos,   &
          &         val = srcElemPos(1:nSrcElems)  )

        ! last position of source elements
        call append(me  = pntDataMapToTree%srcElem%last,         &
          &         val = pntDataMapToTree%srcElem%elemPos%nVals )

        ! weights for srcElements
        call append(me  = pntDataMapToTree%srcElem%weight, &
          &         val = weights(1:nSrcElems)  )
        end if
      end do

      call truncate(pntDataMapToTree%srcElem%first)
      call truncate(pntDataMapToTree%srcElem%last)
      call truncate(pntDataMapToTree%srcElem%elemPos)
      call truncate(pntDataMapToTree%srcElem%weight)

  end subroutine mus_createSrcElemInTreeForGetPoint
  ! ************************************************************************** !

end module mus_varSys_module
! **************************************************************************** !

