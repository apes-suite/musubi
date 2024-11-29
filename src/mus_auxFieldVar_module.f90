! Copyright (c) 2019-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2021-2023 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!! This module contains routine to retrieve auxiliary field variables for
!! getElement, getPoint, setupIndices and getValOfIndex.
!! Auxilary field variables are:
!!    * density and velocity for fluid
!!    * species desity and velocity for multispecies
!!    * potential for poisson
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
module mus_auxFieldVar_module
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,          only: rk, rk_mpi
  use tem_aux_module,      only: tem_abort
  use tem_param_module,    only: rho0, rho0Inv, cs2, cs2inv, cs4inv
  use tem_varSys_module,   only: tem_varSys_type, tem_varSys_op_type
  use tem_time_module,     only: tem_time_type
  use treelmesh_module,    only: treelmesh_type
  use tem_geometry_module, only: tem_CoordOfReal, &
    &                            tem_PosofId, tem_BaryOfId
  use tem_topology_module, only: tem_IdOfCoord
  use tem_topology_module, only: tem_levelOf
  use tem_stencil_module,  only: tem_stencilHeader_type
  use tem_comm_module,     only: tem_commpattern_type, tem_communication_type
  use tem_construction_module, only: tem_levelDesc_type
  use tem_logging_module,  only: logUnit
  use tem_math_module,     only: invert_matrix
  use tem_compileconf_module,        only: vlen
  use tem_debug_module,    only:dbgunit

  use mus_varSys_module,         only: mus_varSys_data_type
  use mus_scheme_type_module,    only: mus_scheme_type
  use mus_scheme_header_module,  only: mus_scheme_header_type
  use mus_derVarPos_module,      only: mus_derVarPos_type
  use mus_physics_module,        only: mus_convertFac_type
  use mus_auxField_module,       only: mus_proc_calcAuxField
  use mus_source_type_module,    only: mus_source_op_type
  use mus_eNRTL_module,          only: mus_calc_thermFactor
  use mus_gradData_module,       only: mus_gradData_type
  use mus_scheme_type_module,    only: mus_scheme_type
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_connectivity_module,   only: mus_intp_getSrcElemPosInTree
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none
  private

  public :: mus_assign_calcAuxField_ptr
  public :: mus_access_auxFieldVar_forElement
  public :: mus_auxFieldVar_forPoint
  public :: mus_auxFieldVar_fromIndex
  public :: mus_addForceToAuxField_fluid
  public :: mus_addForceToAuxField_fluidIncomp
  public :: mus_addForceToAuxField_MSL
  public :: mus_addForceToAuxField_MSL_WTDF
  public :: mus_addElectricToAuxField_MSL
  public :: mus_addElectricToAuxField_MSL_WTDF
  public :: mus_addSrcToAuxField_poisson
  public :: mus_addSponFldToAuxField_fluid
  public :: mus_addDynSponFldToAuxField_fluid
  public :: mus_addTurbChanForceToAuxField_fluid
  public :: mus_addHRRCorrToAuxField_fluid_2D
  public :: mus_addHRRCorrToAuxField_fluid_3D

contains


  ! ************************************************************************* !
  !> This routine assign function pointer to compute auxField var
  subroutine mus_assign_calcAuxField_ptr(schemeHeader, calcAuxField)
    ! --------------------------------------------------------------------- !
    !> scheme defnition
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    !> function pointer to assign
    procedure(mus_proc_calcAuxField), pointer, intent(out) :: calcAuxField
    ! --------------------------------------------------------------------- !
    ! --------------------------------------------------------------------- !
    calcAuxField => mus_calcAuxField_dummy
    select case (trim(schemeHeader%kind))
    case ('fluid', 'fluid_incompressible', 'isotherm_acEq')
      select case (trim(schemeHeader%layout))
      case ('d2q9')
        calcAuxField => mus_calcAuxField_fluid_d2q9
      case ('d3q19')
        calcAuxField => mus_calcAuxField_fluid_d3q19
      case ('d3q27')
        calcAuxField => mus_calcAuxField_fluid_d3q27
      case default
        calcAuxField => mus_calcAuxField_fluid
      end select

    case ('passive_scalar','poisson', 'poisson_boltzmann_linear', &
      &   'poisson_boltzmann_nonlinear')
      calcAuxField => mus_calcAuxField_zerothMoment
    case ('nernst_planck')
      calcAuxField => mus_calcAuxField_nernst_planck
    case ('multispecies_gas', 'multispecies_liquid')
      calcAuxField => mus_calcAuxField_MS
    case default
      calcAuxField => mus_calcAuxField_dummy
    end select
  end subroutine mus_assign_calcAuxField_ptr
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> Return the solver aux variable for a given set of elements
  !!
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_element]].
  !!
  recursive subroutine mus_access_auxFieldVar_forElement(fun, varsys, elempos, &
    &                                                    time, tree, nElems,   &
    &                                                    nDofs, res            )
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
    integer :: statePos, iElem, iComp, iLevel
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
      do iComp = 1, fun%nComponents
        res( (iElem-1)*fun%nComponents+iComp ) =                         &
          & scheme%auxField( iLevel )%val(                               &
          ! position of this aux variable in the aux array
          & (statePos-1)*varSys%nAuxScalars + fun%auxField_varPos(iComp) )
      end do !iComp
    end do !iElem

  end subroutine mus_access_auxFieldVar_forElement
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> Auxilary field variable for a given set of points using linear
  !! interpolation. Unlike mus_deriveVar_forPoint which does not consider
  !! ghost and halo elements, this routine considers them because
  !! auxField vars on ghost elements are interpolated and halo elements are
  !! exchanged.
  !! The interface has to comply to the abstract interface
  !! [[tem_varSys_module:tem_varSys_proc_point]].
  !!
  recursive subroutine mus_auxFieldVar_forPoint(fun, varsys, point, time, &
    &                                           tree, nPnts, res          )
    ! -------------------------------------------------------------------- !
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
    integer :: elemPos, statePos, iPnt, iLevel
    real(kind=rk), allocatable :: srcRes(:), pntVal(:), weights(:)
    integer :: iSrc, iComp, nSrcElems
    integer, allocatable :: srcElemPos(:)
    ! -------------------------------------------------------------------- !
!write(dbgUnit(1),*) 'Derive for point :'//trim(varSys%varname%val(fun%myPos))
    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    allocate(srcRes(scheme%layout%fStencil%QQ*fun%nComponents))
    allocate(pntVal(fun%nComponents))
    allocate(weights(scheme%layout%fStencil%QQ))
    allocate(srcElemPos(scheme%layout%fStencil%QQ))

    res = 0.0_rk
    do iPnt = 1, nPnts
      srcRes = 0.0_rk

      ! get position of source element position in global tree for
      ! interpolation.
      ! Also calculate weights for interpolation using distance between the
      ! point and source element barycenter
      call mus_intp_getSrcElemPosInTree(          &
        & srcElemPos   = srcElemPos,              &
        & weights      = weights,                 &
        & nSrcElems    = nSrcElems,               &
        & point        = point(iPnt,:),           &
        & stencil      = scheme%layout%fStencil,  &
        & tree         = tree,                    &
        & levelPointer = fPtr%solverData%geometry &
        &                    %levelPointer,       &
        & levelDesc    = scheme%levelDesc         )

      ! get source element values
      do iSrc = 1, nSrcElems
        ! position in global tree
        elemPos = srcElemPos(iSrc)
        ! position of element in levelDesc total list
        statePos = fPtr%solverData%geometry%levelPointer( elemPos )
        iLevel = tem_levelOf( tree%treeID( elemPos ) )
        do iComp = 1, fun%nComponents
          srcRes( (iSrc-1)*fun%nComponents + iComp )                       &
            & = scheme%auxField( iLevel )%val(                             &
            & (statePos-1)*varSys%nAuxScalars + fun%auxField_varPos(iComp) )
        end do !iComp
      end do !iSrc

      ! Linear interpolation res = sum(weight_i*phi_i)
      pntVal = 0.0_rk
      do iSrc = 1, nSrcElems
        pntVal(:) = pntVal(:) + weights(iSrc)                         &
          & * srcRes((iSrc-1)*fun%nComponents+1 : iSrc*fun%nComponents)
      end do
      res( (iPnt-1)*fun%nComponents+1 : iPnt*fun%nComponents ) = pntVal

    end do !iPnt
  end subroutine mus_auxFieldVar_forPoint
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> Routine to get the actual value for a given array of indices.
  !! The indices belong to the grwarray of points storing levelwise in
  !! Pointdata%pntLvl(iLevel).
  !! Hence this routines takes the indeices as input, can refer to the pointData
  !! and evaluate the variable and returns the values
  subroutine mus_auxFieldVar_fromIndex( fun, varSys, time, iLevel, idx, &
    &                                   idxLen,  nVals,  res            )
    ! -------------------------------------------------------------------------- !
    !> Description of the method to obtain the variables,
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in) :: time

    !> Level on which values are requested
    integer, intent(in) :: iLevel

    !> Index of points in the growing array and variable val array to
    !! return.
    !! Size: n
    integer, intent(in) :: idx(:)

    !> With idx as start index in contiguous memory,
    !! idxLength defines length of each contiguous memory
    !! Size: dependes on number of first index for contiguous array,
    !! but the sum of all idxLen is equal to nVals
    integer, optional, intent(in) :: idxLen(:)

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nVals

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------------- !
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    integer :: iComp, iSrc, iVal, first, last, loc_level
    integer :: iSrcElems, statePos, nSrcElems
    real(kind=rk), allocatable :: srcRes(:), pntVal(:)
    real(kind=rk) :: weight
    ! -------------------------------------------------------------------------- !

!    write(dbgUnit(4),*) 'get the values of indices for variable ',  &
!      &                  trim(varSys%varname%val(fun%myPos))

    call C_F_POINTER( fun%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    allocate(srcRes(scheme%layout%fStencil%QQ*fun%nComponents))
    allocate(pntVal(fun%nComponents))

    res = 0.0_rk

    ! distinguish if we have an array of index or we have contingous memory
    ! access where index are always first entries!
    if (present(idxLen)) then
      call tem_abort('Error: idxLen is not supported in get_valOfIndex for ' &
        &          //'state variable')
    end if

    ! Get the state value at the specific point
    do iVal = 1, nVals
      if (idx(iVal)>0) then

        ! elemPos in tree
        ! elemPos = fPtr%pointData%pntLvl(iLevel)%elemPos%val(idx(iVal))

        ! Element position in level-wise list
        first = fPtr%pointData%pntLvl(iLevel)%srcElem%first%val(idx(iVal))
        last = fPtr%pointData%pntLvl(iLevel)%srcElem%last%val(idx(iVal))

        ! get pdf's of source elements
        srcRes = 0.0_rk
        iSrcElems = 0
        do iSrc = first, last
          iSrcElems = iSrcElems + 1

          statePos = fPtr%pointData%pntLvl(iLevel)%srcElem%elemPos%val(iSrc)
          loc_level = fPtr%pointData%pntLvl(iLevel)%pntLevel%val(idx(iVal))

          do iComp = 1, fun%nComponents
            srcRes( (iSrcElems-1)*fun%nComponents + iComp )                  &
              & = scheme%auxField(loc_level)%val(                            &
              & (statePos-1)*varSys%nAuxScalars + fun%auxField_varPos(iComp) )
          end do !iComp
        end do !iSrc

        ! Linear interpolation res = sum(weight_i*phi_i)
        pntVal = 0.0_rk
        nSrcElems = 0
        do iSrc = first, last
          weight = fPtr%pointData%pntLvl(iLevel)%srcElem%weight%val(iSrc)
          nSrcElems = nSrcElems + 1
          pntVal(:) = pntVal(:) + weight &
            & * srcRes( (nSrcElems-1)*fun%nComponents+1 &
            &         : nSrcElems*fun%nComponents )
        end do

        ! get the state value for each component of this
        res( (iVal-1)*fun%nComponents+1: iVal*fun%nComponents ) = pntVal
      end if !idx>0
    end do !iVal

  end subroutine mus_auxFieldVar_fromIndex
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> This routine compute auxFields density and velocity for compressible model
  !! for fluid and nGhostFromCoarser elements
  subroutine mus_calcAuxField_fluid(auxField, state, neigh, nSize, nSolve, &
    & iLevel, stencil, varSys, derVarPos, quantities)
    ! -------------------------------------------------------------------------- !
    !> output auxField array
    real(kind=rk), intent(inout) :: auxField(:)
    !> input state array
    real(kind=rk), intent(in) :: state(:)
    !> connectivity array
    integer, intent(in) :: neigh(:)
    !> number of elements in the state array
    integer, intent(in) :: nSize
    !> number of elements excluding halos
    integer, intent(in) :: nSolve
    !> current level
    integer, intent(in) :: iLevel
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! -------------------------------------------------------------------------- !
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: rho(vlen), pdf( stencil%QQ,vlen ), vel(3,vlen)
    integer :: iElem, iDir, QQ, nScalars
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! -------------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    QQ = stencil%QQ
    nScalars = varSys%nScalars

    ! vectorization
    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    !$omp parallel do                                                                            &
    !$omp   & private( rho, pdf, elemPos, elemOff, nChunkElems, iChunks, low_bound, iElem, vel ) &
    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      !NEC$ ivdep
      do iElem = 1, nChunkElems

        elemPos = low_bound + iElem

        !NEC$ shortloop
        do iDir = 1, QQ
          pdf(iDir,iElem) = state( neigh((idir-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        end do

        ! element offset
        elemoff = (elemPos-1)*varSys%nAuxScalars

        ! density
        rho(iElem) = sum( pdf )
        ! store density
        auxField(elemOff+dens_pos) = rho(iElem)

      end do

      ! compute velocity
      vel = quantities%vel_from_pdf_vectorized_ptr(pdf = pdf, dens = rho, &
        &                  cxDirRK = stencil%cxDirRK, nSolve = nChunkElems)

      !NEC$ ivdep
      do iElem = 1, nChunkElems
        ! element offset
        elemPos = low_bound + iElem
        elemoff = (elemPos-1)*varSys%nAuxScalars
        ! store velocity
        auxField(elemOff+vel_pos(1)) = vel(1,iElem)
        auxField(elemOff+vel_pos(2)) = vel(2,iElem)
        auxField(elemOff+vel_pos(3)) = vel(3,iElem)

      end do

    end do
    !$omp end parallel do

  end subroutine mus_calcAuxField_fluid
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> This routine compute auxFields density and velocity for compressible
  !! d2q9 model for fluid and nGhostFromCoarser elements

  subroutine mus_calcAuxField_fluid_d2q9(auxField, state, neigh, nSize, nSolve, &
    & iLevel, stencil, varSys, derVarPos, quantities)
    ! -------------------------------------------------------------------------- !
    !> output auxField array
    real(kind=rk), intent(inout) :: auxField(:)
    !> input state array
    real(kind=rk), intent(in) :: state(:)
    !> connectivity array
    integer, intent(in) :: neigh(:)
    !> number of elements in the state array
    integer, intent(in) :: nSize
    !> number of elements excluding halos
    integer, intent(in) :: nSolve
    !> current level
    integer, intent(in) :: iLevel
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! -------------------------------------------------------------------------- !
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: rho(vlen), pdf(9, vlen), vel(3, vlen)
    integer :: iElem, QQ, nScalars
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! -------------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    QQ = 9
    nScalars = varSys%nScalars

    ! vectorization
    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    !$omp parallel do                                                                            &
    !$omp   & private( rho, pdf, elemPos, elemOff, nChunkElems, iChunks, low_bound, iElem, vel ) &
    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      !NEC$ ivdep
      do iElem = 1, nChunkElems

        elemPos = low_bound + iElem

        pdf(1, iElem) = state( neigh((1-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(2, iElem) = state( neigh((2-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(3, iElem) = state( neigh((3-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(4, iElem) = state( neigh((4-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(5, iElem) = state( neigh((5-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(6, iElem) = state( neigh((6-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(7, iElem) = state( neigh((7-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(8, iElem) = state( neigh((8-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(9, iElem) = state( neigh((9-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)

        ! element offset
        elemoff = (elemPos-1)*varSys%nAuxScalars
        ! density
        rho(iElem) = sum(pdf(:, iElem))
        ! store density
        auxField(elemOff+dens_pos) = rho(iElem)

      end do

      ! compute velocity
      vel = quantities%vel_from_pdf_vectorized_ptr(pdf = pdf, dens = rho, nSolve = nChunkElems)

      !NEC$ ivdep
      do iElem = 1, nChunkElems
        ! element offset
        elemPos = low_bound + iElem
        elemoff = (elemPos-1)*varSys%nAuxScalars
        ! store velocity
        auxField(elemOff+vel_pos(1)) = vel(1, iElem)
        auxField(elemOff+vel_pos(2)) = vel(2, iElem)
        auxField(elemOff+vel_pos(3)) = 0.0_rk
      end do

    end do!iChunks
    !$omp end parallel do

  end subroutine mus_calcAuxField_fluid_d2q9
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> This routine compute auxFields density and velocity for compressible
  !! d3q19 model for fluid and nGhostFromCoarser elements
  subroutine mus_calcAuxField_fluid_d3q19(auxField, state, neigh, nSize, &
    &  nSolve, iLevel, stencil, varSys, derVarPos, quantities)
    ! ----------------------------------------------------------------------- !
    !> output auxField array
    real(kind=rk), intent(inout) :: auxField(:)
    !> input state array
    real(kind=rk), intent(in) :: state(:)
    !> connectivity array
    integer, intent(in) :: neigh(:)
    !> number of elements in the state array
    integer, intent(in) :: nSize
    !> number of elements excluding halos
    integer, intent(in) :: nSolve
    !> current level
    integer, intent(in) :: iLevel
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! ----------------------------------------------------------------------- !
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: rho(vlen), pdf(19, vlen), vel(3, vlen)
    integer :: iElem, QQ, nScalars
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! ----------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    QQ = 19
    nScalars = varSys%nScalars

    ! vectorization
    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    !$omp parallel do                                                                            &
    !$omp   & private( rho, pdf, elemPos, elemOff, nChunkElems, iChunks, low_bound, iElem, vel ) &
    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      !NEC$ ivdep
      do iElem = 1, nChunkElems

        elemPos = low_bound + iElem

        pdf( 1,iElem) = state( neigh(( 1-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf( 2,iElem) = state( neigh(( 2-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf( 3,iElem) = state( neigh(( 3-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf( 4,iElem) = state( neigh(( 4-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf( 5,iElem) = state( neigh(( 5-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf( 6,iElem) = state( neigh(( 6-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf( 7,iElem) = state( neigh(( 7-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf( 8,iElem) = state( neigh(( 8-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf( 9,iElem) = state( neigh(( 9-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(10,iElem) = state( neigh((10-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(11,iElem) = state( neigh((11-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(12,iElem) = state( neigh((12-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(13,iElem) = state( neigh((13-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(14,iElem) = state( neigh((14-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(15,iElem) = state( neigh((15-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(16,iElem) = state( neigh((16-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(17,iElem) = state( neigh((17-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(18,iElem) = state( neigh((18-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(19,iElem) = state( neigh((19-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)

        ! element offset
        elemoff = (elemPos-1)*varSys%nAuxScalars

        ! density
        rho(iElem) = sum(pdf(:,iElem))
        ! store density
        auxField(elemOff+dens_pos) = rho(iElem)

      end do

      ! compute velocity
      vel = quantities%vel_from_pdf_vectorized_ptr(pdf = pdf, dens = rho, nSolve = nChunkElems)

      !NEC$ ivdep
      do iElem = 1, nChunkElems
        ! element offset
        elemPos = low_bound + iElem
        elemoff = (elemPos-1)*varSys%nAuxScalars

        ! store velocity
        auxField(elemOff+vel_pos(1)) = vel(1, iElem)
        auxField(elemOff+vel_pos(2)) = vel(2, iElem)
        auxField(elemOff+vel_pos(3)) = vel(3, iElem)

      end do

    end do!iChunks
    !$omp end parallel do

  end subroutine mus_calcAuxField_fluid_d3q19
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> This routine compute auxFields density and velocity for compressible
  !! d3q27 model for fluid and nGhostFromCoarser elements
  subroutine mus_calcAuxField_fluid_d3q27(auxField, state, neigh, nSize, &
    &  nSolve, iLevel, stencil, varSys, derVarPos, quantities)
    ! ----------------------------------------------------------------------- !
    !> output auxField array
    real(kind=rk), intent(inout) :: auxField(:)
    !> input state array
    real(kind=rk), intent(in) :: state(:)
    !> connectivity array
    integer, intent(in) :: neigh(:)
    !> number of elements in the state array
    integer, intent(in) :: nSize
    !> number of elements excluding halos
    integer, intent(in) :: nSolve
    !> current level
    integer, intent(in) :: iLevel
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! ----------------------------------------------------------------------- !
    integer :: dens_pos, vel_pos(3), elemOff
    real(kind=rk) :: rho(vlen), pdf(27,vlen), vel(3,vlen)
    integer :: iElem, QQ, nScalars
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! ----------------------------------------------------------------------- !
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    QQ = 27
    nScalars = varSys%nScalars

    ! vectorization
    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    !$omp parallel do                                                                            &
    !$omp   & private( rho, pdf, elemPos, elemOff, nChunkElems, iChunks, low_bound, iElem, vel ) &
    do iChunks = 1, nChunks
      ! calculate the end  number of iElem loop
      nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
      low_bound = (iChunks - 1) * vlen

      !NEC$ ivdep
      do iElem = 1, nChunkElems

        elemPos = low_bound + iElem

        pdf( 1,iElem) = state( neigh((1-1)* nsize+ elempos)+(  1-1)* qq+ nscalars*0)
        pdf( 2,iElem) = state( neigh((2-1)* nsize+ elempos)+(  1-1)* qq+ nscalars*0)
        pdf( 3,iElem) = state( neigh((3-1)* nsize+ elempos)+(  1-1)* qq+ nscalars*0)
        pdf( 4,iElem) = state( neigh((4-1)* nsize+ elempos)+(  1-1)* qq+ nscalars*0)
        pdf( 5,iElem) = state( neigh((5-1)* nsize+ elempos)+(  1-1)* qq+ nscalars*0)
        pdf( 6,iElem) = state( neigh((6-1)* nsize+ elempos)+(  1-1)* qq+ nscalars*0)
        pdf( 7,iElem) = state( neigh((7-1)* nsize+ elempos)+(  1-1)* qq+ nscalars*0)
        pdf( 8,iElem) = state( neigh((8-1)* nsize+ elempos)+(  1-1)* qq+ nscalars*0)
        pdf( 9,iElem) = state( neigh((9-1)* nsize+ elempos)+(  1-1)* qq+ nscalars*0)
        pdf(10,iElem) = state( neigh((10-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(11,iElem) = state( neigh((11-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(12,iElem) = state( neigh((12-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(13,iElem) = state( neigh((13-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(14,iElem) = state( neigh((14-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(15,iElem) = state( neigh((15-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(16,iElem) = state( neigh((16-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(17,iElem) = state( neigh((17-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(18,iElem) = state( neigh((18-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(19,iElem) = state( neigh((19-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(20,iElem) = state( neigh((20-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(21,iElem) = state( neigh((21-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(22,iElem) = state( neigh((22-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(23,iElem) = state( neigh((23-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(24,iElem) = state( neigh((24-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(25,iElem) = state( neigh((25-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(26,iElem) = state( neigh((26-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)
        pdf(27,iElem) = state( neigh((27-1)* nsize+ elempos)+( 1-1)* qq+ nscalars*0)

        ! element offset
        elemoff = (elemPos-1)*varSys%nAuxScalars

        ! density
        rho(iElem) = sum(pdf(:,iElem))
        ! store density
        auxField(elemOff+dens_pos) = rho(iElem)

      end do

      ! compute velocity
      vel = quantities%vel_from_pdf_vectorized_ptr(pdf = pdf, dens = rho, nSolve = nChunkElems)

      !NEC$ ivdep
      do iElem = 1, nChunkElems
        ! element offset
        elemPos = low_bound + iElem
        elemoff = (elemPos-1)*varSys%nAuxScalars
        ! store velocity
        auxField(elemOff+vel_pos(1)) = vel(1, iElem)
        auxField(elemOff+vel_pos(2)) = vel(2, iElem)
        auxField(elemOff+vel_pos(3)) = vel(3, iElem)

      end do

    end do
    !$omp end parallel do

  end subroutine mus_calcAuxField_fluid_d3q27
 ! ************************************************************************* !


  ! ************************************************************************* !
  !> This routine compute zeroth moment from state and store in auxField.
  !! use this routine only for models which requires only zeroth-order
  !! moment as auxField
  subroutine mus_calcAuxField_zerothMoment(auxField, state, neigh, nSize, &
    & nSolve, iLevel, stencil, varSys, derVarPos, quantities)
    ! -------------------------------------------------------------------------- !
    !> output auxField array
    real(kind=rk), intent(inout) :: auxField(:)
    !> input state array
    real(kind=rk), intent(in) :: state(:)
    !> connectivity array
    integer, intent(in) :: neigh(:)
    !> number of elements excluding halos
    integer, intent(in) :: nSolve
    !> number of elements in the state array
    integer, intent(in) :: nSize
    !> current level
    integer, intent(in) :: iLevel
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! -------------------------------------------------------------------------- !
    real(kind=rk) :: pdf( stencil%QQ )
    integer :: iElem, iDir
    ! -------------------------------------------------------------------------- !
    ! safety check
    if (varSys%nAuxScalars /= 1) then
      call tem_abort('CalcAuxField_zeroth moment can be used only when' &
        &          //'nAuxScalars = 1')
    end if

    !NEC$ ivdep
    do iElem = 1, nSolve
      !NEC$ shortloop
      do iDir = 1, stencil%QQ
        pdf(iDir) = state(                                                    &
          &  neigh((idir-1)* nsize+ ielem)+( 1-1)* stencil%qq+ varsys%nscalars*0)
      end do

      ! store scalar zeroth moment
      auxField(iElem) = sum(pdf)
    end do
  end subroutine mus_calcAuxField_zerothMoment
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> This routine compute zeroth moment (mole density) from state for each
  !! species and store in auxField.
  subroutine mus_calcAuxField_nernst_planck(auxField, state, neigh, nSize, &
    & nSolve, iLevel, stencil, varSys, derVarPos, quantities)
    ! -------------------------------------------------------------------------- !
    !> output auxField array
    real(kind=rk), intent(inout) :: auxField(:)
    !> input state array
    real(kind=rk), intent(in) :: state(:)
    !> connectivity array
    integer, intent(in) :: neigh(:)
    !> number of elements excluding halos
    integer, intent(in) :: nSolve
    !> number of elements in the state array
    integer, intent(in) :: nSize
    !> current level
    integer, intent(in) :: iLevel
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! -------------------------------------------------------------------------- !
    real(kind=rk) :: pdf( stencil%QQ )
    integer :: iElem, iDir, iFld, nScalars, elemOff
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! -------------------------------------------------------------------------- !
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    nScalars = varSys%nScalars

    !NEC$ ivdep
    do iElem = 1, nSolve
      ! element offset
      elemOff = (iElem-1)*varSys%nAuxScalars
      ! store scalar mole density
      !NEC$ shortloop
      do iFld = 1, scheme%nFields
        do iDir = 1, stencil%QQ
          pdf(iDir) = state(                                                &
            &  neigh((idir-1)* nsize+ ielem)+( ifld-1)* stencil%qq+ nscalars*0)
        end do
        ! Since nernst_planck has only mole density in auxField array.
        ! So, instead of varSys->auxField_varPos, the array can be accessed
        ! using iFld itself
        auxField(elemOff + iFld) = sum(pdf)
      end do !iField
    end do !iElem
  end subroutine mus_calcAuxField_nernst_planck
  ! ************************************************************************* !


  ! ************************************************************************* !
  !> This routine compute auxField density and momentum for each species
  !! for multicomponent models. The momentum computed here is only momentum
  !! of transformed PDF. The momentum of original PDF is computed by solving
  !! linear equation system in compute kernel and the momentum in auxField is
  !! updated there.
  subroutine mus_calcAuxField_MS(auxField, state, neigh, nSize, &
    & nSolve, iLevel, stencil, varSys, derVarPos, quantities)
    ! ------------------------------------------------------------------------ !
    !> output auxField array
    real(kind=rk), intent(inout) :: auxField(:)
    !> input state array
    real(kind=rk), intent(in) :: state(:)
    !> connectivity array
    integer, intent(in) :: neigh(:)
    !> number of elements excluding halos
    integer, intent(in) :: nSolve
    !> number of elements in the state array
    integer, intent(in) :: nSize
    !> current level
    integer, intent(in) :: iLevel
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! ------------------------------------------------------------------------ !
    integer :: iElem, iFld, iDir, elemOff, nScalars
    integer :: QQ, nFields, dens_pos, mom_pos(3)
    real(kind=rk) :: pdf( stencil%QQ )
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! ------------------------------------------------------------------------ !
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme

    nFields = scheme%nFields
    QQ = stencil%QQ
    nScalars = varSys%nScalars

    !NEC$ ivdep
    do iElem = 1, nSolve
      ! element offset
      elemoff = (iElem-1)*varSys%nAuxScalars

      !NEC$ shortloop
      do iFld = 1, scheme%nFields
        do iDir = 1, QQ
          pdf(iDir) = state( neigh((idir-1)* nsize+ ielem)+( ifld-1)* qq+ nscalars*0)
        end do
        ! position of density and velocity of current field in auxField array
        dens_pos = varSys%method%val(derVarPos(iFld)%density)%auxField_varPos(1)
        mom_pos = varSys%method%val(derVarPos(iFld)%momentum)%auxField_varPos(:)

        ! store mass density of species
        auxField(elemOff+dens_pos) = sum(pdf)

        ! store momentum
        auxField(elemOff+mom_pos(1)) = sum(pdf * stencil%cxDirRK(1, :))
        auxField(elemOff+mom_pos(2)) = sum(pdf * stencil%cxDirRK(2, :))
        auxField(elemOff+mom_pos(3)) = sum(pdf * stencil%cxDirRK(3, :))
      end do !iField
    end do

  end subroutine mus_calcAuxField_MS
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> Dummy routine for calcAuxField
  subroutine mus_calcAuxField_dummy(auxField, state, neigh, nSize, nSolve, &
    & iLevel, stencil, varSys, derVarPos, quantities)
    ! -------------------------------------------------------------------------- !
    !> output auxField array
    real(kind=rk), intent(inout) :: auxField(:)
    !> input state array
    real(kind=rk), intent(in) :: state(:)
    !> connectivity array
    integer, intent(in) :: neigh(:)
    !> number of elements excluding halos
    integer, intent(in) :: nSolve
    !> number of elements in the state array
    integer, intent(in) :: nSize
    !> current level
    integer, intent(in) :: iLevel
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! -------------------------------------------------------------------------- !
    call tem_abort('Dummy routine for calcAuxField')
  end subroutine mus_calcAuxField_dummy
  ! ************************************************************************* !

  ! ************************************************************************** !
  !> This routine add body force to velocity in auxField for weakly-compressible
  !! model.
  subroutine mus_addForceToAuxField_fluid(fun, auxField, iLevel, time, varSys, &
    &                                   phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, vel_pos(3)
    real(kind=rk) :: forceTerm(3)
    integer :: iElem, nElems, posInTotal, elemOff
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3), inv_rho
    ! ------------------------------------------------------------------------ !
    ! position of density and velocity field in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems
    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / phyConvFac%body_force

!$omp parallel do schedule(static), private( posInTotal, forceTerm, inv_rho, elemOff )
    !NEC$ ivdep
    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)
      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! inverse of density
      inv_rho = 1.0_rk/auxField(elemOff+dens_pos)
      ! forceterm to add to velocity: F/(2*rho)
      forceTerm = forceField((iElem-1)*3+1 : iElem*3)*0.5_rk*inv_rho
      ! add force to velocity
      auxField(elemOff+vel_pos(1)) = auxField(elemOff+vel_pos(1)) + forceTerm(1)
      auxField(elemOff+vel_pos(2)) = auxField(elemOff+vel_pos(2)) + forceTerm(2)
      auxField(elemOff+vel_pos(3)) = auxField(elemOff+vel_pos(3)) + forceTerm(3)
    end do

  end subroutine mus_addForceToAuxField_fluid
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine add force to velocity in auxField for incompressible model
  subroutine mus_addForceToAuxField_fluidIncomp(fun, auxField, iLevel, time, &
    & varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: vel_pos(3), elemOff
    real(kind=rk) :: forceTerm(3)
    integer :: iElem, nElems, posInTotal
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    ! ------------------------------------------------------------------------ !
    ! position of velocity field in auxField
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems
    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / phyConvFac%body_force

!$omp parallel do schedule(static), private( posInTotal, forceTerm, elemOff )
    !NEC$ ivdep
    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)
      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! forceterm to add to velocity: F/(2*rho)
      forceTerm = forceField((iElem-1)*3+1 : iElem*3)*0.5_rk*rho0Inv
      ! add force to velocity
      auxField(elemOff+vel_pos(1)) = auxField(elemOff+vel_pos(1)) + forceTerm(1)
      auxField(elemOff+vel_pos(2)) = auxField(elemOff+vel_pos(2)) + forceTerm(2)
      auxField(elemOff+vel_pos(3)) = auxField(elemOff+vel_pos(3)) + forceTerm(3)
    end do

  end subroutine mus_addForceToAuxField_fluidIncomp
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine add body force to momentum in auxField for multispecies
  !! liquid model
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  subroutine mus_addForceToAuxField_MSL(fun, auxField, iLevel, time, varSys, &
    &                                   phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, mom_pos(3), depField
    real(kind=rk) :: forceTerm(3)
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, massFrac
    integer :: iElem, nElems, elemOff, nInputStates, iField
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    type(mus_varSys_data_type), pointer :: fPtr
    ! ------------------------------------------------------------------------ !
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! number of pdf states this source depends on
    ! last input is spacetime function so it is neglected
    nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / phyConvFac%body_force

    call c_f_pointer( varSys%method%val(1)%method_data, fPtr )
    associate( scheme => fPtr%solverData%scheme,            &
      &        posInTotal => fun%elemLvl(iLevel)%posInTotal )

!$omp parallel do schedule(static), private( forceTerm, elemOff )
      !NEC$ ivdep
      do iElem = 1, nElems
        ! element offset
        elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars

        ! forceterm to add to momentum: F/2
        forceTerm = forceField((iElem-1)*3+1 : iElem*3)*0.5_rk

        !NEC$ shortloop
        do iField = 1, scheme%nFields
          ! position of density and momentum field in auxField
          dens_pos = varSys%method%val(derVarPos(iField)%density) &
            &                     %auxField_varPos(1)
          ! mass density of species
          mass_dens(iField) = auxField( elemOff + dens_pos )
        end do

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! Update auxField depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        !NEC$ shortloop
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)
          ! Position of depField momentum in auxField array
          mom_pos = varSys%method%val(derVarPos(depField)%momentum) &
            &                    %auxField_varPos(1:3)

          ! add force to momentum
          auxField(elemOff+mom_pos(1)) = auxField(elemOff+mom_pos(1))    &
            &                          + massFrac(depField) * forceTerm(1)
          auxField(elemOff+mom_pos(2)) = auxField(elemOff+mom_pos(2))    &
            &                          + massFrac(depField) * forceTerm(2)
          auxField(elemOff+mom_pos(3)) = auxField(elemOff+mom_pos(3))    &
            &                          + massFrac(depField) * forceTerm(3)
        end do !iField
      end do !iElem
    end associate
  end subroutine mus_addForceToAuxField_MSL
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine add electric force to momentum in auxField for multispecies
  !! liquid model
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  subroutine mus_addElectricToAuxField_MSL(fun, auxField, iLevel, time, &
    &                                      varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, mom_pos(3), depField
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, massFrac
    integer :: iElem, nElems, elemOff, nInputStates, iField
    real(kind=rk) :: electricField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk) :: EF_elem(3)
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce
    real(kind=rk) :: charge_dens, diffForce_cs2inv, minMolWeight
    real(kind=rk), dimension(varSys%nStateVars) :: chargeTerm
    type(mus_varSys_data_type), pointer :: fPtr
    ! ------------------------------------------------------------------------ !
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! number of pdf states this source depends on
    ! last input is spacetime function so it is neglected
    nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

    call c_f_pointer( varSys%method%val(1)%method_data, fPtr )
    associate( scheme => fPtr%solverData%scheme,                            &
      &        posInTotal => fun%elemLvl(iLevel)%posInTotal,                &
      &        physics => fPtr%solverData%physics,                          &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
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

      ! minimum molecular weight
      minMolWeight = minval(species(:)%molWeight)

      ! constant term to multiply forcing term
      diffForce_cs2inv = minMolWeight / ( mixture%gasConst_R_LB &
        &              * mixture%temp0LB )

!$omp parallel do schedule(static), private( elemOff )
      !NEC$ ivdep
      do iElem = 1, nElems
        ! element offset
        elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars

        !NEC$ shortloop
        do iField = 1, scheme%nFields
          ! position of density and momentum field in auxField
          dens_pos = varSys%method%val(derVarPos(iField)%density) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField( elemOff + dens_pos )

          ! chargeTerm for each species: \rho_k z_k Faraday / M_k
          chargeTerm(iField) = mass_dens(iField)            &
            &                * species(iField)%molWeightInv &
            &                * species(iField)%chargeNr     &
            &                * mixture%faradayLB
        end do

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! compute charge density: \sum_k \rho_k z_k Faraday / M_k
        charge_dens = sum(chargeTerm)

        ! electric field for current element
        EF_elem = electricField((iElem-1)*3+1 : iElem*3) * 0.5_rk

        ! compute force on each species
        ! F_k = cs2 min(M) (\rho_k z_k/M_k - y_k \sum_l \rho_l z_l / M_l)
        !      F E / (RT)
        ! forceterm to add to momentum: F_k/2
        !NEC$ shortloop
        do iField = 1, scheme%nFields
          spcForce(:, iField) = EF_elem * cs2 * diffForce_cs2inv      &
            & * ( chargeTerm(iField) - massFrac(iField) * charge_dens )
        end do

        ! Update auxField depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        !NEC$ shortloop
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)
          ! Position of depField momentum in auxField array
          mom_pos = varSys%method%val(derVarPos(depField)%momentum) &
            &                    %auxField_varPos(1:3)

          ! add force to momentum
          auxField(elemOff+mom_pos(1)) = auxField(elemOff+mom_pos(1))    &
            &                          + spcForce(1, depField)
          auxField(elemOff+mom_pos(2)) = auxField(elemOff+mom_pos(2))    &
            &                          + spcForce(2, depField)
          auxField(elemOff+mom_pos(3)) = auxField(elemOff+mom_pos(3))    &
            &                          + spcForce(3, depField)
        end do !iField

      end do !iElem
    end associate

  end subroutine mus_addElectricToAuxField_MSL
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine add body force to momentum in auxField for multispecies
  !! liquid model with thermodynamic factor
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  subroutine mus_addForceToAuxField_MSL_WTDF(fun, auxField, iLevel, time, &
    &                                        varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, mom_pos(3), depField
    real(kind=rk) :: forceTerm(3)
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, massFrac
    real(kind=rk), dimension(varSys%nStateVars) :: num_dens, moleFrac
    integer :: iElem, nElems, elemOff, nInputStates, iField, iField_2
    real(kind=rk) :: forceField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: &
      & thermodynamic_fac, inv_thermodyn_fac
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce, spcForce_WTDF
    type(mus_varSys_data_type), pointer :: fPtr
    ! ------------------------------------------------------------------------ !
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! number of pdf states this source depends on
    ! last input is spacetime function so it is neglected
    nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = forceField                                )

    ! convert physical to lattice
    forceField = forceField / phyConvFac%body_force

    call c_f_pointer( varSys%method%val(1)%method_data, fPtr )
    associate( scheme => fPtr%solverData%scheme,                            &
      &        posInTotal => fun%elemLvl(iLevel)%posInTotal,                &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
      &        species => fPtr%solverData%scheme%field(:)%fieldProp%species )

!$omp parallel do schedule(static), private( forceTerm, elemOff )
      !NEC$ ivdep
      do iElem = 1, nElems
        ! element offset
        elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars

        !NEC$ shortloop
        do iField = 1, scheme%nFields
          ! position of density and momentum field in auxField
          dens_pos = varSys%method%val(derVarPos(iField)%density) &
            &                     %auxField_varPos(1)
          ! mass density of species
          mass_dens(iField) = auxField( elemOff + dens_pos )

          ! number density of species
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv
        end do

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! mole fraction
        moleFrac(:) = num_dens(:)/sum(num_dens)

        ! Thermodynamic factor from C++ code
        call mus_calc_thermFactor( nFields       = scheme%nFields,    &
          &                        temp          = mixture%temp0,     &
          &                        press         = mixture%atm_press, &
          &                        mole_frac     = moleFrac,          &
          &                        therm_factors = thermodynamic_fac  )

        ! invert thermodynamic factor
        inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

        ! forceterm to add to momentum: F/2
        forceTerm = forceField((iElem-1)*3+1 : iElem*3)*0.5_rk

        ! compute force on each species
        ! F_k =  y_k F
        do iField = 1, scheme%nFields
          spcForce(:, iField) = massFrac(iField) * forceTerm(:)
        end do

        ! compute external forcing term
        ! d^m_k = \gamma^{-1}_{k,l} F_k
        ! F_k is diffusive forcing term
        spcForce_WTDF = 0.0_rk
        do iField = 1, scheme%nFields
          do iField_2 = 1, scheme%nFields
            spcForce_WTDF(:, iField ) = spcForce_WTDF(:, iField)           &
              &                      + inv_thermodyn_fac(iField, iField_2) &
              &                      * spcForce(:, iField_2)
          end do
        end do

        ! Update auxField depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        !NEC$ shortloop
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)
          ! Position of depField momentum in auxField array
          mom_pos = varSys%method%val(derVarPos(depField)%momentum) &
            &                    %auxField_varPos(1:3)

          ! add force to momentum
          auxField(elemOff+mom_pos(1)) = auxField(elemOff+mom_pos(1)) &
            &                          + spcForce_WTDF(1, depField)
          auxField(elemOff+mom_pos(2)) = auxField(elemOff+mom_pos(2)) &
            &                          + spcForce_WTDF(2, depField)
          auxField(elemOff+mom_pos(3)) = auxField(elemOff+mom_pos(3)) &
            &                          + spcForce_WTDF(3, depField)
        end do !iField
      end do !iElem
    end associate
  end subroutine mus_addForceToAuxField_MSL_WTDF
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine add electric force to momentum in auxField for multispecies
  !! liquid model with thermodynamic factor
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  subroutine mus_addElectricToAuxField_MSL_WTDF(fun, auxField, iLevel, time, &
    &                                           varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, mom_pos(3), depField
    real(kind=rk), dimension(varSys%nStateVars) :: mass_dens, massFrac
    real(kind=rk), dimension(varSys%nStateVars) :: num_dens, moleFrac
    integer :: iElem, nElems, elemOff, nInputStates, iField, iField_2
    real(kind=rk) :: electricField(fun%elemLvl(iLevel)%nElems*3)
    real(kind=rk), dimension(varSys%nStateVars, varSys%nStateVars) :: &
      & thermodynamic_fac, inv_thermodyn_fac
    real(kind=rk), dimension(3, varSys%nStateVars ) :: spcForce, spcForce_WTDF
    real(kind=rk) :: charge_dens, diffForce_cs2inv, minMolWeight
    real(kind=rk), dimension(varSys%nStateVars) :: chargeTerm
    type(mus_varSys_data_type), pointer :: fPtr
    ! ------------------------------------------------------------------------ !
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! number of pdf states this source depends on
    ! last input is spacetime function so it is neglected
    nInputStates = varSys%method%val(fun%srcTerm_varPos)%nInputs - 1

    call c_f_pointer( varSys%method%val(1)%method_data, fPtr )
    associate( scheme => fPtr%solverData%scheme,                            &
      &        posInTotal => fun%elemLvl(iLevel)%posInTotal,                &
      &        physics => fPtr%solverData%physics,                          &
      &        mixture => fPtr%solverData%scheme%mixture,                   &
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

      ! minimum molecular weight
      minMolWeight = minval(species(:)%molWeight)

      ! constant term to multiply forcing term
      diffForce_cs2inv = minMolWeight / ( mixture%gasConst_R_LB &
        &           * mixture%temp0LB )

!$omp parallel do schedule(static), private( elemOff )
      !NEC$ ivdep
      do iElem = 1, nElems
        ! element offset
        elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars

        !NEC$ shortloop
        do iField = 1, scheme%nFields
          ! position of density and momentum field in auxField
          dens_pos = varSys%method%val(derVarPos(iField)%density) &
            &                     %auxField_varPos(1)

          ! mass density of species
          mass_dens(iField) = auxField( elemOff + dens_pos )

          ! number density of species
          num_dens(iField) = mass_dens(iField) * species(iField)%molWeightInv

          ! chargeTerm for each species: \rho_k z_k Faraday / M_k
          chargeTerm(iField) = num_dens(iField) * species(iField)%chargeNr &
            &                * mixture%faradayLB
        end do

        !mass fraction
        massFrac(:) = mass_dens(:)/sum(mass_dens)

        ! compute charge density: \sum_k \rho_k z_k Faraday / M_k
        charge_dens = sum(chargeTerm)

        ! compute force on each species
        ! F_k = cs2 min(M) (\rho_k z_k/M_k - y_k \sum_l \rho_l z_l / M_l)
        !      F E / (RT)
        ! forceterm to add to momentum: F_k/2
        !NEC$ shortloop
        do iField = 1, scheme%nFields
          spcForce(:, iField) = electricField((iElem-1)*3+1 : iElem*3)  &
            & * ( chargeTerm(iField) - massFrac(iField) * charge_dens ) &
            & * diffForce_cs2inv * 0.5_rk * cs2
        end do

        ! mole fraction
        moleFrac(:) = num_dens(:)/sum(num_dens)

        ! Thermodynamic factor from C++ code
        call mus_calc_thermFactor( nFields       = scheme%nFields,    &
          &                        temp          = mixture%temp0,     &
          &                        press         = mixture%atm_press, &
          &                        mole_frac     = moleFrac,          &
          &                        therm_factors = thermodynamic_fac  )

        ! invert thermodynamic factor
        inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

        ! compute external forcing term
        ! d^m_k = \gamma^{-1}_{k,l} F_k
        ! F_k is diffusive forcing term
        spcForce_WTDF = 0.0_rk
        do iField = 1, scheme%nFields
          do iField_2 = 1, scheme%nFields
            spcForce_WTDF(:, iField ) = spcForce_WTDF(:, iField)           &
              &                      + inv_thermodyn_fac(iField, iField_2) &
              &                      * spcForce(:, iField_2)
          end do
        end do

        ! Update auxField depends on nInputStates
        ! if nInputStates = 1, it is field source else it is global source
        !NEC$ shortloop
        do iField = 1, nInputStates
          depField = varSys%method%val(fun%srcTerm_varPos)%input_varPos(iField)
          ! Position of depField momentum in auxField array
          mom_pos = varSys%method%val(derVarPos(depField)%momentum) &
            &                    %auxField_varPos(1:3)

          ! add force to momentum
          auxField(elemOff+mom_pos(1)) = auxField(elemOff+mom_pos(1))    &
            &                          + spcForce_WTDF(1, depField)
          auxField(elemOff+mom_pos(2)) = auxField(elemOff+mom_pos(2))    &
            &                          + spcForce_WTDF(2, depField)
          auxField(elemOff+mom_pos(3)) = auxField(elemOff+mom_pos(3))    &
            &                          + spcForce_WTDF(3, depField)
        end do !iField

      end do !iElem
    end associate

  end subroutine mus_addElectricToAuxField_MSL_WTDF
  ! ************************************************************************** !


! ************************************************************************** !
  !> This routine add source term with charge density in the Poisson equation
  !! to the potential.
  !! Refer to Appendix in PhD Thesis of K. Masilamani
  !! "Coupled Simulation Framework to Simulate Electrodialysis Process for
  !! Seawater Desalination"
  subroutine mus_addSrcToAuxField_poisson(fun, auxField, iLevel, time, &
    &                                            varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: iElem, nElems
    real(kind=rk) :: rhs(fun%elemLvl(iLevel)%nElems), rhs_Fac
    type(mus_varSys_data_type), pointer :: fPtr
    ! ------------------------------------------------------------------------ !
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    call c_f_pointer( varSys%method%val(fun%srcTerm_varPos)%method_data, fPtr )
    associate( posInTotal => fun%elemLvl(iLevel)%posInTotal,                 &
      &        poisson => fPtr%solverData%scheme%field(1)%fieldProp%poisson, &
      &        physics => fPtr%solverData%physics                            )

      ! factor to multiply rhs
      rhs_Fac = 0.5_rk * poisson%pot_diff / poisson%permittivity

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

!$omp parallel do schedule(static)
      !NEC$ ivdep
      do iElem = 1, nElems
        ! add source term to potential
        auxField(posInTotal(iElem)) = auxField(posInTotal(iElem)) + rhs(iElem)
      end do !iElem
    end associate

  end subroutine mus_addSrcToAuxField_poisson
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine add sponge density and velocity field to density and velocity
  !! in auxField for weakly-compressible model.
  !! Reference:
  !! Jacob, J.; Sagaut, P. (2019): Solid wall and open boundary conditions in
  !! hybrid recursive regularized lattice Boltzmann method for compressible
  !! flows. In Physics of Fluids 31 (12), p. 126103.
  subroutine mus_addSponFldToAuxField_fluid(fun, auxField, iLevel, time, &
    &                                       varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, vel_pos(3)
    real(kind=rk) :: sigma
    real(kind=rk) :: inv_rho_phy, inv_vel_phy
    integer :: iElem, nElems, elemOff
    real(kind=rk) :: spongeField(fun%elemLvl(iLevel)%nElems)
    real(kind=rk) :: dens, vel(3)
    real(kind=rk) :: dens_ref, vel_ref(3), sponDens, sponVel(3)
    ! ------------------------------------------------------------------------ !
    ! position of density and velocity field in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems
    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = spongeField                               )

    inv_rho_phy = 1.0_rk / phyConvFac%press * cs2inv
    inv_vel_phy = 1.0_rk / phyConvFac%vel

    ! target pressure and velocity in lattice unit
    dens_ref = fun%absLayer%config%target_pressure * inv_rho_phy
    vel_ref(1:3) = fun%absLayer%config%target_velocity(1:3) * inv_vel_phy

    associate(posInTotal => fun%elemLvl(iLevel)%posInTotal)
      !NEC$ ivdep
      do iElem = 1, nElems
        ! element offset
        elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars
        ! Local density and velocity
        dens = auxField(elemOff+dens_pos)
        vel(1) = auxField(elemOff+vel_pos(1))
        vel(2) = auxField(elemOff+vel_pos(2))
        vel(3) = auxField(elemOff+vel_pos(3))

        ! SpongeField contains: spongeStrength
        sigma = spongeField(iElem)

        ! Sponge factor for density and velocity field
        sponDens = -sigma*(dens - dens_ref)*0.5_rk
        sponVel(:) = -sigma*(vel - vel_ref)*0.5_rk

        ! add force to velocity
        auxField(elemOff+dens_pos) = dens + sponDens
        auxField(elemOff+vel_pos(1)) = vel(1) + sponVel(1)
        auxField(elemOff+vel_pos(2)) = vel(2) + sponVel(2)
        auxField(elemOff+vel_pos(3)) = vel(3) + sponVel(3)
      end do
    end associate

  end subroutine mus_addSponFldToAuxField_fluid
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine add sponge density and velocity field to density and velocity
  !! in auxField. Density and velocity in far field are computed by time
  !! average.
  !!
  !! Reference:
  !! Jacob, J.; Sagaut, P. (2019): Solid wall and open boundary conditions in
  !! hybrid recursive regularized lattice Boltzmann method for compressible
  !! flows. In Physics of Fluids 31 (12), p. 126103.
  subroutine mus_addDynSponFldToAuxField_fluid(fun, auxField, iLevel, time, &
    &                                          varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, vel_pos(3)
    real(kind=rk) :: sigma
    integer :: iElem, nElems,  elemOff
    real(kind=rk) :: spongeField(fun%elemLvl(iLevel)%nElems)
    real(kind=rk) :: dens, vel(3)
    real(kind=rk) :: densAvgNew, velAvgNew(3)
    real(kind=rk) :: sponDens, sponVel(3)
    ! ------------------------------------------------------------------------ !
    ! position of density and velocity field in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems
    ! Get force which is refered in config file either its
    ! spacetime variable or operation variable
    call varSys%method%val(fun%data_varPos)%get_valOfIndex( &
      & varSys  = varSys,                                   &
      & time    = time,                                     &
      & iLevel  = iLevel,                                   &
      & idx     = fun%elemLvl(iLevel)%idx(1:nElems),        &
      & nVals   = nElems,                                   &
      & res     = spongeField                               )

    associate( dynAvg => fun%elemLvl(iLevel)%dynAvg,         &
      &        posInTotal => fun%elemLvl(iLevel)%posInTotal  )
!$omp parallel do schedule(static), private( elemOff )
      !NEC$ ivdep
      do iElem = 1, nElems
        ! element offset
        elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars
        ! Local density and velocity
        dens = auxField(elemOff+dens_pos)
        vel(1) = auxField(elemOff+vel_pos(1))
        vel(2) = auxField(elemOff+vel_pos(2))
        vel(3) = auxField(elemOff+vel_pos(3))

        !! New avg
        densAvgNew = dynAvg%dens(iElem)
        velAvgNew(1) = dynAvg%velX(iElem)
        velAvgNew(2) = dynAvg%velY(iElem)
        velAvgNew(3) = dynAVg%velZ(iElem)

        ! SpongeField contains: spongeStrength
        sigma = spongeField(iElem)

        ! Sponge factor for density and velocity field
        sponDens = -sigma*(dens - densAvgNew)*0.5_rk
        sponVel(:) = -sigma*(vel - velAvgNew)*0.5_rk

        ! add force to velocity
        auxField(elemOff+dens_pos) = dens + sponDens
        auxField(elemOff+vel_pos(1)) = vel(1) + sponVel(1)
        auxField(elemOff+vel_pos(2)) = vel(2) + sponVel(2)
        auxField(elemOff+vel_pos(3)) = vel(3) + sponVel(3)

      end do
    end associate

  end subroutine mus_addDynSponFldToAuxField_fluid
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine add sponge density and velocity field to density and velocity
  !! in auxField. Density and velocity in far field are computed by time
  !! average.
  !!
  !! Reference:
  !! Jacob, J.; Sagaut, P. (2019): Solid wall and open boundary conditions in
  !! hybrid recursive regularized lattice Boltzmann method for compressible
  !! flows. In Physics of Fluids 31 (12), p. 126103.
  subroutine mus_addHRRCorrToAuxField_fluid_2D(fun, auxField, iLevel, time, &
    &                                          varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, vel_pos(3)
    integer :: nSolve, iElem,  elemOff
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! ------------------------------------------------------------------------ !
    ! Number of elements to apply source terms
    nSolve = fun%elemLvl(iLevel)%nElems

    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    ! number od dimensions
    ! position of density and velocity field in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    ! Calculate phi
    associate( HRR_Corr => fun%elemLvl(iLevel)%HRR_Corr  )

      do iChunks = 1, nChunks
        ! calculate the end  number of iElem loop
        nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
        low_bound = (iChunks - 1) * vlen

        do iElem = 1, nChunkElems

          elemPos = low_bound + iElem

          ! element offset
          elemoff = (elemPos-1)*varSys%nAuxScalars
          ! source term + Local density and velocity
          auxField(elemOff+dens_pos) = 0.5_rk * HRR_Corr%dens(elemPos) + auxField(elemOff+dens_pos)
          auxField(elemOff+vel_pos(1)) = 0.5_rk * HRR_Corr%vel(elemPos,1) + auxField(elemOff+vel_pos(1))
          auxField(elemOff+vel_pos(2)) = 0.5_rk * HRR_Corr%vel(elemPos,2) + auxField(elemOff+vel_pos(2))

        enddo
      enddo

    end associate

  end subroutine mus_addHRRCorrToAuxField_fluid_2D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine add sponge density and velocity field to density and velocity
  !! in auxField. Density and velocity in far field are computed by time
  !! average.
  !!
  !! Reference:
  !! Jacob, J.; Sagaut, P. (2019): Solid wall and open boundary conditions in
  !! hybrid recursive regularized lattice Boltzmann method for compressible
  !! flows. In Physics of Fluids 31 (12), p. 126103.
  subroutine mus_addHRRCorrToAuxField_fluid_3D(fun, auxField, iLevel, time, &
    &                                          varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, vel_pos(3)
    integer :: nSolve, iElem, elemOff
    integer :: nChunks, iChunks, nChunkElems, low_bound, elemPos
    ! ------------------------------------------------------------------------ !
    ! Number of elements to apply source terms
    nSolve = fun%elemLvl(iLevel)%nElems

    nChunks = ceiling(real(nSolve, kind=rk) / real(vlen, kind=rk))

    ! number od dimensions
    ! position of density and velocity field in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    ! Calculate phi
    associate( HRR_Corr => fun%elemLvl(iLevel)%HRR_Corr )

      do iChunks = 1, nChunks
        ! calculate the end  number of iElem loop
        nChunkElems = min(vlen, nSolve - ((iChunks - 1) * vlen))
        low_bound = (iChunks - 1) * vlen

        do iElem = 1, nChunkElems

          elemPos = low_bound + iElem

          ! element offset
          elemoff = (elemPos-1)*varSys%nAuxScalars
          ! source term + Local density and velocity
          auxField(elemOff+dens_pos) = 0.5_rk * HRR_Corr%dens(elemPos) + auxField(elemOff+dens_pos)
          auxField(elemOff+vel_pos(1)) = 0.5_rk * HRR_Corr%vel(elemPos,1) + auxField(elemOff+vel_pos(1))
          auxField(elemOff+vel_pos(2)) = 0.5_rk * HRR_Corr%vel(elemPos,2) + auxField(elemOff+vel_pos(2))
          auxField(elemOff+vel_pos(3)) = 0.5_rk * HRR_Corr%vel(elemPos,3) + auxField(elemOff+vel_pos(3))

        enddo
      enddo

    end associate

  end subroutine mus_addHRRCorrToAuxField_fluid_3D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine add dynamic force to velocity in auxField for
  !! weakly-compressible model for turbulent channel test case.
  !! Force definition:
  !! Force = rho*u_tau^2/H + rho*(u_bulk_ref-uX_bulk_avg)*u_bulk_ref/H
  !! Reference:
  !! 1) https://www.wias-berlin.de/people/john/ELECTRONIC_PAPERS/JR07.IJNMF.pdf
  !! 2) Haussmann, Marc; BARRETO, Alejandro CLARO; KOUYI, Gislain LIPEME;
  !! Rivière, Nicolas; Nirschl, Hermann; Krause, Mathias J. (2019):
  !! Large-eddy simulation coupled with wall models for turbulent channel flows
  !! at high Reynolds numbers with a lattice Boltzmann method — Application to
  !! Coriolis mass flowmeter. In Computers & Mathematics with Applications 78
  !! (10), pp. 3285–3302. DOI: 10.1016/j.camwa.2019.04.033.
  subroutine mus_addTurbChanForceToAuxField_fluid(fun, auxField, iLevel, time, &
    &                                             varSys, phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    integer :: dens_pos, vel_pos(3)
    real(kind=rk) :: forceTerm(3)
    integer :: iElem, nElems, posInTotal, elemOff
    real(kind=rk) :: forceDynL(3)
    ! ------------------------------------------------------------------------ !

    ! position of density and velocity field in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems

    ! Convert dynamic force term in m/s^2 to lattice unit.
    forceDynL = fun%turbChanForce%forceDyn / phyConvFac%accel
    !write(dbgunit(1), *) 'forceDynL: ', forceDynL

!$omp parallel do schedule(static), private( posInTotal, forceTerm, elemOff )
    !NEC$ ivdep
    do iElem = 1, nElems
      posInTotal = fun%elemLvl(iLevel)%posInTotal(iElem)
      ! element offset
      elemoff = (posInTotal-1)*varSys%nAuxScalars
      ! forceterm to add to velocity: F/2
      forceTerm = forceDynL * 0.5_rk
      ! add force to velocity i.e. u(i) = u(i) + F(i)/2
      auxField(elemOff+vel_pos(1)) = auxField(elemOff+vel_pos(1)) + forceTerm(1)
      auxField(elemOff+vel_pos(2)) = auxField(elemOff+vel_pos(2)) + forceTerm(2)
      auxField(elemOff+vel_pos(3)) = auxField(elemOff+vel_pos(3)) + forceTerm(3)
    end do

  end subroutine mus_addTurbChanForceToAuxField_fluid
  ! ************************************************************************** !


end module mus_auxFieldVar_module

