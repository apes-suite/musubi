! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> This module contains data types, function and routines for gradient
!! computation.
!!
!! author: Gregorio Gerardo Spinelli
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
module mus_gradData_module
  ! include treelm modules
  use env_module,                 only: labelLen, rk
  use tem_debug_module,           only: dbgUnit
  use tem_logging_module,         only: logUnit
  use tem_aux_module,             only: tem_abort
  use tem_stencil_module,         only: tem_stencilHeader_type
  use tem_param_module,           only: qN00, q100, q0N0, q010, q00N, q001
  use tem_construction_module,    only: tem_levelDesc_type

  implicit none
  private

  public :: mus_Grad_type
  public :: mus_gradData_type
  public :: mus_init_gradData
  public :: mus_assign_GradCalculation_ptr
  public :: getGradRhoU3
  public :: getGradRhoUVZ
  public :: getGradU
  public :: getGradXXU

  !> collection of properties of the Gradient type
  type mus_Grad_type

    !> function pointer to get Gradient U
    procedure(getGradU), nopass, pointer :: U_ptr => null()
    !> function pointer to get Gradient XX U
    procedure(getGradXXU), nopass, pointer :: XXU_ptr => null()
    !> function pointer to get Gradient rho U^3
    procedure(getGradRhoU3), nopass, pointer :: RhoU3_ptr => null()
    !> function pointer to get Gradient Rho u_x u_y u_z
    procedure(getGradRhoUVZ), nopass, pointer :: RhoUVZ_ptr => null()

  end type mus_Grad_type

  !> Contains information required to compute gradient like position
  !! of six direct neigbors in the state array and coeff for the gradient
  !! operation
  type mus_gradData_type
    !> Stores position of 6 direct face neighbors in the state array for
    !! each element
    !! Size: nSize, stencil%nDims, 2
    !! last index refers to left and right side
    !! 1 is left/negative side and 2 is right/positive side
    integer, allocatable :: neighPos(:,:,:)

    !> coeff to calculate 1st order derivative
    !! use forward difference for element with boundary (coeff = 1.0)
    !! and central difference for inner fluid element (coeff=0.5)
    !! size: nSize, stencil%nDims
    real(kind=rk), allocatable :: FDcoeff(:,:)

  end type mus_gradData_type

  abstract interface
    !> function pointers to obtain full gradient of U along all
    !! directions
    pure function getGradU( auxField, gradData, velPos, nAuxScalars, &
      &                     nDims, nSolve, elemOffset ) result(gradU)
      import :: rk, mus_gradData_type

      !> auxField
      real(kind=rk), intent(in) :: auxField(:)
      !> Number of element to solve in this level
      integer, intent(in) :: nSolve
      !> gradient data
      type(mus_gradData_type), intent(in) :: gradData
      !> Position of velocity field in auxField
      integer, intent(in) :: velPos(3)
      !> Number of scalars in auxField array
      integer, intent(in) :: nAuxScalars
      !> Offset for elements when computing chunkwise
      integer, intent(in) :: elemOffset
      !> Dimensions
      integer, intent(in) :: nDims
      !> output is gradient of velocity
      real(kind=rk) :: gradU(nDims,nDims,nSolve)

    end function getGradU

    !> function pointers to obtain full gradient of U along main
    !! directions derived two times d^2U/dx^2
    pure function getGradXXU( auxField, gradData, velPos, nAuxScalars, &
      &                       nDims, nSolve, elemOffset ) result(gradXXU)
      import :: rk, mus_gradData_type

      !> auxField
      real(kind=rk), intent(in) :: auxField(:)
      !> Number of element to solve in this level
      integer, intent(in) :: nSolve
      !> gradient data
      type(mus_gradData_type), intent(in) :: gradData
      !> Position of velocity field in auxField
      integer, intent(in) :: velPos(3)
      !> Number of scalars in auxField array
      integer, intent(in) :: nAuxScalars
      !> Offset for elements when computing chunkwise
      integer, intent(in) :: elemOffset
      !> Dimensions
      integer, intent(in) :: nDims
      !> output is gradient of velocity:
      ! 1: Dxxu, 2: Dyyv, 3: Dzzw
      real(kind=rk) :: gradXXU(nDims,nSolve)

    end function getGradXXU

    !> function pointers to obtain grad Rho * U^3 along main directions
    pure function getGradRhoU3( auxField, gradData, velPos, densPos, nAuxScalars, &
      &                         nDims, nSolve, elemOffset ) result(gradRhoU3)
      import :: rk, mus_gradData_type

      !> auxField
      real(kind=rk), intent(in) :: auxField(:)
      !> Number of element to solve in this level
      integer, intent(in) :: nSolve
      !> gradient data
      type(mus_gradData_type), intent(in) :: gradData
      !> Position of velocity field in auxField
      integer, intent(in) :: velPos(3), densPos
      !> Number of scalars in auxField array
      integer, intent(in) :: nAuxScalars
      !> Offset for elements when computing chunkwise
      integer, intent(in) :: elemOffset
      !> Dimensions
      integer, intent(in) :: nDims
      !> output is gradient of velocity
      real(kind=rk) :: gradRhoU3(nDims,nSolve)

    end function getGradRhoU3

    !> function pointers to obtain grad Rho * u_x * u_y * u_z
    !! along main directions
    pure function getGradRhoUVZ( auxField, gradData, velPos, densPos, &
      &                          nDims, nAuxScalars, nSolve, elemOffset )    &
      &  result(gradRhoUVZ)
      import :: rk, mus_gradData_type

      !> auxField
      real(kind=rk), intent(in) :: auxField(:)
      !> Number of element to solve in this level
      integer, intent(in) :: nSolve
      !> gradient data
      type(mus_gradData_type), intent(in) :: gradData
      !> Position of velocity field in auxField
      integer, intent(in) :: velPos(3), densPos
      !> Number of scalars in auxField array
      integer, intent(in) :: nAuxScalars
      !> Offset for elements when computing chunkwise
      integer, intent(in) :: elemOffset
      !> Dimensions
      integer, intent(in) :: nDims
      !> output is gradient of velocity
      real(kind=rk) :: gradRhoUVZ(nDims,nSolve)

    end function getGradRhoUVZ
  end interface

contains

  ! ************************************************************************* !
  !> This routine initialize gradData with direct neighbors in state and
  !! finite difference coefficients.
  subroutine mus_init_gradData(me, neigh, & !levelDesc,
    &                          stencil, nSize, nSolve, &
    &                          nScalars)
    !---------------------------------------------------------------------------
    !> Gradient type
    type(mus_gradData_type), intent(out) :: me
    !> neighbor connectivity array
    integer, intent(in) :: neigh(:)
    !!> levelDesc to access communication buffers of state array
    !type(tem_levelDesc_type), intent(in) :: levelDesc
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> Number of elements in state array
    integer, intent(in) :: nSize
    !> Number of elements solved in compute kernel i.e. excluding halos
    integer, intent(in) :: nSolve
    !> number of scalars in state array
    integer, intent(in) :: nScalars
    !---------------------------------------------------------------------------
    integer :: nDims
    integer :: iElem, iDir, iFace, iMeshDir, iD
    integer :: nFaces
    integer :: nghElem
    logical :: hasBnd(3,2)
    !---------------------------------------------------------------------------
    write(logUnit(3),*) 'Allocate gradient data ...'
    nDims = stencil%nDims
    ! allocate gradient data type
    allocate(me%neighPos( nSolve, nDims, 2 ) )
    me%neighPos = -1
    allocate(me%FDcoeff( nSolve, nDims ) )
    me%FDcoeff = -1.0_rk

    write(logUnit(1),*) 'Filling gradient data ...'
    ! fill gradient data with face neighbor information for turbulence model
!write(dbgUnit(1),*) 'Filling gradient data ...'
    ! Number of faces per element
    nFaces = nDims*2

    ! coeff to calculate 1st order derivative
    ! use forward difference for element with boundary (coeff = 1.0)
    ! and central difference for inner fluid element (coeff=0.5)
    me%FDcoeff(:,:) = 0.5_rk

    elemLoop: do iElem = 1, nSolve
!write(dbgUnit(1),*) 'iElem: ', iElem, ' treeID: ', levelDesc%total(iElem)

      iFace = 0
      ! set to true for the direction it has boundary.
      ! This variable is used to set coeff for gradient calculation.
      ! if neighbor element is my element then my neighbor is boundary
      ! use Forward difference to approximate gradient
      hasBnd = .false.

      ! Loop over every link
      neighLoop: do iDir  = 1, stencil%QQN
!write(dbgUnit(1),*) 'iDir: ', iDir

        ! Find neighor using neigh array because for boundary, neigh array
        ! points to current element so no need to check for existence of the
        ! neighbor element
        nghElem = int((neigh(( stencil%cxdirinv( idir)-1)* nsize+ ielem)-1)/ nscalars)+1

        ! map neighbor direction to treelm stencil to identify the face
        iMeshDir = stencil%map(iDir)
        select case(iMeshDir)
        case (qN00) !-x
          iFace = iFace+1
          me%neighPos(iElem,1,1) = nghElem
          if (nghElem == iElem) hasBnd(1,1) = .true.
        case (q100) !+x
          iFace = iFace+1
          me%neighPos(iElem,1,2) = nghElem
          if (nghElem == iElem) hasBnd(1,2) = .true.
        case (q0N0) !-y
          iFace = iFace+1
          me%neighPos(iElem,2,1) = nghElem
          if (nghElem == iElem) hasBnd(2,1) = .true.
        case (q010) !+y
          iFace = iFace+1
          me%neighPos(iElem,2,2) = nghElem
          if (nghElem == iElem) hasBnd(2,2) = .true.
        case (q00N) !-z
          iFace = iFace+1
          me%neighPos(iElem,3,1) = nghElem
          if (nghElem == iElem) hasBnd(3,1) = .true.
        case (q001) !+z
          iFace = iFace+1
          me%neighPos(iElem,3,2) = nghElem
          if (nghElem == iElem) hasBnd(3,2) = .true.
        end select
        ! exit iDir loop
        if (iFace == nFaces) exit neighLoop
      end do neighLoop

      ! Update FDCoeff for the direction it has boundary
      do iD = 1, nDims
        if (any(hasBnd(iD,:))) me%FDcoeff(iElem,iD) = 1.0_rk
      end do
!write(dbgUnit(1),*) 'hasBnd: ', hasBnd
!write(dbgUnit(1),*) 'coeff: ', me%FDcoeff(iElem,:)
    end do elemLoop

!      write(logUnit(1),*) 'Done with fill grad data.'
!      write(dbgUnit(1),*) 'Done with fill grad data.'
!      flush(dbgUnit(1))
!      call tem_abort('Verify fill gradData')
  end subroutine mus_init_gradData
  ! ************************************************************************* !

! ************************************************************************** !
!> This function returns function pointer of nonEquilibrium scaling
!! for interpolation according to scheme definition
function mus_assign_GradCalculation_ptr(label) result(Grad)
  ! --------------------------------------------------------------------------
  !> Scheme header information
  character(len=labelLen), intent(in) :: label
  !> GradRhoU3 function
  type(mus_Grad_type) :: Grad
  ! --------------------------------------------------------------------------
  Grad%RhoU3_ptr => null()
  Grad%RhoUVZ_ptr => null()
  Grad%U_ptr => null()
  Grad%XXU_ptr => null()
  select case (trim(label))
  case ('d1q3')
    Grad%XXU_ptr => getGradXXU_1D
    Grad%RhoU3_ptr => getGradRhoU3_1D
    Grad%U_ptr => getGradU_1D
  case ('d2q9')
    Grad%XXU_ptr => getGradXXU_2D
    Grad%RhoU3_ptr => getGradRhoU3_2D
    Grad%U_ptr => getGradU_2D
  case ('d3q6', 'd3q15', 'd3q19', 'd3q27')
    Grad%XXU_ptr => getGradXXU_3D
    Grad%RhoU3_ptr => getGradRhoU3_3D
    Grad%RhoUVZ_ptr => getGradRhoUVZ_3D
    Grad%U_ptr => getGradU_3D
  case default
    write(*,*) 'stencil label = "', trim(label), '"'
    call tem_abort('Error: getGradAny not implemented for the given stencil')
  end select

end function mus_assign_GradCalculation_ptr
! ************************************************************************** !

  ! ************************************************************************** !
  !> This function computes gradient of velocity from gradient and veleocity
  !! data.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used
  pure function getGradU_1D( auxField, gradData, velPos, nAuxScalars, &
    &                        nDims, nSolve, elemOffset ) result(gradU)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3)
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity
    real(kind=rk) :: gradU(nDims,nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(1), rightngh(1)
    real(kind=rk) :: leftvel(1,1), rightvel(1,1)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
      elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      leftvel(1,1) = auxfield( leftngh(1)+velpos(1) )
      rightvel(1,1) = auxfield( rightngh(1)+velpos(1) )
      gradU(1,1,iElem) = (rightvel(1,1)-leftvel(1,1)) &
        &                 * graddata%fdcoeff(elempos, 1)
    end do

  end function getGradU_1D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This function computes gradient of velocity from gradient and veleocity
  !! data.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used
  pure function getGradU_2D( auxField, gradData, velPos, nAuxScalars, &
    &                        nDims, nSolve, elemOffset ) result(gradU)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3)
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity
    real(kind=rk) :: gradU(nDims,nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(2), rightngh(2)
    real(kind=rk) :: leftvel(2,2), rightvel(2,2)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
      elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      leftngh(2) = (graddata%neighpos(elempos, 2, 1) - 1) * nAuxScalars
      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      rightngh(2) = (graddata%neighpos(elempos, 2, 2) - 1) * nAuxScalars
      leftvel(1,1) = auxfield( leftngh(1)+velPos(1) )
      leftvel(1,2) = auxfield( leftngh(2)+velPos(1) )
      leftvel(2,1) = auxfield( leftngh(1)+velPos(2) )
      leftvel(2,2) = auxfield( leftngh(2)+velPos(2) )
      rightvel(1,1) = auxfield( rightngh(1)+velPos(1) )
      rightvel(1,2) = auxfield( rightngh(2)+velPos(1) )
      rightvel(2,1) = auxfield( rightngh(1)+velPos(2) )
      rightvel(2,2) = auxfield( rightngh(2)+velPos(2) )
      gradU(1,1,iElem) = (rightvel(1,1)-leftvel(1,1)) &
        &                 * graddata%fdcoeff(elempos, 1)
      gradU(1,2,iElem) = (rightvel(1,2)-leftvel(1,2)) &
        &                 * graddata%fdcoeff(elempos, 2)
      gradU(2,1,iElem) = (rightvel(2,1)-leftvel(2,1)) &
        &                 * graddata%fdcoeff(elempos, 1)
      gradU(2,2,iElem) = (rightvel(2,2)-leftvel(2,2)) &
        &                 * graddata%fdcoeff(elempos, 2)
    end do

  end function getGradU_2D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function computes gradient of velocity from gradient and veleocity
  !! data.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used
  pure function getGradU_3D( auxField, gradData, velPos, nAuxScalars, &
    &                        nDims, nSolve, elemOffset ) result(gradU)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3)
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity
    real(kind=rk) :: gradU(nDims,nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(3), rightngh(3)
    real(kind=rk) :: leftvel(3,3), rightvel(3,3)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
      elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      leftngh(2) = (graddata%neighpos(elempos, 2, 1) - 1) * nAuxScalars
      leftngh(3) = (graddata%neighpos(elempos, 3, 1) - 1) * nAuxScalars
      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      rightngh(2) = (graddata%neighpos(elempos, 2, 2) - 1) * nAuxScalars
      rightngh(3) = (graddata%neighpos(elempos, 3, 2) - 1) * nAuxScalars
      leftvel(1,1) = auxfield( leftngh(1)+velPos(1) )
      leftvel(2,1) = auxfield( leftngh(1)+velPos(2) )
      leftvel(3,1) = auxfield( leftngh(1)+velPos(3) )
      leftvel(1,2) = auxfield( leftngh(2)+velPos(1) )
      leftvel(2,2) = auxfield( leftngh(2)+velPos(2) )
      leftvel(3,2) = auxfield( leftngh(2)+velPos(3) )
      leftvel(1,3) = auxfield( leftngh(3)+velPos(1) )
      leftvel(2,3) = auxfield( leftngh(3)+velPos(2) )
      leftvel(3,3) = auxfield( leftngh(3)+velPos(3) )
      rightvel(1,1) = auxfield( rightngh(1)+velPos(1) )
      rightvel(2,1) = auxfield( rightngh(1)+velPos(2) )
      rightvel(3,1) = auxfield( rightngh(1)+velPos(3) )
      rightvel(1,2) = auxfield( rightngh(2)+velPos(1) )
      rightvel(2,2) = auxfield( rightngh(2)+velPos(2) )
      rightvel(3,2) = auxfield( rightngh(2)+velPos(3) )
      rightvel(1,3) = auxfield( rightngh(3)+velPos(1) )
      rightvel(2,3) = auxfield( rightngh(3)+velPos(2) )
      rightvel(3,3) = auxfield( rightngh(3)+velPos(3) )
      gradU(1,1,iElem) = (rightvel(1,1)-leftvel(1,1)) &
        &                 * graddata%fdcoeff(elempos, 1)
      gradU(2,1,iElem) = (rightvel(2,1)-leftvel(2,1)) &
        &                 * graddata%fdcoeff(elempos, 1)
      gradU(3,1,iElem) = (rightvel(3,1)-leftvel(3,1)) &
        &                 * graddata%fdcoeff(elempos, 1)
      gradU(1,2,iElem) = (rightvel(1,2)-leftvel(1,2)) &
        &                 * graddata%fdcoeff(elempos, 2)
      gradU(2,2,iElem) = (rightvel(2,2)-leftvel(2,2)) &
        &                 * graddata%fdcoeff(elempos, 2)
      gradU(3,2,iElem) = (rightvel(3,2)-leftvel(3,2)) &
        &                 * graddata%fdcoeff(elempos, 2)
      gradU(1,3,iElem) = (rightvel(1,3)-leftvel(1,3)) &
        &                 * graddata%fdcoeff(elempos, 3)
      gradU(2,3,iElem) = (rightvel(2,3)-leftvel(2,3)) &
        &                 * graddata%fdcoeff(elempos, 3)
      gradU(3,3,iElem) = (rightvel(3,3)-leftvel(3,3)) &
        &                 * graddata%fdcoeff(elempos, 3)
    end do

  end function getGradU_3D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This function computes gradient of velocity from gradient and veleocity
  !! data.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used
  pure function getGradXXU_1D( auxField, gradData, velPos, nAuxScalars, &
    &                          nDims, nSolve, elemOffset ) result(gradXXU)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3)
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity:
    ! 1: Dxxu, 2: Dyyv, 3: Dzzw
    real(kind=rk) :: gradXXU(nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(1), rightngh(1)
    real(kind=rk) :: leftvel(1), rightvel(1), vel(1)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
      elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      vel(1) = 2._rk * auxField( ElemOffset * nAuxScalars + velpos(1) )
      leftvel(1) = auxfield( leftngh(1)+velpos(1) )
      rightvel(1) = auxfield( rightngh(1)+velpos(1) )
      gradXXU(1,iElem) = (rightvel(1) - vel(1) + leftvel(1)) &
        &                 * graddata%fdcoeff(elempos, 1)**2
    end do

  end function getGradXXU_1D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This function computes gradient of velocity from gradient and veleocity
  !! data.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used
  pure function getGradXXU_2D( auxField, gradData, velPos, nAuxScalars, &
    &                          nDims, nSolve, elemOffset ) result(gradXXU)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3)
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity:
    ! 1: Dxxu, 2: Dyyv, 3: Dzzw
    real(kind=rk) :: gradXXU(nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(2), rightngh(2)
    real(kind=rk) :: leftvel(2), rightvel(2), vel(2)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
        elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      leftngh(2) = (graddata%neighpos(elempos, 2, 1) - 1) * nAuxScalars
      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      rightngh(2) = (graddata%neighpos(elempos, 2, 2) - 1) * nAuxScalars
      vel(1) = 2._rk * auxField( ElemOffset * nAuxScalars + velpos(1) )
      vel(2) = 2._rk * auxField( ElemOffset * nAuxScalars + velpos(2) )
      leftvel(1) = auxfield( leftngh(1)+velPos(1) )
      leftvel(2) = auxfield( leftngh(2)+velPos(2) )
      rightvel(1) = auxfield( rightngh(1)+velPos(1) )
      rightvel(2) = auxfield( rightngh(2)+velPos(2) )
      gradXXU(1,iElem) = (rightvel(1) - vel(1) + leftvel(1)) &
        &                 * graddata%fdcoeff(elempos, 1)**2
      gradXXU(2,iElem) = (rightvel(2) - vel(2) + leftvel(2)) &
        &                 * graddata%fdcoeff(elempos, 2)**2
    end do

  end function getGradXXU_2D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This function computes gradient of velocity from gradient and veleocity
  !! data.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used
  pure function getGradXXU_3D( auxField, gradData, velPos, nAuxScalars, &
    &                          nDims, nSolve, elemOffset ) result(gradXXU)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3)
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity:
    ! 1: Dxxu, 2: Dyyv, 3: Dzzw
    real(kind=rk) :: gradXXU(nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(3), rightngh(3)
    real(kind=rk) :: leftvel(3), rightvel(3), vel(3)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
      elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      leftngh(2) = (graddata%neighpos(elempos, 2, 1) - 1) * nAuxScalars
      leftngh(3) = (graddata%neighpos(elempos, 3, 1) - 1) * nAuxScalars
      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      rightngh(2) = (graddata%neighpos(elempos, 2, 2) - 1) * nAuxScalars
      rightngh(3) = (graddata%neighpos(elempos, 3, 2) - 1) * nAuxScalars
      vel(1) = 2._rk * auxField( ElemOffset * nAuxScalars + velpos(1) )
      vel(2) = 2._rk * auxField( ElemOffset * nAuxScalars + velpos(2) )
      vel(3) = 2._rk * auxField( ElemOffset * nAuxScalars + velpos(3) )
      leftvel(1) = auxfield( leftngh(1)+velPos(1) )
      leftvel(2) = auxfield( leftngh(2)+velPos(2) )
      leftvel(3) = auxfield( leftngh(3)+velPos(3) )
      rightvel(1) = auxfield( rightngh(1)+velPos(1) )
      rightvel(2) = auxfield( rightngh(2)+velPos(2) )
      rightvel(3) = auxfield( rightngh(3)+velPos(3) )
      gradXXU(1,iElem) = (rightvel(1) - vel(1) + leftvel(1)) &
        &                 * graddata%fdcoeff(elempos, 1)**2
      gradXXU(2,iElem) = (rightvel(2) - vel(2) + leftvel(2)) &
        &                 * graddata%fdcoeff(elempos, 2)**2
      gradXXU(3,iElem) = (rightvel(3) - vel(3) + leftvel(3)) &
        &                 * graddata%fdcoeff(elempos, 3)**2
    end do

  end function getGradXXU_3D
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function computes gradient of rho * velocity^3 from gradient, density
  !! and veleocity data. Just derivatives u_x, v_y and w_z.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used.
  pure function getGradRhoU3_1D( auxField, gradData, velPos, densPos, nAuxScalars, &
    &                            nDims, nSolve, elemOffset ) result(gradRhoU3)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3), densPos
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity
    real(kind=rk) :: gradRhoU3(nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(1), rightngh(1)
    real(kind=rk) :: leftvel(1), rightvel(1)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
      elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      leftvel(1) = auxfield( leftngh(1)+denspos ) * (auxfield( leftngh(1)+velpos(1) ) )** 3
      rightvel(1) = auxfield( rightngh(1)+denspos ) * (auxfield( rightngh(1)+velpos(1) ) )** 3
      gradRhoU3(1,iElem) = (rightvel(1)-leftvel(1)) &
        &                 * graddata%fdcoeff(elempos, 1)
    end do

  end function getGradRhoU3_1D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This function computes gradient of rho * velocity^3 from gradient, density
  !! and veleocity data. Just derivatives u_x, v_y and w_z.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used.
  pure function getGradRhoU3_2D( auxField, gradData, velPos, densPos, nAuxScalars, &
    &                            nDims, nSolve, elemOffset ) result(gradRhoU3)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3), densPos
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity
    real(kind=rk) :: gradRhoU3(nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(2), rightngh(2)
    real(kind=rk) :: leftvel(2), rightvel(2)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
      elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      leftngh(2) = (graddata%neighpos(elempos, 2, 1) - 1) * nAuxScalars
      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      rightngh(2) = (graddata%neighpos(elempos, 2, 2) - 1) * nAuxScalars
      leftvel(1) = auxfield( leftngh(1)+denspos ) * (auxfield( leftngh(1)+velPos(1) ) )** 3
      leftvel(2) = auxfield( leftngh(2)+denspos ) * (auxfield( leftngh(2)+velPos(2) ) )** 3
      rightvel(1) = auxfield( rightngh(1)+denspos ) * (auxfield( rightngh(1)+velPos(1) ) )** 3
      rightvel(2) = auxfield( rightngh(2)+denspos ) * (auxfield( rightngh(2)+velPos(2) ) )** 3
      gradRhoU3(1,iElem) = (rightvel(1)-leftvel(1)) &
        &                 * graddata%fdcoeff(elempos, 1)
      gradRhoU3(2,iElem) = (rightvel(2)-leftvel(2)) &
        &                 * graddata%fdcoeff(elempos, 2)
    end do

  end function getGradRhoU3_2D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This function computes gradient of rho * velocity^3 from gradient, density
  !! and veleocity data. Just derivatives u_x, v_y and w_z.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used.
  pure function getGradRhoU3_3D( auxField, gradData, velPos, densPos, nAuxScalars, &
    &                            nDims, nSolve, elemOffset ) result(gradRhoU3)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3), densPos
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity
    real(kind=rk) :: gradRhoU3(nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(3), rightngh(3)
    real(kind=rk) :: leftvel(3), rightvel(3)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
      elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      leftngh(2) = (graddata%neighpos(elempos, 2, 1) - 1) * nAuxScalars
      leftngh(3) = (graddata%neighpos(elempos, 3, 1) - 1) * nAuxScalars
      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      rightngh(2) = (graddata%neighpos(elempos, 2, 2) - 1) * nAuxScalars
      rightngh(3) = (graddata%neighpos(elempos, 3, 2) - 1) * nAuxScalars
      leftvel(1) = auxfield( leftngh(1)+denspos ) * (auxfield( leftngh(1)+velPos(1) ) )** 3
      leftvel(2) = auxfield( leftngh(2)+denspos ) * (auxfield( leftngh(2)+velPos(2) ) )** 3
      leftvel(3) = auxfield( leftngh(3)+denspos ) * (auxfield( leftngh(3)+velPos(3) ) )** 3
      rightvel(1) = auxfield( rightngh(1)+denspos ) * (auxfield( rightngh(1)+velPos(1) ) )** 3
      rightvel(2) = auxfield( rightngh(2)+denspos ) * (auxfield( rightngh(2)+velPos(2) ) )** 3
      rightvel(3) = auxfield( rightngh(3)+denspos ) * (auxfield( rightngh(3)+velPos(3) ) )** 3
      gradRhoU3(1,iElem) = (rightvel(1)-leftvel(1)) &
        &                 * graddata%fdcoeff(elempos, 1)
      gradRhoU3(2,iElem) = (rightvel(2)-leftvel(2)) &
        &                 * graddata%fdcoeff(elempos, 2)
      gradRhoU3(3,iElem) = (rightvel(3)-leftvel(3)) &
        &                 * graddata%fdcoeff(elempos, 3)
    end do

  end function getGradRhoU3_3D
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This function computes gradient of rho * velocity^3 from gradient, density
  !! and veleocity data. Just derivatives u_x, v_y and w_z.
  !! Gradient is computed using central difference.
  !! if an element has an boundary then neighbor refers to current element
  !! then forward difference is used.
  pure function getGradRhoUVZ_3D( auxField, gradData, velPos, densPos, &
    &                             nDims, nAuxScalars, nSolve, elemOffset )    &
    &  result(gradRhoUVZ)
    ! --------------------------------------------------------------------------
    !> auxField
    real(kind=rk), intent(in) :: auxField(:)
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> Position of velocity field in auxField
    integer, intent(in) :: velPos(3), densPos
    !> Number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> Offset for elements when computing chunkwise
    integer, intent(in) :: elemOffset
    !> Dimensions
    integer, intent(in) :: nDims
    !> output is gradient of velocity
    real(kind=rk) :: gradRhoUVZ(nDims,nSolve)
    ! --------------------------------------------------------------------------
    integer :: iElem, elempos
    integer :: leftngh(3), rightngh(3)
    real(kind=rk) :: leftvel(3), rightvel(3)
    ! --------------------------------------------------------------------------

    do iElem = 1, nSolve
      elempos = ElemOffset + iElem
      leftngh(1) = (graddata%neighpos(elempos, 1, 1) - 1) * nAuxScalars
      leftngh(2) = (graddata%neighpos(elempos, 2, 1) - 1) * nAuxScalars
      leftngh(3) = (graddata%neighpos(elempos, 3, 1) - 1) * nAuxScalars

      rightngh(1) = (graddata%neighpos(elempos, 1, 2) - 1) * nAuxScalars
      rightngh(2) = (graddata%neighpos(elempos, 2, 2) - 1) * nAuxScalars
      rightngh(3) = (graddata%neighpos(elempos, 3, 2) - 1) * nAuxScalars

      leftvel(1) = auxfield( leftngh(1)+denspos ) * auxfield( leftngh(1)+velPos(1) ) &
        &       * auxfield( leftngh(1)+velPos(2) )  * auxfield( leftngh(1)+velPos(3) )
      leftvel(2) = auxfield( leftngh(2)+denspos ) * auxfield( leftngh(2)+velPos(1) ) &
        &       * auxfield( leftngh(2)+velPos(2) )  * auxfield( leftngh(2)+velPos(3) )
      leftvel(3) = auxfield( leftngh(3)+denspos ) * auxfield( leftngh(3)+velPos(1) ) &
        &       * auxfield( leftngh(3)+velPos(2) )  * auxfield( leftngh(3)+velPos(3) )

      rightvel(1) = auxfield( rightngh(1)+denspos ) * auxfield( rightngh(1)+velPos(1) ) &
        &         * auxfield( rightngh(1)+velPos(2) ) * auxfield( rightngh(1)+velPos(3) )
      rightvel(2) = auxfield( rightngh(2)+denspos ) * auxfield( rightngh(2)+velPos(1) ) &
        &         * auxfield( rightngh(2)+velPos(2) ) * auxfield( rightngh(2)+velPos(3) )
      rightvel(3) = auxfield( rightngh(3)+denspos ) * auxfield( rightngh(3)+velPos(1) ) &
        &         * auxfield( rightngh(3)+velPos(2) ) * auxfield( rightngh(3)+velPos(3) )

      gradRhoUVZ(1,iElem) = (rightvel(1)-leftvel(1)) &
        &                 * graddata%fdcoeff(elempos, 1)
      gradRhoUVZ(2,iElem) = (rightvel(2)-leftvel(2)) &
        &                 * graddata%fdcoeff(elempos, 2)
      gradRhoUVZ(3,iElem) = (rightvel(3)-leftvel(3)) &
        &                 * graddata%fdcoeff(elempos, 3)

    end do

  end function getGradRhoUVZ_3D
  ! ************************************************************************** !

end module mus_gradData_module
