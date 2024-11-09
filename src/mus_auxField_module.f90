! Copyright (c) 2019-2020, 2022 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
module mus_auxField_module
  use env_module,         only: rk
  use tem_param_module,   only: rho0, rho0Inv
  use tem_varSys_module,  only: tem_varSys_type, tem_varSys_op_type
  use treelmesh_module,   only: treelmesh_type
  use tem_stencil_module, only: tem_stencilHeader_type
  use tem_comm_module,    only: tem_commpattern_type, tem_communication_type, &
    &                           tem_comm_init
  use tem_construction_module, only: tem_levelDesc_type
  use tem_general_module,      only: tem_general_type

  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_scheme_header_module,      only: mus_scheme_header_type
  use mus_interpolate_header_module, only: mus_interpolation_type
  use mus_pdf_module,                only: pdf_data_type
  use mus_source_type_module,        only: mus_source_type
  use mus_field_module,              only: mus_field_type
  use mus_physics_module,            only: mus_convertFac_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none
  private

  public :: mus_init_auxFieldArrays
  public :: mus_auxFieldVar_type
  public :: mus_proc_calcAuxField

  public :: mus_initAuxFieldFluidAndExchange
  public :: mus_calcAuxFieldAndExchange
  public :: mus_intpAuxFieldCoarserAndExchange
  public :: mus_intpAuxFieldFinerAndExchange

  !> Contains auxiliary field variable values per level and communication
  !! buffers
  type mus_auxFieldVar_type
    !> auxiliary field variable values computed from pre-collision PDF
    !! after PDF exchange
    !! Size: nSize*nScalars
    !! Element order is same as state array
    !! Access: (iElem-1)*nScalars + varSys%method%val(iVar)%auxField_varPos
    !! See mus_append_auxField for the name of the variable stored in this
    !! array as it depends on the scheme kind.
    real(kind=rk), allocatable :: val(:)

    !> Local Fluids required by remote processes
    type( tem_communication_type ) :: sendBuffer
    !> Local ghostFromCoarser required by remote processes
    type( tem_communication_type ) :: sendBufferFromCoarser
    !> Local ghostFromFiner required by remote processes
    type( tem_communication_type ) :: sendBufferFromFiner
    !> My halos which are fluids on remote processes
    type( tem_communication_type ) :: recvBuffer
    !> My halos which are ghostFromCoarser on remote processes
    type( tem_communication_type ) :: recvBufferFromCoarser
    !> My halos which are ghostFromFiner on remote processes
    type( tem_communication_type ) :: recvBufferFromFiner
  end type mus_auxFieldVar_type

  abstract interface
    !> Interface to compute auxField vars i.e. conserved macroscopic moments
    !! from pre-collision PDF for fluid and ghostFromCoarser.
    !! auxField on GhostFromFiner elements are interpolated and
    !! halo elements are exchanged
    !! For Multicomponent models: in calcAuxField function, the velocity
    !! is computed on transformed PDF such that force term can be added to it
    !! in addSrcToAuxField routine. The auxField is updated with correct
    !! velocity field in compute kernel
    !! i.e. velocity of original PDF is obtained by solving
    !! linear equation system  in compute kernel
    subroutine mus_proc_calcAuxField(auxField, state, neigh, nSize, nSolve, &
      & iLevel, stencil, varSys, derVarPos, quantities)
      import :: rk, tem_varSys_type, tem_stencilHeader_type, mus_derVarPos_type, &
        &       mus_scheme_derived_quantities_type
      !> output auxField array
      real(kind=rk), intent(inout) :: auxField(:)
      !> input state array
      real(kind=rk), intent(in) :: state(:)
      !> connectivity array
      integer, intent(in) :: neigh(:)
      !> number of elements in the state array
      integer, intent(in) :: nSize
      !> number of fluid elements + ghostFromCoarser
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
    end subroutine mus_proc_calcAuxField

  end interface



contains

  ! ************************************************************************* !
  !> This routine initialize auxField var val array and communication buffers
  subroutine mus_init_auxFieldArrays(me, levelDesc, pattern, nSize, nAuxScalars)
    ! --------------------------------------------------------------------- !
    !> Auxiliary field variable
    type(mus_auxFieldVar_type), intent(out) :: me
    !> levelDesc to access communication buffers of state array
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> communication pattern
    type(tem_commPattern_type), intent(in) :: pattern
    !> Number of elements in state array
    integer, intent(in) :: nSize
    !> Number of scalars in auxiliary variables
    integer, intent(in) :: nAuxScalars
    ! --------------------------------------------------------------------- !
    ! --------------------------------------------------------------------- !
    allocate(me%val(nSize * nAuxScalars))
    me%val(:) = -1000000.0_rk

    ! initialize send buffer
    call init_commBuffer(buffer_aux   = me%sendBuffer,       &
      &                  buffer_State = levelDesc%sendBuffer )

    ! initialize recv buffer
    call init_commBuffer(buffer_aux   = me%recvBuffer,       &
      &                  buffer_State = levelDesc%recvBuffer )

     ! initialize send buffer
    call init_commBuffer(buffer_aux   = me%sendBufferFromCoarser,       &
      &                  buffer_State = levelDesc%sendBufferFromCoarser )

    ! initialize recv buffer
    call init_commBuffer(buffer_aux   = me%recvBufferFromCoarser,       &
      &                  buffer_State = levelDesc%recvBufferFromCoarser )

      ! initialize send buffer
    call init_commBuffer(buffer_aux   = me%sendBufferFromFiner,       &
      &                  buffer_State = levelDesc%sendBufferFromFiner )

    ! initialize recv buffer
    call init_commBuffer(buffer_aux   = me%recvBufferFromFiner,       &
      &                  buffer_State = levelDesc%recvBufferFromFiner )

  contains

    ! ************************************************************************ !
    subroutine init_commBuffer(buffer_aux, buffer_state)
      !-------------------------------------------------------------------------
      !> communication buffer for velocity field
      type(tem_communication_type), intent(out) :: buffer_aux
      !> communication buffer of state array which is already initialized
      !! in tem_construction_module
      type(tem_communication_type), intent(in) :: buffer_state
      !-------------------------------------------------------------------------
      ! positions, in velocity vectors
      integer, allocatable :: pos(:)
      integer :: iProc, iElem, iScal, counter, elemPos
      !-------------------------------------------------------------------------
      ! copy information about target and source procs from pdf sendBuffer to
      ! velocity sendBuffer
      call tem_comm_init(buffer_aux, buffer_state%nProcs)
      buffer_aux%proc = buffer_state%proc
      buffer_aux%nElemsProc = buffer_state%nElemsProc

      allocate(pos(maxval(buffer_aux%nElemsProc*nAuxScalars)))

      do iProc = 1, buffer_aux%nProcs
        counter = 0
        do iElem = 1, buffer_aux%nElemsProc( iProc )
          ! halo element position in the levelDesc total list
          elemPos = buffer_state%elemPos( iProc )%val( iElem )
          do iScal = 1, nAuxScalars
            counter = counter + 1
            pos(counter) = (elemPos-1)*nAuxScalars + iScal
          end do
        end do
        ! copy position array to me%pos, allocate me%val array
        call pattern%initBuf_real( buffer_aux%buf_real( iProc ), pos, counter )
      end do

      deallocate(pos)
    end subroutine init_commBuffer
    ! ************************************************************************ !

  end subroutine mus_init_auxFieldArrays
  ! ************************************************************************* !


  ! ************************************************************************** !
  !> This routine initializes auxField for fluid elements using SAVE access on
  !! PDF initialized by IC
  subroutine mus_initAuxFieldFluidAndExchange(auxField, state, neigh, nElems,  &
    &                                         nSize, nFields, stencil, varSys, &
    &                                         derVarPos, iLevel, general, quantities)
    !---------------------------------------------------------------------------
    !> auxilary field array
    type(mus_auxFieldVar_type), intent(inout) :: auxField
    !> state array
    real(kind=rk), intent(in) :: state(:)
    !> connectivity vector
    integer, intent(in) :: neigh(:)
    !> number of elements to compute auxField
    integer, intent(in) :: nElems
    !> number of elements in state array
    integer, intent(in) :: nSize
    !> number of fields
    integer, intent(in) :: nFields
    !> current level
    integer, intent(in) :: iLevel
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> contains auxField position of all fields in varSys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> contains commPattern, MPI communicator and simControl
    type(tem_general_type), intent(in) :: general
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    !---------------------------------------------------------------------------
    integer :: iField
    !---------------------------------------------------------------------------
    ! calculate auxField from local state
    do iField = 1, nFields
      call derVarPos(iField)%auxFieldFromState( state    = state,        &
        &                                       neigh    = neigh,        &
        &                                       iField   = iField,       &
        &                                       nElems   = nElems,       &
        &                                       nSize    = nSize,        &
        &                                       iLevel   = iLevel,       &
        &                                       stencil  = stencil,      &
        &                                       varSys   = varSys,       &
        &                                       auxField = auxField%val, &
        &                                       quantities = quantities  )
    end do

    ! communicate velocity field. Requires for tubulence to compute ShearRate
    ! from velocity gradient.
    ! exchange velocity halo on current level
    call general%commpattern%exchange_real(   &
      &  send         = auxField%sendBuffer,  &
      &  recv         = auxField%recvBuffer, &
      &  state        = auxField%val(:),      &
      &  message_flag = iLevel+100,           &
      &  comm         = general%proc%comm     )

  end subroutine mus_initAuxFieldFluidAndExchange
  ! ************************************************************************** !


  ! ************************************************************************* !
  !> This routine compute auxField variable from pre-collision pdf and exchange
  !! halos
  subroutine mus_calcAuxFieldAndExchange(auxField, calcAuxField, state, &
    &  pdfData, nFields, field, globSrc, stencil, varSys, derVarPos,    &
    &  phyConvFac, general, iLevel, minLevel, schemeHeader, quantities)
    ! -------------------------------------------------------------------- !
    !> auxilary field array
    type(mus_auxFieldVar_type), intent(inout) :: auxField
    !> function pointer to calculate auxField
    procedure(mus_proc_calcAuxField), pointer, intent(in) :: calcAuxField
    !> state array
    real(kind=rk), intent(in) :: state(:)
    !> contains neigh array and nElems on current level
    type(pdf_data_type), intent(in) :: pdfData
    !> Number of fields
    integer, intent(in) :: nFields
    !> contains sources of all fields
    type(mus_field_type), intent(inout)  :: field(nFields)
    !> global source
    type(mus_source_type), intent(inout)  :: globSrc
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> contains auxField position of all fields in varSys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    !> physics conversion factors for this level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> contains commPattern, MPI communicator and simControl
    type(tem_general_type), intent(in) :: general
    !> current level
    integer, intent(in) :: iLevel
    !> minlevel
    integer, intent(in) :: minLevel
    !> scheme header
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    !> Class that contains pointers to the proper derived quantities functions
    type(mus_scheme_derived_quantities_type), intent(in) :: quantities
    ! -------------------------------------------------------------------- !
    integer :: nSolve, iField, iSrc
    ! -------------------------------------------------------------------- !

    ! calculate auxField only for fluids and ghostfromcoarser (buffer ghost).
    ! ghostFromFiner are interpolated
    nSolve = pdfData%nElems_fluid + pdfData%nElems_ghostFromCoarser
    ! calculate velocity for fluid and ghost elements
    call calcAuxField( auxField   = auxField%val(:), &
      &                state      = state,           &
      &                neigh      = pdfData%neigh,   &
      &                nSize      = pdfData%nSize,   &
      &                nSolve     = nSolve,          &
      &                iLevel     = iLevel,          &
      &                stencil    = stencil,         &
      &                varSys     = varSys,          &
      &                derVarPos  = derVarPos,       &
      &                quantities = quantities       )

    ! update auxField with source term
    ! Field add field source and then add globSrc
    do iField = 1, nFields
      do iSrc = 1, field(iField)%source%varDict%nVals
        call field(iField)%source%method(iSrc)%addSrcToAuxField( &
          & auxField   = auxField%val(:),                        &
          & iLevel     = iLevel,                                 &
          & time       = general%simControl%now,                 &
          & varSys     = varSys,                                 &
          & phyConvFac = phyConvFac,                             &
          & derVarPos  = derVarPos                               )
      end do

      ! Add internal source to auxField. Internal source is added to state
      ! in compute routine
      do iSrc = 1, field(iField)%internalSource%varDict%nVals
        call field(iField)%internalSource%method(iSrc)%addSrcToAuxField( &
          & auxField   = auxField%val(:),                                &
          & iLevel     = iLevel,                                         &
          & time       = general%simControl%now,                         &
          & varSys     = varSys,                                         &
          & phyConvFac = phyConvFac,                                     &
          & derVarPos  = derVarPos                                       )
      end do
    end do

    ! apply global source
    do iSrc = 1, globSrc%varDict%nVals
      call globSrc%method(iSrc)%addSrcToAuxField( &
        & auxField   = auxField%val(:),           &
        & iLevel     = iLevel,                    &
        & time       = general%simControl%now,    &
        & varSys     = varSys,                    &
        & phyConvFac = phyConvFac,                &
        & derVarPos  = derVarPos                  )
    end do

    ! communicate velocity field. Requires for tubulence to compute ShearRate
    ! from velocity gradient.
    ! exchange velocity halo on current level
    call general%commpattern%exchange_real(   &
      &  send         = auxField%sendBuffer,  &
      &  recv         = auxField%recvBuffer , &
      &  state        = auxField%val(:),      &
      &  message_flag = iLevel+100,           &
      &  comm         = general%proc%comm     )

    ! communicate ghost halos from coarser
    if (iLevel > minLevel) then
       call general%commpattern%exchange_real(             &
        &  send         = auxField%sendBufferFromCoarser,  &
        &  recv         = auxField%recvBufferFromCoarser , &
        &  state        = auxField%val(:),                 &
        &  message_flag = iLevel+200,                      &
        &  comm         = general%proc%comm                )
    end if

  end subroutine mus_calcAuxFieldAndExchange
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> This routine interpolate auxField variable for ghostFromFiner and exchange
  !! halos
  subroutine mus_intpAuxFieldCoarserAndExchange(intp, tAuxField, sAuxField,  &
    &                                           tLevelDesc, stencil, iLevel, &
    &                                           nAuxScalars,  general)
    ! -------------------------------------------------------------------- !
    !> Interpolation type
    type(mus_interpolation_type), intent(inout) :: intp
    !> target auxilary field array
    type(mus_auxFieldVar_type), intent(inout) :: tAuxField
    !> source auxilary field array
    type(mus_auxFieldVar_type), intent(in) :: sAuxField
    !> level descriptor on target level
    type(tem_levelDesc_type), intent(in) :: tLevelDesc
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> current level
    integer, intent(in) :: iLevel
    !> number of scalars to interpolate
    integer, intent(in) :: nAuxScalars
    !> contains commPattern, MPI communicator and simControl
    type(tem_general_type), intent(in) :: general
    ! -------------------------------------------------------------------- !
    call intp%fillMineFromFiner%do_intpArbiVal(      &
      & tLevelDesc = tLevelDesc,                     &
      & level      = iLevel,                         &
      & stencil    = stencil,                        &
      & sVal       = sAuxField%val(:),               &
      & tVal       = tAuxField%val(:),               &
      & nTargets   = tLevelDesc%intpFromFiner%nVals, &
      & targetList = tLevelDesc%intpFromFiner%val,   &
      & nScalars   = nAuxScalars                     )

    ! exchange velocity halo fromFiner, required to compute velocity
    ! gradient
    call general%commPattern%exchange_real(            &
      &  send         = tAuxField%sendBufferFromFiner, &
      &  recv         = tAuxField%recvBufferFromFiner, &
      &  state        = tAuxField%val(:),              &
      &  message_flag = iLevel+300,                    &
      &  comm         = general%proc%comm              )

  end subroutine mus_intpAuxFieldCoarserAndExchange
  ! ************************************************************************* !

  ! ************************************************************************* !
  !> This routine interpolate auxField variable for ghostFromCoarser and exchange
  !! halos
  subroutine mus_intpAuxFieldFinerAndExchange(intp, tAuxField, sAuxField,  &
    &                                         tLevelDesc, stencil, iLevel, &
    &                                         nAuxScalars, general)
    ! -------------------------------------------------------------------- !
    !> Interpolation type
    type(mus_interpolation_type), intent(inout) :: intp
    !> target auxilary field array
    type(mus_auxFieldVar_type), intent(inout) :: tAuxField
    !> source auxilary field array
    type(mus_auxFieldVar_type), intent(in) :: sAuxField
    !> level descriptor on target level
    type(tem_levelDesc_type), intent(in) :: tLevelDesc
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> current level
    integer, intent(in) :: iLevel
    !> number of scalars to interpolate
    integer, intent(in) :: nAuxScalars
    !> contains commPattern, MPI communicator and simControl
    type(tem_general_type), intent(in) :: general
    ! -------------------------------------------------------------------- !
    integer :: iOrder
    ! -------------------------------------------------------------------- !
    do iOrder = 0, intp%config%order
      call intp%fillFinerFromMe(iOrder)%do_intpArbiVal(          &
        & tLevelDesc = tLevelDesc,                               &
        & level      = iLevel,                                   &
        & stencil    = stencil,                                  &
        & sVal       = sAuxField%val(:),                         &
        & tVal       = tAuxField%val(:),                         &
        & nTargets   = tLevelDesc%intpFromCoarser(iOrder)%nVals, &
        & targetList = tLevelDesc%intpFromCoarser(iOrder)%val,   &
        & nScalars   = nAuxScalars                               )
     end do

    ! exchange velocity halo fromFiner, required to compute velocity
    ! gradient
    call general%commPattern%exchange_real(              &
      &  send         = tAuxField%sendBufferFromCoarser, &
      &  recv         = tAuxField%recvBufferFromCoarser, &
      &  state        = tAuxField%val(:),                &
      &  message_flag = iLevel+200,                      &
      &  comm         = general%proc%comm                )

  end subroutine mus_intpAuxFieldFinerAndExchange
  ! ************************************************************************* !
end module mus_auxField_module

