! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
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
! *****************************************************************************!
!> author: Kannan Masilamani
!! Module containing subroutines for building MUSUBI specific transport
!! variables to use in compute kernels and source update
!!
module mus_transport_var_module
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer

  ! include treelm modules
  use env_module,               only: rk, long_k, labelLen, solSpecLen, newUnit
  use tem_aux_module,           only: tem_abort
  use tem_varSys_module,        only: tem_varSys_type,               &
    &                                 tem_varSys_append_derVar,      &
    &                                 tem_varSys_proc_point,         &
    &                                 tem_varSys_proc_element,       &
    &                                 tem_varSys_proc_setParams,     &
    &                                 tem_varSys_proc_getParams,     &
    &                                 tem_varSys_proc_setupIndices,  &
    &                                 tem_varSys_proc_getValOfIndex, &
    &                                 tem_varSys_dump
  use tem_varMap_module,        only: tem_possible_variable_type, &
    &                                 init, append, truncate,     &
    &                                 tem_variable_loadMapping
  use tem_stringKeyValuePair_module, only: init, truncate, &
    &                                      grw_stringKeyValuePairArray_type
  use tem_stencil_module,       only: tem_stencilHeader_type
  use treelmesh_module,         only: treelmesh_type
  use tem_geometry_module,      only: tem_BaryOfID
  use tem_logging_module,       only: logUnit
  use tem_operation_module,     only: tem_indexLvl_type
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_dyn_array_module,     only: PositionOfVal

  ! include musubi modules
  use mus_scheme_header_module, only: mus_scheme_header_type

  ! include aotus modules
  use aotus_module,   only: flu_State
  use aot_out_module, only: aot_out_type, aot_out_open, aot_out_close, &
    &                       aot_out_val, aot_out_toChunk

  implicit none
  private

  public :: mus_transport_var_type
  public :: mus_create_poss_transVar
  public :: mus_load_transport_var
  public :: mus_init_transport_var

  ! ***************************************************************************!
  !> Description contains index to access value using variable function
  !! pointer
  type mus_transport_op_type
    !> Position of data variable provided in config file in the varSys
    integer :: data_varpos

    !> Indices for points for nElems_solve (nFluids + nGhostsFromCoarser).
    !! Order of index matches levelDesc%total list, required for setup_index,
    !! getvalof_Index
    type(tem_indexLvl_type) :: pntIndex
  end type mus_transport_op_type
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> Description of musubi transport variable type
  type mus_transport_var_type
    !> Contains variable pntIndex to setup_index and getValOfIndex
    !! Size: varDict%nVals
    type(mus_transport_op_type), allocatable :: method(:)

    !> Dictionary of transport variable with
    !! varDict%val()%key is the name of transport variable and
    !! varDict%val()%value is the name of variable provided for the key
    type(grw_stringKeyValuePairArray_type) :: varDict
  end type mus_transport_var_type
  ! ***************************************************************************!

contains


  ! ***************************************************************************!
  !> Routine initialize possible transport variable depends on scheme kind
  subroutine mus_create_poss_transVar(poss_transVar, schemeHeader)
    ! --------------------------------------------------------------------------!
    !> possible transport variables
    type(tem_possible_variable_type), intent(out) :: poss_transVar

    !> Identifier of the scheme
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------!
    write(logUnit(10),*) 'Creating possible transport variables '
    call init(me = poss_transVar, length = 2 )

    select case(trim(schemeHeader%kind))
    case ('passive_scalar', 'nernst_planck')
      call append(me          = poss_transVar,        &
        &         varName     = 'transport_velocity', &
        &         nComponents = 3                     )
    case default
      write(logUnit(1),*) 'No possible transport variable defined for ' &
        &               //'scheme kind: '//trim(schemeHeader%kind)
    end select
    call truncate(poss_transVar)

  end subroutine mus_create_poss_transVar
  ! ***************************************************************************!


  ! ***************************************************************************!
  !> Routine load musubi transport variables
  subroutine mus_load_transport_var(me, possVars, conf, parent, varSys, &
    &                               schemeHeader)
    ! --------------------------------------------------------------------------!
    !> transport variable type to initialize
    type(mus_transport_var_type), intent(out) :: me
    !> possible transport variables
    type(tem_possible_variable_type), intent(in) :: possVars
    !> flu state
    type( flu_State ) :: conf
    !> parent handle if scheme table is defined
    integer, intent(in), optional :: parent
    !> Global variable system
    type(tem_varSys_type), intent(inout) :: varSys
    !> Identifier of the scheme
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------!
    integer :: iVar
    ! --------------------------------------------------------------------------!
    write(logUnit(1),*) 'Loading transport variables'
    ! initialize growing array stringKeyValuePair
    call init( me = me%varDict )

    ! load the transport variables
    do iVar = 1, possVars%varName%nVals
      call tem_variable_loadMapping(                    &
        & expectedName = possVars%varName%val(iVar),    &
        & conf         = conf,                          &
        & thandle      = parent,                        &
        & varDict      = me%varDict,                    &
        & varSys       = varSys,                        &
        & nComp        = possVars%nComponents%val(iVar) )
    end do

    select case(trim(schemeHeader%kind))
    case ('passive_scalar','nernst_planck')
      if (me%varDict%nVals /= 1) then
        write(logUnit(1),*) 'Error: transport_velocity'         &
          &              // ' variable is not defined for lbm_ps'
        call tem_abort()
      end if
    end select

    call truncate( me = me%varDict )

  end subroutine mus_load_transport_var
  ! ***************************************************************************!

  ! ***************************************************************************!
  !> Initialize transport variable by calling setupIndices for every variable
  !! and store pntIndex
  subroutine mus_init_transport_var(me, varSys, tree, nElems_solve, levelDesc)
    ! --------------------------------------------------------------------------
    !> transport variable to fill in
    type(mus_transport_var_type), intent(inout)   :: me

    !> global variable system
    type(tem_varSys_type), intent(in)      :: varSys

    !> global treelm mesh
    type( treelmesh_type ), intent(in)     :: tree

    !> Number of elements to solve in all levels
    !! nFluids + nGhosts
    integer, intent(in) :: nElems_solve(tree%global%minLevel:)

    !> Level descriptors
    type( tem_levelDesc_type ), intent(in) :: levelDesc(tree%global%minLevel:)
    ! --------------------------------------------------------------------------
    integer :: iLevel, iElem, iVar
    integer :: nSolve, minLevel, maxLevel
    real(kind=rk), allocatable :: bary(:,:)
    integer, allocatable :: idx(:)
    integer :: data_varPos
    ! --------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*)'Initializing transport variables ...'

    ! immediately exit this routine if there are no transport variables
    if ( me%varDict%nVals == 0) then
      write(logUnit(1),*) 'No active transport variables'
      return
    else
      allocate(me%method(me%varDict%nVals))
    end if

    do iVar = 1, me%varDict%nVals
      data_varPos = PositionOfVal( me  = varSys%varName,                  &
        &                          val = trim(me%varDict%val(iVar)%value) )
      if (data_varPos > 0) then
        me%method(iVar)%data_varPos = data_varPos
      else
        write(logUnit(1),*) 'Error: variable '                 &
          &               // trim(me%varDict%val(iVar)%value)  &
          &              // ' is not added to variable system'
        call tem_abort()
      end if
    end do

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel

    write(logunit(10),*) ' setup indices for transport var'
    do ilevel = minlevel, maxlevel
      write(logunit(10),*) 'ilevel: ', ilevel

      nsolve = nelems_solve(ilevel)
      ! gets barycenter of all elements to solve i.e fluid+ghost to
      ! access transport variable inside compute kernel and source update
      allocate(bary(nSolve, 3))
      do iElem = 1, nSolve
        bary(iElem, :) = tem_BaryOfId(tree, levelDesc(iLevel)%total(iElem))
      end do

      allocate(idx(nSolve))

      do iVar = 1, me%varDict%nVals
        idx = 0
        data_varPos = me%method(iVar)%data_varPos
        ! set params
        call varSys%method%val(data_varPos)%set_params( &
          & varSys   = varSys,                          &
          & instring = 'isSurface = false'              )

        call varSys%method%val(data_varPos)%setup_indices( &
          & varSys     = varSys,                           &
          & point      = bary,                             &
          & iLevel     = iLevel,                           &
          & tree       = tree,                             &
          & nPnts      = nSolve,                           &
          & idx        = idx                               )

        call append(me%method(iVar)%pntIndex%indexLvl(iLevel), idx)

      end do !iVar
      deallocate(idx)
      deallocate(bary)
    end do !iLevel

  end subroutine mus_init_transport_var
  ! ***************************************************************************!

end module mus_transport_var_module
! *****************************************************************************!
