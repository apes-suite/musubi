! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
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
! ****************************************************************************** !
!> author: Simon Zimny
!! author: Kannan Masilmani
!! author: Khaled Ibrahim
!! This module provides MUSUBI specific helper functions for extracting
!! information from the solver specific character.
!!
module mus_solSpecHelpers_module

  ! include treelm modules
  use env_module,            only: rk, long_k, PathLen, labelLen
  use tem_variable_module,   only: tem_variable_type
  use tem_stencil_module,    only: tem_stencilHeader_type
  use tem_aux_module,        only: tem_abort
  use tem_time_module,       only: tem_time_type
  use tem_logging_module,    only: logUnit

  ! include aotus modules
  use aotus_module, only: flu_State, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val

  implicit none

  private


  public :: getIdentifyChar
  public :: getNFields
  public :: getFieldPrefixes
  public :: getWeights
  public :: getFieldVariable
  public :: getConversionFac
  public :: getVariable_FromTable

  interface getFieldVariable
    module procedure getFieldVariable_scalar
    module procedure getFieldVariable_array
  end interface getFieldVariable
contains


! ****************************************************************************** !
!             Helper Functions For The Solver Specific Character               !
! ****************************************************************************** !


! ****************************************************************************** !
  !> Get a character from the identify table using a given solver specific
  !! character handle and a given key word.
  !!
  function getIdentifyChar( conf, key )                                        &
    & result( resChar )
    ! ---------------------------------------------------------------------------
    !> handle of the solver specific character
    type( flu_state ), intent(in) :: conf
    !> key to search for
    character(len=*),intent(in) :: key
    !> scheme kind to be returned
    character(len= labelLen) resChar
    ! ---------------------------------------------------------------------------
    ! aotus handle and error variable
    integer :: identifyHandle
    integer :: iError
    ! ---------------------------------------------------------------------------

    call aot_table_open( L = conf, thandle = identifyHandle, key = 'identify' )
    call aot_get_val( val     = resChar,                                       &
      &               ErrCode = iError,                                        &
      &               L       = conf,                                          &
      &               key     = trim(key),                                     &
      &               thandle = identifyHandle,                                &
      &               default = '' )

    call aot_table_close( L = conf, thandle = identifyHandle )

  end function getIdentifyChar
! ****************************************************************************** !

! ****************************************************************************** !
  !> Get the value of variable inside the table name 'key' in scheme
  !! KM: @todo extent it for vectors
  function getVariable_FromTable( conf, varLabel, table )  result( val )
    ! ---------------------------------------------------------------------------
    !> handle of the solver specific character
    type( flu_state ), intent(in) :: conf
    !> complete variable label
    character(len=*), intent(in) :: varLabel
    !> table name
    character(len=*), intent(in) :: table
    !> val to be returned
    real(kind=rk) :: val
    ! ---------------------------------------------------------------------------
    ! aotus handles and error variables
    integer :: iError
    integer :: varHandle
    ! ---------------------------------------------------------------------------
    call aot_table_open( L = conf,                                             &
      &                  thandle = varHandle,                                  &
      &                  key = trim(table) )
    call aot_get_val( val = val,                                               &
      &               ErrCode = iError,                                        &
      &               L = conf,                                                &
      &               key = trim(varLabel),                                    &
      &               thandle = varHandle,                                     &
      &               default = 0._rk )
    ! close the table again
    call aot_table_close( L = conf, thandle = varHandle )

  end function getVariable_FromTable
! ****************************************************************************** !
! ****************************************************************************** !
  !> Get the number of fields from a given solver specific character handle.
  !!
  function getNFields( conf )                                                  &
    & result( nFields )
    ! ---------------------------------------------------------------------------
    !> handle of the solver specific character
    type( flu_state ), intent(in) :: conf
    !> number of fields to be returned
    integer :: nFields
    ! ---------------------------------------------------------------------------
    ! aotus error variable
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! read in nFields
    call aot_get_val( val     = nFields,                                       &
      &               ErrCode = iError,                                        &
      &               L       = conf,                                          &
      &               key     = 'nFields' )

  end function getNFields
! ****************************************************************************** !


! ****************************************************************************** !
  !> Get the right field prefixes from a given solver specific character handle.
  !!
  function getFieldPrefixes( conf, nFields )                                   &
    & result( prefix )
    ! ---------------------------------------------------------------------------
    !> handle of the solver specific character
    type( flu_state ), intent(in) :: conf
    ! total number of fields defined
    integer, intent(in) :: nFields
    !> field prefixes to be returned
    character(len=labelLen) :: prefix(nFields)
    ! ---------------------------------------------------------------------------
    ! counter for field loop
    integer :: iField
    ! aotus handles and error variables
    integer :: prefixHandle
    integer :: iError
    ! ---------------------------------------------------------------------------

    ! open the fieldPrefix table
    call aot_table_open( L       = conf,                                       &
      &                  thandle = prefixHandle,                               &
      &                  key     = 'fieldPrefixes' )
    do iField=1,nFields
      ! read the prefixes from the fieldPrefix table
      call aot_get_val( val     = prefix(iField),                              &
        &               ErrCode = iError,                                      &
        &               L       = conf,                                        &
        &               pos     = iField,                                      &
        &               thandle = prefixHandle,                                &
        &               default = '' )
    end do
    call aot_table_close( L = conf, thandle = prefixHandle )

  end function getFieldPrefixes
! ****************************************************************************** !


! ****************************************************************************** !
  !> Get the field variable name for given field type from a given solver
  !! specific character handle.
  !!
  function getFieldVariable_scalar( conf, varLabel, varName, fieldVar,         &
    &                               fieldProp )  result( val )
    ! ---------------------------------------------------------------------------
    !> handle of the solver specific character
    type( flu_state ), intent(in) :: conf
    !> complete variable label (prefix + pure variable name)
    character(len=*), intent(in) :: varLabel
    !> pure variable name (e.g. density)
    character(len=*), intent(in) :: varName
    !> required name of the field variable
    character(len=*), intent(in) :: fieldVar
    !> Which field type does the field variable belong to.
    !! Example: 'fluid'/'species'
    character(len=*), intent(in) :: fieldProp
    !> val to be returned
    real(kind=rk) :: val
    ! ---------------------------------------------------------------------------
    ! aotus handles and error variables
    integer :: iError
    integer :: fieldHandle, subHandle
    ! starting position of the pure variable name in varLabel
    integer :: pos
    ! prefix of the belonging field
    character(len=labelLen) :: prefix
    ! ---------------------------------------------------------------------------
    ! extract the prefix from the varLabel and open the belonging field table
    pos = INDEX(trim(varLabel), trim(varName))
    prefix = trim(varLabel(1:pos-1))
    call aot_table_open( L       = conf,                                       &
      &                  thandle = fieldHandle,                                &
      &                  key     = trim(prefix)//'field' )
    ! open the species table and extract omega
    call aot_table_open( L       = conf,                                       &
      &                  thandle = subHandle,                                  &
      &                  key     = trim(fieldProp),                            &
      &                  parent  = fieldHandle )
    call aot_get_val( val     = val,                                           &
      &               ErrCode = iError,                                        &
      &               L       = conf,                                          &
      &               key     = trim(fieldVar),                                &
      &               thandle = subHandle,                                     &
      &               default = 0._rk )
    ! close the tables again
    call aot_table_close( L = conf, thandle = subHandle )
    call aot_table_close( L = conf, thandle = fieldHandle )

  end function getFieldVariable_scalar
! ****************************************************************************** !


! ****************************************************************************** !
  !> Get the field variable name for given field type from a given solver
  !! specific character handle.
  !!
  function getFieldVariable_array( conf, varLabel, varName, fieldVar,          &
    &                              fieldProp, nVals ) result( val )
    ! ---------------------------------------------------------------------------
    !> handle of the solver specific character
    type( flu_state ), intent(in) :: conf
    !> complete variable label (prefix + pure variable name)
    character(len=*), intent(in) :: varLabel
    !> pure variable name (e.g. density)
    character(len=*), intent(in) :: varName
    !> required name of the field variable
    character(len=*), intent(in) :: fieldVar
    !> Which field type does the field variable belong to.
    !! Example: 'fluid'/'species'
    character(len=*), intent(in) :: fieldProp
    !> number of entries in the array to read out
    integer, intent(in) :: nVals
    !> val to be returned
    real(kind=rk) :: val( nVals )
    ! ---------------------------------------------------------------------------
    ! aotus handles and error variables
    integer :: iError( nVals )
    integer :: fieldHandle, subHandle
    ! starting position of the pure variable name in varLabel
    integer :: pos
    ! prefix of the belonging field
    character(len=labelLen) :: prefix
    real(kind=rk) :: def0(nVals)
    ! ---------------------------------------------------------------------------
    ! extract the prefix from the varLabel and open the belonging field table
    def0 = 0.0_rk
    pos = INDEX(trim(varLabel), trim(varName))
    prefix = trim(varLabel(1:pos-1))
    call aot_table_open( L       = conf,                                       &
      &                  thandle = fieldHandle,                                &
      &                  key     = trim(prefix)//'field' )
    ! open the species table and extract omega
    call aot_table_open( L       = conf,                                       &
      &                  thandle = subHandle,                                  &
      &                  key     = trim(fieldProp),                            &
      &                  parent  = fieldHandle )
    call aot_get_val( val     = val,                                           &
      &               ErrCode = iError,                                        &
      &               L       = conf,                                          &
      &               default = def0,                                          &
      &               key     = trim(fieldVar),                                &
      &               thandle = subHandle )
    ! close the tables again
    call aot_table_close( L = conf, thandle = subHandle )
    call aot_table_close( L = conf, thandle = fieldHandle )

  end function getFieldVariable_array
! ****************************************************************************** !


! ****************************************************************************** !
  !> Get the the weights of a used stencil from a given solver specific
  !! character handle.
  !!
  function getWeights( conf, stencil )                                         &
    & result( weights )
    ! ---------------------------------------------------------------------------
    !> handle of the solver specific character
    type( flu_state ), intent(in) :: conf
    !> stencil information
    type( tem_stencilHeader_type ), intent(in) :: stencil
    !> weights to be returned
    real(kind=rk) :: weights(stencil%QQ)
    ! ---------------------------------------------------------------------------
    ! aotus error variable
    integer :: errCode(stencil%QQ)
    ! ---------------------------------------------------------------------------

    call aot_get_val( val     = weights,                                       &
      &               ErrCode = errCode,                                       &
      &               L       = conf,                                          &
      &               key     = 'weight' )

  end function getWeights
! ****************************************************************************** !


! ****************************************************************************** !
  !> Get the conversion factor variable from physics table from a given solver
  !! specific character handle.
  !!
  function getConversionFac( conf, facName, nLevels )   &
    & result( val )
    ! ---------------------------------------------------------------------------
    !> handle of the solver specific character
    type( flu_state ), intent(in) :: conf
    !> conversion factor variable label
    character(len=*), intent(in) :: facName
    integer, intent(in) ::  nLevels
    !> val to be returned
    real(kind=rk) :: val(nLevels)
    ! ---------------------------------------------------------------------------
    ! aotus handles and error variables
    integer :: iError
    integer :: phy_Handle, sub_Handle
    !> position of the phyiscs subtable for each level
    integer :: iLevel
    ! ---------------------------------------------------------------------------
    call aot_table_open( L       = conf,                                       &
      &                  thandle = phy_Handle,                                 &
      &                  key     = 'physics' )
    do iLevel = 1, nLevels
      call aot_table_open( L       = conf,                                     &
        &                  thandle = sub_Handle,                               &
        &                  parent  = phy_Handle,                               &
        &                  pos     = iLevel )
      call aot_get_val( val     = val( iLevel ),                               &
        &               ErrCode = iError,                                      &
        &               L       = conf,                                        &
        &               key     = trim(facName),                               &
        &               thandle = sub_Handle  )
      ! close the tables again
      call aot_table_close( L = conf, thandle = sub_Handle )
    end do
    call aot_table_close( L = conf, thandle = phy_Handle )

    if ( minval( val ) < 0.0_rk ) then
      write(logUnit(1),*)'Requested conversion factor variable name: '//       &
        &            trim(facName)
      write(logUnit(1),*)  'by tracking is not well defined'
      call tem_abort()
    end if
  end function getConversionFac
! ****************************************************************************** !


end module mus_solSpecHelpers_module
! ****************************************************************************** !
