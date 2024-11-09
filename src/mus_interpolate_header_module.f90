! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2011-2012 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2015, 2017-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017-2018, 2020 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
! ****************************************************************************** !
!> author: Manuel Hasert
!! author: Jiaxing Qi
!! Interpolation header to load confiugration and type definition
module mus_interpolate_header_module

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_val

  ! include treelm modules
  use env_module,            only: rk, labelLen
  use tem_tools_module,      only: tem_horizontalSpacer
  use tem_aux_module,        only: tem_abort
  use tem_logging_module,    only: logUnit
  use tem_debug_module,      only: dbgUnit
  use tem_varSys_module,     only: tem_varSys_type
  use tem_grow_array_module, only: grw_realArray_type, init
  use tem_construction_module, only: tem_levelDesc_type
  use tem_param_module,      only: qN00, q0N0, q00N, q100, q010, q001, &
    &                              q0NN, q0N1, q01N, q011, qN0N, q10N, &
    &                              qN01, q101, qNN0, qN10, q1N0, q110, &
    &                              qNNN, qNN1, qN1N, qN11, q1NN, q1N1, &
    &                              q11N, q111, q000
  use tem_matrix_module,     only: tem_intpMatrixLSF_type, init
  use tem_time_module,       only: tem_time_type
  use tem_stencil_module,    only: tem_stencilHeader_type

  ! include musubi modules
  use mus_field_prop_module,     only: mus_field_prop_type
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_physics_module,        only: mus_physics_type
  use mus_derVarPos_module,      only: mus_derVarPos_type
  use mus_directions_module,     only: qN0, q0N, q10, q01, qNN, qN1, q1N, q11

  implicit none

  private

  public :: mus_interpolation_type
  public :: mus_interpolation_config_type
  public :: mus_interpolation_method_type
  public :: mus_interpolation_stencil_type
  public :: mus_load_interpolate
  public :: mus_interpolate_out
  public :: mus_set_nSources

  !> Interpolation parameter to choose fillFinerFromMe
  integer, parameter, public ::          no_intp = -1
  integer, parameter, public :: weighted_average = 0
  integer, parameter, public ::           linear = 1
  integer, parameter, public ::        quadratic = 2

  ! D2Q9 directions are imported from mus_directions_module, we only need to
  ! specify q00
  integer,parameter :: q00 = 9  !< rest

  !> This data types contains intpRoutine function pointer for FillFiner
  !! and FillCoarser.
  !! For fillFiner, it build least square fit matrix for linear
  !! quadratic interpolations
  !! For fillCoarser: currently we do simple average
  !!
  !! Why do we need different intpRoutine for fillFinerFromMe?
  !! The order of interpolation to finer depends on available number
  !! of coarser source elements so for every order we use different
  !! interpolation routines.
  !! We start with user defined interpolation order and
  !! If nMinSources for that order is not found then we fall back to lower order
  !! Weighted average will be the lowest level for which nMinSources = 1
  type mus_interpolation_method_type

    !> Routine to interpolate coarse to fine for ghostFromCoarser elements
    !! and interpolate fine to coarse for ghostFromFiner elements.
    !! Sets pdf for ghost elements by f_eq + f_neq
    !! The moments required to compute equilibrium function are obtained
    !! from auxField array and the auxField of ghost elements are interpolate
    !! seperately using do_intpArbitraryField
    procedure(intpRoutine), pointer :: do_intp => null()

    !> Routine to interpolate coarse to fine and fine to coarse for
    !! arbitrary variables
    procedure(intpRoutine_arbitraryVal), pointer :: do_intpArbiVal => null()


    !> Matrix entries for linear/Quadratic interpolation least square fit
    !! ((A^T)A)^-1*(A^T)
    !! Size: (6,9) for D2Q9 stencil
    !! Size: (10,QQ)  for D3Q19 and D3Q27
    type(tem_intpMatrixLSF_type) :: intpMat_forLSF

    !> how many source elements are required by this interpolation order
    integer :: nMinSources

    !> Max number of sources amoung target ghosts
    !! Computed in mus_contruction::mus_intp_complete_coarseDep
    integer :: nMaxSources
  end type mus_interpolation_method_type


  !> Contains stencil for interpolation
  type mus_interpolation_stencil_type
    !> Is active only for specific layouts like d2q9, d3q19, d3q27
    logical :: isActive = .false.
    !> cxDir for interpolation stencil for depFromCoarser
    integer, allocatable :: neighDir(:,:)
  end type mus_interpolation_stencil_type

  !> Contains information loaded from config file
  type mus_interpolation_config_type
    !> name of the order of the interpolation method for fillFinerFromMe
    character(len=labelLen) :: method

    !> Order of the interpolation for fillFinerFromMe
    integer :: order

    !> name of used weighting method
    character(len=labelLen) :: weights_method = 'linear_distance'

    !> Power factor for inverse distance weighting
    integer :: IDW_powerfac = 6

    !> Stencil for linear interpolation.
    !! By default use stencil from weighted average
    logical :: useComputeStencil = .false.

    !> Interpolation test by comparing against the initial condition
    logical :: testInterpolation = .false.
    logical :: testEachElement   = .false.
    logical :: testFluids        = .false.
    logical :: noIntpFromFiner   = .false.
    logical :: noIntpFromCoarser = .false.

  end type mus_interpolation_config_type

  !> definition of the used interpolation method
  type mus_interpolation_type

    !> Information loaded from config file
    type(mus_interpolation_config_type) :: config

    !> Interpolation routines to fillFiner
    !! Size: interpolation order
    type(mus_interpolation_method_type), allocatable :: fillFinerFromME(:)

    !> Interpolation routines to fillCoarser
    type(mus_interpolation_method_type) :: fillMineFromFiner

    !> stencil for weighted average interpolation
    type(mus_interpolation_stencil_type) :: weightedAvgStencil

  end type mus_interpolation_type

! ****************************************************************************** !
  abstract interface
    !> This is the interface for all interpolation methods that
    !  can be called from outside to set state variable for ghost elements
    subroutine intpRoutine( method,  fieldProp, tLevelDesc, level, sState, &
      & sNeigh, snSize, sAuxfield, tState, tNeigh, tnSize,                 &
      & layout, nTargets, targetList, physics, time, varSys, derVarPos     )
      import :: mus_interpolation_method_type, mus_scheme_layout_type,         &
        &       mus_field_prop_type, tem_levelDesc_type, mus_physics_type, rk, &
        &       tem_varSys_type, tem_time_type, mus_derVarPos_type

      class(mus_interpolation_method_type), intent(inout) :: method

      !> Array of field properties (fluid or species)
      type(mus_field_prop_type), target, intent(in) :: fieldProp(:)

      !> my refinement level
      integer,intent(in) :: level

      !> State vector of SOURCE elements
      real(kind=rk),    intent(in) :: sState(:)
      integer,          intent(in) :: sNeigh(:)
      integer,          intent(in) :: snSize

      !> AuxField variable to read rho and vel from source elements
      real(kind=rk), intent(inout) :: sAuxField(:)

      !> State vector of TARGET GHOST elements
      real(kind=rk), intent(inout) :: tState(:)
      integer,          intent(in) :: tNeigh(:)
      integer,          intent(in) :: tnSize

      !> level descriptor on target level
      type( tem_levelDesc_type ), intent(in) :: tLevelDesc

      !> the layout used
      type( mus_scheme_layout_type ), intent(in) :: layout

      !> List of target elements ( their position in depSource list )
      integer, intent(in) :: nTargets
      integer, intent(in) :: targetList(nTargets)

      !> physics type to convert lattice to physics SI unit and vice versa
      type(mus_physics_type), intent(in) :: physics

      !> time required to compute analytical solution for TGV case
      type(tem_time_type), intent(in) :: time

      !> scheme variable system
      type( tem_varSys_type ), intent(in) :: varSys

      !> position of all derive variable in varSys
      type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    end subroutine intpRoutine

    !> This is the interface for all interpolation methods that
    !  can be called from outside to interpolate arbitrary variable
    subroutine intpRoutine_arbitraryVal( method, tLevelDesc, level, stencil, &
      &                                  sVal, tVal, nTargets, targetList,   &
      &                                  nScalars )
      import :: mus_interpolation_method_type, rk, tem_levelDesc_type, &
        &       tem_stencilHeader_type

      class(mus_interpolation_method_type), intent(inout) :: method

      !> level descriptor on target level
      type( tem_levelDesc_type ), intent(in) :: tLevelDesc

      !> my refinement level
      integer,intent(in) :: level

      !> stencil header
      type(tem_stencilHeader_type), intent(in) :: stencil

      !> array of SOURCE elements
      real(kind=rk),    intent(in) :: sVal(:)

      !> array of TARGET GHOST elements
      real(kind=rk), intent(inout) :: tVal(:)

      !> List of target elements ( their position in depSource list )
      integer, intent(in) :: nTargets
      !> position in total list - offset
      integer, intent(in) :: targetList(nTargets)

      !> Number of scalars to interpolate
      integer, intent(in) :: nScalars
    end subroutine intpRoutine_arbitraryVal

  end interface
! ****************************************************************************** !

  contains

! ****************************************************************************** !
  !> Read in the type of interpolation scheme
  !!
  !!```lua
  !! interpolation_method = 'linear'               -- simple definition
  !! interpolation_method = {method ='debug', value = 1.} -- definition in a table
  !!```
  !!
  subroutine mus_load_interpolate( me, conf, parent )
    ! ---------------------------------------------------------------------------
    !> interpolation type to load info to
    type(mus_interpolation_config_type), intent(out)  :: me
    !> lua state to load from
    type(flu_state)                        :: conf
    !> optional parent table to load from
    integer, optional, intent(in)          :: parent
    ! ---------------------------------------------------------------------------
    character(len=32) :: localKey
    integer :: iError, nEntries, intp_handle
    ! ---------------------------------------------------------------------------

    localKey = 'interpolation_method'

    ! if more than one variable is defined in a table then
    ! ref_value should also be a table with reference value for each variable
    ! defined in the variable table.
    if( present( parent )) then
      call aot_table_open( L=conf, thandle=intp_handle, key=localKey,          &
        &                  parent= parent )
    else
      call aot_table_open( L=conf, thandle=intp_handle, key=localKey)
    endif

    ! not a table, just the name of the interpolation method
    if (intp_handle == 0) then
      call aot_get_val(L = conf, key = localKey,               &
        &              val = me%method,                        &
        &              ErrCode = iError, default = 'quadratic' )
    else
      ! interpolation is defined as a table with optional additional information
      ! Get size of table
      nEntries = aot_table_length(L=conf, thandle=intp_handle)
      call aot_get_val(L = conf, thandle = intp_handle,       &
        &              val = me%method, ErrCode = iError, &
        &              key = 'method', default = 'quadratic'  )

      call aot_get_val(L = conf, thandle = intp_handle, &
        &              val     = me%useComputeStencil,  &
        &              ErrCode = iError,                &
        &              key     = 'use_compute_stencil', &
        &              default = .false.                )

      call aot_get_val(L = conf, thandle = intp_handle, &
        &              val     = me%weights_method,     &
        &              ErrCode = iError,                &
        &              key     = 'weights_method',      &
        &              default = 'linear_distance'      )

      call aot_get_val(L = conf, thandle = intp_handle,       &
        &              val     = me%IDW_powerfac,             &
        &              ErrCode = iError,                      &
        &              key     = 'inverse_distance_powerfac', &
        &              default = 6                            )

      call aot_get_val(L=conf, thandle=intp_handle,                            &
        &              val=me%testInterpolation, ErrCode=iError,               &
        &              key = 'test', default=.false.)
      if( me%testInterpolation ) then
        write(logUnit(1),*)' Activated the Interpolation test.'
        call aot_get_val(L=conf, thandle=intp_handle,                          &
          &              val=me%testEachElement, ErrCode=iError,               &
          &              key = 'testEach', default=.false.)
        if( me%testEachElement   ) then
          write(logUnit(1),*)' Testing each element... '
        end if
        call aot_get_val(L=conf, thandle=intp_handle,                          &
          &              val=me%testFluids, ErrCode=iError,                    &
          &              key = 'testFluids', default=.false.)
        if( me%testFluids ) then
          write(logUnit(1),*)' Testing fluid element... '
        end if
      end if
      call aot_get_val(L=conf, thandle=intp_handle,                            &
        &              val=me%noIntpFromFiner, ErrCode=iError,                 &
        &              key = 'noIntpFromFiner', default=.false.)

      call aot_get_val(L=conf, thandle=intp_handle,                            &
        &              val=me%noIntpFromCoarser, ErrCode=iError,               &
        &              key = 'noIntpFromCoarser', default=.false.)

    endif
    call aot_table_close(L=conf, thandle=intp_handle)

    ! set interpolation order for fillFinerFromMe
    select case( trim(me%method) )
    case( 'weighted_average' )
      me%order = weighted_average
    case( 'linear' )
      me%order = linear
    case( 'quadratic' )
      me%order = quadratic
    case( 'none' )
      me%order = no_intp
    case default
      call tem_abort('From mus_load_interpolate: Unknown interpolation method')
    end select

    call interpolate_dump( me, logUnit(1) )

  end subroutine mus_load_interpolate
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine mus_set_nSources( me, nDims, QQ, layout)
    ! ---------------------------------------------------------------------------
    type(mus_interpolation_type), intent(inout) :: me
    integer,     intent(in)  :: nDims, QQ
    character(len=*), intent(in) :: layout
    ! ---------------------------------------------------------------------------
    integer :: iOrder
    ! ---------------------------------------------------------------------------

    ! Set nSources for fillCoarse
    me%fillMineFromFiner%nMinSources = 1
    select case(nDims)
    case(1)
      me%fillMineFromFiner%nMaxSources = 2
    case(2)
      me%fillMineFromFiner%nMaxSources = 4
    case(3)
      me%fillMineFromFiner%nMaxSources = 8
    end select

    ! Set nSources for fillFiner
    if (allocated(me%fillFinerFromMe)) deallocate(me%fillFinerFromMe)
    allocate(me%fillFinerFromMe(0:me%config%order))
    do iOrder = 0, me%config%order
      select case(iOrder)
      case(no_intp)
        ! do nothing
      case(weighted_average)
        me%fillFinerFromMe(iOrder)%nMinSources = 1
        select case(nDims)
        case(1)
          me%fillFinerFromMe(iOrder)%nMaxSources = 2
        case(2)
          me%fillFinerFromMe(iOrder)%nMaxSources = 4
          if (trim(layout) == 'd2q9') me%weightedAvgStencil%isActive = .true.
        case(3)
          select case(trim(layout))
          case('d3q19')
            me%fillFinerFromMe(iOrder)%nMaxSources = 7
            me%weightedAvgStencil%isActive = .true.
          case('d3q27')
            me%fillFinerFromMe(iOrder)%nMaxSources = 8
            me%weightedAvgStencil%isActive = .true.
          case default
            me%fillFinerFromMe(iOrder)%nMaxSources = QQ
          end select
        end select

        ! Use weighted average stencil only if QQ stencil is D2Q9, D3Q19
        ! or D3Q27
        ! else compute weight from all available sources
        if (me%weightedAvgStencil%isActive) then
          if ( allocated(me%weightedAvgStencil%neighDir) ) &
            & deallocate(me%weightedAvgStencil%neighDir)
          allocate(me%weightedAvgStencil%neighDir(me%fillFinerFromMe(iOrder) &
            &                                        %nMaxSources, 8))
          me%weightedAvgStencil%neighDir                              &
            & = init_cxDirWeightedAvg( QQ, me%fillFinerFromMe(iOrder) &
            &                                %nMaxSources             )
        end if

      case(linear)
        ! initialize least square matrix for linear interpolation
        call init( me     = me%fillFinerFromMe(iOrder)%intpMat_forLSF, &
          &        length = 1,                                         &
          &        nDims  = nDims,                                     &
          &        order  = linear                                     )

        me%fillFinerFromMe(iOrder)%nMinSources = me%fillFinerFromMe(iOrder)  &
          &                                         %intpMat_forLSF%nCoeffs

        ! Number of sources in weighted average stencil is enough for linear
        ! interpolation so use this stencil. Using weighted average stencil
        ! showed better result then using compute stencil
        if (me%config%useComputeStencil) then
          me%fillFinerFromMe(iOrder)%nMaxSources = QQ
          me%weightedAvgStencil%isActive = .false.
        else
          select case(nDims)
          case(1)
            me%fillFinerFromMe(iOrder)%nMaxSources = 2
          case(2)
            me%fillFinerFromMe(iOrder)%nMaxSources = 4
            if (trim(layout) == 'd2q9') me%weightedAvgStencil%isActive = .true.
          case(3)
            select case(trim(layout))
            case('d3q19')
              me%fillFinerFromMe(iOrder)%nMaxSources = 7
              me%weightedAvgStencil%isActive = .true.
            case('d3q27')
              me%fillFinerFromMe(iOrder)%nMaxSources = 8
              me%weightedAvgStencil%isActive = .true.
            case default
              me%fillFinerFromMe(iOrder)%nMaxSources = QQ
            end select
          end select
        end if

      case(quadratic)
        ! initialize least square matrix for quadratic interpolation
        call init( me     = me%fillFinerFromMe(iOrder)%intpMat_forLSF, &
          &        length = 1,                                         &
          &        nDims  = nDims,                                     &
          &        order  = quadratic                                  )

        me%fillFinerFromMe(iOrder)%nMinSources = me%fillFinerFromMe(iOrder)  &
          &                                        %intpMat_forLSF%nCoeffs
        me%fillFinerFromMe(iOrder)%nMaxSources = QQ

      case default
        call tem_abort('From mus_set_nSources: Unknown interpolation order')
      end select
    end do

    write(logUnit(1),"(A)") 'Setting for interpolation scheme: ' &
      &                      //trim(me%config%method)
    write(logUnit(3),"(A)") '  Number of sources from coarser: '
    write(logUnit(3),"(A,I0)") '     nMinSources: ', &
      &                        me%fillMineFromFiner%nMinSources
    write(logUnit(3),"(A,I0)") '     nMaxSources: ', &
      &                        me%fillMineFromFiner%nMaxSources
    write(logUnit(3),"(A)") '  Number of sources from finer: '
    do iOrder = 0, me%config%order
      write(logUnit(3),"(A,I0)") '    order: ', iOrder
      write(logUnit(3),"(A,I0)") '      nMinSources: ', &
      &                          me%fillFinerFromMe(iOrder)%nMinSources
      write(logUnit(3),"(A,I0)") '      nMaxSources: ', &
      &                          me%fillFinerFromMe(iOrder)%nMaxSources
    end do
    if (me%weightedAvgStencil%isActive) then
      write(logUnit(3),"(A)") '  use WeightedAvg stencil: T'
    end if

  end subroutine mus_set_nSources
! ****************************************************************************** !


! ****************************************************************************** !
  !> Dump interpolation method to lua
  subroutine mus_interpolate_out( me, conf )
    ! ---------------------------------------------------------------------------
    !> interpolation type to dump info to
    type(mus_interpolation_type), intent(in)  :: me
    !> aotus type handling the output to the file in lua format
    type(aot_out_type), optional, intent(inout) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_val( put_conf = conf, vname = 'interpolation_method', &
      &               val = trim(me%config%method) )

  end subroutine mus_interpolate_out
! ****************************************************************************** !

! ****************************************************************************** !
  !> Dump interpolation method to logUnit
  subroutine interpolate_dump( me, outUnit )
    ! ---------------------------------------------------------------------------
    !> interpolation type to dump info to
    type(mus_interpolation_config_type), intent(in)  :: me
    !> File unit to write to
    integer, intent(in) :: outUnit
    ! ---------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = outUnit)
    write(outUnit,"(A)") 'Interpolation:'
    write(outUnit,"(2A)") '  method: ',trim(me%method)
    write(outUnit,"(A,i2)")  '   order: ', me%order
    write(outUnit,"(2A)") '  weights_method: ',trim(me%weights_method)
    write(outUnit,"(A,i2)") '  inverse distance powefac: ', &
      &                           me%IDW_powerfac
    call tem_horizontalSpacer(fUnit = outUnit)

  end subroutine interpolate_dump
! ****************************************************************************** !

! ****************************************************************************** !
  !> Initialize stencil for weighted average interpolation
  function init_cxDirWeightedAvg( QQ, nSources ) result( me )
    ! ---------------------------------------------------------------------------
    integer, intent(in)  :: QQ
    integer, intent(in)  :: nSources
    integer              :: me(nSources,8)
    ! ---------------------------------------------------------------------------
    integer :: iChild
    integer,parameter :: q000_19 = 19  !< rest density is last for d3q19 stencil
    ! ---------------------------------------------------------------------------

    select case (QQ)
    case (9)
      me(:,:) = reshape( [ qNN, q0N, qN0, q00, & ! A
        &                   q0N, q1N, q00, q10, & ! B
        &                   qN0, q00, qN1, q01, & ! C
        &                   q00, q10, q01, q11, & ! D
        &                   qNN, q0N, qN0, q00, & ! A
        &                   q0N, q1N, q00, q10, & ! B
        &                   qN0, q00, qN1, q01, & ! C
        &                   q00, q10, q01, q11  ], [nSources,8] ) ! D
    case (19)
      me(:,:) = reshape( [ q000_19, qN00, q0N0, q00N, qNN0, qN0N, q0NN, & ! child 1
        &                  q000_19, q100, q0N0, q00N, q1N0, q10N, q0NN, & ! child 2
        &                  q000_19, qN00, q010, q00N, qN10, qN0N, q01N, & ! child 3
        &                  q000_19, q100, q010, q00N, q110, q10N, q01N, & ! child 4
        &                  q000_19, qN00, q0N0, q001, qNN0, qN01, q0N1, & ! child 5
        &                  q000_19, q100, q0N0, q001, q1N0, q101, q0N1, & ! child 6
        &                  q000_19, qN00, q010, q001, qN10, qN01, q011, & ! child 7
        &                  q000_19, q100, q010, q001, q110, q101, q011 ], & ! child 8
        &                                                      [nSources,8])
    case (27)
      me(:,:) = reshape( [ q000, qN00, q0N0, q00N, qNN0, qN0N, q0NN, qNNN, & ! child 1
        &                  q000, q100, q0N0, q00N, q1N0, q10N, q0NN, q1NN, & ! child 2
        &                  q000, qN00, q010, q00N, qN10, qN0N, q01N, qN1N, & ! child 3
        &                  q000, q100, q010, q00N, q110, q10N, q01N, q11N, & ! child 4
        &                  q000, qN00, q0N0, q001, qNN0, qN01, q0N1, qNN1, & ! child 5
        &                  q000, q100, q0N0, q001, q1N0, q101, q0N1, q1N1, & ! child 6
        &                  q000, qN00, q010, q001, qN10, qN01, q011, qN11, & ! child 7
        &                  q000, q100, q010, q001, q110, q101, q011, q111 ],& ! child 8
        &                                                      [nSources,8])
    case default
      write(logUnit(1),"(A)") 'Weighted average interpolation requires D2Q9, D3Q19 '&
        &                    // 'or D3Q27!'
      call tem_abort()
    end select

    write(dbgUnit(2),"(A)") ''
    write(dbgUnit(2),"(A)") ' linear stencil dir '
    do iChild = 1, 8
      ! me(1:3,iNeigh,iChild) = cxDir(1:3, dir(iNeigh, iChild) )
      write(dbgUnit(2), "(A,I0,A,4I3)")  'childNum: ', iChild, &
        &                                ', dir: ', me(1:nSources,iChild)
    end do
    write(dbgUnit(2),"(A)") ''
  end function init_cxDirWeightedAvg
! ****************************************************************************** !


end module mus_interpolate_header_module
! ****************************************************************************** !
