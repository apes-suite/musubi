! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2011-2012 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2013, 2015, 2018-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2019 Jana Gericke <jana.gericke@uni-siegen.de>
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
!! # Interpolation of flow quantities between different grid levels
!!
!! ## Interpolation
!!
!! The routines defined here, fill up the ghost elements with valid data.
!! Ghost elements are employed at grid level interfaces to provide valid
!! pdf values to the neighboring fluid elements. This way, the solvers can
!! act on elements of the same size only, treating the levels successively.
!! Target elements are the ghost elements, which have to be filled with
!! valid values.
!! Source elements are the fluid elements from other levels, from where to
!! take the input values for the interpolation.
!! The target ghost elements on the target level have corresponding source
!! fluid elements on the source level.
!!
!! [[tem_topology_module]] For a detailed description of the grid
!!
!! ## Workflow <a id="interpolation-of"></a>
!!
!! Each interpolation routine acts on a list of ghost elements.
!! This list contains pointers to the position in the total list.
!! For each of these ghost elements, the source elements are identified.
!! Before that, the sourceLevel is identified. However, the code is restricted
!! to work with a level jump of only one level, so the sourceLevel is
!! for sourceLevel = targetLevel+1
!! sourceLevel = targetLevel-1
!!
!! For an overview over implemented interpolation methods, see
!! [Interpolation methods](../page/features/intp_methods.html)
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
module mus_interpolate_module

  ! include treelm modules
  use mpi
  use env_module,            only: rk
  use tem_aux_module,        only: tem_abort
  use tem_debug_module,      only: main_debug, dbgUnit
  use tem_logging_module,    only: logUnit
  use tem_construction_module, only: tem_levelDesc_type, depSource_append
  use tem_dyn_array_module,  only: dyn_intArray_type, init, append, truncate
  use tem_grow_array_module, only: grw_intArray_type, init, append, truncate
  use tem_element_module,    only: eT_ghostFromCoarser, eT_ghostFromFiner
  use tem_stencil_module,    only: tem_stencilHeader_type
  use tem_matrix_module,     only: init, append
  use tem_tools_module,      only: tem_horizontalSpacer
  use tem_property_module,   only: prp_fluid, prp_fineGhostClosestToFluid

  ! include musubi modules
  use mus_field_prop_module,            only: mus_field_prop_type
  use mus_interpolate_header_module,    only: mus_interpolation_type,         &
    &                                         mus_interpolation_config_type,  &
    &                                         mus_interpolation_stencil_type, &
    &                                         weighted_average, linear,       &
    &                                         quadratic
  use mus_interpolate_debug_module,     only: do_nothing, do_nothing_arbi
  use mus_interpolate_average_module,   only:           &
    & fillArbiMyGhostsFromFiner_avg,                    &
    & fillMyGhostsFromFiner_avg_feq_fneq,               &
    & fillMyGhostsFromFiner_avgLES_feq_fneq,            &
    & fillMyGhostsFromFiner_avg2D_feq_fneq,             &
    & fillArbiFinerGhostsFromMe_weighAvg,               &
    & fillFinerGhostsFromMe_weighAvg_feq_fneq,          &
    & fillFinerGhostsFromMe_weighAvgLES_feq_fneq,       &
    & fillFinerGhostsFromMe_weighAvg2D_feq_fneq
  use mus_interpolate_linear_module,    only:         &
    & fillArbiFinerGhostsFromMe_linear,               &
    & fillArbiFinerGhostsFromMe_linear2D,             &
    & fillFinerGhostsFromMe_linear_feq_fneq,          &
    & fillFinerGhostsFromMe_linearLES_feq_fneq,       &
    & fillFinerGhostsFromMe_linear2D_feq_fneq
  use mus_interpolate_quadratic_module, only:       &
    & fillArbiFinerGhostsFromMe_quad,               &
    & fillArbiFinerGhostsFromMe_quad2D,             &
    & fillFinerGhostsFromMe_quad_feq_fneq,          &
    & fillFinerGhostsFromMe_quadLES_feq_fneq,       &
    & fillFinerGhostsFromMe_quad2D_feq_fneq
  use mus_scheme_header_module,     only: mus_scheme_header_type
  use mus_scheme_layout_module,     only: mus_scheme_layout_type

  implicit none

  private

  public :: mus_init_interpolate
  public :: mus_init_levelDescIntpArrays
  public :: mus_dump_levelDescIntp_nElems
  public :: mus_intp_update_depFromCoarser
  public :: mus_intp_update_depFromFiner

  contains

  ! ************************************************************************** !
  !> This subroutine initialzes the interpolation
  !!
  !! - setting the fillMineFromFiner and fillFinerFromMe function pointers
  !!   to the interpolation function chosen by the user and
  !! - pre-calculating the weights based on the distances between source
  !!   and target nodes
  !! - allocate buffers if required
  !!
  subroutine mus_init_interpolate( intp, levelDesc, schemeHeader, stencil, &
    &                              minLevel, maxLevel, fieldProp)
    ! ---------------------------------------------------------------------------
    !> interpolation type
    type(mus_interpolation_type), intent(inout) :: intp
    integer,                        intent(in) :: minLevel, maxLevel
    !> level descriptor is actually used
    type(tem_levelDesc_type), intent(inout) :: levelDesc(minLevel:maxLevel)
    !> the stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> scheme header
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    !> field properties
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    ! ---------------------------------------------------------------------------
    integer :: iOrder
    ! ---------------------------------------------------------------------------

    write(logUnit(1),"(A)") 'Initialize interpolation: '//trim(intp%config%method)
    write(logUnit(1),"(A)") '  with scheme: '//trim(schemeHeader%relaxation)
    flush(logUnit(1))

    select case (trim(schemeHeader%kind))
    case ('fluid', 'fluid_incompressible')
      if (fieldProp(1)%fluid%turbulence%active) then
        call assign_intp_fluidLES( intp  = intp,          &
          &                        nDims = stencil%nDims, &
          &                        QQ    = stencil%QQ     )
      else
        call assign_intp_fluid( intp  = intp,          &
          &                     nDims = stencil%nDims, &
          &                     QQ    = stencil%QQ     )
      end if
    case default
        write(logUnit(1),"(A)") ' Multilevel interpolation is not supported' &
          &                     // 'for scheme: '//trim(schemeHeader%kind)
        call tem_abort()
    end select

    ! Used for debugging. Configurable variable, Default is False
    if (intp%config%noIntpFromFiner) then
      intp%fillMineFromFiner%do_intp => do_nothing
      intp%fillMineFromFiner%do_intpArbiVal => do_nothing_arbi
    end if
    if (intp%config%noIntpFromCoarser) then
      do iOrder = 0, intp%config%order
        intp%fillFinerFromMe(iOrder)%do_intp => do_nothing
        intp%fillFinerFromMe(iOrder)%do_intpArbiVal => do_nothing_arbi
      end do
    end if

    write(logUnit(1),*) 'Done with interpolation method.'
  end subroutine mus_init_interpolate
  ! ************************************************************************** !

! ****************************************************************************** !
  !> Set up interpolation routines for fluid (weakly compressible) scheme
  subroutine assign_intp_fluid( intp, nDims, QQ )
    ! ---------------------------------------------------------------------------
    !> interpolation type
    type(mus_interpolation_type), intent(inout) :: intp
    !> Number of dimensions
    integer,                      intent(in) :: nDims
    !> number of stencil directions
    integer,                      intent(in) :: QQ
    ! ---------------------------------------------------------------------------
    integer :: ii
    ! ---------------------------------------------------------------------------
    write(logUnit(3),*)   ' Assign intp function pointer for fluid:'

    ! assign pointers to fillCoarser
    intp%fillMineFromFiner%do_intpArbiVal => fillArbiMyGhostsFromFiner_avg
    select case(nDims)
    case(2)
      write(logUnit(3),*) 'fillCoarser average 2D'
      intp%fillMineFromFiner%do_intp => fillMyGhostsFromFiner_avg2D_feq_fneq

    case(3)
      write(logUnit(3),*) 'fillCoarser average 3D'
      intp%fillMineFromFiner%do_intp => fillMyGhostsFromFiner_avg_feq_fneq

    case default
      call tem_abort('Unsupported nDims')
    end select

    ! assign pointers to fillFiner
    do ii = 0, intp%config%order
      select case(ii)
      case(weighted_average)
        intp%fillFinerFromMe(ii)%do_intpArbiVal &
          & => fillArbiFinerGhostsFromMe_weighAvg
        select case(nDims)
        case(2)
          write(logUnit(3),*) 'fillFiner weighted average 2D'
          intp%fillFinerFromMe(ii)                           &
            & %do_intp => fillFinerGhostsFromMe_weighAvg2D_feq_fneq

        case(3)
          write(logUnit(3),*) 'fillFiner weighted average 3D'
          intp%fillFinerFromMe(ii)%do_intp &
            & => fillFinerGhostsFromMe_weighAvg_feq_fneq
        end select

      case(linear)
        select case(nDims)
        case(2)
          write(logUnit(3),*) 'fillFiner linear 2D'
          intp%fillFinerFromMe(ii)%do_intpArbiVal    &
            & => fillArbiFinerGhostsFromMe_linear2D
          intp%fillFinerFromMe(ii)                         &
            & %do_intp => fillFinerGhostsFromMe_linear2D_feq_fneq

        case(3)
          write(logUnit(3),*) 'fillFiner linear 3D'
          intp%fillFinerFromMe(ii)%do_intpArbiVal &
            & => fillArbiFinerGhostsFromMe_linear
          intp%fillFinerFromMe(ii)%do_intp => &
            & fillFinerGhostsFromMe_linear_feq_fneq
        end select

      case(quadratic)
        select case(nDims)
        case(2)
          write(logUnit(3),*) 'fillFiner quadratic 2D'
          intp%fillFinerFromMe(ii)%do_intpArbiVal &
            & => fillArbiFinerGhostsFromMe_quad2D
          intp%fillFinerFromMe(ii)                       &
            & %do_intp => fillFinerGhostsFromMe_quad2D_feq_fneq

        case(3)
          write(logUnit(3),*) 'fillFiner quadratic 3D'
          intp%fillFinerFromMe(ii)%do_intpArbiVal &
            & => fillArbiFinerGhostsFromMe_quad
          intp%fillFinerFromMe(ii)%do_intp => fillFinerGhostsFromMe_quad_feq_fneq

        end select

      case default
        call tem_abort('Unknown interpolation method')
      end select

    end do

  end subroutine assign_intp_fluid
! ****************************************************************************** !

  !> Set up interpolation routines for fluid (weakly compressible) scheme and
  !! turbulence active
  subroutine assign_intp_fluidLES( intp, nDims, QQ )
    ! ---------------------------------------------------------------------------
    !> interpolation type
    type(mus_interpolation_type), intent(inout) :: intp
    !> Number of dimensions
    integer,                      intent(in) :: nDims
    !> number of stencil directions
    integer,                      intent(in) :: QQ
    ! ---------------------------------------------------------------------------
    integer :: ii
    ! ---------------------------------------------------------------------------
    write(logUnit(3),*)   ' Assign intp function pointer for fluid LES:'

    ! assign pointers to fillCoarser
    intp%fillMineFromFiner%do_intpArbiVal => fillArbiMyGhostsFromFiner_avg
    select case(nDims)
    case(3)
      write(logUnit(3),*) 'fillCoarser average 3D LES'
      intp%fillMineFromFiner%do_intp => fillMyGhostsFromFiner_avgLES_feq_fneq

    case default
      call tem_abort('Unsupported nDims')
    end select

    ! assign pointers to fillFiner
    do ii = 0, intp%config%order
      select case(ii)
      case(weighted_average)
        intp%fillFinerFromMe(ii)%do_intpArbiVal &
          & => fillArbiFinerGhostsFromMe_weighAvg
        select case(nDims)
        case(3)
          write(logUnit(3),*) 'fillFiner weighted average 3D LES'
          intp%fillFinerFromMe(ii)%do_intp &
            & => fillFinerGhostsFromMe_weighAvgLES_feq_fneq

        end select

      case(linear)
        select case(nDims)
        case(3)
          write(logUnit(3),*) 'fillFiner linear 3D LES'
          intp%fillFinerFromMe(ii)%do_intpArbiVal &
            & => fillArbiFinerGhostsFromMe_linear
          intp%fillFinerFromMe(ii)%do_intp => fillFinerGhostsFromMe_linearLES_feq_fneq

        end select

      case(quadratic)
        select case(nDims)
        case(3)
          intp%fillFinerFromMe(ii)%do_intpArbiVal &
            & => fillArbiFinerGhostsFromMe_quad
          write(logUnit(3),*) 'fillFiner quadratic 3D LES'
          intp%fillFinerFromMe(ii)%do_intp => fillFinerGhostsFromMe_quadLES_feq_fneq

        end select
      case default
        call tem_abort('Unknown interpolation method')
      end select

    end do

  end subroutine assign_intp_fluidLES
! ****************************************************************************** !

! ****************************************************************************** !
  !> Initialize levelwise ghost and source list for interpolation
  subroutine mus_init_levelDescIntpArrays( intp, levelDesc, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    !> interpolation type
    type(mus_interpolation_type), intent(inout)  :: intp
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    !> levelDesc
    type(tem_levelDesc_type), intent(inout) :: levelDesc(minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------
    write(logUnit(10),*) 'Initiale levelDesc intp arrays'
    do iLevel = minLevel, maxLevel
      call init_intpArraysPerLevel(                                   &
        & sourceFromCoarser    = levelDesc(iLevel)%sourceFromCoarser, &
        & intpFromCoarser      = levelDesc(iLevel)%intpFromCoarser,   &
        & sourceFromFiner      = levelDesc(iLevel)%sourceFromFiner,   &
        & intpFromFiner        = levelDesc(iLevel)%intpFromFiner,     &
        & nGhostFromCoarser    = levelDesc(iLevel)%elem               &
        &                        %nElems( eT_ghostFromCoarser ),      &
        & nGhostFromFiner      = levelDesc(iLevel)%elem               &
        &                        %nElems( eT_ghostFromFiner ),        &
        & intp_order           = intp%config%order,                   &
        & nMaxSourcesFromFiner = intp%fillMineFromFiner%nMaxSources,  &
        & iLevel               = iLevel                               )
    end do !iLevel
  end subroutine mus_init_levelDescIntpArrays
! ****************************************************************************** !


! ****************************************************************************** !
  subroutine init_intpArraysPerLevel( sourceFromCoarser, intpFromCoarser,      &
    &                                 sourceFromFiner, intpFromFiner,          &
    &                                 nGhostFromCoarser, nGhostFromFiner,      &
    &                                 intp_order, nMaxSourcesFromFiner, iLevel )
    ! --------------------------------------------------------------------------
    !> dynamic array of source elements from coarser
    type(dyn_intArray_type), intent(out) :: sourceFromCoarser
    !> growing array of ghost elements from coarser for different intp order
    type(grw_intArray_type), allocatable, intent(out) :: intpFromCoarser(:)
    !> dynamic array of source elements from finer
    type(dyn_intArray_type), intent(out) :: sourceFromFiner
    !> growing array of ghost elements from finer
    type(grw_intArray_type), intent(out) :: intpFromFiner
    !> Number of ghost from coarser on this level
    integer, intent(in) :: nGhostFromCoarser
    !> Number of ghost from finer on this level
    integer, intent(in) :: nGhostFromFiner
    !> current level
    integer, intent(in) :: iLevel
    !> interpolation order defined by user
    integer, intent(in) :: intp_order
    !> nMaxSources for fillMineFromFiner
    integer, intent(in) :: nMaxSourcesFromFiner
    ! --------------------------------------------------------------------------
    integer :: iOrder
    ! --------------------------------------------------------------------------
    write(logUnit(10),*) 'Initiale levelDesc intp array on level:', iLevel

    ! fromCoarser
    call init( me = sourceFromCoarser, length = nGhostFromCoarser )

    allocate(intpFromCoarser(0:intp_order))
    ! initialize only user defined order with nGhostFromCoarser since
    ! order below user defined order are fall back which will be small in size
    call init( me     = intpFromCoarser(intp_order), &
      &        length = nGhostFromCoarser            )
    do iOrder = 0, intp_order-1
      call init( me = intpFromCoarser(iOrder) )
    end do

    ! fromFiner
    call init( me     = sourceFromFiner,                       &
      &        length = nGhostFromFiner * nMaxSourcesFromFiner )

    call init( me  = intpFromFiner, length = nGhostFromFiner )

  end subroutine init_intpArraysPerLevel
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine dumps global nElems in intpFromCoarser, intpFromFiner,
  !! sourcesFromCoarser ans sourcesFromFiner
  subroutine mus_dump_levelDescIntp_nElems(intp, levelDesc, minLevel, maxLevel,&
    &                                      root, comm                          )
    ! --------------------------------------------------------------------------
    !> interpolation type
    type(mus_interpolation_type), intent(in)  :: intp
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    !> levelDesc
    type(tem_levelDesc_type), intent(inout) :: levelDesc(minLevel:maxLevel)
    integer, intent(in) :: root !< root process
    integer, intent(in) :: comm !< mpi communicator
    ! --------------------------------------------------------------------------
    integer :: iLevel, iOrder
    ! amount of High and Low Oder Ghost elements form Fine and Coarser
    integer :: nGFC(0:intp%config%order), tSFC, nGFF, tSFF, iErr
    ! ---------------------------------------------------------------------------
    do iLevel = minLevel, maxLevel
      call truncate( me = levelDesc(iLevel)%sourceFromCoarser    )
      do iOrder = 0, intp%config%order
        call truncate( me = levelDesc(iLevel)%intpFromCoarser(iOrder) )
      end do

      call truncate( me = levelDesc(iLevel)%sourceFromFiner    )
      call truncate( me = levelDesc(iLevel)%intpFromFiner )

      ! get ghostFromCoarser for every order ( nGFC )
      do iOrder = 0, intp%config%order
        call MPI_REDUCE( levelDesc(iLevel)%intpFromCoarser(iOrder)%nVals, &
          &              nGFC(iOrder), 1, MPI_INTEGER, mpi_sum,           &
          &              root, comm, iErr)
      end do

      ! get total sourceFromCoarser ( tSFC )
      call MPI_REDUCE( levelDesc( iLevel )%sourceFromCoarser%nVals, &
        &              tSFC, 1, MPI_INTEGER, mpi_sum,               &
        &              root, comm, iErr)

      ! get ghostFromFiner ( nGFF )
      call MPI_REDUCE( levelDesc( iLevel )%intpFromFiner%nVals, &
        &              nGFF, 1, MPI_INTEGER, mpi_sum,           &
        &              root, comm, iErr)

      ! get total sourceFromFiner ( tSFF )
      call MPI_REDUCE( levelDesc( iLevel )%sourceFromFiner%nVals, &
        &              tSFF, 1, MPI_INTEGER, mpi_sum,             &
        &              root, comm, iErr)

      write(logUnit(3),"(A,I0)") '                       level: ', iLevel
      do iOrder = 0, intp%config%order
        write(logUnit(3),"(A,I0)") '         Interpolation order: ', iOrder
        write(logUnit(3),"(A,I0)") '              ghostFromCoarser: ', &
          &                        nGFC(iOrder)
      end do
      write(logUnit(3),"(A,I0)") '     total sourceFromCoarser: ', tSFC
      write(logUnit(3),"(A,I0)") '              ghostFromFiner: ', nGFF
      write(logUnit(3),"(A,I0)") '     total   sourceFromFiner: ', tSFF
    end do ! iLevel = minLevel, maxLevel
    ! ------------------------------------------------------------------------

  end subroutine mus_dump_levelDescIntp_nElems
! ****************************************************************************** !


! ****************************************************************************** !
  !> The required source elements for ghost from coarser elements are
  !! identified in this routine. Moreover, the weights for each sources based
  !! on distance are calculated. \n
  !!
  !! The depdendencies for ghostFromCoarser elements initially only include its
  !! parent element.
  !! Additional dependencies Y for the element x are found which were not
  !! included in the construction routine before hand. However, they are
  !! defined to be inherently part of the stencil and hence exist locally and
  !! can simply be added to the dependencies. Two different ways are
  !! implemented.
  !! Full interpolation stencil
  !!   the complete compute stencil is included into the dependency list
  !!```
  !!  Y       Y       Y
  !!
  !!        x   o
  !!  Y       O       Y
  !!        o   o
  !!
  !!  Y       Y       Y
  !!```
  !!
  !! Note: The view here is according to how the interpolation will be done.
  !!       The coarse level does the interpolation for the finer level.
  !!       Thus, we access the information from the levelDesc with the
  !!       targetLevel and sourceLevel perspective.
  !!
  !! Steps of finding sources:
  !! 1.  loop over all neighbors: iNeigh = 1, QQ
  !! 2.  get neighbor posInTotal: elemPos = neigh(1)%nghElems
  !! 3.  if elemPos > 0,
  !!       add this source into depFromCoarser%elem
  !!       add this source into sourceFromCoarser
  !!       save its position in sourceFromCoarser into depFromCoarser%elemBuffer
  !!       calculate its weight based on relative coordinates in stencil
  !!       weight is inverse proportional to the distance from source to target
  !!       weight is normalized at last, so that sum of weights equals 1
  !! 4.  save its offset on stencil
  !!
  !! After finding sources, append target to low or high interpolation list
  !! based on if all required sources are found.
  !!
  subroutine mus_intp_update_depFromCoarser( intp, levelDesc, stencil, &
    &                                        minLevel, maxLevel )
    ! ---------------------------------------------------------------------------
    !> interpolation type
    type(mus_interpolation_type), intent(inout)  :: intp
    integer,                        intent(in) :: minLevel
    integer,                        intent(in) :: maxLevel
    !> Level Descriptor
    type( tem_levelDesc_type )                 :: levelDesc(minLevel:maxLevel)
    type( tem_stencilHeader_type ), intent(in) :: stencil
    ! ---------------------------------------------------------------------------
    integer :: sourceLevel  ! level of source elements
    integer :: targetLevel  ! level of target elements
    integer :: targetElem   ! position of target element in total list
    integer :: iElem, iNeigh
    integer :: nGhostElems
    integer :: nFoundSources
    integer :: posInTotal, parentPosInTotal
    real(kind=rk) :: targetBary(3) ! barycenter of current target element
    integer :: mySources(stencil%QQ), myNeighDir(stencil%QQ)
    ! ---------------------------------------------------------------------------
    integer :: posInMatArray, intpOrder, childNum
    integer :: iOrder
    logical :: success
    integer :: sourceElem, iSourceElem
    ! ---------------------------------------------------------------------------

    write(logUnit(1),"(A)") 'Updating dependences of ghost from coarser'

! DEBUG output  ----------------------------------------------------------------
call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
write(dbgUnit(3),"(A)") "Inside mus_intp_update_depFromCoarser routine"
write(dbgUnit(3),"(A)") "Up to now, for each ghostFromCoarser we only have its parent."
write(dbgUnit(3),"(A)") "Here we are trying to find all neighbors of this parent."
write(dbgUnit(3),"(A,I0)") ''
! DEBUG output  ----------------------------------------------------------------

    do sourceLevel = minLevel, maxLevel - 1

      targetLevel = sourceLevel + 1

      nGhostElems = levelDesc( targetLevel )%elem%nElems( eT_ghostFromCoarser )
write(dbgUnit(4),"(A,I0)") " source level: ", sourceLevel
write(dbgUnit(4),"(A,I0)") " target level: ", targetLevel
write(dbgUnit(4),"(A,I0)") ' nGhostFromCoarser: ', nGhostElems

      !!suppChild = 1
      !do iElem = 1, nGhostElems
      !  !childNum = levelDesc(targetLevel)%depFromCoarser(iElem)%childNum
      !  !allocate(levelDesc( targetLevel )%depFromCoarser( iElem )%bitmask(stencil%QQ))
      !  !levelDesc( targetLevel )%depFromCoarser( iElem )%bitmask = .true.
      !
      !  !if (suppChild == 1) supp_bitmask = .true.
      !
      !  targetElem = iElem + levelDesc( targetLevel )%             &
      !    &                                offset( 1, eT_ghostFromCoarser)

!write(dbgUnit(5),"(A, I0)") '      iElem: ', iElem
!!write(dbgUnit(5),"(A, I0)") '      suppChild: ', suppChild
!!write(dbgUnit(5),"(A, I0)") '      childNum: ', childNum
!flush(dbgUnit(5))

      !  do iNeigh = 1, stencil%QQ-1 ! no rest position
      !
      !    !if ( iNeigh == stencil%restPosition ) then
      !      !! for rest position
      !      !supp_bitmask(iNeigh) = .true.
      !    !else
      !      ! Find neighor using neigh array because for boundary, neigh array
      !      ! points to current element so no need to check for existence of the
      !      ! neighbor element
      !      ! check for both first and second layer of ghost cells
      !      ! the second condition is true only if the secon neigh is fluid
      !      nghElem = levelDesc( targetLevel )%neigh(1)%nghElems( iNeigh, targetElem )
      !      if (nghElem > 0) then ! nghElem exists
      !        if ( btest(levelDesc( targetLevel )%property( nghElem ), prp_fluid) ) then
      !          ! nghElem is
      !          !invDir = stencil%cxDirInv( iNeigh )
      !          !supp_bitmask(invDir) = .false.
      !          ! then this is the closest fine ghost cell to the fluid cell
      !          levelDesc( targetLevel )%property( targetElem ) = &
      !            & ibset( levelDesc( targetLevel )%property( targetElem ), prp_fineGhostClosestToFluid )
      !         else
      !          ! find the second layer of closest fine ghosts to fluid
      !          nghElem_2 = levelDesc( targetLevel )%neigh(1)%nghElems( iNeigh, nghElem )
      !          if ( nghElem_2 > 0 ) then
      !            if ( btest(levelDesc( targetLevel )%property( nghElem_2 ), prp_fluid) ) then
      !      !        supp_bitmask = .true.
      !              levelDesc( targetLevel )%property( targetElem ) = &
      !            & ibset( levelDesc( targetLevel )%property( targetElem ), prp_fineGhostClosestToFluid )
      !            end if
      !          end if
      !        end if
      !      end if
      !    !end if
      !  end do ! iNeigh
      !
      !
      !  !if (childNum == 8) then
      !  !  do iChild = 1, suppChild
      !  !    levelDesc( targetLevel )%depFromCoarser( iElem - suppChild + iChild )%bitmask = supp_bitmask
      !  !  end do
      !  !  suppChild = 1
      !  !else
      !  !  suppChild = suppChild + 1
      !  !end if
      !
      !end do

      ! Treat all the fromCoarser elements
      do iElem = 1, nGhostElems

        parentPosInTotal = &
          & levelDesc( targetLevel )%depFromCoarser( iElem )%elem%val(1)
        targetBary = levelDesc( targetLevel )%depFromCoarser( iElem )%coord(:)
        childNum = levelDesc(targetLevel)%depFromCoarser(iElem)%childNum

        ! target element treeID
write(dbgUnit(4),"(A,I0)") ''
        ! target element position in total list
        targetElem = iElem + levelDesc( targetLevel )%             &
          &                                offset( 1, eT_ghostFromCoarser)

write(dbgUnit(4),"(A,I0)") ' target treeID: ', leveldesc( targetlevel )%total( targetElem )
write(dbgUnit(4),"(A,I0)") ' child number: ', childNum
write(dbgUnit(4),"(A,3F8.3)") ' coord relative to parent: ', targetBary(1:3)
write(dbgUnit(4),"(A,I0)") ' number of sources (parent): ', &
  &             levelDesc( targetLevel )%depFromCoarser( iElem )%elem%nVals
write(dbgUnit(4),"(A,I0)") ' parent treeID: ', levelDesc( sourceLevel )%total( parentPosInTotal )
write(dbgUnit(4),"(A   )") ''

        ! --------------- treat parent element (first source) --------------
        ! start with treating the parent element, which is located at the
        ! first position in the dependency list
        ! Later on, proceed to siblings of parent in the target element's
        ! directions from its parent element
        ! add the parent to the list of source elements on this level
        nFoundSources = 0
        mySources = 0
        myNeighDir = -1
        ! Now loop over all neighbors and find maximum number of source elements
        do iNeigh = 1, stencil%QQ

write(dbgUnit(5),"(A, I0)") '      iNeigh: ', iNeigh

          if ( iNeigh == stencil%restPosition ) then
            ! normally parent is the last source due to stencil definition
            posInTotal = parentPosInTotal
          else
            ! this is a neihgbor of the parent element ->
            ! access neighbor information.
            posInTotal = levelDesc( sourceLevel )%neigh(1)%       &
              &                   nghElems( iNeigh, parentPosInTotal )
          end if

write(dbgUnit(5),"(A, I0)") '  posInTotal: ', posInTotal
          if ( posInTotal > 0 ) then
write(dbgUnit(5),"(A, I0)") '  sourceID: ', levelDesc(sourceLevel)%total(posInTotal)
            nFoundSources = nFoundSources + 1
            mySources( nFoundSources ) = posInTotal
            myNeighDir( nFoundSources ) = iNeigh
          end if ! posIntotal > 0, i.e. this source is a valid one
        end do ! iNeigh

!write(dbgUnit(3),*) '  bitmask = ', levelDesc( targetLevel )%depFromCoarser( iElem )%bitmask(:)

        ! if sources are found then determine which interpolation to use
        ! depending on available sources
        if ( nFoundSources /= 0 ) then
write(dbgUnit(4),"(A, I0)")   '      nFound: ', nFoundSources
write(dbgUnit(4),*) '   myNeighDir: ', myNeighDir(1:nFoundSources)
write(dbgUnit(4),*) '    mySources: ', mySources(1:nFoundSources)

          ! Find maximum possible interpolation order for current nFoundSources
          ! and update mySources, neighDir and nFoundSources if all sources
          ! in stencil are found.
          call find_possIntpOrderAndUpdateMySources(          &
            &    intp          = intp,                        &
            &    mySources     = mySources(1:nFoundSources),  &
            &    neighDir      = myNeighDir(1:nFoundSources), &
            &    nFoundSources = nFoundSources,               &
            &    childNum      = childNum,                    &
            &    intpOrder     = intpOrder                    )
write(dbgUnit(4),"(A,I0)") '   intpOrder: ', intpOrder
write(dbgUnit(4),"(A, I0)")   '      nFound_new: ', nFoundSources


          ! Compute interpolation matrix for least square fit using stencil
          ! direction of available sources
          ! The parent of target childs coord is 0,0,0 so we could
          ! just use of stencil%cxDir to build up this matrix entries
          ! Every row in matrix is evaluated with coord of source element
          if (intpOrder == quadratic) then
            call append( me       = intp%fillFinerFromMe(intpOrder) &
              &                          %intpMat_forLSF,           &
              &          order    = intpOrder,                      &
              &          QQ       = stencil%QQ,                     &
              &          nDims    = stencil%nDims,                  &
              &          nSources = nFoundSources,                  &
              &          cxDirRK  = stencil%cxDirRK,                &
              &          neighDir = myNeighDir(1:nFoundSources),    &
              &          pos      = posInMatArray,                  &
              &          success  = success                         )
             if (success) then
               levelDesc( targetLevel )%depFromCoarser(iElem)         &
                 &                     %posInIntpMatLSF = posInMatArray
             else
               ! matrix is singular then fall back to linear
               intpOrder = linear
             end if
          end if

          if (intpOrder == linear) then
            call append( me       = intp%fillFinerFromMe(intpOrder) &
              &                          %intpMat_forLSF,           &
              &          order    = intpOrder,                      &
              &          QQ       = stencil%QQ,                     &
              &          nDims    = stencil%nDims,                  &
              &          nSources = nFoundSources,                  &
              &          cxDirRK  = stencil%cxDirRK,                &
              &          neighDir = myNeighDir(1:nFoundSources),    &
              &          pos      = posInMatArray,                  &
              &          success  = success                         )
             if (success) then
               levelDesc( targetLevel )%depFromCoarser(iElem)         &
                 &                     %posInIntpMatLSF = posInMatArray
             else
               ! matrix is singular then fall back to weighted average
               intpOrder = weighted_average
             end if
          end if

          if (intpOrder == weighted_average) then
            ! Compute weights for weighted average interperation.
            ! Weights are computed according to available sources.
            call compute_weight( &
              &    weights       = levelDesc(targetLevel)         &
              &                    %depFromCoarser(iElem)%weight, &
              &    nFoundSources = nFoundSources,                 &
              &    neighDir      = myNeighDir(1:nFoundSources),   &
              &    targetBary    = targetBary,                    &
              &    intp_config   = intp%config,                   &
              &    cxDirRK       = stencil%cxDirRK                )
          end if

write(dbgUnit(4),"(A,I0)") '   Final intpOrder: ', intpOrder
          ! assign current element to correct interpolation scheme according
          ! to available sources
          call append( me  = levelDesc( targetLevel )     &
            &                %intpFromCoarser(intpOrder), &
            &          val = iElem                        )

        else ! no Sources
          write(logUnit(1),*) 'Can NOT find any source for this ghostFromCoarser.'
          write(logUnit(1),*) 'Thus can NOT do interpolation.'
          call tem_abort()
        end if

        call depSource_append(                                                &
          &    me         = levelDesc( targetLevel )%depFromCoarser( iElem ), &
          &    sourceList = levelDesc( targetLevel )%sourceFromCoarser,       &
          &    mySources  = mySources(1:nFoundSources),                       &
          &    n          = nFoundSources                                     )

      end do ! iElem = 1, nGhostElems

      write(dbgUnit(4),"(A,I0)") ' level: ', targetLevel
      write(dbgUnit(4),"(A,I0)") ' number of ghost elements: ', nGhostElems
      do iOrder = 0, intp%config%order
        write(dbgUnit(4),"(A,I0)") ' ghostFromCoarser: ', &
          &      levelDesc( targetLevel )%intpFromCoarser(iOrder)%nVals
      end do
      write(dbgUnit(4),"(A,I0)") ' total sources: ',      &
        &      levelDesc( targetLevel )%sourceFromCoarser%nVals

    end do ! level loop

    do sourceLevel = minLevel, maxLevel-1
      targetLevel = sourceLevel + 1
      write(dbgUnit(4),*) 'sourceIDs for target level: ', targetLevel
      do iSourceElem = 1, levelDesc(targetLevel)%sourceFromCoarser%nVals
        sourceElem = levelDesc(targetLevel)%sourceFromCoarser%val(iSourceElem)
        write(dbgUnit(4),*) iSourceElem, sourceElem, &
          & levelDesc(sourceLevel)%total(sourceElem)
      end do
    end do

    write(dbgUnit(4),"(A)") "Outside of mus_intp_update_depFromCoarser routine"
    call tem_horizontalSpacer( fUnit = dbgUnit(4), after = 1 )
  end subroutine mus_intp_update_depFromCoarser
! ****************************************************************************** !


! **************************************************************************** !
  !> Find maximum possible interpolation order which can be used for
  !! fillFinerFromMe by comparing nFoundSources with nMaxSources of different
  !! interpolation order starting from interpolation order defined by user
  subroutine find_possIntpOrderAndUpdateMySources( intp, mySources, neighDir, &
    &                                              nFoundSources, childNum,   &
    &                                              intpOrder                  )
    ! ---------------------------------------------------------------------------
    !> interpolation type
    type(mus_interpolation_type), intent(in)  :: intp
    !> Number of source elements found
    integer, intent(inout) :: nFoundSources
    !> position of found source elements in total list
    !! Update this list if intpStencil%isActive
    integer, intent(inout) :: mySources(:)
    !> cxDir for found sounce elements
    integer, intent(inout) :: neighDir(:)
    !> Curent finer ghost child number
    integer, intent(in) :: childNum
    !> interpolation order
    integer, intent(out) :: intpOrder
    ! --------------------------------------------------------------------------
    integer :: iOrder, intpOrder_tmp
    logical :: allSrcFound
    ! --------------------------------------------------------------------------
    ! Start with interpolation order defined by user and if requirement
    ! for that interpolation is not satisfied then fall back to lower
    ! order
    intpOrder = intp%config%order

    ! if interpolation order <= quadratic then
    ! Find which order to use for this ghost element depending on
    ! order defined in config and number of sources found
    ! Set to weighted_average if nSources does not satisfy
    ! for any higher order
    intpOrder_tmp = weighted_average
    do iOrder = intpOrder, linear, -1
      if ( nFoundSources >= intp%fillFinerFromMe(iOrder)%nMinSources ) then
        intpOrder_tmp = iOrder
        exit
      end if
    end do
    intpOrder = intpOrder_tmp

    ! if weightedAvgStencil is defined then use only neighDir from
    ! stencilNeighDir as found sources for weighted_average.
    ! If all sources required by the stencil are not found then use found
    ! sources
    if ( (intpOrder == weighted_average  .or. intpOrder == linear ) &
      & .and. intp%weightedAvgStencil%isActive ) then
      call updateMySources(                                         &
        &    mySources     = mySources,                             &
        &    neighDir      = neighDir,                              &
        &    nFoundSources = nFoundSources,                         &
        &    childNum      = childNum,                              &
        &    nMaxSources   = intp%fillFinerFromMe(weighted_average) &
        &                        %nMaxSources,                      &
        &    intpStencil   = intp%weightedAvgStencil,               &
        &    allSrcFound   = allSrcFound                            )
    end if

  contains

    ! ************************************************************************ !
    !> Update mySources if all sources in required by the stencil are found
    subroutine updateMySources(mySources, neighDir, nFoundSources, childNum, &
      &                        nMaxSources, intpStencil, allSrcFound)
      ! ------------------------------------------------------------------------
      !> Number of source elements found
      integer, intent(inout) :: nFoundSources
      !> position of found source elements in total list
      !! Update this list if intpStencil%isActive
      integer, intent(inout) :: mySources(:)
      !> cxDir for found sounce elements
      integer, intent(inout) :: neighDir(:)
      !> ChildNum of current target fine element
      integer, intent(in) :: childNum
      !> Maximum number of sources required for this intpStencil
      integer, intent(in) :: nMaxSources
      !> Interpolation stencil
      type(mus_interpolation_stencil_type), intent(in) :: intpStencil
      !> is true if all sources required by the stencil are found
      logical, intent(out) :: allSrcFound
      ! ------------------------------------------------------------------------
      integer :: iSrc, nFoundSources_new
      integer :: mySources_tmp(nFoundSources), neighDir_tmp(nFoundSources)
      ! ------------------------------------------------------------------------
      allSrcFound = .false.
      nFoundSources_new = 0

      do iSrc = 1, nFoundSources
        if ( any(intpStencil%neighDir(1:nMaxSources, childNum) &
          & == neighDir(iSrc)) ) then
          nFoundSources_new = nFoundSources_new + 1
          neighDir_tmp(nFoundSources_new) = neighDir(iSrc)
          mySources_tmp(nFoundSources_new) = mySources(iSrc)
        end if
      end do
      ! all sources required by this stencil are found. update sources, neighdir
      ! and nFoundSources
      if (nFoundSources_new == nMaxSources) then
        allSrcFound = .true.
        nFoundSources = nFoundSources_new
        mySources(1:nFoundSources) = mySources_tmp(1:nFoundSources)
        neighDir(1:nFoundSources) = neighDir_tmp(1:nFoundSources)
      end if

    end subroutine updateMySources
    ! ************************************************************************ !

  end subroutine find_possIntpOrderAndUpdateMySources
! **************************************************************************** !

! ****************************************************************************** !
  !> This routine computes weights for weighted_average interpolation.
  subroutine compute_weight( weights, nFoundSources, neighDir, targetBary, &
    &                        intp_config, cxDirRK )
    ! ---------------------------------------------------------------------------
    !> computed weight
    real(kind=rk), allocatable, intent(out) :: weights(:)
    !> Number of source elements found
    integer, intent(in) :: nFoundSources
    !> cxDir for found sounce elements
    integer, intent(in) :: neighDir(nFoundSources)
    !> child bary relative to parent
    real(kind=rk), intent(in) :: targetBary(3)
    !> Interpolation config info
    type(mus_interpolation_config_type), intent(in) :: intp_config
    !> cxDir of current source
    real(kind=rk), intent(in) :: cxDirRK(:,:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: scalar_dist
    !> absolute distances
    real(kind=rk) :: distance(3)
    real(kind=rk)  :: sumWeight
    integer :: iSrc, iNeigh
    ! ---------------------------------------------------------------------------
write(dbgUnit(4),"(A)") "Inside compute weight"
    allocate( weights( nFoundSources ) )
    do iSrc = 1, nFoundSources
      iNeigh = neighDir(iSrc)
      ! distance from this source to target
      ! coordinates of source is indeed its offset to parent of target
      ! For 2D: use only X and Y
      distance(1:3) = abs( cxDirRK(1:3, iNeigh) - targetBary )

      scalar_dist = sqrt(dot_product(distance, distance))

      select case ( trim(intp_config%weights_method) )
        case( 'inverse_distance' )
          weights(iSrc) = 1.0_rk / scalar_dist**intp_config%IDW_powerfac
        case( 'modified_shepards' )
          weights(iSrc) = ( max( 0.0_rk, 1.0_rk - scalar_dist ) &
            &           / scalar_dist )**2
        case( 'linear_distance')
          weights(iSrc) = (1.0_rk - distance(1)) * (1.0_rk-distance(2)) &
            &           * (1.0_rk - distance(3))
        case default
          write(*,*) 'Unknown weighting method, check inputs. Aborting ...'
          write(logUnit(0),*) 'weighting_method in lua file should be:'
          write(logUnit(0),*) 'inverse_distance'
          write(logUnit(0),*) 'modified_shepards'
          write(logUnit(0),*) 'Aborting ...'
          call tem_abort()
      end select
write(dbgUnit(4),"(A,3F7.3)") '    distance: ', distance(1:3)
write(dbgUnit(4),"(A,3F7.3)") '      weight: ', weights(iSrc)
    end do !iSrc

    ! Normalize weight, so that sum(weight) = 1
    sumWeight = sum( weights(1:nFoundSources) )
    weights(:) = weights( 1:nFoundSources ) / sumWeight
write(dbgUnit(7),*) '   Weights: ', weights

  end subroutine compute_weight
! ****************************************************************************** !

! ****************************************************************************** !
  !> All sources (children) have been found in treelm,
  !! updated number of sources needed, based on nDims
  !! collect all source elements into sourceFromFiner
  !! assign ghost intp list. Currently only average interpolation is implemented
  !! for fillMineFromFiner
  !!
  subroutine mus_intp_update_depFromFiner( intp, levelDesc, minLevel, maxLevel, &
    & stencil )
    ! ---------------------------------------------------------------------------
    !> interpolation type
    type(mus_interpolation_type), intent(inout)  :: intp
    integer,                        intent(in) :: minLevel
    integer,                        intent(in) :: maxLevel
    !> Level Descriptor
    type( tem_levelDesc_type )                 :: levelDesc(minLevel:maxLevel)
    !> Stencil
    type( tem_stencilHeader_type ), intent(in) :: stencil
    ! ---------------------------------------------------------------------------
    integer :: iLevel, iElem
    integer :: nSources   ! number of source elements actually appended
    integer :: mySources( intp%fillMineFromFiner%nMaxSources )
    integer :: iSourceElem, sourceElem
    !integer :: nghElem, targetElem, iNeigh, invDir
    ! ---------------------------------------------------------------------------

    write(logUnit(1),"(A)") 'Updating dependencies of ghost from finer.'
call tem_horizontalSpacer( fUnit = dbgUnit(4), before = 1 )
    write(dbgUnit(4),"(A)") 'Updating dependencies of ghost from finer.'

    ! all levels except highest level
    do iLevel = minLevel, maxLevel-1

      do iElem = 1, levelDesc( iLevel )%elem%nElems( eT_ghostFromFiner )

        ! Depending of dimension elements are added to source list
        ! nMaxSources in fillMineFromFiner is set according to nDims
        ! for 3D, all 8 elements are considered
        ! for 2D, only first 4 elements are considered
        ! for 1D, only first 2 elements are considered.
        nSources = min( intp%fillMineFromFiner%nMaxSources,                  &
          &             levelDesc( iLevel )%depFromFiner( iElem )%elem%nVals )
        mySources(1:nSources) = levelDesc( iLevel )%depFromFiner( iElem ) &
          &                                        %elem%val(1:nSources)

        call depSource_append(                                       &
          &  me         = levelDesc( iLevel )%depFromFiner( iElem ), &
          &  sourceList = levelDesc( iLevel )%sourceFromFiner,       &
          &  mySources  = mySources(1:nSources),                     &
          &  n          = nSources                                   )

        ! add into interpolation list
        call append( me  = levelDesc( iLevel )%intpFromFiner, &
          &          val = iElem )

        !targetElem = iElem + levelDesc( iLevel )%             &
        !  &                                offset( 1, eT_ghostFromFiner)

        !write(logUnit(1),"(A)") 'Create mask for coarse ghosts.'
        !flush(logunit(1))
        !allocate(levelDesc( iLevel )%depFromFiner( iElem )%bitmask(stencil%QQ))
        !levelDesc( iLevel )%depFromFiner( iElem )%bitmask = .false.

        ! Now loop over all neighbors
        !do iNeigh = 1, stencil%QQ-1 !no rest position
        !  !if ( iNeigh == stencil%restPosition ) then
        !    ! for rest position coarse ghost has the correct value
        !    !levelDesc( iLevel )%depFromFiner( iElem )%bitmask(iNeigh) = .true.
        !  !else if ( .not. iNeigh == stencil%restPosition ) then
        !    ! Find neighor using neigh array because for boundary, neigh array
        !    ! points to current element so no need to check for existence of the
        !    ! neighbor element
        !    ! check for both first and second layer of ghost cells
        !    ! the second condition is true only if the secon neigh is fluid
        !    nghElem = levelDesc( iLevel )%neigh(1)%nghElems( iNeigh, targetElem )
        !    if (nghElem > 0) then
        !      if ( btest(levelDesc( iLevel )%property( nghElem ), prp_fluid) ) then
        !        levelDesc( iLevel )%depFromFiner( iElem )%bitmask(iNeigh) = .true.
        !      end if
        !    end if
        !  !end if
        !end do

!write(dbgUnit(3),*) '  bitmask = ', levelDesc( iLevel )%depFromFiner( iElem )%bitmask(:)


      end do ! elem%nElems( eT_ghostFromFiner )

      write(dbgUnit(4),*) 'SourceIDs on targetlevel:', iLevel
      do iSourceElem = 1, levelDesc(iLevel)%sourceFromFiner%nVals
        sourceElem = levelDesc(iLevel)%sourceFromFiner%val(iSourceElem)
        write(dbgUnit(4),*) iSourceElem, sourceElem, levelDesc(iLevel+1)%total(sourceElem)
      end do

    end do ! iLevel = minLevel, maxLevel-1

    write(dbgUnit(4),"(A)") "Outside of mus_intp_update_depFromFiner routine"
    call tem_horizontalSpacer( fUnit = dbgUnit(4), after = 1 )

  end subroutine mus_intp_update_depFromFiner
! ****************************************************************************** !


end module mus_interpolate_module
! ****************************************************************************** !
!> \page intp_methods Interpolation methods
!!
!! \section intp_average Averaging interpolation
!!\verbatim interpolation_method = 'linear'
!!\endverbatim
!! Simple averaging of all quantities with weighted summation
!! \section intp_linear Linear interpolation
!! Linear interpolation of the first two moments (density and momentum) and
!! separation of equilibrium and non-equilibrium part to account for correct
!! shear stress calculation
!! This fine->coarse interpolation method is based on Dazhi Yu et. al. (2002):
!! "A multi-block lattice Boltzmann method for viscous fluid flows".
!! The basic idea is to split up the source pdfs into equilibrium- and
!! non-equilibrum- parts, then interpolate both parts with average interpolation,
!! and apply a corrective scaling factor to the non-equilibrium part before
!! adding both parts up and writing it to the target element.
!! Values are interpolated in a linear fashion by multiplying the quantities
!! of the source elements by weighting factors
