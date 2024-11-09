! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2016, 2018-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2015, 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2021-2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
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
!! Module containing subroutines for initialize Musubi source
!! variables and update source terms
!!
module mus_source_module
  ! include treelm modules
  use mpi
  use env_module,               only: rk, long_k, long_k_mpi
  use tem_aux_module,           only: tem_abort
  use tem_varSys_module,        only: tem_varSys_type
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_logging_module,       only: logUnit
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_geometry_module,      only: tem_baryOfID
  use tem_timer_module,         only: tem_startTimer, tem_stopTimer
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_debug_module,         only: dbgUnit
  use tem_subTree_module,       only: tem_create_subTree_of
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_stencil_module,       only: tem_stencilHeader_type

  ! include musubi modules
  use mus_derVarPos_module,      only: mus_derVarPos_type
  use mus_pdf_module,            only: pdf_data_type
  use mus_field_module,          only: mus_field_type
  use mus_source_type_module,    only: mus_source_type, mus_source_op_type, &
    &                                  mus_turbChannelForce_type,           &
    &                                  mus_HRRCorrectionTerm_type
  use mus_physics_module,        only: mus_convertFac_type
  use mus_absorbLayer_module,    only: mus_init_absorbLayer
  use mus_timer_module,          only: mus_timerHandles

  implicit none
  private

  public :: mus_init_sourceTerms
  public :: mus_apply_sourceTerms
  public :: mus_update_sourceVars

contains


  ! ************************************************************************** !
  !> This routine does set_params and setupIndices for all sources terms
  !! by gathering points to apply souce term before.
  subroutine mus_init_sourceTerms( field, nFields, globSrc, varSys, tree,    &
    &                              bc_prop, stencil, nElems_solve, levelDesc )
    ! --------------------------------------------------------------------------
    !> Number of fields
    integer, intent(in)                    :: nFields

    !> contains sources of all fields
    type(mus_field_type), intent(inout)   :: field(:)

    !> global source
    type(mus_source_type), intent(inout)   :: globSrc

    !> global variable system
    type(tem_varSys_type), intent(in)      :: varSys

    !> global treelm mesh
    type(treelmesh_type), intent(in)     :: tree

    !> bc property which is used to identify elements belong to certain BCs
    type( tem_bc_prop_type ), intent(in)       :: bc_prop

    !> Number of elements to solve in all levels
    !! nFluids + nGhostFromCoarser
    integer, intent(in) :: nElems_solve(tree%global%minLevel:)

    !> Level descriptors
    type( tem_levelDesc_type ), intent(in) :: levelDesc(tree%global%minLevel:)

    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    ! --------------------------------------------------------------------------
    integer :: iLevel, iField, iSrc
    integer :: nSolve, minLevel, maxLevel
    ! --------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*)'Initializing Source terms ...'
    write(dbgUnit(3),*)'Initializing Source terms ...'

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel

    ! allocate source element array
    do iField = 1, nFields
      do iSrc = 1, field(iField)%internalSource%varDict%nVals
        allocate( field(iField)%internalSource%method(iSrc)%elemLvl(minLevel:maxLevel) )
      end do
    end do

    ! get source parameter to store in data variable for coupling
    write(logUnit(10),*) ' Initializing internal sources'
    do iLevel = minLevel, maxLevel
      write(logUnit(10),*) '  iLevel: ', iLevel
      write(dbgUnit(3),*) '  iLevel: ', iLevel

      ! Use barycenter of all elements to solve i.e fluid+ghost to apply
      ! source terms
      nSolve = nElems_solve(iLevel)

      do iField = 1, nFields
        write(dbgUnit(3),*) '  field source: ', iField
        call mus_init_internalSource(                   &
          &  source     = field(iField)%internalSource, &
          &  varSys     = varSys,                       &
          &  nSolve     = nSolve,                       &
          &  iLevel     = iLevel,                       &
          &  stencil    = stencil                       )

      end do

    end do ! iLevel

    ! immediately exit this routine if there are no source terms to apply
    if ( all(field(:)%source%varDict%nVals == 0) &
      & .and. globSrc%varDict%nVals == 0 ) then
      write(logUnit(1),*) 'No active source terms'
      return
    end if

    ! allocate source element array
    do iField = 1, nFields
      do iSrc = 1, field(iField)%source%varDict%nVals
        allocate( field(iField)%source%method(iSrc)%elemLvl(minLevel:maxLevel) )
      end do
    end do

    do iSrc = 1, globSrc%varDict%nVals
      allocate( globSrc%method(iSrc)%elemLvl(minLevel:maxlevel) )
    end do

    ! get source parameter to store in data variable for coupling
    write(logUnit(10),*) ' Setup indices for sources'
    write(dbgUnit(3),*) ' Setup indices for source'
    do iLevel = minLevel, maxLevel
      write(logUnit(10),*) '  iLevel: ', iLevel
      write(dbgUnit(3),*) '  iLevel: ', iLevel


      ! Use barycenter of all elements to solve i.e fluid+ghost to apply
      ! source terms
      nSolve = nElems_solve(iLevel)

      do iField = 1, nFields
        write(dbgUnit(3),*) '  field source: ', iField
        ! Store idx for active field source variables
        call mus_setupIndices_forSrc( source     = field(iField)%source, &
          &                           varSys     = varSys,               &
          &                           tree       = tree,                 &
          &                           bc_prop    = bc_prop,              &
          &                           stencil    = stencil,              &
          &                           nSolve     = nSolve,               &
          &                           bary       = levelDesc(iLevel)     &
          &                                        %baryOfTotal,         &
          &                           iLevel     = iLevel                )
      end do

      write(dbgUnit(3),*) '  glob source: '
      ! Store idx for active glob source variables
      call mus_setupIndices_forSrc( source     = globSrc,          &
        &                           varSys     = varSys,           &
        &                           tree       = tree,             &
        &                           bc_prop    = bc_prop,          &
        &                           stencil    = stencil,          &
        &                           nSolve     = nSolve,           &
        &                           bary       = levelDesc(iLevel) &
        &                                        %baryOfTotal,     &
        &                           iLevel     = iLevel            )

    end do ! iLevel

    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_init_sourceTerms
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routines does setup indices for given source within a field or
  !! global. Index are stored for points which source term is active
  subroutine mus_setupIndices_forSrc(source, varSys, nSolve, bary, iLevel, &
    &                                tree, bc_prop, stencil)
    ! --------------------------------------------------------------------------
    !> Source term to fill in
    type(mus_source_type), intent(inout) :: source

    !> global variable system
    type(tem_varSys_type), intent(in)      :: varSys

    !> global treelm mesh
    type( treelmesh_type ), intent(in)     :: tree

    !> stencil used to find bcID on certain links
    type( tem_stencilHeader_type ), intent(in) :: stencil

    !> bc property which is used to identify elements belong to certain BCs
    type( tem_bc_prop_type ), intent(in)       :: bc_prop

    !> Number of elements to apply source term  on this level
    integer, intent(in) :: nSolve

    !> Space coordinates to apply source terms
    real(kind=rk), intent(in) :: bary(:,:)

    !> Current level
    integer, intent(in) :: iLevel
    ! --------------------------------------------------------------------------
    integer :: iVar, iElem, counter, src_nElems
    integer :: data_varPos
    integer, allocatable :: idx(:)
    integer(kind=long_k) :: nElems_var(source%varDict%nVals)
    ! --------------------------------------------------------------------------
    allocate(idx(nSolve))

    do iVar = 1, source%varDict%nVals
      idx = 0
      data_varPos = source%method(iVar)%data_varPos
      ! number of components of the source field
      ! nComps = varSys%method%val(data_varPos)%nComponents
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

      ! Store only valid idx to apply source term only on shapes defined
      ! for a data variable
      src_nElems = count(idx>0)
      nElems_var(iVar) = src_nElems
      source%method(iVar)%elemLvl(iLevel)%nElems = src_nElems
      allocate(source%method(iVar)%elemLvl(iLevel)%posInTotal(src_nElems))
      allocate(source%method(iVar)%elemLvl(iLevel)%idx(src_nElems))
!KM!      allocate(source%method(iVar)%elemLvl(iLevel)%val(src_nElems*nComps))

      select case(trim(source%varDict%val(iVar)%key))
      case('absorb_layer', 'absorb_layer_inlet', 'absorb_layer_outlet')
        call mus_init_absorbLayer( absLayer = source%method(iVar)%absLayer, &
          &                        dynAvg   = source%method(iVar)           &
          &                                   %elemLvl(iLevel)%dynAvg,      &
          &                        nElems   = src_nElems                    )
      case('turb_channel_force_accel')
        call mus_init_turbChanForce( turbChanForce = source%method(iVar)   &
          &                                                %turbChanForce, &
          &                          tree          = tree,                 &
          &                          bc_prop       = bc_prop,              &
          &                          stencil       = stencil               )
      end select

      counter = 0
      do iElem = 1, nSolve
        if (idx(iElem) > 0) then
          counter = counter + 1
          source%method(iVar)%elemLvl(iLevel)%posInTotal( counter ) = iElem
          source%method(iVar)%elemLvl(iLevel)%idx( counter ) = idx(iElem)
        end if
      end do
    end do !iVar

    !KM!call MPI_Reduce(nElems_var, glob_nElems_var, source%varDict%nVals, &
    !KM!  &             long_k_mpi, mpi_sum,                              &
    !KM!  &             0, tree%global%comm, ierror                        )
    do iVar = 1, source%varDict%nVals
      data_varPos = source%method(iVar)%data_varPos
      write(dbgUnit(3),*) 'Source iVar: ', iVar, &
        & trim(varSys%varName%val(data_varPos))
      write(dbgUnit(3),*) 'Total source nElems: ', nElems_var(iVar)
      write(logUnit(10),*) 'Source iVar: ', iVar, &
        & trim(varSys%varName%val(data_varPos))
      write(logUnit(10),*) 'Total source nElems: ', nElems_var(iVar)
    end do


  end subroutine mus_setupIndices_forSrc
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routines does setup indices for given source within a field or
  !! global. Index are stored for points which source term is active
  subroutine mus_init_internalSource(source, varSys, nSolve, iLevel, stencil)
    ! --------------------------------------------------------------------------
    !> Source term to fill in
    type(mus_source_type), intent(inout) :: source

    !> global variable system
    type(tem_varSys_type), intent(in)      :: varSys

    !> Number of elements to apply source term  on this level
    integer, intent(in) :: nSolve

    !> Current level
    integer, intent(in) :: iLevel

    !> layout descriptor
    type(tem_stencilHeader_type), intent(in) :: stencil
    ! --------------------------------------------------------------------------
    integer :: iVar, src_nElems
    ! --------------------------------------------------------------------------
    do iVar = 1, source%varDict%nVals
      ! Internal sources are applied to all fluid + ghostFromCoarser
      src_nElems = nSolve
      source%method(iVar)%elemLvl(iLevel)%nElems = nSolve

      select case(trim(source%varDict%val(iVar)%key))
      case('hrr_correction')
        call mus_init_hrrCorrection( HRR_Corr = source%method(iVar)           &
          &                                        %elemLvl(iLevel)%HRR_Corr, &
          &                          nElems   = src_nElems,                   &
          &                          nDim     = stencil%nDims                 )
      end select
    end do !iVar

  end subroutine mus_init_internalSource
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Create subTree and store nElemsGlobal in all proc for turbulent
  !! channel force
  subroutine mus_init_turbChanForce(turbChanForce, tree, bc_prop, stencil)
    ! --------------------------------------------------------------------------
    !> Contains info for turbulent channel force
    type(mus_turbChannelForce_type), intent(inout) :: turbChanForce
    !> global treelm mesh
    type( treelmesh_type ), intent(in)     :: tree
    !> bc property which is used to identify elements belong to certain BCs
    type( tem_bc_prop_type ), intent(in) :: bc_prop
    !> stencil used to find bcID on certain links
    type( tem_stencilHeader_type ), intent(in) :: stencil
    ! --------------------------------------------------------------------------
    integer :: iErr
    ! --------------------------------------------------------------------------
    call tem_create_subTree_of( inTree    = tree,                       &
      &                         bc_prop   = bc_prop,                    &
      &                         stencil   = stencil,                    &
      &                         subTree   = turbChanForce%subTree_utau, &
      &                         inShape   = turbChanForce%geom_utau     )

    call tem_create_subTree_of( inTree    = tree,                        &
      &                         bc_prop   = bc_prop,                     &
      &                         stencil   = stencil,                     &
      &                         subTree   = turbChanForce%subTree_umean, &
      &                         inShape   = turbChanForce%geom_umean     )


    ! Store global nElems in all proc to compute average velocity
    call mpi_allreduce( turbChanForce%subTree_utau%nElems,              &
      &                 turbChanForce%nElemsGlobal_utau,                &
      &                 1, mpi_integer, mpi_sum, tree%global%comm, iErr )

    call mpi_allreduce( turbChanForce%subTree_umean%nElems,             &
      &                 turbChanForce%nElemsGlobal_umean,               &
      &                 1, mpi_integer, mpi_sum, tree%global%comm, iErr )

    turbChanForce%forceDyn = 0.0_rk

  end subroutine mus_init_turbChanForce
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Initialize arrays to store time average density and velocity for
  !! dynamic hrrCorrection.
  !! \todo KM: 20210301 Allocate also for ghost cells!
  subroutine mus_init_hrrCorrection(HRR_Corr, nElems, nDim)
    ! --------------------------------------------------------------------------
    !> HRR correction term type
    type(mus_HRRCorrectionTerm_type), intent(inout) :: HRR_Corr
    !> Number of source elements
    integer, intent(in) :: nElems
    !> number of dimensions
    integer, intent(in) :: nDim
    ! --------------------------------------------------------------------------
    allocate(HRR_Corr%dens(nElems))
    allocate(HRR_Corr%vel(nElems,nDim))
    HRR_Corr%dens(:) = 0.0_rk
    HRR_Corr%vel(:,:) = 0.0_rk

  end subroutine mus_init_hrrCorrection
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Apply all source terms i.e field specific source and global source on
  !! all fields.
  !!
  subroutine mus_apply_sourceTerms( field, nFields, globSrc, pdf, varSys,  &
    &                               iLevel, time, phyConvFac, state,       &
    &                               auxField, derVarPos )
    ! --------------------------------------------------------------------------
    !> Number of fields
    integer, intent(in)                :: nFields

    !> contains sources of all fields
    type(mus_field_type), intent(in)  :: field(nFields)

    !> global source
    type(mus_source_type), intent(in)  :: globSrc

    !> pdf datatype
    type(pdf_data_type), intent(inout) :: pdf

    !> global variable system
    type(tem_varSys_type), intent(in)  :: varSys

    !> current level
    integer, intent(in)                :: iLevel

    !> current timing information
    type(tem_time_type), intent(in)    :: time

    !> state type containing the state vector to update
    real(kind=rk), intent(inout) :: state(:,:)

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! --------------------------------------------------------------------------
    ! counter variables
    integer :: iSrc, iField
    integer :: now, next
    ! --------------------------------------------------------------------------
    call tem_startTimer( timerHandle = mus_timerHandles%source )
    ! \todo create buffer for source to copy state values before applying
    ! source term to avoid using state which overwritten by source variable

    now = pdf%nNow
    next = pdf%nNext

    ! loop over all active sourve variable to apply source.
    ! Note: Internal field sources are added to state in compute routine
    do iField = 1, nFields
      do iSrc = 1, field(iField)%source%varDict%nVals
        call field(iField)%source%method(iSrc)%applySrc( &
          & inState    = state(:, now),                  &
          & outState   = state(:, next),                 &
          & neigh      = pdf%neigh(:),                   &
          & auxField   = auxField,                       &
          & nPdfSize   = pdf%nSize,                      &
          & iLevel     = iLevel,                         &
          & varSys     = varSys,                         &
          & time       = time,                           &
          & derVarPos  = derVarPos,                      &
          & phyConvFac = phyConvFac                      )
      end do
    end do

    ! apply global source
    do iSrc = 1, globSrc%varDict%nVals
      call globSrc%method(iSrc)%applySrc(  &
        & inState    = state(:, now),      &
        & outState   = state(:, next),     &
        & neigh      = pdf%neigh(:),       &
        & auxField   = auxField,           &
        & nPdfSize   = pdf%nSize,          &
        & iLevel     = iLevel,             &
        & varSys     = varSys,             &
        & time       = time,               &
        & derVarPos  = derVarPos,          &
        & phyConvFac = phyConvFac          )
    end do

    call tem_stopTimer( timerHandle = mus_timerHandles%source )
  end subroutine mus_apply_sourceTerms
! **************************************************************************** !

  ! ************************************************************************** !
  !> Updated all source variables i.e field specific source and global source on
  !! all fields.
  !!
  subroutine mus_update_sourceVars( nFields, field, globSrc, varSys, iLevel, &
    &                               auxField, phyConvFac, derVarPos )
    ! --------------------------------------------------------------------------
    !> Number of fields
    integer, intent(in) :: nFields

    !> contains sources of all fields
    type(mus_field_type), intent(inout) :: field(nFields)

    !> global source
    type(mus_source_type), intent(inout) :: globSrc

    !> global variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> current level
    integer, intent(in)                :: iLevel

    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)

    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac

    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! --------------------------------------------------------------------------
    ! counter variables
    integer :: iSrc, iField
    ! --------------------------------------------------------------------------
    call tem_startTimer( timerHandle = mus_timerHandles%source )

    ! update field source variables
    do iField = 1, nFields
      do iSrc = 1, field(iField)%source%varDict%nVals
        call field(iField)%source%method(iSrc)%updateSourceVar( &
          & auxField   = auxField,                              &
          & iLevel     = iLevel,                                &
          & varSys     = varSys,                                &
          & phyConvFac = phyConvFac,                            &
          & derVarPos  = derVarPos                              )
      end do
    end do

    ! update global source variables
    do iSrc = 1, globSrc%varDict%nVals
      call globSrc%method(iSrc)%updateSourceVar( &
        & auxField   = auxField,                 &
        & iLevel     = iLevel,                   &
        & varSys     = varSys,                   &
        & phyConvFac = phyConvFac,               &
        & derVarPos  = derVarPos                 )
    end do

    call tem_stopTimer( timerHandle = mus_timerHandles%source )
  end subroutine mus_update_sourceVars
  ! ************************************************************************** !

end module mus_source_module
! **************************************************************************** !
