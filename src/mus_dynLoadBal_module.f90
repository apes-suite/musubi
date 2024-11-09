! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2016, 2018-2022 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!! This routine serves as a wrapper to call other routines responsible for
!! dynamic load balancing in Musubi. The dynamic load balancing algorithm
!! follows the SPARTA algorithm in an iterative way to perform balancing for
!! Multilevel simulations as well as complex boundary conditions.
!!
module mus_dynLoadBal_module

  ! include treelm modules
  use mpi
  use env_module,              only: PathLen, labelLen, rk, long_k, rk_mpi
  use treelmesh_module,        only: treelmesh_type, &
    &                                exchange_elements,                       &
    &                                tem_dump_weights
  use tem_restart_module,      only: tem_init_restart,                        &
    &                                tem_restart_finalize
  use tem_general_module,      only: tem_general_type
  use tem_bc_prop_module,      only: tem_bc_prop_type
  use tem_topology_module,     only: tem_levelOf
  use tem_aux_module,          only: check_mpi_error , tem_abort
  use tem_tools_module,        only: tem_horizontalSpacer
  use tem_dyn_array_module,    only: append
  use tem_stencil_module,      only: append
  use tem_time_module,         only: tem_time_dump, tem_time_sim_stamp
  use tem_construction_module, only: tem_levelDesc_type
  use tem_sparta_module,       only: tem_balance_sparta, tem_sparta_type,     &
    &                                tem_exchange_sparta,                     &
    &                                tem_destroy_sparta,                      &
    &                                tem_derive_sparta, tem_init_sparta

  use tem_element_module,      only: eT_fluid
  use tem_logging_module,      only: logUnit
  use tem_debug_module,        only: dbgUnit
  use tem_convergence_module,  only: tem_init_convergence
  use tem_varMap_module,       only: tem_varMap_type
  use tem_property_module,     only: prp_hasBnd, prp_hasQVal
  use tem_spacetime_fun_module, only: tem_create_subTree_of_st_funList
  use tem_operation_var_module, only: tem_opVar_reduction_transient_init


  ! include musubi modules
  use mus_bc_general_module,         only: mus_init_boundary
  use mus_bc_header_module,          only: glob_boundary_type
  use mus_bndForce_module,           only: mus_init_bndForce
  use mus_param_module,              only: mus_param_type
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_scheme_layout_module,      only: mus_init_layout, mus_define_layout
  use mus_scheme_module,             only: mus_init_scheme, mus_scheme_cleanup
  use mus_geom_module,               only: mus_geom_type, mus_build_posInProp
  use mus_tools_module,              only: check_density, dump_linear_partition
  use mus_interpolate_module,        only: mus_init_interpolate
  use mus_source_module,             only: mus_init_sourceTerms
  use mus_transport_var_module,      only: mus_init_transport_var
  use mus_construction_module,       only: mus_construct
  use mus_flow_module,               only: fillHelperElementsCoarseToFine, &
    &                                      fillHelperElementsFineToCoarse
  use mus_IBM_module,                only: mus_unload_IBM, mus_reload_IBM, &
    &                                      mus_IBM_globType
  use mus_fluid_module,              only: mus_init_fluid
  use mus_gradData_module,           only: mus_init_gradData
  use mus_time_module,               only: mus_timeControl_homogenize
  use mus_varSys_module,             only: mus_varSys_solverData_type
  use mus_tracking_module,           only: mus_init_tracker
  use mus_timer_module,              only: mus_timerHandles,     &
    &                                      mus_reset_levelTimer, &
    &                                      mus_reset_BCTimer
  use mus_buffer_module,             only: mus_pdf_unserialize,  &
    &                                      mus_pdf_serialize
  use mus_weights_module,            only: mus_getWeights, mus_dumpWeights
  use mus_auxField_module,           only: mus_calcAuxFieldAndExchange,        &
    &                                      mus_intpAuxFieldCoarserAndExchange


  implicit none

  private

  public :: mus_perform_dynLoadBal


contains


! **************************************************************************** !
  !> Wrap up the routines required for dynamic load balancing
  subroutine mus_perform_dynLoadBal( scheme, params, geometry, solverData)
    ! --------------------------------------------------------------------------
    !> scheme type
    type( mus_scheme_type ), intent(inout) :: scheme
    !> Global parameters
    type( mus_param_type ), intent(inout) :: params
    !> Treelmesh data
    type( mus_geom_type ), intent(inout) :: geometry
    !> contains pointer to scheme, physics types
    type( mus_varSys_solverData_type ), target :: solverData
    ! --------------------------------------------------------------------------
    integer :: minLevel, maxLevel
    real(kind=rk), allocatable :: PDFbuffer(:)
    type( tem_sparta_type ) :: sparta
    ! --------------------------------------------------------------------------

    write(logUnit(1),"(A)") 'Performing dynamic load balance ...'

    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel

    ! Cleanup IBM
    write(dbgUnit(3),*) "Cleanup IBM"
    call mus_unload_IBM( me       = geometry%globIBM%IBM, &
      &                  proc     = params%general%proc,  &
      &                  minLevel = minLevel,             &
      &                  maxLevel = maxLevel              )

    ! deallocate( geometry%posInBndID )
    ! deallocate( geometry%posInQVal  )

    ! --------------------------------------------------------------------------
    ! Collect PDF from all levels
    write(logUnit(7),*) "Serialize PDF following the order of treeID list"
    allocate( PDFbuffer( geometry%tree%nElems * scheme%varSys%nScalars ) )
    call mus_pdf_serialize( scheme,                                        &
      &                     geometry%tree%treeID(1:geometry%tree%nElems),  &
      &                     geometry%levelPointer(1:geometry%tree%nElems), &
      &                     scheme%statevarMap, geometry%tree%nElems,      &
      &                     PDFbuffer, minLevel, maxLevel                  )
    ! --------------------------------------------------------------------------

    ! levelPointer not used anymore
    if( allocated( geometry%levelPointer ))then
      deallocate( geometry%levelPointer )
    end if

    ! --------------------------------------------------------------------------
    !         Balance the tree and free all arrays for reconstruction          !
    ! --------------------------------------------------------------------------
    ! Calculate weights for each element,
    ! Dump weight onto disk
    ! Call sparta to obtain new nElems and offset
    call tem_init_sparta( sparta, params%general%proc%comm_size )
    call mus_balance( levelDesc = scheme%levelDesc,      &
      ! &               globIBM   = scheme%field(1)%globIBM,  &
      &               general   = params%general,             &
      &               tree      = geometry%tree,              &
      &               nBCs      = geometry%boundary%nBCtypes, &
      &               globBC    = scheme%globBC,              &
      &               sparta    = sparta,                     &
      &               minLevel  = minLevel,                   &
      &               maxLevel  = maxLevel                    )

    call exchange_tree_bc( sparta, &
      &                    geometry, &
      &                    params%general%proc%comm, &
      &                    params%general%proc%comm_size )

    if ( params%dump_linear_partition ) then
      call dump_linear_partition( treeID = geometry%tree%treeID, &
        &                         nElems = geometry%tree%nElems, &
        &                         offset = geometry%tree%elemOffset, &
        &                         myRank = params%general%proc%rank, &
        &                         iter   = params%general%simControl%now%iter )
    end if
    ! Unload and reload tree and boundary ----------------------------

    ! Destroy the arrays which were allocated for performing new construction
    ! variables deallcoated include:
    !   globBC, levelDesc, pdf, layout%stencil
    call mus_scheme_cleanup( scheme, minLevel, maxLevel, &
      &                      geometry%boundary%nBCtypes  )

    ! --------------------------------------------------------------------------
    !           Reconstruct the levelDesc, reinitialize the simulation         !
    ! --------------------------------------------------------------------------
    !
    ! Re-initialize growing array of stencil as it is destroyed in
    ! mus_finialize_stencil and growing array stencil must be created
    ! for new list of boundary elements after load balancing
    ! Initialize layout growing array
    !call mus_init_layout( scheme%layout )
    ! define fStencil with predefined stencil layouts
    call mus_define_layout( layout      = scheme%layout,        &
      &                     stencilName = scheme%header%layout, &
      &                     nElems      = geometry%tree%nElems  )

    ! Initialize schemes: stencil, interpolation nSources and variable system
    call mus_init_scheme( me         = scheme,        &
      &                   tree       = geometry%tree, &
      &                   solverData = solverData     )

    call mus_construct( scheme    = scheme,   &
      &                 geometry  = geometry, &
      &                 params    = params )

    ! --------------------------------------------------------------------------
    !           exchange PDF and unserialize
    ! --------------------------------------------------------------------------
    ! exchange PDF
    write(logUnit(7),"(A)") "Going to exchange PDF"
    ! call tem_output_sparta( sparta, logUnit(1) )
    call tem_exchange_sparta( sparta, PDFbuffer, &
      &                       scheme%varSys%nScalars, params%general%proc%comm )
    call tem_destroy_sparta( sparta )
    ! copy PDF from buffer to PDF type
    write(logUnit(7),"(A)") "Copy serialized PDF back to level-wise array"
    call mus_pdf_unserialize( scheme,                                                       &
      &                       geometry%tree%treeID(1:geometry%tree%nElems),                 &
      &                       geometry%levelPointer(1:geometry%tree%nElems),                &
      &                       scheme%statevarMap, geometry%tree%nElems, PDFbuffer,          &
      &                       geometry%tree%global%minLevel, geometry%tree%global%maxLevel  )

    deallocate( PDFbuffer )
    ! --------------------------------------------------------------------------

    ! read Restart, reset BC, fill ghost by interpolation, tracking
    call mus_reset_aux(  scheme   = scheme,    &
      &                  geometry = geometry,  &
      &                  params   = params     )

    call mus_reload_IBM( me        = geometry%globIBM%IBM, &
      &                  iField    = 1,                    &
      &                  levelDesc = scheme%levelDesc,     &
      &                  tree      = geometry%tree         )

    write(logUnit(1),"(A)") 'Done with dynamic load balance.'
    write(logUnit(1),"(A)") ''

  end subroutine mus_perform_dynLoadBal
! **************************************************************************** !

! **************************************************************************** !
  !> This subroutine initializes musubi after a dynamic load balancing is
  !! performed.
  !!
  subroutine mus_reset_aux( scheme, params, geometry)
    ! --------------------------------------------------------------------------
    !> scheme type
    type( mus_scheme_type ), intent(inout) :: scheme
    !> Global parameters
    type( mus_param_type ),  intent(inout) :: params
    !> Treelmesh data
    type( mus_geom_type ),   intent(inout) :: geometry
    ! --------------------------------------------------------------------------
    integer :: minLevel, maxLevel, iLevel, ii
    ! real(kind=rk) :: total_density
    ! --------------------------------------------------------------------------

    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel

    !> initialize fluid type which contains relaxation parameter
    !! and function pointers to get mrt paramter and nonEqScaling factor
    !! for interpolation
    select case( trim(scheme%header%kind) )
    case('fluid', 'fluid_incompressible', 'isotherm_acEq')
      if (scheme%nFields > 1) then
        call tem_abort('chosen scheme kind supports only one field')
      end if
      ! initialize fluid viscosity relaxation paramters
      call mus_init_fluid(                                &
        & me           = scheme%field(1)%fieldProp%fluid, &
        & physics      = params%physics,                  &
        & schemeHeader = scheme%header,                   &
        & minLevel     = minLevel,                        &
        & maxLevel     = maxLevel,                        &
        & levelDesc    = scheme%levelDesc(:),             &
        & pdf          = scheme%pdf(:),                   &
        & stencil      = scheme%layout%fStencil,          &
        & general      = params%general,                  &
        & tNow         = params%general%simControl%now    )
    end select

    ! Initialize gradient data. Required for LES tuburbulent and evaluating
    ! gradient of a variable
    allocate(scheme%gradData(minLevel:maxLevel))
    do iLevel = minLevel, maxLevel
      call mus_init_gradData( me        = scheme%gradData(iLevel),         &
        &                     neigh     = scheme%pdf(iLevel)%neigh(:),     &
        !&                     levelDesc = scheme%levelDesc(iLevel),        &
        &                     stencil   = scheme%layout%fStencil,          &
        &                     nSize     = scheme%pdf(iLevel)%nSize,        &
        &                     nSolve    = scheme%pdf(iLevel)%nElems_solve, &
        &                     nScalars  = scheme%varSys%nScalars           )
    end do

    ! create subTree for all spacetime function in the linked list of
    ! spacetime function
    call tem_create_subTree_of_st_funList(     &
      &       me      = scheme%st_funList,     &
      &       tree    = geometry%tree,         &
      &       bc_prop = geometry%boundary,     &
      &       stencil = scheme%layout%fStencil )

    ! initialize the source terms for all fields and global source
    if ( all(scheme%field(:)%source%varDict%nVals /= 0) &
      & .or. scheme%globSrc%varDict%nVals /= 0 ) then
      write(logUnit(0),*) 'Error: In dynamic load balancing while reinitialize '
      write(logUnit(0),*) 'source terms. Source terms must be deallocated first'
      call tem_abort()
      call mus_init_sourceTerms( field        = scheme%field(:),            &
        &                        nFields      = scheme%nFields,             &
        &                        globSrc      = scheme%globSrc,             &
        &                        varSys       = scheme%varSys,              &
        &                        tree         = geometry%tree,              &
        &                        bc_prop      = geometry%boundary,          &
        &                        nElems_solve = scheme%pdf(:)%nElems_solve, &
        &                        levelDesc    = scheme%levelDesc,           &
        &                        stencil      = scheme%layout%fStencil      )
    end if

    ! initialize transport variables like velocity for passive scalar
    if ( scheme%transVar%varDict%nVals /= 0) then
      write(logUnit(0),*) 'Error: In dynamic load balancing while reinitialize '
      write(logUnit(0),*) 'transport var. It must be deallocated first!'
      call tem_abort()
      call mus_init_transport_var( me           = scheme%transVar,            &
        &                          varSys       = scheme%varSys,              &
        &                          tree         = geometry%tree,              &
        &                          nElems_solve = scheme%pdf(:)%nElems_solve, &
        &                          levelDesc    = scheme%levelDesc       )
    end if

    ! --------------------------------------------------------------------------
    !                      Reinitialize the global restart                     !
    ! --------------------------------------------------------------------------
    ! 1. finalize restart by freeing all mpi types
    call tem_restart_finalize( me = params%general%restart )
    ! 2. initialize restart again
    call tem_init_restart( me      = params%general%restart, &
      &                    solver  = params%general%solver,  &
      ! &                    varSys  = scheme%varSys,          &
      &                    varMap  = scheme%stateVarMap,     &
      &                    tree    = geometry%tree           )

    ! --------------------------------------------------------------------------
    !              Initialize tracking, interpolation and boundaries           !
    ! --------------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    !                  Reinitialize the tracking objects                     !
    ! ------------------------------------------------------------------------
    ! initialize tracking objects.
    call mus_init_tracker( scheme    = scheme,   &
      &                    geometry  = geometry, &
      &                    params    = params    )
    ! ------------------------------------------------------------------------
    !                  Reinitialize the tracking objects                     !
    ! ------------------------------------------------------------------------

    ! convergence objects
    if ( params%general%simControl%abortCriteria%steady_state ) then
      write(logUnit(1),"(A)") 'Initializing convergence...'
      write(dbgUnit(1),"(A)") 'init convergence'

      do ii = 1, size( params%general%simControl%abortCriteria%convergence)
        call mus_timeControl_homogenize(                                     &
          &       me = params%general%simControl%abortCriteria               &
          &                          %convergence(ii)%header%timeControl,    &
          &       dt = params%physics%dtLvl( maxLevel ),                     &
          &       reqInt = params%reqInterval                                )
      end do

      call tem_init_convergence( me       = params%general%simControl         &
        &                                   %abortCriteria%convergence,       &
        &                        tree     = geometry%tree,                    &
        &                        bc_prop  = geometry%boundary,                &
        &                        stencil  = scheme%layout%fStencil,           &
        &                        globProc = params%general%proc,              &
        &                        varSys   = scheme%varSys                     )
    end if

    if( minLevel /= maxlevel ) then
      write(logUnit(1),"(A)") 'Initializing interpolation...'
      ! initialize the interpolation
      call mus_init_interpolate(                               &
        &             intp         = scheme%intp,              &
        &             levelDesc    = scheme%levelDesc,         &
        &             schemeHeader = scheme%header,            &
        &             stencil      = scheme%layout%fStencil,   &
        &             minLevel     = minLevel,                 &
        &             maxLevel     = maxLevel,                 &
        &             fieldProp    = scheme%field(:)%fieldProp )
    end if

    ! ------------------------------------------------------------------------
    !                    Reinitialize the interpolation                      !
    ! ------------------------------------------------------------------------
    write(logUnit(1),"(A)") 'Fill ghost element by interpolation and do' &
      &                     //' communication'
    call fillHelperElementsFineToCoarse( scheme     = scheme,         &
      &                                  general    = params%general, &
      &                                  physics    = params%physics, &
      &                                  iLevel     = minLevel,       &
      &                                  maxLevel   = maxLevel        )

    call fillHelperElementsCoarseToFine( scheme     = scheme,         &
      &                                  general    = params%general, &
      &                                  physics    = params%physics, &
      &                                  iLevel     = minLevel,       &
      &                                  minLevel   = minLevel,       &
      &                                  maxLevel   = maxLevel        )

    ! reinitialize auxiliary field variable from state for fluid and ghost
    ! elements. Since set boundary is applied before load balancing, the
    ! auxField for fluid and ghostFromCoarser must be computed state
    ! using FETCH and ghostFromFiner must be interpolated
    do iLevel = minLevel, maxLevel
      call mus_calcAuxFieldAndExchange(                      &
        & auxField          = scheme%auxField(iLevel),       &
        & calcAuxField      = scheme%calcAuxField,           &
        & state             = scheme%state(iLevel)%val(:,    &
        &                       scheme%pdf(iLevel)%nNext),   &
        & pdfData           = scheme%pdf(iLevel),            &
        & nFields           = scheme%nFields,                &
        & field             = scheme%field(:),               &
        & globSrc           = scheme%globSrc,                &
        & stencil           = scheme%layout%fStencil,        &
        & varSys            = scheme%varSys,                 &
        & derVarPos         = scheme%derVarPos,              &
        & general           = params%general,                &
        & phyConvFac        = params%physics%fac(iLevel),    &
        & iLevel            = iLevel,                        &
        & minLevel          = geometry%tree%global%minLevel, &
        & schemeHeader      = scheme%header,                 &
        & quantities        = scheme%layout%quantities       )

      if (iLevel < maxLevel) then
        call mus_intpAuxFieldCoarserAndExchange(     &
          & intp        = scheme%intp,               &
          & tAuxField   = scheme%auxField(iLevel),   &
          & sAuxField   = scheme%auxField(iLevel+1), &
          & tLevelDesc  = scheme%levelDesc(iLevel),  &
          & stencil     = scheme%layout%fStencil,    &
          & iLevel      = iLevel,                    &
          & nAuxScalars = scheme%varSys%nAuxScalars, &
          & general     = params%general             )
       end if
    end do

    ! ------------------------------------------------------------------------
    !                    Reinitialize the boundaries                         !
    ! ------------------------------------------------------------------------
    ! Boundary force calculation is valid only for single field schemes
    ! like fluid and fluid_incompressible so initialize only if nFields = 1
    if (geometry%boundary%nBCtypes > 0 .and. scheme%nFields==1) then
      call mus_init_bndForce(bndForce     = geometry%bndForce,  &
        &                    bndMoment    = geometry%bndMoment, &
        &                    bc_prop      = geometry%boundary,  &
        &                    schemeHeader = scheme%header,      &
        &                    bc           = scheme%field(1)%bc  )
    end if
    write(logUnit(1),*) 'Initialize boundary routines'
    ! initialize each field boundary by looping over the field inside
    ! init_boundary_field
    call mus_init_boundary( field        = scheme%field,      &
      &                     pdf          = scheme%pdf,        &
      &                     state        = scheme%state,      &
      &                     auxField     = scheme%auxField,   &
      &                     tree         = geometry%tree,     &
      &                     levelDesc    = scheme%levelDesc,  &
      &                     layout       = scheme%layout,     &
      &                     schemeHeader = scheme%header,     &
      &                     varSys       = scheme%varSys,     &
      &                     derVarPos    = scheme%derVarPos,  &
      &                     globBC       = scheme%globBC,     &
      &                     bc_prop      = geometry%boundary  )

    call mus_reset_levelTimer()
    call mus_reset_bcTimer()

    ! Reinitialize time reduction operation variable and set last index
    ! with current value.
    call tem_opVar_reduction_transient_init(                &
      &      varSys         = scheme%varSys,                &
      &      tree           = geometry%tree,                &
      &      redTransVarMap = scheme%redTransVarMap,        &
      &      time           = params%general%simControl%now )

  end subroutine mus_reset_aux
! **************************************************************************** !


! **************************************************************************** !
  !> This routine performs the load balancing for multilevel simulations. The
  !! weights are calculated on the basis of levelwise run time, which are then
  !! fed to sparta for calculation of splitting positions. Restart files are
  !! saved and the simulation is restarted with the newly distributed mesh
  subroutine mus_balance( tree, minLevel, maxLevel, levelDesc, nBCs, globBC, &
    &                     general, sparta )
    ! --------------------------------------------------------------------------
    !> geometry infomation
    type(treelmesh_type),intent(inout) :: tree
    !> min level and max level
    integer, intent(in) :: minLevel, maxLevel
    !> Level descriptor
    type( tem_levelDesc_type ), intent(in) :: levelDesc(minLevel:maxLevel)
    !> global IBM type
    ! type( mus_IBM_globType ), intent(in) :: globIBM
    !> Number of boundary conditions
    integer, intent(in) :: nBCs
    !> BC elements information
    type( glob_boundary_type ), intent(in) :: globBC( nBCs )
    !> global parameters
    type( tem_general_type ), intent(in) :: general
    !> Sparta data type
    type( tem_sparta_type ), intent(inout) :: sparta
    ! --------------------------------------------------------------------------
    character(len=PathLen) :: basename
    character(len=labelLen) :: timestamp
    ! weights for fluid elements
    real(kind=rk), allocatable :: weights(:)
    ! --------------------------------------------------------------------------
    ! allocate weights
    allocate( weights( tree%nElems ) )
    ! Calculate weights according to compute, intp and bc routines
    call mus_getWeights( weights   = weights,   &
      &                  tree      = tree,      &
      &                  minLevel  = minLevel,  &
      &                  maxLevel  = maxLevel,  &
      &                  levelDesc = levelDesc, &
      &                  nBCs      = nBCs,      &
      &                  globBC    = globBC     )


    if ( general%balance%weight ) then
      write(logUnit(3),"(A)") "Dump weight file onto disk"
      ! Dump the weights using the current sim time in the basename
      ! Weight file name can be: simulation_name_weight_t
      timestamp = tem_time_sim_stamp( general%simControl%now )
      write(basename,'(a)')                                                 &
        &           './balance/'//trim(general%solver%simName)//'_weight_t' &
        &           //trim(timestamp)
      call mus_dumpWeights( tree     = tree,    &
        &                   weights  = weights, &
        &                   basename = basename )
    end if

    ! Call the sparta to get new myElems and offset and sparta
    call tem_balance_sparta( weight = weights,                 &
      &                      myPart = general%proc%rank,       &
      &                      nParts = general%proc%comm_size,  &
      &                      offset = tree%elemOffset,         &
      &                      myElems= tree%nElems,             &
      &                      sparta = sparta,                  &
      &                      comm   = general%proc%comm        )

    deallocate( weights )

  end subroutine mus_balance
! **************************************************************************** !


! **************************************************************************** !
  subroutine exchange_tree_bc( sparta, geometry, comm, comm_size )
    type( tem_sparta_type ), intent(in) :: sparta
    type( mus_geom_type ), intent(inout) :: geometry
    integer, intent(in) :: comm, comm_size

    type( tem_sparta_type ) :: sparta_bc, sparta_qVal

    ! tree is still OLD, but tree%nElems is already NEW
    if ( geometry%boundary%nBCtypes > 0 ) then
      write(logUnit(7),"(A)") "Set sparta info for BCID exchange"
      call tem_init_sparta( sparta_bc,   comm_size)
      call tem_derive_sparta(                                   &
        &    origin           = sparta,                         &
        &    derived          = sparta_bc,                      &
        &    nElems           = sparta%old_size,                &
        &    elemPropertyBits = geometry%tree%elemPropertyBits, &
        &    prpBit           = prp_hasbnd,                     &
        &    comm             = comm,                           &
        &    nParts           = comm_size                       )
      if ( any( geometry%boundary%hasQval ) ) then
        write(logUnit(7),"(A)") "Set sparta info for qValue exchange"
        call tem_init_sparta( sparta_qVal, comm_size)
        call tem_derive_sparta(                                   &
          &    origin           = sparta,                         &
          &    derived          = sparta_qVal,                    &
          &    nElems           = sparta%old_size,                &
          &    elemPropertyBits = geometry%tree%elemPropertyBits, &
          &    prpBit           = prp_hasqVal,                    &
          &    comm             = comm,                           &
          &    nParts           = comm_size                       )
      end if
    end if

    ! call tem_output_sparta( sparta, logUnit(1) )

    ! Unload and reload tree and boundary ----------------------------
    ! call unload_treelmesh( geometry%tree )
    ! call load_treelmesh( me      = geometry%tree, &
    !   &                  nParts  = params%general%proc%comm_size )
    call exchange_elements( geometry%tree, sparta )

    ! EXCHANGE BCID
    if ( geometry%boundary%nBCtypes > 0 ) then
      write(logUnit(7),"(A)") "Exchange boundaryID by sparta"
      call tem_exchange_sparta( sparta_bc, &
        &                       geometry%boundary%boundary_ID, &
        &                       geometry%boundary%nSides, &
        &                       comm )
      call tem_destroy_sparta( sparta_bc )

      if ( any(geometry%boundary%hasQVal) ) then
        write(logUnit(7),"(A)") "Exchange qValue by sparta"
        call tem_exchange_sparta( sparta_qVal, &
          &                       geometry%boundary%qVal, &
          &                       geometry%boundary%nSides, &
          &                       comm )
        call tem_destroy_sparta( sparta_qVal )
      end if

    end if

    call mus_build_posInProp( geometry )

  end subroutine exchange_tree_bc
! **************************************************************************** !

end module mus_dynLoadBal_module
! **************************************************************************** !
