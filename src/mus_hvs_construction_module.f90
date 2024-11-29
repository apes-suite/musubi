! See copyright notice in the COPYRIGHT file.
!! ****************************************************************************** !
!> author: Kannan Masilamani
!! [[mus_hvs_construct]] "Creation of the data structures"
!! from the information in the configuration and
!! from the mesh read from disk for the musubi harvesting
!!
module mus_hvs_construction_module

  ! include treelm modules
  use mpi
  use env_module,              only: rk, long_k, eps, labelLen, my_status_int
  use treelmesh_module,        only: treelmesh_type
  use tem_comm_module,         only: tem_communication_type,                   &
    &                                tem_commpattern_type
  use tem_tools_module,        only: tem_horizontalSpacer
  use tem_grow_array_module,   only: grw_intArray_type, init, append, destroy, &
    &                                empty, truncate
  use tem_dyn_array_module,    only: dyn_intArray_type, init, append, expand,  &
    &                                destroy, empty, PositionOfVal, truncate,  &
    &                                dyn_labelArray_type
  use tem_construction_module, only: tem_init_elemLevels, tem_find_allElements,&
    &                                tem_build_horizontalDependencies,         &
    &                                tem_debug_horizontalDependencies,         &
    &                                tem_levelDesc_type,                       &
    &                                tem_dumpTreeIDLists
  use tem_element_module,      only: eT_fluid, eT_halo, eT_ghostFromCoarser,   &
    &                                eT_ghostFromFiner, destroy
  use tem_timer_module,        only: tem_startTimer, tem_stopTimer,            &
    &                                tem_getTimerVal
  use tem_debug_module,        only: main_debug, dbgUnit
  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit

  ! include musubi modules
  use mus_param_module,              only: mus_param_type
  use mus_geom_module,               only: mus_geom_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type,             &
    &                                      mus_finalize_layout
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_pdf_module,                only: mus_calc_nElems,                    &
    &                                      mus_pdf_allocate
  use mus_construction_module,       only: calculate_nElems
  use mus_timer_module,              only: mus_timerHandles


  implicit none

  private

  public :: mus_hvs_construct

contains

! ****************************************************************************** !
  !> Initialize Musubi data strucutres based on data provided by Treelm
  !!
  !! Load the mesh and boundary conditions for this process from disk.
  !! Get the level-wise treeID lists and create the required ghost and halo
  !! elements.
  !!
  !! This is achieved by a two-folded identification of elements.
  !!
  !! -# the theoretically required elements are collected based on
  !! [[mus_scheme_layout_module:mus_scheme_layout_type]] stencil information
  !! The [[tem_construction_module:tem_find_allElements]] routine
  !! performs this task for compute fluid elements.
  !! For boundaries which require information from neighbor elements, these
  !! required [[tem_topology_module]] "treeIDs" are collected into the
  !! [[mus_bc_header_module]]
  !! "boundary element type"
  !!
  !! # Additional Tasks
  !!
  !! - receive [[tem_construction_module:tem_build_horizontaldependencies]]
  !! "horizontal"
  !! (within a level for the element updates)
  !! - and [[tem_construction_module:tem_build_verticaldependencies]] "vertical"
  !! dependencies (between levels for ghost-interpolations).
  !! - The main state vector and the neighbor lists on which the kernel then
  !!    acts is created
  !! - The MPI buffers are created.
  !! - For each [[mus_scheme_module]] "Scheme", the
  !! [[tem_construction_module:tem_levelDesc_type]] "Level Descriptor"
  !! is created
  !!
  !! # Result
  !!
  !! After this routine, all data structures for starting the main loop of the
  !! solver are allocated and ready.
  !!
  !! Only difference between this routine and mus_construct is then
  !! creating of boundary elements and its stencil are omitted for harvesting
  !!
  subroutine mus_hvs_construct( scheme, geometry, params)
    ! ---------------------------------------------------------------------------
    !> run-time Parameters
    type( mus_param_type ), intent(inout) :: params
    !> geometric information
    type( mus_geom_type ), intent(inout) :: geometry
    !> scheme information including fluid, boundary and flow information
    type( mus_scheme_type ), intent(inout) :: scheme
    ! ---------------------------------------------------------------------------
    integer :: iLevel      ! counter for level
    integer :: minLevel, maxLevel
    integer :: hwmVal
    integer :: ii
    ! ---------------------------------------------------------------------------
    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(1),*) 'Get into routine: mus_construction'

    call tem_startTimer( timerHandle =  mus_timerHandles%initLvlD )

    ! set up the data structures for the kernel

    write(logUnit(1),*) 'Starting to initialize the geometry'
    call tem_horizontalSpacer(fUnit = logUnit(1))

    ! finalize scheme layout
    ! copy growing array of stencil to allocatable array of stencil
    ! and destroy growing array
    call mus_finalize_layout( layout   = scheme%layout,            &
      &                       nElemsInTree = geometry%tree%nElems, &
      &                       minLevel = minLevel,                 &
      &                       maxLevel = maxLevel,                 &
      &                       proc     = params%general%proc       )
    hwmVal = my_status_int('VmHWM:')
    write(logUnit(10),"(A,I0)") 'After finalize layout, VmHWM: ', hwmVal

    ! 1. store treeIDs which use this stencil in the
    !    levelDesc( elemLevel )%elem
    ! 2. store the neighbors which can be loaded from the disc or accessed
    !    directly
    ! THIS MEANS THAT FOR ALL ELEMENTS (INCL. BOUNDARY ELEMENTS) NEIGHBORS
    ! ARE STORED BUT THEY DO NOT NEED TO EXIST IN THE MESH AT THIS POINT!!!
    call tem_init_elemLevels(                            &
      &                me       = scheme%levelDesc,      &
      &                boundary = geometry%boundary,     &
      &                tree     = geometry%tree,         &
      &                stencils = [scheme%layout%fStencil] )
    hwmVal = my_status_int('VmHWM:')
    write(logUnit(10),"(A,I0)") 'After init elemLevels, VmHWM: ', hwmVal

    ! 1. build all dependencies for fluid, halos and ghosts
    ! 2. build the pointers for each element to its neighbors/stencil elements
    ! All this information is stored in tem_levelDesc_type.
    !
    ! THIS MEANS THAT FOR ALL ELEMENTS FOUND IN em_init_elemLevels THE CORRECT
    ! DEPENDENCIES AND NEIGHBORS ARE SET
    !
    ! This is the heart of the topologic data structure creation and
    ! takes most of the compute time, depending on the mesh you have
    !
    write(logUnit(1),*) 'Creating level descriptor ...'
    write(logUnit(6),*) 'before tem_find_allElements'
    write(dbgUnit(6),*) 'before find_allElements'
    ! The level descriptor is created from the below routine.
    ! the neigh array is created using the LD and communication buffers are filled up.
    call tem_find_allElements(                                         &
      &                  tree            = geometry%tree,              &
      &                  levelDesc       = scheme%levelDesc,           &
      &                  levelPointer    = geometry%levelPointer,      &
      &                  computeStencil  = [scheme%layout%fStencil],   &
      &                  commPattern     = params%general%commPattern, &
      &                  cleanup         = .true.,                     &
      &                  reqNesting      = params%nNesting,            &
      &                  proc            = params%general%proc         )
    hwmVal = my_status_int('VmHWM:')
    write(logUnit(10),"(A,I0)") 'After find allElements, VmHWM: ', hwmVal
    if ( main_debug%dumpTreeIDlists ) then
      call tem_dumpTreeIDlists( minLevel, maxLevel, &
        &                       scheme%levelDesc )
    end if

    write(dbgUnit(6),*) 'before horizontal dep'
    write(logUnit(6),*) 'before tem_build_horizontalDependencies'
    call tem_build_horizontalDependencies(            &
      &       iStencil       = 1,                     &
      &       levelDesc      = scheme%levelDesc,      &
      &       tree           = geometry%tree,         &
      &       computeStencil = scheme%layout%fStencil )
    if( main_debug%checkDependencies ) then
      call tem_debug_HorizontalDependencies( 1,         &
        &         scheme%levelDesc, geometry%tree, &
        &         scheme%layout%fStencil )
    end if
    hwmVal = my_status_int('VmHWM:')
    write(logUnit(10),"(A,I0)") 'After build horizontal, VmHWM: ', hwmVal

    call calculate_nElems( levelDesc = scheme%levelDesc,    &
      &                    proc      = params%general%proc, &
      &                    minLevel  = minLevel,            &
      &                    maxLevel  = maxLevel             )

    ! allocate here since it needs to be deallocated by dynamic load
    ! balancing algorithm
    write(logUnit(6),"(A,I0,A,I0)") 'allocate PDF from Level ', minLevel, &
      &                             ' to ', maxLevel
    allocate( scheme%pdf( minLevel:maxLevel ) )
    allocate( scheme%state( minLevel:maxLevel ) )

    write(logUnit(4),*) 'Allocating PDF state and neighbor array'
    ! Allocate the PDF state array
    do iLevel = minLevel, maxLevel
      write(dbgUnit(1),*) 'calculate nElems on level: ', iLevel
      call mus_calc_nElems(                                              &
        & me                = scheme%pdf( iLevel ),                      &
        & nFluids           = scheme%levelDesc( iLevel )%elem            &
        &                                %nElems( eT_fluid ),            &
        & nGhostFromCoarser = scheme%levelDesc( iLevel )%elem            &
        &                                %nElems( eT_ghostFromCoarser ), &
        & nGhostFromFiner   = scheme%levelDesc( iLevel )%elem            &
        &                                %nElems( eT_ghostFromFiner ),   &
        & nHalos            = scheme%levelDesc( iLevel )%elem            &
        &                                %nElems( eT_halo )              )

      allocate(scheme%state(iLevel)%val( scheme%pdf(iLevel)%nSize    &
        &                                * scheme%varSys%nScalars, 2 ) )
      !For debugging purposes, set complete flow field to invalid
      do ii = 1, scheme%pdf(iLevel)%nSize * scheme%varSys%nScalars
        scheme%state(iLevel)%val(ii, 1) = -1000000.0_rk
        scheme%state(iLevel)%val(ii, 2) = -1000000.0_rk
      end do

      call mus_pdf_allocate( me              = scheme%pdf( iLevel ),      &
        &                    nScalars        = scheme%varSys%nScalars,    &
        &                    QQ              = scheme%layout%fStencil%QQ, &
        &                    nElems_bcBuffer = 0,                         &
        &                    isPDF           = scheme%readVarIsPdf        )
    end do

    hwmVal = my_status_int('VmHWM:')
    write(logUnit(10),"(A,I0)") 'After allocate PDF, VmHWM: ', hwmVal

    write(dbgUnit(6),*) 'clean up elem list in level desctiptor'
    write(logUnit(6),*) 'clean up elem list in level descriptor'
    ! KJ: @todo de-activate for adaptive grid refinement
    do iLevel = minLevel, maxLevel
      call destroy( me = scheme%levelDesc( iLevel )%elem )
    end do

    write(dbgUnit(6),*) 'after clean up'

    call tem_stopTimer( timerHandle =  mus_timerHandles%initLvlD )

    call tem_horizontalSpacer( before = 1, fUnit = logUnit(1) )

    write(logUnit(1),"(A,F10.3)") 'Done initializing geometry. Took [s]', &
      &          tem_getTimerVal( timerHandle =  mus_timerHandles%initLvlD)
    call tem_horizontalSpacer( after  = 1, fUnit = logUnit(1) )

    write(dbgUnit(1),"(A)")  'after the construction'

  end subroutine mus_hvs_construct
! ****************************************************************************** !

end module mus_hvs_construction_module
! ****************************************************************************** !
