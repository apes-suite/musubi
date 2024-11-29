! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2012,2015-2016,2020-2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2015-2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!! ****************************************************************************** !
!> \ref mus_construct "Creation of the data structures"
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
! Make sure loglvl is defined to an integer value.
! Usually this should be defined on the command line with -Dloglvl=
! to some value.
module mus_construction_module

  ! include treelm modules
  use mpi
  use env_module,              only: rk, long_k, long_k_mpi, eps, labelLen, &
    &                                my_status_int
  use treelmesh_module,        only: treelmesh_type
  use tem_comm_module,         only: tem_communication_type,                   &
    &                                tem_commpattern_type
  use tem_param_module,        only: childPosition, tem_rotationMatrix, PI
  use tem_topology_module,     only: tem_LevelOf, tem_ChildNumber,             &
    &                                tem_FirstIdAtLevel, tem_IDofCoord,        &
    &                                tem_CoordOfId, tem_ParentOF
  use tem_geometry_module,     only: tem_determine_discreteVector,             &
    &                                tem_ElemSizeLevel
  use tem_tools_module,        only: tem_horizontalSpacer
  use tem_grow_array_module,   only: init, append, destroy, empty, truncate
  use tem_dyn_array_module,    only: init, append, expand,  &
    &                                destroy, empty, PositionOfVal, truncate,  &
    &                                dyn_labelArray_type, truncate
  use tem_property_module,     only: tem_property_type, prp_solid,             &
    &                                prp_hasQVal, prp_sendHalo, prp_hasBnd
  use tem_bc_prop_module,      only: tem_bc_prop_type, tem_debug_bc_prop
  use tem_construction_module, only: tem_init_elemLevels, tem_find_allElements,&
    &                                tem_build_verticalDependencies,           &
    &                                tem_build_horizontalDependencies,         &
    &                                tem_debug_horizontalDependencies,         &
    &                                tem_levelDesc_type,                       &
    &                                tem_treeIDinTotal,                        &
    &                                tem_dumpTreeIDLists
  use tem_element_module,      only: eT_halo, eT_ghostFromCoarser,      &
    &                                eT_minRelevant, eT_maxRelevant,    &
    &                                eT_ghostFromFiner, eT_fluid, destroy
  use tem_timer_module,        only: tem_startTimer, tem_stopTimer,            &
    &                                tem_getTimerVal
  use tem_stencil_module,      only: tem_stencilHeader_type,              &
    &                                tem_stencil_getLabelForcxDir,        &
    &                                d3q81_cxDir, d3q125_cxDir,           &
    &                                grw_stencilHeaderArray_type, append, &
    &                                destroy, truncate, init
  use tem_debug_module,        only: main_debug, dbgUnit
  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit, tem_toStr
  use tem_subTree_module,      only: tem_write_debugMesh
  use tem_comm_env_module,     only: tem_comm_env_type
  use tem_geometry_module,     only: tem_eligibleChildren

  ! include musubi modules
  use mus_bc_header_module,          only: boundary_type, glob_boundary_type,  &
    &                                      mus_init_bc_elems,                  &
    &                                      mus_set_posInNghElems,              &
    &                                      mus_alloc_fieldBC,                  &
    &                                      debug_glob_boundary_type
  use mus_param_module,              only: mus_param_type
  use mus_geom_module,               only: mus_geom_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type,             &
    &                                      mus_finalize_layout
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_pdf_module,                only: pdf_data_type,   &
    &                                      mus_calc_nElems, mus_pdf_allocate,  &
    &                                      allocate_momBuf
  use mus_statistics_module,         only: mus_statistics_type,                &
    &                                      mus_calc_commAmount
    use mus_debug_module,              only: debug_normals, &
      &                                      debug_dependenciesFromFiner,      &
      &                                      debug_dependenciesFromCoarser,    &
      &                                      debug_connectivity
  use mus_comm_module,               only: mus_init_longBuffers
  use mus_field_module,              only: mus_field_type, remove_solid_in_bc, &
    &                                      mus_check_allWall,                  &
    &                                      mus_field_getSymmetricBCs
  use mus_IBM_module,                only: mus_IBM_globType
  use mus_timer_module,              only: mus_timerHandles
  use mus_interpolate_module,        only: mus_init_levelDescIntpArrays,   &
    &                                      mus_intp_update_depFromCoarser, &
    &                                      mus_intp_update_depFromFiner,   &
    &                                      mus_dump_levelDescIntp_nElems
  use mus_connectivity_module,       only: mus_construct_connectivity, &
    &                                      mus_updateConnectivity_forSymmetricBC
  use mus_auxField_module,           only: mus_init_auxFieldArrays

  implicit none

  private

  public :: mus_construct
  public :: communicate_property
  public :: calculate_nElems
  public :: set_sendHaloBits

  !> Contains bitmask to reduce halo elements communciation links.
  !! This type is used only in init_levelBuffers
  type halo_commBitmask_type
    !> Bitmask of which direction to send,
    !! only directions of fluid neighbor need to be update
    !! Size: nDofs, nElems
    logical, allocatable :: bitmask(:,:)
  end type halo_commBitmask_type

  type logical_array_type
    logical, allocatable :: val(:)
  end type logical_array_type


contains


! **************************************************************************** !
  !> Initialize Musubi data strucutres based on data provided by Treelm
  !!
  !! Load the mesh and boundary conditions for this process from disk.
  !! Get the level-wise treeID lists and create the required ghost and halo
  !! elements.
  !!
  ! !! TODO: Rewrite
  ! !! This is achieved by a two-folded identification of elements.
  ! !! - the theoretically required elements are collected based on
  ! !! [[mus_scheme_layout_module:mus_scheme_layout_type]] "stencil information"
  ! !! The [[tem_construction_module:tem_findneighbor]] "find neighbor routine"
  ! !! performs this task for compute fluid elements.
  ! !! For boundaries which require information from neighbor elements, these
  ! !! required [[tem_topology_module]] "treeIDs" are collected into the
  ! !! [[mus_bc_header_module]] "boundary element type"
  ! !!
  ! !! - All required elements are created in the
  ! !! [[tem_construction_module:tem_createleveldescriptor]]
  ! !! "Level Descriptor creation routine"
  ! !!
  ! !! # Additional Tasks
  ! !!
  ! !! - receive [[tem_construction_module:tem_buildhorizontaldependencies]]
  ! !! "horizontal"
  ! !! (within a level for the element updates)
  ! !! - and [[tem_construction:tem_buildverticaldependencies]] "vertical"
  ! !! dependencies (between levels for ghost-interpolations).
  ! !! - The main state vector and the neighbor lists on which the kernel then
  ! !!    acts is created
  ! !! - The MPI buffers are created.
  ! !! - For each [[mus_scheme_module]] "Scheme", the
  ! !! [[tem_construction_module:tem_levelDesc_type]] "Level Descriptor"
  ! !! is created
  !!
  !! # Result
  !!
  !! After this routine, all data structures for starting the main loop of the
  !! solver are allocated and ready.
  !!
  subroutine mus_construct( scheme, geometry, params)
    ! --------------------------------------------------------------------------
    !> scheme information including fluid, boundary and flow information
    type( mus_scheme_type ), intent(inout) :: scheme
    !> geometric information
    type( mus_geom_type ), intent(inout) :: geometry
    !> run-time Parameters
    type( mus_param_type ), intent(in) :: params
    ! --------------------------------------------------------------------------
    integer :: iLevel      ! counter for level
    integer :: iStencil    ! lists with different stencils
    integer :: iBC         ! counter for different boundary conditions
    integer :: iField      ! counter for different fields
    logical :: requireAll_loc ! are all links required
    integer :: minLevel, maxLevel
    type( mus_statistics_type ) :: stat
    !> mark halo elements which communicate all its links
    type( logical_array_type ), allocatable :: haloRequired(:)
    integer :: symmetricBCs(geometry%boundary%nBCtypes)
    integer :: nSymBCs
    ! --------------------------------------------------------------------------
    minLevel = geometry%tree%global%minLevel
    maxLevel = geometry%tree%global%maxLevel

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(1),*) 'Get into routine: mus_construction'

    call tem_startTimer( timerHandle = mus_timerHandles%initLvlD )


    write(logUnit(1),*) 'Starting to initialize the geometry'
    call tem_horizontalSpacer(fUnit = logUnit(1))

    ! 1. determine the number of elements for each mesh boundary
    ! 2. copy the information to globbc
    ! 3. determine the normal vector and bitmasks for the individual
    !    boundaries
    ! 4. allocate the BC lists in globbc
    ! THIS ROUTINE DOES NOT MODIFY THE BC LISTS IN THE FIELDS!!!
    call build_BClists( globBC    = scheme%globBC,              &
      &                 tree      = geometry%tree,              &
      &                 bc_prop   = geometry%boundary,          &
      &                 minLevel  = minLevel,                   &
      &                 maxLevel  = maxLevel,                   &
      &                 layout    = scheme%layout,              &
      &                 field     = scheme%field(:),            &
      &                 comm      = params%general%proc%comm    )

    ! Treat the boundary lists: in addition to the normal stencil treatment
    ! one has to find all those elements along the element list in the stencil
    ! type, which are boundary elements and give them the special
    ! boundary treatment

    ! 1. define the stencil for the individual boundaries (in the fields)
    !    by setting useAll, requireNeighNeigh, QQ, QQN, nElems, cxDir
    ! 2. fill stencil%elem with the level wise element positions and
    !    stencil%elemLvl with the level wise positions in the primal mesh
    do iField = 1, scheme%nFields
      do iBC = 1, geometry%boundary%nBCtypes
        call mus_alloc_fieldBC( scheme%field( iField )%bc( iBC ),&
          &                     minLevel, maxLevel)
        call mus_build_BCStencils(                                  &
          &      globbc         = scheme%globbc( iBC ),             &
          &      bc             = scheme%field( iField )%bc( iBC ), &
          &      prevailDir     = scheme%layout%prevailDir,         &
          &      prefix         = scheme%field( iField )%label,     &
          &      stencil_labels = scheme%layout%stencil_labels,     &
          &      grwStencil     = scheme%layout%grwStencil,         &
          &      minLevel       = minLevel,                         &
          &      maxLevel       = maxLevel                          )
      enddo ! iBC
    enddo ! iField

    ! Build and append IBM stencils to stencil array
    call mus_build_IBMStencils( globIBM    = geometry%globIBM,         &
      &                         layout     = scheme%layout,            &
      &                         grwStencil = scheme%layout%grwStencil  )

    ! finalize scheme layout
    ! copy growing array of stencil to allocatable array of stencil
    ! and destroy growing array
    call mus_finalize_layout( layout       = scheme%layout,        &
      &                       nElemsInTree = geometry%tree%nElems, &
      &                       minLevel     = minLevel,             &
      &                       maxLevel     = maxLevel,             &
      &                       proc         = params%general%proc   )
    ! peakval = my_status_int('VmPeak:')
    ! write(logUnit(10),"(A,I0)") 'After finalize layout, VmPeak: ', peakval

    ! 1. store treeIDs which use this stencil in the
    !    levelDesc( elemLevel )%elem
    ! 2. store the neighbors which can be loaded from the disc or accessed
    !    directly
    ! THIS MEANS THAT FOR ALL ELEMENTS (INCL. BOUNDARY ELEMENTS) NEIGHBORS
    ! ARE STORED BUT THEY DO NOT NEED TO EXIST IN THE MESH AT THIS POINT!!!
    call tem_init_elemLevels( me       = scheme%levelDesc,   &
      &                       boundary = geometry%boundary,       &
      &                       tree     = geometry%tree,           &
      &                       stencils = scheme%layout%stencil(:) )
    ! peakval = my_status_int('VmPeak:')
    ! write(logUnit(10),"(A,I0)") 'After init elemLevels, VmPeak: ', peakval

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
    !\todo KM: increasing nesting level by 1 if turbulence is active and viscosity
    ! is computed from velocity gradient, valid velocity is required in the
    ! two layers of ghost elements

    ! The level descriptor is created from the below routine.
    ! the neigh array is created using the LD and communication buffers are filled up.
    ! to properly evaluate gradients we need levels of nesting + 1 for fine ghosts
    call tem_find_allElements(                                         &
      &                  tree            = geometry%tree,              &
      &                  levelDesc       = scheme%levelDesc,           &
      &                  levelPointer    = geometry%levelPointer,      &
      &                  computeStencil  = scheme%layout%stencil(:),   &
      &                  commPattern     = params%general%commPattern, &
      &                  cleanup         = .true.,                     &
      &                  reqNesting      = params%nNesting+1,          &
      &                  proc            = params%general%proc         )
    ! peakval = my_status_int('VmPeak:')
    ! write(logUnit(10),"(A,I0)") 'After find allElements, VmPeak: ', peakval

      if ( main_debug%dumpTreeIDlists ) then
        call tem_dumpTreeIDlists( minLevel, maxLevel, &
          &                       scheme%levelDesc )
      end if

    write(dbgUnit(6),*) 'before horizontal dep'
    write(logUnit(6),*) 'before tem_build_horizontalDependencies'
    do iStencil = 1, scheme%layout%nStencils
      call tem_build_horizontalDependencies(                      &
        &       iStencil       = iStencil,                        &
        &       levelDesc      = scheme%levelDesc,                &
        &       tree           = geometry%tree,                   &
        &       computeStencil = scheme%layout%stencil(iStencil)  )
      if( main_debug%checkDependencies ) then
        call tem_debug_HorizontalDependencies( iStencil, &
          &         scheme%levelDesc, geometry%tree, &
          &         scheme%layout%stencil(iStencil) )
      end if
    end do !i Stencil
    ! peakval = my_status_int('VmPeak:')
    ! write(logUnit(10),"(A,I0)") 'After build horizontal, VmPeak: ', peakval

    ! Communicate property after setting missing Neighbor in
    ! tem_build_horizontalDependencies routine

    write(logUnit(6),*) 'before communicate_property'
    do iLevel = minLevel, maxLevel
      call communicate_property(                                  &
        &   send     = scheme%levelDesc( iLevel )%sendbuffer,&
        &   recv     = scheme%levelDesc( iLevel )%recvbuffer,&
        &   property = scheme%levelDesc( iLevel )%property,  &
        &   flag     = iLevel,                                    &
        &   proc     = params%general%proc,                       &
        &   pattern  = params%general%commPattern                 )
      call communicate_property(                                  &
        &   send     = scheme%levelDesc( iLevel )%sendbufferFromCoarser,  &
        &   recv     = scheme%levelDesc( iLevel )%recvbufferFromCoarser,  &
        &   property = scheme%levelDesc( iLevel )%property,       &
        &   flag     = iLevel,                                    &
        &   proc     = params%general%proc,                       &
        &   pattern  = params%general%commPattern                 )
      call communicate_property(                                  &
        &   send     = scheme%levelDesc( iLevel )%sendbufferFromFiner,    &
        &   recv     = scheme%levelDesc( iLevel )%recvbufferFromFiner,    &
        &   property = scheme%levelDesc( iLevel )%property,       &
        &   flag     = iLevel,                                    &
        &   proc     = params%general%proc,                       &
        &   pattern  = params%general%commPattern                 )
    end do ! iLevel
    ! peakval = my_status_int('VmPeak:')
    ! write(logUnit(10),"(A,I0)") 'After comm property, VmPeak: ', peakval

    call tem_horizontalSpacer( after = 1, fUnit = logUnit(1) )

    ! Update the positions stored in the bc%elem lists to the positions in the
    ! levelwise total list
    ! 1. allocate and initialize the pre- and post collision buffers as well
    !    as the bcBuffer for the boundary elements
    ! 2. Update the positions stored in the bc%elems lists to the positions
    !    in the levelwise total list as well as the positions of the neighbors
    !    in bc%elems%neigh
    ! IN THIS ROUTINE THE ALREADY EXISTING NEIGHBORS ARE COPIED AND THE ONES
    ! WHICH DO NOT EXIST ARE SET TO REASONABLE ELEMENTS (CURRENT ELEMENT OR
    ! LAST VALID NEIGHBOR)
    if( main_debug%dumpBoundaries ) then
      write(dbgUnit(1),*) 'Check BC before update_BCLists, i.e. elemPos ', &
        &                 ' with respect to tree'
      call debug_glob_boundary_type( geometry%boundary%nBCtypes, minLevel, maxLevel, &
        &                            scheme%globBC )
    end if

    call update_BClists( scheme       = scheme,                     &
      &                  minLevel     = minLevel,                   &
      &                  maxLevel     = maxLevel,                   &
      &                  nBCs         = geometry%boundary%nBCtypes, &
      &                  LP           = geometry%levelPointer,      &
      &                  remove_solid = params%remove_solid         )

    ! Get vertical dependencies (between levels for interpolation)
    call tem_build_verticalDependencies(          &
      &        levelDesc = scheme%levelDesc, &
      &        minlevel  = minLevel,              &
      &        maxlevel  = maxLevel               )
    ! peakval = my_status_int('VmPeak:')
    ! write(logUnit(10),"(A,I0)") 'After build vertical, VmPeak: ', peakval

    ! --------------------------------------------------------------------------
    ! Only if multilevel
    if ( minLevel /= maxLevel ) then

      ! initialize the levelwise ghost and source list
      call mus_init_levelDescIntpArrays(                    &
        &  intp      = scheme%intp,                         &
        &  levelDesc = scheme%levelDesc(minLevel:maxLevel), &
        &  minLevel  = minLevel,                            &
        &  maxLevel  = maxLevel                             )

      ! Update source elements for fromCoarser ghost elements and compute
      ! weights and least square fit matric depending in interpolation method
      ! and number of source elements found
      call mus_intp_update_depFromCoarser(                  &
        &  intp      = scheme%intp,                         &
        &  levelDesc = scheme%levelDesc(minLevel:maxLevel), &
        &  stencil   = scheme%layout%fStencil,              &
        &  minLevel  = minLevel,                            &
        &  maxLevel  = maxLevel                             )

      ! Update the dependencies fromFiner (in 2d, 4 of 8 get removed)
      call mus_intp_update_depFromFiner(                    &
        &  intp      = scheme%intp,                         &
        &  levelDesc = scheme%levelDesc(minLevel:maxLevel), &
        &  minLevel  = minLevel,                            &
        &  maxLevel  = maxLevel,                            &
        &  stencil   = scheme%layout%fStencil              )

      ! dumps global nElems in intpFromCoarser, intpFromFiner,
      ! sourcesFromCoarser ans sourcesFromFiner
      call mus_dump_levelDescIntp_nElems(                   &
        &  intp      = scheme%intp,                         &
        &  levelDesc = scheme%levelDesc(minLevel:maxLevel), &
        &  minLevel  = minLevel,                            &
        &  maxLevel  = maxLevel,                            &
        &  root      = params%general%proc%root,            &
        &  comm      = params%general%proc%comm             )

    end if ! minLevel /= maxLevel
    ! --------------------------------------------------------------------------
    ! peakval = my_status_int('VmPeak:')
    ! write(logUnit(10),"(A,I0)") 'After update interpolation, VmPeak: ', peakval

    call mus_update_BcghostElem( minLevel  = minLevel,          &
      &                          maxLevel  = maxLevel,          &
      &                          levelDesc = scheme%levelDesc,  &
      &                          layout    = scheme%layout,     &
      &                          bc_prop   = geometry%boundary, &
      &                          globBC    = scheme%globBC      )

    call calculate_nElems( levelDesc = scheme%leveldesc,         &
      &                    proc      = params%general%proc,      &
      &                    minLevel  = minLevel,                 &
      &                    maxLevel  = maxLevel                  )

    call finalize_BClist( scheme       = scheme,                      &
      &                   proc         = params%general%proc,         &
      &                   minLevel     = minLevel,                    &
      &                   maxLevel     = maxLevel,                    &
      &                   nBCs         = geometry%boundary%nBCtypes   )

    ! Build boundary level pointer
    if (geometry%boundary%nBCtypes > 0) then
      call build_bcLevelPointer( bcLevelPointer = geometry%bcLevelPointer, &
        &                        posInBndID     = geometry%posInBndID,     &
        &                        minBcID        = geometry%minBcID,        &
        &                        levelPointer   = geometry%levelPointer,   &
        &                        bc_prop        = geometry%boundary,       &
        &                        tree           = geometry%tree,           &
        &                        globBC         = scheme%globBC            )
    end if

    ! allocate here since it needs to be deallocated by dynamic load
    ! balancing algorithm
    write(logUnit(6),"(A,I0,A,I0)") 'allocate PDF from Level ', minLevel, &
      &                             ' to ', maxLevel
    allocate( scheme%pdf( minLevel:maxLevel ) )
    allocate( scheme%state( minLevel:maxLevel ) )

    write(logUnit(4),*) 'Allocating PDF state, neighbor and bcBuffer array'
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
        &                              * scheme%varSys%nScalars, 2 ) )
      !For debugging purposes, set complete flow field to invalid
      scheme%state(iLevel)%val(:, :) = -1000000.0_rk

      ! allocate neigh array
      call mus_pdf_allocate(                                           &
        & me              = scheme%pdf( iLevel ),                      &
        & nScalars        = scheme%varSys%nScalars,                    &
        & QQ              = scheme%layout%fStencil%QQ,                 &
        & isPDF           = .true.,                                    &
        & nElems_bcBuffer = scheme%levelDesc( iLevel )                 &
        &                                         %bc_elemBuffer%nVals )

!KM \todo momBUF: activate while pre-computing moments for interpolation
!KM!      call allocate_momBuf( scheme%pdf(iLevel),                        &
!KM!        &    max( scheme%levelDesc(iLevel)%sourceFromCoarser%nVals,    &
!KM!        &         scheme%levelDesc(iLevel)%sourceFromFiner%nVals    )  )

    end do
    ! peakval = my_status_int('VmPeak:')
    ! write(logUnit(10),"(A,I0)") 'After allocate PDF, VmPeak: ', peakval

    write(dbgUnit(6),*) 'clean up elem list in level desctiptor'
    write(logUnit(6),*) 'clean up elem list in level descriptor'
    ! KJ: @todo de-activate for adaptive grid refinement
    do iLevel = minLevel, maxLevel
      call destroy( me = scheme%levelDesc( iLevel )%elem )
    end do

    write(dbgUnit(6),*) 'after clean up'

    do iLevel = minLevel, maxLevel
      ! construct connectivity vector pdf( iLevel )%neigh including
      ! bounce back rules
      call mus_construct_connectivity(                   &
        & neigh       = scheme%pdf(iLevel)%neigh,        &
        & nSize       = scheme%pdf(iLevel)%nSize,        &
        & nElems      = scheme%pdf(iLevel)%nElems_local, &
        & levelDesc   = scheme%levelDesc(iLevel),        &
        & stencil     = scheme%layout%fStencil,          &
        & varSys      = scheme%varSys,                   &
        & stateVarMap = scheme%stateVarMap               )

      ! identify boundary ids which are symmetry kind
      call mus_field_getSymmetricBCs(                &
        & symmetricBCs = symmetricBCs,               &
        & nSymBCs      = nSymBCs,                    &
        & nBCs         = geometry%boundary%nBCtypes, &
        & nFields      = scheme%nFields,             &
        & field        = scheme%field                )

      ! Symmetric boundaries are defined. Update connectivity for those BC
      if ( nSymBCS > 0 ) then
        ! update connectivity to treat symmetric boundary condition implicitly
        call mus_updateConnectivity_forSymmetricBC(    &
          & neigh        = scheme%pdf(iLevel)%neigh,   &
          & nSize        = scheme%pdf(iLevel)%nSize,   &
          & iLevel       = iLevel,                     &
          & levelDesc    = scheme%levelDesc(iLevel),   &
          & layout       = scheme%layout,              &
          & varSys       = scheme%varSys,              &
          & stateVarMap  = scheme%stateVarMap,         &
          & nBCs         = geometry%boundary%nBCtypes, &
          & globBC       = scheme%globBC,              &
          & nSymBCs      = nSymBCs,                    &
          & symmetricBCs = symmetricBCs(1:nSymBCs)     )
      end if
    end do


    ! Dump state neighbor array information to debugUnit
    if( main_debug%dumpAuxLists ) then
      call debug_connectivity( scheme   = scheme,   &
        &                      minLevel = minLevel, &
        &                      maxLevel = maxLevel  )
    end if

    ! Set halo elements which require all links for communication.
    call set_halo_commLinks( scheme, minLevel, maxLevel, &
      &                      geometry%boundary%nBCtypes, &
      &                      params%comm_reduced,        &
      &                      haloRequired )

    write(dbgUnit(6),*)  'before init levelBuffers'

    ! Communicate all links if any stencil has requireAll true or
    ! comm_reduced is false
    requireAll_loc = any( scheme%layout%stencil(:)%requireAll ) &
      &             .or. .not. params%comm_reduced
    do iLevel = minLevel, maxLevel
      !    for multi-level, still all the links have to be communicated
      !    so, check for regular Buffer, if the max and min level is the same
      write(dbgUnit(1),*) "Level: ", iLevel
      write(dbgUnit(1),*) "Set send and recv buffer for normal Halo elements"
      call init_levelBuffers(                                 &
        & send       = scheme%levelDesc( iLevel )%sendbuffer, &
        & recv       = scheme%levelDesc( iLevel )%recvbuffer, &
        & pdf        = scheme%pdf(iLevel),                    &
        & requireAll = requireAll_loc,                        &
        & pattern    = params%general%commPattern,            &
        & comm       = params%general%proc%comm,              &
        & scheme     = scheme,                                &
        & stat       = stat,                                  &
        & haloRequired = haloRequired(iLevel)%val,             &
        & offset     = scheme%levelDesc( iLevel )%offset      )


      write(dbgUnit(1),*) "Set send and recv buffer for fromCoarser Halo elements"
      call init_levelBuffers(                                     &
        & send       = scheme%levelDesc( iLevel )%sendbufferFromCoarser, &
        & recv       = scheme%levelDesc( iLevel )%recvbufferFromCoarser, &
        & pdf        = scheme%pdf(iLevel),                        &
        & requireAll = requireAll_loc,                            &
        & pattern    = params%general%commPattern,                &
        & comm       = params%general%proc%comm,                   &
        & scheme     = scheme,                                    &
        & haloRequired = haloRequired(iLevel)%val,             &
        & stat       = stat,                            &
        & offset     = scheme%levelDesc( iLevel )%offset      )

      write(dbgUnit(1),*) "Set MPI buffer for fromFiner Halo elements"
      call init_levelBuffers(                                     &
        & send       = scheme%levelDesc( iLevel )%sendbufferFromFiner, &
        & recv       = scheme%levelDesc( iLevel )%recvbufferFromFiner, &
        & pdf        = scheme%pdf(iLevel),                        &
        & requireAll = requireAll_loc,                            &
        & pattern    = params%general%commPattern,                &
        & comm       = params%general%proc%comm,                  &
        & scheme     = scheme,                                    &
        & stat       = stat,                                      &
        & haloRequired = haloRequired(iLevel)%val,                &
        & offset     = scheme%levelDesc( iLevel )%offset          )

      deallocate( haloRequired(iLevel)%val )

    end do

    deallocate( haloRequired )

    ! Initialize auxField var val array and communication buffers
    allocate(scheme%auxField(minLevel:maxLevel))
    do iLevel = minLevel, maxLevel
      call mus_init_auxFieldArrays( me          = scheme%auxField(iLevel),    &
        &                           levelDesc   = scheme%levelDesc(iLevel),   &
        &                           pattern     = params%general%commPattern, &
        &                           nSize       = scheme%pdf(iLevel)%nSize,   &
        &                           nAuxScalars = scheme%varSys%nAuxScalars   )
    end do

    call tem_write_debugMesh( globtree  = geometry%tree,           &
      &                       levelDesc = scheme%levelDesc,   &
      &                       myPart    = params%general%proc%rank )

    ! Dump Ghost and Source elements information to debugUnit
    if( main_debug%dumpDependencies) then
      call debug_dependenciesFromFiner( &
        &    scheme%levelDesc(minLevel:maxLevel), &
        &    minLevel, maxLevel )
      call debug_dependenciesFromCoarser( &
        &    scheme%levelDesc(minLevel:maxLevel), &
        &    minLevel, maxLevel )
    end if

    ! Dump Boundary debug information to debugUnit
    if( main_debug%dumpBoundaries ) then
      if ( size( geometry%tree%Property ) > 0 ) then

        write(dbgUnit(1),*) 'Check BC AFTER update_BCLists, i.e. elemPos ', &
          &                 ' with respect to total list'
        call debug_glob_boundary_type( geometry%boundary%nBCtypes, minLevel,   &
          &                            maxLevel, scheme%globBC )

        call tem_debug_bc_prop( geometry%boundary )

        call debug_normals( tree   = geometry%tree, &
          &                 scheme = scheme,        &
          &                 nBCs   = geometry%boundary%nBCtypes )

      end if
    end if ! dump Boundaries

    do iLevel = minlevel, maxLevel
      ! set the sendHalo bit for all fluid elements which are send
      ! to remote procs
      call set_sendHaloBits( property = scheme%levelDesc( iLevel )%property, &
        &                    sendBuffer = scheme%levelDesc( iLevel )%sendBuffer)
    end do

    call mus_calc_commAmount( stat = stat,     &
      &                       comm = params%general%proc%comm, &
      &                       nprocs = params%general%proc%comm_size )

    call tem_stopTimer( timerHandle = mus_timerHandles%initLvlD )

    call tem_horizontalSpacer( before = 1, fUnit = logUnit(1) )

    write(logUnit(1),"(A,F10.2)") 'Done initializing geometry. Took [s]',  &
      &           tem_getTimerVal( timerHandle =  mus_timerHandles%initLvlD)
    call tem_horizontalSpacer( after  = 1, fUnit = logUnit(1) )

    write(dbgUnit(1),*)  'after the construction'

  end subroutine mus_construct
! **************************************************************************** !


! **************************************************************************** !
  !> Initialize the communication buffers for a single level
  subroutine communicate_property( send, recv, property, flag, proc, pattern)
    ! --------------------------------------------------------------------------
    !> Communication structure to send and receive property
    type(tem_communication_type), intent(inout) :: send, recv
    !> property to communicate
    integer(long_k), intent(inout)              :: property(:)
    !> communication flag
    integer, intent(in)                         :: flag
    !>
    type(tem_commPattern_type), intent(in)      :: pattern
    !> Process description to use.
    type(tem_comm_env_type), intent(in)         :: proc
    ! --------------------------------------------------------------------------

    call mus_init_longBuffers( comm    = send,                                 &
      &                        pattern = pattern )
    call mus_init_longBuffers( comm    = recv,                                 &
      &                        pattern = pattern )

    call pattern%exchange_long( send         = send,                           &
      &                         recv         = recv,                           &
      &                         state        = property,                       &
      &                         message_flag = flag,                           &
      &                         comm         = proc%comm )

  end subroutine communicate_property
! **************************************************************************** !


! **************************************************************************** !
  !> Initialize the communication buffers for a single level
  subroutine init_levelBuffers( send, recv, pdf, pattern, offset, &
    &                           requireAll, scheme, stat, comm, haloRequired )
    ! --------------------------------------------------------------------------
    !> Communication structure to initialize
    type(tem_communication_type)               :: send, recv
    !> iLevel pdf info with neigh array
    type( pdf_data_type ), intent(in)          :: pdf
    !> communication pattern
    type( tem_commpattern_type ), intent(in)   :: pattern
    !>
    integer, intent(in)         :: offset( 2, eT_minRelevant:eT_maxRelevant )
    !> Statistics
    type( mus_statistics_type ), intent(inout) :: stat
    !> mpi communication enviroment with mpi communicator
    integer, intent(in)                        :: comm
    !> fluid, bnd and flow info
    type( mus_scheme_type ), intent(in)        :: scheme
    !> different place to take values from for interpolation
    logical, intent(in)                        :: requireAll
    !>
    logical, intent(in)                        :: haloRequired(:)
    ! --------------------------------------------------------------------------
    ! myhalos are to be received from remote and remoteHalos are to be
    ! send to remote
    type( halo_commBitmask_type ), allocatable :: myHalos(:), remoteHalos(:)
    ! --------------------------------------------------------------------------

    ! initialize bitmask array
    call init_bitmask( Halos    = myHalos,               &
      &                comm     = recv,                  &
      &                nScalars = scheme%varSys%nScalars )
    call init_bitmask( Halos    = remoteHalos,           &
      &                comm     = send,                  &
      &                nScalars = scheme%varSys%nScalars )

    ! only initialize recvBuffers and then communicate the desired
    ! values to the neighbor processes.
    ! @todo: SZ: think about the nDOFs in case more state variables than
    !            the QQ PDFs exist
    call init_recvBuffers(me         = recv,                 &
      &                   myhalos    = myhalos,              &
      &                   neigh      = pdf%neigh,            &
      &                   nSize      = pdf%nSize,            &
      &                 haloRequired = haloRequired,         &
      &                   requireAll = requireAll,           &
      &                   halo_offset= offset(1,eT_halo),    &
      &                   stat       = stat,                 &
      &                   scheme     = scheme,               &
      &                   pattern    = pattern               )

    ! communicate bitmasks
    call communicate_bitmask()

    call init_sendBuffers(send        = send,                 &
      &                   remoteHalos = remoteHalos,          &
      &                   neigh       = pdf%neigh,            &
      &                   nSize       = pdf%nSize,            &
      &                   scheme      = scheme,               &
      &                   pattern     = pattern               )

    ! destroy send and recv bitmask since its not needed anymore
    call destroy_bitmask( remoteHalos, send%nProcs )
    call destroy_bitmask( myhalos, recv%nProcs )

    contains

    ! ************************************************************************ !
      !> Initialize myHalos and remoteHalos bitmasks.
      !! myHalos bitmask are from my halo elements which are to be
      !! received from remote process and remoteHalos bitmask are
      !! remote process halos which are requested by remote process.
      subroutine init_bitmask( halos, comm, nScalars )
        ! ---------------------------------------------------------------------
        !> contains bitmask of my halos with array size of recv%nProcs
        type(halo_commBitmask_type), allocatable, intent(out) :: halos(:)
        !> Communication structure to initialize
        type(tem_communication_type), intent(in)              :: comm
        !> number of scalars in variable system to communicate
        integer,                                   intent(in) :: nScalars
        ! ----------------------------------------------------------------------
        integer :: iProc
        ! ----------------------------------------------------------------------

        ! number of process my halos depends on
        allocate(halos(comm%nProcs))

        ! Allocate the bitmask for elements to be communicated
        ! from the provider proc
        do iProc = 1, comm%nProcs
          allocate( halos(iProc)%bitmask( nScalars, comm%nElemsProc(iProc) ) )
          halos(iProc)%bitmask = .false.
        end do

      end subroutine init_bitmask
    ! ************************************************************************ !

    ! ************************************************************************ !
      !> Destroy allocated send and recv bitmask buffer
      subroutine destroy_bitmask( halos, nProcs )
        ! ----------------------------------------------------------------------
        !> Communication structure to initialize
        integer, intent(in) :: nProcs
        !> contains bitmask of my halos with array size of recv%nProcs
        type(halo_commBitmask_type), allocatable, intent(inout) :: Halos(:)
        ! ----------------------------------------------------------------------
        integer :: iProc
        ! ----------------------------------------------------------------------

        do iProc = 1, nProcs
          deallocate( halos( iProc )%bitmask )
        end do
        deallocate( halos )

      end subroutine destroy_bitmask
    ! ************************************************************************ !

    ! ************************************************************************ !
      subroutine communicate_bitmask()
        ! ----------------------------------------------------------------------
        integer :: iProc
        ! ------ MPI variables
        ! send%nProcs + recv%nProcs )
        integer,allocatable :: rq_handle(:)
        ! mpi_status_size, recv%nProcs + send%nProcs )
        integer,allocatable :: status(:,:)
        integer  :: iErr
        integer  :: message_flag
        integer  :: nScalars
        ! ----------------------------------------------------------------------

        nScalars = scheme%varSys%nScalars
        ! Choosing random message flag
        message_flag = 1
        allocate( rq_handle( send%nProcs + recv%nProcs))
        allocate( status( mpi_status_size, send%nProcs + recv%nProcs))

        ! First start receive communication of the 'send' bitmasks from the target
        ! Procs. Consumer view: receive from the target Process
        do iProc = 1, send%nProcs
          call mpi_irecv(                              &
            &      remoteHalos( iProc )%bitMask,       & ! me
            &      send%nElemsProc( iProc )*nScalars,  & ! me size
            &      MPI_LOGICAL,                        & ! data type
            &      send%proc( iProc ),                 & ! target me
            &      message_flag,                       & ! flag
            &      comm,                          & ! communicator
            &      rq_handle(iProc),                   & ! handle
            &      iErr  )                               ! error status
        end do

        ! then send the 'recv' bitmask arrays to the neighbor procs
        ! Provider view: send to the source processes
        do iProc = 1, recv%nProcs
          call mpi_isend(                              &
           &      myHalos( iProc )%bitMask,            & ! me
           &      recv%nElemsProc( iProc )*nScalars,   & ! me size
           &      MPI_LOGICAL,                         & ! data type
           &      recv%proc( iProc ),                  & ! target me
           &      message_flag,                        & ! flag
           &      comm,                           & ! communicator
           &      rq_handle( iProc + send%nProcs ),    & ! handle
           &      iErr                    )              ! error status
        end do ! iProc

        ! Wait for above communications to complete
        call mpi_waitall(               &
          & send%nProcs+recv%nProcs,    & ! total number of comm.'s to wait for
          & rq_handle,                  & ! request handles
          & status,                     & ! mpi status
          & iErr         )                ! error status
        ! ----------------------------------------------------------------------
        ! Status: We received the bitmasks from the Provider processes,
        !         from which we can now construct the correct send buffers
        ! ----------------------------------------------------------------------

        deallocate( rq_handle, status )

      end subroutine communicate_bitmask
    ! ************************************************************************ !

  end subroutine init_levelBuffers
! **************************************************************************** !

! **************************************************************************** !
  !> Select the halo elements which require all links
  !!
  !! Determine which pdfs to communicate
  !!
  subroutine set_halo_commLinks( scheme, minLevel, maxLevel, nBCs, comm_reduced, &
    &                            haloRequired )
    ! --------------------------------------------------------------------------
    !> scheme information including fluid, boundary and flow information
    type(mus_scheme_type), intent(in) :: scheme
    !> Global information
    integer, intent(in) :: minLevel, maxLevel, nBCs
    !> reduced communication
    logical, intent(in) :: comm_reduced
    !>
    type( logical_array_type ), allocatable, intent(inout) :: haloRequired(:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iLevel, iSourceElem
    integer :: iBnd, iField, iNeigh, iProc
    integer :: targetElem, sourceElem, targetLevel, sourceLevel
    integer :: nSourceElems
    integer :: elemPos, haloOffset
    ! --------------------------------------------------------------------------

    write(logUnit(5),"(A)") "Set halo elements comm links"

    allocate( haloRequired(minLevel:maxLevel) )
    do iLevel = minLevel, maxLevel
      allocate( haloRequired(iLevel)%val( scheme%pdf(iLevel)%nElems_halo ) )
      haloRequired(iLevel)%val(:) = .false.
    end do

    do iLevel = minLevel, maxLevel

      ! allocate( haloRequired(iLevel)%val( scheme%pdf(iLevel)%nElems_halo ) )

      if( .not. comm_reduced ) then
        ! Communicate all links of all halos
        ! scheme%pdf(iLevel)%haloRequired(:) = .true.
        haloRequired(iLevel)%val(:) = .true.
      else ! Choose which ones are really required
        ! ----------------------------------------------------------------------
        ! GhostFromFiner
        ! check if there are halo deps
        if( iLevel < maxLevel ) then
          targetLevel = iLevel
          sourceLevel = iLevel + 1
          haloOffset  = scheme%levelDesc( sourceLevel )%offset(1,eT_halo)

          do iElem = 1, scheme%pdf( targetLevel )%nElems_ghostFromFiner
            ! position of coarser ghost element in levelDesc total list
            targetElem = iElem + scheme%levelDesc( targetLevel )          &
              &                      %offset( 1, eT_ghostFromFiner)
            nSourceElems = scheme%levelDesc( targetLevel )%               &
              &                                depFromFiner( iElem )%elem%nVals
            ! number of fine elements(childs)
            do iSourceElem = 1, nSourceElems
              ! position of child element in levelDesc total list
              sourceElem = scheme%levelDesc( targetLevel )%               &
                &                 depFromFiner( iElem )%elem%val( iSourceElem )
              call set_requiredLink(                                           &
                &          haloElem  = sourceElem,                             &
                &          offset    = haloOffset,                             &
                &          linkField = haloRequired( sourceLevel )%val )
            end do ! iSourceElem
          end do ! iElem
        end if !iLevel < maxLevel
        ! GhostFromFiner
        ! ----------------------------------------------------------------------

        ! ----------------------------------------------------------------------
        ! GhostFromCoarser
        ! check if there are halo deps
        if( iLevel < maxLevel ) then
          targetLevel = iLevel + 1
          sourceLevel = iLevel
          haloOffset =scheme%levelDesc( sourceLevel )%offset(1,eT_halo)
          do iElem = 1, scheme%pdf( targetLevel )%nElems_ghostFromCoarser
            ! position of ghost finer element in levelDesc total list
            targetElem = iElem + scheme%levelDesc( targetLevel )          &
              &                      %offset( 1, eT_ghostFromCoarser)
            ! number of coarser elements that finer ghost dependent on.
            ! usually this is just 1.
            nSourceElems = scheme%levelDesc( targetLevel )                &
              &            %depFromCoarser( iElem )%elem%nVals
            do iSourceElem = 1, nSourceElems
              sourceElem = scheme%levelDesc( targetLevel )                &
                &          %depFromCoarser( iElem )%elem%val( iSourceElem )
              ! Check the type of the source element.
              ! If it is a halo -> mark as require all links
              call set_requiredLink(                                           &
                &          haloElem  = sourceElem,                             &
                &          offset    = haloOffset,                             &
                &          linkField = haloRequired( sourceLevel )%val )
            end do ! iSourceElem
          end do ! iElem
        end if !level<maxLevel
        ! GhostFromCoarser
        ! ----------------------------------------------------------------------

        haloOffset = scheme%levelDesc( iLevel )%offset(1, eT_halo)
        ! ----------------------------------------------------------------------
        ! Run over recvbufferFromCoarser
        do iProc = 1, scheme%levelDesc(iLevel)%recvBufferFromCoarser%nProcs
          do iElem = 1, scheme%levelDesc( iLevel )%recvBufferFromCoarser% &
              &                                             nElemsProc( iProc )
            elemPos = scheme%levelDesc( iLevel )%recvBufferFromCoarser%   &
              &                                   elemPos( iProc )%val( iElem )
            call set_requiredLink(                                             &
              &                 haloElem  = elemPos,                           &
              &                 offset    = haloOffset,                        &
              &                 linkField = haloRequired( iLevel )%val )
          end do
        end do
        ! ----------------------------------------------------------------------

        ! ----------------------------------------------------------------------
        ! Run over recvbufferFromFiner
        do iProc = 1, scheme%levelDesc( iLevel )%recvBufferFromFiner%nProcs
          do iElem = 1, scheme%levelDesc( iLevel )%recvBufferFromFiner%   &
            &                                               nElemsProc( iProc )
            elemPos = scheme%levelDesc( iLevel )%recvBufferFromFiner%     &
              &                                   elemPos( iProc )%val( iElem )
            call set_requiredLink(                                             &
              &                 haloElem  = elemPos,                           &
              &                 offset    = haloOffset,                        &
              &                 linkField = haloRequired( iLevel )%val )
          end do
        end do
        ! ----------------------------------------------------------------------

        ! ----------------------------------------------------------------------
        ! Run over all the boundary neighbors
        ! If they require neighbor infos, it is usually integral quantities
        ! ( macrosopic)  and thus require all links
        do iBnd = 1, nBCs
          do iElem = 1, scheme%globBC( iBnd )%nElems( iLevel )
            do iField = 1, scheme%nFields
              do iNeigh = 1, scheme%field( iField )%bc( iBnd )%nNeighs
                ! Get the neighbor source element
                sourceElem = scheme%field( iField )%bc( iBnd )                 &
                  &                 %neigh( iLevel )%posInState( iNeigh, iElem )
                call set_requiredLink(                                         &
                  &             haloElem  = sourceElem,                        &
                  &             offset    = haloOffset,                        &
                  &             linkField = haloRequired( iLevel )%val )
              end do ! iNeigh
            end do  ! iField
          end do ! iElem
        end do ! iBnd

        ! Communicate all links of boundary neighbors which are halo.
        ! There are some boundaries which requires all links of neighbors
        ! like wall_linearInterpolation.
        ! Do this for boundaries other than wall
        ! KM: elemPos for wall boundary is not updated in updateBCList
        ! Also wall boundary does not require neighbors
        do iBnd = 1, nBCs
          if( .not. scheme%globBC(iBnd)%isWall ) then
            do iElem = 1, scheme%globBC( iBnd )%nElems( iLevel )
              elemPos = scheme%globBC(iBnd)%elemLvl(iLevel)%elem%val(iElem)
              do iNeigh = 1, scheme%layout%fStencil%QQN
                ! Get the neighbor source element
                sourceElem = scheme%levelDesc(iLevel)%neigh(1)%        &
                  & nghElems( scheme%layout%fstencil%cxdirinv(ineigh),elemPos)
                ! if sourceElem < 0 then this neighbor is boundary
                call set_requiredLink(                                 &
                  &      haloElem  = sourceElem,                       &
                  &      offset    = haloOffset,                       &
                  &      linkField = haloRequired( iLevel )%val )
              end do ! iNeigh
            end do ! iElem
          end if ! not wall
        end do ! iBnd
        ! Boundary neighbors
        ! ----------------------------------------------------------------------
      end if ! comm_reduced?

    end do ! iLevel

  contains

    ! ************************************************************************ !
    !> set haloRequired link to true
    !!
    subroutine set_requiredLink( haloElem, offset, linkField )
      ! ------------------------------------------------------------------------
      !> halo element position in the state array
      integer, intent(in)    :: haloElem
      !> halo element offset in the state array
      !! nr. of elements in state array before halo elements
      integer, intent(in)    :: offset
      !> halo required array used in init_recvBuffers to reduce communication links
      logical, intent(inout) :: linkField(:)
      ! ------------------------------------------------------------------------
      if( haloElem > offset ) then
        ! haloRequired is only nElems_halo long, so remove the preceeding
        ! elem nums
        linkField( haloElem - offset ) = .true.
      end if

    end subroutine set_requiredLink
    ! ************************************************************************ !

  end subroutine set_halo_commLinks
! **************************************************************************** !


! **************************************************************************** !
  !> Create the communication buffers
  !!
  !! Assign the positions in state vectors, where to fetch from or save to
  !! Solver specific part to identify the actual data locations
  !! to send use within the communication buffers.
  !! @todo: this should be pattern specific, with typed exchange we
  !!        do not need to store the pos, as it is done inside MPI
  !!        with the defined datatype
  !!
  !! Access to state array: SAVE
  !! ---------------------------
  !! Why? Because we need to work on the pdfs on which the
  !! compute kernel was working. The kernel stored to SAVE
  !! All actions taking place on the state array after the kernel
  !! need to be performed with SAVE access
  !!
  subroutine init_recvBuffers( me, myHalos, pattern, requireAll, haloRequired, &
    &                          halo_offset, scheme, stat, neigh, nSize )
    ! --------------------------------------------------------------------------
    !> receive communication buffer
    type(tem_communication_type), intent(inout) :: me
    !> contains bitmask of my halos with array size of recv%nProcs
    type(halo_commBitmask_type), intent(inout)  :: myHalos(:)
    !> communication pattern
    type(tem_commPattern_type), intent(in)      :: pattern
    !> requires complete element information or only required links?
    logical, intent(in) :: requireAll
    !> Level descriptor
    integer, intent(in) :: halo_offset
    !> Scheme information on fluid, boundary and flow properties
    type( mus_scheme_type ), intent(in)         :: scheme
    !> Statistics
    type(mus_statistics_type), intent(inout)    :: stat
    !> PDF neighbor array
    integer, intent(in) :: neigh(:)
    !> number of Elements in neigh array
    integer, intent(in) :: nSize
    logical, intent(in) :: haloRequired(:)
    ! --------------------------------------------------------------------------
    ! positions, in state vectors, of links with fluid neighbors
    integer, allocatable :: pos(:)
    integer :: iElem, iDir, iProc, elemPos
    integer :: counter ! number of total links that acturally needs to update
    integer :: iField, nScalars, varPos, neighDir, QQ
    ! Include the current direction for the communciation
    logical :: includeInComm
    integer :: nghElem ! neighbor element position in total list
    integer :: neighPos ! position of neighbor in the neigh array list
    ! --------------------------------------------------------------------------

    nScalars = scheme%varSys%nScalars
    QQ = scheme%layout%fStencil%QQ

    allocate(pos(maxval(me%nElemsProc*nScalars)))

    do iProc = 1, me%nProcs
      counter = 0
      do iElem = 1, me%nElemsProc( iProc )
        ! halo element position in the levelDesc total list
        elemPos = me%elemPos( iProc )%val( iElem )
        ! Here the check for the bitmasks has to be performed
        ! Check only pdf states. No need to communicate omega
        do iDir = 1, QQ
          includeInComm = .false.
          ! Now check all the neighbors of the current halo element
          ! - If we encounter fluid or ghost element, then the direction is
          !   required to exchange for      -> bitmask( iDir, iElem ) = .true.
          ! - If a halo element is ecountered, there is no need for
          !   communication                -> bitmask( iDir, iElem ) = .false.
          ! Get the element position of the neighbor element
          ! By accessing the neighbor array, get the position of the neighbor
          ! element's direction and reconstruct the neighbor element position
          ! in the total List
          neighDir =  scheme%layout%fstencil%cxdirinv( idir)
          neighPos = ( neighdir-1)* nsize + elempos
          nghElem  = int(( neigh( neighpos )-1)/ nscalars)+1
          ! Now check, what kind of element it is. Must be anything but halo.
          ! requireAll is true for IBM stencil or comm_reduced is false
          if( nghElem <= halo_offset .or. requireAll ) then
            includeInComm = .true.
          end if

          ! If it's a ghostHalo or boundary halo, then communicate all
          if( haloRequired( elemPos - halo_offset )) then
            includeInComm = .true.
          end if

          ! If current neighbor is a ghosthalo or boundary halo,
          ! then communicate all
          if( nghElem > halo_offset ) then
            if( haloRequired( nghElem  - halo_offset )) then
              includeInComm = .true.
            end if
          end if

          ! Field loop is inside direction loop since neighbor element
          ! identification is independent of fields
          if( includeInComm ) then
            do iField = 1, scheme%nFields
              counter = counter + 1
              ! ? SAVE ? ... for regular halos
              ! ? IDX ? ... for ghost-halos
              ! - > These two are the same for PULL!! so, no difference.
              ! Why SAVE? Because we need to work on the pdfs on which the
              ! compute kernel was working. The kernel stored to SAVE
              ! All actions taking place on the state array after the kernel
              ! need to be performed with SAVE access
              pos(counter) = (elempos-1)*nscalars+idir+(ifield-1)*qq

              varPos = scheme%varSys%method%val(                      &
                &             scheme%stateVarMap%varPos%val(iField) ) &
                &                   %state_varPos(iDir)
              myHalos( iProc )%bitMask( varPos, iElem ) = .true.
            end do ! field loop
          end if
        end do ! iDir
      end do ! iElem = 1, me%nElemsProc( iProc )

      stat%nLinks_total = stat%nLinks_total + scheme%nFields * QQ * me%nElemsProc(iProc)
      stat%nLinks_comm = stat%nLinks_comm + counter
      ! copy position array to me%pos, allocate me%val array
      write(dbgUnit(1),*) "Set recv buffer, iProc: ", iProc, ", counter: ", counter
      call pattern%initBuf_real( me%buf_real( iProc ), pos, counter )

    end do ! iProc
    deallocate(pos)

  end subroutine init_recvBuffers
! **************************************************************************** !


! **************************************************************************** !
  !> Create the communication buffers
  !!
  !! Receive buffers were before created. Now receive the buffers as
  !! send buffers from the remote processes.
  !! Assign the positions, where to fetch from and where to save to
  !! Solver specific part to identify the actual data locations
  !! to send use within the communication buffers.
  !! @todo: this should be pattern specific, with typed exchange we
  !!        do not need to store the pos, as it is done inside MPI
  !!        with the defined datatype
  !! @todo: does PULL and PUSH give the same results? as SAVE is used here
  !!
  !! Access to state array: SAVE
  !! ---------------------------
  !! Why? Because we need to work on the pdfs on which the
  !! compute kernel was working. The kernel stored to SAVE
  !! All actions taking place on the state array after the kernel
  !! need to be performed with SAVE access
  !!
  subroutine init_sendBuffers( send, remoteHalos, neigh, pattern, scheme, &
    &                          nSize )
    ! --------------------------------------------------------------------------
    !> send communication buffer
    type(tem_communication_type), intent(inout) :: send
    !> contains bitmask of remote halos with array size send%nProcs
    type(halo_commBitmask_type), intent(in)     :: remoteHalos(:)
    !> neighbor array for state array
    integer, intent(in)                         :: neigh(:)
    !> communication pattern
    type(tem_commPattern_type), intent(in)      :: pattern
    !> fluid, boundary and flow information
    type( mus_scheme_type ), intent(in)         :: scheme
    integer, intent(in)                         :: nSize
    ! --------------------------------------------------------------------------
    ! Element positions
    integer, allocatable :: pos(:)
    integer :: iElem ! count variable for elements
    integer :: iDir, iProc, counter, elemPos
    integer :: iField, nScalars, QQ, varPos
    ! --------------------------------------------------------------------------

    nScalars = scheme%varSys%nScalars
    QQ = scheme%layout%fStencil%QQ
    allocate( pos(maxval(send%nElemsProc*nScalars)))
    ! allocate( send%buf_real( send%nProcs ))

    ! Then assign the send buffer positions
    do iProc = 1, send%nProcs
      counter = 0
      do iElem = 1, send%nElemsProc( iProc )
        ! Now assign the exact positions in the sendBuf
        elemPos = send%elemPos( iProc )%val( iElem )
        ! Check only pdf states. No need to communicate omega
        do iDir = 1, QQ
          ! Field loop is inside direction loop to maintain the order of pos
          ! as in init_recvBuffers
          do iField = 1, scheme%nFields
            varPos = scheme%varSys%method%val(                     &
              &            scheme%stateVarMap%varPos%val(iField) ) &
              &                  %state_varPos(iDir)
            if( remoteHalos( iProc )%bitMask( varPos, iElem )) then
              counter = counter + 1
              ! ? SAVE ? ... for regular halos
              ! ? IDX ? ... for ghost-halos
              ! -> These two are the same for PULL!! so, no difference.
              pos( counter ) = (elempos-1)*nscalars+idir+(ifield-1)*qq
            endif
          end do ! field loop
        end do ! Dir loop
      end do ! iElem

      write(dbgUnit(1),*) "Set send buffer, iProc: ", iProc, ", counter: ", counter
      call pattern%initBuf_real( send%buf_real( iProc ), pos, counter )
    end do ! iProc
    write(dbgUnit(1),"(A)" ) ''

    deallocate(pos)

  end subroutine init_sendBuffers
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate global and sum of levels numbers of elements
  !! @todo: This routine can be seperated.
  !!
  subroutine calculate_nElems( levelDesc, proc, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    !> Level range
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    type(tem_levelDesc_type), intent(in) :: levelDesc(minLevel:maxLevel)
    !> mpi communication enviroment with mpi communicator
    type(tem_comm_env_type), intent(in) :: proc
    ! --------------------------------------------------------------------------
    integer :: iLevel, iErr, nLevels
    integer(kind=long_k) :: ElemCount(minLevel:maxLevel)
    integer(kind=long_k) :: nElems_overall
    integer(kind=long_k) :: nElems_allLevel( minLevel:maxLevel )
    integer(kind=long_k) :: nElems_fluidLevel( minLevel:maxLevel )
    integer :: myGhostFromFiner, myGhostFromCoarser
    integer :: nElems_totalGhostFromFiner
    integer :: nElems_totalGhostFromCoarser
    ! --------------------------------------------------------------------------

    nLevels = maxLevel - minLevel + 1
    ! Count the amount of different element types
    ! counting levelwise elements ...
    call tem_horizontalSpacer(fUnit = logUnit(6))
    write(logUnit(6),"(A)") 'Elements per level'

    call tem_horizontalSpacer(fUnit = logUnit(6))
    ! done counting levelwise  elements

    ! Calculate the total number of fluid elements in ALL LEVELS and processes

    ! count level-wise global nElems_total
    ElemCount(minLevel:maxLevel) =   &
      &     levelDesc( minLevel:maxLevel )%elem%nElems( eT_fluid ) &
      &   + levelDesc( minLevel:maxLevel )%elem%nElems( eT_ghostFromCoarser ) &
      &   + levelDesc( minLevel:maxLevel )%elem%nElems( eT_ghostFromFiner )  &
      &   + levelDesc( minLevel:maxLevel )%elem%nElems( eT_halo )
    call mpi_allreduce( ElemCount(minLevel:maxLevel), &
      &                 nElems_allLevel(minLevel:maxLevel),             &
      &                 nLevels, mpi_integer8, mpi_sum, proc%comm, iErr )

    ! count level-wise global nElems_fluid
    ElemCount(minLevel:maxLevel) =   &
      &     levelDesc( minLevel:maxLevel )%elem%nElems( eT_fluid )
    call mpi_allreduce( ElemCount(minLevel:maxLevel), &
      &                 nElems_fluidLevel(minLevel:maxLevel), &
      &                 nLevels, mpi_integer8, mpi_sum, proc%comm, iErr )

    ! count my ghost over all level
    myGhostFromFiner   = 0
    myGhostFromCoarser = 0
    do iLevel = minLevel, maxLevel
      myGhostFromFiner   = myGhostFromFiner &
        &                + levelDesc( iLevel )%elem%nElems( eT_ghostFromFiner )
      myGhostFromCoarser = myGhostFromCoarser &
        &                + levelDesc( iLevel )%elem%nElems( eT_ghostFromCoarser)
    end do
    ! count global ghost
    call mpi_allreduce( myGhostFromFiner, nElems_totalGhostFromFiner, &
      &                 1, mpi_integer, mpi_sum, proc%comm, iErr )
    call mpi_allreduce( myGhostFromCoarser, nElems_totalGhostFromCoarser, &
      &                 1, mpi_integer, mpi_sum, proc%comm, iErr )

    nElems_overall = sum( nElems_allLevel( minLevel:maxLevel ))
    do iLevel = minLevel, maxLevel
      write(logUnit(1),"(A,I0)") ' Total number of elements on Level: ', iLevel
      write(logUnit(1),"(A,I0)") '   Fluid:   ', nElems_fluidLevel( iLevel )
    end do

    write(logUnit(3),"(A   )") ' Global Elements information:'
    write(logUnit(3),"(A,I0)") '   nElems overall: ', nElems_overall
    write(logUnit(3),"(A,I0)") '   nElems ghost from finer: ', &
      &                        nElems_totalGhostFromFiner
    write(logUnit(3),"(A,I0)") '   nElems ghost from coarser: ', &
      &                        nElems_totalGhostFromCoarser

    write(logUnit(1),*) ''

  end subroutine calculate_nElems
! **************************************************************************** !





! **************************************************************************** !
  !> Update the neighbor stencil positions for the Boundaries
  !!
  !! Before, positions in the original treeIDs were assigned to
  !! the boundary element lists. Update these now with level-wise positions in
  !! the [[tem_construction_module:tem_levelDesc_type]] "total list".
  !! When no element is found, the value in nghElems is 0. This has to be
  !! treated in this routine. If an empty neighbor relation is encountered, the
  !! last valid neighbor is used for all other neighbors
  !! Assertion:
  !! Before this routine, we should have:
  !!  scheme%globBC( iBnd )%nElems( iLevel )
  !!  scheme%globBC( iBnd )%elemLvl( iLevel )%elem%val( iElem ))
  !!
  subroutine update_BClists( scheme, minLevel, maxLevel, nBCs, LP, remove_solid)
    ! --------------------------------------------------------------------------
    !> scheme information including fluid, boundary and flow information
    type(mus_scheme_type), intent(inout) :: scheme
    !> Level Pointer
    integer, intent(in) :: LP(:)
    !>
    integer, intent(in) :: minLevel, maxLevel, nBCs
    logical, intent(in) :: remove_solid
    ! --------------------------------------------------------------------------
    integer :: iElem, iLevel, iBnd, iField
    integer :: elemPos, iPos
    integer :: bcBufferSize
    ! --------------------------------------------------------------------------

    ! set posInNghElems in scheme%field(iField)%bc
    do iField = 1, scheme%nFields
      do iBnd = 1, nBCs
        call mus_set_posInNghElems( minLevel, maxLevel, &
          &                         scheme%layout%nStencils, &
          &                         scheme%globBC(iBnd), &
          &                         scheme%field(iField)%bc(iBnd) )
      end do
    end do

    if ( remove_solid .and. (nBCs > 0) ) then
      call remove_solid_in_bc( minLevel, maxLevel, nBCs, &
        &                      scheme%nFields, LP, &
        &                      scheme%levelDesc, scheme%globBC, scheme%field )
    end if

    ! 1. Update the element positions from the treeID to levelDesc%total list
    ! 1a. Create dynamic array of bc_elemBuffer with position of
    ! treeID in the levelDescriptor%total list and store the position
    ! in globBCLvl%elemLvl(iLevel)%posInBcElemBuf
    do iLevel = minLevel, maxLevel
      ! allocate bc_elemBuffer with initial size of all bc nElems
      bcBufferSize = 0
      do iBnd = 1, nBCs
        bcBufferSize = bcBufferSize + scheme%globBC( iBnd )%nElems( iLevel )
      end do ! iBnd
      call init( me     = scheme%levelDesc( iLevel )%bc_elemBuffer,       &
        &        length = bcBufferSize )
      do iBnd = 1, nBCs
        ! if BC kind of all fields are wall then do nothing for this boundary
        ! i.e do not add this boundary to elemBuffer or bc_elemBuffer or
        ! neigh buffers since wall boundaries are treated explicitly as
        ! bounce back via pdf%neigh array
        if ( .not. scheme%globBC(iBnd)%isWall ) then
          ! UPDATE globBD%elemLvl%elem%val and add to bc_elemBuffer
          do iElem = 1, scheme%globBC( iBnd )%nElems( iLevel )
            ! update bc element position in treelm to be pos in total list
            elemPos = LP( scheme%globBC( iBnd )%elemLvl( iLevel )% &
                &                               elem%val( iElem )  )
            scheme%globBC( iBnd )%elemLvl( iLevel )%elem%val( iElem ) = elemPos

            ! store elemPos in bc_elemBuffer
            call append( me  = scheme%levelDesc( iLevel )%bc_elemBuffer,  &
              &          val = elemPos,                                   &
              &          pos = iPos )
            ! store the position in the bc_elemBuffer in the elemBuffer
            call append ( me =  scheme%globBC( iBnd )%elemLvl( iLevel )% &
              &                                       posInBcElemBuf,    &
              &          val = iPos )
          end do ! iElem in globBC%nElems on the current level

          ! elem and elemBuffer for corner nodes
          if (scheme%globBC(iBnd)%treat_corner) then
            do iElem = 1, scheme%globBC( iBnd )%cornerBC%nElems( iLevel )
              elemPos = LP( scheme%globBC( iBnd )%cornerBC%elemLvl( iLevel ) &
                &                                         %elem%val( iElem ))

              scheme%globBC( iBnd )%cornerBC%elemLvl( iLevel )         &
                &                           %elem%val( iElem ) = elemPos

              iPos = PositionOfVal( me  = scheme%levelDesc( iLevel ) &
                &                               %bc_elemBuffer,      &
                &                   val = elemPos                    )
              ! store the position in the bc_elemBuffer in the elemBuffer
              call append ( me =  scheme%globBC( iBnd )%cornerBC            &
                &                        %elemLvl( iLevel )%posInBcElemBuf, &
                &          val = iPos                                       )
            end do
          end if ! treat corner

        end if ! all Wall
      end do ! iBnd
    end do !iLevel

  end subroutine update_BClists
! **************************************************************************** !


! **************************************************************************** !
  subroutine finalize_BClist( scheme, proc, minLevel, maxLevel, nBCs )
    ! --------------------------------------------------------------------------
    !> scheme information including fluid, boundary and flow information
    type(mus_scheme_type), intent(inout) :: scheme
    !>
    integer, intent(in) :: minLevel, maxLevel, nBCs
    !> mpi communication enviroment with mpi communicator
    type(tem_comm_env_type), intent(in) :: proc
    ! --------------------------------------------------------------------------
    integer :: iLevel, iBnd, iField
    integer :: iNeigh, iElem, iDir, neighPos, neighVal
    logical :: inValidNeigh, inValidNeighGlobal
    integer :: iErr
    ! --------------------------------------------------------------------------

    do iLevel = minLevel, maxLevel
      do iBnd = 1, nBCs

        ! if BC kind of all fields are wall then do nothing for this boundary
        ! i.e do not add this boundary to elemBuffer or bc_elemBuffer or
        ! neigh buffers since wall boundaries are treated explicitly as
        ! bounce back via pdf%neigh array
        if ( .not. scheme%globBC(iBnd)%isWall ) then

          ! allocate neighBuffer pre-collision and postcollision for each
          ! boundary in each field and initialize them
          ! @todo: SZ: think about this allocation in case there are more state
          !            variables than the QQ PDFs
          do iField = 1, scheme%nFields
            if (scheme%field( iField )%bc( iBnd )%requireNeighBufPre) then
              allocate( scheme%field( iField )%bc( iBnd )%           &
                & neigh( iLevel )%neighBufferPre(                    &
                &         scheme%field( iField )%bc( iBnd )%nNeighs, &
                &         scheme%globBC( iBnd )%nElems( iLevel )*    &
                &         scheme%layout%fStencil%QQ )                )
              scheme%field(iField)%bc(iBnd)%neigh(iLevel)%neighBufferPre &
                & = -1._rk
            end if
            if (scheme%field( iField )%bc( iBnd )%requireNeighBufPre_nNext) then
              allocate( scheme%field( iField )%bc( iBnd )%           &
                & neigh( iLevel )%neighBufferPre_nNext(              &
                &         scheme%field( iField )%bc( iBnd )%nNeighs, &
                &         scheme%globBC( iBnd )%nElems( iLevel )*    &
                &         scheme%layout%fStencil%QQ )                )
              scheme%field(iField)%bc(iBnd)%neigh(iLevel)%neighBufferPre_nNext &
                & = -1._rk
            end if

            if (scheme%field( iField )%bc( iBnd )%requireNeighBufPost) then
              allocate( scheme%field( iField )%bc( iBnd )%                     &
                & neigh( iLevel )%neighBufferPost(                             &
                &         scheme%field( iField )%bc( iBnd )%nNeighs,           &
                &         scheme%globBC( iBnd )%nElems( iLevel )*              &
                &         scheme%layout%fStencil%QQ ) )
              scheme%field(iField)%bc(iBnd)%neigh(iLevel)%neighBufferPost &
                & = -1._rk
            end if

            if ( scheme%field(iField)%bc(iBnd)%useComputeNeigh ) then
              allocate( scheme%field( iField )%bc( iBnd )%neigh( iLevel )  &
                &                                        %computeNeighBuf( &
                &         scheme%globBC( iBnd )%nElems( iLevel )*          &
                &         scheme%layout%fStencil%QQ ) )
              scheme%field( iField )%bc( iBnd )%neigh( iLevel )        &
                &                              %computeNeighBuf = -1._rk
            end if

          end do ! iField
        end if ! isWall
      end do !nBCs
    end do ! iLevel

    ! build up neigh information for higher order boundaries
    ! loop over the field to update array boundary in each field
    write(logUnit(1),"(A,I0)") 'Start filling field BC neighbor list '
    do iField = 1, scheme%nFields
      do iBnd = 1, nBCs
        ! Copy the neighboring information to the bc%elem%neigh
        if( scheme%field( iField )%bc( iBnd )%nNeighs > 0 ) then
          ! get position of neighbor element in the levelwise list
          call setFieldBCNeigh(                                       &
            &     fieldBC   = scheme%field( iField )%bc( iBnd ),      &
            &     globBC    = scheme%globBC( iBnd ),                  &
            &     levelDesc = scheme%levelDesc(minLevel:maxLevel),&
            &     minLevel  = minLevel,                     &
            &     maxLevel  = maxLevel                      )
        end if ! nNeigh>0
      end do ! iBnd
    end do ! iField

    ! KM: 20220421 Abort if neighbor is ghost from finer and does
    ! not have valid neighbor. It happens when children of
    ! ghostFromFiner is distributed.
    ! These elements cannot be treated if they try to FETCH for Pull
    ! and SAVE for Push.
    inValidNeigh = .false.
    do iLevel = minLevel, maxLevel
      do iBnd = 1, nBCs
        do iField = 1, scheme%nFields
          ! FETCH uses neighbor for Pull and SAVE uses neighbor for Push
          ! so we need to check the valid neighbor depending on BC
          if (scheme%field( iField )%bc( iBnd )%requireNeighBufPre) then
            do iNeigh = 1, scheme%field(iField)%bc(iBnd)%nNeighs
              do iElem = 1, scheme%globBC(iBnd)%nElems(iLevel)
                ! neighbor position in the total list
                neighPos =  scheme%field(iField)%bc(iBnd)%neigh(iLevel)  &
                  &                             %posInState(iNeigh, iElem)
                if (neighPos > scheme%levelDesc(iLevel)            &
                  &                  %offset(1, eT_ghostFromFiner) &
                 & .and. neighPos < scheme%levelDesc(iLevel)%offset(1, eT_halo) ) then
                  do iDir = 1, scheme%layout%fStencil%QQN
                    neighVal = scheme%levelDesc(iLevel)%neigh(1) &
                      &              %nghElems(iDir, neighPos)
                    if (neighVal == 0) then
                      inValidNeigh = .true.
                    end if
                  end do
                end if !ghostFromFiner
              end do !nElems
            end do !nNeigh
          end if
        end do ! iField
      end do ! iBnd
    end do ! iLevel

    call mpi_allreduce( inValidNeigh, inValidNeighGlobal, 1,  &
      &                 mpi_logical, mpi_lor, proc%comm, iErr )
    if (inValidNeighGlobal) then
      call tem_abort('Error: BC neigh is ghostFromFiner and &
        &it does not have valid neighbor')
    end if

  end subroutine finalize_BClist
! **************************************************************************** !


! **************************************************************************** !
  subroutine build_bcLevelPointer( bcLevelPointer, posInBndID, minBcID, &
    &                              levelPointer, bc_prop, tree, globBC )
    ! --------------------------------------------------------------------------
    !> levelPointer from posInBndID to globBC%elemLVl(:)%elem
    integer, allocatable, intent(inout) :: bcLevelPointer(:)
    !> tree element position in boundaryID
    integer, intent(in) :: posInBndID(:)
    !> minBCID for each element
    integer, intent(in) :: minBCID(:)
    !> tree element position in levelDesc total list
    integer, intent(in) :: levelPointer(:)
    !> boundary information from mesh
    type(tem_bc_prop_type), intent(in) :: bc_prop
    !> global boundary information
    type(glob_boundary_type), intent(in) :: globBC(bc_prop%nBCtypes)
    !> fluid tree from mesh
    type(treelmesh_type), intent(in) :: tree
    ! --------------------------------------------------------------------------
    integer :: iElem, iLevel, posInTree
    ! --------------------------------------------------------------------------
    allocate(bcLevelPointer(bc_prop%property%nElems))

    do iElem = 1, bc_prop%property%nElems
      posInTree = bc_prop%property%ElemID(iElem)
      iLevel = tem_levelOf( tree%treeID(posInTree) )
      if (minBcID(iElem) > 0) then
        bcLevelPointer(iElem) = PositionOfVal(                                 &
          &                       me  = globBC(minBcID(iElem))%elemLvl(iLevel) &
          &                                                   %elem,           &
          &                       val = levelPointer(posInTree)                )
      end if
    end do

  end subroutine build_bcLevelPointer
! **************************************************************************** !


! **************************************************************************** !
  !> This routine sets field BC neigh array with position of
  !! neighbor element in the inward normal direction of boundary
  !! in the levelwise list.
  !! if valid 1st neighbor does not exist return current element position.
  !! if valid higher order neighbor does not exist return last valid
  !! neighbor
  subroutine setFieldBCNeigh( fieldBC, globBC, levelDesc, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    !> Level range
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    !> field boundary with boundary neighbor info
    type(boundary_type),      intent(inout) :: fieldBC
    !> boundaries for the elements with bnd property set
    type(glob_boundary_type), intent(in)    :: globBC
    !> Level descriptor
    type(tem_levelDesc_type), intent(in)    :: levelDesc(minLevel:maxLevel)
    ! --------------------------------------------------------------------------
    integer :: iElem, iLevel, iNeigh, x(4)
    integer :: neighVal, normal(3)
    integer :: posInNghElems, stencilPos
    integer(kind=long_k) :: tOffset, neighID, globpos, treeID
    ! --------------------------------------------------------------------------
    ! count up the position where to store the neighborhood
    lvlLoop: do iLevel = minLevel, maxLevel
      write(logUnit(6),"(A,I0)") '  Fill field BC neighbor list on level: ', &
        &                        iLevel
      allocate( fieldBC%neigh( iLevel )%posInState( fieldBC%nNeighs,         &
        &                                           globBC%nElems( iLevel ) ))



      ! loop over all BCelems without ghostelems having BC property
      elemLoop: do iElem = 1, globBC%nElems_Fluid( iLevel )
        stencilPos = fieldBC%elemLvl(iLevel)%stencilPos(iElem)
        ! current element index in nghElems
        ! counter(iLevel, stencilPos) = counter(iLevel, stencilPos) + 1
        posInNghElems = fieldBC%elemLvl(iLevel)%posInNghElems( iElem )


        neighLoop: do iNeigh = 1, fieldBC%nNeighs
          ! neighbor position in the levelwise list
          neighVal = levelDesc( iLevel )%neigh( stencilPos )%                  &
            &                            nghElems( iNeigh, posInNghElems )

          ! No valid neighbor found. set to boundary element pos if iNeigh = 1
          ! else set last valid neigbor
          if ( neighVal < 1 ) then
            if ( iNeigh == 1 ) then

              ! Can not even find the first neighbor, i.e. no neighbors at all
              ! so save simply current bnd elem pos (itself as neighbor)
              fieldBC%neigh( iLevel )%posInState( :, iElem ) = &
                &                     globBC%elemLvl( iLevel )%elem%val( iElem )
              write(logUnit(1),"(A,I5)") 'no neigh found for elem:', iElem
            else

              ! Neighbor fails, so take last neighbor which was valid
              fieldBC%neigh( iLevel )%posInState( iNeigh:, iElem ) = &
                &      fieldBC%neigh( iLevel )%posInState( iNeigh-1, iElem )
              write(logUnit(1),"(A,I0,A,I5,A)") 'neighbor number ',iNeigh,     &
                &    ' not found for elem:', iElem,'. Use last found neighbor.'
            end if ! iNeigh = 1
            ! no need to find further neighbors
            exit neighLoop
          else
            !write(dbgUnit(1),*) 'setting neighVal ', neighVal

            ! regular neighbor: assign
            ! fieldBC%elemLvl( iLevel )%neigh( iNeigh, iElem ) = neighVal
            ! store neighbor position in original treeID list
            fieldBC%neigh( iLevel )%posInState( iNeigh, iElem ) = neighVal
          end if ! neighVal < 1
        end do neighLoop
      end do  elemLoop
    end do  LvlLoop


    ! loop over all remaining ghostelems with BC prop
    do iLevel = minLevel, maxLevel
      write(logUnit(6),"(A,A,I0)") 'Fill field BC neighbor list with ', &
        &                          'BndGhostelem on level: ', iLevel

      do iElem = globBC%nElems_Fluid( iLevel ) +1, globBC%nElems( iLevel )
        globpos = globBC%elemlvl(ilevel)%elem%val(ielem)
        treeID = levelDesc(iLevel)%total( globpos )
        x = tem_CoordOfId( treeID )
        normal = globBC%elemlvl( iLevel )%normal%val(:, iElem)
        neighloop2:do iNeigh = 1, fieldBC%nNeighs
          tOffset = tem_FirstIdAtLevel( x(4) )
          neighID  =  tem_IdOfCoord(                &
            &            [ x(1) + normal(1)*iNeigh, &
            &              x(2) + normal(2)*iNeigh, &
            &              x(3) + normal(3)*iNeigh, &
            &              x(4) ], tOffset)
          neighVal = tem_treeIDinTotal( neighID, levelDesc( iLevel ) )
          ! No valid neighbor found. set to boundary element pos if iNeigh = 1
          ! else set last valid neigbor
          if ( neighVal < 1 ) then
            if ( iNeigh == 1 ) then
              ! Can not even find the first neighbor, i.e. no neighbors at all
              ! so save simply current bnd elem pos (itself as neighbor)
              fieldBC%neigh( iLevel )%posInState( :, iElem ) = &
                &                     globBC%elemLvl( iLevel )%elem%val( iElem )
              write(logUnit(1),"(A,I5)") 'no neigh found for elem:', iElem
            else
              ! Neighbor fails, so take last neighbor which was valid
              fieldBC%neigh( iLevel )%posInState( iNeigh:, iElem ) = &
                &      fieldBC%neigh( iLevel )%posInState( iNeigh-1, iElem )
              write(logUnit(1),"(A,I0,A,I5,A)") 'neighbor number ',ineigh,     &
                &    ' not found for elem:', iElem,'. Use last found neighbor.'
            end if ! iNeigh = 1
            ! no need to find further neighbors
            exit neighLoop2
          else
            ! regular neighbor: assign
            ! fieldBC%elemLvl( iLevel )%neigh( iNeigh, iElem ) = neighVal
            ! store neighbor position in original treeID list
            fieldBC%neigh( iLevel )%posInState( iNeigh, iElem ) = neighVal
          end if ! neighVal < 1
        end do neighLoop2
      end do
    end do

  end subroutine setFieldBCNeigh
! **************************************************************************** !


! **************************************************************************** !
  !> Assemble the level-wise list of elements which adhere to the boundary
  !! conditions.
  !!
  !! The boundaries will then be treated for each level one by one, running over
  !! the list of elements.
  !! The bitmasks are set for the directions pointing into the domain.
  !! In the LBM, the incoming densities have to be updated.
  !! As boundaries are being set before the kernel, the state arrays have to be
  !! stored to the FETCH position
  !!
  subroutine build_BClists( globBC, tree, bc_prop, minLevel, maxLevel, &
    &                       layout, field, comm )
    ! --------------------------------------------------------------------------
    !> fluid tree from mesh
    type( treelmesh_type ), intent(in)                  :: tree
    !> boundary information from mesh
    type( tem_bc_prop_type ), intent(in)                :: bc_prop
    !> contains pdf global information
    integer,  intent(in) :: minLevel, maxLevel
    !> scheme layout
    type( mus_scheme_layout_type ), intent(in)          :: layout
    !> boundaries for the elements with bnd property set
    type(glob_boundary_type), intent(out), allocatable  :: globBC(:)
    !> field type
    type( mus_field_type ), intent(in)                  :: field(:)
    !> mpi communication enviroment with mpi communicator
    integer, intent(in)                                 :: comm
    ! --------------------------------------------------------------------------
    integer :: iBnd, nBCs, iField
    ! --------------------------------------------------------------------------
    write(logUnit(1),*) 'Creating boundary lists ...'

    nBCs = bc_prop%nBCtypes
    if( nBCs > 0 ) then

      ! allocate bc states for each boundary with size of nFields
      ! to collect bc states of all fields of this particular boundary
      ! into one global list. Needed when boundary condition is
      ! dependent on all fields.
      ! contains element information about each boundary type
      allocate( globBC( nBCs ) )
      do iBnd = 1, nBCs
        allocate( globBC(iBnd)%elemLvl(           minLevel:maxLevel ) )
        allocate( globBC(iBnd)%nElems(            minLevel:maxLevel ) )
        allocate( globBC(iBnd)%nElems_totalLevel( minLevel:maxLevel ) )
        allocate( globBC(iBnd)%nElems_Fluid(      minLevel:maxLevel ) )
        ! check if all field BCs are wall
        globBC(iBnd)%isWall = mus_check_allWall( size(field), field, iBnd )

        ! collect corner elements if any field bc requires it for moments BC
        do iField = 1, size(field)
          globBC(iBnd)%treat_corner = field(iField)%bc(iBnd)%treat_corner
          if (globBC(iBnd)%treat_corner) exit
        end do
        if (globBC(iBnd)%treat_corner) then
          allocate( globBC(iBnd)%cornerBC%nElems(   minLevel:maxLevel ) )
          allocate( globBC(iBnd)%cornerBC%elemLvl(  minLevel:maxLevel ) )
        end if
      end do
      ! copy bc label from mesh boundary to solver globBC
      globBC(:)%label = bc_prop%BC_label(:)

      ! boundary%hasQval is allocated only if prp_hasQval is set in seeder
      ! so initialize them first and use globBC(:)%hasQval in rest of the
      ! routine
      globBC(:)%hasQVal = .false.
      if (allocated(bc_prop%hasQval)) globBC(:)%hasQVal = bc_prop%hasQval(:)

      ! 1. Identify the number of boundary condition elements of each BC type
      call countnBnds( globBC      = globBC,              &
        &              boundaryID  = bc_prop%boundary_ID, &
        &              tree        = tree,                &
        &              stencil     = layout%fStencil,     &
        &              nBCs        = nBCs,                &
        &              comm        = comm                 )

      ! 2. Allocate the BC lists
      call allocateBCList( globBC   = globBC,             &
        &                  nBCs     = nBCs,               &
        &                  minLevel = minLevel,           &
        &                  maxLevel = maxLevel,           &
        &                  QQN      = layout%fStencil%QQN )

      ! 3. Assign the BC lists
      ! Run over all the elements with the property boundary
      ! and check each direction.
      ! Assign all common boundaries to the level-wise representation
      call assignBCList( tree        = tree,            &
        &                stencil     = layout%fStencil, &
        &                bc_prop     = bc_prop,         &
        &                globBC      = globBC           )

      ! 4 Normalize the normal vectors for regular and corner elements
      !   and choose the corresponding prevailing direction from the stencil
      call normalizeBC( globBC   = globBC,   &
        &               field    = field,    &
        &               nBCs     = nBCs,     &
        &               minLevel = minLevel, &
        &               maxLevel = maxLevel, &
        &               layout   = layout,   &
        &               comm     = comm      )

    else
      allocate( globBC(0) )
    end if ! nBCs > 0 ?

  end subroutine build_BClists
! **************************************************************************** !


! **************************************************************************** !
  !> Identify the number of boundary condition elements of each BC type
  !! and number of elements with multiple BC types.
  !! This is 1st step in Build_BClists
  subroutine countnBnds( globBC, boundaryID, tree, stencil, nBCs, comm )
    ! --------------------------------------------------------------------------
    !> boundary information from mesh
    integer(kind=long_k),           intent(in) :: boundaryID(:,:)
    !> fluid tree from mesh
    type( treelmesh_type ),         intent(in) :: tree
    !> stencil
    type( tem_stencilHeader_type ), intent(in) :: stencil
    !> contains pdf global information
    integer,                        intent(in) :: nBCs
    !> boundaries for the elements with bnd property set
    type(glob_boundary_type), intent(inout)    :: globBC(nBCs)
    !> mpi communication enviroment with mpi communicator
    integer,                        intent(in) :: comm
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, iBnd, level, iMeshDir, iErr, iLevel
    integer :: minLevel, maxLevel
    logical :: found( nBCs )
    integer(kind=long_k) :: bID
    !> number of local elements of each boundary type in each level
    integer, allocatable :: nBnds(:,:)
    !> number of local elements with multiple boundary type in each level
    integer, allocatable :: nCornerBnds(:,:)
    integer, allocatable :: nBnds_total(:,:)
    ! --------------------------------------------------------------------------

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel

    ! number of local elements of each boundary type in each level
    allocate( nBnds( minLevel:maxLevel, nBCs ))
    ! number of total elements(all processes) of each boundary type in each
    ! level
    allocate( nBnds_total( minLevel:maxLevel, nBCs ))
    ! number of elements with multiple boundaries in each level
    allocate( nCornerBnds( minLevel:maxLevel, nBCs ))

    nBnds       = 0
    nCornerBnds = 0

    do iElem = 1, tree%property(1)%nElems
      found = .false.
      level = tem_LevelOf( tree%treeID( tree%property(1)%ElemID(iElem)))
      do iDir = 1, stencil%QQN
        iMeshDir = stencil%map( iDir )
        ! the index in boundary_ID does not correspond to the stencil!
        bID = boundaryID( iMeshDir, iElem )
        if( bID > 0  ) then
          ! checks if number of boundary IDs dumped by seeder is not
          ! greater than number of boundary types in bnd.lua
          if( bID > nBCs ) then
            write(logUnit(1),*) ' found BID:     ', bID
            write(logUnit(1),*) ' where the highest should be ', nBCs
            write(logUnit(1),*) ' at element     ', iElem
            write(logUnit(1),*) ' found higher bID than nr of bcs exist.'
          end if

          if ( .not. found ( bID ) ) then
            ! Increase the counter only, if this bID has not been
            ! identified before at iElem
            nBnds( level, bID ) = nBnds( level, bID )+1
            found( bID ) = .true.
          end if ! found ( bID ) == .false.
        end if ! bid > 0 ?
      end do ! iDir

      ! if element is intersected by more than one boundary
      ! also count it as a corner bc element
      if ( count(found) > 1 ) then
        do iBnd = 1, nBCs
          if (found(iBnd)) then
            nCornerBnds(level, iBnd) = nCornerBnds(level, iBnd) + 1
          endif
        enddo
      endif !count > 1 for corner nodes

    end do ! iElem

    ! Counting number of boundary elements on local partition and global
    ! Also copies field bcstates to globBC bcstates to access all
    ! field bcstates in each boundary routines
    do iBnd = 1, nBCs
      ! number of elements for each iBnd on each level
      globBC( iBnd )%nElems(minLevel:maxLevel) = nBnds(minLevel:maxLevel, iBnd )
      globBC( iBnd )%nElems_fluid(minLevel:maxLevel) &
        & = nBnds(minLevel:maxLevel, iBnd )
      ! Sum up the elements (all levels) of the local partition for each iBnd
      globBC( iBnd )%nElems_local = sum( nBnds(minLevel:maxLevel, iBnd) )

      ! number of elements with multiple boundaries
      if (globBC(iBnd)%treat_corner) then
        globBC( iBnd )%cornerBC%nElems(minLevel:maxLevel) = &
          &                  nCornerBnds(minLevel:maxLevel, iBnd )
      end if
    end do ! iBnd

    ! calculate total number of boundary elements in all levels and processes
    call mpi_allreduce( nBnds, nBnds_total, (maxLevel-minLevel+1)*nBCs,  &
      &                 mpi_integer, mpi_sum, comm, iErr      )

    ! Avoid division by 0.
    nBnds_total = max(nBnds_total,1)

    call tem_horizontalSpacer( fUnit = dbgUnit(1), before = 1 )
    write(dbgUnit(1),"(A)") "Number of boundary elements"
    do iBnd = 1, nBCs
      globBC( iBnd )%nElems_totalLevel(:) = nBnds_total( :, iBnd )
      globBC( iBnd )%nElems_total = sum(globBC( iBnd )%nElems_totalLevel(:))
      write(dbgUnit(1), "(A,I0,A)") "iBC: ", iBnd, ", label: " &
        &               //trim(globBC(iBnd)%label)
      do iLevel = minLevel, maxLevel
        write(dbgUnit(1),"(A,I2,A,I8,A,F5.2)") "level: ", iLevel, &
          &                 ", nElems: ", nBnds(iLevel,iBnd), &
          &                 ", % of nGlobal: ", &
          &     dble(nBnds(iLevel,iBnd))/dble(nBnds_total(iLevel,iBnd)) * 100._rk
      end do
      write(dbgUnit(1),"(A)") ''
    end do
    call tem_horizontalSpacer( fUnit = dbgUnit(1), after = 1 )

   end subroutine countnBnds
! **************************************************************************** !


! **************************************************************************** !
  !> Allocate BC lists, 2nd step in Build_BClists
  subroutine allocateBCList( globBC, nBCs, minLevel, maxLevel, QQN )
    ! --------------------------------------------------------------------------
    integer,                  intent(in)    :: nBCs, minLevel, maxLevel, QQN
    !> boundaries for the elements with bnd property set
    type(glob_boundary_type), intent(inout) :: globBC(nBCs)
    ! --------------------------------------------------------------------------
    integer :: iBC, iLevel
    ! --------------------------------------------------------------------------

    do iBC = 1, nBCs
      do iLevel = minLevel, maxLevel

        call mus_init_bc_elems( me      = globBC(iBC)%elemLvl(iLevel), &
          &                     nElems  = globBC(iBC)%nElems(iLevel),  &
          &                     QQN     = QQN,                         &
          &                     hasQVal = globBC(iBC)%hasQVal          )

        if( globBC(iBC)%nElems(iLevel) > 0) then
          write(logUnit(6),"(A,A,I2,A,I0)") ' BC: '//trim(globBC(iBC)%label),&
            &  ', level: ', iLevel, ', nElems: ', globBC( iBC )%nElems( iLevel )
          globBC( iBC )%elemLvl( iLevel )%elem%val(:) = 0
          globBC( iBC )%elemLvl( iLevel )%bitmask%val = .false.
          globBC( iBC )%elemLvl( iLevel )%normal%val  = 0
        end if

        if (globBC(iBC)%treat_corner) then
          ! elements with multiple boundaries
          call mus_init_bc_elems(                             &
            & me      = globBC(iBC)%cornerBC%elemLvl(iLevel), &
            & nElems  = globBC(iBC)%cornerBC%nElems(iLevel),  &
            & QQN     = QQN,                                  &
            & hasQVal = .false.                               )

          if ( globBC(iBC)%cornerBC%nElems(iLevel) > 0 ) then
            write(logUnit(6),"(A,A,I2,A,I0)") ' BC: '//trim(globBC(iBC)%label),&
              &  ', level: ', iLevel, &
              &  ', corner nElems: ', globBC( iBC )%cornerBC%nElems( iLevel )
            globBC(iBC)%cornerBC%elemLvl( iLevel )%elem%val(:) = 0
            globBC(iBC)%cornerBC%elemLvl( iLevel )%bitmask%val = .false.
            globBC(iBC)%cornerBC%elemLvl( iLevel )%normal%val  = 0
          end if
        end if

      end do !iLevel
    end do !iBnd
  end subroutine allocateBCList
! **************************************************************************** !


! **************************************************************************** !
  !> This routine assigns the BC lists
  !! Run over all the elements with the property boundary
  !! and check each direction.
  !! Assign all common boundaries to the level-wise representation
  !! 3rd step in build_BCLists
  subroutine assignBCList( tree, bc_prop, stencil, globBC )
    ! --------------------------------------------------------------------------
    !> fluid tree from mesh
    type( treelmesh_type ), intent(in)         :: tree
    !> boundary information from mesh
    type( tem_bc_prop_type ), intent(in)       :: bc_prop
    !> stencil
    type( tem_stencilHeader_type ), intent(in) :: stencil
    !> boundaries for the elements with bnd property set
    type( glob_boundary_type ), intent(inout)  :: globBC(:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, iBnd, iMeshDir
    integer(kind=long_k) :: bID
    integer :: normal(3), level
    logical :: corner_bitmask( stencil%QQN ), bitmask(stencil%QQN), wasadded
    logical :: found( bc_prop%nBCtypes )
    integer :: iElem_qVal
    integer :: weight(stencil%QQN), length, posInTree
    real(kind=rk) :: qval(stencil%QQN)
    ! number of elem for each BC
    integer, allocatable :: nBnds(:)
    real(kind=rk) :: minQVal, qValWeight(stencil%QQN), qValiDir
    ! --------------------------------------------------------------------------

    ! number of local elements of each boundary type
    allocate( nBnds( bc_prop%nBCtypes ))
    !Reset the counters
    nBnds       = 0
    iElem_qVal  = 0
    qval        = 0

    ! calculate weights
    ! Get length of currently treated direction and determine
    ! If it is one of the main coordinate axis
    do iDir = 1, stencil%QQN
      length   = stencil%cxDir( 1, iDir )**2 &
        &      + stencil%cxDir( 2, iDir )**2 &
        &      + stencil%cxDir( 3, iDir )**2

      ! For the normals, we weight the orthogonal directions
      ! of the stencil more heavily
      if( length == 1) then
        weight(iDir) = 4  ! axis orthogonal weighs twice as much as all others
        qValWeight(iDir) = 1.0_rk
      else if ( length == 2 ) then
        weight(iDir) = 2  ! non-axis orthogonal
        qValWeight(iDir) = sqrt(2.0_rk)
      else ! length == 3
        weight(iDir) = 1  ! non-axis orthogonal
        qValWeight(iDir) = sqrt(3.0_rk)
      end if
    end do

    do iElem = 1, tree%property(1)%nElems

      posInTree = tree%Property(1)%ElemID(iElem)
      ! Count the elem with qVal property to access the qVal from
      ! bc_prop%qVal
      if( btest( tree%elemPropertyBits( posInTree ), prp_hasQVal) ) &
        & iElem_qVal = iElem_qVal + 1

      found = .false.
      level = tem_LevelOf( tree%treeID( posInTree ) )
      minQVal = huge(minQVal)
      ! run over each direction in the stencil
      do iDir = 1, stencil%QQN
        ! Stencil order of directions might differ from
        ! the treelm numbering. -> look up correct treelm direction number
        iMeshDir = stencil%map( iDir )
        ! Get the boundary ID as stored in the mesh
        bID      = bc_prop%boundary_ID( iMeshDir, iElem )
        ! KM: 20170825 \todo WAll boundary can be skipped here since it is treated
        ! with simple boundary back
        ! Exclude no-boundary or periodic boundaries
        if( bID > 0 ) then
          ! If this boundary ID was not found for current element yet,
          ! assign to list

          if( .not. allocated( globBC(bID)%elemLvl(level)%elem%val))then
            write(dbgUnit(1),*) 'Error: bc not allctd lvl ', level,        &
              &                     ' bid ', bID
          end if

          ! append the position of the treeID in tree%treeID list
          ! to the boundary elem lists
          call append( me = globBC( bID )%elemLvl( level )%elem,             &
            &         val = posInTree,                                       &
            &         pos = nBnds( bID ),                                    &
            &    wasadded = wasadded )

          if ( wasadded ) then
            bitmask     = .false.
            normal      = 0
            !append the bitmask for bID
            call append( me = globBC( bID )%elemLvl( level )%bitmask, &
              &          val = bitmask  )

            !append normal
            call append( me = globBC( bID )%elemLvl( level )%normal, &
              &          val = normal )
            if ( globBC( bID )%hasQVal ) then
              !append qVal
              call append( me = globBC( bID )%elemLvl( level )%qVal, &
                &         val = qVal  )
            end if ! haQval

            found( bID ) = .true.
          end if

          ! Set the bitmask for the incoming directions which have to be
          ! updated, which is opposed to BCid direction.
          globBC( bID )%elemLvl( level )%bitmask%val(       &
            & stencil%cxDirInv(iDir), nBnds( bID ) ) = .true.

          if( globBC( bID )%hasQVal ) then
            ! Set the qvalue from the boundary array
            ! q-values stored in outgoing direction
            globBC( bID )%elemLvl( level )%qVal%val(iDir, nBnds( bID ) ) &
              & = bc_prop%qVal(iMeshDir, iElem_qVal)
            if (bc_prop%qVal(iMeshDir, iElem_qVal) > 0.0_rk) then
              qValiDir = qValWeight(iDir) * bc_prop%qVal(iMeshDir, iElem_qVal)
              minQVal = min(minQVal, qValiDir)
              if (minQVal == qValiDir ) then
                globBC( bID )%elemLvl( level )%normal%val(:, nBnds( bID ) ) &
                  & = - stencil%cxDir( :, iDir )
              end if
            end if
          else
            ! Add this boundary direction to the sum for the normal vector
            ! use the minus, as we want to go into the flow domain, and the
            ! stencil is pointing outside
            globBC( bID )%elemLvl( level )%normal%val(:, nBnds( bID ) )        &
              &  =  globBC( bID )%elemLvl( level )%normal%val(:, nBnds( bID ) )&
              &     -  weight(iDir) * stencil%cxDir( :, iDir )
          end if ! hasQVal
        end if ! bID > 0
      end do  ! iDir

      ! if element is intersected by more than one boundary
      if ( count(found) > 1 ) then
        corner_bitmask = .false.
        normal = 0
        do iBnd = 1, bc_prop%nBCtypes
          if (found(iBnd) .and. globBC(iBnd)%treat_corner ) then
            ! Assign the position of the treeID to the corner elem list
            call append( me = globBC( iBnd )%cornerBC%elemLvl( level )%elem,   &
              &          val = posInTree )

            corner_bitmask = corner_bitmask .or. &
              & globBC( iBnd )%elemLvl(level)%bitmask%val(:, nBnds( iBnd ))
             normal = normal + &
              & globBC( iBnd )%elemLvl(level)%normal%val(:, nBnds( iBnd ))
          endif
        end do

        do iBnd = 1, bc_prop%nBCtypes
          if (found(iBnd) .and. globBC(iBnd)%treat_corner ) then
            call append( me = globBC(iBnd)%cornerBC%elemLvl(level)%bitmask,    &
              &          val = corner_bitmask )

            call append( me = globBC( iBnd )%cornerBC%elemLvl( level )%normal, &
              &          val = normal )
          endif
        end do !iBnd
      endif ! count(found) > 1

    end do ! iElem = 1, tree%property(1)%nElems
  end subroutine assignBCList
! **************************************************************************** !


! **************************************************************************** !
  !> This routine normalizes the normal vectors of boundary elements including
  !! the corner elements as well as assigns the corresponding prevailing
  !! direction from the stencil
  subroutine normalizeBC( nBCs, minLevel, maxLevel, globBC, layout, field, &
    &                     comm )
    ! --------------------------------------------------------------------------
    !> number of boundaries
    integer,                        intent(in) :: nBCs, minLevel, maxLevel
    !> boundaries for the elements with bnd property set
    type(glob_boundary_type),    intent(inout) :: globBC(:)
    !> scheme layout
    type( mus_scheme_layout_type ), intent(in) :: layout
    !> field type
    type( mus_field_type ),         intent(in) :: field(:)
    !> mpi communication enviroment with mpi communicator
    integer,                        intent(in) :: comm
    ! --------------------------------------------------------------------------
    integer :: iBnd, iLevel, iElem, iDir, iField, counter, iErr
    real(kind=rk) :: angle, angleMax, oneDeginRad
    real(kind=rk) :: min_nz_comp, max_comp
    integer :: elem_normal(3), bc_normal(3), bc_globNormal(3)
    integer(kind=long_k) :: bc_normlong(3), bc_globNormlong(3)
    logical :: check_angle(nBCs), curved(nBCs)
    ! --------------------------------------------------------------------------
    write(dbgUnit(5),*) 'Normalize BC'
    ! one degree in radian
    oneDegInRad = PI/180.0_rk

    ! check which boundary needs to check for angle
    check_angle = .false.
    curved = .false.
    do iBnd = 1, nBCs
      do iField = 1, size(field)
        if ( field(iField)%bc(iBnd)%nNeighs > 0 ) then
          check_angle(iBnd) = .true.
        end if
        if ( field(iField)%bc(iBnd)%curved ) then
          curved(iBnd) = .true.
        end if
      end do
    end do

    ! normalize all boundary normals
    do iBnd = 1, nBCs
      angleMax = 0.0_rk
      counter = 0
      bc_normal = 0.0_rk
      do iLevel = minLevel, maxLevel
        do iElem = 1, globBC( iBnd )%nElems( iLevel )
          call tem_determine_discreteVector(                               &
            &   globBC( iBnd )%elemLvl( iLevel )%normal%val( :, iElem ),   &
            &   layout%prevailDir, angle )
          ! element normal
          elem_normal = globBC( iBnd )%elemLvl( iLevel )%normal%val( :, iElem )
          ! average boundary normal
          bc_normal = bc_normal + elem_normal
          ! find maximum angle between vectors
          angleMax = max(angleMax, abs(angle))
          ! count number of elements with angle > 1deg
          if (abs(angle) > oneDegInRad) counter = counter + 1
          ! Find the index in the stencil corresponding to the normal
          ! direction
          do iDir = 1, layout%fStencil%QQ
            if( layout%fStencil%cxDir(1, iDir ) == elem_normal(1)         &
              & .and. layout%fStencil%cxDir(2, iDir ) == elem_normal(2)   &
              & .and. layout%fStencil%cxDir(3, iDir ) == elem_normal(3) ) then
              call append (                                                &
                &    me = globBC(iBnd)%elemLvl(iLevel)%normalInd,          &
                &   val = iDir )
            end if
          end do
        end do  ! iElem
      end do  ! iLevel

      bc_normlong = int(bc_normal, kind=long_k)

      ! calculate global normal for given boundary
      call mpi_allreduce( bc_normlong, bc_globNormlong, 3, long_k_mpi, &
        &                 mpi_sum, comm, iErr                          )

      if (maxval(bc_globNormlong) < huge(bc_globNormal(1))) then
        bc_globNormal = int(bc_globNormlong)
      else
        ! Largest Component does not fit into regular integer, rescale to
        ! integer range:
        fit_range: do
          ! Repeat neglecting components until we found a representation that
          ! fits into a default integer.
          min_nz_comp = real( minval( pack( abs(bc_globNormlong),            &
            &                               bc_globNormlong /= 0_long_k ) ), &
            &                 kind=rk                                        )
          max_comp = real( maxval(abs(bc_globNormlong)), kind=rk )
          if ( ceiling(max_comp/min_nz_comp)         &
            & < real(huge(bc_globNormal(1)),kind=rk) ) then
            bc_globNormal = nint( real(bc_globNormlong, kind=rk) / min_nz_comp )
            EXIT fit_range
          else
            ! Smallest nonzero component is too small in comparison to the largest
            ! one to fit into a default integer. To fit into integer the integer
            ! range, we can neglect it.
            where (abs(bc_globNormlong) <= min_nz_comp) bc_globNormlong = 0_long_k
          end if
        end do fit_range
      end if

      ! Compute global normal of this boundary
      call tem_determine_discreteVector( vector        = bc_globNormal,     &
        &                                compareVector = layout%prevailDir, &
        &                                angle         = angle              )

      globBC(iBnd)%normal = bc_globNormal
      do iDir = 1, layout%fStencil%QQ
        if( layout%fStencil%cxDir(1, iDir ) == bc_globNormal(1)         &
          & .and. layout%fStencil%cxDir(2, iDir ) == bc_globNormal(2)   &
          & .and. layout%fStencil%cxDir(3, iDir ) == bc_globNormal(3) ) then
          globBC(iBnd)%normalInd = iDir
        end if
      end do

      ! For straight boundaries use angle from bc_normal to check_angle
      if (.not. curved(iBnd)) angleMax = angle

      write(logUnit(5),'(A,I0,A,3I3,A,I0)') 'BC: '//trim(globBC(iBnd)%label) &
        & //', iBnd: ', iBnd, ', normal: ',  bc_globNormal, ', Ind: ',       &
        & globBC(iBnd)%normalInd

      ! Normal direction is required only to create stencil for
      ! boundary elements which requires neighbors so check for
      ! angle only if any field bc has nNeigh > 0
      if ( check_angle(iBnd)  ) then
        ! convert angle in radian to degree
        angleMax = angleMax / oneDegInRad
        ! if angle is > 1 deg then print warning message
        if (angleMax > 1.0_rk) then
          write(logUnit(6),'(a,i0,a)') ' WARNING: Normal direction of ',&
            &  counter, ' elements of "'// trim(globBC(iBnd)%label)//'" boundary'
          write(logUnit(6),*) 'deviates from stencil discrete vectors.'
          write(logUnit(6),'(a,f5.2)') ' Max. angle between bnd normal and '//&
            & 'discrete vector is:', angleMax
          write(logUnit(6),*) 'This might cause some problems in '// &
            & 'outlet BC.'
          write(logUnit(6),*) 'Solution: use outlet_pab'
        end if
      end if
    end do ! iBnd


    ! normalize normal of elements with multiple boundaries
    do iBnd = 1, nBCs
      do iLevel = minLevel, maxLevel
        if (globBC(iBnd)%treat_corner) then
          do iElem = 1, globBC(iBnd)%cornerBC%nElems( iLevel )
            call tem_determine_discreteVector(                                 &
             &   globBC( iBnd)%cornerBC%elemLvl( iLevel )%                     &
             &                                          normal%val( :, iElem ),&
             &   layout%prevailDir)
            ! Find the index in the stencil corresponding to the normal
            ! direction
            do iDir = 1, layout%fStencil%QQ
              if( layout%fStencil%cxDir(1, iDir) ==                            &
                & globBC(iBnd)%cornerBC%elemLvl(iLevel)%normal%val(1, iElem)   &
                & .and.                                                        &
                & layout%fStencil%cxDir(2, iDir) ==                            &
                & globBC(iBnd)%cornerBC%elemLvl(iLevel)%normal%val(2, iElem)   &
                & .and.                                                        &
                & layout%fStencil%cxDir(3, iDir) ==                            &
                & globBC(iBnd)%cornerBC%elemLvl(iLevel)%normal%val(3, iElem))  &
                & then
              call append (                                                    &
                &      me = globBC(iBnd)%cornerBC%elemLvl(iLevel)%normalInd,   &
                &     val = iDir )
              end if
            end do
          end do ! iElem
        end if ! treat corner
      end do ! iLevel
    end do ! iBnd

  end subroutine normalizeBC
! **************************************************************************** !



! **************************************************************************** !
  !> subroutine to find neighbours of element with individual (for each element)
  !! stencil definitions.
  !! Unique stencil label for boundary stencils are created with boundary label
  !! and stencil%cxDir therefore each stencil is limited to one boundary type
  subroutine mus_build_BCStencils( globBC, bc, prevailDir, prefix, minLevel, &
    &                              maxLevel, stencil_labels, grwStencil )
    ! --------------------------------------------------------------------------
    !> boundaries for the elements with bnd property set
    type(glob_boundary_type), intent(in)             :: globBC
    !> field boundary with boundary neighbor info
    type(boundary_type), intent(inout)               :: bc
    !> scheme layout
    real(kind=rk),                  intent(in)       :: prevailDir(:,:)
    !> field label
    character(len=*), intent(in)                     :: prefix
    !> min and max level
    integer,                        intent(in)       :: minLevel, maxLevel
    !> dynamic array of stencil labels
    type(dyn_labelArray_type), intent(inout)         :: stencil_labels
    !> growing array of stencils
    type(grw_stencilHeaderArray_type), intent(inout) :: grwStencil
    ! --------------------------------------------------------------------------
    integer :: iElem, elemPos, iLevel, iNeigh, stencilPos
    integer :: normal(3)
    type(tem_stencilHeader_type) :: stencil
    logical :: wasAdded
    character(len=labelLen) :: stnLabel
    ! --------------------------------------------------------------------------

    ! only need stencil when BC require neighbors
    if ( bc%nNeighs > 0 ) then

      stencil%label = trim(prefix)//trim(bc%label)
      stencil%useAll = .false.
      ! only true for nNeighs == 2?
      stencil%requireNeighNeigh = .true.
      stencil%QQ     = bc%nNeighs
      stencil%QQN    = bc%nNeighs
      allocate( stencil%cxDir( 3, stencil%QQN ))
      stencil%cxDir = 0

      write(logUnit(5),"(A)")    ' Set BC stencil: '//trim(stencil%label)
      write(logUnit(6),"(A,I0)") '             QQ: ', stencil%QQ

      allocate( stencil%elemLvl( minLevel:maxLevel ) )

      do iLevel = minLevel, maxLevel

        ! allocate field bc stencilPos to size of globBC nElems of that level
        allocate( bc%elemLvl(iLevel)%stencilPos(    globBC%nElems( iLevel ) ) )
        allocate( bc%elemLvl(iLevel)%posInNghElems( globBC%nElems( iLevel ) ) )

        do iElem = 1, globBC%nElems( iLevel )
          ! Get the normal direction of this boundary
          ! only for higher order boundaries.
          ! For curved boundaries, use normal direction of each element
          ! for straight boundaries, use boundary direction
          if (bc%curved) then
            normal(:) = globBC%elemLvl( iLevel )%normal%val( :, iElem )
          else
            normal(:) = globBC%normal(:)
          end if
          ! And store the normal direction into the stencil
          ! QQN stands for the number of neighbors
          ! Here it means, how many neighbors in the normal direction
          do iNeigh = 1, stencil%QQN
            ! ?: what is the purpose of multiplying iStencilElem with normal
            !    direction?
            ! !: We are going along the normal direction for more than
            ! one neighbor in one direction (required for higher order boundary)
            stencil%cxDir(:, iNeigh) = normal * iNeigh
          end do !iNeigh

          ! Get unique stencil label from cxDir
          stnLabel = tem_stencil_getLabelForcxDir(  me         = stencil,   &
            &                                       prevailDir = prevailDir )

          ! extend stencil label with bc label such that
          ! a stencil is restricted strictly to one boundary type
          stnLabel = trim(stencil%label)//'_'//trim(stnLabel)

write(dbgUnit(1),*) 'stencil label ', trim(stnLabel)
          ! append stencil label to dynamic array of stencil labels,
          ! if label is new then wasAdded is true else it is false
          call append( me       = stencil_labels, &
            &          val      = stnlabel,       &
            &          pos      = stencilPos,     &
            &          wasAdded = wasAdded        )
write(dbgUnit(1),*) 'wasAdded ', wasAdded, ' stencilPos ', stencilPos
          if ( wasAdded ) then
            write(logUnit(10),"(A)") "stencil: "//trim(stnlabel)//" was added"
            write(logUnit(10),"(A,I0)") "stencilPos: ", stencilPos
          endif

          ! new stencil, append this stencil to stencil array
          ! and initialize elemLvl and elem in stencil type
          if ( wasAdded ) then
            write(logUnit(10),"(A)") "append stencil to grwStencil"
            call append( me  = grwStencil, &
              &          val = stencil     )

            ! initialize elemLvl and elem list for new stencil
            write(logUnit(10),"(A)") "init grwStencil%elemLvl"
            if ( .not. allocated(grwStencil%val(stencilPos)%elemLvl) ) then
              write(logUnit(10),"(A)") "grwStencil%elemlvl not allocated"
              ! @todo: not allocated for dynamic load balance, needs check!
              allocate( grwStencil%val(stencilPos)%elemLvl( minLevel:maxLevel ) )
            end if
            call init( me     = grwStencil%val(stencilPos)%elemLvl( iLevel ), &
              &        length = 8                                              )
            write(logUnit(10),"(A)") "init grwStencil%elem"
            call init( me     = grwStencil%val(stencilPos)%elem, &
              &        length = 8                                 )

          end if

          ! Get the position of the treeID in the primal mesh
          elemPos = globBC%elemLvl( iLevel )%elem%val( iElem )

          ! append element position in treeID list
          ! to level wise list and level independent list
          call append( me  = grwStencil%val(stencilPos)%elemLvl( iLevel ), &
            &          val = elemPos                                        )
          call append( me  = grwStencil%val(stencilPos)%elem, &
            &          val = elemPos                           )

          ! store stencil position current boundary element
          bc%elemLvl(iLevel)%stencilPos(iElem) = stencilPos
        end do !iElem
      end do !iLevel

      deallocate( stencil%elemLvl )

    end if ! bc%nNeighs > 0


  end subroutine mus_build_BCStencils
! **************************************************************************** !


! **************************************************************************** !
  !> This routine build and append IBM stencils to scheme stencil array
  subroutine mus_build_IBMStencils( globIBM, layout, grwStencil )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_globType ), intent(inout) :: globIBM
    !> scheme stencil layout
    type( mus_scheme_layout_type ), intent(inout) :: layout
    !> contains array of stencils
    type( grw_stencilHeaderArray_type ), intent(inout) :: grwStencil
    ! --------------------------------------------------------------------------
    type( tem_stencilHeader_type ) :: IBM_stencil
    character(len=labelLen) :: label
    integer :: iIBM, stencilPos
    logical :: wasAdded
    ! --------------------------------------------------------------------------
    if (globIBM%nIBMs > 0) write(logUnit(10),*) 'Build and append IBM stencil'

    do iIBM = 1, globIBM%nIBMs
      if (globIBM%IBM(iIBM)%movPredef) then
        call init( me     = IBM_stencil,  &
          &        QQN    = 80,           &
          &        QQ     = 81,           &
          &        useAll = .true.,       &
          &        nDims  = 3,            &
          &        label  = 'd3q81',      &
          &        cxDir  = d3q81_cxDir() )
      else
        call init( me     = IBM_stencil,   &
          &        QQN    = 124,           &
          &        QQ     = 125,           &
          &        useAll = .true.,        &
          &        nDims  = 3,             &
          &        label  = 'd3q125',      &
          &        cxDir  = d3q125_cxDir() )
      end if
      IBM_stencil%requireAll = .true.

      label = IBM_stencil%label

      ! append stencil label to dynamic array of stencil labels,
      ! if label is new then wasAdded is true else it is false
      call append( me       = layout%stencil_labels, &
        &          val      = label,                 &
        &          pos      = stencilPos,            &
        &          wasAdded = wasAdded               )

      ! if stencil is already there then just store the position
      ! else append stencil to stencil array and store its position
      if ( wasAdded ) then
        call append( me  = grwStencil, &
          &          val = IBM_stencil )
      end if
      ! store position of IBM stencil in stencil array
      globIBM%IBM(iIBM)%stencilPos = stencilPos

      write(logUnit(10),*) 'Stored stencil for IBM ', iIBM, ' at position ',   &
        &                 globIBM%IBM(iIbm)%stencilPos
    end do
  end subroutine mus_build_IBMStencils
! **************************************************************************** !


! **************************************************************************** !
  !> set the sendHalo bit for all fluid elements which are send
  !! to remote procs
  subroutine set_sendHaloBits( property, sendBuffer )
    ! --------------------------------------------------------------------------
    type( tem_communication_type ), intent(in) :: sendBuffer
    integer(kind=long_k), intent(inout) :: property(:)
    ! --------------------------------------------------------------------------
    integer :: iProc, iElem, pos
    ! --------------------------------------------------------------------------

    ! loop over the integer send buffers ...
    do iProc = 1, sendBuffer%nProcs
      ! ... and over the corresponding elements ...
      do iElem = 1, sendBuffer%elemPos(iProc)%nVals
        pos = sendBuffer%elemPos(iProc)%val(iElem)
        ! ... and set the property bits
        property(pos) = ibset(property(pos), prp_sendHalo )
      end do
    end do

  end subroutine set_sendHaloBits
! **************************************************************************** !



! **************************************************************************** !
  ! adds ghost elements from coarser with boundary property to bc list
  subroutine mus_update_BcghostElem( minLevel, maxLevel, levelDesc, layout, &
    &                                bc_prop, globBC                        )
    ! --------------------------------------------------------------------------
    !> Min and max level
    integer, intent(in)                        :: minLevel, maxLevel
    !> Level Descriptor
    type( tem_levelDesc_type ),intent(inout)   :: levelDesc(minLevel:maxLevel)
    !> scheme layout
    type( mus_scheme_layout_type ), intent(in) :: layout
    !> boundary information from mesh
    type( tem_bc_prop_type ), intent(in)       :: bc_prop
    !> boundaries for the elements with bnd property set
    type( glob_boundary_type ), intent(inout)  :: globBC(:)
    ! --------------------------------------------------------------------------
    integer :: nGhostFromCoarser!, nGhostFromFiner
    integer :: iLevel, iDir, offset, length
    integer :: weight(layout%fstencil%QQ)
    integer :: iBnd
    ! --------------------------------------------------------------------------
    ! if no boundary leave routine
    if ( bc_prop%nBCtypes == 0 ) return

    ! calculate weights
    ! Get length of currently treated direction and determine
    ! If it is one of the main coordinate axis
    do iDir = 1, layout%fstencil%QQN
      length =   layout%fstencil%cxDir( 1, iDir )**2 &
        &      + layout%fstencil%cxDir( 2, iDir )**2 &
        &      + layout%fstencil%cxDir( 3, iDir )**2
      ! For the normals, we weight the orthogonal directions
      ! of the stencil more heavily
      if( length == 1) then
        weight(iDir) = 4  ! axis orthogonal weighs twice as much as all others
      else if ( length == 2 ) then
        weight(iDir) = 2  ! non-axis orthogonal
      else ! length == 3
        weight(iDir) = 1  ! non-axis orthogonal
      end if
    end do
    ! Level loop
    do iLevel = minLevel, maxLevel
      ! offset for ghost from Coarser Elements
      offset = levelDesc( iLevel )%offset(1, eT_ghostFromCoarser )
      ! nGhost from Coarser Elements
      nGhostFromCoarser = levelDesc( iLevel )%offset(2, eT_ghostFromCoarser ) &
        & - levelDesc( iLevel )%offset(1, eT_ghostFromCoarser )
      ! add all ghost from Coarser ghostelements to bc list
      call  mus_add_BcghostElem(                                             &
        &   levelDesc = levelDesc(iLevel),                                   &
        &   stencil   = layout%fstencil,                                     &
        &   bc_prop   = bc_prop,                                             &
        &   globBC    = globBc,                                              &
        &   nGhosts   = nGhostFromCoarser,                                   &
        &   offset    = offset,                                              &
        &   weight    = weight,                                              &
        &   iLevel    = iLevel                                               )

      do iBnd = 1 , bc_prop%nBCtypes
        ! store number of Bc Elems without GhostBoundaryElems
        globBC( iBnd )%nElems_Fluid(iLevel) = globBC( iBnd )%nElems( ilevel )
        ! replace nElems without GhostBoundaryElems
        ! by nElems incl. GhostFromCoarser boundary elem
        globBC( iBnd )%nElems( ilevel ) =                                      &
         &                        globBC( iBnd )%elemLvl( ilevel )%elem%nVals

        if ( globBC( iBnd )%nElems( ilevel ) -                                 &
         &                     globBC( iBnd )%nElems_Fluid( ilevel ) > 0)  then
          write(logUnit(5),"(A,I5,A,I2,A,A)") 'Added ',                        &
           &globBC( iBnd )%nElems( ilevel ) -                                  &
           &globBC( iBnd )%nElems_Fluid( ilevel ), ' Ghostelements on Level: ',&
           &iLevel,' to Boundary: ', globBC( iBND )%label
        end if

      end do !iBnd
    end do !iLevel

    call set_normalIndGhost(                                  &
      &                 nBCs     = bc_prop%nBCtypes,          &
      &                 minLevel = minLevel,                  &
      &                 maxLevel = maxLevel,                  &
      &                 globBC   = globBC,                    &
      &                 layout   = layout                     )

  contains

  subroutine set_normalIndGhost( nBCs, minLevel, maxLevel, globBC, layout )
    ! --------------------------------------------------------------------------
    !> number of boundaries
    integer,                        intent(in) :: nBCs, minLevel, maxLevel
    !> boundaries for the elements with bnd property set
    type(glob_boundary_type),    intent(inout) :: globBC(:)
    !> scheme layout
    type( mus_scheme_layout_type ), intent(in) :: layout
    ! --------------------------------------------------------------------------
    integer :: iBnd, iLevel, iElem, iDir
    real(kind=rk) :: angle
    ! --------------------------------------------------------------------------

    ! Set normalInd for all Ghost Elements only, since they are allready set
    ! for all Fluid elems
    do iBnd = 1 , nBCs
      do iLevel = minLevel, maxLevel
        do iElem = globBC( iBnd )%nElems_fluid( iLevel ) +1,                  &
                 & globBC( iBnd )%nElems( iLevel )
          call tem_determine_discreteVector(                                  &
            &   globBC( iBnd )%elemLvl( iLevel )%normal%val( :, iElem ),      &
            &   layout%prevailDir, angle )
          ! loop over all directions
          do iDir = 1, layout%fstencil%QQN
            if( layout%fstencil%cxDir(1, iDir ) ==                            &
              & globBC( iBnd )%elemLvl( iLevel )%normal%val( 1, iElem )       &
              & .and.                                                         &
              & layout%fstencil%cxDir(2, iDir ) ==                            &
              & globBC( iBnd )%elemLvl( iLevel )%normal%val( 2, iElem )       &
              & .and. layout%fstencil%cxDir(3, iDir ) ==                      &
              & globBC( iBnd )%elemLvl( iLevel )%normal%val( 3, iElem ) )     &
              & then
              call append (                                                   &
                &    me = globBC( iBnd )%elemLvl( iLevel )%normalInd,         &
                &   val = iDir                                                )
            end if ! iDir
          end do ! iDir
        end do !iElem

        ! truncate growing arrays, since no more Elems will be added
        call truncate ( me = globBC( iBnd )%elemLvl( ilevel )%elem )
        call truncate ( me = globBC( iBnd )%elemLvl( ilevel )%posInBcElemBuf )
        call truncate ( me = globBC( iBnd )%elemLvl( ilevel )%normalInd )
        call truncate ( me = levelDesc( iLevel )%bc_elemBuffer )
        write(dbgUnit(6),*) 'bc_elemBuffer nVals: ', &
          &                  levelDesc( iLevel )%bc_elemBuffer%nVals
      end do ! iLevel
    end do !iBnd
      end subroutine set_normalIndGhost
  end subroutine mus_update_BcghostElem
! **************************************************************************** !


! **************************************************************************** !
  ! This routine appends all ghostelement on this level to BC lists
  ! Run over all the ghostelements with boundary property and check each
  ! direction. Also sets normal, bitmask and qVal for ghostelements.
  subroutine mus_add_BcghostElem(levelDesc, stencil, bc_prop, globBC, nGhosts, &
      &                          offset, weight, iLevel )
    ! --------------------------------------------------------------------------
    !> Level Descriptor for iLevel
    type( tem_levelDesc_type ),intent(inout)    :: levelDesc
    !> stencil
    type(tem_stencilHeader_type),intent(in)     :: stencil
    !> boundary information from mesh
    type( tem_bc_prop_type ), intent(in)        :: bc_prop
    !> boundaries for the elements with bnd property set
    type( glob_boundary_type ), intent(inout)  :: globBC(:)
    integer,intent(in)                   :: iLevel, offset
    integer, intent(in)                  :: weight(stencil%QQN)
    ! --------------------------------------------------------------------------
    integer              :: iElem, iDir, iMeshDir
    integer              :: iPos, nGhosts, ghostpos
    integer              :: Elempos(bc_prop%nBCtypes), normal(3)
    integer(kind=long_k) :: bID, ghostID
    logical              :: added(bc_prop%nBCtypes), bitmask(stencil%QQN)
    real(kind=rk)        :: qval(stencil%QQN)
    ! --------------------------------------------------------------------------
    ! Looping over all Ghosts
    do iElem = 1, nGhosts
      ! set initial values to false and 0
      added      = .false.
      bitmask    = .false.
      Elempos    = 0
      normal     = 0
      qval       = 0.0_rk
      ! position in total of ghost Element
      ghostpos = offset + iElem
      ! TreeID of Ghost
      ghostID = levelDesc%total( ghostpos )
      do iDir = 1, stencil%QQN
        if (levelDesc%neigh(1)%nghElems(iDir,ghostPos) < 0) then
          bID = ABS( levelDesc%neigh(1)%nghElems(iDir,ghostPos) )
          iMeshDir = stencil%map( iDir )
          if ( .not. globBC( bID )%isWall ) then
            if ( globBC( bID )%hasQVal ) then
              write(logUnit(0),*) 'WARNING: q-values are not implemented for '   &
                &                 // 'ghost elements. Might cause problems for ' &
                &                 // 'higher order wall boundaries'
            end if

            if ( .not. added(bID) ) then
              ! storem elemPos in globBc
              call append( me  = globBC( bID )%elemLvl( iLevel )%elem, &
                &          val = ghostpos,                             &
                &          pos = Elempos( bID )                        )
              ! store elemPos in bc_elemBuffer
              call append( me  = levelDesc%bc_elemBuffer, &
                &          val = ghostpos,                &
                &          pos = iPos                     )
              ! store the position in the bc_elemBuffer in the elemBuffer
              call append ( me  = globBC( bID )%elemLvl( iLevel ) &
                  &                            %posInBcElemBuf,   &
                  &         val = iPos                            )
              !append empty bitmask
              call append ( me  = globBC( bID )%elemLvl( ilevel )%bitmask, &
                &           val = bitmask                                  )
              !append empty normal
              call append ( me = globBC( bID )%elemLvl( ilevel )%normal, &
                &           val = normal                                 )
              ! check for qVal
              if( globBC( bID )%hasQVal ) then
                !append qVal
                call append( me  = globBC( bID )%elemLvl( ilevel )%qVal, &
                  &          val = qVal  )
              end if ! haQval
              ! update normal
               globBC( bID )%elemLvl( ilevel )%normal%val(:, Elempos( bID ) ) &
                 &  =  globBC( bID )%elemLvl( ilevel )                        &
                 &                  %normal%val(:,Elempos( bID ))             &
                 &     -  weight(iDir) * stencil%cxDir( :, iDir )

              !update bitmask
              globBC( bID )%elemLvl( ilevel )%bitmask%                   &
               & val( stencil%cxDirInv(iDir), Elempos( bID ) ) = .true.
              added(bID) = .true.
            else
              ! update bitmask
              globBC( bID )%elemLvl( ilevel )%bitmask%                   &
                & val( stencil%cxDirInv(iDir), Elempos( bID ) ) = .true.
              ! update normal
              globBC( bID )%elemLvl( ilevel )%normal%val(:, Elempos( bID ) ) &
                &  =  globBC( bID )%elemLvl( ilevel )                        &
                &                  %normal%val(:,Elempos( bID ))             &
                &     -  weight(iDir) * stencil%cxDir( :, iDir )
            end if ! added
            ! update qVal
            if( globBC( bID )%hasQVal ) then
              ! Set the qvalue from the boundary array
              ! q-values stored in outgoing direction
              ! qVal are initiated with 0.5 for ALL directions
              ! since we check the bitmask and only take the qVal for
              ! bitmask == true, there is no problem with that.
              globBC(bID)%elemLvl(iLevel)%qVal%val(iDir, ElemPos( bID ) ) &
               & = 0.5_rk
            end if ! hasQVal
          end if ! isWall
        end if ! BC
      end do !iDir
    end do !iElem
  end subroutine mus_add_BcghostElem
! **************************************************************************** !

end module mus_construction_module
! **************************************************************************** !
