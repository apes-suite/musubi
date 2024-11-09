! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013, 2015-2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
! ***************************************************************************** !
!> author: Simon Zimny
!! scheme_layout module, providing the scheme_layout datatype and the
!! functionality to read the lua files and to set the predefined stencils.
!!
!!
!! # Scheme Layout
!!
!! The layout includes all information about the stencil, the inverse
!! directions, the prevailing directions and the weights.
!!
!! In MUSUBI one can choose among the following predefined layouts:
!!
!! -  'predefined_d3q19':
!!    The D3Q19 model. The advection relaxation kernel for this model is highly
!!    optimized.
!! -  'predefined_d3q7':
!!    The more simple D3Q7 model. The advection relaxation kernel for this
!!    model is highly optimized.
!! -  'Flekkoy':
!!    The Flekkoy model consists of a D3Q6 stencil and is used for passive
!!    scalar transport.
!!
!! In MUSUBI it is also possible to define a new layout in the lua file. This
!! feature is implemented to test new layouts. In case one wants to run
!! multiple simulations using this layout, it is highly recommended to
!! implement a new kernel (the following files have to be extended:
!! [[mus_flow_module]] and related modules in source/compute.
!!
!! Defining a new Layout
!! To define a new layout one has to set `layout = 'new_stencil'`. To define the
!! stencil the following information has to be provided:
!!
!!  -  the number of offsets (QQ)
!!  -  the offsets (disc_vel)
!!  -  the weights (weight)
!!  -  the inverse directions (inv_dir) as integers pointing on the belonging
!!     offsets
!!  -  the prevailing directions (prev_dir) as reals giving the distance in the
!!     different directions
!!
!! The following example shows how to define the standard D3Q19 as a new
!! layout:
!!
!!```lua
!! -- Scheme
!! scheme = {
!!  label = 'test',
!!  layout = 'new_stencil',
!!  -- Initial condition
!!  initial_condition = {
!!                       density = 1.0,
!!                       velocityX = 0.0,
!!                       velocityY = 0.0,
!!                       velocityZ = 0.0 },
!!  -- Boundary conditions
!!  boundary_condition = {
!!  { label = 'wall',
!!    kind = 'velocity_bounceback',
!!    velocityX = 0.03, velocityY = 0.0, velocityZ = 0.0 },
!!  { label = 'wall',
!!    kind = 'wall',
!!    velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }
!!  },
!!  fluid = { relaxation_scheme = 'BGK',
!!            omega = 1.8,
!!            rho0 = 1.0 },
!!  -- Defining the new Stencil
!!  stencil = {
!!             QQ = 19,
!!             disc_vel = {
!!              {-1,0,0},{0,-1,0},{0,0,-1},{1,0,0},{0,1,0},{0,0,1},{0,-1,-1},{0,-1,1},{0,1,-1},
!!              {0,1,1},{-1,0,-1},{1,0,-1},{-1,0,1},{1,0,1},{-1,-1,0},{-1,1,0},{1,-1,0},{1,1,0},
!!              {0,0,0}
!!             },
!!             weight = {
!!              (1.0/18.0),(1.0/18.0),(1.0/18.0),(1.0/18.0),(1.0/18.0),(1.0/18.0),(1.0/36.0),(1.0/36.0),(1.0/36.0),
!!              (1.0/36.0),(1.0/36.0),(1.0/36.0),(1.0/36.0),(1.0/36.0),(1.0/36.0),(1.0/36.0),(1.0/36.0),(1.0/36.0),
!!              (1.0/3.0)
!!             },
!!             inv_dir = {
!!                         4,  5,  6,  1,  2,  3, 10,  9,  8,
!!                         7, 14, 13, 12, 11, 18, 17, 16, 15,
!!                        19
!!             },
!!             prev_dir = {
!!              {1.,0.,0.},{0.,1.,0.},{0.,0.,1.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.},
!!              {0.,0.5*math.sqrt(2.),0.5*math.sqrt(2.)},{0.,0.5*math.sqrt(2.),0.5*math.sqrt(2.)},
!!              {0.,0.5*math.sqrt(2.),0.5*math.sqrt(2.)},{0.,0.5*math.sqrt(2.),0.5*math.sqrt(2.)},
!!              {0.5*math.sqrt(2.),0.,0.5*math.sqrt(2.)},{0.5*math.sqrt(2.),0.,0.5*math.sqrt(2.)},
!!              {0.5*math.sqrt(2.),0.,0.5*math.sqrt(2.)},{0.5*math.sqrt(2.),0.,0.5*math.sqrt(2.)},
!!              {0.5*math.sqrt(2.),0.5*math.sqrt(2.),0.},{0.5*math.sqrt(2.),0.5*math.sqrt(2.),0.},
!!              {0.5*math.sqrt(2.),0.5*math.sqrt(2.),0.},{0.5*math.sqrt(2.),0.5*math.sqrt(2.),0.}
!!             }
!!  }
!! }
!!```
!!
module mus_scheme_layout_module

  ! include treelm modules
  use mpi
  use env_module,            only: rk, long_k, labelLen
  use tem_aux_module,        only: tem_abort
  use tem_stencil_module,    only: tem_stencilHeader_type, tem_loadStencil,   &
    &                              tem_identify_inverseDirections,            &
    &                              tem_identify_prevailDirections,            &
    &                              tem_stencil_zeroPos,                       &
    &                              tem_create_stencil,                        &
    &                              grw_stencilHeaderArray_type, init,         &
    &                              destroy, truncate
  use tem_logging_module,    only: logUnit
  use tem_param_module,      only: div1_3, div1_18, div1_36, div1_54, div1_216,&
    &                              div2_27, div8_27, div1_4, div1_8, div1_6, &
    &                              div2_3, div1_9, div4_9, div1_24, div1_2
  use tem_dyn_array_module,  only: dyn_labelArray_type, init, append, destroy
  use tem_grow_array_module, only: grw_intArray_type, init, append, destroy,   &
    &                              truncate, grw_longArray_type
  use tem_comm_env_module,   only: tem_comm_env_type
  use tem_tools_module,      only: tem_horizontalSpacer

  ! include aotus modules
  use aotus_module,     only: flu_state, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_val

  ! include musubi modules
  use mus_moments_type_module, only: mus_moment_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none

  private

  public :: mus_scheme_layout_type
  public :: mus_load_newLayout
  public :: mus_init_layout
  public :: mus_finalize_layout
  public :: mus_destroy_stencil
  public :: mus_define_layout
  public :: mus_define_d3q19
  public :: mus_define_d3q27
  public :: mus_define_d3q7
  public :: mus_define_d3q6
  public :: mus_define_d2q9
  public :: mus_define_d1q3
  public :: mus_weights_out
  public :: mus_set_weights_d1q3
  public :: mus_set_weights_d2q9
  public :: mus_set_weights_d3q6
  public :: mus_set_weights_d3q7
  public :: mus_set_weights_d3q19
  public :: mus_set_weights_d3q27

  !> data structure containing all information related to the
  !! compute stencil. Several stencils can be defined.
  !! [[mus_moments_module]] Moments are directly related to the
  !! stencil layout and are therefore defined here
  type mus_scheme_layout_type

    !> fluid stencil same as stencil(1)
    type( tem_stencilHeader_type ) :: fStencil

    !> number of stencils used in this scheme
    integer :: nStencils

    !> Temporary growing array of stencil
    !! It is copied to stencil(:) and destroyed, where this is destroyed?
    type(grw_stencilHeaderArray_type) :: grwStencil

    !> The list of stencil types, the stencils for the individual schemes
    !! is ordered as follows:
    !! --------------------------------------------------------------------------------------
    !! | flSt | bcSt1_field1 ... bcStN_field1 ... bcSt1_fieldM ... bcStN_fieldM | addSt ... |
    !! --------------------------------------------------------------------------------------
    !! Unique stencil label for boundary stencils are created with boundary label
    !! and stencil%cxDir therefore each stencil is limited to one boundary type
    type( tem_stencilHeader_type ), allocatable :: stencil(:)

    !> dynamic array of labels created from stencil directions to create unique
    !! growing array of grwStencil
    type(dyn_labelArray_type) :: stencil_labels

    !> position of fluid stencil in grwStencil
    integer :: fStencil_pos

    !> The weights for the different discrete velocities
    real(kind=rk), allocatable :: weight(:)

    !> Lattice speed of sound for fStencil
    !! $\sum_i (weight_i*cx_i*cx_i) = c_s^2 I
    real(kind=rk) :: cs

    !> Prevailing directions
    real(kind=rk), allocatable :: prevailDir(:,:)

    !> Moment space definition
    type(mus_moment_type) :: moment

    !> New stencil definition loaded from config
    logical :: new_stencil = .false.

    !> derive quantities that depends on the layout such as velocity, pdf_eq, etc..
    type(mus_scheme_derived_quantities_type) :: quantities

  end type mus_scheme_layout_type

contains

! ****************************************************************************** !
  !> load a new stencil definition from the lua file
  !!
  !! - a label
  !! - the stencil (predefined or new)
  !! - for a new defined stencil ( weights, inverse directions, prevailing
  !!   directions )
  !!
  subroutine mus_load_newLayout( me, parent_handle, conf )
    ! ---------------------------------------------------------------------------
    integer, intent(in), optional :: parent_handle
    type(mus_scheme_layout_type) :: me
    type(flu_State) :: conf
    ! ---------------------------------------------------------------------------
    ! defining local variables
    integer :: stencil_handle
    integer :: weight_handle
    integer :: nWeights, iWeight
    integer :: invDir_handle
    integer :: iInvDir, nInvDirs
    integer :: prevDir_handle
    integer :: iPrevDir, nPrevDirs
    integer :: prevDirCoo_handle
    integer :: iPrevDirCoo
    integer :: iError
    ! ---------------------------------------------------------------------------

    if( present( parent_handle))then
      ! open stencil table
      call aot_table_open( L       = conf,                                     &
        &                  parent  = parent_handle,                            &
        &                  thandle = stencil_handle,                           &
        &                  key     = 'stencil' )
    else
      ! open stencil table
      call aot_table_open( L=conf, thandle=stencil_handle, key = 'stencil' )
    end if

    ! If a new stencil is defined, the flag is set to true to
    ! differentiate between predefined and new stencil in initialization
    if (stencil_handle > 0) then
      me%new_stencil = .true.

      write(logUnit(1),*)'Reading the new Layout...'
      me%fStencil%label = 'new_stencil'

      ! Get the number of dimensions in which the stencil is defined
      call aot_get_val( L       = conf,                                        &
        &               thandle = stencil_handle,                              &
        &               val     = me%fStencil%nDims,                          &
        &               ErrCode = iError,                                      &
        &               key     = 'nDims',                                     &
        &               default = -1 )

      ! load the stencil
      call tem_loadStencil( stencil       = me%fStencil,                      &
        &                   parent_handle = stencil_handle,                    &
        &                   conf          = conf )

      ! load additional information like weights, inv. directions, etc.
      call aot_table_open( L       = conf,                                     &
        &                  parent  = stencil_handle,                           &
        &                  thandle = weight_handle,                            &
        &                  key     = 'weight' )
      nWeights = aot_table_length( L = conf, thandle = weight_handle )
      if ( nWeights == me%fStencil%QQ ) then
        allocate( me%weight( nWeights ))

        do iWeight = 1, nWeights
          call aot_get_val( L       = conf,                                    &
            &               thandle = weight_handle,                           &
            &               val     = me%weight( iWeight ),                    &
            &               ErrCode = iError,                                  &
            &               pos     = iWeight )
        end do
        call aot_table_close( L=conf, thandle=weight_handle )
      else
        write(logUnit(1),*) 'The number of defined weights does not '//        &
          &                 'match the number of discrete velocities.'
        write(logUnit(1),*) nWeights, ' vs. ', me%fStencil%QQ
        call tem_abort()
      end if
      call aot_table_open( L       = conf,                                     &
        &                  parent  = stencil_handle,                           &
        &                  thandle = invDir_handle,                            &
        &                  key     = 'inv_dir' )
      nInvDirs = aot_table_length( L = conf, thandle = invDir_handle )
      if ( nInvDirs == me%fStencil%QQ ) then
        allocate( me%fStencil%cxDirInv( nInvDirs ))

        do iInvDir = 1, nInvDirs
          call aot_get_val( L       = conf,                                    &
            &               thandle = invDir_handle,                           &
            &               val     = me%fStencil%cxDirInv( iInvDir ),       &
            &               ErrCode = iError,                                  &
            &               pos     = iInvDir )
        end do
        call aot_table_close( L = conf, thandle = invDir_handle )
      else
        call tem_identify_inverseDirections( me%fStencil%cxDirInv,           &
          &                                  me%fStencil%cxDir )
      end if

      call aot_table_open( L       = conf,                                     &
        &                  parent  = stencil_handle,                           &
        &                  thandle = prevDir_handle,                           &
        &                  key     = 'prev_dir' )
      nPrevDirs = aot_table_length( L = conf, thandle = prevDir_handle )
      if(nPrevDirs > 0)then
        allocate( me%prevailDir( 3, nPrevDirs ))
        do iPrevDir = 1, nPrevDirs
          call aot_table_open( L       = conf,                                 &
            &                  parent  = prevDir_handle,                       &
            &                  thandle = prevDirCoo_handle,                    &
            &                  pos     = iPrevDir )
          do iPrevDirCoo = 1, 3
            call aot_get_val( L       = conf,                                  &
              &               thandle = prevDirCoo_handle,                     &
              &               val     = me%prevailDir( iPrevDirCoo, iPrevDir ),&
              &               ErrCode = iError,                                &
              &               pos     = iPrevDirCoo )
          end do
          call aot_table_close( L = conf, thandle = prevDirCoo_handle )
        end do
        call aot_table_close( L=conf, thandle=prevDir_handle )
      else
        call tem_identify_prevailDirections( me%prevailDir,                    &
          &                                  me%fStencil%cxDir )
      end if

      me%fStencil%restPosition = tem_stencil_zeroPos( me%fStencil )

      write(logUnit(1),*) 'A new stencil has been defined successfully.'
      if( me%fStencil%nDims <= 0 ) then
        write(logUnit(1),*)'Error: number of dimensions is not given for '//   &
          &            'the stencil'
        write(logUnit(1),*)'Please define in the configuration file as '
        write(logUnit(1),*)'  stencil = { nDims = 3, weights = ... }'
        call tem_abort
      end if
    end if

  end subroutine mus_load_newLayout
! ****************************************************************************** !


! ****************************************************************************** !
  !> Dump the weights in lua format.
  !!
  !! The style is:
  !! `weights = {w1, w2, ... , wN}`
  !!
  subroutine mus_weights_out( me, conf )
    ! ---------------------------------------------------------------------------
    !> weights
    real(kind=rk), intent(in) :: me(:)
    type(aot_out_type) :: conf
    ! ---------------------------------------------------------------------------
    call aot_out_val( put_conf     = conf,                                     &
      &               vname        = 'weight',                                 &
      &               val          = me,                                       &
      &               max_per_line = 2 )

  end subroutine mus_weights_out
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initialize growing array of stencils
  !!
  subroutine mus_init_layout( layout )
    ! ---------------------------------------------------------------------------
    !> musubi schemes layout type
    type( mus_scheme_layout_type ), intent(inout) :: layout
    ! ---------------------------------------------------------------------------
    write(logUnit(1),*) 'Initializing stencil array: '
    ! Init is explicitly called here to clear the eventuellay used arrays
    call init( me = layout%grwStencil )
    call init( me = layout%stencil_labels)

  end subroutine mus_init_layout
! ****************************************************************************** !


! ****************************************************************************** !
  !> This routine finialize grwStencil by truncating stencil elem arrays and
  !! set stencil%nElems
  subroutine mus_finalize_layout( layout, nElemsInTree, minLevel, maxLevel, proc )
    ! ---------------------------------------------------------------------------
    !> scheme layout
    type( mus_scheme_layout_type ), intent(inout) :: layout
    !> fluid tree from mesh
    integer,                           intent(in) :: nElemsInTree
    !> min and max level
    integer,                           intent(in) :: minLevel, maxLevel
    !> mpi communication type
    type(tem_comm_env_type),           intent(in) :: proc
    ! ---------------------------------------------------------------------------
    integer :: iStencil, iLevel, iProc
    integer, allocatable :: nElems_totalStencil(:)
    integer :: nStencils_all(proc%comm_size), nStencils_total
    integer :: offset(proc%comm_size)
    integer :: iErr, charType, stencilPos
    logical :: wasAdded
    type(dyn_labelArray_type) :: stencil_labels
    type(grw_longArray_type) :: nElems
    character(len=labelLen), allocatable :: stencil_labels_total(:)
    ! ---------------------------------------------------------------------------
    ! truncate stencil array
    call truncate( me = layout%grwStencil )
    write(logUnit(5),'(a)') ' Finalizing stencils...'

    allocate( layout%stencil(layout%grwStencil%nVals) )
    layout%nStencils = layout%grwStencil%nVals

    ! copy stencil from growing array to allocatable array
    layout%stencil(:) = layout%grwStencil%val(1:layout%nStencils)
    ! update fluid stencil
    layout%fStencil = layout%stencil(1)

    ! detroy growing array
    call destroy( me = layout%grwStencil )

    do iStencil = 1, layout%nStencils
      ! if stencil use all = .true. use all elements in tree
      ! else use stencil elements
      if ( layout%stencil(iStencil)%useAll ) then
        layout%stencil(iStencil)%nElems = nElemsInTree
      else
        call truncate( me = layout%stencil(iStencil)%elem )
        layout%stencil(iStencil)%nElems = layout%stencil(iStencil)%elem%nVals

        do iLevel = minLevel, maxLevel
          call truncate( me = layout%stencil(iStencil)%elemLvl(iLevel) )
        end do
      end if
    end do

    ! gather number of stencils across process
    call MPI_GATHER( layout%nStencils, 1, MPI_INTEGER, nStencils_all, 1, &
      &              MPI_INTEGER, proc%root, proc%comm, iErr )

    if (proc%rank == proc%root) then
      ! total number of stencils
      nStencils_total = sum(nStencils_all)
      ! number of elements in each stencil in each process
      allocate( nElems_totalStencil(nStencils_total) )
      ! stencil labels
      allocate( stencil_labels_total(nStencils_total) )
      ! displacement to gather nElems vector
      offset(1) = 0
      do iProc = 2, proc%comm_size
        offset(iProc) = offset(iProc-1) + nStencils_all(iProc-1)
      end do
    else
      allocate( nElems_totalStencil(0) )
      allocate( stencil_labels_total(0) )
    end if

    ! gather number of elements per stencil on each process
    call MPI_GATHERV( layout%stencil(:)%nElems, layout%nStencils, MPI_INTEGER, &
      & nElems_totalStencil, nStencils_all, offset, MPI_INTEGER, proc%root,    &
      & proc%comm, iErr )

    ! cast character array type to gather stencil labels
    call MPI_TYPE_CONTIGUOUS( labelLen, MPI_CHARACTER, charType, iErr )
    call MPI_TYPE_COMMIT( charType, iErr )

    ! gather stencil label to create unique stencil list
    call MPI_GATHERV( layout%stencil_labels%val(:), layout%nStencils, charType,&
      & stencil_labels_total, nStencils_all, offset, charType, proc%root,      &
      & proc%comm, iErr )

    if (proc%rank == proc%root) then
      ! create unique stencil labels and reduce nElems on each stencil
      do iProc = 1, proc%comm_size
        do iStencil = 1, nStencils_all(iProc)
          call append( me       = stencil_labels,                              &
            &          val      = stencil_labels_total(offset(iProc)+iStencil),&
            &          pos      = stencilPos,                                  &
            &          wasAdded = wasAdded                                     )
          if (wasAdded) then
            call append( me  = nElems,                                         &
              &          val = int(nElems_totalStencil(offset(iProc)+iStencil),&
              &                    kind=long_k)                                )
          else
            nElems%val(stencilPos) = nElems%val(stencilPos) +                  &
              & int(nElems_totalStencil(offset(iProc)+iStencil), kind=long_k)
          end if
        end do
      end do

      write(logUnit(5),'(a,i0)') ' Total Number of stencils: ', stencil_labels%nVals
      do iStencil = 1, stencil_labels%nVals
        write(logUnit(5),'(a,i2,2a)') ' iStencil: ', iStencil, &
          & ', label: ', trim(stencil_labels%val(iStencil))
        write(logUnit(5),'(a,i0)') ' nElems: ', nElems%val(iStencil)
      end do
      call destroy(stencil_labels)
      call destroy(nElems)
    end if

  end subroutine mus_finalize_layout
! ****************************************************************************** !


! ****************************************************************************** !
  !> Destroy the stencil
  !!
  subroutine mus_destroy_stencil( stencil )
    ! ---------------------------------------------------------------------------
    !>musubi schemes stencil type
    type( tem_stencilHeader_type ), allocatable, intent(out) :: stencil(:)
    ! ---------------------------------------------------------------------------
    write(logUnit(5),*)'Deallocating stencil layout...'
    if( allocated(stencil) ) deallocate( stencil )

  end subroutine mus_destroy_stencil
! ****************************************************************************** !


! **************************************************************************** !
  !> This routine defines layout for predefined stencils
  subroutine mus_define_layout( layout, stencilName, nElems )
    ! ------------------------------------------------------------------------ !
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    !> Name of the stencil to create
    character(len=*), intent(in) :: stencilName
    !> number of elements use this layout
    integer, intent(in) :: nElems
    ! ------------------------------------------------------------------------ !
    ! create fStencil
    select case ( trim( stencilName ) )
    case ( 'd2q9' )
      call mus_define_d2q9( layout = layout, &
        &                   nElems = nElems  )
    case ( 'd2q5' )
      call mus_define_d2q5( layout = layout, &
        &                   nElems = nElems  )
    case ( 'd3q19' )
      call mus_define_d3q19( layout = layout, &
        &                    nElems = nElems  )
    case ( 'd3q13' )
      call mus_define_d3q13( layout = layout, &
        &                    nElems = nElems  )
    case ( 'd3q27' )
      call mus_define_d3q27( layout = layout, &
        &                    nElems = nElems  )
    case ( 'd3q7' )
      call mus_define_d3q7( layout = layout, &
        &                   nElems = nElems  )
    case ( 'd3q6', 'flekkoy' )
      call mus_define_d3q6( layout = layout, &
        &                   nElems = nElems  )
    case ( 'd1q3' )
      call mus_define_d1q3( layout = layout, &
        &                   nElems = nElems  )
    case default
      call tem_horizontalSpacer(fUnit = logUnit(1))
      write(logUnit(1),*) 'The chosen scheme layout is not available. '
      write(logUnit(1),*) 'STOPPING!'
      write(logUnit(1),*) 'Choose one of:  '
      write(logUnit(1),*) '  - d2q9 '
      write(logUnit(1),*) '  - d2q5 '
      write(logUnit(1),*) '  - d3q13 '
      write(logUnit(1),*) '  - d3q19 '
      write(logUnit(1),*) '  - d3q27 '
      write(logUnit(1),*) '  - d3q7  '
      write(logUnit(1),*) '  - d3q6  '
      write(logUnit(1),*) '  - d1q3  '
      call tem_abort()
    end select

    ! if ( .not. allocated( layout%weight ) &
      ! &  allocate( layout%weight( layout%fStencil%QQ ) )

  end subroutine mus_define_layout
! **************************************************************************** !

! ****************************************************************************** !
  !> This subroutine sets the parameters for the predefined d3q13 stencil.
  !!
  subroutine mus_define_d3q13( layout, nElems )
    ! ---------------------------------------------------------------------------
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    !> number of elements use this layout
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------

    ! fill the tem_stencilHeader_type
    call tem_create_stencil( layout%fStencil, 'd3q13' )

    ! setting the prevailing directions (normalized major directions without
    ! sign)
    call tem_identify_prevailDirections( layout%prevailDir,      &
      &                                  layout%fStencil%cxDir )

    if ( .not. allocated(layout%weight) ) &
      & allocate( layout%weight( layout%fStencil%QQ ))
    call mus_set_weights_d3q13( layout%weight )

    layout%fStencil%useAll = .true.
    layout%fStencil%nElems = nElems

    layout%fStencil%restPosition = tem_stencil_zeroPos( layout%fStencil )

    ! speed of sound
    layout%cs = mus_calculate_speed_of_sound( layout )

  end subroutine mus_define_d3q13
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine mus_set_weights_d3q13( weights )
    real( kind=rk ) :: weights(13)

    weights(1:12) = div1_24
    weights(13)   = div1_2
  end subroutine mus_set_weights_d3q13
! ****************************************************************************** !

! ****************************************************************************** !
  !> This subroutine sets the parameters for the predefined d3q19 stencil.
  !!
  subroutine mus_define_d3q19( layout, nElems )
    ! ---------------------------------------------------------------------------
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    !> number of elements use this layout
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------

    ! fill the tem_stencilHeader_type
    call tem_create_stencil( layout%fStencil, 'd3q19' )

    ! setting the prevailing directions (normalized major directions without
    ! sign)
    call tem_identify_prevailDirections( layout%prevailDir,      &
      &                                  layout%fStencil%cxDir )

    if ( .not. allocated(layout%weight) ) &
      & allocate( layout%weight( layout%fStencil%QQ ))
    call mus_set_weights_d3q19( layout%weight )

    layout%fStencil%useAll = .true.
    layout%fStencil%nElems = nElems

    layout%fStencil%restPosition = tem_stencil_zeroPos( layout%fStencil )

    ! speed of sound
    layout%cs = mus_calculate_speed_of_sound( layout )

  end subroutine mus_define_d3q19
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine mus_set_weights_d3q19( weights )
    real( kind=rk ) :: weights(19)

    weights(1:6 ) = div1_18
    weights(7:18) = div1_36
    weights(19)   = div1_3
  end subroutine mus_set_weights_d3q19
! ****************************************************************************** !

! ****************************************************************************** !
  !> This subroutine sets the parameters for the predefined d3q27 stencil.
  !!
  subroutine mus_define_d3q27( layout, nElems )
    ! ---------------------------------------------------------------------------
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    !> number of elements use this layout
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------

    ! fill the tem_stencilHeader_type
    call tem_create_stencil( layout%fStencil, 'd3q27' )

    ! setting the prevailing directions (normalized major directions without
    ! sign)
    call tem_identify_prevailDirections( layout%prevailDir,                    &
      &                                  layout%fStencil%cxDir )

    if ( .not. allocated(layout%weight) ) &
      & allocate( layout%weight( layout%fStencil%QQ ))
    call mus_set_weights_d3q27( layout%weight )

    layout%fStencil%useAll = .true.
    layout%fStencil%nElems = nElems

    layout%fStencil%restPosition = tem_stencil_zeroPos( layout%fStencil )

    ! speed of sound
    layout%cs = mus_calculate_speed_of_sound( layout )

  end subroutine mus_define_d3q27
! ****************************************************************************** !

  subroutine mus_set_weights_d3q27( weights )
    real( kind=rk ) :: weights(27)

    weights( 1: 6) = div2_27
    weights( 7:18) = div1_54
    weights(19:26) = div1_216
    weights(27)    = div8_27
  end subroutine mus_set_weights_d3q27

! ****************************************************************************** !
  !> This subroutine sets the parameters for the predefined d3q7 stencil.
  !!
  subroutine mus_define_d3q7( layout, nElems )
    ! ---------------------------------------------------------------------------
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    !> number of elements use this layout
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------

    ! fill the tem_stencilHeader_type
    call tem_create_stencil( layout%fStencil, 'd3q7' )

    ! setting the prevailing directions (normalized major directions without
    ! sign)
    call tem_identify_prevailDirections( layout%prevailDir,                    &
      &                                  layout%fStencil%cxDir )

    allocate( layout%weight( layout%fStencil%QQ ))
    call mus_set_weights_d3q7( layout%weight )

    layout%fStencil%useAll = .true.
    layout%fStencil%nElems = nElems

    layout%fStencil%restPosition = tem_stencil_zeroPos( layout%fStencil )

    ! speed of sound
    layout%cs = mus_calculate_speed_of_sound( layout )

  end subroutine mus_define_d3q7
! ****************************************************************************** !

  subroutine mus_set_weights_d3q7( weights )
    real( kind=rk ) :: weights(7)

    weights(1:6) = div1_8
    weights(7)   = div1_4
  end subroutine mus_set_weights_d3q7

! ****************************************************************************** !
  !> This subroutine sets the parameters for the predefined d3q6
  !! layout%fStencil, used by the Flekkoy model of passive scalar transport.
  !!
  subroutine mus_define_d3q6( layout, nElems )
    ! ---------------------------------------------------------------------------
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    !> number of elements use this layout
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------

    ! fill the tem_stencilHeader_type
    call tem_create_stencil( layout%fStencil, 'd3q6' )

    ! setting the prevailing directions (normalized major directions without
    ! sign)
    call tem_identify_prevailDirections( layout%prevailDir,                    &
      &                                  layout%fStencil%cxDir )

    allocate( layout%weight(layout%fStencil%QQ) )
    call mus_set_weights_d3q6( layout%weight )

    layout%fStencil%useAll = .true.
    layout%fStencil%nElems = nElems

    layout%fStencil%restPosition = tem_stencil_zeroPos( layout%fStencil )

    ! speed of sound
    layout%cs = mus_calculate_speed_of_sound( layout )

  end subroutine mus_define_d3q6
! ****************************************************************************** !

  subroutine mus_set_weights_d3q6( weights )
    real( kind=rk ) :: weights(6)

    weights(1:6) = div1_6
  end subroutine mus_set_weights_d3q6

! ****************************************************************************** !
  !> This subroutine sets the parameters for the predefined d2q9 stencil.
  !!
  subroutine mus_define_d2q9( layout, nElems )
    ! ---------------------------------------------------------------------------
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    !> number of elements use this layout
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------

    ! fill the tem_stencilHeader_type
    call tem_create_stencil( layout%fStencil, 'd2q9' )

    ! setting the prevailing directions (normalized major directions without
    ! sign)
    call tem_identify_prevailDirections( layout%prevailDir,    &
      &                                  layout%fStencil%cxDir )

    if ( .not. allocated(layout%weight) ) &
      & allocate( layout%weight( layout%fStencil%QQ ))
    call mus_set_weights_d2q9( layout%weight )

    layout%fStencil%useAll = .true.
    layout%fStencil%nElems = nElems

    layout%fStencil%restPosition = tem_stencil_zeroPos( layout%fStencil )

    ! speed of sound
    layout%cs = mus_calculate_speed_of_sound( layout )

  end subroutine mus_define_d2q9
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine mus_set_weights_d2q9( weights )
    real( kind=rk ) :: weights(9)

    weights(1:4) = div1_9
    weights(5:8) = div1_36
    weights(9)   = div4_9
  end subroutine mus_set_weights_d2q9
! ****************************************************************************** !

! ****************************************************************************** !
  !> This subroutine sets the parameters for the predefined d2q5 stencil.
  !!
  subroutine mus_define_d2q5( layout, nElems )
    ! ---------------------------------------------------------------------------
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    !> number of elements use this layout
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------

    ! fill the tem_stencilHeader_type
    call tem_create_stencil( layout%fStencil, 'd2q5' )

    ! setting the prevailing directions (normalized major directions without
    ! sign)
    call tem_identify_prevailDirections( layout%prevailDir,    &
      &                                  layout%fStencil%cxDir )

    if ( .not. allocated(layout%weight) ) &
      & allocate( layout%weight( layout%fStencil%QQ ))
    call mus_set_weights_d2q5( layout%weight )

    layout%fStencil%useAll = .true.
    layout%fStencil%nElems = nElems

    layout%fStencil%restPosition = tem_stencil_zeroPos( layout%fStencil )

    ! speed of sound
    layout%cs = mus_calculate_speed_of_sound( layout )

  end subroutine mus_define_d2q5
! ****************************************************************************** !

! ****************************************************************************** !
  subroutine mus_set_weights_d2q5( weights )
    real( kind=rk ) :: weights(5)

    weights(1:4) = div1_4
    weights(5)   = 0.0_rk
  end subroutine mus_set_weights_d2q5
! ****************************************************************************** !

! ****************************************************************************** !
  !> This subroutine sets the parameters for the predefined d2q9 stencil.
  !!
  subroutine mus_define_d1q3( layout, nElems )
    ! ---------------------------------------------------------------------------
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    !> number of elements use this layout
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------

    ! fill the tem_stencilHeader_type
    call tem_create_stencil( layout%fStencil, 'd1q3' )

    ! setting the prevailing directions (normalized major directions without
    ! sign)
    call tem_identify_prevailDirections( layout%prevailDir,                    &
      &                                  layout%fStencil%cxDir )

    allocate( layout%weight( layout%fStencil%QQ ))
    call mus_set_weights_d1q3( layout%weight )

    layout%fStencil%useAll = .true.
    layout%fStencil%nElems = nElems

    layout%fStencil%restPosition = tem_stencil_zeroPos( layout%fStencil )

    ! speed of sound
    layout%cs = mus_calculate_speed_of_sound( layout )

  end subroutine mus_define_d1q3
! ****************************************************************************** !

  subroutine mus_set_weights_d1q3( weights )
    real( kind=rk ) :: weights(3)

    ! Weights from [1]  E. M. Viggen, The lattice Boltzmann method in
    ! acoustics, Scandinavian Symposium on Physical Acoustics, pp. 1–5, Mar.
    ! 2010.
    weights(1:2) = div1_6
    weights(3)   = div2_3
  end subroutine mus_set_weights_d1q3
! ****************************************************************************** !


! ****************************************************************************** !
  !> Calculate lattice speed of sound for given stencil
  function mus_calculate_speed_of_sound(layout) result(c_sound)
    ! ---------------------------------------------------------------------------
    !> scheme layout for pdf state
    type( mus_scheme_layout_type ), intent(inout)   :: layout
    real(kind=rk) :: c_sound
    ! ---------------------------------------------------------------------------
    c_sound = sqrt(sum( layout%fStencil%cxcx(1,:) * layout%weight(:) ))

  end function mus_calculate_speed_of_sound
! ****************************************************************************** !

end module mus_scheme_layout_module
! ****************************************************************************** !

