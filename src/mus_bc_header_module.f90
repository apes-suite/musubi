! Copyright (c) 2012-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2015, 2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2015-2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2018,2020 Jana Gericke <jana.gericke@uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! **************************************************************************** !
!> This module contains boundary type definitions, different boundary treatment
!! types, abtract interface and routine to load boundary table from config file
!!
!! The boundary conditions are stored in boundary type definitions for each
!! boundary which is defined in the main lua file.
!!```lua
!! boundary_condition = {  }
!!```
!!
!! The definitions there must equal the boundary definition within the mesh on
!! disk, in terms of the amount of boundary conditions and the labels.
!! A detailed description on the implementation details are given in
!! [[tem_bc_module]] "boundary condition implementation"
!!
!! This module will be used as a header to use boundary definitions in other
!! modules.
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
module mus_bc_header_module

  ! include treelm modules
  use env_module,               only: rk, LabelLen
  use tem_aux_module,           only: tem_abort
  use tem_float_module,         only: operator(.feq.)
  use tem_property_module,      only: prp_solid
  use tem_bc_module,            only: tem_bc_state_type, tem_load_bc_state
  use tem_bc_header_module,     only: tem_load_bc_header, tem_bc_header_type
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_dyn_array_module,     only: init, destroy, dyn_intArray_type, &
    &                                 dyn_longArray_type
  use tem_grow_array_module,    only: init, destroy, grw_int2darray_type, &
    &                                 grw_logical2darray_type,            &
    &                                 grw_intarray_type,                  &
    &                                 grw_real2darray_type
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_logging_module,       only: logUnit
  use tem_debug_module,         only: dbgUnit
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_geometry_module,      only: tem_BaryOfId,tem_ElemSizeLevel
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_stringKeyValuePair_module, only: init, truncate,                &
    &                                      grw_stringKeyValuePairArray_type

  ! include aotus modules
  use aotus_module,     only: flu_State, aoterr_Fatal, aoterr_NonExistent, &
    &                         aoterr_WrongType, aot_top_get_val, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close

  ! include musubi modules
  use mus_field_prop_module,     only: mus_field_prop_type
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_derVarPos_module,      only: mus_derVarPos_type
  ! use mus_param_module,          only: mus_param_type
  use mus_physics_module,        only: mus_physics_type, mus_convertFac_type
  use mus_mixture_module,        only: mus_mixture_type
  use mus_pdf_module,            only: pdf_data_type
  use mus_turb_wallFunc_module,  only: mus_turb_wallFunc_type, &
    &                                  mus_load_turb_wallFunc

  implicit none

  private

  public :: boundary_type
  public :: glob_boundary_type
  public :: mus_load_bc
  public :: mus_init_bc_elems
  public :: check_solid_in_bc
  public :: rearrange_bc_elems
  public :: debug_glob_boundary_type
  public :: mus_set_posInNghElems
  public :: mus_set_bcLinks
  public :: mus_set_bouzidi
  public :: mus_alloc_bouzidi
  public :: mus_alloc_fieldBC
  public :: mus_set_inletUbb
  public :: mus_set_inletBfl
  public :: mus_set_nonEqExpol
  public :: mus_set_outletExpol
  public :: mus_fieldBC_cleanup

  !> information needed for moments based boundary condition
  type bc_moments_type
    !> Number of unknown state links to update
    integer :: nUnKnownPdfs
    !> known moments position in moments vector
    integer, allocatable :: knownMom_pos(:)
    !> inverse matrix of unknown pdf matrix
    real(kind=rk), allocatable :: unKnownPdfs_MatInv(:,:)

  end type bc_moments_type


  !> Level wise boundary elements information
  type bc_elems_type

    !> Positions in levelDesc total list
    !! Its purpose is to get the treeID of BC elements
    !! size is globBC%nElems
    !! to use: levelDesc(iLevel)%total( globBC%elemLvl(iLevel)%elem%val(iElem) )
    type( dyn_intArray_type ) :: elem

    !> Position of this boundary in the bc_elemBuffer
    !! bc_elemBuffer is growing array
    !! initial size: globBC%nElems
    !! It is initiated in mus_init_bc_elems of mus_bc_header_module
    !! Only non-wall BC needs elemBuffer
    type ( grw_intarray_type ) :: posInBcElemBuf

    !> Normal vector for each boundary element pointing into the domain.
    !! normal is a growing array with size: normal%val(3, nElems)
    type ( grw_int2darray_type ) :: normal

    !> which is the index in the stencil corresponding to the normal vector?
    type ( grw_intarray_type ) :: normalInd

    !> bit mask for each node holding directions which have to be updated.
    !! The bitmask points into the incoming direction into the flow domain,
    !!  which actually we want to update
    !! * For PUSH, we write to the incoming position,
    !!   as the kernel reads it from there without propagation.
    !! * For PULL, we need to write to the inverse direction, as the kernel
    !!   performs a bounce back before reading it.
    !!   However, this bounced back direction actually comes from the
    !!   non-existent boundary element and would point into the incoming
    !!   direction, so the value has to be treated and set as if it points
    !!   really into the incoming direction.
    !! bitmask is a growing array, the values are in bitmask%val(:,:)
    !! 1st index size is QQN since center is not treated for bcElems
    type ( grw_logical2darray_type ) :: bitmask

    !> The q-Values for the exact wall positions
    !! Its size: QQN, nElemsBC
    !! It is allocated in routine: allocateBCList
    !! assigned in routine: assignBCList
    type( grw_real2darray_type ) ::  qVal

  end type bc_elems_type


  !> contains information needed to treat corner nodes or nodes intersected
  !! by more than one boundary type
  type mus_bc_corner_type
    !> number of local boundary elements per level
    integer, allocatable :: nElems(:)
    !> boundary neighbor in different refinement level
    type( bc_elems_type ), allocatable :: elemLvl(:)
    !> boundary id in each direction for each corner node
    !! size layout%stencil%QQ, nElems
    integer, allocatable :: bcid(:,:)
  end type mus_bc_corner_type

  !> Wall bouzidi data for one level
  type bc_wall_bouzidi_type

    !> size is links(iLevel)%nVals
    real(kind=rk), allocatable :: qVal(:)
    real(kind=rk), allocatable ::  cIn(:)
    real(kind=rk), allocatable :: cOut(:)
    real(kind=rk), allocatable :: cNgh(:)

    !> size is links(iLevel)%nVals
    !! set in routine: mus_set_bouzidi
    integer, allocatable ::  inPos(:) ! position in bcBuffer
    integer, allocatable :: outPos(:) ! position in bcBuffer
    integer, allocatable :: nghPos(:) ! position in computeNeighBuf

  end type bc_wall_bouzidi_type


  !> Provides bouzidi coefficients and q-values needed for link-wise
  !! implementation of certain inlet boundary conditions.
  type bc_inlet_bouzidi_type

    !> size is links(iLevel)%nVals
    real(kind=rk), allocatable :: qVal(:) ! q-Values
    real(kind=rk), allocatable ::  cIn(:) ! coefficient for incoming directions
    real(kind=rk), allocatable :: cOut(:) ! coefficient for outgoing directions
    real(kind=rk), allocatable :: cNgh(:) ! coefficient for neighbours

    !> size is links(iLevel)%nVals
    !! set in routine: mus_set_bouzidi
    integer, allocatable ::  inPos(:) ! position in bcBuffer
    integer, allocatable :: outPos(:) ! position in bcBuffer
    integer, allocatable :: nghPos(:) ! position in computeNeighBuf

    !> size is links(iLevel)%nVals
    real(kind=rk), allocatable ::    cVel(:) ! coefficient for eqPlus

    !> size is links(iLevel)%nVals
    !! set in routine: mus_set_bouzidi
    integer, allocatable :: iDir(:) ! incoming direction from boundary

    integer, allocatable ::   posInBuffer(:) !< used for position in bcBuffer
  end type bc_inlet_bouzidi_type


  !> Provides coefficients needed for link-wise implementation of
  !! nonEquilibrium extrapolation scheme for boundary conditions.
  type bc_nonEqExpol_type

    !> size is links(iLevel)%nVals
    real(kind=rk), allocatable ::      c_w(:) ! coefficient for surface
    real(kind=rk), allocatable ::      c_f(:) ! coefficient for fluid
    real(kind=rk), allocatable ::     c_ff(:) ! coefficient for overnext fluid
    real(kind=rk), allocatable ::  c_neq_f(:) ! nonEq coefficient for fluid
    real(kind=rk), allocatable :: c_neq_ff(:) ! nonEq coeff. for overnext fluid

    !> size is links(iLevel)%nVals
    !! set in routine: mus_set_bouzidi
    integer, allocatable ::          iDir(:) ! incoming direction
    integer, allocatable ::   posInBuffer(:) ! used for position in bcBuffer
    ! Position of fluid neigh in incoming direction in neighBuffer
    integer, allocatable :: posInNeighBuf(:)
    integer, allocatable :: posInBCelems(:)

  end type bc_nonEqExpol_type


  !> Provides needed positions for link-wise implementation of outlet boundary
  !> condition outlet_expol
  type bc_outlet_type

    !> size is links(iLevel)%nVals
    integer, allocatable ::  statePos(:) ! position of state
    integer, allocatable ::  iElem(:) ! position of state

  end type bc_outlet_type


  !> variable definition for non-reflective type of boundary conditions\n
  !! These boundary condition is taken from the paper:
  !!```lua
  !!   S. Izquierdo and N. Fueyo,
  !!   "Characteristic nonreflecting boundary
  !!   conditions for open boundaries in lattice Boltzmann methods,"
  !!   Physical Review E, vol. 78, no. 46707, 2008.
  !!```
  type bc_nrbc_type

    !> specific heat ratio
    !! appears in eq 3: cs = sqrt( kappa * R * T )
    real(kind=rk) :: kappa

    !> from eq 17: k = sigma * ( 1 - Ma^2 ) * cs / L
    real(kind=rk) :: K_mod

    !> used in eq 17
    real(kind=rk) :: cs_mod

    !> used in eq 17
    real(kind=rk) :: sigma

    !> length between inlet and outlet, represented by L in the paper
    real(kind=rk) :: lodi_length

    !> Lattice Mach number characterised by lattice flow velocity
    real(kind=rk) :: Ma_L

  end type bc_nrbc_type


  !> Contains variables for black-box membrane boundary condition
  type bc_blackbox_membrane_type
    !> type of membrance (AEM or CEM)
    character(len=labelLen) :: label
    !> transference number defining transfer property
    !! of ionic species on the membrance
    !! KM @todo check transNr constrain
    !!          sumOfIonicSpecies(transNr) = 1
    real(kind=rk) :: transNr
    !> osmotic solvent transport coefficient
    real(kind=rk) :: solvTrans_coeff
  end type bc_blackbox_membrane_type


  !> contains all possible bc state variables needed for boundary
  type bc_states_type

    !> Velocity at boundary
    type(tem_bc_state_type) :: velocity

    !> mass flow rate at boundary node for mfr (mass flow rate)
    !! boundary condition
    type(tem_bc_state_type) :: massFlowRate

    !> Pressure at boundary
    type(tem_bc_state_type) :: pressure

    !> mole fraction for species boundary
    type(tem_bc_state_type) :: moleFrac

    !> molar flux for species (c_i*u_i)
    type(tem_bc_state_type) :: moleFlux

    !> mole density for species (c_i)
    type(tem_bc_state_type) :: moleDens

    !> molar diffusion flux for species
    type(tem_bc_state_type) :: moleDiff_flux

    !> probability density function for boundary
    type(tem_bc_state_type) :: pdf

    !> potential for Poisson boundary
    type(tem_bc_state_type) :: potential

    !> surface charge density for Poisson Neumann boundary
    type(tem_bc_state_type) :: surChargeDens

  end type bc_states_type


  !> Simple array type which include array and its size
  type array_type

    !> array
    integer, allocatable :: val(:)

    !> size
    integer :: nVals

  end type array_type


  !> Information about boundary elements
  !!
  !! This derived type encapsulates the definition of a boundary
  !! with all necessary parameters and data fields
  type glob_boundary_type

    !> name of this boundary
    character(len=LabelLen) :: label

    !> elements list of this boundary in different refinement level
    !! size is minLevel -> maxLevel
    type( bc_elems_type ), allocatable :: elemLvl(:)

    !> Collect corner BC i.e elements which are interesected by multiple
    !! boundaries only for moments based BC
    logical :: treat_corner

    !> list corner boundary elements
    type( mus_bc_corner_type ) :: cornerBC

    !> number of local (this process) boundary elements per level
    !! including ghostFromCoarser.
    !! GhostFromFiner are interpolated but ghostFromCoarser needs
    !! correct boundary value in fine level sub-iteration.
    integer, allocatable :: nElems(:)

    !> number of local (this process) boundary elements per level
    !! wihtout Ghost elements
    integer, allocatable :: nElems_Fluid(:)

    !> number of total (total processes) boundary elements per level
    integer, allocatable :: nElems_totalLevel(:)

    !> number of local boundary elements
    integer :: nElems_local = 0

    !> number of total boundary elements
    integer :: nElems_total = 0

    !> whether this boundary requires qVal
    logical :: hasQVal = .false.

    !> has qVal initialized with default qVal = 0.5
    logical :: qValInitialized = .false.

    !> if this is a wall or symmetry BC (implicit treated during streaming)
    logical :: isWall = .false.

    !> Average normal direction of this boundary along layout%prevailDir
    !! Normal in elemLvl is normal direction per element and this
    !! is normal per boundary
    integer :: normal(3)

    !> which is the index in the stencil corresponding to the normal vector?
    integer :: normalInd
  end type glob_boundary_type


  !> Level wise information about neighbor of boundary elements
  !! for a single field
  type bc_neigh_type

    !> Neighbor position in level-wise state array
    !! size is nNeighs, nElems
    !! allocated in routine: setFieldBCNeigh, mus_construction_module
    integer, allocatable :: posInState(:,:)

    !> Pre-collision state values of neighbors on normal direction on next
    !! time step
    !! size is nNeighs, nElems*stencil%QQ
    !! It is allocated in routine update_BClists
    real(kind=rk), allocatable :: neighBufferPre_nNext(:,:)

    !> Pre-collision state values of neighbors on normal direction on previous
    !! time step
    !! size is nNeighs, nElems*stencil%QQ
    !! It is allocated in routine update_BClists
    real(kind=rk), allocatable :: neighBufferPre(:,:)

    !> Post-collision state values of neighbors on normal direction
    !! size is nNeighs, nElems*stencil%QQ
    !! It is allocated in routine update_BClists
    !! It is filled in fill_neighBuffer
    real(kind=rk), allocatable :: neighBufferPost(:,:)

    !> Post-collision state values of neighbors on all directions at
    !! current time step.
    !! Its also a Pre-collision state values of an element next to boundary
    !! at next time step
    !! size is  nElems * computeStencil%QQ
    !! It is allocated in routine update_BClists (mus_construction_module),
    !! filled in routine: fill_computeNeighBuf (mus_bc_general_module).
    !! it use AOS layout always!
    real(kind=rk), allocatable :: computeNeighBuf(:)

  end type bc_neigh_type


  !> Contains field specific element information for each level
  type bc_field_elems_type
    !> Position of stencil of each element in array of stencils
    !! in layout.
    !! Unique stencil label for boundary stencils are created with boundary label
    !! and stencil%cxDir therefore each stencil is limited to one boundary type
    !! Dimension: globBC%nElems(iLevel)
    !! allocated in mus_build_BCStencils (mus_construction_module)
    integer, allocatable :: stencilPos(:)

    !> LODI array for the NRBC
    !! Dimension: 4, 3, globBC%nElems(iLevel)
    !! allocated in init_nrbc routine
    real(kind=rk), allocatable :: lodi(:,:,:)

    !> defined moments position in moments matrix and
    !! inverse of unknown pdf matrix
    !! Dimension: bc%nElems(iLevel)
    type( bc_moments_type ), allocatable :: moments(:)

    !> Position in LevelDesc%neigh%NghElems
    !! size is globBC%nElems(iLevel)
    !! allocated in mus_build_BCStencils routine, mus_construction_module
    integer, allocatable :: posInNghElems(:)

  end type bc_field_elems_type


  !> Boundary treatment data type.
  !! It include boundary treatment for only ONE field.
  !! For the same boundary elements, each field may have a different treatment.
  !! It includes function pointers and neighbor elements required by this
  !! boundary treatment
  type boundary_type

    !> Boundary ID
    integer :: bc_id

    !> kind of this boundary
    character(len=LabelLen) :: BC_kind

    !> name of this boundary
    character(len=LabelLen) :: label

    !> number of neighbors fluid element needed for extrapolation to
    !! get data at boundary element
    integer :: nNeighs =  0

    !> Curved boundaries are treated with qValues and extrapolation of
    !! nonEq from 2nd neighbor
    logical :: curved = .false.

    !> reference point for boundary moment calculation
    real(kind=rk) :: bndMomRefPnt(3)

    !> Level wise neighbor information of boundary elements
    !! Each field may have a different handlement on the same boundary
    !! Thus requires a different neighborhood
    !! allocated in mus_alloc_fieldBC
    type(bc_neigh_type), allocatable :: neigh(:)

    !> Level wise element information of boundary elements
    !! like stencil position in array of stencils,
    !! moments type for momentBC and LODI array for NRBC
    !! allocated in mus_alloc_fieldBC
    type(bc_field_elems_type), allocatable :: elemLvl(:)

    !> function pointer to boundary routine to call
    procedure(boundaryRoutine), pointer :: fnct => null()

    !> function pointer to compute bndForce on wall
    procedure(mus_proc_calcBndForce), pointer :: calcBndForce => null()

    !> musubi bc state variables
    type(bc_states_type) :: BC_states

    !> Dictionary of boundary variables with
    !! varDict%val()%key is the name of boundary variable and
    !! varDict%val()%value is the name of spacetime function variable
    type(grw_stringKeyValuePairArray_type) :: varDict

    !> order of boundary which determines number of neighbor need for this boundary
    integer :: order

    !> standard q-value for higher order boundaries without any definition
    !! in seeder
    real(kind=rk) :: qVal

    !> Non-reflective boundary type
    type(bc_nrbc_type) :: nrbc

    !> slip wall factor
    real(kind=rk) :: slip_fac

    !> black box membrane boundary
    type(bc_blackbox_membrane_type) :: blackbox_mem

    !> whether this boundary needs neighbors on compute stencil
    logical :: useComputeNeigh = .false.

    !> Is true if this boundary requires boundary variables
    !! on boundary surface linkwise.
    !! Required to compute real coordinates points on boundary surface
    logical :: evalBcVar_link = .false.

    !> Is true only if boundary requires neighbuf Pre-collision values from
    !! next time step
    logical :: requireNeighBufPre_nNext = .false.

    !> Is true only if boundary requires neighbuf Pre-collision values from
    !! current time step
    logical :: requireNeighBufPre = .false.

    !> Is true only if boundary requires neighbuf Post-collision values
    logical :: requireNeighBufPost = .false.

    !> Collect corner BC i.e elements which are interesected by multiple
    !! boundaries only for moments based BC
    logical :: treat_corner = .false.

    !> Level-wise links that are to be updated
    !! size is minLevel:maxLevel
    !! allocated in mus_alloc_fieldBC
    type(array_type), allocatable :: links(:)

    !> Level wise wall bouzidi data
    !! size is minLevel:maxLevel
    type(bc_wall_bouzidi_type), allocatable :: bouzidi(:)

    !> Level wise inletUbbQVal data
    !! size is minLevel:maxLevel
    type(bc_inlet_bouzidi_type), allocatable :: inletUbbQVal(:)

    !> Level wise inletBfl data
    !! size is minLevel:maxLevel
    type(bc_inlet_bouzidi_type), allocatable :: inletBFL(:)

    !> Level wise nonEquilibrium extrapolation data for link-wise implementation
    !! size is minLevel:maxLevel
    type(bc_nonEqExpol_type), allocatable :: nonEqExpol(:)

    !> Level wise outletExpol data
    !! size is minLevel:maxLevel
    type(bc_outlet_type), allocatable :: outletExpol(:)

    !> Turbulent wall model type contains function pointers to
    !! compute friction and stream-wise velocity. Also stores friction
    !! velocity and turbulent viscosity on boundary elements
    type(mus_turb_wallFunc_type) :: turbwallFunc
  end type boundary_type


  !> Interface definition for the boundary condition subroutines
  abstract interface
    !>  This abstract interface defines the interface for boundary routines. All
    !! boundary routines that need to be called via [[boundary_type:fnct]]
    !! function pointer need to implement this interface.
    !! In case of new boundary routines, mark them with a comment reffering to
    !! [[boundary_type:fnct]] in order to find the routine in case the interface
    !! needs to be changed.
    subroutine boundaryRoutine( me, state, bcBuffer, globBC, levelDesc, tree, &
      &                         nSize, iLevel, sim_time, neigh, layout,       &
      &                         fieldProp, varPos, nScalars, varSys,          &
      &                         derVarPos, physics, iField, mixture           )
      import :: boundary_type, rk, mus_field_prop_type,                    &
        &       mus_scheme_layout_type, tem_time_type, tem_leveldesc_type, &
        &       treelmesh_type, mus_derVarPos_type, tem_varSys_type,       &
        &       glob_boundary_type, mus_physics_type, mus_mixture_type
      !> field boundary type
      class( boundary_type ) :: me
      !> State array
      real(kind=rk), intent(inout) :: state(:)
      !> size of state array ( in terms of elements )
      integer, intent(in) :: nSize
      !> state values of boundary elements of all fields of iLevel
      real(kind=rk),intent(in) :: bcBuffer(:)
      !> global treelm mesh
      type( treelmesh_type ), intent(in) ::tree
      !> iLevel descriptor
      type(tem_levelDesc_type), intent(in) :: levelDesc
      !> level which invokes boundary
      integer,intent(in) :: iLevel
      !> global time information
      type( tem_time_type ), intent(in)  :: sim_time
      !> global parameters
      ! type(mus_param_type),intent(in) :: params
      integer,intent(in)  :: neigh(:)  !< connectivity array
      !> scheme layout
      type( mus_scheme_layout_type ),intent(in) :: layout
      !> Fluid property
      type( mus_field_prop_type ), intent(in)  :: fieldProp
      !> pointer to field variable in the state vector
      integer, intent(in) :: varPos(:)
      !> number of Scalars in the scheme var system
      integer, intent(in) :: nScalars
      !> scheme variable system
      type( tem_varSys_type ), intent(in) :: varSys
      !> position of derived quantities in varsys
      type( mus_derVarPos_type ), intent(in) :: derVarPos
      !> scheme global boundary type
      type( glob_boundary_type ), intent(in) :: globBC
      !> scheme global boundary type
      type( mus_physics_type ), intent(in) :: physics
      !> current field
      integer, intent(in) :: iField
      !> mixture info
      type(mus_mixture_type), intent(in) :: mixture
    end subroutine boundaryRoutine

    !>  This abstract interface defines the interface for bndForce calculation.
    !! [[boundary_type:calcBndForce]] in order to find the routine in case the
    !! interface needs to be changed.
    subroutine mus_proc_calcBndForce( me, bndForce, bndMoment, posInBndID, &
      &                               globBC, currState, levelDesc, nSize, &
      &                               iLevel, neigh, layout, nScalars,     &
      &                               phyConvFac)
      import :: rk, mus_scheme_layout_type, tem_leveldesc_type, &
        &       glob_boundary_type, boundary_type, mus_convertFac_type
      !> field boundary type
      class( boundary_type ), intent(in) :: me
      !> bndForce to fill
      real(kind=rk), intent(inout) :: bndForce(:,:)
      !> bndMom to fill
      real(kind=rk), intent(inout) :: bndMoment(:,:)
      ! position of boundary element in boundary%bcID
      integer, intent(in) :: posInBndID(:)
      !> scheme global boundary type
      type( glob_boundary_type ), intent(in) :: globBC
      !> current state array to access post-collision values
      real(kind=rk), intent(in) :: currState(:)
      !> size of state array ( in terms of elements )
      integer, intent(in) :: nSize
      !> iLevel descriptor
      type(tem_levelDesc_type), intent(in) :: levelDesc
      !> level which invokes boundary
      integer,intent(in) :: iLevel
      !> global parameters
      ! type(mus_param_type),intent(in) :: params
      integer,intent(in)  :: neigh(:)  !< connectivity array
      !> scheme layout
      type( mus_scheme_layout_type ),intent(in) :: layout
      !> number of Scalars in the scheme var system
      integer, intent(in) :: nScalars
      !> physics conversion factor
      type(mus_convertFac_type), intent(in) :: phyConvFac
    end subroutine mus_proc_calcBndForce
  end interface


contains


  ! ------------------------------------------------------------------------ !
  !> Read in the boundary conditions from the LUA parameter file.\n
  !!
  subroutine mus_load_bc( me, bc_prop, conf, parent, varSys, stencil )
    ! -------------------------------------------------------------------- !
    !> field boundary type
    type( boundary_type ), allocatable   :: me(:)
    !> boundary data from mesh
    type( tem_bc_prop_type ), intent(in) :: bc_prop
    !> lua flu state
    type( flu_state )               :: conf
    !> bc parent handle
    integer, intent(in), optional   :: parent
    !> Global variable system required to append annoymous boundary variables
    type(tem_varSys_type), intent(inout) :: varSys
    !> stencil information
    type(tem_stencilHeader_type) :: stencil
    ! -------------------------------------------------------------------- !
    integer :: bc_handle
    integer :: sub_handle
    integer :: iError, vError(3)
    integer :: errFatal(3)
    ! loop variable
    integer :: iBC, myBCID
    ! boundary header type
    type( tem_bc_header_type ) :: bc_header
    ! -------------------------------------------------------------------- !

    errFatal = aoterr_Fatal

    ! if fluid informations in scheme table parentHandle /= 0
    ! load bounary labels  in lua config, and match them with those of
    ! loaded from mesh (information in tem_bc_prop_type)
    ! If a BC in bc_prop can not be match by lua, then code abort
    ! If a BC in lua can not be match by bc_prop, it gives a warning.
    call tem_load_bc_header( me           = bc_header, &
      &                      conf         = conf,      &
      &                      parentHandle = parent,    &
      &                      BC_prop      = bc_prop    )
    call aot_table_open( L       = conf,                &
      &                  parent  = parent,              &
      &                  thandle = bc_handle,           &
      &                  key     = 'boundary_condition' )

    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),'(A,I0)') ' Loading boundary definition for nBCs: ', &
      &                        bc_header%nBCs

    ! load tables inside the boundary table for bc kind and their state
    ! variables boundaries defined in lua are match with bc label read from
    ! mesh using bc_header
    allocate( me( bc_header%nBCs ))

    do iBC = 1, bc_header%nBCs
      myBCID = bc_header%BC_ID(iBC)
      call aot_table_open( L       = conf,       &
        &                  parent  = bc_handle,  &
        &                  thandle = sub_handle, &
        &                  pos     = iBC         )
      ! Skip BC definitions that are not present in the mesh.
      if ( (sub_handle /= 0) .and. (myBCID > 0) ) then
        write(logUnit(1),"(A,I0)") '   Boundary: ',iBC
        ! matching the boundary label from mesh to boundary label read for lua
        ! file bc_header%BC_ID(iBC) returns the position of boundary label iBC
        ! in the mesh.
        me( myBCID )%bc_id = myBCID
        me( myBCID )%label   = bc_header%label(iBC)
        me( myBCID )%BC_kind = bc_header%BC_kind(iBC)

        ! Initialize varDict for current boundary
        call init( me = me(myBCID)%varDict )

        ! Load qVal for the boundary, if certain boundary kind requires
        ! default qVal then they are set in the select case
        call aot_get_val( L       = conf,              &
          &               thandle = sub_handle,        &
          &               val     = me( myBCID )%qVal, &
          &               ErrCode = iError,            &
          &               key     = 'qVal',            &
          &               default = 0.5_rk             )

        ! Flag to treat curved boundaries. Especially nonEq expolation
        ! boundary kind. For curved boundaries, neighbors are taken along
        ! the link direction
        call aot_get_val( L       = conf,                &
          &               thandle = sub_handle,          &
          &               val     = me( myBCID )%curved, &
          &               ErrCode = iError,              &
          &               key     = 'curved',            &
          &               default = .false.              )
        write(logUnit(1),*)'   Curved: ', me( myBCID )%curved

        ! Load reference point for boundary moment calculation
        call aot_get_val( L       = conf,                         &
          &               thandle = sub_handle,                   &
          &               val     = me( myBCID )%bndMomRefPnt,    &
          &               key     = 'moment_ref_point',           &
          &               default = (/ 0.0_rk, 0.0_rk, 0.0_rk /), &
          &               ErrCode = vError )
        if ( any(btest(vError, errFatal)) ) then
          write(logUnit(1), *) 'Error reading moment_ref_point for boundary ' &
            &                // trim(me( myBCID )%label)
          call tem_abort()
        end if

        select case( trim( bc_header%BC_kind(iBC) ))
        case( 'wall' )
          write(logUnit(1),*) '    Simple bounce back treatment for '// &
            &  'no-slip wall '// trim(me( myBCID )%label)
          me( myBCID )%nNeighs = 0

        case( 'symmetry' )
          write(logUnit(1),*) '    Symmetry treatment for slip wall ' &
            & // trim(me( myBCID )%label)
          me( myBCID )%nNeighs = 0

        case( 'spc_bb_wall' )
          write(logUnit(1),*) '    Species bounce back treatment for '// &
            &  'no-slip wall '// trim(me( myBCID )%label)
          me( myBCID )%nNeighs = 1

        case( 'moments_wall', 'spc_moments_wall' )
          write(logUnit(1),*) '    Moments wall treatment for '// &
            &  'no-slip wall '// trim(me( myBCID )%label)
          me( myBCID )%nNeighs = 0

        case( 'spc_outflow', 'spc_solvent_outflow' )
          write(logUnit(1),*) '    '//trim(bc_header%BC_kind(iBC))        &
            &                 //' treatment for '//trim(me( myBCID )%label)
          if (.not. me( myBCID )%curved) me( myBCID )%qVal = 0.0_rk
          me( myBCID )%nNeighs = 2

        case( 'spc_moments_outflow' )
          write(logUnit(1),*) '    Spc moments outflow treatment for '// &
            &  trim(me( myBCID )%label)
          me( myBCID )%nNeighs = 1

        case( 'wall_multiReflection','wall_mr')
          write(logUnit(1),*) '    Wall multi reflection BC does NOT work!'
          call tem_abort()

          bc_header%BC_kind(iBC) = 'wall_multiReflection'
          me( myBCID )%BC_kind = bc_header%BC_kind(iBC)
          write(logUnit(1),*) '    Wall multi reflection for no-slip wall ' &
            &             //trim(me( myBCID )%label)
          me( myBCID )%nNeighs = 2

        case( 'wall_linearInterpolation',    &
          &   'wall_libb', 'wall_libb_link', &
          &   'wall_bouzidi'                 )
          bc_header%BC_kind(iBC) = 'wall_libb'
          me( myBCID )%BC_kind = bc_header%BC_kind(iBC)
          write(logUnit(1),*) '    Wall LIBB treatment: '// &
            &                 trim(me( myBCID )%label)
          me( myBCID )%nNeighs = 0
          me( myBCID )%useComputeNeigh = .true.

        case( 'slip_wall', 'spc_slip_wall' )
          write(logUnit(1),*) '   '//trim( bc_header%BC_kind(iBC) ) // &
            &    ' treatment for '//trim(me( myBCID )%label)
          call aot_get_val( L       = conf,                  &
            &               thandle = sub_handle,            &
            &               val     = me( myBCID )%slip_fac, &
            &               ErrCode = iError,                &
            &               key     = 'fac',                 &
            &               default = 1.0_rk                 )
          me( myBCID )%nNeighs = 0

        case( 'turbulent_wall', 'turbulent_wall_noneq_expol', &
          & 'turbulent_wall_eq')
          write(logUnit(1),*) '    Wall function ' //                   &
            &                 trim( bc_header%BC_kind(iBC) ) //         &
            &                 ' treatment for '//trim( me(myBCID)%label )
          me( myBCID )%qVal = 0.5_rk
          ! Use upto 2 neighs because turbulent wall requires precollision
          ! from 1st neigh which uses fetch.
          ! Also, this BC works only if 1st neigh is local fluid not halo.
          me( myBCID )%nNeighs = 2
          me( myBCID )%requireNeighBufPre_nNext = .true.
          select case( trim( bc_header%BC_kind(iBC) ))
          case ('turbulent_wall')
            me( myBCID )%useComputeNeigh = .true.
          case ('turbulent_wall_noneq_expol', 'turbulent_wall_eq')
            ! Use compute neigh buffer to compute density locally for
            ! straight boundary
            if(.not. me( myBCID )%curved) me( myBCID )%useComputeNeigh = .true.
          end select
          ! load wall model type
          call mus_load_turb_wallFunc(me     = me(myBCID)%turbwallFunc, &
            &                         conf   = conf,                    &
            &                         parent = sub_handle               )

        case( 'bc_pdf' )
          write(logUnit(1),*) '   '//trim( bc_header%BC_kind(iBC) ) //     &
            &    ' treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states%pdf, &
            &                     state_name = 'pdf',                    &
            &                     nComp      = stencil%QQ,               &
            &                     conf       = conf,                     &
            &                     bc_handle  = sub_handle,               &
            &                     varDict    = me(myBCID)%varDict,       &
            &                     varSys     = varSys                    )

        case( 'velocity_bounceback', 'velocity_eq', 'spc_inlet_eq',     &
          &   'spc_outlet_vel', 'spc_vel_bb', 'velocity_momentsbased',  &
          &   'spc_inlet', 'spc_moments_vel', 'spc_bb_vel_test',        &
          &   'velocity_bfl', 'vel_neq', 'velocity_noneq_expol',        &
          &   'spc_inflow', 'spc_velocity_noneq_expol',                 &
          &   'spc_solvent_inflow', 'inlet_nrbc'                        )
          write(logUnit(1),*) '   '//trim( bc_header%BC_kind(iBC) ) &
            &                 //' treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states &
            &                                            %velocity, &
            &                     state_name = 'velocity',          &
            &                     nComp      = 3,                   &
            &                     conf       = conf,                &
            &                     bc_handle  = sub_handle,          &
            &                     varDict    = me(myBCID)%varDict,  &
            &                     varSys     = varSys               )

          call aot_get_val( L       = conf,               &
            &               thandle = sub_handle,         &
            &               val     = me( myBCID )%order, &
            &               ErrCode = iError,             &
            &               key     = 'order',            &
            &               default = 1                   )

          me(myBCID)%nNeighs = me(myBCID)%order

          ! load molefraction in addition for spc_inlet_eq
          select case( trim( bc_header%BC_kind(iBC) ))
          case( 'spc_velocity_noneq_expol', 'spc_solvent_inflow' )
            if (.not. me( myBCID )%curved) me( myBCID )%qVal = 1.0_rk
            me( myBCID )%nNeighs = 1
          case( 'spc_inlet_eq', 'spc_inlet', 'spc_inflow' )
            call tem_load_bc_state( bc         = me(myBCID)%BC_states &
              &                                            %moleFrac, &
              &                     state_name = 'mole_fraction',     &
              &                     nComp      = 1,                   &
              &                     conf       = conf,                &
              &                     bc_handle  = sub_handle,          &
              &                     varDict    = me(myBCID)%varDict,  &
              &                     varSys     = varSys               )

            if ( trim( bc_header%BC_kind(iBC) ) == 'spc_inflow') then
              if (.not. me( myBCID )%curved) me( myBCID )%qVal = 1.0_rk
              me( myBCID )%requireNeighBufPre_nNext = .true.
              me( myBCID )%nNeighs = 1
            else
              if (.not. me( myBCID )%curved) me( myBCID )%qVal = 0.5_rk
              me( myBCID )%requireNeighBufPost = .true.
              me( myBCID )%nNeighs = 2
            end if

          case( 'moments_inflow' )
            call tem_load_bc_state( bc         = me(myBCID)%BC_states &
              &                                            %pressure, &
              &                     state_name = 'pressure',          &
              &                     nComp      = 1,                   &
              &                     conf       = conf,                &
              &                     bc_handle  = sub_handle,          &
              &                     varDict    = me(myBCID)%varDict,  &
              &                     varSys     = varSys               )

          !case ('inlet_bouzidi', 'velocity_bounceback_qval', 'velocity_bfl')
          case ('velocity_bounceback', 'velocity_bfl')
            write(logUnit(1),*) '    Inlet higher order with linear'   //  &
              &             ' interpolation treatment for velocity bc '//  &
              &            trim(me( myBCID )%label)
            if (.not. me( myBCID )%curved) me( myBCID )%qVal = 0.5_rk
            me( myBCID )%nNeighs = 0
            me( myBCID )%useComputeNeigh = .true.

          ! Load additionals for velocity_noneq_expol
          case( 'velocity_noneq_expol' )
            write(logUnit(1),*) '    Non-equilibrium extrapolation ' &
              &               //' treatment for '// trim(me( myBCID )%label)

            if (.not. me( myBCID )%curved) then
               me( myBCID )%qVal = 0.0_rk
               me( myBCID )%requireNeighBufPre_nNext = .true.
            else
               me( myBCID )%requireNeighBufPost = .true.
            end if
            me( myBCID )%nNeighs = 1
          case( 'inlet_nrbc' )
            write(logUnit(1),*) '    Characteristic based non-reflective ' &
              &                 //'inflow velocity boundary condition '// &
              &                 trim(me( myBCID )%label)
            me( myBCID )%requireNeighBufPre_nNext = .true.
            me( myBCID )%nNeighs = 2
          end select

        ! inlet mass flow rate boundary condition
        case( 'mfr_bounceback', 'mfr_eq' )
          write(logUnit(1),*) '   '//trim( bc_header%BC_kind(iBC) ) //  &
            &                 ' treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states     &
            &                                            %massFlowRate, &
            &                     state_name = 'mass_flowrate',         &
            &                     nComp      = 1,                       &
            &                     conf       = conf,                    &
            &                     bc_handle  = sub_handle,              &
            &                     varDict    = me(myBCID)%varDict,      &
            &                     varSys     = varSys                   )

          call aot_get_val( L       = conf,               &
            &               thandle = sub_handle,         &
            &               val     = me( myBCID )%order, &
            &               ErrCode = iError,             &
            &               key     = 'order',            &
            &               default = 1                   )

          me(myBCID)%nNeighs = me(myBCID)%order

        case( 'pressure_eq', 'pressure_antibounceback', 'outlet_bouzidi',      &
          &   'pressure_expol', 'pressure_expol_slow', 'pressure_noneq_expol', &
          &   'pressure_momentsbased', 'press_neq'                             )
          select case( trim( bc_header%BC_kind(iBC) ))
          case( 'pressure_eq' )
            write(logUnit(1),*) '    Outlet pressure equilibrium type for ' &
              &                 // trim(me( myBCID )%label)
            me( myBCID )%requireNeighBufPre_nNext = .true.
            me( myBCID )%nNeighs = 2

          case( 'pressure_antibounceback' )
            write(logUnit(1),*) '    Outlet pressure anti bounce back for ' &
              &                 // trim(me( myBCID )%label)
            me( myBCID )%requireNeighBufPost = .true.
            me( myBCID )%nNeighs = 1

          case( 'outlet_bouzidi' )
            write(logUnit(1),*) '    Outlet bouzidi boundary condition for' &
              &                 // ' ' // trim(me( myBCID )%label)
            me( myBCID )%nNeighs = 1

          case( 'press_neq' )
            write(logUnit(1),*) '    Pressure non-equilibrium boundary ' &
              &            //'condition for '//trim(me( myBCID )%label)
            me( myBCID )%nNeighs = 1

          case( 'pressure_noneq_expol' )
            write(logUnit(1),*) '    Pressure Non-equilibrium extrapolation' &
              &               //' treatment for '//trim(me( myBCID )%label)

            me( myBCID )%qVal = 0.0_rk
            me( myBCID )%requireNeighBufPre_nNext = .true.
            me( myBCID )%nNeighs = 1

          case( 'pressure_expol', 'pressure_expol_slow' )
            write(logUnit(1),*) '  Outlet extrapolation type for '//  &
              &                      trim(me( myBCID )%label)
            me( myBCID )%qVal = 0.0_rk
            me( myBCID )%requireNeighBufPre_nNext = .true.
            me( myBCID )%nNeighs = 2

          case( 'pressure_momentsbased' )
            write(logUnit(1),*) '    Moments based pressure boundary '// &
              &                 'type for '//                            &
              &                 trim(me( myBCID )%label)
            me( myBCID )%nNeighs = 0
          end select

          call tem_load_bc_state( bc         = me(myBCID)%BC_states &
            &                                            %pressure, &
            &                     state_name = 'pressure',          &
            &                     conf       = conf,                &
            &                     bc_handle  = sub_handle,          &
            &                     varDict    = me(myBCID)%varDict,  &
            &                     varSys     = varSys               )
        case( 'spc_outlet_eq', 'spc_outlet_expol' )
          select case( trim( bc_header%BC_kind(iBC) ))
          case( 'spc_outlet_eq' )
            write(logUnit(1),*) '    Species outlet pressure '// &
              &                 'equilibrium type for '//        &
              &                  trim(me( myBCID )%label)
            me( myBCID )%qVal = 0.5_rk
            me( myBCID )%nNeighs = 1

          case( 'spc_outlet_expol' )
            write(logUnit(1),*) '    Species outlet extrapolation type for ' &
              &                 // trim(me( myBCID )%label)
            me( myBCID )%qVal = 0.5_rk
            me( myBCID )%requireNeighBufPre_nNext = .true.
            me( myBCID )%nNeighs = 2
          end select

        case( 'outlet_nrbc', 'outlet_nrbc_eq' )
          me( myBCID )%requireNeighBufPre_nNext = .true.
          me( myBCID )%nNeighs = 2
          select case( trim( bc_header%BC_kind(iBC) ))
          case( 'outlet_nrbc' )
            write(logUnit(1),*) '    Characteristic based non-reflective ' &
              &                 //'outflow pressure boundary condition '// &
              &                 trim(me( myBCID )%label)
            call tem_load_bc_state( bc         = me(myBCID)%BC_states &
              &                                            %pressure, &
              &                     state_name = 'pressure',          &
              &                     conf       = conf,                &
              &                     bc_handle  = sub_handle,          &
              &                     varDict    = me(myBCID)%varDict,  &
              &                     varSys     = varSys               )

          case( 'outlet_nrbc_eq' )
            write(logUnit(1),*) '    Eq-type Characteristic based ' // &
              &             'non-reflective boundary condition '    // &
              &             trim(me( myBCID )%label)
            call tem_load_bc_state( bc         = me(myBCID)%BC_states &
              &                                            %pressure, &
              &                     state_name = 'pressure',          &
              &                     conf       = conf,                &
              &                     bc_handle  = sub_handle,          &
              &                     varDict    = me(myBCID)%varDict,  &
              &                     varSys     = varSys               )
          end select

          call aot_get_val( L       = conf,                    &
            &               thandle = sub_handle,              &
            &               val     = me( myBCID )%nrbc%kappa, &
            &               ErrCode = iError,                  &
            &               key     = 'kappa',                 &
            &               default = 1.0_rk                   )
          call aot_get_val( L       = conf,                    &
            &               thandle = sub_handle,              &
            &               val     = me( myBCID )%nrbc%sigma, &
            &               ErrCode = iError,                  &
            &               key     = 'sigma',                 &
            &               default = 0.58_rk                  )
          call aot_get_val(                                    &
            &         L       = conf,                          &
            &         thandle = sub_handle,                    &
            &         val     = me( myBCID )%nrbc%lodi_length, &
            &         ErrCode = iError,                        &
            &         key     = 'length'                       )

          if ( btest(iError, aoterr_Fatal) ) then
            write(logUnit(1),*)'FATAL Error occured in loading LODI length.'
            write(logUnit(1),*)'HINT:It is nElements between inlet and outlet'
            write(logUnit(1),*)'STOPPING'
            call tem_abort()
          end if

          if ( me(myBCID)%nrbc%lodi_length .feq. 0._rk ) then
            write(logUnit(1),*)'LODI length is zero.'
            write(logUnit(1),*)'HINT:It is nElements between inlet and outlet'
            write(logUnit(1),*)'STOPPING'
            call tem_abort()
          end if

          call aot_get_val(                             &
            &         L       = conf,                   &
            &         thandle = sub_handle,             &
            &         val     = me( myBCID )%nrbc%Ma_L, &
            &         ErrCode = iError,                 &
            &         key     = 'mach_lat'              )
          if ( btest(iError, aoterr_Fatal) ) then
            write(logUnit(1),*) 'FATAL Error occured in loading mach_lat.'
            write(logUnit(1),*) 'HINT: Ma_L = u_L/cs_L'
            write(logUnit(1),*) 'STOPPING'
            call tem_abort()
          end if

          write(logUnit(1),"(A,F9.2)") ' kappa: ', me( myBCID )%nrbc%kappa
          write(logUnit(1),"(A,F9.5)") ' sigma: ', me( myBCID )%nrbc%sigma
          write(logUnit(1),"(A,F9.2)") ' lodi L:', me( myBCID )%nrbc%lodi_length
          write(logUnit(1),"(A,F9.2)") ' mach_lat:', me( myBCID )%nrbc%Ma_L

        case( 'outlet_zero_prsgrd', 'spc_outlet_zero_prsgrd' )
          write(logUnit(1),*) '    Outlet zero pressure gradient for '  // &
            &                      trim(me( myBCID )%label)
          me( myBCID )%nNeighs = 1

        case( 'outlet_dnt' )
          write(logUnit(1),*) '    Outlet do-nothing type for '// &
            &                      trim(me( myBCID )%label)
          me( myBCID )%requireNeighBufPost = .true.
          me( myBCID )%nNeighs = 1

        case( 'moments_outflow' )
          write(logUnit(1),*) '    Outlet moments type for '// &
            &                      trim(me( myBCID )%label)
          me( myBCID )%requireNeighBufPre_nNext = .true.
          me( myBCID )%nNeighs = 2

        case( 'spc_molefrac', 'spc_molefrac_wtdf', 'spc_moments_molefrac', &
          &   'spc_molefrac_eq' )
          write(logUnit(1),*) '    Species '//trim(bc_header%BC_kind(iBC)) &
            &                 //'  treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states &
            &                                            %moleFrac, &
            &                     state_name = 'mole_fraction',     &
            &                     conf       = conf,                &
            &                     bc_handle  = sub_handle,          &
            &                     varDict    = me(myBCID)%varDict,  &
            &                     varSys     = varSys               )

          me( myBCID )%nNeighs = 0

        case( 'spc_mole_fraction_noneq_expol' )
          write(logUnit(1),*) '    Species '//trim(bc_header%BC_kind(iBC)) &
            &                 //'  treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states &
            &                                            %moleFrac, &
            &                     state_name = 'mole_fraction',     &
            &                     conf       = conf,                &
            &                     bc_handle  = sub_handle,          &
            &                     varDict    = me(myBCID)%varDict,  &
            &                     varSys     = varSys               )
          if (.not. me( myBCID )%curved) me( myBCID )%qVal = 0.0_rk
          me( myBCID )%nNeighs = 1

        case( 'spc_moledens_eq' )
          write(logUnit(1),*) '    Species '//trim(bc_header%BC_kind(iBC)) &
            &                 //'  treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states &
            &                                            %moleDens, &
            &                     state_name = 'mole_density',      &
            &                     conf       = conf,                &
            &                     bc_handle  = sub_handle,          &
            &                     varDict    = me(myBCID)%varDict,  &
            &                     varSys     = varSys               )

          me( myBCID )%nNeighs = 0


        case( 'spc_moleflux', 'spc_moleflux_eq', 'spc_moments_moleflux' )
          write(logUnit(1),*) '    Species molar flux treatment for ' // &
            &                 trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states &
            &                                            %moleFlux, &
            &                     state_name = 'mole_flux',         &
            &                     nComp      = 3,                   &
            &                     conf       = conf,                &
            &                     bc_handle  = sub_handle,          &
            &                     varDict    = me(myBCID)%varDict,  &
            &                     varSys     = varSys               )

          me( myBCID )%nNeighs = 0

        case( 'spc_blackbox_mem_ion', 'spc_blackbox_mem_solvent' )
          write(logUnit(1),*) '    Species black box membrance ' //    &
            &                 'treatment for '//trim(me( myBCID )%label)

          call aot_get_val( L       = conf,                                &
            &               thandle = sub_handle,                          &
            &               val     = me( myBCID )%blackbox_mem%transNr,   &
            &               ErrCode = iError,                              &
            &               key     = 'transference_number'                )
          if( btest( iError, aoterr_Fatal ))then
            write(logUnit(1),*)'FATAL Error occured in definition of '// &
              &            'black box membrance boundary condition'
            write(logUnit(1),*)'while retrieving transference number'
            write(logUnit(1),*)'STOPPING'
            call tem_abort()
          end if

          select case( trim( bc_header%BC_kind(iBC) ))
            case( 'spc_blackbox_mem_solvent' )
              call aot_get_val( L       = conf,                            &
                &               thandle = sub_handle,                      &
                &               val     = me( myBCID )                     &
                &                         %blackbox_mem%solvTrans_coeff,   &
                &               ErrCode = iError,                          &
                &               key     = 'transference_number'            )
              if( btest( iError, aoterr_Fatal ))then
                write(logUnit(1),*)'FATAL Error occured in definition of ' &
                  &            //'black box '//                            &
                  &            'membrance solvent boundary condition'
                write(logUnit(1),*)'while retrieving solv_trans_coeff'
                write(logUnit(1),*)'STOPPING'
                call tem_abort()
              end if
          end select
        case( 'spc_molediff_flux' )
          write(logUnit(1),*) '    Species molar diffusion flux '//    &
            &                 'treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states &
            &                                            %moleFrac, &
            &                     state_name = 'mole_fraction',     &
            &                     conf       = conf,                &
            &                     bc_handle  = sub_handle,          &
            &                     varDict    = me(myBCID)%varDict,  &
            &                     varSys     = varSys               )

          call tem_load_bc_state( bc         = me(myBCID)%BC_states      &
            &                                            %moleDiff_flux, &
            &                     state_name = 'mole_diff_flux',         &
            &                     nComp      = 3,                        &
            &                     conf       = conf,                     &
            &                     bc_handle  = sub_handle,               &
            &                     varDict    = me(myBCID)%varDict,       &
            &                     varSys     = varSys                    )

          me( myBCID )%nNeighs = 0

        case( 'potential_noneq_expol' )
          write(logUnit(1),*) '    Non-equilibrium extrapolation '//   &
            &                 'treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states  &
            &                                            %potential, &
            &                     state_name = 'potential',          &
            &                     conf       = conf,                 &
            &                     bc_handle  = sub_handle,           &
            &                     varDict    = me(myBCID)%varDict,   &
            &                     varSys     = varSys                )

          if (.not. me( myBCID )%curved) me( myBCID )%qVal = 0.0_rk
          me( myBCID )%requireNeighBufPost = .true.
          me( myBCID )%nNeighs = 1

        case( 'potential_neumann' )
          write(logUnit(1),*) '    Non-equilibrium extrapolation  '//  &
            &                 'treatment for '//trim(me( myBCID )%label)
          call tem_load_bc_state( bc         = me(myBCID)%BC_states      &
            &                                            %surChargeDens, &
            &                     state_name = 'surface_charge_density', &
            &                     conf       = conf,                     &
            &                     bc_handle  = sub_handle,               &
            &                     varDict    = me(myBCID)%varDict,       &
            &                     varSys     = varSys                    )

          !KM: For Neumann BC, qVal should be 0.5 to be consistent with second
          ! order extrapolation of potential field
          if (.not. me( myBCID )%curved) me( myBCID )%qVal = 0.5_rk
          me( myBCID )%requireNeighBufPost = .true.
          me( myBCID )%nNeighs = 1

        case( 'moledens_noneq_expol' )
          write(logUnit(1),*) '    Non-equilibrium extrapolation '//  &
            &                 'treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states &
            &                                            %velocity, &
            &                     state_name = 'velocity',          &
            &                     nComp      = 3,                   &
            &                     conf       = conf,                &
            &                     bc_handle  = sub_handle,          &
            &                     varDict    = me(myBCID)%varDict,  &
            &                     varSys     = varSys               )

          call tem_load_bc_state( bc         = me(myBCID)%BC_states  &
            &                                            %moleDens,  &
            &                     state_name = 'mole_density',       &
            &                     conf       = conf,                 &
            &                     bc_handle  = sub_handle,           &
            &                     varDict    = me(myBCID)%varDict,   &
            &                     varSys     = varSys                )

          if (.not. me( myBCID )%curved) me( myBCID )%qVal = 1.0_rk
          me( myBCID )%nNeighs = 1

        case( 'moledens_neumann' )
          write(logUnit(1),*) '    Non-equilibrium extrapolation' &
            &               //' treatment for '//trim(me( myBCID )%label)

          call tem_load_bc_state( bc         = me(myBCID)%BC_states &
            &                                            %velocity, &
            &                     state_name = 'velocity',          &
            &                     nComp      = 3,                   &
            &                     conf       = conf,                &
            &                     bc_handle  = sub_handle,          &
            &                     varDict    = me(myBCID)%varDict,  &
            &                     varSys     = varSys               )

          call tem_load_bc_state( bc         = me(myBCID)%BC_states  &
            &                                            %moleFlux,  &
            &                     state_name = 'mole_flux',          &
            &                     conf       = conf,                 &
            &                     bc_handle  = sub_handle,           &
            &                     varDict    = me(myBCID)%varDict,   &
            &                     varSys     = varSys                )

          if (.not. me( myBCID )%curved) me( myBCID )%qVal = 0.5_rk
          me( myBCID )%requireNeighBufPost = .true.
          me( myBCID )%nNeighs = 1

        case( 'flekkoy_outlet' )
          write(logUnit(1),*) '   Outlet BCs for the flekkoy ' // &
            &                 'diffusion model'
          me( myBCID )%nNeighs = 2

        case( 'flekkoy_inlet' )
          write(logUnit(1),*) '   Inlet BCs for the flekkoy ' // &
            &                 'diffusion model'
          me( myBCID )%nNeighs = 0

        case( 'test' )
          write(logUnit(1),*) '   WARNING!! '
          write(logUnit(1),*) '   Test boundary conditions for testing ' // &
            &                 'neighbor info exchange'
          me( myBCID )%nNeighs = 20

        case default
          write(logUnit(1),*) 'Unknown Boundary condition kind!'
          write(logUnit(1),*) 'Stopping ...'
          call tem_abort()
        end select

        write(logUnit(3),*)'   Default qVal: ', me( myBCID )%qVal

        if ( me(myBCID)%curved ) then
           me( myBCID )%requireNeighBufPost = .true.
           !KM! \todo activate this line after implementing nonEq expol boundary
           ! condition for curved and straight boundary seperately
           ! me( myBCID )%requireNeighBufPre = .false.

           ! Set nNeighs only if its set to zero before else use already set
           ! nNeighs
           if ( me( myBCID )%nNeighs == 0 ) then
             me( myBCID )%nNeighs = 1
           end if
        end if

        write(logUnit(1),"(A,I0)") '   required neighbors: ', &
          &                        me( myBCID )%nNeighs
        ! Do special treatment for corner elements which are intersected
        ! by multiple boundaries only for moments based BC
        select case( trim( bc_header%BC_kind(iBC) ))
        case( 'pressure_momentsbased', 'moments_wall', 'velocity_momentsbased',&
          &   'moments_inflow', 'moments_outflow',                             &
          &   'spc_moments_molefrac', 'spc_moments_moleflux',                  &
          &   'spc_moments_wall', 'spc_moments_vel', 'spc_moments_outflow'     )
          me( myBCID )%treat_corner = .true.
        case default
          me( myBCID )%treat_corner = .false.
        end select

        call truncate( me = me(myBCID)%varDict )
      end if  ! sub_handle/=0
      call aot_table_close( L=conf, thandle=sub_handle )
    end do ! iBC
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_load_bc
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_init_bc_elems( me, nElems, QQN, hasQVal )
    ! -------------------------------------------------------------------- !
    type( bc_elems_type ), intent(inout) :: me
    integer, intent(in) :: nElems
    !> Number of direction excluding center since center is never treated for
    !! bc elems
    integer, intent(in) :: QQN
    logical, intent(in) :: hasQval
    ! -------------------------------------------------------------------- !

    if ( allocated( me%elem%val ) ) then
      call destroy( me%normalInd  )
      call destroy( me%posInBcElemBuf )
      call destroy(me%elem)
      call destroy(me%normal)
      call destroy(me%bitMask)
    end if

    call init( me = me%normalInd, length=nElems )
    call init( me = me%posInBcElemBuf, length=nElems)
    call init( me = me%elem, length=nElems )
    call init( me = me%bitMask, width = QQN, length=nElems )
    call init( me = me%normal, width = 3, length=nElems )

    if ( hasQVal ) then
      if ( allocated( me%qVal%val ) ) then
        call destroy( me%qVal )
      end if

      call init( me = me%qVal, width=QQN )
    end if

  end subroutine mus_init_bc_elems
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> It removes non-valid elements while still maintaining the origianl order.
  !! The given bc elements (elems) contains both valid and non-valid elements.
  !! Position of valid elements are given by posInBCElem.
  !! Valid elements are moved towards the start of the elems so that they become
  !! continuous in the elems.
  !!
  subroutine rearrange_bc_elems( minLevel, maxLevel, nValid, posInBCElem, &
    &                            nElems, elemLvl )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: minLevel, maxLevel
    !> number of valid BC elements
    integer, intent(in) :: nValid(minLevel:maxLevel)
    !> BC elements information
    integer, intent(in) :: nElems(minLevel:maxLevel)
    !> their position in bc_elems_type
    integer, intent(in) :: posInBCElem(maxval(nElems),minLevel:maxLevel)
    !> BC elements information
    type( bc_elems_type ), intent(inout) :: elemLvl(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iLevel, oldPos
    ! -------------------------------------------------------------------- !

    write(dbgUnit(1), *) 'Inside remove solid in bcElems routine ...'

    do iLevel = minLevel, maxLevel

      ! Re arrange bc elements only when any solid are found
      if ( nValid(iLevel) < nElems(iLevel) ) then
        do iElem = 1, nValid( iLevel )
          oldPos = posInBCElem(iElem,iLevel)
          elemLvl(iLevel)%elem%val( iElem ) = elemLvl(iLevel)%elem%val( oldPos )
          elemLvl(iLevel)%bitMask%val( :,iElem ) =                             &
            &                            elemLvl(iLevel)%bitmask%val( :,oldPos )
          elemLvl(iLevel)%normal%val( :,iElem )  =                             &
            &                            elemLvl(iLevel)%normal%val( :,oldPos )
          elemLvl(iLevel)%normalInd%val( iElem ) =                             &
            &                            elemLvl(iLevel)%normalInd%val( oldPos )
          write(dbgUnit(6), *) 'elem: ', elemLvl( iLevel )%elem%val( iElem )
        end do

        ! update growing array nVals
        elemLvl(iLevel)%elem%nVals = nValid(iLevel)
        elemLvl(iLevel)%bitMask%nVals = nValid(iLevel)
        elemLvl(iLevel)%normal%nVals = nValid(iLevel)
        elemLvl(iLevel)%normalInd%nVals = nValid(iLevel)

      end if ! nValid(iLevel) == nElems(iLevel)

    end do ! iLevel

    write(dbgUnit(1), *) 'Leave  remove solid in bcElems routine ...'
    write(dbgUnit(1), *) ''

  end subroutine rearrange_bc_elems
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> It count valid (non-solid) elements in BC elements list.
  !! Input:
  !!   minLevel, maxLevel
  !!   LevelPointer
  !!   LevelDesc
  !!   nElems - number of BC elements
  !!   elems - positions of BC elements in tree or levelPointer
  !! Output:
  !!   nValid - number of valid BC elements
  !!   posInBCElem - positions of valid elements in BC elements list
  !!
  subroutine check_solid_in_bc( minLevel, maxLevel, levelPointer, levelDesc, &
    &                           nElems, elemLvl, nValid, posInBCElem )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: minLevel, maxLevel
    !> Level Pointer ( position of a element in level desc )
    integer, intent(in) :: levelPointer(:)
    !> Level Descriptor
    type( tem_levelDesc_type ), intent(in) :: levelDesc(minLevel:maxLevel)
    !> number of BC elements
    integer, intent(in) :: nElems(minLevel:maxLevel)
    !> BC elements list that contains their position in levelPointer
    type( bc_elems_type ), intent(in) :: elemLvl(minLevel:maxLevel)
    !> number of valid (non-solid) elements
    integer, intent(out) :: nValid(minLevel:maxLevel)
    !> positions of valid elements in globBC elements list
    integer, intent(out) :: posInBCElem(maxval(nElems),minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iLevel, posInTotal
    ! -------------------------------------------------------------------- !

    write(dbgUnit(1), *) 'Inside check solid in bc routine ...'
    write(dbgUnit(1), *) 'If any solid are found, they are listed as below.'

    posInBCElem(:,:) = 0
    nValid(:) = 0

    do iLevel = minLevel, maxLevel

      do iElem = 1, nElems( iLevel )
        posInTotal = levelPointer( elemLvl( iLevel )%elem%val( iElem ) )
        if ( btest( levelDesc(iLevel)%property(posInTotal), prp_solid ) ) then
          write(dbgUnit(3), "(A,I2,2(A,I0))")    &
            &  'found a Solid! level: ', iLevel, &
            &  ', posInTotal: ', posInTotal,     &
            &  ', treeID: ', levelDesc(iLevel)%total(posInTotal)
        else
          nValid(iLevel) = nValid(iLevel) + 1
          posInBCElem( nValid(iLevel), iLevel ) = iElem
        end if
        write(dbgUnit(6), *) 'elem: ', elemLvl( iLevel)%elem%val( iElem )

      end do ! iElem

      if ( nElems(iLevel) /= nValid( iLevel) ) then
        write(logUnit(5),"(2(A,I0))") &
          &  ' Removed ', nElems(iLevel) - nValid(iLevel), &
          &  ' elements on level ', iLevel
      end if
      write(dbgUnit(3), *) 'nElems: ', nElems( iLevel )
      write(dbgUnit(3), *) 'nValid: ', nValid( iLevel )

    end do ! iLevel

    write(dbgUnit(1), *) 'Leave  check solid in bc routine ...'
    write(dbgUnit(1), *) ''

  end subroutine check_solid_in_bc
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine debug_glob_boundary_type( nBCs, minLevel, maxLevel, me )
    ! -------------------------------------------------------------------- !
    integer,                    intent(in) :: nBCs, minLevel, maxLevel
    type( glob_boundary_type ), intent(in) :: me( nBCs )
    ! -------------------------------------------------------------------- !
    integer :: iBC, iLevel, iElem
    ! -------------------------------------------------------------------- !

    write(dbgUnit(1),"(A)")'-- DEBUG glob_boundary_type  ----------------------'
    do iBC = 1, nBCs
      write(dbgUnit(5), "(I0,A)") iBC, ', BC label: '//trim(me(iBC)%label)
      write(dbgUnit(5), "(A,L5)")  'hasQVal: ', me(iBC)%hasQVal
      write(dbgUnit(5), "(A,L5)")  'qValInitialized: ', me(iBC)%qValInitialized
      write(dbgUnit(5), "(A,L5)")  'isWall: ', me(iBC)%isWall
      write(dbgUnit(5), "(A,I0)") 'nElems_local: ', me(iBC)%nElems_local
      write(dbgUnit(5), "(A,I0)") 'nElems_total: ', me(iBC)%nElems_total
      if ( .not. me( iBC )%isWall ) then
        do iLevel = minLevel, maxLevel
          write(dbgUnit(5), "(A,I0)") 'iLevel: ', iLevel
          write(dbgUnit(5), "(A,I0)") 'nElems: ', me(iBC)%nElems(iLevel)
          write(dbgUnit(5), "(A,I0)") 'nElems_totalLevel: ', &
            &                         me(iBC)%nElems_totalLevel(iLevel)
          write(dbgUnit(5), "(4A10)") 'elem', 'elemBuf', 'normalInd', 'normal'
          do iElem = 1, me(iBC)%nElems(iLevel)
            write(dbgUnit(5),"(6I10)")                              &
              &  me(iBC)%elemLvl(iLevel)%elem%val(iElem),           &
              &  me(iBC)%elemLvl(iLevel)%posInBcElemBuf%val(iElem), &
              &  me(iBC)%elemLvl(iLevel)%normalInd%val(iElem),      &
              &  me(iBC)%elemLvl(iLevel)%normal%val(1:3,iElem)
          end do ! iElem
          write(dbgUnit(1), "(A)") ''
        end do ! iLevel
      else ! isWall == true
        write(dbgUnit(5),"(A)") 'This is wall. No need to check bc elems'
      end if ! isWall
      write(dbgUnit(1), "(A)") ''
    end do ! iBC

  end subroutine debug_glob_boundary_type
  ! ------------------------------------------------------------------------ !

  ! ------------------------------------------------------------------------ !
  !> Set BC elements positions in LevelDesc%neigh%nghElems
  !!
  subroutine mus_set_posInNghElems( minLevel, maxLevel, nStencils, &
      &                             globBC, BC )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: minLevel, maxLevel, nStencils
    type( glob_boundary_type ), intent(in) :: globBC
    type( boundary_type ),   intent(inout) :: BC
    ! -------------------------------------------------------------------- !
    integer :: iLevel, iElem
    integer :: counter(minLevel:maxLevel,nStencils), stencilPos
    ! -------------------------------------------------------------------- !

    if ( BC%nNeighs > 0 ) then
      counter = 0
      do iLevel = minLevel, maxLevel
        do iElem = 1, globBC%nElems(iLevel)

          stencilPos = BC%elemLvl(iLevel)%stencilPos(iElem)
          counter(iLevel, stencilPos) = counter(iLevel, stencilPos) + 1
          BC%elemLvl(iLevel)%posInNghElems(iElem) = counter(iLevel, stencilPos)

        end do ! iElem
      end do ! iLevel
    end if

  end subroutine mus_set_posInNghElems
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_set_bcLinks( iField, QQ, QQN, nScalars, nElems, elemLvl, &
    &                         nSize, neigh, links )
    ! -------------------------------------------------------------------- !
    integer,                intent(in) :: iField, QQ, QQN, nScalars
    !> number of BC elements
    integer,                intent(in) :: nElems, nSize
    type( bc_elems_type ),  intent(in) :: elemLvl
    integer,                intent(in) :: neigh(:)
    type( array_type ),    intent(out) :: links
    ! -------------------------------------------------------------------- !
    integer :: iElem, iDir, iLink, posInTotal
    integer, allocatable :: dirs(:)
    ! -------------------------------------------------------------------- !

    allocate( dirs( nElems * QQ ) )
    iLink = 0

    write(dbgUnit(7), "(A,I0)") 'nSize: ', nSize

    do iElem = 1, nElems
      posInTotal = elemLvl%elem%val( iElem )
      write(dbgUnit(7), "(2(A,I0))") 'iElem: ', iElem, &
        &                            ', posInTotal: ', posInTotal
      do iDir = 1, QQN
        if( elemLvl%bitmask%val( iDir, iElem )) then
          iLink = iLink + 1
          dirs(iLink) =  neigh((idir-1)* nsize+ posintotal)+( ifield-1)* qq+ nscalars*0
        end if
      end do ! iDir
    end do ! iElem

    allocate( links%val( iLink ) )
    links%val(1:iLink) = dirs(1:iLink)
    links%nVals        = iLink

    deallocate( dirs )

  end subroutine mus_set_bcLinks
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Set the coefficients of bouzidi linear interpolation boundary condition.
  !!
  subroutine set_bouzidi_coeff( qVal, cIn, cOut, cNgh, cVel )
    ! -------------------------------------------------------------------- !
    real(kind=rk),  intent(in) :: qVal
    real(kind=rk), intent(out) :: cIn, cOut, cNgh, cVel
    ! -------------------------------------------------------------------- !
    if ( qVal >= 0.5_rk ) then
      cIn  = 1._rk - 0.5_rk / qVal
      cOut =         0.5_rk / qVal
      cNgh = 0.0_rk
      cVel = 0.5_rk / qVal
    else
      cIn  = 0.0_rk
      cOut =         2.0_rk * qVal
      cNgh = 1._rk - 2.0_rk * qVal
      cVel = 1._rk
    end if

  end subroutine set_bouzidi_coeff
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Set necessary data for Wall Bouzidi BC
  !! bouzidi should be allocated beforehand
  subroutine mus_set_bouzidi( iLevel, QQ, QQN, nScalars, globBC, cxDirInv, &
    &                         varPos, bouzidi )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: iLevel, QQ, QQN, nScalars
    type( glob_boundary_type ), intent(in) :: globBC
    integer, intent(in)  :: cxDirInv(QQ)
    integer, intent(in)  :: varPos(:)
    !> cIn, cOut, cNgh should have size of nLinks
    type( bc_wall_bouzidi_type ), intent(inout) :: bouzidi
    ! -------------------------------------------------------------------- !
    integer :: iLink, iDir, invDir, iElem, posInBCBuf
    real(kind=rk) :: cVelTmp ! not required for wall_libb
    ! -------------------------------------------------------------------- !

    iLink = 0

    do iElem = 1, globBC%nElems(iLevel)

      posInBCBuf = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )

      do iDir = 1, QQN
        if( globBC%elemLvl(iLevel)%bitmask%val( iDir, iElem )) then
          invDir = cxDirInv(iDir)
          iLink = iLink + 1
          bouzidi%qVal( iLink ) = globBC%elemLvl( iLevel )       &
            &                           %qVal%val( invDir, iElem )
          call set_bouzidi_coeff(            &
            &    qVal = bouzidi%qVal(iLink), &
            &    cIn  = bouzidi%cIn(iLink),  &
            &    cOut = bouzidi%cOut(iLink), &
            &    cNgh = bouzidi%cNgh(iLink), &
            &    cVel = cVelTmp              )

          ! inPos and outPos are used in bcBuffer
          bouzidi%inPos( iLink ) = varPos(iDir) + ( posInBCBuf-1 ) * nScalars
          bouzidi%outPos( iLink ) = varPos(invDir) + ( posInBCBuf-1 ) * nScalars
          ! nghPos is used in computeNeighBuf
          bouzidi%nghPos( iLink ) = invDir + (iElem-1)*QQ
        end if
      end do ! iDir

    end do ! iElem

  end subroutine mus_set_bouzidi
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_alloc_bouzidi( me, nVals, minLevel, maxLevel )
    ! -------------------------------------------------------------------- !
    type( bc_wall_bouzidi_type ), allocatable :: me(:)
    integer, intent(in) :: minLevel, maxLevel
    integer, intent(in) :: nVals(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    ! -------------------------------------------------------------------- !

    allocate( me( minLevel:maxLevel ) )
    write(logUnit(7), "(2(A,I0))") 'Allocated wall libb type level: ', &
      &                             minLevel, ' - ', maxLevel

    do iLevel = minLevel, maxLevel
      if ( nVals(iLevel) > 0 ) then
        write(logUnit(9), "(2(A,I0))") 'level: ', iLevel, ', nVal: ', &
          &                            nVals(iLevel)
      end if
      if ( nVals(iLevel) >= 0 ) then
        allocate( me(iLevel)%  qVal( nVals(iLevel) ) )
        allocate( me(iLevel)%   cIn( nVals(iLevel) ) )
        allocate( me(iLevel)%  cOut( nVals(iLevel) ) )
        allocate( me(iLevel)%  cNgh( nVals(iLevel) ) )
        allocate( me(iLevel)% inPos( nVals(iLevel) ) )
        allocate( me(iLevel)%outPos( nVals(iLevel) ) )
        allocate( me(iLevel)%nghPos( nVals(iLevel) ) )
        ! Intialize with zero
        me(iLevel)%  qVal=0.0_rk
        me(iLevel)%   cIn=0.0_rk
        me(iLevel)%  cOut=0.0_rk
        me(iLevel)%  cNgh=0.0_rk
        me(iLevel)% inPos=0
        me(iLevel)%outPos=0
        me(iLevel)%nghPos=0
      end if
    end do ! iLevel

  end subroutine mus_alloc_bouzidi
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_alloc_fieldBC( me, minLevel, maxLevel )
    ! -------------------------------------------------------------------- !
    type( boundary_type ) :: me
    integer, intent(in) :: minLevel, maxLevel
    ! -------------------------------------------------------------------- !

    if ( .not. allocated( me%elemLvl )) allocate( me%elemLvl(minLevel:maxLevel))
    if ( .not. allocated( me%neigh ) ) allocate( me%neigh(minLevel:maxLevel) )
    if ( .not. allocated( me%links ) ) allocate( me%links(minLevel:maxLevel) )

  end subroutine mus_alloc_fieldBC
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Set necessary data for BC velocity_bounceback_qval
  subroutine mus_set_inletUbb( me, tree, stencil, nScalars,       &
    &                          globBC, levelDesc, varPos, nLinks, &
    &                          minLevel, maxLevel                 )
    ! -------------------------------------------------------------------- !
    !> setting type for bc
    type(bc_inlet_bouzidi_type), allocatable, intent(out) :: me(:)
    !> using mesh information
    type(treelmesh_type), intent(in)         :: tree
    !> for directions
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> number of scalars
    integer, intent(in)                      :: nScalars
    !> for number of elements in boundary and position in buffer
    type(glob_boundary_type), intent(in)  :: globBC
    !> minimum and maximum level
    integer, intent(in)                   :: minLevel, maxLevel
    !> level descriptor
    type(tem_levelDesc_type), intent(in)  :: levelDesc(minLevel:maxLevel)
    !> for position of outgoing elements
    integer, intent(in)                   :: varPos(:)
    !> for linkwise treatment
    integer, intent(in)                   :: nLinks(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iLink, iDir, invDir, iElem, posInBCBuf
    integer :: iLevel
    real(kind=rk) :: dx, qVal
    real(kind=rk), dimension(3) :: bary
    ! -------------------------------------------------------------------- !

    write(logUnit(10), "(A)") 'Setting inletUbb type ...'

    allocate( me(minLevel:maxLevel) )

    ! loop over levels
    do iLevel = minLevel, maxLevel
      if (nLinks(iLevel) > 0) then
        write(logUnit(10), "(2(A,I0))") 'level: ', iLevel, ', nVal: ', &
          &                             nLinks(iLevel)
      end if
      if (nLinks(iLevel) >= 0) then
        ! allocate arrays
        allocate( me(iLevel)%  qVal(nLinks(iLevel)) )
        allocate( me(iLevel)%outPos(nLinks(iLevel)) )
        allocate( me(iLevel)%  iDir(nLinks(iLevel)) )
        allocate( me(iLevel)%  posInBuffer(nLinks(iLevel)) )
        ! Intialize arrays with zero
        do iLink = 1, nLinks(iLevel)
          me(iLevel)%  qVal(iLink) = 0.0_rk
          me(iLevel)%outPos(iLink) = 0
          me(iLevel)%  iDir(iLink) = 0
          me(iLevel)%  posInBuffer(iLink) = -99
        end do
      end if

      ! fill arrays
      iLink = 0
      dx    = tem_ElemSizeLevel(tree, iLevel)

      ! loop over elements
      do iElem = 1, globBC%nElems(iLevel)
        ! position in buffer and position of barycenter
        posInBCBuf = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)
        bary(:)    = tem_BaryOfId(tree, levelDesc(iLevel)          &
          &                        %total(globBC%elemLvl(iLevel)   &
          &                        %elem%val(iElem) ) )

        ! loop over directions
        do iDir = 1, stencil%QQN
          !> Bitmask is true for incoming direction
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then

            invDir  = stencil%cxDirInv(iDir)
            iLink   = iLink + 1
            ! posInBcElemBuf is needed for bcBuffer
            me(iLevel)%posInBuffer(iLink) = posInBCBuf
            ! qValues
            qVal    = globBC%elemLvl(iLevel)%qVal%val(invDir, iElem)
            !me(iLevel)%   qVal(iLink) = qVal needed here?
            ! position for outgoing elememts
            me(iLevel)% outPos(iLink) = varPos(invDir)            &
              &                                   + (posInBCBuf-1) * nScalars
            ! store directions to access stencil weight and cxDirRK
            me(iLevel)%   iDir(iLink) = iDir

          end if
        end do ! iDir

      end do ! iElem

    end do ! iLevel

  end subroutine mus_set_inletUbb
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
 subroutine mus_set_inletBfl( me, tree, stencil, nScalars, globBC,           &
    &                         levelDesc, varPos, nLinks, minLevel, maxLevel  )
    ! -------------------------------------------------------------------- !
    !> setting type for bc
    type(bc_inlet_bouzidi_type), allocatable, intent(out) :: me(:)
    !> using mesh information
    type(treelmesh_type), intent(in)         :: tree
    !> for directions
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> number of scalars
    integer, intent(in)                      :: nScalars
    !> for number of elements in boundary and position in buffer
    type(glob_boundary_type), intent(in)     :: globBC
    !> minimum and maximum level
    integer, intent(in)                      :: minLevel, maxLevel
    !> level descriptor
    type(tem_levelDesc_type), intent(in)     :: levelDesc(minLevel:maxLevel)
    !> for position of elements
    integer, intent(in)                      :: varPos(:)
    !> for linkwise treatment
    integer, intent(in)          :: nLinks(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iLink, iDir, invDir, iElem, posInBCBuf
    integer :: iLevel
    real(kind=rk) :: dx, qVal
    real(kind=rk), dimension(3) :: bary
    ! -------------------------------------------------------------------- !

    write(logUnit(10), "(A)") 'Setting inletBFL ...'

    allocate( me(minLevel:maxLevel) )

    ! loop over levels
    do iLevel = minLevel, maxLevel
      if (nLinks(iLevel) > 0) then
        write(logUnit(10),"(2(A,I0))") 'level: ', iLevel, &
          &                            ', nVal: ', nLinks(iLevel)
      end if
      if (nLinks(iLevel) >= 0) then
        ! allocate arrays
        ! allocate( me(iLevel)%  qVal(nLinks(iLevel)) )
        allocate( me(iLevel)%   cIn(nLinks(iLevel)) )
        allocate( me(iLevel)%  cOut(nLinks(iLevel)) )
        allocate( me(iLevel)%  cNgh(nLinks(iLevel)) )
        allocate( me(iLevel)%  cVel(nLinks(iLevel)) )
        allocate( me(iLevel)% inPos(nLinks(iLevel)) )
        allocate( me(iLevel)%outPos(nLinks(iLevel)) )
        allocate( me(iLevel)%nghPos(nLinks(iLevel)) )
        allocate( me(iLevel)%  iDir(nLinks(iLevel)) )
        allocate( me(iLevel)%  posInBuffer(nLinks(iLevel)) )
        ! Intialize arrays with zero
        ! me(iLevel)%  qVal = 0.0_rk
        do iLink = 1, nLinks(iLevel)
          me(iLevel)%   cIn(iLink) = 0.0_rk
          me(iLevel)%  cOut(iLink) = 0.0_rk
          me(iLevel)%  cNgh(iLink) = 0.0_rk
          me(iLevel)% inPos(iLink) = 0
          me(iLevel)%outPos(iLink) = 0
          me(iLevel)%nghPos(iLink) = 0
          me(iLevel)%  iDir(iLink) = 0
          me(iLevel)%  posInBuffer(iLink) = -99
        end do
      end if

      ! fill arrays
      iLink = 0
      dx    = tem_ElemSizeLevel( tree, iLevel )

      ! loop over elements
      do iElem = 1, globBC%nElems(iLevel)

        ! position in buffer and position of barycenter
        posInBCBuf = globBC%elemLvl( iLevel )%posInBcElemBuf%val( iElem )
        bary(:)    = tem_BaryOfId( tree, levelDesc( iLevel )        &
          &                        %total(globBC%elemLvl ( iLevel ) &
          &                        %elem%val( iElem ) ) )

        ! loop over directions
        do iDir = 1, stencil%QQN
          !> Bitmask is true for incoming direction
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            invDir  = stencil%cxDirInv(iDir)
            iLink   = iLink + 1
            ! posInBcElemBuf is needed for bcBuffer
            me(iLevel)%posInBuffer(iLink) = posInBCBuf
            ! qValues
            qVal = globBC%elemLvl( iLevel )%qVal%val( invDir, iElem )
            ! coefficients for state values
            call set_bouzidi_coeff(               &
              &    qVal = qVal,                   &
              &    cIn  = me(iLevel)%cIn(iLink),  &
              &    cOut = me(iLevel)%cOut(iLink), &
              &    cNgh = me(iLevel)%cNgh(iLink), &
              &    cVel = me(iLevel)%cVel(iLink)  )
            ! positions of incoming, outgoing and neighbour elemets
            ! inPos and outPos are used in bcBuffer
            me(iLevel)%  inPos(iLink) = varPos(iDir) &
              &                        + (posInBCBuf-1) * nScalars
            me(iLevel)% outPos(iLink) = varPos(invDir) &
              &                        + (posInBCBuf-1) * nScalars
            ! nghPos is used in computeNeighBuf
            me(iLevel)% nghPos(iLink) = invDir + (iElem-1)*stencil%QQ
            ! store directions to access stencil weight and cxDirRK
            me(iLevel)%   iDir(iLink) = iDir

          end if
        end do ! iDir

      end do ! iElem

    end do ! iLevel

  end subroutine mus_set_inletBfl
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Linkwise non-equilibrium extrapolation (can handle curved walls)
  !!
  !! All the coefficients are pre-calculated here and then used for the specific
  !! boundary conditions
  !! - potential_noneq_expol & potential_neumann
  !! - velocity_noneq_expol & pressure_noneq_expol
  !! - spc_moledens_noneq_expol & spc_solvent_noneq_expol
  !!
  !! The pdf is decomposed into equilibrium (eq) and non-equilibrium (neq) part:
  !! f = f_eq + f_neq
  !! - f_eq is calculated by weighting a fictitious X, which is obtained
  !!   by an extrapolation using the fluid neighbor(s)
  !! - f_neq is approximated by second-order extrapolation using the fluid
  !!   neighbor(s)
  !! - for qVal < 0.75 even the second neighbor is used for the extrapolations
  !!
  subroutine mus_set_nonEqExpol( me, curved, tree, stencil, nScalars, globBC, &
    &                            bc_neigh, pdf, levelDesc, varPos, nLinks,    &
    &                            minLevel, maxLevel  )
    ! -------------------------------------------------------------------- !
    !> minimum and maximum level
    integer, intent(in)                      :: minLevel, maxLevel
    !> setting type for bc
    type(bc_nonEqExpol_type), allocatable, intent(out) :: me(:)
    !> Curved or straight boundary
    logical, intent(in) :: curved
    !> using mesh information
    type(treelmesh_type), intent(in)         :: tree
    !> for directions
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> number of scalars
    integer, intent(in)                      :: nScalars
    !> scheme global boundary type
    type(glob_boundary_type), intent(in)     :: globBC
    !> boundary neighbor
    type(bc_neigh_type), intent(in)          :: bc_neigh(minLevel:maxLevel)
    !> contains global state vector
    type( pdf_data_type ), intent(in)        :: pdf( tree%global%minlevel   &
      &                                              : tree%global%maxlevel )
    !> level descriptor
    type(tem_levelDesc_type), intent(in)     :: levelDesc(minLevel:maxLevel)
    !> for position of elements
    integer, intent(in)                      :: varPos(:)
    !> for linkwise treatment
    integer, intent(in)                      :: nLinks(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iLink, invDir, iElem
    integer :: iLevel, iDir, elemPos, QQ, QQN
    real(kind=rk) :: qVal
    ! -------------------------------------------------------------------- !

    write(logUnit(10), "(A)") 'Setting nonEqExpol ...'

    QQ  = stencil%QQ
    QQN = stencil%QQN

    allocate( me(minLevel:maxLevel) )

    ! loop over levels
    do iLevel = minLevel, maxLevel
      if (nLinks(iLevel) > 0) then
        write(logUnit(10),"(2(A,I0))") 'level: ', iLevel, ', nVal: ', &
          &                            nLinks(iLevel)
      end if
      if (nLinks(iLevel) >= 0) then
        ! allocate arrays
        if (curved) then
          allocate( me(iLevel)%          c_w(nLinks(iLevel)) )
          allocate( me(iLevel)%          c_f(nLinks(iLevel)) )
          allocate( me(iLevel)%         c_ff(nLinks(iLevel)) )
          allocate( me(iLevel)%      c_neq_f(nLinks(iLevel)) )
          allocate( me(iLevel)%     c_neq_ff(nLinks(iLevel)) )
        end if
        allocate( me(iLevel)%  posInBuffer(nLinks(iLevel)) )
        allocate( me(iLevel)%posInNeighBuf(nLinks(iLevel)) )
        allocate( me(iLevel)%         iDir(nLinks(iLevel)) )
        allocate( me(iLevel)%  posInBCelems(nLinks(iLevel)) )

        ! Intialize arrays
        if (curved) then
          me(iLevel)%          c_w(:) = -99.0_rk
          me(iLevel)%          c_f(:) = -99.0_rk
          me(iLevel)%         c_ff(:) = -99.0_rk
          me(iLevel)%      c_neq_f(:) = -99.0_rk
          me(iLevel)%     c_neq_ff(:) = -99.0_rk
        end if
        me(iLevel)%posInNeighBuf(:) = -99
        me(iLevel)%         iDir(:) = -99
        me(iLevel)%  posInBuffer(:) = -99
        me(iLevel)%  posInBCelems(:) = -99
      end if

      ! fill arrays
      iLink = 1

      ! loop over elements
      do iElem = 1, globBC%nElems(iLevel)
        ! loop over directions
        do iDir = 1, QQN
          !> Bitmask is true for incoming direction
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            ! posInBcElemBuf is needed for bcBuffer
            me(iLevel)%posInBuffer(iLink) = globBC%elemLvl( iLevel ) &
              &                                   %posInBcElemBuf%val( iElem )

            !save iDir for current link
            me(iLevel)%iDir(iLink) = iDir

            ! elemPos is needed for overnext neighbor later on
            elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
            ! determine neighbor along incoming direction
            invDir  = stencil%cxDirInv(iDir)

            ! position of the boundary element of this link in globBC elems list
            me(iLevel)%posInBCelems(iLink) = iElem

            ! Position of fluid Neigbor in neighBufferPre/Post
            ! Use neighbor along the normal direction of boundary
            ! if neighbor along direction of link is not found or
            ! boundary is defined as straight boundary
            me(iLevel)%posInNeighBuf(iLink) = iElem

            ! For curved boundary find the neighbor along the incoming link
            if (curved) then
              ! qValues
              qVal = globBC%elemLvl( iLevel )%qVal%val( invDir, iElem )

              ! coefficients for state values depending on qVal
              ! compute coeffient for physical boundary (_w)
              if ( (qVal - 0.75_rk) >= 0.0_rk ) then
                me(iLevel)%c_w(iLink) = 1.0_rk / qVal
                ! compute coeffient for fluid c_f
                me(iLevel)%c_f(iLink) = ( qVal - 1.0_rk ) / qVal
                ! compute coeffient for overnext fluid c_ff
                me(iLevel)%c_ff(iLink) = 0.0_rk
                ! compute nonEq coefficient for fluid c_nEq_f
                me(iLevel)%c_nEq_f(iLink) = 1.0_rk
                ! compute nonEq coefficient for overnext fluid c_nEq_ff
                me(iLevel)%c_nEq_ff(iLink) = 0.0_rk
              else
                me(iLevel)%c_w(iLink) = 1.0_rk + ( 2.0_rk * (1.0_rk - qVal) ) &
                  &                           / ( 1.0_rk + qVal )
                ! compute coeffient for fluid c_f
                me(iLevel)%c_f(iLink) = qVal - 1.0_rk
                ! compute coeffient for overnext fluid c_ff
                me(iLevel)%c_ff(iLink) = ( (1.0_rk - qVal)    &
                  &                       * (qVal - 1.0_rk) ) &
                  &                      / (1.0_rk + qVal)

                me(iLevel)%c_nEq_f(iLink) = qVal
                me(iLevel)%c_nEq_ff(iLink) = 1.0_rk - qVal
              end if ! qVal > 0.75
            end if ! curved

            ! increase link counter
            iLink   = iLink + 1

          end if ! bitmask
        end do ! iDir
      end do ! iElem
    end do ! iLevel

  end subroutine mus_set_nonEqExpol
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine mus_set_outletExpol( me, stencil, globBC,      &
    &                             nLinks, minLevel, maxLevel )
    ! -------------------------------------------------------------------- !
    !> setting type for bc
    type(bc_outlet_type), allocatable, intent(out) :: me(:)
    !> for directions
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> number of scalars
    ! integer, intent(in)                      :: nScalars
    !> for number of elements in boundary and position in buffer
    type(glob_boundary_type), intent(in)     :: globBC
    !> minimum and maximum level
    integer, intent(in) :: minLevel, maxLevel
    !> for position of elements
    ! integer, intent(in)                      :: varPos(:)
    !> Number of links on each level
    integer, intent(in) :: nLinks(minLevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iLevel, iLink, iElem, iDir
    ! -------------------------------------------------------------------- !

    write(logUnit(10), "(A)") 'Setting outlet expol ...'

    allocate( me(minLevel:maxLevel) )

    ! loop over levels
    do iLevel = minLevel, maxLevel
      if ( nLinks(iLevel) > 0) then
        write(logUnit(10), '(2(A,I0))') 'level: ', iLevel, &
          &                             ', nLinks: ', nLinks(iLevel)
      end if
      if ( nLinks(iLevel) >= 0) then
        ! allocate arrays
        allocate( me(iLevel)%statePos(nLinks(iLevel)) )
        allocate( me(iLevel)%iElem(nLinks(iLevel)) )
        ! Intialize arrays with zero
        me(iLevel)%statePos = 0
      end if

      iLink = 0

      ! loop over elements
      do iElem = 1, globBC%nElems(iLevel)
        ! position in buffer is needed afterwards
        ! posInBCBuf = globBC%elemLvl( iLevel )%posInBcElemBuf%val(iElem)

        ! loop over directions
        do iDir = 1, stencil%QQN
          ! check bitmask
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            ! invDir  = stencil%cxDirInv(iDir)
            iLink   = iLink + 1
            ! store position
            me(iLevel)%statePos(iLink) = iDir + (iElem-1)*stencil%QQ
            me(iLevel)%iElem(iLink) = iElem
          end if ! check bitmask

        end do ! iDir

      end do ! iElem

    end do ! iLevel

  end subroutine mus_set_outletExpol
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routines deallocates allocatables in field%bc boundary_type for
  !! dynamic load balancing
  subroutine mus_fieldBC_cleanup( me, nBCs, minLevel, maxLevel )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: nBCs, minLevel, maxLevel
    type( boundary_type ), intent(inout) :: me(nBCs)
    ! -------------------------------------------------------------------- !
    integer :: iBC, iLevel
    ! -------------------------------------------------------------------- !

    write(dbgUnit(1),*) "Enter mus_fieldBC_cleanup"

    do iBC = 1, nBCs
      if( me(iBC)%nNeighs > 0 ) then
        do iLevel = minLevel, maxLevel
          deallocate( me(iBC)%neigh( iLevel )%posInState )
          if( me(iBC)%requireNeighBufPre_nNext  ) &
            & deallocate( me(iBC)%neigh( iLevel )%neighBufferPre_nNext )
          if( me(iBC)%requireNeighBufPre  ) &
            & deallocate( me(iBC)%neigh( iLevel )%neighBufferPre )
          if( me(iBC)%requireNeighBufPost ) &
            & deallocate( me(iBC)%neigh( iLevel )%neighBufferPost )
          deallocate( me(iBC)%elemLvl( iLevel )%stencilPos )
          deallocate( me(iBC)%elemLvl( iLevel )%posInNghElems )
        end do
      end if

      if( me(iBC)%useComputeNeigh ) then
        do iLevel = minLevel, maxLevel
          deallocate( me(iBC)%neigh( iLevel )%computeNeighBuf )
        end do
      end if

      do iLevel = minLevel, maxLevel

        ! @todo: for pure wall, following array is not allocated
        !        we should have better way to detect this.
        if ( allocated( me(iBC)%links( iLevel )%val ) ) then
          deallocate( me(iBC)%links( iLevel )%val )
        end if

        if ( allocated( me(iBC)%elemLvl( iLevel )%lodi ) ) then
          deallocate( me(iBC)%elemLvl( iLevel )%lodi )
        end if

        if ( allocated( me(iBC)%elemLvl( iLevel )%moments ) ) then
          deallocate( me(iBC)%elemLvl( iLevel )%moments )
        end if

        if ( allocated( me(iBC)%bouzidi ) ) then
          deallocate( me(iBC)%bouzidi( iLevel )%qVal )
          deallocate( me(iBC)%bouzidi( iLevel )%cIn  )
          deallocate( me(iBC)%bouzidi( iLevel )%cOut )
          deallocate( me(iBC)%bouzidi( iLevel )%cNgh )
          deallocate( me(iBC)%bouzidi( iLevel )%inPos )
          deallocate( me(iBC)%bouzidi( iLevel )%outPos )
          deallocate( me(iBC)%bouzidi( iLevel )%nghPos )
        end if

        if ( allocated( me(iBC)%inletUbbQVal ) ) then
          deallocate( me(iBC)%inletUbbQVal( iLevel )%qVal )
          deallocate( me(iBC)%inletUbbQVal( iLevel )%outPos )
          deallocate( me(iBC)%inletUbbQVal( iLevel )%iDir )
          deallocate( me(iBC)%inletUbbQVal( iLevel )%posInBuffer )
        end if

        if ( allocated( me(iBC)%inletBFL ) ) then
          deallocate( me(iBC)%inletBFL( iLevel )%cIn )
          deallocate( me(iBC)%inletBFL( iLevel )%cOut )
          deallocate( me(iBC)%inletBFL( iLevel )%cNgh )
          deallocate( me(iBC)%inletBFL( iLevel )%inPos )
          deallocate( me(iBC)%inletBFL( iLevel )%outPos )
          deallocate( me(iBC)%inletBFL( iLevel )%nghPos )
          deallocate( me(iBC)%inletBFL( iLevel )%cVel )
          deallocate( me(iBC)%inletBFL( iLevel )%iDir )
          deallocate( me(iBC)%inletBFL( iLevel )%posInBuffer )
        end if

        if ( allocated( me(iBC)%nonEqExpol ) ) then
          if ( me(iBC)%curved ) then
            deallocate( me(iBC)%nonEqExpol( iLevel )%c_w )
            deallocate( me(iBC)%nonEqExpol( iLevel )%c_f)
            deallocate( me(iBC)%nonEqExpol( iLevel )%c_ff )
            deallocate( me(iBC)%nonEqExpol( iLevel )%c_neq_f )
            deallocate( me(iBC)%nonEqExpol( iLevel )%c_neq_ff )
          end if
          deallocate( me(iBC)%nonEqExpol( iLevel )%posInNeighBuf )
          deallocate( me(iBC)%nonEqExpol( iLevel )%iDir )
          deallocate( me(iBC)%nonEqExpol( iLevel )%posInBuffer )
          deallocate( me(iBC)%nonEqExpol( iLevel )%posInBCelems )
        end if

        if ( allocated( me(iBC)%outletExpol ) ) then
          deallocate( me(iBC)%outletExpol( iLevel )%statePos )
          deallocate( me(iBC)%outletExpol( iLevel )%iElem )
        end if

        if ( allocated( me(iBC)%turbWallFunc%dataOnLvl ) ) then
          deallocate( me(iBC)%turbWallFunc%dataOnLvl( iLevel )%tVisc )
          deallocate( me(iBC)%turbWallFunc%dataOnLvl( iLevel )%velTau )
          deallocate( me(iBC)%turbWallFunc%dataOnLvl( iLevel )%distToBnd )
          deallocate( me(iBC)%turbWallFunc%dataOnLvl( iLevel )%neighDistToBnd )
          deallocate( me(iBC)%turbWallFunc%dataOnLvl( iLevel )%unitNormal )
        end if

      end do ! iLevel = minLevel, maxLevel

      if ( allocated( me(iBC)%links ) ) then
        deallocate( me(iBC)%neigh )
        deallocate( me(iBC)%elemLvl )
        deallocate( me(iBC)%links )
      end if

      if ( allocated( me(iBC)%bouzidi ) ) then
        deallocate( me(iBC)%bouzidi )
      end if
      if ( allocated( me(iBC)%inletUBBqVal ) ) then
        deallocate( me(iBC)%inletUBBqVal )
      end if
      if ( allocated( me(iBC)%inletBFL ) ) then
        deallocate( me(iBC)%inletBFL )
      end if
      if ( allocated( me(iBC)%nonEqExpol ) ) then
        deallocate( me(iBC)%nonEqExpol )
      end if
      if ( allocated( me(iBC)%outletExpol ) ) then
        deallocate( me(iBC)%outletExpol )
      end if
      if ( allocated( me(iBC)%turbWallFunc%dataOnLvl ) ) then
        deallocate( me(iBC)%turbWallFunc%dataOnLvl )
      end if
    end do

    write(dbgUnit(1),*) "Done!"

  end subroutine mus_fieldBC_cleanup
  ! ------------------------------------------------------------------------ !

end module mus_bc_header_module
! **************************************************************************** !
