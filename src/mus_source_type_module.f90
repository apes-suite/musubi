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
! ***************************************************************************** !
!> author: Kannan Masilamani
!! Module containing subroutines for building MUSUBI specific source
!! variables
!!
module mus_source_type_module
  ! include treelm modules
  use env_module,               only: rk, long_k, labelLen
  use tem_tools_module,         only: upper_to_lower
  use tem_aux_module,           only: tem_abort
  use tem_varMap_module,        only: tem_possible_variable_type, &
    &                                 init, append, truncate,     &
    &                                 tem_variable_loadMapping
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_stringKeyValuePair_module, only: grw_stringKeyValuePairArray_type, &
    &                                      init, append, truncate
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_logging_module,       only: logUnit
  use tem_shape_module,         only: tem_shape_type, tem_load_shape
  use tem_subTree_type_module,  only: tem_subTree_type

  ! include musubi modules
  use mus_scheme_header_module, only: mus_scheme_header_type
  use mus_physics_module,       only: mus_convertFac_type
  use mus_derVarPos_module,     only: mus_derVarPos_type
  use mus_absorbLayer_module,   only: mus_absorbLayer_type,        &
    &                                 mus_absorbLayer_dynAvg_type, &
    &                                 mus_load_absorbLayer

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_Fatal, &
    &                         aoterr_NonExistent, aoterr_WrongType
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val

  implicit none
  private

  public :: mus_source_type
  public :: mus_source_op_type
  public :: mus_turbChannelForce_type
  public :: mus_HRRCorrectionTerm_type

  public :: mus_create_poss_srcVar
  public :: mus_load_source_var
  public :: mus_source_cleanup

  public :: mus_applySrc_dummy
  public :: mus_addSrcToAuxField_dummy
  public :: mus_updateSrcVar_dummy

  ! ************************************************************************** !
  !> Stores correction term for HRR_bgk
  type mus_HRRCorrectionTerm_type
    !> density
    real(kind=rk), allocatable :: dens(:)
    !> velocity
    real(kind=rk), allocatable :: vel(:,:)
  end type mus_HRRCorrectionTerm_type
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Contains information to compute average bulk velocity for dynamic_force.
  !! In turb_channel_force, the force term is adapted according to difference
  !! between reference bulk velocity and simulated plane average bulk velocity
  !! to avoid linear increase in simulated bulk velocits.
  !! For more information:
  !! https://www.wias-berlin.de/people/john/ELECTRONIC_PAPERS/JR07.IJNMF.pdf
  type mus_turbChannelForce_type
    !> tracking shapes
    type(tem_shape_type)  :: geom_utau(1) !1
    type(tem_shape_type)  :: geom_umean(1) !1

    !> sub-tree resulting from the elements within the tracking shape
    !! The sub-tree also holds the sub-communicator
    !! This data needs to be UPDATED after balance
    type(tem_subTree_type) :: subTree_utau !one for utau and umean
    type(tem_subTree_type) :: subTree_umean !one for utau and umean

    !> Reference bulk velocity in physical unit
    real(kind=rk) :: refVelBulk

    !> Characteristic height in physical unit
    real(kind=rk) :: refHeight

    !> Stream-wise direction to compute average velocity
    !! x=1, y=2, z=3
    integer :: flow_direction

    !> Dynamic Force term for turbulent channel in physical unit [m/s^2]
    !! F_dyn = (refVelBulk-avgVelXBulk) * refVelBulk / refHeight
    real(kind=rk) :: forceDyn(3)

    !> Global number of elements in defined shape
    integer :: nElemsGlobal_utau
    integer :: nElemsGlobal_umean
  end type mus_turbChannelForce_type
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Contains source elements position in state array and idx to access
  !! data variable refered in config file.
  !! This type is defined for each level
  type mus_source_elems_type
    !> Number of source elements on this level.
    !! nFluids + nGhosts
    integer :: nElems

    !> Position of elements in state array to apply source terms.
    !! Position in state array is same as position in total list
    !! Size: nElems
    integer, allocatable :: posInTotal(:)

    !> Index to access point data type to retrieve values from variable
    !! refered for source variable
    integer, allocatable :: idx(:)

    !> source field value obtained from ST_fun data variable.
    !! Filled only for elements where source is active i.e. elements in
    !! posInTotal.
    !! size: nElems*nComponents
    !! \todo KM: might be not neccessary
!KM!    real(kind=rk), allocatable :: val(:)

    !> Contains time average values of density and velocity for dynamic
    !! absorblayer.
    !! \todo KM: 02042021 Introduce method_data c_ptr and point to
    !! dynAvg for absorbLayer and change intent(inout) to intent(in) in
    !! proc_addSrcToAuxField.
    type(mus_absorbLayer_dynAvg_type) :: dynAvg

    ! source term for HRR_bgk
    type(mus_HRRCorrectionTerm_type) :: HRR_Corr
  end type mus_source_elems_type
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Description contains list of elements on which source is active and
  !! function pointer to update source
  type mus_source_op_type
    !> Position of this source term variable in the varSys
    integer :: srcTerm_varPos

    !> Position of data variable provided in config file in the varSys
    integer :: data_varPos

    !> Contains source elements position in state array/total list for
    !! each level
    type(mus_source_elems_type), allocatable :: elemLvl(:)

    !> Function to update state with source term
    procedure(proc_apply_source), pointer :: applySrc => null()

    !> to use source field array
    !KM!logical :: useSrcFieldVal = .false.

    !> name of the source variable
    character(len=labelLen) :: varname

    !> Function pointer to append source field to auxilary variable
    procedure(proc_addSrcToAuxField), pointer :: addSrcToAuxField => null()

    !> Function pointer to update source variable which are dependent on
    !! auxField.
    procedure(proc_updateSourceVar), pointer :: updateSourceVar => null()

    !> Order of approximation for source like force, electric_field,
    !! charge_density.
    !! Order = 1, uses force term in BE approximated by forward Euler method
    !! Order = 2, uses force term in BE approximated by Trapezoidal method.
    !! For order 2, macroscopic source is also added to auxField.
    !! For fluid, fluid_incompressible, multispecies_liquid: source is added
    !! to momentum and for poisson: source is added to potential.
    !! Default: order = 2.
    integer :: order

    !> Additional config information for absorbLayer
    type(mus_absorbLayer_type) :: absLayer

    !> Contains information to compute average bulk velocity for dynamic_force
    type(mus_turbChannelForce_type) :: turbChanForce
  end type mus_source_op_type
  ! *************************************************************************** !


  ! *************************************************************************** !
  !> Description of musubi source type
  type mus_source_type
    !> Contains source elements position in tree%treeID and
    !! function pointer to update source
    !! Size: varDict%nVals
    type(mus_source_op_type), allocatable :: method(:)

    !> Dictionary of source variable with
    !! varDict%val()%key is the name of source variable and
    !! varDict%val()%value is the name of variable provided for the key
    type(grw_stringKeyValuePairArray_type) :: varDict
  end type mus_source_type
  ! *************************************************************************** !


  ! *************************************************************************** !
  abstract interface
    !> Abstract interface to update state with source terms
    subroutine proc_apply_source( fun, inState, outState, neigh, auxField, &
      & nPdfSize, iLevel, varSys, time, phyConvFac, derVarPos )
      import :: rk, mus_source_op_type, tem_varSys_type, tem_time_type, &
        &       mus_convertFac_type, mus_derVarPos_type

      !> Description of method to update source
      class(mus_source_op_type), intent(in) :: fun

      !> input  pdf vector
      !! \todo KM: instate is passed to compute auxField.
      !! Since auxField is precomputed from instate and passed to this routine.
      !! instate can be removed
      real(kind=rk), intent(in)  ::  inState(:)

      !> output pdf vector
      real(kind=rk), intent(inout) :: outState(:)

      !> connectivity Array corresponding to state vector
      integer,intent(in) :: neigh(:)

      !> auxField array
      real(kind=rk), intent(in) :: auxField(:)

      !> number of elements in state Array
      integer, intent(in) :: nPdfSize

      !> current level
      integer, intent(in) :: iLevel

      !> variable system
      type(tem_varSys_type), intent(in) :: varSys

      !> Point in time at which to evaluate the variable.
      type(tem_time_type), intent(in)  :: time

      !> Physics conversion factor for current level
      type(mus_convertFac_type), intent(in) :: phyConvFac

      !> position of derived quantities in varsys
      type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    end subroutine proc_apply_source

    !> Interface to add source to auxField vars in source_op_type for
    !! all nSolve elements (nFluids+nGhostFromCoarser+nGhostFromFiner).
    !! Halo elements are exchanged
    subroutine proc_addSrcToAuxField(fun, auxField, iLevel, time, varSys, &
      & phyConvFac, derVarPos)
      import :: rk, tem_varSys_type, mus_derVarPos_type, mus_convertFac_type, &
        &       tem_time_type, mus_source_op_type

      !> Description of method to update source
      class(mus_source_op_type), intent(inout) :: fun
      !> output auxField array
      real(kind=rk), intent(inout)          :: auxField(:)
      !> current level
      integer, intent(in)                   :: iLevel
      !> current timing information
      type(tem_time_type), intent(in)       :: time
      !> variable system definition
      type(tem_varSys_type), intent(in)     :: varSys
      !> Physics conversion factor for current level
      type(mus_convertFac_type), intent(in) :: phyConvFac
      !> position of derived quantities in varsys
      type(mus_derVarPos_type), intent(in)  :: derVarPos(:)
    end subroutine proc_addSrcToAuxField

    !> Interface to update source variable which has dependency on auxField.
    !! Applied on all nSolve elements (nFluids+nGhostFromCoarser+nGhostFromFiner).
    !! Halo elements are exchanged
    !! This should be called after adding sorce term to state so that both
    !! auxField and apply_source uses same source value in one multilevel cycle.
    subroutine proc_updateSourceVar(fun, auxField, iLevel, varSys, phyConvFac, &
      & derVarPos)
      import :: rk, tem_varSys_type, mus_derVarPos_type, mus_convertFac_type, &
        &       treelmesh_type, mus_source_op_type

      !> Description of method to update source
      class(mus_source_op_type), intent(inout) :: fun
      !> input auxField array on current level
      real(kind=rk), intent(in)          :: auxField(:)
      !> current level
      integer, intent(in) :: iLevel
      !> variable system definition
      type(tem_varSys_type), intent(in)     :: varSys
      !> Physics conversion factor on current level
      type(mus_convertFac_type), intent(in) :: phyConvFac
      !> position of derived quantities in varsys
      type(mus_derVarPos_type), intent(in)  :: derVarPos(:)
    end subroutine proc_updateSourceVar


  end interface
  ! *************************************************************************** !

contains

  ! *************************************************************************** !
  !> Routine initialize possible source variable depends on scheme kind
  subroutine mus_create_poss_srcVar(poss_srcVar, schemeHeader)
    ! --------------------------------------------------------------------------!
    !> possible source variables
    type(tem_possible_variable_type), intent(out) :: poss_srcVar

    !> Identifier of the scheme
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------!
    integer :: QQ
    ! --------------------------------------------------------------------------!
    ! used only for HRR Correction term
    select case (trim(schemeHeader%layout))
    case ('d2q9')
      QQ = 9
    case ('d3q19')
      QQ = 19
    case ('d3q27')
      QQ = 27
    case default
      ! Musubi will be aborted later, other layouts not supported!!!!
      QQ = -1
    end select

    write(logUnit(10),*) 'Creating possible source terms '
    call init(me = poss_srcVar, length = 2 )

    select case(trim(schemeHeader%kind))
    case ('fluid', 'fluid_incompressible')
      ! body force
      call append(me          = poss_srcVar, &
        &         varName     = 'force',     &
        &         nComponents = 3            )
      ! dynamic body force for turbulent channel
      call append(me          = poss_srcVar,                &
        &         varName     = 'turb_channel_force_accel', &
        &         nComponents = 3                           )
      ! absorb layer, STfun should return sponge_strength
      ! Target pressure and velocity are defined under absorb_layer_target
      call append(me          = poss_srcVar,    &
        &         varName     = 'absorb_layer', &
        &         nComponents = 1               )

      ! Absorn layer for inlet boundary with dynamic pressure and
      ! constant velocity.
      call append(me          = poss_srcVar,          &
        &         varName     = 'absorb_layer_inlet', &
        &         nComponents = 1                     )
      ! Absorn layer for outlet boundary with constant pressure and
      ! dynamic velocity.
      call append(me          = poss_srcVar,           &
        &         varName     = 'absorb_layer_outlet', &
        &         nComponents = 1                     )
      ! HRR correction term, might be usefull also for other schemes??
      call append(me          = poss_srcVar,      &
        &         varName     = 'hrr_correction', &
        &         nComponents = QQ                )
      ! In the Brinkman term \( -\nu / K \mathbf{u} \),
      ! this is its coefficient \nu / K
      call append(me          = poss_srcVar, &
        &         varName     = 'brinkman',  &
        &         nComponents = 1            )

    case ('fluid_GNS', 'fluid_incompressible_GNS')
      ! body force
      call append(me          = poss_srcVar, &
        &         varName     = 'force',     &
        &         nComponents = 3            )

    case ('passive_scalar')
      call append(me          = poss_srcVar,       &
        &         varName     = 'equal_injection', &
        &         nComponents = 1                  )

      call append(me          = poss_srcVar, &
        &         varName     = 'injection', &
        &         nComponents = 1            )

      ! source term \( \alpha C \) for passive scalar C
      ! where \( \alpha \) is ps_sourceCoeff
      call append(me          = poss_srcVar,      &
        &         varName     = 'ps_sourceCoeff', &
        &         nComponents = 1                 )

    case ('nernst_planck')
      call append(me          = poss_srcVar,      &
        &         varName     = 'electric_field', &
        &         nComponents = 3                 )

    case ('multispecies_liquid')
      call append(me          = poss_srcVar,      &
        &         varName     = 'electric_field', &
        &         nComponents = 3                 )

      call append(me          = poss_srcVar, &
        &         varName     = 'force',     &
        &         nComponents = 3            )

    case ('poisson')
      call append(me          = poss_srcVar,      &
        &         varName     = 'charge_density', &
        &         nComponents = 1                 )

    case default
      write(logUnit(1),*) 'No possible source term defined for scheme kind:' &
        &                 //trim(schemeHeader%kind)
    end select
    call truncate(poss_srcVar)

  end subroutine mus_create_poss_srcVar
  ! *************************************************************************** !



  ! ***************************************************************************!
  !> Routine load musubi source terms for given key.
  !! key is glob_source or source
  subroutine mus_load_source_var(me, possVars, conf, parent, key, varSys)
    ! --------------------------------------------------------------------------!
    !> Source variable type to initialize
    type(mus_source_type), intent(out) :: me
    !> possible source variables
    type(tem_possible_variable_type), intent(in) :: possVars
    !> flu state
    type( flu_State ) :: conf
    !> parent handle if scheme table is defined
    integer, intent(in), optional :: parent
    !> key to load source term
    character(len=*), intent(in) :: key
    !> Global variable system
    type(tem_varSys_type), intent(inout) :: varSys
    ! --------------------------------------------------------------------------!
    integer :: iSrc, iError, srchandle
    character(len=labelLen) :: varname_order
    ! --------------------------------------------------------------------------!
    ! initialize growing array stringKeyValuePair
    call init( me = me%varDict )

    ! load the global source variables
    call tem_variable_loadMapping( possVars = possVars,   &
      &                            conf     = conf,       &
      &                            parent   = parent,     &
      &                            key      = key,        &
      &                            varDict  = me%varDict, &
      &                            varSys   = varSys      )

    ! truncate varDict
    call truncate( me = me%varDict )

    ! If source is defined then reopen source table to load order of
    ! approximatation to use for defined source variable
    if (me%varDict%nVals > 0) then
      ! allocate source method
      allocate(me%method(me%varDict%nVals))

      call aot_table_open( L       = conf,      &
        &                  parent  = parent,    &
        &                  thandle = srchandle, &
        &                  key     = key        )

      do iSrc = 1, me%varDict%nVals
         varname_order = trim(me%varDict%val(iSrc)%key)//'_order'
         call aot_get_val( L       = conf,                  &
           &               thandle = srchandle,             &
           &               key     = varname_order,         &
           &               val     = me%method(iSrc)%order, &
           &               default = 2,                     &
           &               ErrCode = iError                 )

         if (btest(iError, aoterr_Fatal)) then
           write(logUnit(1),*)'FATAL Error occured, while retrieving source' &
             &                //'order'
           if (btest(iError, aoterr_WrongType)) then
             write(logUnit(1),*)'Variable has wrong type!'
             write(logUnit(1),*)'STOPPING'
             call tem_abort()
           endif
         end if
         write(logUnit(1),'(A,I0)') ' Order of approximation for '          &
           &                        //trim(me%varDict%val(iSrc)%key)//': ', &
           &                        me%method(iSrc)%order

         ! Load additional information for absorblayer
         select case (trim(me%varDict%val(iSrc)%key))
         case ('turb_channel_force_accel')
           call load_turbChanForce( me     = me%method(iSrc)%turbChanForce, &
             &                      conf   = conf,                          &
             &                      key    = 'turb_channel_force_dynamic',  &
             &                      parent = srchandle                      )
         case ('absorb_layer')
           call mus_load_absorbLayer( me       = me%method(iSrc)%absLayer &
             &                                     %config,               &
             &                        conf     = conf,                    &
             &                        key      = 'absorb_layer_target',   &
             &                        parent   = srchandle,               &
             &                        loadPres = .true.,                  &
             &                        loadVel  = .true.                   )
         case ('absorb_layer_inlet')
           call mus_load_absorbLayer( me       = me%method(iSrc)%absLayer     &
             &                                     %config,                   &
             &                        conf     = conf,                        &
             &                        key      = 'absorb_layer_inlet_target', &
             &                        parent   = srchandle,                   &
             &                        loadPres = .false.,                     &
             &                        loadVel  = .true.                       )
         case ('absorb_layer_outlet')
           call mus_load_absorbLayer( me       = me%method(iSrc)%absLayer      &
             &                                     %config,                    &
             &                        conf     = conf,                         &
             &                        key      = 'absorb_layer_outlet_target', &
             &                        parent   = srchandle,                    &
             &                        loadPres = .true.,                       &
             &                        loadVel  = .false.                       )
         end select
      end do
    else
      allocate(me%method(0))
    end if

  end subroutine mus_load_source_var
  ! ***************************************************************************!

  ! ************************************************************************** !
  !> Load shape, bulk velocity and height for turbulent channel force
  subroutine load_turbChanForce(me, conf, key, parent)
    ! -------------------------------------------------------------------------!
    !> Turbulent channel force
    type(mus_turbChannelForce_type), intent(out) :: me
    !> flu state
    type( flu_State ) :: conf
    !> Table name to load target states
    character(len=*), intent(in) :: key
    !> parent source handle
    integer, intent(in) :: parent
    ! -------------------------------------------------------------------------!
    integer :: turbForce_handle, iError
    character(len=1) :: flow_direction
    ! -------------------------------------------------------------------------!
    ! -------------------------------------------------------------------------!

    call aot_table_open( L       = conf,             &
      &                  parent  = parent,           &
      &                  thandle = turbForce_handle, &
      &                  key     = trim(key)         )

    ! Load Shape, bulk velocity and height for turbulent channel force
    write(logUnit(1),*) ' * Turbulent channel force:'
    call aot_get_val( L       = conf,                &
      &               thandle = turbForce_handle,    &
      &               key     = 'ref_velocity_bulk', &
      &               val     = me%refVelBulk,       &
      &               ErrCode = iError               )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1), *) 'Error loading ref_velocity_bulk!'
      call tem_abort()
    end if
    write(logUnit(1),*) '    ref_velocity_bulk =', me%refVelBulk

    call aot_get_val( L       = conf,             &
      &               thandle = turbForce_handle, &
      &               key     = 'ref_height',     &
      &               val     = me%refHeight,     &
      &               ErrCode = iError            )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1), *) 'Error loading ref_height!'
      call tem_abort()
    end if

    call aot_get_val( L       = conf,              &
      &               thandle = turbForce_handle,  &
      &               key     = 'flow_direction', &
      &               val     = flow_direction,   &
      &               ErrCode = iError            )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1), *) 'Error loading flow_direction!'
      call tem_abort()
    end if

    write(logUnit(1),*) '    flow_direction ='//flow_direction
    select case(flow_direction)
    case ('x', 'X')
      me%flow_direction = 1
    case ('y', 'Y')
      me%flow_direction = 2
    case ('z', 'Z')
      me%flow_direction = 3
    case default
      call tem_abort('Unknown flow_direction')
    end select

    ! load geometry
    call tem_load_shape( conf    = conf,             &
      &                  parent  = turbForce_handle, &
      &                  key     = 'shape_utau',     &
      &                  me      = me%geom_utau(1)   )

    call tem_load_shape( conf    = conf,             &
      &                  parent  = turbForce_handle, &
      &                  key     = 'shape_umean',    &
      &                  me      = me%geom_umean(1)  )

    if (.not. allocated(me%geom_utau(1)%canoND)) then
      allocate(me%geom_utau(1)%canoND(0))
    end if
    if (.not. allocated(me%geom_utau(1)%bcLabels)) then
      allocate(me%geom_utau(1)%bcLabels(0))
    end if

    if (size(me%geom_utau(1)%canoND) == 0        &
      & .and. size(me%geom_utau(1)%bcLabels) == 0) then
      write(logUnit(1),*) 'Error: Requires canond or boundary shape for u_tau'
      call tem_abort()
    end if

    if (.not. allocated(me%geom_umean(1)%canoND)) then
      allocate(me%geom_umean(1)%canoND(0))
    end if

    if (size(me%geom_umean(1)%canoND) == 0         &
      & .and. trim(me%geom_umean(1)%kind) /= 'all' ) then
      write(logUnit(1),*) 'Error: Requires canond or all shape for u_mean'
      call tem_abort()
    end if


  end subroutine load_turbChanForce
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine act as a destructor for source type.
  !! The arrays allocated in mus_init_sourceTerms are destroyed here
  subroutine mus_source_cleanup(me)
    ! --------------------------------------------------------------------------
    type(mus_source_type), intent(inout)   :: me
    ! --------------------------------------------------------------------------
    integer :: iSrc
    ! --------------------------------------------------------------------------
    ! KM: DO NOT DESTROY VARDICT IN SOURCE_TYPE AS IT CONTAINS CONFIG INFO
    do iSrc = 1, me%varDict%nVals
      deallocate(me%method(iSrc)%elemLvl)
    end do
  end subroutine mus_source_cleanup
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Dummy routine for apply source
  subroutine mus_applySrc_dummy( fun, inState, outState, neigh, auxField, &
    & nPdfSize, iLevel, varSys, time, phyConvFac, derVarPos )
    ! --------------------------------------------------------------------------
    !> Description of method to update source
    class(mus_source_op_type), intent(in) :: fun
    !> input  pdf vector
    !! \todo KM: instate is passed to compute auxField.
    !! Since auxField is precomputed from instate and passed to this routine.
    !! instate can be removed
    real(kind=rk), intent(in)  ::  inState(:)
    !> output pdf vector
    real(kind=rk), intent(inout) :: outState(:)
    !> connectivity Array corresponding to state vector
    integer,intent(in) :: neigh(:)
    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)
    !> number of elements in state Array
    integer, intent(in) :: nPdfSize
    !> current level
    integer, intent(in) :: iLevel
    !> variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! --------------------------------------------------------------------------
    ! abort only if source is active and function pointer is not assinged
    write(logUnit(6),*) 'WARNING: Dummy routine for applySrc'
  end subroutine mus_applySrc_dummy
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Dummy routine for add source to auxField
  subroutine mus_addSrcToAuxField_dummy(fun, auxField, iLevel, time, varSys, &
    &                                   phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> output auxField array
    real(kind=rk), intent(inout)         :: auxField(:)
    !> current level
    integer, intent(in)                :: iLevel
    !> current timing information
    type(tem_time_type), intent(in)    :: time
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> Physics conversion factor for current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos(:)
    ! ------------------------------------------------------------------------ !
    ! abort only if source is active and function pointer is not assinged
    write(logUnit(6),*) 'WARNING: Dummy routine for addSrcToAuxField'
  end subroutine mus_addSrcToAuxField_dummy
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Dummy routine for update source variable
  subroutine mus_updateSrcVar_dummy(fun, auxField, iLevel, varSys, phyConvFac, &
    &                               derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> input auxField array on current level
    real(kind=rk), intent(in)          :: auxField(:)
    !> current level
    integer, intent(in) :: iLevel
    !> variable system definition
    type(tem_varSys_type), intent(in)     :: varSys
    !> Physics conversion factor on current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in)  :: derVarPos(:)
    ! --------------------------------------------------------------------------
    write(logUnit(6),*) 'WARNING: Dummy routine for updateSourceVar'
  end subroutine mus_updateSrcVar_dummy
  ! ************************************************************************** !


end module mus_source_type_module
! **************************************************************************** !
