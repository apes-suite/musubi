! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2017, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014, 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015-2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!> author: Simon Zimny
!! This module contains datatypes and subroutines for the immersed boundary
!! method (IBM).
!!
!! @todo: IBM: some more information on IBM
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
!!
module mus_IBM_module

!$ use omp_lib

  ! include treelm modules
  use env_module,               only: rk, long_k, LabelLen, globalMaxLevels,   &
    &                                 tem_connect_toNull, newunit, PathLen
  use tem_spacetime_fun_module, only: tem_spacetime_fun_type,                  &
    &                                 tem_load_spacetime, tem_spacetime_for
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: tem_logging_type, logUnit, tem_last_lu,  &
    &                                 tem_logging_load
  use tem_surfaceData_module,   only: tem_surfData_type,                       &
    &                                 tem_load_surfData, tem_init_surfData,    &
    &                                 tem_readAndUnify_surfData,               &
    &                                 tem_freeSurfData,                        &
    &                                 tem_update_surfPos, tem_calcTriaAreas
  use tem_timer_module,         only: tem_labeledtimer_type, tem_startTimer,   &
    &                                 tem_stopTimer, tem_addTimer,             &
    &                                 tem_getMaxTimerVal, tem_resetTimer
  use tem_stlb_io_module,       only: tem_dump_stlb
  use tem_stencil_module,       only: tem_stencilHeader_type
  use treelmesh_module,         only: treelmesh_type
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_comm_module,          only: tem_commPattern_type,                    &
    &                                 tem_communication_type
  use tem_comm_env_module,      only: tem_comm_env_type
  use tem_param_module,         only: q__W, q__S, q__B, q__E, q__N, q__T
  use tem_grow_array_module,    only: grw_intArray_type, grw_realArray_type,   &
    &                                 grw_longArray_type,                      &
    &                                 init, append, destroy, empty
  use tem_dyn_array_module,     only: dyn_intArray_type, dyn_longArray_type,   &
    &                                 PositionOfVal, init, append, expand,     &
    &                                 destroy
  use tem_construction_module,  only: tem_levelDesc_type, tem_treeIDinTotal,   &
    &                                 tem_updateTree_properties
  use tem_math_module,          only: inamuroDelta3D
  use tem_geometry_module,      only: tem_baryOfId
  use tem_varSys_module,        only: tem_varSys_type
  use tem_general_module,       only: tem_general_type
  use tem_property_module,      only: prp_solid, prp_sendHalo
  use tem_time_module,          only: tem_time_type, tem_time_sim_stamp,       &
    &                                 tem_time_dump, tem_time_advance
  use tem_topology_module,      only: tem_CoordOfId, tem_IdOfCoord,            &
    &                                 tem_TreeIDComparison
  use tem_element_module,       only: eT_minRelevant, eT_maxRelevant,          &
    &                                 eT_fluid, eT_halo
  use tem_timeControl_module,   only: tem_timeControl_check
  use tem_simControl_module,    only: tem_simControl_type
  use tem_balance_module,       only: tem_balance_type

  ! include musubi modules
  use mus_param_module,              only: mus_param_type
  use mus_pdf_module,                only: pdf_data_type
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_derivedQuantities_module2, only: getVelocity, getDensity
  use mus_physics_module,            only: mus_convertFac_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_time_module,               only: mus_time_modulo
  use mus_timer_module,              only: mus_timerHandles

  ! include aotus modules
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aotus_module,     only: flu_State, aot_get_val, aoterr_Fatal,            &
    &                         aoterr_NonExistent, aoterr_WrongType

  implicit none

  private

  public :: mus_IBM_type
  public :: mus_IBM_globType
  public :: mus_load_IBM
  public :: mus_init_IBM
  public :: mus_inamuro_IBM
  public :: mus_IBM_setParentIDs
  public :: mus_finishIBM
  public :: mus_buildBuffIBM
  public :: mus_unload_IBM
  public :: mus_reload_IBM

  !> log units for the IBM
  integer, parameter :: nTimers = 14
  integer, save :: IBM_logUnit(0:tem_last_lu)

  integer, save :: mess_velX = 22
  integer, save :: mess_forceXk = 23
  integer, save :: mess_amountXk = 24
  integer, save :: mess_posXk = 25
  integer, save :: mess_pointsXk = 26
  integer, save :: mess_velXk = 27
  integer, save :: mess_posXk2 = 28
  integer, save :: mess_velXk2 = 29

  integer, save :: fileID = 1

  !> This datatype includes all necessary variables needed for 1 IBM step.
  !! It will be filled at every timestep.
  type mus_IBM_tmpData_type
    !> tmp array holding the body force at the lagrangian surface points Xk
    !! (Inamuro paper: \f$ g_{l}(X_{k}, t+\delta t) \f$)
    !! size: nPoints*3
    real(kind=rk), allocatable :: force_Xk(:)
    !> tmp array holding the body force at the eulerian grid elements X
    !! (Inamuro paper: \f$ g_{l}(x, t+\delta t) \f$)
    !! size: nPoints*stencil%QQ, 3
    real(kind=rk), allocatable :: force_X(:)
    !> tmp array holding the velocity at the lagrangian surface points Xk
    !! (Inamuro paper: \f$ u_{l}(X_{k}, t+\delta t) \f$)
    !! size: nPoints*3
    real(kind=rk), allocatable :: vel_Xk(:)
    !> tmp array holding the velocity at the eulerian grid elements X
    !! (Inamuro paper: \f$ u_{l}(x, t+\delta t) \f$)
    !! size: nPoints*stencil%QQ, 3
    real(kind=rk), allocatable :: vel_X(:)
    !> tmp array holding the initial velocity at the lagrangian points Xk
    !! (Inamuro paper: \f$ u*(X_{k}, t+\delta t) \f$)
    !! size: nPoints*3
    real(kind=rk), allocatable :: vel_Xk_ini(:)
    !> tmp array holding the initial velocity at the eulerian grid elements X
    !! (Inamuro paper: \f$ u*(x, t+\delta t) \f$)
    !! size: nPoints*stencil%QQ, 3
    real(kind=rk), allocatable :: vel_X_ini(:)
    !> tmp array holding the predef. velocity at the lagrangian surface points X
    !! (Inamuro paper: \f$ U_{k}(t+\delta t) \f$)
    !! size: nPoints*3
    real(kind=rk), allocatable :: vel_Xk_surf(:)
    !> temporary array storing pointers to the dynamic array of neighbor
    !! positions in the total treeID list of the parent elements of Xk
    !! size: nPoints * stencil%QQ
    integer, allocatable :: ptrs_neighs_Xk(:)
    !> temporary dynamic array for storing the actual neighbor positions of
    !! the Xk
!    type( dyn_intArray_type ) :: neighs_Xk
    !> tmp array for storing the neighbor positions in the total treeID list of
    !! the eulerian grid points
    type( grw_intArray_type ), allocatable :: neighs_x(:)
    !> tmp array for the result of the inamura delta function for the lagrangian
    !! points X
    !! size: nPoints * stencil%QQ
    real(kind=rk), allocatable :: inaDelta_Xk(:)
    !> tmp array for the result of the inamura delta function for the eulerian
    !! points X
    !! size: nElems (fluid, ghosts, halos)
    type( grw_realArray_type ), allocatable :: inaDelta_X(:)
    !> sendBuffer to communicate the force values on Xk (on this process)
    type( tem_communication_type ) :: IBMSend_Xk
    !> recvBuffer to communicate the force values on Xk (on this process)
    type( tem_communication_type ) :: IBMRecv_Xk
    !> tmp array of growing arrays storing the Xk positions
    type( grw_intArray_type ), allocatable :: posXk(:)
    !> sendBuffer to communicate the force values on X (on this process)
    type( tem_communication_type ) :: IBMSend_X
    !> recvBuffer to communicate the force values on X (on this process)
    type( tem_communication_type ) :: IBMRecv_X
    !> map from the global send proc array to the local one
    integer, allocatable :: map2globSend(:)
    !> map from the global recv proc array to the local one
    integer, allocatable :: map2globRecv(:)
    !> sendBuffer to communicate the pdf values on X (on this process)
    type( tem_communication_type ) :: IBMSend_X_pdf
    !> recvBuffer to communicate the pdf values on X (on this process)
    type( tem_communication_type ) :: IBMRecv_X_pdf
    !> map from the global send proc array to the local one
    integer, allocatable :: map2globSend_pdf(:)
    !> map from the global recv proc array to the local one
    integer, allocatable :: map2globRecv_pdf(:)
    !> tmp array of growing arrays storing the treeID position at the send proc
    type( grw_intArray_type ), allocatable :: treeIDs(:)
  end type mus_IBM_tmpData_type

  !> Datatype containing information on the immersed boundary method
  type mus_IBM_type
    !> is this IBM active?
    logical :: active = .false.
    !> surface data information incl. the filenames, point coordinates and
    !! corresponding triangle data
    type( tem_surfData_type ) :: surfData
    !> label for indentifying the type of IBM
    character(len=LabelLen) :: label
    !> position of the stencil in layout%stencil array
    integer :: stencilPos
    !> use the initial positions in the movement function
    !! or use the updated values
    logical :: useInitPos
    !> is the motion predefined?
    !! If a movement and velocity spacetime function is provided
    !! the new positions can be caluclated locally. This reduced the
    !! amount of communication.
    logical :: movPredef
    !> spacetime function type describing the movement of the points
    type( tem_spacetime_fun_type ) :: movement
    !> spacetime function type describing the velocity of movement
    type( tem_spacetime_fun_type ) :: velocity
    !> number of iterations for calculating the force
    integer :: nInaIters
    !> temporary dynamic array for storing the actual neighbor positions of
    !! the Xk
    type( dyn_intArray_type ) :: neighs_Xk
    !> number of local neighbor elements
    integer :: locNeighs_Xk = 0
    !> timer type for evaluating runtime in different routines
    type( tem_labeledtimer_type ) :: timings
    !> array of timer handles (definition in mus_init_IBM)
    integer, allocatable :: timerHandles(:)
    !> temporary data used
    type( mus_IBM_tmpData_type ) :: IBMData
  end type mus_IBM_type

  !> This datatype is a wrapper for the immersed boundary information
  !! and a possible logging type for debugging
  type mus_IBM_globType
    !> the immersed boundary data
    type( mus_IBM_type ), allocatable :: IBM(:)
    !> number of immersed boundaries
    integer :: nIBMs
    !> logging unit for IBM debug output to file
    type( tem_logging_type ) :: logIBM
  end type mus_IBM_globType


contains


! **************************************************************************** !
  !> Load the IBM information from the lua config file.
  !!
  !! The information has to be provided in the \a immersed_boundary \a table as
  !! a single table
  !!```lua
  !!  immersed_boundary = {
  !!                        logging = {
  !!                          level = 20,   -- level of logging
  !!                          filename = 'IBMdeb/test', -- file to write the log
  !!                          root_only = false -- should only root write msgs?
  !!                        },
  !!                        label = 'surface',
  !!                        kind = 'IBM',
  !!                        stlfiles = {'stl/surface1.stl', 'stl/surface2.stl'},
  !!                        dump_stl = {
  !!                          time_control = {
  !!                            min = {iter=0},
  !!                            max = {iter=tmax},
  !!                            interval = {iter=interval}
  !!                          },
  !!                          outprefix = 'test_',
  !!                        },
  !!                        inaIters = 4,
  !!                        movement = mov_pulse,
  !!                        velocity = vel_pulse
  !!  }
  !!```
  !! or as multiple tables
  !!```lua
  !!  immersed_boundary = {
  !!                        logging = {
  !!                          level = 20,   -- level of logging
  !!                          filename = 'IBMdeb/test', -- file to write the log
  !!                          root_only = false -- should only root write msgs?
  !!                        },
  !!                       {
  !!                        label = 'surface',
  !!                        kind = 'IBM',
  !!                        stlfiles = {'stl/surface1.stl', 'stl/surface2.stl'},
  !!                        dump_stl = {
  !!                          time_control = {
  !!                            min = {iter=0},
  !!                            max = {iter=tmax},
  !!                            interval = {iter=interval}
  !!                          },
  !!                          outprefix = 'test_',
  !!                        },
  !!                        inaIters = 4,
  !!                        movement = mov_pulse,
  !!                        velocity = vel_pulse
  !!                       }
  !!  }
  !!```
  !! in the mus_field table.
  subroutine mus_load_IBM( me, conf, rank )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_globType ), intent(inout) :: me
    !> handle of the lua config file
    type( flu_state ), intent(in) :: conf
    !> the current rank
    integer, intent(in) :: rank
    ! --------------------------------------------------------------------------
    ! handles for the ibm table and subtable
    integer :: ibm_handle, ibm_sub_handle, log_handle
    ! counter
    integer :: iIBM
    ! offset for the number of IBM in the immersed boundary table
    integer :: offset
    ! --------------------------------------------------------------------------

    ! open the IBM table
    call aot_table_open( L       = conf,                                       &
        &                thandle = ibm_handle,                                 &
        &                key     = 'immersed_boundary' )

    ! check wether ibm table is existing
    if( ibm_handle == 0 )then
      write(logUnit(1),*) 'No IBM found!'
      ! ... not existing close it again
      call aot_table_close( L       = conf,                                    &
        &                   thandle = ibm_handle )
      call tem_horizontalSpacer(fUnit=logUnit(10))
      me%nIBMs = 0
      allocate( me%IBM( me%nIBMs ))
      return
    else

      write(logUnit(3),*) 'Loading IBM ...'

      offset = 0

      ! open the logging table if available
      call aot_table_open( L       = conf,                                     &
        &                  parent  = ibm_handle,                               &
        &                  thandle = log_handle,                               &
        &                  key     = 'logging' )

      if( log_handle > 0 )then

        ! load the logging type for the IBM
        call tem_logging_load( conf    = conf,                                 &
          &                    thandle = log_handle,                           &
          &                    rank    = rank,                                 &
          &                    me      = me%logIBM )
        offset = 1
        ! now set the module variable IBM_logUnits for the IBM logging type
        IBM_logUnit = me%logIBM%fUnit
      else
        ! connect the module variable to the null device
        call tem_connect_toNull(IBM_logUnit(0))
        IBM_logUnit(1:tem_last_lu) = IBM_logUnit(0)
      end if

      call aot_table_close( L = conf, thandle = log_handle )

      ! ... an ibm table exists now check wether it is a single ibm or
      !     a table of ibms
      ! try to open the first subtable ...
      call aot_table_open( L       = conf,                                     &
          &                parent  = ibm_handle,                               &
          &                thandle = ibm_sub_handle,                           &
          &                pos     = 1 )

      if( ibm_sub_handle == 0 )then
        ! ... if the subtable does not exist close it again
        call aot_table_close( L       = conf,                                  &
          &                   thandle = ibm_sub_handle )

        ! allocate the ibm array with one
        me%nIBMs = 1
        allocate( me%IBM( me%nIBMs ))

        write(logUnit(5),*) 'IBM is a single table'

        ! ... read the data
        call mus_load_IBM_single( me         = me%IBM(1), &
          &                       conf       = conf,      &
          &                       sub_handle = ibm_handle )
      else
        ! ... subtable is existing close the table
        call aot_table_close( L       = conf,                                  &
          &                   thandle = ibm_sub_handle )
        ! ... get the table length
        me%nIBMs = aot_table_length( L=conf, thandle=ibm_handle ) - offset

        ! ... allocate the ibm array
        allocate( me%IBM( me%nIBMs ))

        write(logUnit(3),*) 'Number of IBMs: ', me%nIBMs

        ! ... loop over the entries and ...
        do iIBM = 1, me%nIBMs
          write(logUnit(5),*) 'Loading IBM ', iIBM
          ! ... open the table for the individual positions
          call aot_table_open( L       = conf,                                 &
              &                parent  = ibm_handle,                           &
              &                thandle = ibm_sub_handle,                       &
              &                pos     = iIBM )

          ! ... read the data
          call mus_load_IBM_single( me         = me%IBM(iIBM),  &
            &                       conf       = conf,          &
            &                       sub_handle = ibm_sub_handle )

          ! ... and close the table again
          call aot_table_close( L       = conf,                                &
            &                   thandle = ibm_sub_handle )

          write(logUnit(5),*) 'Done IBM ', iIBM
        end do
      end if
    end if
    call aot_table_close( L       = conf,                                      &
      &                   thandle = ibm_handle )
    call tem_horizontalSpacer(fUnit=logUnit(10))

  end subroutine mus_load_IBM
! **************************************************************************** !


! **************************************************************************** !
  !> Load a single IBM table from the config file.
  !!
  subroutine mus_load_IBM_single( me, conf, sub_handle )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(inout) :: me
    !> handle of the lua config file
    type( flu_state ), intent(in) :: conf
    !> handle for the surfaceData table
    integer, intent(in) :: sub_handle
    ! --------------------------------------------------------------------------
    ! error variable
    integer :: iError
    ! --------------------------------------------------------------------------

    ! read the surface data
    call tem_load_surfData( me        = me%surfData,                           &
      &                     conf      = conf,                                  &
      &                     sd_handle = sub_handle )

    ! load the number of iterations
    call aot_get_val( L       = conf,                                          &
      &               thandle = sub_handle,                                    &
      &               val     = me%nInaIters,                                  &
      &               ErrCode = iError,                                        &
      &               key     = 'inaIters',                                    &
      &               default = 1)

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(0),*)'FATAL Error occured, while retrieving inaIters:'
      if (btest(iError, aoterr_NonExistent))then
        write(logUnit(0),*)'Variable not existent!'
        write(logUnit(0),*)' Using the default value 1'
      else if (btest(iError, aoterr_WrongType))then
        write(logUnit(0),*)'Variable has wrong type!'
        write(logUnit(0),*)'STOPPING'
        call tem_abort()
      end if
    end if

    ! load wether to use the initial or updated positions when calculating
    ! the new positions
    call aot_get_val( L       = conf,                                          &
      &               thandle = sub_handle,                                    &
      &               val     = me%useInitPos,                                 &
      &               ErrCode = iError,                                        &
      &               key     = 'useInitPos',                                  &
      &               default = .false.)

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(0),*)'FATAL Error occured, while retrieving useInitPos'
      if (btest(iError, aoterr_NonExistent))then
        write(logUnit(0),*)'Variable not existent!'
        write(logUnit(0),*)' Using the default value 1'
      else if (btest(iError, aoterr_WrongType))then
        write(logUnit(0),*)'Variable has wrong type!'
        write(logUnit(0),*)'STOPPING'
        call tem_abort()
      end if
    end if

    ! load wether to use the initial or updated positions when calculating
    ! the new positions
    call aot_get_val( L       = conf,                                          &
      &               thandle = sub_handle,                                    &
      &               val     = me%movPredef,                                  &
      &               ErrCode = iError,                                        &
      &               key     = 'movPredef',                                   &
      &               default = .true.)

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(0),*)'FATAL Error occured, while retrieving movPredef:'
      if (btest(iError, aoterr_NonExistent))then
        write(logUnit(0),*)'Variable not existent!'
        write(logUnit(0),*)' Using the default value 1'
      else if (btest(iError, aoterr_WrongType))then
        write(logUnit(0),*)'Variable has wrong type!'
        write(logUnit(0),*)'STOPPING'
        call tem_abort()
      end if
    end if

    ! not predefined motion is still buggy at the moment!!!
    if( .not. me%movPredef )then
      write(logUnit(0),*)'ATTENTION: Not predefining the motion is NOT '//     &
        &                'working as expected. To fix this the surface '//     &
        &                'areas need to be communicated. ABORTING!!!!'
      call tem_abort()
    end if

    ! Load the spacetime function for the movement (xyz-component)
    call tem_load_spacetime( me     = me%movement, &
      &                      conf   = conf,        &
      &                      parent = sub_handle,  &
      &                      key    = 'movement',  &
      &                      nComp  = 3            )

    ! @todo: SZ: a spacetime function is not a must
    ! check if the spacetime function has been read
    if(trim(me%movement%fun_kind)=='none')then
      write(logUnit(0),*) 'The spacetime function for the movement ' //        &
        &                 'has not been loaded properly.'
      call tem_abort()
    end if

    ! Load the spacetime function for the velocity (xyz-component)
    call tem_load_spacetime( me     = me%velocity, &
      &                      conf   = conf,        &
      &                      parent = sub_handle,  &
      &                      key    = 'velocity',  &
      &                      nComp  = 3            )

    ! @todo: SZ: a spacetime function is not a must
    ! check if the spacetime function has been read
    if(trim(me%velocity%fun_kind)=='none')then
      write(logUnit(0),*) 'The spacetime function for the velocity ' //        &
        &                 'has not been loaded properly.'
      call tem_abort()
    end if

    allocate(me%timerHandles(nTimers))

    ! set this IBM to be active at all
    me%active = .true.

  end subroutine mus_load_IBM_single
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine initializes the IBM data incl. reading the stl, unifying
  !! the coordinates, storing the connectivity, allocating the parentIDs array
  !! and initializing the stencil used.
  !!
  subroutine mus_init_IBM( me, globTree )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_globType ), intent(inout) :: me
    !> global tree information
    type( treelmesh_type ), intent(in) :: globTree
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iIBM
    ! --------------------------------------------------------------------------

    ! loop over the IBMs
    do iIBM = 1, me%nIBMs
      ! ... read and unify the surface data
      call tem_readAndUnify_surfData( me         = me%IBM(iIBM)%surfData,      &
        &                             useInitPos = me%IBM(iIBM)%useInitPos )
      ! ... allocate the array of parent IDs with the number of levels
      allocate( me%IBM( iIBM )%surfData%parentIDs( globTree%global%minLevel:   &
        &                                          globTree%global%maxLevel ))

      ! add timers for debugging:
      ! 1. initialization timer incl. allocation
      call tem_addTimer( me          = me%IBM( iIBM )%timings,         &
        &                timerHandle = me%IBM( iIBM )%timerHandles(1), &
        &                timerName   = 'IBMinit' )
      ! 2. build_Xk communicator timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,         &
        &                timerHandle = me%IBM( iIBM )%timerHandles(2), &
        &                timerName   = 'buildXk' )
      ! 3. build_X communicator timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,         &
        &                timerHandle = me%IBM( iIBM )%timerHandles(3), &
        &                timerName   = 'buildX' )
      ! 4. IBM loop timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,         &
        &                timerHandle = me%IBM( iIBM )%timerHandles(4), &
        &                timerName   = 'IBMloop' )
      ! 5. post loop timer incl. deallocation
      call tem_addTimer( me          = me%IBM( iIBM )%timings,         &
        &                timerHandle = me%IBM( iIBM )%timerHandles(5), &
        &                timerName   = 'IBMfinish' )

      ! 6. init pre comm timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,         &
        &                timerHandle = me%IBM( iIBM )%timerHandles(6), &
        &                timerName   = 'IBMInit_pre' )
      ! 7. init post comm timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,         &
        &                timerHandle = me%IBM( iIBM )%timerHandles(7), &
        &                timerName   = 'IBMInit_post' )
      ! 8. buidXk pre comm timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,         &
        &                timerHandle = me%IBM( iIBM )%timerHandles(8), &
        &                timerName   = 'buildXk_pre' )
      ! 9. buildXk post comm timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,         &
        &                timerHandle = me%IBM( iIBM )%timerHandles(9), &
        &                timerName   = 'buildXk_post' )
      ! 10. buidX pre comm timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,          &
        &                timerHandle = me%IBM( iIBM )%timerHandles(10), &
        &                timerName   = 'buildX_pre' )
      ! 11. buildX search in buffers timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,          &
        &                timerHandle = me%IBM( iIBM )%timerHandles(11), &
        &                timerName   = 'buildX_buff' )
      ! 12. buildX post comm timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,          &
        &                timerHandle = me%IBM( iIBM )%timerHandles(12), &
        &                timerName   = 'buildX_post' )
      ! 13. initPre update timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,          &
        &                timerHandle = me%IBM( iIBM )%timerHandles(13), &
        &                timerName   = 'IBMInit_Pre_update' )
      ! 14. initPre surf vel init IBMData timer
      call tem_addTimer( me          = me%IBM( iIBM )%timings,          &
        &                timerHandle = me%IBM( iIBM )%timerHandles(14), &
        &                timerName   = 'IBMInit_Pre_initVel' )

    end do

  end subroutine mus_init_IBM
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine sets the positions of the parent IDs in the level
  !! descriptor.
  !!
  subroutine mus_IBM_setParentIDs( nIBMs, me, levelDesc, tree )
    ! --------------------------------------------------------------------------
    !> number of IBM types
    integer, intent(in) :: nIBMs
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(inout) :: me(:)
    !> the level descriptor incl. ghost and halo elements as well as the
    !! communicator information on the level iLevel
    type( tem_levelDesc_type ), intent(inout) :: levelDesc(:)
    !> global Tree information
    type( treelmesh_type ), intent(inout) :: tree
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iIBM
    integer :: iLevel
    ! --------------------------------------------------------------------------

    ! loop over the IBMs
    do iIBM = 1, nIBMs
      write(logUnit(1),*) 'Initializing IBM data ...'
      ! loop over the levels
      do iLevel = tree%global%minLevel, tree%global%maxLevel
        ! initialize the parent ID positions on the different levels
        call tem_init_surfData( me        = me(iIBM)%surfData,            &
          &                     levelDesc = levelDesc(iLevel),            &
          &                     globtree  = tree,                         &
          &                     iLevel    = tree%global%minLevel-1+iLevel )
        ! and update the properties in the tree
        call tem_updateTree_properties( levelDesc = levelDesc(iLevel), &
          &                             tree      = tree )
      end do
      write(logUnit(1),*) 'Done initializing IBM data.'
    end do

  end subroutine mus_IBM_setParentIDs
! **************************************************************************** !


! **************************************************************************** !
  !>
  !!
  subroutine mus_buildBuffIBM( me, commPattern, globTree, params, layout,      &
    &                          levelDesc, iLevel )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(inout) :: me(:)
    !> communication pattern
    type( tem_commPattern_type ), intent(inout) :: commPattern
    !> global tree information
    type( treelmesh_type ) :: globTree
    !> global parameters
    type( mus_param_type ) :: params
    !> scheme layout of the current scheme incl. array of stencils
    type( mus_scheme_layout_type ) :: layout
    !> the level descriptor incl. ghost and halo elements as well as the
    !! communicator information on the level iLevel
    type( tem_levelDesc_type ), intent(inout) :: levelDesc
    !> the current level
    integer, intent(in) :: iLevel
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iIBM
    ! tmp variable for the total number of neighbors: nUniquePoints*QQ
    integer :: tmp_totNeighs
    ! local sim control
    type( tem_simControl_type ) :: loc_simControl
    integer :: loc_level
    ! --------------------------------------------------------------------------
    loc_level = iLevel

    ! advance the local time control type
    call tem_time_advance( me     = loc_simControl%now,                        &
      &                    sim_dt = params%physics%dtLvl(                      &
      &                                     globTree%global%maxLevel ))

    loc_simControl = params%general%simControl

    write(IBM_logUnit(3),*)'Starting to build communicators for all IBM ' //   &
      &                    'for time: '                                   //   &
      &                    trim(tem_time_sim_stamp(loc_simControl%now))

    ! loop over the IBMs in this field
    do iIBM = 1, size(me)

      ! start the IBMinit timer
      call tem_startTimer( me          = me( iIBM )%timings%timedat, &
        &                  timerHandle = me( iIBM )%timerHandles(1)  )

      ! start the IBMinit timer
      call tem_startTimer( me          = me( iIBM )%timings%timedat, &
        &                  timerHandle = me( iIBM )%timerHandles(6)  )

      ! start the IBMinit timer
      call tem_startTimer( me          = me( iIBM )%timings%timedat, &
        &                  timerHandle = me( iIBM )%timerHandles(14) )

      call destroy( me( iIBM )%neighs_Xk )

      write(IBM_logUnit(3),*)'Starting to process IBM ', iIBM, '...'

      ! set the total number of neighbors
      tmp_totNeighs = me(iIBM)%surfData%nUniquePoints_total                    &
        &           * layout%stencil( me(iIBM)%stencilPos )%QQ

      call mus_init_IBMData( me           = me( iIBM )%IBMData,                &
        &                    IBM          = me( iIBM ),                        &
        &                    totNeighs    = tmp_totNeighs,                     &
        &                    nElems_fluid = levelDesc%elem%nElems( eT_fluid ) )

      write(IBM_logUnit(3),*)'   Calculating the new surface position and '//  &
        &                 'updating the surface Areas...'

      ! calculate the surface velocity based on the old Xk
      call mus_IBM_getSurfVel( me        = me(iIBM),                           &
        &                      IBMData   = me(iIBM)%IBMData,                   &
        &                      levelDesc = levelDesc,                          &
        &                      general   = params%general,                     &
        &                      iLevel    = iLevel )

      ! stop the IBMinit timer
      call tem_stopTimer( me          = me( iIBM )%timings%timedat, &
        &                 timerHandle = me( iIBM )%timerHandles(14) )

      ! start the IBMinit timer
      call tem_startTimer( me          = me( iIBM )%timings%timedat, &
        &                  timerHandle = me( iIBM )%timerHandles(13) )

      ! compute the new positions of the surface data points, assign the correct
      ! properties and update the pointers to the parent IDs
      call tem_update_surfPos( me         = me(iIBM)%surfData,                 &
        &                      levelDesc  = levelDesc,                         &
        &                      globTree   = globTree,                          &
        &                      movement   = me(iIBM)%movement,                 &
        &                      time       = loc_simControl%now,                &
        &                      iLevel     = loc_level,                         &
        &                      IBMUnit    = IBM_logUnit,                       &
        &                      useInitPos = me(iIBM)%useInitPos,               &
        &                      movPredef  = me(iIBM)%movPredef )

      ! stop the IBMinit timer
      call tem_stopTimer( me          = me( iIBM )%timings%timedat, &
        &                 timerHandle = me( iIBM )%timerHandles(13) )

      ! stop the IBMinit timer
      call tem_stopTimer( me          = me( iIBM )%timings%timedat, &
        &                 timerHandle = me( iIBM )%timerHandles(1)  )

      ! stop the IBMinit timer
      call tem_stopTimer( me          = me( iIBM )%timings%timedat, &
        &                 timerHandle = me( iIBM )%timerHandles(6)  )

      ! build up the send and receive buffers for the new coordinates and
      ! force terms
      ! For this:
      ! @todo: IBM: Add description
      !
      ! -> buffers are ready for this timestep and can be filled with
      !    coordinates of force terms on Xk (3 reals per Xk which needs to be
      !    communicated)
      !
      call mus_IBM_buildSendRecv_Xk( me          = me(iIBM),                   &
        &                            IBMData     = me(iIBM)%IBMData,           &
        &                            levelDesc   = levelDesc,                  &
        &                            commPattern = commPattern,                &
        &                            globTree    = globTree,                   &
        &                            iLevel      = iLevel,                     &
        &                            params      = params )

      write(IBM_logUnit(3),*)'   Done building the communicator and '//        &
        &                    'communicating the new positions!'

      write(IBM_logUnit(3),*)'   Filling the IBM neighbor array AND '//        &
        &                    ' building the send and recv buffers for X ...'

      ! start the IBMbuildX timer
      call tem_startTimer( me          = me( iIBM )%timings%timedat, &
        &                  timerHandle = me( iIBM )%timerHandles(3)  )

      ! prepare the buffers for the eulerian elements
      call mus_IBM_prepareSendRecv_X( IBMData = me( iIBM )%IBMData )

      ! stop the IBMbuildX timer
      call tem_stopTimer( me          = me( iIBM )%timings%timedat, &
        &                 timerHandle = me( iIBM )%timerHandles(3)  )

      ! build the neighbor array filled with all neighbors (unique) of all
      ! parentIDs
      ! AND initialize the send and receive buffers for the eulerian elements X
      call mus_fillNeigh_surfData( me          = me(iIBM),                     &
        &                          IBMData     = me( iIBM )%IBMData,           &
        &                          stencil     = layout%stencil,               &
        &                          levelDesc   = levelDesc,                    &
        &                          parentIDs   = me(iIBM)%surfData%            &
        &                                          parentIDs(iLevel)%ptrs,     &
        &                          globTree    = globTree,                     &
        &                          commPattern = commPattern,                  &
        &                          params      = params )

    end do

  end subroutine mus_buildBuffIBM
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine modifies the state vector according to the method described
  !! in the paper \a Lift generation by a two-dimensional symmetric flapping
  !! wing: immersed boundary-lattice Boltzmann simulations \a by Inamuro et al.
  !! @cite Ota:2012bx .
  !!
  subroutine mus_inamuro_IBM( me, commPattern, globTree, general, pdf, layout, &
    &                         levelDesc, globSys, stateVarMap, convFac,        &
    &                         iField, iLevel, state                            )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(inout) :: me(:)
    !> communication pattern
    type( tem_commPattern_type ), intent(inout) :: commPattern
    !> global tree information
    type( treelmesh_type ) :: globTree
    !> general data
    type( tem_general_type ), intent(in) :: general
    !> pdf_data_type incl. connectivity array on all levels
    type( pdf_data_type ), intent(inout) :: pdf
    !> scheme layout of the current scheme incl. array of stencils
    type( mus_scheme_layout_type ) :: layout
    !> the level descriptor incl. ghost and halo elements as well as the
    !! communicator information on the level iLevel
    type( tem_levelDesc_type ), intent(inout) :: levelDesc
    !> global variable system of the current scheme
    type(tem_varSys_type) :: globSys
    !> Position of state variables in globSys
    integer, intent(in) :: stateVarMap(:)
    !> conversion factors
    type(mus_convertFac_type), intent(in) :: convFac
    !> the current field
    integer, intent(in) :: iField
    !> the current level
    integer, intent(in) :: iLevel
    !> state_data type
    real(kind=rk), intent(inout) :: state(:,:)
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iIBM, iIter, iProc
    ! shall the stl be dumped?
    logical :: triggered

    integer :: forceunit
    character( len=LabelLen ) :: forceName
    character(len=7) :: rankstamp
    ! timestamp for the filename
    character(len=16)     :: timeStamp
    logical :: writeForce = .false.
    real(kind=rk) :: maxForce
    real(kind=rk) :: maxVel_X
    real(kind=rk) :: maxVel_X_ini
    real(kind=rk) :: maxVel_Xk
    real(kind=rk) :: maxForce_xk
    ! --------------------------------------------------------------------------

    forceunit = newunit()
!! !$omp single

    write(IBM_logUnit(3),*)'Starting to process all IBM for time: '//          &
      &                  trim(tem_time_sim_stamp(general%simControl%now))

!! !$omp end single

    ! loop over the IBMs in this field
    do iIBM = 1, size(me)

!! !$omp single

      call mus_IBMFinishBuff( me          = me(iIBM),          &
          &                   IBMData     = me(iIBM)%IBMData,  &
          &                   levelDesc   = levelDesc,         &
          &                   commPattern = commPattern,       &
          &                   globTree    = globTree,          &
          &                   iLevel      = iLevel,            &
          &                   comm        = general%proc%comm, &
          &                   stencil     = layout%stencil)

      ! start the dyn load balancing timer after the communication
      call tem_startTimer( timerHandle = mus_timerHandles%doIBM(iLevel))

      ! start the IBMinit timer
      call tem_startTimer( me          = me( iIBM )%timings%timedat, &
        &                  timerHandle = me( iIBM )%timerHandles(1)  )

      ! start the IBMinit timer
      call tem_startTimer( me          = me( iIBM )%timings%timedat, &
        &                  timerHandle = me( iIBM )%timerHandles(7)  )

      write(IBM_logUnit(3),*)'   Updating the surface areas ...'

      ! calculate the area for the different surface points (1/3 of the sum of
      ! attached triangle areas)
      call tem_calcTriaAreas( me = me(iIBM)%surfData )

      write(IBM_logUnit(3),*)'   Done updating the surface areas!'


      call tem_timeControl_check(                                              &
        &               me     = me( iIBM )%surfData%timeControl,              &
        &               now    = general%simControl%now,                &
        &               comm   = general%proc%comm,                     &
        &               triggered = triggered )

      ! @todo: IBM: use a sim_control table to set timings%timedat
      !             in which the stl should be dumped
      ! check convergence if timeControl is triggered
      if( triggered &
       ! ... then make sure it is really the last iteration due to recursiveness
        !& mus_time_modulo( general%simControl%now, reqInterval ) &
        & )  then
        writeForce = me(iIBM)%surfData%dumpForce
        if( me(iIBM)%surfData%dumpForce )then
          write(rankstamp, '(a1,I6.6)') '.',general%proc%rank
          write(forceName,'(a)')'IBMdeb/force'
          timestamp = trim(tem_time_sim_stamp(general%simControl%now))
          open( unit = forceunit,                                              &
            &   file = trim(forceName)//trim(rankstamp)//'_'//trim(timestamp))
        end if
        call tem_dump_stlb( outprefix = trim(me(iIBM)%surfData%outprefix),     &
          &                 nodes     = me(iIBM)%surfData%pointCoords,         &
          &                 triangles = me(iIBM)%surfData%trias,               &
          &                 proc      = general%proc,                          &
          &                 time      = general%simControl%now )
      end if

      write(IBM_logUnit(3),*)'   Filling the initial force and velocity '//    &
        &                    'arrays ...'

!! !$omp end single

      ! calculate the initial force, initial velocity and the surface velocity
      ! at the lagrangian surface points Xk from the local information
      call mus_inamuroIni( me          = me(iIBM),                             &
        &                  IBMData     = me(iIBM)%IBMData,                     &
        &                  state       = state(:, pdf%nNext),                  &
        &                  globtree    = globTree,                             &
        &                  stencil     = layout%stencil,                       &
        &                  stencilPos  = me(iIBM)%stencilPos,                  &
        &                  levelDesc   = levelDesc,                            &
        &                  globSys     = globSys,                              &
        &                  pdfPos      = stateVarMap(iField),                  &
        &                  iLevel      = iLevel,                               &
        &                  convFac     = convFac                               )

!! !$omp single

      write(IBM_logUnit(3),*)'   Done calculating the initial velocity and '// &
        &                    'force!'

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!
!write(IBM_logUnit(3),*)'inaDelta_Xk after ini'
!do iProc = 1, tmp_totNeighs
!  write(IBM_logUnit(3),*)IBMData%inaDelta_Xk( iProc)
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logunit(1),*)'IBMSend_Xk: '
!write(IBM_logunit(1),*)'nProcs: ', IBMData%IBMSend_Xk%nProcs
!do iProc = 1, IBMData%IBMSend_Xk%nProcs
!  write(IBM_logunit(1),*)'  proc   ', IBMData%IBMSend_Xk%proc(iProc)
!  write(IBM_logunit(1),*)'  nVals: ', &
!    &                    IBMData%IBMSend_Xk%buf_real(iProc)%nVals
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logunit(1),*)'IBMRecv_Xk: '
!write(IBM_logunit(1),*)'nProcs: ', IBMData%IBMRecv_Xk%nProcs
!do iProc = 1, IBMData%IBMRecv_Xk%nProcs
!  write(IBM_logunit(1),*)'  proc   ', IBMData%IBMRecv_Xk%proc(iProc)
!  write(IBM_logunit(1),*)'  nVals: ', &
!    &                    IBMData%IBMRecv_Xk%buf_real(iProc)%nVals
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logunit(1),*)'IBMSend_X: '
!write(IBM_logunit(1),*)'nProcs: ', IBMData%IBMSend_X%nProcs
!do iProc = 1, IBMData%IBMSend_X%nProcs
!  write(IBM_logunit(1),*)'  proc   ', IBMData%IBMSend_X%proc(iProc)
!  write(IBM_logunit(1),*)'  nVals: ', &
!    &                    IBMData%IBMSend_X%buf_real(iProc)%nVals
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logunit(1),*)'IBMRecv_X: '
!write(IBM_logunit(1),*)'nProcs: ', IBMData%IBMRecv_X%nProcs
!do iProc = 1, IBMData%IBMRecv_X%nProcs
!  write(IBM_logunit(1),*)'  proc   ', IBMData%IBMRecv_X%proc(iProc)
!  write(IBM_logunit(1),*)'  nVals: ', &
!    &                    IBMData%IBMRecv_X%buf_real(iProc)%nVals
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logunit(1),*)'GlobSend: '
!write(IBM_logunit(1),*)'nProcs: ', levelDesc%sendbuffer%nProcs
!do iProc = 1, levelDesc%sendbuffer%nProcs
!  write(IBM_logunit(1),*)'  proc   ', levelDesc%sendbuffer%proc(iProc)
!  write(IBM_logunit(1),*)'  nVals: ', &
!    &                    levelDesc%sendbuffer%buf_real(iProc)%nVals
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logunit(1),*)'Globrecv: '
!write(IBM_logunit(1),*)'nProcs: ', levelDesc%recvbuffer%nProcs
!do iProc = 1, levelDesc%recvbuffer%nProcs
!  write(IBM_logunit(1),*)'  proc   ', levelDesc%recvbuffer%proc(iProc)
!  write(IBM_logunit(1),*)'  nVals: ', &
!    &                    levelDesc%recvbuffer%buf_real(iProc)%nVals
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

      ! -----------------------------------------------------------------------!
      !      S T A R T I N G   T H E  I N A M U R O  I T E R A T I O N S       !
      ! -----------------------------------------------------------------------!

      ! stop the IBMinit timer
      call tem_stopTimer( me          = me( iIBM )%timings%timedat, &
        &                 timerHandle = me( iIBM )%timerHandles(1)  )

      ! stop the IBMinit timer
      call tem_stopTimer( me          = me( iIBM )%timings%timedat, &
        &                 timerHandle = me( iIBM )%timerHandles(7)  )

      ! start the IBMloop timer
      call tem_startTimer( me          = me( iIBM )%timings%timedat, &
        &                  timerHandle = me( iIBM )%timerHandles(4)  )

      write(IBM_logUnit(3),*) '   Starting the Inamuro iterations ',           &
        &                   me(iIBM)%nInaIters, ' times ...'

!! !$omp end single

      ! loop over the inamuro iterations
      do iIter = 1, me(iIBM)%nInaIters

!! !$omp single

    maxForce = 0._rk
    maxVel_X = 0._rk
    maxVel_X_ini = 0._rk
    maxVel_Xk = 0._rk
    maxforce_Xk = 0._rk

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'force_Xk before comm'
!do iProc = 1, me(iIBM)%surfData%nUniquePoints_total
!  write(IBM_logUnit(3),*)IBMData%force_Xk( (iProc-1)*3+1:(iProc-1)*3+3 )
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
        write(IBM_logUnit(3),*) '      Starting with Inamuro iteration ', iIter

        ! stop the dyn load balancing timer since communication spoils the
        ! result
        call tem_stopTimer( timerHandle = mus_timerHandles%doIBM(iLevel) )

        write(IBM_logUnit(3),*) '      Communicating the force values ...'

        ! ... communicate the force values for Xk
        call commPattern%exchange_real(                                    &
          &                    send         = me(iIBM)%IBMData%IBMSend_Xk, &
          &                    recv         = me(iIBM)%IBMData%IBMRecv_Xk, &
          &                    state        = me(iIBM)%IBMData%force_Xk,   &
          &                    message_flag = mess_forceXk,                &
          &                    comm         = general%proc%comm            )

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'exchange force_Xk'
!do iProc = 1, IBMData%IBMSend_Xk%nProcs
!  write(IBM_logUnit(3),*)'send: ', iProc
!  write(IBM_logUnit(3),*)IBMData%IBMSend_Xk%buf_real(iProc)%val(              &
!    &                                1:IBMData%IBMSend_Xk%elemPos(iProc)%nVals)
!end do
!do iProc = 1, IBMData%IBMrecv_Xk%nProcs
!  write(IBM_logUnit(3),*)'recv: ', iProc
!  write(IBM_logUnit(3),*)IBMData%IBMrecv_Xk%buf_real(iProc)%val(              &
!    &                                1:IBMData%IBMrecv_Xk%elemPos(iProc)%nVals)
!end do
!write(IBM_logUnit(3),*)'force_Xk after comm'
!do iProc = 1, me(iIBM)%surfData%nUniquePoints_total
!  write(IBM_logUnit(3),*)IBMData%force_Xk( (iProc-1)*3+1:(iProc-1)*3+3 )
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!
!write(IBM_logUnit(3),*)'force_X before calc'
!do iProc = 1, me(iIBM)%neighs_Xk%nVals
!  write(IBM_logUnit(3),*)IBMData%force_X( (iProc-1)*3+1:(iProc-1)*3+3 )
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

        write(IBM_logUnit(3),*) '      Done communicating the force values!'

        ! start the dyn load balancing timer
        call tem_startTimer( timerHandle = mus_timerHandles%doIBM(iLevel) )

        write(IBM_logUnit(3),*) '      Calculating the new forces at the '//   &
          &                     'grid points ...'
!! !$omp end single
        ! calculate the new forces at the grid points
        call mus_calcForce_X( me           = me(iIBM),                         &
          &                   IBMData      = me(iIBM)%IBMData,                 &
          &                   nElems_fluid = levelDesc%elem%nElems( eT_fluid ),&
          &                   convFac      = convFac )

!write(IBM_logUnit(3),*)'force_X after calc'
!do iProc = 1, me(iIBM)%neighs_Xk%nVals
!  write(IBM_logUnit(3),*)IBMData%force_X( (iProc-1)*3+1:(iProc-1)*3+3 )
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

!! !$omp single
        write(IBM_logUnit(3),*) '      Done calculating the new forces '//     &
          &                     'at the grid points!'


        write(IBM_logUnit(3),*) '      Correct the velocity at the grid '//    &
          &                     'points ...'

!write(IBM_logUnit(3),*)'vel_X before correct: '
!do iProc = 1, me(iIBM)%neighs_Xk%nVals
!  write(IBM_logUnit(3),*)IBMData%vel_X( (iProc-1)*3+1:(iProc-1)*3+3)
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

!! !$omp end single
        ! correct the velocity at the grid points
        call mus_corrVel_X( me        = me( iIBM ),                            &
          &                 IBMData   = me(iIBM)%IBMData,                      &
          &                 convFac   = convFac )

!! !$omp single
        write(IBM_logUnit(3),*) '      Done correcting the velocity!'

!write(IBM_logUnit(3),*)'vel_X after correct: '
!do iProc = 1, me(iIBM)%neighs_Xk%nVals
!  write(IBM_logUnit(3),*)IBMData%vel_X( (iProc-1)*3+1:(iProc-1)*3+3)
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'exchange vel_X'
!write(IBM_logUnit(3),*)'nProcs: ', IBMData%IBMSend_X%nProcs
!do iProc = 1, IBMData%IBMSend_X%nProcs
!  write(IBM_logUnit(3),*)'send:  ', iProc
!  write(IBM_logUnit(3),*)'proc:  ', IBMData%IBMSend_X%proc( iProc )
!  write(IBM_logUnit(3),*)'nVals: ', IBMData%IBMSend_X%elemPos(iProc)%nVals
!end do
!write(IBM_logUnit(3),*)'nProcs: ', IBMData%IBMRecv_X%nProcs
!do iProc = 1, IBMData%IBMrecv_X%nProcs
!  write(IBM_logUnit(3),*)'recv:  ', iProc
!  write(IBM_logUnit(3),*)'proc:  ', IBMData%IBMrecv_X%proc( iProc )
!  write(IBM_logUnit(3),*)'nVals: ',IBMData%IBMrecv_X%elemPos(iProc)%nVals
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

!flush( IBM_logUnit(3) )
!!
!call mpi_barrier( general%proc%comm, ierr )
!
!call tem_abort()


        ! stop the dyn load balancing timer
        call tem_stopTimer( timerHandle = mus_timerHandles%doIBM(iLevel) )

        write(IBM_logUnit(3),*) '      Communicating the velocity on the '//   &
          &                     'eulerian elements ...'

        ! ... communicate the force values for X
        call commPattern%exchange_real(                                        &
          &                    send         = me(iIBM)%IBMData%IBMSend_X,      &
          &                    recv         = me(iIBM)%IBMData%IBMRecv_X,      &
          &                    state        = me(iIBM)%IBMData%vel_X,          &
          &                    message_flag = mess_velX,                       &
          &                    comm         = general%proc%comm )

        ! start the dyn load balancing timer
        call tem_startTimer( timerHandle = mus_timerHandles%doIBM(iLevel) )



        write(IBM_logUnit(3),*) '      Done communicating the velocity X!'

        write(IBM_logUnit(3),*) '      Interpolate the velocity at the '//     &
          &                     'surface points ...'

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'vel_xk before intp'
!do iProc = 1, me(iIBM)%surfData%nUniquePoints_total
!  write(IBM_logUnit(3),*)IBMData%vel_xk( (iProc-1)*3+1:(iProc-1)*3+3 )
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

!! !$omp end single
        ! interpolate the velocity at the surface points
        call mus_intpVel_Xk( IBMData      = me(iIBM)%IBMData,                  &
          &                  nNeighs      = layout%stencil(                    &
          &                                           me(iIBM)%stencilPos)%QQ, &
          &                  convFac      = convFac,                           &
          &                  nPoints      = me(iIBM)%surfData%                 &
          &                                               nUniquePoints_total, &
          &                  parentIDs    = me(iIBM)%surfData%                 &
          &                                            parentIDs(iLevel)%ptrs, &
          &                  nElems_fluid = levelDesc%elem%nElems( eT_fluid ) )
!! !$omp single

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'vel_xk after intp'
!do iProc = 1, me(iIBM)%surfData%nUniquePoints_total
!  write(IBM_logUnit(3),*)IBMData%vel_xk( (iProc-1)*3+1:(iProc-1)*3+3 )
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

        write(IBM_logUnit(3),*) '      Done interpolating the velocity!'

        write(IBM_logUnit(3),*) '      Correct the force terms at the '//      &
          &                     'surface points ...'

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'force_Xk before correct'
!do iProc = 1, me(iIBM)%surfData%nUniquePoints_total
!  write(IBM_logUnit(3),*)IBMData%force_Xk( (iProc-1)*3+1:(iProc-1)*3+3 )
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

!! !$omp end single
        ! correct the forces
        call mus_corrForce_Xk( IBMData      = me(iIBM)%IBMData,                &
          &                    convFac      = convFac,                         &
          &                    nPoints      = me(iIBM)%surfData%               &
          &                                               nUniquePoints_total, &
          &                    parentIDs    = me(iIBM)%surfData%               &
          &                                            parentIDs(iLevel)%ptrs, &
          &                   nElems_fluid = levelDesc%elem%nElems( eT_fluid ) )
!! !$omp single

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'force_Xk after correct'
!do iProc = 1, me(iIBM)%surfData%nUniquePoints_total
!  write(IBM_logUnit(3),*)IBMData%force_Xk( (iProc-1)*3+1:(iProc-1)*3+3 )
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

        write(IBM_logUnit(3),*) '      Done correcting the forces!'

        if( writeForce )  then
          do iProc = 1, me(iIBM)%surfData%nUniquePoints_total
            if( sqrt( sum( me(iIBM)%IBMData%vel_Xk( (iProc-1)*3                &
              &                       +1:(iProc-1)*3+3 )**2 )) > maxVel_Xk)then
              maxVel_Xk = sqrt( sum( me(iIBM)%IBMData%vel_Xk( (iProc-1)*3      &
                &                                      +1:(iProc-1)*3+3 )**2 ))
            end if
            if( sqrt( sum( me(iIBM)%IBMData%force_xk( (iProc-1)*3              &
              &                     +1:(iProc-1)*3+3 )**2 )) > maxforce_xk)then
              maxforce_xk = sqrt( sum( me(iIBM)%IBMData%force_xk( (iProc-1)*3  &
                &                                      +1:(iProc-1)*3+3 )**2 ))
            end if
          end do

!          do iProc = 1, IBMData%neighs_Xk%nVals
          do iProc = 1, me( iIBM )%neighs_Xk%nVals
            if( sqrt( sum( me(iIBM)%IBMData%vel_X( (iProc-1)*3                 &
              &                        +1:(iProc-1)*3+3 )**2 )) > maxVel_X)then
              maxVel_X = sqrt( sum( me(iIBM)%IBMData%vel_X( (iProc-1)*3        &
                &                                      +1:(iProc-1)*3+3 )**2 ))
            end if
            if( sqrt( sum( me(iIBM)%IBMData%vel_X_ini( (iProc-1)*3             &
              &                    +1:(iProc-1)*3+3 )**2 )) > maxVel_X_ini)then
              maxVel_X_ini = sqrt(sum( me(iIBM)%IBMData%vel_X_ini( (iProc-1)*3 &
                &                                      +1:(iProc-1)*3+3 )**2 ))
            end if
            if( sqrt( sum( me(iIBM)%IBMData%force_X( (iProc-1)*3               &
              &                        +1:(iProc-1)*3+3 )**2 )) > maxForce)then
              maxForce = sqrt( sum( me(iIBM)%IBMData%force_X( (iProc-1)*3      &
                &                                      +1:(iProc-1)*3+3 )**2 ))
            end if
          end do
          write(forceUnit,*) iIter, maxForce, maxVel_X, maxVel_X_ini,          &
            &                maxForce_xk, abs(0._rk - maxVel_Xk)
        end if
!! !$omp end single
      end do ! iIter

!! !$omp single

      ! stop the IBMloop timer
      call tem_stopTimer( me          = me( iIBM )%timings%timedat, &
        &                 timerHandle = me( iIBM )%timerHandles(4)  )

      ! start the IBMfinish timer
      call tem_startTimer( me          = me( iIBM )%timings%timedat, &
        &                  timerHandle = me( iIBM )%timerHandles(5)  )

      ! -----------------------------------------------------------------------!
      !      D O N E  W I T H   T H E  I N A M U R O  I T E R A T I O N S      !
      ! -----------------------------------------------------------------------!

      write(IBM_logUnit(3),*) '   Done with Inamuro iterations!'

      ! stop the dyn load balancing timer
      call tem_stopTimer( timerHandle = mus_timerHandles%doIBM(iLevel) )

      write(IBM_logUnit(3),*) '      Communicating the force values ...'
      ! ... communicate the force values for Xk
      call commPattern%exchange_real(                                          &
        &                    send         = me(iIBM)%IBMData%IBMSend_Xk,       &
        &                    recv         = me(iIBM)%IBMData%IBMRecv_Xk,       &
        &                    state        = me(iIBM)%IBMData%force_Xk,         &
        &                    message_flag = mess_forceXk,                      &
        &                    comm         = general%proc%comm )

      ! start the dyn load balancing timer
      call tem_startTimer( timerHandle = mus_timerHandles%doIBM(iLevel) )


      write(IBM_logUnit(3),*) '      Done communicating the force values!'

      write(IBM_logUnit(3),*) '   Calculate the new forces at the grid '//     &
        &                     'points ...'

!! !$omp end single

      ! calculate the new forces at the grid points
      call mus_calcForce_X( me           = me(iIBM),                           &
        &                   IBMData      = me(iIBM)%IBMData,                   &
        &                   nElems_fluid = levelDesc%elem%nElems( eT_fluid ),  &
        &                   convFac      = convFac )

!! !$omp single


      write(IBM_logUnit(3),*) '   Done calculating the new forces!'
      write(IBM_logUnit(3),*) '   Apply the new forces on the grid points ...'

!! !$omp end single

      ! apply the force on the grid points
      call mus_applyForce_X( me           = me( iIBM ),                        &
        &                    state        = state(:, pdf%nNext),           &
        &                    IBMData      = me(iIBM)%IBMData,                  &
        &                    layout       = layout,                            &
        &                    varPos       = globSys%method%val(                &
        &                                     stateVarMap(iField))             &
        &                                       %state_varPos,                 &
        &                    nScalars     = globSys%nScalars,                  &
        &                    nElems_fluid = levelDesc%elem%nElems( eT_fluid ), &
        &                    convFac      = convFac,                           &
        &                    time         = general%simControl%now,     &
        &                    levelDesc    = levelDesc                          )

!! !$omp single

      write(IBM_logUnit(3),*) '   Done applying the new forces!'

      ! stop the dyn load balancing timer
      call tem_stopTimer( timerHandle = mus_timerHandles%doIBM(iLevel) )

      ! communicate the new pdfs on the halo elements
      call commPattern%exchange_real(                                          &
        &                    send         = me(iIBM)%IBMData%IBMSend_X_pdf,    &
        &                    recv         = me(iIBM)%IBMData%IBMRecv_X_pdf,    &
        &                    state        = state(:, pdf%nNext),           &
        &                    message_flag = iLevel,                            &
        &                    comm         = general%proc%comm )

      ! start the dyn load balancing timer
      call tem_startTimer( timerHandle = mus_timerHandles%doIBM(iLevel) )

      write(IBM_logUnit(3),*)'Done processing IBM ', iIBM

      ! free all temporary variables
      call mus_free_IBMData( me           = me(iIBM)%IBMData,                  &
        &                    commPattern  = commPattern,                       &
        &                    nElems_fluid = levelDesc%elem%nElems( eT_fluid ) )

      ! stop the IBMfinish timer
      call tem_stopTimer( me          = me( iIBM )%timings%timedat, &
        &                 timerHandle = me( iIBM )%timerHandles(5)  )

!! !$omp end single

    end do ! iIBM

!! !$omp single

    if( writeForce )  then
      close(unit = forceunit)
    end if

    ! stop the dyn load balancing timer
    call tem_stopTimer( timerHandle = mus_timerHandles%doIBM(iLevel) )

!! !$omp end single

  end subroutine mus_inamuro_IBM
! **************************************************************************** !


! **************************************************************************** !
  !> This routine builds the neighbor lists for Xk -> x (neigh_Xk) and
  !! x->Xk (neigh_x) as well as the send and receive buffers for the eulerian
  !! elements.
  !!
  subroutine mus_fillNeigh_surfData( me, IBMData, stencil, levelDesc, globTree,&
    &                                parentIDs, commPattern, params )
    ! --------------------------------------------------------------------------
    !> IBM data including the surface data
    type( mus_IBM_type ), intent(inout) :: me
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> array of stencils (1 is the fluid stencil)
    type( tem_stencilHeader_type ), intent(in) :: stencil(:)
    !> the level descriptor incl. ghost and halo elements on the level iLevel
    type( tem_levelDesc_type ), intent(inout) :: levelDesc
    !> global tree information
    type( treelmesh_type ), intent(in) :: globTree
    !> array of parentID positions hosting the lagrangian points
    integer, intent(in) :: parentIDs(:)
    !> communication pattern to be used
    type( tem_commPattern_type ), intent(inout) :: commPattern
    !> global parameters
    type(mus_param_type), intent(inout)    :: params
    ! --------------------------------------------------------------------------
    ! tmp element position
    integer :: treeIDPos
    ! integer coordinates of the parent element
    integer :: parCoord(4)
    ! temporary integer coordinates to search for
    integer :: neighCoord(4)
    ! treeID of the neighbor element
    integer( kind=long_k) :: neighID
    ! counters
    integer :: iPoint
    integer :: iQQ
    integer :: iType
    integer :: iProc
    integer :: iVal
    ! tmp variable for the position in the send buffer
    integer :: send_pos
    ! tmp variable for the position in the recv buffer
    integer :: recv_pos
    ! tmp variable for the position in the dynamic neighbor array of the Xk
    integer :: neighPos
    ! tmp variable if the position was added to the dynamic neighbor array of
    ! the Xk
    logical :: wasAdded
    ! tmp aray of logicals for marking elements that have been added for
    ! communication for each participating process
    ! size: nElems, nProcsSend
    logical, allocatable :: toComm(:,:)
    logical :: match
    integer :: firstHaloPos
    ! has the parent ID of this Xk been processed before?
    logical, allocatable :: parentUsed(:)
    ! counter for the number of processes attached to this process
    integer :: nProcs_send
    ! adjacent processes to the current process
    integer, allocatable :: adjProcs_send(:)
    ! counter for the number of processes attached to this process
    integer :: nProcs_recv
    ! adjacent processes to the current process
    integer, allocatable :: adjProcs_recv(:)
    ! tmp array for storing the number of Xk to be received by this process
    integer, allocatable :: nElems_recv(:)
    ! does the process participate in the communication?
    logical, allocatable :: procPart(:)
    ! --------------------------------------------------------------------------

    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(3)  )

    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(10) )

    allocate( procPart( levelDesc%sendBuffer%nProcs ))
    procPart = .false.
    allocate( IBMData%map2globSend_pdf( levelDesc%sendBuffer%nProcs ))
    IBMData%map2globSend_pdf = 0
    allocate( IBMData%map2globRecv_pdf( levelDesc%recvBuffer%nProcs ))
    IBMData%map2globRecv_pdf = 0

    nProcs_send = 0
    nProcs_recv = 0

    allocate( adjProcs_send( levelDesc%sendBuffer%nProcs ))
    allocate( IBMData%treeIDs( levelDesc%sendBuffer%nProcs ))

    ! allocate the array for marking elements for sending and receiving
    allocate( toComm( levelDesc%nElems, params%general%proc%comm_size ))
    ! initialize toSend array to false
    toComm = .false.

    allocate( parentUsed( levelDesc%nElems ))
    parentUsed = .true.

    firstHaloPos = levelDesc%nElems-levelDesc%elem%nElems( eT_halo ) + 1

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
    ! initialize the wasAdded with false
    wasAdded = .false.
    ! initialize the array of pointers to the dynamic array of neighbors of Xk
    IBMData%ptrs_neighs_Xk = 0

    ! loop over the surface coordinates and ...
    do iPoint = 1, me%surfData%nUniquePoints_total
      ! ... check if the parent ID of the element is actually existing on
      !     this process
      if( parentIDs( iPoint ) > 0 )then
        if( parentUsed( parentIDs( iPoint )))then
          ! 1. get the integer coordinates of the current parentID
          parCoord = tem_CoordOfId( levelDesc%total( parentIDs( iPoint )))
          ! 2. loop over all possible neighbors and ...
          do iQQ = 1, stencil(me%stencilPos)%QQ
            ! ... update the neighbor coordinates
            neighCoord(1:3) = parCoord(1:3) - stencil(me%stencilPos) &
              &                                 %cxDir(:, iQQ)
            neighCoord(4)   = parCoord(4)
            ! ... get the neighboring treeID
            neighID = tem_IdOfCoord( neighCoord )
            ! ... search for the treeID in the total list
            eT_loop: do iType = eT_minRelevant, eT_maxRelevant
              treeIDPos = tem_treeIDinTotal( neighID, levelDesc, iType )
              if( treeIDPos > 0 )then
                ! ... and store the resulting postition in the dynamic array of
                !     neighbors
                call append( me       = me%neighs_Xk,                          &
                  &          val      = treeIDPos,                             &
                  &          pos      = neighPos,                              &
                  &          wasAdded = wasAdded )
  !call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
  !write(IBM_logUnit(3),*)'added neigh: ', wasAdded
  !write(IBM_logUnit(3),*)'iPoint:      ', iPoint
  !write(IBM_logUnit(3),*)'iQQ:         ', iQQ
  !write(IBM_logUnit(3),*)'neighPos:    ', treeIDpos
  !write(IBM_logUnit(3),*)'neighCoord:  ', neighCoord
  !write(IBM_logUnit(3),*)'neighID:     ', neighID
  !write(IBM_logUnit(3),*)'parentID:    ', levelDesc%total( parentIDs( iPoint ))
  !write(IBM_logUnit(3),*)'parentPos:   ', parentIDs( iPoint )
  !write(IBM_logUnit(3),*)'parentCoord: ', parCoord
  !call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
  !flush( IBM_logUnit(3) )
                ! ... and append the added position to the pointer array for
                !     each Xk
                IBMData%ptrs_neighs_Xk(                                        &
                  &       (iPoint-1)*stencil(me%stencilPos)%QQ+iQQ ) = neighPos
                ! ... if the position is an actual fluid element on this
                !     partition ...
                if( treeIDpos <= levelDesc%elem%nElems( eT_fluid ) )then
                  ! ... store the point position for the neighbor elements
                  call append( me  = IBMData%neighs_X( treeIDPos ),            &
                    &          val = iPoint )
                  if( wasAdded )then
                    me%locNeighs_Xk = me%locNeighs_Xk + 1
                  end if
                  ! ... check if the parentID of the Xk was a halo and the
                  !     neighbor was not added to the communication before ...
                  if( parentIDs( iPoint ) >= firstHaloPos )then
                    ! ... set the position of the element in the velocity array
                    !     which will be send
                    send_pos = (neighPos-1)*3

    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(11) )
                    ! ... fill the elemPos array in the send buffer
                    !     pass the map2globRecv here since it was build
                    !     depending on Xk which is the inverse of the X
                    !     for sending
                    call mus_IBM_fillSendPos_X(                                &
                      &               IBMData      = IBMData,                  &
                      &               globSend     = levelDesc%sendBuffer,     &
                      &               treeIDPos    = treeIDPos,                &
                      &               startPos     = send_pos,                 &
                      &               parentID     = levelDesc%total(          &
                      &                                parentIDs( iPoint )),   &
                      &               globTree     = globTree,                 &
                      &               added        = toComm( treeIDPos,: ),    &
                      &               match        = match )
    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(11) )
  !if( match )then
  !write(IBM_logUnit(3),*)'added vel for send: '
  !write(IBM_logUnit(3),*)'iPoint:      ', iPoint
  !write(IBM_logUnit(3),*)'iQQ:         ', iQQ
  !write(IBM_logUnit(3),*)'neighPos:    ', treeIDpos
  !write(IBM_logUnit(3),*)'treeID:      ', levelDesc%total(treeIDPos)
  !write(IBM_logUnit(3),*)'parent:      ', levelDesc%total(parentIDs( iPoint ))
  !write(IBM_logUnit(3),*)'parentPos:   ', parentIDs( iPoint )
  !write(IBM_logUnit(3),*)'send_pos:    ', send_pos
  !write(IBM_logUnit(3),*)'treeIDPos:   ', treeIDPos
  !write(IBM_logUnit(3),*)'toComm:      ', toComm( treeIDPos,: )
  !call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
  !end if
  !flush( IBM_logUnit(3) )
                  end if
                ! ... else if the neighbor is a halo element ...
                else if( treeIDPos >= firstHaloPos)then
                  ! ... check if the parentID of the Xk was a fluid which equals
                  !     the neighbor was added to the dynamic array of neighbor
                  !     positions of the Xk ...
                  if ( parentIDs(iPoint)                  &
                    &  <= levelDesc%elem%nElems(eT_fluid) ) then
                    ! ... set the position of the element in the velocity array
                    !     which will be receive
                    recv_pos = (neighPos-1)*3
    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(11) )
                    ! ... fill the elemPos array in the receive buffer
                    !     pass the map2globSend here since it was build
                    !     depending on Xk which is the inverse of the X
                    !     for receiving
                    call mus_IBM_fillRecvPos_X(                             &
                      &              IBMData       = IBMData,               &
                      &              globRecv      = levelDesc%recvBuffer,  &
                      &              treeIDPos     = treeIDPos,             &
                      &              startPos      = recv_pos,              &
                      &              added         = toComm( treeIDPos,: ), &
                      &              match         = match                  )
    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(11) )
  !if( match )then
  !write(IBM_logUnit(3),*)'added vel for recv: '
  !write(IBM_logUnit(3),*)'iPoint:      ', iPoint
  !write(IBM_logUnit(3),*)'iQQ:         ', iQQ
  !write(IBM_logUnit(3),*)'neighPos:    ', treeIDpos
  !write(IBM_logUnit(3),*)'treeID:      ', levelDesc%total(treeIDPos)
  !write(IBM_logUnit(3),*)'parent:      ', levelDesc%total(parentIDs( iPoint ))
  !write(IBM_logUnit(3),*)'parentPos:   ', parentIDs( iPoint )
  !write(IBM_logUnit(3),*)'recv_pos:    ', recv_pos
  !write(IBM_logUnit(3),*)'treeIDPos:   ', treeIDPos
  !write(IBM_logUnit(3),*)'toComm:      ', toComm( treeIDPos,: )
  !call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
  !end if
  !flush( IBM_logUnit(3) )
                  end if
                end if
                ! in case the sendHalo bit is set this neighbor has to be
                ! communicated ...
                if( btest( levelDesc%property( treeIDPos ), prp_sendHalo ))then
                  ! ... loop over the global number of processes and ...
                  do iProc = 1, levelDesc%sendBuffer%nProcs
                    ! ... search in the global sendBuffer%elemPos array for the
                    !     process to send to
                    do iVal = 1, levelDesc%sendBuffer%elemPos( iProc )%nVals
                      ! ... if the position in the sendBuffer%elemPos equals to
                      !     the position in the parentID ...
                      if( levelDesc%sendBuffer%elemPos(iProc)%val(iVal) ==     &
                        &                                      treeIDPos .and. &
                        & wasAdded )then
                        ! ... and append the corresponding treeIDs to the
                        !     growing array
                        call append( me  = IBMData%treeIDs( iProc ),           &
                          &          val = treeIDPos )
                        ! this process participates in the communication
                        procPart( iProc ) = .true.
                      end if
                    end do
                  end do
                end if
                exit eT_loop
              else
                IBMData%ptrs_neighs_Xk(                                        &
                  &            (iPoint-1)*stencil(me%stencilPos)%QQ + iQQ) = 0
              end if ! treeIDPos > 0
            end do eT_loop
          end do ! iQQ
          parentUsed( parentIDs( iPoint )) = .false.
        end if ! parentUsed??
      end if ! parentID > 0
    end do ! iPoint

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!do iProc = 1, levelDesc%sendBuffer%nProcs
!  write(IBM_logUnit(3),*)'treeIDs of iProc ', iProc
!  if( treeIDs(iProc)%nVals > 0 )then
!    write(IBM_logUnit(3),*) treeIDs(iProc)%val(1:treeIDs(iProc)%nVals)
!  end if
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

    ! --------------------------------------------------------------------------
    ! Now the neighbor array and the buffers for the velocity on the eulerian
    ! elements X are available. Now we need to build the buffers for the
    ! pdf on the eulerian elements
    !
    ! IMPORTANT: The buffer for the velocity on the eulerian elements DOES NOT
    !            equal the one for the pdfs !!!!!

    ! fill the helper arrays for building the subcommunicators
    do iProc = 1, levelDesc%sendBuffer%nProcs
      ! if this process participates in the communication ...
      if( procPart( iProc ))then
        ! ... increase the process counter by 1
        nProcs_send = nProcs_send + 1
        ! ... add the corresponding globsend process to the array
        adjProcs_send( nProcs_send ) = levelDesc%sendBuffer%proc( iProc )
        ! ... add the process position in the global buffer
        !     to map2globSend
        IBMData%map2globSend_pdf( nProcs_send ) = iProc
      end if
    end do

    allocate( nElems_recv( levelDesc%recvBuffer%nProcs ))
    nElems_recv = 0
    ! initialize the global integer buffer
    do iProc = 1, levelDesc%sendBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%sendBuffer%buf_int( iProc ))
      call commPattern%initBuf_int( me    = levelDesc%sendBuffer%              &
        &                                                    buf_int( iProc ), &
        &                           pos   = (/iProc/),                         &
        &                           nVals = 1 )
    end do

    do iProc = 1, levelDesc%recvBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%recvBuffer%buf_int( iProc ))
      call commPattern%initBuf_int( me    = levelDesc%recvBuffer%              &
        &                                                    buf_int( iProc ), &
        &                           pos   = (/iProc/),                         &
        &                           nVals = 1 )
    end do

    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(3)  )

    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(10) )

    ! ... exchange the amount of Xk to be communicated
    call commPattern%exchange_int( send         = levelDesc%sendBuffer,        &
      &                            recv         = levelDesc%recvBuffer,        &
      &                            state        = nElems_recv,                 &
      &                            message_flag = mess_amountXk,               &
      &                            send_state   = IBMData%treeIDs( 1:levelDesc%&
      &                                              sendBuffer%nProcs )%nVals,&
      &                            comm         = params%general%proc%comm )

    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(3)  )

    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(10) )

    ! ... free the buffers
    do iProc = 1, levelDesc%sendBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%sendBuffer%buf_int( iProc ))
    end do
    do iProc = 1, levelDesc%recvBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%recvBuffer%buf_int( iProc ))
    end do

    ! set the number of procs
    IBMData%IBMSend_X_pdf%nProcs = nProcs_send
    ! allocate the array of processes to send to and set them
    allocate( IBMData%IBMSend_X_pdf%proc( nProcs_send ))
    IBMData%IBMSend_X_pdf%proc = adjProcs_send(1:nProcs_send)
    allocate( IBMData%IBMSend_X_pdf%buf_long( nProcs_send ))
    ! allocate the array of nElemsProc with the number of processes
    allocate( IBMData%IBMSend_X_pdf%nElemsProc( nProcs_send ))
    allocate( IBMData%IBMSend_X_pdf%elemPos( nProcs_send ))

    allocate( adjProcs_recv( levelDesc%recvBuffer%nProcs ))

    ! ... IBMRecv_X_pdf buffers
    nProcs_recv = 0
    ! loop over the procs in the global receive buffer ...
    do iProc = 1, levelDesc%recvBuffer%nProcs
      ! ... if elements will be received by this proc ...
      if( nElems_recv( iProc ) > 0 )then
        ! ... increase the recv counter by 1
        nProcs_recv = nProcs_recv + 1
        ! add the process to the adjacent recv array
        adjProcs_recv( nProcs_recv ) = levelDesc%recvBuffer%proc( iProc )
        ! add this process to the map
        IBMData%map2globRecv_pdf( nProcs_recv ) = iProc
      end if
    end do

    ! set the number of procs
    IBMData%IBMRecv_X_pdf%nProcs = nProcs_recv
    ! allocate the array of processes to recv to and set them
    allocate( IBMData%IBMRecv_X_pdf%proc( nProcs_recv ))
    IBMData%IBMRecv_X_pdf%proc = adjProcs_recv(1:nProcs_recv)
    allocate( IBMData%IBMRecv_X_pdf%buf_long( nProcs_recv ))
    ! allocate the array of nElemsProc with the number of processes
    allocate( IBMData%IBMRecv_X_pdf%nElemsProc( nProcs_recv ))
    allocate( IBMData%IBMRecv_X_pdf%elemPos( nProcs_recv ))

    ! loop over all receive ranks and ...
    do iProc = 1, nProcs_recv
      ! ... copy the number of elements to be received for each rank
      IBMData%IBMRecv_X_pdf%nElemsProc( iProc ) =                              &
        &                        nElems_recv( IBMData%map2globRecv_pdf( iProc ))
    end do

    deallocate( nElems_recv )
    deallocate( adjProcs_send )
    deallocate( adjProcs_recv )
    deallocate( procPart )

    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(3)  )

    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(10) )
  end subroutine mus_fillNeigh_surfData
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine fills the initial force, initial velocity, predef. velocity
  !! array for the surface points Xk as well as the velocity array for the
  !! neighbors.
  !!
  subroutine mus_inamuroIni( me, IBMData, state, globTree, stencil, stencilPos,&
    &                        levelDesc, globSys, pdfPos, iLevel, convFac )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(in) :: me
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> the state array holding the pdfs
    real(kind=rk), intent(in) :: state(:)
    !> global tree information
    type( treelmesh_type ), intent(in) :: globTree
    !> stencil used by the IBM
    type( tem_stencilHeader_type ), intent(in) :: stencil(:)
    !> position of the IBM stencil
    integer, intent(in) :: stencilPos
    !> level descriptor incl. ghost and fluid elements
    type( tem_levelDesc_type ), intent(in) :: levelDesc
    !> global variable system of the current scheme
    type(tem_varSys_type) :: globSys
    !> position of the velocity in the global variable system
    integer :: pdfPos
    !> the current level
    integer :: iLevel
    !> conversion factors
    type(mus_convertFac_type), intent(in) :: convFac
    ! --------------------------------------------------------------------------
    ! counters
    integer :: iPoint
    integer :: iNeigh
    integer :: iDir
    ! tmp variable to store the point coordinates to
    real(kind=rk) :: pos(1,3)
    ! tmp variable to store the velocity
    real(kind=rk) :: bary(3)
    ! neighboring position
    integer :: neighPos
    ! neighboring position of the Xk to access the unique neighbor array
    integer :: neighPos_X
    ! min and max position in the pointCoords array
    integer :: minPos, maxPos
    ! min and max position for the neighbor
    integer :: minNeigh, maxNeigh
    ! --------------------------------------------------------------------------

    ! reset the arrays to 0
    IBMData%inaDelta_Xk = 0._rk

    ! initialize the array with 0
    IBMData%vel_X_ini = 0._rk

    ! initialize the array with 0
    IBMData%vel_Xk_ini = 0._rk

    ! initialize the array with 0
    IBMData%force_Xk = 0._rk

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'InaDelta_Xk'
!write(IBM_logUnit(3),*)"dx: ", convFac%length
!write(IBM_logUnit(3),*)"dt: ", convFac%time

!$omp parallel
    ! 1. fill the array of inamura delta functions for the Xk
    !
    ! loop over the surface points and ...
!$omp do private( minPos, maxPos, pos, bary, neighPos_X )
    do iPoint = 1, me%surfData%nUniquePoints_total
      ! ... check if the parent element is a fluid element on this process
      if ( me%surfData%parentIDs(iLevel)%ptrs(iPoint) > 0     &
         & .and. ( me%surfData%parentIDs(iLevel)%ptrs(iPoint) &
         &         <= levelDesc%elem%nElems(eT_fluid) )       ) then
        ! ... set the min and max pos
        minPos = (iPoint-1)*3+1
        maxPos = (iPoint-1)*3+3
        ! ... get the x,y,z coordinates
        pos(1,1:3) = me%surfData%pointCoords(minPos:maxPos)
        ! ... loop over the neighbors and ...
        do iDir = 1, stencil(stencilPos)%QQ
          neighPos_X = IBMData%ptrs_neighs_Xk(                                 &
            &                        (iPoint-1)*stencil(stencilPos)%QQ + iDir )
          if( neighPos_X > 0 )then
            ! ... calculate the barycenter
            bary = tem_baryOfId( globTree,                                     &
              &              levelDesc%total(me%neighs_Xk%val(neighPos_X)))
            ! ... and the delta function
            IBMData%inaDelta_Xk((iPoint-1)*stencil(stencilPos)%QQ + iDir) =    &
              &                             inamuroDelta3D( pos(1,1:3)-bary,   &
              &                                             convFac%length)
          end if
        end do
      end if ! parentID is fluid
    end do ! iPoint

!$omp end do

!$omp end parallel

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

!! !$omp single

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'InaDelta_X'

    ! 2. fill the array of inamura delta functions for the X
    !
    ! loop over the elements and ...
    do iNeigh = 1, levelDesc%elem%nElems( eT_fluid )
      do iPoint = 1, IBMData%neighs_X(iNeigh)%nVals
        ! ... set the min and max pos
        minPos = (IBMData%neighs_X(iNeigh)%val(iPoint)-1)*3+1
        maxPos = (IBMData%neighs_X(iNeigh)%val(iPoint)-1)*3+3
        pos(1,1:3)=me%surfData%pointCoords(minPos:maxPos)
        bary = tem_baryOfId( globTree, levelDesc%total( iNeigh ))
        call append( me  = IBMData%inaDelta_X( iNeigh ),                       &
          &          val = inamuroDelta3D( pos(1,1:3)-bary, convFac%length ))
      end do
    end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

!! !$omp end single



    ! 3. fill the initial velocity array for the neighbors
    !    (Inamuro paper: u*(X))
    !
    ! loop over the unique neighbors and ...
!$omp parallel
!    do iNeigh = 1, IBMData%neighs_Xk%nVals
!$omp do private( minPos, maxPos )
    do iNeigh = 1, me%neighs_Xk%nVals
      ! ... set min and max position
      minPos = (iNeigh - 1)*3 + 1
      maxPos = (iNeigh - 1)*3 + 3
      ! ... calculate the velocity in LB units
      IBMData%vel_X_ini( minPos:maxPos ) = getVelocity(                        &
        &                  state    = state,                                   &
        &                  elem     = me%neighs_Xk%val( iNeigh ),              &
        &                  stencil  = stencil(1),                              &
        &                  varPos   = globSys%method%val(pdfPos)%state_varPos, &
        &                  nScalars = globSys%nScalars )
      ! rescale the velocity from LB to physical units
      IBMData%vel_X_ini( minPos:maxPos ) =                                   &
        &                    IBMData%vel_X_ini( minPos:maxPos ) * convFac%vel
    end do
!$omp end do


    ! 4. fill the velocity array for the interpolated velocity
    !    (Inamuro paper: u*(X_k))
    !
    ! loop over the surface points and ...
!$omp do private( minPos, maxPos, neighPos, neighPos_X, minNeigh, maxNeigh, iDir )
    do iPoint = 1, me%surfData%nUniquePoints_total
      ! ... check if the parent element is a fluid element on this process
      if ( me%surfData%parentIDs(iLevel)%ptrs(iPoint) > 0     &
        &  .and. ( me%surfData%parentIDs(iLevel)%ptrs(iPoint) &
        &          <= levelDesc%elem%nElems(eT_fluid) )       ) then
        ! ... set the min and max pos
        minPos = (iPoint-1)*3+1
        maxPos = (iPoint-1)*3+3
        ! ... loop over the neighbors and ...
        do iDir = 1, stencil(stencilPos)%QQ
          ! ... set min and max neighbor pos in vel_X_ini array and neighpos
          !     in inaDelta_Xk array
          neighPos = (iPoint-1)*stencil(stencilPos)%QQ + iDir
          neighPos_X = IBMData%ptrs_neighs_Xk( neighPos )
          if( neighPos_X > 0 )then
            minNeigh = (neighPos_X-1)*3+1
            maxNeigh = (neighPos_X-1)*3+3
            ! ... calculate the initial velocity
            IBMData%vel_Xk_ini( minPos:maxPos ) =                             &
              &                          IBMData%vel_Xk_ini( minPos:maxPos )  &
              &                        + IBMData%vel_X_ini(minNeigh:maxNeigh) &
              &                        * IBMData%inaDelta_Xk( neighPos )      &
              &                        * convFac%length**3
          end if
        end do
      end if ! parentID is fluid
    end do ! iPoint

!$omp end do

!write(IBM_logUnit(3),*)'vel_Xk_ini in inamuroIni:'
!do iPoint = 1, me%surfData%nUniquePoints_total
!  minPos = (iPoint-1)*3+1
!  maxPos = (iPoint-1)*3+3
!  write(IBM_logUnit(3),*)IBMData%vel_Xk_ini( minPos:maxPos)
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

    ! 5. fill the force array for the Xk
    !    (Inamuro paper: g_0(X_k))
    !
    ! loop over the elements in the neighbor array and ...
!$omp do private( minPos, maxPos )
    do iPoint = 1, me%surfData%nUniquePoints_total
      ! ... check if the parent element is a fluid element on this process
      if ( me%surfData%parentIDs(iLevel)%ptrs(iPoint) > 0     &
         & .and. ( me%surfData%parentIDs(iLevel)%ptrs(iPoint) &
         &         <= levelDesc%elem%nElems(eT_fluid) )       ) then
        minPos = (iPoint-1)*3+1
        maxPos = (iPoint-1)*3+3
        ! ... calculate the force
        IBMData%force_Xk( minPos:maxPos ) =                                    &
          &                         ( IBMData%vel_Xk_surf( minPos:maxPos )     &
          &                         - IBMData%vel_Xk_ini( minPos:maxPos ))     &
          &                         / convFac%time
      end if ! parentID is fluid
    end do ! iPoint

!$omp end do


!$omp end parallel

  end subroutine mus_inamuroIni
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine fills the force array for the X (neighbors).
  !! (Inamuro paper: step 1, fill g_l(X))
  !!
  subroutine mus_calcForce_X( me, IBMData, nElems_fluid, convFac )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(in) :: me
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> number of fluid elements on this process
    integer, intent(in) :: nElems_fluid
    !> conversion factors
    type(mus_convertFac_type), intent(in) :: convFac
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iNeigh
    integer :: iPoint
    ! neighboring position
    integer :: neighPos
    ! min and max position in the pointCoords array
    integer :: minPos, maxPos
    ! min and max position in the force_X array
    integer :: minForce, maxForce
    ! --------------------------------------------------------------------------

!! !$omp single

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)"calcForce_X"
    ! reset the array to 0
    IBMData%force_X = 0._rk
!write(IBM_logUnit(3),*)"loopTo: ", me%neighs_Xk%nVals

!! !$omp end single

    ! loop over the force array and ...
!    do iNeigh = 1, IBMData%neighs_Xk%nVals

!$omp parallel

!$omp do private( neighPos, iPoint, minPos, maxPos, minForce, maxForce )

    do iNeigh = 1, me%neighs_Xk%nVals
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)"iNeigh: ", iNeigh
      ! ... get the position in the neigbor array of the X
!      neighPos = IBMData%neighs_Xk%val( iNeigh )
      neighPos = me%neighs_Xk%val( iNeigh )
      ! ... loop over all Xk relevant in the neighborhood and ...
      if( neighPos <= nElems_fluid )then
        do iPoint = 1, IBMData%neighs_X( neighPos )%nVals
!write(IBM_logUnit(3),*)"iPoint: ", iPoint
          if( IBMData%inaDelta_X( neighPos )%val( iPoint ) > 0._rk )then
            ! ... set the min and max pos
            minPos = (IBMData%neighs_X(neighPos)%val(iPoint)-1)*3+1
            maxPos = (IBMData%neighs_X(neighPos)%val(iPoint)-1)*3+3
            minForce = (iNeigh-1)*3+1
            maxForce = (iNeigh-1)*3+3
            ! ... calculate the force and multiply it with the area times dx
            !     by this the units of deltaFct and area*dx cancel out
            IBMData%force_X( minForce:maxForce ) =                             &
              &             IBMData%force_X( minForce:MaxForce )               &
              &           + IBMData%force_Xk( minPos:maxPos )                  &
              &           * IBMData%inaDelta_X( neighPos )%val(iPoint)         &
              &           * me%surfData%surfArea(IBMData%neighs_X(neighPos)%   &
              &                                                    val(iPoint))&
              &           * convFac%length
          end if
        end do
      end if
    end do

!$omp end do

!$omp end parallel

  end subroutine mus_calcForce_X
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine corrects the velocity values according to the force on X
  !! (neighbors).
  !! (Inamuro paper: step 2, correct u_l(X))
  !!
  subroutine mus_corrVel_X( me, IBMData, convFac )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(in) :: me
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> conversion factors
    type(mus_convertFac_type), intent(in) :: convFac
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iNeigh
    ! minimum and maximum position for the velocity at one neighbor
    integer :: minPos, maxPos
    ! --------------------------------------------------------------------------

!! !$omp single

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'corr vel_X'
    ! reset the array to 0
    IBMData%vel_X = 0._rk

!! !$omp end single

    ! loop over the velocity array and ...
!    do iNeigh = 1, IBMData%neighs_Xk%nVals

!$omp parallel

!$omp do private( minPos, maxPos )

    do iNeigh = 1, me%neighs_Xk%nVals
      ! ... set the min and max pos
      minPos = (iNeigh-1)*3+1
      maxPos = (iNeigh-1)*3+3
      ! ... calculate the velocities
      IBMData%vel_X( minPos:maxPos ) = IBMData%vel_X_ini( minPos:maxPos )      &
        &                            + IBMData%force_X( minPos:maxPos )        &
        &                            * convFac%time
    end do
!$omp end do

!$omp end parallel

  end subroutine mus_corrVel_X
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine interpolates the velocity values using the velocity on Xk.
  !! (neighbors).
  !! (Inamuro paper: step 3, correct u_l(X_k))
  !!
  subroutine mus_intpVel_Xk( IBMData, nNeighs, convFac, nPoints, parentIDs,    &
    &                        nElems_fluid )
    ! --------------------------------------------------------------------------
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> the number of  neighbors for the surface points
    integer, intent(in) :: nNeighs
    !> conversion factors
    type(mus_convertFac_type), intent(in) :: convFac
    !> number of points
    integer, intent(in) :: nPoints
    !> array of parentID positions hosting the lagrangian points
    integer, intent(in) :: parentIDs(:)
    !> number of fluid elements on this process
    integer, intent(in) :: nElems_fluid
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iPoint
    integer :: iNeigh
    ! min and max neighbor position
    integer :: minNeigh, maxNeigh
    ! min and max position in the pointCoords array
    integer :: minPos, maxPos
    ! neighboring position of the Xk
    integer :: neighPos_X
    ! --------------------------------------------------------------------------

!! !$omp single

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'intp vel_Xk:'
    ! reset the array to 0
    IBMData%vel_Xk = 0._rk

!! !$omp end single

    ! loop over the surface points and ...

!$omp parallel

!$omp do private( minPos, maxPos, neighPos_X, minNeigh, maxNeigh )

    do iPoint = 1, nPoints
      ! ... check if the parent is a fluid element on this process
      if( parentIDs( iPoint ) > 0 .and.                                     &
        & parentIDs( iPoint ) <= nElems_fluid )then
        ! ... set the min and max pos
        minPos = (iPoint-1)*3+1
        maxPos = (iPoint-1)*3+3
        ! ... loop over the neighbors of the Xk ...
        do iNeigh = 1, nNeighs
          neighPos_X = IBMData%ptrs_neighs_Xk( (iPoint-1)*nNeighs+iNeigh )
          if( neighPos_X > 0 )then
            minNeigh = (neighPos_X-1)*3+1
            maxNeigh = (neighPos_X-1)*3+3
            ! ... interpolate the velocity at the Xk
            IBMData%vel_Xk(minPos:maxPos)                &
              &  = IBMData%vel_Xk( minPos:maxPos )       &
              &    + IBMData%vel_X( minNeigh: maxNeigh ) &
              &    * IBMData%inaDelta_Xk(                &
              &           (iPoint-1)*nNeighs + iNeigh )  &
              &    * convFac%length**3
          end if
        end do
      end if
    end do

!$omp end do

!$omp end parallel

  end subroutine mus_intpVel_Xk
! **************************************************************************** !


! **************************************************************************** !
  !>
  !!
  subroutine mus_corrForce_Xk( IBMData, convFac, nPoints, parentIDs,           &
    &                          nElems_fluid )
    ! --------------------------------------------------------------------------
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> conversion factors
    type(mus_convertFac_type), intent(in) :: convFac
    !> number of points
    integer, intent(in) :: nPoints
    !> array of parentID positions hosting the lagrangian points
    integer, intent(in) :: parentIDs(:)
    !> number of fluid elements on this process
    integer, intent(in) :: nElems_fluid
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iPoint
    ! min and max position in the pointCoords array
    integer :: minPos, maxPos
    ! --------------------------------------------------------------------------

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'corr force_Xk:'
    ! loop over the surface points and ...

!$omp parallel

!$omp do private( minPos, maxPos )

    do iPoint = 1, nPoints
      ! ... check if the parent is a fluid element
      if( parentIDs( iPoint ) > 0 .and.                                        &
        & parentIDs( iPoint ) <= nElems_fluid )then
        ! ... set the min and max pos
        minPos = (iPoint-1)*3+1
        maxPos = (iPoint-1)*3+3
        ! ... correct the force term
        IBMData%force_Xk(minPos:maxPos) = IBMData%force_Xk(minPos:maxPos)      &
          &                             + (IBMData%vel_Xk_surf(minPos:maxPos)  &
          &                             -  IBMData%vel_Xk(minPos:maxPos))      &
          &                             / convFac%time
!write(IBM_logUnit(3),*)'result: ', IBMData%force_Xk(minPos:maxPos)
      end if
    end do

!$omp end do

!$omp end parallel

  end subroutine mus_corrForce_Xk
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine applies the force calculated to the eulerian elements.
  !!
  subroutine mus_applyForce_X( me, state, IBMData, layout, varPos, nScalars,   &
    &                          nElems_fluid, convFac, time, levelDesc )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(in) :: me
    !> the state array holding the pdfs
    real(kind=rk), intent(inout) :: state(:)
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> scheme layout of the current scheme incl. array of stencils and weights
    type( mus_scheme_layout_type ), intent(in) :: layout
    !> variable positions of the state variables in the levelDesc/state vector
    integer, intent(in) :: varPos(:)
    !> number of scalars in the global variable system
    integer, intent(in) :: nScalars
    !> number of fluid elements on this partition and level
    integer, intent(in) :: nElems_fluid
    !> conversion factors
    type(mus_convertFac_type), intent(in) :: convFac
    !> current time
    type(tem_time_type), intent(in) :: time
    !> level descriptor incl. ghost and fluid elements
    type( tem_levelDesc_type ), intent(in) :: levelDesc
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iNeigh, iDir
    ! tmp density
    real(kind=rk) :: tmp_dens
    ! element position in the state vector/levelDesc
    integer :: elem
    ! tmp variable storing force * c_i
    real(kind=rk) :: tmp_forceDir
    ! position of the force values for the elements
    integer :: forcePos
    ! number of elements inside state array
    integer :: nElems
    ! --------------------------------------------------------------------------
    nElems = size( state ) / nScalars
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'APPLYFORCE'
    ! loop over all neighbors of the Xk ...
!    do iNeigh = 1, IBMData%neighs_Xk%nVals

!$omp parallel

!$omp do private( elem, tmp_dens, forcePos, tmp_forceDir )

    do iNeigh = 1, me%neighs_Xk%nVals
      ! ... if the element is existing ...
!      if( IBMData%neighs_Xk%val( iNeigh ) <= nElems_fluid )then
      if( me%neighs_Xk%val( iNeigh ) <= nElems_fluid )then
!        elem = IBMData%neighs_Xk%val( iNeigh )
        elem = me%neighs_Xk%val( iNeigh )
        ! ... calculate the density at this element
        tmp_dens = getDensity( state    = state,           &
          &                    elem     = elem,            &
          &                    stencil  = layout%fStencil, &
          &                    varPos   = varPos,          &
          &                    nScalars = nScalars         )
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'iElem:         ', levelDesc%total( elem )
        ! ... and loop over the directions of the fluid stencil ...
        do iDir = 1, layout%fStencil%QQN
          forcePos = (iNeigh-1)*3
          ! .. calculate c_i*g and scale it to lattice units
          tmp_forceDir =                                                     &
            &    ( layout%fStencil%cxDir(1,iDir)*IBMData%force_X(forcePos+1) &
            &    + layout%fStencil%cxDir(2,iDir)*IBMData%force_X(forcePos+2) &
            &    + layout%fStencil%cxDir(3,iDir)*IBMData%force_X(forcePos+3))&
            &    / convFac%accel*3*layout%weight( iDir )*tmp_dens
          if( abs(tmp_forceDir) > 0.01 )then
            write(*,*)'Force is to strong!!!'
            write(*,*)'iter:   ', time%iter
            write(*,*)'iNeigh: ', iNeigh
            write(*,*)'treeID: ', levelDesc%total(me%neighs_Xk%val(iNeigh))
            write(*,*)'force:  ', tmp_forceDir
          end if

          ! ... apply the force on the fluid elements
          if( .not. btest( levelDesc%property( elem ), prp_solid ))then
            state( ( elem-1)* nscalars+ varpos(idir)) =     &
              &   state( ( elem-1)* nscalars+varpos(idir) ) &
              &   + tmp_forceDir
          end if
        end do
      end if
    end do

!$omp end do

!$omp end parallel

  end subroutine mus_applyForce_X
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine builds the communication types for the lagrangian elements
  !! Xk.
  !!
  subroutine mus_IBM_buildSendRecv_Xk( me, IBMData, levelDesc, commPattern, &
    &                                  globTree, iLevel, params             )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(inout) :: me
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> the level descriptor incl. the global send and receive buffers
    type( tem_levelDesc_type ), intent(inout) :: levelDesc
    !> communication pattern to be used
    type( tem_commPattern_type ), intent(inout) :: commPattern
    !> global tree information
    type( treelmesh_type ), intent(inout) :: globTree
    !> current level
    integer, intent(in) :: iLevel
    !> global parameters
    type(mus_param_type), intent(inout)    :: params
    ! --------------------------------------------------------------------------
    ! counters
    integer :: iProc, iVal, iPoint
    ! counter for the number of processes attached to this process
    integer :: nProcs_send
    ! adjacent processes to the current process
    integer, allocatable :: adjProcs_send(:)
    ! counter for the number of processes attached to this process
    integer :: nProcs_recv
    ! adjacent processes to the current process
    integer, allocatable :: adjProcs_recv(:)
    ! tmp array for storing the number of Xk to be received by this process
    integer, allocatable :: nElems_recv(:)
    ! does the process participate in the communication?
    logical, allocatable :: procPart(:)
    ! --------------------------------------------------------------------------

    ! in case the motion is not predefined communicate all Xk which just moved
    ! to a halo on the local proc
    if( .not. me%movPredef )then
      call mus_IBM_commNewPos( IBMData     = IBMData,                 &
        &                      levelDesc   = levelDesc,               &
        &                      commPattern = commPattern,             &
        &                      globTree    = globTree,                &
        &                      surfData    = me%surfData,             &
        &                      iLevel      = iLevel,                  &
        &                      comm        = params%general%proc%comm )
    end if

    ! start the IBMbuildXk timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(2)  )

    ! start the IBMbuildXk timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(8)  )

    ! --------------------------------------------------------------------------
    ! 1. fill the array of positions and helper arrays
    !
    allocate( IBMData%map2globSend( levelDesc%sendBuffer%nProcs ))
    IBMData%map2globSend = 0
    allocate( IBMData%map2globRecv( levelDesc%recvBuffer%nProcs ))
    IBMData%map2globRecv = 0
    allocate( procPart( levelDesc%sendBuffer%nProcs ))
    procPart = .false.

    nProcs_send = 0
    nProcs_recv = 0

    allocate( adjProcs_send( levelDesc%sendBuffer%nProcs ))
    allocate( IBMData%posXk( levelDesc%sendBuffer%nProcs ))

    ! loop over the Xk
    do iPoint = 1, me%surfData%nUniquePoints_total
      ! ... if the parent is a fluid element search in my send buffer
      !     for the corresponding positions
      if ( me%surfData%parentIDs(iLevel)%ptrs(iPoint) > 0     &
         & .and. ( me%surfData%parentIDs(iLevel)%ptrs(iPoint) &
         &         <= levelDesc%elem%nElems(eT_fluid) )       ) then
        ! ... loop over the global number of processes and ...
        do iProc = 1, levelDesc%sendBuffer%nProcs
          ! ... search in the global sendBuffer%elemPos array for the process to
          !     send to
          do iVal = 1, levelDesc%sendBuffer%elemPos( iProc )%nVals
            ! ... if the position in the sendBuffer%elemPos equals to the
            !     position in the parentID ...
            if( levelDesc%sendBuffer%elemPos(iProc)%val(iVal) ==               &
              &              me%surfData%parentIDs(iLevel)%ptrs( iPoint ))then
              ! ... and append them to the growing array
              call append( me  = IBMData%posXk( iProc ),                       &
                &          val = iPoint )
              ! this process participates in the communication
              procPart( iProc ) = .true.
            end if
          end do
        end do
      end if
    end do
    ! fill the helper arrays for building the subcommunicators
    do iProc = 1, levelDesc%sendBuffer%nProcs
      ! if this process participates in the communication ...
      if( procPart( iProc ))then
        ! ... increase the process counter by 1
        nProcs_send = nProcs_send + 1
        ! ... add the corresponding globsend process to the array
        adjProcs_send( nProcs_send ) = levelDesc%sendBuffer%proc( iProc )
        ! ... add the process position in the global buffer
        !     to map2globSend
        IBMData%map2globSend( nProcs_send ) = iProc
      end if
    end do

    ! --------------------------------------------------------------------------
    ! 2. send the number of elements their corresponding positions as well as
    !    actual new coordinates to the right processes
    !
    ! therefor use the global communicator and communicate the number of Xk
    ! that will be send
    !
    allocate( nElems_recv( levelDesc%recvBuffer%nProcs ))
    nElems_recv = 0
    ! initialize the global integer buffer
    do iProc = 1, levelDesc%sendBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%sendBuffer%buf_int( iProc ))
      call commPattern%initBuf_int( me    = levelDesc%sendBuffer%              &
        &                                                    buf_int( iProc ), &
        &                           pos   = (/iProc/),                         &
        &                           nVals = 1 )
    end do

    do iProc = 1, levelDesc%recvBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%recvBuffer%buf_int( iProc ))
      call commPattern%initBuf_int( me    = levelDesc%recvBuffer%              &
        &                                                    buf_int( iProc ), &
        &                           pos   = (/iProc/),                         &
        &                           nVals = 1 )
    end do

    ! stop the IBMbuildXk timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(2) )

    ! stop the IBMbuildXk timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(8)  )

    ! ... exchange the amount of Xk to be communicated
    call commPattern%exchange_int( send         = levelDesc%sendBuffer,        &
      &                            recv         = levelDesc%recvBuffer,        &
      &                            state        = nElems_recv,                 &
      &                            message_flag = mess_amountXk,               &
      &                            send_state   = IBMData%posXk( 1:levelDesc%  &
      &                                              sendBuffer%nProcs )%nVals,&
      &                            comm         = params%general%proc%comm )

    ! start the IBMbuildXk timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(2)  )

    ! start the IBMbuildXk timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(8)  )

    ! ... free the buffers
    do iProc = 1, levelDesc%sendBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%sendBuffer%buf_int( iProc ))
    end do
    do iProc = 1, levelDesc%recvBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%recvBuffer%buf_int( iProc ))
    end do


    ! --------------------------------------------------------------------------
    ! 3. Now each proc knows how many elements will be communicated.
    !    Initialize the IBMSend and IBMRecv buffers.
    !
    ! therefor build the IBMData%IBMSend_Xk and ...
    ! set the number of procs
    IBMData%IBMSend_Xk%nProcs = nProcs_send
    ! allocate the array of processes to send to and set them
    if( allocated( IBMData%IBMSend_Xk%proc ))                                  &
      & deallocate( IBMData%IBMSend_Xk%proc )
    allocate( IBMData%IBMSend_Xk%proc( nProcs_send ))
    IBMData%IBMSend_Xk%proc = adjProcs_send(1:nProcs_send)
    ! allocate the array of integer buffers
    if( allocated( IBMData%IBMSend_Xk%buf_int ))                               &
      & deallocate( IBMData%IBMSend_Xk%buf_int )
    allocate( IBMData%IBMSend_Xk%buf_int( nProcs_send ))
    ! allocate the array of nElemsProc with the number of processes
    if( allocated( IBMData%IBMSend_Xk%nElemsProc ))                            &
      & deallocate( IBMData%IBMSend_Xk%nElemsProc )
    allocate( IBMData%IBMSend_Xk%nElemsProc( nProcs_send ))
    if( allocated( IBMData%IBMSend_Xk%elemPos ))                               &
      & deallocate( IBMData%IBMSend_Xk%elemPos )
    allocate( IBMData%IBMSend_Xk%elemPos( nProcs_send ))

    allocate( adjProcs_recv( levelDesc%recvBuffer%nProcs ))

    ! ... IBMRecv_Xk buffers
    nProcs_recv = 0
    ! loop over the procs in the global receive buffer ...
    do iProc = 1, levelDesc%recvBuffer%nProcs
      ! ... if elements will be received by this proc ...
      if( nElems_recv( iProc ) > 0 )then
        ! ... increase the recv counter by 1
        nProcs_recv = nProcs_recv + 1
        ! add the process to the adjacent recv array
        adjProcs_recv( nProcs_recv ) = levelDesc%recvBuffer%proc( iProc )
        ! add this process to the map
        IBMData%map2globRecv( nProcs_recv ) = iProc
      end if
    end do
    ! set the number of procs
    IBMData%IBMRecv_Xk%nProcs = nProcs_recv
    ! allocate the array of processes to recv to and set them
    if( allocated( IBMData%IBMRecv_Xk%proc ))                                  &
      & deallocate( IBMData%IBMRecv_Xk%proc )
    allocate( IBMData%IBMRecv_Xk%proc( nProcs_recv ))
    IBMData%IBMRecv_Xk%proc = adjProcs_recv(1:nProcs_recv)
    ! allocate the array of integer buffers
    if( allocated( IBMData%IBMRecv_Xk%buf_int ))                               &
      & deallocate( IBMData%IBMRecv_Xk%buf_int )
    allocate( IBMData%IBMRecv_Xk%buf_int( nProcs_recv ))
    ! allocate the array of nElemsProc with the number of processes
    if( allocated( IBMData%IBMRecv_Xk%nElemsProc ))                            &
      & deallocate( IBMData%IBMRecv_Xk%nElemsProc )
    allocate( IBMData%IBMRecv_Xk%nElemsProc( nProcs_recv ))
    if( allocated( IBMData%IBMRecv_Xk%elemPos ))                               &
      & deallocate( IBMData%IBMRecv_Xk%elemPos )
    allocate( IBMData%IBMRecv_Xk%elemPos( nProcs_recv ))
    if( allocated( IBMData%IBMRecv_Xk%nElemsProc ))                            &
      & deallocate( IBMData%IBMRecv_Xk%nElemsProc )
    allocate( IBMData%IBMRecv_Xk%nElemsProc( nProcs_recv ))

    ! loop over all receive ranks and ...
    do iProc = 1, nProcs_recv
      ! ... copy the number of elements to be received for each rank
      IBMData%IBMRecv_Xk%nElemsProc( iProc ) =                                 &
        &                        nElems_recv( IBMData%map2globRecv( iProc ))
    end do

    deallocate( nElems_recv )
    deallocate( adjProcs_send )
    deallocate( adjProcs_recv )
    deallocate( procPart )

    ! stop the IBMbuildXk timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(2)  )

    ! stop the IBMbuildXk timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(8)  )

  end subroutine mus_IBM_buildSendRecv_Xk
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine prepares the send and receive buffers for the eulerian
  !! elements by copying information from the send and receive buffers for the
  !! lagrangian elements.
  !!
  subroutine mus_IBM_prepareSendRecv_X( IBMData )
    ! --------------------------------------------------------------------------
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    ! --------------------------------------------------------------------------
    ! counters
    integer :: iProc
    ! --------------------------------------------------------------------------

    ! Copy the information from the Xk buffers to the X buffers
    !
    ! 1. send buffer
    !
    ! set number of procs
    IBMData%IBMSend_X%nProcs = IBMData%IBMRecv_Xk%nProcs
    ! allocate and set array of procs
    allocate( IBMData%IBMSend_X%proc( IBMData%IBMSend_X%nProcs ))
    IBMData%IBMSend_X%proc = IBMData%IBMRecv_Xk%proc
    ! allocate array of integer buffers
    allocate( IBMData%IBMSend_X%buf_int( IBMData%IBMSend_X%nProcs ))
    ! allocate the array of nElemsProc and the array of element positions per
    ! process
    allocate( IBMData%IBMSend_X%nElemsProc( IBMData%IBMSend_X%nProcs ))
    allocate( IBMData%IBMSend_X%elemPos( IBMData%IBMSend_X%nProcs ))
    ! initialize the arrays for communicating the number of halos to be send
    do iProc = 1, IBMData%IBMSend_X%nProcs
      call init( IBMData%IBMSend_X%elemPos( iProc ), 1)
    end do
    !
    ! 2. receive buffer
    !
    ! set number of procs
    IBMData%IBMRecv_X%nProcs = IBMData%IBMSend_Xk%nProcs
    ! allocate and set array of procs
    allocate( IBMData%IBMRecv_X%proc( IBMData%IBMRecv_X%nProcs ))
    IBMData%IBMRecv_X%proc = IBMData%IBMSend_Xk%proc
    ! allocate array of integer buffers
    allocate( IBMData%IBMRecv_X%buf_int( IBMData%IBMRecv_X%nProcs ))
    ! allocate the array of nElemsProc and the array of element positions per
    ! process
    allocate( IBMData%IBMRecv_X%nElemsProc( IBMData%IBMRecv_X%nProcs ))
    allocate( IBMData%IBMRecv_X%elemPos( IBMData%IBMRecv_X%nProcs ))
    ! initialize the arrays for communicating the number of halos to be received
    do iProc = 1, IBMData%IBMRecv_X%nProcs
      call init( IBMData%IBMRecv_X%elemPos( iProc ), 1)
    end do
  end subroutine mus_IBM_prepareSendRecv_X
! **************************************************************************** !


! **************************************************************************** !
  !>
  !!
  !!
  subroutine mus_IBM_fillSendPos_X( IBMData, globSend, treeIDPos, startPos,    &
    &                               parentID, globTree, added, match )
    ! --------------------------------------------------------------------------
    !> IBM temporary datatype incl. map2glob and communicator for send and recv
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> global send communicator
    type( tem_communication_type ), intent(in) :: globSend
    !> element position in the level desc total list
    integer, intent(in) :: treeIDPos
    !> element position in the level desc total list
    integer(kind=long_k), intent(in) :: parentID
    !> global tree information
    type( treelmesh_type ), intent(in) :: globTree
    !> starting position of what to send as elemPos
    integer, intent(in) :: startPos
    !> is this element added to the communication
    logical, intent(inout) :: added(:)
    logical :: match
    ! --------------------------------------------------------------------------
    ! counters
    integer :: iProc
    integer :: iCoord
    ! tmp variable storing the process position from the map2glob
    integer :: procPos
    ! logical to determine wether the parent was found on a certain proc
    logical :: foundParent( size( added ))
    ! tmp variable for the absolut process
    integer :: proc
    integer :: nVals
    ! --------------------------------------------------------------------------

    foundParent  = .false.
    match = .false.

    globProc_loop: do iProc = 1, globTree%global%nParts
      if(tem_TreeIDComparison(parentID, globTree%Part_First(iProc)) > -1 .and. &
        &tem_TreeIDComparison(parentID, globTree%Part_Last(iProc)) < 1 )then
        foundParent(iProc) = .true.
        exit globProc_loop
      end if
    end do globProc_loop

    proc_loop: do iProc = 1, IBMData%IBMSend_X%nProcs
      procPos = IBMData%map2globRecv( iProc )
      ! set the absolut proc position (+1 due to mpi numbering)
      proc = globSend%proc( procPos )+1
      nVals = globSend%elemPos( procPos )%nVals
      if(any(globSend%elemPos( procPos )%val(1:nVals) == treeIDPos ) .and.     &
        & foundParent(proc) .and. .not. added(proc) )then
        ! ... append the position of the velocity coordinates
        !     in the velocity array to the element positions in
        !     the recv buffer
        match = .true.
        do iCoord = 1, 3
          call append( me  = IBMData%IBMSend_X%elemPos( iProc ),               &
            &          val = startPos+iCoord )
!write(IBM_logUnit(3),*)startPos+iCoord
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
        end do
        ! set the flag for the corresponding proc to true
        added(proc) = .true.
!        ! exit the proc loop
        exit proc_loop
      end if
    end do proc_loop
  end subroutine mus_IBM_fillSendPos_X
! **************************************************************************** !


! **************************************************************************** !
  !>
  !!
  !!
  subroutine mus_IBM_fillRecvPos_X( IBMData, globRecv, treeIDPos, startPos,    &
    &                               added, match )
    ! --------------------------------------------------------------------------
    !> IBM temporary datatype incl. map2glob and communicator for send and recv
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> global recv communicator
    type( tem_communication_type ), intent(in) :: globRecv
    !> element position in the level desc total list
    integer, intent(in) :: treeIDPos
    !> starting position of what to recv as elemPos
    integer, intent(in) :: startPos
    !> is this element added to the communication
    logical, intent(inout) :: added(:)
    logical :: match
    ! --------------------------------------------------------------------------
    ! counters
    integer :: iProc
    integer :: iCoord
    ! tmp variable storing the process position
    integer :: procPos
    ! tmp variable for the absolut process
    integer :: proc
    integer :: nVals
    ! --------------------------------------------------------------------------
    match = .false.

    proc_loop: do iProc = 1, IBMData%IBMRecv_X%nProcs
      procPos = IBMData%map2globSend( iProc )
      ! set the absolut proc position (+1 due to mpi numbering)
      proc = globRecv%proc( procPos )+1
      nVals = globRecv%elemPos( procPos )%nVals
      if(any(globRecv%elemPos( procPos )%val(1:nVals) == treeIDPos ) .and.     &
        & .not. added(proc) )then
        ! ... append the position of the velocity coordinates
        !     in the velocity array to the element positions in
        !     the recv buffer
        match = .true.
        do iCoord = 1, 3
          call append( me  = IBMData%IBMRecv_X%elemPos( iProc ),               &
            &          val = startPos+iCoord )
!write(IBM_logUnit(3),*)startPos+iCoord
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
        end do
        ! set the flag for the corresponding proc to true
        added(proc) = .true.
        ! exit the proc loop
        exit proc_loop
      end if
    end do proc_loop
  end subroutine mus_IBM_fillRecvPos_X
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine initializes all arrays in the mus_IBM_tmpData_type.
  !!
  subroutine mus_init_IBMData( me, IBM, totNeighs, nElems_fluid )
    ! --------------------------------------------------------------------------
    !> tmpData type to be initialized
    type( mus_IBM_tmpData_type ), intent(inout) :: me
    !> IBM datatype
    type( mus_IBM_type ), intent(inout) :: IBM
    !> total number of neighbors: number of surface points*stencil%QQ
    integer, intent(in) :: totNeighs
    !> number of fluid nodes on this process
    integer, intent(in) :: nElems_fluid
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iElem
    ! --------------------------------------------------------------------------

    ! allocate the arrays depending on the number of surface points
    allocate( me%vel_Xk( IBM%surfData%nUniquePoints_total*3 ))
    allocate( me%vel_Xk_ini( IBM%surfData%nUniquePoints_total*3 ))
    allocate( me%vel_Xk_surf( IBM%surfData%nUniquePoints_total*3 ))
    allocate( me%force_Xk( IBM%surfData%nUniquePoints_total*3 ))
    ! allocate the arrays depending on the number of surface points*stencil%QQ
    allocate( me%ptrs_neighs_Xk( totNeighs ))
    allocate( me%inaDelta_Xk( totNeighs ))
    ! and set them to 0
    me%vel_Xk      = 0._rk
    me%vel_Xk_ini  = 0._rk
    me%vel_Xk_surf = 0._rk
    me%force_Xk    = 0._rk
    me%inaDelta_Xk = 0._rk
!    call init( me = me%neighs_Xk, unique = .true. )

    ! allocate the array of growing arrays for all neighbors on this partition
    allocate( me%neighs_X( nElems_fluid ))
    ! allocate the array of growing arrays for the inamuro delta function
    ! results
    allocate( me%inaDelta_X( nElems_fluid ))

    ! initialize the growing arrays with size 0
    do iElem = 1, nElems_fluid
      call init( me%neighs_X( iElem ))
      call init( me%inaDelta_X( iElem ))
    end do
  end subroutine mus_init_IBMData
! **************************************************************************** !


! **************************************************************************** !
  !> This routine frees all temporary variables and destroys growing arrays
  !! as well as the communicators.
  !!
  subroutine mus_free_IBMData( me, commPattern, nElems_fluid )
    ! --------------------------------------------------------------------------
    !> tmpData type to be initialized
    type( mus_IBM_tmpData_type ), intent(inout) :: me
    !> communication pattern to be used
    type( tem_commPattern_type ), intent(inout) :: commPattern
    !> number of fluid nodes on this process
    integer, intent(in) :: nElems_fluid
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iElem
    integer :: iProc
    ! --------------------------------------------------------------------------

    ! deallocate the temporary arrays
    deallocate( me%vel_Xk )
    deallocate( me%vel_Xk_ini )
    deallocate( me%vel_Xk_surf )
    deallocate( me%force_Xk )
    deallocate( me%vel_X )
    deallocate( me%vel_X_ini )
    deallocate( me%force_X )
    deallocate( me%ptrs_neighs_Xk )
    deallocate( me%inaDelta_Xk )

    ! destroy the growing arrays ...
    do iElem = 1, nElems_fluid
      call destroy( me%neighs_X( iElem ))
      call destroy( me%inaDelta_X( iElem ))
    end do

    ! and deallocate the array of growing arrays ...
    deallocate( me%neighs_X )
    deallocate( me%inaDelta_X )

    ! destroy the dynamic array
!    call destroy( me%neighs_Xk )

    ! free all additional  buffers
    do iProc = 1, me%IBMSend_Xk%nProcs
      call commPattern%finBuf_real( me = me%IBMSend_Xk%buf_real( iProc ))
      call commPattern%finBuf_int( me = me%IBMSend_Xk%buf_int( iProc ))
      call destroy( me%IBMSend_Xk%elemPos( iProc ))
    end do
    do iProc = 1, me%IBMRecv_Xk%nProcs
      call commPattern%finBuf_real( me = me%IBMRecv_Xk%buf_real( iProc ))
      call commPattern%finBuf_int( me = me%IBMRecv_Xk%buf_int( iProc ))
      call destroy( me%IBMRecv_Xk%elemPos( iProc ))
    end do

    do iProc = 1, me%IBMSend_X%nProcs
      call commPattern%finBuf_real( me = me%IBMSend_X%buf_real( iProc ))
      call commPattern%finBuf_int( me = me%IBMSend_X%buf_int( iProc ))
      call destroy( me%IBMSend_X%elemPos( iProc ))
    end do
    do iProc = 1, me%IBMRecv_X%nProcs
      call commPattern%finBuf_real( me = me%IBMRecv_X%buf_real( iProc ))
      call commPattern%finBuf_int( me = me%IBMRecv_X%buf_int( iProc ))
      call destroy( me%IBMRecv_X%elemPos( iProc ))
    end do

    do iProc = 1, me%IBMSend_X_pdf%nProcs
      call commPattern%finBuf_real( me = me%IBMSend_X_pdf%buf_real( iProc ))
      call commPattern%finBuf_long( me = me%IBMSend_X_pdf%buf_long( iProc ))
      call destroy( me%IBMSend_X_pdf%elemPos( iProc ))
    end do
    do iProc = 1, me%IBMRecv_X_pdf%nProcs
      call commPattern%finBuf_real( me = me%IBMRecv_X_pdf%buf_real( iProc ))
      call commPattern%finBuf_long( me = me%IBMRecv_X_pdf%buf_long( iProc ))
      call destroy( me%IBMRecv_X_pdf%elemPos( iProc ))
    end do

    deallocate( me%IBMSend_Xk%buf_real )
    deallocate( me%IBMRecv_Xk%buf_real )
    deallocate( me%IBMSend_Xk%buf_int )
    deallocate( me%IBMRecv_Xk%buf_int )
    deallocate( me%IBMSend_Xk%elemPos )
    deallocate( me%IBMRecv_Xk%elemPos )
    deallocate( me%IBMSend_Xk%proc )
    deallocate( me%IBMRecv_Xk%proc )
    deallocate( me%IBMSend_Xk%nElemsProc )
    deallocate( me%IBMRecv_Xk%nElemsProc )

    deallocate( me%IBMSend_X%buf_real )
    deallocate( me%IBMRecv_X%buf_real )
    deallocate( me%IBMSend_X%buf_int )
    deallocate( me%IBMRecv_X%buf_int )
    deallocate( me%IBMSend_X%elemPos )
    deallocate( me%IBMRecv_X%elemPos )
    deallocate( me%IBMSend_X%proc )
    deallocate( me%IBMRecv_X%proc )
    deallocate( me%IBMSend_X%nElemsProc )
    deallocate( me%IBMRecv_X%nElemsProc )

    deallocate( me%IBMSend_X_pdf%buf_real )
    deallocate( me%IBMRecv_X_pdf%buf_real )
    deallocate( me%IBMSend_X_pdf%buf_long )
    deallocate( me%IBMRecv_X_pdf%buf_long )
    deallocate( me%IBMSend_X_pdf%elemPos )
    deallocate( me%IBMRecv_X_pdf%elemPos )
    deallocate( me%IBMSend_X_pdf%proc )
    deallocate( me%IBMRecv_X_pdf%proc )
    deallocate( me%IBMSend_X_pdf%nElemsProc )
    deallocate( me%IBMRecv_X_pdf%nElemsProc )

    deallocate( me%map2globSend )
    deallocate( me%map2globRecv )

    deallocate( me%map2globSend_pdf )
    deallocate( me%map2globRecv_pdf )
  end subroutine mus_free_IBMData
! **************************************************************************** !


! **************************************************************************** !
  !> This routine calculates the surface velocity for all local xk.
  !!
  subroutine mus_IBM_getSurfVel( me, IBMData, levelDesc, general, iLevel )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(in) :: me
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> level descriptor incl. ghost and fluid elements
    type( tem_levelDesc_type ), intent(in) :: levelDesc
    !> general info
    type(tem_general_type), intent(in) :: general
    !> the current level
    integer :: iLevel
    ! --------------------------------------------------------------------------
    ! counter
    integer :: iPoint
    ! tmp variable to store the point coordinates to
    real(kind=rk) :: pos(1,3)
    ! min and max position in the pointCoords array
    integer :: minPos, maxPos
    ! tmp variable to store the velocity
    real(kind=rk) :: vel(1,3)
    ! --------------------------------------------------------------------------

    IBMData%vel_Xk_surf = 0._rk

    do iPoint = 1, me%surfData%nUniquePoints_total
      if( me%movPredef .and.                                                   &
        &             me%surfData%parentIDs( iLevel )%ptrs( iPoint ) > 0 )then
        ! ... set the min and max pos
        minPos = (iPoint-1)*3+1
        maxPos = (iPoint-1)*3+3
        ! ... calculate the velocity at lagrangian surface points in
        !     physical units
        ! check wether the initial point coordinates shall be used ...
        if( me%useInitPos )then
          ! ... yes: use array backPointCoords and ...
          ! ... store the coordinates in a temporary variable
          pos(1,1:3) = me%surfData%backPointCoords(minPos:maxPos)
        else
          ! ... yes: use array pointCoords and ...
          ! ... store the coordinates in a temporary variable
          pos(1,1:3) = me%surfData%pointCoords(minPos:maxPos)
        end if
        vel = tem_spacetime_for( me    = me%velocity,            &
          &                      coord = pos,                    &
          &                      time  = general%simControl%now, &
          &                      n     = 1,                      &
          &                      nComp = 3                       )
        IBMData%vel_Xk_surf(minPos:maxPos) = vel(1,1:3)! /convFac%vel
      ! ... check if the parent element is a fluid element on this process
      else if( me%surfData%parentIDs( iLevel )%ptrs( iPoint ) > 0 .and.        &
        &      me%surfData%parentIDs( iLevel )%ptrs( iPoint ) <=               &
        &                                levelDesc%elem%nElems( eT_fluid ).and.&
        &      .not. me%movPredef )then
        ! ... set the min and max pos
        minPos = (iPoint-1)*3+1
        maxPos = (iPoint-1)*3+3
        ! ... calculate the velocity at lagrangian surface points in
        !     physical units
        ! check wether the initial point coordinates shall be used ...
        if( me%useInitPos )then
          ! ... yes: use array backPointCoords and ...
          ! ... store the coordinates in a temporary variable
          pos(1,1:3) = me%surfData%backPointCoords(minPos:maxPos)
        else
          ! ... yes: use array pointCoords and ...
          ! ... store the coordinates in a temporary variable
          pos(1,1:3) = me%surfData%pointCoords(minPos:maxPos)
        end if
        vel = tem_spacetime_for( me    = me%velocity,            &
          &                      coord = pos,                    &
          &                      time  = general%simControl%now, &
          &                      n     = 1,                      &
          &                      nComp = 3                       )
        IBMData%vel_Xk_surf(minPos:maxPos) = vel(1,1:3)! /convFac%vel
      end if
    end do
  end subroutine mus_IBM_getSurfVel
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine communicates all elements which just moved from the fluids
  !! to the halo elements.
  !!
  subroutine mus_IBM_commNewPos( IBMData, levelDesc, commPattern,              &
    &                            globTree, surfData, iLevel, comm )
    ! --------------------------------------------------------------------------
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> the level descriptor incl. the global send and receive buffers
    type( tem_levelDesc_type ), intent(inout) :: levelDesc
    !> communication pattern to be used
    type( tem_commPattern_type ), intent(inout) :: commPattern
    !> global tree information
    type( treelmesh_type ), intent(inout) :: globTree
    !> the surface data incl. the coordinates for the Xk
    type( tem_surfData_type ), intent(inout) :: surfData
    !> current level
    integer, intent(in) :: iLevel
    !> process description to use
    integer, intent(in) :: comm
    ! --------------------------------------------------------------------------
    ! counters
    integer :: iProc, iElem, iVal, iCoord, iPoint, iRecvProc
    ! tmp variable to store a position
    integer :: pos
    ! tmp variable for counting elements written before
    integer :: nElemPos
    ! counter for the number of processes attached to this process
    integer :: nProcs_send
    integer :: nProcs_send2
    ! adjacent processes to the current process
    integer, allocatable :: adjProcs_send(:)
    ! tmp array of growing arrays storing the Xk positions
    ! size: nProcs_send
    type( grw_intArray_type ), allocatable :: posXk(:)
    ! counter for the number of processes attached to this process
    integer :: nProcs_recv
    ! adjacent processes to the current process
    integer, allocatable :: adjProcs_recv(:)
    ! tmp array for storing the number of Xk to be received by this process
    integer, allocatable :: nElems_recv(:)
    ! array storing the Xk positions to be send in a linearized way
    integer, allocatable :: posXk_recv(:)
    ! array storing the Xk positions to be send in a linearized way
    integer, allocatable :: posXk_send(:)
    ! position of the first halo element in the levelDesc totallist
    integer :: firstHaloPos
    ! does the process participate in the communication?
    logical, allocatable :: procPart(:)
    ! --------------------------------------------------------------------------
    ! initialize nProcs_send with 0 and ...
    nProcs_send = 0
    nProcs_send2 = 0
    ! ... allocate the array of processes in the sendBuffer with the maximum
    !     number of processes from the global sendBuffer
    allocate( adjProcs_send( levelDesc%recvBuffer%nProcs ))
    ! ... allocate the array of growing arrays with the maximum number of
    !     processes from the global sendBuffer
    allocate( posXk( levelDesc%recvBuffer%nProcs ))
    ! ... allocate the map from the local to the global proc arrays
    allocate( IBMData%map2globSend( levelDesc%recvBuffer%nProcs ))
    ! ... set the position of the first halo in the levelDesc%total
    firstHaloPos = levelDesc%nElems-levelDesc%elem%nElems( eT_halo ) + 1
    ! ... allocate the array of logicals indicating wether a proc participates
    !     in the communication and initialize it
    allocate( procPart( levelDesc%recvBuffer%nProcs ))
    procPart = .false.


    ! --------------------------------------------------------------------------
    ! 1a. search for those Xk which just moved from a fluid element to a halo
    !     element on this proc and store their proc, the number of total procs
    !     and their position in the pointCoords array
    !
    ! therefor loop over the Xk ...
    do iPoint = 1, surfData%nUniquePoints_total
      ! ... if the parent is a halo element search in my receive buffer
      !     for the corresponding positions
      if( surfData%parentIDs(iLevel)%ptrs( iPoint ) >= firstHaloPos )then
        ! ... loop over the global number of processes and ...
        do iProc = 1, levelDesc%recvBuffer%nProcs
          ! ... search in the global recvBuffer%elemPos array for the process to
          !     send to
          pos_loop: do iVal = 1, levelDesc%recvBuffer%elemPos( iProc )%nVals
            ! ... if the position in the recvBuffer%elemPos equals to the
            !     position in the parentID ...
            if( levelDesc%recvBuffer%elemPos(iProc)%val(iVal) ==             &
              &                  surfData%parentIDs(iLevel)%ptrs( iPoint ))then
              write(IBM_logUnit(1),*) 'proc in old: ', &
                &                     levelDesc%recvBuffer%proc(iProc)
              ! ... and append them to the growing array
              call append( me  = posXk( iProc ),                               &
                &          val = iPoint )
              ! this process participates in the communication
              procPart( iProc ) = .true.
              ! this element is only received by one process so we can
              ! exit the loop in case we found the position
              exit pos_loop
            end if
          end do pos_loop
        end do
      end if
    end do
    ! fill the helper arrays for building the subcommunicators
    do iProc = 1, levelDesc%recvBuffer%nProcs
      ! if this process participates in the communication ...
      if( procPart( iProc ))then
        ! ... increase the process counter by 1
        nProcs_send = nProcs_send + 1
        ! ... add the corresponding globRecv process to the array
        adjProcs_send( nProcs_send ) = levelDesc%recvBuffer%proc( iProc )
        ! ... add the process position in the global buffer
        !     to map2globSend
        IBMData%map2globSend( nProcs_send ) = iProc
      end if
    end do

!write(IBM_logUnit(3),*)'after setting posXk'
!write(IBM_logUnit(3),*)'partFirst: ', globTree%Part_First
!write(IBM_logUnit(3),*)'partLast:  ', globTree%Part_Last

    ! maybe new implementation
    ! therefor loop over the Xk ...
    do iPoint = 1, surfData%nUniquePoints_total
      ! ... if the parent is a halo element search in my receive buffer
      !     for the corresponding positions
      if( surfData%parentIDs(iLevel)%ptrs( iPoint ) >= firstHaloPos )then
        ! ... loop over all participating procs and ...
        proc_loop: do iProc = 1, globTree%global%nParts
          ! ... check for the parent proc ...
          if( levelDesc%total( surfData%parentIDs(iLevel)%ptrs( iPoint ))      &
            & >= globTree%Part_First(iProc) .and.                              &
            & levelDesc%total( surfData%parentIDs(iLevel)%ptrs( iPoint ))      &
            & <= globTree%Part_Last(iProc))then
            ! ... loop over all procs in the levelDesc receive buffer and ...
            recvProc_loop: do iRecvProc = 1, levelDesc%recvBuffer%nProcs
              ! ... check for the corresponding iProc
              if( iProc == levelDesc%recvBuffer%proc(iRecvProc) )then
                write(IBM_logUnit(1),*)'proc in new: ', iRecvProc
                call append( me  = posXk( iRecvProc ), &
                  &          val = iPoint )
                procPart( iRecvProc ) = .true.
                exit recvProc_loop
              end if
            end do recvProc_loop
            exit proc_loop
          end if
        end do proc_loop
      end if
    end do
    ! fill the helper arrays for building the subcommunicators
    do iRecvProc = 1, levelDesc%recvBuffer%nProcs
      ! if this process participates in the communication ...
      if( procPart( iRecvProc ))then
        ! ... increase the proc counter by 1
        nProcs_send2 = nProcs_send2 + 1
        ! ... add the corresponding globRecv process to the array
        adjProcs_send( nProcs_send2 ) = levelDesc%recvBuffer%proc( iRecvProc )
        ! ... add the process position in the global buffer
        !     to map2globSend
        IBMData%map2globSend( nProcs_send2 ) = iRecvProc
      end if
    end do

    ! --------------------------------------------------------------------------
    ! 1b. send the number of elements their corresponding positions as well as
    !     actual new coordinates to the right processes
    !
    ! therefor use the global communicator and communicate the number of Xk
    ! that will be send

    ! allocate the array of number of elements to be recv with the maximum
    ! number of processes from the global receive buffer
    allocate( nElems_recv( levelDesc%sendBuffer%nProcs ))
    nElems_recv = 0

!  write(IBM_logUnit(3),*)'posXk: '
!  write(IBM_logUnit(3),*)posXk(1:levelDesc%recvBuffer%nProcs)%nVals
!
!  write(IBM_logUnit(3),*)'SendBuffer amount Xk: '
    ! initialize the global integer buffer
    do iProc = 1, levelDesc%sendBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%sendBuffer%buf_int( iProc ))
      call commPattern%initBuf_int( me    = levelDesc%sendBuffer%              &
        &                                                    buf_int( iProc ), &
        &                           pos   = (/iProc/),                         &
        &                           nVals = 1 )
!write(IBM_logUnit(3),*)'Proc:  ', levelDesc%SendBuffer%proc(iProc)
!write(IBM_logUnit(3),*)'nVals: ', levelDesc%SendBuffer%buf_int( iProc )%nVals
    end do

!write(IBM_logUnit(3),*)'RecvBuffer amount Xk: '
    do iProc = 1, levelDesc%recvBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%recvBuffer%buf_int( iProc ))
      call commPattern%initBuf_int( me    = levelDesc%recvBuffer%              &
        &                                                    buf_int( iProc ), &
        &                           pos   = (/iProc/),                         &
        &                           nVals = 1 )
!write(IBM_logUnit(3),*)'Proc:  ', levelDesc%recvBuffer%proc(iProc)
!write(IBM_logUnit(3),*)'nVals: ', levelDesc%recvBuffer%buf_int( iProc )%nVals
    end do

    ! ... exchange the amount of Xk to be communicated
    call commPattern%exchange_int( send         = levelDesc%sendBuffer,        &
      &                            recv         = levelDesc%recvBuffer,        &
      &                            state        = nElems_recv,                 &
      &                            message_flag = mess_amountXk,               &
      &                            send_state   = posXk( 1:levelDesc%          &
      &                                              recvBuffer%nProcs )%nVals,&
      &                            comm         = comm )

    ! ... free the buffers
    do iProc = 1, levelDesc%sendBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%sendBuffer%buf_int( iProc ))
    end do
    do iProc = 1, levelDesc%recvBuffer%nProcs
      call commPattern%finBuf_int( me = levelDesc%recvBuffer%buf_int( iProc ))
    end do

    ! --------------------------------------------------------------------------
    ! 2. Now we need to update the communication types such that the actual
    !    positions of the elements can be exchanged
    !
    ! therefor build the IBMSend_Xk and ...
    ! set the number of procs
    IBMData%IBMSend_Xk%nProcs = nProcs_send
    ! allocate the array of processes to send to and set them
    allocate( IBMData%IBMSend_Xk%proc( nProcs_send ))
    IBMData%IBMSend_Xk%proc = adjProcs_send(1:nProcs_send)
    ! allocate the array of integer buffers
    allocate( IBMData%IBMSend_Xk%buf_int( nProcs_send ))
    ! allocate the array of nElemsProc with the number of processes
    allocate( IBMData%IBMSend_Xk%nElemsProc( nProcs_send ))
    allocate( IBMData%IBMSend_Xk%elemPos( nProcs_send ))

    ! ... IBMRecv_Xk buffers
    allocate( adjProcs_recv( levelDesc%sendBuffer%nProcs ))
    allocate( IBMData%map2globRecv( levelDesc%sendBuffer%nProcs ))
    nProcs_recv = 0
    ! loop over the procs in the global receive buffer ...
    do iProc = 1, levelDesc%sendBuffer%nProcs
      ! ... if elements will be received by this proc ...
      if( nElems_recv( iProc ) > 0 )then
        ! ... increase the recv counter by 1
        nProcs_recv = nProcs_recv + 1
        ! add the process to the adjacent recv array
        adjProcs_recv( nProcs_recv ) = levelDesc%sendBuffer%proc( iProc )
        ! add this process to the map
        IBMData%map2globRecv( nProcs_recv ) = iProc
      end if
    end do
    ! set the number of procs
    IBMData%IBMRecv_Xk%nProcs = nProcs_recv
    ! allocate the array of processes to recv to and set them
    allocate( IBMData%IBMRecv_Xk%proc( nProcs_recv ))
    IBMData%IBMRecv_Xk%proc = adjProcs_recv(1:nProcs_recv)
    ! allocate the array of integer buffers
    allocate( IBMData%IBMRecv_Xk%buf_int( nProcs_recv ))
    ! allocate the array of nElemsProc with the number of processes
    allocate( IBMData%IBMRecv_Xk%nElemsProc( nProcs_recv ))
    allocate( IBMData%IBMRecv_Xk%elemPos( nProcs_recv ))

    ! ... allocate the array of Xk positions with the total number of elements
    !     to be received
    allocate( posXk_recv( sum( nElems_recv( : ))))
    ! ... and the one to be send
    allocate( posXk_send( sum( posXk( IBMData%                                 &
      &                                map2globSend( 1:nProcs_send ))%nVals )))
    ! ... loop over the processes and empty the growing arrays of positions in
    !     the send buffer
    ! ... fill the posXk_send array with the correct positions process wise
    !     and store the position in IBMSend_Xk%elemPos(iProc)
    ! ... initialize the integer buffer for the 2nd communication
    pos = 0
    do iProc = 1, nProcs_send
      call init( IBMData%IBMSend_Xk%elemPos( iProc ), 1)
      do iElem = 1, posXk( IBMData%map2globSend( iProc ))%nVals
        pos = pos + 1
        posXk_send( pos ) = posXk( IBMData%map2globSend( iProc ))%val(iElem)
        call append( me  = IBMData%IBMSend_Xk%elemPos( iProc ),                &
          &          val = pos )
      end do
      call commPattern%initBuf_int(                                            &
        &                   me    = IBMData%IBMSend_Xk%buf_int( iProc ),       &
        &                   pos   = IBMData%IBMSend_Xk%elemPos( iProc )%val,   &
        &                   nVals = IBMData%IBMSend_Xk%elemPos( iProc )%nVals )
    end do
    ! ... loop over the processes and empty the growing arrays of positions in
    !     the receive buffer
    ! ... set the correct positions in IBMRecv_Xk%elemPos(iProc)
    ! ... initialize the integer buffer for the 2nd communication
    pos = 0
    do iProc = 1, nProcs_recv
      call init( IBMData%IBMRecv_Xk%elemPos( iProc ), 1)
      do iElem = 1, nElems_recv( IBMData%map2globRecv( iProc ))
        pos = pos + 1
        call append( me  = IBMData%IBMRecv_Xk%elemPos( iProc ),                &
          &          val = pos )
      end do
      call commPattern%initBuf_int(                                            &
        &                   me    = IBMData%IBMRecv_Xk%buf_int( iProc ),       &
        &                   pos   = IBMData%IBMRecv_Xk%elemPos( iProc )%val,   &
        &                   nVals = IBMData%IBMRecv_Xk%elemPos( iProc )%nVals )
    end do

!write(IBM_logUnit(3),*)'Before exchanging positions1'

    ! Now exchange the positions in the Xk array
    call commPattern%exchange_int( send         = IBMData%IBMSend_Xk,          &
      &                            recv         = IBMData%IBMRecv_Xk,          &
      &                            state        = posXk_recv,                  &
      &                            message_flag = mess_posXk,                  &
      &                            send_state   = posXk_send,                  &
      &                            comm         = comm )

    ! ... free the integer buffers
    do iProc = 1, nProcs_send
      call commPattern%finBuf_int( me    = IBMData%IBMSend_Xk%buf_int( iProc ))
    end do
    do iProc = 1, nProcs_recv
      call commPattern%finBuf_int( me    = IBMData%IBMRecv_Xk%buf_int( iProc ))
    end do

    ! --------------------------------------------------------------------------
    ! 3. Now we need to update the communication types such that the buffer for
    !    reals is initialized and the right elemPos for 3 reals per Xk are set
    !    the right way
    !
    ! therefor ...
    ! ... allocate the real buffers in the IBMSend_Xk and IBMRecv_Xk with the
    !     corresponding number of processes
    allocate( IBMData%IBMSend_Xk%buf_real( nProcs_send ))
    allocate( IBMData%IBMRecv_Xk%buf_real( nProcs_recv ))
    ! ... loop over the send processes and set the new element positions and
    !     initialize the real buffer
    pos = 0
    do iProc = 1, nProcs_send
      call empty( IBMData%IBMSend_Xk%elemPos( iProc ))
      do iElem = 1, posXk( IBMData%map2globSend( iProc ))%nVals
        do iCoord = 1, 3
          call append( me  = IBMData%IBMSend_Xk%elemPos( iProc ),              &
            &          val = ( posXk( IBMData%map2globSend(iProc) )%val(iElem) &
            &                - 1)*3 + iCoord )
        end do
      end do
      ! ... initialize the real buffer
      call commPattern%initBuf_real(                                           &
        &                    me    = IBMData%IBMSend_Xk%buf_real( iProc ),     &
        &                    pos   = IBMData%IBMSend_Xk%elemPos( iProc )%val,  &
        &                    nVals = IBMData%IBMSend_Xk%elemPos( iProc )%nVals )
    end do
    ! ... loop over the receive processes and set the new element positions and
    !     initialize the real buffer
    nElemPos = 0
    do iProc = 1, nProcs_recv
      call empty( IBMData%IBMRecv_Xk%elemPos( iProc ))
      do iElem = 1, nElems_recv( IBMData%map2globRecv( iProc ))
        do iCoord = 1, 3
          call append( me  = IBMData%IBMRecv_Xk%elemPos( iProc ),              &
            &          val = (posXk_recv(nElemPos+iElem)-1)*3+iCoord )
        end do
      end do
      ! ... update the position counter
      nElemPos = nElemPos + nElems_recv( IBMData%map2globRecv( iProc ))
      ! ... initialize the real buffer
      call commPattern%initBuf_real(                                           &
        &                    me    = IBMData%IBMRecv_Xk%buf_real( iProc ),     &
        &                    pos   = IBMData%IBMRecv_Xk%elemPos( iProc )%val,  &
        &                    nVals = IBMData%IBMRecv_Xk%elemPos( iProc )%nVals )
    end do

    ! Now communicate the new coordinates ...
    call commPattern%exchange_real( send         = IBMData%IBMSend_Xk,         &
      &                             recv         = IBMData%IBMRecv_Xk,         &
      &                             state        = surfData%pointCoords,       &
      &                             message_flag = mess_pointsXk,              &
      &                             comm         = comm )

!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))
!write(IBM_logUnit(3),*)'surf vel before first commu'
!do iProc =1, surfData%nUniquePoints_total
!  write(IBM_logUnit(3),*)IBMData%vel_Xk_surf( (iProc-1)*3+1:(iProc-1)*3+3)
!end do

    ! ... and exchange the surface velocity based on the old Xk ...
    call commPattern%exchange_real( send         = IBMData%IBMSend_Xk,         &
      &                             recv         = IBMData%IBMRecv_Xk,         &
      &                             state        = IBMData%vel_Xk_surf,        &
      &                             message_flag = mess_velXk,                 &
      &                             comm         = comm )

!write(IBM_logUnit(3),*)'surf vel after first commu'
!do iProc =1, surfData%nUniquePoints_total
!  write(IBM_logUnit(3),*)IBMData%vel_Xk_surf( (iProc-1)*3+1:(iProc-1)*3+3)
!end do
!call tem_horizontalSpacer(fUnit=IBM_logUnit(3))

    ! ... and initialize the parent elements
    call tem_init_surfData( me        = surfData,                              &
      &                     levelDesc = levelDesc,                             &
      &                     globTree  = globTree,                              &
      &                     iLevel    = iLevel )

    ! --------------------------------------------------------------------------
    ! 4a. Now each Xk has a fluid parentID on a single proc.
    !     Search for all these Xk and store their proc, the number of total
    !     procs and their position in the pointCoords array
    !
    ! therefor reset the necessary arrays and variables and free the growing
    ! array of positions
    do iProc = 1, levelDesc%recvBuffer%nProcs
      call empty( posXk( iProc ))
    end do

    do iProc = 1, IBMData%IBMSend_Xk%nProcs
      call empty( IBMData%IBMSend_Xk%elemPos( iProc ))
      call commPattern%finBuf_int( me  = IBMData%IBMSend_Xk%buf_int( iProc ))
      call commPattern%finBuf_real( me = IBMData%IBMSend_Xk%buf_real( iProc ))
    end do

    do iProc = 1, IBMData%IBMRecv_Xk%nProcs
      call empty( IBMData%IBMRecv_Xk%elemPos( iProc ))
      call commPattern%finBuf_int( me = IBMData%IBMRecv_Xk%buf_int( iProc ))
      call commPattern%finBuf_real( me = IBMData%IBMRecv_Xk%buf_real( iProc ))
    end do

    ! deallocate maps
    deallocate( IBMData%map2globSend )
    deallocate( IBMData%map2globRecv )


  end subroutine mus_IBM_commNewPos
! **************************************************************************** !


! **************************************************************************** !
  !> This routine dumps the timings%timedat to disc
  !!
  subroutine mus_finishIBM( me, params, useTime )
    ! --------------------------------------------------------------------------
    !> global IBM datatype incl. array of IBM datatypes
    type( mus_IBM_globType ), intent(inout) :: me
    !> global parameters
    type( mus_param_type ), intent(in)   :: params
    !> use the timestamps when dumping the info?
    logical, intent(in) :: useTime
    ! --------------------------------------------------------------------------
    integer       :: iIBM, iTimer
    integer       :: iUnit
    real(kind=rk) :: cTime
    character(len=PathLen) :: header
    character(len=PathLen) :: output
    character(len=PathLen) :: filename
    ! timestamp for the filename
    character(len=16)     :: timeStamp
    ! --------------------------------------------------------------------------

    do iIBM = 1, me%nIBMs
      header = ''
      output = ''
      do iTimer = 1, me%IBM(iIBM)%timings%timedat%nTimers
        cTime = tem_getMaxTimerVal( me          = me%IBM(iIBM)%timings%timedat,&
          &                         timerHandle = iTimer,                      &
          &                         comm        = params%general%proc%comm     )
        write(header,'(a,a16)') trim(header),                               &
          &                     trim( me%IBM(iIBM)%timings%label%val(iTimer))
        write(output, '(a,en16.2)') trim(output), cTime
        ! reset the timer
        call tem_resetTimer( timerHandle = iTimer,                      &
          &                  me          = me%IBM(iIBM)%timings%timedat )
      end do
      ! now the root process writes the data to disc
      if( params%general%proc%rank == params%general%proc%root ) then
        ! get a new unit and ...
        iUnit = newunit()
        if( useTime )then
          timestamp = trim(tem_time_sim_stamp( params%general%simControl%now ))
          ! ... construct the filename using the timestamp and ...
          write( filename, '(i4.4,a,i2.2,a)')params%general%proc%comm_size,    &
            &               '_',iIBM,'_IBM_MAXTimings_'//trim(timestamp)//'.res'
        else
          ! ... construct the filename and ...
          write( filename, '(i4.4,a,i2.2,a)')params%general%proc%comm_size,'_',&
            &                                         iIBM,'_IBM_MAXTimings.res'
        end if
        ! ... open the output file
        open( unit = iUnit, file = 'timings/'//trim( filename ))
        write( iUnit, '(a)' )trim(header)
        write( iUnit, '(a)' )trim(output)
        close( unit = iUnit )
      end if
    end do

  end subroutine mus_finishIBM
! **************************************************************************** !


! **************************************************************************** !
  !> This routine finishes the buffers for Xk and X_pdf. This is moved to
  !! a seperate routine since both buffers depend on a local communication
  !! which should be done nearby the global synchronization point (mus_exchange)
  !!
  subroutine mus_IBMFinishBuff( me, IBMData, levelDesc, commPattern, globTree, &
    &                           iLevel, comm, stencil )
    ! --------------------------------------------------------------------------
    !> datatype to store the surface information
    type( mus_IBM_type ), intent(inout) :: me
    !> tmp IBMData type to be filled
    type( mus_IBM_tmpData_type ), intent(inout) :: IBMData
    !> the level descriptor incl. the global send and receive buffers
    type( tem_levelDesc_type ), intent(inout) :: levelDesc
    !> communication pattern to be used
    type( tem_commPattern_type ), intent(inout) :: commPattern
    !> global tree information
    type( treelmesh_type ), intent(inout) :: globTree
    !> current level
    integer, intent(in) :: iLevel
    !> MPI communicator
    integer, intent(in) :: comm
    !> array of stencils (1 is the fluid stencil)
    type( tem_stencilHeader_type ), intent(in) :: stencil(:)
    ! --------------------------------------------------------------------------
    ! counters
    integer :: iProc, iElem, iCoord
    ! tmp variable to store a position
    integer :: pos
    ! tmp variable for counting elements written before
    integer :: nElemPos
    ! array storing the Xk positions to be send in a linearized way
    integer, allocatable :: posXk_recv(:)
    ! array storing the Xk positions to be send in a linearized way
    integer, allocatable :: posXk_send(:)
    ! array storing the treeIDs to be send in a linearized way
    integer(kind=long_k), allocatable :: treeIDs_recv(:)
    ! array storing the treeIDs to be received in a linearized way
    integer(kind=long_k), allocatable :: treeIDs_send(:)
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    !       F I N I S H     T H E     B U F F E R    F O R    T H E    Xk      !
    ! --------------------------------------------------------------------------

    ! start the IBMbuildXk timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(2)  )

    ! start the IBMbuildXk timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(9)  )

    ! ... allocate the array of Xk positions with the total number of
    !     elements to be received
    if( allocated( posXk_recv )) &
      & deallocate( posXk_recv )
    allocate( posXk_recv( sum( IBMData%IBMRecv_Xk%nElemsProc( : ))))
    ! ... and the one to be send
    if( allocated( posXk_send )) &
      & deallocate( posXk_send )
    allocate( posXk_send( sum( IBMData%posXk( IBMData%map2globSend( 1:         &
      &                                   IBMData%IBMSend_Xk%nProcs ))%nVals)))
    ! ... loop over the processes and empty the growing arrays of positions in
    !     the send buffer
    ! ... fill the posXk_send array with the correct positions process wise
    !     and store the position in IBMSend_Xk%elemPos(iProc)
    ! ... initialize the integer buffer for the 2nd communication
    pos = 0
    do iProc = 1, IBMData%IBMSend_Xk%nProcs
      call init( IBMData%IBMSend_Xk%elemPos( iProc ), 1)
      do iElem = 1, IBMData%posXk( IBMData%map2globSend( iProc ))%nVals
        pos = pos + 1
        posXk_send(pos) = IBMData%posXk(IBMData%map2globSend(iProc))%val(iElem)
        call append( me  = IBMData%IBMSend_Xk%elemPos( iProc ),                &
          &          val = pos )
      end do
      call commPattern%initBuf_int(                                            &
        &                   me    = IBMData%IBMSend_Xk%buf_int( iProc ),       &
        &                   pos   = IBMData%IBMSend_Xk%elemPos( iProc )%val,   &
        &                   nVals = IBMData%IBMSend_Xk%elemPos( iProc )%nVals )
    end do
    ! ... loop over the processes and empty the growing arrays of positions in
    !     the receive buffer
    ! ... set the correct positions in IBMRecv_Xk%elemPos(iProc)
    ! ... initialize the integer buffer for the 2nd communication
    pos = 0
    do iProc = 1, IBMData%IBMRecv_Xk%nProcs
      call init( IBMData%IBMRecv_Xk%elemPos( iProc ), 1)
      do iElem = 1, IBMData%IBMRecv_Xk%nElemsProc(iProc)
        pos = pos + 1
        call append( me  = IBMData%IBMRecv_Xk%elemPos( iProc ),                &
          &          val = pos )
      end do
      call commPattern%initBuf_int(                                            &
        &                   me    = IBMData%IBMRecv_Xk%buf_int( iProc ),       &
        &                   pos   = IBMData%IBMRecv_Xk%elemPos( iProc )%val,   &
        &                   nVals = IBMData%IBMRecv_Xk%elemPos( iProc )%nVals )
    end do

    ! stop the IBMbuildXk timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(2)  )

    ! stop the IBMbuildXk timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(9)  )

    ! Now exchange the positions in the Xk array
    call commPattern%exchange_int( send         = IBMData%IBMSend_Xk, &
      &                            recv         = IBMData%IBMRecv_Xk, &
      &                            state        = posXk_recv,         &
      &                            message_flag = mess_posXk2,        &
      &                            send_state   = posXk_send,         &
      &                            comm         = comm                )

    ! start the IBMbuildXk timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(2)  )

    ! start the IBMbuildXk timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(9)  )

    ! ... free the integer buffers
    do iProc = 1, IBMData%IBMSend_Xk%nProcs
      call commPattern%finBuf_int( me    = IBMData%IBMSend_Xk%buf_int( iProc ))
    end do
    do iProc = 1, IBMData%IBMRecv_Xk%nProcs
      call commPattern%finBuf_int( me    = IBMData%IBMRecv_Xk%buf_int( iProc ))
    end do

    ! --------------------------------------------------------------------------
    ! Now we need to update the communication types such that the buffer for
    ! reals is initialized and the right elemPos for 3 reals per Xk are set
    ! the right way
    !
    ! therefor ...
    ! ... reallocate the real buffers in the IBMSend_Xk and IBMRecv_Xk with the
    !     corresponding number of processes
    if( allocated( IBMData%IBMSend_Xk%buf_real )) &
      & deallocate( IBMData%IBMSend_Xk%buf_real )
    allocate( IBMData%IBMSend_Xk%buf_real( IBMData%IBMSend_Xk%nProcs ))
    if( allocated( IBMData%IBMRecv_Xk%buf_real )) &
      & deallocate( IBMData%IBMRecv_Xk%buf_real )
    allocate( IBMData%IBMRecv_Xk%buf_real( IBMData%IBMRecv_Xk%nProcs ))
    ! ... loop over the send processes and set the new element positions and
    !     initialize the real buffer
    pos = 0
    do iProc = 1, IBMData%IBMSend_Xk%nProcs
      call empty( IBMData%IBMSend_Xk%elemPos( iProc ))
      do iElem = 1, IBMData%posXk( IBMData%map2globSend( iProc ))%nVals
        do iCoord = 1, 3
          call append( me  = IBMData%IBMSend_Xk%elemPos( iProc ),              &
            &          val = ( IBMData%posXk( IBMData%map2globSend( iProc ))%  &
            &                  val(iElem) - 1)*3 + iCoord )
        end do
      end do
      ! ... initialize the real buffer
      call commPattern%initBuf_real(                                           &
        &                    me    = IBMData%IBMSend_Xk%buf_real( iProc ),     &
        &                    pos   = IBMData%IBMSend_Xk%elemPos( iProc )%val,  &
        &                    nVals = IBMData%IBMSend_Xk%elemPos( iProc )%nVals )
    end do
    ! ... loop over the receive processes and set the new element positions and
    !     initialize the real buffer
    nElemPos = 0
    do iProc = 1, IBMData%IBMRecv_Xk%nProcs
      call empty( IBMData%IBMRecv_Xk%elemPos( iProc ))
      do iElem = 1, IBMData%IBMRecv_Xk%nElemsProc(iProc)
        do iCoord = 1, 3
          call append( me  = IBMData%IBMRecv_Xk%elemPos( iProc ),              &
            &          val = (posXk_recv(nElemPos+iElem)-1)*3+iCoord )
        end do
      end do
      ! ... update the position counter
      nElemPos = nElemPos + IBMData%IBMRecv_Xk%nElemsProc( iProc )
      ! ... initialize the real buffer
      call commPattern%initBuf_real(                                           &
        &                    me    = IBMData%IBMRecv_Xk%buf_real( iProc ),     &
        &                    pos   = IBMData%IBMRecv_Xk%elemPos( iProc )%val,  &
        &                    nVals = IBMData%IBMRecv_Xk%elemPos( iProc )%nVals )
    end do

    ! --------------------------------------------------------------------------
    ! Now we the buffers are ready to communicate all infromation on the
    ! surface points Xk (e.g. velocity, force, coordinates, ...)

    ! in case of a non predefined motion ...
    if( .not. me%movPredef )then
      ! ... communicate the new coordinates ...
      call commPattern%exchange_real( send         = IBMData%IBMSend_Xk,       &
        &                             recv         = IBMData%IBMRecv_Xk,       &
        &                             state        = me%surfData%pointCoords,  &
        &                             message_flag = iLevel,                   &
        &                             comm         = comm )

      ! ... and exchange the surface velocity based on the old Xk
      call commPattern%exchange_real( send         = IBMData%IBMSend_Xk,       &
        &                             recv         = IBMData%IBMRecv_Xk,       &
        &                             state        = IBMData%vel_Xk_surf,      &
        &                             message_flag = mess_velXk2,              &
        &                             comm         = comm )

      ! ... and initialize the parent elements as well as ...
      call tem_init_surfData( me        = me%surfData,                         &
        &                     levelDesc = levelDesc,                           &
        &                     globTree  = globTree,                            &
        &                     iLevel    = iLevel )

      ! ... update the properties from the levelDesc to the global tree
      call tem_updateTree_properties( levelDesc = levelDesc,                   &
        &                             Tree  = globTree )
    end if ! .not. movPredef

    ! deallocate temporary arrays and destroy growing arrays
    deallocate( posXk_recv )
    deallocate( posXk_send )

    do iProc = 1, levelDesc%sendBuffer%nProcs
      call destroy( IBMData%posXk( iProc ))
    end do

    deallocate( IBMData%posXk )

    ! stop the IBMbuildXk timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(2)  )

    ! stop the IBMbuildXk timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(9)  )


    ! --------------------------------------------------------------------------
    !  F I N I S H     T H E     B U F F E R    F O R    T H E    P D F S      !
    ! --------------------------------------------------------------------------

    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(3)  )

    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(12) )

    ! ... allocate the array of Xk positions with the total number of
    !     elements to be received
    allocate( treeIDs_recv( sum( IBMData%IBMRecv_X_pdf%nElemsProc( : ))))
    ! ... and the one to be send
    allocate( treeIDs_send( sum( IBMData%treeIDs( IBMData%map2globSend_pdf(    &
      &                             1:IBMData%IBMSend_X_pdf%nProcs ))%nVals )))
    ! ... loop over the processes and empty the growing arrays of positions in
    !     the send buffer
    ! ... fill the treeIDs_send array with the correct positions process wise
    !     and store the position in IBMSend_X_pdf%elemPos(iProc)
    ! ... initialize the integer buffer for the 2nd communication
    pos = 0
    do iProc = 1, IBMData%IBMSend_X_pdf%nProcs
      call init( IBMData%IBMSend_X_pdf%elemPos( iProc ), 1)
      do iElem = 1, IBMData%treeIDs( IBMData%map2globSend_pdf( iProc ))%nVals
        pos = pos + 1
        treeIDs_send( pos ) = levelDesc%total( IBMData%treeIDs( IBMData%       &
          &                            map2globSend_pdf( iProc ))%val( iElem ))
        call append( me  = IBMData%IBMSend_X_pdf%elemPos( iProc ),             &
          &          val = pos )
      end do
      call commPattern%initBuf_long(                                           &
        &               me    = IBMData%IBMSend_X_pdf%buf_long( iProc ),       &
        &               pos   = IBMData%IBMSend_X_pdf%elemPos( iProc )%val,    &
        &               nVals = IBMData%IBMSend_X_pdf%elemPos( iProc )%nVals )
    end do
    ! ... loop over the processes and empty the growing arrays of positions in
    !     the receive buffer
    ! ... set the correct positions in IBMRecv_X_pdf%elemPos(iProc)
    ! ... initialize the integer buffer for the 2nd communication
    pos = 0
    do iProc = 1, IBMData%IBMRecv_X_pdf%nProcs
      call init( IBMData%IBMRecv_X_pdf%elemPos( iProc ), 1)
      do iElem = 1, IBMData%IBMRecv_X_pdf%nElemsProc( iProc )
        pos = pos + 1
        call append( me  = IBMData%IBMRecv_X_pdf%elemPos( iProc ),             &
          &          val = pos )
      end do
      call commPattern%initBuf_long(                                           &
        &               me    = IBMData%IBMRecv_X_pdf%buf_long( iProc ),       &
        &               pos   = IBMData%IBMRecv_X_pdf%elemPos( iProc )%val,    &
        &               nVals = IBMData%IBMRecv_X_pdf%elemPos( iProc )%nVals )
    end do

    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(3)  )

    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(12) )

    ! Now exchange the treeIDs to be send
    call commPattern%exchange_long( send         = IBMData%IBMSend_X_pdf, &
      &                             recv         = IBMData%IBMRecv_X_pdf, &
      &                             state        = treeIDs_recv,          &
      &                             message_flag = mess_posXk2,           &
      &                             send_state   = treeIDs_send,          &
      &                             comm         = comm                   )

    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(3)  )

    ! start the IBMbuildX timer
    call tem_startTimer( me          = me%timings%timedat, &
      &                  timerHandle = me%timerHandles(12) )

    ! ... free the integer buffers
    do iProc = 1, IBMData%IBMSend_X_pdf%nProcs
      call commPattern%finBuf_long( me = IBMData%IBMSend_X_pdf%buf_long(iProc))
    end do
    do iProc = 1, IBMData%IBMRecv_X_pdf%nProcs
      call commPattern%finBuf_long( me = IBMData%IBMRecv_X_pdf%buf_long(iProc))
    end do

    ! ... reallocate the real buffers in the IBMSend_X_pdf and IBMRecv_X_pdf
    !     with the corresponding number of processes
    allocate( IBMData%IBMSend_X_pdf%buf_real( IBMData%IBMSend_X_pdf%nProcs ))
    allocate( IBMData%IBMRecv_X_pdf%buf_real( IBMData%IBMRecv_X_pdf%nProcs ))
    ! ... loop over the send processes and set the new element positions and
    !     initialize the real buffer
    do iProc = 1, IBMData%IBMSend_X_pdf%nProcs
      call empty( IBMData%IBMSend_X_pdf%elemPos( iProc ))
      do iElem = 1, IBMData%treeIDs( IBMData%map2globSend_pdf( iProc ))%nVals
        do iCoord = 1, stencil(1)%QQ
          call append( me  = IBMData%IBMSend_X_pdf%elemPos( iProc ),           &
            &          val = ( IBMData%treeIDs( IBMData%map2globSend_pdf(      &
            &                  iProc ))%val(iElem)-1)*stencil(1)%QQ + iCoord )
        end do
      end do
      ! ... initialize the real buffer
      call commPattern%initBuf_real(                                           &
        &                 me    = IBMData%IBMSend_X_pdf%buf_real( iProc ),     &
        &                 pos   = IBMData%IBMSend_X_pdf%elemPos( iProc )%val,  &
        &                 nVals = IBMData%IBMSend_X_pdf%elemPos( iProc )%nVals )
    end do
    ! ... loop over the receive processes and set the new element positions and
    !     initialize the real buffer
    nElemPos = 0
    do iProc = 1, IBMData%IBMRecv_X_pdf%nProcs
      call empty( IBMData%IBMRecv_X_pdf%elemPos( iProc ))
      do iElem = 1, IBMData%IBMRecv_X_pdf%nElemsProc( iProc )
        ! search in the halos of the corresponding proc for the position and ...
        pos = tem_treeIDinTotal( treeIDs_recv(nElemPos+iElem), levelDesc,      &
          &                                                        eT_halo)
        ! ... append the coordinates of the state vector to the recv buffer
        do iCoord = 1, stencil(1)%QQ
          call append( me  = IBMData%IBMRecv_X_pdf%elemPos( iProc ),           &
            &          val = (pos-1)*stencil(1)%QQ+iCoord )
        end do
      end do
      nElemPos = nElemPos + IBMData%IBMRecv_X_pdf%nElemsProc( iProc )
      ! ... initialize the real buffer
      call commPattern%initBuf_real(                                           &
        &                 me    = IBMData%IBMRecv_X_pdf%buf_real( iProc ),     &
        &                 pos   = IBMData%IBMRecv_X_pdf%elemPos( iProc )%val,  &
        &                 nVals = IBMData%IBMRecv_X_pdf%elemPos( iProc )%nVals )
    end do

    ! allocate the arrays for the force, velocity and initial velocity on the
    ! eulerian elements with neighs_Xk%nVals * 3
    allocate( IBMData%vel_X( me%neighs_Xk%nVals*3 ))
    allocate( IBMData%vel_X_ini( me%neighs_Xk%nVals*3 ))
    allocate( IBMData%force_X( me%neighs_Xk%nVals*3 ))

    ! initialize the arrays with 0
    IBMData%vel_X = 0._rk
    IBMData%vel_X_ini = 0._rk
    IBMData%force_X = 0._rk

    ! allocate the array of real buffers with the number of procs
    allocate( IBMData%IBMSend_X%buf_real( IBMData%IBMSend_X%nProcs ))
    allocate( IBMData%IBMRecv_X%buf_real( IBMData%IBMRecv_X%nProcs ))

    ! Now the actual send and receive buffers are initialized
    do iProc = 1, IBMData%IBMSend_X%nProcs
      call commPattern%initBuf_real(                                           &
        &                    me    = IBMData%IBMSend_X%buf_real( iProc ),      &
        &                    pos   = IBMData%IBMSend_X%elemPos( iProc )%val(1: &
        &                           IBMData%IBMSend_X%elemPos( iProc )%nVals), &
        &                    nVals = IBMData%IBMSend_X%elemPos( iProc )%nVals )
    end do

    do iProc = 1, IBMData%IBMRecv_X%nProcs
      call commPattern%initBuf_real(                                           &
        &                    me    = IBMData%IBMRecv_X%buf_real( iProc ),      &
        &                    pos   = IBMData%IBMRecv_X%elemPos( iProc )%val(1: &
        &                           IBMData%IBMRecv_X%elemPos( iProc )%nVals), &
        &                    nVals = IBMData%IBMRecv_X%elemPos( iProc )%nVals )
    end do

    ! deallocate temporary arrays and destroy growing arrays
    deallocate( treeIDs_recv )
    deallocate( treeIDs_send )

    do iProc = 1, levelDesc%sendBuffer%nProcs
      call destroy( IBMData%treeIDs( iProc ))
    end do

    deallocate( IBMData%treeIDs )

    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(3)  )

    ! stop the IBMbuildX timer
    call tem_stopTimer( me          = me%timings%timedat, &
      &                 timerHandle = me%timerHandles(12) )
  end subroutine mus_IBMFinishBuff
! **************************************************************************** !


! **************************************************************************** !
  subroutine mus_unload_IBM( me, proc, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    !> IBM data type
    type( mus_IBM_type), intent(inout) :: me(:)
    !> Global parameters
    type( tem_comm_env_type ), intent(in) :: proc
    !> Level range
    integer, intent(in) :: minLevel, maxLevel
    ! --------------------------------------------------------------------------
    integer :: iIBM
    character(len=PathLen) :: stlFilename
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    !          Dump the stl files and free the arryas in case of IBM           !
    ! --------------------------------------------------------------------------
    !
    ! This part has to be called BEFORE the mus_balance call since
    ! the tree is unloaded in mus_balance. params%general%balance%fileID
    ! will be updated in mus_balance as well, so only temporarily update
    ! it here.
    !
    if ( any( me(:)%active )) then
      do iIBM = 1, size( me )
        ! construct the filename for dumping the stl
        ! structure for Scheme 1, Field 1, fileID 1, iIBM 1:
        !           -> balanceDirectory/surface_01_01_01_1
        ! file extension '.stl' will be added in the call tem_dump_stlb
        ! fileID = mod(params%general%balance%fileID,2)+1
        write( stlFilename,'(a,i2.2,a,i2.2,a,i1)')                  &
          ! &   trim(params%general%balance%restart%controller%writePrefix)//&
          &   './balance/surface_',iIBM,'_',fileID
        ! dump the stl files for each scheme and field
        call tem_dump_stlb( outprefix = trim(stlFilename),              &
          &                 nodes     = me(iIBM)%surfData%pointCoords, &
          &                 triangles = me(iIBM)%surfData%trias,       &
          &                 proc      = proc )
        call tem_freeSurfData( me       = me(iIBM)%surfData, &
          &                    minLevel = minLevel, maxLevel = maxLevel )
      end do
    end if

  end subroutine mus_unload_IBM
! **************************************************************************** !


! **************************************************************************** !
  subroutine mus_reload_IBM( me, iField, LevelDesc, tree )
    ! --------------------------------------------------------------------------
    !> IBM data type
    type( mus_IBM_type), intent(inout) :: me(:)
    !> scheme and field number
    integer, intent(in) :: iField
    !> Global parameters
    ! type( tem_balance_type ), intent(in) :: balance
    !> Level desc
    type( tem_levelDesc_type ), intent(inout) :: LevelDesc(:)
    !> tree
    type( treelmesh_type ), intent(inout) :: tree
    ! --------------------------------------------------------------------------
    integer :: iIBM, nIBMs
    character(len=PathLen) :: stlFilename
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    !    Read the stl files and identify the parent elements in case of IBM    !
    ! --------------------------------------------------------------------------

    if( any( me(:)%active ))then
      nIBMs = size( me )
      do iIBM = 1, nIBMs
        ! construct the filename for dumping the stl
        ! structure for Scheme 1, Field 1, fileID 1, iIBM 1:
        !           -> balanceDirectory_surface_01_01_01_1
        ! file extension '.stl' has to be added beforehand
        write( stlFilename,'(a,i2.2,a,i2.2,a,i1,a)')                &
          ! &  trim(balance%restart%controller%writePrefix)// &
          &     './balance/surface_',iField,'_',iIBM,'_',fileID,'.stl'
        ! allocate the stlHead in the tem_surfaceData_type with 1
        allocate( me(iIBM)%surfData%stlHead(1))
        ! and set the filename
        me(iIBM)%surfData%stlHead(1)%filename = trim(stlFilename)
        ! read the stl data
        call tem_readAndUnify_surfData( me = me(iIBM)%surfData )
      end do ! iIBM
      ! set the parentIDs
      call mus_IBM_setParentIDs(  nIBMs = nIBMs,     &
        &                     me        = me,        &
        &                     levelDesc = levelDesc, &
        &                     tree      = tree )
    end if ! IBM is active
    ! --------------------------------------------------------------------------

  end subroutine mus_reload_IBM
! **************************************************************************** !

end module mus_IBM_module
! **************************************************************************** !
