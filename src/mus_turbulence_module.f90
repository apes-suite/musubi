! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2021 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! ****************************************************************************** !
!> This module contains data types, function and routines for turbulence
!! model. Especially the routine to compute turbulent eddy viscosity for
!! different turbulence model
!!
!! author: Kannan Masilamani
module mus_turbulence_module
  ! include treelm modules
  use env_module,                only: rk, labelLen
  use tem_tools_module,          only: upper_to_lower, &
    &                                  tem_horizontalSpacer
  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logUnit
  use tem_comm_module,           only: tem_communication_type, &
    &                                  tem_commPattern_type,   &
    &                                  tem_comm_init
  use tem_construction_module,   only: tem_levelDesc_type

  ! include aotus modules
  use aotus_module,              only: aoterr_Fatal,       &
    &                                  aoterr_NonExistent, &
    &                                  flu_State,          &
    &                                  aoterr_WrongType,   &
    &                                  aot_get_val
  use aot_table_module,          only: aot_table_open, &
    &                                  aot_table_close

  ! include musubi modules
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_gradData_module,       only: mus_gradData_type, mus_Grad_type

  implicit none
  private

  public :: mus_turbulence_type
  public :: mus_turbulence_config_type
  public :: mus_turbulence_data_type
  public :: mus_init_turbulenceData
  public :: mus_load_turbulence
  public :: mus_turb_calcVisc
  public :: mus_turb_updateViscOfTurbWall

  !> Contains large Eddy Turbulence (LES) model coefficients
  type les_coeff_type
    !> Smagorinsky constant.
    !! C_s = sqrt(C_k sqrt(C_k/C_e) ) = 0.16778594 = 0.17
    real(kind=rk) :: C_s
    !> Model constant for WALE (Wall-Adapting Local Eddy-Viscosity)
    !! default: 0.5
    real(kind=rk) :: C_w
    !> Model constant for Vreman model
    !! In literature: C_v = sqrt(2.5) C_s = 0.27
    real(kind=rk) :: C_v
    !> Modal constant for turbulent kinetic energy dissipation rate
    !! default: 1.048
    real(kind=rk) :: C_e
    !> Model constant for eddy-viscosity coefficient
    !! default: 0.094_rk
    !! https://caefn.com/openfoam/smagorinsky-sgs-model
    real(kind=rk) :: C_k
  end type les_coeff_type

  !> Contains turbulence information loaded from config file
  type mus_turbulence_config_type
    !> turbulence model type
    character(len=labelLen) :: model

    !> les model coefficients
    type(les_coeff_type) :: coeff

    !> To compute strain-rate from non-equilibrium PDF for Smagorinsky les model.
    !! If true then velocity and grad data are not required
    logical :: compSR_fromPDF = .false.

    !> Use Kolmogorov scale for interpolation turbulent viscosity for multilevel
    logical :: useKolmogorovScale = .true.
  end type mus_turbulence_config_type

  !> Contains velocity and gradient data to compute eddy viscosity
  type mus_turbulence_data_type
    !> Communication buffers to communicate visoscity field
    !! Local Fluids required by remote processes
    type( tem_communication_type ) :: sendBuffer
    !> My halos which are fluids on remote processes
    type( tem_communication_type ) :: recvBuffer
    !> Local ghostFromCoarser required by remote processes
    type( tem_communication_type ) :: sendBufferFromCoarser
    !> Local ghostFromFiner required by remote processes
    type( tem_communication_type ) :: sendBufferFromFiner
    !> My halos which are ghostFromCoarser on remote processes
    type( tem_communication_type ) :: recvBufferFromCoarser
    !> My halos which are ghostFromFiner on remote processes
    type( tem_communication_type ) :: recvBufferFromFiner

    !> Normalized turbulence viscosity
    !! i.e. viscosity scaled to current level i.e. visc/dtL
    !! Size: nSize (nFluids+nGhosts+nHalos)
    !! Used gradData to compute viscosity for nFluids and nGhosts.
    !! This viscosity is interpolated and scaled for setting nonEq term
    !! interpolation routines. The source element of interpolation might be
    !! halo so they are communicated.
    !!
    !! Simple scaling assumping norm of strainrate tensor |S| in different
    !! level is small:
    !! Schneider, A. (2015). A Consistent Large Eddy Approach for
    !! Lattice Boltzmann Methods and its Application to Complex Flows.
    !! Technical University Kaiserslautern.
    !! v_c = 4 v_f. Scaled visc: v^s = v/dt.
    !! => v^s_c dtL_c = 4 v^s_f dtL_f => v^s_c = 2 v^s_f
    !!
    !! Kolmogorov scaling:
    !! Touil, H., Ricot, D., & Lévêque, E. (2014). Direct and large-eddy
    !! simulation of turbulent flows on composite multi-resolution grids by
    !! the lattice Boltzmann method. Journal of Computational Physics, 256,
    !! 220–233.
    !! v^s_c = 2^(1/3) v^s_f
    real(kind=rk), allocatable :: visc(:)
  end type mus_turbulence_data_type


  !> Contains function pointers to obtain normalized turbulence viscosity.
  !! Viscosity is normalized to current level i.e. v_s = v/dt
  type mus_turbulence_visc_proc_type
    !> this procedure compute eddy viscosity from velocity field depending
    !! turbulence and lbm (compressible/incompressible) models
    procedure(proc_calc_turb_visc_fromGradU), pointer, nopass &
      & :: fromGradU => null()

    !> this procedure compute eddy viscosity from preCollision PDF.
    !! It is used for Smagorinsky model which depends only on strain rate
    !! that can be calculated using local nonEquilibrium.
    !! Is assigned when compSR_fromPDF is .true.
    procedure(proc_calc_turb_visc_fromPreColPDF), pointer, nopass &
      & :: fromPreColPDF => null()
  end type mus_turbulence_visc_proc_type

  !> Contains information required to compute eddy viscosity
  type mus_turbulence_type
    !> is true if turbulence table is defined
    logical :: active
    !> information loaded from config file
    type(mus_turbulence_config_type) :: config
    !> contains level-wise turbulence data to compute eddy viscosity
    !! size: minlevel:maxLevel
    type(mus_turbulence_data_type), allocatable :: dataOnLvl(:)

    !> contains turbulence viscosity function pointers
    type(mus_turbulence_visc_proc_type) :: calcVisc

    !> Factor to scale normalized turbulent viscosity from coarse to fine
    !! depending on whether useKolmogorovScale true or false
    !! if useKolmogorovScale fac_c2f = 1/2^(1/3) else fac_c2f = 1/2
    !! How to use: v^s_f = fac_c2f v^s_c
    real(kind=rk) :: fac_c2f

    !> Factor to scale normalized turbulent viscosity from fine to coarse
    !! depending on whether useKolmogorovScale true or false
    !! if useKolmogorovScale fac_f2c = 2^(1/3) else fac_f2c = 2
    !! How to use: v^s_c = fac_f2c v^s_f
    real(kind=rk) :: fac_f2c

  end type mus_turbulence_type

  !> interface to calculate subgrid scale turbulent eddy viscosity
  abstract interface
    !> This function computes turbulent viscosity from gradient U
    subroutine proc_calc_turb_visc_fromGradU(turbVisc, turbConfig, gradData, &
      & auxField, velPos, nSolve, nAuxScalars, dxL, dtL, Grad)
      import :: rk, mus_turbulence_config_type, mus_gradData_type, mus_Grad_type

      !> output: turbulent viscosity
      real(kind=rk), intent(out) :: turbVisc(:)

      !> turbulence config contains oefficients
      type(mus_turbulence_config_type), intent(in) :: turbConfig

      !> gradient data
      type(mus_gradData_type), intent(in) :: gradData

      !> Auxiliary field variable array
      real(kind=rk), intent(in) :: auxField(:)

      !> position of velocity components in auxField
      integer, intent(in) :: velPos(3)

      !> Number of element to solve in this level
      integer, intent(in) :: nSolve

      !> number of scalars in auxField array
      integer, intent(in) :: nAuxScalars

      !> turbulence coefficients
      !> current level lattice element size
      real(kind=rk), intent(in) :: dxL

      !> current level lattice time step size
      real(kind=rk), intent(in) :: dtL

      !> Object that contains pointers to calculate gradients
      type(mus_Grad_type), intent(in) :: Grad
    end subroutine proc_calc_turb_visc_fromGradU

    !> This function compute turbulent viscosity from pre-collision PDF
    subroutine proc_calc_turb_visc_fromPreColPDF(turbVisc, turbConfig, state,  &
      & neigh, auxField, densPos, velPos, nSize, nSolve, nScalars, nAuxScalars,&
      & layout, dxL, dtL, viscKine)
      import :: rk, mus_turbulence_config_type, mus_scheme_layout_type

      !> output: turbulent viscosity
      real(kind=rk), intent(out) :: turbVisc(:)

      !> turbulence type is implicitly passed to access turbulence coefficients
      type(mus_turbulence_config_type), intent(in) :: turbConfig

      !> state array
      real(kind=rk), intent(in) :: state(:)

      !> neigh array to obtain precollision pdf
      integer, intent(in) :: neigh(:)

      !> Auxiliary field variable array
      real(kind=rk), intent(in) :: auxField(:)

      !> position of density in auxField
      integer, intent(in) :: densPos

      !> position of velocity components in auxField
      integer, intent(in) :: velPos(3)

      !> number of elements in state array
      integer, intent(in) :: nSize

      !> Number of element to solve in this level
      integer, intent(in) :: nSolve

      !> number of scalars in state array
      integer, intent(in) :: nScalars

      !> number of scalars in auxField array
      integer, intent(in) :: nAuxScalars

      !> scheme layout
      type(mus_scheme_layout_type), intent(in) :: layout

      !> current level lattice element size
      real(kind=rk), intent(in) :: dxL

      !> current level lattice time step size
      real(kind=rk), intent(in) :: dtL

      !> Background kinematic viscosity divided by dtL
      real(kind=rk), intent(in) :: viscKine(:)
    end subroutine proc_calc_turb_visc_fromPreColPDF
  end interface

contains

  ! ************************************************************************** !
  !> load turbulence table
  subroutine mus_load_turbulence(me, conf, parent)
    !--------------------------------------------------------------------------
    !> fluid type
    type( mus_turbulence_type ),intent(out) :: me
    !> lua state
    type( flu_state )              :: conf
    !> parent handle
    integer, intent(in), optional  :: parent
    !--------------------------------------------------------------------------
    integer :: turb_table, iError
    !--------------------------------------------------------------------------

    ! if fluid informations in scheme table parentHandle /= 0
    if( present(parent) ) then
      call aot_table_open( L       = conf,        &
        &                  parent  = parent,      &
        &                  thandle = turb_table,  &
        &                  key     = 'turbulence' )
    else
      call aot_table_open( L=conf, thandle = turb_table, key = 'turbulence' )
    end if

    if ( turb_table == 0 ) then
      write(logUnit(1),*)'No turbulence table defined'
      me%config%compSR_fromPDF = .false.
      me%active = .false.
      return
    endif

    write(logUnit(1),*) 'Loading turbulence informations'
    me%active = .true.
    ! default is set to false for other models, its optional only for
    ! Smagorinsky model
    me%config%compSR_fromPDF = .false.

    ! load turbulence model
    call aot_get_val(L       = conf,            &
      &              thandle = turb_table,      &
      &              key     = 'model',         &
      &              val     = me%config%model, &
      &              default = 'smagorinsky',   &
      &              ErrCode = iError           )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving model :'
      if (btest(iError, aoterr_NonExistent)) then
        write(logUnit(1),*)'Variable not existent!'
        write(logUnit(1),*) 'Using default model: smagorinsky'
      end if

      if (btest(iError, aoterr_WrongType)) then
         write(logUnit(1),*)'Variable has wrong type!'
         write(logUnit(1),*)'STOPPING'
         call tem_abort()
      end if
    end if

    write(logUnit(1),*) '  model: ', trim(me%config%model)
    select case(upper_to_lower(trim(me%config%model)))
    case('smagorinsky')
      call aot_get_val(L       = conf,                         &
        &              thandle = turb_table,                   &
        &              key     = 'compute_strainrate_fromPDF', &
        &              val     = me%config%compSR_fromPDF,     &
        &              default = .true.,                       &
        &              ErrCode = iError                        )

      write(logUnit(1),*) '  use preCol PDF to calc strainrate: ', &
        &                 me%config%compSR_fromPDF
      ! use Smagorinsky constant when computing eddy viscosity from PDF
      ! else use C_k and C_e
      if (me%config%compSR_fromPDF) then
        call aot_get_val(L       = conf,                &
          &              thandle = turb_table,          &
          &              key     = 'c_s',               &
          &              val     = me%config%coeff%C_s, &
          &              default = 0.17_rk,             &
          &              ErrCode = iError               )
        write(logUnit(1),*) '  C_s: ', me%config%coeff%C_s

      else
        call aot_get_val(L       = conf,                &
          &              thandle = turb_table,          &
          &              key     = 'c_e',               &
          &              val     = me%config%coeff%C_e, &
          &              default = 1.048_rk,            &
          &              ErrCode = iError               )
        write(logUnit(1),*) '  C_e: ', me%config%coeff%C_e

        call aot_get_val(L       = conf,                &
          &              thandle = turb_table,          &
          &              key     = 'c_k',               &
          &              val     = me%config%coeff%C_k, &
          &              default = 0.094_rk,            &
          &              ErrCode = iError               )
        write(logUnit(1),*) '  C_k: ', me%config%coeff%C_k
      end if
    case('wale')
      call aot_get_val(L       = conf,                &
        &              thandle = turb_table,          &
        &              key     = 'c_w',               &
        &              val     = me%config%coeff%C_w, &
        &              default = 0.5_rk,              &
        &              ErrCode = iError               )


      write(logUnit(1),*) '  C_w: ', me%config%coeff%C_w
    case('vreman')
      ! load model coefficients.
      ! Vreman model constant is related to Smagorinsky constant: c_v = 2.5*c_s^2
      ! For c_s=0.17, c_v is approximately 0.07
      call aot_get_val(L       = conf,                &
        &              thandle = turb_table,          &
        &              key     = 'c_v',               &
        &              val     = me%config%coeff%C_v, &
        &              default = 0.07_rk,              &
        &              ErrCode = iError               )

      write(logUnit(1),*) '  C_v: ', me%config%coeff%C_v
    case default
        call tem_abort('Error: Unknown turbulence model')
    end select

    ! to use Kolmogorov scale for multilevel
    call aot_get_val(L       = conf,                         &
      &              thandle = turb_table,                   &
      &              key     = 'use_kolmogorov_scale',       &
      &              val     = me%config%useKolmogorovScale, &
      &              default = .true.,                       &
      &              ErrCode = iError                        )


    write(logUnit(1),*) '  Use Kolmogorov scale: ', me%config%useKolmogorovScale
    ! set scaling factor to convert turbulent viscosity from coarse to fine
    ! and vice versa
    if(me%config%useKolmogorovScale) then
      me%fac_f2c = 2.0_rk**(1.0_rk/3.0_rk)
      me%fac_c2f = 1.0_rk/me%fac_f2c
    else
      me%fac_f2c = 2.0_rk
      me%fac_c2f = 1.0_rk/me%fac_f2c
    end if

    call aot_table_close( L=conf, thandle=turb_table )
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_load_turbulence
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This initialize turbulence data type which includes velocity array
  !! and communication buffer
  subroutine mus_init_turbulenceData(me, & !turbConfig,
    & levelDesc, pattern, nSize)
    !---------------------------------------------------------------------------
    !> turbulence data
    type(mus_turbulence_data_type), intent(out) :: me
    !!> tubulence configuration
    !type(mus_turbulence_config_type), intent(in) :: turbConfig
    !> levelDesc to access communication buffers of state array
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> communication pattern
    type(tem_commPattern_type), intent(in)      :: pattern
    !> Number of elements in state array
    integer, intent(in) :: nSize
    !---------------------------------------------------------------------------
    ! allocate tubulence viscosity array
    allocate(me%visc(nSize))

    ! initialize turbulent viscosity with zeros to fill ghost elements.
    ! Actual turbulent viscosity are computed before stream-collision after
    ! calculating auxField in control_routine
    me%visc(:) = 0.0_rk

    ! allocate send and recv buffer to exchange turbulenc viscosity
    !if (.not. turbConfig%compSR_fromPDF) then

      ! initialize communication buffers for viscosity
      call init_commBuffer(buffer_visc  = me%sendBuffer,       &
        &                  buffer_State = levelDesc%sendBuffer )
      call init_commBuffer(buffer_visc  = me%recvBuffer,       &
        &                  buffer_State = levelDesc%recvBuffer )

      call init_commBuffer(buffer_visc  = me%sendBufferFromCoarser,      &
        &                  buffer_State = levelDesc%SendBufferFromCoarser)
      call init_commBuffer(buffer_visc  = me%recvBufferFromCoarser,      &
        &                  buffer_State = levelDesc%recvBufferFromCoarser)

      call init_commBuffer(buffer_visc  = me%sendBufferFromFiner,      &
        &                  buffer_State = levelDesc%sendBufferFromFiner)
      call init_commBuffer(buffer_visc  = me%recvBufferFromFiner,      &
        &                  buffer_State = levelDesc%recvBufferFromFiner)
    !end if

    contains

    ! ************************************************************************ !
    subroutine init_commBuffer(buffer_visc, buffer_state)
      !-------------------------------------------------------------------------
      !> communication buffer for velocity field
      type(tem_communication_type), intent(out) :: buffer_visc
      !> communication buffer of state array which is already initialized
      !! in tem_construction_module
      type(tem_communication_type), intent(in) :: buffer_state
      !-------------------------------------------------------------------------
      integer :: iProc
      !-------------------------------------------------------------------------
      ! copy information about target and source procs from pdf sendBuffer to
      ! velocity sendBuffer
      call tem_comm_init(buffer_visc, buffer_state%nProcs)
      buffer_visc%proc = buffer_state%proc
      buffer_visc%nElemsProc = buffer_state%nElemsProc
      do iProc = 1, buffer_visc%nProcs
        call pattern%initBuf_real( me    = buffer_visc%buf_real( iProc ),   &
          &                        pos   = buffer_state%elemPos(iProc)%val, &
          &                        nVals = buffer_state%nElemsProc(iProc)   )
      end do

    end subroutine init_commBuffer
    ! ************************************************************************ !

  end subroutine mus_init_turbulenceData
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine compute turbulence viscosity and stores in turbulence data type
  subroutine mus_turb_calcVisc( turbData, turbConfig, calcTurbVisc, state,    &
    &                           neigh, auxField, gradData, densPos, velPos,   &
    &                           nSize, nSolve, nScalars, nAuxScalars, layout, &
    &                           dxL, dtL, viscKine, Grad                      )
    ! --------------------------------------------------------------------------
    !> turbulence data type
    type(mus_turbulence_data_type), intent(inout) :: turbData
    !> turbulence configuration
    type(mus_turbulence_config_type), intent(in) :: turbConfig
    !> turbulence function
    type(mus_turbulence_visc_proc_type), intent(in) :: calcTurbVisc
    !> state array
    real(kind=rk), intent(in) :: state(:)
    !> neigh array to obtain precollision pdf
    integer, intent(in) :: neigh(:)
    !> Auxiliary field variable array
    real(kind=rk), intent(in) :: auxField(:)
    !> gradient data
    type(mus_gradData_type), intent(in) :: gradData
    !> position of density in auxField
    integer, intent(in) :: densPos
    !> position of velocity components in auxField
    integer, intent(in) :: velPos(3)
    !> number of elements in state array
    integer, intent(in) :: nSize
    !> Number of element to solve in this level
    integer, intent(in) :: nSolve
    !> number of scalars in state array
    integer, intent(in) :: nScalars
    !> number of scalars in auxField array
    integer, intent(in) :: nAuxScalars
    !> scheme layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> lattice element size in current level
    real(kind=rk), intent(in) :: dxL
    !> current level lattice time step size
    real(kind=rk), intent(in) :: dtL
    !> Background kinematic viscosity divided by dtL
    real(kind=rk), intent(in) :: viscKine(:)
    !> Object that contains pointers to calculate gradients
    type(mus_Grad_type), intent(in) :: Grad
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! calculate turbulent viscosity
    if (turbConfig%compSR_fromPDF) then
      call calcTurbVisc%fromPreColPDF( turbVisc    = turbData%visc(1:nSolve), &
        &                              turbConfig  = turbConfig,              &
        &                              state       = state,                   &
        &                              neigh       = neigh,                   &
        &                              auxField    = auxField,                &
        &                              densPos     = densPos,                 &
        &                              velPos      = velPos,                  &
        &                              nSize       = nSize,                   &
        &                              nSolve      = nSolve,                  &
        &                              nScalars    = nScalars,                &
        &                              nAuxScalars = nAuxScalars,             &
        &                              layout      = layout,                  &
        &                              dxL         = dxL,                     &
        &                              dtL         = dtL,                     &
        &                              viscKine    = viscKine                 )
    else
      call calcTurbVisc%fromGradU( turbVisc    = turbData%visc(1:nSolve), &
        &                          turbConfig  = turbConfig,              &
        &                          gradData    = gradData,                &
        &                          auxField    = auxField,                &
        &                          velPos      = velPos,                  &
        &                          nSolve      = nSolve,                  &
        &                          nAuxScalars = nAuxScalars,             &
        &                          dxL         = dxL,                     &
        &                          dtL         = dtL,                     &
        &                          Grad        = Grad                     )
    end if

  end subroutine mus_turb_calcVisc
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine update turbulent viscosity of boundary elements from RANS
  !! viscosity computed in turbulent_wall boundary.
  subroutine mus_turb_updateViscOfTurbWall(turbData, viscTurbWall, nElems_bnd, &
    &                                      elemPos)
    ! --------------------------------------------------------------------------
    !> turbulence data type
    type(mus_turbulence_data_type), intent(inout) :: turbData
    !> Turbulent viscosity on turbulent wall boundary computed in set boundary
    real(kind=rk), intent(in) :: viscTurbWall(:)
    !> Number of elements in turbulent_wall boundary
    integer, intent(in) :: nElems_bnd
    !> Position of boundary element in levelwise total list or state array
    integer, intent(in) :: elemPos(:)
    ! --------------------------------------------------------------------------
    integer :: iElem
    ! --------------------------------------------------------------------------
    do iElem = 1, nElems_bnd
      turbData%visc(elemPos(iElem)) = viscTurbWall(iElem)
    end do
  end subroutine mus_turb_updateViscOfTurbWall
  ! ************************************************************************** !

end module mus_turbulence_module
