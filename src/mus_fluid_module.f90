! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012-2014, 2016-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012, 2021 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Philipp Otte <otte@mathcces.rwth-aachen.de>
! Copyright (c) 2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2019, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> This module keeps all information about the fluid
!!
!! In this module, all infos about the fluid is collected.
!! The file type [[mus_fluid_type]] contains all relevant information.
!! This includes physical parameters such as the viscosity and reference
!! density. Also, LBM-specific parameters such as the relaxation rates are
!! defined.
!!
module mus_fluid_module
  use mpi
  ! include treelm modules
  use env_module,           only: rk
  use tem_tools_module,     only: tem_horizontalSpacer
  use tem_aux_module,       only: tem_abort
  use tem_param_module,     only: cs2, cs2inv, rho0, div2_3
  use tem_logging_module,   only: logUnit
  use tem_time_module,      only: tem_time_type
  use tem_stencil_module,   only: tem_stencilHeader_type
  use tem_general_module,   only: tem_general_type
  use tem_spacetime_fun_module, only: tem_load_spacetime
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_grow_array_module,    only: init

  ! include aotus modules
  use aotus_module,     only: flu_State,          &
    &                         aot_get_val,        &
    &                         aoterr_NonExistent, &
    &                         aoterr_Fatal
  use aot_table_module, only: aot_table_open, aot_table_close
  use aot_out_module,   only: aot_out_type, aot_out_val, aot_out_open_table,   &
    &                         aot_out_close_table

  ! include musubi modules
  use mus_physics_module,       only: mus_physics_type
  use mus_nonNewtonian_module,  only: mus_nNwtn_type, mus_nNwtn_load, &
    &                                 mus_assign_nNwtnVisc_ptr,       &
    &                                 mus_nNwtn_dump2outUnit, mus_nNwtn_save2lua
  use mus_pdf_module,            only: pdf_data_type
  use mus_turbulence_module,     only: mus_turbulence_type,  &
    &                                  mus_load_turbulence,  &
    &                                  mus_init_turbulenceData
  use mus_turb_viscosity_module, only: mus_assign_turbVisc_ptr
  use mus_scheme_header_module,  only: mus_scheme_header_type
  use mus_relaxationParam_module, only: mus_viscosity_type,              &
    &                                   mus_init_relaxParam,             &
    &                                   mus_calcOmegaFromVisc,           &
    &                                   mus_update_relaxParamFromViscSTfun
  use mus_mrtRelaxation_module,   only: mus_proc_mrt,                      &
    &                                   mus_assign_mrt_ptr
  use mus_cumulantInit_module,    only: cumulant_omega_check

  implicit none

  private

  public :: mus_fluid_type
  public :: mus_load_fluid
  public :: mus_fluid_save2lua
  public :: mus_init_fluid
  public :: mus_fluid_cleanup

  !> collection of properties of the fluid
  type mus_fluid_type

    logical :: active = .false.       !< is this object a fluid?

    !> Magic value for TRT collision model
    !! Lambda = ( 1/omega_+ - 0.5 ) * ( 1/omega_- - 0.5 )
    real(kind=rk) :: lambda = 0.25_rk

    !> level-wise bulk omegas, used as relaxation in mrt model
    !! allocated in mus_init_fluid
    real(kind=rk), allocatable :: omegaBulkLvl(:)
    !> Level wise bulk viscosity in lattice
    real(kind=rk), allocatable :: viscBulkLvl(:)

    !> Contains information for turbulence model
    type(mus_turbulence_type) :: turbulence

    !> nonNewtonian fluid parameter
    type(mus_nNwtn_type) :: nNwtn

    !> function pointer to get MRT diagonal relaxation matrix
    procedure(mus_proc_mrt), nopass, pointer :: mrtPtr => null()

    !> kinematic viscosity
    !! \todo KM: implement interpolation routine for constant viscosity
    type(mus_viscosity_type) :: viscKine

    real(kind=rk) :: viscBulk_phy  !< physical bulk viscosity

    real(kind=rk) :: force(3)        !< forcing term

    ! also works around ICE in Intel 19.1
    ! HRR
    real(kind=rk) :: HRR_sigma = 0.98_rk

    ! DRT
    real(kind=rk) :: DRT_tauN = 0.70_rk

    ! Cumulant extended
    real(kind=rk) :: omega_Cum(10)
    real(kind=rk) :: omega_Lim(3)

  end type mus_fluid_type


contains


! **************************************************************************** !
  !> Read in the fluid property from the LUA parameter file
  subroutine mus_load_fluid( me, conf, parent, minLevel, physics, schemeHeader )
    !--------------------------------------------------------------------------
    !> fluid type
    type( mus_fluid_type ),intent(out) :: me
    !> lua state
    type( flu_state )              :: conf
    !> parent handle
    integer, intent(in), optional  :: parent
    !> global pdf info
    integer, intent(in) :: minLevel
    !> physics type to convert physics to lattice unit or vice versa
    type( mus_physics_type ), intent(in) :: physics
    !> identifier of the scheme
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------
    integer :: fluid_table, ii
    integer :: iError
    integer :: vecError(3)
    real(kind=rk) :: omegaBulk_def ! default omegaBulk for incompressible model
    character(len=6) :: omega_name = "omega_"
    character(len=8) :: omega_string
    character(len=10) :: omega_lim_name = "omega_lim_"
    character(len=11) :: omega_lim_string
    real(kind=rk) :: omegaCum_def
    ! --------------------------------------------------------------------------


    ! if fluid informations in scheme table parentHandle /= 0
    if( present(parent) ) then
      call aot_table_open( L       = conf,                                     &
        &                  parent  = parent,                                   &
        &                  thandle = fluid_table,                              &
        &                  key     = 'fluid' )
    else
      call aot_table_open( L=conf, thandle = fluid_table, key = 'fluid' )
    end if

    if ( fluid_table == 0 ) then
      write(logUnit(1),*)'No fluid table defined'
      call tem_abort()
    endif

    write(logUnit(1),*) 'Loading fluid informations'
    ! Activating the current fluid object
    me%active = .true.

    ! load kinematic viscosity as spacetime function
    call tem_load_spacetime( me      = me%viscKine%STfun,      &
      &                      conf    = conf,                   &
      &                      parent  = fluid_table,            &
      &                      key     = 'kinematic_viscosity',  &
      &                      errCode = iError                  )

    if (iError /= 0) then
      call tem_abort('ERROR: Unable to load "kinematic_viscosity"' &
        &          //'from fluid table')
    end if

    ! Load bulk viscosity as scalar and convert it to LB
    call aot_get_val( L       = conf,             &
      &               thandle = fluid_table,      &
      &               key     = 'bulk_viscosity', &
      &               val     = me%viscBulk_phy,  &
      &               ErrCode = iError            )

    ! If not provided, use default value only for incompressible model
    if (btest(iError, aoterr_NonExistent)) then
      if (trim(schemeHeader%kind) == 'fluid') then
        if (trim(schemeHeader%layout) == 'd2q9') then
          omegaBulk_def = 1.63_rk
          write(logUnit(1),*) 'Warning: Using default omegaBulk=1.63'
        else if (trim(schemeHeader%layout) == 'd3q27') then
          omegaBulk_def = 1.54_rk
          write(logUnit(1),*) 'Warning: Using default omegaBulk=1.54'
        else
          call tem_abort('ERROR: "bulk_viscosity" is not provided in '// &
            &            'fluid table')
        end if
      else
        ! using default value from
        ! D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P.,
        ! & Luo, L.-S. (2002).
        ! Multiple-relaxation-time lattice Boltzmann models in three
        ! dimensions. Philosophical Transactions. Series A, Mathematical,
        ! Physical, and Engineering Sciences, 360(1792), 437–51.
        select case(trim(schemeHeader%layout))
        case ('d2q9')
          ! default value for d2q9 is from
          ! Chen, S., Peng, C., Teng, Y., Wang, L. P., & Zhang, K. (2016).
          ! Improving lattice Boltzmann simulation of moving particles in a
          ! viscous flow using local grid refinement. Computers and Fluids,
          omegaBulk_def = 1.63_rk
        case ('d3q15')
          omegaBulk_def = 1.60_rk
        case ('d3q27')
          omegaBulk_def = 1.54_rk
        case default
          omegaBulk_def = 1.19_rk
        end select

        write(logUnit(1),*) 'Warning: Using default omegaBulk:', omegaBulk_def
      end if

      if (trim(schemeHeader%layout) == 'd2q9') then
        !! "Theory of the lattice Boltzmann method: Dispersion, dissipation,
        !! isotropy, Galilean invariance, and stability", Pierre Lallemand and
        !! Li-Shi Luo, Phys. Rev. E 61, 2000.
        me%viscBulk_phy = ((8._rk - 12._rk * cs2) / 24._rk) &
          &             * (1.0_rk/omegaBulk_def - 0.5_rk )  &
          &             * physics%fac(minLevel)%visc
      else
        me%viscBulk_phy = ((5._rk - 9._rk * cs2) / 9._rk)  &
          &             * (1.0_rk/omegaBulk_def - 0.5_rk ) &
          &             * physics%fac(minLevel)%visc
      endif

    endif

    ! load forces
    call aot_get_val( L       = conf,                &
      &               thandle = fluid_table,         &
      &               key     = 'force',             &
      &               val     = me%force,            &
      &               default = [0._rk,0._rk,0._rk], &
      &               ErrCode = vecError             )

    ! convert force to lattice
    me%force = me%force / physics%fac(minLevel)%force

    select case(trim(schemeHeader%relaxation))

    case('trt')
      ! load lambda (magic value for trt model)
      call aot_get_val(L       = conf,        &
        &              thandle = fluid_table, &
        &              key     = 'lambda',    &
        &              val     = me%lambda,   &
        &              default = 0.25_rk,     &
        &              ErrCode = iError       )

    case('hrr_bgk', 'hrr_bgk_corrected')
      ! load sigma for HRR_bgk
      call aot_get_val(L       = conf,        &
        &              thandle = fluid_table, &
        &              key     = 'hrr_sigma',     &
        &              val     = me%HRR_sigma,&
        &              default = 0.98_rk,     &
        &              ErrCode = iError       )

      if (btest(iError, aoterr_fatal)) then
        call tem_abort('ERROR: wrong data type for sigma_HRR ' &
            &          // 'fluid table')
      end if

    case('prr_bgk_corrected')
      me%HRR_sigma = 0._rk

    case('rr_bgk_corrected')
      me%HRR_sigma = 1._rk

    case('drt_bgk')
      ! load tauN for DRT_bgk
      call aot_get_val(L       = conf,        &
        &              thandle = fluid_table, &
        &              key     = 'drt_taun',     &
        &              val     = me%DRT_tauN,&
        &              default = 0.70_rk,     &
        &              ErrCode = iError       )

      if (btest(iError, aoterr_fatal)) then
        call tem_abort('ERROR: wrong data type for DRT tauN ' &
            &          // 'fluid table')
      end if

    case('cumulant_extended', 'cumulant_extended_generic')
      ! load omega vec for cumulant and limiter
      me%omega_Cum = -1._rk

      call aot_get_val(L       = conf,            &
        &              thandle = fluid_table,     &
        &              key     = 'omega_2',       &
        &              val     = me%omega_Cum(2), &
        &              default = -1._rk,          &
        &              ErrCode = iError           )

      if (me%omega_Cum(2) < 0._rk) then
        write(logUnit(1),*) 'Cumulant_extended is active'
        write(logUnit(1),*) '  Warning: Using omega_2 = omegaBulk'
      end if

      if (btest(iError, aoterr_fatal)) then
        call tem_abort('ERROR: wrong data type for omega_2 '// &
            &            'fluid table')
      end if

      do ii = 3, 10
        if (ii < 10) then
          write (omega_string, "(A6,I1)") omega_name, ii
        else
          write (omega_string, "(A6,I2)") omega_name, ii
        endif

        if (ii < 6) then
          omegaCum_def = -1._rk
        else
          omegaCum_def = 1._rk
        endif

        call aot_get_val(L       = conf,               &
          &              thandle = fluid_table,        &
          &              key     = trim(omega_string), &
          &              val     = me%omega_Cum(ii),   &
          &              default = omegaCum_def,       &
          &              ErrCode = iError              )

        if (btest(iError, aoterr_fatal)) then
          call tem_abort('ERROR: wrong data type for ' // trim(omega_string) &
              & // ' cumulant fluid table')
        end if

      enddo

      do ii = 1, 3
        write (omega_lim_string, "(A10,I1)") omega_lim_name, ii

        call aot_get_val(L       = conf,                       &
          &              thandle = fluid_table,                &
          &              key     = trim(omega_lim_string),     &
          &              val     = me%omega_Lim(ii),           &
          &              default = 1.e-2_rk,                   &
          &              ErrCode = iError                      )

        if (btest(iError, aoterr_fatal)) then
          call tem_abort('ERROR: wrong data type for ' &
            &            // trim(omega_lim_string) &
            &            // ' cumulant fluid table')
        end if

      enddo

    end select

    ! load nonNewtonian fluid feature
    call mus_nNwtn_load( me     = me%nNwtn,   &
      &                  conf   = conf,       &
      &                  parent = fluid_table )

    ! load turbulence
    call mus_load_turbulence( me     = me%turbulence, &
      &                       conf   = conf,          &
      &                       parent = fluid_table    )

    call aot_table_close( L=conf, thandle=fluid_table )
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_load_fluid
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This routine initilizes fluid visocity and relaxation paramters for each
  !! level
  subroutine mus_init_fluid(me, physics, schemeHeader, minLevel, maxLevel,    &
    &                       levelDesc, pdf, stencil, general, tNow)
    ! --------------------------------------------------------------------------
    !> fluid type
    type(mus_fluid_type), intent(inout) :: me
    !> physics type to convert physics to lattice unit or vice versa
    type(mus_physics_type), intent(in) :: physics
    !> scheme header
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    !> min and max level
    integer, intent(in) :: minLevel, maxLevel
    !> level descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc(minLevel:maxLevel)
    !> pdf info with neigh array for all levels
    type(pdf_data_type), intent(in) :: pdf(minLevel:maxLevel)
    !> stencil header
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> general type contains communication pattern and proc info
    type(tem_general_type), intent(in)      :: general
    !> current simulation time
    type(tem_time_type),intent(in) :: tNow
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------
    write(logUnit(1),'(A)') 'Initialize fluid'
    ! allocate array to store kinematic viscosity and bulk viscosity
    allocate(me%viscKine%dataOnLvl(minLevel:maxLevel))
    do iLevel = minLevel, maxLevel
      call init(me%viscKine%dataOnLvl(iLevel), pdf(iLevel)%nElems_local)
      me%viscKine%dataOnLvl(iLevel)%val = -100.0_rk
    end do
    allocate(me%viscBulkLvl(minLevel:maxLevel))
    allocate(me%omegaBulkLvl(minLevel:maxLevel))

    ! New initialization
    ! allocate relaxation parameter array to size nSolve
    call mus_init_relaxParam(omLvl    = me%viscKine%omLvl,  &
      &                      minLevel = minLevel,           &
      &                      maxLevel = maxLevel,           &
      &                      nElems   = pdf(:)%nElems_local )

    ! assign function pointer
    call mus_assign_mrt_ptr(me%mrtPtr, schemeHeader)

    ! intialize kinematic omega
    write(logUnit(3),'(A)') '  Set relaxation parameter from kinematic ' &
      &                     //'viscosity'
    do iLevel = minLevel, maxLevel
      call mus_update_relaxParamFromViscSTfun(             &
        & omega       = me%viscKine%omLvl(iLevel)%val,     &
        & visc        = me%viscKine%dataOnLvl(iLevel)%val, &
        & viscSTfun   = me%viscKine%STfun,                 &
        & viscRef     = physics%fac(iLevel)%visc,          &
        & nElems      = pdf(iLevel)%nElems_local,          &
        & baryOfTotal = levelDesc(iLevel)%baryOfTotal,     &
        & tNow        = tNow                               )

      ! Initialize bulk omega
      me%viscBulkLvl(iLevel) = me%viscBulk_phy / physics%fac(iLevel)%visc
      ! for 2D the formula is different!!!!!!
      !write(logUnit(3),'(A,I0)') '   nDims: ', stencil%nDims
      if (stencil%nDims == 2) then
        ! for 2d the formula of omega from bulk is different!
        ! see ref: "Theory of the Lattice Boltzmann Method: Dispersion,
        ! Dissipation, Isotropy, Galilean Invariance, and Stability"
        ! Luo 2000, PHYSICAL REVIEW E
        me%omegaBulkLvl(iLevel) = 1.0_rk / (24.0_rk * me%viscBulkLvl(iLevel) &
          & / (8._rk - 12._rk * cs2) + 0.5_rk )
      else
        ! see ref: "A D3Q27 multiple-relaxation-time lattice Boltzmann method
        ! for turbulent flows", Suga 2015, Computers&Fluids
        !me%omegaBulkLvl(iLevel) = 2.0_rk                                    &
        !  &                     / (9.0_rk * me%viscBulkLvl(iLevel) + 1.0_rk )
        me%omegaBulkLvl(iLevel) = 1.0_rk / (9.0_rk * me%viscBulkLvl(iLevel) &
          & / (5._rk - 9._rk * cs2) + 0.5_rk )
      endif

    end do

    ! Assign function pointer for nonNewtonian model
    if (me%nNwtn%active) then
      write(logUnit(3),'(A)') '  Assign function to compute viscosity ' &
        &                     //'for non-Newtonian model'
      call mus_assign_nNwtnVisc_ptr(me%nNwtn, schemeHeader)
    end if

    if (me%turbulence%active) then
      write(logUnit(3),'(A)') '  Assign function to compute viscosity ' &
        &                     //'for turbulence model'
      ! assign function pointer to compute turbulence viscosity
      call mus_assign_turbVisc_ptr(me%turbulence, schemeHeader)

      ! Initialize communication buffer for viscosity to compute
      ! velocity gradient for turbulence model
      ! This step must be done after construct_connectivity and
      ! init_levelBuffers
      allocate( me%turbulence%dataOnLvl( minLevel:maxLevel ) )
      do iLevel = minLevel, maxLevel
        call mus_init_turbulenceData(                     &
          & me         = me%turbulence%dataOnLvl(iLevel), &
          !& turbConfig = me%turbulence%config,            &
          & levelDesc  = levelDesc(iLevel),               &
          & pattern    = general%commPattern,             &
          & nSize      = pdf(iLevel)%nSize                )
      end do !iLevel

    end if

    ! Dump fluid information
    call mus_fluid_dump( me, minLevel, maxLevel, physics, pdf(:)%nElems_solve, &
      &                  general, logUnit(1), schemeHeader )

  end subroutine mus_init_fluid
  ! ************************************************************************** !

  ! ************************************************************************** !
  subroutine mus_fluid_dump( me, minLevel, maxLevel, physics, nSolve, general, &
    &                        outUnit, schemeHeader )
    ! --------------------------------------------------------------------------
    type( mus_fluid_type ), intent(inout) :: me !< fluid type
    !> minlevel and maxlevel
    integer, intent(in) :: minLevel, maxLevel
    !> physics type to convert physics to lattice unit or vice versa
    type( mus_physics_type ), intent(in) :: physics
    !> number of elements to solve per level (fluid+ghost)
    integer, intent(in) :: nSolve(minLevel:maxLevel)
    !> general type
    type(tem_general_type), intent(in) :: general
    integer, intent(in) :: outUnit
    !> scheme header
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------
    integer :: iLevel, iError
    real(kind=rk) :: visc_min, visc_max, viscRef, om_min, om_max
    real(kind=rk) :: glob_visc_min, glob_visc_max
    ! --------------------------------------------------------------------------

    write(outUnit,"(A)")       'Fluid properties:'
    write(outUnit,"(A,F10.7)") '  Reference density rho0: ', rho0
    write(outUnit,"(A)") '  Kinematic viscosity STfun kind: ' // &
      &                        trim(me%viscKine%STfun%fun_kind)

    do iLevel = minLevel, maxLevel
      ! get global min and max viscosity
      visc_min = minval(me%viscKine%dataOnLvl(iLevel)%val(1:nSolve(iLevel)))
      visc_max = maxval(me%viscKine%dataOnLvl(iLevel)%val(1:nSolve(iLevel)))

      call mpi_reduce( visc_min, glob_visc_min, 1, mpi_double_precision, &
                       mpi_min, 0, general%proc%comm, ierror    )
      call mpi_reduce( visc_max, glob_visc_max, 1, mpi_double_precision, &
                       mpi_max, 0, general%proc%comm, ierror    )
      viscRef = physics%fac(iLevel)%visc
      om_min = mus_calcOmegaFromVisc(glob_visc_min)
      om_max = mus_calcOmegaFromVisc(glob_visc_max)
      write(outUnit,'(A,I0)') 'Kinematic viscosity, (min, max) on level: ', &
        &                        iLevel
      write(outUnit,'(A,F10.7,A,F10.7,A)') '  Physical: (', &
        & glob_visc_min*viscRef, ',', glob_visc_max*viscRef, ')'
      write(outUnit,'(A,F10.7,A,F10.7,A)') '  Lattice: (', glob_visc_min, &
        &                                     ',', glob_visc_max, ')'
      write(outUnit,'(A,F10.7,A,F10.7,A)') '  Kinematic omega (min,max):(', &
        &                                     om_min,',', om_max, ')'

      write(outUnit,"(A)") ' viscBulk lattice and omega: '
      write(outUnit,"(A,F10.7)") '  viscBulk:  ', me%viscBulkLvl(iLevel)
      write(outUnit,"(A,F10.7)") '  omegaBulk: ', me%omegaBulkLvl(iLevel)

    end do

    write(outUnit, "(A,F10.7)") 'Magic lambda for TRT: ', me%lambda

    write(outUnit, "(A,F10.7)") 'HRR sigma: ', me%HRR_sigma

    write(outUnit, "(A,F10.7)") 'DRT tauN: ', me%DRT_tauN

    if (trim(schemeHeader%relaxation) == 'cumulant_extended' .or. &
      & trim(schemeHeader%relaxation) == 'cumulant_extended_generic') then
      write(outUnit, "(A)") 'omega Cumulant: -1.000 means adjust during runtime'
      do iLevel = 1, 10
        write(outUnit, "(A,I2,A,F10.7)") '  ', iLevel,' = ', me%omega_Cum(iLevel)
      enddo
      write(outUnit, "(A)") 'omega limiter Cumulant: 1.e10 means unlimited'
      do iLevel = 1, 3
        write(outUnit, "(A,I2,A,F10.7)") '  ', iLevel,' = ', me%omega_Lim(iLevel)
      enddo

      ! check Cumulant_extended range of admissibility for omegas
      write(logUnit(1),*) 'Checking admissibility of omegas for parametrized Cumulant scheme'

      do iLevel = minLevel, maxLevel
        call cumulant_omega_check( omegaVisc = me%viscKine%omLvl(iLevel)%val(:), &
          &                        omegaBulk = me%omegaBulkLvl(iLevel),          &
          &                        omegaIn   = me%omega_Cum(:),                  &
          &                        nSolve    = nSolve(iLevel),                   &
          &                        level     = iLevel                            )
      enddo
    end if

    if( maxval( abs( me%force )) > 0._rk )  &
      & write(outUnit,*) ' Forcing:             ', real(me%force)

    ! Dump nonNewtonian parameters to outUnit
    call mus_nNwtn_dump2outUnit( me%nNwtn, outUnit )

  end subroutine mus_fluid_dump
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> write fluid prop into a lua file
  !!
  subroutine mus_fluid_save2lua( me, conf )
    ! ---------------------------------------------------------------------------
    !> single fluid type
    type( mus_fluid_type ), intent(in) :: me
    type( aot_out_type ) :: conf
    ! ---------------------------------------------------------------------------

    call aot_out_open_table( put_conf = conf, tname = 'fluid' )

    call aot_out_val( put_conf = conf,             &
      &               vname    = 'bulk_viscosity', &
      &               val      = me%viscBulk_phy   )
    if (trim(me%viscKine%STfun%fun_kind)=='const') then
      call aot_out_val( put_conf = conf,                      &
        &               vname    = 'kine_viscosity',          &
        &               val      = me%viscKine%STfun%const(1) )
    end if

    if ( me%nNwtn%active ) call mus_nNwtn_save2lua( me%nNwtn, conf )

    call aot_out_close_table( put_conf = conf )

  end subroutine mus_fluid_save2lua
  ! **************************************************************************** !

  ! ************************************************************************** !
  !> This routines act as a destructor for fluid type
  subroutine mus_fluid_cleanup( me )
    ! ---------------------------------------------------------------------------
    !> single fluid type
    type( mus_fluid_type ), intent(inout) :: me
    ! ---------------------------------------------------------------------------

    if (allocated(me%viscKine%dataOnLvl)) then
      deallocate(me%viscKine%dataOnLvl)
    end if
    if (allocated(me%viscKine%omLvl)) then
      deallocate(me%viscKine%omLvl)
    end if
    if (allocated(me%viscBulkLvl)) then
      deallocate(me%viscBulkLvl)
    end if
    if (allocated(me%omegaBulkLvl)) then
      deallocate(me%omegaBulkLvl)
    end if
    if (me%turbulence%active) then
      if (allocated(me%turbulence%dataOnLvl)) then
        deallocate(me%turbulence%dataOnLvl)
      end if
    end if

  end subroutine mus_fluid_cleanup
  ! ************************************************************************** !

end module mus_fluid_module
! ****************************************************************************** !
