! Copyright (c) 2013-2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
!> This module contains all information about fluid mixture
module mus_mixture_module
  use env_module,               only: rk, globalMaxLevels, single_k,           &
    &                                 eps_single_k, labelLen
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_aux_module,           only: tem_abort
  use tem_param_module,         only: cs2, cs2inv
  use tem_ini_condition_module, only: tem_ini_condition_type, tem_load_ic
  use tem_spatial_module,       only: tem_spatial_type, tem_load_spatial
  use tem_temporal_module,      only: tem_temporal_type, tem_load_temporal
  use tem_spacetime_fun_module, only: tem_spacetime_fun_type
  use tem_logging_module,       only: logUnit

  use aot_table_module, only: aot_table_open, aot_table_close
  use aotus_module,     only: flu_State, aot_get_val,                          &
    &                         aoterr_Fatal, aoterr_NonExistent,                &
    &                         aoterr_WrongType, flu_State
  use aot_out_module,   only: aot_out_type, aot_out_val, aot_out_open_table,   &
    &                         aot_out_close_table

  use mus_physics_module,       only: mus_physics_type
  use mus_scheme_header_module, only: mus_scheme_header_type
  use mus_eNRTL_module,         only: mus_init_eNRTL

  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char

  implicit none

  private

  public :: mus_mixture_type
  public :: mus_load_mixture
  public :: mus_mixture_out

  !> relaxation paraemters multispecies required for each level
  type mixRelaxation_type
    !> lattice shear viscosity of the mixture for different level
    real(kind=rk) :: visc
    !> lattice free parameter B for different level
    real(kind=rk) :: paramB
    !> lattice bulk viscosity for different level
    real(kind=rk) :: bulkVisc
    !> relaxation parameter,
    !! \( \omega = K*B/\rho \), B - free parameter unit same as resistivity
    real(kind=rk) :: omega_diff
    real(kind=rk) :: omega_kine !< kinematic viscosity relaxation parameter
  end type mixRelaxation_type

  !> This type contains mixture information
  !! @todo KM: implement parameters which depend on dx and dt for all levels
  type mus_mixture_type
    !> initialization case, initial condition of the mixture
    type( tem_ini_condition_type ) :: ic
    !> mass density of the mixture
    !! SI unit: kg/m^3.
    !! Physics to lattice conversion: rho0/physics%rho0
    real(kind=rk) :: rho0
    !> lattice mass density of the mixture
    real(kind=rk) :: rho0LB
    !> number density of the mixture or total mixture molar density
    !! SI unit: mol/m^3.
    !! physics to lattice conversion:
    !! nT0/physics%mol*physics%fac(minlevel)%length^3.
    !! mixture molar density is required if initial condition
    !! are defined by molar fraction.
    !! set only if initial molefraction is space independent
    real(kind=rk) :: moleDens0
    !> lattice number density of the mixture or total mixture molar density
    !! physics to lattice conversion:
    !! moleDens0/physics%mol*physics%fac(minlevel)%length^3.
    real(kind=rk) :: moleDens0LB
    !> lattice kinematic shear viscosity of the mixture
    real(kind=rk) :: kine_viscosityLB
    !> physical kinematic shear viscosity of the mixture
    real(kind=rk) :: kine_viscosity
    !> lattice bulk viscosity
    real(kind=rk) :: bulk_viscosityLB
    !> physical bulk viscosity
    real(kind=rk) :: bulk_viscosity
    !> lattice bulk modulus of the liquid mixture
    !! \( K = c_s^2*\rho \),
    !! \( c_s \) - speed of sound (in lattice unit: \f$ 1/\sqrt{3} \f$
    real(kind=rk) :: bulk_modulusLB
    !> relaxation parameter,
    !! \( \omega = K*B/\rho \), B - free parameter unit same as resistivity
    real(kind=rk) :: omega_diff
    real(kind=rk) :: omega_kine !< kinematic viscosity relaxation parameter
    !> omega for higher order moments
    real(kind=rk) :: omega_hom
    !> relaxation parameters for each level
    type(mixRelaxation_type) :: relaxLvl(globalMaxLevels)
    !> temporal omega for ramping etc.
    type(tem_temporal_type) :: omega_ramping
    !> free parameter
    !! \f$ B = \omega*\rho/K \f$
    real(kind=rk) :: paramB
    !> temperature
    real(kind=rk) :: temp0
    !> temperature
    real(kind=rk) :: temp0LB
    !> equilibrium theta to choose between mixture velocity
    !! and equilibrium species velocity in the quadratic term
    !! equilibrium function.
    !! theta = 0 -> mixture velocity
    !! theta = 1 -> equilibrium species velocity
    !! \todo KM: remove theta_eq and use
    !! mixture velocity in quadratic term of equilibrium function
    real(kind=rk) :: theta_eq
    !> spatial omega definition, e.g. for sponge layers
    type(tem_spatial_type) :: viscSpatial
    !> external electrical force
    !@todo use source term to define external force
    real(kind=rk) :: electricField(3)
    !> gravitational force
    real(kind=rk) :: gravityField(3)
    !> faraday constant(C/mol)
    real(kind=rk) :: faraday
    !> faraday constant in lattice
    real(kind=rk) :: faradayLB
    !> gas constant R (Nm/(mol *K))
    real(kind=rk) :: gasConst_R
    !> gas constant R in lattice
    real(kind=rk) :: gasConst_R_LB
    !> eNRTL file with species properties
    character(kind=c_char, len=labelLen) :: prop_file
    !> atmospheric pressure
    real(kind=rk) :: atm_press
    !> atmospheric pressure
    real(kind=rk) :: atm_pressLB
  end type mus_mixture_type


contains


! **************************************************************************** !
  !> This routine load mixture table from scheme table.
  !! Define either mass density or number density.
  !! If mass density is specified, number density can be computed at runtime
  !! or vice versa.
  !! KM: @todo Currently, the simulation is initialized by density, extend
  !! it to initialize from mixture number density/volume fraction
  !! and mole fraction
  !! \verbatim
  !! mixture = { rho0 = 1.0, omega }
  !! \endverbatim
  subroutine mus_load_mixture( me, conf, parent, minLevel, maxLevel, physics,  &
    &                          schemeHeader, nFields )
    ! --------------------------------------------------------------------------
    !> contains mixture information
    type( mus_mixture_type ), intent(out) :: me
    type( flu_State ) :: conf !< flu state
    integer, intent(in), optional :: parent !< parent lua handle
    !> identifier of the scheme
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    integer, intent(in) :: minLevel, maxLevel
    !> physics type to convert physics to lattice unit or vice versa
    type( mus_physics_type ), intent(in) :: physics
    !> number of fields defined in lua file
    integer, intent(in) :: nFields
    ! --------------------------------------------------------------------------
    !local variables
    integer :: mix_handle
    integer :: iError(1)
    integer :: vError(3), errFatal(3)
    logical :: prop_read
    integer :: nFields_loc
    ! --------------------------------------------------------------------------
    errFatal = aotErr_Fatal

    call aot_table_open( L = conf, parent = parent, thandle = mix_handle,    &
      &                  key = 'mixture')

    ! if mixture handle is not defined
    if ( mix_handle == 0 ) then
      write(logUnit(1),*)' ERROR: No mixture table defined'
      write(logUnit(1),*)'        Mixture properties are neccessary for'
      write(logUnit(1),*)'        multispecies scheme'
      call tem_abort()
    endif

    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*)' Loading mixture information'

    ! load initial condition for mixture
    call tem_load_ic( me        = me%ic,               &
      &               conf      = conf,                &
      &               parent    = mix_handle,          &
      &               key       = 'initial_condition', &
      &               errCode   = iError,              &
      &               StateName = ['pressure']         )

    ! load species property file for thermodynamic factor only if
    ! relaxation kind is bgk_withthermodynfac
    if (trim(schemeHeader%relaxation) == 'bgk_withthermodynfac' .or. &
      & trim(schemeHeader%relaxation) == 'mrt_withthermodynfac') then
      ! load species property from eNRTL file
      call aot_get_val( L = conf, thandle = mix_handle, key = 'prop_file',     &
        &               val = me%prop_file, default = '', ErrCode=iError(1) )

      ! if property file is defined, initialize eNRTL parameters using
      ! c-code
      if(trim(me%prop_file) /= '') then
        me%prop_file = trim(me%prop_file)//C_NULL_CHAR
        write(logUnit(1),*) 'Thermodynamic property file name ' &
          &                 // trim(me%prop_file)
        prop_read = mus_init_eNRTL( me%prop_file, nFields_loc )
        if (.not. prop_read) then
          write(logUnit(1),*) 'ERROR: loading prop_file need to compute ' &
            &                 // 'thermodynamic factors'
          write(logUnit(1),*) 'Solution: Compile musubi with "--with_ext_tdf"' &
            &                 //'to compile musubi with external c-code which'
          write(logUnit(1),*) '          computes thermodynamic factor'
          call tem_abort()
        else
          if (nFields_loc /= nFields) then
            write(logUnit(1),*) 'nFields in config file: ', nFields
            write(logUnit(1),*) 'nFields in property file: ', nFields_loc
            call tem_abort("Error: nFields in config file /= "&
              &          //"nFields in thermodynamic property file")
          end if
        endif
      endif
    endif

    !get mixture density
    call aot_get_val( L = conf, thandle = mix_handle, key = 'rho0',            &
      &               val = me%rho0, default = 1.0_rk, ErrCode=iError(1) )

    if (btest(iError(1), aoterr_Fatal)) then
      write(logUnit(1),*) ' FATAL Error occured, while retrieving mixture '//  &
        &             'density :'
      if (btest(iError(1), aoterr_WrongType)) then
        write(logUnit(1),*)' Variable has wrong type!'
        call tem_abort()
      end if
    end if

    !convert to lattice
    me%rho0LB = me%rho0/physics%rho0

    !bulk modulus K = cs2 * rho
    me%bulk_modulusLB = cs2 * me%rho0

    !get mixture number density
    call aot_get_val( L = conf, thandle = mix_handle, key = 'moleDens0',       &
      &               val = me%moleDens0, default = 1.0_rk, ErrCode=iError(1) )
    !KM: @todo check whether we need moleDens for different level
    me%moleDens0LB = me%moleDens0 / physics%moleDens0

    if (btest(iError(1), aoterr_Fatal)) then
      write(logUnit(1),*) ' FATAL Error occured, while retrieving mixture '//  &
        &             'density :'
      if (btest(iError(1), aoterr_WrongType)) then
        write(logUnit(1),*)' Variable has wrong type!'
        call tem_abort()
      end if
    end if

    !get equilibrium theta
    call aot_get_val( L = conf, thandle = mix_handle, key = 'theta_eq',        &
      &               val = me%theta_eq, ErrCode=iError(1) )
    if (btest(iError(1), aoterr_Fatal)) then
      write(logUnit(1),*) ' FATAL Error occured, while retrieving theta_eq '
      if (btest(iError(1), aoterr_NonExistent)) then
        write(logUnit(1),*)' ATTENTION: Setting theta_eq=1.0 i.e using mass '
        write(logUnit(1),*)'averaged mixture velocity in quadratic part of '
        write(logUnit(1),*)'equilibrium function'
        me%theta_eq = 1.0_rk
      end if
      if (btest(iError(1), aoterr_WrongType)) then
        write(logUnit(1),*)' Variable has wrong type!'
        call tem_abort()
      end if
    end if

    ! Omega ramping
    call tem_load_temporal( me     = me%omega_ramping, &
      &                     conf   = conf,             &
      &                     parent = mix_handle,       &
      &                     key    = 'omega_ramping'   )

    me%viscSpatial%isStored = .false.
    call tem_load_spatial( me     = me%viscSpatial, &
      &                    conf   = conf,           &
      &                    parent = mix_handle,     &
      &                    key    = 'visc_spatial'  )

    !get diffusivity relaxation parameter
    call aot_get_val(L = conf, thandle = mix_handle, key = 'omega_diff',       &
      &              val = me%omega_diff, default = 2.0_rk, ErrCode = iError(1))
    if (btest(iError(1), aoterr_Fatal)) then
      write(logUnit(1),*) ' FATAL Error occured, while retrieving relaxation ' &
        &                 // 'parameter omega_diff :'
      if (btest(iError(1), aoterr_WrongType)) then
        write(logUnit(1),*)' Variable has wrong type!'
        call tem_abort()
      end if
    end if

    ! if omega is not provided. load free parameter B and compute omega
    ! from paramB
    if ( btest( iError(1), aotErr_NonExistent )) then
      write(logUnit(1),*) 'omega_diff not defined. Load paramB '
      call aot_get_val(L = conf, thandle = mix_handle, key = 'paramB',      &
        &              val = me%paramB, ErrCode = iError(1))
      if (btest(iError(1), aoterr_NonExistent)) then
        write(logUnit(1),*) 'ATTENTION: neither omega_diff nor paramB'
        write(logUnit(1),*) '           Setting default value to     '
        write(logUnit(1),*) '           omega_diff = 2.0'
        me%omega_diff = 2.0
        me%paramB = me%omega_diff / cs2
      else
        !compute omega from paramB
        ! omega_diff = cs2 * B
        ! cs2 = p/rho for multispecies gas
        ! cs2 = K/rho for multispecies liquid
        me%omega_diff = me%paramB * cs2
      end if
    else
      me%paramB = me%omega_diff / cs2
    endif

    !get kinematic shear viscosity relaxation parameter
    call aot_get_val(L = conf, thandle = mix_handle, key = 'omega_kine',       &
      &              val = me%omega_kine, ErrCode = iError(1))
    if (btest(iError(1), aoterr_Fatal)) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving relaxation '//&
        &             'parameter omega_kine:'
      if (btest(iError(1), aoterr_WrongType)) then
        write(logUnit(1),*)'Variable has wrong type!'
        call tem_abort()
      end if
    end if
    ! if omega_kine is not provided. load kinematic viscosity in physical unit
    ! if physics table is active
    ! from paramB
    if ( btest( iError(1), aotErr_NonExistent )) then
      write(logUnit(1),*) 'omega_kine is not defined.'//                       &
        &                 'Load kinematic shear viscosity'
      call aot_get_val(L = conf, thandle = mix_handle,             &
        &              key = 'kinematic_viscosity',                &
        &              val = me%kine_viscosity, ErrCode = iError(1))
      if (btest(iError(1), aoterr_NonExistent)) then
        write(logUnit(1),*) 'ATTENTION: neither omega_kine nor'//              &
          &             'kine_shear_viscosity is defined'
        write(logUnit(1),*) '           Setting default value to lattice'//    &
          &                 'kine_shear_vis c= 1/(3*omega_diff)'
        me%kine_viscosityLB = cs2 / me%omega_diff
        me%kine_viscosity = me%kine_viscosityLB &
          &               * physics%fac(minLevel)%visc
        me%omega_kine = me%omega_diff
      else
        ! convert physical visc to lattice
        me%kine_viscosityLB = me%kine_viscosity                                &
          &                 / physics%fac(minLevel)%visc
        !compute omega_kine from kinematic shear viscosity
        ! omega_kine = cs2 / viscosity
        me%omega_kine = cs2 / me%kine_viscosityLB
      end if
    else
    ! omega = cs^2/nu -> nu = cs^2/omega
      me%kine_viscosityLB = cs2 / me%omega_kine
      !physical
      me%kine_viscosity = me%kine_viscosityLB*physics%fac(minLevel)%visc
    endif

    call aot_get_val(L = conf, thandle = mix_handle, key = 'bulk_viscosity',   &
      &              val = me%bulk_viscosity,                                  &
      &              ErrCode = iError(1))

    !convert to lattice
    me%bulk_viscosityLB = me%bulk_viscosity/physics%fac(minLevel)%visc

    if (btest(iError(1), aoterr_NonExistent)) then
      ! formula for bulk viscosity.
      ! http://scienceworld.wolfram.com/physics/BulkViscosity.html
      ! require second viscosity coeff which is missing here
      write(logUnit(1),*)' ATTENTION: Bulk viscosity is not provided.'
      write(logUnit(1),*)'            It is computed from shear viscosity as ' &
        &                //'bulk_visc = 2*nu/3'
      me%bulk_viscosityLB = 2._rk/3._rk*me%kine_viscosityLB
      !me%bulk_viscosityLB = me%kine_viscosityLB
      me%bulk_viscosity = me%bulk_viscosityLB                                  &
        &               * physics%fac(minLevel)%visc
    endif

    !get relaxation parameter for higher order moments
    call aot_get_val(L = conf, thandle = mix_handle, key = 'omega_hom',        &
      &              val = me%omega_hom, default = 2.0_rk, ErrCode = iError(1))

    !! load paramters for forces
    !electric field
    call aot_get_val( L = conf, thandle = mix_handle, key = 'electricField',   &
      &               val = me%electricField, default=(/0.0_rk,0.0_rk,0.0_rk/),&
      &               ErrCode=vError )

    if ( any(btest(vError, errFatal))) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving electricField '
      call tem_abort()
    end if

    ! convert electricfield into lattice unit (N/C)
    me%electricField = me%electricField * physics%coulomb0                     &
      &              / physics%fac(minLevel)%force

    !gravity field
    call aot_get_val( L = conf, thandle = mix_handle, key = 'gravityField',    &
      &               val = me%gravityField, default=(/0.0_rk,0.0_rk,0.0_rk/), &
      &               ErrCode=vError )

    if ( any(btest(vError, errFatal))) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving gravityField '
      call tem_abort()
    end if

    ! convert gravityfield into lattice unit (m/s^2)
    me%gravityField = me%gravityField                                          &
      &              * physics%dtLvl(minLevel)**2                    &
      &              / physics%dxLvl(minLevel)

    !temperature
    call aot_get_val( L = conf, thandle = mix_handle, key = 'temp',            &
      &               val = me%temp0, default=273._rk, ErrCode=iError(1) )

    if ( btest(iError(1), aoterr_Fatal)) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving ' //          &
        &                 'temp (temperature)'
      call tem_abort()
    end if
    me%temp0LB = me%temp0/physics%temp0

    ! atmospheric pressure
    call aot_get_val( L = conf, thandle = mix_handle, key = 'atm_press', &
      &               val = me%atm_press, default=1.01325e5_rk,          &
      &               ErrCode=iError(1)                                  )
    me%atm_pressLB = me%atm_press / physics%fac(minLevel)%press

    !set constant parameter faraday and gasconstant in lattice units
    me%faraday = 96485.3365 ![C/mol]
    me%faradayLB = me%faraday/physics%fac(minLevel)%faraday

    !gas constant
    me%gasConst_R = 8.3144621 ![N m /(mol * K)]
    me%gasConst_R_LB = me%gasConst_R / physics%fac(minLevel)%gasConst

    call aot_table_close(L=conf, thandle=mix_handle)

    write(logUnit(1),*) ' Mixture properties '
    write(logUnit(1),*) '   theta_eq:                        ',                &
      & real(me%theta_eq)
    write(logUnit(1),*) '   Phy.reference mass density:      ', real(me%rho0)
    write(logUnit(1),*) '   Lat.reference mass density:      ', real(me%rho0LB)
    write(logUnit(1),*) '   Phy.reference mole density:      ',                &
      & real(me%moleDens0)
    write(logUnit(1),*) '   Lat.reference mole density:      ',                &
      & real(me%moleDens0LB)
    write(logUnit(1),*) '   Bulk modulus lattice B:          ',                &
      & real(me%bulk_modulusLB)
    write(logUnit(1),*) '   Diffusivity relaxation parameter:',                &
      & real(me%omega_diff)
    write(logUnit(1),*) '   Free parameter:                  ', real(me%paramB)
    write(logUnit(1),*) '   Shear viscosity Physical:        ', &
      & real(me%kine_viscosity)
    write(logUnit(1),*) '   Shear viscosityLB:               ', &
      & real(me%kine_viscosityLB)
    write(logUnit(1),*) '   Kinematic shear viscosity omega: ',                &
      & real(me%omega_kine)
    write(logUnit(1),*) '   Bulk viscosity Physical:         ',                &
      & real(me%bulk_viscosity)
    write(logUnit(1),*) '   Bulk viscosityLB:                ',                &
      & real(me%bulk_viscosityLB)
    write(logUnit(1),*) '   Omega of Higher order moments: ', real(me%omega_hom)

    write(logUnit(1),*)('   Setting omegas')
    call set_omegasLvl( mixture  = me,       &
      &                 minLevel = minLevel, &
      &                 maxLevel = maxLevel, &
      &                 physics  = physics   )

    write(logUnit(1),*) ' Forcing parameters'
    write(logUnit(1),*) '   ElectricField:               ',                    &
      & real(me%electricField)
    write(logUnit(1),*) '   GravityField:                ',                    &
      & real(me%gravityField)
    write(logUnit(1),*) '   Physical Temperature:        ', real(me%temp0)
    write(logUnit(1),*) '   Lattice Temperature:         ', real(me%temp0LB)
    write(logUnit(1),*) '   Atmospheric pressure:        ', real(me%atm_press)
    write(logUnit(1),*) '   Atm pressure lattice:        ', real(me%atm_pressLB)
    write(logUnit(1),*) '   Faraday constant:            ', real(me%faraday)
    write(logUnit(1),*) '   Faraday constant lattice:    ', real(me%faradayLB)
    write(logUnit(1),*) '   Ideal Gas constant:          ', real(me%gasconst_R)
    write(logUnit(1),*) '   Ideal Gas constant lattice:  ',                    &
      & real(me%gasconst_R_LB)

  end subroutine mus_load_mixture
! **************************************************************************** !


! **************************************************************************** !
  !> Set the omegas according to the time step setting
  subroutine set_omegasLvl( mixture, minLevel, maxLevel, physics )
    ! --------------------------------------------------------------------------
    type( mus_mixture_type ),intent(inout) :: mixture !< mixture type
    integer, intent(in) :: minLevel, maxLevel
    !> physics type to convert physics to lattice unit or vice versa
    type( mus_physics_type ), intent(in) :: physics
    ! --------------------------------------------------------------------------
    integer :: iLevel
    real(kind=rk) :: viscPhy_loc, error_rel
    ! --------------------------------------------------------------------------
    mixture%relaxLvl(:)%paramB = 0._rk
    mixture%relaxLvl(:)%omega_diff = 0._rk
    mixture%relaxLvl(:)%visc = 0._rk
    mixture%relaxLvl(:)%omega_kine = 0._rk
    mixture%relaxLvl(:)%bulkvisc = 0._rk

    do iLevel = minLevel, maxLevel
      write(logUnit(5),*) '    level: ', iLevel
      ! paramB is inverse of resistivity unit [s/m^2]
      mixture%relaxLvl( iLevel )%paramB = mixture%paramB                       &
        &                               * physics%fac( iLevel )%visc           &
        &                               / physics%fac( minLevel )%visc

      mixture%relaxLvl( iLevel )%omega_diff = mixture%relaxLvl(iLevel)%paramB  &
        &                                   * cs2
      write(logUnit(5),*)'    omega_diff ',                                    &
        &                real(mixture%relaxLvl( iLevel )%omega_diff)

      ! kinematic viscosity
      mixture%relaxLvl( iLevel )%visc = mixture%kine_viscosityLB               &
        &                             * physics%fac( minLevel )%visc           &
        &                             / physics%fac( iLevel )%visc
      mixture%relaxLvl(iLevel)%omega_kine = cs2/mixture%relaxLvl( iLevel )%visc
      write(logUnit(5),*)'    omega_kine ',                                    &
        &                real(mixture%relaxLvl( iLevel )%omega_kine)

      !bulk omega
      mixture%relaxLvl( iLevel )%bulkVisc = mixture%bulk_viscosityLB           &
        &                                 * physics%fac( minLevel )%visc       &
        &                                 / physics%fac( iLevel )%visc
    end do

    !> cross check whether omega at each level is set correctly
    !! by computing physical viscosity from omega and check it with
    !! specified physical kinematic viscosity
    !! \( \nu = (1/ \omega - 0.5)*cs2 \)
    do iLevel = minLevel, maxLevel
      viscPhy_loc = (1.0_rk/mixture%relaxLvl(iLevel)%omega_kine) * cs2 &
        &         * physics%fac(iLevel)%visc
      error_rel = (mixture%kine_viscosity-viscPhy_loc)/mixture%kine_viscosity
      write(logUnit(3),*) ' ATTENTION: Relative error between defined '//      &
        &                 'physical viscosity to'
      write(logUnit(3),*) ' computed physical viscosity on level: ', ilevel
      write(logUnit(3),*) '                                  is : ', error_rel
      if( abs(real(error_rel, kind=single_k)) > eps_single_k) then
        write(logUnit(1),*) 'Error: Physical kinematic viscosity computed '//  &
          &                 'from omega'
        write(logUnit(1),*) '       does not match with specified physical '// &
          &                 'kinematic viscosity'
        write(logUnit(1),*) ' Physical kinematic shear viscosity: ' ,          &
          &                 real(mixture%kine_viscosity)
        write(logUnit(1),*) ' Calculated phy kinematic shear viscosity: ',     &
          &                 real(viscPhy_loc)
        write(logUnit(1),*) ' Defined omega: ',                                &
          &                 real(mixture%relaxLvl(iLevel)%omega_kine)
        call tem_abort()
      end if
    end do

  end subroutine set_omegasLvl
! **************************************************************************** !


! **************************************************************************** !
  !> This routine write mixture properties into lua file
  subroutine mus_mixture_out( me, conf, schemeHeader )
    ! --------------------------------------------------------------------------
    !> mixture info
    type( mus_mixture_type ), intent(in) :: me
    type(aot_out_type) :: conf
    !> identifier of the scheme
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------
    select case( trim(schemeHeader%kind) )
    case('multispecies_gas', 'multispecies_liquid')
      call aot_out_open_table( put_conf = conf, tname = 'mixture')

      call aot_out_val( put_conf = conf, vname = 'rho0', val = me%rho0 )
      call aot_out_val( put_conf = conf, vname = 'omega_diff',                 &
        &               val = me%omega_diff )
      call aot_out_val( put_conf = conf, vname = 'omega_kine',                 &
        &               val = me%omega_kine )
      call aot_out_val( put_conf = conf, vname = 'omega_hom',                  &
        &               val = me%omega_hom )
      call aot_out_val( put_conf = conf, vname = 'bulk_viscosity',             &
        &               val = me%bulk_viscosity )
      call aot_out_val( put_conf = conf, vname = 'paramB', val = me%paramB )
      call aot_out_val( put_conf = conf, vname = 'theta_eq',                   &
        &               val = me%theta_eq )
      call aot_out_val( put_conf = conf, vname = 'moleDens0',                  &
        &               val = me%moleDens0 )
      call aot_out_val( put_conf = conf, vname = 'temp', val = me%temp0 )
      call aot_out_val( put_conf = conf, vname = 'atm_press',                  &
        &               val = me%atm_pressLB )
      call aot_out_close_table( put_conf = conf )
    end select

  end subroutine mus_mixture_out
end module mus_mixture_module
