! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2017-2018 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2018 Jana Gericke <jana.gericke@uni-siegen.de>
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
! **********************************************************************!
!! In this module the poission equation is considered.
!! The file type contains all relevant information.
!! It is a one dimensional case.
!> author: sindhuja

module mus_poisson_module

  use env_module,         only: rk
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit, tem_toStr

  use aotus_module,     only: flu_State, aot_get_val, aoterr_Fatal, &
  &                           aoterr_NonExistent, aoterr_WrongType
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_val, aot_out_open_table,   &
  &                           aot_out_close_table

  use mus_physics_module, only: mus_physics_type, faraday, gasConst_R

  implicit none

  private

  public :: mus_poisson_type
  public :: mus_load_poisson
  public :: mus_load_poissonBoltzmann

  !> Contains information to calculate rhs for poisson boltzmann equation.
  !! Definition of linear and non-linear poisson boltzmann equation
  !! can be found in
  !! Masilamani, K. (2010). WaLBerla : Investigation of Electrostatic Effects
  !! in Particulate and Electro-Osmotic Flows. Master Thesis.
  !! FRIEDRICH-ALEXANDER-UNIVERSITÄT ERLANGEN-NÜRNBERG.
  type mus_poisson_boltzmann_type
    !> Neccesary if source term is poisson_boltzmann
    logical       :: active
    !> abosulte temperature in Kelvin
    real(kind=rk) :: temp
    !> Number of ions
    integer :: nIons
    !> valence of the ion
    !! size: nions
    integer, allocatable :: valence(:)
    !> Mole density of ions in the bulk
    real(kind=rk) :: moleDens0
    !> RHS coeff for linear poisson boltzmann equation
    real(kind=rk) :: RHS_coeff
    !> Gas constant in lattice unit
    real(kind=rk) :: gasConst_R_LB
    !> Faraday constant in lattice unit
    real(kind=rk) :: faradayLB
    !!> Boltzmann constant in lattice unit
    !real(kind=rk) :: k_b
    !!> Fundamental charge in lattice unit
    !real(kind=rk) :: charge
  end type mus_poisson_boltzmann_type

  !> Contains information to solve the poission equation
  type mus_poisson_type
    !> Potential diffusivty to tune omega and stability
    real(kind=rk) :: pot_diff
    !> relaxation parameter
    real(kind=rk) :: omega
    !> the dielectric constant C^2 J^-1 m^-1
    real(kind=rk) :: permittivity
    !> information of poisson boltmann equation
    type(mus_poisson_boltzmann_type) :: PB
  end type mus_poisson_type

contains

  ! ****************************************************************************
  !> load input to solve poisson equation
  subroutine mus_load_poisson( me, conf, parent, minLevel, cs_lattice, physics,&
    &                          schemeKind )
    !-------------------------------------------------------------------------
    !> poisson type
    type(mus_poisson_type ), intent(out) :: me
    !> flu state
    type( flu_state )             :: conf
    !> parent handle
    integer, intent(in), optional  :: parent
    !> minlevel
    integer, intent(in) :: minLevel
    !> physics type to convert physics to lattice unit or vice versa
    type( mus_physics_type ), intent(in) :: physics
    !> lattice speed of sound calculated for defined stencil layout
    !! required to compute omega from potential diffusivity
    real(kind=rk), intent(in) :: cs_lattice
    !> scheme kind
    character(len=*), intent(in) :: schemeKind
    ! -------------------------------------------------------------------------
    integer :: poisson_handle
    integer :: iError, iIon
    real(kind=rk) :: debye_length_inv_sqr, valence_sqr_sum, debye_length
    ! -------------------------------------------------------------------------

    ! if poisson informations in scheme table parentHandle /= 0
    if ( present(parent) ) then
      call aot_table_open( L       = conf,           &
        &                  parent  = parent,         &
        &                  thandle = poisson_handle, &
        &                  key     = 'poisson'       )
    else
      call aot_table_open( L=conf, thandle = poisson_handle, &
        &                  key = 'poisson' )
    end if
    if ( poisson_handle == 0 ) then
       write(logUnit(1),*)'No poisson table defined'
       call tem_abort()
    end if

    write(logUnit(1),*) 'Loading poisson informations'
    ! Load potential diffusivity to control evolution speed and stability
    call aot_get_val(L       = conf,                    &
      &              thandle = poisson_handle,          &
      &              key     = 'potential_diffusivity', &
      &              val     = me%pot_diff,             &
      &              ErrCode = iError                   )
    me%pot_diff = me%pot_diff / physics%fac(minLevel)%visc

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving '// &
        &                'potentail_diffusivity:'
      if (btest(iError, aoterr_NonExistent))        &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))            &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if

    ! Relaxation parameter omega is compted from potential diffusivity
    me%omega = 1.0_rk/(me%pot_diff/cs_lattice**2 + 0.5_rk)

    ! load permittivity
    call aot_get_val(L       = conf,            &
      &              thandle = poisson_handle,  &
      &              key     = 'permittivity',  &
      &              val     = me%permittivity, &
      &              ErrCode = iError           )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*)'FATAL Error occured, while retrieving '// &
        &                'permittivity:'
      if (btest(iError, aoterr_NonExistent))        &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))            &
        & write(logUnit(1),*)'Variable has wrong type!'
      write(logUnit(1),*)'STOPPING'
      call tem_abort()
    end if
    me%permittivity = me%permittivity / ( physics%coulomb0**2                &
      &             / (physics%fac(minLevel)%energy*physics%dxLvl(minLevel)) )


    select case (trim(schemeKind))
    case('poisson_boltzmann_linear', 'poisson_boltzmann_nonlinear')
      call mus_load_poissonBoltzmann( me         = me%PB,          &
        &                             conf       = conf,           &
        &                             parent     = poisson_handle, &
        &                             physics    = physics         )
    case default
      me%PB%active = .false.
    end select

    call aot_table_close( L=conf, thandle=poisson_handle )

    write(logUnit(1),"(A)") 'Poisson properties:'
    write(logUnit(1),"(A)") '  potential_diffusivty: ' &
      &                   //trim(tem_toStr(me%pot_diff))
    write(logUnit(1),"(A)") '  omega: '//trim(tem_toStr(me%omega))
    write(logUnit(1),"(A)") '  permittivity: '//trim(tem_toStr(me%permittivity))
    if (me%PB%active) then
      ! Convert boltzmann constant into lattice units
      ! me%PB%k_b = k_b / ( physics%fac(minLevel)%energy / physics%temp0 )

      ! Convert Faraday constant into lattice unit
      me%PB%faradayLB = faraday/physics%fac(minLevel)%faraday

      ! convert gas constant into lattice unit
      me%PB%gasConst_R_LB = gasConst_R / physics%fac(minLevel)%gasConst

      ! electric charge in lattice units
      ! me%PB%charge = coulomb_ref / physics%coulomb0

      ! calculate debye length in lattice units
      valence_sqr_sum = sum(me%PB%valence**2)
      !debye_length_inv_sqr = valence_sqr_sum * me%PB%moleDens0    &
      !  &                   * me%PB%charge**2                     &
      !  &                  / (me%permittivity*me%PB%k_b*me%PB%temp)
      debye_length_inv_sqr = valence_sqr_sum * me%PB%moleDens0              &
        &                   * me%PB%faradayLB**2                            &
        &                  / (me%permittivity*me%PB%gasConst_R_LB*me%PB%temp)
      me%PB%RHS_coeff = debye_length_inv_sqr
      debye_length = 1.0_rk/sqrt(debye_length_inv_sqr) &
        &          * physics%dxLvl(minLevel)

      write(logUnit(1),"(A)") '  Paramters for poisson boltzmann Eq:'
      write(logUnit(1),"(A)") '  temp: '//trim(tem_toStr(me%PB%temp))
      write(logUnit(1),"(A)") '  moleDens: '//trim(tem_toStr(me%PB%moleDens0))
      write(logUnit(1),"(A)") '  faradayLB: '//trim(tem_toStr(me%PB%faradayLB))
      write(logUnit(1),"(A)") '  Gas Const R: ' &
        &                     //trim(tem_toStr(me%PB%gasConst_R_LB))
      write(logUnit(1),"(A)") '  valence : '
      do iIon = 1, me%PB%nIons
        write(logUnit(1),"(A)") '     '//trim(tem_toStr(me%PB%valence(iIon)))
      end do
      write(logUnit(1),"(A)") '  RHS_coeff: '//trim(tem_toStr(me%PB%RHS_coeff))
      write(logUnit(1),"(A)") '  Debye length: '//trim(tem_toStr(debye_length))
    end if

  end subroutine mus_load_poisson
  ! ***************************************************************************!

  ! ****************************************************************************
  !> Load input to solve poisson boltzmann equation
  subroutine mus_load_poissonBoltzmann( me, conf, parent, physics )
    !-------------------------------------------------------------------------
    !> poisson bolztmann type
    type(mus_poisson_boltzmann_type), intent(out) :: me
    !> flu state
    type( flu_state )             :: conf
    !> parent handle
    integer, intent(in), optional  :: parent
    !> physics type to convert physics to lattice unit or vice versa
    type( mus_physics_type ), intent(in) :: physics
    ! -------------------------------------------------------------------------
    integer :: PB_handle, valence_handle
    integer :: iError
    integer, allocatable :: vError(:), errFatal(:)
    ! -------------------------------------------------------------------------
    ! Load information for poisson boltzmann equation
    call aot_table_open( L       = conf,           &
      &                  parent  = parent,         &
      &                  thandle = PB_handle,      &
      &                  key = 'poisson_boltzmann' )
    if ( PB_handle /= 0 ) then
       me%active = .true.

      ! load absolute temperature
      call aot_get_val(L       = conf,       &
        &              thandle = PB_handle,  &
        &              key     = 'temp',     &
        &              val     = me%temp,    &
        &              default = 273._rk,    &
        &              ErrCode = iError      )
      me%temp = me%temp/physics%temp0

      ! load mole density at bulk
      call aot_get_val( L       = conf,         &
        &               thandle = PB_handle,    &
        &               key     = 'moleDens0',  &
        &               val     = me%moleDens0, &
        &               ErrCode = iError        )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*)'FATAL Error occured, while retrieving '// &
          &                'moleDens0:'
        if (btest(iError, aoterr_NonExistent))        &
          & write(logUnit(1),*)'Variable not existent!'
        if (btest(iError, aoterr_WrongType))            &
          & write(logUnit(1),*)'Variable has wrong type!'
        write(logUnit(1),*)'STOPPING'
        call tem_abort()
      end if

      me%moleDens0 = me%moleDens0 / physics%moleDens0

      ! load valence
      call aot_table_open( L       = conf,           &
        &                  parent  = PB_handle,      &
        &                  thandle = valence_handle, &
        &                  key     = 'valence'       )

      me%nIons = aot_table_length( L=conf, thandle = valence_handle )
      allocate( errFatal(me%nIons) )
      errFatal = aotErr_Fatal
      call aot_get_val( L         = conf,       &
        &               thandle   = PB_handle,  &
        &               key       = 'valence',  &
        &               val       = me%valence, &
        &               maxlength = me%nIons,   &
        &               ErrCode   = vError      )

      if ( any(btest(vError, errFatal)) ) then
         write(logUnit(1),*) 'FATAL Error occured, while retrieving '//      &
            &            'resi_coeff table'
         call tem_abort()
      endif

!      call aot_get_val(L       = conf,       &
!        &              thandle = PB_handle,  &
!        &              key     = 'valence ', &
!        &              val     = me%valence, &
!        &              default = 1,          &
!        &              ErrCode = iError      )
    else
       write(logUnit(3),*)'No poisson boltzmann table defined'
       me%active = .false.
    end if

    call aot_table_close( L=conf, thandle=PB_handle )

  end subroutine mus_load_poissonBoltzmann
  ! ***************************************************************************!

end module mus_poisson_module
