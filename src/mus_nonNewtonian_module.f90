! Copyright (c) 2013-2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013, 2016, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2020 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> author: Jiaxing Qi, Kannan Masilamani
!! This module keeps all information about the nonNewtonian models.
!! Contains routines which calculates non-Newtonian kinematic
!! viscosity according to non-Newtonian model.
!!
!! It supports three non-Newtonian models:
!! Power law, Carrear Yasuda and Casson.
!! All these models are described in
!! Ashrafizaadeh, M., & Bakhshaei, H. (2009). A comparison of non-Newtonian
!! models for lattice Boltzmann blood flow simulations.
!! Computers and Mathematics with Applications, 58(5), 1045–1054.
!!
!! For further information about the theory visit the
!! [non-newtonian theory page](../page/mus_nonNewtonianTheory.html).
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
module mus_nonNewtonian_module

  ! include treelm modules
  use env_module,           only: rk, labellen
  use tem_tools_module,     only: tem_horizontalSpacer, upper_to_lower
  use tem_logging_module,   only: logUnit
  use tem_aux_module,       only: tem_abort
  use tem_param_module,     only: rho0, rho0Inv, cs2, cs2inv, div1_2

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val
  use aot_table_module, only: aot_table_open, aot_table_close
  use aot_out_module,   only: aot_out_type, aot_out_val, aot_out_open_table,   &
    &                         aot_out_close_table

  ! include musubi modules
  use mus_physics_module,       only: mus_convertFac_type
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_scheme_header_module,  only: mus_scheme_header_type
  use mus_derivedQuantities_module2, only: secondMom_minus_cs2_3D

  implicit none

  private

  public :: mus_nNwtn_type
  public :: mus_nNwtn_load
  public :: mus_nNwtn_save2lua
  public :: mus_nNwtn_dump2outUnit
  public :: mus_assign_nNwtnVisc_ptr
  public :: calcVisc_CY

  !> Identifier for Power-Law model
  integer, parameter, public :: powerLaw = 1

  !> Identifier for Casson model
  integer, parameter, public :: casson = 2

  !> Identifier for Carreau-Yasuda model
  integer, parameter, public :: carreauYasuda = 3

  ! ---------------------------------------------------------------------------
  !> The nonNewtonian power law model parameter
  !!
  !! This date type gathers parameters of a power law (PL) model.
  !! It is encapsulated in mus_nNwtn_type.
  type mus_nNwtn_PL_type

    !> exponentiation parameter
    real(kind=rk) :: n = 0.5_rk

    !> Unit consistency index.
    !! Dynamic viscosity parameter when shear rate equals to 1.
    !! i.e. backgroud viscosity for power-law model
    real(kind=rk) :: visc0 = 0.0035_rk

    !> parameter for computation
    real(kind=rk) :: nMinus1

  end type mus_nNwtn_PL_type
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  !> The nonNewtonian power law model parameter
  !!
  !! This date type gathers parameters of the Carreau-Yasuda (CY) model
  !! It is encapsulated in mus_nNwtn_type.
  type mus_nNwtn_CY_type

    !> model parameter
    real(kind=rk) :: n = 0.2128_rk

    !> model parameter
    real(kind=rk) :: a = 0.64_rk

    !> model parameter
    real(kind=rk) :: lambda = 8.2_rk

    !> model parameter, dynamic viscosity at zero shear-rate
    real(kind=rk) :: visc0 = 0.16_rk

    !> model parameter, dynamic viscosity in infinity shear-rate
    real(kind=rk) :: viscInf = 0.0035_rk

    !> calculated parameter for later usage, nMinus1Div_a = (n-1)/a
    real(kind=rk) :: nMinus1Div_a = (0.2128_rk - 1._rk) / 0.64_rk

  end type mus_nNwtn_CY_type
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  !> The nonNewtonian power law model parameter
  !!
  !! This date type gathers parameters of the Casson model
  !! It is encapsulated in mus_nNwtn_type.
  type mus_nNwtn_CS_type

    !> model parameter
    real(kind=rk) :: k0 = 0.1937_rk

    !> model parameter
    real(kind=rk) :: k1 = 0.055_rk

  end type mus_nNwtn_CS_type
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  !> The nonNewtonian fluid feature description type
  !!
  !! This date type gathers related parameters of a nonNewtonian fluid.
  type mus_nNwtn_type
    !> Indicator whether nonNewtonian feature is active
    !! maybe not useful. schemeHeader%kind can used to check if nNwtn is active
    logical :: active = .false.

    !> nonNewtonian fluid model label
    character(len=labellen) :: label

    !> nonNewtonian fluid model identifier
    integer :: model

    !> Power Law (PL) model type
    type( mus_nNwtn_PL_type ) :: PL

    !> Carreau-Yasuda (CY) model type
    type( mus_nNwtn_CY_type ) :: CY

    !> Casson model type
    type( mus_nNwtn_CS_type ) :: CS

    !> this procedure compute kinematic viscosity in lattice unit on current
    !! level from preCollision PDF based on non-Newtonian model.
    !! It uses shear-rate to compute viscosity.
    !! Non-newtonian model is given in dynamic viscosity in physical unit
    !! so it is dimensionalized using viscDyna in physics%fac and
    !! lattice kinematic viscosity = lattice dynamic viscosity
    !! / rho (local density)  for compressible model and
    !! lattice kinematic viscosity = lattice dynamic viscosity / rho0
    !! for incompressible model.
    procedure(proc_calc_nNwtn_visc_fromPreColPDF), pointer :: calcVisc => null()

  end type mus_nNwtn_type
  ! ----------------------------------------------------------------------------

  !> Interface to calculate kinematic viscosity for non-Newtonian model.
  !! Viscosity is computed from shear rate which is computed from strain rate
  !! which is computed from nonEquilibrium PDF which in turn is computed from
  !! pre-collision PDF
  abstract interface
    subroutine proc_calc_nNwtn_visc_fromPreColPDF(nNwtn, viscKine, omega, &
      & state, neigh, auxField, densPos, velPos, nSize, nSolve, nScalars, &
      & nAuxScalars, layout, convFac)
      import :: rk, mus_nNwtn_type, mus_scheme_layout_type, mus_convertFac_type

      !> contains non-Newtonian model parameters loaded from config file
      class(mus_nNwtn_type), intent(in) :: nNwtn

      !> output is kinematic viscosity from non-Netonian model
      real(kind=rk), intent(inout) :: viscKine(:)

      !> Kinematic viscosity omega from last timestep
      real(kind=rk), intent(in) :: omega(:)

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

      !> conversion factor to convert lattice to physical units
      type(mus_convertFac_type), intent(in) :: convFac
    end subroutine proc_calc_nNwtn_visc_fromPreColPDF
  end interface


contains


  ! ************************************************************************** !
  !> Read in the nonNewtonian table from Lua file and dump parameters to logUnit
  !! Specificly, this routine calls each model parameter loader.
  !!
  subroutine mus_nNwtn_load( me, conf, parent )
    ! --------------------------------------------------------------------------
    !> nonNewtonian type
    type( mus_nNwtn_type ), intent(out) :: me
    !> lua state
    type( flu_state ), intent(inout) :: conf
    !> parent handle
    integer, intent(in), optional  :: parent
    ! --------------------------------------------------------------------------
    integer :: nonNwt_table
    integer :: iError
    CHARACTER(LEN=12) :: nonNwt_table_str
    ! --------------------------------------------------------------------------

    nonNwt_table_str = "nonNewtonian"

    ! if nonNewtonian informations in scheme table parentHandle /= 0
    call aot_table_open( L       = conf,            &
      &                  parent  = parent,          &
      &                  thandle = nonNwt_table,    &
      &                  key     = nonNwt_table_str )

    if ( nonNwt_table == 0 ) then
      write(logUnit(1),*)'No nonNewtonian table defined'
      me%active = .false.
      return
    endif

    ! when table exist, read in parameters from table
    write(logUnit(1),*) 'Loading nonNewtonian informations'
    ! Set nonNewtonian feature in on
    me%active = .true.

    ! load model label name
    call aot_get_val(L       = conf,        &
      &              tHandle = nonNwt_table,&
      &              key     = 'model',     &
      &              val     = me%label,    &
      &              default = 'power_law', &
      &              ErrCode = iError       )

    ! load model parameters by calling model loader
    ! set model identifier
    select case( trim(upper_to_lower(me%label)) )
    case ( 'power_law' )
      me%model = powerLaw
      call mus_nNwtn_PL_load( me%PL, conf, nonNwt_table )
    case ( 'carreau_yasuda')
      me%model = carreauYasuda
      call mus_nNwtn_CY_load( me%CY, conf, nonNwt_table )
    case ( 'casson' )
      me%model = casson
      call mus_nNwtn_CS_load( me%CS, conf, nonNwt_table )
    case default
      call tem_abort('Error: Unknown non-Newtonian model')
    end select

    call aot_table_close( L=conf, thandle = nonNwt_table )

  end subroutine mus_nNwtn_load
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> Read in the nonNewtonian Power Law (PL) model parameters from Lua file
  subroutine mus_nNwtn_PL_load( me, conf, tHandle )
    ! --------------------------------------------------------------------------
    !> nonNewtonian type
    type( mus_nNwtn_PL_type ), intent(out) :: me
    !> lua state
    type( flu_state ), intent(inout) :: conf
    !> nonNewtonian table handle
    integer, intent(inout) :: tHandle
    ! --------------------------------------------------------------------------
    integer :: iError
    ! --------------------------------------------------------------------------

    ! load n
    call aot_get_val( L       = conf,         &
      &               thandle = tHandle,      &
      &               key     = 'n',          &
      &               val     = me%n,         &
      &               default = 0.5_rk,       &
      &               ErrCode = iError        )

    ! load k
    call aot_get_val( L       = conf,                  &
      &               thandle = tHandle,               &
      &               key     = 'dynamic_viscosity_0', &
      &               val     = me%visc0,              &
      &               default = 0.0000035_rk,          &
      &               ErrCode = iError                 )

    me%nMinus1 = me%n - 1._rk

  end subroutine mus_nNwtn_PL_load
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Read in the nonNewtonian Carreau-Yasuda (CY) model parameters from Lua file
  subroutine mus_nNwtn_CY_load( me, conf, tHandle )
    ! --------------------------------------------------------------------------
    !> nonNewtonian type
    type( mus_nNwtn_CY_type ), intent(out) :: me
    !> lua state
    type( flu_state ), intent(inout) :: conf
    !> nonNewtonian table handle
    integer, intent(inout) :: tHandle
    ! --------------------------------------------------------------------------
    integer :: iError
    ! --------------------------------------------------------------------------

    ! load visc0
    call aot_get_val(L       = conf,                  &
      &              thandle = tHandle,               &
      &              key     = 'dynamic_viscosity_0', &
      &              val     = me%visc0,              &
      &              default = 0.16_rk,               &
      &              ErrCode = iError                 )

    ! load viscInf
    call aot_get_val( L       = conf,                      &
      &               thandle = tHandle,                   &
      &               key     = 'dynamic_viscosity_infty', &
      &               val     = me%viscInf,                &
      &               default = 0.0035_rk,                 &
      &               ErrCode = iError                     )

    ! load lambda
    call aot_get_val( L       = conf,         &
      &               thandle = tHandle,      &
      &               key     = 'lambda',     &
      &               val     = me%lambda,    &
      &               default = 8.2_rk,       &
      &               ErrCode = iError        )

    ! load a
    call aot_get_val( L       = conf,         &
      &               thandle = tHandle,      &
      &               key     = 'a',          &
      &               val     = me%a,         &
      &               default = 0.64_rk,      &
      &               ErrCode = iError        )

    ! load n
    call aot_get_val( L       = conf,         &
      &               thandle = tHandle,      &
      &               key     = 'n',          &
      &               val     = me%n,         &
      &               default = 0.2128_rk,    &
      &               ErrCode = iError        )

    ! calculate intermediate parameter
    me%nMinus1Div_a = ( me%n - 1._rk ) / me%a

  end subroutine mus_nNwtn_CY_load
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Read in the nonNewtonian Casson model parameters from Lua file
  subroutine mus_nNwtn_CS_load( me, conf, tHandle )
    ! --------------------------------------------------------------------------
    !> nonNewtonian type
    type( mus_nNwtn_CS_type ), intent(out) :: me
    !> lua state
    type( flu_state ), intent(inout) :: conf
    !> nonNewtonian table handle
    integer, intent(inout) :: tHandle
    ! --------------------------------------------------------------------------
    integer :: iError
    ! --------------------------------------------------------------------------

    ! load k0
    call aot_get_val(L       = conf,      &
      &              tHandle = tHandle,   &
      &              key     = 'k0',      &
      &              val     = me%k0,     &
      &              default = 0.1937_rk, &
      &              ErrCode = iError     )

    ! load k1
    call aot_get_val( L       = conf,     &
      &               tHandle = tHandle,  &
      &               key     = 'k1',     &
      &               val     = me%k1,    &
      &               default = 0.055_rk, &
      &               ErrCode = iError    )

  end subroutine mus_nNwtn_CS_load
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> write nonNewtonian fluid parameters into a lua file
  !!
  subroutine mus_nNwtn_save2lua( me, conf )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_type ), intent(in) :: me
    type( aot_out_type ) :: conf
    ! --------------------------------------------------------------------------

    call aot_out_open_table( put_conf = conf, tname = 'nonNewtonian' )

    call aot_out_val( put_conf = conf,          &
      &               vname    = 'model',       &
      &               val      = trim(me%label) )

    select case( me%model )
      case ( powerLaw )
        call mus_nNwtn_PL_save( me%PL, conf )
      case ( carreauYasuda )
        call mus_nNwtn_CY_save( me%CY, conf )
      case ( casson )
        call mus_nNwtn_CS_save( me%CS, conf )
    end select

    call aot_out_close_table( put_conf = conf )

  end subroutine mus_nNwtn_save2lua
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> write nonNewtonian Power Law (PL) parameters into a lua file
  !!
  subroutine mus_nNwtn_PL_save( me, conf )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_PL_type ), intent(in) :: me
    type( aot_out_type ) :: conf
    ! --------------------------------------------------------------------------

    call aot_out_val( put_conf = conf, &
      &               vname    = 'n',  &
      &               val      = me%n  )
    call aot_out_val( put_conf = conf,                  &
      &               vname    = 'dynamic_viscosity_0', &
      &               val      = me%visc0               )

  end subroutine mus_nNwtn_PL_save
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> write nonNewtonian (CY) parameters into a lua file
  !!
  subroutine mus_nNwtn_CY_save( me, conf )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_CY_type ), intent(in) :: me
    type( aot_out_type ) :: conf
    ! --------------------------------------------------------------------------

    call aot_out_val( put_conf = conf,  &
      &               vname    = 'n',   &
      &               val      = me%n   )
    call aot_out_val( put_conf = conf,  &
      &               vname    = 'a',   &
      &               val      = me%a   )
    call aot_out_val( put_conf = conf,      &
      &               vname    = 'lambda',  &
      &               val      = me%lambda  )
    call aot_out_val( put_conf = conf,     &
      &               vname    = 'dynamic_viscosity_0', &
      &               val      = me%visc0               )
    call aot_out_val( put_conf = conf,                      &
      &               vname    = 'dynamic_viscosity_infty', &
      &               val      = me%viscInf                 )

  end subroutine mus_nNwtn_CY_save
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> write nonNewtonian Casson parameters into a lua file
  !!
  subroutine mus_nNwtn_CS_save( me, conf )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_CS_type ), intent(in) :: me
    type( aot_out_type ) :: conf
    ! --------------------------------------------------------------------------

    call aot_out_val( put_conf = conf, &
      &               vname    = 'k0', &
      &               val      = me%k0 )
    call aot_out_val( put_conf = conf, &
      &               vname    = 'k1', &
      &               val      = me%k1 )

  end subroutine mus_nNwtn_CS_save
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Dump nonNewtonian fluid parameters to outUnit
  !!
  subroutine mus_nNwtn_dump2outUnit( me, outUnit )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_type ), intent(in) :: me
    integer, intent(in) :: outUnit
    ! --------------------------------------------------------------------------

    if ( me%active ) then
      write(outUnit,'(A)') 'nonNewtonian fluid parameters:'
      write(outUnit,'(A)') '  model label: ', trim(me%label)

      ! dump model parameters by calling model dumper
      select case( me%model )
        case ( powerLaw )
          call mus_nNwtn_PL_dump( me%PL, outUnit )
        case ( carreauYasuda )
          call mus_nNwtn_CY_dump( me%CY, outUnit )
        case ( casson )
          call mus_nNwtn_CS_dump( me%CS, outUnit )
      end select

    else
      write(outUnit,'(A)') 'No nonNewtonian table defined.'
    end if

    call tem_horizontalSpacer( fUnit = outUnit )

  end subroutine mus_nNwtn_dump2outUnit
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Dump nonNewtonian Power Law (PL) parameters to outUnit
  !!
  subroutine mus_nNwtn_PL_dump( me, outUnit )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_PL_type ), intent(in) :: me
    integer, intent(in) :: outUnit
    ! --------------------------------------------------------------------------

    write(outUnit,"( '  n = ', F8.4)") me%n
    write(outUnit,"( '  dynamic_viscosity_0 = ', F8.4)") me%visc0

  end subroutine mus_nNwtn_PL_dump
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Dump nonNewtonian (CY) parameters to outUnit
  !!
  subroutine mus_nNwtn_CY_dump( me, outUnit )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_CY_type ), intent(in) :: me
    integer, intent(in) :: outUnit
    ! --------------------------------------------------------------------------

    write(outUnit, "('  n       = ', F8.4)") me%n
    write(outUnit, "('  a       = ', F8.4)") me%a
    write(outUnit, "('  lambda  = ', F8.4)") me%lambda
    write(outUnit, "('  dynamic_viscosity_0     = ', F8.4)") me%visc0
    write(outUnit, "('  dynamic_viscosity_infty = ', F8.4)") me%viscInf

  end subroutine mus_nNwtn_CY_dump
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Dump nonNewtonian (CY) parameters to outUnit
  !!
  subroutine mus_nNwtn_CS_dump( me, outUnit )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_CS_type ), intent(in) :: me
    integer, intent(in) :: outUnit
    ! --------------------------------------------------------------------------

    write(outUnit,"('  k0 = ', F8.4)") me%k0
    write(outUnit,"('  k1 = ', F8.4)") me%k1

  end subroutine mus_nNwtn_CS_dump
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> nonNewtonian power-law model
  !!
  real(kind=rk) function viscPhy_PL( me, shearRate )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_type ), intent(in)  :: me
    real(kind=rk),          intent(in)  :: shearRate
    ! --------------------------------------------------------------------------

    viscPhy_PL = ( shearRate ** me%PL%nMinus1 ) * me%PL%visc0

  end function viscPhy_PL
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> nonNewtonian Casson model
  !!
  real(kind=rk) function viscPhy_CS( me, shearRate )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_type ), intent(in)  :: me
    real(kind=rk),           intent(in)  :: shearRate
    ! --------------------------------------------------------------------------
    real(kind=rk) :: t

    t = ( me%CS%k0 + me%CS%k1 * sqrt(shearRate) )
    viscPhy_CS = t * t / shearRate

  end function viscPhy_CS
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> nonNewtonian Carreau-Yasuda model
  !!
  real(kind=rk) function viscPhy_CY( me, shearRate )
    ! --------------------------------------------------------------------------
    !> nonNewtonian parameters
    type( mus_nNwtn_type ), intent(in)  :: me
    real(kind=rk),           intent(in)  :: shearRate
    ! --------------------------------------------------------------------------
    real(kind=rk) :: t

    t = ( 1._rk + (me%CY%lambda*shearRate) ** me%CY%a ) ** me%CY%nMinus1Div_a
    viscPhy_CY = me%CY%viscInf + ( me%CY%visc0 - me%CY%viscInf ) * t

  end function viscPhy_CY
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine assigns function pointer to compute non-Newtonian viscosity
  subroutine mus_assign_nNwtnVisc_ptr(nNwtn, schemeHeader)
    ! --------------------------------------------------------------------------
    !> non-Newtonian type
    type(mus_nNwtn_type), intent(inout) :: nNwtn
    !> scheme header
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------
    select case(trim(schemeHeader%kind))
    case('fluid')
      select case(nNwtn%model)
      case (powerLaw)
        nNwtn%calcVisc => calcVisc_PL
      case (carreauYasuda)
        nNwtn%calcVisc => calcVisc_CY
      case (casson)
        nNwtn%calcVisc => calcVisc_CS
      end select
    case('fluid_incompressible')
      select case(nNwtn%model)
      case (powerLaw)
        nNwtn%calcVisc => calcVisc_incomp_PL
      case (carreauYasuda)
        nNwtn%calcVisc => calcVisc_incomp_CY
      case (casson)
        nNwtn%calcVisc => calcVisc_incomp_CS
      end select
    case default
      call tem_abort('Unknown scheme kind for non-Newtonian model')
    end select

  end subroutine mus_assign_nNwtnVisc_ptr
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate kinematic viscosity from nonNewtonian power-law model.
  !! $\mu = K shearRate^(n-1)$.
  !! Shear rate is computed from strain rate which is computed from
  !! nonEquilibrium PDF which in turn computed from pre-collision PDF
  subroutine calcVisc_PL(nNwtn, viscKine, omega, state, neigh, auxField, &
    & densPos, velPos, nSize, nSolve, nScalars, nAuxScalars, layout, convFac)
    ! --------------------------------------------------------------------------
    !> contains non-Newtonian model parameters loaded from config file
    class(mus_nNwtn_type), intent(in) :: nNwtn
    !> output is physical kinematic viscosity will be overwritten by
    !! non-Netonian model
    real(kind=rk), intent(inout) :: viscKine(:)
    !> Kinematic viscosity omega from last timestep
    real(kind=rk), intent(in) :: omega(:)
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
    !> conversion factor to convert lattice to physical units
    type(mus_convertFac_type), intent(in) :: convFac
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: rho, inv_rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(6), nEqTensMag
    real(kind=rk) :: shearRate, strainRate, viscDynaPhy, coeffSR
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    ! constant coefficients in strainRate computation
    coeffSR = div1_2 * cs2inv * convFac%strainRate

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
        f_preCol(iDir) = state(                                &
          &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      ! Access density and velocity from auxField
      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      inv_rho = 1.0_rk/rho
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = auxField( elemOff + velPos(3) )

      ! Calculate the equilibrium distribution function
      fEq = layout%quantities%pdfEq_ptr( rho = rho,               &
        &                                vel = vel,               &
        &                                QQ = QQ                  )

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric deviatoric second-order tensor of
      ! nonEquilibrium part
      ! the static part cs2 I is usually neglected for weakly compressible flows
      ! however, in current implementation it is considered
      nEqTens = secondMom_minus_cs2_3D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      !nEqTens = nEqTens * (-1.5_rk) * omega(iElem)*convFac%strainRate*inv_rho
      ! compute strain
      ! magnitude of second-order tensor
      nEqTensMag = sqrt(nEqTens(1)**2 + nEqTens(2)**2 + nEqTens(3)**2 &
        &        + 2.0_rk*(nEqTens(4)**2 + nEqTens(5)**2 + nEqTens(6)**2) )

      ! omega from last time step
      ! convert shear-rate into physical unit because only
      ! non-Newtonian model requies it.
      ! physical unit conversion factor is pre-multiplied in coeffSR
      ! KM: Actual formula to compute strainrate has negative but since we
      ! calculating a magnitude, it is not used
      strainRate = coeffSR * omega(iElem) * inv_rho * nEqTensMag
      !strainRate = nEqTensMag

      ! compute shearRate
      shearRate = 2.0_rk * strainRate

      ! compute physical dynamic viscosity from non-Newtonian powerlaw model
      viscDynaPhy = (shearRate ** nNwtn%PL%nMinus1) * nNwtn%PL%visc0
      ! viscKine_L = viscDyna_L / rho
      viscKine(iElem) = (viscDynaPhy / convFac%viscDyna) * inv_rho

    end do

  end subroutine calcVisc_PL
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate kinematic viscosity from nonNewtonian Casson model.
  !! $\mu = (k0 + k1 * sqrt(shearRate))^2/shearRate$
  !!
  !! Shear rate is computed from strain rate which is computed from
  !! nonEquilibrium PDF which in turn computed from pre-collision PDF
  subroutine calcVisc_CS(nNwtn, viscKine, omega, state, neigh, auxField, &
    & densPos, velPos, nSize, nSolve, nScalars, nAuxScalars, layout, convFac)
    ! --------------------------------------------------------------------------
    !> contains non-Newtonian model parameters loaded from config file
    class(mus_nNwtn_type), intent(in) :: nNwtn
    !> output is physical kinematic viscosity will be overwritten by
    !! non-Netonian model
    real(kind=rk), intent(inout) :: viscKine(:)
    !> Kinematic viscosity omega from last timestep
    real(kind=rk), intent(in) :: omega(:)
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
    !> conversion factor to convert lattice to physical units
    type(mus_convertFac_type), intent(in) :: convFac
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: rho, inv_rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(6), nEqTensMag
    real(kind=rk) :: shearRate, strainRate, viscTerm, coeffSR, viscDynaPhy
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    ! constant coefficients in strainRate computation
    coeffSR = div1_2 * cs2inv * convFac%strainRate

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
        f_preCol(iDir) = state (                               &
          &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      ! Access density and velocity from auxField
      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      inv_rho = 1.0_rk/rho
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = auxField( elemOff + velPos(3) )

      ! Calculate the equilibrium distribution function
      fEq = layout%quantities%pdfEq_ptr( rho = rho,               &
        &                                vel = vel,               &
        &                                QQ = QQ                  )

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric deviatoric second-order tensor of
      ! nonEquilibrium part
      ! the static part cs2 I is usually neglected for weakly compressible flows
      ! however, in current implementation it is considered
      nEqTens = secondMom_minus_cs2_3D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      ! compute strain
      ! magnitude of second-order tensor
      nEqTensMag = sqrt(nEqTens(1)**2 + nEqTens(2)**2 + nEqTens(3)**2 &
        &        + 2.0_rk*(nEqTens(4)**2 + nEqTens(5)**2 + nEqTens(6)**2) )

      ! omega from last time step
      ! convert shear-rate into physical unit because only
      ! non-Newtonian model requies it
      ! physical unit conversion factor is pre-multiplied in coeffSR
      strainRate = coeffSR * omega(iElem) * inv_rho * nEqTensMag

      ! compute shearRate
      shearRate = 2.0_rk * strainRate

      ! compute dynamic viscosity from non-Newtonian Casson model
      ! mu = (k0 + k1 * sqrt(shearRate))**2/shearRate
      viscTerm = ( nNwtn%CS%k0 + nNwtn%CS%k1 * sqrt(shearRate) )
      viscDynaPhy = viscTerm * viscTerm / shearRate
      ! convert to lattice kinematic viscosity
      viscKine(iElem) = ( viscDynaPhy / convFac%viscDyna) * inv_rho

    end do

  end subroutine calcVisc_CS
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate kinematic viscosity from nonNewtonian Carreau-Yasuda model.
  !! $\mu = \mu_\inf + (\mu_0-\mu_\inf)(1+(\lambda*shearRate)*a)^((n-1)/a)$
  !!
  !! Shear rate is computed from strain rate which is computed from
  !! nonEquilibrium PDF which in turn computed from pre-collision PDF
  subroutine calcVisc_CY(nNwtn, viscKine, omega, state, neigh, auxField, &
    & densPos, velPos, nSize, nSolve, nScalars, nAuxScalars, layout, convFac)
    ! --------------------------------------------------------------------------
    !> contains non-Newtonian model parameters loaded from config file
    class(mus_nNwtn_type), intent(in) :: nNwtn
    !> output is physical kinematic viscosity will be overwritten by
    !! non-Netonian model
    real(kind=rk), intent(inout) :: viscKine(:)
    !> Kinematic viscosity omega from last timestep
    real(kind=rk), intent(in) :: omega(:)
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
    !> conversion factor to convert lattice to physical units
    type(mus_convertFac_type), intent(in) :: convFac
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: rho, inv_rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(6), nEqTensMag
    real(kind=rk) :: shearRate, strainRate, v0_vInf, coeffSR
    real(kind=rk) :: viscDynaPhy, viscTerm
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    v0_vInf = nNwtn%CY%visc0 - nNwtn%CY%viscInf

    ! constant coefficients in strainRate computation
    coeffSR = div1_2 * cs2inv * convFac%strainRate

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
        f_preCol(iDir) = state (                               &
          &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      ! Access density and velocity from auxField
      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      inv_rho = 1.0_rk/rho
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = auxField( elemOff + velPos(3) )

      ! Calculate the equilibrium distribution function
      fEq = layout%quantities%pdfEq_ptr( rho = rho,               &
        &                                vel = vel,               &
        &                                QQ = QQ                  )

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric deviatoric second-order tensor of
      ! nonEquilibrium part
      ! the static part cs2 I is usually neglected for weakly compressible flows
      ! however, in current implementation it is considered
      nEqTens = secondMom_minus_cs2_3D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      ! compute strain
      ! magnitude of second-order tensor
      nEqTensMag = sqrt(nEqTens(1)**2 + nEqTens(2)**2 + nEqTens(3)**2 &
        &        + 2.0_rk*(nEqTens(4)**2 + nEqTens(5)**2 + nEqTens(6)**2) )

      ! omega from last time step
      ! convert shear-rate into physical unit because only
      ! non-Newtonian model requies it.
      ! physical unit conversion factor is pre-multiplied in coeffSR
      strainRate = coeffSR * omega(iElem) * inv_rho * nEqTensMag

      ! compute shearRate = 2*strainRate
      shearRate = 2.0_rk * strainRate

      ! compute dynamic viscosity from non-Newtonian Casson model
      ! mu = (k0 + k1 * sqrt(shearRate))**2/shearRate
      viscTerm = 1.0_rk + (nNwtn%CY%lambda*shearRate)**nNwtn%CY%a
      viscDynaPhy = nNwtn%CY%viscInf + v0_vInf &
        &                              * (viscTerm**nNwtn%CY%nMinus1Div_a)
      ! viscKine_L = viscDyna_L / rho
      viscKine(iElem) = (viscDynaPhy / convFac%viscDyna) * inv_rho

    end do

  end subroutine calcVisc_CY
  ! ************************************************************************** !

  ! ************************************************************************** !
  ! Incompressible model                                                       !
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate kinematic viscosity from nonNewtonian power-law model for
  !! incompressible model
  !! $\mu = K shearRate^(n-1)$
  !!
  !! Shear rate is computed from strain rate which is computed from
  !! nonEquilibrium PDF which in turn computed from pre-collision PDF
  subroutine calcVisc_incomp_PL(nNwtn, viscKine, omega, state, neigh,      &
    & auxField, densPos, velPos, nSize, nSolve, nScalars, nAuxScalars, layout, &
    & convFac)
    ! --------------------------------------------------------------------------
    !> contains non-Newtonian model parameters loaded from config file
    class(mus_nNwtn_type), intent(in) :: nNwtn
    !> output is physical kinematic viscosity will be overwritten by
    !! non-Netonian model
    real(kind=rk), intent(inout) :: viscKine(:)
    !> Kinematic viscosity omega from last timestep
    real(kind=rk), intent(in) :: omega(:)
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
    !> conversion factor to convert lattice to physical units
    type(mus_convertFac_type), intent(in) :: convFac
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(6), nEqTensMag
    real(kind=rk) :: shearRate, strainRate, viscDynaPhy, coeffSR
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    ! constant coefficients in strainRate computation
    coeffSR = div1_2 * cs2inv * convFac%strainRate

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
       f_preCol(iDir) = state (                               &
         &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      ! Access density and velocity from auxField
      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = auxField( elemOff + velPos(3) )

      ! Calculate the equilibrium distribution function
      fEq = layout%quantities%pdfEq_ptr( rho = rho, &
        &                                vel = vel, &
        &                                QQ = QQ    )

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric deviatoric second-order tensor of
      ! nonEquilibrium part
      ! the static part cs2 I is usually neglected for weakly compressible flows
      ! however, in current implementation it is considered
      nEqTens = secondMom_minus_cs2_3D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      ! compute strain
      ! magnitude of second-order tensor
      nEqTensMag = sqrt(nEqTens(1)**2 + nEqTens(2)**2 + nEqTens(3)**2 &
        &        + 2.0_rk*(nEqTens(4)**2 + nEqTens(5)**2 + nEqTens(6)**2) )

      ! omega from last time step
      ! convert shear-rate into physical unit because only
      ! non-Newtonian model requies it.
      ! physical unit conversion factor is pre-multiplied in coeffSR
      strainRate = coeffSR * omega(iElem) * rho0Inv * nEqTensMag

      ! compute shearRate
      shearRate = 2.0_rk * strainRate

      ! compute physical dynamic viscosity from non-Newtonian powerlaw model
      viscDynaPhy = (shearRate ** nNwtn%PL%nMinus1) * nNwtn%PL%visc0
      ! viscKine_L = viscDyna_L / rho
      viscKine(iElem) = (viscDynaPhy / convFac%viscDyna) * rho0Inv

    end do

  end subroutine calcVisc_incomp_PL
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate kinematic viscosity from nonNewtonian Casson model for
  !! incompressible model.
  !! $\mu = (k0 + k1 * sqrt(shearRate))^2/shearRate$
  !!
  !! Shear rate is computed from strain rate which is computed from
  !! nonEquilibrium PDF which in turn computed from pre-collision PDF
  subroutine calcVisc_incomp_CS(nNwtn, viscKine, omega, state, neigh,      &
    & auxField, densPos, velPos, nSize, nSolve, nScalars, nAuxScalars, layout, &
    & convFac)
    ! --------------------------------------------------------------------------
    !> contains non-Newtonian model parameters loaded from config file
    class(mus_nNwtn_type), intent(in) :: nNwtn
    !> output is physical kinematic viscosity will be overwritten by
    !! non-Netonian model
    real(kind=rk), intent(inout) :: viscKine(:)
    !> Kinematic viscosity omega from last timestep
    real(kind=rk), intent(in) :: omega(:)
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
    !> conversion factor to convert lattice to physical units
    type(mus_convertFac_type), intent(in) :: convFac
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(6), nEqTensMag
    real(kind=rk) :: shearRate, strainRate, viscTerm, coeffSR, viscDynaPhy
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    ! constant coefficients in strainRate computation
    coeffSR = div1_2 * cs2inv * convFac%strainRate

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
       f_preCol(iDir) = state (                               &
         &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      ! Access density and velocity from auxField
      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = auxField( elemOff + velPos(3) )

      ! Calculate the equilibrium distribution function
      fEq = layout%quantities%pdfEq_ptr( rho = rho, &
        &                                vel = vel, &
        &                                QQ = QQ    )

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric deviatoric second-order tensor of
      ! nonEquilibrium part
      ! the static part cs2 I is usually neglected for weakly compressible flows
      ! however, in current implementation it is considered
      nEqTens = secondMom_minus_cs2_3D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      ! compute strain
      ! magnitude of second-order tensor
      nEqTensMag = sqrt(nEqTens(1)**2 + nEqTens(2)**2 + nEqTens(3)**2 &
        &        + 2.0_rk*(nEqTens(4)**2 + nEqTens(5)**2 + nEqTens(6)**2) )

      ! omega from last time step
      ! convert shear-rate into physical unit because only
      ! non-Newtonian model requies it
      ! physical unit conversion factor is pre-multiplied in coeffSR
      strainRate = coeffSR * omega(iElem) * rho0Inv * nEqTensMag

      ! compute shearRate
      shearRate = 2.0_rk * strainRate

      ! compute dynamic viscosity from non-Newtonian Casson model
      ! mu = (k0 + k1 * sqrt(shearRate))**2/shearRate
      viscTerm = ( nNwtn%CS%k0 + nNwtn%CS%k1 * sqrt(shearRate) )
      viscDynaPhy = viscTerm * viscTerm / shearRate
      ! convert to lattice kinematic viscosity
      viscKine(iElem) = ( viscDynaPhy / convFac%viscDyna) * rho0Inv

    end do

  end subroutine calcVisc_incomp_CS
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Calculate kinematic viscosity from nonNewtonian Carreau-Yasuda model for
  !! incompressible model.
  !! $\mu = \mu_\inf + (\mu_0-\mu_\inf)(1+(\lambda*shearRate)*a)^((n-1)/a)$
  !!
  !! Shear rate is computed from strain rate which is computed from
  !! nonEquilibrium PDF which in turn computed from pre-collision PDF
  subroutine calcVisc_incomp_CY(nNwtn, viscKine, omega, state, neigh,      &
    & auxField, densPos, velPos, nSize, nSolve, nScalars, nAuxScalars, layout, &
    & convFac)
    ! --------------------------------------------------------------------------
    !> contains non-Newtonian model parameters loaded from config file
    class(mus_nNwtn_type), intent(in) :: nNwtn
    !> output is physical kinematic viscosity will be overwritten by
    !! non-Netonian model
    real(kind=rk), intent(inout) :: viscKine(:)
    !> Kinematic viscosity omega from last timestep
    real(kind=rk), intent(in) :: omega(:)
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
    !> conversion factor to convert lattice to physical units
    type(mus_convertFac_type), intent(in) :: convFac
    ! --------------------------------------------------------------------------
    integer :: iElem, iDir, QQ, elemOff
    real(kind=rk) :: rho, vel(3)
    !> precollision PDF
    real(kind=rk) :: f_preCol(layout%fStencil%QQ)
    real(kind=rk) :: fEq(layout%fStencil%QQ), nEq(layout%fStencil%QQ)
    real(kind=rk) :: nEqTens(6), nEqTensMag
    real(kind=rk) :: shearRate, strainRate, v0_vInf, coeffSR
    real(kind=rk) :: viscDynaPhy, viscTerm
    ! --------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    v0_vInf = nNwtn%CY%visc0 - nNwtn%CY%viscInf

    ! constant coefficients in strainRate computation
    coeffSR = div1_2 * cs2inv * convFac%strainRate

    do iElem = 1, nSolve
      ! Get pre-collisiton PDF
      do iDir = 1, QQ
       f_preCol(iDir) = state (                               &
         &  neigh((idir-1)* nsize+ ielem)+( 1-1)* qq+ nscalars*0)
      end do

      ! Access density and velocity from auxField
      elemOff = (iElem-1)*nAuxScalars
      ! density
      rho = auxField( elemOff + densPos)
      ! velocity
      vel(1) = auxField( elemOff + velPos(1) )
      vel(2) = auxField( elemOff + velPos(2) )
      vel(3) = auxField( elemOff + velPos(3) )

      ! Calculate the equilibrium distribution function
      fEq = layout%quantities%pdfEq_ptr( rho = rho, &
        &                                vel = vel, &
        &                                QQ = QQ    )

      ! Calculate the non-equilibrium part
      nEq(:) = f_preCol(:) - fEq(:)

      ! Now calculate the symmetric deviatoric second-order tensor of
      ! nonEquilibrium part
      ! the static part cs2 I is usually neglected for weakly compressible flows
      ! however, in current implementation it is considered
      nEqTens = secondMom_minus_cs2_3D(layout%fStencil%cxcx, nEq, layout%fStencil%QQ)

      ! compute strain
      ! magnitude of second-order tensor
      nEqTensMag = sqrt(nEqTens(1)**2 + nEqTens(2)**2 + nEqTens(3)**2 &
        &        + 2.0_rk*(nEqTens(4)**2 + nEqTens(5)**2 + nEqTens(6)**2) )

      ! omega from last time step
      ! convert shear-rate into physical unit because only
      ! non-Newtonian model requies it.
      ! physical unit conversion factor is pre-multiplied in coeffSR
      strainRate = coeffSR * omega(iElem) * rho0Inv * nEqTensMag

      ! compute shearRate = 2*strainRate
      shearRate = 2.0_rk * strainRate

      ! compute dynamic viscosity from non-Newtonian Casson model
      ! mu = (k0 + k1 * sqrt(shearRate))**2/shearRate
      viscTerm = 1.0_rk + (nNwtn%CY%lambda*shearRate)**nNwtn%CY%a
      viscDynaPhy = nNwtn%CY%viscInf + v0_vInf &
        &                              * (viscTerm**nNwtn%CY%nMinus1Div_a)
      ! viscKine_L = viscDyna_L / rho
      viscKine(iElem) = (viscDynaPhy / convFac%viscDyna) * rho0Inv

    end do

  end subroutine calcVisc_incomp_CY
  ! ************************************************************************** !

end module mus_nonNewtonian_module
! **************************************************************************** !
