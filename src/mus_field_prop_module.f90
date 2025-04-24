! Copyright (c) 2012-2013, 2016-2017, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Philipp Otte <otte@mathcces.rwth-aachen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! See copyright notice in he COPYRIGHT file.
! **************************************************************************** !
!> This module contains mus_field_prop_type and modules related
!! to fiels properties.
!!
!! \author Kannan Masilamani, Simon Zimny
!!
module mus_field_prop_module

  ! include treelm modules
  use env_module,         only: rk, labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  ! include aotus modules
  use flu_binding,     only: flu_State
  use aot_out_module,  only: aot_out_type

  ! include musubi modules
  use mus_poisson_module,       only: mus_load_poisson, &
    &                                 mus_poisson_type
  use mus_scheme_header_module, only: mus_scheme_header_type
  use mus_fluid_module,         only: mus_fluid_type, mus_load_fluid, &
    &                                 mus_fluid_save2lua
  use mus_species_module,       only: mus_species_type, mus_load_species, &
    &                                 mus_species_out
  use mus_physics_module,       only: mus_physics_type

  implicit none

  private

  public :: mus_field_prop_type
  public :: mus_load_field_prop
  public :: mus_field_prop_out

  !> This type contains parameter needed for field
  type mus_field_prop_type
    !> contains fluid information
    type( mus_fluid_type ) :: fluid
    !> contains species information
    type( mus_species_type ) :: species
    !> Contains information for poisson equation
    type( mus_poisson_type) :: poisson
  end type mus_field_prop_type


contains


! **************************************************************************** !
  !> load fluid properties like fluid and species table from
  !! lua file based on the scheme kind
  subroutine mus_load_field_prop( me, conf, parent, minLevel,    &
    &                             schemeHeader, nFields, physics, &
    &                             cs_lattice )
    ! --------------------------------------------------------------------------
    !> field property type
    type( mus_field_prop_type ), intent(out) :: me
    !> number of fields defined in lua file
    integer, intent(in) :: nFields
    !> flu state
    type( flu_State ), intent(inout) :: conf
    !> parent lua handle
    integer, intent(in), optional :: parent
    integer, intent(in) :: minLevel
    !> identifier of the scheme
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    !> physics type to convert physics to lattice unit or vice versa
    type( mus_physics_type ), intent(in) :: physics
    !> lattice speed of sound calculated for defined stencil layout
    real(kind=rk), intent(in) :: cs_lattice
    ! --------------------------------------------------------------------------
    ! load fluid info
    select case( trim(schemeHeader%kind) )
      case('fluid', 'fluid_incompressible', 'fluid_GNS', &
          & 'fluid_incompressible_GNS','isotherm_acEq'   )
        write(logUnit(1),"(A)") ' Loading the fluid properties.'
        call mus_load_fluid( me           = me%fluid,    &
          &                  conf         = conf,        &
          &                  parent       = parent,      &
          &                  minLevel     = minLevel,    &
          &                  physics      = physics,     &
          &                  schemeHeader = schemeHeader )
      case('poisson', 'poisson_boltzmann_linear', 'poisson_boltzmann_nonlinear')
        write(logUnit(1),"(A)") ' Loading properties for poisson Eq:'
        call mus_load_poisson( me         = me%poisson,             &
          &                    conf       = conf,                   &
          &                    parent     = parent,                 &
          &                    minLevel   = minLevel,               &
          &                    cs_lattice = cs_lattice,             &
          &                    physics    = physics,                &
          &                    schemeKind = trim(schemeHeader%kind) )
      case('passive_scalar')
        write(logUnit(1),"(A)") ' Loading properties for passive scalar:'
        call mus_load_species( me =  me%species, conf = conf, parent = parent, &
          &                    minLevel = minLevel, nFields = nFields,       &
          &                    physics = physics, cs_lattice=cs_lattice )
      case('multispecies_gas','multispecies_liquid','nernst_planck')
        write(logUnit(1),"(A)") ' Loading properties for multispecies:'
        call mus_load_species( me =  me%species, conf = conf, parent = parent, &
          &                    minLevel = minLevel, nFields = nFields,         &
          &                    physics = physics, cs_lattice=cs_lattice )
      case default
        write(logUnit(1),*) 'The selected scheme kind is unknown '//           &
          &             trim(schemeHeader%kind)
        call tem_abort()
    end select

  end subroutine mus_load_field_prop
! **************************************************************************** !


! **************************************************************************** !
  !> write field prop into a lua file
  !!
  subroutine mus_field_prop_out( me, conf, schemeHeader )
    ! --------------------------------------------------------------------------
    !> single field type
    type( mus_field_prop_type ), intent(in) :: me
    !> identifier of the scheme
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    !> aotus out type
    type(aot_out_type), intent(inout) :: conf
    ! --------------------------------------------------------------------------

    select case( trim(schemeHeader%kind) )
      case( 'fluid', 'fluid_incompressible' )
        call mus_fluid_save2lua( me = me%fluid, conf = conf )
      case('passive_scalar', 'multispecies_gas', &
        &  'multispecies_liquid')
        call mus_species_out( me = me%species, conf = conf )
    end select

  end subroutine mus_field_prop_out
! **************************************************************************** !

end module mus_field_prop_module
! **************************************************************************** !
