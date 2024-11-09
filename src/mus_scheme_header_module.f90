! Copyright (c) 2012-2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2024 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
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
! **************************************************************************** !
!> This module contains scheme property type and module related to scheme prop
module mus_scheme_header_module

  ! include treelm modules
  use env_module,                    only: LabelLen, rk
  use tem_tools_module,              only: tem_horizontalSpacer
  use tem_aux_module,                only: tem_abort
  use tem_logging_module,            only: logUnit

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_NonExistent
  use aot_table_module, only: aot_table_open, aot_table_close
  use aot_out_module,   only: aot_out_type, aot_out_val, aot_out_open_table,   &
    &                         aot_out_close_table

  implicit none

  private

  public :: mus_scheme_header_type, mus_scheme_header_out
  public :: mus_load_scheme_header
  public :: mus_relaxation_header_type


  !> Datatype containing additional options for the relaxation like variant and
  !! other variant specific parameters
  type mus_relaxation_header_type
    !> Varaint name of the relaxation. Set to "default" to select default
    !! relaxation
    character(len=labelLen) :: variant
    !> Addtional information to load for regularited bgk like
    !! "regularized", "recursive_regularited" and "hybrid_recursive_regularized"
    !! variant.
    real(kind=rk) :: regularization_omega

    !> todo move omega_Cum, omega_Lim, DRT_tauN, lambda from mus_fluid_type
    !! to here
  end type mus_relaxation_header_type


  !> Datatype containing information to identify the scheme
  !!
  !! Combination of scheme kind, relaxation and layout%stencilKind
  !! are used to choose the correct compute kernel for the
  !! scheme
  !!
  !!> | type | options |
  !!> |:-----------------|:--------------|
  !!> | kind | **fluid** (default)             |
  !!> |      | **fluid_incompressible**        |
  !!> |      | **isotherm_acEq**               |
  !!> |      | **multispecies_gas**            |
  !!> |      | **multispecies_liquid**         |
  !!> |      | **nernst_planck**               |
  !!> |      | **passive_scalar**              |
  !!> |      | **poisson**                     |
  !!> |      | **poisson_boltzmann_linear**    |
  !!> |      | **poisson_boltzmann_nonlinear** |
  !!> |------|---------------------------------|
  !!> | layout | **d2q9**           |
  !!> |        | **d3q19** (default)|
  !!> |        | **d3q27**          |
  !!> |        | _d1q3_             |
  !!> |        | _d2q5_             |
  !!> |        | _d3q6_             |
  !!> |        | _d3q7_             |
  !!> |        | _d3q13_            |
  !!> |        | _d3q15_            |
  !!> |        | _flekkoy_          |
  !!> |--------|--------------------|
  !!> | relaxation | **bgk**  (default)       |
  !!> |            | **mrt**                  |
  !!> |            | **trt**                  |
  !!> |            | **bgk_withthermodynfac** |
  !!> |            | **mrt_withthermodynfac** |
  !!> |            | _cumulant_               |
  !!> |            | _cascaded_               |
  !!> |            | _vec_fma_                |
  !!> |            | _test_                   |
  !!> |            | _bgk_noFluid_            |
  !!> |------------|--------------------------|
  !!> | variant for bgk relaxation | **standard**  (default)       |
  !!> |                            | **improved**                  |
  !!> |                            | **block**                     |
  !!> |                            | **mrt**                       |
  !!> |------------------------------------------------------------|
  type mus_scheme_header_type
    !> scheme kind, Ex: fluid, fluid_incompressible, multispecies_gas,
    !! multispecies_liquid, poisson, poisson_boltzmann_linear,
    !! poisson_boltzmann_nonlinear, nernst_planck, isotherm_acEq
    character(len=labelLen) :: kind
    !> scheme layout, Ex: d3q19
    character(len=labelLen) :: layout
    !> scheme relaxation type Ex: BGK, MRT, bgk_pl, bgk_cy, bgk_cs...
    character(len=labelLen) :: relaxation
    !> Variant and additional options for a relaxation
    type(mus_relaxation_header_type) :: relaxHeader
  end type mus_scheme_header_type

contains

  ! ************************************************************************** !
  !> load scheme header info from lua file identify table or from scheme table
  !! or from config
  !!
  !! Load scheme label, kind, layoutKind and relaxation
  !!```lua
  !! identify = { kind = 'simType',
  !!              layout = 'stencilLayout',
  !!              relaxation = 'relaxationType' }
  !!```
  !! For a possible selection of the other parameters
  !! - simType: fluid, fluid_incompressible, multispecies_liquid
  !! - [[mus_scheme_layout_module]]: d2q9, d3q19, ...
  !! - relaxationType: bgk, mrt, ...
  !!
  subroutine mus_load_scheme_header(me, conf, parent, scaling)
    !---------------------------------------------------------------------------
    !> returns scheme identify information
    type(mus_scheme_header_type), intent(out) :: me
    type(flu_State) :: conf !< flu state
    !> parent handle if scheme table is defined
    integer, intent(in), optional :: parent
    !> scaling: diffusive or acoustic?
    character(len=*), intent(in) :: scaling
    ! --------------------------------------------------------------------------
    integer :: thandle !< handle for scheme identify table
    integer :: relax_handle
    integer :: iError
    ! --------------------------------------------------------------------------

    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1), '(A)') 'Loading Scheme identify table: '
    !default values
    me%kind = 'fluid'
    me%layout = 'd3q19'
    me%relaxation = 'bgk'

    call aot_table_open( L       = conf,      &
      &                  parent  = parent,    &
      &                  thandle = thandle,   &
      &                  key     = 'identify' )

    if (thandle > 0) then
      ! get schemekind
      call aot_get_val( L       = conf,    &
        &               thandle = thandle, &
        &               key     = 'kind',  &
        &               val     = me%kind, &
        &               default = 'fluid', &
        &               ErrCode = iError   )

      ! get layoutkind
      call aot_get_val( L       = conf,      &
        &               thandle = thandle,   &
        &               key     = 'layout',  &
        &               val     = me%layout, &
        &               default = 'd3q19',   &
        &               ErrCode = iError     )

      ! Load relaxation as table to load additional information for relaxation.
      ! if not a table then variant is set to default.
      call aot_table_open( L       = conf,         &
        &                  parent  = thandle,      &
        &                  thandle = relax_handle, &
        &                  key     = 'relaxation'  )

      if (relax_handle ==  0) then
        ! get relaxation
        call aot_get_val( L       = conf,          &
          &               thandle = thandle,       &
          &               key     = 'relaxation',  &
          &               val     = me%relaxation, &
          &               default = 'bgk',         &
          &               ErrCode = iError         )
        me%relaxHeader%variant = 'standard'
      else ! load relaxation options from a table
        ! get relaxation
        call aot_get_val( L       = conf,          &
          &               thandle = relax_handle,  &
          &               key     = 'name',        &
          &               val     = me%relaxation, &
          &               default = 'bgk',         &
          &               ErrCode = iError         )
        call load_relaxation_header( me      = me%relaxHeader, &
          &                          conf    = conf,           &
          &                          thandle = relax_handle    )
      end if
      call aot_table_close(L=conf, thandle=relax_handle)
    else
      write(logUnit(1), '(A)') 'Scheme Identify table not defined.'
      write(logUnit(1), '(A)') 'Setting default values for scheme..'
    end if

    call aot_table_close(L=conf, thandle=thandle)

    write(logUnit(1), '(A)') 'kind: '// trim(me%kind)
    write(logUnit(1), '(A)') 'Layout: '// trim(me%layout)
    write(logUnit(1), '(A)') 'relaxation: '// trim(me%relaxation)
    write(logUnit(1), '(A)') '  variant: ' //           &
      &                      trim(me%relaxHeader%variant)
    call tem_horizontalSpacer(fUnit = logUnit(1))

    ! Both multispeciees and poisson equation must have diffusive scaling
    ! since diffusive scaling is used to recover macroscopic equations
    ! from asymptotic analysis
    select case(trim(me%kind))
    case ('fluid', 'fluid_incompressible', 'isotherm_acEq' )
      if (trim(scaling) /= 'acoustic') then
         call tem_abort('ERROR: Choose scaling = "acoustic" for ' &
           &            // trim(me%kind))
      end if
    case ( 'multispecies_gas', 'multispecies_liquid', 'nernst_planck', &
      &    'passive_scalar, ''poisson', 'poisson_boltzmann_linear',    &
      &    'poisson_boltzmann_nonlinear'                               )
      if(trim(scaling) /= 'diffusive') then
         call tem_abort('ERROR: Choose scaling = "diffusive" for ' &
           &            // trim(me%kind))
      end if
    end select

  end subroutine mus_load_scheme_header
  ! ************************************************************************** !

  ! ***************************************************************************!
  !> Load relaxation options from a table
  subroutine load_relaxation_header(me, conf, thandle)
    ! --------------------------------------------------------------------------
    type(mus_relaxation_header_type), intent(out) :: me
    type(flu_State) :: conf !< flu state
    !> relaxation handle
    integer, intent(in) :: thandle
    ! --------------------------------------------------------------------------
    integer :: iError
    ! --------------------------------------------------------------------------
    call aot_get_val( L       = conf,       &
      &               thandle = thandle,    &
      &               key     = 'variant',  &
      &               val     = me%variant, &
      &               default = 'standard', &
      &               ErrCode = iError      )


    select case (trim(me%variant))
    case ('regularized', 'recursive_regularied', 'hybrid_recursive_regularized')
      call aot_get_val( L       = conf,                    &
        &               thandle = thandle,                 &
        &               key     = 'regularization_omega',  &
        &               val     = me%regularization_omega, &
        &               ErrCode = iError                   )
      if (btest(iError, aoterr_NonExistent)) then
        call tem_abort('Error: regularization_omega is not specified for ' &
          &            // 'regularized variant')
      end if
      write(logUnit(1),'(A, F10.7)') '  regularization_omega: ', &
        &                            me%regularization_omega

    end select

  end subroutine load_relaxation_header
  ! ***************************************************************************!

  ! ************************************************************************** !
  !> Dumps scheme header
  subroutine mus_scheme_header_out(me, conf)
  ! ----------------------------------------------------------------------------
    !> returns scheme identify information
    type(mus_scheme_header_type), intent(in) :: me
    type(aot_out_type) :: conf
    ! --------------------------------------------------------------------------

    ! the label does not have to be outputted.
    ! because it is the name of the outputted varSys
    call aot_out_open_table(put_conf = conf, tname = 'identify')
    call aot_out_val( put_conf = conf,          &
      &               vname    = 'kind',        &
      &               val      = trim( me%kind ))
    call aot_out_val( put_conf = conf,                &
      &               vname    = 'relaxation',        &
      &               val      = trim( me%relaxation ))
    call aot_out_val( put_conf = conf,            &
      &               vname    = 'layout',        &
      &               val      = trim( me%layout ))
    call aot_out_close_table(put_conf = conf)
  end subroutine mus_scheme_header_out
  ! ************************************************************************** !

end module mus_scheme_header_module
! **************************************************************************** !
