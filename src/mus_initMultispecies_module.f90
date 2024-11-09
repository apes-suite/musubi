! Copyright (c) 2013-2014, 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! ****************************************************************************** !
!> This module contains routines which initiliaze advection relaxation and
!! flow field for multispecies lbm gas model and liquid model.
module mus_initMultispecies_module

  ! include treelm modules
  use env_module,         only: labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  ! include musubi modules
  use mus_MSGas_module, only: bgk_advRel_MSGas_generic, &
    &                         bgk_advRel_d3q19f3_MSGas
  use mus_MSLiquid_module, only: bgk_advRel_MSLiquid_generic,         &
    &                            bgk_advRel_MSLiquid_generic_WTDF,    &
    &                            bgk_advRel_d3q19f3_MSLiquid,         &
    &                            bgk_advRel_d3q19f3_MSLiquid_WTDF,    &
    &                            bgk_forcing_advRel_MSLiquid_generic, &
    &                            mrt_advRel_MSLiquid_generic,         &
    &                            mrt_advRel_MSLiquid_generic_WTDF,    &
    &                            mrt_advRel_d3q19f3_MSLiquid,         &
    &                            mrt_advRel_d3q19f3_MSLiquid_WTDF
  use mus_scheme_type_module, only: kernel

  implicit none

  private

  public :: mus_init_advRel_multispecies_gas
  public :: mus_init_advRel_multispecies_liquid

contains

! ****************************************************************************** !
  !> Initialize the relaxation model for multispecies gas model
  subroutine mus_init_advRel_multispecies_gas( relaxation, layout, nFields, &
    &                                          compute )
    ! ---------------------------------------------------------------------------
    character(len=labelLen), intent(in) :: relaxation
    character(len=labelLen), intent(in) :: layout
    integer,                 intent(in) :: nFields
    procedure( kernel ), pointer, intent(out) :: compute
    ! ---------------------------------------------------------------------------

    select case ( trim(relaxation) )
    case ('bgk')
      select case ( trim(layout) )
      case ('d3q19')
        if ( nFields == 3 ) then
          write(logUnit(1),*) 'Choosen optimized d3q19f3 multispecies '//      &
            &                 'gas kernel'
          compute => bgk_advRel_d3q19f3_MSGas
        else
          write(logUnit(1),*) 'Choosen testing multispecies gas kernel'
          compute => bgk_advRel_MSGas_generic
        end if
      case ('d2q9')
        write(logUnit(1),*) 'Choosen testing multispecies gas kernel'
        compute => bgk_advRel_MSGas_generic
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select ! layout
    case default
      write(logUnit(1),*) 'Relaxation '//trim(relaxation)//' is not supported!'
      call tem_abort()
    end select ! relaxation

  end subroutine mus_init_advRel_multispecies_gas
! ****************************************************************************** !

! ****************************************************************************** !
  !> Initialize the relaxation model for multispecies liquid model
  subroutine mus_init_advRel_multispecies_liquid( relaxation, layout, nFields, &
    &                                             compute                      )
    ! ---------------------------------------------------------------------------
    character(len=labelLen), intent(in) :: relaxation
    character(len=labelLen), intent(in) :: layout
    integer,                 intent(in) :: nFields
    procedure( kernel ), pointer, intent(out) :: compute
    ! ---------------------------------------------------------------------------

    select case ( trim(relaxation) )
    case ('bgk')
      if ( (trim(layout) == 'd3q19') .and. (nFields == 3) ) then
        write(logUnit(1),*) 'Choosen optimized d3q19f3 multispecies ' &
          &                 // 'liquid kernel'
        compute => bgk_advRel_d3q19f3_MSLiquid
      else
        write(logUnit(1),*) 'Choosen generic generic multispecies ' &
          &                 // 'liquid kernel'
        compute => bgk_advRel_MSLiquid_generic
      end if

    case( 'mrt' )
      if ( (trim(layout) == 'd3q19') .and. (nFields == 3) ) then
        write(logUnit(1),*) 'Choosen optimized mrt d3q19f3 multispecies '//  &
          &                 'liquid kernel'
        compute => mrt_advRel_d3q19f3_MSLiquid
      else
        write(logUnit(1),*) 'Choosen generic mrt multispecies liquid '// &
          &                 'kernel'
        compute => mrt_advRel_MSLiquid_generic
      end if

    case ('bgk_forcing')
      write(logUnit(1),*) 'Choosen generic bgk forcing multispecies '//  &
        &                 'liquid kernel'
      compute => bgk_forcing_advRel_MSLiquid_generic

    case ('bgk_withthermodynfac')
      if ( (trim(layout) == 'd3q19') .and. (nFields == 3) ) then
        write(logUnit(1),*) 'Choosen optimized d3q19f3 multispecies ' &
          &                 // 'liquid kernel'
        compute => bgk_advRel_d3q19f3_MSLiquid_WTDF
      else
        write(logUnit(1),*) 'Choosen generic bgk with thermodynamic '//  &
          &             'factor multispecies liquid kernel'
        compute => bgk_advRel_MSLiquid_generic_WTDF
      end if

    case ('mrt_withthermodynfac')
      if ( (trim(layout) == 'd3q19') .and. (nFields == 3) ) then
        write(logUnit(1),*) 'Choosen optimized d3q19f3 multispecies ' &
          &                 // 'liquid kernel'
        compute => mrt_advRel_d3q19f3_MSLiquid_WTDF
      else
        write(logUnit(1),*) 'Choosen generic mrt multispecies liquid '// &
          &                 'kernel with thermodynamic factor'
        compute => mrt_advRel_MSLiquid_generic_WTDF
      end if

    case default
      write(logUnit(1),*) 'Relaxation '//trim(relaxation)//' is not supported!'
      call tem_abort()

    end select ! relaxation

  end subroutine mus_init_advRel_multispecies_liquid
! ****************************************************************************** !

end module mus_initMultispecies_module
