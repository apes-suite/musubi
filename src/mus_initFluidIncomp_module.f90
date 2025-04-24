! Copyright (c) 2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2024 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2013-2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!> This module contains routines which initiliaze advection relaxation and
!! flow field for lbm incompressible model.
module mus_initFluidIncomp_module
  ! include treelm modules
  use env_module,         only: labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  ! include musubi modules
  use mus_bgk_module,   only: mus_advRel_kCFD_rBGK_vStdNoOpt_l
  use mus_d3q19_module, only: mus_advRel_kFluidIncomp_rBGK_vStd_lD3Q19, &
    &                         mus_advRel_kFluidIncomp_rTRT_vStd_lD3Q19, &
    &                         mus_advRel_kFluidIncompGNS_rBGK_vStd_lD3Q19
  use mus_d3q27_module, only: mus_advRel_kCFD_rBGK_vStd_lD3Q27
  use mus_d2q9_module,  only: mus_advRel_kFluidIncomp_rBGK_vStd_lD2Q9, &
    &                         mus_advRel_kFluidIncomp_rMRT_vStd_lD2Q9
  use mus_mrt_d3q19_module, only: mus_advRel_kCFD_rMRT_vStdNoOpt_l,            &
    &                             mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q19,    &
    &                             mus_advRel_kFluidIncomp_rMRT_vStdNoOpt_lD3Q19
  use mus_mrt_d3q27_module, only: mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q27, &
    &                             mus_advRel_kCFD_rMRT_vStdNoOpt_lD3Q27
  use mus_scheme_type_module, only: kernel

  implicit none

  private

  public :: mus_init_advRel_fluidIncomp
  public :: mus_init_advRel_fluidIncomp_GNS

contains

  ! ************************************************************************** !
  !> Initialize the relaxation model for lbm incompressible model
  subroutine mus_init_advRel_fluidIncomp(relaxation, variant, layout, compute)
    ! --------------------------------------------------------------------------
    character(len=labelLen), intent(inout) :: relaxation
    character(len=labelLen), intent(in) :: variant
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! --------------------------------------------------------------------------
    write(logUnit(1),*) 'Choosing fluid_incompressible relaxation model: ' &
      &                 // trim(relaxation) // ' for layout ' // trim(layout)

    select case (trim(relaxation))
    case ('bgk')
      call mus_init_advRel_fluidIncomp_bgk( variant, layout, compute )

    case ('mrt')
      call mus_init_advRel_fluidIncomp_mrt( variant, layout, compute )

    case ('trt')
      select case (trim(variant))
      case ('standard')
        select case (trim(layout))
        case ('d3q19')
          compute => mus_advRel_kFluidIncomp_rTRT_vStd_lD3Q19
        case default
          call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
            &            // 'relaxation "'//trim(relaxation)//'"')
        end select
      case default
        call tem_abort('Unsupported variant "'//trim(variant)//'" for ' &
          &            // 'relaxation "trt"')
      end select

    case default
      call tem_abort('Unsupported relaxation "'//trim(relaxation)//'" for ' &
        &            // 'kind "fluid_incompressible"')
    end select

  end subroutine mus_init_advRel_fluidIncomp
  ! ************************************************************************** !
  ! ************************************************************************** !
  !> Assigning compute kernel routine by scheme relaxation type for fluid GNS 
  !! kind for unresolved LBM-DEM particulate flow simulations
  subroutine mus_init_advRel_fluidIncomp_GNS( relaxation, layout, compute )
    ! --------------------------------------------------------------------------
    character(len=labelLen), intent(inout) :: relaxation
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! --------------------------------------------------------------------------

    write(logUnit(1),*) 'Choosing fluid relaxation model: '                 &
      &                 // trim(relaxation) // ' for layout ' // trim(layout)

    select case (trim(relaxation))

    case ('bgk')
      select case (trim(layout))
      case ('d3q19')
        compute => mus_advRel_kFluidIncompGNS_rBGK_vStd_lD3Q19
      case default
        write(logUnit(1),*) 'ERROR: layout not supported!'
        write(logUnit(1),*) 'Currently only d3q19 and d2q9 layouts are supported for fluid_GNS'
        call tem_abort()
      end select

    case default
      write(logUnit(1),*) 'Relaxation '//trim(relaxation)//' is not supported!'
      write(logUnit(1),*) 'Currently only relaxation bgk is supported for fluid_GNS'
      call tem_abort()
    end select

  end subroutine mus_init_advRel_fluidIncomp_GNS
  ! **************************************************************************** !

  ! ************************************************************************** !
  !> This routine assigns compute routine for bgk relaxation
  !!
  !! Supported variants are:
  !!   * standard        - Optimized routines for specifc layouts.
  !!   * standard_no_opt - no optimized routines for any layouts.
  subroutine mus_init_advRel_fluidIncomp_bgk(variant, layout, compute)
    ! --------------------------------------------------------------------------
    character(len=labelLen), intent(in) :: variant
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! --------------------------------------------------------------------------
    select case (trim(variant))
    case ('standard')
      select case (trim(layout))
      case ('d3q27')
        compute => mus_advRel_kCFD_rBGK_vStd_lD3Q27
      case ('d3q19')
        compute => mus_advRel_kFluidIncomp_rBGK_vStd_lD3Q19
      case ('d2q9')
        compute => mus_advRel_kFluidIncomp_rBGK_vStd_lD2Q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for '         &
          &            //' relaxation "bgk" for variant "'//trim(variant)//'"')
      end select

    case ('standard_no_opt')
      compute => mus_advRel_kCFD_rBGK_vStdNoOpt_l

    case default
      call tem_abort('Unsupported variant "'//trim(variant)//'" for ' &
        &            //'relaxation "bgk"')
    end select

  end subroutine mus_init_advRel_fluidIncomp_bgk
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine assigns compute routine for mrt relaxation
  !!
  !! Supported variants are:
  !!   * standard        - Optimized routines for specifc layouts.
  !!   * standard_no_opt - no optimized routines for any layouts.
  !!   * bgk             - Uses no optimized routine but
  !!                       in mus_mrtRelaxation_module all relaxation parameters
  !!                       are set to same omega to recover bgk.
  subroutine mus_init_advRel_fluidIncomp_mrt(variant, layout, compute)
    ! --------------------------------------------------------------------------
    character(len=labelLen), intent(in) :: variant
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! --------------------------------------------------------------------------
    select case (trim(variant))
    case ('standard')
      select case (trim(layout))
      case ('d3q27')
        compute => mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q27
      case ('d3q19')
        compute => mus_advRel_kFluidIncomp_rMRT_vStd_lD3Q19
      case ('d2q9')
        compute => mus_advRel_kFluidIncomp_rMRT_vStd_lD2Q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for '        &
          &            //'relaxation "mrt" for variant "'//trim(variant)//'"')
      end select

    case ('standard_no_opt')
      select case (trim(layout))
      case ('d3q27')
        compute => mus_advRel_kCFD_rMRT_vStdNoOpt_lD3Q27
      case ('d3q19')
        compute => mus_advRel_kFluidIncomp_rMRT_vStdNoOpt_lD3Q19
      case default
        compute => mus_advRel_kCFD_rMRT_vStdNoOpt_l
      end select

    case ('bgk')
      compute => mus_advRel_kCFD_rMRT_vStdNoOpt_l

    case default
      call tem_abort('Unsupported variant "'//trim(variant)//'" for ' &
        &            //'relaxation "mrt"')
    end select

  end subroutine mus_init_advRel_fluidIncomp_mrt
  ! ************************************************************************** !

end module mus_initFluidIncomp_module
