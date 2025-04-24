! Copyright (c) 2013, 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2024 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2021 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!! flow field for lbm model.
module mus_initFluid_module
  use env_module,         only: labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  use mus_bgk_module,       only: mus_advRel_kCFD_rBGK_vStdNoOpt_l
  use mus_compute_cumulant_module,  only: cumulant_d3q27, cascaded_d3q27,  &
    &                                     cumulant_d3q27_extended_generic, &
    &                                     cumulant_d3q27_extended_fast
  use mus_d3q27_module,     only: mus_advRel_kCFD_rBGK_vStd_lD3Q27,        &
    &                             mus_advRel_kFluid_rTRT_vStd_lD3Q27,      &
    &                             mus_advRel_kFluid_rBGK_vImproved_lD3Q27, &
    &                             bgk_Regularized_d3q27,                   &
    &                             bgk_RecursiveRegularized_d3q27,          &
    &                             bgk_ProjectedRecursiveRegularized_d3q27, &
    &                             bgk_HybridRecursiveRegularized_d3q27,    &
    &                             bgk_DualRelaxationTime_RR_d3q27,         &
    &                             bgk_HybridRecursiveRegularizedCorr_d3q27
  use mus_d3q19_module,     only: mus_advRel_kFluid_rBGK_vStd_lD3Q19,      &
    &                             mus_advRel_kFluid_rBGK_vBlock_lD3Q19,    &
    &                             mus_advRel_kFluid_rTRT_vStd_lD3Q19,      &
    &                             bgk_Regularized_d3q19,                   &
    &                             bgk_RecursiveRegularized_d3q19,          &
    &                             bgk_ProjectedRecursiveRegularized_d3q19, &
    &                             bgk_HybridRecursiveRegularized_d3q19,    &
    &                             bgk_DualRelaxationTime_RR_d3q19,         &
    &                             bgk_HybridRecursiveRegularizedCorr_d3q19,&
    &                             bgk_advRel_d3q19_GNS
  use mus_d2q9_module,      only: mus_advRel_kFluid_rMRT_vStd_lD2Q9,       &
    &                             mus_advRel_kFluid_rBGK_vStd_lD2Q9,       &
    &                             mus_advRel_kFluid_rBGK_vImproved_lD2Q9,  &
    &                             bgk_Regularized_d2q9,                    &
    &                             bgk_RecursiveRegularized_d2q9,           &
    &                             bgk_ProjectedRecursiveRegularized_d2q9,  &
    &                             bgk_HybridRecursiveRegularized_d2q9,     &
    &                             bgk_DualRelaxationTime_RR_d2q9,          &
    &                             bgk_HybridRecursiveRegularizedCorr_d2q9, &
    &                             mus_advRel_kFluidGNS_rBGK_vStd_lD2Q9                    
  use mus_mrt_d3q19_module, only: mus_advRel_kFluid_rMRT_vStd_lD3Q19,      &
    &                             mus_advRel_kFluid_rMRT_vStdNoOpt_lD3Q19, &
    &                             mus_advRel_kCFD_rMRT_vStdNoOpt_l
  use mus_mrt_d3q27_module, only: mus_advRel_kFluid_rMRT_vStd_lD3Q27, &
    &                             mus_advRel_kCFD_rMRT_vStdNoOpt_lD3Q27
  use mus_test_module,      only: vec_fma
  use mus_scheme_type_module, only: kernel

  implicit none

  private

  public :: mus_init_advRel_fluid
  public :: mus_init_advRel_fluid_GNS

contains

  ! ************************************************************************** !
  !> Assigning compute kernel routine by scheme relaxation type for fluid kind.
  !!
  subroutine mus_init_advRel_fluid(relaxation, variant, layout, compute)
    ! --------------------------------------------------------------------------
    character(len=labelLen), intent(in) :: relaxation
    character(len=labelLen), intent(in) :: variant
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! --------------------------------------------------------------------------

    write(logUnit(1),*) 'Choosing fluid relaxation model: '                 &
      &                 // trim(relaxation) // ' for layout ' // trim(layout)

    select case (trim(relaxation))
    case ('bgk')
      call mus_init_advRel_fluid_bgk(variant, layout, compute)

    case ('mrt')
      call mus_init_advRel_fluid_mrt(variant, layout, compute)

    case ('trt')
      call mus_init_advRel_fluid_trt(variant, layout, compute)

    case ('cumulant')
      select case (trim(layout))
      case ('d3q27')
        compute => cumulant_d3q27
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case ('cumulant_extended')
      select case (trim(layout))
      case ('d3q27')
        compute => cumulant_d3q27_extended_fast
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case ('cumulant_extended_generic')
      select case (trim(layout))
      case ('d3q27')
        compute => cumulant_d3q27_extended_generic
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case ('cascaded')
      select case (trim(layout))
      case ('d3q27')
        compute => cascaded_d3q27
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case ('vec_fma', 'test')
      compute => vec_fma

    case ('hrr_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_HybridRecursiveRegularized_d3q27
      case ('d3q19')
        compute => bgk_HybridRecursiveRegularized_d3q19
      case ('d2q9')
        compute => bgk_HybridRecursiveRegularized_d2q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case ('hrr_bgk_corrected', 'prr_bgk_corrected', 'rr_bgk_corrected')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_HybridRecursiveRegularizedCorr_d3q27
      case ('d3q19')
        compute => bgk_HybridRecursiveRegularizedCorr_d3q19
      case ('d2q9')
        compute => bgk_HybridRecursiveRegularizedCorr_d2q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case ('drt_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_DualRelaxationTime_RR_d3q27
      case ('d3q19')
        compute => bgk_DualRelaxationTime_RR_d3q19
      case ('d2q9')
        compute => bgk_DualRelaxationTime_RR_d2q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case ('rr_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_RecursiveRegularized_d3q27
      case ('d3q19')
        compute => bgk_RecursiveRegularized_d3q19
      case ('d2q9')
        compute => bgk_RecursiveRegularized_d2q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case ('prr_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_ProjectedRecursiveRegularized_d3q27
      case ('d3q19')
        compute => bgk_ProjectedRecursiveRegularized_d3q19
      case ('d2q9')
        compute => bgk_ProjectedRecursiveRegularized_d2q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case ('r_bgk')
      select case (trim(layout))
      case ('d3q27')
        compute => bgk_Regularized_d3q27
      case ('d3q19')
        compute => bgk_Regularized_d3q19
      case ('d2q9')
        compute => bgk_Regularized_d2q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for ' &
          &            //'relaxation "'//trim(relaxation)//'"')
      end select

    case default
      call tem_abort('Unsupported relaxation "'//trim(relaxation)//'" for ' &
        &            //'kind "fluid"')
    end select

  end subroutine mus_init_advRel_fluid
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> Assigning compute kernel routine by scheme relaxation type for fluid GNS kind.
  !!
  subroutine mus_init_advRel_fluid_GNS( relaxation, layout, compute )
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
        compute => bgk_advRel_d3q19_GNS
      case('d2q9')
        compute => mus_advRel_kFluidGNS_rBGK_vStd_lD2Q9
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

  end subroutine mus_init_advRel_fluid_GNS

  ! ************************************************************************** !
  !> This routine assigns compute routine for bgk relaxation.
  !!
  !! Supported variants are:
  !!   * standard        - Optimized routines for specifc layouts.
  !!   * standard_no_opt - Semi or no optimized routines for any layouts.
  !!   * improved        - improved BGK with Galilean correction term for
  !!                       specific layouts.
  !!   * block           - routines for vector machine. Implemented only for
  !!                       D3Q19
  subroutine mus_init_advRel_fluid_bgk(variant, layout, compute)
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
        compute => mus_advRel_kFluid_rBGK_vStd_lD3Q19
      case ('d2q9')
        compute => mus_advRel_kFluid_rBGK_vStd_lD2Q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for '        &
          &            //'relaxation "bgk" for variant "'//trim(variant)//'"')
      end select

    case ('standard_no_opt')
      compute => mus_advRel_kCFD_rBGK_vStdNoOpt_l

    case ('improved')
      select case (trim(layout))
      case ('d3q27')
        compute => mus_advRel_kFluid_rBGK_vImproved_lD3Q27
      case ('d2q9')
        compute => mus_advRel_kFluid_rBGK_vImproved_lD2Q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for '        &
          &            //'relaxation "bgk" for variant "'//trim(variant)//'"')
      end select

    case ('block')
      select case (trim(layout))
      case ('d3q19')
        compute => mus_advRel_kFluid_rBGK_vBlock_lD3Q19
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for '        &
          &            //'relaxation "bgk" for variant "'//trim(variant)//'"')
      end select

    case default
      call tem_abort('Unsupported variant "'//trim(variant)//'" for '        &
        &            //'relaxation "bgk"')
    end select

  end subroutine mus_init_advRel_fluid_bgk
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
  subroutine mus_init_advRel_fluid_mrt(variant, layout, compute)
    ! --------------------------------------------------------------------------
    character(len=labelLen), intent(in) :: variant
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! --------------------------------------------------------------------------
    select case (trim(variant))
    case ('standard')
      select case (trim(layout))
      case ('d3q27')
        compute => mus_advRel_kFluid_rMRT_vStd_lD3Q27
      case ('d3q19')
        compute => mus_advRel_kFluid_rMRT_vStd_lD3Q19
      case ('d2q9')
        compute => mus_advRel_kFluid_rMRT_vStd_lD2Q9
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for '        &
          &            //'relaxation "mrt" for variant "'//trim(variant)//'"')
      end select

    case ('standard_no_opt')
      select case (trim(layout))
      case ('d3q27')
        compute => mus_advRel_kCFD_rMRT_vStdNoOpt_lD3Q27
      case ('d3q19')
        compute => mus_advRel_kFluid_rMRT_vStdNoOpt_lD3Q19
      case default
        compute => mus_advRel_kCFD_rMRT_vStdNoOpt_l
      end select


    case ('bgk')
      compute => mus_advRel_kCFD_rMRT_vStdNoOpt_l

    case default
      call tem_abort('Unsupported variant "'//trim(variant)//'" for ' &
        &            //'relaxation "mrt"')
    end select

  end subroutine mus_init_advRel_fluid_mrt
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This routine assigns compute routine for trt relaxation
  !!
  !! Supported variants are:
  !!   * standard        - Optimized routines for specifc layouts.
  subroutine mus_init_advRel_fluid_trt(variant, layout, compute)
    ! --------------------------------------------------------------------------
    character(len=labelLen), intent(in) :: variant
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! --------------------------------------------------------------------------
    select case (trim(variant))
    case ('standard')
      select case (trim(layout))
      case ('d3q19')
        compute => mus_advRel_kFluid_rTRT_vStd_lD3Q19
      case ('d3q27')
        compute => mus_advRel_kFluid_rTRT_vStd_lD3Q27
      case default
        call tem_abort('Unsupported layout "'//trim(layout)//'" for '        &
          &            //'relaxation "trt" for variant "'//trim(variant)//'"')
      end select
    case default
      call tem_abort('Unsupported variant "'//trim(variant)//'" for '        &
        &            //'relaxation "trt"')
    end select

  end subroutine mus_init_advRel_fluid_trt
  ! ************************************************************************** !

end module mus_initFluid_module
! **************************************************************************** !
