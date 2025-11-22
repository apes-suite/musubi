! Copyright (c) 2013 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2025 Mengyu Wang <m.wang-2@utwente.nl>
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
!> This module contains routines which assigns compute kernel for passive scalar
!! model
module mus_initLBMPS_module

  ! include treelm modules
  use env_module,         only: labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  ! include musubi modules
  use mus_compute_passiveScalar_module, only: mus_advRel_kPS_rBGK_v1st_l,      &
    &                                         mus_advRel_kPS_rBGK_v2nd_l,      &
    &                                         mus_advRel_kPS_rTRT_vStdNoOpt_l
  use mus_d3q19_module,                 only: mus_advRel_kPS_rBGK_vEmodel_lD3Q19, &
    &                                         mus_advRel_kPS_rBGK_vEmodelCorr_lD3Q19, &
    &                                         mus_advRel_kPS_rMRT_vEmodelCorr_lD3Q19, &
    &                                         mus_advRel_kPS_rTRT_vEmodel_lD3Q19, &
    &                                         mus_advRel_kPS_rTRT_vEmodelCorr_lD3Q19, &
    &                                         mus_advRel_kPS_rTRT_vLmodel_lD3Q19

  use mus_scheme_type_module,           only: kernel

  implicit none

  private

  public :: mus_init_advRel_lbm_ps

contains

! ****************************************************************************** !
  !> Initialize the relaxation model for lbm passive scalar scheme kind
  subroutine mus_init_advRel_lbm_ps( relaxation, layout, relaxation_variant, compute )
    ! ---------------------------------------------------------------------------
    character(len=labelLen), intent(in) :: relaxation
    character(len=labelLen), intent(in) :: layout
    character(len=labelLen), intent(in) :: relaxation_variant
    procedure( kernel ), pointer, intent(out) :: compute
    ! ---------------------------------------------------------------------------

    write(logUnit(1),*) 'Choosing LBM Passive Scalar relaxation model: '//     &
      &                 trim(relaxation)

    select case( trim(relaxation) )
    case( 'bgk' )
      select case( trim(relaxation_variant) )
      case( 'first' )
        write(logUnit(1), *) 'Using bgk_advRel scheme with first order.'
        compute => mus_advRel_kPS_rBGK_v1st_l
      case( 'second' )
        write(logUnit(1), *) 'Using bgk_advRel scheme with second order.'
        compute => mus_advRel_kPS_rBGK_v2nd_l
      case( 'Emodel' )
        write(logUnit(1), *) 'Using bgk_advRel scheme with E-model.'
        select case( trim(layout) )
        case( 'd3q19' )
          compute => mus_advRel_kPS_rBGK_vEmodel_lD3Q19
        case default
          write(logUnit(1),*) 'The selected layout is not supported: ' &
            &                  // trim(layout)
          call tem_abort()
        end select
      case( 'EmodelCorr' )
        write(logUnit(1), *) 'Using bgk_advRel scheme with E-model correction.'
        select case( trim(layout) )
        case( 'd3q19' )
          compute => mus_advRel_kPS_rBGK_vEmodelCorr_lD3Q19
        case default
          write(logUnit(1),*) 'The selected layout is not supported: ' &
            &                  // trim(layout)
          call tem_abort()
        end select
      case default
        write(logUnit(1),*) 'relaxation_variant '//trim(relaxation_variant)//   &
          &                 ' is not supported yet!'
        call tem_abort()
      end select
    case( 'trt' )
      select case( trim(relaxation_variant) )
      case( 'Emodel' )
        write(logUnit(1), *) 'Using trt_advRel scheme with E-model.'
        select case( trim(layout) )
        case( 'd3q19' )
          compute => mus_advRel_kPS_rTRT_vEmodel_lD3Q19
        case default
          write(logUnit(1),*) 'The selected layout is not supported: ' &
            &                  // trim(layout)
          call tem_abort()
        end select
      case( 'EmodelCorr' )
        write(logUnit(1), *) 'Using trt_advRel scheme with E-model correction.'
        select case( trim(layout) )
        case( 'd3q19' )
          compute => mus_advRel_kPS_rTRT_vEmodelCorr_lD3Q19
        case default
          write(logUnit(1),*) 'The selected layout is not supported: ' &
            &                  // trim(layout)
          call tem_abort()
        end select
      case( 'Lmodel' )
        write(logUnit(1), *) 'Using trt_advRel scheme with L-model.'
        select case( trim(layout) )
        case( 'd3q19' )
          compute => mus_advRel_kPS_rTRT_vLmodel_lD3Q19
        case default
          write(logUnit(1),*) 'The selected layout is not supported: ' &
            &                  // trim(layout)
          call tem_abort()
        end select
      case default
        write(logUnit(1),*) 'Using trt_advRel scheme with standard no optimization.'
        compute => mus_advRel_kPS_rTRT_vStdNoOpt_l
      end select
    case( 'mrt' )
      select case( trim(layout) )
      case( 'd3q19' )
        write(logUnit(1), *) 'Using mrt_advRel scheme with E-model correction.'
        compute => mus_advRel_kPS_rMRT_vEmodelCorr_lD3Q19
      case default
        write(logUnit(1),*) 'The selected layout is not supported: ' &
          &                  // trim(layout)
        call tem_abort()
      end select
    case default
      write(logUnit(1),*) 'The selected relaxation model is not supported: '// &
        &                  trim(relaxation)
      call tem_abort()
    end select

  end subroutine mus_init_advRel_lbm_ps
! ****************************************************************************** !

end module mus_initLBMPS_module
