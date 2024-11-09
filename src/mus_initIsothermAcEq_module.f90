! Copyright (c) 2016 Philipp Otte <otte@mathcces.rwth-aachen.de>
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
!! flow field for lbm model for acoustic equations.
module mus_initIsothermAcEq_module
  use env_module,         only: labelLen
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  use mus_isotherm_acEq_module, only: bgk_advRel_isotherm_acEq_d3q19
  use mus_scheme_type_module,   only: kernel

  implicit none

  private

  public :: mus_init_advRel_isotherm_acEq

contains

  ! **************************************************************************** !
  !> Assigning compute kernel routine by scheme relaxation type for
  !> isotherm_acEq kind.
  !!
  subroutine mus_init_advRel_isotherm_acEq( relaxation, layout, compute )
    ! ---------------------------------------------------------------------------
    character(len=labelLen), intent(inout) :: relaxation
    character(len=labelLen), intent(in) :: layout
    procedure( kernel ), pointer, intent(out) :: compute
    ! ---------------------------------------------------------------------------

    write(logUnit(1),*) 'Choosing relaxation model: '//trim( relaxation )
    select case( trim( relaxation ))
    case( 'bgk' )
      select case( trim(layout) )
      case( 'd3q19' )
        compute => bgk_advRel_isotherm_acEq_d3q19
      case default
        write(logUnit(1),*) 'Stencil '//trim(layout)//' is not supported yet!'
        call tem_abort()
      end select
    case default
      write(logUnit(1),*) 'Relaxation '//trim(relaxation)//' is not supported!'
      call tem_abort()
    end select

  end subroutine mus_init_advRel_isotherm_acEq
  ! **************************************************************************** !

end module mus_initIsothermAcEq_module
! ****************************************************************************** !
