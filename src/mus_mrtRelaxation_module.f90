! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> This module contains functions for MRT relaxation paramater for different
!! stencil layouts.
!! NOTE: The order of relaxation entries is consistent with moments
!! transformation matrix used in compute kernel.
!! author: Kannan Masilamani
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
module mus_mrtRelaxation_module
  use env_module,       only: rk
  use tem_aux_module,   only: tem_abort
  use tem_float_module, only: operator(.fne.)

  use mus_moments_type_module,  only: mus_moment_type
  use mus_scheme_header_module,      only: mus_scheme_header_type

  implicit none
  private

  public :: mus_assign_mrt_ptr
  public :: mus_proc_mrt

  abstract interface
    !> function pointers to obtain relaxation matrix for MRT collision operator
    pure function mus_proc_mrt(omegaKine, omegaBulk, QQ) result (s_mrt)
      import :: rk

      !> omega related to kinematic viscosity
      real(kind=rk), intent(in) :: omegaKine

      !> omega related to bulk viscosity for compressible model
      real(kind=rk), intent(in) :: omegaBulk

      !> number of directions
      integer, intent(in) :: QQ

      !> output: MRT diagonal relaxation parameter
      real(kind=rk) :: s_mrt(QQ)
    end function mus_proc_mrt

  end interface

contains

  ! ************************************************************************** !
  !> This function returns mrt function pointer according to scheme definition.
  !! In Jonas Toelke paper (2006) about MRT, the following notions are used:\n
  !!  s(a) = s(2)
  !!  s(b) = s(3)
  !!  s(c) = s(5) = s(7) = s(9)
  !!  s(d) = s(11) = s(13
  !!  s(e) = s(17) = s(18) = s(19)
  !!  s(w) = s(10) = s(12) = s(14) = s(15) = s(16)
  !! It is suggested that, for D3Q19,
  !!  s(a) = s(b) = s(c) = s(d) = s(e) = max( s(w), -1.0 )
  !! Notice that the collision matrix S used in this papar corresponds to
  !! -omega in BGK model, because it express the LB equation is slightly
  !! different way.
  !!
  subroutine mus_assign_mrt_ptr(mrtPtr, schemeHeader)
    ! --------------------------------------------------------------------------
    !> mrt function pointer
    procedure(mus_proc_mrt), pointer :: mrtPtr
    !> Scheme header information
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------
    mrtPtr => null()
    select case (trim(schemeHeader%relaxation))
    case ( 'mrt' )
      select case (trim(schemeHeader%relaxHeader%variant))
      case ('standard', 'standard_no_opt')
        select case (trim(schemeHeader%layout))
        case ('d2q9')
          if (trim(schemeHeader%kind) == 'fluid_incompressible') then
            mrtPtr => mrt_d2q9_incomp
          else
            mrtPtr => mrt_d2q9
          end if
        case ('d3q15')
          mrtPtr => mrt_d3q15
        case ('d3q19')
          mrtPtr => mrt_d3q19
        case ('d3q27')
          mrtPtr => mrt_d3q27
        case default
          call tem_abort('Error: Unknown layout for mrt relaxation')
        end select
      case ('bgk')
        ! set all relaxation paramter to same omega value for bgk variant
        mrtPtr => mrt_bgk
      case default
        call tem_abort('Error: Unknown variant for mrt relaxation')
      end select
    case default
      ! set all relaxation paramter to same omega value as default
      mrtPtr => mrt_bgk
    end select
  end subroutine mus_assign_mrt_ptr
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns mrt relaxation diagonal matrix for d2q9 layout
  !! Parameters are taken from:
  !! Lallemand, P., & Luo, L. (2000).
  !! Theory of the lattice boltzmann method: dispersion, dissipation,
  !! isotropy, galilean invariance, and stability. Physical Review. E,
  !! Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary
  !! Topics, 61(6 Pt A), 6546–62.
  !!
  !! Another paper for d2q9 with different set of parameters
  !! Chen, S., Peng, C., Teng, Y., Wang, L. P., & Zhang, K. (2016).
  !! Improving lattice Boltzmann simulation of moving particles in a
  !! viscosity flow using local grid refinement. Computers and Fluids,
  !! 136, 228–246.
  pure function mrt_d2q9_incomp(omegaKine, omegaBulk, QQ) result(s_mrt)
    ! --------------------------------------------------------------------------
    !> omega related to kinematic viscosity
    real(kind=rk), intent(in) :: omegaKine
    !> omega related to bulk viscosity
    real(kind=rk), intent(in) :: omegaBulk
    !> number of directions
    integer, intent(in) :: QQ
    !> output mrt diagonal matrix
    real(kind=rk) :: s_mrt(QQ)
    ! --------------------------------------------------------------------------
    s_mrt = 0.0_rk
    ! set relaxation by bulk viscostiy
    ! set to 1.63 for both compressible and incompressible model
    ! KM: Tried omegaBulk for compressible model and it is unstable
    s_mrt(2)   = 1.63_rk
    s_mrt(3)   = 1.14_rk
    s_mrt(5)   = 1.92_rk
    s_mrt(7)   = 1.92_rk
    s_mrt(8:9) = omegaKine

  end function mrt_d2q9_incomp
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns mrt relaxation diagonal matrix for d2q9 layout
  !! Parameters are taken from:
  !! Lallemand, P., & Luo, L. (2000).
  !! Theory of the lattice boltzmann method: dispersion, dissipation,
  !! isotropy, galilean invariance, and stability. Physical Review. E,
  !! Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary
  !! Topics, 61(6 Pt A), 6546–62.
  !!
  !! Another paper for d2q9 with different set of parameters
  !! Chen, S., Peng, C., Teng, Y., Wang, L. P., & Zhang, K. (2016).
  !! Improving lattice Boltzmann simulation of moving particles in a
  !! viscosity flow using local grid refinement. Computers and Fluids,
  !! 136, 228–246.
  !! modified by GGSpinelli. All V^n moments are relaxed with omegaBulk
  pure function mrt_d2q9(omegaKine, omegaBulk, QQ) result(s_mrt)
    ! --------------------------------------------------------------------------
    !> omega related to kinematic viscosity
    real(kind=rk), intent(in) :: omegaKine
    !> omega related to bulk viscosity
    real(kind=rk), intent(in) :: omegaBulk
    !> number of directions
    integer, intent(in) :: QQ
    !> output mrt diagonal matrix
    real(kind=rk) :: s_mrt(QQ)
    ! --------------------------------------------------------------------------
    s_mrt = 0.0_rk
    ! set relaxation by bulk viscostiy
    ! set to 1.63 for both compressible and incompressible model
    ! KM: Tried omegaBulk for compressible model and it is unstable
    ! GGS: omega2 = omegaBulk & omega3 = 1.14 => not stable!
    ! GGS: omega2 = omega3 = omegaBulk => works!
    s_mrt(2:3) = omegaBulk
    s_mrt(5)   = 1.92_rk
    s_mrt(7)   = 1.92_rk
    s_mrt(8:9) = omegaKine

  end function mrt_d2q9
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns mrt relaxation diagonal matrix for d3q15 layout
  !! Parameters are taken from:
  !! D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S.
  !! (2002). Multiple-relaxation-time lattice Boltzmann models in three
  !! dimensions. Philosophical Transactions. Series A, Mathematical,
  !! Physical, and Engineering Sciences, 360(1792), 437–51.
  pure function mrt_d3q15(omegaKine, omegaBulk, QQ) result(s_mrt)
    ! --------------------------------------------------------------------------
    !> omega related to kinematic viscosity
    real(kind=rk), intent(in) :: omegaKine
    !> omega related to bulk viscosity
    real(kind=rk), intent(in) :: omegaBulk
    !> number of directions
    integer, intent(in) :: QQ
    !> output mrt diagonal matrix
    real(kind=rk) :: s_mrt(QQ)
    ! --------------------------------------------------------------------------
    s_mrt = 0.0_rk
    s_mrt(2)     = omegaBulk
    s_mrt(3)     = 1.20_rk
    s_mrt(5)     = 1.60_rk
    s_mrt(7)     = 1.60_rk
    s_mrt(9)     = 1.60_rk
    s_mrt(10:14) = omegaKine
    s_mrt(15)    = 1.20_rk

  end function mrt_d3q15
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns mrt relaxation diagonal matrix for d3q19 layout
  !! Parameters are taken from:
  !! D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S.
  !! (2002). Multiple-relaxation-time lattice Boltzmann models in three
  !! dimensions. Philosophical Transactions. Series A, Mathematical,
  !! Physical, and Engineering Sciences, 360(1792), 437–51.
  pure function mrt_d3q19(omegaKine, omegaBulk, QQ) result(s_mrt)
    ! --------------------------------------------------------------------------
    !> omega related to kinematic viscosity
    real(kind=rk), intent(in) :: omegaKine
    !> omega related to bulk viscosity
    real(kind=rk), intent(in) :: omegaBulk
    !> number of directions
    integer, intent(in) :: QQ
    !> output mrt diagonal matrix
    real(kind=rk) :: s_mrt(QQ)
    ! --------------------------------------------------------------------------
    s_mrt = 0.0_rk
    s_mrt(2)     = omegaBulk
    s_mrt(3)     = 1.40_rk
    s_mrt(5)     = 1.20_rk
    s_mrt(7)     = 1.20_rk
    s_mrt(9)     = 1.20_rk
    s_mrt(10)    = omegaKine
    s_mrt(11)    = 1.40_rk
    s_mrt(12)    = omegaKine
    s_mrt(13)    = 1.40_rk
    s_mrt(14:16) = omegaKine
    s_mrt(17:19) = 1.98_rk

  end function mrt_d3q19
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> This function returns mrt relaxation diagonal matrix for d3q27 layout
  pure function mrt_d3q27(omegaKine, omegaBulk, QQ) result(s_mrt)
    ! --------------------------------------------------------------------------
    !> omega related to kinematic viscosity
    real(kind=rk), intent(in) :: omegaKine
    !> omega related to bulk viscosity
    real(kind=rk), intent(in) :: omegaBulk
    !> number of directions
    integer, intent(in) :: QQ
    !> output mrt diagonal matrix
    real(kind=rk) :: s_mrt(QQ)
    ! --------------------------------------------------------------------------
    s_mrt = 0.0_rk
    s_mrt(5:9)   = omegaKine
    s_mrt(10)    = omegaBulk
    s_mrt(11:13) = 1.50_rk
    s_mrt(14:17) = 1.74_rk
    s_mrt(18)    = 1.4_rk
    s_mrt(19:23) = 1.98_rk
    s_mrt(24:26) = 1.83_rk
    s_mrt(27)    = 1.61_rk

  end function mrt_d3q27
  ! ************************************************************************** !

  ! ************************************************************************** !
  !> set all relaxation parameter to same omega, results in bgk collision
  pure function mrt_bgk(omegaKine, omegaBulk, QQ) result(s_mrt)
    ! --------------------------------------------------------------------------
    !> omega related to kinematic viscosity
    real(kind=rk), intent(in) :: omegaKine
    !> omega related to bulk viscosity
    real(kind=rk), intent(in) :: omegaBulk
    !> number of directions
    integer, intent(in) :: QQ
    !> output mrt diagonal matrix
    real(kind=rk) :: s_mrt(QQ)
    ! --------------------------------------------------------------------------
    s_mrt(:) = omegaKine
  end function mrt_bgk
  ! ************************************************************************** !

end module mus_mrtRelaxation_module
