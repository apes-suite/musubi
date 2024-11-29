! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2011 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2015, 2018-2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012, 2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2018 Jana Gericke <jana.gericke@uni-siegen.de>
! Copyright (c) 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2021-2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> author: Jiaxing Qi
!! This module provides functions for calculating macroscopic quantities
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
module mus_derivedQuantities_module2

  ! include treelm modules
  use env_module,               only: rk
  use tem_stencil_module,       only: tem_stencilHeader_type
  use tem_param_module,         only: div1_3, div3_4, cs2inv, div1_9, div1_54, &
    &                                 cs2, cs4inv, div1_2, sqrt3, rho0
  use tem_logging_module,       only: logUnit
  use tem_float_module,         only: operator(.fne.)
  use tem_aux_module,           only: tem_abort
  !use tem_property_module,      only: prp_fluid
  !use tem_varSys_module,        only: tem_varSys_type
  !use tem_construction_module,  only: tem_levelDesc_type

  ! include musubi modules
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_moments_module,       only: get_moment
  !use mus_varSys_module,        only: mus_varSys_data_type

  implicit none
  private

  public :: getDensity, getVelocity, getVelocity_incomp
  public :: getEquilibriumIncomp, getEquilibrium
  public :: getNEq_diffusive
  public :: getNEq_acoustic
  public :: convPrePost
  public :: secondMom_2D
  public :: secondMom_3D
  public :: secondMom_minus_cs2_2D
  public :: secondMom_minus_cs2_3D
  public :: getShearRate
  public :: getNonEqFac_intp
  public :: geteqbydensvel
  public :: getNonEqFac_intp_coarse_to_fine
  public :: getNonEqFac_intp_fine_to_coarse

  interface getEquilibrium
    module procedure getEquilibrium_forElemfromState
    module procedure getEqByDensVel
    module procedure getEquilibrium_forPdfSubset
  end interface

  interface getDensity
    module procedure getDensity_forElemfromState
    module procedure getDensity_forPdfSubset
  end interface

  interface getVelocity
    module procedure getVelocity_forElemFromState_noForce
    module procedure getVelocity_forPdfSubset
  end interface

  interface getVelocity_incomp
    module procedure getVelocity_forPdfSubset_incomp
  end interface getVelocity_incomp

contains

! **************************************************************************** !
  !> Calculate the density of a given subset of pdfs
  !!        vector (sum up all values)
  !!
  pure function getDensity_forPdfSubset( subset, stencil, varPos ) result( res )
    ! --------------------------------------------------------------------------
    type(tem_stencilHeader_type), intent(in) :: stencil
    real(kind=rk), intent(in) :: subset(:)
    integer, intent(in) :: varPos(:) !< varPos of current field variable
    real(kind=rk)             :: res !< return value
    ! --------------------------------------------------------------------------
    ! local variables
    integer :: iDir
    ! --------------------------------------------------------------------------

    res = 0._rk
    do iDir = 1,stencil%QQ
      res = res + subset( varPos( iDir ) )
    enddo

  end function getDensity_forPdfSubset
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the density of a given element number with the given state
  !!        vector (sum up all values)
  !!
  pure function getDensity_forElemFromState( state, elem, stencil, varPos,     &
    &                                        nScalars ) result( res )
    ! --------------------------------------------------------------------------
    type(tem_stencilHeader_type), intent(in) :: stencil
    real(kind=rk), intent(in) :: state(:)
    integer, intent(in)       :: elem
    integer, intent(in) :: varPos(:) !< varPos of current field variable
    integer, intent(in) :: nScalars !< number of scalars in global system
    real(kind=rk)             :: res !< return value
    ! --------------------------------------------------------------------------
    ! local variables
    integer :: iDir
    integer :: nElems
    ! --------------------------------------------------------------------------
    nElems = size( state ) / nScalars

    res = 0._rk
    do iDir = 1, stencil%QQ
      res = res + state( ( elem-1)* nscalars+ varpos(idir))
    enddo

  end function getDensity_forElemFromState
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the velocity in all 3 directions
  !!        from a subset given, ordered according to the stencil
  !!
  pure function getVelocity_forPdfSubset( subset, stencil, varPos ) &
    &                             result( vel )
    ! --------------------------------------------------------------------------
    type(tem_stencilHeader_type), intent(in) :: stencil !< stencil information
    real(kind=rk), intent(in) :: subset(:) !< complete state of one level
    integer,       intent(in) :: varPos(:) !< varPos of current field variable
    real(kind=rk)             :: vel(3)    !< return value
    ! --------------------------------------------------------------------------
    real(kind=rk) :: dens
    integer       :: iDir
    ! --------------------------------------------------------------------------

    vel = 0._rk
    dens = 0._rk
    do iDir = 1,stencil%QQ
      vel(:) = vel(:) + subset(varPos(iDir)) * stencil%cxDirRK(:,iDir)
      dens   = dens   + subset(varPos(iDir))
    enddo
    vel = vel / dens

  end function getVelocity_forPdfSubset
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the velocity in all 3 directions
  !!        from a subset given, ordered according to the stencil
  !!
  pure function getVelocity_forPdfSubset_incomp( subset, stencil, varPos )     &
    &                             result( vel )
    ! --------------------------------------------------------------------------
    type(tem_stencilHeader_type), intent(in) :: stencil !< stencil information
    real(kind=rk), intent(in) :: subset(:) !< complete state of one level
    integer,       intent(in) :: varPos(:) !< varPos of current field variable
    real(kind=rk)             :: vel(3)    !< return value
    ! --------------------------------------------------------------------------
    integer       :: iDir
    ! --------------------------------------------------------------------------

    vel = 0._rk
    do iDir = 1,stencil%QQ
      vel(:) = vel(:) + subset(varPos(iDir)) * stencil%cxDirRK(:,iDir)
    enddo
    vel = vel / rho0

  end function getVelocity_forPdfSubset_incomp
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the velocity in all 3 directions
  !! from the element indicated (elem),
  !! reading the pdf (state information) from the state array.
  !! state array includes all the pdfs of all elements.
  !! The access to the state array has to be done via the generic
  !! access macro IDX, as we want to access post-collision values.
  !!
  pure function getVelocity_forElemFromState_noForce( state, elem, stencil,    &
    &                                           varPos, nScalars ) result( vel )
    ! --------------------------------------------------------------------------
    type(tem_stencilHeader_type), intent(in) :: stencil !< stencil information
    real(kind=rk), intent(in) :: state(:)  !< complete state of one level
    !> element index, for which to calc velocity
    integer, intent(in) :: elem
    integer, intent(in) :: varPos(:) !< varPos of current field variable
    integer, intent(in) :: nScalars !< number of scalars in global system
    real(kind=rk)       :: vel(3)    !< return value
    ! --------------------------------------------------------------------------
    real(kind=rk) :: dens
    integer       :: iDir
    integer :: nElems
    ! --------------------------------------------------------------------------
    nElems = size( state ) / nScalars

    vel = 0._rk
    dens = 0._rk
    do iDir = 1,stencil%QQ
      vel(:) = vel(:) + state( ( elem-1)* nscalars+ varpos(idir))&
        &               * stencil%cxDirRK(:,iDir)
      dens   = dens   + state( ( elem-1)* nscalars+ varpos(idir))
    enddo

    vel = vel/dens

  end function getVelocity_forElemFromState_noForce
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the equilibrium distribution function in all directions
  !!
  !! The equilibrium distribution function is:\n
  !! \[ f^{eq}_i = w_i \rho ( 1 + \frac{\vec c_i \cdot \vec u}{c^2_s}
  !!                      + \frac{ {(\vec c_i \cdot \vec u)}^2}{2c^4_s}
  !!                      - \frac{\vec u \cdot \vec u}{2c^2_s}) \]
  !!
  !! where \(w_i\) is the weight in each direction,\n
  !! \(\rho\) is the macroscopic value of density,\n
  !! \(c_s\) is the speed of sound,\n
  !! \(\vec c_i\) is the lattice unit velocity in each direction,\n
  !! \(\vec u\) is the macroscopic value of velocity.
  !!
  pure function getEqByDensVel( dens, vel, layout ) result( equil )
    ! --------------------------------------------------------------------------
    type(mus_scheme_layout_type), intent(in) :: layout !scheme layout
    real(kind=rk), intent(in) :: dens
    real(kind=rk), intent(in) :: vel(3)
    real(kind=rk)             :: equil(layout%fStencil%QQ) !< return value
    ! --------------------------------------------------------------------------
    real(kind=rk) :: ucx, usq
    integer :: iDir
    ! --------------------------------------------------------------------------

    ! square of velocity
    usq = vel(1)*vel(1) + vel(2)*vel(2) + vel(3)*vel(3)

    do iDir = 1, layout%fStencil%QQ

      ! velocity times lattice unit velocity
      ucx =   layout%fStencil%cxDirRK(1, iDir) * vel(1) &
        &   + layout%fStencil%cxDirRK(2, iDir) * vel(2) &
        &   + layout%fStencil%cxDirRK(3, iDir) * vel(3)

      ! calculate equilibrium density
      equil(iDir) =   layout%weight(iDir) * dens * ( 1._rk + ucx*cs2inv &
        &           + ucx*ucx*cs2inv*cs2inv*div1_2                     &
        &           - usq*cs2inv*div1_2 )

    enddo

  end function getEqByDensVel
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the equilibrium distribution function in all directions
  !!
  !! The equilibrim distribution function is:\n
  !! \[ f^{eq}_i = w_i \rho ( 1 + \frac{\vec c_i \cdot \vec u}{c^2_s}
  !!                      + \frac{ {(\vec c_i \cdot \vec u)}^2}{2c^4_s}
  !!                      - \frac{\vec u \cdot \vec u}{2c^2_s}) \]\n
  !! where \(w_i\) is the weight in each direction,\n
  !! \(\rho\) is the macroscopic value of density,\n
  !! \(c_s\) is the speed of sound,\n
  !! \(\vec c_i\) is the lattice unit velocity in each direction,\n
  !! \(\vec u\) is the macroscopic value of velocity.
  !!
  pure function getEquilibrium_forPdfSubset( subset, layout, varPos )          &
    &                                                            result( equil )
    ! --------------------------------------------------------------------------
    type(mus_scheme_layout_type), intent(in) :: layout !scheme layout
    real(kind=rk), intent(in) :: subset(:)  !< pdf array
    integer, intent(in) :: varPos(:) !< varPos of current field variable
    real(kind=rk) :: equil(layout%fStencil%QQ) !< return value
    ! --------------------------------------------------------------------------
    real(kind=rk) :: rho, vel(3)
    real(kind=rk) :: ucx, usq
    integer :: iDir
    ! --------------------------------------------------------------------------

    rho = getDensity_forPdfSubset(  subset, layout%fStencil, varPos )
    vel = getVelocity_forPdfSubset( subset, layout%fStencil, varPos )
    ! square of velocity
    usq = vel(1)*vel(1) + vel(2)*vel(2) + vel(3)*vel(3)

    do iDir = 1, layout%fStencil%QQ

      ! velocity times lattice unit velocity
      ucx =   layout%fStencil%cxDirRK(1, iDir)*vel(1) &
        &   + layout%fStencil%cxDirRK(2, iDir)*vel(2) &
        &   + layout%fStencil%cxDirRK(3, iDir)*vel(3)

      ! calculate equilibrium density
      equil( iDir ) = layout%weight( iDir ) * rho * ( 1._rk + ucx*cs2inv &
        &           + ucx*ucx*cs2inv*cs2inv*0.5_rk                       &
        &           - usq*cs2inv*0.5_rk )

    enddo

  end function getEquilibrium_forPdfSubset
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the equilibrium distribution function in all directions
  !!
  !! The equilibrim distribution function is:\n
  !! \[ f^{eq}_i = w_i \rho ( 1 + \frac{\vec c_i \cdot \vec u}{c^2_s}
  !!                      + \frac{ {(\vec c_i \cdot \vec u)}^2}{2c^4_s}
  !!                      - \frac{\vec u \cdot \vec u}{2c^2_s}) \]\n
  !! where \(w_i\) is the weight in each direction,\n
  !! \(\rho\) is the macroscopic value of density,\n
  !! \(c_s\) is the speed of sound,\n
  !! \(\vec c_i\) is the lattice unit velocity in each direction,\n
  !! \(\vec u\) is the macroscopic value of velocity.
  !!
  pure function getEquilibrium_forElemfromState( state, elem, layout,          &
    &                                  varPos, nScalars, neigh ) result( equil )
    ! --------------------------------------------------------------------------
    type(mus_scheme_layout_type), intent(in) :: layout !scheme layout
    real(kind=rk), intent(in) :: state(:)   !< pdf array
    integer, intent(in)       :: elem       !< treeID of the target element
    integer, intent(in) :: varPos(:) !< varPos of current field variable
    integer, intent(in) :: nScalars !< number of scalars in global system
    integer, intent(in) :: neigh(:)   !< connectivity vector
    real(kind=rk) :: equil(layout%fStencil%QQ) !< return value
    ! --------------------------------------------------------------------------
    real(kind=rk) :: rho, vel(3)
    real(kind=rk) :: ucx, usq
    integer :: iDir
    ! --------------------------------------------------------------------------
    rho = getDensity( state, elem, layout%fStencil, varPos, nScalars )
    vel(:) = getVelocity( state, elem, layout%fStencil, varPos, nScalars )

    ! square of velocity
    usq = vel(1)*vel(1) + vel(2)*vel(2) + vel(3)*vel(3)

    do iDir = 1, layout%fStencil%QQ

      ! velocity times lattice unit velocity
      ucx =   layout%fStencil%cxDirRK( 1, iDir )*vel(1)          &
        &   + layout%fStencil%cxDirRK( 2, iDir )*vel(2)          &
        &   + layout%fStencil%cxDirRK( 3, iDir )*vel(3)

      ! calculate equilibrium density
      equil( iDir ) = layout%weight( iDir ) * rho * ( 1._rk + ucx*cs2inv        &
        &           + ucx*ucx*cs2inv*cs2inv*0.5_rk                              &
        &           - usq*cs2inv*0.5_rk )

    enddo

  end function getEquilibrium_forElemfromState
! **************************************************************************** !


! **************************************************************************** !
  !> author: Jiaxing Qi
  !! Calculate the Shear Rate
  !!
  !! The Shear Rate is defined as
  !! \[
  !!  \dot{\gamma} = 2\sqrt{ D_{II} }
  !! \]
  !! where \( D_{II} \) is the second invariant of the strain rate tensor and
  !! defined as
  !! \[
  !!    D_{II} = \sum^{l}_{\alpha,\beta=l} S_{\alpha\beta} S_{\alpha\beta}
  !! \]
  !! where \( S_{\alpha\beta} \) is the strain rate tensor.
  !!
  pure function getShearRate( strain ) result( shear )
    ! --------------------------------------------------------------------------
    !> strain rate tensor: xx, yy, zz, xy, yz, zx
    real(kind=rk), intent(in) :: strain(:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: shear !< shear rate
    ! --------------------------------------------------------------------------

    shear = sqrt( sum( strain(:) * strain(:) ) ) * 2._rk

  end function getShearRate
! **************************************************************************** !


! **************************************************************************** !
  !> author: Jiaxing Qi
  !! Setting the non-equilibrium part based on the acoustic scaling
  !!
  !! The non-equilibirium part of pdf is computed from strain rate tensor
  !! by \cite Latt:2011vr
  !! \[
  !!  f_i^{neq} = - \frac{t_i \rho_0}{c_{s}^2 \omega)}
  !!                            Q_{i\alpha\beta}:S_{\alpha\beta}
  !! \]
  !! where \( \boldsymbol{A} : \boldsymbol{B} \) is Frobenius inner product.
  !! \[
  !!   \boldsymbol{A} : \boldsymbol{B} = \sum_{i,j} A_{ij}B_{ij}
  !! \]
  !! \[
  !!  Q_{i\alpha\beta} = c_{i\alpha}c_{i\beta} - c_s^2 \delta_{\alpha\beta}
  !! \]
  !! and \( S \) is the strain rate tensor
  !! \[
  !!  S_{\alpha\beta} = -\frac{1}{2}
  !!          ( \partial_{\alpha}u_{\beta} + \partial_{\alpha}u_{\beta} )
  !! \]
  !!
  function getNEq_acoustic( layout, omega, Sxx ) result( nEq )
    ! --------------------------------------------------------------------------
    type( mus_scheme_layout_type ), intent(in) :: layout
    real(kind=rk), intent(in) :: omega
    !> strain rate tensor
    real(kind=rk), intent(in) :: Sxx(3,3)
    real(kind=rk) :: nEq(layout%fStencil%QQ)
    ! --------------------------------------------------------------------------
    integer :: iVal, jVal, iDir
    real(kind=rk) :: tau(3,3)
    real(kind=rk) :: nu, coeff
    ! --------------------------------------------------------------------------

    nu = cs2 * ( 1._rk / omega - div1_2 )
    ! convert strain rate to stress
    tau(:,:) = 2._rk * nu * Sxx(:,:)
    coeff = cs4inv / ( 2._rk - omega )

    ! Recover the non-equilibrium part from stress (tau)
    nEq(:) = 0._rk
    do iDir = 1, layout%fStencil%QQ
      do jVal = 1, layout%fStencil%nDims
        do iVal = 1, layout%fStencil%nDims
          Neq( iDir ) = Neq( iDir ) + &
          &   tau(iVal,jVal) * layout%fStencil%cxDirRK(iVal,iDir) &
          &                   *layout%fStencil%cxDirRK(jVal,iDir)
        end do ! iVal
        neq( iDir ) = neq( iDir ) - cs2 * tau( jVal, jVal )
      end do ! jVal
    end do ! iDir
    neq(:) = -layout%weight(:) * neq(:) * coeff

    ! Convert from post to pre-collision
    ! KM: convert factor is zero for omega = 1.0 and dividing with 0.0
    ! leads to NaN so divide conv factor only if omega is not equal to 1
    if (omega .fne. 1.0_rk) Neq = Neq / convPrePost( omega )

  end function getNEq_acoustic
! **************************************************************************** !


! **************************************************************************** !
  !> author: Jiaxing Qi
  !! Calculate the non-equilibrium part of pdf from strain rate tensor
  !! based on the diffusive scaling
  !!
  !! According to \cite Junk:2005cr \n
  !! The non-equilibrium part of pdf \( f_i^{neq} \) is set by \n
  !! \[
  !!   f_i^{neq} = -\frac{t_i}{2\kappa c_{s}^4} \nu' S^{(1)}:\Lambda
  !!             = -\frac{t_i}{2\kappa c_{s}^4} \frac{\kappa c_s^2}{\omega} S^{(1)}:\Lambda
  !!             = -\frac{t_i}{2 c_{s}^2} S^{(1)}:\Lambda
  !!             = -\frac{t_i}{c_{s}^2} \bm S:\Lambda
  !! \]
  !! where \( \nu' = \frac{\kappa c_s^2}{\omega} \) is the viscosity,
  !! \( \bm S = \frac{1}{2}S^{(1)} \) is the strain rate tensor and
  !! \[
  !! \Lambda_{i\alpha\beta} =
  !!    c_{i\alpha} c_{i\beta} - \frac{1}{D}\sum_{\gamma}(c_{i\gamma}c_{i\gamma})
  !!                                                      \delta_{\alpha\beta}
  !! \]
  !! and \( D\) is the number of dimension. \n
  !! Notice here that strain rate tensor above has to be a traceless tensor, i.e.
  !! \( Tr(S) = 0 \).
  !! In current implementation,
  !! the above equation is slightly modified so that
  !! the strain rate tensor is not required to be traceless anymore.
  !! In this way, \( f_i^{neq} \) calculated by this routine can recover the
  !! input strain rate tensor no matter it is traceless or not.\n
  !! Specificly the \( \Lambda \) in above equation is modified slightly, i.e.
  !! \[
  !! \Lambda_{i\alpha\beta} =
  !!                       c_{i\alpha} c_{i\beta} - c_s^2 \delta_{\alpha\beta}
  !! \]
  !! This routine has a unit test program utest/mus_fNeq_diffusive_test
  !!
  function getNEq_diffusive( layout, omega, Sxx ) result( nEq )
    ! --------------------------------------------------------------------------
    type( mus_scheme_layout_type ), intent(in) :: layout
    real(kind=rk), intent(in) :: omega
    !> Strain rate tensor. It is a symmetric 3x3 matrix
    real(kind=rk), intent(in) :: Sxx(3,3)
    real(kind=rk) :: nEq(layout%fStencil%QQ)
    ! --------------------------------------------------------------------------
    integer :: iVal, jVal, iDir
    real(kind=rk) :: strain(3,3)
    ! --------------------------------------------------------------------------

    strain = Sxx
    ! First calculate the part of f_neq = strain : Lambda
    nEq(:) = 0._rk
    do iDir = 1, layout%fStencil%QQ
      ! cx2Sum = 0
      ! do iVal = 1, layout%fStencil%nDims
      !   cx2Sum = cx2Sum + layout%fStencil%cxDir(iVal,iDir)**2
      ! end do
      ! cx2Sum = cx2Sum / nDims
      do jVal = 1, layout%fStencil%nDims
        do iVal = 1, layout%fStencil%nDims
          Neq( iDir ) = Neq( iDir ) + strain( iVal, jVal )     &
            &           *layout%fStencil%cxDirRK( iVal,iDir) &
            &           *layout%fStencil%cxDirRK( jVal,iDir)
        end do

        ! By this equation, strain rate tensor is required to be traceless
        ! Neq( iDir ) = Neq( iDir ) - strain( jVal, jVal ) * real(cx2sum, rk)

        ! By the following, strain rate tensor is NOT required to be traceless
        Neq( iDir ) = Neq( iDir ) - strain( jVal, jVal ) * cs2

      end do
    end do
    ! Then multiply the pre-factor
    ! f_neq = -t_i/cs^2/omega * (strain:Lambda)
    Neq(:) = -layout%weight(:) / omega / cs2 * Neq(:)

    ! Convert from post to pre-collision
    ! KM: convert factor is zero for omega = 1.0 and dividing with 0.0
    ! leads to NaN so divide conv factor only if omega is not equal to 1
    if (omega .fne. 1.0_rk) Neq = Neq / convPrePost( omega )

  end function getNEq_diffusive
! **************************************************************************** !


! **************************************************************************** !
  !> author: Jiaxing Qi
  !! Conversion factor betwen the pre- and post-collision quantity for the shear
  !! stress.
  !!
  !! Shear stress calculation requires the non-equilibirium value of pdf before
  !! collision. However that value not may be accessable directly when PULL
  !! scheme is uitilized, as only pdf after collision is available.
  !! So this conversion factor is introduced to help
  !! recover fNeq before collision from fNeq after collision as long as
  !! relaxation parameter (omage) does not equal to 1.0. When omage equals to 1,
  !! this conversion factor is set to be 0.
  !!
  !! How to use this pre-factor?
  !!```
  !!  shearstress = convPrePost(omega) * omega * cs2inv * shearLB_postColl
  !!```
  !!
  pure function convPrePost( omega ) result( conv )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: omega !< relaxation parameter
    real(kind=rk) :: conv !< conversion factor
    ! --------------------------------------------------------------------------

      conv = 1._rk / ( 1._rk - omega )
  end function convPrePost
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the conversion factor to convert nonEq moments
  !! between fine and coarser.
  pure function getNonEqFac_intp( omegaS, omegaT ) result( fac )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: omegaS  !< omega value on source level
    real(kind=rk), intent(in) :: omegaT  !< omage value on target level
    real(kind=rk) :: fac
    ! --------------------------------------------------------------------------
      fac =   omegaS * ( 1._rk - omegaT ) &
        &   / ( ( 1._rk - omegaS ) * omegaT )

  end function getNonEqFac_intp
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the conversion factor to convert nonEq pdfs
  !! from coarse to fine.
  pure function getNonEqFac_intp_coarse_to_fine( omegaC, omegaF ) result( fac )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: omegaC  !< omega value on coarse level
    real(kind=rk), intent(in) :: omegaF  !< omage value on target level
    real(kind=rk) :: fac
    ! --------------------------------------------------------------------------
    fac = 0.5_rk * getNonEqFac_intp(omegaS = omegaC, omegaT = omegaF)

  end function getNonEqFac_intp_coarse_to_fine
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the conversion factor to convert nonEq pdfs
  !! from fine to coarse.
  pure function getNonEqFac_intp_fine_to_coarse( omegaC, omegaF ) result( fac )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: omegaC  !< omega value on coarse level
    real(kind=rk), intent(in) :: omegaF  !< omage value on target level
    real(kind=rk) :: fac
    ! --------------------------------------------------------------------------
    fac = 2._rk * getNonEqFac_intp(omegaS = omegaF, omegaT = omegaC)

  end function getNonEqFac_intp_fine_to_coarse
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate the equilibrium distribution function in all directions
  !! This is the incompressible formulation with reference density rho0
  !!
  !! The equilibrium distribution function is:\n
  !! \[ f^{eq}_i = w_i  ( \rho + \frac{\vec c_i \cdot \vec u}{c^2_s}
  !!                      + \frac{ {(\vec c_i \cdot \vec u)}^2}{2c^4_s}
  !!                      - \frac{\vec u \cdot \vec u}{2c^2_s}) \]\n
  !! where \(w_i\) is the weight in each direction,\n
  !! \(\rho = \sum_i f_i\) is the macroscopic density,\n
  !! \(c_s\) is the speed of sound,\n
  !! \(\vec c_i\) is the lattice unit velocity in each direction,\n
  !! \(\vec u = \sum_i c_i f_i\) is the macroscopic value of velocity.
  !!
  pure function getEquilibriumIncomp( dens, vel, layout, rho0 )  result( equil )
    ! --------------------------------------------------------------------------
    type(mus_scheme_layout_type), intent(in) :: layout !scheme layout
    real(kind=rk), intent(in) :: dens
    real(kind=rk), intent(in) :: rho0
    real(kind=rk), intent(in) :: vel(3)
    real(kind=rk)             :: equil(layout%fStencil%QQ) !< return value
    ! --------------------------------------------------------------------------
    ! local variables
    real(kind=rk) :: ucx, usq
    integer :: iDir
    ! --------------------------------------------------------------------------

    ! square of velocity
    usq = vel(1)*vel(1) + vel(2)*vel(2) + vel(3)*vel(3)

    do iDir = 1, layout%fStencil%QQ

      ! velocity times lattice unit velocity
      ucx =   layout%fStencil%cxDirRK(1,iDir) * vel(1) &
        &   + layout%fStencil%cxDirRK(2,iDir) * vel(2) &
        &   + layout%fStencil%cxDirRK(3,iDir) * vel(3)

      ! calculate equilibrium density
      equil( iDir ) =   layout%weight( iDir )         &
        &             * ( dens + rho0*( ucx*cs2inv    &
        &             + ucx*ucx*cs2inv*cs2inv*div1_2  &
        &             - usq*cs2inv*div1_2 ))

    end do

  end function getEquilibriumIncomp
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate second moments of some quantity \( f \)
  !! \[
  !!    M_{\alpha\beta} = \sum_{i=1}^{Q} c_{i\alpha} c_{i\beta} f_i
  !! \]
  !! where Q is number of discrete velocity.\n
  !! The output is 1 dimentional array which has 6 componenents.\n
  !! Specifically,
  !! \[ m_1 = \sum_{i=1}^{Q} c_{i1} c_{i1} f_i \]
  !! \[ m_2 = \sum_{i=1}^{Q} c_{i2} c_{i2} f_i \]
  !! \[ m_3 = \sum_{i=1}^{Q} c_{i1} c_{i2} f_i \]
  !! This function is used by shear stress and strain rate.
  !! 1=xx, 2=yy, 3=xy
  !!
  pure function secondMom_2D( cxcx, f, QQ ) result ( m )
    ! --------------------------------------------------------------------------
    integer,       intent(in) :: QQ !< number of discrete directions (=QQ)
    real(kind=rk), intent(in) :: cxcx(6,QQ)   !< discrete velocity of stencil
    !> quantity to which second moment is applied
    real(kind=rk), intent(in) :: f(QQ)
    real(kind=rk)             :: m(3) !< output array
    ! --------------------------------------------------------------------------

    m(1) = sum( cxcx(1,:) * f(:) )
    m(2) = sum( cxcx(2,:) * f(:) )
    m(3) = sum( cxcx(3,:) * f(:) )

  end function secondMom_2D
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate second moments of some quantity \( f \)
  !! \[
  !!    M_{\alpha\beta} = \sum_{i=1}^{Q} c_{i\alpha} c_{i\beta} f_i
  !! \]
  !! where Q is number of discrete velocity.\n
  !! The output is 1 dimentional array which has 6 componenents.\n
  !! Specifically,
  !! \[ m_1 = \sum_{i=1}^{Q} c_{i1} c_{i1} f_i \]
  !! \[ m_2 = \sum_{i=1}^{Q} c_{i2} c_{i2} f_i \]
  !! \[ m_3 = \sum_{i=1}^{Q} c_{i3} c_{i3} f_i \]
  !! \[ m_4 = \sum_{i=1}^{Q} c_{i1} c_{i2} f_i \]
  !! \[ m_5 = \sum_{i=1}^{Q} c_{i2} c_{i3} f_i \]
  !! \[ m_6 = \sum_{i=1}^{Q} c_{i3} c_{i1} f_i \]
  !! This function is used by shear stress and strain rate.
  !! 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=xz
  !! in 2D:
  !! 1=xx, 2=yy, 3=xy
  !!
  pure function secondMom_3D( cxcx, f, QQ ) result ( m )
    ! --------------------------------------------------------------------------
    integer,       intent(in) :: QQ !< number of discrete directions (=QQ)
    real(kind=rk), intent(in) :: cxcx(6,QQ)   !< discrete velocity of stencil
    !> quantity to which second moment is applied
    real(kind=rk), intent(in) :: f(QQ)
    real(kind=rk)             :: m(6) !< output array
    ! --------------------------------------------------------------------------

    m(1) = sum( cxcx(1,:) * f(:) )
    m(2) = sum( cxcx(2,:) * f(:) )
    m(3) = sum( cxcx(3,:) * f(:) )
    m(4) = sum( cxcx(4,:) * f(:) )
    m(5) = sum( cxcx(5,:) * f(:) )
    m(6) = sum( cxcx(6,:) * f(:) )

  end function secondMom_3D
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate second moments of some quantity \( f \)
  !! \[
  !!    M_{\alpha\beta} = \sum_{i=1}^{Q} c_{i\alpha} c_{i\beta} f_i
  !! \]
  !! where Q is number of discrete velocity.\n
  !! The output is 1 dimentional array which has 6 componenents.\n
  !! Specifically,
  !! \[ m_1 = \sum_{i=1}^{Q} c_{i1} c_{i1} f_i \]
  !! \[ m_2 = \sum_{i=1}^{Q} c_{i2} c_{i2} f_i \]
  !! \[ m_3 = \sum_{i=1}^{Q} c_{i1} c_{i2} f_i \]
  !! This function is used by shear stress and strain rate.
  !! 1=xx, 2=yy, 3=xy
  !!
  pure function secondMom_minus_cs2_2D( cxcx, f, QQ ) result ( m )
    ! --------------------------------------------------------------------------
    integer,       intent(in) :: QQ !< number of discrete directions (=QQ)
    real(kind=rk), intent(in) :: cxcx(6,QQ)   !< discrete velocity of stencil
    !> quantity to which second moment is applied
    real(kind=rk), intent(in) :: f(QQ)
    real(kind=rk)             :: m(3) !< output array
    ! --------------------------------------------------------------------------

    m(1) = sum( (cxcx(1,:) - cs2) * f(:) )
    m(2) = sum( (cxcx(2,:) - cs2) * f(:) )
    m(3) = sum( cxcx(3,:) * f(:) )

  end function secondMom_minus_cs2_2D
! **************************************************************************** !

! **************************************************************************** !
  !> Calculate second moments of some quantity \( f \)
  !! \[
  !!    M_{\alpha\beta} = \sum_{i=1}^{Q} c_{i\alpha} c_{i\beta} f_i
  !! \]
  !! where Q is number of discrete velocity.\n
  !! The output is 1 dimentional array which has 6 componenents.\n
  !! Specifically,
  !! \[ m_1 = \sum_{i=1}^{Q} c_{i1} c_{i1} f_i \]
  !! \[ m_2 = \sum_{i=1}^{Q} c_{i2} c_{i2} f_i \]
  !! \[ m_3 = \sum_{i=1}^{Q} c_{i3} c_{i3} f_i \]
  !! \[ m_4 = \sum_{i=1}^{Q} c_{i1} c_{i2} f_i \]
  !! \[ m_5 = \sum_{i=1}^{Q} c_{i2} c_{i3} f_i \]
  !! \[ m_6 = \sum_{i=1}^{Q} c_{i3} c_{i1} f_i \]
  !! This function is used by shear stress and strain rate.
  !! 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=xz
  !!
  pure function secondMom_minus_cs2_3D( cxcx, f, QQ ) result ( m )
    ! --------------------------------------------------------------------------
    integer,       intent(in) :: QQ !< number of discrete directions (=QQ)
    real(kind=rk), intent(in) :: cxcx(6,QQ)   !< discrete velocity of stencil
    !> quantity to which second moment is applied
    real(kind=rk), intent(in) :: f(QQ)
    real(kind=rk)             :: m(6) !< output array
    ! --------------------------------------------------------------------------

    m(1) = sum( (cxcx(1,:) - cs2) * f(:) )
    m(2) = sum( (cxcx(2,:) - cs2) * f(:) )
    m(3) = sum( (cxcx(3,:) - cs2) * f(:) )
    m(4) = sum( cxcx(4,:) * f(:) )
    m(5) = sum( cxcx(5,:) * f(:) )
    m(6) = sum( cxcx(6,:) * f(:) )

  end function secondMom_minus_cs2_3D
! **************************************************************************** !

end module mus_derivedQuantities_module2
! **************************************************************************** !
