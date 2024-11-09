! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2017-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2018 Jana Gericke <jana.gericke@uni-siegen.de>
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
! ****************************************************************************** !
!> Boundary condition treatment routines for Poisson equation
!!
!! A detailed description on the implementation details are given in
!! boundary_implementation
!!
module mus_bc_poisson_module

  ! include treelm modules
  use env_module,               only: rk, eps
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_debug_module,         only: dbgUnit

  ! include musubi modules
  use mus_bc_header_module,      only: boundary_type, glob_boundary_type
  use mus_scheme_layout_module,  only: mus_scheme_layout_type
  use mus_field_prop_module,     only: mus_field_prop_type
  use mus_derVarPos_module,      only: mus_derVarPos_type
  use mus_param_module,          only: mus_param_type
  use mus_physics_module,        only: mus_physics_type
  use mus_mixture_module,        only: mus_mixture_type

  implicit none

  private

  public :: potential_nonEqExpol
  public :: potential_nonEqExpol_curved
  public :: potential_neumann
  public :: potential_neumann_curved

contains

! ****************************************************************************** !
  !> Linkwise Dirichlet potential non-equilibrium boundary condition for curved
  !> and straight walls. For straight wall, physical boundary overlaps with
  !! boundary node i.e. qVal=0.0.
  !!
  !! The pdf is decomposed into equilibrium (eq) and non-equilibrium (neq) part:
  !! f = f_eq(x_b,t) + f_neq(x_f,t)
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'inner',
  !!    kind = 'potential_noneq_expol',
  !!    potential = pot_inner,
  !!  }
  !!}
  !!```
  !! For straight boundaries, qVal=1.0:
  !! Luo, K., Wu, J., Yi, H., & Tan, H. (2016). Lattice Boltzmann model for
  !! Coulomb-driven flows in dielectric liquids dielectric liquids.
  !! Physical Review E, 23309(93), 1–11.
  !! http://doi.org/10.1103/PhysRevE.93.023309
  ! subroutine potential_nonEqExpol
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine potential_nonEqExpol( me, state, bcBuffer, globBC, levelDesc,     &
    &                              tree, nSize, iLevel, sim_time, neigh,       &
    &                              layout, fieldProp, varPos, nScalars,        &
    &                              varSys, derVarPos, physics, iField, mixture )
    ! -------------------------------------------------------------------- !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    ! also determined in mus_bc_header_module
    integer :: iDir
    ! variables for fictious boundary element
    real(kind=rk) :: feq_b
    ! variables for fluid element
    real(kind=rk) :: feq_ff, pot_f, pot_ff, pot_b
    real(kind=rk) :: inv_pot
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdf_f(layout%fStencil%QQ)
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    ! potential on surface (link-wise)
    real(kind=rk) :: pot_w(globBC%nElems(iLevel))
    integer :: bcPot_pos, iLink, QQ, QQN
    ! ---------------------------------------------------------------------------
    QQ = layout%fstencil%QQ
    QQN = layout%fstencil%QQN
    inv_pot = 1._rk / physics%fac(iLevel)%potential

    ! position of boundary potential in varSys
    bcPot_pos = me%bc_states%potential%varPos
    ! Get potential_phy on boundary
    call varSys%method%val(bcPot_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%potential              &
      &           %pntIndex%indexLvl(ilevel)          &
      &           %val(1:globBC%nElems(iLevel)),      &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = pot_w                               )

    ! convert physical potential into LB potential
    pot_w = pot_w * inv_pot

    associate( neighBufferPost => me%neigh(iLevel)%neighBufferPost,  &
      &        posInNeighBuf => me%nonEqExpol(iLevel)%posInNeighBuf, &
      &        posInBCelems => me%nonEqExpol(iLevel)%posInBCelems,   &
      &        posInBuffer => me%nonEqExpol(iLevel)%posInBuffer      )

      do iLink = 1, me%links(iLevel)%nVals
        ! link-wise direction
        iDir = me%nonEqExpol(iLevel)%iDir(iLink)

        ! calulate potential of the current element (_f)
        pdf_f(1:QQ) = bcBuffer( (posInBuffer(iLink)-1)*nScalars+varPos(1)  &
          &                   : (posInBuffer(iLink)-1)*nScalars+varPos(QQ) )
        pot_f = sum(pdf_f)

        ! determine needed quantities of the overnext fluid neighbor element
        ! x_ff
        pdf_ff = neighBufferPost(1, (posInNeighBuf(iLink)-1)*QQ+1  &
          &                       : (posInNeighBuf(iLink)-1)*QQ+QQ )
        pot_ff = sum(pdf_ff)

        ! fEq_loc is local element value
        feq_ff = layout%weight(iDir)*pot_ff

        ! compute potential on boundary (eq.18)
        pot_b = 3.0_rk*pot_w(posInBCelems(iLink)) - pot_f - pot_ff

        ! compute equlibrium (according to eq.17)
        ! fEq_b is on boundary node
        feq_b = layout%weight(iDir)*pot_b

        ! write into state
        ! feq_b + fneq_ff
        ! fneq is computed from post-collision pdf so no need to multiply
        ! (1-omega) as in the literature
        state( me%links(iLevel)%val(iLink) ) = feq_b + (pdf_ff(iDir) - feq_ff)

      end do !iLink
    end associate

  end subroutine potential_nonEqExpol
! ****************************************************************************** !

! ****************************************************************************** !
  !> Linkwise Dirichlet potential non-equilibrium boundary condition for curved
  !! wall
  !!
  !! The pdf is decomposed into equilibrium (eq) and non-equilibrium (neq) part:
  !! f = f_eq + f_neq
  !! - f_eq is calculated by weighting a fictitious potential, which is obtained
  !!   by an extrapolation using the fluid neighbor(s)
  !! - f_neq is approximated by second-order extrapolation using the fluid
  !!   neighbor(s)
  !! - for qVal < 0.75 even the second neighbor is used for the extrapolations
  !! - Dirichlet: potential on wall is directly known
  !! - Linkwise: as much as possible is outsourced to the subroutine
  !!             "mus_set_nonEqExpol"
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'inner',
  !!    kind = 'potential_noneq_expol',
  !!    potential = pot_inner,
  !!    curved = true
  !!  }
  !!}
  !!```
  !! This is described in the paper:
  !! Luo, K.; Wu, J.; Yi HL. & Tan HP. (2016). A lattice Boltzmann method for
  !! electric field-space charge coupled problems.
  !! Proceedings of the 2016 Electrostatics Joint Conference (June 2016).
  !!
  ! subroutine potential_nonEqExpol
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine potential_nonEqExpol_curved( me, state, bcBuffer, globBC,         &
    &                                     levelDesc, tree, nSize, iLevel,      &
    &                                     sim_time, neigh, layout, fieldProp,  &
    &                                     varPos, nScalars, varSys, derVarPos, &
    &                                     physics, iField, mixture             )
    ! -------------------------------------------------------------------- !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    ! coefficients which are calculated in mus_bc_header_module
    real(kind=rk) :: c_w, c_f, c_ff, c_neq_f, c_neq_ff
    ! also determined in mus_bc_header_module
    integer :: iDir, posInBuffer, posInNeighBuf
    ! variables for fictious boundary element
    real(kind=rk) :: feq_b, fneq_b, pot_b
    ! variables for fluid element
    real(kind=rk) :: feq_f, fneq_f, pot_f
    ! variables for overnext fluid element
    real(kind=rk) :: feq_ff, fneq_ff, pot_ff
    real(kind=rk) :: inv_pot
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdf_f(layout%fStencil%QQ)
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    ! potential on surface (link-wise)
    real(kind=rk) :: pot_w(me%links(iLevel)%nVals)
    integer :: bcPot_pos, iLink, QQ, QQN
    ! ---------------------------------------------------------------------------
    QQ = layout%fstencil%QQ
    QQN = layout%fstencil%QQN
    inv_pot = 1._rk / physics%fac(iLevel)%potential

    ! position of boundary potential in varSys
    bcPot_pos = me%bc_states%potential%varPos
    ! Get potential_phy on boundary
    call varSys%method%val(bcPot_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%potential              &
      &           %pntIndex%indexLvl(ilevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = me%links(iLevel)%nVals,             &
      & res     = pot_w                               )

    ! convert physical potential into LB potential
    pot_w = pot_w * inv_pot

    do iLink = 1, me%links(iLevel)%nVals

      ! load coefficients
      c_w      = me%nonEqExpol(iLevel)%     c_w(iLink)
      c_f      = me%nonEqExpol(iLevel)%     c_f(iLink)
      c_ff     = me%nonEqExpol(iLevel)%    c_ff(iLink)
      c_neq_f  = me%nonEqExpol(iLevel)% c_neq_f(iLink)
      c_neq_ff = me%nonEqExpol(iLevel)%c_neq_ff(iLink)

      ! link-wise direction
      iDir = me%nonEqExpol(iLevel)%iDir(iLink)

      ! calulate potential of the current element (_f)
      posInBuffer  = me%nonEqExpol(iLevel)%posInBuffer(iLink)
      pdf_f(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                    ( posInBuffer-1)*nScalars+varPos(QQ) )
      pot_f = sum(pdf_f)

      ! calulate potential of the overnext fluid neighbor element (_ff)
      posInNeighBuf = me%nonEqExpol(iLevel)%posInNeighBuf(iLink)
      pdf_ff = me%neigh(iLevel)%neighBufferPost(1,                    &
        &             (posInNeighBuf-1)*QQ+1: (posInNeighBuf-1)*QQ+QQ )
      pot_ff = sum(pdf_ff)

      ! compute potential on boundary (eq.18)
      pot_b = c_w*pot_w(iLink) + c_f*pot_f + c_ff*pot_ff

      ! compute equlibrium (according to eq.17)
      ! fEq_b is on boundary node
      feq_b = layout%weight(iDir)*pot_b
      ! fEq_loc is local element value
      feq_f = layout%weight(iDir)*pot_f
      ! feq_ff is of overnext fluid
      feq_ff = layout%weight(iDir)*pot_ff

      ! use pdf_b = fEq_b + fneq_b (eq.16) to determine
      ! non-equilibrium components (eq.19)
      fneq_f  = pdf_f(iDir) - feq_f
      fneq_ff = pdf_ff(iDir) - feq_ff
      fneq_b  = c_nEq_f * fneq_f + c_nEq_ff * fneq_ff


      state( me%links(iLevel)%val(iLink) ) = feq_b + fneq_b

    end do !iLink

  end subroutine potential_nonEqExpol_curved
! ****************************************************************************** !

! ****************************************************************************** !
  !> Linkwise neumann potential non-equilibrium boundary condition for curved
  !> and straight walls (zero gradient).
  !! For straight wall, values are extrapolated along
  !! boundary normal instead of along the link. The accuracy of straight wall
  !! depends on the qVal defined in config file and default is set to 0.5
  !!
  !! The pdf is decomposed into equilibrium (eq) and non-equilibrium (neq) part:
  !! - f_eq is calculated by weighting a fictitious potential, which is obtained
  !!   by an extrapolation using the fluid neighbor(s)
  !! - f_neq is approximated by second-order extrapolation using the fluid
  !!   neighbor(s)
  !! - for qVal < 0.75 even the second neighbor is used for the extrapolations
  !! - boundary potential (pot_b) has to be extrapolated using two fluid neighb:
  !!   pot_b = (4*pot_f - pot_ff)/3 (typo in paper)
  !! - Linkwise: as much as possible is outsourced to the subroutine
  !!             "mus_set_nonEqExpol"
  !!
  !! NOTE: - possibility to extend the equation to extraplotate pot_b by a
  !!         gradient if a non-zero one is desired
  !!         (see: Huang H, Lee T S, and Shu C. “Thermal curved boundary
  !!               treatment for the thermal lattice Boltzmann equation,”
  !!               Int. J. Mod. Phys. C, vol. 17(05), pp. 631-643, 2006)
  !!       - more accurate schemes to determine pot_b can be found in:
  !!         Chen Q, Zhang X, and Zhang J. “Improved treatments for general
  !!         boundary conditions in the lattice Boltzmann method for convection-
  !!         diffusion and heat transfer processes,” Phys Rev E, vol. 88(3),
  !!         033304, 2013
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'outer',
  !!    kind = 'potential_neumann_link',
  !!    curved = true
  !!  }
  !!}
  !!```
  !! This is described in the paper:
  !! Luo, K.; Wu, J.; Yi HL. & Tan HP. (2016). A lattice Boltzmann method for
  !! electric field-space charge coupled problems.
  !! Proceedings of the 2016 Electrostatics Joint Conference (June 2016).
  !!
  ! subroutine potential_neumann
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine potential_neumann( me, state, bcBuffer, globBC, levelDesc, tree, &
    &                           nSize, iLevel, sim_time, neigh, layout,       &
    &                           fieldProp, varPos, nScalars, varSys,          &
    &                           derVarPos, physics, iField, mixture           )
    ! -------------------------------------------------------------------- !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    ! also determined in mus_bc_header_module
    integer :: iDir
    ! variables for fictious boundary element
    real(kind=rk) :: feq_b, fneq_b, pot_b
    ! variables for fluid element
    real(kind=rk) :: feq_f, fneq_f, pot_f
    ! variables for overnext fluid element
    real(kind=rk) :: feq_ff, fneq_ff, pot_ff
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdf_f(layout%fStencil%QQ)
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    ! surface charge density on surface (link-wise)
    real(kind=rk) :: surChargeDens_w(me%links(iLevel)%nVals)
    integer :: iLink, QQ
    integer :: bcSCD_pos
    real(kind=rk) :: inv_permittivity
    real(kind=rk) :: normal(3), surChargeDens_fac
    ! ---------------------------------------------------------------------------
    !write(*,*) 'bclabel ', trim(me%label)
    QQ = layout%fStencil%QQ

    ! position of boundary surface charge density variable in varSys
    bcSCD_pos = me%bc_states%surChargeDens%varPos
    ! Get potential_phy on boundary
    call varSys%method%val(bcSCD_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%surChargeDens          &
      &           %pntIndex%indexLvl(ilevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = me%links(iLevel)%nVals,             &
      & res     = surChargeDens_w                     )

    ! convert physical surface charge density into Lattice unit
    surChargeDens_w = surChargeDens_w * physics%dxLvl(iLevel)**2 &
      &             / physics%coulomb0

    inv_permittivity = 1.0_rk / fieldProp%poisson%permittivity

    associate( neighBufferPost => me%neigh(iLevel)%neighBufferPost,  &
      &        posInNeighBuf => me%nonEqExpol(iLevel)%posInNeighBuf, &
      &        posInBCelems => me%nonEqExpol(iLevel)%posInBCelems,   &
      &        posInBuffer => me%nonEqExpol(iLevel)%posInBuffer      )

      do iLink = 1, me%links(iLevel)%nVals

        ! link-wise direction
        iDir = me%nonEqExpol(iLevel)%iDir(iLink)

        ! calulate potential of the current element (_f)
        pdf_f(1:QQ) = bcBuffer( (posInBuffer(iLink)-1)*nScalars+varPos(1)  &
          &                   : (posInBuffer(iLink)-1)*nScalars+varPos(QQ) )
        pot_f = sum(pdf_f)

        ! calulate potential of the overnext fluid neighbor element (_ff)
        pdf_ff = neighBufferPost(1, (posInNeighBuf(iLink)-1)*QQ+1  &
          &                       : (posInNeighBuf(iLink)-1)*QQ+QQ )
        pot_ff = sum(pdf_ff)

        normal = globBC%elemLvl(iLevel)%normal                          &
          &            %val(:, me%nonEqExpol(iLevel)%posInBCelems(iLink))
        surChargeDens_fac = -surChargeDens_w(iLink)                       &
          &               * inv_permittivity                              &
          &               * dot_product(layout%fStencil%cxDirRK(:, iDir), &
          &                             normal)

        ! calulate potential on boundary
        pot_b = ( 4.0_rk*pot_f - pot_ff - 2.0_rk*surChargeDens_fac ) / 3.0_rk

        ! compute equlibrium (according to eq.17)
        ! fEq_b is on boundary
        feq_b = layout%weight(iDir)*pot_b
        ! fEq_loc is local element value
        feq_f = layout%weight(iDir)*pot_f
        ! feq is of overnext fluid
        fEq_ff = layout%weight(iDir)*pot_ff

        ! use pdf_b = fEq_b + fneq_b (eq.16) to determine
        ! non-equilibrium components (eq.19)
        fneq_f  = pdf_f(iDir) - feq_f
        fneq_ff = pdf_ff(iDir) - feq_ff
        ! KM: For straight boundary, qVal is 0.5 in consistent with
        ! second order extrapolation of potential
        fneq_b  = 0.5_rk * fneq_f + 0.5_rk * fneq_ff

        state( me%links(iLevel)%val(iLink) ) = feq_b + fneq_b

      end do !iLink
    end associate

  end subroutine potential_neumann
! ****************************************************************************** !

! ****************************************************************************** !
  !> No comment yet!
  !!
  !! @TODO add comment
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine potential_neumann_curved( me, state, bcBuffer, globBC, levelDesc, &
    &                                  tree, nSize, iLevel, sim_time, neigh,   &
    &                                  layout, fieldProp, varPos, nScalars,    &
    &                                  varSys, derVarPos, physics, iField,     &
    &                                  mixture                                 )
    ! -------------------------------------------------------------------- !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> iLevel descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) :: layout
    !> the level On which this boundary was invoked
    integer, intent(in) :: iLevel
    !> connectivity array corresponding to state vector
    integer, intent(in) :: neigh(:)
    !> global time information
    type(tem_time_type), intent(in)  :: sim_time
    !> pointer to field variable in the state vector
    integer, intent(in) :: varPos(:)
    !> number of Scalars in the scheme var system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type(tem_varSys_type), intent(in) :: varSys
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme global boundary type
    type(glob_boundary_type), intent(in) :: globBC
    !> scheme global boundary type
    type(mus_physics_type), intent(in) :: physics
    !> current field
    integer, intent(in) :: iField
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! -------------------------------------------------------------------- !
    ! coefficients which are calculated in mus_bc_header_module
    real(kind=rk) :: c_neq_f, c_neq_ff
    ! also determined in mus_bc_header_module
    integer :: iDir, posInBuffer, posInNeighBuf
    ! variables for fictious boundary element
    real(kind=rk) :: feq_b, fneq_b, pot_b
    ! variables for fluid element
    real(kind=rk) :: feq_f, fneq_f, pot_f
    ! variables for overnext fluid element
    real(kind=rk) :: feq_ff, fneq_ff, pot_ff
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdf_f(layout%fStencil%QQ)
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    ! surface charge density on surface (link-wise)
    real(kind=rk) :: surChargeDens_w(me%links(iLevel)%nVals)
    integer :: iLink, QQ
    integer :: bcSCD_pos
    real(kind=rk) :: inv_permittivity
    real(kind=rk) :: normal(3), surChargeDens_fac
    ! ---------------------------------------------------------------------------
    !write(*,*) 'bclabel ', trim(me%label)
    QQ = layout%fStencil%QQ

    ! position of boundary surface charge density variable in varSys
    bcSCD_pos = me%bc_states%surChargeDens%varPos
    ! Get potential_phy on boundary
    call varSys%method%val(bcSCD_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%surChargeDens          &
      &           %pntIndex%indexLvl(ilevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = me%links(iLevel)%nVals,             &
      & res     = surChargeDens_w                     )

    ! convert physical surface charge density into Lattice unit
    surChargeDens_w = surChargeDens_w * physics%dxLvl(iLevel)**2 &
      &             / physics%coulomb0
    inv_permittivity = 1.0_rk / fieldProp%poisson%permittivity

    do iLink = 1, me%links(iLevel)%nVals

      ! load coefficients
      c_neq_f  = me%nonEqExpol(iLevel)% c_neq_f(iLink)
      c_neq_ff = me%nonEqExpol(iLevel)%c_neq_ff(iLink)

      ! link-wise direction
      iDir = me%nonEqExpol(iLevel)%iDir(iLink)

      ! calulate potential of the current element (_f)
      posInBuffer  = me%nonEqExpol(iLevel)%posInBuffer(iLink)
      pdf_f(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                    ( posInBuffer-1)*nScalars+varPos(QQ) )
      pot_f = sum(pdf_f)

      ! calulate potential of the overnext fluid neighbor element (_ff)
      posInNeighBuf = me%nonEqExpol(iLevel)%posInNeighBuf(iLink)
      pdf_ff = me%neigh(iLevel)%neighBufferPost(1,                    &
        &             (posInNeighBuf-1)*QQ+1: (posInNeighBuf-1)*QQ+QQ )
      pot_ff = sum(pdf_ff)

      !write(*,*) 'normal ', globBC%elemLvl(iLevel)%normal%val(:, me%nonEqExpol(iLevel)%posInBCelems(iLink))
      !write(*,*) 'normal Ind', globBC%elemLvl(iLevel)%normalInd%val(me%nonEqExpol(iLevel)%posInBCelems(iLink))
      normal = globBC%elemLvl(iLevel)%normal                          &
        &            %val(:, me%nonEqExpol(iLevel)%posInBCelems(iLink))
      surChargeDens_fac = -surChargeDens_w(iLink)*inv_permittivity &
        & * dot_product(layout%fStencil%cxDirRK(:, iDir), normal)

      ! calulate potential on boundary
      pot_b = ( 4.0_rk*pot_f - pot_ff - 2.0_rk*surChargeDens_fac ) / 3.0_rk

      ! compute equlibrium (according to eq.17)
      ! fEq_b is on boundary
      fEq_b = layout%weight(iDir)*pot_b
      ! fEq_loc is local element value
      fEq_f = layout%weight(iDir)*pot_f
      ! feq is of overnext fluid
      fEq_ff = layout%weight(iDir)*pot_ff

      ! use pdf_b = fEq_b + fneq_b (eq.16) to determine
      ! non-equilibrium components (eq.19)
      fneq_f  = pdf_f(iDir) - feq_f
      fneq_ff = pdf_ff(iDir) - feq_ff
      fneq_b  = c_nEq_f * fneq_f + c_nEq_ff * fneq_ff

      state( me%links(iLevel)%val(iLink) ) = feq_b + fneq_b

    end do !iLink

  end subroutine potential_neumann_curved
! ****************************************************************************** !

end module mus_bc_poisson_module
! ****************************************************************************** !
