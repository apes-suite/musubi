! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
! Copyright (c) 2019, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!> Boundary condition treatment routines for Nernst-Planck equation
!!
!! A detailed description on the implementation details are given in
!! boundary_implementation
!!
module mus_bc_nernstPlanck_module

  ! include treelm modules
  use env_module,               only: rk, eps
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_construction_module,  only: tem_levelDesc_type
  use tem_debug_module,         only: dbgUnit
  use tem_param_module,         only: cs2inv
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

  public :: moleDens_nonEqExpol_curved
  public :: moleDens_neumann_curved
  public :: moleDens_nonEqExpol
  public :: moleDens_neumann

contains

! **************************************************************************** !
  !> No comment yet!
  !!
  !! @TODO add comment
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine moleDens_nonEqExpol_curved( me, state, bcBuffer, globBC,         &
    &                                    levelDesc, tree, nSize, iLevel,      &
    &                                    sim_time, neigh, layout, fieldProp,  &
    &                                    varPos, nScalars, varSys, derVarPos, &
    &                                    physics, iField, mixture             )
    ! -------------------------------------------------------------------- !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> Level descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) ::layout
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
    real(kind=rk) :: feq_b, fneq_b, moleDens_b
    ! variables for fluid element
    real(kind=rk) :: feq_f, fneq_f, moleDens_f
    ! variables for overnext fluid element
    real(kind=rk) :: feq_ff, fneq_ff, moleDens_ff
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdfTmp(layout%fStencil%QQ)
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    ! potential on surface (link-wise)
    real(kind=rk) :: moleDens_w(me%links(iLevel)%nVals)
    real(kind=rk) :: vel_w(me%links(iLevel)%nVals*3), vel_b(3)
    integer :: bcVel_pos, bcMoleDens_pos, iLink, QQ, QQN
    real(kind=rk) :: ucx
    ! --------------------------------------------------------------------------
    QQ = layout%fstencil%QQ
    QQN = layout%fstencil%QQN

    ! position of boundary moleDensity in varSys
    bcMoleDens_pos = me%bc_states%moleDens%varPos
    ! Get moleDensity on boundary
    call varSys%method%val(bcMoleDens_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleDens                    &
      &           %pntIndex%indexLvl(ilevel)               &
      &           %val(1:me%links(iLevel)%nVals),          &
      & nVals   = me%links(iLevel)%nVals,                  &
      & res     = moleDens_w                               )

    ! convert physical moleDensity into LB moleDensity
    moleDens_w = moleDens_w / physics%moleDens0

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_w                               )

    ! convert physical velocity into LB velocity
    vel_w = vel_w / physics%fac(iLevel)%vel

    do iLink = 1, me%links(iLevel)%nVals

      ! load coefficients
      c_w      = me%nonEqExpol(iLevel)%     c_w(iLink)
      c_f      = me%nonEqExpol(iLevel)%     c_f(iLink)
      c_ff     = me%nonEqExpol(iLevel)%    c_ff(iLink)
      c_neq_f  = me%nonEqExpol(iLevel)% c_neq_f(iLink)
      c_neq_ff = me%nonEqExpol(iLevel)%c_neq_ff(iLink)

      ! link-wise direction
      iDir = me%nonEqExpol(iLevel)%iDir(iLink)

      ! calulate moleDensity of the current element (_f)
      posInBuffer  = me%nonEqExpol(iLevel)%posInBuffer(iLink)
      pdfTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                    ( posInBuffer-1)*nScalars+varPos(QQ) )
      moleDens_f = sum(pdfTmp)

      ! velocity at boundary node
      vel_b = vel_w((iLink-1)*3+1 : iLink*3)

      !> Pre-calculate velocitiy terms
      ucx = dble(layout%fStencil%cxDirRK( 1, iDir ))*vel_b(1) &
        & + dble(layout%fStencil%cxDirRK( 2, iDir ))*vel_b(2) &
        & + dble(layout%fStencil%cxDirRK( 3, iDir ))*vel_b(3)

      ! calulate moleDensity of the overnext fluid neighbor element (_ff)
      posInNeighBuf = me%nonEqExpol(iLevel)%posInNeighBuf(iLink)
      pdf_ff = me%neigh(iLevel)%neighBufferPost(1,                    &
        &             (posInNeighBuf-1)*QQ+1: (posInNeighBuf-1)*QQ+QQ )
      moleDens_ff = sum(pdf_ff)

      ! compute moleDensity on boundary (eq.18)
      moleDens_b = c_w*moleDens_w(iLink) + c_f*moleDens_f + c_ff*moleDens_ff

      ! compute equlibrium (according to eq.17)
      ! fEq_b is on boundary node
      feq_b = layout%weight(iDir)*moleDens_b*(1.0_rk + ucx * cs2inv)
      ! fEq_loc is local element value
      feq_f = layout%weight(iDir)*moleDens_f*(1.0_rk + ucx * cs2inv)
      ! feq_ff is of overnext fluid
      feq_ff = layout%weight(iDir)*moleDens_ff*(1.0_rk + ucx * cs2inv)

      ! use pdf_b = fEq_b + fneq_b (eq.16) to determine
      ! non-equilibrium components (eq.19)
      fneq_f  = pdfTmp(iDir) - feq_f
      fneq_ff = pdf_ff(iDir) - feq_ff
      fneq_b  = c_nEq_f * fneq_f + c_nEq_ff * fneq_ff

      state( me%links(iLevel)%val(iLink) ) = feq_b + fneq_b

    end do !iLink


  end subroutine moleDens_nonEqExpol_curved
! **************************************************************************** !

! **************************************************************************** !
  !> No comment yet!
  !!
  !! @TODO add comment
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine moleDens_neumann_curved( me, state, bcBuffer, globBC, levelDesc, &
    &                                 tree, nSize, iLevel, sim_time, neigh,   &
    &                                 layout, fieldProp, varPos, nScalars,    &
    &                                 varSys, derVarPos, physics, iField,     &
    &                                 mixture                                 )
    ! -------------------------------------------------------------------- !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> Level descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) ::layout
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
    real(kind=rk) :: feq_b, fneq_b, moleDens_b
    ! variables for fluid element
    real(kind=rk) :: feq_f, fneq_f, moleDens_f
    ! variables for overnext fluid element
    real(kind=rk) :: feq_ff, fneq_ff, moleDens_ff
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdfTmp(layout%fStencil%QQ)
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    ! mole flux on surface (link-wise)
    real(kind=rk) :: moleFlux_w(me%links(iLevel)%nVals)
    real(kind=rk) :: vel_w(me%links(iLevel)%nVals*3), vel_b(3)
    integer :: bcVel_pos,iLink, QQ
    integer :: bcMoleFlux_pos
    real(kind=rk) :: normal(3), moleFlux_b
    real(kind=rk) :: ucx
    ! --------------------------------------------------------------------------
    !write(*,*) 'bclabel ', trim(me%label)
    QQ = layout%fStencil%QQ

    ! position of boundary mole flux variable in varSys
    bcMoleFlux_pos = me%bc_states%moleFlux%varPos
    ! Get moleDensity on boundary
    call varSys%method%val(bcmoleFlux_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%moleFLux               &
      &           %pntIndex%indexLvl(ilevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = me%links(iLevel)%nVals,             &
      & res     = moleFlux_w                          )

    ! convert physical moleflux into Lattice unit
    moleFlux_w = moleFlux_w / physics%fac(iLevel)%moleflux

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_w                               )

    ! convert physical velocity into LB velocity
    vel_w = vel_w / physics%fac(iLevel)%vel

    do iLink = 1, me%links(iLevel)%nVals

      ! load coefficients
      c_neq_f  = me%nonEqExpol(iLevel)% c_neq_f(iLink)
      c_neq_ff = me%nonEqExpol(iLevel)%c_neq_ff(iLink)

      ! link-wise direction
      iDir = me%nonEqExpol(iLevel)%iDir(iLink)

      ! calulate potential of the current element (_f)
      posInBuffer  = me%nonEqExpol(iLevel)%posInBuffer(iLink)
      pdfTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                    ( posInBuffer-1)*nScalars+varPos(QQ) )
      moleDens_f = sum(pdfTmp)

      ! velocity at boundary node
      vel_b = vel_w((iLink-1)*3+1 : iLink*3)

     !> Pre-calculate velocitiy terms
      ucx = dble(layout%fStencil%cxDirRK( 1, iDir ))*vel_b(1) &
        & + dble(layout%fStencil%cxDirRK( 2, iDir ))*vel_b(2) &
        & + dble(layout%fStencil%cxDirRK( 3, iDir ))*vel_b(3)


      ! calulate potential of the overnext fluid neighbor element (_ff)
      posInNeighBuf = me%nonEqExpol(iLevel)%posInNeighBuf(iLink)
      pdf_ff = me%neigh(iLevel)%neighBufferPost(1,                    &
        &             (posInNeighBuf-1)*QQ+1: (posInNeighBuf-1)*QQ+QQ )
      moleDens_ff = sum(pdf_ff)

      !write(*,*) 'normal ', globBC%elemLvl(iLevel)%normal%val(:, me%nonEqExpol(iLevel)%posInBCelems(iLink))
      !write(*,*) 'normal Ind', globBC%elemLvl(iLevel)%normalInd%val(me%nonEqExpol(iLevel)%posInBCelems(iLink))
      normal = globBC%elemLvl(iLevel)%normal                          &
        &            %val(:, me%nonEqExpol(iLevel)%posInBCelems(iLink))
      moleFlux_b = moleFlux_w(iLink)                          &
        & * dot_product(layout%fStencil%cxDirRK(:, iDir), normal)

      ! calulate moleDensity on boundary
      moleDens_b = ( 4.0_rk*moleDens_f - moleDens_ff - 2.0_rk*moleFlux_b ) &
        &        / 3.0_rk

      ! compute equlibrium (according to eq.17)
      ! fEq_b is on boundary
      fEq_b = layout%weight(iDir)*moleDens_b*(1.0_rk + ucx * cs2inv)
      ! fEq_loc is local element value
      fEq_f = layout%weight(iDir)*moleDens_f*(1.0_rk + ucx * cs2inv)
      ! feq is of overnext fluid
      fEq_ff = layout%weight(iDir)*moleDens_ff*(1.0_rk + ucx * cs2inv)

      ! use pdf_b = fEq_b + fneq_b (eq.16) to determine
      ! non-equilibrium components (eq.19)
      fneq_f  = pdfTmp(iDir) - feq_f
      fneq_ff = pdf_ff(iDir) - feq_ff
      fneq_b  = c_nEq_f * fneq_f + c_nEq_ff * fneq_ff

      state( me%links(iLevel)%val(iLink) ) = feq_b + fneq_b

    end do !iLink


  end subroutine moleDens_neumann_curved
! **************************************************************************** !

! **************************************************************************** !
  !> No comment yet!
  !!
  !! @TODO add comment
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine moleDens_nonEqExpol( me, state, bcBuffer, globBC, levelDesc,      &
    &                             tree, nSize, iLevel, sim_time, neigh,        &
    &                             layout, fieldProp, varPos, nScalars, varSys, &
    &                             derVarPos, physics, iField, mixture          )
    ! -------------------------------------------------------------------- !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> Level descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) ::layout
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
    integer :: iDir, posInBuffer
    ! variables for fictious boundary element
    real(kind=rk) :: feq_b, moleDens_b
    ! variables for overnext fluid element
    real(kind=rk) :: feq_f, moleDens_f
    real(kind=rk) :: inv_moleDens
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdfTmp(layout%fStencil%QQ)
    ! potential on surface (link-wise)
    real(kind=rk) :: moleDens_w(me%links(iLevel)%nVals)
    real(kind=rk) :: vel_w(me%links(iLevel)%nVals*3), vel_b(3)
    integer :: bcVel_pos, bcMoleDens_pos, iLink, QQ, QQN
    real(kind=rk) :: ucx
    ! --------------------------------------------------------------------------
    QQ = layout%fstencil%QQ
    QQN = layout%fstencil%QQN
    inv_moleDens = 1._rk / mixture%moleDens0LB

    ! position of boundary moleDensity in varSys
    bcMoleDens_pos = me%bc_states%moleDens%varPos
    ! Get potential_phy on boundary
    call varSys%method%val(bcMoleDens_pos)%get_valOfIndex( &
      & varSys  = varSys,                                  &
      & time    = sim_time,                                &
      & iLevel  = iLevel,                                  &
      & idx     = me%bc_states%moleDens                    &
      &           %pntIndex%indexLvl(ilevel)               &
      &           %val(1:me%links(iLevel)%nVals),          &
      & nVals   = me%links(iLevel)%nVals,                  &
      & res     = moleDens_w                               )

    ! convert physical moleDensity into LB potential
    moleDens_w = moleDens_w * inv_moleDens

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_w                               )

    ! convert physical velocity into LB velocity
    vel_w = vel_w / physics%fac(iLevel)%vel


    do iLink = 1, me%links(iLevel)%nVals

      ! link-wise direction
      iDir = me%nonEqExpol(iLevel)%iDir(iLink)

      ! calulate moleDensity of the current element (_f)
      posInBuffer  = me%nonEqExpol(iLevel)%posInBuffer(iLink)
      pdfTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                    ( posInBuffer-1)*nScalars+varPos(QQ) )
      moleDens_f = sum(pdfTmp)

      ! velocity at boundary node
      vel_b = vel_w((iLink-1)*3+1 : iLink*3)

      !> Pre-calculate velocitiy terms
      ucx = dble(layout%fStencil%cxDirRK( 1, iDir ))*vel_b(1) &
        & + dble(layout%fStencil%cxDirRK( 2, iDir ))*vel_b(2) &
        & + dble(layout%fStencil%cxDirRK( 3, iDir ))*vel_b(3)

      ! compute moleDensity on boundary (eq.18)
      moleDens_b = moleDens_w(iLink)

      ! velocity at boundary node
      vel_b = vel_w((iLink-1)*3+1 : iLink*3)

      ! compute equlibrium (according to eq.17)
      ! fEq_b is on boundary node
      feq_b = layout%weight(iDir)*moleDens_b*(1.0_rk + ucx * cs2inv)

      ! feq_ff is of overnext fluid
      feq_f = layout%weight(iDir)*moleDens_f*(1.0_rk + ucx * cs2inv)

      ! use pdf_b = fEq_b + fneq_b (eq.16) to determine
      ! non-equilibrium components (eq.19)
      state( me%links(iLevel)%val(iLink) ) = feq_b + (pdfTmp(iDir) - feq_f)

    end do !iLink


  end subroutine moleDens_nonEqExpol
! **************************************************************************** !

! **************************************************************************** !
  !> No comment yet!
  !!
  !! @TODO add comment
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine moleDens_neumann( me, state, bcBuffer, globBC, levelDesc, tree,   &
    &                          nSize, iLevel, sim_time, neigh, layout,         &
    &                          fieldProp, varPos, nScalars, varSys, derVarPos, &
    &                          physics, iField, mixture                        )
    ! -------------------------------------------------------------------- !
    !> global boundary type
    class(boundary_type) :: me
    !> Current state vector of iLevel
    real(kind=rk), intent(inout) :: state(:)
    !> size of state array ( in terms of elements )
    integer, intent(in) :: nSize
    !> state values of boundary elements of all fields of iLevel
    real(kind=rk), intent(in) :: bcBuffer(:)
    !> Level descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc
    !> Treelm Mesh
    type(treelmesh_type), intent(in) :: tree
    !> fluid parameters and properties
    type(mus_field_prop_type), intent(in) :: fieldProp
    !> stencil layout information
    type(mus_scheme_layout_type), intent(in) ::layout
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
    integer :: iDir, posInBuffer, posInNeighBuf
    ! variables for fictious boundary element
    real(kind=rk) :: feq_b, fneq_b, moleDens_b
    real(kind=rk) :: feq_f, fneq_f, moleDens_f
    ! variables for overnext fluid element
    real(kind=rk) :: feq_ff, fneq_ff, moleDens_ff
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdfTmp(layout%fStencil%QQ)
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    ! molefluxy on surface (link-wise)
    real(kind=rk) :: moleFlux_w(me%links(iLevel)%nVals)
    integer :: bcVel_pos,iLink, QQ
    integer :: bcMoleFlux_pos
    real(kind=rk) :: normal(3), moleFlux_b
    real(kind=rk) :: vel_w(me%links(iLevel)%nVals*3), vel_b(3)
    real(kind=rk) ucx
    ! --------------------------------------------------------------------------
    !write(*,*) 'bclabel ', trim(me%label)
    QQ = layout%fStencil%QQ

    ! position of boundary moleflux variable in varSys
    bcMoleFlux_pos = me%bc_states%moleFlux%varPos
    ! Get potential_phy on boundary
    call varSys%method%val(bcMoleFlux_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%moleFlux               &
      &           %pntIndex%indexLvl(ilevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = me%links(iLevel)%nVals,             &
      & res     = moleFlux_w                          )

    ! convert physical moleflux into Lattice unit
    moleFlux_w = moleFlux_w / physics%fac(iLevel)%moleflux

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! Get velocity
    call varSys%method%val(bcVel_pos)%get_valOfIndex( &
      & varSys  = varSys,                             &
      & time    = sim_time,                           &
      & iLevel  = iLevel,                             &
      & idx     = me%bc_states%velocity               &
      &           %pntIndex%indexLvl(iLevel)          &
      &           %val(1:me%links(iLevel)%nVals),     &
      & nVals   = globBC%nElems(iLevel),              &
      & res     = vel_w                               )

    ! convert physical velocity into LB velocity
    vel_w = vel_w / physics%fac(iLevel)%vel


    do iLink = 1, me%links(iLevel)%nVals

      ! link-wise direction
      iDir = me%nonEqExpol(iLevel)%iDir(iLink)

      ! calulate potential of the current element (_f)
      posInBuffer  = me%nonEqExpol(iLevel)%posInBuffer(iLink)
      pdfTmp(1:QQ) = bcBuffer( (posInBuffer-1)*nScalars+varPos(1) :  &
        &                    ( posInBuffer-1)*nScalars+varPos(QQ) )
      moleDens_f = sum(pdfTmp)

      ! velocity at boundary node
      vel_b = vel_w((iLink-1)*3+1 : iLink*3)

      !> Pre-calculate velocitiy terms
      ucx = dble(layout%fStencil%cxDirRK( 1, iDir ))*vel_b(1) &
        & + dble(layout%fStencil%cxDirRK( 2, iDir ))*vel_b(2) &
        & + dble(layout%fStencil%cxDirRK( 3, iDir ))*vel_b(3)

      ! calulate potential of the overnext fluid neighbor element (_ff)
      posInNeighBuf = me%nonEqExpol(iLevel)%posInNeighBuf(iLink)
      pdf_ff = me%neigh(iLevel)%neighBufferPost(1,                    &
        &             (posInNeighBuf-1)*QQ+1: (posInNeighBuf-1)*QQ+QQ )
      moleDens_ff = sum(pdf_ff)

      !write(*,*) 'normal ', globBC%elemLvl(iLevel)%normal%val(:, me%nonEqExpol(iLevel)%posInBCelems(iLink))
      !write(*,*) 'normal Ind', globBC%elemLvl(iLevel)%normalInd%val(me%nonEqExpol(iLevel)%posInBCelems(iLink))
      normal = globBC%elemLvl(iLevel)%normal                          &
        &            %val(:, me%nonEqExpol(iLevel)%posInBCelems(iLink))
      moleFlux_b = moleFlux_w(iLink) &
        & * dot_product(layout%fStencil%cxDirRK(:, iDir), normal)

      ! calulate moleDens on boundary
      moleDens_b = ( 4.0_rk*moleDens_f - moleDens_ff - 2.0_rk*moleFlux_b) &
        &        / 3.0_rk

      ! compute equlibrium (according to eq.17)
      ! fEq_b is on boundary
      fEq_b = layout%weight(iDir)*moleDens_b*(1.0_rk + ucx * cs2inv)

      ! feq is of fluid
      fEq_f = layout%weight(iDir)*moleDens_f*(1.0_rk + ucx * cs2inv)

      ! feq is of overnext fluid
      fEq_ff = layout%weight(iDir)*moleDens_ff*(1.0_rk + ucx * cs2inv)

      ! use pdf_b = fEq_b + fneq_b (eq.16) to determine
      ! non-equilibrium components (eq.19)
      fneq_f  = pdfTmp(iDir) - feq_f
      fneq_ff = pdf_ff(iDir) - feq_ff
      ! KM: For straight boundary, qVal is 0.5 in consistent with
      ! second order extrapolation of potential
      fneq_b  = 0.5_rk * fneq_f + 0.5_rk * fneq_ff

      state( me%links(iLevel)%val(iLink) ) = feq_b + fneq_b

    end do !iLink

  end subroutine moleDens_neumann
! **************************************************************************** !

end module mus_bc_nernstPlanck_module
! **************************************************************************** !
