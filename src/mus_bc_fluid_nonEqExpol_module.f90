! Copyright (c) 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2021-2023 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
!> author: Kannan Masilamani
!! Non-equilibrium extrapolatioon boundary condition treatment routines for
!! fluid simulation. The boundary condition for straight boundary with qVal=0
!! is based on: Guo, Z., & Shi, B. (2002). "Non-equilibrium extrapolation method
!! for velocity and pressure boundary conditions in the lattice Boltzmann
!! method." Chinese Physics.
!! The boundary condition for curved boundary is based on
!! Guo, Z., Shi, B., & Zheng, C. (2002). An extrapolation method for boundary
!! conditions in lattice Boltzmann method. Physics of Fluids, 14(May), 10–14.
!!
!! Details on implementation can be found in [[tem_bc_module]]

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
module mus_bc_fluid_nonEqExpol_module

  ! include treelm modules
  use env_module,               only: rk
  use tem_param_module,         only: div1_3, div4_3, csInv, cs, cs2, cs2inv, &
    &                                 rho0, rho0Inv
  use tem_time_module,          only: tem_time_type
  use treelmesh_module,         only: treelmesh_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_construction_module,  only: tem_levelDesc_type

  ! include musubi modules
  use mus_bc_header_module,          only: boundary_type, &
    &                                      glob_boundary_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_field_prop_module,         only: mus_field_prop_type
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_physics_module,            only: mus_physics_type
  use mus_mixture_module,            only: mus_mixture_type

  implicit none

  private

  public :: velocity_nonEqExpol
  public :: velocity_nonEqExpol_curved
  public :: pressure_nonEqExpol

contains

! **************************************************************************** !
  !> Element-wise Dirichlet velocity non-equilibrium boundary condition for
  !! straight boundary to update all directions.
  !! For straight wall, values are extrapolated along
  !! boundary normal instead of along the link and qVal =0.0 for straight wall.
  !!
  !! Notation: b (fictious boundary)   w  (physical bolundary or surface)
  !!           f (local element) ff (first neighbor fluid element)
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'inlet',
  !!    kind = 'velocity_nonEqExpol',
  !!    velocity = 'inlet_vel,
  !!  }
  !!}
  !!variable = {
  !!  name = 'inlet_vel',
  !!  ncomponents = 3,
  !!  vartype = 'st_fun',
  !!  st_fun = {0.06, 0.0, 0.0}
  !!}
  !!```
  !! This is described in the paper:
  !! Non-equilibrium extrapolation method
  !! for velocity and pressure boundary conditions in the lattice Boltzmann
  !! method." Chinese Physics.
  !!
  !! More informations concerning the nonEqExpol can be found in the
  !! corresponding subroutine.
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine velocity_nonEqExpol( me, state, bcBuffer, globBC, levelDesc,      &
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
    real(kind=rk) :: feq_b(layout%fStencil%QQ), vel_b(3)
    ! variables for surface (density and link-wise velocity)
    real(kind=rk) :: vel_w(globBC%nElems(iLevel)*3)
    ! variables for overnext fluid element
    real(kind=rk) :: feq_ff(layout%fStencil%QQ), rho_ff, vel_ff(3)
    real(kind=rk) :: inv_vel
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    integer :: bcVel_pos, QQ, iElem, elemPos
    ! -------------------------------------------------------------------- !
    QQ = layout%fstencil%QQ
    inv_vel = 1._rk / physics%fac(iLevel)%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! get velocity_phy on boundary
    call varSys%method%val(bcVel_pos)%get_valOfIndex(                    &
      &    varSys  = varSys,                                             &
      &    time    = sim_time,                                           &
      &    iLevel  = iLevel,                                             &
      &    idx     = me%bc_states%velocity%pntIndex                      &
      &                                   %indexLvl(iLevel)              &
      &                                   %val(1:globBC%nElems(iLevel)), &
      &    nVals   = globBC%nElems(iLevel),                              &
      &    res     = vel_w                                               )

    ! convert physical velocity into LB velocity
    vel_w = vel_w * inv_vel

    ! Update all directions of boundary elements
    do iElem = 1, globBC%nElems(iLevel)
      ! determine needed quantities of the overnext fluid neighbor element
      ! x_ff
      pdf_ff = me%neigh(iLevel)%neighBufferPre_nNext(1, (iElem-1) * QQ + 1 &
        &                                            : (iElem-1) * QQ + QQ )

      ! calculate density
      rho_ff = sum(pdf_ff)
      ! calulate velocity
      vel_ff = layout%quantities%vel_from_pdf_ptr( pdf = pdf_ff, dens = rho_ff )

      ! determine needed quantities of fictious boundary (_b)
      ! compute velocity on boundary (eq.18)
      ! qVal = 0.0, boundary node is one lattice away from fluid node!
      vel_b = vel_w( (iElem-1)*3+1:iElem*3)

      ! position of this boundary element in state array
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! derive feq from velocity and density for all direction
      fEq_ff = layout%quantities%pdfEq_ptr( rho = rho_ff,            &
        &                                   vel = vel_ff,            &
        &                                   QQ = QQ                  )
      ! extrapolate velocity on boundary
      ! derive feq from velocity and density per direction
      fEq_b = layout%quantities%pdfEq_ptr( rho = rho_ff,             &
        &                                   vel = vel_b,             &
        &                                   QQ = QQ                  )
      do iDir = 1, layout%fStencil%QQ
        ! Update pre-collision PDF
        state(  neigh(( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
          & = fEq_b(iDir) + (pdf_ff(iDir) - fEq_ff(iDir))
      end do

    end do !iElem

  end subroutine velocity_nonEqExpol
! **************************************************************************** !

! **************************************************************************** !
  !> Linkwise Dirichlet velocity non-equilibrium boundary condition for
  !! curved using the subroutine "mus_set_nonEqExpol".
  !! For curved wall, values are extrapolated along
  !! element normal
  !!
  !! Notation: b (fictious boundary)   w  (physical bolundary or surface)
  !!           f (local element)      ff  (overnext fluidneighbor)
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'inlet',
  !!    kind = 'velocity_nonEqExpol',
  !!    velocity = 'inlet_vel,
  !!    curved = true
  !!  }
  !!}
  !!variable = {
  !!  name = 'inlet_vel',
  !!  ncomponents = 3,
  !!  vartype = 'st_fun',
  !!  st_fun = {0.06, 0.0, 0.0}
  !!}
  !!```
  !! This is described in the paper:
  !! Guo, Z.; Zheng, C. & Shi B. (2002). An extrapolation method for boundary
  !! conditions in lattice Boltzmann method.
  !! Physics of Fluids 14, 2007 (2002); https://doi.org/10.1063/1.1471914
  !!
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine velocity_nonEqExpol_curved( me, state, bcBuffer, globBC,         &
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
    real(kind=rk) :: feq_b, fneq_b, rho_b, vel_b(3)
    ! variables for surface (density and link-wise velocity)
    real(kind=rk) :: vel_w(me%links(iLevel)%nVals*3)
    ! variables for fluid element
    real(kind=rk) :: feq_f, fneq_f, rho_f, vel_f(3)
    ! variables for overnext fluid element
    real(kind=rk) :: feq_ff, fneq_ff, rho_ff, vel_ff(3)
    real(kind=rk) :: inv_vel
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdfTmp(layout%fStencil%QQ)
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    integer :: bcVel_pos, iLink, QQ, offset
    ! relaxation parameter
    !real(kind=rk) :: omega
    ! -------------------------------------------------------------------- !
    QQ = layout%fstencil%QQ
    inv_vel = 1._rk / physics%fac(iLevel)%vel

    ! position of boundary velocity in varSys
    bcVel_pos = me%bc_states%velocity%varPos
    ! get velocity_phy on boundary (surface)
    call varSys%method%val(bcVel_pos)%get_valOfIndex(                     &
      &    varSys  = varSys,                                              &
      &    time    = sim_time,                                            &
      &    iLevel  = iLevel,                                              &
      &    idx     = me%bc_states%velocity%pntIndex                       &
      &                                   %indexLvl(ilevel)               &
      &                                   %val(1:me%links(iLevel)%nVals), &
      &    nVals   = me%links(iLevel)%nVals,                              &
      &    res     = vel_w                                                )

    ! convert physical velocity into LB velocity
    vel_w = vel_w * inv_vel

    ! loop over links
    do iLink = 1, me%links(iLevel)%nVals

      ! load pre-calculated coefficients from subroutine
      c_w      = me%nonEqExpol(iLevel)%     c_w(iLink)
      c_f      = me%nonEqExpol(iLevel)%     c_f(iLink)
      c_ff     = me%nonEqExpol(iLevel)%    c_ff(iLink)
      c_neq_f  = me%nonEqExpol(iLevel)% c_neq_f(iLink)
      c_neq_ff = me%nonEqExpol(iLevel)%c_neq_ff(iLink)

      ! set link-wise direction
      iDir = me%nonEqExpol(iLevel)%iDir(iLink)

      ! determine needed quantities of the current element (_f)
      ! calculate density
      posInBuffer  = me%nonEqExpol(iLevel)%posInBuffer(iLink)
      pdfTmp(1:QQ) = bcBuffer( (posInBuffer - 1) * nScalars + varPos(1)     &
        &                      : ( posInBuffer - 1) * nScalars + varPos(QQ) )
      rho_f = sum(pdfTmp)
      ! calulate velocity (facilitate it by using a pure function)
      vel_f = layout%quantities%vel_from_pdf_ptr( pdf = pdfTmp, dens = rho_f )
      ! derive feq from velocity and density per direction
      feq_f = layout%quantities%pdfEq_iDir_ptr( rho = rho_f,                               &
        &                                       vel = vel_f,                               &
        &                                       iDir = iDir,                               &
        &                                       cxDirRK = layout%fStencil%cxDirRK(:,iDir), &
        &                                       weight = layout%weight(iDir)               )

      ! determine needed quantities of the overnext fluid neighbor element (_ff)
      ! calculate density
      posInNeighBuf = me%nonEqExpol(iLevel)%posInNeighBuf(iLink)
      pdf_ff = me%neigh(iLevel)%neighBufferPost(                            &
        &        1,                                                         &
        &        (posInNeighBuf - 1) * QQ + 1:(posInNeighBuf - 1) * QQ + QQ )
      rho_ff = sum(pdf_ff)
      ! calulate velocity (facilitate it by using a pure function)
      vel_ff = layout%quantities%vel_from_pdf_ptr( pdf = pdf_ff, dens = rho_ff )
      ! derive feq from velocity and density per direction
      feq_ff = layout%quantities%pdfEq_iDir_ptr( rho = rho_ff,                              &
        &                                        vel = vel_ff,                              &
        &                                        iDir = iDir,                               &
        &                                        cxDirRK = layout%fStencil%cxDirRK(:,iDir), &
        &                                        weight = layout%weight(iDir)               )

      ! need offset for velocity
      offset = (iLink-1)*3

      ! determine needed quantities of fictious boundary (_b)
      ! compute velocity on boundary (eq.18)
      vel_b = (c_w * vel_w(offset+1:offset+3)) + (c_f * vel_f) + (c_ff * vel_ff)
      ! extrapolate density on boundary
      !rho_b = div4_3*rho_f - div1_3*rho_ff
      rho_b = rho_f
      ! derive feq from velocity and density per direction
      feq_b = layout%quantities%pdfEq_iDir_ptr( rho = rho_b,                               &
        &                                       vel = vel_b,                               &
        &                                       iDir = iDir,                               &
        &                                       cxDirRK = layout%fStencil%cxDirRK(:,iDir), &
        &                                       weight = layout%weight(iDir)               )

      ! use pdf_b = fEq_b + fneq_b (eq.16) to determine
      ! non-equilibrium components (eq.19)
      fneq_f  = pdfTmp(iDir) - feq_f
      fneq_ff = pdf_ff(iDir) - feq_ff
      fneq_b  = c_nEq_f * fneq_f + c_nEq_ff * fneq_ff

      ! write into state
      ! fneq is computed from post-collision pdf so no need to multiply
      ! (1-omega) as in the literature
      state( me%links(iLevel)%val(iLink) ) = feq_b + fneq_b

    end do !iLink

  end subroutine velocity_nonEqExpol_curved
! **************************************************************************** !


! **************************************************************************** !
  !> Element-wise Dirichlet pressure non-equilibrium boundary condition for
  !! straight boundary and updates all directions.
  !! For straight wall, values are extrapolated along
  !! boundary normal instead of along the link and qVal=0.
  !!
  !! Notation: b (fictious boundary)   w  (physical bolundary or surface)
  !!           f (local element)      ff  (overnext fluidneighbor)
  !!
  !! Usage
  !! -----
  !!```lua
  !!boundary_condition = {
  !!  { label = 'outlet',
  !!    kind = 'pressure_nonEqExpol',
  !!    pressure = 'press,
  !!    curved = true
  !!  }
  !!}
  !!variable = {
  !!  name = 'press',
  !!  ncomponents = 1,
  !!  vartype = 'st_fun',
  !!  st_fun = 0.0
  !!}
  !!```
  !! It is based on the paper:
  !! Guo, Z.; Zheng, C. & Shi B. (2002). An extrapolation method for boundary
  !! conditions in lattice Boltzmann method.
  !! Physics of Fluids 14, 2007 (2002); https://doi.org/10.1063/1.1471914
  !!
  !! More informations concerning the nonEqExpol can be found in the
  !! corresponding subroutine.
  !!
  !!
  !! This subroutine's interface must match the abstract interface definition
  !! [[boundaryRoutine]] in bc/[[mus_bc_header_module]].f90 in order to be
  !! callable via [[boundary_type:fnct]] function pointer.
  subroutine pressure_nonEqExpol( me, state, bcBuffer, globBC, levelDesc,      &
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
    real(kind=rk) :: feq_b(layout%fStencil%QQ), rho_b
    ! variables for surface (density and link-wise velocity)
    real(kind=rk) :: rho_w(globBC%nElems(iLevel))
    ! variables for overnext fluid element
    real(kind=rk) :: feq_ff(layout%fStencil%QQ), rho_ff, rhoInv, vel_ff(3)
    ! temporary local pdf values and the ones of overnext fluid
    real(kind=rk) :: pdf_ff(layout%fStencil%QQ)
    integer :: bcPress_pos, QQ, iElem, elemPos
    ! relaxation parameter
    ! real(kind=rk) :: omega
    real(kind=rk) :: inv_rho_phy
    ! -------------------------------------------------------------------- !
    QQ = layout%fstencil%QQ
    inv_rho_phy = 1.0_rk / physics%fac(iLevel)%press * cs2inv

    ! position of boundary velocity in varSys
    bcPress_pos = me%bc_states%pressure%varPos
    ! get pressure_phy on boundary (surface) and store it to rho_w
    call varSys%method%val(bcPress_pos)%get_valOfIndex(                  &
      &    varSys  = varSys,                                             &
      &    time    = sim_time,                                           &
      &    iLevel  = iLevel,                                             &
      &    idx     = me%bc_states%pressure%pntIndex                      &
      &                                   %indexLvl(ilevel)              &
      &                                   %val(1:globBC%nElems(iLevel)), &
      &    nVals   = globBC%nElems(iLevel),                              &
      &    res     = rho_w                                               )

    ! convert physical physical pressure into LB density
    rho_w = rho_w * inv_rho_phy

    ! Update all directions of boundary elements
    do iElem = 1, globBC%nElems(iLevel)
      ! determine needed quantities of the overnext fluid neighbor element
      ! x_ff
      pdf_ff = me%neigh(iLevel)%neighBufferPre_nNext(1, (iElem-1) * QQ + 1 &
        &                                            : (iElem-1) * QQ + QQ )

      ! calculate density
      rho_ff = sum(pdf_ff)
      rhoInv = 1.0_rk / rho_ff
      ! calulate velocity
      vel_ff = layout%quantities%vel_from_pdf_ptr( pdf = pdf_ff, dens = rho_ff )

      ! density on boundary, qVal=0: boundary and fluid node overlaps
      rho_b = rho_w(iElem)

      ! position of this boundary element in state array
      elemPos = globBC%elemLvl(iLevel)%elem%val( iElem )
      ! derive feq from velocity and density for all direction
      fEq_ff = layout%quantities%pdfEq_ptr( rho = rho_ff,            &
        &                                   vel = vel_ff,            &
        &                                   QQ = QQ                  )
      ! extrapolate velocity on boundary
      ! derive feq from velocity and density per direction
      fEq_b = layout%quantities%pdfEq_ptr( rho = rho_b,             &
        &                                  vel = vel_ff,            &
        &                                  QQ = QQ                  )
      do iDir = 1, layout%fStencil%QQ
        state(  neigh(( idir-1)*nsize+ elempos)+( ifield-1)* qq+ nscalars*0) &
          & = fEq_b(iDir) + (pdf_ff(iDir) - fEq_ff(iDir))
      end do

    end do !iElem

  end subroutine pressure_nonEqExpol
! **************************************************************************** !

end module mus_bc_fluid_nonEqExpol_module
