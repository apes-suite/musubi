! Copyright (c) 2012-2022 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2012, 2014-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013, 2015, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2015-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2018, 2020 Jana Gericke <jana.gericke@uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
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
!  See copyright notice in the COPYRIGHT file.
! ****************************************************************************** !
!> author: Manuel Hasert
!! This module contains general boundary routines
!!
!! like initializing boundary, setting boundary at every time step for each
!! field and other general boundary routines
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
module mus_bc_general_module

  ! include treelm modules
  use env_module,              only: rk, long_k, labelLen, pathLen, newUnit
  use tem_param_module,        only: cs, qOffset_inChar, q000
  use treelmesh_module,        only: treelmesh_type
  use tem_aux_module,          only: tem_abort
  use tem_float_module,        only: operator(.feq.), operator(.fne.)
  use tem_varSys_module,       only: tem_varSys_type
  use tem_varMap_module,       only: tem_varMap_type
  use tem_bc_prop_module,      only: tem_bc_prop_type
  use tem_topology_module,     only: tem_levelOF
  use tem_debug_module,        only: dbgUnit
  use tem_logging_module,      only: logUnit
  use tem_math_module,         only: invert_matrix
  use tem_geometry_module,     only: tem_BaryOfId,tem_ElemSizeLevel
  use tem_stencil_module,      only: tem_stencilHeader_type
  use tem_bc_module,           only: tem_bc_state_type
  use tem_grow_array_module,   only: init, append, destroy
  use tem_construction_module, only: tem_levelDesc_type
  use tem_timer_module,        only: tem_startTimer, tem_stopTimer

  ! include musubi modules
  use mus_field_module,                 only: mus_field_type
  use mus_field_prop_module,            only: mus_field_prop_type
  use mus_scheme_layout_module,         only: mus_scheme_layout_type
  use mus_scheme_header_module,         only: mus_scheme_header_type
  use mus_scheme_type_module,           only: array2D_type
  use mus_pdf_module,                   only: pdf_data_type
  use mus_bc_header_module,             only: boundary_type,                   &
    &                                         glob_boundary_type,              &
    &                                         mus_set_bcLinks, mus_set_bouzidi,&
    &                                         mus_alloc_bouzidi,               &
    &                                         mus_set_inletUbb,                &
    &                                         mus_set_inletBfl,                &
    &                                         mus_set_nonEqExpol,              &
    &                                         mus_set_outletExpol
  use mus_bc_fluid_module,              only: velocity_eq, vel_neq,            &
    &                                         outlet_nrbc, outlet_nrbc_eq,     &
    &                                         outlet_nrbc_incomp,              &
    &                                         inlet_nrbc, inlet_nrbc_incomp,   &
    &                                         outlet_dnt, pressure_eq,         &
    &                                         pressure_expol, press_neq,       &
    &                                         pressure_antiBounceBack,         &
    &                                         outlet_zero_prsgrd,              &
    &                                         mfr_bounceback, mfr_eq,          &
    &                                         velocity_bounceback,             &
    &                                         velocity_bounceback_incomp,      &
    &                                         velocity_bfl, bc_pdf,            &
    &                                         velocity_bfl_incomp
  use mus_bc_fluid_nonEqExpol_module,   only: velocity_nonEqExpol,             &
    &                                         pressure_nonEqExpol,             &
    &                                         velocity_nonEqExpol_curved
  use mus_bc_fluid_wall_module,         only: slip_wall,                &
    &                                         wall_libb, do_nothing,    &
    &                                         spc_slip_wall
  use mus_bc_fluid_turbulent_module,    only: turbulent_wall,                  &
    &                                  turbulent_wall_libb,                    &
    &                                  turbulent_wall_eq,                      &
    &                                  turbulent_wall_eq_curved,               &
    &                                  turbulent_wall_noneq_expol,             &
    &                                  turbulent_wall_noneq_expol_curved
  use mus_bc_fluid_experimental_module, only: pressure_expol_slow
  use mus_bc_fluid_module,              only: pressure_momentsbased,           &
    &                                         moments_wall,                    &
    &                                         velocity_momentsbased,           &
    &                                         velocity_momentsbased_incomp,    &
    &                                         pressure_momentsbased_incomp
  use mus_bc_fluid_experimental_module, only: moments_inflow, moments_outflow, &
    &                                         spc_moments_outflow,             &
    &                                         spc_bb_wall, spc_bb_vel_test
  use mus_bc_passiveScalar_module,      only: outlet_pasScal, inlet_pasScal,   &
    &                                         pressure_antiBounceBack_pasScal
  use mus_bc_species_module,            only: spc_outlet_zero_prsgrd,          &
    &                                         spc_moleFrac, spc_moleFlux,      &
    &                                         spc_moleFrac_wtdf,               &
    &                                         spc_moleFrac_eq, spc_inlet,      &
    &                                         spc_moleDens_eq,                 &
    &                                         spc_mole_fraction_noneq_expol,   &
    &                                         spc_velocity_noneq_expol,        &
    &                                         spc_moleDiff_Flux, spc_inlet_eq, &
    &                                         spc_outlet_eq, spc_outlet_vel,   &
    &                                         spc_outlet_expol,                &
    &                                         spc_outflow,                     &
    &                                         spc_inflow,                      &
    &                                         spc_solvent_inflow,              &
    &                                         spc_solvent_outflow,             &
    &                                         spc_moleFlux_eq, spc_vel_bb,     &
    &                                         spc_blackbox_mem_ion,            &
    &                                         spc_blackbox_mem_solvent,        &
    &                                         spc_moments_moleFrac,            &
    &                                         spc_moments_moleFlux,            &
    &                                         spc_moments_wall,                &
    &                                         spc_moments_vel
  use mus_bc_poisson_module,            only: potential_nonEqExpol,        &
    &                                         potential_nonEqExpol_curved, &
    &                                         potential_neumann,           &
    &                                         potential_neumann_curved
  use mus_bc_nernstPlanck_module,       only: moleDens_nonEqExpol,        &
    &                                         moleDens_nonEqExpol_curved, &
    &                                         moleDens_neumann,           &
    &                                         moleDens_neumann_curved

  use mus_relaxationParam_module,       only: mus_viscosity_type
  use mus_auxField_module,              only: mus_auxFieldVar_type
  use mus_derVarPos_module,             only: mus_derVarPos_type
  use mus_param_module,                 only: mus_param_type
  use mus_physics_module,               only: mus_physics_type
  use mus_mixture_module,               only: mus_mixture_type
  use mus_timer_module,                 only: mus_timerHandles

  implicit none

  private

  public :: set_boundary
  public :: mus_init_boundary
  public :: mus_get_points_fromBC

contains


! ****************************************************************************** !
  !> Call the functions associated with each boundary condition
  !!
  !! Loop over each field and Run over all the boundary conditions for
  !! current iLevel and call the function pointer.
  !! The function pointer was before assigned in init_boundary
  !! This routine is being called from the
  !! [[mus_control_module:do_recursive_multilevel]] "main control routine"
  !! before the compute (advection_relaxation) kernel call
  subroutine set_boundary( field, pdf, levelDesc, tree, iLevel, nBCs, params, &
    &                      layout, varSys, derVarPos, globBC, mixture,        &
    &                      physics, state )
    ! ---------------------------------------------------------------------------
    !> fluid parameters and properties
    type( mus_field_type ), intent(inout) :: field(:)
    !> contains global state vector
    type( pdf_data_type ), intent(inout) :: pdf
    !> state arrays fo current iLevel both now and next
    real(kind=rk), intent(inout) :: state(:,:)
    !> global type contains iLevel descriptor
    type( tem_levelDesc_type ), intent(in) :: levelDesc
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> the iLevel on which this boundary was invoked
    integer, intent(in) :: iLevel
    !> number of BC
    integer, intent(in) :: nBCs
    !> global parameters
    type(mus_param_type),intent(in) :: params
    !> stencil layout information
    type( mus_scheme_layout_type ), intent(in) ::layout
    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys
    !> position of derived quantities in varsys
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    !> scheme global boundary type
    type( glob_boundary_type ), intent(in) :: globBC(:)
    !> scheme global boundary type
    type( mus_physics_type ), intent(in) :: physics
    !> mixture info
    type(mus_mixture_type), intent(in) :: mixture
    ! --------------------------------------------------------------------------
    integer :: iField, nFields, iBnd
    ! --------------------------------------------------------------------------

    nFields = size( field )

    if ( nBCs > 0 ) then
      call tem_startTimer( timerHandle = mus_timerHandles%bcBuffer(iLevel) )
      call fill_bcBuffer(                                 &
        &    currState = state( :, pdf%nNext ),           &
        &    neigh     = pdf%neigh,                       &
        &    bcBuffer  = pdf%bcBuffer,                    &
        &    nSize     = pdf%nSize,                       &
        &    nElems_bc = levelDesc%bc_elemBuffer%nVals,   &
        &    posInTotal= levelDesc%bc_elemBuffer%val,     &
        &    nFields   = nFields,                         &
        &    QQ        = layout%fStencil%QQ,              &
        &    varSys    = varSys                           )
      call tem_stopTimer( timerHandle = mus_timerHandles%bcBuffer(iLevel) )

      ! @todo: As neigh buffer is inside
      !        field( iField )%bc( iBnd )%neigh( iLevel ),
      !        we may fill it in the following fileds and BCs loop.
      call fill_neighBuffer(                       &
        &    prevstate  = state( :, pdf%nNow ),    &
        &    currstate  = state( :, pdf%nNext ),   &
        &    neigh  = pdf%neigh,                   &
        &    globBC = globBC,                      &
        &    nBCs   = nBCs,                        &
        &    field  = field,                       &
        &    varSys = varSys,                      &
        &    QQ     = layout%fStencil%QQ,          &
        &    nSize  = pdf%nSize,                   &
        &    iLevel = iLevel                       )

      ! loop over fields
      do iField = 1, nFields

        ! Treat all boundary conditions
        do iBnd = 1, nBCs

          ! write(dbgUnit(10),*) 'Do boundary condition: '// trim( globBC( iBnd )%label )

          call tem_startTimer( timerHandle =  mus_timerHandles%setBnd(iBnd) )
          call field(iField)%bc( iBnd )%fnct(                   &
            &      state       = state( :, pdf%nNext ),         &
            &      nSize       = pdf%nSize,                     &
            ! we want to write into the current time step
            ! the field pointers have been inverted, so the newest
            ! information is stored in nNow
            &      bcBuffer    = pdf%bcBuffer,                  &
            &      globBC      = globBC( iBnd ),                &
            &      levelDesc   = levelDesc,                     &
            &      tree        = tree,                          &
            &      iLevel      = iLevel,                        &
            &      sim_time    = params%general%simControl%now, &
            ! &      params      = params,                        &
            &      physics     = physics,                       &
            &      mixture     = mixture,                       &
            &      iField      = iField,                        &
            &      neigh       = pdf%neigh,                     &
            &      layout      = layout,                        &
            &      fieldProp   = field(iField)%fieldProp,       &
            &      varSys      = varSys,                        &
            &      derVarPos   = derVarPos(iField),             &
            ! each field may has multiple state variables
            &      varPos      = varSys%method%val(iField)      &
            &                      %state_varPos,               &
            &      nScalars    = varSys%nScalars                )

          call tem_stopTimer( timerHandle =  mus_timerHandles%setBnd(iBnd) )
        end do  ! iBnd
      end do  ! iField
    end if
  end subroutine set_boundary
! ****************************************************************************** !


! ****************************************************************************** !
  ! > This routine call init_boundary to initialize boundary for each field
  !!
  subroutine mus_init_boundary( field, pdf, tree, levelDesc, layout,     &
    &                           schemeHeader, varSys, derVarPos, globBC, &
    &                           bc_prop, state, auxField )
    !---------------------------------------------------------------------------
    !> fluid parameters and properties
    type( mus_field_type ), intent(inout) :: field(:)
    !> global treelm mesh
    type( treelmesh_type ), intent(in) :: tree
    !> contains global state vector
    type( pdf_data_type ), intent(inout) :: pdf(tree%global%minLevel &
      &                                         :tree%global%maxLevel)
    !> state array
    type( array2D_type ), intent(inout) :: state(tree%global%minLevel &
      &                                          :tree%global%maxLevel)
    !> AuxField array
    type( mus_auxFieldVar_type), intent(in) :: auxField(tree%global%minLevel &
      &                                                 :tree%global%maxLevel)
    !> scheme layout type
    type( mus_scheme_layout_type ), intent(in) ::layout
    !> scheme header info
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    !> Level Descriptor
    type( tem_leveldesc_type ), intent(in) :: levelDesc(tree%global%minLevel:tree%global%maxLevel)
    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys
    !> position of derived quantities in varsys
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    !> scheme global boundary type
    type( glob_boundary_type ), intent(inout) :: globBC(:)
    !> boundary property type
    type( tem_bc_prop_type ) :: bc_prop
    ! ---------------------------------------------------------------------------
    integer :: iField, iBC, iLevel, minLevel, maxLevel
    ! ---------------------------------------------------------------------------

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel

    if ( bc_prop%nBCtypes > 0 ) then
      !> Check prerequisite for multi-species boundary conditions
      if (size(field)>1) then
        call check_BCs_preRequisite_MS( field, bc_prop%nBCtypes )
      end if

      do iField = 1, size(field)
        ! Set up links which are needed to be updated
        write(logUnit(7), *) 'Counting actual number of links that need '//&
          &                  'to be update:'
        do iBC = 1, bc_prop%nBCtypes
          if ( .not. globBC( iBC )%isWall ) then
            do iLevel = minLevel, maxLevel
              if ( globBC(iBC)%nElems(iLevel) >= 0 ) then
                write(logUnit(7), "(2(A,I0))") &
                  &  trim(globBC(iBC)%label)//', level: ', iLevel, &
                  &  ', local nElems: ', globBC(iBC)%nElems(iLevel)
                call mus_set_bcLinks( iField   = iField,                     &
                  &                   QQ       = layout%fStencil%QQ,         &
                  &                   QQN      = layout%fStencil%QQN,        &
                  &                   nScalars = varSys%nScalars,            &
                  &                   nElems   = globBC(iBC)%nElems(iLevel), &
                  &                   nSize    = pdf(iLevel)%nSize,          &
                  &                   elemLvl  = globBC(iBC)%elemLvl(iLevel),&
                  &                   neigh    = pdf(iLevel)%neigh,          &
                  &                   links    = field(iField)%bc(iBC)       &
                  &                                           %links(iLevel) )
                write(logUnit(7), "(A,I0)") &
                  &     ' nLinks: ', field(iField)%bc(iBC)%links(iLevel)%nVals
              end if
            end do ! iLevel
          end if
        end do ! iBC

        call init_boundary_single(                                     &
          &    bc           = field(iField)%bc,                        &
          &    pdf          = pdf(minLevel:maxLevel),                  &
          &    state        = state(minLevel:maxLevel),                &
          &    auxField     = auxField(minLevel:maxLevel),             &
          &    tree         = tree,                                    &
          &    leveldesc    = levelDesc,                               &
          &    layout       = layout,                                  &
          &    schemeHeader = schemeHeader,                            &
          &    varPos       = varSys%method%val(iField)%state_varPos,  &
          &    varSys       = varSys,                                  &
          &    dervarPos    = derVarPos(iField),                       &
          &    globBC       = globBC,                                  &
          &    bc_prop      = bc_prop,                                 &
          &    fieldProp    = field(iField)%fieldProp                  )
      end do ! iField
    end if

  end subroutine mus_init_boundary
! ****************************************************************************** !


! ****************************************************************************** !
  !> This subroutine sets the right boundary conditions for the different
  !! boundaries.
  !!
  subroutine init_boundary_single(bc, pdf, tree, levelDesc, layout,          &
    &                             schemeHeader, varPos, varSys, derVarPos,   &
    &                             globBC,bc_prop, state, auxField, fieldProp )
    ! ---------------------------------------------------------------------------
    !> global array boundary type
    type( boundary_type ) :: bc(:)
    !> global treelm mesh
    type( treelmesh_type ), intent(in) ::tree
    !> contains global state vector
    type( pdf_data_type ), intent(in) :: pdf( tree%global%minLevel &
      &                                     : tree%global%maxLevel )
    !> contains global state vector
    type( array2D_type ), intent(in) :: state( tree%global%minLevel &
      &                                      : tree%global%maxLevel )
    !> AuxField array
    type( mus_auxFieldVar_type), intent(in) :: auxField(tree%global%minLevel &
      &                                                 :tree%global%maxLevel)
    !> scheme layout type
    type( mus_scheme_layout_type ), intent(in) ::layout
    !> scheme header info
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    !> global pdf type
    type( tem_leveldesc_type ), intent(in) :: levelDesc( tree%global%minLevel &
      &                                                : tree%global%maxLevel )
    !> varPos of current field variable
    integer, intent(in) :: varPos(:)
    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys
    !> position of derived quantities in varsys
    type( mus_derVarPos_type ), intent(in) :: derVarPos
    !> scheme global boundary type
    type( glob_boundary_type ), intent(inout) :: globBC(:)
    !> boundary property type
    type( tem_bc_prop_type ), intent(in) :: bc_prop
    !> fluid parameters and properties
    type( mus_field_prop_type ), intent(in) :: fieldProp
    ! ---------------------------------------------------------------------------
    integer :: iBnd, iLevel, minLevel, maxLevel, nBCs
    logical :: isWall
    ! ---------------------------------------------------------------------------

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel
    nBCs     = bc_prop%nBCtypes

    !> The boundary conditions were read in the order of config.lua
    !! the order in the boundary condition description file however
    !! might be different.
    !! We have to match the labels.
    write(logUnit(1),'(A,I0,A)') 'Initialize ', nBCs, ' boundary conditions'

    do iBnd = 1, nBCs
      if( trim(bc( iBnd )%label ) /= trim( globBC( iBnd )%label )) then
        write(logUnit(1),*) 'Error: The boundary with label from the mesh:'    &
          &                 //trim( globBC(iBnd)%label )
        write(logUnit(1),*) 'does not match the one specified in the lua file:'&
          &                 //trim(bc(iBnd)%label)
        write(logUnit(1),*) 'Stopping ...'
        call tem_abort()
      end if

      write(logUnit(1),*) 'Initializing BC: '//trim(bc(iBnd)%label)

      isWall = .false.
      select case (trim(bc(iBnd)%BC_kind))
      case('wall', 'symmetry')
        isWall = .true.
        bc( iBnd )%fnct => do_nothing
      case('wall_libb', 'velocity_bounceback', 'velocity_bfl')
        bc( iBnd )%evalBcVar_link = .true.

        ! Here we have to allocate and set the q-values if it is not
        ! provided by seeder and set the qVal from musubi.lua
        if( .not. globBC(iBnd)%hasQVal .and.        &
          & .not. globBC(iBnd)%qValInitialized ) then
          call init_qVals( refQval   = bc( iBnd )%qVal, &
            &              minLevel  = minLevel,        &
            &              maxLevel  = maxLevel,        &
            &              layout    = layout,          &
            &              globBC    = globBC(iBnd)     )
        end if

        select case (trim(bc( iBnd )%BC_kind))
        case('wall_libb')
          isWall = .true.
          bc( iBnd )%fnct => wall_libb
          ! set link-wise data
          call mus_alloc_bouzidi( me       = bc(iBnd)%bouzidi,       &
            &    nVals    = bc(iBnd)%links(minLevel:maxLevel)%nVals, &
            &    minLevel = minLevel,                                &
            &    maxLevel = maxLevel                                 )
          do iLevel = minLevel, maxLevel
            call mus_set_bouzidi( iLevel   = iLevel,                    &
              &                   QQ       = layout%fStencil%QQ,        &
              &                   QQN      = layout%fStencil%QQN,       &
              &                   nScalars = varSys%nScalars,           &
              &                   globBC   = globBC(iBnd),              &
              &                   cxDirInv = layout%fStencil%cxDirInv,  &
              &                   varPos   = varPos,                    &
              &                   bouzidi  = bc( iBnd )%bouzidi(iLevel) )
          end do
        case('velocity_bounceback')
          select case(trim(schemeHeader%kind))
          case('fluid')
            bc( iBnd )%fnct => velocity_bounceback
          case('fluid_incompressible')
            bc( iBnd )%fnct => velocity_bounceback_incomp
          case default
            call tem_abort('Unknown scheme kind for velocity_bounceback')
          end select
          ! set link-wise data
          ! bc(iBnd)%links is built in mus_set_bcLinks
          call mus_set_inletUbb( me           = bc(iBnd)%inletUbbQVal, &
            &                    tree         = tree,                  &
            &                    stencil      = layout%fStencil,       &
            &                    nScalars     = varSys%nScalars,       &
            &                    globBC       = globBC(iBnd),          &
            &                    levelDesc    = levelDesc(minLevel:maxLevel), &
            &                    varPos       = varPos,                &
            &        nLinks = bc(iBnd)%links(minLevel:maxLevel)%nVals, &
            &                    minLevel     = minLevel,              &
            &                    maxLevel     = maxLevel               )
        case('velocity_bfl')
          select case(trim(schemeHeader%kind))
          case('fluid')
            bc( iBnd )%fnct => velocity_bfl
          case('fluid_incompressible')
            bc( iBnd )%fnct => velocity_bfl_incomp
          case default
            call tem_abort('Unknown scheme kind for velocity_bfl')
          end select
          ! set link-wise data
          call mus_set_inletBfl(                                       &
            &     me        = bc(iBnd)%inletBfl,                       &
            &     tree      = tree,                                    &
            &     stencil   = layout%fStencil,                         &
            &     nScalars  = varSys%nScalars,                         &
            &     globBC    = globBC(iBnd),                            &
            &     levelDesc = levelDesc(minLevel:maxLevel),            &
            &     varPos    = varPos,                                  &
            &     nLinks    = bc(iBnd)%links(minLevel:maxLevel)%nVals, &
            &     minLevel  = minLevel,                                &
            &     maxLevel  = maxLevel                                 )
        end select

      case('turbulent_wall', 'turbulent_wall_noneq_expol', 'turbulent_wall_eq')
        isWall = .true.
        ! Here we have to allocate and set the q-values if it is not
        ! provided by seeder and set the qVal from musubi.lua
        if (.not. globBC(iBnd)%hasQVal .and.       &
          & .not. globBC(iBnd)%qValInitialized) then
          call init_qVals( refQval   = bc( iBnd )%qVal, &
            &              minLevel  = minLevel,        &
            &              maxLevel  = maxLevel,        &
            &              layout    = layout,          &
            &              globBC    = globBC(iBnd)     )
        end if

        ! allocate turbulent viscosity on boundary elements,
        ! and friction velocity and normal distance to boundary on first neigbor
        ! of boundary elements.
        ! Initialize friction velocity from stream-wise
        ! velocity component computed on first neighbor and normal distance to
        ! to boundary.
        ! \todo KM: Calculate friction velocity from wall shear
        ! stress computed on the first neighbor using non-equilibrium pdf
        write(logUnit(10), "(A)") 'Initializing turbulent_wall ...'
        allocate(bc(iBnd)%turbwallFunc%dataOnLvl(minLevel: maxLevel))
        do iLevel = minLevel, maxLevel
          call mus_init_turb_wallFunc( bc        = bc(iBnd),                 &
            &                          globBC    = globBC(iBnd),             &
            &                          auxField  = auxField(iLevel)%val(:),  &
            &                          viscKine  = fieldProp%fluid%viscKine, &
            &                          derVarPos = derVarPos,                &
            &                          varSys    = varSys,                   &
            &                          stencil   = layout%fStencil,          &
            &                          iLevel    = iLevel                    )
        end do

        if (bc(iBnd)%curved) then
          ! set link-wise data for bouzidi wall linear interpolation
          call mus_alloc_bouzidi( me       = bc(iBnd)%bouzidi,       &
            &    nVals    = bc(iBnd)%links(minLevel:maxLevel)%nVals, &
            &    minLevel = minLevel,                                &
            &    maxLevel = maxLevel                                 )
          do iLevel = minLevel, maxLevel
            call mus_set_bouzidi( iLevel   = iLevel,                    &
              &                   QQ       = layout%fStencil%QQ,        &
              &                   QQN      = layout%fStencil%QQN,       &
              &                   nScalars = varSys%nScalars,           &
              &                   globBC   = globBC(iBnd),              &
              &                   cxDirInv = layout%fStencil%cxDirInv,  &
              &                   varPos   = varPos,                    &
              &                   bouzidi  = bc( iBnd )%bouzidi(iLevel) )
          end do
        end if


        select case (trim(bc(iBnd)%BC_kind))
        case ('turbulent_wall_noneq_expol')
          select case(trim(schemeHeader%kind))
          case('fluid', 'fluid_incompressible')
            if (bc(iBnd)%curved) then
              bc( iBnd )%fnct => turbulent_wall_noneq_expol_curved
            else
              bc( iBnd )%fnct => turbulent_wall_noneq_expol
            end if
          case default
            call tem_abort('Unknown scheme kind for '//trim(bc(iBnd)%BC_kind))
          end select
        case ('turbulent_wall_eq')
          select case(trim(schemeHeader%kind))
          case('fluid', 'fluid_incompressible')
            if (bc(iBnd)%curved) then
              bc( iBnd )%fnct => turbulent_wall_eq_curved
            else
              bc( iBnd )%fnct => turbulent_wall_eq
            end if
          case default
            call tem_abort('Unknown scheme kind for '//trim(bc(iBnd)%BC_kind))
          end select
        case ('turbulent_wall')
          select case(trim(schemeHeader%kind))
          case('fluid')
            if (bc(iBnd)%curved) then
              bc( iBnd )%fnct => turbulent_wall_libb
            else
              bc( iBnd )%fnct => turbulent_wall
            end if
          case default
            call tem_abort('Unknown scheme kind for '//trim(bc(iBnd)%BC_kind))
          end select
        end select

      case('velocity_eq')
        bc( iBnd )%fnct => velocity_eq
      case('mfr_bounceback')
        bc( iBnd )%fnct => mfr_bounceback
      case('mfr_eq')
        bc( iBnd )%fnct => mfr_eq
      case('vel_neq')
        bc( iBnd )%fnct => vel_neq
      case('bc_pdf')
        bc( iBnd )%fnct => bc_pdf
      case('outlet_nrbc','outlet_nrbc_eq', 'inlet_nrbc')
        select case (trim(bc( iBnd )%BC_kind))
          case( 'outlet_nrbc' )
            select case(trim(schemeHeader%kind))
            case ('fluid')
              bc( iBnd )%fnct => outlet_nrbc
            case('fluid_incompressible')
              bc( iBnd )%fnct => outlet_nrbc_incomp
            case default
              call tem_abort('Unknown scheme kind for outlet_nrbc')
            end select
          case( 'outlet_nrbc_eq' )
            bc( iBnd )%fnct => outlet_nrbc_eq
          case( 'inlet_nrbc' )
            select case(trim(schemeHeader%kind))
            case ('fluid')
              bc( iBnd )%fnct => inlet_nrbc
            case('fluid_incompressible')
              bc( iBnd )%fnct => inlet_nrbc_incomp
            case default
              call tem_abort('Unknown scheme kind for inlet_nrbc')
            end select
        end select

        bc( iBnd )%nrbc%cs_mod = cs * sqrt( bc( iBnd )%nrbc%kappa )

        write(logUnit(5),"(A)")       'Initializing NRBC with parameters:'
        do iLevel = minLevel, maxLevel
          call init_nrbc( bc        = bc(iBnd),                                 &
            &             state     = state(iLevel)%val(:,pdf(iLevel)%nNext),   &
            &             nSize     = pdf(iLevel)%nSize,                        &
            &             neigh     = pdf(iLevel)%neigh(:),                     &
            &             layout    = layout,                                   &
            &             level     = iLevel,                                   &
            &             nScalars  = varSys%nScalars,                          &
            &             varSys    = varSys,                                   &
            &             derVarPos = derVarPos,                                &
            &             elemPos   = globBC(iBnd)%elemLvl(iLevel)%elem%val(:), &
            &             nElems    = globBC(iBnd)%nElems(iLevel)               )
        end do
      case('outlet_dnt')
        bc( iBnd )%fnct => outlet_dnt
      case('pressure_expol', 'pressure_expol_slow')
        ! set link-wise data
        call mus_set_outletExpol(                                             &
          &            me          = bc(iBnd)%outletExpol,                    &
          &            stencil     = layout%fStencil,                         &
          &            globBC      = globBC(iBnd),                            &
          &            nLinks      = bc(iBnd)%links(minLevel:maxLevel)%nVals, &
          &            minLevel    = minLevel,                                &
          &            maxLevel    = maxLevel                                 )
        if( trim(bc(iBnd)%BC_kind) == 'pressure_expol' ) then
          bc( iBnd )%fnct => pressure_expol
        else
          bc( iBnd )%fnct => pressure_expol_slow
        end if
      case('press_neq')
        bc( iBnd )%fnct => press_neq
      case('pressure_eq')
        bc( iBnd )%fnct => pressure_eq
      case('pressure_antibounceback')
        select case(trim(schemeHeader%kind))
        case('fluid')
           bc( iBnd )%fnct => pressure_antiBounceBack
        case('passive_scalar')
           bc( iBnd )%fnct => pressure_antiBounceBack_pasScal
        case default
          call tem_abort('Unknown scheme kind for pressure_antiBounceBack')
        end select
      case('outlet_zero_prsgrd')
        bc( iBnd )%fnct => outlet_zero_prsgrd
      ! passive scalar boundary conditions
      case('flekkoy_inlet')
        bc( iBnd )%fnct => inlet_pasScal
      case('flekkoy_outlet')
        bc( iBnd )%fnct => outlet_pasScal
      ! multispecies boundary conditions
      case('spc_slip_wall')
        isWall = .true.
        bc( iBnd )%fnct => spc_slip_wall
      case('slip_wall')
        isWall = .true.
        bc( iBnd )%fnct => slip_wall
      case('spc_bb_wall')
        bc( iBnd )%fnct => spc_bb_wall
      case('spc_bb_vel_test')
        bc( iBnd )%fnct => spc_bb_vel_test
      case('spc_outlet_zero_prsgrd')
        bc( iBnd )%fnct => spc_outlet_zero_prsgrd
      case('spc_inlet_eq')
        bc( iBnd )%fnct => spc_inlet_eq
      case('spc_inlet')
        bc( iBnd )%fnct => spc_inlet
      case('spc_inflow')
        bc( iBnd )%fnct => spc_inflow
      case('spc_solvent_inflow')
        bc( iBnd )%fnct => spc_solvent_inflow
      case('spc_velocity_noneq_expol')
        bc( iBnd )%fnct => spc_velocity_noneq_expol
      case('spc_mole_fraction_noneq_expol')
        bc( iBnd )%fnct => spc_mole_fraction_noneq_expol
      case('spc_vel_bb')
        bc( iBnd )%fnct => spc_vel_bb
      case('spc_outlet_eq')
        bc( iBnd )%fnct => spc_outlet_eq
      case('spc_outlet_expol')
        bc( iBnd )%fnct => spc_outlet_expol
      case('spc_outflow')
        bc( iBnd )%fnct => spc_outflow
      case('spc_solvent_outflow')
        bc( iBnd )%fnct => spc_solvent_outflow
      case('spc_outlet_vel')
        bc( iBnd )%fnct => spc_outlet_vel
      case('spc_molefrac')
        bc( iBnd )%fnct => spc_moleFrac
      case('spc_molefrac_eq')
        bc( iBnd )%fnct => spc_moleFrac_eq
      case('spc_moledens_eq')
        bc( iBnd )%fnct => spc_moleDens_eq
      case('spc_molefrac_wtdf')
        bc( iBnd )%fnct => spc_moleFrac_wtdf
      case('spc_moleflux')
        bc( iBnd )%fnct => spc_moleFlux
      case('spc_moleflux_eq')
        bc( iBnd )%fnct => spc_moleFlux_eq
      case('spc_molediff_flux')
        bc( iBnd )%fnct => spc_moleDiff_Flux
      case( 'spc_blackbox_mem_ion')
        bc( iBnd )%fnct => spc_blackbox_mem_ion
      case( 'spc_blackbox_mem_solvent')
        bc( iBnd )%fnct => spc_blackbox_mem_solvent

      ! cases which use the non-equilibrium extrapolation (nonEqExpol)
      case('potential_noneq_expol', 'potential_neumann',    &
        &  'velocity_noneq_expol', 'pressure_noneq_expol' )
        bc( iBnd )%evalBcVar_link = .true.

        ! Here we have to allocate and set the q-values if it is not
        ! provided by seeder and set the qVal from musubi.lua
        if( .not. globBC(iBnd)%hasQVal .and.        &
          & .not. globBC(iBnd)%qValInitialized ) then
          call init_qVals( refQval   = bc(iBnd)%qVal, &
            &              minLevel  = minLevel,      &
            &              maxLevel  = maxLevel,      &
            &              layout    = layout,        &
            &              globBC    = globBC(iBnd)   )
        end if

        ! set link-wise data
        call mus_set_nonEqExpol(                                     &
          &     me        = bc(iBnd)%nonEqExpol,                     &
          &     curved    = bc(iBnd)%curved,                         &
          &     tree      = tree,                                    &
          &     stencil   = layout%fStencil,                         &
          &     nScalars  = varSys%nScalars,                         &
          &     globBC    = globBC(iBnd),                            &
          &     bc_neigh  = bc(iBnd)%neigh(minLevel:maxLevel),       &
          &     pdf       = pdf(minLevel:maxLevel),                  &
          &     levelDesc = levelDesc(minLevel:maxLevel),            &
          &     varPos    = varPos,                                  &
          &     nLinks    = bc(iBnd)%links(minLevel:maxLevel)%nVals, &
          &     minLevel  = minLevel,                                &
          &     maxLevel  = maxLevel                                 )

        select case ( trim(bc(iBnd)%BC_kind ) )
        case('potential_noneq_expol')
          if (bc(iBnd)%curved) then
            bc( iBnd )%fnct => potential_nonEqExpol_curved
          else
            bc( iBnd )%fnct => potential_nonEqExpol
          end if
        case('potential_neumann')
          if (bc(iBnd)%curved) then
            bc( iBnd )%fnct => potential_neumann_curved
          else
            bc( iBnd )%fnct => potential_neumann
          end if
        case('velocity_noneq_expol')
          if (trim(schemeHeader%relaxation(1:3)) == 'trt') then
            call tem_abort('velocity_noneq_expol is not supported for trt!')
          end if

          if (bc(iBnd)%curved) then
            bc( iBnd )%fnct => velocity_nonEqExpol_curved
          else
            bc( iBnd )%fnct => velocity_nonEqExpol
          end if

        case('pressure_noneq_expol')
          if (trim(schemeHeader%relaxation(1:3)) == 'trt') then
            call tem_abort('velocity_noneq_expol is not supported for trt!')
          end if
          bc( iBnd )%fnct => pressure_nonEqExpol

        end select

      ! cases which use the nernst_planck
      case('moledens_noneq_expol', 'moledens_neumann')
        bc( iBnd )%evalBcVar_link = .true.

        ! Here we have to allocate and set the q-values if it is not
        ! provided by seeder and set the qVal from musubi.lua
        if( .not. globBC(iBnd)%hasQVal .and.        &
          & .not. globBC(iBnd)%qValInitialized ) then
          call init_qVals( refQval   = bc( iBnd )%qVal, &
            &              minLevel  = minLevel,        &
            &              maxLevel  = maxLevel,        &
            &              layout    = layout,          &
            &              globBC    = globBC(iBnd)     )
        end if

        ! set link-wise data
        call mus_set_nonEqExpol(                                     &
          &     me        = bc(iBnd)%nonEqExpol,                     &
          &     curved    = bc(iBnd)%curved,                         &
          &     tree      = tree,                                    &
          &     stencil   = layout%fStencil,                         &
          &     nScalars  = varSys%nScalars,                         &
          &     globBC    = globBC(iBnd),                            &
          &     bc_neigh  = bc(iBnd)%neigh(minLevel:maxLevel),       &
          &     pdf       = pdf(minLevel:maxLevel),                  &
          &     levelDesc = levelDesc(minLevel:maxLevel),            &
          &     varPos    = varPos,                                  &
          &     nLinks    = bc(iBnd)%links(minLevel:maxLevel)%nVals, &
          &     minLevel  = minLevel,                                &
          &     maxLevel  = maxLevel                                 )

        select case ( trim(bc(iBnd)%BC_kind ) )
        case('moleDens_noneq_expol')
          if (bc(iBnd)%curved) then
            bc( iBnd )%fnct => moleDens_nonEqExpol_curved
          else
            bc( iBnd )%fnct => moleDens_nonEqExpol
          end if
        case('moledens_neumann')
          if (bc(iBnd)%curved) then
            bc( iBnd )%fnct => moleDens_neumann_curved
          else
            bc( iBnd )%fnct => moleDens_neumann
          end if
        end select

      case( 'pressure_momentsbased', 'moments_wall', 'velocity_momentsbased', &
        &   'moments_inflow', 'moments_outflow',                              &
        &   'spc_moments_molefrac', 'spc_moments_moleflux',                   &
        &   'spc_moments_wall', 'spc_moments_vel', 'spc_moments_outflow'      )
        select case (trim(bc( iBnd  )%BC_kind))
        case( 'pressure_momentsbased')
          select case (trim(schemeHeader%kind))
          case ('fluid')
            bc( iBnd )%fnct => pressure_momentsbased
          case ('fluid_incompressible')
            bc( iBnd )%fnct => pressure_momentsbased_incomp
          case default
            write(logUnit(1),*) 'Chosen boundary kind '                        &
              &                 //trim(bc(iBnd)%BC_kind)//' is not supported'  &
              &                //'scheme kind'//trim(schemeHeader%kind)
            call tem_abort()
          end select
        case( 'moments_wall')
          bc( iBnd )%fnct => moments_wall
        case( 'velocity_momentsbased')
          select case (trim(schemeHeader%kind))
          case ('fluid')
            bc( iBnd )%fnct => velocity_momentsbased
          case('fluid_incompressible')
            bc( iBnd )%fnct => velocity_momentsbased_incomp
          case default
            write(logUnit(1),*) 'Chosen boundary kind '                        &
              &                 //trim(bc(iBnd)%BC_kind)//' is not supported'  &
              &                //'scheme kind'//trim(schemeHeader%kind)
            call tem_abort()
          end select
        case( 'moments_inflow')
          bc( iBnd )%fnct => moments_inflow
        case( 'moments_outflow')
          bc( iBnd )%fnct => moments_outflow
        case( 'spc_moments_molefrac')
          bc( iBnd )%fnct => spc_moments_moleFrac
        case( 'spc_moments_moleflux')
          bc( iBnd )%fnct => spc_moments_moleFlux
        case( 'spc_moments_outflow')
          bc( iBnd )%fnct => spc_moments_outflow
        case( 'spc_moments_wall')
          bc( iBnd )%fnct => spc_moments_wall
        case( 'spc_moments_vel')
          bc( iBnd )%fnct => spc_moments_vel
        end select
        call init_momentsBC( bc        = bc(iBnd),     &
          &                  leveldesc = levelDesc,    &
          &                  layout    = layout,       &
          &                  globBC    = globBC(iBnd), &
          &                  minLevel  = minLevel,     &
          &                  maxLevel  = maxLevel      )
      case default
        write(logUnit(1),*)'BC type '//trim(bc( iBnd  )%bc_kind)// &
          & ' for label '//trim(bc( iBnd  )%label)// 'not found'
        call tem_abort()
      end select

      if ( .not. isWall ) then
        ! setup indices for each boundaries
        call mus_setupIndices_forBC( bc       = bc(iBnd),        &
          &                          globBC   = globBC(iBnd),    &
          &                          tree     = tree,            &
          &                          stencil  = layout%fStencil, &
          &                          levelDesc= levelDesc,       &
          &                          varSys   = varSys,          &
          &                          minLevel = minLevel,        &
          &                          maxLevel = maxLevel         )
      end if

    end do ! iBnd

  end subroutine init_boundary_single
! ****************************************************************************** !

  ! ************************************************************************** !
  !> Check prerequisite for some boundary conditions
  subroutine check_BCs_preRequisite_MS(field, nBCs)
    ! --------------------------------------------------------------------------
    !> fluid parameters and properties
    type( mus_field_type ), intent(inout) :: field(:)
    !> Number of boundary types
    integer, intent(in) :: nBCs
    ! --------------------------------------------------------------------------
    integer :: iField, iBnd, counter(nBCs)
    logical :: checkBnd(nBCs)
    character(len=labelLen) :: checkKind(nBCs)
    ! --------------------------------------------------------------------------
    counter = 0
    checkBnd = .false.
    do iField = 1, size(field)
      do iBnd = 1, nBCs
        select case (trim(field(iField)%bc(iBnd)%BC_kind))
        case ('spc_inlet_eq', 'spc_inlet', 'spc_outlet_vel')
          ! this boundaries can be applied only if all species has this kind
          counter(iBnd) = counter(iBnd) + 1
          checkBnd(iBnd) = .true.
          checkKind(iBnd) = field(iField)%bc(iBnd)%BC_kind
        end select
      end do !iBnd
    end do !iField

    do iBnd = 1, nBCs
      if (checkBnd(iBnd)) then
        if (counter(iBnd) /= size(field)) then
           call tem_abort( 'Error: Not all species has kind: ' &
             &            //trim(checkKind(iBnd)) )
        end if
      end if
    end do
  end subroutine check_BCs_preRequisite_MS
  ! **************************************************************************** !

  ! **************************************************************************** !
  !> Initialize the values required for the moments BC
  subroutine init_momentsBC( bc, leveldesc, layout, globBC, minLevel, maxLevel )
    ! ---------------------------------------------------------------------------
    !> Level range
    integer, intent(in) :: minLevel, maxLevel
    !> global array boundary type
    type( boundary_type ), intent(inout) :: bc
    !> Level descriptor
    type( tem_leveldesc_type ), intent(in) :: levelDesc(minLevel:maxLevel)
    !> Layout
    type( mus_scheme_layout_type), intent(in) :: layout
    !> scheme global boundary type
    type( glob_boundary_type ), intent(inout) :: globBC
    ! ---------------------------------------------------------------------------
    integer :: iElem, iLevel, iDir
    integer :: d2q9_xNormal(3), d2q9_yNormal(3)
    integer :: d3q19_xNormal(5), d3q19_yNormal(5), d3q19_zNormal(5)
    integer, allocatable, dimension(:) :: xNormal_mom, yNormal_mom,            &
      &                                   zNormal_mom,                         &
      &                                   xyNormal_mom, yzNormal_mom,          &
      &                                   xzNormal_mom,                        &
      &                                   xyzNormal_mom
    integer, allocatable, dimension(:,:) :: xNorm_links, yNorm_links,          &
      &                                     zNorm_links, xyNorm_links,         &
      &                                     yzNorm_links, xzNorm_links,        &
      &                                     xyzNorm_links
    character(len=labelLen) :: normalIndex
    real(kind=rk) :: normal(3)
    integer :: nLinks, iLink, normal_nLinks, edge_nLinks, corner_nLinks
    integer, allocatable :: missing_links(:)
    real(kind=rk), allocatable :: unKnown_Mat(:,:)
    integer :: elemPos
    logical :: bitmask( layout%fStencil%QQN )
    integer(kind=long_k), allocatable :: corner_elems(:)
    integer(kind=long_k) :: treeID
    logical :: corner_node, update_allMoments
    integer :: iCorner
    ! ---------------------------------------------------------------------------
!    write(dbgUnit(1),*) 'Boundary label ', trim(bc%label)

    ! @todo KM: move this generic info to tem_stencil_module or tem_param_module
    ! known moments positions in moments array
    select case (trim(bc%BC_kind))
    case('pressure_momentsbased', 'moments_outflow', 'spc_moments_molefrac')
      d2q9_xNormal = (/ 1, 5, 3/) !m0, mY, mYY
      d2q9_yNormal = (/ 1, 4, 2/) !m0, mX, mXX
      d3q19_xNormal = (/ 1, 3, 4, 6, 7/) !m0, mY, mZ, mYY, mZZ
      d3q19_yNormal = (/ 1, 2, 4, 5, 7/) !m0, mX, mZ, mXX, mZZ
      d3q19_zNormal = (/ 1, 2, 3, 5, 6/) !m0, mX, mY, mXX, mYY
    case('moments_wall','velocity_momentsbased','moments_inflow', &
      &  'spc_moments_moleflux','spc_moments_wall','spc_moments_vel', &
      & 'spc_moments_outflow' )
      ! ordering is based on such that 1st knownMoments pos is
      ! moment in normal direction and rest follows
      d2q9_xNormal = (/ 2, 3, 5/) !mX, mY, mYY
      d2q9_yNormal = (/ 3, 2, 4/) !mX, mY, mXX
      d3q19_xNormal = (/ 2, 3, 4, 6, 7/) !mX, mY, mZ, mYY, mZZ
      d3q19_yNormal = (/ 3, 2, 4, 5, 7/) !mX, mY, mZ, mXX, mZZ
      d3q19_zNormal = (/ 4, 2, 3, 5, 6/) !mX, mY, mZ, mXX, mYY
    case default
      d2q9_xNormal = 0
      d2q9_yNormal = 0
      d3q19_xNormal = 0
      d3q19_yNormal = 0
      d3q19_zNormal = 0
      write(logUnit(1),*)'This boundary kind is not supported'
      call tem_abort()
    end select

    select case (trim(layout%fStencil%label))
    case('d2q9')
      normal_nLinks = 3
      ! in 2d edge and corner or same
      edge_nLinks = 5
      corner_nLinks = 5
      !axis-normal
      allocate(xNormal_mom(normal_nLinks))
      allocate(yNormal_mom(normal_nLinks))
      !edge normal
      allocate(xyNormal_mom(edge_nLinks))
      xNormal_mom = d2q9_xNormal
      yNormal_mom = d2q9_yNormal
      xyNormal_mom = (/ 1, 2, 3, 4, 5/) !m0, mX, mY, mXX, mYY

      ! normal links
      allocate(xNorm_links(normal_nLinks,2))
      xNorm_links = reshape((/ 1, 5, 6,    & !x-
        &                      3, 7, 8 /), & !x+
        &                  (/normal_nLinks,2/))

      allocate(yNorm_links(normal_nLinks,2))
      yNorm_links = reshape((/ 2, 5, 7,    & !y-
        &                      4, 6, 8 /), & !y+
        &                  (/normal_nLinks,2/))

      allocate(xyNorm_links(edge_nLinks,4))
      xyNorm_links = reshape((/ 1, 2, 5, 6, 7,   & !x-,y-
        &                       1, 4, 5, 6, 8,   & !x-.y+
        &                       2, 3, 5, 7, 8,   & !x+,y-
        &                       3, 4, 6, 7, 8 /),& !x+,y+
        &                   (/edge_nLinks,4/))
      allocate(zNorm_links(normal_nLinks,2))
      zNorm_links = 0
      allocate(yzNorm_links(edge_nLinks,4))
      yzNorm_links = 0
      allocate(xzNorm_links(edge_nLinks,4))
      xzNorm_links = 0
      allocate(xyzNorm_links(corner_nLinks,8))
      xyzNorm_links = 0
    case('d3q19')
      normal_nLinks = 5
      edge_nLinks = 9
      corner_nLinks = 12
      !axis-normal
      allocate(xNormal_mom(normal_nLinks))
      allocate(yNormal_mom(normal_nLinks))
      allocate(zNormal_mom(normal_nLinks))
      !edge normal
      allocate(xyNormal_mom(edge_nLinks))
      allocate(yzNormal_mom(edge_nLinks))
      allocate(xzNormal_mom(edge_nLinks))
      !corner normal
      allocate(xyzNormal_mom(corner_nLinks))
      xNormal_mom = d3q19_xNormal
      yNormal_mom = d3q19_yNormal
      zNormal_mom = d3q19_zNormal
      ! m0, mX, mY, mZ, mXX, mYY, mZZ, mYZ, mZZX
      xyNormal_mom = (/ 1, 2, 3, 4, 5, 6, 7, 9, 15/)
      ! m0, mX, mY, mZ, mXX, mYY, mZZ, mXZ, mXXY
      yzNormal_mom = (/ 1, 2, 3, 4, 5, 6, 7, 10, 11/)
      ! m0, mX, mY, mZ, mXX, mYY, mZZ, mXY, mYYZ
      xzNormal_mom = (/ 1, 2, 3, 4, 5, 6, 7, 8, 14/)
      ! m0, mX, mY, mZ, mXX, mYY, mZZ, mXY, mYZ, mXXY, mYYX, mZZY
      xyzNormal_mom = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 16/)

      ! normal links
      allocate(xNorm_links(normal_nLinks,2))
      xNorm_links = reshape((/ 1, 11, 13, 15, 16,    & !x-
        &                      4, 12, 14, 17, 18 /), & !x+
        &                  (/normal_nLinks,2/))

      allocate(yNorm_links(normal_nLinks,2))
      yNorm_links = reshape((/ 2,  7,  8, 15, 17,    & !y-
        &                      5,  9, 10, 16, 18 /), & !y+
        &                  (/normal_nLinks,2/) )

      allocate(zNorm_links(normal_nLinks,2))
      zNorm_links = reshape((/ 3,  7,  9, 11, 12,    & !z-
        &                      6,  8, 10, 13, 14 /), & !z+
        &                  (/normal_nLinks,2/))

      allocate(xyNorm_links(edge_nLinks,4))
      xyNorm_links = reshape((/ 1,  2,  7,  8, 11, 13, 15, 16, 17,   & !x-,y-
        &                       1,  5,  9, 10, 11, 13, 15, 16, 18,   & !x-,y+
        &                       2,  4,  7,  8, 12, 14, 15, 17, 18,   & !x+,y-
        &                       4,  5,  9, 10, 12, 14, 16, 17, 18 /),& !x+,y+
        &                   (/edge_nLinks,4/))

      allocate(yzNorm_links(edge_nLinks,4))
      yzNorm_links = reshape((/ 2,  3,  7,  8,  9, 11, 12, 15, 17,   & !y-,z-
        &                       2,  6,  7,  8, 10, 13, 14, 15, 17,   & !y-,z+
        &                       3,  5,  7,  9, 10, 11, 12, 16, 18,   & !y+,z-
        &                       5,  6,  8,  9, 10, 13, 14, 16, 18 /),& !y+,z+
        &                   (/edge_nLinks,4/))

      allocate(xzNorm_links(edge_nLinks,4))
      xzNorm_links = reshape((/ 1,  3,  7,  9, 11, 12, 13, 15, 16,   & !x-,z-
        &                       1,  6,  8, 10, 11, 13, 14, 15, 16,   & !x-,z+
        &                       3,  4,  7,  9, 11, 12, 14, 17, 18,   & !x+,z-
        &                       4,  6,  8, 10, 12, 13, 14, 17, 18 /),& !x+,z+
        &                   (/edge_nLinks,4/))

      allocate(xyzNorm_links(corner_nLinks,8))
      xyzNorm_links = reshape((/ &
        &        1,  2,  3,  7,  8,  9, 11, 12, 13, 15, 16, 17,   & !x-,y-,z-
        &        1,  2,  6,  7,  8, 10, 11, 13, 14, 15, 16, 17,   & !x-,y-,z+
        &        1,  3,  5,  7,  9, 10, 11, 12, 13, 15, 16, 18,   & !x-,y+,z-
        &        1,  5,  6,  8,  9, 10, 11, 13, 14, 15, 16, 18,   & !x-,y+,z+
        &        2,  3,  4,  7,  8,  9, 11, 12, 14, 15, 17, 18,   & !x+,y-,z-
        &        2,  4,  6,  7,  8, 10, 12, 13, 14, 15, 17, 18,   & !x+,y-,z+
        &        3,  4,  5,  7,  9, 10, 11, 12, 14, 16, 17, 18,   & !x+,y+,z-
        &        4,  5,  6,  8,  9, 10, 12, 13, 14, 16, 17, 18 /),& !x+,y+,z+
        &                    (/corner_nLinks,8/))
    case default
      normal_nLinks = 0
      edge_nLinks   = 0
      corner_nLinks = 0
      allocate(xNormal_mom(0))
      allocate(yNormal_mom(0))
      allocate(zNormal_mom(0))
      allocate(xyNormal_mom(0))
      allocate(yzNormal_mom(0))
      allocate(xzNormal_mom(0))
      write(logUnit(1),*) trim(layout%fStencil%label)//&
        &' layout is not supported for this boundary'
      call tem_abort()
    end select


    do iLevel = minLevel, maxLevel
      allocate( bc%elemLvl( iLevel )%moments( globBC%nElems( iLevel )))

      ! list of corner elements treeIDs
      allocate(corner_elems(globBC%cornerBC%nElems(iLevel)))
      corner_elems =                                                           &
       & levelDesc(iLevel)%total(globBC%cornerBC%elemLvl(iLevel)%elem%val(:))

      iCorner = 0
      do iElem = 1, globBC%nElems( iLevel )
        ! when moments combination are not defined properly
        ! update all moments
        update_allMoments = .false.
        ! if corrent node is corner node i.e intersected by multiple boundaries
        corner_node = .false.

!write(dbgUnit(1),*) 'iElem ', iElem
        treeID = levelDesc(iLevel)%total(globBC%elemLvl(iLevel)%elem%val(iElem))
        if ( any(corner_elems == treeID) ) corner_node = .true.
!write(dbgUnit(1),*) 'corner_node ', corner_node

        ! update bitmask, normal and normalInd with cornerBC information
        ! if corrent node is corner node
        if (corner_node) then
          iCorner = iCorner + 1
          globBC%elemLvl(iLevel)%bitmask%val(:, iElem) = &
            & globBC%cornerBC%elemLvl(iLevel)%bitmask%val(:, iCorner)
          globBC%elemLvl(iLevel)%normal%val(:, iElem) =  &
            & globBC%cornerBC%elemLvl(iLevel)%normal%val(:, iCorner)
          globBC%elemLvl(iLevel)%normalInd%val(iElem) =  &
            & globBC%cornerBC%elemLvl(iLevel)%normalInd%val(iCorner)
        end if ! corner node

        ! element position in state array and treeID list
        elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
        ! normal of this corrent node
        normal = globBC%elemLvl(iLevel)%normal%val(:,iElem)
        bitmask = globBC%elemLvl(iLevel)%bitmask%val(:, iElem)

        nLinks = count(globBC%elemLvl(iLevel)%bitmask%val(:, iElem))
        bc%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs = nLinks

        allocate(bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos(nLinks))
        allocate(missing_links(nLinks))
        iLink = 0
        do iDir = 1,layout%fStencil%QQN
          if (globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem)) then
            iLink = iLink+1
            missing_links(iLink) = iDir
          endif
        end do

        normalIndex = 'default'

        ! axis normals
        if (nLinks == normal_nLinks) then
          do iDir = 1, 2
            if ( all(missing_links == xNorm_links(:, iDir)) ) &
              & normalIndex = 'xNormal'
            if ( all(missing_links == yNorm_links(:, iDir)) ) &
              & normalIndex = 'yNormal'
            if ( all(missing_links == zNorm_links(:, iDir)) ) &
              & normalIndex = 'zNormal'
          end do
        else if (nLinks == edge_nLinks) then
        ! edge direction
          do iDir = 1, 4
            if ( all(missing_links == xyNorm_links(:, iDir)) ) &
              & normalIndex = 'xyNormal'

            if ( all(missing_links == yzNorm_links(:, iDir)) ) &
              & normalIndex = 'yzNormal'

            if ( all(missing_links == xzNorm_links(:, iDir)) ) &
              & normalIndex = 'xzNormal'
          end do
        else if (nLinks == corner_nLinks) then
        ! corner direction
          do iDir = 1, 8
            if ( all(missing_links == xyzNorm_links(:,iDir)) ) &
              & normalIndex = 'xyzNormal'
          end do
        end if

!write(dbgUnit(1),*) 'normal ', globBC%elemLvl(iLevel)%normal%val(:, iElem)
!write(dbgUnit(1),*) 'normalIndex ', normalIndex
!write(dbgUnit(1),*) 'missing_links ', missing_links
        select case(trim(normalIndex))
        case ('xNormal')
          bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos = xNormal_mom
        case ('yNormal')
          bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos = yNormal_mom
        case ('zNormal')
          bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos = zNormal_mom
        case ('xyNormal')
          bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos = xyNormal_mom
        case ('yzNormal')
          bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos = yzNormal_mom
        case ('xzNormal')
          bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos = xzNormal_mom
        case ('xyzNormal')
          bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos = xyzNormal_mom
        case default
          ! undefined normal so update all links for this boundary
          globBC%elemLvl(iLevel)%bitmask%val(:, iElem) = .true.
          nLinks = layout%fStencil%QQ
          deallocate(bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos)
          allocate(bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos(nLinks))
          do iLink = 1,nLinks
            bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos(iLink) = iLink
          enddo

          deallocate(missing_links)
          allocate(missing_links(nLinks))
          iLink = 0
          do iDir = 1,layout%fStencil%QQN
            if (globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem)) then
              iLink = iLink+1
              missing_links(iLink) = iDir
            endif
          end do
        end select

        bc%elemLvl(iLevel)%moments(iElem)%nUnKnownPdfs = nLinks
        allocate(bc%elemLvl(iLevel)%moments(iElem)                   &
          &                      %unKnownPdfs_MatInv(nLinks,nLinks))
        allocate(unKnown_Mat(nLinks,nLinks))

        do iLink = 1, nLinks
          unKnown_Mat(:,iLink) = layout%moment%toMoments%A( &
            & bc%elemLvl(iLevel)%moments(iElem)%knownMom_pos, &
            & missing_links(iLink) )
        end do

        bc%elemLvl(iLevel)%moments(iElem)%unKnownPdfs_MatInv &
          & = invert_matrix(unKnown_Mat)

        deallocate(unKnown_Mat)
        deallocate(missing_links)

      end do ! iElem
    end do ! level

  end subroutine init_momentsBC
  ! **************************************************************************** !


! ****************************************************************************** !
  !> assign qVal to corresponding BC and level-wise if qVal not provided by
  !! seeder. qVal from seeder are assigned in assignBCList in
  !! mus_construction_module, So set qVal from config only when it is not
  !! provided by seeder.
  subroutine init_qVals(refQVal, layout, globBC, minLevel, maxLevel)
    ! ---------------------------------------------------------------------------
    real(kind=rk),                 intent(in) :: refQVal
    type( mus_scheme_layout_type), intent(in) :: layout
    type( glob_boundary_type ), intent(inout) :: globBC
    integer,                       intent(in) :: minLevel, maxLevel
    ! ---------------------------------------------------------------------------
    integer :: iElem, iLevel, iDir, invDir, QQN
    ! ---------------------------------------------------------------------------

    globBC%qValInitialized = .true.
    QQN = layout%fStencil%QQN

    ! 2. allocate qVal array level-wise
    lvlLoop: do iLevel = minLevel, maxLevel

      ! if qVal is not available from seeder then initialize qVal
      call init( me     = globBC%elemLvl(iLevel)%qVal, &
        &        width  = QQN,                         &
        &        length = globBC%nElems(iLevel)        )
        globBC%elemLvl( iLevel )%qVal%val = 0.0_rk

      ! loop over elements
      do iElem = 1, globBC%nElems(iLevel)
        ! loop over directions
        do iDir = 1, QQN
          ! Bitmask is true for incoming direction.
          ! so use invDir to access qVal and cxDirRK
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            invDir  = layout%fStencil%cxDirInv(iDir)
            ! if no qVal from seeder then use refQVal
            globBC%elemLvl(iLevel)%qVal%val(invDir, iElem) = refQVal
          end if
        end do !iDir
      end do !iElem
    end do lvlLoop

  end subroutine init_qVals
! ****************************************************************************** !


! ****************************************************************************** !
  !> Initialize the values required for the characteristic boundary conditions
  !!
  !! Simply calculate the macroscopic values from the current pdf distributions
  !! in the boundary element and its neighbor elements and store it to the
  !! %nrbc%lodi field
  !! NOTE:
  !! To be consistent, the same calculation procedure of the LODI values need to
  !! be performed as in the
  !! [[mus_bc_fluid_module:outlet_nrbc]] "NRBC routine" itself
  !!
  subroutine init_nrbc( bc, state, nSize, neigh, layout, level, &
    &                   nScalars, varSys, derVarPos, elemPos, nElems )
    ! ---------------------------------------------------------------------------
    !> global array boundary type
    type( boundary_type ) :: bc
    !> level
    integer, intent(in) :: level
    !> State array
    real(kind=rk), intent(in) :: state(:)
    !> nSize
    integer, intent(in) :: nSize
    !> neighbor array
    integer, intent(in) :: neigh(:)
    !> fluid parameters
    type( mus_scheme_layout_type), intent(in) :: layout
    !> number of scalars in global system
    integer, intent(in) :: nScalars
    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys
    !> position of derived quantities in varsys
    type( mus_derVarPos_type ), intent(in) :: derVarPos
    !> BC elements positions in total list
    integer, intent(in) :: elemPos(:)
    !> number of BC elements
    integer, intent(in) :: nElems
    ! ---------------------------------------------------------------------------
    integer :: iElem, iDir, iField, QQ
    integer :: iNeigh, neighPos, posInTotal
    real(kind=rk) :: ff( layout%fStencil%QQ )
    ! ---------------------------------------------------------------------------

    QQ = layout%fStencil%QQ
    iField = 1

    ! Assign the initial values to the lodi variables
    allocate( bc%elemLvl( level )%lodi( 4, 3, nElems))
    bc%elemLvl( level )%lodi = 0._rk

    ! reset
    do iElem = 1, nElems
      posInTotal = elemPos( iElem )

      ! adopt initial density from the fluid cells
      do iDir = 1, QQ
        ff(iDir) = state( (posintotal-1)*nscalars+idir+(ifield-1)*qq )
      end do !iDir
      bc%elemLvl( level )%lodi( 1, 1, iElem ) = sum( ff )

      ! adapt initial velocities from the fluid cells
      call derVarPos%velFromState(                                &
        &           state  = ff,                                  &
        &           iField = iField,                              &
        &           nElems = 1,                                   &
        &           varSys = varSys,                              &
        &           layout = layout,                              &
        &           res    = bc%elemLvl(level)%lodi(2:4, 1, iElem) )

      do iNeigh = 1,2
        ! Get the values for the neighbors as well
        neighPos = bc%neigh( level )%posInState( iNeigh, iElem )
        if( neighPos > 0 ) then
          do iDir = 1,QQ
            ff(iDir) = state( (neighpos-1)*nscalars+idir+(ifield-1)*qq)
          end do
          bc%elemLvl( level )%lodi( 1, 1+iNeigh, iElem ) = sum( ff )

          ! velocities
          call derVarPos%velFromState(                                  &
            &       state  = ff,                                        &
            &       iField = iField,                                    &
            &       nElems = 1,                                         &
            &       varSys = varSys,                                    &
            &       layout = layout,                                    &
            &       res    = bc%elemLvl(level)%lodi(2:4, 1+iNeigh, iElem) )
        end if
      end do ! iNeigh
    end do ! iElem

  end subroutine init_nrbc
! ****************************************************************************** !


  ! -------------------------------------------------------------------------- !
  !> This routine allocates turbulent viscosity and friction velocity on
  !! boundary elements. It also initialize friction velocity from stream-wise
  !! velocity component on first neighbor in normal direction.
  subroutine mus_init_turb_wallFunc(bc, globBC, auxField, viscKine, derVarPos,&
    &                                varSys, stencil, iLevel)
    ! -------------------------------------------------------------------- !
    !> Field bc which contains turbwallFunc type and neighbor info
    type(boundary_type), intent(inout) :: bc
    !> global bc of current boundary with elemPos and normal info
    type(glob_boundary_type), intent(inout) :: globBC
    !> auxField array
    real(kind=rk), intent(in) :: auxField(:)
    !> Kinematic viscosity
    type(mus_viscosity_type) :: viscKine
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in) :: derVarPos
    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys
    !> stencil info
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> current level
    integer, intent(in) :: iLevel
    ! -------------------------------------------------------------------- !
    integer :: iElem, neighPos, elemOff, elemPos, normalInd_inv, vel_pos(3)
    integer :: nElems
    real(kind=rk) :: velNeigh(3), normal(3), streamwise(3), unitSW(3)
    real(kind=rk) :: velSW, normalMag, qVal, distToBnd, streamwise_mag
    real(kind=rk) :: unitNormal(3)
    ! -------------------------------------------------------------------- !
    nElems = globBC%nElems(iLevel)
    allocate( bc%turbwallFunc%dataOnLvl(iLevel)%tVisc(nElems) )
    bc%turbwallFunc%dataOnLvl(iLevel)%tVisc = 0.0_rk

    allocate( bc%turbwallFunc%dataOnLvl(iLevel)%velTau(nElems) )
    bc%turbwallFunc%dataOnLvl(iLevel)%velTau = 0.0_rk

    allocate( bc%turbwallFunc%dataOnLvl(iLevel)%distToBnd(nElems) )
    bc%turbwallFunc%dataOnLvl(iLevel)%distToBnd = 0.0_rk

    allocate( bc%turbwallFunc%dataOnLvl(iLevel)%neighDistToBnd(nElems) )
    bc%turbwallFunc%dataOnLvl(iLevel)%neighDistToBnd = 0.0_rk

    allocate( bc%turbwallFunc%dataOnLvl(iLevel)%unitNormal(3, nElems) )
    bc%turbwallFunc%dataOnLvl(iLevel)%unitNormal = 0.0_rk

    allocate( bc%turbwallFunc%dataOnLvl(iLevel)%bndForce(3, nElems) )
    bc%turbwallFunc%dataOnLvl(iLevel)%bndForce = 0.0_rk

    allocate( bc%turbwallFunc%dataOnLvl(iLevel)%bndMoment(3, nElems) )
    bc%turbwallFunc%dataOnLvl(iLevel)%bndMoment = 0.0_rk

    ! Position of velocity in auxField array
    vel_pos = varSys%method%val(derVarPos%velocity)%auxField_varPos(1:3)

    do iElem = 1, nElems
      ! Get velocity on the first neighbor element from auxField array
      neighPos = bc%neigh(iLevel)%posInState(1, iElem)
      elemOff = (neighPos-1)*varSys%nAuxScalars
      velNeigh(1) = auxField(elemOff + vel_pos(1))
      velNeigh(2) = auxField(elemOff + vel_pos(2))
      velNeigh(3) = auxField(elemOff + vel_pos(3))

      ! Calculate local stream-wise unit vector
      ! e' = (u - (u . n) . n) / |(u - (u . n) . n|
      !
      ! NormalInd in bc elem is pointing into the domain so we invert it
      normalInd_inv = stencil%cxDirInv( globBC%elemLvl(iLevel)  &
        &                                 %normalInd%val(iElem) )
      normal = stencil%cxDirRK(:, normalInd_inv)
      normalMag = sqrt(dot_product(normal, normal))
      unitNormal = normal / normalMag
      bc%turbwallFunc%dataOnLvl(iLevel)%unitNormal(:,iElem) = unitNormal
      streamwise = velNeigh - dot_product( velNeigh, unitNormal ) * unitNormal
      streamwise_mag = sqrt(dot_product(streamwise, streamwise))
      unitSW = 0.0_rk
      ! if velocity is zero then compute streamwise unit vector from bc normal
      if (streamwise_mag .fne. 0.0_rk) then
        unitSW = streamwise / streamwise_mag
      end if
      ! Unit vector

      ! stream-wise velocity component
      velSW = dot_product(velNeigh, unitSW)

      ! use qValue in normal direction as distance of first node from boundary
      qVal = globBC%elemLvl(iLevel)%qVal%val(normalInd_inv, iElem)
      distToBnd = qVal + normalMag
      bc%turbwallFunc%dataOnLvl(iLevel)%distToBnd(iElem) = qVal
      bc%turbwallFunc%dataOnLvl(iLevel)%neighDistToBnd(iElem) = distToBnd

      ! Initialize friction velocity on first neighbor with werner and wengle
      ! wall model in linear profile i.e. u+=y+ -> u_tau = sqrt(u * nu / y)
      elemPos = globBC%elemLvl(iLevel)%elem%val(iElem)
      bc%turbwallFunc%dataOnLvl(iLevel)%velTau(iElem)              &
        & = sqrt( velSW * viscKine%dataOnLvl(iLevel)%val( elemPos ) &
        &        / distToBnd )
    end do
  end subroutine mus_init_turb_wallFunc
  ! -------------------------------------------------------------------------- !


! ****************************************************************************** !
  !> Transfer pre- and post-collision PDF of neighbors of boundary elements
  !! into neighBufferPre and neighBufferPost.
  !! Access to state array
  !! ---------------------
  !!
  !! * NeighBufferPre: FETCH
  !! * NeighBufferPost: SAVE
  !! @todo: for PUSH, post-collision is not correct yet.
  !!
  subroutine fill_neighBuffer( prevstate, currstate, neigh, globBC, nBCs, &
    &                          field, varSys, QQ, nSize, iLevel )
    !---------------------------------------------------------------------------
    !> Previous state vector of iLevel
    real(kind=rk),intent(in) :: prevstate(:)
    !> Current state vector of iLevel
    real(kind=rk),intent(in) :: currstate(:)
    !> connectivity array corresponding to state vector
    integer,      intent(in) :: neigh(:)
    !> scheme global boundary type
    type( glob_boundary_type ), intent(in) :: globBC(:)
    !> number of total BC
    integer, intent(in) :: nBCs
    !> fluid parameters and properties
    type( mus_field_type ), intent(inout) :: field(:)
    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys
    !> number of total links
    integer, intent(in) :: QQ
    !> number of elements in state vector
    integer, intent(in) :: nSize
    !> the iLevel on which this boundary was invoked
    integer, intent(in) :: iLevel
    ! ---------------------------------------------------------------------------
    integer :: iField, iElem, iDir, iBnd, iNeigh
    integer :: pdfVarPos(QQ), neighPos, myPos
    ! ---------------------------------------------------------------------------
    ! write(dbgUnit(5),*) 'Fill neighBufferPre and neighBufferPost'

    do iField = 1, size( field )

      pdfVarPos(:) = varSys%method%val(iField)%state_varPos(1:QQ)

      ! fill boundary neighbor pre- and post-collision state values
      do iBnd = 1, nBCs
        call tem_startTimer( timerHandle = mus_timerHandles%setBnd(iBnd) )
        if (field(iField)%bc( iBnd )%requireNeighBufPre) then
! write(dbgUnit(5),"(A,I0)") 'iBC: ', iBnd
          do iNeigh = 1, field(iField)%bc( iBnd )%nNeighs
! write(dbgUnit(5),"(A,I0)") 'iNeigh: ', iNeigh
            do iElem = 1, globBC( iBnd )%nElems( iLevel )
              ! neighbor position in total list
              neighPos =  field(iField)%bc( iBnd )%neigh( iLevel )          &
                &                                 %posInState( iNeigh, iElem )
! write(dbgUnit(5),"(2(A,I0))") 'iElem: ', iElem, ', neighPos: ', neighPos
              do iDir = 1, QQ
                ! pre-collision neigh
                field(iField)%bc( iBnd )%neigh( iLevel )%neighBufferPre( &
  & iNeigh, (iElem-1)*QQ+iDir ) = prevstate(                             &
  &  neigh (( idir-1)* nsize+ neighpos)+( ifield-1)* qq+ varsys%nscalars*0 )
! write(dbgUnit(5),"(A,I0)") 'iDir: ', iDir
! write(dbgUnit(5),"(3(A,I0),A)") &
!   & 'neighBufferPre(', iNeigh, ',',(iElem-1)*QQ+iDir, ') = state(', &
!   &  neigh (( idir-1)* nsize+ neighpos)+( ifield-1)* qq+ varsys%nscalars*0, &
!   & ')'
              end do ! iDir
            end do ! iElem
          end do ! iNeigh
        end if ! neighBufferPre

        if (field(iField)%bc( iBnd )%requireNeighBufPre_nNext) then
! write(dbgUnit(5),"(A,I0)") 'iBC: ', iBnd
          do iNeigh = 1, field(iField)%bc( iBnd )%nNeighs
! write(dbgUnit(5),"(A,I0)") 'iNeigh: ', iNeigh
            do iElem = 1, globBC( iBnd )%nElems( iLevel )
              ! neighbor position in total list
              neighPos =  field(iField)%bc( iBnd )%neigh( iLevel )          &
                &                                 %posInState( iNeigh, iElem )
! write(dbgUnit(5),"(2(A,I0))") 'iElem: ', iElem, ', neighPos: ', neighPos
              do iDir = 1, QQ
                ! pre-collision neigh
                field(iField)%bc( iBnd )%neigh( iLevel )%neighBufferPre_nNext( &
  & iNeigh, (iElem-1)*QQ+iDir ) = currstate(                                   &
  &  neigh (( idir-1)* nsize+ neighpos)+( ifield-1)* qq+ varsys%nscalars*0 )
! write(dbgUnit(5),"(A,I0)") 'iDir: ', iDir
! write(dbgUnit(5),"(3(A,I0),A)") &
!   & 'neighBufferPre_nNext(', iNeigh, ',',(iElem-1)*QQ+iDir, ') = state(', &
!   &  neigh (( idir-1)* nsize+ neighpos)+( ifield-1)* qq+ varsys%nscalars*0, &
!   & ')'
              end do ! iDir
            end do ! iElem
          end do ! iNeigh
        end if ! neighBufferPre_nNext

        if (field(iField)%bc( iBnd )%requireNeighBufPost) then
! write(dbgUnit(5),"(A,I0)") 'iBC: ', iBnd
          do iNeigh = 1, field(iField)%bc( iBnd )%nNeighs
! write(dbgUnit(5),"(A,I0)") 'iNeigh: ', iNeigh
            do iElem = 1, globBC( iBnd )%nElems( iLevel )
              ! neighbor position in total list
              neighPos =  field(iField)%bc( iBnd )%neigh( iLevel )          &
                &                                 %posInState( iNeigh, iElem )
! write(dbgUnit(5),"(2(A,I0))") 'iElem: ', iElem, ', neighPos: ', neighPos
              do iDir = 1, QQ
                ! for PULL, this is post-collision value
                ! for PUSH, this is  pre-collision value
                field(iField)%bc( iBnd )%neigh( iLevel )%neighBufferPost( &
  & iNeigh, (iElem-1)*QQ+iDir ) = currstate(          &
  & ( neighpos-1)* varsys%nscalars+ idir+( ifield-1)* qq )
!  & ( neighpos-1)* varsys%nscalars+ pdfvarpos(idir))

! write(dbgUnit(5),"(A,I0)") 'iDir: ', iDir
! write(dbgUnit(5),"(3(A,I0),A)") &
!   & 'neighBufferPost(', iNeigh, ',',(iElem-1)*QQ+iDir, ') = state(', &
!   & ( neighpos-1)* varsys%nscalars+ pdfvarpos(idir), &
!   & ')'
              end do ! iDir
            end do ! iElem
          end do ! iNeigh
        end if !neighBufferPost

        if ( field(iField)%bc( iBnd )%useComputeNeigh ) then
          do iElem = 1, globBC( iBnd )%nElems( iLevel )
            ! my (bnd element) position in total list
            myPos = globBC( iBnd )%elemLvl( iLevel )%elem%val( iElem )
            do iDir = 1, QQ
              ! get post-collision state values of my neighbors
              ! computeNeighBuf always uses AOS layout
field(iField)%bc(iBnd)%neigh(iLevel)%computeNeighBuf( (iElem-1)*QQ+iDir ) =    &
& currstate(  neigh (( idir-1)* nsize+ mypos)+( ifield-1)* qq+ varsys%nscalars*0 )
            end do ! iDir
          end do ! iElem
        end if ! useComputeNeigh

        call tem_stopTimer( timerHandle = mus_timerHandles%setBnd(iBnd) )
      end do ! iBnd
    end do ! iField

  end subroutine fill_neighBuffer
! ****************************************************************************** !


! ****************************************************************************** !
  !> Transfer pdf of boundary elements into bcBuffer which is used by all
  !! boundary routines.
  ! for PULL, this is post-collision value
  ! for PUSH, this is  pre-collision value
  subroutine fill_bcBuffer( bcBuffer, currState, neigh, nSize, nElems_bc, &
    &                       posInTotal, nFields, QQ, varSys )
    ! ---------------------------------------------------------------------------
    !> state values of all boundary elements
    real(kind=rk),intent(out) :: bcBuffer(:)
    !> Current state vector
    real(kind=rk),intent(in) :: currState(:)
    !> connectivity array corresponding to state vector
    integer,      intent(in) :: neigh(:)
    !> nSize
    integer, intent(in) :: nSize
    !> number of boundary elements
    integer, intent(in) :: nElems_bc
    !> positions in total list of boundary elements
    integer, intent(in) :: posInTotal(nElems_bc)
    !> Number of fields
    integer, intent(in) :: nFields
    !> number of total links
    integer, intent(in) :: QQ
    !> scheme variable system
    type( tem_varSys_type ), intent(in) :: varSys
    ! ---------------------------------------------------------------------------
    integer :: iElem, iField, varPos, iDir
    ! ---------------------------------------------------------------------------

    ! write(dbgUnit(10), "(A)") ' Fill bcBuffer. '

    ! Loop over fields
    do iField = 1, nFields
      !NEC$ ivdep
      !dir$ ivdep
      !ibm* independent
      do iElem = 1, nElems_bc
        do iDir = 1, QQ
          ! bcBuffer always uses AOS data structure
          varPos = varSys%method%val(iField)%state_varPos(iDir)
          bcBuffer( varPos+(iElem-1)*varSys%nScalars ) = currState( &
& ( posintotal(ielem)-1)* varsys%nscalars+ idir+( ifield-1)* qq )
        end do
      end do ! iElem
    end do

  end subroutine fill_bcBuffer
! ****************************************************************************** !

! ***************************************************************************** !
  !> Get Surface points on boundary elements.
  !! For boundary state variable which are evaluated linkwise, extract surface
  !! points for each link and for non-link based variables project barycenter
  !! on the boundary surface.
  !! Return real coordinates on boundary surface and offset bit which encodes
  !! direction.
  subroutine mus_get_points_fromBC( bc, globBC, tree, stencil, total, iLevel, &
    &                               nPoints, points, offset_bit )
    !---------------------------------------------------------------------------
    !> Field boundary type
    type(boundary_type), intent(in) :: bc
    !> for number of elements in boundary and position in buffer
    type(glob_boundary_type), intent(in)     :: globBC
    !> global treelm mesh
    type(treelmesh_type), intent(in) ::tree
    !> global pdf type
    ! type(tem_levelDesc_type), intent(in) :: levelDesc
    integer(kind=long_k), intent(in) :: total(:)
    !> for directions
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> Current level
    integer, intent(in) :: iLevel
    !> Number of points
    integer, intent(out) :: nPoints
    !> 3-d real coordinates on which boundary variables are evaluated
    real(kind=rk), allocatable, intent(out) :: points(:,:)
    !> Offset bit encodes direction of boundary.
    !! used by apesmate to translate space coordinate in the offset direction
    !! to determine the treeID in remote domain
    character, allocatable, intent(out) :: offset_bit(:)
    !---------------------------------------------------------------------------
    real(kind=rk) :: dx, bary(3)
    real(kind=rk), allocatable :: qVal(:)
    integer :: iElem, iDir, iPnt, invDir
    logical :: qValIsZero
    !---------------------------------------------------------------------------
    write(logUnit(10),*) ' Get points from BC '
    ! Element size on current level
    dx = tem_ElemSizeLevel(tree, iLevel)
    allocate(qVal(stencil%QQN))

    ! if qVal is not provided by seeder and default is zero then no need to
    ! extract point on boundary surface. Just extract on barycenter
    if ((bc%qVal .feq. 0.0_rk) .and. .not. globBC%hasQVal &
      & .and. .not. bc%curved) then
      qValIsZero = .true.
    else
      qValIsZero = .false.
    end if

    ! Get qValue if boundary variables are evaulated link wise barycenter
    ! projection on boundary surface only if qVal is not zero
    if (bc%evalBcVar_link .and. .not. qValIsZero) then
      nPoints = bc%links(iLevel)%nVals
    else
      ! boundary variables are evaluated at bary center projection on
      ! boundary surface along boundary normal
      nPoints = globBC%nElems(iLevel)
    end if

    ! if no points on this level for this boundary, leave this routine
    ! if (nPoints==0) return

    allocate(points(nPoints,3))
    allocate(offset_bit(nPoints))

    iPnt = 0
    ! loop over elements
    do iElem = 1, globBC%nElems(ilevel)
      bary(:) = tem_BaryOfId( tree, total(globbc%elemLvl(iLevel)  &
        &                                       %elem%val(iElem)) )

      ! Set qValue to 0.5 if qValue is not defined.
      ! Mostly qVal is already set for link based boundary conditions
      if( .not. globBC%hasQVal .and.        &
        & .not. globBC%qValInitialized ) then
        qVal = bc%qVal
      else
        qVal = globBC%elemLvl(iLevel)%qVal%val(:, iElem)
      end if

      ! boundary variables are evaluated at link wise projection of barycenter
      ! on boundary surface
      if (bc%evalBcVar_link .and. .not. qValIsZero) then

        ! loop over directions
        do iDir = 1, stencil%QQN
          ! Bitmask is true for incoming direction.
          ! so use invDir to access qVal and cxDirRK
          if ( globBC%elemLvl(iLevel)%bitmask%val(iDir, iElem) ) then
            iPnt   = iPnt + 1
            invDir  = stencil%cxDirInv(iDir)
            ! position of boundary surface
            points(iPnt, :) = bary(:) + dx * qVal(invDir)  &
              &             * stencil%cxDirRK(:, invDir)

            offset_bit(iPnt) = qOffset_inChar( stencil%map( invDir ) )
          end if
        end do !iDir
      elseif (bc%BC_kind == 'bc_pdf') then
        ! For this boundary, boundary elements must overlap between two domains
        ! so return barycenter and offset 0
        iPnt = iPnt + 1
        points(iPnt, :) = bary(:)
        offset_bit(iPnt) = qOffset_inChar(q000)
      else

        ! boundary variables are evaluated at bary center projection on
        ! boundary surface along boundary normal
        iPnt   = iPnt + 1
        invDir  = stencil%cxDirInv(globBC%elemLvl(iLevel)%normalInd%val(iElem))
        if (invDir<=stencil%QQN) then
          points(iPnt, :) = bary(:) + dx * qVal(invDir)   &
            &             * stencil%cxDirRK(:, invDir)
          offset_bit(iPnt) = qOffset_inChar( stencil%map( invDir ) )
        else
          points(iPnt, :) = bary(:)
          offset_bit(iPnt) = qOffset_inChar(q000)
        end if

      end if !bc%evalBcVar_link
    end do !iElem

  end subroutine mus_get_points_fromBC
! ***************************************************************************** !


! ***************************************************************************** !
  !> This routine setup indices for boundary variables in bc_State_type
  !! pntIndex for the points on which boundaries are treated.
  subroutine mus_setupIndices_forBC( bc, globBC, tree, stencil, levelDesc, &
    &                                varSys, minLevel, maxLevel)
    ! --------------------------------------------------------------------------
    !> Field boundary type
    type(boundary_type), target, intent(inout) :: bc
    !> for number of elements in boundary and position in buffer
    type(glob_boundary_type), intent(in)     :: globBC
    !> global treelm mesh
    type(treelmesh_type), intent(in) ::tree
    !> Min and Max level in mesh
    integer, intent(in) :: minLevel, maxLevel
    !> global pdf type
    type(tem_levelDesc_type), intent(in) :: levelDesc(minLevel:maxLevel)
    !> for directions
    type(tem_stencilHeader_type), intent(in) :: stencil
    !> Global variable system
    type(tem_varSys_type), intent(in) :: varSys
    ! --------------------------------------------------------------------------
    !> Number of points
    integer :: nPoints
    !> 3-d real coordinates on which boundary variables are evaluated
    real(kind=rk), allocatable :: points(:,:)
    !> Offset bit encodes direction of boundary.
    !! used by apesmate to translate space coordinate in the offset direction
    !! to determine the treeID in remote domain
    character, allocatable :: offset_bit(:)
    integer, allocatable :: idx(:)
    integer :: iLevel, iVar
    character(len=labelLen) :: bc_varName
    type(tem_bc_state_type), pointer :: bc_state => NULL()
    character(len=pathLen) :: isSurface
    ! --------------------------------------------------------------------------


    write(logUnit(10),*) ' Setup indices for BC: ', trim(bc%label)
    do iLevel = minLevel, maxLevel
      write(logUnit(10),*) 'iLevel: ', iLevel
      call mus_get_points_fromBC( bc         = bc,                      &
        &                         globBC     = globBC,                  &
        &                         tree       = tree,                    &
        &                         stencil    = stencil,                 &
        &                         total      = levelDesc(iLevel)%total, &
        &                         iLevel     = iLevel,                  &
        &                         nPoints    = nPoints,                 &
        &                         points     = points,                  &
        &                         offset_bit = offset_bit               )

      ! if no points on this level for this boundary, goto next level
      ! KM: Skipping append idx results in not allocated indexLvl%val array
      ! which causes problems in calling get_valOfIndex for intel compiler
      !if (nPoints==0) cycle

      allocate(idx(nPoints))
      ! store points in spacetime function and store indices of points
      ! in bc_state
      do iVar = 1, bc%varDict%nVals
        bc_varName = trim(bc%varDict%val(iVar)%key)
        write(logUnit(10),*) ' Storing index for variable: ', trim(bc_varName)
        isSurface = 'isSurface = true'

        bc_state => NULL()
        select case (trim(bc_varName))
        case ('velocity')
          bc_state => bc%bc_states%velocity
        case ('pdf')
          bc_state => bc%bc_states%pdf
          ! For pdf boundary, both domain boundary layer must overlap
          isSurface = 'isSurface = false'
        case ('pressure')
          bc_state => bc%bc_states%pressure
        case ('mass_flowrate')
          bc_state => bc%bc_states%massFlowRate
        case ('mole_fraction')
          bc_state => bc%bc_states%moleFrac
        case ('mole_density')
          bc_state => bc%bc_states%moleDens
        case ('mole_flux')
          bc_state => bc%bc_states%moleFlux
        case ('mole_diff_flux')
          bc_state => bc%bc_states%moleDiff_flux
        case ('potential')
          bc_state => bc%bc_states%potential
        case ('surface_charge_density')
          bc_state => bc%bc_states%surChargeDens
        case default
          write(logUnit(1),*) 'Error: Unknown boundary variable: "'// &
            & trim(bc_varName)//'"'
          call tem_abort()
        end select

        if (.not. associated(bc_state)) then
          call tem_abort('Error: bc_state is not assosiated')
        end if

        ! set params
        write(logUnit(10),*) ' Set params in bc variable: '          &
          &                //trim(varSys%varName%val(bc_state%varPos))
        call varSys%method%val(bc_state%varPos)%set_params( &
          & varSys   = varSys,                              &
          & instring = isSurface                            )

        call varSys%method%val(bc_state%varPos)%setup_indices( &
          & varSys     = varSys,                               &
          & point      = points,                               &
          & offset_bit = offset_bit,                           &
          & iLevel     = iLevel,                               &
          & tree       = tree,                                 &
          & nPnts      = nPoints,                              &
          & idx        = idx                                   )

        if (nPoints == 0) then
          ! KM: Intel compiler fails if indexLvl is not allocated and accessed
          ! in get_valOfIndex
          ! initialize array with size zero
          call init(bc_state%pntIndex%indexLvl(iLevel))
        else
          call append(bc_state%pntIndex%indexLvl(iLevel), idx)
        end if
      end do !iVar
      deallocate(idx)
    end do !iLevel

  end subroutine mus_setupIndices_forBC
! ***************************************************************************** !


end module mus_bc_general_module
! ****************************************************************************** !
