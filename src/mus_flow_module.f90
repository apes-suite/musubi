! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2011-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2011-2013, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012-2015 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2016-2017 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Philipp Otte <otte@mathcces.rwth-aachen.de>
! Copyright (c) 2016-2018 Raphael Haupt <raphael.haupt@uni-siegen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
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
! **************************************************************************** !
!> Flow-related routines such as initialization of the flow field and parameters
!!
module mus_flow_module

  ! include treelm modules
  use env_module,               only: rk, io_buffer_size, long_k, &
    &                                 single_k, eps_single_k
  use tem_param_module,         only: cs2inv, cs2, t2cs2inv, t2cs4inv
  use treelmesh_module,         only: treelmesh_type
  use tem_aux_module,           only: tem_abort
  use tem_global_module,        only: tem_global_type
  use tem_geometry_module,      only: tem_baryOfId
  use tem_spatial_module,       only: tem_spatial_for
  use tem_subTree_type_module,  only: tem_subTree_type, tem_local_subTree_from,&
    &                                 tem_destroy_subTree
  use tem_logging_module,       only: logUnit
  use tem_debug_module,         only: dbgUnit
  use tem_tools_module,         only: tem_horizontalSpacer
  use tem_math_module,          only: invert_matrix
  use tem_ini_condition_module, only: tem_ini_condition_type
  use tem_varSys_module,        only: tem_varSys_type
  use tem_varMap_module,        only: tem_create_varMap
  use tem_general_module,       only: tem_general_type
  use tem_stencil_module,       only: tem_stencilHeader_type

  ! include musubi modules
  ! use mus_comm_module,               only: mus_exchange
  use mus_param_module,              only: mus_param_type
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_derivedQuantities_module2, only: getNEq_diffusive,                   &
    &                                      getNEq_acoustic
  use mus_field_module,              only: mus_field_type
  use mus_fluid_module,              only: mus_fluid_type
  use mus_mixture_module,            only: mus_mixture_type
  use mus_initfluid_module,          only: mus_init_advRel_fluid, &
                                         & mus_init_advRel_fluid_GNS
  use mus_initfluidIncomp_module,    only: mus_init_advRel_fluidIncomp, &
                                         & mus_init_advRel_fluidIncomp_GNS
  use mus_initLBMPS_module,          only: mus_init_advRel_LBM_PS
  use mus_initMultispecies_module,   only: mus_init_advRel_multispecies_gas,   &
    &                                      mus_init_advRel_multispecies_liquid
  use mus_initIsothermAcEq_module,   only: mus_init_advRel_isotherm_acEq
  use mus_initPoisson_module,        only: mus_init_advRel_Poisson,  &
    &                                      mus_init_advRel_PBLinear, &
    &                                      mus_init_advRel_PBnonLinear
  use mus_initNernstPlanck_module,   only: mus_init_advRel_nernstPlanck
  use mus_eNRTL_module,              only: mus_calc_thermFactor,               &
    &                                      mus_calc_MS_DiffMatrix
  use mus_restart_module,            only: mus_readRestart
  use mus_interpolate_verify_module, only: mus_testInterpolation
  use mus_physics_module,            only: mus_convertFac_type, mus_physics_type
  use mus_nernstPlanck_module,       only: mus_nernstPlanck_type
  use mus_derVarPos_module,          only: mus_derVarPos_type
!  use mus_gradData_module,      only: mus_gradData_type
!  use mus_derivedQuantities_module2, only: getGradU
  use mus_auxField_module,           only: mus_initAuxFieldFluidAndExchange,   &
    &                                      mus_intpAuxFieldCoarserAndExchange, &
    &                                      mus_intpAuxFieldFinerAndExchange

  implicit none

  private

  public :: mus_init_flow
  public :: fillHelperElementsCoarseToFine
  public :: fillHelperElementsFineToCoarse
  public :: mus_initAuxField

contains

! **************************************************************************** !
  !> Choose the relaxation model
  !!
  subroutine init_advRel( scheme )
    ! --------------------------------------------------------------------------
    !> scheme type
    type(mus_scheme_type) :: scheme
    ! --------------------------------------------------------------------------

    write(logUnit(1),*) 'Initialize compute kernel for scheme kind: '&
      &                   //trim( scheme%header%kind )

    select case( trim(scheme%header%kind) )
    case ('fluid')
      call mus_init_advRel_fluid( relaxation = scheme%header%relaxation, &
        &                         variant    = scheme%header%relaxHeader &
        &                                                   %variant,    &
        &                         layout     = scheme%header%layout,     &
        &                         compute    = scheme%compute            )
    case ('fluid_incompressible')
      call mus_init_advRel_fluidIncomp( relaxation = scheme%header%relaxation, &
        &                               variant    = scheme%header%relaxHeader &
        &                                                         %variant,    &
        &                               layout     = scheme%header%layout,     &
        &                               compute    = scheme%compute            )
    case ('fluid_GNS')
      call mus_init_advRel_fluid_GNS( relaxation = scheme%header%relaxation, &
        &                             layout     = scheme%header%layout,     &
        &                             compute    = scheme%compute            )
    case ('fluid_incompressible_GNS')
      call mus_init_advRel_fluidIncomp_GNS( relaxation = scheme%header%relaxation, &
        &                                   layout     = scheme%header%layout,     &
        &                                   compute    = scheme%compute            )
    case ('multispecies_gas')
      call mus_init_advRel_multispecies_gas(        &
        &    relaxation = scheme%header%relaxation, &
        &    layout     = scheme%header%layout,     &
        &    nFields    = scheme%nFields,           &
        &    compute    = scheme%compute            )
    case ('multispecies_liquid')
      call mus_init_advRel_multispecies_liquid(     &
        &    relaxation = scheme%header%relaxation, &
        &    layout     = scheme%header%layout,     &
        &    nFields    = scheme%nFields,           &
        &    compute    = scheme%compute            )
    case ('passive_scalar')
      ! lattice boltzmann passive scalar
      call mus_init_advRel_lbm_ps( relaxation         = scheme%header%relaxation,         &
        &                          layout             = scheme%header%layout,             &
        &                          relaxation_variant = scheme%header%relaxHeader%variant,&
        &                          compute            = scheme%compute )
    case ('nernst_planck')
      call mus_init_advRel_nernstPlanck(            &
        &    relaxation = scheme%header%relaxation, &
        &    layout     = scheme%header%layout,     &
        &    compute    = scheme%compute            )
    case ('poisson')
      call mus_init_advRel_Poisson( relaxation = scheme%header%relaxation, &
        &                           layout     = scheme%header%layout,     &
        &                           compute    = scheme%compute            )
    case ('poisson_boltzmann_linear')
      call mus_init_advRel_PBLinear( relaxation = scheme%header%relaxation, &
        &                            layout     = scheme%header%layout,     &
        &                            compute    = scheme%compute            )
    case ('poisson_boltzmann_nonlinear')
      call mus_init_advRel_PBnonLinear( relaxation = scheme%header%relaxation, &
        &                               layout     = scheme%header%layout,     &
        &                               compute    = scheme%compute            )
    case ('isotherm_acEq')
      ! lattice boltzmann model for the isothermal acoustic equations
      call mus_init_advRel_isotherm_acEq(           &
        &    relaxation = scheme%header%relaxation, &
        &    layout     = scheme%header%layout,     &
        &    compute    = scheme%compute            )
    case default
      write(logUnit(1),*) 'The selected scheme kind model is not supported: ' &
        &                 //trim( scheme%header%kind )
      call tem_abort()
    end select

  end subroutine init_advRel
! **************************************************************************** !

! **************************************************************************** !
  !> Initialize flow field depends on read restart or initial condition
  !!
  subroutine mus_init_flow( scheme, tree, general, physics, scaling, &
    &                       levelPointer)
    ! --------------------------------------------------------------------------
    !> Scheme type
    type(mus_scheme_type), intent(inout) :: scheme
    type( tem_general_type ), intent(inout) :: general
    type( mus_physics_type ), intent(in) :: physics
    character(len=*), intent(in) :: scaling
    !> tree
    type(treelmesh_type), intent(in) :: tree
    !> global info type
    integer, intent(in) :: levelPointer(:)
    ! --------------------------------------------------------------------------
    integer :: iLevel, minLevel, maxLevel!, velPos(3), grad_pos(9), iElem
    ! --------------------------------------------------------------------------

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel


    ! Fill state vector from restart file
    if ( general%restart%controller%readRestart ) then
      !@todo: no general read option available in restart read
      call mus_readRestart( levelPointer = levelPointer,           &
        &                   restart      = general%restart,        &
        &                   scheme       = scheme,                 &
        &                   tree         = tree                    )
    else
      ! initialize state vector depends on scheme type
      call mus_init_byIC( scheme = scheme,       &
        &                 tree   = tree,         &
        &                 scaling= scaling,      &
        &                 fac    = physics%fac  )
    end if !restart?

    ! init auxiliary field variable from state for fluid elements in state and
    ! interpolate for ghostFromFiner elements in do_intp routine.
    call mus_initAuxField(scheme, general, minLevel, maxLevel)

    ! Fill all elements (ghost, halo) with valid values from fluid elements
    call fillHelperElementsFineToCoarse( scheme     = scheme,   &
      &                                  general    = general,  &
      &                                  physics    = physics,  &
      &                                  iLevel     = minLevel, &
      &                                  maxLevel   = maxLevel  )

    call fillHelperElementsCoarseToFine( scheme     = scheme,   &
      &                                  general    = general,  &
      &                                  physics    = physics,  &
      &                                  iLevel     = minLevel, &
      &                                  minLevel   = minLevel, &
      &                                  maxLevel   = maxLevel  )

    if ( scheme%intp%config%testInterpolation ) then
      do iLevel = minLevel, maxLevel
        call mus_testInterpolation( scheme     = scheme,                &
          &                         tree       = tree,                  &
          &                         general    = general,               &
          &                         fac        = physics%fac(iLevel),   &
          &                         iLevel     = iLevel,                &
          &                         minLevel   = minLevel,              &
          &                         maxLevel   = maxLevel,              &
          &                         pdf        = scheme%pdf             )
      end do
    end if

    ! Choose the advection relaxation scheme
    call init_advRel( scheme = scheme )

  end subroutine mus_init_flow
! **************************************************************************** !


! **************************************************************************** !
  !> Initialize flow field by calling corresponding routine according to scheme
  !! kind.
  !!
  subroutine mus_init_byIC( scheme, tree, fac, scaling )
    ! --------------------------------------------------------------------------
    !> Scheme type
    type(mus_scheme_type), intent(inout) :: scheme
    character(len=*), intent(in) :: scaling
    !> tree
    type(treelmesh_type), intent(in) :: tree
    !> Global parameters
    type(mus_convertFac_type), intent(in) :: fac(tree%global%minLevel &
      &                                          :tree%global%maxLevel)
    ! --------------------------------------------------------------------------
    integer :: nElems, nSize, iField, iLevel, minLevel, maxLevel
    ! --------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = logUnit(1))

    minLevel = tree%global%minLevel
    maxLevel = tree%global%maxLevel

    write(logUnit(1),"(A)") "Initialize flow"

    do iLevel = minLevel, maxlevel

      write(logUnit(7),"(A,I0)") '  on level: ', iLevel

      nElems = scheme%pdf( iLevel )%nElems_fluid
      nSize  = scheme%pdf( iLevel )%nSize

      select case (trim(scheme%header%kind))
      case ('fluid', 'fluid_incompressible',        &
           & 'fluid_GNS', 'fluid_incompressible_GNS')
        do iField = 1, scheme%nFields
          call mus_init_pdf( me     = scheme,                    &
            &    tree   = tree,                                  &
            &    fac    = fac(iLevel),                           &
            &    scaling= scaling,                               &
            &    nElems = nElems,                                &
            &    nSize  = nSize,                                 &
            &    state  = scheme%state(iLevel)%                  &
            &             val(:,scheme%pdf(iLevel)%nNext ),      &
            &    neigh  = scheme%pdf(iLevel)%neigh(:),           &
            &    iField = iField,                                &
            &    iLevel = iLevel,                                &
            &    field  = scheme%field(iField)                   )
        end do ! ifield
      case ('poisson','poisson_boltzmann_linear', 'poisson_boltzmann_nonlinear')
        if (scheme%nFields== 1) then
          call mus_init_poisson(                            &
             &    me     = scheme,                          &
             &    tree   = tree,                            &
             &    fac    = fac(iLevel),                     &
             &    scaling= scaling,                         &
             &    nElems = nElems,                          &
             &    nSize  = nSize,                           &
             &    state  = scheme%state(iLevel)%            &
             &             val(:,scheme%pdf(iLevel)%nNext), &
             &    neigh  = scheme%pdf(iLevel)%neigh(:),     &
             &    iLevel = iLevel,                          &
             &    field  = scheme%field(1)                  )
         else
           call tem_abort('nFields>1 for poisson')
        end if
      case ('nernst_planck')
        do iField = 1, scheme%nFields
          call mus_init_nernst_planck(                           &
            &    me           = scheme,                          &
            &    tree         = tree,                            &
            &    fac          = fac(iLevel),                     &
            &    nElems       = nElems,                          &
            &    nSize        = nSize,                           &
            &    state        = scheme%state(iLevel)%            &
            &                   val(:,scheme%pdf(iLevel)%nNext), &
            &    neigh        = scheme%pdf(iLevel)%neigh(:),     &
            &    iField       = iField,                          &
            &    iLevel       = iLevel,                          &
            &    nernstPlanck = scheme%nernstPlanck,             &
            &    field        = scheme%field(iField)             )
        end do
      case ('passive_scalar')
        ! \todo KM: 20161206 Implement compute kernel for multifield
        !                    passive scalar
        do iField = 1, scheme%nFields
          call mus_init_passiveScalar(                          &
            &    me     = scheme,                               &
            &    tree   = tree,                                 &
            &    fac    = fac(iLevel),                          &
            &    scaling= scaling,                              &
            &    nElems = nElems,                               &
            &    nSize  = nSize,                                &
            &    state  = scheme%state(iLevel)%                 &
            &             val(:,scheme%pdf(iLevel)%nNext),      &
            &    neigh  = scheme%pdf(iLevel)%neigh(:),          &
            &    iField = iField,                               &
            &    iLevel = iLevel,                               &
            &    field  = scheme%field(iField)                  )
        end do
      case ('multispecies_liquid')
        call mus_init_MSLiquid(                                 &
          &    me      = scheme,                                &
          &    tree    = tree,                                  &
          &    fac     = fac(iLevel),                           &
          &    nElems  = nElems,                                &
          &    nSize   = nSize,                                 &
          &    state   = scheme%state(iLevel)%                  &
          &              val(:,scheme%pdf(iLevel)%nNext ),      &
          &    neigh   = scheme%pdf(iLevel)%neigh(:),           &
          &    mixture = scheme%mixture,                        &
          &    iLevel  = iLevel,                                &
          &    field   = scheme%field                           )
      case ('multispecies_gas')
        call mus_init_MSGas(                                   &
          &    me     = scheme,                                &
          &    tree   = tree,                                  &
          &    fac    = fac(iLevel),                           &
          &    nElems = nElems,                                &
          &    nSize  = nSize,                                 &
          &    state  = scheme%state(iLevel)%                  &
          &             val(:,scheme%pdf(iLevel)%nNext ),      &
          &    neigh  = scheme%pdf(iLevel)%neigh(:),           &
          &    iLevel = iLevel,                                &
          &    field  = scheme%field                           )
      case ('isotherm_acEq')
        do iField = 1, scheme%nFields
          call mus_init_isotherm_acEq( me     = scheme,          &
            &    tree   = tree,                                  &
            &    fac    = fac(iLevel),                           &
            &    nElems = nElems,                                &
            &    nSize  = nSize,                                 &
            &    state  = scheme%state(iLevel)%                  &
            &             val(:,scheme%pdf(iLevel)%nNext ),      &
            &    neigh  = scheme%pdf(iLevel)%neigh(:),           &
            &    iField = iField,                                &
            &    iLevel = iLevel,                                &
            &    field  = scheme%field(iField)                   )
        end do ! ifield
      case default
        write(logUnit(1),"(A)") ' Scheme kind '//trim(scheme%header%kind)&
          &                     //'is wrong! Can NOT do IC!'
        call tem_abort()
      end select
    end do !iLevel

  end subroutine mus_init_byIC
! **************************************************************************** !


! **************************************************************************** !
  !> Initialize the flow from pressure, velocity and strain rate.\n
  !! First equilibirium pdf (fEq) is calculated from pressure and velocity.
  !! Then non-equilibirium (fnEq) is calculated from strain rate.
  !! At last set the pdf of each element by sum up these two parts (fEq+fnEq).
  !!
  subroutine mus_init_pdf(me, tree, fac, scaling, Field, iField, state, neigh, &
    &                     nElems, nSize, iLevel )
    ! --------------------------------------------------------------------------
    !> Scheme type
    type(mus_scheme_type), intent(in) :: me
    !> Global parameters
    type( mus_convertFac_type ), intent(in) :: fac
    !> scaling
    character(len=*), intent(in) :: scaling
    !> tree type
    type( treelmesh_type ), intent(in) :: tree
    !> Field type
    type(mus_field_type), intent(inout) :: field
    !> Field index
    integer, intent(in)                 :: iField
    !> Number of local elements
    integer, intent(in)                 :: nElems
    !> number of elements as size
    integer, intent(in)                 :: nSize
    !> Level index
    integer, intent(in)                 :: iLevel
    !> PDF
    real(kind=rk), intent(inout)        :: state(:)
    !> Connectivity array
    integer, intent(in) :: neigh(:)
    ! --------------------------------------------------------------------------
    integer :: iDir, iElem
    real(kind=rk), allocatable :: fEq(:), fnEq(:), rho(:)
    real(kind=rk), allocatable :: xc(:,:), vel(:,:)
    real(kind=rk), allocatable :: Sxx(:,:) ! Sxx, Syy, Szz, Sxy, Syz, Sxz
    integer :: iChunk, nChunks, chunkSize, nChunkElems, elemOff, elemPos, QQ
    integer :: offset! , nElems_local
    real(kind=rk) :: inv_p, inv_v, inv_s
    integer :: nScalars
    ! --------------------------------------------------------------------------

    ! when AOS, nSize is not used in this routine
    ! by this we can avoid compiler warning
    iDir = nSize
    QQ = me%layout%fStencil%QQ
    nScalars = me%varSys%nScalars

    inv_p = 1._rk / fac%press
    inv_v = 1._rk / fac%vel
    inv_s = 1._rk / fac%strainRate

    ! find chunksize and number of chunks required for initialzation
    chunkSize = io_buffer_size / QQ
    nChunks = ceiling( dble(nElems)/dble(chunkSize) )
    allocate(xc(chunkSize, 3))
    allocate(rho(  chunkSize))
    allocate(vel(3,chunkSize))
    allocate(Sxx(6,chunkSize))  ! Sxx, Syy, Szz, Sxy, Syz, Sxz
    ! use AOS layout for fEq and fnEq
    allocate(fEq( QQ*chunkSize))
    allocate(fnEq(QQ*chunkSize))

    do iChunk = 1, nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ( (iChunk-1)*chunksize )
      nChunkElems = min(chunkSize, nElems - elemOff)

      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        ! Calculate the coordinates
        xc(iElem,1:3) = tem_BaryOfId( tree,     &
          &                           me%levelDesc(iLevel)%total(elemPos))
      end do

      rho(1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(1),  &
        &                                   coord = xc(1:nChunkElems,1:3),  &
        &                                   n     = nChunkElems             )
      vel(1,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(2), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      vel(2,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(3), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      vel(3,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(4), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      ! Read in the shear rate tensor
      ! This corresponds to S = grad(u) + grad(u)^T
      Sxx(1,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(5), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      Sxx(2,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(6), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      Sxx(3,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(7), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      Sxx(4,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(8), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      Sxx(5,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(9), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      Sxx(6,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(10), &
        &                                     coord = xc(1:nChunkElems,1:3),  &
        &                                     n     = nChunkElems             )

      ! convert these quantities from physics to LB
!cdir nodep
!ibm* novector
!dir$ novector
      do iElem = 1, nChunkElems
        rho(iElem) = rho(iElem) * cs2inv * inv_p
        vel(1,iElem) = vel(1,iElem) * inv_v
        vel(2,iElem) = vel(2,iElem) * inv_v
        vel(3,iElem) = vel(3,iElem) * inv_v
        Sxx(1,iElem) = Sxx(1,iElem) * inv_s
        Sxx(2,iElem) = Sxx(2,iElem) * inv_s
        Sxx(3,iElem) = Sxx(3,iElem) * inv_s
        Sxx(4,iElem) = Sxx(4,iElem) * inv_s
        Sxx(5,iElem) = Sxx(5,iElem) * inv_s
        Sxx(6,iElem) = Sxx(6,iElem) * inv_s
      end do

      call me%derVarPos(iField)%equilFromMacro( &
        &    density  = rho(1:nChunkElems),     &
        &    velocity = vel(1:3,1:nChunkElems), &
        &    iField   = iField,                 &
        &    nElems   = nChunkElems,            &
        &    varSys   = me%varSys,              &
        &    layout   = me%layout,              &
        &    res      = fEq                     )

      select case ( trim(scaling) )
      case ('acoustic')
        ! Acoustic fNeq with defining the stress tensor instead of the
        ! shear rate tensor in the lua file
        do iElem = 1, nChunkElems
          fNeq( (iElem-1)*QQ+1 : iElem*QQ ) =                     &
            &   getnEq_acoustic(                                  &
            &     layout = me%layout,                             &
            &     omega  = field%fieldProp%fluid%viscKine         &
            &                  %omLvl(iLevel)%val(elemOff+iElem), &
            &     Sxx    = [ Sxx(1,iElem), Sxx(4,iElem),          &
            &                Sxx(6,iElem), Sxx(4,iElem),          &
            &                Sxx(2,iElem), Sxx(5,iElem),          &
            &                Sxx(6,iElem ), Sxx(5,iElem),         &
            &                Sxx(3,iElem) ] )
        end do ! iElem
      case ('diffusive')
        ! Diffusive non-equilibrium part with the shear rate tensor
        ! given as an input
        do iElem = 1, nChunkElems
          fNeq( (iElem-1)*QQ+1 : iElem*QQ ) =                   &
          &  getnEq_diffusive(                                  &
          &     layout = me%layout,                             &
          &     omega  = field%fieldProp%fluid%viscKine         &
          &                  %omLvl(iLevel)%val(elemOff+iElem), &
          &     Sxx    = [ Sxx(1,iElem), Sxx(4,iElem),          &
          &                Sxx(6,iElem), Sxx(4,iElem),          &
          &                Sxx(2,iElem), Sxx(5,iElem),          &
          &                Sxx(6,iElem), Sxx(5,iElem),          &
          &                Sxx(3,iElem) ] )
        end do ! iElem
      end select ! scaling

      ! assign pdf by fEq + fNeq
      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        offset = (iElem-1)*QQ
        do iDir = 1, QQ
          state( ( elempos-1)* nscalars+idir+( 1-1)* qq ) &
            & = fEq(offset + iDir) + fNeq(offset + iDir)
        end do ! iDir
      end do ! iElem

    end do ! chunk

    ! deallocate memory for next chunk
    deallocate(xc)
    deallocate(rho)
    deallocate(vel)
    deallocate(Sxx)
    deallocate(fEq)
    deallocate(fnEq)

  end subroutine mus_init_pdf
! **************************************************************************** !

!***************************************************************************** !
  !> Initialize passive scalar from pressure and velocity.\n
  !! Equilibirium pdf (fEq) is calculated from pressure and velocity.
  !!
  subroutine mus_init_passiveScalar(me, tree, fac, scaling, Field, iField, &
    &                               state, neigh, nElems, nSize, iLevel)
    ! --------------------------------------------------------------------------
    !> Scheme type
    type(mus_scheme_type), intent(in) :: me
    !> Global parameters
    type( mus_convertFac_type ), intent(in) :: fac
    !> scaling
    character(len=*), intent(in) :: scaling
    !> tree type
    type( treelmesh_type ), intent(in) :: tree
    !> Field type
    type(mus_field_type), intent(inout) :: field
    !> Field index
    integer, intent(in)                 :: iField
    !> Number of local elements
    integer, intent(in)                 :: nElems
    !> number of elements as size
    integer, intent(in)                 :: nSize
    !> Level index
    integer, intent(in)                 :: iLevel
    !> PDF
    real(kind=rk), intent(inout)        :: state(:)
    !> Connectivity array
    integer, intent(in) :: neigh(:)
    ! --------------------------------------------------------------------------
    integer :: iDir, iElem
    real(kind=rk), allocatable :: fEq(:), rho(:)
    real(kind=rk), allocatable :: xc(:,:), vel(:,:)
    integer :: iChunk, nChunks, chunkSize, nChunkElems, elemOff, elemPos, QQ
    integer :: offset
    real(kind=rk) :: inv_p, inv_v
    integer :: nScalars
    ! --------------------------------------------------------------------------

    ! when AOS, nSize is not used in this routine
    ! by this we can avoid compiler warning
    iDir = nSize
    QQ = me%layout%fStencil%QQ
    nScalars = me%varSys%nScalars

    inv_p = 1._rk / fac%press
    inv_v = 1._rk / fac%vel

    ! find chunksize and number of chunks required for initialzation
    chunkSize = io_buffer_size / QQ
    nChunks = ceiling( dble(nElems)/dble(chunkSize) )
    allocate(xc(chunkSize, 3))
    allocate(rho(  chunkSize))
    allocate(vel(3,chunkSize))
    ! use AOS layout for fEq and fnEq
    allocate(fEq( QQ*chunkSize))

    do iChunk = 1, nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ( (iChunk-1)*chunksize )
      nChunkElems = min(chunkSize, nElems - elemOff)

      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        ! Calculate the coordinates
        xc(iElem,1:3) = tem_BaryOfId( tree,     &
          &                           me%levelDesc(iLevel)%total(elemPos))
      end do

      rho(1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(1),  &
        &                                   coord = xc(1:nChunkElems,1:3),  &
        &                                   n     = nChunkElems             )
      vel(1,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(2), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      vel(2,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(3), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      vel(3,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(4), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )

      ! convert these quantities from physics to LB
!cdir nodep
!ibm* novector
!dir$ novector
      do iElem = 1, nChunkElems
        rho(iElem) = rho(iElem) * cs2inv *inv_p
        vel(1,iElem) = vel(1,iElem) * inv_v
        vel(2,iElem) = vel(2,iElem) * inv_v
        vel(3,iElem) = vel(3,iElem) * inv_v
      end do

      call me%derVarPos(iField)%equilFromMacro( &
        &    density  = rho(1:nChunkElems),     &
        &    velocity = vel(1:3,1:nChunkElems), &
        &    iField   = iField,                 &
        &    nElems   = nChunkElems,            &
        &    varSys   = me%varSys,              &
        &    layout   = me%layout,              &
        &    res      = fEq                     )

      ! fNeq = zero for passive_scalar
      ! assign pdf = fEq
      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        offset = (iElem-1)*QQ
        do iDir = 1, QQ
          state( ( elempos-1)* nscalars+idir+( 1-1)* qq ) &
            & = fEq(offset + iDir)
        end do ! iDir
      end do ! iElem

    end do ! chunk

    ! deallocate memory for next chunk
    deallocate(xc)
    deallocate(rho)
    deallocate(vel)
    deallocate(fEq)

  end subroutine mus_init_passiveScalar
! **************************************************************************** !


! ****************************************************************************!
  !> Initialize poisson lbm from potential
  !! Equilibirium pdf (fEq) is calculated from potential.
  !!
  subroutine mus_init_poisson(me, tree, fac, scaling, Field, state, neigh, &
    &                         nElems, nSize, iLevel)
    ! -------------------------------------------------------------------------
    !> Scheme type
    type(mus_scheme_type), intent(in) :: me
    !> Global parameters
    type( mus_convertFac_type ), intent(in) :: fac
    !> scaling
    character(len=*), intent(in) :: scaling
    !> tree type
    type( treelmesh_type ), intent(in) :: tree
    !> Field type
    type(mus_field_type), intent(inout) :: field
    !> Number of local elements
    integer, intent(in)                 :: nElems
    !> number of elements as size
    integer, intent(in)                 :: nSize
    !> Level index
    integer, intent(in)                 :: iLevel
    !> PDF
    real(kind=rk), intent(inout)        :: state(:)
    !> Connectivity array
    integer, intent(in) :: neigh(:)
    ! -------------------------------------------------------------------------
    integer :: iDir, iElem
    real(kind=rk), allocatable :: potential(:)
    real(kind=rk), allocatable :: xc(:,:)
    integer :: iChunk, nChunks, chunkSize, nChunkElems, elemOff, elemPos, QQ
    real(kind=rk) :: inv_potential, fEq
    integer :: nScalars
    ! --------------------------------------------------------------------------

    ! when AOS, nSize is not used in this routine
    ! by this we can avoid compiler warning
    iDir = nSize
    QQ = me%layout%fStencil%QQ
    nScalars = me%varSys%nScalars


    inv_potential = 1._rk/ fac%potential

    ! find chunksize and number of chunks required for initialzation
    chunkSize = io_buffer_size / QQ
    nChunks = ceiling( dble(nElems)/dble(chunkSize) )
    allocate(xc(chunkSize, 3))
    allocate(potential(chunkSize))

    do iChunk = 1, nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ( (iChunk-1)*chunksize )
      nChunkElems = min(chunkSize, nElems - elemOff)

      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        ! Calculate the coordinates
        xc(iElem,1:3) = tem_BaryOfId( tree,                              &
          &                           me%levelDesc(iLevel)%total(elemPos))
      end do

      potential(1:nChunkElems)                              &
        & = tem_spatial_for( me    = field%ic%ini_state(1), &
        &                    coord = xc(1:nChunkElems,1:3), &
        &                    n     = nChunkElems            )

      ! convert these quantities from physics to LB
      potential = potential * inv_potential

      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        do iDir = 1,QQ
          !> Calculate equilibrium distribution functions fEq
          fEq = me%layout%weight(iDir)*potential(iElem)

          ! assign pdf = fEq
          state( ( elempos-1)* nscalars+idir+( 1-1)* qq ) &
            & = fEq

        end do ! iDir
      end do !iElem
    end do !chunk

    ! deallocate memory for next chunk
    deallocate(xc)
    deallocate(potential)

  end subroutine mus_init_poisson
 ! **************************************************************************!


! **************************************************************************** !
  !> Initialize nernst planck from           and           .\n
  !! Equilibirium pdf (fEq) is calculated from          and         .
  !!
  subroutine mus_init_nernst_planck(me, tree, fac, Field, iField, state,   &
    &                               neigh, nElems, nSize, iLevel, nernstPlanck)
    ! --------------------------------------------------------------------------
    !> Scheme type
    type(mus_scheme_type), intent(in) :: me
    !> Global parameters
    type( mus_convertFac_type ), intent(in) :: fac
    !> tree type
    type( treelmesh_type ), intent(in) :: tree
    !> Field type
    type(mus_field_type), intent(inout) :: field
    !> Field index
    integer, intent(in)                 :: iField
    !> Number of local elements
    integer, intent(in)                 :: nElems
    !> number of elements as size
    integer, intent(in)                 :: nSize
    !> Level index
    integer, intent(in)                 :: iLevel
    !> PDF
    real(kind=rk), intent(inout)        :: state(:)
    !> Connectivity array
    integer, intent(in) :: neigh(:)
    !> Contins solvent information
    type(mus_nernstPlanck_type), intent(in) :: nernstPlanck
    ! --------------------------------------------------------------------------
    integer :: iDir, iElem
    real(kind=rk), allocatable :: moleDens(:)
    real(kind=rk), allocatable :: xc(:,:), vel(:,:)
    integer :: iChunk, nChunks, chunkSize, nChunkElems, elemOff, elemPos, QQ
    integer :: offset
    real(kind=rk) :: inv_v, ucx
    integer :: nScalars
    ! --------------------------------------------------------------------------

    ! when AOS, nSize is not used in this routine
    ! by this we can avoid compiler warning
    iDir = nSize
    QQ = me%layout%fStencil%QQ
    nScalars = me%varSys%nScalars

    inv_v = 1._rk / fac%vel

    ! find chunksize and number of chunks required for initialzation
    chunkSize = io_buffer_size / QQ
    nChunks = ceiling( dble(nElems)/dble(chunkSize) )
    allocate(xc(chunkSize, 3))
    allocate(moleDens(chunkSize))
    allocate(vel(3,chunkSize))

    do iChunk = 1, nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ( (iChunk-1)*chunksize )
      nChunkElems = min(chunkSize, nElems - elemOff)

      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        ! Calculate the coordinates
        xc(iElem,1:3) = tem_BaryOfId( tree,     &
          &                           me%levelDesc(iLevel)%total(elemPos))
      end do

      ! mole fraction
      moleDens(1:nChunkElems) = tem_spatial_for(                  &
        &                         me    = field%ic%ini_state(1),  &
        &                         coord = xc(1:nChunkElems,1:3),  &
        &                         n     = nChunkElems             )
      vel(1,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(2), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      vel(2,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(3), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      vel(3,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(4), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )

      ! convert these quantities from physics to LB
!cdir nodep
!ibm* novector
!dir$ novector
      do iElem = 1, nChunkElems
        moleDens(iElem) = moleDens(iElem) * nernstPlanck%moleDens
        vel(1,iElem) = vel(1,iElem) * inv_v
        vel(2,iElem) = vel(2,iElem) * inv_v
        vel(3,iElem) = vel(3,iElem) * inv_v
      end do

      ! fNeq = zero for passive_scalar
      ! assign pdf = fEq
      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        offset = (iElem-1)*QQ
        do iDir = 1, QQ
          ucx = dot_product(me%layout%fStencil%cxDir(:, iDir), &
            &               vel(:,iElem))
          state( ( elempos-1)* nscalars+idir+( 1-1)* qq ) &
            & = me%layout%weight(iDir)*moleDens(iElem)*(1.0_rk + ucx*cs2inv)
        end do ! iDir
      end do ! iElem

    end do ! chunk

    ! deallocate memory for next chunk
    deallocate(xc)
    deallocate(moleDens)
    deallocate(vel)

  end subroutine mus_init_nernst_planck
! **************************************************************************** !

! **************************************************************************** !
  !> Initialize the flow from calculated quantitites like density, velocity etc.
  !! for multispecies lbm
  subroutine mus_init_MSLiquid( me, tree, fac, state, neigh, &
    &                           Field, mixture, nElems, nSize, iLevel )
    ! --------------------------------------------------------------------------
    type(mus_scheme_type), intent(inout) :: me !< Scheme type
    type(mus_convertFac_type), intent(in) :: fac  !< Global parameters
    type(treelmesh_type), intent(in) :: tree !< tree
    type(mus_field_type), intent(inout) :: field(:) !< Field type
    type(mus_mixture_type), intent(inout) :: mixture !< mixture type
    integer, intent(in)                 :: nElems !< Number of elements
    integer, intent(in)                 :: nSize  !< Number of elements as size
    integer, intent(in)                 :: iLevel !< Level index
    !> PDF
    real(kind=rk), intent(inout)        :: state(:)
    !> Connectivity array
    integer, intent(in) :: neigh(:)
    ! --------------------------------------------------------------------------
    integer :: iDir, iElem, iField, nFields, QQ
    real(kind=rk), allocatable :: fEqStar(:), fEq(:), press(:)
    real(kind=rk), allocatable :: xc(:,:), ux(:,:), uy(:,:), uz(:,:)
    real(kind=rk), allocatable :: rho(:,:), moleFrac(:,:)
    real(kind=rk) :: phi(me%nFields), tot_massDens
    integer :: iChunk, nChunks, chunkSize, nChunkElems, elemPos, elemOff
    real(kind=rk), dimension(3) :: eqVel, velAvg, velQuadStar, velQuad
    real(kind=rk) :: vel(3, me%nFields)
    integer :: iField_2, iField_3
    real(kind=rk) :: resi_coeff(me%nFields, me%nFields)
    real(kind=rk) :: diff_coeff(me%nFields, me%nFields)
    real(kind=rk) :: ucx, ucxStar, usq, usqStar, ucxQuad, ucxQuadStar
    real(kind=rk), dimension(me%nFields, me%nFields) :: &
      & thermodynamic_fac, inv_thermodyn_fac
    integer :: restPosition
    integer :: nScalars
    ! --------------------------------------------------------------------------

    nScalars = me%varSys%nScalars
    QQ = me%layout%fStencil%QQ
    nFields = me%nFields
    ! molecular weight ratios
    phi = field(:)%fieldProp%species%molWeigRatio
    ! pdf stencil at rest (center)
    restPosition = me%layout%fStencil%restPosition

    !resistivities
    do iField = 1, nFields
      resi_coeff(iField,:) = field(iField)%fieldProp%species%resi_coeff
    end do

    ! find chunksize and number of chunks required for initialzation
    chunkSize = io_buffer_size/QQ
    nChunks = ceiling( real(nElems, kind=rk)/real(chunkSize, kind=rk))

    do iChunk = 1, nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ( (iChunk-1)*chunksize )
      nChunkElems = min(chunkSize, nElems - elemOff)

      !allocate memory for density, vel, eq, coord
      allocate(xc(nChunkElems, 3))
      allocate(moleFrac(nFields, nChunkElems))
      allocate(rho(nFields, nChunkElems))
      allocate(ux(nFields, nChunkElems))
      allocate(uy(nFields, nChunkElems))
      allocate(uz(nFields, nChunkElems))
      allocate(Press(nChunkElems))
      ! Equilibrium part of the pdf
      ! use AOS layout for fEqStart and fEq
      allocate(fEqStar(QQ*nChunkElems))
      allocate(fEq(    QQ*nChunkElems))

      do iElem = 1, nChunkElems
        ! Calculate the coordinates
        elemPos = elemOff + iElem
        xc(iElem,1:3) = tem_BaryOfId( tree, &
          & me%levelDesc( iLevel )%total( elemPos ) )
      end do

      do iField = 1, nFields
        ! read velocity
        ux(iField, :) = tem_spatial_for(                         &
          &               me    = field(iField)%ic%ini_state(2), &
          &               coord = xc,                            &
          &               n     = nChunkElems                    )
        ux(iField, :) = ux(iField, :)/fac%vel

        uy(iField, :) = tem_spatial_for(                         &
          &               me    = field(iField)%ic%ini_state(3), &
          &               coord = xc,                            &
          &               n     = nChunkElems                    )
        uy(iField, :) = uy(iField, :)/fac%vel

        uz(iField, :) = tem_spatial_for(                         &
          &               me    = field(iField)%ic%ini_state(4), &
          &               coord = xc,                            &
          &               n     = nChunkElems                    )
        uz(iField, :) = uz(iField, :)/fac%vel

        ! mole fraction of each species
        moleFrac(iField, :) = tem_spatial_for(                         &
          &                     me    = field(iField)%ic%ini_state(1), &
          &                     coord = xc,                            &
          &                     n     = nChunkElems                    )
      end do

      ! read mixture hydrodynamic pressure
      Press = tem_spatial_for( me    = mixture%ic%ini_state(1), &
        &                      coord = xc,                      &
        &                      n     = nChunkElems              )

      !convert to lattice unit
      press = press/fac%press

      !compute species density from moleFraction and mixture number density
      !\rho_i = n0 * \chi_i * m_i + mass_fraction_i*rho0*KinePress/(cs^2*phi_i)
      !\rho_i = n0 * \chi_i * m_i + mass_fraction_i*rho0*Press/rho/(cs^2*phi_i)
      !\rho_i = n0 * \chi_i * m_i + mass_fraction_i*Press/(cs^2*phi_i)
      do iElem = 1, nChunkElems
        ! check if physical constraint sum(MolFrac) = 1
        ! check upto single precision
        if( abs( 1.0_single_k - real( sum(moleFrac(:,iElem)), kind=single_k ) )&
          & > eps_single_k ) then
          write(logUnit(1),*)'Error: Initial sum(molefraction) of all fields is'
          write(logUnit(1),*)'       not equal to 1 for iElem', iElem
          write(logUnit(1),*)' totMolFrac:', sum(moleFrac(:,iElem))
          write(logUnit(1),*)' MoleFracs: ', moleFrac(:,iElem)
          call tem_abort()
        end if
        !mixMolWeight = sum( chi_i * m_i )
        !mixMolWeight = sum( moleFrac(:, iElem)                                &
        !  &        * field(:)%fieldProp%species%molWeight )
        do iField = 1, nFields
!         rho(iField, iElem) = mixture%rho0LB * moleFrac(iField, iElem)        &
!           & * field(iField)%fieldProp%species%molWeight / mixMolWeight
          rho(iField, iElem) = mixture%moleDens0LB * moleFrac(iField, iElem)   &
            & * field(iField)%fieldProp%species%molWeight
        end do
        ! add kinematic pressure term to density
        tot_massDens = sum(rho(:,iElem))
        do iField = 1, nFields
          rho(iField, iElem) = rho(iField, iElem) &
            & + ( cs2inv*(rho(iField,iElem)/tot_massDens)*Press(iElem)  &
            & / phi(iField) )
        end do

      end do !iElem

      do iField = 1, nFields
        ! Calculate the equilibrium distribution
        fEqStar = 0.0_rk
        fEq = 0.0_rk

        do iElem = 1, nChunkElems
          do iField_2 = 1, nFields
           vel(1, iField_2) = ux(iField_2, iElem)
           vel(2, iField_2) = uy(iField_2, iElem)
           vel(3, iField_2) = uz(iField_2, iElem)
          end do
          ! if scheme relaxation is bgk_withthermodynfac then
          ! calculate thermodynamic factor from c++ code and
          ! compute equilibrium velocity from thermodynamic factor
          if (trim(me%header%relaxation) == 'bgk_withthermodynfac' .or. &
            & trim(me%header%relaxation) == 'mrt_withthermodynfac') then

            call mus_calc_MS_DiffMatrix(nFields, me%mixture%temp0,             &
              &                        me%mixture%atm_press,                   &
              &                        moleFrac(:, iElem)*me%mixture%moleDens0,&
              &                        diff_coeff )
            call mus_calc_thermFactor( nFields, me%mixture%temp0,              &
              &                        me%mixture%atm_press,                   &
              &                        moleFrac(:, iElem),                     &
              &                        thermodynamic_fac )

            ! calculate resi_coeff from diff_coeff from C-code
            resi_coeff = fac%diffusivity/diff_coeff
            do iField_2 = 1, nFields
              resi_coeff(iField_2, iField_2) = mixture%paramB
            end do

            inv_thermodyn_fac = invert_matrix( thermodynamic_fac )

            eqVel = rho(iField, iElem)*vel(:, iField)
            do iField_2 = 1, nFields
              do iField_3 = 1, nFields
                eqVel(:) = eqVel(:)                                            &
                  &      + inv_thermodyn_fac(iField, iField_2)                 &
                  &      * rho(iField_2, iElem)                                &
                  &      * resi_coeff(iField_2, iField_3)                      &
                  &      * phi(iField_2) * moleFrac(iField_3, iElem)           &
                  &      * ( vel(:, iField_3) - vel( :, iField_2 ) )           &
                  &      / me%mixture%paramB
              end do
            end do
          else
            eqVel = vel(:, iField)
            do iField_2 = 1, nFields
              eqVel(:) = eqVel(:) + resi_coeff(iField, iField_2)               &
                &      * phi(iField) * moleFrac(iField_2, iElem)               &
                &      * ( vel(:, iField_2) - vel( :, iField ) )               &
                &      / me%mixture%paramB
            end do
            eqVel = rho(iField, iElem)*eqVel
          end if !Thermodynamic factor

          tot_massDens = sum(rho(:,iElem))
          ! mass averaged mixture velocity
          velAvg(1) = dot_product( rho(:, iElem), vel(1, :) )/tot_massDens
          velAvg(2) = dot_product( rho(:, iElem), vel(2, :) )/tot_massDens
          velAvg(3) = dot_product( rho(:, iElem), vel(3, :) )/tot_massDens

          ! velocity in quadratic term of equilibrium
          ! eqVel is div by rho since rho is multiplied in the
          ! calculate of eqVel with thermodynamic factor
          ! due to this reason Bilinear part in fEq is not multiplied
          ! by rho in fEq calculation.
          velQuadStar(:) = me%mixture%theta_eq*velAvg(:)                       &
            &            + (1.0_rk-me%mixture%theta_eq)                        &
            &            * ( eqVel(:) / rho(iField, iElem) )
          velQuad(:) = me%mixture%theta_eq*velAvg(:)                           &
            &            + (1.0_rk-me%mixture%theta_eq)*vel(:, iField)

          ! compute Eq with eqVel and initial velocity
          usqStar = dot_product(velQuadStar, velQuadStar)*t2cs2inv
          usq = dot_product(velQuad, velQuad)*t2cs2inv

          do iDir =1,QQ
            ucx = dot_product(me%layout%fStencil%cxDir(:, iDir),             &
              &               vel(:,iField))
            ucxQuad = dot_product(me%layout%fStencil%cxDir(:, iDir),         &
              &               velQuad)

            ucxStar = dot_product(me%layout%fStencil%cxDir(:, iDir),         &
              &                   eqVel)
            ucxQuadStar = dot_product(me%layout%fStencil%cxDir(:, iDir),     &
              &                   velQuadStar)

            ! eqVel is actually is rho_i*eqVel so ucxStar is not multiplied
            ! with rho in below equation
            fEqStar( (iElem-1)*QQ+iDir ) &
              & = me%layout%weight(iDir) * ( rho(iField, iElem) * ( phi(iField)&
              & + ucxQuadStar * ucxQuadStar * t2cs4inv - usqStar )             &
              & + ucxStar * cs2inv )

            fEq( (iElem-1)*QQ+iDir ) &
              & = me%layout%weight(iDir) * rho(iField, iElem)                  &
              & * ( phi(iField) + ucx * cs2inv + ucxQuad * ucxQuad * t2cs4inv  &
              & - usq )
          end do ! iDir

          ! equilibrium at rest
          select case( trim(me%header%layout) )
          case('d2q9')
            fEqStar( (iElem-1)*QQ+restPosition ) &
              & = me%layout%weight(restPosition) * rho(iField, iElem)          &
              & * ( (9._rk - 5._rk*phi(iField)) / 4.0_rk - usqStar )
            fEq( (iElem-1)*QQ+restPosition ) &
              & = me%layout%weight(restPosition) * rho(iField, iElem)          &
              & * ( (9._rk - 5._rk*phi(iField)) /4.0_rk - usq )
          case('d3q19')
            fEqStar((iElem-1)*QQ+restPosition ) &
              & = me%layout%weight(restPosition) * rho(iField, iElem)          &
              & * ( (3._rk - 2._rk*phi(iField)) - usqStar )
            fEq((iElem-1)*QQ+restPosition ) &
              & = me%layout%weight(restPosition) * rho(iField, iElem)          &
              & * ( (3._rk - 2._rk*phi(iField)) - usq )
          case default
            write(logUnit(1),*)'ERROR: initializing Multispecies LBM '
            write(logUnit(1),*)'Unknown stencil layout'
            call tem_abort()
          end select

        end do !iElem

        ! set initial pdf
        do iElem = 1, nChunkElems
          elemPos = elemOff + iElem
          do iDir = 1, QQ
            ! In general we want to set for each element a valid entry.
            ! This was  IDX, but for PULL, it doesnt change anything
            ! for PUSH however, we need to set it to save!
            state( ( elempos-1)* nscalars+idir+( ifield-1)* qq ) &
              & =   fEq((iElem-1)*QQ+iDir) &
              &   + me%mixture%omega_diff*0.5_rk                         &
              &   * ( fEq((iElem-1)*QQ+iDir) - fEqStar((iElem-1)*QQ+iDir))
          end do !idir
        end do ! iElem
      end do ! iField

      ! deallocate memory for next chunk
      deallocate(xc)
      deallocate(rho)
      deallocate(moleFrac)
      deallocate(ux)
      deallocate(uy)
      deallocate(uz)
      deallocate(press)
      deallocate(fEqStar)
      deallocate(fEq)

    end do ! chunk

  end subroutine mus_init_MSLiquid
! **************************************************************************** !

! **************************************************************************** !
  !> Initialize the flow from calculated quantitites like density, velocity etc.
  !! for multispecies lbm
  subroutine mus_init_MSGas( me, tree, fac, state, neigh, Field, nElems, &
    &                        nSize, iLevel )
    ! --------------------------------------------------------------------------
    type(mus_scheme_type), intent(inout) :: me !< Scheme type
    type(mus_convertFac_type), intent(in) :: fac  !< Global parameters
    type(treelmesh_type), intent(in) :: tree
    type(mus_field_type), intent(inout) :: field(:) !< Field type
    integer, intent(in)                 :: nElems !< Number of elements
    integer, intent(in)                 :: nSize  !< Number of elements as size
    integer, intent(in)                 :: iLevel !< Level index
    !> PDF
    real(kind=rk), intent(inout)        :: state(:)
    !> Connectivity array
    integer, intent(in) :: neigh(:)
    ! --------------------------------------------------------------------------
    integer :: iDir, iElem, iField, nFields, offset, QQ
    real(kind=rk), allocatable :: fEq(:)
    real(kind=rk), allocatable :: xc(:,:), ux(:,:), uy(:,:), uz(:,:)
    real(kind=rk), allocatable :: rho(:,:), rhoAll(:), velAll(:,:)
    real(kind=rk) :: phi(me%nFields)
    integer :: iChunk, nChunks, chunkSize, nChunkElems, elemPos, elemOff
    integer :: nScalars
    ! --------------------------------------------------------------------------

    nScalars = me%varSys%nScalars
    QQ = me%layout%fStencil%QQ
    nFields = me%nFields
    ! molecular weight ratios
    phi = field(:)%fieldProp%species%molWeigRatio
!write(dbgUnit,*) 'theta_eq ', theta_eq

    ! find chunksize and number of chunks required for initialzation
    chunkSize = io_buffer_size/QQ
    nChunks = ceiling( real(nElems, kind=rk)/real(chunkSize, kind=rk))

    do iChunk = 1, nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ( (iChunk-1)*chunksize )
      nChunkElems = min(chunkSize, nElems - elemOff)

      !allocate memory for density, vel, eq, coord
      allocate(xc(nChunkElems, 3))
      allocate(rho(nFields, nChunkElems))
      allocate(rhoAll(nFields*nChunkElems))
      allocate(ux(nFields, nChunkElems))
      allocate(uy(nFields, nChunkElems))
      allocate(uz(nFields, nChunkElems))
      allocate(velAll(3, nFields*nChunkElems))
      ! use AOS layout for fEq
      allocate(fEq(QQ*nChunkElems))

      do iElem = 1, nChunkElems
        ! Calculate the coordinates
        elemPos = elemOff + iElem
        xc(iElem,1:3) = tem_BaryOfId( tree, &
          & me%levelDesc( iLevel )%total( elemPos ) )
      end do

      do iField = 1, nFields
        !partial pressure
        rho(iField, :) = tem_spatial_for(                         &
          &                me    = field(iField)%ic%ini_state(1), &
          &                coord = xc,                            &
          &                n     = nChunkElems                    )
        ! read velocity
        ux(iField, :) = tem_spatial_for( me = field(iField)%ic%ini_state(2), &
          &                              coord = xc,                         &
          &                              n     = nChunkElems                 )
        ux(iField, :) = ux(iField, :)/fac%vel

        uy(iField, :) = tem_spatial_for( me = field(iField)%ic%ini_state(3), &
          &                              coord = xc,                         &
          &                              n     = nChunkElems                 )
        uy(iField, :) = uy(iField, :)/fac%vel

        uz(iField, :) = tem_spatial_for( me = field(iField)%ic%ini_state(4), &
          &                              coord = xc,                         &
          &                              n     = nChunkElems                 )
        uz(iField, :) = uz(iField, :)/fac%vel
      end do

      !partial pressure
      do iField = 1, nFields
        !p_sigma = r_sigma * cs2 * phi_sigma
        !rho = press * cs2inv / phi
        rho(iField, :) = ( rho(iField, :) * cs2inv / phi(iField) ) &
          &            / fac%press
      end do

!write(dbgUnit,*) 'iField ', iField
      !linearized the input for equilibrium function pointer
      !with following order
      !|dens_1 | dens_2 | .. | dens_nElems |
      !!| ux_1 | uy_1 | uz_1 | ..
      !!.. | ux_nelems| uy_nelems| uz_nelems|
      do iElem = 1, nChunkElems
        do iField = 1, nFields
          ! MH: Fix for the intel compiler:
          ! entries have to be assigned one by one instead of
          ! assigning them as an array of (/ux, uy, uz/)
          offset = (iElem-1)*nFields + iField
          rhoAll(offset) = rho(iField, iElem)
          velAll(:,offset) = (/ ux( iField, iElem ), &
            &                   uy( iField, iElem ), &
            &                   uz( iField, iElem ) /)
        end do
      end do

      do iField = 1, nFields
        ! Calculate the equilibrium distribution
        fEq = 0.0_rk
        call me%derVarPos(iField)%equilFromMacro( density  = rhoAll,      &
          &                                       velocity = velAll,      &
          &                                       iField   = iField,      &
          &                                       nElems   = nChunkElems, &
          &                                       varSys   = me%varSys,   &
          &                                       layout   = me%layout,   &
          &                                       res      = fEq          )

        do iElem = 1, nChunkElems
          elemPos = elemOff + iElem
          do iDir = 1, QQ
            state( ( elempos-1)* nscalars+idir+( ifield-1)* qq ) &
              & = fEq((iElem-1)*QQ + iDir)
          end do ! iDir = 1, QQ
        end do ! iElem = 1, nChunkElems
      end do ! iField

      ! deallocate memory for next chunk
      deallocate(xc)
      deallocate(rho)
      deallocate(rhoAll)
      deallocate(ux)
      deallocate(uy)
      deallocate(uz)
      deallocate(velAll)
      deallocate(fEq)

    end do ! chunk

  end subroutine mus_init_MSGas
! **************************************************************************** !




! **************************************************************************** !
  !> Initialize the isothermal acEq flow from density and velocity\n
  !! equilibrium pdf (fEq) is calculated from density and velocity
  !!
  subroutine mus_init_isotherm_acEq(me, tree, fac, Field, iField, state, &
    &                               neigh, nElems, nSize, iLevel )
    ! --------------------------------------------------------------------------
    !> Scheme type
    type(mus_scheme_type), intent(inout) :: me
    !> Global parameters
    type(mus_convertFac_type), intent(in) :: fac
    !> tree type
    type( treelmesh_type ), intent(in) :: tree
    !> Field type
    type(mus_field_type), intent(inout) :: field
    !> Field index
    integer, intent(in)                 :: iField
    !> Number of local elements
    integer, intent(in)                 :: nElems
    !> number of elements as size
    integer, intent(in)                 :: nSize
    !> Level index
    integer, intent(in)                 :: iLevel
    !> PDF
    real(kind=rk), intent(inout)        :: state(:)
    !> Connectivity array
    integer, intent(in) :: neigh(:)
    ! --------------------------------------------------------------------------
    integer :: iDir, iElem
    real(kind=rk), allocatable :: fEq(:), rho(:)
    real(kind=rk), allocatable :: xc(:,:), vel(:,:)
    integer :: iChunk, nChunks, chunkSize, nChunkElems, elemOff, elemPos, QQ
    integer :: offset
    real(kind=rk) :: inv_p, inv_v, inv_s
    integer :: nScalars
    ! --------------------------------------------------------------------------

    ! when AOS, nSize is not used in this routine
    ! by this we can avoid compiler warning
    iDir = nSize
    QQ = me%layout%fStencil%QQ
    nScalars = me%varSys%nScalars

    inv_p = 1._rk / fac%press
    inv_v = 1._rk / fac%vel
    inv_s = 1._rk / fac%strainRate

    ! find chunksize and number of chunks required for initialzation
    chunkSize = io_buffer_size / QQ
    nChunks = ceiling( dble(nElems)/dble(chunkSize) )
    allocate(xc(chunkSize, 3))
    allocate(rho(  chunkSize))
    allocate(vel(3,chunkSize))
    ! use AOS layout for fEq and fnEq
    allocate(fEq( QQ*chunkSize))

    do iChunk = 1, nChunks
      ! Number of elements read so far in previous chunks.
      elemOff = ( (iChunk-1)*chunksize )
      nChunkElems = min(chunkSize, nElems - elemOff)

      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        ! Calculate the coordinates
        xc(iElem,1:3) = tem_BaryOfId( tree,     &
          &                           me%levelDesc(iLevel)%total(elemPos))
      end do

      rho(1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(1),  &
        &                                   coord = xc(1:nChunkElems,1:3),  &
        &                                   n     = nChunkElems             )
      vel(1,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(2), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      vel(2,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(3), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )
      vel(3,1:nChunkElems) = tem_spatial_for( me    = field%ic%ini_state(4), &
        &                                     coord = xc(1:nChunkElems,1:3), &
        &                                     n     = nChunkElems            )

      ! convert these quantities from physics to LB
!cdir nodep
!ibm* novector
!dir$ novector
      do iElem = 1, nChunkElems
        rho(iElem) = rho(iElem) * cs2inv *inv_p
        vel(1,iElem) = vel(1,iElem) * inv_v
        vel(2,iElem) = vel(2,iElem) * inv_v
        vel(3,iElem) = vel(3,iElem) * inv_v
      end do

      call me%derVarPos(iField)%equilFromMacro( &
        &    density  = rho(1:nChunkElems),     &
        &    velocity = vel(1:3,1:nChunkElems), &
        &    iField   = iField,                 &
        &    nElems   = nChunkElems,            &
        &    varSys   = me%varSys,              &
        &    layout   = me%layout,              &
        &    res      = fEq                     )

      ! assign pdf by fEq + fNeq
      do iElem = 1, nChunkElems
        elemPos = elemOff + iElem
        offset = (iElem-1)*QQ
        do iDir = 1, QQ
          state( ( elempos-1)* nscalars+idir+( 1-1)* qq ) &
            & = fEq(offset + iDir)
        end do ! iDir
      end do ! iElem

    end do ! chunk

    ! deallocate memory for next chunk
    deallocate(xc)
    deallocate(rho)
    deallocate(vel)
    deallocate(fEq)

  end subroutine mus_init_isotherm_acEq
! **************************************************************************** !


! **************************************************************************** !
  !> Recursively fill all the helper elements (i.e. ghost, halo) with valid
  !! information from the fluid elements.
  !!
  !! This step is required before each run of the simulation. It would be
  !! possible to fill also the helper elements with the initial conditions.
  !! However, we are only able to fill the fluid elements with valid data
  !! (restart files have no information about the helper elements)
  !!
  recursive subroutine fillHelperElementsFineToCoarse( scheme, general, &
    &                                                  physics, iLevel, &
    &                                                  maxLevel         )
    ! --------------------------------------------------------------------------
    !> containers for the scheme
    !! contains interpolation type
    type(mus_scheme_type), intent(inout) :: scheme
    !> global parameters
    type(tem_general_type), intent(in) :: general
    type(mus_physics_type), intent(in) :: physics
    !> global flow quantity properties
    integer, intent(in) :: maxLevel
    !> Level counter variable
    integer, intent(in) :: iLevel
    ! --------------------------------------------------------------------------
    integer :: sNext, tNext ! time layer to use for source and target
    ! --------------------------------------------------------------------------

    tNext = scheme%pdf( iLevel )%nNext

write(dbgUnit(5), "(A)") ''
write(dbgUnit(5), "(A)") '---- Enter fillHelperElementsFineToCoarse -----------'
write(dbgUnit(5), "(A,I0)") 'level: ', iLevel
write(dbgUnit(5), "(A,I0)") ' nNow: ', scheme%pdf( iLevel )%nNow
write(dbgUnit(5), "(A,I0)") 'nNext: ', scheme%pdf( iLevel )%nNext

    if( iLevel < maxLevel ) then

      sNext = scheme%pdf( iLevel+1 )%nNext

      call fillHelperElementsFineToCoarse( scheme, general, physics,  &
        &                                  iLevel+1, maxLevel         )

      ! interpolate state variable PDF
      call scheme%intp%fillMineFromFiner%do_intp(                          &
        &    fieldProp   = scheme%field(:)%fieldProp,                      &
        &    sState      = scheme%state(iLevel+1)%val(:,sNext),            &
        &    sNeigh      = scheme%pdf(iLevel+1)%Neigh,                     &
        &    snSize      = scheme%pdf(iLevel+1)%nSize,                     &
        &    sAuxField   = scheme%auxField(iLevel+1)%val(:),               &
        &    tState      = scheme%state(iLevel)%val(:,tNext),              &
        &    tNeigh      = scheme%pdf(iLevel)%Neigh,                       &
        &    tnSize      = scheme%pdf(iLevel)%nSize,                       &
        &    tLevelDesc  = scheme%levelDesc(iLevel),                       &
        &    level       = iLevel,                                         &
        &    layout      = scheme%layout,                                  &
        &    nTargets    = scheme%levelDesc( iLevel )%intpFromFiner%nVals, &
        &    targetList  = scheme%levelDesc( iLevel )%intpFromFiner%Val,   &
        &    physics     = physics,                                        &
        &    varSys      = scheme%varSys,                                  &
        &    derVarPos   = scheme%derVarPos(:),                            &
        &    time        = general%simControl%now                          )

      call general%commPattern%exchange_real(                           &
        &     recv    = scheme%levelDesc( iLevel )%recvbufferFromFiner, &
        &     send    = scheme%levelDesc( iLevel )%sendbufferFromFiner, &
        &     state   = scheme%state( iLevel )%val( :, tNext ),         &
        &     message_flag   = iLevel,                                  &
        &     comm    = general%proc%comm )
    end if

    call general%commPattern%exchange_real(                  &
      &    recv    = scheme%levelDesc( iLevel )%recvbuffer,  &
      &    send    = scheme%levelDesc( iLevel )%sendbuffer,  &
      &    state   = scheme%state( iLevel )%val( :, tNext ), &
      &    message_flag   = iLevel,                          &
      &    comm    = general%proc%comm )

write(dbgUnit(5), *) '---- Leave fillHelperElementsFineToCoarse --------------'
write(dbgUnit(5), *) ''

  end subroutine fillHelperElementsFineToCoarse
! **************************************************************************** !


! **************************************************************************** !
  !> Recursively fill all the helper elements (i.e. ghost, halo) with valid
  !! information from the fluid elements.
  !!
  !! This step is required before each run of the simulation. It would be
  !! possible to fill also the helper elements with the initial conditions.
  !! However, we are only able to fill the fluid elements with valid data
  !! (restart files have no information about the helper elements)
  !!
  recursive subroutine fillHelperElementsCoarseToFine( scheme, general,   &
    &                                                  physics, iLevel,   &
    &                                                  minLevel, maxLevel )
    ! -------------------------------------------------------------------------
    !> containers for the scheme
    !! contains interpolation type
    type(mus_scheme_type), intent(inout) :: scheme
    !> global parameters
    type(tem_general_type), intent(in) :: general
    type(mus_physics_type), intent(in) :: physics
    !> level range
    integer, intent(in)       :: minLevel, maxLevel
    !> Level counter variable
    integer, intent(in)       :: iLevel
    ! --------------------------------------------------------------------------
    integer :: sNext, tNext ! time layer to use for source and target
    integer :: iOrder
    ! --------------------------------------------------------------------------

    sNext = scheme%pdf( iLevel )%nNext

write(dbgUnit(5), "(A)") ''
write(dbgUnit(5), "(A)") '---- Enter fillHelperElementsCoarseToFine -----------'
write(dbgUnit(5), "(A,I0)") 'level: ', iLevel
write(dbgUnit(5), "(A,I0)") ' nNow: ', scheme%pdf( iLevel )%nNow
write(dbgUnit(5), "(A,I0)") 'nNext: ', scheme%pdf( iLevel )%nNext

    if (iLevel > minLevel) then
      call general%commPattern%exchange_real(                             &
        &     recv    = scheme%levelDesc( iLevel )%recvbufferFromCoarser, &
        &     send    = scheme%levelDesc( iLevel )%sendbufferFromCoarser, &
        &     state   = scheme%state( iLevel )%val(:,sNext),              &
        &     message_flag   = iLevel,                                    &
        &     comm    = general%proc%comm )
    end if

    if( iLevel < maxLevel ) then

      tNext = scheme%pdf( iLevel+1 )%nNext

      do iOrder = 0, scheme%intp%config%order
        ! interpolate state variable
        call scheme%intp%fillFinerFromMe(iOrder)%do_intp(               &
          &    fieldProp   = scheme%field(:)%fieldProp,                 &
          &    sState      = scheme%state(iLevel)%val(:,sNext),         &
          &    sNeigh      = scheme%pdf(iLevel)%Neigh,                  &
          &    snSize      = scheme%pdf(iLevel)%nSize,                  &
          &    sAuxField   = scheme%auxField(iLevel)%val(:),            &
          &    tState      = scheme%state(iLevel+1)%val(:,tNext),       &
          &    tNeigh      = scheme%pdf(iLevel+1)%Neigh,                &
          &    tnSize      = scheme%pdf(iLevel+1)%nSize,                &
          &    tLevelDesc  = scheme%levelDesc( iLevel+1 ),              &
          &    level       = iLevel,                                    &
          &    nTargets    = scheme%levelDesc(iLevel+1)                 &
          &                             %intpFromCoarser(iOrder)%nVals, &
          &    targetList  = scheme%levelDesc(iLevel+1)                 &
          &                             %intpFromCoarser(iOrder)%Val,   &
          &    layout      = scheme%layout,                             &
          &    physics     = physics,                                   &
          &    varSys      = scheme%varSys,                             &
          &    derVarPos   = scheme%derVarPos(:),                       &
          &    time        = general%simControl%now                     )
      end do

      call fillHelperElementsCoarseToFine( scheme, general, physics,  &
        &                         iLevel+1, minLevel, maxLevel )

    end if

write(dbgUnit(5), *) '---- Leave fillHelperElementsCoarseToFine ---------------'
write(dbgUnit(5), *) ''

  end subroutine fillHelperElementsCoarseToFine
! **************************************************************************** !

  ! ************************************************************************** !
  !> This routine initialize auxField variable from PDF values initialized by
  !! initial condition. AuxField is computed from state using SAVE access for
  !! fluid elements and interpolated for ghost elements
  subroutine mus_initAuxField(scheme, general, minLevel, maxLevel)
    !---------------------------------------------------------------------------
    !> containers for the scheme
    !! contains interpolation type
    type(mus_scheme_type), intent(inout) :: scheme
    !> contains commPattern, MPI communicator and simControl
    type(tem_general_type), intent(in) :: general
    !> level range
    integer, intent(in)       :: minLevel, maxLevel
    !---------------------------------------------------------------------------
    integer :: iLevel
    !---------------------------------------------------------------------------

    do iLevel = minLevel, maxLevel
      call mus_initAuxFieldFluidAndExchange(             &
        & auxField   = scheme%auxField(iLevel),          &
        & state      = scheme%state(iLevel)%val(:,       &
        &                     scheme%pdf(iLevel)%nNext), &
        & neigh      = scheme%pdf(iLevel)%neigh(:),      &
        & nElems     = scheme%pdf(iLevel)%nElems_fluid,  &
        & nSize      = scheme%pdf(iLevel)%nSize,         &
        & nFields    = scheme%nFields,                   &
        & iLevel     = iLevel,                           &
        & stencil    = scheme%layout%fStencil,           &
        & varSys     = scheme%varSys,                    &
        & derVarPos  = scheme%derVarPos,                 &
        & general    = general,                          &
        & quantities = scheme%layout%quantities          )
    end do

    ! Initilialize auxField ghostFromFiner and ghostFromCoarser with
    ! interpolation for init with PDF
    do iLevel = maxLevel-1, minLevel,-1
      call mus_intpAuxFieldCoarserAndExchange(     &
        & intp        = scheme%intp,               &
        & tAuxField   = scheme%auxField(iLevel),   &
        & sAuxField   = scheme%auxField(iLevel+1), &
        & tLevelDesc  = scheme%levelDesc(iLevel),  &
        & stencil     = scheme%layout%fStencil,    &
        & iLevel      = iLevel,                    &
        & nAuxScalars = scheme%varSys%nAuxScalars, &
        & general     = general                    )
    end do

    do iLevel = minLevel+1, maxLevel
      call mus_intpAuxFieldFinerAndExchange(        &
        & intp        = scheme%intp,                &
        & tAuxField   = scheme%auxField(iLevel),    &
        & sAuxField   = scheme%auxField(iLevel-1),  &
        & tLevelDesc  = scheme%levelDesc(iLevel),   &
        & stencil     = scheme%layout%fStencil,     &
        & iLevel      = iLevel,                     &
        & nAuxScalars = scheme%varSys%nAuxScalars,  &
        & general     = general                     )
    end do

  end subroutine mus_initAuxField
  ! ************************************************************************** !

end module  mus_flow_module
! **************************************************************************** !
