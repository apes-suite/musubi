! Copyright (c) 2012-2016, 2018, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
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
! **************************************************************************** !
!> author: Kannan Masilmani
!! This module contains mus_species_type and routines to load species table
!! from config file.
!!
module mus_species_module

  ! include treelm modules
  use env_module,         only: rk, globalMaxLevels
  use tem_aux_module,     only: tem_abort
  use tem_tools_module,   only: tem_horizontalSpacer
  use tem_param_module,   only: cs2
  use tem_logging_module, only: logUnit

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_NonExistent,      &
    &                         aoterr_Fatal, aoterr_WrongType
  use aot_table_module, only: aot_table_open, aot_table_close, aot_table_length
  use aot_out_module,   only: aot_out_type, aot_out_val, aot_out_open_table,   &
    &                         aot_out_close_table

  ! include musubi modules
  use mus_physics_module, only: mus_physics_type

  implicit none

  private

  public :: mus_species_type
  public :: mus_load_species, compute_molWeightRatio, compute_bulkViscOmega
  public :: mus_species_out
  public :: Dxx, Dyy, Dzz, Dxy, Dxz, Dyz

  !> MRT species type
  type mrt_species_type
    !> relaxation matrix for mrt
    !! size of this matrix is (layout%QQ, layout%QQ)
    real(kind=rk), allocatable :: s_mrt(:,:)
    !> transformed relaxation matrix-moments factor
    !! omegaMoments = (Moments^-1.s_mrt.Moments)
    !!               .(I+(Moments^-1.s_mrt.Moments)/2.0)^-1
    real(kind=rk), allocatable :: omegaMoments(:,:)
    !> Omega factor for 2nd order force term
    !! omegaMomForce = (I+(Moments^-1.s_mrt.Moments)/2.0)^-1
    real(kind=rk), allocatable :: omegaMomForce(:,:)
  end type mrt_species_type

  !> diffusion tensor index
  integer, parameter :: Dxx = 1
  integer, parameter :: Dyy = 2
  integer, parameter :: Dzz = 3
  integer, parameter :: Dxy = 4
  integer, parameter :: Dxz = 5
  integer, parameter :: Dyz = 6

  !> this type contains species parameters
  !! @todo KM: extent level dependent parameter for multilevel
  type mus_species_type
    !> molecular weight of the species
    real(kind=rk) :: molWeight
    !> Inverse of molecular weight of the species.
    !! This parameter is required to convert mass density to mole density
    real(kind=rk) :: molWeightInv
    !> ratio of molecular weight  \f$ \phi_\sigma = min(M)/M_\sigma i \f$
    real(kind=rk) :: molWeigRatio
    !> coefficient of diffusivity  of the species (size of nspecies)
    real(kind=rk), allocatable :: diff_coeff(:)
    !> coefficient of resisivity of species which is
    !! reciprocal of diffusivity of the species
    real(kind=rk), allocatable :: resi_coeff(:)
    !> full diffusion tensor of the species
    real(kind=rk) :: diff_tensor(6)
    !! KM:@todo set diffusivity and resistivity for multilevel
    !> molar fraction of this species in the mixture
!    real(kind=rk) :: molarFrac
    !> Volume fraction of is species in the mixture
!    real(kind=rk) :: volFrac
    !> mrt relaxation for each level
    type(mrt_species_type) :: mrt(globalMaxLevels)
    !> bulk relaxation parameter
    !! omBulk_k = (2-phi_k)/3*bulkViscosity
    real(kind=rk) :: omBulk
    !> bulk relaxation parameter for each level
    real(kind=rk) :: ombulkLvl(globalMaxLevels)
    !> relaxation parameter for Nernst-Planck equation
    real(kind=rk) :: omega
    !> relaxation parameter for trt scheme
    real(kind=rk) :: lambda
    !> charge number of the species
    real(kind=rk) :: chargeNr
  end type mus_species_type


contains


! **************************************************************************** !
  !> this routines load species table from config file
  !!
  !!```
  !! species = { molweight = 1.0, diff_coeff = { 0.5,0.3,0.1 } }
  !!```
  subroutine mus_load_species( me, conf, parent, minLevel, nFields, physics, &
    &                          cs_lattice )
    ! --------------------------------------------------------------------------
    type( mus_species_type ), intent(out) :: me !< contains species information
    type( flu_State ) :: conf !< flu state
    integer, intent(in), optional :: parent !< parent lua handle
    integer, intent(in) :: minLevel
    integer, intent(in) :: nFields !< number of fields defined in lua file
    !> physics type to convert physics to lattice unit or vice versa
    type( mus_physics_type ), intent(in) :: physics
    !> lattice speed of sound calculated for defined stencil layout
    !! required to compute omega from potential diffusivity
    real(kind=rk), intent(in) :: cs_lattice
    ! --------------------------------------------------------------------------
    !local variables
    integer :: spc_handle, sub_handle
    integer :: iError
    integer, allocatable :: vError(:), errFatal(:)
    integer :: nCoeff
    ! --------------------------------------------------------------------------

    call aot_table_open( L = conf,                                             &
      &                  parent = parent,                                      &
      &                  thandle = spc_handle,                                 &
      &                  key = 'species')
    !> if species handle is not defined
    if ( spc_handle == 0 ) then
      write(logUnit(1),*)' No species table defined'
      call tem_abort()
    endif

    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*)' Loading species information'

    !get molecular weight
    call aot_get_val( L       = conf,                                          &
      &               thandle = spc_handle,                                    &
      &               key     = 'molweight',                                   &
      &               val     = me%molWeight,                                  &
      &               ErrCode = iError,                                        &
      &               default = 1.0_rk )
    me%molWeight = me%molWeight/physics%molWeight0
    me%molWeightInv = 1.0_rk / me%molWeight

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving molecular '// &
        &             'weight of species :'
      if ( btest( iError, aotErr_NonExistent ))                                &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
    end if

    !get trt relaxation parameter
    call aot_get_val( L       = conf,                                          &
      &               thandle = spc_handle,                                    &
      &               key     = 'lambda',                                      &
      &               val     = me%lambda,                                     &
      &               ErrCode = iError,                                        &
      &               default = 0.25_rk )

    !get specific charge
    call aot_get_val( L = conf, thandle = spc_handle, key = 'charge_nr',       &
      &               val = me%chargeNr, ErrCode = iError,                     &
      &               default = 0.0_rk )

    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(1),*) 'FATAL Error occured, while retrieving charge_nr '// &
        &             'of species :'
      if ( btest( iError, aotErr_NonExistent ))                                &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                                     &
        & write(logUnit(1),*)'Variable has wrong type!'
    end if

    !get diffusivities
    call aot_table_open( L = conf,                                             &
      &                  parent = spc_handle,                                  &
      &                  thandle = sub_handle,                                 &
      &                  key = 'diff_coeff')

    if ( sub_handle == 0 ) then
      nCoeff = 1
      !check whether resistivity table is defined
      call aot_table_open( L = conf, parent = spc_handle, thandle = sub_handle,&
        &                 key = 'resi_coeff')
      if ( sub_handle == 0 ) then
        ! coefficients are not defined as a table. try to load single constant
        ! value
        ! In the case of isotropic diffusion, there is no diffusion tensor
        ! In the case of anisotropic diffusion, diff_coeff is usually set 
        ! as the averaged diagonal diffusivity:
        ! \(D = (D_xx + D_yy + D_zz) / 3 \)
        allocate(me%diff_coeff(nCoeff))
        allocate(me%resi_coeff(nCoeff))
        ! diff_coeff may be single constant value
        call aot_get_val(L = conf, thandle = spc_handle, key = 'diff_coeff',   &
          &              val = me%diff_coeff(1), ErrCode = iError )

        if (btest(iError, aoterr_Fatal)) then
          ! check whether resi_coeff is defined
          call aot_get_val(L = conf, thandle = spc_handle, key = 'resi_coeff', &
            &              val = me%resi_coeff(1), ErrCode = iError )
          if (btest(iError, aoterr_Fatal)) then
            write(logUnit(1),*) 'FATAL Error occured, while retrieving '//     &
              &             'diff_coeff/resi_coeff of species :'
            if ( btest( iError, aotErr_NonExistent ))                          &
              & write(logUnit(1),*)'Variable not existent!'
            if (btest(iError, aoterr_WrongType))                               &
              & write(logUnit(1),*)'Variable has wrong type!'
            call tem_abort()
          endif
          me%diff_coeff = 1._rk/me%resi_coeff
        endif
        me%resi_coeff = 1._rk/me%diff_coeff
      else
        ! resisivity coeff is defined as table
        nCoeff = aot_table_length( L=conf, thandle = sub_handle )
        allocate( errFatal(nCoeff) )
        errFatal = aotErr_Fatal
        call aot_get_val( L = conf,                                            &
          &               thandle = spc_handle,                                &
          &               key = 'resi_coeff',                                  &
          &               val = me%resi_coeff,                                 &
          &               maxlength = nCoeff,                                  &
          &               ErrCode = vError )
        if ( any(btest(vError, errFatal)) ) then
           write(logUnit(1),*) 'FATAL Error occured, while retrieving '//      &
              &            'resi_coeff table'
           call tem_abort()
        endif
        allocate(me%diff_coeff(nCoeff))
        me%diff_coeff = 1._rk/me%resi_coeff
      endif
      call aot_table_close( L = conf, thandle = sub_handle )
    else
    ! diff_coeff is defined as a table
      nCoeff = aot_table_length( L=conf, thandle = sub_handle )
      allocate( errFatal(nCoeff) )
      errFatal = aotErr_Fatal
      call aot_get_val( L = conf,                                              &
        &               thandle = spc_handle,                                  &
        &               key = 'diff_coeff',                                    &
        &               val = me%diff_coeff,                                   &
        &               maxlength = nCoeff,                                    &
        &               ErrCode = vError )
      if ( any(btest(vError, errFatal)) ) then
         write(logUnit(1),*) 'FATAL Error occured, while retrieving '//        &
            &                'diff_coeff table'
         call tem_abort()
      endif
      allocate(me%resi_coeff(nCoeff))
      me%resi_coeff = 1._rk/me%diff_coeff
    endif
    call aot_table_close( L = conf, thandle = sub_handle )

    if(nCoeff /= nFields) then
      write(logUnit(1),*) 'ERROR: In loading diff_coeff or resi_coeff'
      write(logUnit(1),*) '       nCoeff does not match nFields'
      call tem_abort()
    endif

    !> Get diffusivity tensor if defined
    call aot_table_open( L = conf,             &
      &                  parent = spc_handle,  &
      &                  thandle = sub_handle, &
      &                  key = 'diff_tensor'   )

    if ( sub_handle /= 0 ) then
      call aot_get_val( L = conf,                  &
        &               thandle = sub_handle,      &
        &               key = 'Dxx',               &
        &               val = me%diff_tensor(Dxx), &
        &               ErrCode = iError,          &
        &               default = 0.0_rk           )

      call aot_get_val( L = conf,                  &
        &               thandle = sub_handle,      &
        &               key = 'Dyy',               &
        &               val = me%diff_tensor(Dyy), &
        &               ErrCode = iError,          &
        &               default = 0.0_rk           )
      
      call aot_get_val( L = conf,                  &
        &               thandle = sub_handle,      &
        &               key = 'Dzz',               &
        &               val = me%diff_tensor(Dzz), &
        &               ErrCode = iError,          &
        &               default = 0.0_rk           )

      call aot_get_val( L = conf,                  &
        &               thandle = sub_handle,      &
        &               key = 'Dxy',               &
        &               val = me%diff_tensor(Dxy), &
        &               ErrCode = iError,          &
        &               default = 0.0_rk           )

      call aot_get_val( L = conf,                  &
        &               thandle = sub_handle,      &
        &               key = 'Dxz',               &
        &               val = me%diff_tensor(Dxz), &
        &               ErrCode = iError,          &
        &               default = 0.0_rk           )

      call aot_get_val( L = conf,                  &
        &               thandle = sub_handle,      &
        &               key = 'Dyz',               &
        &               val = me%diff_tensor(Dyz), &
        &               ErrCode = iError,          &
        &               default = 0.0_rk           )

      if (any(me%diff_tensor /= 0.0_rk)) then
        write(logUnit(1),*) ' Species diffusion is anisotropic.'
        write(logUnit(1),*) '   Diffusion tensor components:'
        write(logUnit(1),*) '    Dxx: ', me%diff_tensor(Dxx)
        write(logUnit(1),*) '    Dyy: ', me%diff_tensor(Dyy)
        write(logUnit(1),*) '    Dzz: ', me%diff_tensor(Dzz)
        write(logUnit(1),*) '    Dxy: ', me%diff_tensor(Dxy)
        write(logUnit(1),*) '    Dxz: ', me%diff_tensor(Dxz)
        write(logUnit(1),*) '    Dyz: ', me%diff_tensor(Dyz)

        me%diff_tensor = me%diff_tensor / physics%fac(minLevel)%diffusivity
      else
        write(logUnit(1),*) ' Error: diffusion tensor is defined but all ' &
          &                 // 'components are zero!'
        call tem_abort()       
      end if
    end if

    call aot_table_close( L = conf, thandle = sub_handle )

    !> convert physics to lattice unit
    me%diff_coeff = me%diff_coeff/physics%fac(minLevel)%diffusivity
    me%resi_coeff = 1.0_rk/me%diff_coeff

    ! Relaxation parameter omega is compted from diffusivity coefficient.
    ! Used for Nernst-Planck equation
    ! For Emodel(Corr) used in passive scalar transport,
    ! omega is used as the free parameter to adjust stability 
    !> \todo KM: Compute omega for each level
    me%omega = 1.0_rk/(me%diff_coeff(1)/cs_lattice**2 + 0.5_rk)

    write(logUnit(1),*) ' Species properties          '
    write(logUnit(1),*) '   Molecular weight:         ', real(me%molWeight)
    write(logUnit(1),*) '   Inverse molecular weight: ', real(me%molWeightInv)
    write(logUnit(1),*) '   Charge number:            ', real(me%chargeNr)
    write(logUnit(1),*) '   Diffusivity coefficients: ', real(me%diff_coeff)
    write(logUnit(1),*) '   Resitivity coefficients:  ', real(me%resi_coeff)
    write(logUnit(1),*) '   Relaxation parameter:  ', real(me%omega)

    call aot_table_close( L=conf, thandle = spc_handle )
    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_load_species
! **************************************************************************** !


! **************************************************************************** !
  !> This routine computes the molecular weight ratio for all species
  !! based asinari model
  !!
  !! "Lattice Boltzmann scheme for mixture modeling: Analysis of the continuum
  !! diffusion regimes recovering Maxwell-Stefan model and incompressible
  !! Navier-Stokes equations. Pietro Asinari(2009)"
  !! \f$ m_\sigma = \frac{min_\varsigma(m_\varsigma)}{m_\sigma} \le 1 \f$
  subroutine compute_molWeightRatio( molWeights, molWeigRatios )
    ! --------------------------------------------------------------------------
    !> molecular weight of the species
    real(kind=rk), intent(in) :: molWeights(:)
    !> ratio of molecular weight  \f$ \phi_\sigma = min(M)/M_\sigma \f$
    real(kind=rk), intent(out) :: molWeigRatios(:)
    ! --------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*)' Compute species molecular weight ratio'
    molWeigRatios(:) = minval(MolWeights)/molWeights(:)

    write(logUnit(1),*)'  Molecular weight ratio:',real(molWeigRatios(:))
    call tem_horizontalSpacer(fUnit = logUnit(1))
  end subroutine compute_molWeightRatio
! **************************************************************************** !


! **************************************************************************** !
  !> This routine compute bulk viscosity omega for species for all levels
  !! omega_bulk = (2-molWeigRatio_k)/(3*bulk_visc)
  subroutine compute_bulkViscOmega( species, bulkVisc, bulkViscLvl, &
    &                               minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    !> contains species information
    type( mus_species_type ), intent(inout) :: species
    !> bulk viscosity of the mixture
    real(kind=rk), intent(in) :: bulkvisc
    real(kind=rk), intent(in) :: bulkviscLvl(globalMaxLevels)
    integer,       intent(in) :: minLevel, maxLevel
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------
    species%omBulk = (2.0_rk - species%molWeigRatio) * cs2 / bulkvisc
    write(logUnit(1),*) '   Bulk omega: ', real(species%omBulk)

    do iLevel = minLevel, maxLevel
      write(logUnit(1),*)'    level ', iLevel
      species%omBulkLvl(iLevel) = ( 2.0_rk - species%molWeigRatio ) * cs2      &
        &                       / bulkViscLvl(iLevel)
      write(logUnit(1),*) '   Bulk omega ', real(species%omBulkLvl(iLevel))
    end do
  end subroutine compute_bulkViscOmega
! **************************************************************************** !

  !> writes species propertries into a lua file
  !!
  subroutine mus_species_out(me, conf)
    ! --------------------------------------------------------------------------
    type( mus_species_type ), intent(in) :: me
    type(aot_out_type) :: conf
    ! --------------------------------------------------------------------------

    call aot_out_open_table( put_conf = conf, tname = 'species')

    call aot_out_val( put_conf = conf,                                         &
      &               vname    = 'molweight',                                  &
      &               val      = me%molWeight )
    call aot_out_val( put_conf = conf,                                         &
      &               vname    = 'diff_coeff',                                 &
      &               val      = me%diff_coeff )
    call aot_out_val( put_conf = conf,                                         &
      &               vname    = 'resi_coeff',                                 &
      &               val      = me%resi_coeff )
    call aot_out_val( put_conf = conf,                                         &
      &               vname    = 'charge_nr',                                  &
      &               val      = me%chargeNr )
    call aot_out_close_table( put_conf = conf )

  end subroutine mus_species_out
! **************************************************************************** !


end module mus_species_module
! **************************************************************************** !
