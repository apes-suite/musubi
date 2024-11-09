! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013-2015, 2017-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
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
! *************************************************************************** !
!> This module deals with the calculation of moments from pdfs
!!
module mus_moments_module

  ! include treelm modules
  use env_module,         only: rk, labelLen
  use tem_aux_module,     only: tem_abort
  use tem_math_module,    only: invert_matrix
  use tem_logging_module, only: logUnit
  use tem_debug_module,   only: dbgUnit
  use tem_matrix_module,  only: tem_matrix_dump

  ! include musubi modules
  use mus_moments_type_module,  only: mus_moment_type
  use mus_scheme_header_module, only: mus_scheme_header_type
  use mus_mrtInit_module,       only: MMtrD3Q19, MMivD3Q19, WMMIvD3Q27, &
    &                                 WMMtrD3Q27

  implicit none

  private

  public :: get_moment
  public :: get_momentVector
  public :: mus_init_moments
  public :: set_momentIndices
  public :: mus_dump_moments


 contains


! **************************************************************************** !
  !> Initialize the moment space
  !!
  subroutine mus_init_moments( me, QQ, cxDir, label, schemeHeader )
    ! --------------------------------------------------------------------------
    type( mus_moment_type ), intent(inout) :: me
    !>
    integer, intent(in) :: QQ
    integer, intent(in) :: cxDir(3, QQ)
    character(len=labelLen) :: label
    !> identifier of the scheme
    type( mus_scheme_header_type ), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------
    logical :: momDefined
    ! --------------------------------------------------------------------------
    ! if moments are already defined then do nothing
    if (me%mom_ready) return
    me%mom_ready = .true.

    write(logUnit(1),"(A)") 'Initializing the moment matrix by stencil: '&
      &                     //trim(label)

    ! Set both entries for the matrix to the stencil size (quadratic matrix)
    write(logUnit(10),"(A,I0)") 'Set nEntries: ', QQ
    me%toMoments%nEntries = QQ
    me%toPdf%nEntries = QQ

    write(logUnit(10),"(A,2I4)") 'allocate toMoment: ', &
      &                         me%toPdf%nEntries(1), me%toPdf%nEntries(2)
    if ( allocated( me%toMoments%A ) ) deallocate( me%toMoments%A )
    if ( allocated( me%toPdf%A     ) ) deallocate( me%toPdf%A     )
    allocate( me%toMoments%A( me%toMoments%nEntries(1), &
      &                       me%toMoments%nEntries(2) ))
    allocate( me%toPdf%A( me%toPdf%nEntries(1), me%toPdf%nEntries(2) ))

    ! Allocate and reset the moment labels
    write(logUnit(10),"(A,I0)") 'allocate momLabel: ', QQ
    if(allocated( me%momLabel )) deallocate( me%momLabel )
    allocate( me%momLabel(QQ))

    me%momLabel = ''
    momDefined = .true.
    select case(trim(schemeHeader%kind))
    case( 'fluid', 'fluid_incompressible' )
      call init_transformation_matrix_fluid( QQ, cxDir, label, me,      &
        &                                    me%toMoments%A, me%toPdf%A )
    case('multispecies_gas','multispecies_liquid')
      call init_transformation_matrix_MS( QQ, cxDir, label, me, me%toMoments%A,&
        &                                 me%toPdf%A )
    case default
      write(logUnit(1),*)'WARNING: Moments matrix is not defined for scheme ' &
        &              //'kind: '//trim(schemeHeader%kind)
      momDefined = .false.
    end select

    ! Dump moments transformation matrix
    if (momDefined) call mus_dump_moments(me, dbgUnit(10))
  end subroutine mus_init_moments
! **************************************************************************** !

! **************************************************************************** !
  !> Initialize Moments transformation matrix for LBM compressible and
  !! incompressible fluid model. This matrix must be consistent with the
  !! relaxation matrix used in compute kernel and interpolation routines
  subroutine init_transformation_matrix_fluid( QQ, cxDir, label, me, toMoment, &
    &                                          toPdf )
    ! --------------------------------------------------------------------------
    !>
    integer, intent(in) :: QQ
    integer, intent(in) :: cxDir(3, QQ)
    character(len=labelLen) :: label
    type( mus_moment_type ), intent(inout) :: me
    real(kind=rk), intent(inout) :: toMoment( me%toMoments%nEntries(1),&
      &                                       me%toMoments%nEntries(2) )
    real(kind=rk), intent(inout) :: toPdf( me%toPDF%nEntries(1), &
      &                                    me%toPDF%nEntries(2)  )
    ! --------------------------------------------------------------------------
    real(kind=rk), dimension(QQ) :: uV, cx, cy, cz, cxsqr, cysqr, czsqr, csqr
    real(kind=rk) :: invMat( me%toPDF%nEntries(1), me%toPDF%nEntries(2) )
    real(kind=rk) :: transMat( me%toPDF%nEntries(1), me%toPDF%nEntries(2) )
    ! --------------------------------------------------------------------------

    write(logUnit(10),"(A)") 'init transformation matrix'

    uV = 1.0_rk ! unit vector
    cx = get_momentVector(QQ, cxDir, (/1,0,0/))
    cy = get_momentVector(QQ, cxDir, (/0,1,0/))
    cz = get_momentVector(QQ, cxDir, (/0,0,1/))
    ! cx^2 is cx.*cx element wise multiplication
    cxsqr = get_momentVector(QQ, cxDir, (/2,0,0/))
    cysqr = get_momentVector(QQ, cxDir, (/0,2,0/))
    czsqr = get_momentVector(QQ, cxDir, (/0,0,2/))

    ! Set the moments one after another on the matrix.
    ! The moment description might come from some place else.
    ! Lua? -> check mus_scheme_type_module the type(mus_momentSpace_type)
    select case( trim(label) )
    case('d2q9')
      ! This moments are according to the paper
      ! Chen, S., Peng, C., Teng, Y., Wang, L. P., & Zhang, K. (2016).
      ! Improving lattice Boltzmann simulation of moving particles in a
      ! viscous flow using local grid refinement. Computers and Fluids,
      ! 136, 228–246. http://doi.org/10.1016/j.compfluid.2016.06.009
      !
      ! lattice velocity square
      csqr = cxsqr + cysqr

      me%momLabel(1) = 'density      '  ! rho = uV
      toMoment( 1,:) = uV
      me%momLabel(2) = 'energy       '  ! e = -4 + 3(cx^2+ cy^2)
      toMoment( 2,:) = -4._rk*uV + 3._rk*csqr
      me%momLabel(3) = 'energy_square'  ! e^2 = 4 - 21/2 (cx^2+cy^2)
                                        !         +  9/2 (cx^2+cy^2)^2
      toMoment( 3,:) = 4._rk*uV - 10.5_rk*csqr + 4.5_rk*csqr**2
      me%momLabel(4) = 'momX         '  ! j_x = cx
      toMoment( 4,:) = cx
      me%momLabel(5) = 'energy_fluxX '  ! q_x = (-5+3*(cx^2+cy^2))*cx
      toMoment( 5,:) = (-5._rk*uV + 3._rk*csqr)*cx
      me%momLabel(6) = 'momY         '  ! j_y = cy
      toMoment( 6,:) = cy
      me%momLabel(7) = 'energy_fluxY '  ! q_y = (-5+3*(cx^2+cy^2))*cy
      toMoment( 7,:) = (-5._rk*uV + 3._rk*csqr)*cy
      me%momLabel(8) = 'normal_stress'  ! p_xx = cx^2-cy^2
      toMoment( 8,:) = cxsqr - cysqr
      me%momLabel(9) = 'shear_stress '  ! p_xy = cx*cy
      toMoment( 9,:) = cx*cy ! p_xy
      ! set moments positions
      !position of each moments pos in D2Q9, QQ,cxDir
      allocate(me%first_moments(2))
      me%first_moments = (/ 4, 6 /)
      allocate(me%second_moments(4))
      me%second_moments = (/ 2, 3, 8, 9 /)
      allocate(me%third_moments(2))
      me%third_moments = (/ 5, 7 /)
      allocate(me%fourth_moments(0))

      transMat = transpose(toMoment)
      invMat = invert_matrix( matmul(toMoment, transMat) )
      toPDF = matmul(transMat, invMat)

    case('d3q15')
      ! D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S.
      ! (2002). Multiple-relaxation-time lattice Boltzmann models in three
      ! dimensions. Philosophical Transactions. Series A, Mathematical,
      ! Physical, and Engineering Sciences, 360(1792), 437–51.
      ! lattice velocity square
      csqr = cxsqr + cysqr + czsqr

      me%momLabel(1)  = 'density        ' ! rho = 1.0_rk
      toMoment( 1,:)  =  uV
      me%momLabel(2)  = 'energy         '
      toMoment( 2,:)  = csqr - 2.0_rk*uV
      me%momLabel(3)  = 'enegry_square  '
      toMoment( 3,:)  = (15._rk*csqr**2 - 55._rk*csqr + 32._rk*uV)*0.5_rk
      me%momLabel(4)  = 'momX           '
      toMoment( 4,:)  = cx
      me%momLabel(5)  = 'energy_fluxX   '
      toMoment( 5,:)  = (5._rk*csqr - 13._rk*uV)*cx*0.5_rk
      me%momLabel(6)  = 'momY           '
      toMoment( 6,:)  = cy
      me%momLabel(7)  = 'energy_fluxY   '
      toMoment( 7,:)  = (5._rk*csqr - 13._rk)*cy*0.5_rk
      me%momLabel(8)  = 'momZ           '
      toMoment( 8,:)  = cz
      me%momLabel(9)  = 'energy_fluxZ'
      toMoment( 9,:)  = (5._rk*csqr - 13._rk)*cz*0.5_rk
      me%momLabel(10) = 'normal_stress_1'
      toMoment(10,:)  = 3._rk*cxsqr - csqr
      me%momLabel(11) = 'normal_stress_2'
      toMoment(11,:)  = cysqr - czsqr
      me%momLabel(12) = 'shear_stress_XY'
      toMoment(12,:)  = cx * cy
      me%momLabel(13) = 'shear_stress_YZ'
      toMoment(13,:)  = cy * cz
      me%momLabel(14) = 'shear_stress_XZ'
      toMoment(14,:)  = cx * cz
      me%momLabel(15) = 'mXYZ'
      toMoment(15,:)  = cx * cy * cz
      ! set moments positions
      !position of each moments pos in D3Q15 layout
      allocate(me%first_moments(3))
      me%first_moments = (/ 4, 6, 8 /)
      allocate(me%second_moments(5))
      me%second_moments = (/ 10, 11, 12, 13, 14 /)
      allocate(me%third_moments(5))
      me%third_moments = (/ 2, 3, 5, 7, 9 /)
      allocate(me%fourth_moments(1))
      me%fourth_moments = (/ 15 /)

      transMat = transpose(toMoment)
      invMat = invert_matrix( matmul(toMoment, transMat) )
      toPDF = matmul(transMat, invMat)

    case('d3q19')
      ! This moments are according to the paper
      ! Tölke, J., Freudiger, S., & Krafczyk, M. (2006).
      ! An adaptive scheme using hierarchical grids for lattice Boltzmann
      ! multi-phase flow simulations. Computers & Fluids, 35(8–9), 820–830.
      toMoment = MMtrD3Q19
      toPDF = MMIvD3Q19
      ! c2 = cx^2+cy^2+cz^2
      me%momLabel(1)  = 'density        ' ! rho = uV
      me%momLabel(2)  = 'energy         ' ! csqr - 1
      me%momLabel(3)  = 'enegry_square  ' ! 3*c2^2 - 6*c2 + 1
      me%momLabel(4)  = 'momX           ' ! cx
      me%momLabel(5)  = 'energy_fluxX   ' ! (3*c2 - 5)*cx
      me%momLabel(6)  = 'momY           ' ! cy
      me%momLabel(7)  = 'energy_fluxY   ' ! (3*c2 - 5)*cy
      me%momLabel(8)  = 'momZ           ' ! cz
      me%momLabel(9)  = 'energy_fluxZ'    ! (3*c2 - 5)*cz
      me%momLabel(10) = 'normal_stress_1' ! 3*cx^2 - c2
      me%momLabel(11) = '2ndOrd_energy_1' ! (2*c2 - 3*uV)*(3*cx^2-c2)
      me%momLabel(12) = 'normal_stress_2' ! cy^2 - cz^2
      me%momLabel(13) = '2ndOrd_energy_2' ! (2*c2 - 3*uV)*(cy^2-cz^2)
      me%momLabel(14) = 'shear_stress_XY' ! cx * cy
      me%momLabel(15) = 'shear_stress_YZ' ! cy * cz
      me%momLabel(16) = 'shear_stress_XZ' ! cx * cz
      me%momLabel(17) = '3rdOrd_tensor_1' ! (cy^2 - cz^2)*cx
      me%momLabel(18) = '3rdOrd_tensor_2' ! (cz^2 - cx^2)*cy
      me%momLabel(19) = '3rdOrd_tensor_3' ! (cx^2 - cy^2)*cz
      ! set moments positions
      !position of each moments pos in D3Q19 layout
      allocate(me%first_moments(3))
      me%first_moments = (/ 4, 6, 8 /)
      allocate(me%second_moments(7))
      me%second_moments = (/ 2, 3, 10, 12, 14, 15, 16 /)
      allocate(me%third_moments(6))
      me%third_moments = (/ 5, 7, 9, 17, 18, 19 /)
      allocate(me%fourth_moments(2))
      me%fourth_moments = (/ 11, 13 /)

    case('d3q27')
      toMoment = WMMtrD3Q27
      toPDF = WMMIvD3Q27
      ! c2 = cx^2+cy^2+cz^2
      ! density
      me%momLabel(1)  = 'density        ' ! rho = uV
      ! momentum
      me%momLabel(2)  = 'momX           ' ! j_x = cx
      me%momLabel(3)  = 'momY           ' ! j_y = cy
      me%momLabel(4)  = 'momZ           ' ! j_z = cz
      ! second order tensor
      me%momLabel(5)  = 'shear_stress_XY' ! p_xy = cx*cy
      me%momLabel(6)  = 'shear_stress_YZ' ! p_yz = cy*cz
      me%momLabel(7)  = 'shear_stress_XZ' ! p_xz = cx*cz
      me%momLabel(8)  = 'normal_stress_1' ! p_xx = 3*cx^2-c2
      me%momLabel(9)  = 'normal_stress_2' ! p_ww = cy^2-cz^2
      ! kinetic energy
      me%momLabel(10) = 'kinetic_energy ' ! e = c2-1
      ! fluxes of energy
      me%momLabel(11) = 'energy_flux_X  ' ! q_x = (3*c2-5)*cx
      me%momLabel(12) = 'energy_flux_Y  ' ! q_y = (3*c2-5)*cy
      me%momLabel(13) = 'energy_flux_Z  ' ! q_z = (3*c2-5)*cz
      ! third order pseudo-vector
      me%momLabel(14) = '3rdOrd_tensor_1' ! tau_x = (cy^2-cz^2)*cx
      me%momLabel(15) = '3rdOrd_tensor_2' ! tau_y = (cz^2-cx^2)*cy
      me%momLabel(16) = '3rdOrd_tensor_3' ! tau_z = (cx^2-cy^2)*cz
      ! antisymmetric tensor
      me%momLabel(17) = 'mXYZ' ! XYZ = cx*cy*cz
      ! square energy
      me%momLabel(18) = 'square_energy  ' ! eps = 0.5*(3*c2^2-7*c2+2)
      ! product of the second order tensor and the energy
      me%momLabel(19) = '2ndOrd_tensor_1' ! XXe = (3*c2-4)*(3*cx^2-c2)
      me%momLabel(20) = '2ndOrd_tensor_2' ! WWe = (3*c2-4)*(cy^2-cz)
      me%momLabel(21) = '2ndOrd_tensor_3' ! XYe = (cx*cy)*(3*c2-7)
      me%momLabel(22) = '2ndOrd_tensor_4' ! YZe = (cy*cz)*(3*c2-7)
      me%momLabel(23) = '2ndOrd_tensor_5' ! XZe = (cz*cx)*(3*c2-7)
      ! fluxes of square of energy
      me%momLabel(24) = '2ndOrd_energy_X' ! qsqr_x = 0.5*cx*(9*c2^2-33*c2+26)
      me%momLabel(25) = '2ndOrd_energy_Y' ! qsqr_y = 0.5*cy*(9*c2^2-33*c2+26)
      me%momLabel(26) = '2ndOrd_energy_Z' ! qsqr_z = 0.5*cz*(9*c2^2-33*c2+26)
      ! cube energy
      me%momLabel(27) = 'cube_energy    ' ! e3 = 0.5*(9*c2^3-36*c2^2+33*c2-2)
      ! set moments positions
      !position of each moments pos in D3Q19 layout
      allocate(me%first_moments(3))
      me%first_moments = (/ 2, 3, 4 /)
      allocate(me%second_moments(6))
      me%second_moments = (/ 5, 6, 7, 8, 9, 10 /)
      allocate(me%third_moments(6))
      me%third_moments = (/ 11, 12, 13, 24, 25, 26 /)
      allocate(me%fourth_moments(11))
      me%fourth_moments = (/ 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 27 /)

    end select

  end subroutine init_transformation_matrix_fluid
! **************************************************************************** !

! **************************************************************************** !
  !> Intialize the moment transformation matrix for multispecies.
  !! This matrix must be consistent with relaxation matrix used for
  !! multispecies MRT collision routines
  subroutine init_transformation_matrix_MS( QQ, cxDir, label, me, toMoment, &
    &                                       toPdf )
    ! --------------------------------------------------------------------------
    !>
    integer, intent(in) :: QQ
    integer, intent(in) :: cxDir(3, QQ)
    character(len=labelLen) :: label
    type( mus_moment_type ), intent(inout) :: me
    real(kind=rk), intent(inout) :: toMoment( me%toMoments%nEntries(1),&
      &                                       me%toMoments%nEntries(2) )
    real(kind=rk), intent(inout) :: toPdf( me%toPDF%nEntries(1), &
      &                                    me%toPDF%nEntries(2)  )
    ! --------------------------------------------------------------------------
    logical :: requireMom
    ! --------------------------------------------------------------------------

    write(logUnit(10),"(A)") 'init transformation matrix'
    ! require moments for stencil
    requireMom = .true.

    ! Set the moments one after another on the matrix.
    ! The moment description might come from some place else.
    ! Lua? -> check mus_scheme_type_module the type(mus_momentSpace_type)
    select case( trim(label) )
    case('d2q9')
      me%momLabel(1) = 'density'
      toMoment( 1,:) = get_momentVector( QQ,cxDir, (/0,0,0/)) ! m0
      me%momLabel(2) = 'momX'
      toMoment( 2,:) = get_momentVector( QQ,cxDir, (/1,0,0/)) ! mx
      me%momLabel(3) = 'momY'
      toMoment( 3,:) = get_momentVector( QQ,cxDir, (/0,1,0/)) ! my
      me%momLabel(4) = 'momXX'
      toMoment( 4,:) = get_momentVector( QQ,cxDir, (/2,0,0/)) ! mxx
      me%momLabel(5) = 'momYY'
      toMoment( 5,:) = get_momentVector( QQ,cxDir, (/0,2,0/)) ! myy
      me%momLabel(6) = 'momXY'
      toMoment( 6,:) = get_momentVector( QQ,cxDir, (/1,1,0/)) ! mxy
      me%momLabel(7) = 'momXXY'
      toMoment( 7,:) = get_momentVector( QQ,cxDir, (/2,1,0/)) ! mxxy
      me%momLabel(8) = 'momXYY'
      toMoment( 8,:) = get_momentVector( QQ,cxDir, (/1,2,0/)) ! mxyy
      me%momLabel(9) = 'momXXYY'
      toMoment( 9,:) = get_momentVector( QQ,cxDir, (/2,2,0/)) ! mxxyy
      ! set moments positions
      !position of each moments pos in D3Q19 QQ,cxDir
      allocate(me%first_moments(2))
      me%first_moments = (/ 1, 2 /)
      allocate(me%second_moments(3))
      me%second_moments = (/ 1, 2, 4 /)
      allocate(me%third_moments(2))
      me%third_moments = (/ 1, 3 /)
      allocate(me%fourth_moments(1))
      me%fourth_moments = (/ 1 /)

    case('d2q9_reference')

      toMoment(1,:) = (/ 1._rk,  1._rk,  1._rk,  1._rk,  1._rk, &
        &                        1._rk,  1._rk,  1._rk,  1._rk /)

      toMoment(2,:) = (/-1._rk,  0._rk,  1._rk,  0._rk, -1._rk, &
        &                       -1._rk,  1._rk,  1._rk,  0._rk /)

      toMoment(3,:) = (/ 0._rk, -1._rk,  0._rk,  1._rk, -1._rk, &
        &                        1._rk, -1._rk,  1._rk,  0._rk /)

      toMoment(4,:) = (/ 1._rk,  0._rk,  1._rk,  0._rk,  1._rk, &
        &                        1._rk,  1._rk,  1._rk,  0._rk /)

      toMoment(5,:) = (/ 0._rk,  1._rk,  0._rk,  1._rk,  1._rk, &
        &                        1._rk,  1._rk,  1._rk,  0._rk /)

      toMoment(6,:) = (/ 0._rk,  0._rk,  0._rk,  0._rk,  1._rk, &
        &                       -1._rk, -1._rk,  1._rk,  0._rk /)

      toMoment(7,:) = (/ 0._rk,  0._rk,  0._rk,  0._rk, -1._rk, &
        &                        1._rk, -1._rk,  1._rk,  0._rk /)

      toMoment(8,:) = (/ 0._rk,  0._rk,  0._rk,  0._rk, -1._rk, &
        &                       -1._rk,  1._rk,  1._rk,  0._rk /)

      toMoment(9,:) = (/ 0._rk,  0._rk,  0._rk,  0._rk,  1._rk, &
        &                        1._rk,  1._rk,  1._rk,  0._rk /)

      ! set moments positions
      !position of each moments pos in D3Q19 QQ,cxDir
      allocate(me%first_moments( 2))
      allocate(me%second_moments(3))
      allocate(me%third_moments( 2))
      allocate(me%fourth_moments(1))
      me%first_moments  = (/ 1, 2 /)
      me%second_moments = (/ 1, 2, 4 /)
      me%third_moments  = (/ 1, 3 /)
      me%fourth_moments = (/ 1 /)

    case('d3q19')

      me%momLabel(1) = 'density'
      toMoment( 1,:) = get_momentVector( QQ,cxDir, (/0,0,0/)) ! m0
      me%momLabel(2) = 'momX   '
      toMoment( 2,:) = get_momentVector( QQ,cxDir, (/1,0,0/)) ! mx
      me%momLabel(3) = 'momY   '
      toMoment( 3,:) = get_momentVector( QQ,cxDir, (/0,1,0/)) ! my
      me%momLabel(4) = 'momZ   '
      toMoment( 4,:) = get_momentVector( QQ,cxDir, (/0,0,1/)) ! mz
      me%momLabel(5) = 'momXX  '
      toMoment( 5,:) = get_momentVector( QQ,cxDir, (/2,0,0/)) ! mxx
      me%momLabel(6) = 'momYY  '
      toMoment( 6,:) = get_momentVector( QQ,cxDir, (/0,2,0/)) ! myy
      me%momLabel(7) = 'momZZ  '
      toMoment( 7,:) = get_momentVector( QQ,cxDir, (/0,0,2/)) ! mzz
      me%momLabel(8) = 'momXY  '
      toMoment( 8,:) = get_momentVector( QQ,cxDir, (/1,1,0/)) ! mxy
      me%momLabel(9) = 'momYZ  '
      toMoment( 9,:) = get_momentVector( QQ,cxDir, (/0,1,1/)) ! myz
      me%momLabel(10) = 'momXZ'
      toMoment(10,:) = get_momentVector( QQ,cxDir, (/1,0,1/)) ! mzx
      me%momLabel(11) = 'momXXY '
      toMoment(11,:) = get_momentVector( QQ,cxDir, (/2,1,0/)) ! mxxy
      me%momLabel(12) = 'momXXZ '
      toMoment(12,:) = get_momentVector( QQ,cxDir, (/2,0,1/)) ! mxxz
      me%momLabel(13) = 'momYYX '
      toMoment(13,:) = get_momentVector( QQ,cxDir, (/1,2,0/)) ! myyx
      me%momLabel(14) = 'momYYZ '
      toMoment(14,:) = get_momentVector( QQ,cxDir, (/0,2,1/)) ! myyz
      me%momLabel(15) = 'momZZX '
      toMoment(15,:) = get_momentVector( QQ,cxDir, (/1,0,2/)) ! mzzx
      me%momLabel(16) = 'momZZY '
      toMoment(16,:) = get_momentVector( QQ,cxDir, (/0,1,2/)) ! mzzy
      me%momLabel(17) = 'momXXYY'
      toMoment(17,:) = get_momentVector( QQ,cxDir, (/2,2,0/)) ! mxxyy
      me%momLabel(18) = 'momYYZZ'
      toMoment(18,:) = get_momentVector( QQ,cxDir, (/0,2,2/)) ! myyzz
      me%momLabel(19) = 'momZZXX'
      toMoment(19,:) = get_momentVector( QQ,cxDir, (/2,0,2/)) ! mzzxx
      ! set moments positions
      !position of each moments pos in D3Q19 layout
      allocate(me%first_moments(3))
      me%first_moments = (/ 1, 2, 3 /)
      allocate(me%second_moments(6))
      me%second_moments = (/ 1, 2, 3, 4, 5, 6/)
      allocate(me%third_moments(6))
      me%third_moments = (/ 1, 2, 3, 4, 5, 6/)
      allocate(me%fourth_moments(3))
      me%fourth_moments = (/ 1, 2, 3 /)

    case('d3q27')

      me%momLabel(1) = 'density'
      toMoment( 1,:) = get_momentVector( QQ,cxDir, [0,0,0]) ! m0
      me%momLabel(2) = 'momX   '
      toMoment( 2,:) = get_momentVector( QQ,cxDir, [1,0,0]) ! mx
      me%momLabel(3) = 'momY   '
      toMoment( 3,:) = get_momentVector( QQ,cxDir, [0,1,0]) ! my
      me%momLabel(4) = 'momZ   '
      toMoment( 4,:) = get_momentVector( QQ,cxDir, [0,0,1]) ! mz
      me%momLabel(5) = 'momXX  '
      toMoment( 5,:) = get_momentVector( QQ,cxDir, [2,0,0]) ! mxx
      me%momLabel(6) = 'momYY  '
      toMoment( 6,:) = get_momentVector( QQ,cxDir, [0,2,0]) ! myy
      me%momLabel(7) = 'momZZ  '
      toMoment( 7,:) = get_momentVector( QQ,cxDir, [0,0,2]) ! mzz
      me%momLabel(8) = 'momXY  '
      toMoment( 8,:) = get_momentVector( QQ,cxDir, [1,1,0]) ! mxy
      me%momLabel(9) = 'momYZ  '
      toMoment( 9,:) = get_momentVector( QQ,cxDir, [0,1,1]) ! myz
      me%momLabel(10) = 'momXZ'
      toMoment(10,:) = get_momentVector( QQ,cxDir, [1,0,1]) ! mzx

      toMoment(11,:) = get_momentVector( QQ,cxDir, [2,1,0]) ! mxxy
      toMoment(12,:) = get_momentVector( QQ,cxDir, [2,0,1]) ! mxxz
      toMoment(13,:) = get_momentVector( QQ,cxDir, [1,2,0]) ! mxyy
      toMoment(14,:) = get_momentVector( QQ,cxDir, [0,2,1]) ! myyz
      toMoment(15,:) = get_momentVector( QQ,cxDir, [1,0,2]) ! mxzz
      toMoment(16,:) = get_momentVector( QQ,cxDir, [0,1,2]) ! myzz
      toMoment(17,:) = get_momentVector( QQ,cxDir, [1,1,1]) ! mxyz

      toMoment(18,:) = get_momentVector( QQ,cxDir, [2,2,0]) ! mxxyy
      toMoment(19,:) = get_momentVector( QQ,cxDir, [2,0,2]) ! mxxzz
      toMoment(20,:) = get_momentVector( QQ,cxDir, [0,2,2]) ! myyzz
      toMoment(21,:) = get_momentVector( QQ,cxDir, [2,1,1]) ! mxxyz
      toMoment(22,:) = get_momentVector( QQ,cxDir, [1,2,1]) ! mxyyz
      toMoment(23,:) = get_momentVector( QQ,cxDir, [1,1,2]) ! mxyzz

      toMoment(24,:) = get_momentVector( QQ,cxDir, [2,2,1]) ! mxxyyz
      toMoment(25,:) = get_momentVector( QQ,cxDir, [2,1,2]) ! mxxyzz
      toMoment(26,:) = get_momentVector( QQ,cxDir, [1,2,2]) ! mxyyzz

      toMoment(27,:) = get_momentVector( QQ,cxDir, [2,2,2]) ! mxxyyzz
    case default
      requireMom = .false.
    end select

    if ( requireMom ) then
      toPdf = invert_matrix( toMoment )
    end if

  end subroutine init_transformation_matrix_MS
! **************************************************************************** !


! **************************************************************************** !
  !> The integer moment vector for a given cxDir and order.
  !!
  !! Assuming 0**0 = 1 here.
  !!
  pure function mus_iMomVector( cxDir, expX, QQ ) result(iMom)
    ! --------------------------------------------------------------------------
    !> number of velocity channels (include rest)
    integer, intent(in) :: QQ
    !> discrete velocity
    integer, intent(in) :: cxDir(:,:)
    !> order in each direction
    integer, intent(in) :: expX(3)
    !>
    integer :: iMom(QQ)
    ! --------------------------------------------------------------------------
    integer :: nIndices
    integer :: iVal
    integer :: iX
    integer :: X_idx(3)
    integer :: xx(3)
    ! --------------------------------------------------------------------------

    nIndices = 0
    do iX=1,3
      if (expX(iX) > 0) then
        nIndices = nIndices + 1
        X_idx(nIndices) = iX
        xx(nIndices) = expX(iX)
      end if
    end do

    select case(nIndices)
    case(0)
      imom = 1
    case(1)
      do iVal=1,QQ
        imom(iVal) = cxDir(X_idx(1),iVal)**xx(1)
      end do
    case(2)
      do iVal=1,QQ
        imom(iVal) = cxDir(X_idx(1),iVal)**xx(1) &
          &        * cxDir(X_idx(2),iVal)**xx(2)
      end do
    case(3)
      do iVal=1,QQ
        imom(iVal) = cxDir(1,iVal)**xx(1) &
          &        * cxDir(2,iVal)**xx(2) &
          &        * cxDir(3,iVal)**xx(3)
      end do
    end select

  end function mus_iMomVector
! **************************************************************************** !


! **************************************************************************** !
  !> Calculate the moment of a centain order
  !! The moment of a distribution \( f_i \) is defined as:\n
  !! \[
  !!    m_{x^{p}y^{q}z^{r}} = \sum_{i}^{Q} c^{p}_{xi} c^{q}_{yi} c^{r}_{zi} f_i
  !! \]
  !! The fucntion argument `expX` is array of size 3,
  !! which contains the values of \f$p, q, r\f$
  !!
  pure function get_moment( QQ, cxDir, expX, pdf ) result( mom )
    ! --------------------------------------------------------------------------
    integer, intent(in) :: QQ
    !> distribution value
    real(kind=rk), intent(in) :: pdf(QQ)
    !>
    integer, intent(in) :: cxDir(3, QQ)
    integer, intent(in) :: expX(3) !< exponents of the moments
    real(kind=rk) :: mom     !< moment
    ! --------------------------------------------------------------------------

    !HK: assuming 0**0 should be 1 here
    mom = sum( pdf * get_momentVector( QQ, cxDir, expX ) )

  end function get_moment
! **************************************************************************** !


! **************************************************************************** !
  !> get the moment vector to calculate the moment from the pdf
  !!
  pure function get_momentVector( QQ, cxDir, expX ) result( mom )
    ! --------------------------------------------------------------------------
    !>
    integer, intent(in) :: QQ
    integer, intent(in) :: cxDir(3, QQ)
    !> exponents of the moments
    integer, intent(in) :: expX(3)
    !> moment vector
    real(kind=rk) :: mom(QQ)
    ! --------------------------------------------------------------------------

    !HK: assuming 0**0 should be 1 here
    Mom = real( mus_iMomVector(cxDir = cxDir,  &
      &                        expX  = expX,   &
      &                        QQ    = QQ),    kind = rk )

  end function get_momentVector
! **************************************************************************** !


! **************************************************************************** !
  !> set indices for accessing the pressure, velocity and the shear from a 1d
  !! vector
  !!
  subroutine set_momentIndices( nDims, iPress, iVelMin, iVelMax, iSMin, iSMax )
    ! --------------------------------------------------------------------------
    !> number of dimensions
    integer, intent(in) :: nDims
    !> index for the pressure / density
    integer, intent(out) :: iPress
    !> starting index for velocity
    integer, intent(out) :: iVelMin
    !> ending index for velocity
    integer, intent(out) :: iVelMax
    !> starting index for shear
    integer, intent(out) :: iSMin
    !> ending index for shear
    integer, intent(out) :: iSMax
    ! --------------------------------------------------------------------------

    ! Pressure / density is always first entry
    iPress = 1
    ! Indices from where / to where to read / store the velocity quantities
    iVelMin = 2
    iVelMax = nDims + 1
    ! Indices from where / to where to read / store the stress   quantities
    iSMin = iVelMax + 1
    iSMax = iSMin   + max(3*(nDims-1)-1, 0)

  end subroutine set_momentIndices
! **************************************************************************** !


! **************************************************************************** !
  !> Dump moments matrix: toPDF and toMoment
  subroutine mus_dump_moments(me, outUnit)
    ! --------------------------------------------------------------------------
    type( mus_moment_type ), intent(in) :: me
    integer, intent(in) :: outUnit
    ! --------------------------------------------------------------------------
    write(outUnit, "(A)") ' toMoment: '
    call tem_matrix_dump(me%toMoments, outUnit)
    write(outUnit,*)
    write(outUnit, "(A)") ' toPDF: '
    call tem_matrix_dump(me%toPDF, outUnit)

  end subroutine mus_dump_moments
! **************************************************************************** !

end module mus_moments_module
! **************************************************************************** !
