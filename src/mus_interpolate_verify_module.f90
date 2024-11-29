! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013-2016, 2018, 2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <raphael.haupt@uni-siegen.de>
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
!> author: Manuel Hasert
!! Interpolation scheme tools
!!
!! For an overview over implemented interpolation methods, see
!! [Interpolation methods](../page/features/intp_methods.html)
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
module mus_interpolate_verify_module

  ! include aotus modules
  use aotus_module,     only: flu_State

  ! include treelm modules
  use mpi
  use env_module,              only: rk, labelLen, long_k, newUnit, pathLen,   &
    &                                rk_mpi
  use treelmesh_module,        only: treelmesh_type
  use tem_aux_module,          only: tem_abort
  use tem_element_module,      only: eT_GhostFromCoarser,                  &
    &                                eT_ghostFromFiner, eT_fluid
  use tem_geometry_module,     only: tem_baryOfId
  use tem_grow_array_module,   only: grw_intArray_type, init, append, destroy
  use tem_param_module,        only: cs2inv, cs2
  use tem_varSys_module,       only: tem_varSys_type
  use tem_subTree_type_module, only: tem_subTree_type
  use tem_subTree_module,      only: tem_subTree_from
  use tem_logging_module,      only: logUnit
  use tem_spatial_module,      only: tem_spatial_for
  use tem_construction_module, only: tem_levelDesc_type
  use tem_general_module,      only: tem_general_type

  ! include musubi modules
  use mus_config_module,             only: mus_open_config
  use mus_param_module,              only: mus_param_type
  use mus_interpolate_header_module, only: mus_interpolation_type
  use mus_pdf_module,                only: pdf_data_type
  use mus_scheme_layout_module,      only: mus_scheme_layout_type
  use mus_scheme_type_module,        only: mus_scheme_type
  use mus_moments_module,            only: set_momentIndices
  use mus_derVarPos_module,          only: mus_derVarPos_type
  use mus_physics_module,            only: mus_convertFac_type

  implicit none

  private

  public :: mus_testInterpolation

 contains

! ****************************************************************************** !
  !> Call tests to determine the actual error from the interpolation routines on
  !! the ghost elements. Compare against the analytical solution, which is
  !! given in terms of the initial conditions.
  !! Call this routine after the initial values are set and the ghost elements
  !! have been filled once, but no computation was started
  !! -> after fillHelperElements in the mus_aux_module
  !!
  subroutine mus_testInterpolation( scheme, tree, general, fac, iLevel, minLevel,&
    &                               maxLevel, pdf )
    ! ---------------------------------------------------------------------------
    !>
    type(mus_scheme_type), intent(inout)  :: scheme
    !> treelm tree
    type( treelmesh_type ), intent(in) :: tree
    !>
    type( pdf_data_type ), intent(inout) :: pdf(tree%global%minlevel:tree%global%maxlevel)
    !>
    type(tem_general_type), intent(in) :: general
    type(mus_convertFac_type), intent(in) :: fac
    !>
    integer, intent(in) :: minLevel, maxLevel
    !> Level counter variable
    integer, intent(in) :: iLevel
    ! ---------------------------------------------------------------------------
    integer :: iElem     ! element counter
    type( grw_intArray_type ) :: local_ind
    ! ---------------------------------------------------------------------------

    write(logUnit(3),'(A,I0)') ' Check interpolation error on level:', iLevel

    ! Check fluid error
    if( scheme%intp%config%testFluids ) then
      call init( me = local_ind, length = pdf( iLevel )%nElems_fluid )
      do iElem = 1, pdf( iLevel )%nElems_fluid
        call append( me = local_ind, val = iElem )
      end do

      call mus_intp_error( pdf       = pdf(iLevel),     &
        &                  scheme    = scheme,          &
        &                  general   = general,         &
        &                  fac       = fac,             &
        &                  tree      = tree,            &
        &                  level     = iLevel,          &
        &                  eType     = eT_fluid,        &
        &                  intp      = scheme%intp,     &
        &                  ind       = local_ind,       &
        &                  layout    = scheme%layout,   &
        &                  varSys    = scheme%varSys,   &
        &                  derVarPos = scheme%derVarPos )

      call destroy( me = local_ind )
    end if ! checkFluid

    ! Check FromFiner Error
    if( iLevel < maxLevel ) then
      write(logUnit(1),*)' Checking fromFiner interpolation error...'
      call mus_intp_error( pdf   = pdf(iLevel),                 &
          &                scheme    = scheme,                  &
          &                general   = general,                 &
          &                fac       = fac,                     &
          &                tree      = tree,                    &
          &                level     = iLevel,                  &
          &                eType     = eT_ghostFromFiner,       &
          &                intp      = scheme%intp,             &
          &                ind       = scheme%levelDesc(iLevel) &
          &                                  %intpFromFiner,    &
          &                layout    = scheme%layout,           &
          &                varSys    = scheme%varSys,           &
          &                derVarPos = scheme%derVarPos         )
    end if ! check FF

    ! check FromCoarser error
    if( iLevel > minLevel ) then
      write(logUnit(1),*)' Checking fromCoarser interpolation error...'
      call mus_intp_error(                                                   &
        &      pdf       = pdf(iLevel),                                      &
        &      scheme    = scheme,                                           &
        &      general   = general,                                          &
        &      fac       = fac,                                              &
        &      tree      = tree,                                             &
        &      level     = iLevel,                                           &
        &      eType     = eT_ghostFromCoarser,                              &
        &      intp      = scheme%intp,                                      &
        &      ind       = scheme%levelDesc(iLevel)                          &
        &                        %intpFromCoarser(scheme%intp%config%order), &
        &      layout    = scheme%layout,                                    &
        &      varSys    = scheme%varSys,                                    &
        &      derVarPos = scheme%derVarPos                                  )
    end if

  end subroutine mus_testInterpolation
! ****************************************************************************** !


! ****************************************************************************** !
  !> Determine the numerical error of the interpolated quantities to the given
  !! initial conditions
  !!
  subroutine mus_intp_error( pdf, level, tree, general, fac, scheme, layout,   &
    &                        intp, ind, varSys, derVarPos, eType )
    ! ---------------------------------------------------------------------------
    !> global fluid parameters
    type( pdf_data_type ), intent(in) :: pdf
    !>
    type(tem_general_type), intent(in) :: general
    type(mus_convertFac_type), intent(in) :: fac
    !> treelm tree
    type( treelmesh_type ), intent(in) :: tree
    !> the current scheme
    type( mus_scheme_type ), intent(in)   :: scheme
    !> level
    integer, intent(in)   :: level
    !> element type, eT_ghostFromFiner, eT_ghostFromCoarser
    integer, intent(in)   :: eType
    !> the layout used
    type( mus_scheme_layout_type), intent(in) :: layout
    !> interpolation method info
    type( mus_interpolation_type ), intent(in)   :: intp
    !> indirection list
    type( grw_intArray_type ), intent(in)   :: ind
    !> global variable system
    type( tem_varSys_type ), intent(in) :: varSys
    !> required variable system which maps to global varsys
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    ! ---------------------------------------------------------------------------
    integer :: targetElem       ! treeId of current source element
    integer :: iVal, iPos       ! value counter
    integer :: iElem, nElems, globnElems    ! element counter for outer loop
    integer :: nStresses        ! number of stress variables
    integer :: indElem          ! element index in from the index list
    integer :: fUnit            ! unit for file write access
    real(kind=rk) :: targetPdf( layout%fStencil%QQ ) !< pdf to reconstruct
    real(kind=rk) :: targetMom( layout%fStencil%QQ )
    integer :: t_target ! which time layer to write to
    integer :: iErr ! mpi error
    integer :: nChunkELems, nDims, iVelMin, iVelMax
    integer :: iStressMin, iStressMax, iPress
    integer(kind=long_k) :: targetID
    integer(kind=long_k), allocatable :: minErrorID(:), maxErrorID(:)
    real(kind=rk) :: targetVel(3), refVel(3), targetPress(1), refPress
    real(kind=rk) :: targetStrain(6)
    ! temporary handle for the config file
    type(flu_state), allocatable :: conf(:)
    real(kind=rk), allocatable :: xc(:,:), ux(:), uy(:), uz(:), rho(:)
    real(kind=rk), allocatable :: Sxx(:,:)
    real(kind=rk), allocatable :: maxError(:), minError(:), error(:)
    real(kind=rk), allocatable :: errorNorm(:)
    real(kind=rk), allocatable :: globMaxError(:), globminError(:)
    real(kind=rk), allocatable :: globErrorNorm(:)
    character(len=labelLen) :: intpString
    character(len=PathLen) :: buffer, filename, header
    logical :: file_exists
    ! local subTree needed to use the function pointers
    !type(tem_subTree_type) :: subTree
    integer :: bound(2)
    integer :: nSize, QQ
    ! ---------------------------------------------------------------------------

    ! Allocate even though no test might been done, for parallel runs
    ! Other processes might have elements to check
    QQ = layout%fStencil%QQ
    allocate(Error   (QQ))
    allocate(minError(QQ))
    allocate(maxError(QQ))
    allocate(minErrorID(QQ))
    allocate(maxErrorID(QQ))
    allocate(errorNorm(QQ))
    allocate(globerrorNorm(QQ))
    allocate(globminError(QQ))
    allocate(globmaxError(QQ))
    minError = huge(minError)
    maxError = tiny(minError)
    error = 0._rk
    errorNorm = 0._rk
    nElems = ind%nVals

    ! copy the global information of the global tree to the subTree
    ! subTree%global = tree%global

    ! Only compute errors if there are any elements to test.
    if( nElems  > 0 ) then
      if( eType == eT_ghostFromCoarser) then
        intpString = 'intpFromCoarser'
      elseif( eType == eT_ghostFromFiner) then
        intpString = 'intpFromFiner'
      else
        intpString = 'noIntp'
      endif
      nDims = layout%fStencil%nDims
      call set_momentIndices( nDims, iPress, iVelMin, iVelMax,                 &
        &                     iStressMin, iStressMax )
      nStresses = iStressMax - iStressMin + 1

      call mus_open_config( conf     = conf,                                   &
        &                   filename = general%solver%configFile,       &
        &                   proc     = general%proc )
      nSize = pdf%nSize
      t_target = pdf%nNext
      nChunkElems = 1
      allocate(xc(nChunkElems, 3))
      allocate(rho(nChunkElems))
      allocate(ux(nChunkElems))
      allocate(uy(nChunkElems))
      allocate(uz(nChunkElems))
      allocate(Sxx(6,nChunkElems))
      minErrorID = 0_long_k
      maxErrorID = 0_long_k
      errorNorm = 0._rk

      maxError = tiny(1._rk)
      minError = huge(1._rk)
      do indElem = 1, nElems
        iElem = ind%val( indElem )
        ! Read the target element treeId
        targetElem = iElem + scheme%levelDesc( level )%offset( 1, eType )
        targetID   = scheme%levelDesc( level )%total( targetElem )
        ! get the target coordinates
        xc(1,:) = tem_baryOfId( tree, targetID )

        ! initialize the lower and upper bound
        bound(1)=indElem
        bound(2)=indElem
        ! call tem_subTree_from( me = subTree, treeID = (/targetID/))

        do iVal = 1, QQ
          targetPdf(iVal) =  &
          scheme%state(level)%val( ( targetelem-1)* qq+ ival,t_target  )
        enddo
        targetMom = matmul( layout%moment%toMoments%A, targetPdf)

        ! @todo: calculate pressure
        targetPress = sum(targetPDF)
        targetPress = targetPress * fac%press / 3

        call derVarPos(1)%velFromState(  &
          &           state  = targetPDF,   &
          &           iField = 1,        &
          &           nElems = 1,           &
          &           varSys = varSys,      &
          &           layout = layout,      &
          &           res    = targetVel    )
        targetVel = targetVel * fac%vel

        ! @todo: calculate strain rate
        targetStrain = targetStrain * fac%strainRate

        rho = tem_spatial_for( me    = scheme%field(1)%ic%ini_state(1), &
          &                    coord = xc,                                 &
          &                    n     = nChunkElems                         )
        ux =  tem_spatial_for( me    = scheme%field(1)%ic%ini_state(2), &
          &                    coord = xc,                                 &
          &                    n     = nChunkElems                         )
        uy =  tem_spatial_for( me    = scheme%field(1)%ic%ini_state(3), &
          &                    coord = xc,                                 &
          &                    n     = nChunkElems                         )
        uz =  tem_spatial_for( me    = scheme%field(1)%ic%ini_state(4), &
          &                    coord = xc,                                 &
          &                    n     = nChunkElems                         )
        ! Read in the shear rate tensor
        ! This corresponds to S = grad(u) + grad(u)^T
        Sxx = 0._rk
        Sxx(1,:) =  tem_spatial_for( me    = scheme%field(1)%ic%ini_state(5), &
          &                          coord = xc,                                 &
          &                          n     = nChunkElems                         )
        Sxx(2,:) =  tem_spatial_for( me    = scheme%field(1)%ic%ini_state(6), &
          &                          coord = xc,                                 &
          &                          n     = nChunkElems                         )
        Sxx(3,:) =  tem_spatial_for( me    = scheme%field(1)%ic%ini_state(7), &
          &                          coord = xc,                                 &
          &                          n     = nChunkElems                         )
        Sxx(4,:) =  tem_spatial_for( me    = scheme%field(1)%ic%ini_state(8), &
          &                          coord = xc,                                 &
          &                          n     = nChunkElems                         )
        Sxx(5,:) =  tem_spatial_for( me    = scheme%field(1)%ic%ini_state(9), &
          &                          coord = xc,                                 &
          &                          n     = nChunkElems                         )
        Sxx(6,:) =  tem_spatial_for( me    = scheme%field(1)%ic%ini_state(10), &
          &                          coord = xc,                                  &
          &                          n     = nChunkElems                          )
        if( scheme%layout%fStencil%nDims == 2 ) then
          Sxx(3,:) = Sxx(4,:)
          targetStrain(3) = targetStrain(4)
        end if

        refVel(:) = [ux(1), uy(1), uz(1)]
        refPress  = rho(1)
        ! Error computation
        error = 0._rk
        ! Absolute error
        error(1) = abs( (targetPress(1)-refPress) )
        error(iVelMin:iVelMax) = abs( (targetVel(1:nDims)-refVel(1:nDims)) )
        error(iStressMin:iStressMax ) = abs( (targetStrain(1:nStresses)      &
          &                                 - Sxx(1:nStresses,1)) )
        ! Relative error for pressure
        if( abs(refPress) > 0.000000001_rk ) then
          error(1) = error(1)/refPress
        end if

        ! Relative error for velocity
        do iVal = iVelMin, iVelMax
          if( abs(refVel(iVal-1)) > 0.000000001_rk ) then
            error( iVal ) = abs(error( iVal)/refVel(iVal-1 ))
          end if
        end do
        ! Relative error for shear stress
        do iVal = iStressMin, iStressMax
          iPos = iVal - iStressMin + 1
          if( abs(Sxx(iPos, 1)) > 0.000000001_rk ) then
            error( iVal ) = abs(error( iVal)/Sxx(iPos, 1 ))
          end if
        end do
        ! Determine min and max error and where it occurred
        do iVal = 1, size( error )
          errorNorm( iVal ) = errorNorm( iVal ) + error( iVal )**2
          if( abs(error( iVal ))> maxError( iVal)) then
            maxError( iVal ) = abs(error( iVal ))
            maxErrorID( iVal ) = targetID
          end if
          if( abs(error( iVal)) < minError( iVal)) then
            minError( iVal ) = abs(error( iVal ))
            minErrorID( iVal ) = targetID
          end if
        end do
        if( intp%config%testEachElement ) then
          write(logUnit(1),*) targetID
          write(logUnit(1),*) 'ref  ',refPress, refVel(1:nDims),             &
            &             Sxx(1:nStresses, 1)
          write(logUnit(1),*) 'val  ',targetPress(1), targetVel(1:nDims),    &
            &             targetStrain(1:nStresses)
          write(logUnit(1),*) 'error ',error(1:iStressMax)
        end if
      enddo
      ! Get rid of the max and min values
      do iVal = 1, size(error)
        if( minError( iVal ) < 1.E-90_rk ) then
          minError( iVal ) = 0._rk
        end if
        if( maxError( iVal ) < 1.E-90_rk ) then
          maxError( iVal ) = 0._rk
        end if
      end do
    end if ! nVals > 0


    if( general%proc%comm_size > 1 ) then
      ! Communicate error to root process
      call mpi_reduce( errorNorm, globErrorNorm, size(errorNorm), rk_mpi,      &
        &              mpi_sum, 0, general%proc%comm, iErr )
      call mpi_reduce( maxError, globMaxError, size(maxError), rk_mpi,         &
        &              mpi_max, 0, general%proc%comm, iErr )
      call mpi_reduce( minError, globMinError, size(minError), rk_mpi,         &
        &              mpi_min, 0, general%proc%comm, iErr )
      call mpi_reduce( nElems,   globnELems,  1, mpi_integer,                  &
        &              mpi_sum, 0, general%proc%comm, iErr )
    else
      globErrorNorm = errorNorm
      globMaxError  = maxError
      globMinError  = minError
      globnElems    = nElems
    end if

    if( general%proc%rank == 0 ) then
      ! Build the L2 norm error = 1/n*(sqrt( sum (phi - phi_error)^2 )
      globErrorNorm = sqrt( globErrorNorm ) / real( globnElems, kind=rk )
      write(buffer,'(2i4)') tree%global%minLevel, tree%global%maxlevel
      ! write(buffer,'(3a)') trim( buffer ), '  ', trim( params%scaling )
      write(buffer,'(3a)') trim( buffer ), '  ', trim( intp%config%method )
      do iVal = 1, size(error)
        write(logUnit(1),*)trim(layout%moment%momLabel( iVal )),               &
          &                   ' E max ', globmaxError( iVal ),                 &
          &                   ' E min ', globminError( iVal ),                 &
          &                   ' norm ', globErrorNorm( iVal )
        write(buffer, '(2a, e15.5)') trim( buffer ), '  ',                     &
          &     maxError( iVal )
      end do
      do iVal = 1, size( error )
        write(buffer, '(a, e16.5)') trim( buffer ),                            &
          &     globErrorNorm( iVal )
      end do

      write(filename,'(4a)') 'error_', trim(intpString), '.res'
      header = ''
      write(header,'(a,a1)') trim(header), '#'
      write(header,'(a,a)') trim(header), 'Lmin Lmax'
      write(header,'(3a)') trim(header), '  scaling'
      write(header,'(3a)') trim(header), '  intp '
      do iVal = 1, size( error )
        write(header,'(4a)') trim(header), ' max',                             &
          &         trim(adjustl(layout%moment%momLabel( iVal )))
      end do
      do iVal = 1, size( error )
        write(header,'(4a)') trim(header), ' rms',                             &
          &         trim(adjustl(layout%moment%momLabel( iVal )))
      end do

      inquire(file=filename,exist=file_exists)
      fUnit = newunit()
      open(unit=fUnit,file=trim(filename),position='append' )
      if( .not. file_exists ) then
         write(fUnit,'(a)') trim(header)
      end if
      write( fUnit, *) trim( buffer )
      close( fUnit )
    end if

  end subroutine mus_intp_error
! ****************************************************************************** !


end module mus_interpolate_verify_module
! ****************************************************************************** !
