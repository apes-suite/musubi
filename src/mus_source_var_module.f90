! Copyright (c) 2016-2021 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
! Copyright (c) 2021-2022 Gregorio Gerardo Spinelli <gregoriogerardo.spinelli@dlr.de>
! Copyright (c) 2022 Kannan Masilamani <kannan.masilamani@dlr.de>
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
! ***************************************************************************** !
!> author: Kannan Masilamani
!! Module containing subroutines for building MUSUBI specific source
!! variables
!!
module mus_source_var_module

  ! include treelm modules
  use mpi
  use env_module,               only: rk
  use tem_param_module,         only: rho0, cs2, cs2inv
  use tem_varSys_module,        only: tem_varSys_type
  use tem_varMap_module,        only: tem_possible_variable_type, &
    &                                 init, append, truncate
  use tem_stringKeyValuePair_module, only: init, truncate, &
    &                                      tem_stringKeyValuePair_type
  use tem_time_module,          only: tem_time_type
  use tem_spacetime_fun_module, only: tem_spacetime_fun_type
  use tem_spacetime_var_module, only: tem_varSys_append_stFun
  use tem_dyn_array_module,    only: PositionOfVal

  ! include musubi modules
  use mus_physics_module,       only: mus_convertFac_type
  use mus_derVarPos_module,     only: mus_derVarPos_type
  use mus_source_type_module,   only: mus_source_op_type, &
                                      mus_source_type
  use mus_scheme_header_module,      only: mus_scheme_header_type

  implicit none
  private

  ! update auxField dependent source variables
  public :: mus_updateSrcVar_dynSponFld
  public :: mus_add_internal_source_var

  contains

  ! ************************************************************************** !
  !> Compute density and velocity in sponge layer for dynamic sponge
  subroutine mus_updateSrcVar_dynSponFld(fun, auxField, iLevel, varSys, &
    &                                    phyConvFac, derVarPos)
    ! ------------------------------------------------------------------------ !
    !> Description of method to update source
    class(mus_source_op_type), intent(inout) :: fun
    !> input auxField array on current level
    real(kind=rk), intent(in)          :: auxField(:)
    !> current level
    integer, intent(in) :: iLevel
    !> variable system definition
    type(tem_varSys_type), intent(in)     :: varSys
    !> Physics conversion factor on current level
    type(mus_convertFac_type), intent(in) :: phyConvFac
    !> position of derived quantities in varsys
    type(mus_derVarPos_type), intent(in)  :: derVarPos(:)
    ! --------------------------------------------------------------------------
    integer :: dens_pos, vel_pos(3)
    integer :: iElem, nElems,  elemOff
    real(kind=rk) :: dens, vel(3)
    real(kind=rk) :: smoothFac
    real(kind=rk) :: inv_rho_phy, inv_vel_phy
    real(kind=rk) :: dens_ref, vel_ref(3)
    ! --------------------------------------------------------------------------
    ! position of density and velocity field in auxField
    dens_pos = varSys%method%val(derVarPos(1)%density)%auxField_varPos(1)
    vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

    inv_rho_phy = 1.0_rk / phyConvFac%press * cs2inv
    inv_vel_phy = 1.0_rk / phyConvFac%vel
    ! Number of elements to apply source terms
    nElems = fun%elemLvl(iLevel)%nElems
    ! Calculate time average value
    ! Set initial density and velocity if no initialized
    associate( dynAvg => fun%elemLvl(iLevel)%dynAvg,         &
      &        posInTotal => fun%elemLvl(iLevel)%posInTotal  )
      ! Do Exponential moving average over nRecord window
      ! P_n = P_{n-1} + smoothFac * (P - P_{n-1})
      smoothFac = fun%absLayer%smoothFac
      ! DO time average for density
      if (fun%absLayer%config%isPressDyn) then
        if (dynAvg%isInitDens) then
          dynAvg%isInitDens = .false.
          ! Initialize dynAvg%dens with initialize density
          do iElem = 1, nElems
            ! element offset
            elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars
            dynAvg%dens(iElem) = auxField(elemOff+dens_pos)
          end do
        else
          do iElem = 1, nElems
            ! element offset
            elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars
            ! Local density and velocity
            dens = auxField(elemOff+dens_pos)
            ! New average
            dynAvg%dens(iElem) = dynAvg%dens(iElem)                     &
              &                + smoothFac * (dens - dynAvg%dens(iElem) )
          end do
        end if
      else
        ! target pressure in lattice unit
        dens_ref = fun%absLayer%config%target_pressure * inv_rho_phy
        dynAvg%dens(:) = dens_ref
      end if

      ! DO time average for velocity
      if (fun%absLayer%config%isVelDyn) then
        if (dynAvg%isInitVel) then
          dynAvg%isInitVel = .false.
          ! Initialize dynAvg%vel with initialize velocity
          do iElem = 1, nElems
            ! element offset
            elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars
            dynAvg%velX(iElem) = auxField(elemOff+vel_pos(1))
            dynAvg%velY(iElem) = auxField(elemOff+vel_pos(2))
            dynAvg%velZ(iElem) = auxField(elemOff+vel_pos(3))
          end do
        else
          do iElem = 1, nElems
            ! element offset
            elemoff = (posInTotal(iElem)-1)*varSys%nAuxScalars
            ! Local velocity
            vel(1) = auxField(elemOff+vel_pos(1))
            vel(2) = auxField(elemOff+vel_pos(2))
            vel(3) = auxField(elemOff+vel_pos(3))
            ! New average
            dynAvg%velX(iElem) = dynAvg%velX(iElem)                       &
              &                + smoothFac * (vel(1) - dynAvg%velX(iElem) )
            dynAvg%velY(iElem) = dynAvg%velY(iElem)                       &
              &                + smoothFac * (vel(2) - dynAvg%velY(iElem) )
            dynAvg%velZ(iElem) = dynAvg%velZ(iElem)                       &
              &                + smoothFac * (vel(3) - dynAvg%velZ(iElem) )
          end do
        end if
      else
        ! target velocity in lattice unit
        vel_ref(1:3) = fun%absLayer%config%target_velocity(1:3) * inv_vel_phy
        dynAvg%velX(:) = vel_ref(1)
        dynAvg%velY(:) = vel_ref(2)
        dynAvg%velZ(:) = vel_ref(3)
      end if
    end associate
  end subroutine mus_updateSrcVar_dynSponFld
  ! ************************************************************************** !


  ! ***************************************************************************!
  !> Routine load musubi source terms for given key.
  !! key is glob_source or source
  subroutine mus_add_internal_source_var(me, possVars, varSys, schemeHeader)
    ! --------------------------------------------------------------------------!
    !> Source variable type to initialize
    type(mus_source_type), intent(out) :: me
    !> possible source variables
    type(tem_possible_variable_type), intent(in) :: possVars
    !> Global variable system
    type(tem_varSys_type), intent(inout) :: varSys
    !> Identifier of the scheme
    type(mus_scheme_header_type), intent(in) :: schemeHeader
    ! --------------------------------------------------------------------------!
    integer :: hrrCorr_varPOS
    type(tem_stringKeyValuePair_type) :: kvp
    type(tem_spacetime_fun_type), pointer :: stfun(:)
    ! --------------------------------------------------------------------------!
    ! initialize growing array stringKeyValuePair
    call init( me = me%varDict )

    if (trim(schemeHeader%relaxation) == 'hrr_bgk_corrected' .or. &
      & trim(schemeHeader%relaxation) == 'prr_bgk_corrected' .or. &
      & trim (schemeHeader%relaxation) == 'rr_bgk_corrected'      ) then
      kvp%value = 'hrr_correction'
      kvp%key = 'hrr_correction'
      allocate(stfun(1))
      ! use global mesh for anonymous variable
      stfun(1)%subTree%useGlobalMesh = .true.
      hrrCorr_varPOS = PositionOfVal( me  = possVars%varName, &
        &                             val = trim(kvp%key)     )
      stfun(1)%nComps = possVars%nComponents%val(hrrCorr_varPOS)
      allocate(stfun(1)%const(stfun(1)%nComps))
      stfun(1)%const(:) = 0._rk
      call tem_varSys_append_stFun(                       &
            & varSys              = varSys,               &
            & stfun               = stfun,                &
            & varname             = kvp%value,            &
            & nComp               = stfun(1)%nComps,      &
            & evalType            = 'firstonly_asglobal'  )
      call append(me = me%varDict, val = kvp)
    endif

    ! truncate varDict
    call truncate( me = me%varDict )

    if (me%varDict%nVals > 0) then
      ! allocate source method
      allocate(me%method(me%varDict%nVals))
    else
      allocate(me%method(0))
    end if
  end subroutine mus_add_internal_source_var
  ! ***************************************************************************!

end module mus_source_var_module
! **************************************************************************** !
