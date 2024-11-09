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
!! variables for turbulent channel flow. To avoid cyclic inclusions
!!
module mus_source_var_turbChanForce_module
    use, intrinsic :: iso_c_binding, only: c_f_pointer

    ! include treelm modules
    use mpi
    use env_module,               only: rk, rk_mpi
    use tem_param_module,         only: rho0, cs2
    use tem_varSys_module,        only: tem_varSys_type
    use tem_varMap_module,        only: init, append, truncate
    use tem_stringKeyValuePair_module, only: init, truncate
    use tem_topology_module,      only: tem_levelOf

    ! include musubi modules
    use mus_physics_module,       only: mus_convertFac_type
    use mus_derVarPos_module,     only: mus_derVarPos_type
    use mus_source_type_module,   only: mus_source_op_type
    use mus_varSys_module,        only: mus_varSys_data_type

    implicit none
    private

    ! update auxField dependent source variables
    public :: mus_updateSrcVar_turbChanForce

    contains

    ! ************************************************************************** !
    !> Compute dynamic force term using auxField for turbulent channel
    !! force.
    subroutine mus_updateSrcVar_turbChanForce(fun, auxField, iLevel, varSys, &
      &                                       phyConvFac, derVarPos)
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
      integer :: vel_pos(3)
      integer :: iElem, nElemsGlobal_uTau, posInTotal, elemOff
      type(mus_varSys_data_type), pointer :: fPtr
      integer :: iErr, iBC, iLvlLoc
      real(kind=rk) :: vel_bulk, vel_tau
      real(kind=rk) :: avgVelBulk, avgVelTau
      ! index 1 for vel_tau and 2 for vel_bulk to get global avg with one mpi call
      real(kind=rk) :: avgVel(2), avgVelGlobal(2)
      real(kind=rk) :: forceDyn, gradU(3,3,1)
      logical :: velTauFromBnd
      ! --------------------------------------------------------------------------

      ! position of velocity field in auxField
      vel_pos = varSys%method%val(derVarPos(1)%velocity)%auxField_varPos(1:3)

      ! First variable is always pdf so just use to access method_data
      call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
      ! Calculate average bulk velocityX over defined shape.
      associate( levelPointer => fPtr%solverData%geometry%levelPointer,          &
        &        auxFieldTot => fPtr%solverData%scheme%auxField,                 &
        &        gradDataTot => fPtr%solverData%scheme%gradData,                 &
        &        viscKineTot => fPtr%solverData%scheme%field(1)%fieldProp        &
        &                        %fluid%viscKine%dataOnLvl,                      &
        &        map2Global_umean => fun%turbChanForce%subTree_umean%map2global, &
        &        map2Global_utau => fun%turbChanForce%subTree_utau%map2global,   &
        &        physics => fPtr%solverData%physics,                             &
        &        tree => fPtr%solverData%geometry%tree,                          &
        &        grad => fPtr%solverData%scheme%Grad                             )
        ! Do average only at every multilevel cycle i.e. when iLevel=minLevel
        if (iLevel == tree%global%minLevel) then
         ! Initialize variables
          vel_bulk = 0.0_rk
          vel_tau = 0.0_rk
          velTauFromBnd = .false.
          nElemsGlobal_uTau = 0

          !> If wall model BC is applied than compute the friction velocity
          ! that is obtained directly from the wall model and perform the spatial
          ! averaging.
          do iBC = 1, fPtr%solverData%geometry%boundary%nBCtypes
            select case(trim(fPtr%solverData%scheme%field(1)%BC(iBC)%BC_kind))
            case('turbulent_wall', 'turbulent_wall_noneq_expol', &
              & 'turbulent_wall_eq')
              if (allocated(fPtr%solverData%scheme%field(1)%BC(iBC) &
                & %turbWallFunc%dataOnLvl)) then
                velTauFromBnd = .true.
                do iLvlLoc = tree%global%minLevel, tree%global%maxLevel
                  ! vel_tau in turbWallFunc is in lattice unit so convert it
                  vel_tau = vel_tau + sum( fPtr%solverData%scheme%field(1) &
                    &                      %BC(iBC)%turbWallFunc           &
                    &                      %dataOnLvl(iLvlLoc)%velTau(:) ) &
                    &               * physics%fac(iLvlLoc)%vel
                end do
                nElemsGlobal_uTau = nElemsGlobal_uTau               &
                  & + fPtr%solverData%scheme%globBC(iBC)%nElems_Total
              end if
            end select
          end do

          !> If wall model BC is not used than compute the friction velocity from
          !! the single sided finite difference's and perform the spatial averaging.
          !! Friction velocity is computed only on elements intersected by
          !! shape_utau defined in musubi.lua.
          if (.not. velTauFromBnd) then
            nElemsGlobal_uTau = fun%turbChanForce%nElemsGlobal_utau
            select case(fun%turbChanForce%flow_direction)
            case (1) ! x-direction
              do iElem = 1, fun%turbChanForce%subTree_utau%nElems
                ! map2Global refers to position in treeid list
                ! levelPointer refers to position in level wise total list
                posInTotal = levelPointer( map2Global_utau(iElem) )
                iLvlLoc = tem_levelOf( tree%treeID( map2Global_utau(iElem) ) )
                ! elemOffset refers to the previous number of elements processed
                gradU(:,:,:) = grad%U_ptr(                     &
                  & auxField    = auxFieldTot(iLvlLoc)%val(:), &
                  & gradData    = gradDataTot(iLvlLoc),        &
                  & velPos      = vel_Pos,                     &
                  & nAuxScalars = varSys%nAuxScalars,          &
                  & nDims       = 3,                           &
                  & nSolve      = 1,                           &
                  & elemOffset  = posInTotal - 1               )
                ! compute uTau using the dudy component of velocity gradient tensor
                vel_tau = vel_tau + sqrt( viscKineTot(iLvlLoc)%val(posInTotal) &
                  &                      * abs(gradU(1,2,1)) )                &
                  &               * physics%fac(iLvlLoc)%vel
                !write(dbgunit(1), *) 'posInTotal: ', posInTotal, 'gradU: ', gradU(1,2,1)
              end do
            case (2) ! y-direction
              do iElem = 1, fun%turbChanForce%subTree_utau%nElems
                ! map2Global refers to position in treeid list
                ! levelPointer refers to position in level wise total list
                posInTotal = levelPointer( map2Global_utau(iElem) )
                iLvlLoc = tem_levelOf( tree%treeID( map2Global_utau(iElem) ) )
                gradU(:,:,:) = grad%U_ptr(                     &
                  & auxField    = auxFieldTot(iLvlLoc)%val(:), &
                  & gradData    = gradDataTot(iLvlLoc),        &
                  & velPos      = vel_Pos,                     &
                  & nAuxScalars = varSys%nAuxScalars,          &
                  & nDims       = 3,                           &
                  & nSolve      = 1,                           &
                  & elemOffset  = posInTotal - 1               )
                ! compute uTau using the dvdx component of velocity gradient tensor
                vel_tau = vel_tau + sqrt( viscKineTot(iLvlLoc)%val(posInTotal) &
                  &                      * abs(gradU(2,1,1)) )                &
                  &               * physics%fac(iLvlLoc)%vel
              end do
            case (3) ! z-direction
              do iElem = 1, fun%turbChanForce%subTree_utau%nElems
                ! map2Global refers to position in treeid list
                ! levelPointer refers to position in level wise total list
                posInTotal = levelPointer( map2Global_utau(iElem) )
                iLvlLoc = tem_levelOf( tree%treeID( map2Global_utau(iElem) ) )
                gradU(:,:,:) = grad%U_ptr(                     &
                  & auxField    = auxFieldTot(iLvlLoc)%val(:), &
                  & gradData    = gradDataTot(iLvlLoc),        &
                  & velPos      = vel_Pos,                     &
                  & nAuxScalars = varSys%nAuxScalars,          &
                  & nDims       = 3,                           &
                  & nSolve      = 1,                           &
                  & elemOffset  = posInTotal - 1               )
                ! compute uTau using the dwdy component of velocity gradient tensor
                vel_tau = vel_tau + sqrt( viscKineTot(iLvlLoc)%val(posInTotal) &
                  &                      * abs(gradU(3,2,1)) )                &
                  &               * physics%fac(iLvlLoc)%vel
              end do
            end select
          end if

          !> Bulk mean velocity part of forcing is independent whether a wall
          ! modelBC is used or not
          select case(fun%turbChanForce%flow_direction)
          case (1) ! x-direction
            do iElem = 1, fun%turbChanForce%subTree_umean%nElems
              ! map2Global refers to position in treeid list
              ! levelPointer refers to position in level wise total list
              posInTotal = levelPointer( map2Global_umean(iElem) )
              ! elemoffset for auxField
              elemoff = (posInTotal-1)*varSys%nAuxScalars
              ! velocity X in defined shape to compute average
              iLvlLoc = tem_levelOf( tree%treeID( map2Global_umean(iElem) ) )
              vel_bulk = vel_bulk + auxFieldTot(iLvlLoc)%val(elemOff+vel_pos(1)) &
                &                 * physics%fac(iLvlLoc)%vel
            end do
          case (2) ! y-direction
            do iElem = 1, fun%turbChanForce%subTree_umean%nElems
              ! map2Global refers to position in treeid list
              ! levelPointer refers to position in level wise total list
              posInTotal = levelPointer( map2Global_umean(iElem) )
              ! elemoffset for auxField
              elemoff = (posInTotal-1)*varSys%nAuxScalars
              ! velocity X in defined shape to compute average
              iLvlLoc = tem_levelOf( tree%treeID( map2Global_umean(iElem) ) )
              vel_bulk = vel_bulk + auxFieldTot(iLvlLoc)%val(elemOff+vel_pos(2)) &
                &                 * physics%fac(iLvlLoc)%vel
            end do
            case (3) ! z-direction
            do iElem = 1, fun%turbChanForce%subTree_umean%nElems
              ! map2Global refers to position in treeid list
              ! levelPointer refers to position in level wise total list
              posInTotal = levelPointer( map2Global_umean(iElem) )
              ! elemoffset for auxField
              elemoff = (posInTotal-1)*varSys%nAuxScalars
              ! velocity X in defined shape to compute average
              iLvlLoc = tem_levelOf( tree%treeID( map2Global_umean(iElem) ) )
              vel_bulk = vel_bulk + auxFieldTot(iLvlLoc)%val(elemOff+vel_pos(3)) &
                &                 * physics%fac(iLvlLoc)%vel
            end do
          end select

          avgVel(1) = vel_tau
          avgVel(2) = vel_bulk
          ! Calculate total friction and bulk velocity
          call mpi_allreduce( avgVel, avgVelGlobal,                      &
            &                 2, rk_mpi, mpi_sum, tree%global%comm, iErr )

          ! compute average friction velocity
          avgVelTau = avgVelGlobal(1) / nElemsGlobal_uTau
          ! compute average bulk velocity in physical unit
          avgVelBulk = avgVelGlobal(2) / fun%turbChanForce%nElemsGlobal_umean


          ! Dynamic force term for turbulent channel
          ! F_dyn = (refVelBulk-avgVelBulk) * refVelBulk / refHeight
          forceDyn = avgVelTau**2 / fun%turbChanForce%refHeight               &
            &      + ( fun%turbChanForce%refVelBulk - avgVelBulk )            &
            &      * fun%turbChanForce%refVelBulk / fun%turbChanForce%refHeight
          fun%turbChanForce%forceDyn = 0.0_rk
          !write(dbgunit(1), *) 'avgVelBulk: ', avgVelBulk , 'avgVelTau:', avgVelTau
          !write(dbgunit(1), *) 'vel_tau: ', vel_tau
          !flush(dbgunit(1))
          !write(dbgunit(1), *) 'forceDyn_phy: ', forceDyn

          select case(fun%turbChanForce%flow_direction)
          case(1) ! x-direction
            fun%turbChanForce%forceDyn(1) = forceDyn
          case(2) ! y-direction
            fun%turbChanForce%forceDyn(2) = forceDyn
          case(3) ! z-direction
            fun%turbChanForce%forceDyn(3) = forceDyn
          end select
        end if
      end associate

    end subroutine mus_updateSrcVar_turbChanForce

  end module mus_source_var_turbChanForce_module
  ! **************************************************************************** !
