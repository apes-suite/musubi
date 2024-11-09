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
! **************************************************************************** !
!> author: Kannan Masilamani
!! This module include the routine required for element wie dumping weight
!! for better load balancing.
!! Dump weights at end of simulation if 'write_weights' is defined in the config
!! file. Weights are based on element wise time measurements
module mus_weights_module
  ! include treelm modules
  use env_module,              only: PathLen, rk, tem_create_EndianSuffix
  use treelmesh_module,        only: treelmesh_type, tem_dump_weights
  use tem_construction_module, only: tem_levelDesc_type
  use tem_topology_module,     only: tem_levelOf
  use tem_element_module,      only: eT_fluid
  use tem_logging_module,      only: logUnit
  use tem_debug_module,        only: dbgUnit

  ! include musubi modules
  use mus_bc_header_module,          only: glob_boundary_type
  use mus_timer_module,              only: mus_timerHandles, get_computeTime, &
    &                                      get_intpFromCoarserTime,           &
    &                                      get_intpFromFinerTime,             &
    &                                      get_bcBufferTime,                  &
    &                                      get_boundaryTime


  implicit none

  private

  public :: mus_getWeights
  public :: mus_dumpWeights


contains


  ! ************************************************************************** !
  !> Calculate weights using timing from compute kernel, interpolation and
  !! boundary routines
  subroutine mus_getWeights(weights, tree, minLevel, maxLevel, levelDesc, &
    &                       nBCs, globBC)
    ! --------------------------------------------------------------------------
    ! weights for fluid elements
    real(kind=rk), intent(out) :: weights(:)
    !> geometry infomation
    type(treelmesh_type),intent(in) :: tree
    !> min level and max level
    integer, intent(in) :: minLevel, maxLevel
    !> Level descriptor
    type(tem_levelDesc_type), intent(in) :: levelDesc(minLevel:maxLevel)
    !> global IBM type
    ! type( mus_IBM_globType ), intent(in) :: globIBM
    !> Number of boundary conditions
    integer, intent(in) :: nBCs
    !> BC elements information
    type(glob_boundary_type), intent(in) :: globBC( nBCs )
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! Assign weights based on compute and intp time
    call assign_level_weights(weights)

    ! Assign weights based on boundary condition
    call assign_bc_weights(weights)

  contains

    ! ************************************************************************ !
    !> Assign weights based on compute kernel and interpolation
    subroutine assign_level_weights(weights)
      ! ------------------------------------------------------------------------
      ! weights for fluid elements
      real(kind=rk), intent(out) :: weights(:)
      ! ------------------------------------------------------------------------
      real(kind=rk) :: level_weight(minLevel:maxLevel)
      integer :: iLevel, iElem
      ! ------------------------------------------------------------------------

      write(logUnit(7),"(A)") 'Assign weights by Computation and Interpolation'
      write(dbgUnit(1),"(A)") 'Assign weights by Computation and Interpolation'
      write(dbgUnit(1),"(A5,4A10)") 'level', 'computeT', 'intpT', 'nElems', &
        &                           'weight'

      do iLevel = minLevel, maxLevel

        if( levelDesc(iLevel)%elem%nElems( eT_fluid ) > 0) then
          ! include compute and intp time
          level_weight(iLevel) = ( get_computeTime( iLevel ) &
            &                      + get_intpFromCoarserTime(iLevel) &
            &                      + get_intpFromFinerTime( iLevel ) &
            &                    ) / levelDesc(iLevel)%elem%nElems( eT_fluid )
        else
          level_weight(iLevel) = 0._rk
        end if
        write(dbgUnit(1),"(I5,E10.2,E10.2,I10,E10.2)") iLevel,   &
          &                           get_computeTime( iLevel ), &
          &    get_intpFromCoarserTime(iLevel)                   &
          &    + get_intpFromFinerTime(iLevel),                  &
          &    levelDesc(iLevel)%elem%nElems( eT_fluid ),        &
          &    level_weight(iLevel)

      end do ! iLevel = minLevel, maxLevel

      do iElem = 1, tree%nElems
        iLevel = tem_levelOf( tree%treeID( iElem ))
        weights( iElem ) = level_weight( iLevel )
      end do

      write(dbgUnit(1),"(A)") ''

      ! Debug output
!      write(dbgUnit(1),*) "Calculate level weight based on intp and " &
!        &                 // "compute time"
!      write(dbgUnit(1),"(A8, 4A15)") "level", "nFluids", "t intp", &
!        &                            "t compute", "weight"
!      write(dbgUnit(1),"(I8, I15, 3F15.10)") iLevel,                        &
!        &  levelDesc(iLevel)%elem%nElems( eT_fluid ),                       &
!        &  tem_getTimerVal(timerHandle = mus_timerHandles%compute(iLevel)), &
!        &  tem_getTimerVal(timerHandle = mus_timerHandles                   &
!        &                                %intpFromCoarser(iLevel)),         &
!        &  level_weight( iLevel )
    end subroutine assign_level_weights
    !************************************************************************* !

    ! ************************************************************************ !
    !> Append weights based on boundary condition
    subroutine assign_bc_weights(weights)
      ! ------------------------------------------------------------------------
      ! weights for fluid elements
      real(kind=rk), intent(inout) :: weights(:)
      ! -----------------------------------------------------------------------
      integer :: iBC, iLevel, iElem, posInTotal, posInTree
      real(kind=rk) :: bc_w, buffer_w
      ! ------------------------------------------------------------------------

      write(logUnit(7),"(A)") 'Assign weights by Boundary Conditions'
      write(dbgUnit(1),"(A)") 'Assign weights by Boundary Conditions'

      ! Assume all bc elements shared the same cost of filling bc buffer
      if ( levelDesc(minLevel)%bc_elemBuffer%nVals > 0 ) then
        buffer_w = get_bcBufferTime() / dble(levelDesc(minLevel)%bc_elemBuffer &
          &                                                     %nVals         )
      else
        buffer_w = 0.0_rk
      end if

      write(dbgUnit(1),"(A,ES10.2)") 'buffer weight: ', buffer_w
      write(dbgUnit(1),"(A2,3A10)") 'BC', 'Time', 'nElems', 'weight'
      do iBC = 1, nBCs
        if ( globBC( iBC )%nElems_local > 0 ) then
          ! calc BC weight: bc_w(ii) = bc_time(ii) / bc_nElems(ii)
          bc_w = get_boundaryTime( iBC ) / dble(globBC( iBC )%nElems_local) &
            &    + buffer_w
          do iLevel = minLevel, maxLevel
            do iElem = 1, globBC( iBC )%nElems_Fluid( iLevel )
              posInTotal = globBC(iBC)%elemLvl( iLevel )%elem%val( iElem )
              posInTree  = levelDesc(iLevel)%pntTID( posInTotal )
              weights( posInTree ) = weights( posInTree ) + bc_w
            end do
          end do

        write(dbgUnit(1),"(I2,ES10.2,I10,ES10.2)") iBC, &
          &                     get_boundaryTime( iBC ), &
          &                     globBC( iBC )%nElems_local, &
          &                     bc_w
        end if ! globBC( iBC )%nElems_local > 0
      end do
      write(dbgUnit(1),"(A)") ''

      end subroutine assign_bc_weights
      ! ********************************************************************** !

    end subroutine mus_getWeights
    ! ************************************************************************ !

    ! ************************************************************************ !
    !> Dump weights to a file.
    subroutine mus_dumpWeights(tree, weights, basename)
      ! ------------------------------------------------------------------------
      type(treelmesh_type), intent(in)   :: tree
      real(kind=rk), intent(in)          :: weights(:)
      character(len=pathLen), intent(in) :: basename
      ! ------------------------------------------------------------------------
      integer                            :: iError
      character(len=PathLen)             :: filename
      character(len=4)                   :: EndianSuffix
      ! ------------------------------------------------------------------------

      iError  = 0
      EndianSuffix = tem_create_EndianSuffix()
      filename = trim(basename)//EndianSuffix

      call tem_dump_weights( me        = tree,     &
        &                    filename  = filename, &
        &                    weights   = weights   )

    end subroutine mus_dumpWeights
  ! ************************************************************************** !

end module mus_weights_module
