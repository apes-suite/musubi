! Copyright (c) 2012-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012-2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2015 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
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
! **************************************************************************** !
!> author: Simon Zimny
!! Module containing subroutines to prepare buffers for transformation of
!! quantities (tracking) and the restart/tracking IO.
!!
!! To save memory, the transformation of quantities and the IO-functionality of
!! MUSUBI is done in linear chunks. The size of the chunk can be defined in the
!! lua input file and is defined by the function [[tem_load_env_params]].
!! The default size is set to be 8MB.
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

module mus_buffer_module

  ! include treelm modules
  use mpi
  use env_module,            only: rk, long_k
  use tem_topology_module,   only: tem_LevelOf
  use tem_varSys_module,     only: tem_varSys_type
  use tem_varMap_module,     only: tem_varMap_type

  ! include musubi modules
  use mus_scheme_type_module, only: mus_scheme_type

  implicit none
  private

  public :: mus_pdf_serialize
  public :: mus_pdf_unserialize


contains


  ! ************************************************************************** !
  !> Preparation of the serialize PDF data
  !!
  !! Serialize the information from the (level-wise) state vector (scheme%state)
  !! into chunks for writing it in original treeIDlist order to disk (which is
  !! not sorted by levels, but by the space-filling curve)
  !! The data is stored like follows (e.g. varSys, 2 elements):
  !!
  !!```
  !!   Elem = 1     Elem = 2
  !! ---------------------------
  !! | stateVars | stateVars
  !! ---------------------------
  !!          varSys
  !!```
  !!
  subroutine mus_pdf_serialize( scheme, treeID, levelPointer, varMap, &
    &                            nElems, buffer,minLevel, maxLevel    )
    ! ---------------------------------------------------------------------- !
    !> scheme type containing the different state vectors
    type(mus_scheme_type), intent(in) :: scheme
    !> number of valid elements in this buffer
    integer, intent(in) :: nElems
    !> Partial treeID list
    integer(kind=long_k), intent(in) :: treeID(nElems)
    !> Partial Level pointer
    integer, intent(in) :: levelPointer(nElems)
    !> varaible map information
    type(tem_varMap_type), intent(in) :: varMap
    !> Data buffer
    real(kind=rk), intent(inout) :: buffer(:)
    integer :: minLevel, maxLevel
    ! ---------------------------------------------------------------------- !
    integer :: iVar, iElem, iLevel, iComp, iIndex
    integer :: nComp   ! max number of components
    integer :: nScalars, varPos, nSize, QQ, elemPos
    ! ---------------------------------------------------------------------- !

    ! Global counter for all words put into the buffer.
    iIndex = 0
    ! Run over all variable systems.
    ! @todo: for performance, loop order should be different for SOA or AOS

    QQ = scheme%layout%fStencil%QQ
    nScalars = varMap%nScalars

    ! For the current chunk (given from the caller), treat the corresponding
    ! elements from the tree ID list
    do iElem = 1, nElems
      ! On which level lies the current element
      iLevel = tem_LevelOf( treeID( iElem ))
      nSize  = scheme%pdf(iLevel)%nSize
      elemPos = levelPointer(iElem)

      ! nVars gives the number of non derived variables in the var system
      ! Run over the variables defined in the var system
      do iVar = 1, varMap%varPos%nVals
        varPos = varMap%varPos%val(iVar)
        nComp = scheme%varSys%method%val( varPos )%nComponents

        ! ... treat each component
        do iComp = 1, nComp
          iIndex = iIndex + 1
          buffer(iIndex) = scheme%state( iLevel )%val(                         &
& ( elempos-1)* nscalars+icomp+( ivar-1)* qq, &
            & scheme%pdf( iLevel )%nNext)
        end do  ! iComp

      end do  ! iVar
    end do  ! iElem

  end subroutine mus_pdf_serialize
  ! ************************************************************************** !


  ! ************************************************************************** !
  !> This subroutine unserializes the read data and stores it in the state-
  !! vector to perform a restart.
  !!
  subroutine mus_pdf_unserialize( scheme, treeID, levelPointer, varMap, &
    &                              nElems, buffer, minlevel, maxlevel)
    ! ---------------------------------------------------------------------- !
    !> scheme type containing the different state vectors
    type(mus_scheme_type), intent(inout) :: scheme
    integer, intent(in) :: nElems
    !> global tree data type
    integer(kind=long_k), intent(in) :: treeID(nElems)
    !> Level Pointer - from treelm to level descriptor
    integer,  intent(in) :: levelPointer(nElems)
    !> variable map
    type(tem_varMap_type), intent(in) :: varMap
    real(kind=rk), intent(in) :: buffer(:)
    integer :: minLevel, maxLevel
    ! ---------------------------------------------------------------------- !
    integer :: iVar, iElem, iLevel, iComp ! counter variables
    integer :: iIndex       ! amount of written variable components per varsys
    integer :: nComp        ! max number of components
    integer :: nScalars, varPos, nSize, QQ, elemPos
    ! ---------------------------------------------------------------------- !

    ! For all variable systems...
    iIndex = 0

    ! Actually this loop over the elems should be moved into and an
    ! array of treeIDs should be passed to the derive subroutine
    QQ = scheme%layout%fStencil%QQ
    nScalars = varMap%nScalars

    ! fill the state vector with the values from the chunk
    do iElem =  1, nElems

      iLevel = tem_LevelOf( treeID( iElem ))
      nSize  = scheme%pdf(iLevel)%nSize
      elemPos = levelPointer(iElem)

      ! nVars gives the number of state variables in the var system
      do iVar = 1, varMap%varPos%nVals
        varPos = varMap%varPos%val(iVar)
        nComp = scheme%varSys%method%val( varPos )%nComponents

        !NEC$ ivdep
        !IBM* INDEPENDENT
        !DIR$ IVDEP
        do iComp = 1, nComp
          iIndex = iIndex + 1
            ! MH: Only set nNext (is being switched to nNow before the compute
            ! kernel, which reads from nNow)
          scheme%state(iLevel)%val(                                  &
& ( elempos-1)* nscalars+icomp+( ivar-1)* qq, &
            & scheme%pdf(iLevel)%nNext ) = buffer(iIndex)

        end do  ! iComp
      end do  ! iVar
    end do  ! iElem

  end subroutine mus_pdf_unserialize
  ! ************************************************************************** !

end module mus_buffer_module
! **************************************************************************** !
