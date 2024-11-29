! Copyright (c) 2012-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2013, 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! **************************************************************************** !
!> author: Manuel Hasert
!! Collection of debugging functions for writing stuff to files
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
module mus_debug_tools_module

  ! include treelm modules
  use env_module,          only: long_k, rk
  use tem_element_module,  only: eT_fluid, eT_ghostFromFiner,          &
    &                            eT_ghostFromCoarser, eT_halo
  use tem_construction_module, only: tem_levelDesc_type
  use tem_debug_module,    only: dbgUnit

  ! include musubi modules
  use mus_pdf_module,      only: pdf_data_type

  implicit none

  private

  public :: dumpPdfAll


contains

! **************************************************************************** !
  !> write all pdf entries of a all elements of the pdf array for a given level
  !! to the debug unit
  !!
  subroutine dumpPdfAll( pdf, levelDesc, QQ, text, iTime, nFields,      &
    &                    nScalars, dumpHalos, state )
    ! --------------------------------------------------------------------------
    type( pdf_data_type ), intent(in) :: pdf
    type( tem_levelDesc_type ), intent(in) :: levelDesc
    real(kind=rk), intent(in) :: state(:,:)
    integer, intent(in) :: QQ
    logical, optional, intent(in)   :: dumpHalos
    character(len=*) :: text
    integer, intent(in), optional :: iTime      !< Level counter variable
    integer, intent(in) :: nFields !< number of fields
    integer, intent(in) :: nScalars !< total number of scalars in state array
    ! --------------------------------------------------------------------------
    integer :: iElem, elemPos, iType, nElems, iDir, iTLayer, iField
    integer :: nTypes
    character(len=1024) :: buffer
    logical :: loc_dumpHalos
    ! --------------------------------------------------------------------------

    nElems = 0
    iTLayer = pdf%nNext
    if( present( iTime ) ) then
      iTLayer = iTime
    endif
    if( present( dumpHalos ) ) then
      loc_dumpHalos = dumpHalos
    else
      loc_dumpHalos = .true.
    endif
    if( loc_dumpHalos ) then
      nTypes = 4
    else
      nTypes = 3
    end if

    write(dbgUnit(1),*) trim( text )
    do iType = 1, nTypes
      select case( iType )
      case( eT_fluid )
        nElems = pdf%nElems_Fluid
        write(dbgUnit(1),*) 'FLUID'
      case( eT_ghostFromCoarser)
        nElems = pdf%nElems_ghostFromCoarser
        write(dbgUnit(1),*) 'GHOSTFROMCOARSER'
      case( eT_ghostFromFiner)
        nElems = pdf%nElems_ghostFromFiner
        write(dbgUnit(1),*) 'GHOSTFROMFINER'
      case( eT_halo)
        nElems = pdf%nElems_halo
        write(dbgUnit(1),*) 'HALO'
      end select

      do iElem = 1, nElems
        elemPos = iElem + levelDesc%offset( 1, iType )
        write(buffer,'(i8)')  levelDesc%total( elemPos )
        do iField = 1, nFields
          do iDir = 1, QQ
            write(buffer,'(2a,en17.7)') trim(buffer),' ',state(           &
& ( elempos-1)* nscalars+ idir+( ifield-1)* qq, iTLayer )
          end do ! QQ
          write(dbgUnit(1),*) trim( buffer )
        end do
      end do
    end do
  end subroutine dumpPdfAll
! **************************************************************************** !


end module mus_debug_tools_module
! **************************************************************************** !
