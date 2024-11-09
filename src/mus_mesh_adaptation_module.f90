! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Julia Moos <julia.moos@student.uni-siegen.de>
! Copyright (c) 2015-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!! author: Kartik Jain
!> This module contains the code responsible for adaptively refining the mesh
!! during run time.
!!
module mus_mesh_adaptation_module

  ! include treelm modules
  use mpi
  use env_module,               only: rk, long_k
  use tem_adaptation_module,    only: tem_adapt_dump_newMesh
  use tem_comm_env_module,      only: tem_comm_env_type
  use tem_adaptation_config_module, only: tem_adapt_type

  ! include musubi modules
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_geom_module,          only: mus_geom_type

  implicit none

  private

  public :: mus_adapt_refine


contains


! **************************************************************************** !
  !> Wrap up the routines required for dynamic load balancing
  subroutine mus_adapt_refine( geometry, scheme, proc, adapt )
    ! --------------------------------------------------------------------------
    !> Treelmesh data
    type( mus_geom_type )       , intent(inout)     :: geometry
    !> scheme type
    type( mus_scheme_type )     , intent(inout)     :: scheme
    type(tem_comm_env_type)     , intent(in)        :: proc
    !> mesh adaptation
    type(tem_adapt_type)        ,  intent(in)       :: adapt
    ! --------------------------------------------------------------------------

    ! Now you get the adaptively refined mesh back
    ! Dump it to the disk at the current moment
    call tem_adapt_dump_newMesh( levelDesc    = scheme%levelDesc, &
      &                          tree         = geometry%tree,    &
      &                          proc         = proc              )

  end subroutine mus_adapt_refine

end module mus_mesh_adaptation_module
