! Copyright (c) 2021 Harald Klimach <harald.klimach@uni-siegen.de>
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
!
!> Musubi specific criteria to abort the simulation.
!>
!! This module provies the [[mus_abortCriteria_type]] that holds parameters for
!! simulation aborts.
!! The configuration parameters are read from the `abort_criteria` table in
!! `sim_control`, see also [[tem_simControl_module]].
!!
!! Additional parameters for Musubi in the `abort_criteria` table:
!!
!! * velocity_lat_max: The maximal lattice velocity that may appear in the
!!                     domain. Defaults to 0.15.
!!                     Lattice velocities larger than 0.1 are considered to
!!                     produce large errors.
!!
!! A simple complete example for the Musubi `abort_criteria` table within the
!! `sim_control` table is:
!!
!!```lua
!!    abort_criteria = {
!!      stop_file = 'stop',
!!      velocity_lat_max = 0.15 -- Maximum lattice velocity for Musubi
!!    }
!!```
!!
module mus_abortCriteria_module
  use aotus_module, only: flu_state, aot_get_val

  use env_module, only: rk
  use tem_abortCriteria_module, only: tem_solverAborts_type

  implicit none

  private

  !> Musubi specific abort criteria.
  type, public, extends(tem_solverAborts_type) :: mus_abortCriteria_type

    !> Maximal lattice velocity that will be tolerated in the simulation.
    !! The lattice velocity should usually be smaller than 0.1.
    real(kind=rk) :: velLat_max = 0.15_rk

  contains

    procedure :: load => mus_abortCriteria_load

  end type mus_abortCriteria_type


contains


  ! ------------------------------------------------------------------------ !
  !> Loading Musubi specific abort criteria from the `abort_criteria` table.
  !!
  !! This routine is passed to the loading of abort criteria in
  !! [[tem_abortcriteria_module:tem_abortcriteria_load]].
  subroutine mus_abortCriteria_load(me, conf, abort_table)
    ! -------------------------------------------------------------------- !
    !> Object to hold the solver specific configuration parameters.
    class(mus_abortCriteria_type), intent(inout) :: me

    !> Handle to the Lua script with the configuration.
    type(flu_state), intent(in) :: conf

    !> Handle to the opened `abort_criteria` table that holds the
    !! abort parameters to load.
    integer, intent(in) :: abort_table
    ! -------------------------------------------------------------------- !
    integer :: iErr
    ! -------------------------------------------------------------------- !

    call aot_get_val( L       = conf,               &
      &               thandle = abort_table,        &
      &               val     = me%velLat_max,      &
      &               key     = 'velocity_lat_max', &
      &               default = 0.15_rk,            &
      &               ErrCode = iErr                )

  end subroutine mus_abortCriteria_load
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module mus_abortCriteria_module
