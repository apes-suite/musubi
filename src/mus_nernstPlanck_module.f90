! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
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
! **********************************************************************!
!! In this module the nernst_planck equation is considered.
!! The file type contains all relevant information.
!! It is a one dimensional case.
!> author: seyfettin bilgi
module mus_nernstPlanck_module
  use env_module,               only: rk
  use tem_aux_module,           only: tem_abort
  use tem_logging_module,       only: logUnit, tem_toStr

  use aotus_module,             only: flu_State, aot_get_val, aoterr_Fatal, &
    &                                 aoterr_NonExistent, aoterr_WrongType
  use aot_table_module,         only: aot_table_open, aot_table_close, &
    &                                 aot_table_length
  use aot_out_module,           only: aot_out_type, aot_out_val, &
    &                                 aot_out_open_table, &
    &                                 aot_out_close_table

  use mus_physics_module,       only: mus_physics_type, faraday, gasConst_R
  implicit none

  private

  public :: mus_nernstPlanck_type
  public :: mus_load_nernstPlanck

  !> Contains configuration to calculate nernst_planck equation.
  !! More information can be found in
  !! A Coupled Lattice Boltzmann Method to Solve Nernst_Planck Model for
  !! simulating Electro-Osmatic Flows by Yang x., Shi B., Chai Z., Guo Z.
  type mus_nernstPlanck_type
    !> abosulte temperature in Kelvin
    real(kind=rk) :: temp
    !> Molar density of ions
    real(kind=rk) :: moleDens
  end type mus_nernstPlanck_type


contains


  ! ****************************************************************************
  !> load input to solve nernst_planck equation
  subroutine mus_load_nernstPlanck( me, conf, parent, physics )
    !-------------------------------------------------------------------------
    !> nernst_planck type
    type(mus_nernstPlanck_type), intent(out) :: me
    !> flu state
    type(flu_state) :: conf
    !> parent handle
    integer, intent(in), optional :: parent
    !> physics type to convert physics to lattice unit or vice versa
    type(mus_physics_type), intent(in) :: physics
    !-------------------------------------------------------------------------
    integer :: NP_handle
    integer :: iError
    !-------------------------------------------------------------------------

    ! if nernst planck informations in scheme table parentHandle /= 0
    if ( present(parent) ) then
      call aot_table_open( L       = conf,           &
        &                  parent  = parent,         &
        &                  thandle = NP_handle,      &
        &                  key     = 'nernst_planck' )
    else
      call aot_table_open( L=conf, thandle = NP_handle, &
        &                  key = 'nernst_planck'        )
    end if

    if ( NP_handle == 0 ) then
       write(logUnit(1),*)'No nernst_planck table defined'
       call tem_abort()
    end if

    write(logUnit(1),*) 'Loading nernst planck informations'

    ! load absolute temperature
    call aot_get_val( L       = conf,           &
      &               thandle = NP_handle,      &
      &               key     = 'mole_density', &
      &               val     = me%moleDens,    &
      &               ErrCode = iError          )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(3),*)'WARNING: FATAL Error occured, while retrieving '//   &
        &                'mole_density:'
      if (btest(iError, aoterr_NonExistent))                            &
        & write(logUnit(1),*)'Variable not existent!'
      if (btest(iError, aoterr_WrongType))                              &
        & write(logUnit(1),*)'Variable has wrong type!'
        write(logUnit(1),*)'STOPPING'
        call tem_abort()
    end if
    me%moleDens = me%moleDens / physics%moleDens0

    ! load absolute temperature
    call aot_get_val( L       = conf,           &
      &               thandle = NP_handle,      &
      &               key     = 'temp',         &
      &               val     = me%temp, &
      &               default = 273._rk,        &
      &               ErrCode = iError          )
    if (btest(iError, aoterr_Fatal)) then
      write(logUnit(3),*)'WARNING: FATAL Error occured, while retrieving '//   &
        &                'temp:'
      write(logUnit(3),*) 'WARNING: Using default value 273 K'
    end if
    me%temp = me%temp / physics%temp0

    call aot_table_close( L=conf, thandle=NP_handle )

    write(logUnit(1),"(A)") 'Nernst_Planck properties:'
    write(logUnit(1),"(A)") '  mole_density of solvent: ' &
      &                     // trim(tem_toStr(me%moleDens))
    write(logUnit(1),"(A)") '  temp: '//trim(tem_toStr(me%temp))

  end subroutine mus_load_nernstPlanck
  ! ****************************************************************************

end module mus_nernstPlanck_module
