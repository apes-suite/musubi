! Copyright (c) 2019 Seyfettin Bilgi <seyfettin.bilgi@student.uni-siegen.de>
! Copyright (c) 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! Make sure loglvl is defined to an integer value.
! Usually this should be defined on the command line with -Dloglvl=
! to some value.
! ****************************************************************************** !
!> This module provides the definition and methods for
!! BGK advection relaxation scheme.
module mus_compute_nernstPlanck_module
  use iso_c_binding, only: c_f_pointer
  ! include treelm modules
  use env_module,               only: rk
  use tem_varSys_module,        only: tem_varSys_type !, tem_varSys_op_type
  use tem_param_module,         only: div1_3, div1_36, div1_8, div3_4h, div1_4,&
    &                                 div3_8, div9_16, div3_16, cs2inv, cs4inv,&
    &                                 rho0

  ! include musubi modules
  use mus_field_prop_module,    only: mus_field_prop_type
  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_scheme_type_module,   only: mus_scheme_type
  use mus_param_module,         only: mus_param_type
  use mus_varSys_module,        only: mus_varSys_data_type
  use mus_derVarPos_module,     only: mus_derVarPos_type

  implicit none

  private

  public :: mus_nernstPlanck_advRel_generic

contains

! ****************************************************************************** !
  !> Advection relaxation routine for the
  !! nernst planvk model with an explicit calculation of all equilibrium
  !! quantities. Slow and simple. This routine should only be
  !! used for testing purposes
  !!
  !! This subroutine interface must match the abstract interface definition
  !! [[kernel]] in scheme/[[mus_scheme_type_module]].f90 in order to be callable
  !! via [[mus_scheme_type:compute]] function pointer.
  subroutine mus_nernstPlanck_advRel_generic( fieldProp, inState, outState,    &
    &                                         auxField, neigh, nElems, nSolve, &
    &                                         level, layout, params, varSys,   &
    &                                         derVarPos                        )
    ! -------------------------------------------------------------------- !
    !> Array of field properties (fluid or species)
    type(mus_field_prop_type), intent(in) :: fieldProp(:)
    !> variable system definition
    type(tem_varSys_type), intent(in) :: varSys
    !> current layout
    type(mus_scheme_layout_type), intent(in) :: layout
    !> number of elements in state Array
    integer, intent(in) :: nElems
    !> input  pdf vector
    real(kind=rk), intent(in)  ::  inState(nElems * varSys%nScalars)
    !> output pdf vector
    real(kind=rk), intent(out) :: outState(nElems * varSys%nScalars)
    !> Auxiliary field computed from pre-collision state
    !! Is updated with correct velocity field for multicomponent models
    real(kind=rk), intent(inout) :: auxField(nElems * varSys%nAuxScalars)
    !> connectivity vector
    integer, intent(in) :: neigh(nElems * layout%fStencil%QQ)
    !> number of elements solved in kernel
    integer, intent(in) :: nSolve
    !> current level
    integer,intent(in) :: level
    !> global parameters
    type(mus_param_type),intent(in) :: params
    !> position of derived quantities in varsys for all fields
    type( mus_derVarPos_type ), intent(in) :: derVarPos(:)
    ! -------------------------------------------------------------------- !
    integer :: iElem, iDir, iField                       ! voxel element counter
    integer :: QQ, nScalars, nFields, elemOff
    ! temporary distribution variables
    real(kind=rk) pdfTmp
    real(kind=rk) inv_vel,transVel( nsolve*3 ), vel_fluid(3)   ! local velocit
    integer :: vel_varPos
    real(kind=rk) moleDens     ! local density
    ! derived constants
    ! equilibrium calculation variables
    real(kind=rk) ucx(layout%fStencil%QQ)
    real(kind=rk) eqState
    real(kind=rk) omega
    type(mus_varSys_data_type), pointer :: fPtr
    type(mus_scheme_type), pointer :: scheme
    ! ---------------------------------------------------------------------------
    QQ = layout%fStencil%QQ
    ! nElems = size(neigh)/QQ
    nScalars = varSys%nScalars

    ! access scheme via 1st variable method data which is a state variable
    call C_F_POINTER( varSys%method%val(1)%method_Data, fPtr )
    scheme => fPtr%solverData%scheme
    nFields = scheme%nFields

    ! passive scalar has only one transport Variable
    vel_varPos = scheme%transVar%method(1)%data_varPos
    ! Get velocity field
    call varSys%method%val(vel_varPos)%get_valOfIndex( &
      & varSys  = varSys,                              &
      & time    = params%general%simControl%now,       &
      & iLevel  = level,                               &
      & idx     = scheme%transVar%method(1)            &
      &           %pntIndex%indexLvl(level)            &
      &           %val(1:nSolve),                      &
      & nVals   = nSolve,                              &
      & res     = transVel                             )

    ! convert physical velocity into LB velocity
    inv_vel = 1.0_rk / params%physics%fac( level )%vel
    transVel = transVel * inv_vel


    nodeloop: do iElem=1, nSolve
      elemOff = (iElem-1)*varSys%nAuxScalars

      ! x-, y- and z-velocity from transport field
      vel_fluid = transVel( (iElem-1)*3+1 : iElem*3 )

      ! since both species are transported by same velocity field.
      ! this can be computed just one
      do iDir = 1, QQ

        !> Pre-calculate velocitiy terms
        ucx(iDir) = dble(layout%fStencil%cxDirRK( 1, iDir ))*vel_fluid(1) &
          &       + dble(layout%fStencil%cxDirRK( 2, iDir ))*vel_fluid(2) &
          &       + dble(layout%fStencil%cxDirRK( 3, iDir ))*vel_fluid(3)
      end do

      fieldloop: do iField = 1, nFields
        !> relaxation parameter
        omega = fieldProp(iField)%species%omega

        ! mole density
        moleDens = auxField(elemOff + iField)

        do iDir = 1, QQ

          !> Calculate equilibrium distribution functions fEq
          eqState = layout%weight(iDir)*moleDens*( 1.0_rk + ucx(iDir)*cs2inv )

          !> Generic fetching step:
          !! Streaming for pull
          !! Local copy for push
          pdfTmp = inState(                                            &
            & neigh(( idir-1)* nelems+ ielem)+( ifield-1)* qq+ nscalars*0 )

          !> Relaxation
          outState( ( ielem-1)* nscalars+ idir+( ifield-1)* qq) &
            & = pdfTmp - omega*( pdfTmp - eqState )
        end do  !< iDir

      end do fieldloop
    end do nodeloop

  end subroutine mus_nernstPlanck_advRel_generic
! ****************************************************************************** !

end module mus_compute_nernstPlanck_module
! ****************************************************************************** !

