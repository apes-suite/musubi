! Copyright (c) 2014-2016, 2019-2020 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! ****************************************************************************** !
!> author: Kannan Masilamani
!! This module contains the type defintion which stores position of all
!! derive variables in global varSys.
!!
module mus_derVarPos_module
  use env_module,        only: rk
  use tem_aux_module,    only: tem_abort
  use tem_varSys_module, only: tem_varSys_type
  use tem_stencil_module, only: tem_stencilHeader_type

  use mus_scheme_layout_module, only: mus_scheme_layout_type
  use mus_scheme_derived_quantities_module, only: mus_scheme_derived_quantities_type

  implicit none

  private

  public :: mus_derVarPos_type
  public :: mus_derive_FromMacro_dummy
  public :: mus_derive_FromState_dummy
  public :: mus_derive_FromPreColState_dummy

  !> This type stores the position of each variable in the global sys
  type mus_derVarPos_type
    integer :: pdf               = -1
    integer :: fetch_pdf         = -1
    integer :: omega             = -1
    integer :: density           = -1
    integer :: moleDensity       = -1
    integer :: pressure          = -1
    integer :: kinePress         = -1
    integer :: velocity          = -1
    !integer :: grad_velocity     = -1
    ! all species velocities
    integer :: spc_velocities    = -1
    integer :: velMag            = -1
    integer :: momentum          = -1
    ! all species momentums
    integer :: spc_momenta       = -1
    integer :: shearStress       = -1
    integer :: wss               = -1
    integer :: shearMag          = -1
    integer :: strainRate        = -1
    integer :: shearRate         = -1
    integer :: kineticEnergy     = -1
    integer :: temperature       = -1
    integer :: moleFrac          = -1
    integer :: massFrac          = -1
    integer :: moleflux          = -1
    integer :: equilibrium       = -1
    integer :: nonEquilibrium    = -1
    integer :: equilibriumVel    = -1
    integer :: potential         = -1
    integer :: vol_frac          = -1
    procedure( derive_FromMacro ), pointer, nopass :: equilFromMacro => null()
    ! return velocity of single species
    procedure( derive_FromState ), pointer, nopass :: velFromState   => null()
    procedure( derive_FromState ), pointer, nopass :: equilFromState => null()
    ! return momentum of single species
    procedure( derive_FromState ), pointer, nopass :: momFromState => null()
    ! return velocity of all species
    procedure( derive_FromState ), pointer, nopass :: velocitiesFromState => null()
    ! return momentum of all species
    procedure( derive_FromState ), pointer, nopass :: momentaFromState => null()
    ! return velocity of single species from precollision state
    procedure( derive_FromPreColState ), pointer, nopass &
      & :: velFromPreColState => null()

    !> return auxField from local pdf state
    procedure(derive_auxFromState), pointer :: auxFieldFromState => null()
    ! KM: \todo use this function pointer in interpolation
    !> return equilibrium from auxilary variable for given nElems
    procedure(derive_equilFromAux), pointer, nopass :: equilFromAux => null()
    !> function to compute equilFromAux for single element
  end type mus_derVarPos_type

  abstract interface
    !> interface to derive equilibrium from macro.
    !! Mainly used in initial condition and boundary condition routines
    !! to avoid dublication of routines between different scheme kinds.
    !! In this interface, solver definition can be access via
    !! varSys%method%val(1)%method_data c_ptr
    subroutine derive_FromMacro( density, velocity, iField, nElems, varSys,    &
      &                          layout, res )
      import :: rk, tem_varSys_type, mus_scheme_layout_type
      !> Array of density.
      !! Single species: dens_1, dens_2 .. dens_n
      !! multispecies: dens_1_sp1, dens_1_sp2, dens_2_sp1, dens_2_sp2 ...
      !!                dens_n_sp1, dens_n_sp2
      !! Access: (iElem-1)*nFields + iField
      real(kind=rk), intent(in) :: density(:)

      !> Array of velocity.
      !! Size: (3, n*nFields)
      !! Access: ( iComp, (iElem-1)*nFields + iField )
      real(kind=rk), intent(in) :: velocity(:, :)

      !> Current field
      integer, intent(in) :: iField

      !> number of elements
      integer, intent(in) :: nElems

      !> variable system which is required to access fieldProp
      !! information via variable method data c_ptr
      type(tem_varSys_type), intent(in) :: varSys

      !> scheme layout contains stencil definition and lattice weights
      type(mus_scheme_layout_type), intent(in) :: layout

      !> Output of this routine
      !! Dimension: n*nComponents of res
      real(kind=rk), intent(out) :: res(:)
    end subroutine derive_FromMacro

    !> Interface that takes state array as input
    !! calculate density, velocity or eq as output
    !! State should be AOS layout
    subroutine derive_FromState( state, iField, nElems, varSys, layout, res )
      import :: rk, tem_varSys_type, mus_scheme_layout_type
      !> Array of state
      !! n * layout%stencil(1)%QQ * nFields
      real(kind=rk), intent(in) :: state(:)

      !> Current field
      integer, intent(in) :: iField

      !> number of elements
      integer, intent(in) :: nElems

      !> variable system which is required to access fieldProp
      !! information via variable method data c_ptr
      type(tem_varSys_type), intent(in) :: varSys

      !> scheme layout contains stencil definition and lattice weights
      type(mus_scheme_layout_type), intent(in) :: layout

      !> Output of this routine
      !! Dimension: n * nComponents of res
      !! Access: (iElem-1)*nComp + iComp
      !! To derive velocities of all species, dimension: n*nFields*nComp
      !! Access: (iElem-1)*nFields*nComp + (iField-1)*nComp + iComp
      real(kind=rk), intent(out) :: res(:)
    end subroutine derive_FromState

    !> Interface that takes state array as input
    !! calculate density, velocity or eq as output from FETCH state i.e.
    !! precollision state
    !! State should be AOS layout
    subroutine derive_FromPreColState( state, neigh, iField, nSize, nElems, &
      &                                varSys, layout, res )
      import :: rk, tem_varSys_type, mus_scheme_layout_type
      !> Array of state
      !! n * layout%stencil(1)%QQ * nFields
      real(kind=rk), intent(in) :: state(:)

      !> connectivity array
      integer, intent(in) :: neigh(:)

      !> Current field
      integer, intent(in) :: iField

      !> number of elements in state array
      integer, intent(in) :: nSize

      !> number of elements
      integer, intent(in) :: nElems

      !> variable system which is required to access fieldProp
      !! information via variable method data c_ptr
      type(tem_varSys_type), intent(in) :: varSys

      !> scheme layout contains stencil definition and lattice weights
      type(mus_scheme_layout_type), intent(in) :: layout

      !> Output of this routine
      !! Dimension: n * nComponents of res
      !! Access: (iElem-1)*nComp + iComp
      !! To derive velocities of all species, dimension: n*nFields*nComp
      !! Access: (iElem-1)*nFields*nComp + (iField-1)*nComp + iComp
      real(kind=rk), intent(out) :: res(:)
    end subroutine derive_FromPreColState

    !> Derive equilibrium from auxField for given nelems
    subroutine derive_equilFromAux( derVarPos, auxField, iField, nElems, &
      &                             varSys, layout, res )
      import :: rk, tem_varSys_type, mus_scheme_layout_type, mus_derVarPos_type
      !> Position of current field derive variable in variable system
      class(mus_derVarPos_type), intent(in) :: derVarPos
      !> Array of auxField.
      !! Single species: dens_1, vel_1, dens_2, vel_2, .. dens_n, vel_n
      !! multispecies: dens_1_sp1, vel_1_spc1, dens_1_sp2, vel_1_spc2,
      !!                dens_2_sp1, vel_2_spc2, dens_2_sp2, vel_2_spc2 ...
      !!                dens_n_sp1, vel_n_sp1, dens_n_sp2, vel_n_spc2
      !! Access: (iElem-1)*nFields + iField
      real(kind=rk), intent(in) :: auxField(:)

      !> Current field
      integer, intent(in) :: iField

      !> number of elements
      integer, intent(in) :: nElems

      !> variable system which is required to access fieldProp
      !! information via variable method data c_ptr
      type(tem_varSys_type), intent(in) :: varSys

      !> scheme layout contains stencil definition and lattice weights
      type(mus_scheme_layout_type), intent(in) :: layout

      !> Output of this routine
      !! Dimension: n*QQ of res
      real(kind=rk), intent(out) :: res(:)
    end subroutine derive_equilFromAux

    !> Derive auxField from local state.
    !! \todo KM: pass external force to add to auxField
    subroutine derive_auxFromState( derVarPos, state, neigh, iField, nElems, &
      &                             nSize, iLevel, stencil, varSys, auxField, quantities )
      import :: rk, tem_varSys_type, mus_derVarPos_type, tem_stencilHeader_type, &
        & mus_scheme_derived_quantities_type
      !> Position of current field derive variable in variable system
      class(mus_derVarPos_type), intent(in) :: derVarPos
      !> Array of state
      !! nSize * layout%stencil(1)%QQ * nFields
      !! use IDX macro to access this state
      real(kind=rk), intent(in) :: state(:)

      !> connectivity vector
      integer, intent(in) :: neigh(:)

      !> Current field
      integer, intent(in) :: iField

      !> number of elements to compute
      integer, intent(in) :: nElems

      !> number of elements in state array
      integer, intent(in) :: nSize

      !> current level
      integer, intent(in) :: iLevel

      !> stencil header contains discrete velocity vectors
      type(tem_stencilHeader_type), intent(in) :: stencil

      !> variable system which is required to access fieldProp
      !! information via variable method data c_ptr
      type(tem_varSys_type), intent(in) :: varSys

      !> Class that contains pointers to the proper derived quantities functions
      type(mus_scheme_derived_quantities_type), intent(in) :: quantities

      !> Output of this routine
      !! auxField is inout to allow storing auxField for each species
      !! seperately
      !! Size: nElems*nAuxScalars
      real(kind=rk), intent(inout) :: auxField(:)
    end subroutine derive_auxFromState

    !> Derive equilibrium from auxField for single element
    pure function derive_equilFromAuxFunc( derVarPos, auxField, iField, &
      &                                    varSys, layout) result( res )
      import :: rk, tem_varSys_type, mus_scheme_layout_type, mus_derVarPos_type
      !> Position of current field derive variable in variable system
      class(mus_derVarPos_type), intent(in) :: derVarPos
      !> Array of auxField of single element.
      !! Single species: dens_1, vel_1
      !! multispecies: dens_1_sp1, vel_1_spc1, dens_1_sp2, vel_1_spc2,
      !! Access: (iElem-1)*nFields + iField
      real(kind=rk), intent(in) :: auxField(:)

      !> Current field
      integer, intent(in) :: iField

      !> variable system which is required to access fieldProp
      !! information via variable method data c_ptr
      type(tem_varSys_type), intent(in) :: varSys

      !> scheme layout contains stencil definition and lattice weights
      type(mus_scheme_layout_type), intent(in) :: layout

      !> Output of this routine
      !! Dimension: QQ of res
      real(kind=rk) :: res(layout%fStencil%QQ)
    end function derive_equilFromAuxFunc

    !> Derive auxField from state for single element
    pure function derive_auxFromStateFunc( derVarPos, state, iField, stencil, &
      &                                    varSys ) result( res )
      import :: rk, tem_varSys_type, mus_derVarPos_type, tem_stencilHeader_type
      !> Position of derive variable in variable system
      class(mus_derVarPos_type), intent(in) :: derVarPos
      !> Array of state
      !! layout%stencil(1)%QQ * nFields
      real(kind=rk), intent(in) :: state(:)

      !> Current field
      integer, intent(in) :: iField

      !> stencil header contains discrete velocity vectors
      type(tem_stencilHeader_type), intent(in) :: stencil

      !> variable system which is required to access fieldProp
      !! information via variable method data c_ptr
      type(tem_varSys_type), intent(in) :: varSys

      !> Output of this routine
      !! Size: nAuxScalars
      real(kind=rk) :: res(varSys%nAuxScalars)
    end function derive_auxFromStateFunc
  end interface

contains
  ! ************************************************************************** !
  subroutine mus_derive_FromMacro_dummy( density, velocity, iField, nElems, &
    &                                    varSys, layout, res )
    !> Array of density.
    !! Single species: dens_1, dens_2 .. dens_n
    !! multispecies: dens_1_sp1, dens_1_sp2, dens_2_sp1, dens_2_sp2 ...
    !!                dens_n_sp1, dens_n_sp2
    !! Access: (iElem-1)*nFields + iField
    real(kind=rk), intent(in) :: density(:)

    !> Array of velocity.
    !! Size: (3, n*nFields)
    !! Access: ( iComp, (iElem-1)*nFields + iField )
    real(kind=rk), intent(in) :: velocity(:, :)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n*nComponents of res
    real(kind=rk), intent(out) :: res(:)

    call tem_abort('Dummy routine for derive_FromMacro')
  end subroutine mus_derive_FromMacro_dummy
  ! ************************************************************************** !

  ! ************************************************************************** !
  ! Interface that takes state array as input
  ! calculate density, velocity or eq as output
  ! State should be AOS layout
  subroutine mus_derive_FromState_dummy(state, iField, nElems, varSys, layout, &
    &                                   res)
    !> Array of state
    !! n * layout%stencil(1)%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    !! Access: (iElem-1)*nComp + iComp
    !! To derive velocities of all species, dimension: n*nFields*nComp
    !! Access: (iElem-1)*nFields*nComp + (iField-1)*nComp + iComp
    real(kind=rk), intent(out) :: res(:)

    call tem_abort('Dummy routine for derive_FromState')
  end subroutine mus_derive_FromState_dummy
  ! ************************************************************************** !

  ! ************************************************************************** !
  ! Interface that takes state array as input
  ! calculate density, velocity or eq as output from FETCH state i.e.
  ! precollision state
  ! State should be AOS layout
  subroutine mus_derive_FromPreColState_dummy( state, neigh, iField, nSize, &
    &                                          nElems, varSys, layout, res )
    !> Array of state
    !! n * layout%stencil(1)%QQ * nFields
    real(kind=rk), intent(in) :: state(:)

    !> connectivity array
    integer, intent(in) :: neigh(:)

    !> Current field
    integer, intent(in) :: iField

    !> number of elements in state array
    integer, intent(in) :: nSize

    !> number of elements
    integer, intent(in) :: nElems

    !> variable system which is required to access fieldProp
    !! information via variable method data c_ptr
    type(tem_varSys_type), intent(in) :: varSys

    !> scheme layout contains stencil definition and lattice weights
    type(mus_scheme_layout_type), intent(in) :: layout

    !> Output of this routine
    !! Dimension: n * nComponents of res
    !! Access: (iElem-1)*nComp + iComp
    !! To derive velocities of all species, dimension: n*nFields*nComp
    !! Access: (iElem-1)*nFields*nComp + (iField-1)*nComp + iComp
    real(kind=rk), intent(out) :: res(:)

    call tem_abort('Dummy routine for derive_FromPreColState')
  end subroutine mus_derive_FromPreColState_dummy
  ! ************************************************************************** !


end module mus_derVarPos_module

