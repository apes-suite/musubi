! Copyright (c) 2013,2024 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2013-2020,2024 Kannan Masilamani <kannan.masilamani@dlr.de>
! Copyright (c) 2013-2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Sindhuja Budaraju <nagasai.budaraju@student.uni-siegen.de>
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
!!
!! This module contains data type and modules related to musubi
!! lattice to physical unit convertion and vice versa.
!! physics data type is global for all schemes, it is defined
!! in the following format:
!!
!!```lua
!! physics = { dt = dt_phy, -- physical time step size
!!             rho0 = rho0_phy, -- reference density
!!             temp0 = t_phy -- reference temperature
!!}
!!```
!!
!! Of these quantities, dt is mandetory for conversion.
!! Others can be omitted, thus use default values.
!!
!! To add a new conversion factor, one has to do the following:
!!
!! 1. add this new factor into [mus_convertFac_type]
!! 2. add the defination of factor inside routine [mus_set_convFac]
!! 3. add this new factor into routine [mus_physics_out]
!!
module mus_physics_module

  ! include treelm modules
  use env_module,          only: rk
  use tem_param_module,    only: cs
  use tem_aux_module,      only: tem_abort
  use treelmesh_module,    only: treelmesh_type
  use tem_geometry_module, only: tem_ElemSizeLevel
  use tem_logging_module,  only: logUnit, tem_toStr
  use tem_tools_module,    only: tem_horizontalSpacer

  ! include aotus modules
  use aotus_module,     only: flu_State, aot_get_val, aoterr_Fatal, &
    &                         aoterr_NonExistent, aoterr_WrongType
  use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val
  use aot_out_module,   only: aot_out_type, aot_out_open_table, &
    &                         aot_out_close_table, aot_out_val

  implicit none

  private

  !> Reference coulomb is set to fundamental electrical charge
  real(kind=rk), parameter, public :: coulomb_ref = 1.60217657e-19_rk

  !> Reference mole is set to inverse of Avogadro's constant
  real(kind=rk), parameter, public :: mole_ref = 1e-23_rk / 6.02214129

  !> the boltzmann constant J K^-1
  real(kind=rk), parameter, public :: k_b = 1.38064852e-23_rk

  !> Faraday constant C/mol
  real(kind=rk), parameter, public :: faraday = 96485.3365_rk

  !> Ideal gas constant N m / (mol K)
  real(kind=rk), parameter, public :: gasConst_R = 8.3144621_rk

  public :: mus_physics_type
  public :: mus_load_physics
  public :: mus_create_funcStr
  public :: mus_physics_out
  public :: mus_physics_out_conv
  public :: mus_convertFac_type
  public :: mus_set_convFac
  public :: set_values_by_levels
  public :: mus_physics_dump2outUnit
  public :: mus_set_scaleFac

  !> This type contains the converstion factor for derived variables
  !! from lattice to physical.
  !!
  !! Their inverses can be used to convert physical to lattice units
  !! use reference density to parmeterize kg and reference mole density
  !! to parmeterize mol
  type mus_convertFac_type
    !> length (m) = dx
    real(kind=rk) :: length
    !> time (s) = dt
    real(kind=rk) :: time
    !> velocity(m/s) = dx/dt
    real(kind=rk) :: vel
    !> kinematic viscosity(m^2/s) = dx^2/dt
    real(kind=rk) :: visc
    !> Dynamic viscosity (Pa s) = kg/m/s
    real(kind=rk) :: viscDyna
    !> acceleration(m/s^2) = dx/dt^2
    real(kind=rk) :: accel
    !> Force(N)(kg m/s^2) = rho0*dx^4/dt^2
    real(kind=rk) :: force
    !> Force per unit volume (N/m^3)(kg/s^2/m^2) = rho0*dx/dt^2
    real(kind=rk) :: body_force
    !> Pressure(N/m^2)(kg/m/s^2) = rho0*dx^2/dt^2
    real(kind=rk) :: press
    !> Strain Rate (1/s) = 1/dt
    real(kind=rk) :: strainRate
    !> Energy (N-m) (kg*m^2/s^2) = rho0*dx^5/dt^2
    real(kind=rk) :: energy
    !> mole density(mol/m^3) = mole0/dx^3
    !real(kind=rk) :: moleDens
    !> Charge density (C/m^3) = Coulomb0/dx^3
    real(kind=rk) :: chargeDens
    !> Current density (C/s/m^2) = Coulomb0/dt/dx^2
    real(kind=rk) :: currentDens
    !> mole flux(mol/m^2/s) = moleDen0*dx/dt
    real(kind=rk) :: moleFlux
    !> mass flux(kg/m^2/s) = rho0*dx/dt
    real(kind=rk) :: flux
    !> diffusivity(m^2/s) = dx^2/dt
    real(kind=rk) :: diffusivity
    !> faraday (C/mol) = coulomb0/moleDens0/dx^3
    real(kind=rk) :: faraday
    !> gas constant (J/mol/K) (N m/mol/K) (kg m^2/s^2/mol/K)
    !! = rho0*dx^5/dt^2/mole0/temp0
    real(kind=rk) :: gasConst
    !> Potential (V) (kg m^2/(C*S^2))
    real(kind=rk) :: potential
    !> sourceCoeff (1/s) = 1/dt
    real(kind=rk) :: sourceCoeff
  end type mus_convertFac_type

  !> This type contains the reference values as defined in the physics
  !! table by the user
  !!
  !! They are used to convert lattice to physical unit and vice versa.
  !! keep reference mass density, mole density and molecular weight
  !! same for all levels.
  !! Use [mus_load_physics] to fill the datatype properly.
  type mus_physics_type
    !> needed to check if physics table is defined
    logical :: active = .false.
    !> reference length - discretization size of the coarsest level
    !! SI unit - meter
    real(kind=rk) :: dx = -1.0_rk
    real(kind=rk), allocatable :: dxLvl(:)
    !> reference time - time discretization for discretization size
    !! of the coarsest level
    !! SI unit - seconds
    real(kind=rk) :: dt = -1.0_rk
    real(kind=rk), allocatable :: dtLvl(:)
    !> reference physical mass density
    !! SI unit - kg/m^3
    real(kind=rk) :: rho0 = -1.0_rk
    !> reference physical mole density
    !! SI unit - mol/m^3
    real(kind=rk) :: moleDens0 = -1.0_rk
    !> reference molecular weight
    !! SI unit - kg/mol
    real(kind=rk) :: molWeight0 = -1.0_rk
    !> reference temperature
    !! SI unit - Kelvin
    real(kind=rk) :: temp0 = -1.0_rk
    !> reference fundamental electrical charge
    !! SI unit - Coulomb
    real(kind=rk) :: coulomb0 = -1.0_rk
    !> mole is defined by inverse of Avogadro Constant
    !! Avogadro Constant = 6.02214129e23 [1/mol]
    real(kind=rk) :: mole0 = -1.0_rk
    !> reference mass in kg derived from density or moleweight
    !! SI unit :: kg
    real(kind=rk) :: mass0 = -1.0_rk

    !> Level-wise conversion factor for derived variables
    !! size: minLevel:maxLevel
    !! allocated in mus_load_physics
    !! \todo KM: conversion factor should not be level-dependent.
    !! it should be same for all levels, the lattice dx and dt for each level
    !! must be considered to scale variables in multilevel.
    !! Implemented force, visc, etc using dtL according to formula
    !! Introduced lattice speed variable: dx/dt for each level
    !! it should be same for all level for acoustic scaling and different
    !! for diffusive scaling
    type( mus_convertFac_type ), allocatable :: fac(:)

    !> Pressure (strain rate) over level scale factor.
    !! This factor is meant to convert pressure in LB unit on source level to
    !! the required pressure on target level.
    !! It is mainly used in interpolation routine.
    !! It is allocated as: allocate(pFac( minLevel:maxLevel, minLevel:maxLevel))
    !! It is allocated and initialized in routine: mus_set_scaleFac
    !! How to use it in the code:
    !! pTargetLevel = pSourceLevel * pFac( sourceLevel, targetLevel )
    real(kind=rk), allocatable :: pFac(:,:)

    !> Velocity over level scale factor
    !! Its usage is the same as pressure scale factor
    real(kind=rk), allocatable :: vFac(:,:)

    !> Strain rate over level scale factor
    !! Its usage is the same as pressure scale factor
    real(kind=rk), allocatable :: sFac(:,:)

  end type mus_physics_type


contains


  ! ************************************************************************** !
  !> This routine loads the physics table from musubi config file
  !!
  !! If no physics table is provided, the conversion factors default to
  !! 1, resulting in the lattice units being directly used i.e.
  !! dx and dt are set default to 1.
  !! See the [mus_physics_type] for a description of the various factors that
  !! can be set here.
  subroutine mus_load_physics( me, conf, tree, scaleFactor )
    ! --------------------------------------------------------------------------
    !> physics type
    type( mus_physics_type ), intent(out) :: me
    !> flu state
    type( flu_State ) :: conf
    !> global treelm mesh
    type( treelmesh_type), intent(in) :: tree
    !> scaling factor: diffusive -> 4; acoustic -> 2
    integer, intent(in) :: scaleFactor
    ! --------------------------------------------------------------------------
    integer :: thandle
    integer :: iError
    real(kind=rk) :: cs_phy ! physical speed of sound
    ! --------------------------------------------------------------------------
    call tem_horizontalSpacer(fUnit = logUnit(1))
    write(logUnit(1),*) ' Loading physics table ...'
    call aot_table_open( L=conf, thandle=thandle, key='physics' )

    if (thandle > 0) then
      me%active = .true.
      ! Found a physics table.
      ! Activate the conversion
      write(logUnit(1), *) ' Physics table is defined.'
      write(logUnit(5), *) ' All values are expressed in physical units.'

      ! reference dx is always the coarsest level in the tree
      me%dx = tem_ElemSizeLevel( tree, tree%global%minLevel )

      ! load cs
      call aot_get_val( L = conf, thandle = thandle, key = 'cs', &
        &               val = cs_phy, ErrCode = iError           )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(7),*)'FATAL Error occured, while retrieving cs.'
        if (btest(iError, aoterr_WrongType)) then
          write(logUnit(1),*)'Variable has wrong type!'
          write(logUnit(1),*)'STOPPING'
          call tem_abort()
        endif
      end if

      if (btest(iError, aoterr_NonExistent)) then
        write(logUnit(1),*) 'Speed of sound (cs) is not defined.'
        write(logUnit(1),*) 'Attempting to load time step dt on coarsest level:'

        ! load dt
        call aot_get_val( L = conf, thandle = thandle, key = 'dt', &
          &               val = me%dt, ErrCode = iError            )
        if (btest(iError, aoterr_Fatal)) then
          write(logUnit(7),*)'FATAL Error occured, while retrieving dt.'
          if (btest(iError, aoterr_WrongType)) then
            write(logUnit(1),*)'Variable has wrong type!'
            write(logUnit(1),*)'STOPPING'
            call tem_abort()
          endif
        end if
      
        ! load time step size if speed of sound is not defined
        if (btest(iError, aoterr_NonExistent)) then
          write(logUnit(1),*) 'ERROR: Neither speed or sound (cs), nor'
          write(logUnit(1),*) '       time step length (dt) is defined'
          write(logUnit(1),*) "Solution: Provide speed of sound 'cs' " &
            &               //"in m/s"
          write(logUnit(1),*) "          or time step 'dt' in s"
          call tem_abort()
        end if
        if (btest(iError, aoterr_Fatal)) then
          write(logUnit(7),*) 'FATAL Error occured, while retrieving dt.'
          if (btest(iError, aoterr_WrongType)) then
            write(logUnit(1),*)'Variable has wrong type!'
            call tem_abort()
          endif
        end if
        cs_phy = me%dx*cs/me%dt
      else
        me%dt = me%dx * cs / cs_phy
      end if

      write(logUnit(1),*) '  cs = '//trim(tem_toStr(cs_phy))
      write(logUnit(1),*) '  dt = '//trim(tem_toStr(me%dt))

      ! define mole before loading density because if molecular weight is
      ! defined to compute reference mass than we need reference mole.
      ! try to load reference mole density, if not defined then set reference
      ! mole to inverse of avogadro's constant
      call aot_get_val( L = conf, thandle = thandle, key = 'moleDens0', &
        &               val = me%moleDens0, ErrCode = iError            )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(7),*) 'No value given for moleDens0.'
        if (btest(iError, aoterr_WrongType)) then
          write(logUnit(1),*)'Error! moleDens0 has wrong type!'
          call tem_abort()
        end if
      end if

      ! reference mole density not defined so try to load reference mole
      if (btest(iError, aoterr_NonExistent)) then
        write(logUnit(7),*)'WARNING: Reference moleDens0 is not found. '
        write(logUnit(7),*)'Loading reference mole0(mol):'
        call aot_get_val( L = conf, thandle = thandle, key = 'mole0', &
          &               val = me%mole0, ErrCode = iError            )
        if (btest(iError, aoterr_Fatal)) then
          write(logUnit(7),*)'FATAL Error occured, while retrieving mole0.'
          if (btest(iError, aoterr_WrongType)) then
            write(logUnit(1),*)'Variable has wrong type!'
            call tem_abort()
          end if
        end if
        ! reference mole is also not defined set inverse of Avogadro's constant
        if (btest(iError, aoterr_NonExistent)) then
          write(logUnit(7),*)'WARNING: Reference mole0 is not found. '
          write(logUnit(7),*)"Setting reference mole to inverse of Avogadro's" &
            &                //"constant."
          me%mole0 = mole_ref
        endif
        ! derive moleDens0 from mole0
        me%moleDens0 = me%mole0 / me%dx**3
      else ! if defined moleDens0 derive mole from moledensity
        me%mole0 = me%moleDens0 * me%dx**3
      endif

      ! try to load reference density, If not defined then try to load reference
      ! molecular weight, if not defined then try to load reference mass in 'kg'
      ! if that also is not defined prompt an error message
      ! load rho0
      call aot_get_val( L = conf, thandle = thandle, key = 'rho0', &
        &               val = me%rho0, ErrCode = iError            )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(3),*) 'No value given for rho0.'
        if (btest(iError, aoterr_WrongType)) then
          write(logUnit(1),*) 'rho0 has wrong type!'
          call tem_abort()
        end if
      end if

      ! reference mass density not defined so try to load reference
      ! molecular weight
      if (btest(iError, aoterr_NonExistent)) then
        write(logUnit(1),*)'WARNING: Reference mass density rho0 is not found.'
        write(logUnit(1),*)'Loading reference molecular weight ' &
          &              //'molWeight0(kg/mol):'
        call aot_get_val( L = conf, thandle = thandle, key = 'molWeight0', &
          &               val = me%molWeight0, ErrCode = iError            )
        if (btest(iError, aoterr_Fatal)) then
          write(logUnit(3),*)'No value given for molWeight0.'
          if (btest(iError, aoterr_WrongType)) then
            write(logUnit(1),*) 'molWeight has wrong type!'
            call tem_abort()
          end if
        end if
        ! if reference molecular weight is not defined than load reference mass
        if (btest(iError, aoterr_NonExistent)) then
          write(logUnit(5),*)'WARNING: Reference molecular weight molWeight0 ' &
          &                //'is not found. '
          write(logUnit(5),*)'Loading reference mass0(kg):'
          call aot_get_val( L = conf, thandle = thandle, key = 'mass0', &
            &               val = me%mass0, ErrCode = iError            )
          if (btest(iError, aoterr_Fatal)) then
            write(logUnit(3),*)'No value given for mass0'
            if (btest(iError, aoterr_WrongType)) then
              write(logUnit(1),*)'Variable has wrong type!'
              call tem_abort()
            end if
          end if
          ! if reference mass is also not defined prompt an error message
          if (btest(iError, aoterr_NonExistent)) then
            write(logUnit(1),*) 'ERROR: Unable to obtain reference mass0(kg)'
            write(logUnit(1),*) "Solution: Provide reference density 'rho0' " &
              &               //"in kg/m^3"
            write(logUnit(1),*) "or reference molecular weight 'molWeight0' " &
              &               //"in kg/mol"
            write(logUnit(1),*) "or reference mass 'mass0' in kg"
            call tem_abort()
          else
            ! reference mass is defined derive density and molWeight0
            me%rho0 = me%mass0/me%dx**3
            ! same as rho0/moleDens0
            ! molecular Weight in independent of dx
            me%molWeight0 = me%mass0/me%mole0
          endif
        else
        ! reference molecular weight is defined derive mass and density
          me%mass0 = me%molWeight0 * me%mole0
          me%rho0 = me%mass0/me%dx**3
        endif
      else
        ! reference density is defined derive mass and molWeigh
        me%mass0 = me%rho0 * me%dx**3
        ! same as rho0/moleDens0
        ! molecular Weight in independent of dx
        me%molWeight0 = me%mass0 / me%mole0
      endif

      ! try to load coloumb, if not defined set to fundamental electrical charge
      call aot_get_val( L = conf, thandle = thandle, key = 'coulomb0', &
        &               val = me%coulomb0, ErrCode = iError            )
      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(3),*)'WARNING: No value given for coulomb0'
        if (btest(iError, aoterr_WrongType)) then
          write(logUnit(1),*) 'coulomb0 has wrong type!'
          call tem_abort()
        end if
      end if
      if (btest(iError, aoterr_NonExistent)) then
        write(logUnit(4),*) 'WARNING: Reference coulomb0 is not found. '
        write(logUnit(4),*) "Setting reference coulomb0 to fundamental" &
          &                 //" electrical charge."
        me%coulomb0 = coulomb_ref
        !me%coulomb0 = me%rho0*me%dx**5 / me%dt**2
      endif

      ! load temp0
      call aot_get_val(L = conf, thandle = thandle, key = 'temp0',        &
        &              val = me%temp0, default = 1.0_rk, ErrCode = iError )
      if (btest(iError, aoterr_WrongType)) then
        write(logUnit(1),*)'temp0 has wrong type! Stopping!'
        call tem_abort()
      end if
    else
      ! No physics table defined.
      me%dt         = 1._rk
      me%dx         = 1._rk
      me%rho0       = 1._rk
      me%moleDens0  = 1._rk
      me%molWeight0 = 1._rk
      me%mass0      = 1._rk
      me%mole0      = 1._rk
      me%temp0      = 1._rk
      me%coulomb0   = 1._rk
      write(logUnit(1),*) ' Physics table is not defined.'
      write(logUnit(5),*) ' All values are expressed in Lattice Units.'
    end if

    call aot_table_close( L=conf, thandle=thandle )

    ! Assign dx and dt for each level
    allocate(me%dxLvl( tree%global%minLevel:tree%global%maxLevel ))
    allocate(me%dtLvl( tree%global%minLevel:tree%global%maxLevel ))
    allocate(  me%fac( tree%global%minLevel:tree%global%maxLevel ))

    ! dx and dt is set according scaling type
    me%dxLvl( tree%global%minLevel:tree%global%maxLevel ) =  &
      & set_values_by_levels( me%dx, tree%global%minLevel,   &
      &                              tree%global%maxLevel, 2 )

    me%dtLvl( tree%global%minLevel:tree%global%maxLevel ) =            &
      & set_values_by_levels( me%dt, tree%global%minLevel,             &
      &                              tree%global%maxLevel, scaleFactor )

    ! compute and store conversion factors in converstionFac type
    call mus_set_convFac( me = me, minLevel = tree%global%minLevel, &
      &                            maxLevel = tree%global%maxLevel  )
    ! set scale factor
    call mus_set_scaleFac( me, tree%global%minLevel, tree%global%maxLevel )

    call mus_physics_dump2outUnit( me, logUnit(4), tree%global%minLevel, &
      &                           tree%global%maxLevel                   )

    call tem_horizontalSpacer(fUnit = logUnit(1))

  end subroutine mus_load_physics
! **************************************************************************** !


! **************************************************************************** !
  !> This routine computed conversion factors for lattice to physical units.
  !! inverse of this factors can be used to convert from physical to lattice
  !! units.\n
  !! use reference density to parmeterize kg and reference mole density
  !! to parmeterize mol.\n
  !! Multiply these factors with the LB quantity to get the physical quantity
  !! Divide the physical quantity by these factors to get the LB units.
  subroutine mus_set_convFac( me, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    type( mus_physics_type ), intent(inout) :: me !< physics type
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------

    do iLevel = minLevel, maxLevel
      ! length
      me%fac( iLevel )%length = me%dxLvl( iLevel )
      ! time
      me%fac( iLevel )%time = me%dtLvl( iLevel )
      ! velocity
      me%fac( iLevel )%vel = me%dxLvl( iLevel )/me%dtLvl( iLevel )
      ! kinematic viscosity
      me%fac( iLevel )%visc = me%dxLvl( iLevel )**2/me%dtLvl( iLevel )
      ! dynamic viscosity
      me%fac( iLevel )%viscDyna = me%rho0 * me%fac( iLevel )%visc
      ! acceleration
      me%fac( iLevel )%accel = me%dxLvl( iLevel )/me%dtLvl( iLevel )**2
      ! force
      me%fac( iLevel )%force = me%rho0*me%dxLvl( iLevel )**4                   &
        &                    / me%dtLvl( iLevel )**2
      ! body_force
      me%fac( iLevel )%body_force = me%rho0*me%dxLvl( iLevel )                 &
        &                         / me%dtLvl( iLevel )**2
      ! pressure or shear stress
      me%fac( iLevel )%press = me%rho0 * me%dxLvl( iLevel )**2                 &
        &                    / me%dtLvl( iLevel )**2
      ! strain rate
      me%fac( iLevel )%strainRate = 1._rk / me%dtLvl( iLevel )
      ! Energy
      me%fac( iLevel )%energy = me%rho0*me%dxLvl( iLevel )**5                  &
        &                     / me%dtLvl( iLevel )**2
      ! Mass
      !me%fac( iLevel )%mass = me%rho0*me%dxLvl( iLevel )**3
      ! Molar mass or Molecular weight
      !me%fac( iLevel )%molWeigh = me%fac( iLevel )%mass/me%mole0
      ! mole density
      !me%fac( iLevel )%moleDens = me%mole0/me%dxLvl( iLevel )**3
      ! mole flux (mol/m^2/s) (moleDens*velocity)
      me%fac( iLevel )%moleFlux = me%moleDens0*me%dxLvl( iLevel ) &
        &                       / me%dtLvl( iLevel )
      ! Charge density (C/m^3) = Coulomb0/dx^3
      me%fac( iLevel )%chargeDens = me%coulomb0 / me%dxLvl( iLevel )**3
      ! Current density (C/s/m^2) = Coulomb0/dt/dx^2
      me%fac( iLevel )%currentDens = me%coulomb0 / me%dxLvl( iLevel )**2 &
        &                          / me%dtLvl( iLevel )
      ! mass flux (kg/m^2/s)  (density*velocity)
      me%fac( iLevel )%flux = me%rho0*me%dxLvl( iLevel )/me%dtLvl( iLevel )
      !diffusivity
      me%fac( iLevel )%diffusivity = me%dxLvl( iLevel )**2/me%dtLvl( iLevel )
      !faraday
      me%fac( iLevel )%faraday = me%coulomb0/me%mole0
      !gas constant
      me%fac( iLevel )%gasConst = me%rho0*me%dxLvl( iLevel )**5     &
        &                       / me%dtLvl( iLevel )**2 / me%mole0  &
        &                       / me%temp0
      ! Potential
      me%fac( iLevel )%potential = me%rho0*me%dxLvl( iLevel )**5  &
        &                     / me%dtLvl( iLevel )**2 / me%coulomb0
      ! SourceCoeff
      me%fac( iLevel )%sourceCoeff = 1.0_rk / me%dtLvl( iLevel )

    end do

  end subroutine mus_set_convFac
! **************************************************************************** !


! **************************************************************************** !
  !> This routine creates musubi specific lua function to compute dx and dt.
  !!
  subroutine mus_create_funcStr( fun_str )
    ! --------------------------------------------------------------------------
    !> This string contains lua functions to compute dt from visocosity or
    !! velocity
    character(len=*) :: fun_str
    ! --------------------------------------------------------------------------
    character(len=256) :: dtVel, dtOmega, dtVisc, dxLevel, dxnL, omega_dt,     &
      &                   vel_dt
    ! --------------------------------------------------------------------------
    write(dxLevel,*) 'function getdxFromLevel(t) '//                           &
      &              ' return t.len_bnd/2^t.level end'

    write(dxnL,*) 'function getdxFromnL(t)'//                                  &
      &           ' return t.len_p/t.nL end'

    write(dtVel,*) 'function getdtFromVel(t) return t.dx*t.u_l/t.u_p end'

    write(dtOmega,*) 'function getdtFromOmega(t) '//                           &
      &              ' return ((1.0/t.omega - 0.5)/3.0)*t.dx*t.dx/t.nu_p end'

    write(dtVisc,*) ' function getdtFromVisc(t)'//                             &
      &             ' return t.dx*t.dx*t.nu_l/t.nu_p end'

    write(omega_dt,*) 'function getOmegaFromdt(t) '//                          &
      &               ' return 1.0/(3.0*t.nu_p*t.dt/(t.dx*t.dx) + 0.5) end'

    write(vel_dt, *) 'function getVelFromdt(t) '//                             &
      &              ' return t.u_p*t.dt/t.dx end'

    write(fun_str,*) trim(dtVel)//' '//trim(dtOmega)//' '//trim(dtVisc)        &
      & //' '//trim(dxLevel)//' '//trim(dxnL)//' '//trim(omega_dt)             &
      & //' '//trim(vel_dt)

  end subroutine mus_create_funcStr
! **************************************************************************** !


! **************************************************************************** !
  !> This routine write reference physics parameters into solver specific
  !! string in lua format.
  !!
  !! This dumped table is loaded back using mus_load_physics
  subroutine mus_physics_out( me, conf )
    ! --------------------------------------------------------------------------
    type( mus_physics_type ), intent(in) :: me !< physics type
    type( aot_out_type ) :: conf
    ! --------------------------------------------------------------------------

    call aot_out_open_table( put_conf = conf, tname = 'physics' )

    call aot_out_val( put_conf = conf, vname = 'dt', val = me%dt )
    call aot_out_val( put_conf = conf, vname = 'rho0', val = me%rho0 )
    call aot_out_val( put_conf = conf, vname = 'moleDens0', val = me%moleDens0 )
    call aot_out_val( put_conf = conf, vname = 'molWeight0', &
      &               val = me%molWeight0 )
    call aot_out_val( put_conf = conf, vname = 'temp0', val = me%temp0 )
    call aot_out_val( put_conf = conf, vname = 'coulomb0', val = me%coulomb0 )

    call aot_out_close_table( put_conf = conf )

  end subroutine mus_physics_out
! **************************************************************************** !


! **************************************************************************** !
  !> This routine write physics convert factor into solver specific string in
  !! lua format.
  !! use reference density to parmeterize kg and reference mole density
  !! to parmeterize mol.
  subroutine mus_physics_out_conv( me, conf, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    type( mus_physics_type ), intent(in) :: me !< physics type
    type( aot_out_type ) :: conf
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------

    call aot_out_open_table( put_conf = conf, tname = 'physics' )
    do iLevel = minLevel, maxLevel
      call aot_out_open_table( put_conf = conf )
      call aot_out_val( put_conf = conf, vname = 'length',                     &
        &               val = me%fac(iLevel)%length )
      call aot_out_val( put_conf = conf, vname = 'time',                       &
        &               val = me%fac(iLevel)%time )
      call aot_out_val( put_conf = conf, vname = 'density',                    &
        &               val = me%rho0 )
      call aot_out_val( put_conf = conf, vname = 'vel',                        &
        &               val = me%fac(iLevel)%vel )
      call aot_out_val( put_conf = conf, vname = 'visc',                       &
        &               val = me%fac(iLevel)%visc )
      call aot_out_val( put_conf = conf, vname = 'viscDyna',                   &
        &               val = me%fac(iLevel)%viscDyna )
      call aot_out_val( put_conf = conf, vname = 'accel',                      &
        &               val = me%fac(iLevel)%accel )
      call aot_out_val( put_conf = conf, vname = 'force',                      &
        &               val = me%fac(iLevel)%force )
      call aot_out_val( put_conf = conf, vname = 'body_force',                 &
        &               val = me%fac(iLevel)%body_force )
      call aot_out_val( put_conf = conf, vname = 'press',                      &
        &               val = me%fac(iLevel)%press )
      call aot_out_val( put_conf = conf, vname = 'strainRate',                 &
        &               val = me%fac(iLevel)%strainRate )
      call aot_out_val( put_conf = conf, vname = 'energy',                     &
        &               val = me%fac(iLevel)%energy )
      call aot_out_val( put_conf = conf, vname = 'temp',                       &
        &               val = me%temp0 )
      call aot_out_val( put_conf = conf, vname = 'moleDens',                   &
        &               val = me%moleDens0 )
      call aot_out_val( put_conf = conf, vname = 'molWeight',                  &
        &               val = me%molWeight0 )
      call aot_out_val( put_conf = conf, vname = 'coulomb',                    &
        &               val = me%coulomb0 )
      call aot_out_val( put_conf = conf, vname = 'moleFlux',                   &
        &               val = me%fac(iLevel)%moleFlux )
      call aot_out_val( put_conf = conf, vname = 'chargeDens',                 &
        &               val = me%fac(iLevel)%chargeDens )
      call aot_out_val( put_conf = conf, vname = 'currentDens',                &
        &               val = me%fac(iLevel)%currentDens )
      call aot_out_val( put_conf = conf, vname = 'flux',                       &
        &               val = me%fac(iLevel)%flux )
      call aot_out_val( put_conf = conf, vname = 'diffusivity',                &
        &               val = me%fac(iLevel)%diffusivity )
      call aot_out_val( put_conf = conf, vname = 'faraday',                    &
        &               val = me%fac(iLevel)%faraday )
      call aot_out_val( put_conf = conf, vname = 'gasConst',                   &
        &               val = me%fac(iLevel)%gasConst )
      call aot_out_close_table( put_conf = conf )
    end do

    call aot_out_close_table( put_conf = conf )

  end subroutine mus_physics_out_conv
! **************************************************************************** !

  ! ************************************************************************** !
  pure function set_values_by_levels( valMinLevel, minLevel, maxLevel, &
    &                                 scaleFac )                       &
    &           result( values )
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: valMinLevel !< value at min level
    integer,       intent(in) :: minLevel
    integer,       intent(in) :: maxLevel
    integer,       intent(in) :: scaleFac !< scale factor between levels
    real(kind=rk) :: values(minLevel:maxLevel) !< return value
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------

    values( minLevel:maxLevel ) = valMinLevel / &
      & real([(scaleFac**(iLevel-minLevel), iLevel=minLevel, maxLevel)],rk)

  end function set_values_by_levels
  ! ************************************************************************** !

  ! ************************************************************************** !
  subroutine mus_physics_dump2outUnit( me, outUnit, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    type( mus_physics_type ), intent(in) :: me
    integer, intent(in) :: outUnit
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    ! --------------------------------------------------------------------------
    integer :: iLevel

    write(outUnit,'(A)')       'Reference physical quantities:'
    write(outUnit,"(A)") '  rho0       = '//trim(tem_toStr(me%rho0))
    write(outUnit,"(A)") '  mass0      = '//trim(tem_toStr(me%mass0))
    write(outUnit,"(A)") '  mole0      = '//trim(tem_toStr(me%mole0))
    write(outUnit,"(A)") '  moleDens0  = '//trim(tem_toStr(me%moleDens0))
    write(outUnit,"(A)") '  molWeight0 = '//trim(tem_toStr(me%molWeight0))
    write(outUnit,"(A)") '  Temp0      = '//trim(tem_toStr(me%temp0))
    write(outUnit,"(A)") '  coulomb0   = '//trim(tem_toStr(me%coulomb0))

    do iLevel = minLevel, maxLevel
      write(outUnit,"(A,I0)")    'Conversion factors on level: ', iLevel
      write(outUnit,"(A)") '  dx            = ' &
        &                  // trim(tem_toStr(me%dxLvl( iLevel )))
      write(outUnit,"(A)") '  dt            = ' &
        &                  // trim(tem_toStr(me%dtLvl( iLevel )))
      write(outUnit,"(A)") '  pressure      = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%press))
      write(outUnit,"(A)") '  velocity      = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%vel))
      write(outUnit,"(A)") '  kine. visc.   = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%visc))
      write(outUnit,"(A)") '  dyn. visc.    = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%viscDyna))
      write(outUnit,"(A)") '  acceleration  = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%accel))
      write(outUnit,"(A)") '  force         = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%force))
      write(outUnit,"(A)") '  body force    = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%body_force))
      write(outUnit,"(A)") '  strain rate   = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%strainRate))
      write(outUnit,"(A)") '  energy        = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%energy))
      write(outUnit,"(A)") '  moleFlux      = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%moleFlux))
      write(outUnit,"(A)") '  chargeDens    = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%chargeDens))
      write(outUnit,"(A)") '  currentDens   = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%currentDens))
      write(outUnit,"(A)") '  massflux      = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%flux))
      write(outUnit,"(A)") '  diffusivity   = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%diffusivity))
      write(outUnit,"(A)") '  faraday       = ' &
        &                  // trim(tem_toStr(me%fac( iLevel )%faraday))
      write(outUnit,"(A)") '  gasConst      = ' &
        &                  //trim(tem_toStr(me%fac( iLevel )%gasConst))
      write(outUnit,"(A)") '  potential     = ' &
        &                  //trim(tem_toStr(me%fac( iLevel )%potential))
      write(outUnit,"(A)") ''
    end do

    write(outUnit,"(A)")   'Over level scale factor: '
    do iLevel = minLevel, maxLevel
      write(outUnit,"(A)") '    pressure: '                                    &
        &                  // trim( tem_toStr(me%pFac(iLevel,                  &
        &                                             minLevel:maxLevel), ',') )
      write(outUnit,"(A)") '    velocity: '                                    &
        &                  // trim( tem_toStr(me%vFac(iLevel,                  &
        &                                             minLevel:maxLevel), ',') )
      write(outUnit,"(A)") ' strain rate: '                                    &
        &                  // trim( tem_toStr(me%sFac(iLevel,                  &
        &                                             minLevel:maxLevel), ',') )
    end do

  end subroutine mus_physics_dump2outUnit
  ! ************************************************************************** !

  ! ************************************************************************** !
  subroutine mus_set_scaleFac( me, minLevel, maxLevel )
    ! --------------------------------------------------------------------------
    type( mus_physics_type ), intent(inout) :: me
    integer, intent(in) :: minLevel
    integer, intent(in) :: maxLevel
    ! --------------------------------------------------------------------------
    integer :: iLevel
    ! --------------------------------------------------------------------------

    allocate( me%pFac( minLevel:maxLevel, minLevel:maxLevel ) )
    allocate( me%vFac( minLevel:maxLevel, minLevel:maxLevel ) )
    allocate( me%sFac( minLevel:maxLevel, minLevel:maxLevel ) )

    do iLevel = minLevel, maxLevel
      me%pFac( :, iLevel ) = me%fac(:)%press      / me%fac(iLevel)%press
      me%vFac( :, iLevel ) = me%fac(:)%vel        / me%fac(iLevel)%vel
      me%sFac( :, iLevel ) = me%fac(:)%strainRate / me%fac(iLevel)%strainRate
    end do

  end subroutine mus_set_scaleFac
  ! ************************************************************************** !

end module mus_physics_module
! **************************************************************************** !
