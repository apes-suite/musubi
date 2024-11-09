require 'seeder'
require 'musubi'

print('-----------------------------------------------------------------------')
print('Simulation name: ', simulation_name)
print('-----------------------------------------------------------------------')
print('---Mesh parameters---')
print('Inner dia.   =', dia_inner)
print('Outer dia.   =', dia_outer)
print('Inner nDia.  =', nDia_inner)
print('Outer nDia.  =', nDia_outer)
print('length_bnd   =', length_bnd)
print('level        =', level)
print('-----------------------------------------------------------------------')
print('---Flow parameters---')
print('---In physical units---')
print('Re                =', Re)
print('Vel.              =', vel_phy, '[m/s]')
print('Angular Vel.      =', angular_vel, '[rad/s]')
print('Kinematic visc.   =', nu_phy, '[m^s/2]')
print('Density           =', rho0_phy, '[kg/m^3]')
print('Press. ambient    =', press_ambient, '[N/m^2]')
print('Element size (dx) =', dx, '[m]')
print('Time step (dt)    =', dt, '[s]')
print('--- In lattice units---')
print('Vel.              =', vel_lat)
print('Ma                =', Ma)
-- Lattice viscosity
nu_lat = nu_phy * dt / dx^2
print('Kinematic visc.   =', nu_lat)
-- Lattice relaxation parameter
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5 )
print('Relaxation param. =', omega)
print('-----------------------------------------------------------------------')
