require 'seeder'
require 'musubi'

print('-----------------------------------------------------------------------')
print('Simulation name: ', simulation_name)
print('-----------------------------------------------------------------------')
print('---Mesh parameters---')
print('length       =', length)
print('height       =', length)
print('nLength      =', nLength)
print('length_bnd   =', length_bnd)
print('level        =', level)
print('-----------------------------------------------------------------------')
print('---Flow parameters---')
print('---In physical units---')
print('Re                =', Re)
print('Vel.              =', vel_phy, '[m/s]')
print('Kinematic visc.   =', nu_phy, '[m^s/2]')
print('Density           =', rho0_phy, '[kg/m^3]')
print('Press. ambient    =', press_ambient, '[N/m^2]')
print('Element size (dx) =', dx, '[m]')
print('Time step (dt)    =', dt, '[s]')
print('Press. drop       =', press_drop, '[N/m^2]')
print('Press. west       =', press_ambient+press_drop, '[N/m^2]')
print('Press. east       =', press_ambient, '[N/m^2]')
print('Press. gradient   =', press_grad, '[N/m^3]')
print('--- In lattice units---')
print('Vel.              =', vel_lat)
print('Ma                =', Ma_lat)
-- Lattice viscosity
nu_lat = nu_phy * dt / dx^2
print('Kinematic visc.   =', nu_lat)
-- Lattice relaxation parameter
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5 )
print('Relaxation param. =', omega)
print('-----------------------------------------------------------------------')
