require 'seeder'
require 'musubi'

print('-----------------------------------------------------------------------')
print('Simulation name: ', simulation_name)
print('-----------------------------------------------------------------------')
print('---Mesh parameters---')
print('length       =', length)
print('height       =', height)
print('nLength      =', nLength)
print('nHeight      =', nHeight)
print('length_bnd   =', length_bnd)
print('level        =', level)
print('diameter     =', dia_cyl)
print('X pos cyl    =', pos_cyl[1])
print('Y pos cyl    =', pos_cyl[2])
print('Z pos cyl    =', pos_cyl[3])
print('-----------------------------------------------------------------------')
print('---Flow parameters---')
print('---In physical units---')
print('Re                =', Re)
print('Ma                =', Ma)
print('Vel. mean         =', vel_mean_phy, '[m/s]')
print('Vel. max.         =', vel_max_phy, '[m/s]')
print('Kinematic visc.   =', nu_phy, '[m^s/2]')
print('Speed of sound    =', cs_phy, '[m^s/2]')
print('Density           =', rho0_phy, '[kg/m^3]')
print('Press. ambient    =', press_ambient, '[N/m^2]')
print('Element size (dx) =', dx, '[m]')
print('Time step (dt)    =', dt, '[s]')
print('--- In lattice units---')
print('Vel.              =', vel_lat)
-- Lattice viscosity
nu_lat = nu_phy * dt / dx^2
print('Kinematic visc.   =', nu_lat)
-- Lattice relaxation parameter
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5 )
print('Relaxation param. =', omega)
print('-----------------------------------------------------------------------')
