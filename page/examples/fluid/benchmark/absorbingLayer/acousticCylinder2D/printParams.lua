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
print('-----------------------------------------------------------------------')
print('---Flow parameters---')
print('---In physical units---')
print('Ma                =', Ma)
print('Re                =', Re)
print('Ma                =', vel_mean_phy/cs_phy)
print('Vel. mean         =', vel_mean_phy, '[m/s]')
print('Kinematic visc.   =', nu_phy, '[m^s/2]')
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
print('Relaxation param:')
for iL = minlevel, cylLevel do
  print('level:',iL,' omLvl: ', 1.0 / ( 0.5 + (1.0/omega-0.5)*2^(iL-minlevel) ))
end 
print('-----------------------------------------------------------------------')
