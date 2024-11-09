require 'seeder'
require 'musubi'

print('-----------------------------------------------------------------------')
print('Simulation name: ', simulation_name)
print('-----------------------------------------------------------------------')
print('---Mesh parameters---')
print('length       =', length)
print('height       =', height)
print('width        =', width)
print('diameter     =', diameter)
print('nLength      =', nLength)
print('nHeight      =', nHeight)
print('nWidth       =', nWidth)
print('nDiameter    =', nDiameter*2^(sphere_level-minlevel))
print('length_bnd   =', length_bnd)
print('level        =', level)
print('-----------------------------------------------------------------------')
print('---Flow parameters---')
print('---In physical units---')
print('Re                =', Re)
print('Vel. mean         =', vel_phy, '[m/s]')
print('Kinematic visc.   =', nu_phy, '[m^s/2]')
print('Density           =', rho0_phy, '[kg/m^3]')
print('Press. ambient    =', press_ambient, '[N/m^2]')
print('Bnd layer thick.  =', 1.3*diameter/math.sqrt(Re)) 
print('Element size (dx) =', dx, '[m]')
print('Time step (dt)    =', dt, '[s]')
print('Finest elem. size =', dx/2^(sphere_level-minlevel))
print('Ramping time      =', t_ramp, '[s]')
print('Sim. end time     =', tmax_phy, '[s]')
print('--- In lattice units---')
print('Vel.              =', vel_lat)
print('Ma                =', Ma_lat)
-- Lattice viscosity
nu_lat = nu_phy * dt / dx^2
print('Kinematic visc.   =', nu_lat)
-- Lattice relaxation parameter
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5 )
print('Relaxation param:')
for iL = minlevel, sphere_level do
  print('level:',iL,' omLvl: ', 1.0 / ( 0.5 + (1.0/omega-0.5)*2^(iL-minlevel) ))
end 
print('-----------------------------------------------------------------------')
