require 'seeder'
require 'musubi'

print('-----------------------------------------------------------------------')
print('Simulation name: ', simulation_name)
print('-----------------------------------------------------------------------')
print('---Mesh parameters---')
print('length       =', l_ch, '[m]')
print('height       =', h_ch, '[m]')
print('l_nozzle     =', l_nozzle, '[m]')
print('IDtoOD       =', inner_to_outer_ratio)
print('length_bnd   =', length_bnd, '[m]')
print('nLength      =', nLength)
print('nHeight      =', nElems_h_ch)
print('level        =', level)
print('nozzle_level =', nozzleLevel)
print('-----------------------------------------------------------------------')
print('---Flow parameters---')
print('---In physical units---')
print('Re                =', Re)
print('Vel. inflow       =', vel_phy, '[m/s]')
print('Kinematic visc.   =', nu_phy, '[m^s/2]')
print('Density           =', rho0_phy, '[kg/m^3]')
print('Press. ambient    =', press_ambient, '[N/m^2]')
print('Element size (dx) =', dx, '[m]')
print('Time step (dt)    =', dt, '[s]')
print('--- In lattice units---')
print('Vel.              =', vel_lat)
print('Ma                =', Ma_lat)
print('Exp. Ma at neck   =', vel_lat/cs_lat/inner_to_outer_ratio)
-- Lattice viscosity
nu_lat = nu_phy * dt / dx^2
print('Kinematic visc.   =', nu_lat)
-- Lattice relaxation parameter
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5 )
print('Relaxation param:')
for iL = minlevel, nozzleLevel do
  print('level:',iL,' omLvl: ', 1.0 / ( 0.5 + (1.0/omega-0.5)*2^(iL-minlevel) ))
end 

