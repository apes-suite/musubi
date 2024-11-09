require 'seeder'
require 'musubi'

print('-----------------------------------------------------------------------')
print('Simulation name: ', simulation_name)
if useRecheck then
  print('-----------------------------------------------------------------------')
  print("Attention! You are using a truncated version of the test case. To\n"
        .."reproduce the results presented in the documentation\n"
        .."please use the following option:\n"
        .."useRecheck = false")
  print('-----------------------------------------------------------------------')
end
print('-----------------------------------------------------------------------')
print('------Mesh parameters------')
print('height       =', height     )
print('width        =', width      )
print('length       =', length     )
print('in length    =', nLength    )
print('length_bnd   =', length_bnd )
print('level        =', level      )
print('-----Number of elements----')
print('in height    =', nHeight    )
print('in width     =', nHeight    )
print('in length    =', nLength    )
print('--------Resolution---------')
print('spatial    =', dx           )
print('temporal   =', dt           )
print('-----------------------------------------------------------------------')
print('------Flow parameters------'                  )
print('------In physical units----'                  )
print('Re                =', Re                      )
print('Vel. mean         =', vel_mean_phy, '[m/s]'   )
print('Vel. max.         =', vel_max_phy,  '[m/s]'   )
print('Kinematic visc.   =', nu_phy,       '[m^s/2]' )
print('Density           =', rho0_phy,     '[kg/m^3]')
print('Press. ambient    =', press_ambient,'[N/m^2]' )
print('Speed of sound    =', cs_phy,  '[m/s]'   )
print('Element size (dx) =', dx,   '[m]'     )
print('Time step (dt)    =', dt,   '[s]'     )
print('------In lattice units-----' )
print('Vel.              =', vel_lat)
print('Ma                =', Ma)
-- Lattice viscosity
nu_lat = nu_phy * dt / dx^2
print('Kinematic lattice visc. =', nu_lat)
-- Lattice relaxation parameter
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5  )
print('Relaxation param. =', omega)
print('-----------------------------------------------------------------------')
