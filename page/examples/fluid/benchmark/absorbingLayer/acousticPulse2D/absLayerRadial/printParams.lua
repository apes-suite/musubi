require "seeder"
require "musubi"
print('-----------------------------------------------------------------------')
print('Simulation name: ', simulation_name )
print('-----------------------------------------------------------------------')
print('-------Mesh parameters [m]-----')
print('Length of domain (x) =', length )
print('Height of domain (y) =', height )
print('-------Number of elements------')
print('in height    =', nLength        )
print('in length    =', nLength        )
print('----------Resolution-----------')
print('dx (spatial)  =', dx            )
print('dt (temporal) =', dt            )
print('Level for dx  =', level         )

print('-----------Absorbing layer--------------')
print('Thickness of abs.layer =', abs_thickness )
print('Damping factor (min)   =', 4.0/2-0.001   )
print('Damping factor         =', damp_factor   )
print('Damping factor (max)   =', 4.0/0.5-0.001 )

print('------------Acoustic source--------------')
print('Background density   =', background       )
print('Background pressure  =', press_ambient    )
print('Amplitude of density =', amplitude        )

print('-------------Timings---------------------')
print('Max. simulation time         =', tmax     )

print('-----------------Flow parameters----------------')
print('Kinematic visc.  =', nu_phy,   '[m^s/2]'         )
print('Ma               =', Ma                          )
print('---------------In physical units----------------')
print('Mean velocity      = ', vel_mean_phy,  '[m/s]'   )
print('Density            = ', rho0_phy,   '[kg/m^3]'   )
print('------------In lattice units--------------')
-- Lattice viscosity
nu_lat = nu_phy * dt / dx^2
print('Kinematic lattice visc. =', nu_lat         )
print('Lattice velocity        =', vel_lat        )
-- Lattice relaxation parameter
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5           )
print('Relaxation param.       =', omega          )
print('-----------------------------------------------------------------------')


--print("Parameter File ")
--print(" Shear visc.:   "..nu_phy)
--print(" cs Phys ref:   "..csPhys)
--print("       Level:   "..level)
--print("       refdx:   "..dx)
--print("          dt:   "..dt)
--print("          T0:   "..T0)
--print("          Tp:   "..Tp)
--print("                ")
--print("   cs LB ref:   "..cs_lat)
--nuLB = nu_phy * dt / dx^2
--print("   nu LB ref:   "..nuLB)
--omega = 1.0 / ( nuLB/cs_lat^2.0 + 0.5 )
--print(" Resulting om   "..omega)


