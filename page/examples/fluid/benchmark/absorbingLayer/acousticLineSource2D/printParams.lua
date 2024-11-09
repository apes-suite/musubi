require "seeder"
require "musubi"
print("Parameter File ")
print(" Shear visc.:   "..nuPhys)
print(" cs Phys ref:   "..csPhys)
print("       Level:   "..level)
print("       refdx:   "..dx)
print("          dt:   "..dt)
print("          T0:   "..T0)
print("          Tp:   "..Tp)
print("                ")
print("   cs LB ref:   "..csLB)
nuLB = nuPhys * dt / dx^2
print("   nu LB ref:   "..nuLB)
omega = 1.0 / ( nuLB/csLB^2.0 + 0.5 )
print(" Resulting om   "..omega)


