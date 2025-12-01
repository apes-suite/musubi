resolution = 2 ^ 1
umax = 0.7 / 4
omega = 1.0 / 0.6
simulation_name = 'brinkman'
tmax               = 2000 * resolution^2
interval           = 200 * resolution^2

--Dimensions and refinement level
length      = 20.0               --cm
domainlen   = 16.0
dx          = 0.25 / resolution
dx_eps = dx/2^20
nLength     = length/dx 
level       = math.ceil(math.log(nLength)/math.log(2))
length_bnd  = (2^level)*dx
Dia         = 2                  --cm
Dia_L       = Dia/dx
radius      = Dia/2.0            --cm

--Physical parameters
rho_phy     = 1000e-6     --Kg/cm^3
mu_phy      = 2.8e-5      --Kg/cm/s (0.035 Poise)
nu_phy      = mu_phy/rho_phy
F0          = 10
Re          = (Dia*umax)/nu_phy

--Lattice parameters with diffusive scaling
nu_L        = (1.0/omega-0.5)/3.0
dt          = nu_L*dx^2/nu_phy
u_in_L    = umax*dt/dx

--Lattice pressure is cs^2*rho_L = 1./3.
press_ref   = rho_phy*(dx^2)/(3.*dt^2)
press_phy   = 0.+press_ref
Re_L        = (Dia_L*u_in_L)/nu_L

--BC
bc_origin = { -dx, -dx, -dx/2.0-2*dx } 
seed_orig = { length/2, Dia/2, 0.0 }

--Brinkman term
Fc = F0 * nu_phy
k1 = rho_phy * umax / (math.exp(math.sqrt(F0) * Dia) - math.exp(-math.sqrt(F0) * Dia))
k2 = -k1