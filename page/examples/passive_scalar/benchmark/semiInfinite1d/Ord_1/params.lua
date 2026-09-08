-- simulation params
t_ini = 0.
t_total = 0.5

-- physical params
D = 5e-2
C_0 = 1
dx = 0.016
L = 2

-- lattice params
-- source coefficient in lattice unit
lambda_L = 0.1
tau = 1.2
D_L = (tau-0.5)/3.0
v_L = 0.1 

-- derived params
dt = D_L*dx^2/D
v = v_L*dx/dt
-- source coefficient in physical unit
lambda = lambda_L / dt
level_usr = math.ceil(math.log(L/dx) / math.log(2))
length_bnd = dx * 2^level_usr
-- maximum time step
mStep = math.ceil((t_total - t_ini) / dt)

