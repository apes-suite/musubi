----------------------- PLEASE READ THIS ---------------------------!!!
-- This input file is set up to run for regression check
-- Please make sure you DO NOT MODIFY AND PUSH it to the repository
-- 
-- This test case is based on the paper:
-- Peng, C., Geneva, N., Guo, Z., & Wang, L. P. (2018). 
-- Direct numerical simulation of turbulent pipe flow using the lattice 
-- Boltzmann method. Journal of Computational Physics, 357, 16–42. 
--------------------------------------------------------------------!!!
-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

run_initial_period = true

--! [Local variables]
-- Flow parameters
-- Mach number
Ma = 0.1
-- Frictional Reynolds number 
Re_tau = 180
-- bulk Reynolds number according to the reference paper (see top)
Re_bulk = 5300
-- mean bulk velocity
vel_bulk_phy = 1.0
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = vel_bulk_phy * diameter / Re_bulk
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Friction velocity [m/s]
vel_fric_phy = Re_tau * nu_phy / radius
-- Characteristics length scale or wall unit [m]
y_star = nu_phy / vel_fric_phy
-- Characteristics time time or large eddy-turnover time [s]
T_star = radius / vel_fric_phy
-- Acceleration for driving force [m/s^2]
g_acce = 2 * vel_fric_phy * vel_fric_phy / radius
-- Normalized distance from wall to center
y_plus_radius = radius / y_star
-- Speed of sound
cs_phy = vel_bulk_phy / Ma
-- Ambient pressure
press_ambient = rho0_phy * cs_phy^2

------------ Compute physical time step from lattice Mach number ---------------
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Physical time step computed from physical and lattice velocity
dt = dx * cs_lat / cs_phy
-- Lattice viscosity
nu_lat = nu_phy*dt/dx^2
-- Lattice relaxation parameter
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5 )
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
if shepherd then
  tmax_phy = T_star/8
elseif run_initial_period then
  tmax_phy = 3*T_star
else
  tmax_phy = 100*T_star
end
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = T_star/8
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = tmax_phy
-- Termination wall clock time [s]
wall_clock = 02*60*60-5*60
-- Sampling average iteration
samp_avg_iter = math.min(math.ceil(50*T_star/dt-5), tmax_iter-5)
------------------------- End of time settings ---------------------------------

---------------------------Output directories-----------------------------------
-- tracking folder
tracking_fol = './tracking/'
-- restart folder
restart_fol = 'restart/'
---------------------End of output directories----------------------------------
--! [Local variables]

---------------------------- Lua functions -------------------------------------
-- Initial velocities
function vel_x_initial(x, y, z)
  r = math.sqrt( y*y + z*z)
  y_plus = (radius - r) / y_star
  if y_plus <= 10.8 then
    return vel_fric_phy * y_plus
  else
    return vel_fric_phy * ( math.log(y_plus)/0.4 + 5.0 )
  end
end
-- Driving force till 3-eddy turn over times
-- Weighting parameter that distributes the perturbation in radial and azimuthal direction
kappa = 0.5 
k_x = 3.0 -- Wavenumber of the perturbation force in streamwise direction
k_t = 2.0 -- Wavenumber of the perturbation force in azimuthal direction
T_f = 2000 -- Forcing period
B_0 = 50.0 -- Forcing magnitude
pi = math.pi
-- Region to add forcing perturbation
min_radius = 0.2*radius
max_radius = 0.4*radius
function force_till3T_star(x, y, z, t)
  r = radius-math.sqrt(y*y + z*z)
  theta = math.tan(y/z)
  if r >= min_radius and r <= min_radius+max_radius then
    term_1 = 2*pi*t/T_f
    term_2 = 2*pi*(radius - r - min_radius)/max_radius
    term_3 = k_x*2*pi*x/length
    term_4 = k_t*theta
    f_x_pert = - g_acce * B_0 * radius/r * math.sin(term_1)
             * math.sin(term_2) * math.sin(term_3) * math.cos(term_4)
    f_r_pert = - g_acce * kappa * B_0 * radius/r * k_x * max_radius/length
             * math.sin(term_1) * (1.0-math.cos(term_2)) 
             * math.cos(term_3) * math.cos(term_4)
    f_t_pert = g_acce * (1.0 - kappa) * B_0 * k_x/k_t * 2 * pi * radius/length
             * math.sin(term_1) * math.sin(term_2) 
             * math.cos(term_3) * math.sin(term_4)
    f_y_pert = f_r_pert * math.cos(f_t_pert)
    f_z_pert = f_r_pert * math.sin(f_t_pert)
    return {rho0_phy * g_acce + f_x_pert, f_y_pert, f_z_pert}
  else
    return {rho0_phy * g_acce, 0.0, 0.0}
  end
end
-----------------------End of Lua functions ------------------------------------

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'pipeLES'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Location of mesh files
mesh = './mesh/'
-- Logging output from simulation
logging = {
  level = 3,
  --filename = 'log' -- filename to write logging output
}
-- Scaling for multilevel simulation
scaling = 'acoustic'
-- Interpolation method for multilevel simulation
interpolation_method = 'quadratic'

NOdebug = {logging ={level=10, filename='dbg', root_only=false}}

--! [Simulation control]
sim_control = {
  time_control = {
    max = { sim=tmax_phy, clock=wall_clock },
    interval = interval_phy
  },
  abort_criteria = {
    stop_file = 'stop',
  }
}

--! [Physics parameters]
-- Required to convert physical unit to lattice unit
physics = {
  dt = dt,
  rho0 = rho0_phy
}
--! [Physics parameters]

--! [Scheme identifier]
identify = {
  label = '3D',
  kind = 'fluid_incompressible', -- Physics
  layout = 'd3q19',              -- Stencil
  relaxation = 'mrt'             -- Collision
}
--! [Scheme identifier]

--! [Fluid]
-- Fluid properties
fluid = {
  kinematic_viscosity = nu_phy,
  turbulence = {
    model = 'wale',
    c_w = 0.5 
  }
}
--! [Fluid]

--! [Initial condition]
initial_condition = {
  pressure = press_ambient,--press_analy,
  velocityX = vel_x_initial,
  velocityY = 0.0,
  velocityZ = 0.0,
}
--! [Initial condition]

--! [Boundary conditions]
-- Label is a boundary identifier and it should be same as in seeder
-- configuration
boundary_condition = {
  {
    label = 'pipe',
    kind = 'wall_libb', -- wall with q-values
  }
}
--! [Boundary conditions]

--! [Source]
if run_initial_period then
  glob_source = {
    force = force_till3T_star,
    force_order = 2
  }
else
  glob_source = {
    force = {rho0_phy * g_acce, 0.0, 0.0},
    force_order = 2
  }
end
--! [Source]

--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  {
    name = 'vel_ref',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = vel_x_initial
  },
  {
    name = 'press_avg',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'pressure_phy'},
      reduction_transient = {
        kind = 'average',
        nrecord = samp_avg_iter
      }
    }
  },
  {
    name = 'vel_avg',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'velocity_phy'},
      reduction_transient = {
        kind = 'average',
        nrecord = samp_avg_iter
      }
    }
  },
}
--! [User defined variables]

--! [Tracking]
tracking = {
  {
    label = 'vtk',
    folder = tracking_fol,
    variable = { 'pressure_phy', 'velocity_phy', 'vel_avg',
                 'turb_viscosity_phy' },
    shape = {kind = 'all'},
    time_control = {
      min = 0,
      max = tmax_phy,
      interval = T_star
    },
    output = {format = 'vtk'}
  },
--  {
--    label = 'hvs',
--    folder = tracking_fol,
--    variable = { 'press_avg', 'vel_avg', 'pressure_phy', 'velocity_phy' },
--    shape = {kind = 'all'},
--    time_control = {
--      min = 0,
--      max = tmax_phy,
--      interval = T_star
--    },
--    output = {format = 'harvester'}
--  },
  {
    label = 'wssBnd',
    folder = tracking_fol,
    variable = {'wss_phy'},
    shape = {
      kind = 'boundary',
      boundary = {'pipe'}
    },
    reduction = {'average'},
    time_control = {
      min= trac_start,
      max = tmax_phy,
      interval = 10*dt
    },
    output = {format = 'ascii'}
  },
  {
    label = 'probeAtCenter',
    folder = tracking_fol,
    variable = {'pressure_phy','velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5, pipe_center[2], pipe_center[3] }
      }
    },
    time_control = {
      min= trac_start,
      max = tmax_phy,
      interval = 10*dt
    },
    output = {format = 'ascii'}
  },
  {
    label = 'planeCenter',
    variable = { 'velocity_phy', 'vel_avg' },
    folder = tracking_fol,
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5, pipe_origin[2], pipe_origin[3] },
        vec = {
          { 0.0, diameter, 0.0 },
          { 0.0, 0.0, diameter }
        }
      }
    },
    time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = tmax_phy
    },
    output = {format = 'asciiSpatial'}
  },
}
--! [Tracking]

--! [Restart]
restart = {
  NOread = restart_fol..'pipeLES_lastHeader.lua',
  write = restart_fol,
  -- without timeControl restart will be dumped by default at end
  -- of simulation when write restart is set
  time_control = {
    min = rest_start,
    max = tmax_phy,
    interval = tmax_phy
  }
}
--! [Restart]
--------------------------- Musubi configuration -------------------------------
