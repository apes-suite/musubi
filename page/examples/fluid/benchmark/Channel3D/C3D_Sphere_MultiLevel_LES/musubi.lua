-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Reynolds number of the flow
Re = 1000
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = 1e-3
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Mean inflow velocity computed from Reynolds number [m/s]
vel_phy = Re * nu_phy / diameter

-- Ambient pressure
press_ambient = 0.0

-- Drag and lift coefficient factors
-- Geier, M., Schönherr, M., Pasquali, A., & Krafczyk, M. (2015). 
-- The cumulant lattice Boltzmann equation in three dimensions: Theory and 
-- validation. Computers and Mathematics with Applications, 70(4), 507–547.
coeff_fac = 2 / (rho0_phy * vel_phy * vel_phy * math.pi * radius^2)

------------ Compute physical time step from lattice Mach number ---------------
-- Lattice Mach number
Ma_lat = 0.06
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Lattice maximum velocity
vel_lat = Ma_lat * cs_lat
-- Physical time step computed from physical and lattice velocity
dt = dx * vel_lat / vel_phy
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Ramping time for velocity
t_ramp = 5*dt--length/vel_phy/10.0
-- Physical simulation end time [s]
tmax_phy = 2*t_ramp
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = dt
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
------------------------- End of time settings ---------------------------------

-- Define velocity at west, north, south, top and bottom as space-time function
velocity_func = {
  predefined = 'combined',
  temporal  = {
    predefined = 'smooth',
    min_factor = 0,
    max_factor = 1.0,  -- Factor to multiply with spatial
    from_time = 0,
    to_time = t_ramp -- ramping time
  },
  spatial = {vel_phy, 0.0, 0.0}
}

-- Define viscosity as space-time function to apply sponge at outlet
sponge_length = 5*diameter
sponge_start_outlet = origin_east[1] - sponge_length - 5*dx
sponge_end_inlet = origin_west[1] + sponge_length
function viscosity_func(x,y,z,t)
  if x < sponge_end_inlet then
    return nu_phy + (sponge_end_inlet-x)/sponge_length*(50*nu_phy-nu_phy)
  elseif x > sponge_start_outlet then
    return nu_phy + (x-sponge_start_outlet)/(sponge_length-5*dx)*(50*nu_phy-nu_phy)
  elseif x >= origin_east[1] - 5*dx then
    return 50*nu_phy
  else
    return nu_phy
  end
end

--! [Local variables]

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'channel3D_sphere'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Location of mesh files
mesh = './mesh/'
--weights = 'mesh/weight'--'balance/channel3D_sphere_weight_t_'
--write_weights = ''--'mesh/weight'
-- Logging output from simulation
logging = {
  level = 3,
  --filename = 'log' -- filename to write logging output
}
-- Scaling for multilevel simulation
scaling = 'acoustic'
-- Interpolation method for multilevel simulation
interpolation_method = 'quadratic'
-- Debug outputs to write additional information
NOdebug = {
  logging = {
    level = 1,
    filename = 'dbg',
    root_only = false -- all involved MPI processes writes output
  }
}

-- Dynamic load balancing
NObalance = {
  dynamic = true,
  --weight = true,
  time_control = {
    min = t_ramp,
    max = t_ramp,
    interval = t_ramp/2,
    check_iter = 1,
  },
}

--! [Simulation control]
sim_control = {
  time_control = {
    max = { sim = tmax_phy },
    interval = { sim = interval_phy }
   },
  -- Abort conditions for simulation
  abort_criteria = {
    -- Abort if file with the name 'stop' exist
    stop_file = 'stop',
  }
}
--! [Simulation control]

--! [Physics parameters]
-- Required to convert physical unit to lattice unit
physics = {
  dt = dt,
  rho0 = rho0_phy
}
--! [Physics parameters]

--! [Scheme identifier]
identify = {
  layout = 'd3q19',    -- Stencil
  relaxation = 'mrt',  -- Collision
  kind = 'fluid'       -- Physics
}
--! [Scheme identifier]

--! [Fluid]
fluid = {
  kinematic_viscosity = viscosity_func,
  bulk_viscosity = 2.0*nu_phy/3.0,
  turbulence = {
    model = 'wale',
    c_w = 0.5
  }
}
--! [Fluid]

--! [Initial condition]
initial_condition = {
  pressure = press_ambient,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0
}
--! [Initial condition]

--! [Boundary conditions]
-- Label is a boundary identifier and it should be same as in seeder
-- configuration
boundary_condition = {
  {
    label = 'west',
    kind = 'velocity_noneq_expol',
    velocity = velocity_func
  },
  {
    label = 'east',
    kind = 'pressure_noneq_expol',
    pressure = press_ambient,
  },
  {
    label = 'north',
    kind = 'symmetry',
    velocity = velocity_func
  },
  {
    label = 'south',
    kind = 'symmetry',
    velocity = velocity_func,
  },
  {
    label = 'top',
    kind = 'symmetry',
    velocity = velocity_func,
  },
  {
    label = 'bottom',
    kind = 'symmetry',
    velocity = velocity_func,
  },
  {
    label = 'sphere',
    kind = 'wall_libb',
  },
}
--! [Boundary conditions]

--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  -- Lift and drag coefficient factors
  {
    name = 'coeff_fac',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = coeff_fac
  },
  -- Multiple force on boundary with coefficient factors
  {
    name  = 'coeff',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'multiply_scalar_times_vector',
      input_varname = { 'coeff_fac', 'bnd_force_phy' }
    }
  }
}
--! [User defined variables]

--! [Tracking]
tracking = {
  -- Output file to visualize in Paraview.
  {
    label = 'vtk',
    folder = 'tracking/',
    variable = { 'pressure_phy','velocity_phy', 'vorticity_phy',
                 'q_criterion_phy', 'turb_viscosity_phy',
                 'kine_viscosity_phy', 'process'
    },
    --shape = { kind = 'all' },
    shape = {
      kind = 'canoND',
      object = {
        origin = {origin_west[1], origin_west[2], 0.0},
        vec = {
          {length+2*dx, 0.0, 0.0},
          {0.0, height+2*dx, 0.0, 0.0}
        }
      }
    },
    time_control = { min = 0, max = tmax_phy, interval = 1.0/10 },
    output = { format = 'vtk' }
  },
  -- Track pressure and velocity at the center of the channel over time.
--  {
--    label = 'probeAtCenter',
--    folder = 'tracking/',
--    variable = { 'pressure_phy', 'velocity_phy' },
--    shape = {
--      kind = 'canoND',
--      object = {
--        origin = { origin_west[1]+length*0.5 , 0.0, 0.0 }
--      }
--    },
--    time_control = {
--      min = trac_start,
--      max = tmax_phy,
--      interval = 10*dt
--    },
--    -- ascii format writes variable values and simulation time into a single
--    -- file at every tracking interval. It is usally used for point tracking.
--    output = {format = 'ascii'}
--  },
--  -- Track pressure, velocity along the center axis of the channel.
--  -- Write this output only at the end of the simulation.
--  {
--    label = 'centerLine',
--    folder = 'tracking/',
--    variable = { 'pressure_phy','velocity_phy'},
--    shape = {
--      kind = 'canoND',
--      object = {
--        origin = { origin_west[1]+length*0.5 , 0.0, 0.0 },
--        vec = { length, 0.0, 0.0 },
--      }
--    },
--    time_control = {
--      min = tmax_phy,
--      max = tmax_phy,
--      interval = tmax_phy
--    },
--    -- asciiSpatial format writes variable values and space coordinates into a
--    -- seperate file for every tracking interval.
--    -- It is usally used for line tracking.
--    output = { format = 'asciiSpatial' }
--  },
  -- Track drag and lift coefficient on cylinder over time
  {
    label = 'force',
    folder = 'tracking/',
    variable = { 'coeff', 'bnd_force_phy' },
    shape = {
      kind = 'boundary',
      boundary = {'sphere'}
    },
    time_control = {min = trac_start, max = tmax_p, interval = 2*dt},
    reduction = {'min', 'min'},
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
}
--! [Tracking]

--! [Restart]
-- Without timeControl restart will be dumped by default at end
-- of simulation when write restart is set
restart = {
  NOread = 'restart/channel3D_sphere_lastHeader.lua',
  write = 'restart/',
  time_control = {min = trac_start, max = tmax_p, interval = 1.0},
}
--! [Restart]

--------------------------- Musubi configuration -------------------------------
