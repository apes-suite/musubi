-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Mach number
Ma = 0.05
-- Reynolds number of the flow
Re = 100
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Speed of sound [m/s]
cs_phy = 343.0

-- Mean inflow velocity computed from Reynolds number [m/s]
vel_mean_phy = Ma * cs_phy
-- Maximum inflow velocity for parabolic profile
vel_max_phy = 3.0 * vel_mean_phy / 2.0
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = vel_mean_phy * height / Re
-- Ambient pressure
press_ambient = rho0_phy * cs_phy^2

-- mass flow rate rho*vel_mag*A [kg/s] 
mfr= rho0_phy * vel_mean_phy * height * dx

-- Pressure drop across the channel length for Poiseuille flow
press_drop = 8 * vel_max_phy * rho0_phy * nu_phy * length / height^2
-- Pressure gradient
press_grad = press_drop / length

----------------- Compute physical time step from omega ------------------------
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Lattice velocity
vel_lat = Ma * cs_lat
-- Physical time step
dt = cs_lat / cs_phy * dx
-- Lattice viscosity
nu_lat = nu_phy*dt /dx^2
-- Relaxation parameter
omega = 1.0/(nu_lat/cs_lat^2 + 0.5)

--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
tmax_phy = 10000
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = 1
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
------------------------- End of time settings ---------------------------------
--! [Global variables]


---------------------------- Lua functions -------------------------------------
--! [Analytical functions]
-- Analytical solution for velocity-x profile for Poiseuille flow
function vel_analy(x,y,z)
  return ( 0.5 / (rho0_phy*nu_phy) ) * press_grad * y * (height-y)
end

-- Analytical solution for wall shear stress across the height
function wss_analy(x, y, z)
  return nu_phy * 2 * math.abs(y-height*0.5) * vel_max_phy / (height*0.5)^2.0
end

-- Analytical solution for strain rate across the height
function strainRate_analy(x, y, z)
  return -1*(y) * vel_max_phy / R / R
end

-- Analytical solution for pressure profile across the length
function press_analy(x,y,z,t)
  return press_ambient + press_drop * ( length-x ) / length
end
--! [Analytical functions]

--! [Inflow velocity profile]
-- Velocity profile for inflow boundary condition
function vel_inflow(x,y,z,t)
  return {vel_analy(x,y,z), 0.0, 0.0} -- parabolic profile for inflow
  --return {vel_mean_phy, 0.0,0.0}    -- block profile for inflow
end
--! [Inflow velocity profile]

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'channel'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Location of mesh files
mesh = './mesh/' 
-- Scaling for multilevel simulation
scaling = 'acoustic'
-- Debug outputs to write additional information
NOdebug = {
  logging = {
    level = 1,
    filename = 'dbg',
    root_only = false -- all involved MPI processes writes output
  }
}

-- Interpolation method for multilevel simulation
interpolation_method = 'quadratic'

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
    -- Abort if steady state is reached with condition defined in convergence
    -- table
    steady_state = true,
    -- Convergence condition for steady state check
    convergence = {
      -- Variables to check
      variable = {'pressure_phy','vel_mag_phy'},
      -- Check entire domain
      shape = {kind = 'all'},
      -- Reduce variables values in space to single average value
      reduction = {'average','average'},
      -- How often to do convergence?
      time_control = {
        min = 1,          -- Start convergence check after 1 [s]
        max = tmax_phy,   -- DO convergence until end of simulation
        interval = 10*dt  -- Do convernce check every 10*dt [s]
      },
      --
      norm='average',
      nvals = 100,
      absolute = true,
      -- Condition to statisfy to every variable
      condition = {
        { threshold = 1.e-10, operator = '<=' },
        { threshold = 1.e-10, operator = '<=' }
      }
    },
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
  label = '2D',
  layout = 'd2q9',    -- Stencil
  relaxation = {
    name = 'mrt', -- Collision
    variant = 'standard_no_opt', -- a variant of collision
  },
  kind = 'fluid'      -- Physics
}
--! [Scheme identifier]


--! [Fluid]
fluid = {
  kinematic_viscosity = nu_phy,
  bulk_viscosity = 2./3.*nu_phy
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
    kind = 'mfr_eq',
    mass_flowrate = mfr,
  },
  {
    label = 'east',
    kind = 'pressure_eq',
    pressure = 'press_analy',
  },
  {
    label = 'north',
    kind = 'wall',
  },
  {
    label = 'south',
    kind = 'wall',
  },
}
--! [Boundary conditions]

--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  -- Analytical velocity profile across the channel height
  {
    name = 'vel_analy',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = vel_analy
  },
  -- Analytical wall shear stress across the channel height
  {
    name = 'wss_analy',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = wss_analy
  },
  -- Analytical pressure profile across the channel length
  {
    name = 'press_analy',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = press_analy
  },
  -- Difference between numerical pressure and analytical pressure
  {
    name = 'press_diff',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'pressure_phy', 'press_analy' }
    }
  },
  -- Difference between numerical velocity and analytical velocity
  {
    name = 'vel_diff',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'vel_mag_phy', 'vel_analy' }
    }
  },
  -- Difference between numerical WSS and analytical WSS
  {
    name = 'wss_diff',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'wss_phy', 'wss_analy' }
    }
  },
}
--! [User defined variables]

--! [Tracking]
tracking = {
  -- Output file to visualize in Paraview.
  {
    label = 'DiffAlongHeight',
    folder = './tracking/',
    variable = { 'vel_diff', 'press_diff', 'wss_diff' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.85, 0.0, dx/2.0 },
        vec = { 0.0, height, 0.0},
      }
    },
    time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = tmax_phy
    },
    output = {format = 'asciiSpatial'}
  },
--  {
--    label = 'vtk',
--    folder = 'tracking/',
--    variable = { 'pressure_phy','velocity_phy', 'press_diff', 'vel_diff' },
--    shape = { kind = 'all' },
--    time_control = { min = 0, max = tmax_phy, interval = 4 },
--    output = { format = 'vtk' }
--  },
  -- Track pressure and velocity at the center of the channel over time.
  {
    label = 'probeAtCenter',
    folder = 'tracking/',
    variable = { 'pressure_phy', 'velocity_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5 + x_offset, height*0.5 + y_offset - dx/2, dx/2.0 }
      }
    },
    time_control = {
      min = 0,
      max = tmax_phy,
      interval = 10*dt
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  -- Track pressure profile along the center axis of the channel.
  -- Write this output only at the end of the simulation.
  {
    label = 'pressAlongLength',
    folder = 'tracking/',
    variable = { 'pressure_phy', 'press_analy', 'press_diff' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, height*0.5 - dx/2, dx/2.0 },
        vec = { length, 0.0, 0.0 },
      }
    },
    time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = tmax_phy
    },
    -- asciiSpatial format writes variable values and space coordinates into a
    -- seperate file for every tracking interval.
    -- It is usally used for line tracking.
    output = { format = 'asciiSpatial' }
  },
  -- Track velocity along the center axis of the channel.
  -- Write this output only at the end of the simulation.
  {
    label = 'velAlongLength',
    folder = 'tracking/',
    variable = {'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, height*0.5 - dx/2, dx/2.0 },
        vec = { length, 0.0, 0.0 },
      }
    },
    time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = tmax_phy
    },
    output = { format = 'asciiSpatial' }
  },
  -- Track velocity along the height of the channel at the middle of the channel
  -- length. Write this output only at the end of the simulation.
  {
    label = 'velAlongHeight',
    folder = './tracking/',
    variable = { 'velocity_phy','vel_mag_phy','vel_analy', 'vel_diff' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5, 0.0, dx/2.0 },
        vec = { 0.0, height, 0.0},
      }
    },
    time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = tmax_phy
    },
    output = {format = 'asciiSpatial'}
  },

  -- Track WSS along the height of the channel at the middle of the channel
  -- length. Write this output only at the end of the simulation.
  {
    label = 'wss_spatial',
    folder = './tracking/',
    variable = { 'wss_phy', 'wss_analy', 'wss_diff'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5, 0.0, dx/2.0 },
        vec = { 0.0, height, 0.0 },
      }
    },
    time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = tmax_phy
    },
    output = {format = 'asciiSpatial'}
  },
  -- Calculate L2-norm of press_diff variable at the end of the simulation
  {
    label = 'press_l2norm',
    folder = './tracking/',
    variable = {'press_diff','press_analy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, height*0.5 - dx/2, dx/2.0 },
        vec = { length, 0.0, 0.0 },
        segments = {200},
        distribution='equal'
      }
    },
    time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = tmax_phy
    },
    reduction = { 'l2norm', 'l2norm' },
    output = { format = 'asciiSpatial', use_get_point = true }
  },
  -- Calculate L2-norm of wss_diff variable at the end of the simulation
  {
    label = 'wss_l2norm',
    folder = './tracking/',
    variable = { 'wss_diff', 'wss_analy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5, 0.0, dx/2.0 },
        vec = { 0.0, height, 0.0 },
        segments = {200},
        distribution='equal'
      }
    },
    time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = tmax_phy
    },
    reduction = { 'l2norm', 'l2norm' },
    output = { format = 'asciiSpatial', use_get_point = true }
  },
  -- Calculate L2-norm of vel_diff variable at the end of the simulation
  {
    label = 'vel_l2norm',
    folder = './tracking/',
    variable = { 'vel_diff', 'vel_analy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5, 0.0, dx/2.0 },
        vec = { 0.0, height, 0.0 },
        segments = {200},
        distribution='equal'
      }
    },
    time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = tmax_phy
    },
    reduction = { 'l2norm', 'l2norm' },
    output = { format = 'asciiSpatial', use_get_point = true }
  },
}
--! [Tracking]

--! [Restart]
Norestart = {
  NOread = 'restart/channel2D_lastHeader.lua',
  write = 'restart/',
}
--! [Restart]


--------------------------- Musubi configuration -------------------------------
