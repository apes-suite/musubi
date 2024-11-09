-- Geometry information like length, width, height, dx are loaded from seeder
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

-- Note that we have reduced the physical runtime of the case for our
-- recheck. Therefore, 'useRecheck' is set to true. If you do not want to
-- use the case for that, please deactivate this option by setting it to false.
-- This case reaches steady state at 32.645s.
useRecheck = true

--! [Local variables]
-- Flow parameters
-- Reynolds number of the flow relative to the height of the channel:
-- (Re = vel_mean_phy * height / nu_phy)
Re = 20  -- steady
-- Mach number
Ma = 0.1
-- Density of the fluid:
rho0_phy = 1.0 -- [kg/m^3]

-- Computed quantities:
-- Maximum inflow velocity for parabolic profile
vel_max_phy = 0.45
-- Mean velocity used for inflow [m/s]
vel_mean_phy = 4/9 * vel_max_phy
-- Speed of sound in the fluid:
cs_phy = vel_max_phy / Ma
-- Kinematic viscosity of the fluid:
nu_phy = (vel_mean_phy*height)/Re
-- Ambient pressure
press_ambient = rho0_phy*cs_phy^2

------------ Compute physical time step from lattice Mach number ---------------
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Lattice maximum velocity
vel_lat = Ma * cs_lat
-- Physical time step computed from physical and lattice velocity
dt = dx * cs_lat / cs_phy
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
if useRecheck then
  tmax_phy = 2.0
else
  tmax_phy = 50.0
end
-- Interval to check status of the simulation [s]
interval_phy = 1
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
------------------------- End of time settings ---------------------------------
--! [Global variables]

---------------------------- Lua functions -------------------------------------
--! [Inflow velocity profile]
-- Velocity profile in x for inflow boundary condition
-- Ref: "Benchmark Computations of Laminar Flow Around a Cylinder" - M. Schaefer
-- and S. Turek
function vel_profile(x,y,z)
  return 16 * vel_max_phy * y * z * ( height - y ) * ( height - z ) / height^4
end
-- Velocity for boundary condition
function vel_inflow(x,y,z,t)
  return {vel_profile(x,y,z), 0.0, 0.0} -- parabolic profile for inflow
end
--! [Inflow velocity profile]

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'C3D_Simple'
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
      nvals = 50,
      absolute = true,
      -- Condition to statisfy to every variable
      condition = {
        { threshold = 1.e-8, operator = '<=' },
        { threshold = 1.e-8, operator = '<=' }
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
  label = '3D',
  layout = 'd3q19',   -- Stencil
  relaxation = {
    name = 'bgk', -- Collision
    variant = 'standard_no_opt', -- a variant of collision
  },
  kind = 'fluid_incompressible' -- Physics
}
--! [Scheme identifier]


--! [Fluid]
fluid = {
  kinematic_viscosity = nu_phy
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
    label = 'west', --inlet
    kind = 'velocity_noneq_expol',
    velocity = vel_inflow
  },
  {
    label = 'east', --outlet
    kind = 'pressure_noneq_expol',
    pressure = press_ambient
  },
  {
    label = 'north',
    kind = 'wall',
  },
  {
    label = 'south',
    kind = 'wall',
  },
  {
    label = 'top',
    kind = 'wall',
  },
  {
    label = 'bottom',
    kind = 'wall',
  },
}
--! [Boundary conditions]


--! [Tracking]
tracking = {
  -- Output file to visualize in Paraview.
--  {
--    label = 'vtk',
--    folder = 'tracking/',
--    variable = { 'pressure_phy','velocity_phy' },
--    shape = { kind = 'all' },
--    time_control = { min = tmax_phy, max = tmax_phy, interval = tmax_phy },
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
        origin = {
          length*0.5 + x_offset,
          height*0.5 + y_offset,
          width*0.5  + z_offset
        }
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
    variable = { 'pressure_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, height*0.5, width*0.5 + z_offset },
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
    variable = {'velocity_phy', 'vel_mag_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, height*0.5, width*0.5  + z_offset },
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
    variable = { 'velocity_phy','vel_mag_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5, 0.0, width*0.5 + z_offset },
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
    --label = 'wssAlongHeight',
    folder = './tracking/',
    variable = { 'wss_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5, 0.0, width*0.5 + z_offset },
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
}
--! [Tracking]

--! [Restart]
restart = {
  NOread = 'restart/'..simulation_name..'_lastHeader.lua',
  write = 'restart/'----, add
}
--! [Restart]


--------------------------- Musubi configuration -------------------------------
