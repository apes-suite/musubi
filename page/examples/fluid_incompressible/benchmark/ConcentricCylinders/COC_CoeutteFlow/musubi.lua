-- Geometry information like rad_outer, rad_inner and dx are loaded from seeder
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters --
-- Reynolds number of the flow:
Re = 10
-- Mach number:
Ma = 0.01

-- Fluid properties --
-- Speed of sound:
cs_phy = 343 -- [m/s]
-- Density of the fluid:
rho0_phy = 1.0 -- [kg/m^3]


-- Computed quantities --
-- Velocity for rotating inner cylinder [m/s]
vel_phy = Ma*cs_phy

-- Kinematic viscosity of the fluid:
nu_phy = vel_phy/Re * (rad_outer-rad_inner)

-- Angular velocity for rotating inner cylinder [rad/s]
angular_vel = vel_phy / rad_inner

-- Ambient pressure
press_ambient = rho0_phy * cs_phy^2

------------ Compute physical time step from lattice Mach number ---------------
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Lattice maximum velocity
vel_lat = Ma * cs_lat
-- Physical time step computed from physical and lattice velocity
dt = dx * vel_lat / vel_phy
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Physical simulation end time:
tmax_phy = 10 -- [s]

-- Interval to check status of the simulation [s]
interval_phy = 5
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
------------------------- End of time settings ---------------------------------
--! [Local variables]

---------------------------- Lua functions -------------------------------------
function vel_analy(x,y,z)
  r = math.sqrt( x^2 + y^2 )
  return vel_phy * ratio * (rad_outer/r - r/rad_outer)/(1-ratio^2)
end

function vel_inflow(x,y,z,t)
  return { -y*angular_vel, x*angular_vel, 0.0 }
end
------------------------- End of Lua functions ---------------------------------

---------------------------- Regions of interest -------------------------------
track_line = {
  kind = 'canoND',
  object = {
    origin = {rad_inner, 0.0, 0.0},
    vec = {rad_outer-rad_inner, 0.0, 0.0}
  }
}
--------------------- End of Regions of interest -------------------------------

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'concentricCylinder'
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
-- Debug outputs to write additional information
NOdebug = {
  logging = {
    level = 1,
    filename = 'dbg',
    root_only = false -- all involved MPI processes writes output
  }
}

--! [Simulation control]
sim_control = {
  time_control = {
    max = tmax_phy,
    interval = interval_phy
  },
  abort_criteria = {
    stop_file = 'stop',
    steady_state = true,
    convergence = {
      variable = { 'vel_mag_phy' },
      shape = track_line,
      time_control = { min = 0, max = tmax_phy, interval = 10*dt },
      reduction = { 'average' },
      norm='average',
      nvals = 50,
      absolute = true,
      condition = {
        { threshold = 1.e-8, operator = '<=' }
      }
    }
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
  kind = 'fluid_incompressible',     -- Physics
  layout = 'd2q9',                   -- Stencil
  relaxation = 'mrt',                -- Collision
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
    label = 'inner',
    kind = 'velocity_noneq_expol',
    velocity = vel_inflow,
    curved = true,
  },
  {
    label = 'outer',
    kind = 'velocity_noneq_expol',
    velocity = {0.0,0.0,0.0},
    curved =true,
  },
--  {
--    label = 'z_wall',
--    kind = 'wall'
--  }
}
--! [Boundary conditions]

--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be referred to as variable in boundary condition and source
variable = {
  {
    name = 'vel_an',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = vel_analy
  }
}
--! [User defined variables]

--! [Tracking]
tracking = {
  {
    label = 'probeAvg',
    folder = 'tracking/',
    --variable = {'pressure_phy','velocity_phy'},
    variable = {'vel_mag_phy','vel_an'},
    shape = track_line,
    time_control = {
      min = trac_start,
      max = tmax_phy,
      interval = 10*dt
    },
    reduction = {'average','average'},
    output = {format = 'ascii'}
  },
  {
    label = 'line',
    folder = 'tracking/',
    variable = {'velocity_phy','vel_mag_phy','vel_an'},
    shape = track_line,
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
--    variable = { 'pressure_phy', 'velocity_phy', 'vel_an' },
--    shape = { kind = 'all' },
--    time_control = { min= 0, max = tmax_phy, interval = interval_phy },
--    output = { format = 'vtk' }
--  },
}
--! [Tracking]

--! [Restart]
-- Without timeControl restart will be dumped by default at end
-- of simulation when write restart is set
restart = {
  NOread = 'restart/cirularChannel_lastHeader.lua',
  write = 'restart/',
}
--! [Restart]
