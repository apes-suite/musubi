----------------------- PLEASE READ THIS ---------------------------!!!

-- This input file is set up to run for regression check
-- Please make sure you DO NOT MODIFY AND PUSH it to the repository
                                                                          
--------------------------------------------------------------------!!!
-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Reynolds number
Re = 400
-- Speed of sound [m/s]
cs_phy = 343.0
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Inflow velocity of the cavity [m/s]
vel_phy = 25.0
-- Mach number computed from inflow velocity
Ma = vel_phy / cs_phy
-- Kinematic viscosity of the fluid [m^2/s] computed from Reynolds number
nu_phy = vel_phy * length / Re
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
-- Physical simulation end time [s]
tmax_phy = length/vel_phy * 10
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = 10000*dt
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = tmax_phy
------------------------- End of time settings ---------------------------------

---------------------------Output directories-----------------------------------
-- tracking folder
tracking_fol = './tracking/'
-- restart folder
restart_fol = 'restart/'
---------------------End of output directories----------------------------------
--! [Local variables]

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'lidcavity'
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

--! [Simulation control]
sim_control = {
  time_control = {
    max = tmax_phy,
    interval = interval_phy,
    check_iter = 10
  },
  -- Abort conditions for simulation
  abort_criteria = {
    -- Abort if file with the name 'stop' exist
    stop_file = 'stop',
    -- Abort if steady state is reached with condition defined in convergence
    -- table
    steady_state = true,
    convergence = {
      variable = { 'pressure_phy', 'vel_mag_phy' },
      shape = {
        kind = 'canoND',
        object = {
          origin = { length/2.0, length/2.0, 0.0 }
        }
      },
      time_control = { min = 0, max = tmax_p, interval = 10*dt },
      norm='average',
      nvals = 100,
      absolute = true,
      reduction = {'average','average'},
      condition = {
        { threshold = 1.e-8, operator = '<=' },
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
  kind = 'fluid_incompressible', -- Physics
  layout = 'd2q9',               -- Stencil
  relaxation = 'bgk'             -- Collision
}
--! [Scheme identifier]

--! [Fluid]
-- Fluid properties
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
    label = 'north',
    kind = 'velocity_bounceback',
    velocity = {
      predefined = 'combined',
      temporal  = {
        predefined = 'smooth',
        min_factor = 0,
        max_factor = 1.0,  -- Factor to multiply with spatial
        from_time = 0,
        to_time = 1.0 -- ramp up to 10 sec
      },
      spatial = { vel_phy, 0.0, 0.0 }
    }
  },
  {
    label = 'south',
    kind = 'wall'
  },
  {
    label = 'west',
    kind = 'wall'
  },
  {
    label = 'east',
    kind = 'wall'
  }
}
--! [Boundary conditions]

--! [Tracking]
tracking = {
  {
    label = 'vtk',
    folder = tracking_fol,
    variable = { 'pressure_phy','velocity_phy' },
    shape = { kind='all' },
    time_control = { min = tmax_phy, max = tmax_phy, interval = tmax_phy },
    output = { format = 'vtk' }
  },
  {
    label = 'probe',
    folder = tracking_fol,
    variable = { 'pressure_phy','velocity_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { length/2.0, length/2.0, width/2.0 }
      }
    },
    time_control = { min = trac_start, max = tmax_phy, interval = 100*dt },
    output = { format = 'ascii' }
  },
  {
    label = 'verticalLine',
    folder = tracking_fol,
    variable = {'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { length/2.0, 0., width/2.0 },
        vec = { 0., length, 0. },
        -- define segments to extract values using exact point.
        -- use_get_point in output table must be true.
        segments = nLength
      },
    },
    time_control = { min = tmax_phy, max = tmax_phy, interval = tmax_phy },
    output = { format = 'asciiSpatial', use_get_point = true },
  },
  {
    label = 'horizondalLine',
    folder = tracking_fol,
    variable = {'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, length/2.0, width/2.0},
        vec = { length, 0., 0.0},
        -- define segments to extract values using exact point.
        -- use_get_point in output table must be true.
        segments = nLength
      }
    },
    time_control = { min = tmax_phy, max = tmax_phy, interval = tmax_phy },
    output = { format = 'asciiSpatial', use_get_point = false },
  },
}
--! [Tracking]

--! [Restart]
restart = {
  ead = restart_fol..'lidcavity_lastHeader.lua',
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
