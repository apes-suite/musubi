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
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = 1e-3
-- Speed of sound [m/s]
cs_phy = 343.0
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Ambient pressure
press_ambient = rho0_phy * cs_phy^2
-- Pressure drop in [Pa]
press_drop = 0.1
-- pressure at west boundary
press_west = press_ambient + press_drop
-- pressure at north east boundary
press_north_east = press_ambient
-- pressure at south east boundary
press_south_east = press_ambient

------------ Compute physical time step from fixed omega ------- ---------------
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Relaxation parameter omega
omega = 1.95
-- Lattice viscosity
nu_lat = ( 1.0 / omega - 0.5) * cs_lat^2
-- Physical time step computed from physical and lattice viscosity
dt = dx^2 * nu_lat / nu_phy
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
if recheck then
  tmax_phy = 2
else
  tmax_phy = 100
end
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = 1
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = tmax_phy
-- Termination wall clock time [s]
wall_clock = 4*60*60-10*60
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
simulation_name = 'pipe'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Location of mesh files
mesh = './mesh/'
-- Logging output from simulation
logging = {
  level = 5,
  --filename = 'log' -- filename to write logging output
}
-- Scaling for multilevel simulation
scaling = 'acoustic'
-- Interpolation method for multilevel simulation
interpolation_method = 'quadratic'

NOdebug = {logging ={level=1, filename='dbg', root_only=false}}

--! [Simulation control]
sim_control = {
  time_control = {
    max = { sim=tmax_phy, clock=wall_clock },
    interval = interval_phy
  },
  abort_criteria = {
    stop_file = 'stop',
    velocity_lat_max = 0.2
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
  layout = 'd3q27',              -- Stencil
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
  velocityZ = 0.0,
}
--! [Initial condition]

--! [Boundary conditions]
-- Label is a boundary identifier and it should be same as in seeder
-- configuration
boundary_condition = {
  {
    label = 'west',
    kind = 'pressure_noneq_expol',
    pressure = press_west,
  },
  {
    label = 'north_east',
    kind = 'pressure_noneq_expol',
    pressure = press_north_east
  },
  {
    label = 'south_east',
    kind = 'pressure_noneq_expol',
    pressure = press_south_east
  },
  {
    label = 'pipe',
    kind = 'wall_libb', -- wall with q-values
  }
}
--! [Boundary conditions]

--! [Tracking]
tracking = {
  {
    label = 'vtk',
    folder = tracking_fol,
    variable = { 'pressure_phy', 'velocity_phy' },
    shape = {kind = 'all'},
    time_control = {
      min = trac_start,
      max = tmax_phy,
      interval = interval_phy
    },
    output = {format = 'vtk'}
  },
  {
    label = 'probeAtCenter',
    folder = tracking_fol,
    variable = {'pressure_phy', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {0.0, 0.0, 0.0}
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
    label = 'probeAtWest',
    folder = tracking_fol,
    variable = {'pressure_phy', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {
          pipe_west_origin[1]+dx/2.0,
          pipe_west_origin[2],
          pipe_west_origin[3]
        },
        vec = {
          { 0.0, diameter, 0.0 },
          { 0.0, 0.0, diameter }
        }
      }
    },
    reduction = {'average', 'average'},
    time_control = {
      min= trac_start,
      max = tmax_phy,
      interval = 10*dt
    },
    output = {format = 'ascii'}
  },
  {
    label = 'probeAtNorthEast',
    folder = tracking_fol,
    variable = {'pressure_phy', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {
          pipe_northeast_origin[1]-dx/2.0,
          pipe_northeast_origin[2],
          pipe_northeast_origin[3]
        },
        vec = {
          { 0.0, diameter, 0.0 },
          { 0.0, 0.0, diameter }
        }
      }
    },
    reduction = {'average', 'average'},
    time_control = {
      min= trac_start,
      max = tmax_phy,
      interval = 10*dt
    },
    output = {format = 'ascii'}
  },
  {
    label = 'probeAtSouthEast',
    folder = tracking_fol,
    variable = {'pressure_phy', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {
          pipe_southeast_origin[1]-dx/2.0,
          pipe_southeast_origin[2],
          pipe_southeast_origin[3]
        },
        vec = {
          { 0.0, diameter, 0.0 },
          { 0.0, 0.0, diameter }
        }
      }
    },
    reduction = {'average', 'average'},
    time_control = {
      min= trac_start,
      max = tmax_phy,
      interval = 10*dt
    },
    output = {format = 'ascii'}
  },
}
--! [Tracking]

--! [Restart]
restart = {
  NOread = restart_fol..'pipe_lastHeader.lua',
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
