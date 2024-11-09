-- Geometry information like length, width, height, dia_cyl and dx are loaded
-- from seeder configuration file because they are required by Seeder to
-- generate the mesh and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Reynolds number of the flow:
Re = 100
-- Mach number:
Ma = 0.1

-- Fluid properties --
-- Density of the fluid
rho0_phy = 1.0 -- [kg/m^3]

-- Maximum inflow velocity for parabolic profile:
vel_max_phy = 1.5 -- [m/s]

-- Computed values --
-- Mean velocity of the inflow
vel_mean_phy = vel_max_phy * 2/3

-- Mean inflow velocity computed from Ma number [m/s]
cs_phy = vel_mean_phy/Ma
-- Kinematic viscosity of the fluid according to Re [m^2/s]
nu_phy = dia_cyl * vel_mean_phy / Re

-- Ambient pressure
press_ambient = rho0_phy*cs_phy^2

-- Pressure drop across the channel length for Poiseuille flow
press_drop = 8 * vel_max_phy * rho0_phy * nu_phy * length / height^2
-- Pressure gradient
press_grad = press_drop / length

-- Drag and lift coefficient factors
force_coeff = 2 / (rho0_phy * vel_mean_phy * vel_mean_phy * dia_cyl * dx)
-- pressure coefficient
pressure_coeff = 2. / (rho0_phy * vel_mean_phy^2)

------------ Compute physical time step from lattice Mach number ---------------
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Lattice maximum velocity
vel_lat = Ma * cs_lat
-- Physical time step computed from physical and lattice velocity
dt = dx * vel_lat / vel_max_phy
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Chacterisitc time [s]
T_c = dia_cyl/vel_mean_phy
-- Physical simulation end time [s]
tmax_phy = 7.5
-- Ramping time for velocity [s]
t_ramp = 1.0
-- Physical simulation time to write output vtk files [s]
t_output = tmax_phy
-- Physical simulation time to write restart files [s]
t_restart = tmax_phy
-- Iterations of time to average [integer]
t_avg_interval = 10*math.ceil(T_c/dt) - 1
-- Physical sampling time with probe [s]
t_sampling = 10*dt
-- Interval to check status of the simulation [s]
interval_phy = dt*5000
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
------------------------- End of time settings ---------------------------------
--! [Local variables]

---------------------------- Lua functions -------------------------------------
--! [Analytical functions]
-- Analytical solutions are taken from 
-- L.D. LANDAU, E.M. LIFSHITZ, Fluid Mechanics (Second Edition),
-- CHAPTER II - VISCOUS FLUIDS,
-- Editor(s): L.D. LANDAU, E.M. LIFSHITZ,  1987
-- https://doi.org/10.1016/B978-0-08-033933-7.50010-6.
-- Analytical solution for velocity-x profile for Poiseuille flow
function vel_analy(x,y,z)
  return ( 0.5 / (rho0_phy*nu_phy) ) * press_grad * y * (height-y)
end
--! [Analytical functions]

--! [Inflow velocity profile]
-- Velocity profile for inflow boundary condition
function vel_inflow(x,y,z)
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
  label = '2D',
  layout = 'd2q9',               -- Stencil
  relaxation = 'bgk',            -- Collision
  kind = 'fluid_incompressible'  -- Physics
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
    label = 'west',
    kind = 'velocity_bounceback',
    velocity = {
      predefined = 'combined',
      temporal  = {
        predefined = 'smooth',
        min_factor = 0,
        max_factor = 1.0,  -- Factor to multiply with spatial
        from_time = 0,
        to_time = t_ramp -- ramp up to 1.0 sec
      },
      spatial = vel_inflow
    }
  },
  {
    label = 'east',
    kind = 'pressure_expol',
    pressure = press_ambient,
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
    label = 'cylinder',
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
    name = 'coeff_fac_force',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun =  force_coeff
  },
  -- Pressure coefficient factor
  {
    name = 'coeff_fac_pressure',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = pressure_coeff
  },

  {
    name = 'ambntPressure',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = press_ambient
  },
  -- Multiple force on boundary with coefficient factors
  {
    name  = 'coeff',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'multiply_scalar_times_vector',
      input_varname = { 'coeff_fac_force', 'bnd_force_phy' }
    }
  },
	{
    name  = 'press_diff',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'pressure_phy', 'ambntPressure' }
    }
  },
  {
    name  = 'coeffPressure', --Cp static
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'multiplication',
      input_varname = { 'coeff_fac_pressure', 'press_diff' }
    }
  },
  {
    name = 'reynolds_stress_diagonals',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'multiplication',
      input_varname = {'velocity_phy','velocity_phy'},
    }
  },
  {
    name = 'coeffPressureAvg',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'coeffPressure'},
      reduction_transient = {
        kind = 'average',
        nrecord = t_avg_interval --math.ceil(tmax_phy/dt)/4
      }
    }
  },
  { name='press_avg',
    ncomponents=1,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'pressure_phy'},
      reduction_transient = {
        kind = 'average',
        nrecord = t_avg_interval
      }
    }
  },
  { name='velocity_avg',
    ncomponents=3,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'velocity_phy'},
      reduction_transient = {
        kind = 'average',
        nrecord = t_avg_interval
      }
    }
  },
  { name='Re_stress_diag_avg', -- Reynolds stresses diagonals average,
    ncomponents=3,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'reynolds_stress_diagonals'},
      reduction_transient = {
        kind = 'average',
        nrecord = t_avg_interval
      }
    }
  },
}
--! [User defined variables]

--! [Tracking]
tracking = {
  -- Output file to visualize in Paraview.
  {
    label = 'vtk',
    folder = 'tracking/',
    variable = { 'pressure_phy', 'velocity_phy', 
                 'q_criterion_phy', 'vorticity_phy',
                 'velocity_avg', 'press_avg', 'Re_stress_diag_avg'
    },
    shape = { kind = 'all' },
    time_control = { min = tmax_phy, max = tmax_phy, interval = t_output },
    output = { format = 'vtk' }
  },
  -- Track pressure and velocity at the center of the channel over time.
  {
    label = 'probeAtCenter',
    folder = 'tracking/',
    variable = { 'pressure_phy', 'velocity_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { length*0.5 + x_offset, height*0.5 + y_offset, dx/2.0 }
      }
    },
    time_control = {
      min = trac_start,
      max = tmax_phy,
      interval = t_sampling
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  -- Track pressure before cylinder over time.
  {
    label = 'probePressAtCylBack',
    folder = 'tracking/',
    variable = { 'pressure_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { x_offset + pos_cyl[1]+rad_cyl+dx, y_offset + 0.2, dx/2.0 }
      }
    },
    time_control = {
      min = trac_start,
      max = tmax_phy,
      interval = t_sampling
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  -- Track pressure after cylinder over time.
  {
    label = 'probePressAtCylFront',
    folder = 'tracking/',
    variable = { 'pressure_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { x_offset + pos_cyl[1]-rad_cyl-dx, y_offset + 0.2, dx/2.0 }
      }
    },
    time_control = {
      min = trac_start,
      max = tmax_phy,
      interval = t_sampling
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  -- Track pressure, velocity along the center axis of the channel.
  -- Write this output only at the end of the simulation.
  {
    label = 'centerLine',
    folder = 'tracking/',
    variable = { 'pressure_phy','velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { x_offset, height*0.5 + y_offset, dx/2.0 },
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
  -- Track drag and lift coefficient on cylinder over time
  {
    label = 'cyl_force',
    folder = 'tracking/',
    variable = { 'coeff', 'bnd_force_phy' },
    shape = {
      kind = 'boundary',
      boundary = {'cylinder'}
    },
    time_control = {min = trac_start, max = tmax_p, interval = t_sampling},
    reduction = {'sum', 'sum'},
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  {
    label = 'Cp',
    folder = 'tracking/',
    variable = { 'coeffPressureAvg'},
    shape = {
      kind = 'boundary',
      boundary = {'cylinder'},
      cutoff_qvalue = 0.75 -- Cutoff boundary elements above this qvalue
    },
    time_control = { min = tmax_phy, max = tmax_phy, interval = 10*T_c },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'asciiSpatial', use_get_point=true}
  }
}
--! [Tracking]

--! [Restart]
-- Without timeControl restart will be dumped by default at end
-- of simulation when write restart is set
NOrestart = {
  NOread = 'restart/channel_lastHeader.lua',
  write = 'restart/',
}
--! [Restart]

--------------------------- Musubi configuration -------------------------------
