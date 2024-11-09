-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Mach number
Ma = 0.2
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = 1.49e-3
-- Physical speed of sound [m/s]
cs_phy = 343.2
-- Lattice speed of sound
cs_lat = 1./math.sqrt(3.)
-- Density of the fluid [kg/m^3]
rho0_phy = 1.2043

-- Reynolds number of the flow
Re = 150
-- Mean inflow velocity computed from Reynolds number [m/s]
vel_mean_phy = Re * nu_phy / dia_cyl

-- Ambient pressure
-- For outlet do nothing BC, the background pressure must be zero.
press_ambient = rho0_phy * cs_phy^2

-- Drag and lift coefficient factors
coeff_fac = 2 / (rho0_phy * vel_mean_phy * vel_mean_phy * dia_cyl * dx)

------------ Compute physical time step from speed of sound ---- ---------------
-- Physical time step
dt = cs_lat / cs_phy * dx
-- Lattice maximum velocity
vel_lat = Ma * cs_lat
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Ramping time
t_ramp = length/vel_mean_phy
-- Physical simulation end time [s]
if shepherd then
  tmax_phy = t_ramp
else
  tmax_phy = 5*t_ramp
end
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = dt*500
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
-- Sampling time press coeff

------------------------- End of time settings ---------------------------------
--! [Local variables]


---------------------------- Lua functions -------------------------------------
--! [Analytical functions]
-- Analytical solutions are taken from 
-- L.D. LANDAU, E.M. LIFSHITZ, Fluid Mechanics (Second Edition),
-- Circular sponge
sponge_start = height/2.0 - abs_thickness
sponge_length = abs_thickness

nu_lat = nu_phy * dt / dx^2
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5 )
damp_factor = 1.0 --maximum value: 4.0/omega - 0.001
thickness = abs_thickness
absorbStart = length/2.0-thickness-dx/2.

-- Velocity boundary condition
velocity_bc = {
  predefined = 'combined',
  NOtemporal  = {
    predefined = 'smooth',
    min_factor = 0,
    max_factor = 1.0,  -- Factor to multiply with spatial
    from_time = 0,
    to_time = t_ramp/2.0 -- ramp up to 1.0 sec
  },
  --spatial = vel_inflow
  temporal = 1.0,
  spatial = {vel_mean_phy, 0.0, 0.0}
}

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'channel'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Location of mesh files
mesh = './mesh/'
write_weights = './mesh/weight'
weights = write_weights
-- Logging output from simulation
logging = {
  level = 5, -- 10 max
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
    max = { sim = tmax_phy, clock = 24*60*60 - 10*60 },
    interval = { sim = interval_phy }
   },
  -- Abort conditions for simulation
  abort_criteria = {
    -- Abort if file with the name 'stop' exist
    stop_file = 'stop',
    velocity_lat_max = 0.25
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
  layout = 'd2q9',    -- Stencil --d2q9
  relaxation = 'mrt', -- Collision (bgk)
  kind = 'fluid'      -- Physics
}
--! [Scheme identifier]

--! [Fluid]
fluid = {
  kinematic_viscosity = {
    predefined = 'combined',
    NOtemporal  = {
      predefined = 'smooth',
      min_factor = 2.0,
      max_factor = 1.0,  -- Factor to multiply with spatial
      from_time = 0.0,
      to_time = t_ramp/4.0 -- Ramp until initial wave hits the boundary
    },
    --spatial = nu_phy,
    temporal = 1.0,
    spatial = {
      predefined = 'viscous_spongelayer_radial_2d',
      origin = {0.0, 0.0, 0.0},
      radius = sponge_start,
      thickness = sponge_length,
      damp_factor = 10.0,
      damp_exponent = 1.0,
      target_state = {
        viscosity = nu_phy
      }
    }
  },
  NOkinematic_viscosity = nu_phy,
  bulk_viscosity = 2 * nu_phy / 3
}
--! [Fluid]

--! [Initial condition]
initial_condition = {
  pressure = press_ambient,
  velocityX = vel_mean_phy,
  velocityY = 0.0,
  velocityZ = 0.0
}
--! [Initial condition]

--! [Boundary conditions]
-- Label is a boundary identifier and it should be same as in seeder
-- configuration
boundary_condition = {
  {
    label = 'inlet',
    kind = 'velocity_noneq_expol',
    velocity = velocity_bc
  },
  {
    label = 'outlet',
    kind = 'pressure_noneq_expol',
    pressure = press_ambient,
  },
  {
    label = 'north',
    kind = 'symmetry',
  },
  {
    label = 'south',
    kind = 'symmetry',
  },
  {
    label = 'cylinder',
    kind = 'wall_libb',
  },
}
--! [Boundary conditions]

--! [Source terms]
source = {
  absorb_layer = 'absorblayer_box',
  absorb_layer_target = {
    pressure = press_ambient,
    --velocity = 'dynamic',
    velocity = {vel_mean_phy, 0., 0.},
    nrecord = math.ceil(thickness/dx)
  },
}
--! [Source terms]

--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  {
    name = 'absorblayer_box',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      temporal = 1.0,
      spatial = {
        predefined = 'spongelayer_box_2d',
        origin = {-length/2.0+thickness+dx/2.0,
                  -height/2.0+thickness+dx/2.0,
                  0.0025-dx
        },
        extent = { length-2*(thickness+dx/2.0),
                   height-2*(thickness+dx/2.0),
                   2*dx
        },
        damp_profile = 'polynomial_n5',
        rounded_corner = true,
        corner_radius = 5*thickness,
        thickness = thickness,
        damp_factor = damp_factor,
      },
      shape = {
        inverted = true,
        kind = 'cylinder',
        object = {
          origin = {0.0, 0.0, 0.0025 - dx_half - dx_eps},
          radius = absorbStart - dx,
          vec = { 0.0, 0.0,  dx + 2*dx_eps}
        }
      },
    },
  },
--  {
--    name = 'absorblayer_radial',
--    ncomponents = 1,
--    vartype = 'st_fun',
--    st_fun = {
--      predefined = 'combined',
--      temporal = 1.0,
--      spatial = {
--        predefined = 'spongelayer_radial_2d',
--        origin = {0.0, 0.0, 0.0025},
--        radius = absorbStart,
--        damp_profile = 'polynomial_n5',
--        thickness = thickness,
--        damp_factor = damp_factor,
--      },
--      shape = {
--        inverted = true,
--        kind = 'cylinder',
--        object = {
--          origin = {0.0, 0.0, 0.0025 - dx_half - dx_eps},
--          radius = absorbStart - dx,
--          vec = { 0.0, 0.0,  dx + 2*dx_eps}
--        }
--      },
--    }
--  },
  -- Lift and drag coefficient factors
  {
    name = 'coeff_fac',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = coeff_fac
  },
  {
    name = 'pressure_ambient',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = press_ambient
  },
  {
    name = 'pressure_coeff_fac',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 1.0/(0.5*rho0_phy*vel_mean_phy*vel_mean_phy)
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
  },
  {
    name = 'pressure_diff',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'pressure_phy', 'pressure_ambient' }
    }
  },
  {
    name  = 'press_coeff',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'multiplication',
      input_varname = { 'pressure_coeff_fac', 'pressure_diff' }
    }
  },
  {
    name = 'press_coeff_avg',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'press_coeff'},
      reduction_transient = {
        kind = 'average',
        nrecord = 10000--math.ceil(t_ramp/dt)-10
      }
    }
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
        nrecord = 10000--math.ceil(t_ramp/dt)-10
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
        nrecord = 10000--math.ceil(t_ramp/dt)-10
      }
    }
  },
  {
    name = 'pfluc',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = { 'press_avg', 'pressure_phy' }
    }
  }
}
--! [User defined variables]

--! [Tracking]
tracking = {
  -- Output file to visualize in Paraview.
  {
    label = 'vtk',
    folder = 'tracking_vtk/',
    variable = { 'pressure_phy','velocity_phy', 'vorticity_phy',
                 'kine_viscosity_phy',
                 'absorblayer_box',
                 'press_avg', 'vel_avg', 'pfluc'},
    NOshape = { kind = 'all' },
    shape = { kind = 'canoND',
      object = {
        origin = {inlet_origin[1], inlet_origin[2], 0.0025 + 2*dx_eps },
        vec = {
          {length+2*dx, 0.0, 0.0},
          {0.0, height+2*dx, 0.0}
        }
      }
    },
    time_control = { min = 0, max = tmax_phy, interval = t_ramp/5.0 },
    --time_control = { min = 0, max = tmax_phy, interval = 1*dt },
    output = { format = 'vtk' }
  },
  -- Track pressure and velocity at the center of the channel over time.
  {
    label = 'probeAt90deg',
    folder = 'tracking/',
    variable = { 'pressure_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0, dia_cyl/2 + dx_cyl/2, 0.0025 }
      }
    },
    time_control = {
      min = trac_start,
      max = tmax_phy,
      interval = 10*dt
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  {
    label = 'probeAt0deg',
    folder = 'tracking/',
    variable = { 'pressure_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { -dia_cyl/2 - dx_cyl/2, 0, 0.0025 }
      }
    },
    time_control = {
      min = trac_start,
      max = tmax_phy,
      interval = 10*dt
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
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
    time_control = {min = trac_start, max = tmax_phy, interval = 10*dt},
    reduction = {'sum', 'sum'},
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },

  {
    label = 'Cp',
    folder = 'tracking/',
    variable = { 'press_coeff_avg'},
    shape = {
      kind = 'boundary',
      boundary = {'cylinder'}
    },
    time_control = {min = tmax_phy, max = tmax_phy, interval = tmax_phy},
    --reduction = {'sum', 'sum'},
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'asciiSpatial'}
  },
}
--! [Tracking]

--! [Restart]
-- Without timeControl restart will be dumped by default at end
-- of simulation when write restart is set
restart = {
  NOread = 'restart/channel_lastHeader.lua',
  write = 'restart/',
}
--! [Restart]

--------------------------- Musubi configuration -------------------------------
