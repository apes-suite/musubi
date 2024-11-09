-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Properties of air: https://en.wikipedia.org/wiki/Standard_sea_level
--
-- Dynamic viscosity of the fluid [Pa s]
mu = 1.789e-5 --Pa s
-- Density of the fluid [kg/m^3]
rho0_phy = 1.225 --kg/m^3
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = mu / rho0_phy
-- Inflow velocity [m/s]
vel_phy = 0.5

-- Reynolds number of the flow
Re = outer_dia_nozzle*vel_phy/nu_phy

-- Ambient pressure: 100kPa or 1atm
press_ambient = 1e5

------------ Compute physical time step from lattice Mach number ---------------
-- Lattice Mach number
Ma_lat = 0.04
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Lattice maximum velocity
vel_lat = Ma_lat * cs_lat
-- Physical time step computed from physical and lattice velocity
dt = dx * vel_lat / vel_phy
-- Musubi lua functions to get dx, dt and omega
--dx = getdxFromLevel( {len_bnd=length_bnd, level=level})
--dt = getdtFromVel( {dx = dx, u_p = vel_phy, u_l = vel_lat } )
--omega = getOmegaFromdt( {dx = dx, dt = dt, nu_p = nu_phy } )
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
tmax_phy = 5.0
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Simulation end time interms of CPU wallclock [s]
wallclock = 2*60*60-5*60
-- Interval to check status of the simulation [s]
interval_phy = 1000*dt
-- Ramping time
t_ramp = 2*l_ch / vel_phy
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
--tracking time_control for ascii output
trac_time = {min = trac_start, max= tmax_phy, interval = {iter=100}}
------------------------- End of time settings ---------------------------------
--! [Local variables]

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'nozzle'
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

-- Time step settings
--! [Simulation control]
sim_control = {
  time_control = {
    max = {sim = tmax_phy, clock = wallclock},
    interval = {sim = interval_phy}
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
      variable = {'normalized_pressure', 'velocity_phy'},
      -- Check center-axis of nozzle
      shape = {
        kind = 'canoND',
        object= {
          origin = { -l_nozzle/2.0, 0.0, z_pos },
          vec = { l_nozzle, 0.0, 0.0 }
        }
      },
      -- How often to do convergence check?
      time_control = {
        min = t_ramp,      -- Start convergence check after t_ramp
        max = tmax_phy,    -- DO convergence until end of simulation
        interval = 10*dt   -- Do convergence check every 10*dt [s]
      },
      -- Average values in center-axis
      reduction = { 'average', 'average' },
      norm ='average',
      nvals = 50,
      absolute = true,
      -- Condition to statisfy to every variable
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
  layout = 'd2q9',              -- Stencil
  relaxation = 'mrt',           -- Collision
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
    label = 'west',                  -- Boundary label defined in seeder.lua
    kind = 'velocity_noneq_expol',   -- Boundary condition kind for solver
    velocity = {
      predefined = 'combined',
      -- Ramp velocity from zero to vel_phy from 0 s to 0.25 s
      temporal = {
        predefined ='smooth',
        min_factor = 0.0,
        max_factor = 1.0,
        from_time = 0,
        to_time = t_ramp
      },
      spatial = { vel_phy, 0.0, 0.0}
    },
  },
  {
    label = 'east',
    kind = 'pressure_noneq_expol',
    pressure = press_ambient,
  },
  {
    label = 'nozzle',
    kind = 'wall_libb'
  }
}

--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  {
    name = 'ref_pressure',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = press_ambient
  },
  {
    name = 'normalized_pressure',
    ncomponents = 1,
    vartype = 'operation',
    operation = { 
      kind = 'difference',
      input_varname = {'pressure_phy','ref_pressure'}
    }
  },
  {
    name = 'cs_inv',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 1.0/cs_lat
  },
  {
    name = 'mach_nr',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'multiplication',
      input_varname = {'vel_mag', 'cs_inv'},
    }
  }
}
--! [User defined variables]

--! [Tracking]
tracking = {
  -- Harvester format to post-process after simulation.
  {
    label = 'hvs',
    folder = 'tracking/',
    variable = {
      'normalized_pressure',
      'velocity_phy',
      'shear_stress_phy',
      'vorticity_phy',
      'mach_nr'
    },
    shape = {kind = 'all'},
    time_control = {
      min = 0,
      max= tmax_phy,
      interval=tmax_phy/100,
    },
    output = {format = 'harvester'}
  },
  {
    label = 'hline',
    folder = 'tracking/',
    variable = {'normalized_pressure', 'velocity_phy'},
    shape = {
      kind='canoND',
      object = {
        origin ={ -inlet_2_nozzleCenter, 0.0, dx_c_half},
        vec = { l_ch, 0.0, 0.0 },
      }
    },
    time_control = {
      min = tmax_phy,
      max = {
        sim = tmax_phy,
        clock = wallclock,
      },
      interval = tmax_phy
    },
    -- asciiSpatial format writes variable values and space coordinates into a
    -- seperate file for every tracking interval.
    -- It is usally used for line tracking.
    output = {format = 'asciiSpatial'}
  },
  {
    label = 'vline',
    folder = 'tracking/',
    variable = {'normalized_pressure', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { nozzle_inner_dia_X, -inner_rad_nozzle, dx_c_half },
        vec = { inner_dia_nozzle, 0.0, 0.0 },
      }
    },
    time_control = {
      min = tmax_phy,
      max = {
        sim = tmax_phy,
        clock = wallclock,
      },
      interval = tmax_phy
    },
    -- asciiSpatial format writes variable values and space coordinates into a
    -- seperate file for every tracking interval.
    -- It is usally used for line tracking.
    output = {format = 'asciiSpatial'}
  },
  {
    label = 'probe_neck',
    folder = 'tracking/',
    variable = {'normalized_pressure', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { nozzle_inner_dia_X, 0.0, dx_c_half }
      }
    },
    time_control = trac_time,
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  {
    label = 'probe_inlet',
    folder = 'tracking/',
    variable = {'normalized_pressure', 'velocity_phy'},
    shape = {
      kind='canoND',
      object = {
        origin = { -inlet_2_nozzleCenter+3*dx_c_half, 0.0, dx_c_half }
      }
    },
    time_control = trac_time,
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  {
    label = 'probe_outlet',
    folder = 'tracking/',
    variable = {'normalized_pressure', 'velocity_phy'},
    shape = {
      kind='canoND',
      object={
        origin={ outlet_2_nozzleCenter/2.0, 0.0, dx_c_half }
      }
    },
    time_control = trac_time,
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  {
    label = 'probe_center',
    folder = 'tracking/',
    variable = {'normalized_pressure', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin={ 0.0, 0.0, dx_c_half }
      }
    },
    time_control = trac_time,
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  }
}
--! [Tracking]

--! [Restart]
-- Without timeControl restart will be dumped by default at end
-- of simulation when write restart is set
restart = {
   --read = restart_fol..'/nozzle_lastHeader.lua',
   write = restart_fol,
   time_control = { min = tmax_phy}
}
--! [Restart]

--------------------------- Musubi configuration -------------------------------
