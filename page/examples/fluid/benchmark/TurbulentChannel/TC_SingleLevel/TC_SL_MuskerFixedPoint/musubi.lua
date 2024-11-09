-- Geometry information like length, width, height, dx are loaded from seeder
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'
--! [Local variables]
-------------------------- Lattice quantities --------------------------------
-- Lattice density
rho0_L = 1.0
-- Lattice speed of sound
cs_L = math.sqrt(1.0/3.0)
-- Lattice bulk velocity
vel_bulk_L = Ma*cs_L
--------------------------------------------------------------------------------

------------ Compute physical time step from speed of sound --------------------
---- Physical timestep computed from physical and lattice speed of sound
dt = cs_L * dx / cs_phy
-- Lattice viscosity
nu_L = nu_phy*dt/dx^2
-- Relaxation parameter
omega = 1.0/(3.0*nu_L+0.5)
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- characteristics time or channel passage time
T_c = length/vel_bulk_phy
-- Maximum wall clock time
wallclock = 24*60*60 - 10*60
-- Physical simulation end time [s]
if shepherd then
  tmax_phy = 100*T_c --s
else
  tmax_phy = 600*T_c --s
end
-- Interval to check status of the simulation [s]
interval_phy = 5*T_c -- s
-- Sampling average iteration
samp_avg_iter = math.ceil(100*T_c/dt-5)
--------------------------------------------------------------------------------

bc_kind = 'turbulent_wall_noneq_expol'
curved = false
wall_function = 'musker' -- options: 'musker', 'reichardt', 'power_law'
nonlinear_solver = 'fixed_point' -- options: 'fixed_point', 'newton'
-- Boundary condition for north boundary
if useSymmetry then
  bc_north = 'symmetry'
else
  bc_north = bc_kind
end

---------------------------------- Lua functions ------------------------------
-- Paper: Simons, T., & Pletcher, R. (2013). Large eddy simulation of turbulent
-- flows using unstructured grids, (October). https://doi.org/10.2514/6.1998-3314
-- magnitude of perturbation
mag_pert = 1e-8*vel_bulk_phy--3.65e-5
PI = math.pi
-- perturbation in x-velocity
function u_x_perturbation(x,y,z)
  -- normalized distance
  x_star = x/length
  y_star = y/height
  z_star = z/width
  term1 = math.cos(2*PI*x_star)*math.sin(2*PI*y_star)
  term2 = 0.5*math.cos(4*PI*x_star)*math.sin(2*PI*y_star)
  term3 = math.cos(2*PI*x_star)*math.sin(4*PI*y_star)
  return mag_pert*length*math.sin(PI*z_star)*(term1+term2+term3)
end

-- perturbation in y-velocity
function u_y_perturbation(x,y,z)
  -- normalized distance
  x_star = x/length
  y_star = y/height
  z_star = z/width
  term1 = 0.5*math.sin(2*PI*x_star)*math.cos(2*PI*y_star)
  term2 = 0.5*math.sin(4*PI*x_star)*math.cos(2*PI*y_star)
  term3 = 0.25*math.sin(2*PI*x_star)*math.cos(4*PI*y_star)
  return -mag_pert*height*math.sin(PI*z_star)*(term1+term2+term3)
end

-- perturbation in z-velocity
function u_z_perturbation(x,y,z)
  -- normalized distance
  x_star = x/length
  y_star = y/height
  z_star = z/width
  term1 = 0.5*math.sin(2*PI*x_star)*math.sin(2*PI*y_star)
  term2 = 0.5*math.sin(4*PI*x_star)*math.sin(2*PI*y_star)
  term3 = 0.25*math.sin(2*PI*x_star)*math.sin(4*PI*y_star)
  return -mag_pert*width*(math.cos(PI*z_star))*(term1+term2+term3)
end

--karman constant
k = 0.41
function u_x_initial(x,y,z)
  if y<=halfHeight then
    y_plus = vel_fric_phy*y/nu_phy
  else
    y_plus = vel_fric_phy*(height-y)/nu_phy
  end
  if y_plus < 11.81 then
    u_plus = y_plus
  else
    u_plus = C_m * y_plus^m
  end
  --if y > halfHeight/4.0 and y < height-halfHeight/4.0 then
    return vel_fric_phy*u_plus + 0.15*(math.random()*2.0-1.0)*vel_bulk_phy
  --else
  --  return vel_fric_phy*u_plus
  --end
end

function u_initial_rand(x,y,z)
  --if y > halfHeight/4.0 and y < height-halfHeight/4.0 then
    return (math.random()*2.0-1.0)*vel_bulk_phy*0.15
  --else
  --  return 0
  --end
end
--! [Local variables]

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'Channel'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
mesh = 'mesh/' -- Mesh information
-- Scaling for multilevel simulation
scaling = 'acoustic'
-- Interpolation method for multilevel simulation
interpolation_method = 'quadratic'
-- Debug outputs to write additional information
NOdebug = {logging = {level=5, filename='dbg', root_only=false}}
-- Logging output from simulation
logging = {level=5}
comm_reduced = true

--! [Simulation control]
sim_control = {
  time_control = {
    max = { sim = tmax_phy, clock = wallclock },
    interval = { sim = interval_phy },
    check_iter = 10
   },
  abort_criteria = {
    stop_file = 'stop',
    velocity_lat_max = 0.2
  }
}
--! [Simulation control]

--! [Physics parameters]
-- Required to convert physical unit to lattice unit
physics = {
  dt = dt,
  rho0 = rho_phy
}
--! [Physics parameters]

--! [Scheme identifier]
identify = {
  kind = 'fluid',     -- Physics
  relaxation='mrt',   -- Collision
  layout='d3q19'      -- Stencil
}
--! [Scheme identifier]

--! [Fluid]
if identify.relaxation == "cumulant_extended" then
  fluid = {
    kinematic_viscosity = nu_phy,
    bulk_viscosity = 2*nu_phy,
    omega2 = -1., -- it uses omegaBulk. For a fix value use 1.! 
    -- if it fails reduce the limiter value by a factor of 10 till the code starts working
    omega_lim_1 = 0.01, -- this is the limiter for omega3, by default = 0.001
    omega_lim_2 = 0.01, -- this is the limiter for omega4, by default = 0.001
    omega_lim_3 = 0.01  -- this is the limiter for omega5, by default = 0.001 
  }
else
  fluid = {
    kinematic_viscosity = nu_phy,
    bulk_viscosity = 2*nu_phy/3.0,
    turbulence = {
      model = 'vreman',
      c_s = 0.17,--065, -- 0.17, -- Smagornisky constant
      c_v = 0.07, -- Constant for vreman model
      c_w = 0.5  -- Constant for WALE model
    },
    hrr_sigma = 0.98 -- dafualt = 0.98
  }
end


--! [Initial condition]
initial_condition = {
  pressure = press_ambient,
  velocityX = u_x_initial,
  velocityY = u_initial_rand,
  velocityZ = u_initial_rand
  --velocityY = 0.0,
  --velocityZ = 0.0,
  --velocityY = u_y_perturbation,
  --velocityZ = u_z_perturbation
}
--! [Initial condition]

--! [Boundary conditions]
-- Label is a boundary identifier and it should be same as in seeder
-- configuration
boundary_condition = {
  {
    label = 'north',
    kind = bc_north,
    curved = curved,
    wall_function = wall_function,
    nonlinear_solver = nonlinear_solver
  },
  {
    label = 'south',
    kind = bc_kind,
    curved = curved,
    wall_function = wall_function,
    nonlinear_solver = nonlinear_solver
  },
}
--! [Boundary conditions]

--! [Source term]
source = {
  turb_channel_force_accel = {accel, 0.0, 0.0},
  turb_channel_force_dynamic = {
    -- Shape to compute average bulk velocity from simulation
    shape_umean = {
      kind = 'canoND',
      object = {
        origin = {0.0, 0.0, 0.0}, -- In case of volume avg
        vec = {
          {length, 0.0, 0.0},
          {0.0, height, 0.0},
          {0, 0, width}
        }
      },
    },
    shape_utau = {
      kind = 'canoND',
      object = {
        origin = {0, dx/2.0, 0.0},
        vec = {
          {length,0.0,0.0},
          {0, 0, width}
        }
      },
    },
    ref_velocity_bulk = vel_bulk_phy,
    ref_height = halfHeight,
    flow_direction = 'x'
  }
}
--! [Source term]

variable = {
  {
    name = 'velX',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = {'velocity_phy'},
      input_varindex = {1}
    }
  },
  {
    name = 'velY',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = {'velocity_phy'},
      input_varindex = {2}
    }
  },
  {
    name = 'velZ',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = {'velocity_phy'},
      input_varindex = {3}
    }
  },
  {
    name = 'reynolds_stress_X',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'multiply_scalar_times_vector',
      input_varname = {'velX','velocity_phy'},
    }
  },
  {
    name = 'reynolds_stress_Y',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'multiply_scalar_times_vector',
      input_varname = {'velY','velocity_phy'},
    }
  },
  {
    name = 'reynolds_stress_Z',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'multiply_scalar_times_vector',
      input_varname = {'velZ','velocity_phy'},
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
  {
    name = 're_stress_avgX',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'reynolds_stress_X'},
      reduction_transient = {
        kind = 'average',
        nrecord = samp_avg_iter
      }
    }
  },
  {
    name = 're_stress_avgY',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'reynolds_stress_Y'},
      reduction_transient = {
        kind = 'average',
        nrecord = samp_avg_iter
      }
    }
  },
  {
    name = 're_stress_avgZ',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'reduction_transient',
      input_varname = {'reynolds_stress_Z'},
      reduction_transient = {
        kind = 'average',
        nrecord = samp_avg_iter
      }
    }
  },
}

-- Tracking
tracking = {
  {
    label = 'vtk',
    folder = 'tracking/',
    variable = {'pressure_phy','velocity_phy', 'press_avg','vel_avg'},
    shape = {
      kind = 'all'
    },
    time_control = {min= 0, max = tmax_phy, interval = 100*T_c},
    --time_control = {min= 0, max = tmax_phy, interval = 100*dt},
    output = {format = 'vtk'}
  },
--  {
--    label = 'hvs',
--    folder = 'tracking/',
--    variable = {'pressure_phy','velocity_phy', 'press_avg','vel_avg',
--                --'treeid'
--    },
--    shape = {
--      kind = 'all'
--    },
--    time_control = {min= tmax_phy, max = tmax_phy, interval = 20*T_c},
--    --time_control = {min= 0, max = tmax_phy, interval = 10*dt},
--    output = {format = 'harvester'}
--  },
  {
    label = 'centerPlaneAvg',
    folder = 'tracking/',
    variable = {'velocity_phy','vel_avg', 'pressure_phy', 'press_avg'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {0.0, halfHeight-dx/2.0, 0.0},
        vec = {
          {length, 0.0, 0.0},
          {0.0, 0.0, width}
        }
      }
    },
    reduction = {'average','average','average','average'},
    time_control = {min= 0, max = tmax_phy, interval = 100*dt},
    output = {format = 'ascii'}
  },
  {
    label = 'meanVel',
    folder = 'tracking/',
    variable = { 'velocity_phy', 'vel_avg'},
    shape = {
      kind = 'all',
    },
    reduction = {'average','average'},
    time_control = {min= 0, max = tmax_phy, interval = 100*dt},
    output = {format = 'ascii'}
  },
  {
    label = 'bcTurbWall',
    folder = 'tracking/',
    variable = { 'bc_fric_velocity_phy', 'bc_y_plus'},
    shape = {
      kind = 'boundary',
      boundary = {'south'}
    },
    reduction = {'average','average'},
    time_control = {min= 0, max = tmax_phy, interval = 100*dt},
    output = {format = 'ascii'}
  },
  {
    label = 'bndPlaneAvg',
    folder = 'tracking/',
    variable = {'velocity_phy','vel_avg', 'pressure_phy', 'press_avg','wss_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, dx/2.0, 0.0},
        vec = {
          {length, 0.0, 0.0},
          {0.0,0.0,width}
        }
      }
    },
    reduction = {'average','average','average','average','average'},
    time_control = {min= 0, max = tmax_phy, interval = 100*dt},
    output = {format = 'ascii'}
  },
  {
    label = 'half',
    folder = 'tracking/',
    variable = {'pressure_phy', 'velocity_phy', 'press_avg','vel_avg',
                're_stress_avgX', 're_stress_avgY', 're_stress_avgZ'
    },
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, 0.0, 0.0},
        vec = {
          {length, 0.0, 0.0},
          {0.0, halfHeight-dx/2.0, 0.0},
          {0.0,0.0,width}
        }
      }
    },
    --time_control = {min= 100*T_c, max = tmax_phy, interval = 10*T_c},
    time_control = {min = tmax_phy, max = tmax_phy, interval = tmax_phy},
    output = {format = 'asciiSpatial'}
  }
}
--! [Tracking]

if not shepherd and identify.relaxation ~= "cumulant_extended" then
  table.insert(tracking,
    {
      label = 'turbVisc',
      folder = 'tracking/',
      variable = {'turb_viscosity_phy'},
      shape = {
        kind = 'all'
      },
      time_control = {min= tmax_phy, max = tmax_phy, interval = 100*T_c},
      output = {format = 'vtk'}
    }
  )
end

--! [Restart]
-- Without timeControl restart will be dumped by default at end
-- of simulation when write restart is set
restart = {
  NOread = 'restart/Channel_lastHeader.lua',
  write = 'restart/',
}
--! [Restart]
--------------------------- Musubi configuration -------------------------------
