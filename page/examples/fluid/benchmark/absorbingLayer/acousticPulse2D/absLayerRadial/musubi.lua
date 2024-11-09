-- Geometry information like length, width, height, dx are loaded from seeder
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = 1.49e-5
-- Physical speed of sound [m/s]
cs_phy = 343.0
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Lattice speed of sound
cs_lat = 1./math.sqrt(3.)
-- Ambient pressure
press_ambient = rho0_phy * cs_phy^2

------------ Compute physical time step from speed of sound ---- ---------------
dt = cs_lat / cs_phy * dx
-- Background mean flow velocity
vel_mean_phy = 0.0 * dx/dt
--------------------------------------------------------------------------------

------------ Compute Mach number from velocity and speed of sound --------------
Ma = vel_mean_phy / cs_phy
-- Lattice maximum velocity
vel_lat = Ma * cs_lat
--------------------------------------------------------------------------------

nu_lat = nu_phy * dt / dx^2
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5 )
-- Damping factor for absorbing layer
damp_factor = 1.5 -- Maximum limit: 4.0/omega-0.001
----------------------------- Time settings ------------------------------------
tmax = 250*dt
------------------------- End of time settings ---------------------------------
--! [Local variables]

----------------------- Parameters for acoustic Pulse ..........................
centerX = 0.0
centerY = 0.
halfwidth = 1.0/20.
amplitude = 1e-3
background = rho0_phy

-- Function for 1D acoustic pulse
function ic_1Dgauss_pulse(x, y, z, t)
  r = ( x - centerX )^2
  return (background + amplitude * math.exp(-0.5/(halfwidth^2)*r)) * cs_phy^2
end

-- Function for 2D acoustic pulse
function ic_2Dgauss_pulse(x, y, z, t)
  r = ( x - centerX )^2+( y - centerY )^2
  return (background + amplitude * math.exp(-0.5/(halfwidth^2)*r)) * cs_phy^2
end
-------------------------------------------------------------------------------

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'pulse2D'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Location of mesh files
mesh = './mesh/'
-- Logging output from simulation
logging = {
  --level = 4,
  level = 3,
  --filename = 'log' -- filename to write logging output
}

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
    max = tmax,
    interval = tmax/10
  },
  abort_criteria = {
    velocity_lat_max = 0.2 -- Maximum lattice velocity permitted
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
  kind = 'fluid',     -- Physics
  relaxation = 'bgk', -- Collision
  layout = 'd2q9'     -- Stencil
}
--! [Scheme identifier]

--! [Fluid]
fluid = {
  kinematic_viscosity = nu_phy
}
--! [Fluid]


--! [Initial condition]
initial_condition = {
  pressure = ic_2Dgauss_pulse,
  velocityX = vel_mean_phy,
  velocityY = 0.0,
  velocityZ = 0.0
}
--! [Initial condition]

--! [Absorbing layer is defined as source term]
-- Absorb_layer variable should return absorbing strength profile.
-- Three shapes are provided as predefined spatial function in TreElm: plane,
-- box and radial.
source = {
  absorb_layer = 'absorblayer_radial',
  -- Define pressure and velocity in the absorbing layer in absorb_layer_target.
  -- Both variables can be defined either as constant or as string "dynamic".
  -- If defined as dynamic then nrecord must be defined to compute the
  -- smoothfactor time averaging using exponential moving average function.
  absorb_layer_target = {
    pressure = press_ambient,
    velocity = {vel_mean_phy,0.0,0.0},
    --velocity = 'dynamic',
    --nrecord = 4*abs_thickness/dx
  }
}
variable = {
  {
    name = 'absorblayer_radial',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      temporal = 1.0,
      spatial = {
        predefined = 'spongelayer_radial_2d',
        origin = {0.0, 0.0, dx/2.0},
        radius = length/2.0-abs_thickness-dx/2.0,
        damp_profile = 'polynomial_n5',
        thickness = abs_thickness,
        damp_factor = damp_factor,
      },
      shape = {
        inverted = true,
        kind = 'cylinder',
        object = {
          origin = {0.0,0.0,0.0},
          radius = length/2.0-abs_thickness-dx,
          vec = {0.0, 0.0, dx},
        }
      },
    },
  },
}

tracking = {
  {
    label = 'probe',
    variable = {'pressure_phy', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {-length/4.0, 0.0, dx/2.0},
      },
    },
    output = {format='ascii'},
    folder='tracking/',
    time_control = {interval=10*dt, min={iter= 0}, max = tmax}
  },
--  {
--    label = 'vtk',
--    variable = {'pressure_phy', 'absorblayer_radial',
--                'density_phy', 'velocity_phy'},
--    --shape= { kind = 'all'},
--    shape = {
--      kind = 'canoND',
--      object = {
--        origin = {-length/2.0, -length/2.0, dx/2.0},
--        vec = {
--          {length, 0.0, 0.0},
--          {0.0, length, 0.0}
--        }
--      },
--    },
--    output = {format='vtk'},
--    folder='tracking/',
--    time_control = {interval=tmax/10, min={iter= 0}, max = tmax}
--  }
}
