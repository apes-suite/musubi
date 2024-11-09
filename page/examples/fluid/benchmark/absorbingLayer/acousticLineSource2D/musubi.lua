-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = 1.49e-5
-- Physical speed of sound [m/s]
cs_phy = 343.2
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Lattice speed of sound
cs_lat = 1./math.sqrt(3.)
-- Damping factor for absorbing layer
damp_factor = 1.5 -- Maximum limit: 4.0/omega-0.001
-- Ambient pressure
press_ambient = rho0_phy * cs_phy^2

------------ Compute physical time step from speed of sound ---- ---------------
dt = cs_lat / cs_phy * dx
-- Background mean flow velocity
vel_mean_phy = 0.1 * dx/dt
--------------------------------------------------------------------------------

------------ Compute Mach number from velocity and speed of sound --------------
Ma = vel_mean_phy / cs_phy
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Time to reach outlet
T0 = length / cs_phy
-- Time of wave propagation of one wavelength
tmax = T0
------------------------- End of time settings ---------------------------------
--! [Local variables]

----------------------- Parameters for acoustic source ..........................
amplitude = 1e-3
background = rho0_phy
frequency = 2.0
Tp = 1.0/frequency

-- Pressure as line source
function bc_acousticLineSrc(x, y, z, t)
  return (background+amplitude*math.sin(2*math.pi*frequency*t*8/T0)) * cs_phy^2
end
-------------------------------------------------------------------------------

------------------ Absorbing layer plane as a lua function ---------------------
-- It is not used. Just provided as an exampe
westStart = -length/2.0+abs_thickness+dx/2.0
eastStart = length-abs_thickness-dx/2.0
southStart = -length/2.0+abs_thickness+dx/2.0
northStart = length/2.0-abs_thickness-dx/2.0
function absorbLayer_fun(x,y,z,t)
  fac = 0.0
  if x > eastStart then
    fac = 3125*(abs_thickness+eastStart-x)*(x-eastStart)^4/(256*(abs_thickness)^5)
  end

  sigma_s = damp_factor*fac
  return sigma_s
end
-------------------------------------------------------------------------------

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'acousticLineSrc'
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
      interval = tmax/50
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
  pressure = press_ambient, 
  velocityX = vel_mean_phy,
  velocityY = 0.0,
  velocityZ = 0.0
}
--! [Initial condition]


--! [Boundary condition]
boundary_condition = {
  {
    label = 'west',
    kind = 'pressure_noneq_expol',
    pressure = bc_acousticLineSrc,
  },
  {
    label = 'east',
    kind = 'pressure_noneq_expol',
    pressure = press_ambient, 
  },
}
--! [Boundary condition]

--! [Absorb layer as source term]
source = {
  --absorb_layer = absorbLayer_fun
  absorb_layer = 'absorblayer',
  absorb_layer_target = {
    pressure = press_ambient, 
    velocity = {vel_mean_phy, 0.,0.},
    --velocity = 'dynamic',
    --nrecord = math.ceil(abs_thickness/dx)
  }
}
--! [Absorb layer as source term]

variable = {
  {
    name = 'absorblayer',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      --fun = absorbLayer_fun,
      predefined = 'combined',
      temporal = 1.0,
      spatial = {
        predefined = 'spongelayer_plane',
        origin = {length-abs_thickness-dx/2.0, 0.0, dx/2.0},
        normal = {1.0,0.0,0.0},
        damp_profile = 'polynomial_n5',
        thickness = abs_thickness,
        damp_factor = damp_factor
      },
      shape = {
        inverted = false,
        kind = 'canoND',
        object = {
          origin = {length-abs_thickness-dx/2.0, 0.0, dx/2.0},
          vec = {
            {abs_thickness+dx, 0.0, 0.0},
            {0.0, height, 0.0}
          }
        }
      },
    }
  }
}

tracking = {
--  {
--    label = 'vtk',
--    variable = {'pressure_phy', 'absorblayer', 'density_phy', 'velocity_phy'},
--    --shape= { kind = 'all'},
--    shape = {
--      kind = 'canoND',
--      object = {
--        origin = {0.0, 0.0, dx/2.0},
--        vec = {
--          {length, 0.0, 0.0},
--          {0.0, height, 0.0}
--        }
--      },
--    },
--    output = {format='vtk'},
--    folder='tracking/',
--    time_control = {interval=tmax/20, min={iter= 0}, max =tmax}
--  },
  {
    label = 'probe',
    variable = {'pressure_phy', 'density_phy'},
    --shape= { kind = 'all'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {dx/2.0, height/2.0, dx/2.0},
      },
    },
    output = {format='ascii'},
    folder='tracking/',
    time_control = {interval=10*dt, min={iter= 0}, max =tmax}
  },
  {
    label = 'line',
    variable = {'pressure_phy', 'density_phy'},
    --shape= { kind = 'all'},
    shape = {
      kind = 'canoND',
      object = {
        origin = {dx/2.0, height/2.0, dx/2.0},
        vec = {length, 0.0, 0.0}
      },
    },
    output = {format='asciiSpatial'},
    folder='tracking/',
    time_control = {interval=tmax, min={iter= 0}, max =tmax}
  }

}
