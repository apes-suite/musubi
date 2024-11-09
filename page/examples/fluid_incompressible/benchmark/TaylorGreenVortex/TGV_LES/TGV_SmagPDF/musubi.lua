-- This is the configuration file for musubi. The geometry is a simple cube,
-- therefore the mesh is created internally (during the musubi run).

--! [Mesh information]
----------------------------- Create mesh internally ----------------------------
-- Length of the cube
length = 2.0*math.pi
-- Refinement level
refinementLevel = 6

-- Here it is a simple cube, which can be defined internally
mesh = { predefined = 'cube',
         origin = {-length/2.0, -length/2.0, -length/2.0},
         length = length,
         refinementLevel = refinementLevel
}

-- Calculate dx based on length and refinementLevel
dx = length / 2.0^refinementLevel
----------------------------------------------------------------------------------
--! [Mesh information]

--! [Local variables]
------------------------------- Flow parameters --------------------------------
-- Now the main or actual musubi file starts.
--! [Local variables]
-- Flow parameters
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Reynolds number of the flow
Re = 1600
-- Velocity of the flow [m/s]
vel_max_phy = 1.0
-- Kinematic viscosity of the fluid [m^2/s] is calculated using the Re.
-- Characteristic velocity and length are both 1, such that Re reduces to the
-- inverse of kinematic viscosity:
nu_phy = 1.0 / Re
-- Maximum inflow velocity [m/s]
vel_max_phy = 1.0
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
--------------------------------------------------------------------------------

------------ Compute physical time step from lattice velocity ------------------
-- Lattice maximum velocity
vel_lat = 0.1
-- Physical time step computed from physical and lattice velocity
dt = dx * vel_lat / vel_max_phy
--------------------------------------------------------------------------------

------------------------ Compute remaining quantities --------------------------
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Mach number
Ma = vel_lat / cs_lat
-- Speed of sound [m/s]
cs_phy = vel_max_phy / Ma
-- Ambient pressure
press_ambient =  rho0_phy * cs_phy^2
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
tmax_phy = 10.0
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = tmax_phy / 10.0
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
-- wallclock is set to time which is used in sbatch-script
nhours     = 2
wall_clock = nhours*60*60 - 10*60
------------------------- End of time settings ---------------------------------
--! [Global variables]


---------------------------- Lua functions -------------------------------------
--! [Analytical functions for Taylor-Green Vortex]
-- Parameters needed for analytic solution of Taylor-Green Vortex
-- Velocity for analytic solution of Taylor-Green Vortex
u0 = vel_max_phy
-- Pressure for analytic solution of Taylor-Green Vortex
p0 = press_ambient
-- Characteristic time
T = 1.0 / ( 2* nu_phy ) -- equal to Re / 2.0

-- Analytical solution for velocity-x
function velocityX(x,y,z,t)
  return u0 * math.sin(x) * math.cos(y) * math.cos(z) * math.exp(-t/T)
end

-- Analytical solution for velocity-y
function velocityY(x,y,z,t)
  return -u0 * math.cos(x) * math.sin(y) * math.cos(z) * math.exp(-t/T)
end

-- Analytical solution for velocity-z
function velocityZ(x,y,z,t)
  return 0
end

-- Analytical solution for pressure
function pressure(x,y,z,t)
  p1 = math.cos(2*x) * (math.cos(2*z) + 2.0)
  p2 = math.cos(2*y) * (math.cos(2*z) + 2.0)
  return (p0 + (p1+p2)/16.0) * math.exp(-2.0*t/T)
end
--! [Analytical functions for Taylor-Green Vortex]

--! [Initial conditions based on analytical functions for Taylor-Green Vortex]
-- Initial condition for velocity-x
function ic_velocityX(x,y,z)
  return velocityX(x,y,z,0)
end

-- Initial condition for velocity-y
function ic_velocityY(x,y,z)
  return velocityY(x,y,z,0)
end

-- Initial condition for velocity-z
function ic_velocityZ(x,y,z)
  return velocityZ(x,y,z,0)
end

-- Initial condition for pressure
function ic_pressure(x,y,z)
  return pressure(x,y,z,0)
end

-- Additional IC
-- Sxx
function Sxx(x,y,z,t)
  sxx = math.cos(x) * math.cos(y) * math.cos(z)
  return sxx * math.exp(-t/T)
end
function ic_Sxx(x,y,z)
  return Sxx(x,y,z,0)
end
-- Syy
function Syy(x,y,z,t)
  return -Sxx(x,y,z,t)
end
function ic_Syy(x,y,z)
  return Syy(x,y,z,0)
end
-- Sxy
function Sxy(x,y,z,t)
  return 0.0  -- when kx == ky
end
function ic_Sxy(x,y,z)
  return Sxy(x,y,z,0)
end
-- Sxz
function Sxz(x,y,z,t)
  sxz = -0.5 * math.sin(x) * math.cos(y) * math.sin(z)
  return sxz * math.exp(-t/T)
end
function ic_Sxz(x,y,z)
  return Sxz(x,y,z,0)
end
-- Syz
function Syz(x,y,z,t)
  syz =  0.5 * math.cos(x) * math.sin(y) * math.sin(z)
  return syz * math.exp(-t/T)
end
function ic_Syz(x,y,z)
  return Syz(x,y,z,0)
end
------------ End of additional ones

--! [Initial conditions based on analytical functions for Taylor-Green Vortex]

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'TGV_LES_SmagPDF'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Scaling for multilevel simulation
scaling = 'acoustic'
-- Choose logging level
logging = { level = 3 }
-- Set IO-buffer size
io_buffer_size = 16

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
    max = { sim = tmax_phy, clock=wall_clock },
    interval = { sim = interval_phy }
  },
  abort_criteria = {
    stop_file = './stop'
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
  layout = 'd3q19',    -- Stencil
  relaxation = 'mrt', -- Collision
  kind = 'fluid_incompressible' -- Physics
}
--! [Scheme identifier]


--! [Fluid]
fluid = {
  kinematic_viscosity = nu_phy,
  bulk_viscosity = 2*nu_phy/3.0,
  turbulence = {
    model = 'smagorinsky', -- LES Smagorinsky model
    compute_strainrate_fromPDF = true, -- if false, then from velocity gradient
    c_s = 0.17, -- Smagorinsky constant, default value
    use_kolmogorov_scale = false -- for interpolation of turb eddy viscosity (ML)
  }
}
--! [Fluid]

--! [Initial condition]
initial_condition = {
  pressure = ic_pressure,
  velocityX = ic_velocityX,
  velocityY = ic_velocityY,
  velocityZ = ic_velocityZ,
  Sxx = ic_Sxx,
  Syy = ic_Syy,
  Szz = 0.0,
  Sxy = 0.0,
  Syz = ic_Syz,
  Sxz = ic_Sxz,
}
--! [Initial condition]

--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  { -- Square of strain rate
    name = 'strainSq',
    ncomponents = 6,
    vartype = 'operation',
    operation = {
      kind = 'multiplication',
      input_varname = {'strain_rate_phy', 'strain_rate_phy'}
    }
  },
  { -- Magnitude of strain rate squared
    name = 'strainSq_mag',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'magnitude',
      input_varname = {'strainSq'}
    }
  },
}
--! [User defined variables]

--! [Tracking]
tracking = {
--  -- Output file to visualize in Paraview.
--  {
--    label = 'vtk',
--    folder = 'vtkfiles/',
--    variable = { 'pressure_phy','velocity_phy', 'kinetic_energy_phy',
--                 'turb_viscosity_phy', 'q_criterion_phy' },
--    shape = { kind = 'all' },
--    time_control = { min= 0, max = tmax_phy, interval = tmax_phy/20.0 },
--    output = { format = 'vtk' }
--  },
  -- Track spatial average kinetic energy (KE)
  {
    label = 'kE_all',
    folder = 'tracking/',
    variable = { 'kinetic_energy_phy'},
    reduction = { 'sum' },
    shape = { kind = 'all' },
    time_control = {
      min = 0,
      max = tmax_phy,
      interval = 1/100
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
  -- Tracking for PSD
  {
    label = 'probeAtCenter',
    folder = 'tracking/',
    variable = { 'velocity_phy', 'pressure_phy', 'strainSq_mag' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0, 0, 0 }
      }
    },
    time_control = {
      min = 0,
      max = tmax_phy,
      interval = 1/100
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
}
--! [Tracking]

--! [Restart]
restart = {
  NOread = 'restart/'..simulation_name..'_lastHeader.lua',
  write = 'restart/',
  time_control = {
      min = tmax_phy,
      max = tmax_phy,
      interval = 1.0
    },
}
--! [Restart]

--------------------------- Musubi configuration -------------------------------
