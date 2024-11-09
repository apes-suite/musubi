-- This is the configuration file for musubi. As it is a simple case, the mesh
-- is defined internally such that seeder is not needed.

-- To run the recheck, shepherd is set to true
shepherd = true

--------------------------- Create mesh internally -----------------------------
length = 10.0

if shepherd then
  refinementLevel = 4
else
  refinementLevel = 6
end

-- Here it is a simple cube, which can be defined internally
mesh = {
  predefined = 'cube',
  origin = {0.0, 0.0, 0.0},
  length = length,
  refinementLevel = refinementLevel
}

-- Calculate dx based on length and refinementLevel
dx = length / 2.0^refinementLevel
--------------------------------------------------------------------------------

--! [Local variables]
-- Flow parameters
-- Kinematic viscosity of the fluid [m^2/s]
--nu_phy = 1e-3
nu_phy = 0.01
-- Physical speed of sound [m/s]
cs_phy = 343.0
-- Density of the fluid [kg/m^3]
rho0_phy = 1.0
-- Lattice speed of sound
cs_lat = 1./math.sqrt(3.)
-- Ambient pressure
press_ambient = rho0_phy * cs_phy^2

------------ Compute physical time step from speed of sound ---- ---------------
dt = cs_lat /cs_phy * dx
-- Background mean flow velocity
vel_mean_phy = 0.0 * dx/dt
--------------------------------------------------------------------------------

-- Kinematic lattice viscosity
nu_lat = nu_phy * dt / dx^2
-- Relaxation parameter
omega = 1.0 / ( nu_lat/cs_lat^2.0 + 0.5 )

------------ Compute Mach number from velocity and speed of sound --------------
Ma = vel_mean_phy / cs_phy
-- Lattice maximum velocity
vel_lat = Ma * cs_lat
--------------------------------------------------------------------------------

--------------- Parameters for gaussian pulse in pressure ----------------------
center = { 5.0, 5.0, 5.0 }
centerX = center[1]
centerY = center[2]
centerZ = center[3]
-- Background pressure
background = press_ambient
halfwidth = 1.0
amplitude = 1.20

-- Function for 3D acoustic pulse
function ic_3Dgauss_pulse(x, y, z, t)
  r = ( x - centerX )^2 + ( y - centerY )^2 + ( z - centerZ )^2
  return background + amplitude * math.exp((-math.log(2.0)/(halfwidth^2))*r)
end

----------------------------- Time settings ------------------------------------
-- Physical simulation end time [s]
tmax_phy = 10.0
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Interval to check status of the simulation [s]
interval_phy = 1
---- Number of iterations used to run the simulation.
---- tmax_iter is also number of lattice iterations
--tmax_iter =  50
-- Starting time for tracking output
trac_start = 0
-- Starting time for restart output
rest_start = 0
------------------------- End of time settings ---------------------------------
--! [Global variables]

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'gaussianPulse'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- File to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Scaling for multilevel simulation
scaling = 'acoustic'
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
  layout = 'd3q19',   -- Stencil
  relaxation = 'bgk', -- Collision
  kind = 'fluid'      -- Physics
}
--! [Scheme identifier]

--! [Fluid]
fluid = {
  kinematic_viscosity = nu_phy,
  bulk_viscosity = nu_phy
}
--! [Fluid]

--! [Initial condition]
initial_condition = {
  pressure = ic_3Dgauss_pulse,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0
}
--! [Initial condition]

--! [Tracking]
tracking = {
  -- Track density, pressure and velocity over line at center.
  {
    label = 'pressAlongLength',
    folder = 'tracking/',
    variable = { 'density_phy', 'pressure_phy', 'velocity_phy'},
    shape = {
      kind = 'canoND',
      object = {
        origin = { 0.0, 0.5*length, 0.5*length },
        vec = { length, 0.0, 0.0 },
      }
     },
    time_control = {
      min = trac_start,
      max = tmax_phy,
      interval = tmax_phy
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'asciiSpatial'}
  },
  -- Track density, pressure and velocity at the center of the channel over time
  {
    label = 'probeAt1',
    folder = 'tracking/',
    variable = { 'density_phy', 'pressure_phy', 'velocity_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { 1.0, 5.0, 5.0 }
      }
     },
    time_control = {
      min = trac_start,
      max = tmax_phy,
      interval = {iter=1}
    },
    -- ascii format writes variable values and simulation time into a single
    -- file at every tracking interval. It is usally used for point tracking.
    output = {format = 'ascii'}
  },
}
--! [Tracking]
--
--! [Restart]
restart = {
  NOread = 'restart/gaussianPulse_lastHeader.lua',
  write = 'restart/',
}
--! [Restart]

--------------------------- Musubi configuration -------------------------------
