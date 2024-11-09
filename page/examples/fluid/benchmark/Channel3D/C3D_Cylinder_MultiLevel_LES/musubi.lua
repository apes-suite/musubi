-- Geometry information like length, width, height, dx are loaded from seeder 
-- configuration file because they are required by Seeder to generate the mesh
-- and they are also required for Musubi configuration.
require 'seeder'

--! [Local variables]
-- Flow parameters
-- Reynolds number of the flow
Re = 3900
-- Mach number
Ma = 0.10
-- speed of sound [m/s]
cs_phy = 343.20
-- Inflow velocity [m/s]
vel_phy = Ma * cs_phy
-- Kinematic viscosity of the fluid [m^2/s]
nu_phy = vel_phy * diameter / Re
-- Density of the fluid [kg/m^3]
rho_phy = 1.225

---------------------------- Lattice quantities --------------------------------
-- Lattice density
rho_lat = 1.0
-- Lattice speed of sound
cs_lat = math.sqrt(1.0/3.0)
-- Lattice reference pressure
p0_lat = rho_lat*cs_lat^2
--------------------------------------------------------------------------------

------------ Compute physical time step from lattice Mach number ---------------
-- Lattice inflow velocity
vel_lat = Ma * cs_lat
-- Physical time step computed from physical and lattice velocity
dt = vel_lat * dx / vel_phy
-- reference pressure
ref_press = rho_phy*cs_phy^2
--------------------------------------------------------------------------------

----------------------------- Time settings ------------------------------------
-- Chacterisitc time
T_c = diameter/vel_phy
-- Ramping time for velocity
t_ramp = 1.0
-- Physical simulation end time [s]
tmax_phy = 20.0
-- Number of iterations required to reach physical simulation end time.
-- tmax_iter is also number of lattice iterations
tmax_iter =  math.ceil(tmax_phy/dt)
-- Maximum wallclock time
tmax_clock = 24*60*60 - 10*60
-- Interval to check status of the simulation [s]
interval_phy = 0.01
-- Starting time for tracking output [s]
trac_start = 0
-- Starting time for restart output [s]
rest_start = 0
------------------------- End of time settings ---------------------------------

-------------------------- Lua functions ---------------------------------------
-- Sponge is defined as a viscosity function
sponge_length = 5*diameter
function viscosity_func(x,y,z)
  if ( x < 5*diameter ) then
    return nu_phy + (50*nu_phy-nu_phy)*(sponge_length-x)/sponge_length
  elseif ( x > (length-5*diameter) ) then
    return nu_phy + (50*nu_phy-nu_phy)*(x-(length-sponge_length))/sponge_length
  else
    return nu_phy
  end 
end
--! [Local variables]

--------------------------- Musubi configuration -------------------------------
-- Simulation name used in tracking and restart outputs
simulation_name = 'C3D_Cylinder'
-- Print runtime information like memory usage at the end of the simulation
printRuntimeInfo = false
-- file to write measurement time of simulation for each solver step
timing_file = 'mus_timing.res'
-- Location of mesh files
mesh = 'mesh/' -- Mesh information
--weights = 'mesh/weight'--'balance/channel3D_sphere_weight_t_'
--write_weights = ''--'mesh/weight'
-- Logging output from simulation
logging = {
  level = 5,
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
    max = { sim = tmax_phy, clock = tmax_clock },
    interval = { sim = interval_phy }
   },
  -- Abort conditions for simulation
  abort_criteria = {
    -- Abort if file with the name 'stop' exist
    stop_file = 'stop',
  }
}

--! [Physics parameters]
-- Required to convert physical unit to lattice unit
physics = {
  dt = dt,
  rho0 = rho_phy
}
--! [Physics parameters]

--! [Scheme identifier]
identify = {
  layout = 'd3q19',    -- Stencil
  relaxation = 'mrt',  -- Collision
  kind = 'fluid'         -- Physics
}
--! [Scheme identifier]


--! [Fluid]
fluid = {
  kinematic_viscosity = viscosity_func,
  bulk_viscosity = 2.0*nu_phy/3.0,
  turbulence = {
    model = 'vreman',
    c_s = 0.17,  -- Smagronisky constant
    c_v = 0.07,   -- Constant for vreman model
    c_w = 0.5    -- Constant for WALE model
  }
}
--! [Fluid]

--! [Initial condition]
initial_condition = {
  pressure = ref_press,
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
    kind = 'velocity_noneq_expol',
    velocity = {
      predefined = 'combined',
      temporal  = {
        predefined = 'smooth',
        min_factor = 0.0, max_factor = 1.0,
        from_time = 0, to_time = t_ramp
      },
      spatial = { vel_phy, 0.0, 0.0 }
    }
  },
  {
    label = 'east',
    kind = 'pressure_noneq_expol',
    pressure = ref_press
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
    label = 'front',
    kind = 'symmetry',
  },
  {
    label = 'back',
    kind = 'symmetry',
  },
  {
    label = 'cylinder',
    kind = 'wall_libb',
  }
}
--! [Boundary conditions]

--! [User defined variables]
-- Mainly used for tracking.
-- This variable can be refered to as variable in boundary condition and source
variable = {
  -- Lift and drag coefficient factors
  {
    name = 'coeff_fac',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun =  2.0 / (rho_phy * vel_phy * vel_phy * width * diameter)
  },
  -- Multiple force on boundary with coefficient factors
  {
    name = 'coeff',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'multiply_scalar_times_vector',
      input_varname = {'coeff_fac','bnd_force_phy'}
    }
  }
--  { name='press_avg',
--    ncomponents=1,
--    vartype = 'operation',
--    operation = {
--      kind = 'reduction_transient',
--      input_varname = {'pressure_phy'},
--      reduction_transient = {
--        kind = 'average',
--        nrecord = 14*T_c
--      }
--    }
--  },
--  { name='velocity_avg',
--    ncomponents=3,
--    vartype = 'operation',
--    operation = {
--      kind = 'reduction_transient',
--      input_varname = {'velocity_phy'},
--      reduction_transient = {
--        kind = 'average',
--        nrecord = 14*T_c
--      }
--    }
--  },
}
--! [User defined variables]

--! [Tracking]
tracking = {
  {
    label = 'vtk',
    folder = 'tracking/',
    variable = { 'pressure_phy', 'velocity_phy', 'turb_viscosity_phy',
                 'kine_viscosity_phy'
    },
    shape = { kind = 'all' },
    time_control = { min= 0, max = tmax_phy, interval = 0.25 },
    output = { format = 'vtk' }
  },
  {
    label = 'probe',
    folder = 'tracking/',
    variable = { 'pressure_phy', 'velocity_phy' },
    shape = {
      kind = 'canoND',
      object = {
        origin = { cylinder_x + 2*diameter, height/2.0 - dx/2.0, 0.0 },
      }
    },
    time_control = { min= 0, max = tmax_phy, interval = 10*dt },
    output = { format = 'ascii' }
  },
  {
    label = 'forceOnCyl',
    folder = 'tracking/',
    variable = { 'coeff' },
    shape = { kind = 'boundary', boundary = {'cylinder'} },
    time_control = {min= 0, max = tmax_phy, interval = 10*dt },
    reduction = { 'sum' },
    output = { format = 'ascii' },
  }
}
--! [Tracking]

--! [Restart]
-- Without timeControl restart will be dumped by default at end
-- of simulation when write restart is set
restart = {
  NOread = 'restart/C3D_Cylinder_lastHeader.lua',
  write = 'restart/',
}
--! [Restart]

--------------------------- Musubi configuration -------------------------------
