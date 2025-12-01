require 'params'
require 'func'
-------------------------------------------------------------------------------
mesh               = './mesh/'  
-------------------------------------------------------------------------------
sim_control        = {
  time_control     = {
    min      = { iter = 0        },
    max      = { iter = tmax     },
    interval = { iter = interval }
  }
}
-------------------------------------------------------------------------------
-- The best magic number for this case is given analytically in
--  Ginzburg, I.. doi: 10.1103/PhysRevE.77.066704.
-- Readers can try by uncommenting the following two lines.
-- B = math.sqrt(F0) / 2
-- magic = 3*(1-B^2*coth(B)^2) / (4*B^2*(1-3))

-- In general, we use magic = 0.25.
magic = 0.25
physics  = { dt    = dt,    rho0 = rho_phy }
fluid    = { 
  rho0 = rho_phy, 
  kinematic_viscosity = nu_phy, 
  lambda = magic 
}
identify = {
  label      = 'fluid_2D',
  kind       = 'fluid_incompressible',
  relaxation = 'trt',
  layout     = 'd3q19'
}
initial_condition = {
  pressure  = press_phy,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0
}
boundary_condition = {
  { label = 'bottom', kind = 'wall'},
  { 
    label = 'top', 
    kind      = 'velocity_bounceback',
    velocity  = utop,
  },
  { label     = 'inlet',
    kind      = 'velocity_bounceback',
    velocity  = uanasol,
  },
  { label     = 'outlet',
    kind      = 'velocity_bounceback',
    velocity  = uanasol,
  },
}
glob_source = {
  varname = 'channelForce',
  brinkman = Fc
}

variable = {
  {
    name = 'var_uanasol',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = uxanaval
  },
  {
    name = 'u_solve',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind='extract',
      input_varname={'velocity_phy'},
      input_varindex = {1}
    }
  },
  {
    name = 'u_diff',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind='difference',
      input_varname={'var_uanasol', 'u_solve'}
    }
  }
}
-------------------------------------------------------------------------------
tracking = {
  label = "F0_"..tostring(F0),
  variable = {
    'u_solve',
    'var_uanasol'
  },
  shape = {
    kind    = 'canoND',
    object  = {
      origin = { domainlen/2, 0, 0 },
      vec = { {0, Dia, 0.0}
        },
    },
  },
  folder = 'tracking/',
  output = {format = 'asciiSpatial'},
  time_control = {
    min = {iter = tmax},
    max = {iter = tmax},
    interval = {iter = interval}
  }
}

