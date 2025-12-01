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
physics  = { dt = dt, rho0 = rho_phy }
fluid    = { 
  rho0 = rho_phy, 
  kinematic_viscosity = nu_phy
}
identify = {
  label      = 'fluid_2D',
  kind       = 'fluid_incompressible',
  relaxation = 'bgk',
  layout     = 'd2q9'
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
    interval = {iter = tmax}
  }
}

