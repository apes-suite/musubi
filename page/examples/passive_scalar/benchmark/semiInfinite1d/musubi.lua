require 'params'
require 'func'
-------------------------------------------------------------------------------
mesh               = './mesh/'  
cs2 = 1./3.

-- initial dirichlet function profile (t_ini = 0)
function initial_concentration(x, y, z)
  if (bearSol(x, y, z, t_ini) > 1e-300) then
    return bearSol(x, y, z, t_ini) * cs2 * dx^2 / dt^2
  else
    return 1e-300 * cs2 * dx^2 / dt^2
  end
end

physics  = { dt = dt,  rho0 = 1. }

-------------------------------------------------------------------------------
sim_control        = {
  time_control     = {
    min      = { iter = 0        },
    max      = { iter = mStep     },
    interval = { iter = interval }
  }
}
-------------------------------------------------------------------------------
identify  = {
  label = 'species', 
  kind = 'passive_scalar', 
  relaxation={
    name = 'bgk',
    variant = 'first'
  }, 
  layout='d1q3',
}

transport_velocity = 'velocity_fluid'

glob_source = {
  varname = 'ps_sink',
  ps_sourceCoeff = 'lambda'
}

variable = {
  {
    name = 'velocity_fluid',
    ncomponents = 3,
    vartype = 'st_fun',
    st_fun = {v, 0.0, 0.0},
  },
  {
    name = 'concentration_real',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = bearSol,
  },
  {
    name = 'c_diff',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'difference',
      input_varname = {
        'concentration_real', 'spc1_density'
      },
    },
  },
  {
    name = 'lambda',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = -lambda
  },
  {
    name = 'sink',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'multiplication',
      input_varname = {
        'lambda', 'spc1_density'
      }
    }
  }
}

tracking  = { 
  { label     = 'spc1',
    variable  = {'spc1_density', 'concentration_real'},
    shape = {
      kind = 'all'
    },
    folder    = 'tracking/',
    output    = {format = 'asciispatial'},  
    time_control     = { 
      min = { iter = mStep }, max = { iter = mStep }, interval = { iter = mStep } }
  }
}

field = { 
  label   = 'spc1',
  species = {diff_coeff = D},
  initial_condition = { pressure  = initial_concentration,
                        velocityX = v,
                        velocityY = 0.0,
                        velocityZ = 0.0 
                      },
  boundary_condition = { 
                        { 
                          label = 'inlet', 
                          kind = 'pressure_antibounceback',
                          pressure = cs2 * dx^2 / dt^2,
                        },
                        { 
                          label = 'outlet', 
                          kind = 'pressure_antibounceback',
                          pressure = 1e-300 * cs2 * dx^2 / dt^2
                        }
                       }
}