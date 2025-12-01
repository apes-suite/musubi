require 'args'
require 'func'
-------------------------------------------------------------------------------
-- import mesh from ./mesh folder
mesh               = './mesh/'
cs2 = 1./3.

-- initial condition for the gaussian hill
function gauss_pulse(x, y, z)
  c = c0*math.exp(-0.5*(( x-x_origin )^2+( y-y_origin )^2) / sigma0^2) + c_add
  return c*cs2
end
-- print("amplititude of pulse: "..gauss_pulse(-50,0,0)/cs2)
-------------------------------------------------------------------------------
sim_control        = {
  time_control     = {
    min      = { iter = 0        },
    max      = { iter = t_total     },
    interval = { iter = t_total / 10 }
  }
}
-------------------------------------------------------------------------------
-- select collison scheme between bgk and trt by selecting "relaxation"
-- select order between first and second for different equilibria
identify  = {
  label = 'species',
  kind = 'passive_scalar',
  relaxation = {
    name = collision,
    variant = order
  },
  layout='d2q9'
}

transport_velocity = 'velocity_fluid'

variable = {
  {
    name = 'velocity_fluid',
    ncomponents = 3,
    vartype = 'st_fun',
    -- define background velocity for advection
    st_fun = {u_field, 0.0, 0.0},
  },
  {
    name = 'concentration_real',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = gauss_pulse_real,
  }
}
tracking  = {
  { label     = 'spc1',
    variable  = {'spc1_density', 'concentration_real'},
    shape = {
      -- kind = 'all'
      kind  = 'canoND',
      object  = {
        origin = {-nelem-1, 0, 0},
        vec = { {2.0*(nelem+1), 0., 0.0} }
      }
    },
    folder    = 'tracking/',
    output    = {format = 'asciispatial'},
    time_control     = {
      min = { iter = t_total }, max = { iter = t_total }, interval = { iter = t_total } }
  }
}

field = {
  label   = 'spc1',
  -- define a diffusion parameter
  species = {diff_coeff = {(tau-0.5)/3} },
  initial_condition = { pressure  = gauss_pulse,
                        velocityX = 0.0,
                        velocityY = 0.0,
                        velocityZ = 0.0
                      },
}

restart = {
  write = 'restart/',
  time_control = {
    min      = { iter = t_total },
    max      = { iter = t_total },
    interval = { iter = t_total }
  },
}
