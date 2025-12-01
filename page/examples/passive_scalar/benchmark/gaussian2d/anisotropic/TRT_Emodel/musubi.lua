require 'args'
require 'func'
-------------------------------------------------------------------------------
mesh               = './mesh/'  
cs2 = 1./3.

function gauss_pulse(x, y, z)
  c = c0*math.exp(-0.5*(( x-x_origin )^2+( y-y_origin )^2) / sigma0^2) + c_add
  return c*cs2
end
-------------------------------------------------------------------------------
sim_control        = {
  time_control     = {
    min      = { iter = 0        },
    max      = { iter = t_total     },
    interval = { iter = t_total / 10 }
  }
}
-------------------------------------------------------------------------------
identify  = {
  label = 'species', 
  kind = 'passive_scalar', 
  relaxation= {
    -- bgk, trt, mrt are supported for anisotropic diffusion
    name = 'trt',
    -- There are some models for anisotropic diffusion, e.g.,
    -- Emodel, EmodelCorr, Lmodel
    variant = 'Emodel'
  },
  -- d3q19 is the only layout supporting anisotropic diffusion
  layout='d3q19'
}

transport_velocity = 'velocity_fluid'

variable = {
  {
    name = 'velocity_fluid',
    ncomponents = 3,
    vartype = 'st_fun',
    st_fun = {u_field[1], u_field[2], 0.0}
  },
  {
    name = 'concentration_real',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = gauss_pulse_real,
  },
  {
    name = 'c_diff',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'difference',
      input_varname = {
        'concentration_real', 
        'spc1_density'
      },
    },
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
      min = { iter = t_total }, max = { iter = t_total }, interval = { iter = interval } 
    }
  }
}

field = { 
  label   = 'spc1',
  species = {
    -- diff_coeff controls the free parameter of diffusion, 
    -- which is set to the average of the tensor components
    -- if test fails e.g. with bgk model, try setting diff_coeff with the given tau
    -- i.e. diff_coeff = (tau - 0.5) / 3
    diff_coeff = (Dxx+Dyy)/3,
    -- diff_tensor sets the anisotropic diffusion tensor
    diff_tensor = {
      Dxx = Dxx,
      Dyy = Dyy,
      Dzz = Dzz,
      Dxy = Dxy
    }
  },
  initial_condition = { 
    pressure  = gauss_pulse,
    velocityX = u_field[1],
    velocityY = u_field[2],
    velocityZ = 0.0 
  }
}
