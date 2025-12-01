-- Diffusion process inside a 2D cylinder. The concentration profile is compared
-- to the analytical solution. In this simulation, anti-bounce back boundary
-- condition is used.

require 'args'
-------------------------------------------------------------------------------
mesh               = './mesh/'
cs2 = 1./3.
c_init = 1.
maxN = 150.0
tau = 0.516
diff = (tau - 0.5) / 3
-- initial condition impacts the time-dependent result
tIni = 0  

function factorial(n)
  if n < 0.1 and n > -0.1 then
      return 1
  else
      return n * factorial(n - 1)
  end
end

-- 0th and 1st Bessel function is computed for the concentration profile
function bessel(param, order)
  local sum = 0.
  for iN = 0., maxN - 1 do
    sum = sum + (-1)^iN / (factorial(iN) * factorial(iN + order)) * (
      param / 2) ^ (2 * iN + order)
  end
  return sum
end

-- Analytical solution of the time evolution of the concentration profile
function concentration_real(x, y, z, t)
  local mu = {2.4048, 5.5201, 8.6537, 11.7915, 14.9309}
  local sum = 0
  local r = math.sqrt((x-nelem)^2 + (y-nelem)^2)

  if (r > nelem) then
    return cs2
  end

  for i = 1, #mu do
    sum = sum + 2 / (mu[i] * bessel(mu[i], 1)) * math.exp(-mu[i] ^ 2
      * diff * t / nelem ^ 2) * bessel(mu[i] * r / nelem, 0)
  end

  local c = c_init * (1 - sum)
  if c < 0.001 then
    return 0.001 * cs2
  end
  return c*cs2
end

-- Initial boundary condition can be assigned with analytical solution
-- at a certain time step.
function initial_concentration(x, y, z)
  return concentration_real(x, y, z, tIni)
end

function concentration_afterIni(x, y, z, t)
  return concentration_real(x, y, z, t + tIni) / cs2
end

-------------------------------------------------------------------------------
sim_control        = {
  time_control     = {
    min      = { iter = 0        },
    max      = { iter = t_total     },
    interval = { iter = interval }
  }
}
-------------------------------------------------------------------------------
identify  = {
  label = 'species',
  kind = 'passive_scalar',
  relaxation='trt', -- trt is used for a converged solution
  layout='d2q9',
  variant = 'second'
}

transport_velocity = 'velocity_fluid'

variable = {
  {
    name = 'velocity_fluid',
    ncomponents = 3,
    vartype = 'st_fun',
    st_fun = {0., 0.0, 0.0},
  },
  {
    name = 'concentration_real',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = concentration_afterIni,
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
  }
}
tracking  = {
  { label     = 'spc1',
    variable  = {'spc1_density', 'concentration_real'},
    shape = {
      -- kind = 'all'
      kind  = 'canoND',
      object  = {
        origin = {nelem, nelem, 5},
        vec = { {nelem, 0., 0.0} }
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
  species = {diff_coeff = {(tau-0.5)/3} },
  initial_condition = { pressure  = initial_concentration,
                        velocityX = 0.0,
                        velocityY = 0.0,
                        velocityZ = 0.0
                      },
  boundary_condition = {
                        {
                          label = 'circle',
                          -- anti-bounce back boundary condition is assigned
                          kind = 'pressure_antibounceback',
                          pressure = c_init * cs2,
                        }
                       }
}
