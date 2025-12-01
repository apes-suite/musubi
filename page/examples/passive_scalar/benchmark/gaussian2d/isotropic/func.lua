require 'args'

x_origin = 0.
y_origin = 0.

-- The analytical solution for gaussian hill
function gauss_pulse_real(x, y, z, t)
  D = (tau - 0.5) / 3
  c = c0*sigma0^2/(sigma0^2+2.*D*t)*math.exp(
    -0.5*(( x - x_origin - u_field*t )^2+( y - y_origin )^2) / (sigma0^2 + 2*D*t)
  ) + c_add
  return c
end