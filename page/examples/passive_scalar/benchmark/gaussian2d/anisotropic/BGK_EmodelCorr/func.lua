require 'args'

x_origin = 0.
y_origin = 0.

function gauss_pulse_real(x, y, z, t)

  sigma_xx = sigma0^2+2.*Dxx*t
  sigma_yy = sigma0^2+2.*Dyy*t
  sigma_xy = 2.*Dxy*t
  sigma = sigma_xx * sigma_yy - sigma_xy^2
  x_dev = x - x_origin - u_field[1]*t
  y_dev = y - y_origin - u_field[2]*t
  c = c0*sigma0^2/math.sqrt(sigma)*math.exp(
    -0.5*( x_dev^2 * sigma_yy / sigma + 
            y_dev^2 *sigma_xx / sigma - 
            2 * x_dev * y_dev * sigma_xy / sigma )
  ) + c_add
  return c
end