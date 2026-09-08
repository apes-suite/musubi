require 'params'
-- Analytical solutions 

function uxanaval(x,y,z,t)
  j = k1 * math.exp(math.sqrt(F0) * y) + k2 * math.exp(-math.sqrt(F0) * y)
  return j / rho_phy
end

function uanasol(x, y, z, t)
  return {uxanaval(x, y, z, t), 0, 0}
end

function utop(x, y, z, t)
  return {umax, 0, 0}
end

function coth(x)
  return (math.exp(x) + math.exp(-x)) / (math.exp(x) - math.exp(-x))
end