-- This setup illustrates the use of non-reflecting boundary conditions.
--
-- This kind of boundary condition is selected by 'outlet_nrbc'.
-- It requires the definition of the following parameters:
--
-- * kappa: TODO, add description
-- * sigma: TODO, add description
-- * length: TODO, add description
-- * mach_lat: TODO, add description

simulation_name = 'nrbc'
require "seeder"
omega = 1.8
model = 'fluid'
Re = 200
amplitude = 0.1
cs2LB = 1./3.
rho0LB = 1.
bcKind='outlet_nrbc'
verbose = true

--Scaling-----------------------------------------------------------------------
scaling = 'acoustic'

-- at this point in time, the reflection from the nrbc outlet has propagated
-- back to the tracking point to yield a maximum there.
tEnd = 5
u0 = 1.

-- length and level is defined in seeder.lua
viscosity = u0*length/Re
dx = length/2.^level
u0LB = 0.05

uLB = u0LB
dt = uLB/u0*dx
viscLB = viscosity*dt/dx/dx
omega = 1./(3.*viscLB + 0.5)

nIter = math.ceil(tEnd/dt)
p0 = 0
MaLB = uLB / math.sqrt(cs2LB)

-------------------------------------------------------------------------------
kind={}
kind['xM'] = 'wall'
kind['xP'] = 'wall'
kind['yM'] = 'wall'
kind['yP'] = 'wall'
kind['zM'] = 'wall'
kind['zP'] = 'wall'

-- checkDir is defined in seeder.lua
kind[checkDir] = bcKind

--
factorTrack = {}
factorTrack['xM'] = { 0.25, 0.5, 0.5 }
factorTrack['xP'] = { 0.75, 0.5, 0.5 }
factorTrack['yM'] = { 0.5, 0.25, 0.5 }
factorTrack['yP'] = { 0.5, 0.75, 0.5 }
factorTrack['zM'] = { 0.5, 0.5, 0.25 }
factorTrack['zP'] = { 0.5, 0.5, 0.75 }
trackPos = {
  origin[1] + factorTrack[checkDir][1]*length,
  origin[2] + factorTrack[checkDir][2]*length,
  origin[3] + factorTrack[checkDir][3]*length
}
print( 'tracking Position:  ' .. trackPos[1] .. ' '
                              .. trackPos[2] .. ' '
                              .. trackPos[3]
     )

rho0 =  1.0
kappa = 1.0
sigma = 0.1
LcharLB = 2^level

timing_file = 'mus_timing.res'
-- Simulation name
mesh = 'mesh/'-- Mesh information
-- Time step settigs
tmax   = nIter
interval = tmax/20
sim_control = {
  time_control = {max = tEnd}
}
interpolation_method = 'quadratic'
physics = { dt = dt, rho0 = 1.  }

identify = {
  kind       = 'fluid',      -- simulation type of this scheme
  relaxation = 'bgk', -- relaxation type (bgk, mrt, ...)
  layout     = 'd3q19'
}

originX = origin[1] + length/2.
originY = origin[2] + length/2.
originZ = origin[3] + length/2.
halfwidth = length/50.

function ic_1Dgauss_pulseX(x, y, z, t)
  return p0 + amplitude*math.exp(-0.5/(halfwidth^2)*( x - originX )^2)
end
function ic_1Dgauss_pulseY(x, y, z, t)
  return p0 + amplitude*math.exp(-0.5/(halfwidth^2)*( y - originY )^2)
end
function ic_1Dgauss_pulseZ(x, y, z, t)
  return p0 + amplitude*math.exp(-0.5/(halfwidth^2)*( z - originZ )^2)
end
function ic_2Dgauss_pulse(x, y, z, t)
  return p0 + amplitude*math.exp( -0.5 / (halfwidth^2)
                                       * ( (x - originX)^2 + (y - originY)^2)
                                )
end
function ic_3Dgauss_pulse(x, y, z, t)
  return p0 + amplitude*math.exp( -0.5 / (halfwidth^2)
                                       * ( (x - originX)^2 + (y - originY)^2
                                                           + (z - originZ)^2)
                                )
end

function gausspulse(x,y,z, t)
  if checkDir == 'xM' or checkDir == 'xP'  then
    val = ic_1Dgauss_pulseX(x,y,z)
  elseif checkDir == 'yM' or checkDir == 'yP'  then
    val = ic_1Dgauss_pulseY(x,y,z)
  elseif checkDir == 'zM' or checkDir == 'zP'  then
    val = ic_1Dgauss_pulseZ(x,y,z)
  end
  return val
end

-- Initial condition
initial_condition = {
  pressure = gausspulse,
  velocityX = 0.0,
  velocityY = 0.0,
  velocityZ = 0.0,
  Sxx = 0, Syy = 0, Szz = 0,
  Sxy = 0, Syz = 0, Sxz = 0,
}
-- outer omega cutoff ratio
w_min = 10.
-- inner omega cutoff ratio
w_max = 1.
cutoff_min = 0.45*length*0.5
cutoff_max = 0.5*length*0.5


function spatialFunction(x,y,z)
  if (x < cutoff_min) then
    res = w_max
  elseif (x  >= cutoff_max) then
    res =  w_min
  else
    slope = (w_max-w_min)/(cutoff_min-cutoff_max)
    res =  x*slope + w_max - cutoff_min*slope
  end
  return res
end


fluid = {
  kinematic_viscosity = viscosity,
  bulk_viscosity = 2*viscosity/3.0
}

boundary_condition = {
  {
    label = 'wall_xM',
    kind = kind['xM'],
    pressure = 'p0',
    kappa = kappa,
    sigma = sigma,
    length = LcharLB,
    mach_lat = MaLB
  },
  {
    label = 'wall_xP',
    kind = kind['xP'],
    pressure = 'p0',
    kappa = kappa,
    sigma = sigma,
    length = LcharLB,
    mach_lat = MaLB
  }
}

logging = {level =10}


-- user variables
variable = {
  {
    name = 'p0',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = p0
  }
}

-- Tracking
tracking = {
  {
    label    = 'probe_'..identify.kind,
    variable = {'pressure_phy','density',},
    folder   = 'tracking/',
    shape    = { kind   = 'canoND',
                 object = {origin = trackPos}
               },
    time_control = {min = 0, max = tmax, interval = dt},
    output       = {format = 'ascii'},
  },
  {
    label    = 'hvs',
    variable = {'pressure_phy','density'},
    folder   = 'tracking/',
    shape    = { kind   = 'canoND',
                 object = {origin = trackPos}
               },
    time_control = {min = 0, max = tmax, interval = interval},
    output       = {format = 'harvester'} }
}

if verbose then
  print('Scaling  '..scaling..' on level '..level)
  print('  Re       '..Re)
  print('  omega    '..omega)
  print('  uLB      '..uLB)
  print('  dx       '..dx)
  print('  dt       '..dt)
  print('  nIter    '..nIter)
end

restart = {
  write        = 'restart/',
  time_control = {min = 0., max = nIter, interval = interval}
}
