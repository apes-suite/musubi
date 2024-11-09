--! [Global variables]
------------------------GEOMETRY PARAMETERS---------------------------
shepherd = true
if shepherd then
  -- length to height ratio
  l_h = 2.0
  -- width to height ratio
  w_h = 2.0
else
  -- length to height ratio
  l_h = 2.0*math.pi
  -- width to height ratio
  w_h = 2.0*math.pi
end
-- Height of the channel
halfHeight =  1.0 --m
height = 2*halfHeight
-- length of the channel
length = l_h*halfHeight
-- width of the channel
width = w_h*halfHeight
-- Set to true to simulate only half of channel height with symmetry BC
-- for north
useSymmetry = false
if useSymmetry then
  simHeight = halfHeight
else
  simHeight = height
end

-----------------------FLOW PARAMETERS---------------------------
-- Reference papers:
-- Malaspinas, O.; Sagaut, P. (2014): Wall model for large-eddy simulation based
-- on the lattice Boltzmann method. In Journal of Computational Physics 275,
-- pp. 25–40. DOI: 10.1016/j.jcp.2014.06.020.
--
-- Haussmann, Marc; BARRETO, Alejandro CLARO; KOUYI, Gislain LIPEME; Rivière,
-- Nicolas; Nirschl, Hermann; Krause, Mathias J. (2019): Large-eddy simulation
-- coupled with wall models for turbulent channel flows at high Reynolds numbers
-- with a lattice Boltzmann method — Application to Coriolis mass flowmeter.
-- In Computers & Mathematics with Applications 78 (10), pp. 3285–3302.
-- DOI: 10.1016/j.camwa.2019.04.033.
-- Velocity profile constants
C_m = 8.3
m = 1.0/7.0
-- Friction Reynolds number
Re_tau = 1000
-- bulk Reynolds number computing Eq.68 from
-- Malaspinas, O.; Sagaut, P. (2014): Wall model for large-eddy simulation based
-- on the lattice Boltzmann method. In Journal of Computational Physics 275, 
-- pp. 25–40. DOI: 10.1016/j.jcp.2014.06.020.
Re_bulk = (8.0/0.073)^(4.0/7.0)*Re_tau^(8.0/7.0)
-- mean bulk velocity
vel_bulk_phy = 1.0
-- Kinematic viscosity [m^2/s]
nu_phy = vel_bulk_phy * (2*halfHeight) / Re_bulk
---- Mach number
Ma = 0.1
-- density
rho_phy = 1.0 --kg/m^3
-- friction velocity
vel_fric_phy = Re_tau * nu_phy / halfHeight
-- characteristics length
delta_nu = nu_phy / vel_fric_phy
-- driving force acceleration
accel = vel_fric_phy * vel_fric_phy / halfHeight
-- Speed of sound
cs_phy = vel_bulk_phy / Ma
-- reference pressure
press_ambient = rho_phy*cs_phy^2

------------------------MESH PARAMETERS---------------------------
-- Normalized distance from boundary to first cell
y_plus = 50
-- element size
dx = y_plus * nu_phy / vel_fric_phy * 2
-- number of elements in height
nHeight = math.ceil(height/dx)
nHalfHeight = nHeight/2
-- number of elements in length
nLength = nHalfHeight*l_h
-- length of the channel
length = nLength*dx
-- number of elements in bounding cube
--=number of elements in channel + inlet + outlet
nLength_bnd = nLength+2
-- level required to reach computed dx
level = math.ceil(math.log(nLength_bnd)/math.log(2))
-- length of the bounding cube
length_bnd = (2^level)*dx
-- level of wall
level_wall = level + 0

-- smallest possible element size
dx_eps = length_bnd/2^20
dx_half = dx*0.5
zpos = dx_half
--! [Global variables]

---------------------------Seeder configurations -------------------------------
-- directory to write mesh
folder = 'mesh/'

-- File to write measurement time of mesh generation each step
timing_file = 'sdr_timing.res'

-- How detailed should the output be?
-- The higher the level, the more output you'll get from seeder.
logging = { level = 6 }

NOdebug = {debugMode = true, debugFiles=true}

-- bounding_cube: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {origin = {-dx/1.,-dx/1.,-dx/1.},
               length = length_bnd}

-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0.
minlevel = level

-- *********************** Table of spatial objects *********************** --
-- Each spatial object is defined by an attribute and some geometric entity
-- attached to this attribute. Attributes might be defined multiple times.
-- Attributes are described by a kind (boundary, seed or refinement), a level
-- and maybe further kind specific values, like a label for the boundary.
-- Geometric objects might by right now:
-- - canoND (point, line, plane or box)
-- - STL
-- - Sphere
-- - Cylinder
--
-- Periodic boundaries are special, spatial objects of this kind can only make
-- use of geometric objects of the kind 'periodic'.
spatial_object = {
-- Seed point
  {
    attribute = {
      kind = 'seed',  ------seed
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { dx_eps, height/4.0+dx_eps, dx_eps },
        vec = {
          {length-2*dx_eps, 0.0, 0.0},
          {0.0, 0.0, width-2*dx_eps}
        }
      } --- object
    }
  },
  {
    attribute = {
      kind = 'boundary',  ---kind in attribute is seed/boundary/refinement/periodic
      calc_dist = false,
      label = 'north',     -- for north
    },
    geometry = {
      kind = 'canoND',    -- kind in geometry is canoND/sphere/stl 
      object = {
        origin = { -dx_eps, simHeight+dx_eps, -dx_eps },
        vec = {
          {length+2*dx_eps, 0.0, 0.0},
          {0.0, 0.0, width+2.*dx_eps}
        }
      }
    }
  },
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false,
      label = 'south',   --- for south
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {-dx_eps, -dx_eps, -dx_eps},
        vec = {
          {length+2*dx_eps, 0.0, 0.0},
          {0.0, 0.0, width+2.*dx_eps}
        }
      }
    }
  },

  {
    attribute = {
      kind = 'periodic',
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1={
          origin = {-dx_eps, -dx_eps, -dx_eps},
          vec = {
            {0.0, 0.0, width+2*dx_eps},
            {0.0, height+2*dx_eps, 0.0}
          }
        },       
        plane2 = {
          origin = {length+dx_eps, -dx_eps, -dx_eps},
          vec = {
            {0.0,height+2*dx_eps,0.0},
            {0.0,0.0,width+2.*dx_eps}
          }
        }
      }  
    }
  },
  {
    attribute = {
      kind = 'periodic', --kind is periodic
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = {-dx_eps, -dx_eps, width+dx_eps},
          vec = {
            {length+2*dx_eps, 0.0, 0.0},
            {0.0, height+2*dx_eps, 0.0}
          }
        },  --- plane1
        plane2 = {
          origin = {-dx_eps, -dx_eps, -dx_eps},
          vec = {
            {0.0, height+2*dx_eps, 0.0},
            {length+2*dx_eps,0.0,0.0}
          }
        } --- plane2        
      }  
    }
  },
}
------------------- End of seeder configurations -------------------------------
