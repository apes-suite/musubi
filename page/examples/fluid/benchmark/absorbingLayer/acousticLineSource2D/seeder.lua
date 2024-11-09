--! [Global variables]
-- Thickness of the absorbing layer at outlet
abs_thickness = 0.8 --m
-- Height of the channel
height = 1.0
-- Length of the domain with absorbing layer at outlet
length = 4*height + abs_thickness
-- Number of elements in height
nHeight = 200
-- Element size
dx = height/nHeight
-- nElements in length
nLength = math.ceil(length/dx)
-- Level required to reach computed dx including one element for boundary at
-- both ends
level = math.ceil(math.log(nLength+2)/math.log(2))
-- Length of the bounding cube
length_bnd = (2^level)*dx
-- Smallest possible element size
dx_eps = length_bnd/2.0^20
--! [Global variables]

---------------------------Seeder configurations -------------------------------
-- Directory to write mesh
folder = 'mesh/'

-- file to write measurement time of mesh generation each step
timing_file = 'sdr_timing.res'

-- How detailed should the output be?
-- The higher the level, the more output you'll get from seeder.
logging = { level = 10 }

-- Debug outputs
NOdebug = {
  debugMode = true,
  debugFiles = true
}


-- Bounding_cube: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
  origin = {-dx, -dx, -dx},
  length = length_bnd
}

-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0.
minlevel = level

-- spatial object is composed  by the attribute and geometry

spatial_object = {
  {
   attribute = {
      kind = 'seed',        --- kind is seed/boundary/refinement/periodic
    },
    geometry = {
      kind = 'canoND',     --- canoND is nothing but the line/plane/point/box
      object = {
        origin = {length/2.0,height/2.0,dx/2.0}
        } ---object
      } --- geometry
  }, --- attribute
  {
    attribute = {
      kind = 'boundary',
      label = 'west',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {-dx_eps, -dx_eps, -dx_eps},
        vec = {
          {0.0, height + 2*dx_eps, 0.0},
          {0.0, 0.0, dx + 2*dx_eps}
        }
      }
    }
  },
  {
    attribute = {
      kind = 'boundary',
      label = 'east',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {length+dx_eps, -dx_eps, -dx_eps},
        vec = {
          {0.0, height + 2*dx_eps, 0.0},
          {0.0, 0.0, dx + 2*dx_eps}
        }
      }
    }
  },
  -- periodic at top and bottom
  {
    attribute = {
    kind = 'periodic',
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = {-dx_eps, -dx_eps, -dx_eps},
          vec = {
            {length+2*dx_eps, 0.0, 0.0},
            {0.0, 0.0, dx + 2*dx_eps}
          }
        },
        plane2 = {
          origin = {-dx_eps, height+dx_eps, -dx_eps },
          vec = {
            {0.0, 0.0, dx + 2*dx_eps},
            {length+dx_eps, 0.0, 0.0}
          }
        }
      } --- object
    } -- geometry
  }, --- attribute
  -- periodic at front and back
  {
    attribute = {
      kind = 'periodic',
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {          --- plane is composed by the a origin and two vectors
          origin = {-dx_eps, -dx_eps, dx+dx_eps },
          vec = {
            {length+2*dx_eps,0.0,0.0},
            {0.0,height+2*dx_eps,0.0}
          }
        },
        plane2 = {
          origin = {-dx_eps, -dx_eps,-dx_eps},
          vec = {
            {0.0,height+2*dx_eps,0.0},
            {length+2*dx_eps,0.0,0.0}
          }
        }
      } --- object
    } -- geometry
  } --- attribute
}  --- spatial object
