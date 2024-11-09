--! [Global variables]
-- Height of the channel [m]
height =  0.41
-- Width of the channel [m]
width = height
-- Number of elements in height
nHeight = 32
-- Element size
dx = height/nHeight
-- Half of element size
dx_half = dx*0.5
---- Length of the channel [m]
--length = nLength*dx
---- Length of the channel [m]
length = 2.5
-- Number of elements in length
nLength = math.ceil(length/dx)
-- Number of elements in bounding cube
-- = number of elements in channel + inlet + outlet
nLength_bnd = nLength + 2
-- Level required to reach computed dx
level = math.ceil(math.log(nLength_bnd)/math.log(2))
-- Length of the bounding cube
length_bnd = (2^level)*dx
-- Smallest possible element size
dx_eps = length_bnd/2^20
-- Offset in x-direction to control the origin and end of the channel length
x_offset = 0.0 --dx/2.0
-- Offset in y-direction to control the origin and end of the channel height
y_offset = 0.0 --dx/2.0
-- Offset in z-direction to control the origin and end of the channel width
z_offset = 0.0 --dx/2.0
--! [Global variables]

---------------------------Seeder configurations -------------------------------
-- Directory to write mesh
folder = 'mesh/'

-- File to write measurement time of mesh generation each step
timing_file = 'sdr_timing.res'

-- How detailed should the output be?
-- The higher the level, the more output you'll get from seeder.
logging = { level = 3 }

-- Debug outputs
NOdebug = {
  debugMode = true,
  debugFiles = true
}

-- Bounding_cube: two entries: origin and length in this order, if no keys are
-- used
bounding_cube = {
  -- Center is in left back corner
  origin = { -dx-x_offset, -dx-y_offset, -dx },
  length = length_bnd
}

-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0.
minlevel = level

-- *********************** Table of spatial objects ************************* --
-- Each spatial object is defined by an attribute and some geometric entity
-- attached to this attribute. Attributes might be defined multiple times.
-- Attributes are described by a kind (boundary, periodic, seed or refinement),
-- a level and maybe further kind specific values, like a label for the
-- boundary. Geometric objects might by right now:
-- - canoND (point, line, plane or box)
-- - stl
-- - sphere
-- - cylinder
-- - Spacer
-- - triangle
-- - ellipsoid
-- - periodic
--
-- Periodic boundaries are special spatial objects.
-- Attribute with kind = 'periodic' can only make use of geometric objects of
-- kind = 'periodic'.

spatial_object = {
--! [Seed point]
  {
    attribute = {
      kind = 'seed',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {
          length*0.5 + dx_eps + x_offset,
          height*0.5 + dx_eps + y_offset,
          width*0.5  + dx_eps + z_offset
        } -- origin
      } --- object
    }
  },
--! [Seed point]

--! [West boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false,
      label = 'west'     -- label for the boundary at x = -x_offset

    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -x_offset-dx_eps, -y_offset-dx_eps, -z_offset-dx_eps },
        vec = {
          { 0.0, height + 2*y_offset + 2*dx_eps, 0.0 },
          { 0.0, 0.0, width + 2*z_offset + 2*dx_eps }
        }
      }
    }
  },
--! [West boundary]

--! [East boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false,
      label = 'east'     -- label for the boundary at x = length + x_offset
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { length+x_offset+dx_eps, -y_offset-dx_eps, -z_offset-dx_eps },
        vec = {
          { 0.0, height + 2*y_offset + 2*dx_eps, 0.0 },
          { 0.0, 0.0, width + 2*z_offset + 2*dx_eps }
        }
      }
    }
  },
--! [East boundary]


--! [North boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false, -- to calculate q-values
      label = 'north',   -- label for the boundary at y = height + y_offset
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -x_offset-dx_eps, height+y_offset+dx_eps, -z_offset-dx_eps },
        vec = {
          { length + 2*x_offset + 2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, width + 2*z_offset + 2*dx_eps }
        }
      }
    }
  },
--! [North boundary]

--! [South boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false,
      label = 'south',   -- label for the boundary at y = -y_offset
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -x_offset-dx_eps, -y_offset-dx_eps, -z_offset-dx_eps },
        vec = {
          { length + 2*x_offset + 2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, width + 2*z_offset + 2*dx_eps }
        }
      }
    }
  },
--! [South boundary]

--! [Bottom boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false,
      label = 'bottom'     -- boundary at z = -z_offset, parallel to y

    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -x_offset-dx_eps, -y_offset-dx_eps, -z_offset-dx_eps },
        vec = {
          { length + 2*x_offset + 2*dx_eps, 0.0, 0.0 },
          { 0.0, height + 2*y_offset + 2*dx_eps, 0.0 }
        }
      }
    }
  },
--! [Bottom boundary]

--! [Top boundary]
  {
    attribute = {
      kind = 'boundary',
      calc_dist = false,
      label = 'top'     -- boundary at z = width + z_offset, parallel to y
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {
          -x_offset-dx_eps,
          -y_offset-dx_eps,
           width + 2*z_offset + 2*dx_eps
        },
        vec = {
          { length + 2*x_offset + 2*dx_eps, 0.0, 0.0 },
          { 0.0, height + 2*y_offset + 2*dx_eps, 0.0 }
        }
      }
    }
  },
} -- spatial object
------------------- End of seeder configurations -------------------------------
