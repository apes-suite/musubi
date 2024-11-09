--! [Global variables]
-- cylinder radius [m]
radius = 0.05
-- cylinder diameter [m]
diameter = radius*2.0
-- Channel height [m]
height =  10*diameter
-- Channel width [m]
width = 4.0*diameter
-- Channel length [m]
length = 30*diameter
-- coarsest element size is fixed to contain 10 elements in z-direction
nDiameter = 10
-- Element size
dx = diameter/nDiameter
dx_half = dx*0.5
-- Number of elements in height
nHeight = nDiameter * 10
-- Number of elements in length
nLength = nDiameter * 30
-- Number of elements in length + inlet boundary + outlet boundary
nLength_bnd = nLength+2
-- Level required to reach computed dx
level = math.ceil(math.log(nLength_bnd)/math.log(2))
-- Length of the bounding cube
length_bnd = (2^level)*dx
-- Smallest possible element size
dx_eps = length_bnd/2^20
-- Refinement level near cylinder
level_cyl = level + 3

-- cylinder offset from origin
cylinder_x = 6.5*diameter
cylinder_y = height/2.0

---------------------------Seeder configurations -------------------------------
-- Directory to write mesh
folder = 'mesh/'

-- File to write measurement time of mesh generation each step
timing_file = 'sdr_timing.res'

-- How detailed should the output be?
-- The higher the level, the more output you'll get from seeder.
logging = { level = 6 }

-- Debug outputs
NOdebug = {
  debugMode = true,
  debugFiles = true
}

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
  origin = { -dx, -dx, -dx },
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
--! [Seed box]
  {
    attribute = {
      kind = 'seed',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { cylinder_x + diameter, dx_eps, dx_eps },
        vec = {
          { length - cylinder_x - 2*diameter, 0.0, 0.0 },
          { 0.0, height - dx_eps, 0.0 },
          { 0.0, 0.0, width - dx_eps }
        }
      }
    }
  },
--! [Seed box]

--! [Refinement box 1]
  {
    attribute = {
      kind = 'refinement',
      level = level + 1
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { cylinder_x - 5.0*diameter, height/8.0, -dx_eps },
        vec = {
          { length - 5.0*diameter, 0.0, 0.0 },
          { 0.0, height - height/4.0, 0.0 },
          { 0.0, 0.0, width + 2.*dx_eps }
        }
      }
    }
  },
--! [Refinement box 1]

--! [Refinement box 2]
  {
    attribute = {
      kind = 'refinement',
      level = level + 2
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { cylinder_x - 2.0*diameter, height/4.0, -dx_eps },
        vec = {
          { length/2.0 + 2*diameter, 0.0, 0.0 },
          { 0.0, height - height/2.0, 0.0 },
          { 0.0, 0.0, width + 2.*dx_eps }
        }
      }
    }
  },
--! [Refinement box 2]

--! [East boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'east'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { length + dx_eps, -dx_eps, -dx_eps },
        vec = {
          { 0.0, height + 2*dx_eps, 0.0 },
          { 0.0, 0.0, width + 2.*dx_eps }
        }
      }
    }
  },
--! [East boundary]

--! [West boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'west'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -dx_eps, -dx_eps, -dx_eps },
        vec = {
          { 0.0, height + 2*dx_eps, 0.0 },
          { 0.0, 0.0, width + 2.*dx_eps }
        }
      }
    }
  },
--! [West boundary]

--! [North boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'north'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -dx_eps, height + dx_eps, -dx_eps },
        vec = {
          { length + 2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, width + 2.*dx_eps }
        }
      }
    }
  },
--! [North boundary]

--! [South boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'south'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -dx_eps, -dx_eps, -dx_eps },
        vec = {
          { length + 2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, width + 2.*dx_eps }
        }
      }
    }
  },
--! [South boundary]

--! [Front boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'front'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -dx_eps, -dx_eps, width + dx_eps },
        vec = {
          { length + 2*dx_eps, 0.0, 0.0 },
          { 0.0, height + 2.*dx_eps, 0.0 }
        }
      }
    }
  },
--! [Front boundary]

--! [Bottom boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'back'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -dx_eps, -dx_eps, -dx_eps },
        vec = {
          { length + 2*dx_eps, 0.0, 0.0 },
          { 0.0, height + 2.*dx_eps, 0.0 }
        }
      }
    }
  },
--! [Bottom boundary]

--! [Cylinder boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'cylinder',
      calc_dist = true,
      level = level_cyl,
      distance_refine = {
        {
          radius = radius/2.0,
          level_offset = 0,
        },
        {
          radius = radius/2,
          level_offset = -1
        },
      }
    },
    geometry = {
      kind = 'stl',
      object = {
        filename = 'cylinder.stl'
      }
    },
    transformation ={
      deformation = { 1, 1, 20.0},
      translation = { cylinder_x, cylinder_y, -diameter }
    }
  },
--! [Cylinder boundary]
}
------------------- End of seeder configurations -------------------------------
