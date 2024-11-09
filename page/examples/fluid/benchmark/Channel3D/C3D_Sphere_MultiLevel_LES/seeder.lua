-- Path to stl file
shepherd = true
if shepherd then
  stlFolder = '$!stl_path!$'
else
  stlFolder = './'
end
--! [Global variables]
-- Geometry information of sphere.stl
-- Center of sphere
sphere_center = {0.0, 0.0, 0.0}
-- Diameter of sphere [m]
diameter = 0.1
-- Radius of sphere [m]
radius = diameter/2.0
-- Height to diameter ratop
h_d = 10.0
-- Height of the channel [m]
height = h_d*diameter
-- Width of the channel [m]
width = height
-- length to height ratio
l_h = 3.0
-- Length of the channel [m]
length = l_h*height
-- Distance between sphere and west boundary
westbnd_to_sphere = 6.5*diameter
-- Distance between sphere and east boundary
sphere_to_eastbnd = length - westbnd_to_sphere
-- Number of elements in diameter
nDiameter = 5
-- Number of elements in height
nHeight = h_d*nDiameter
-- Number of elements in width
nWidth = nHeight
-- Number of elements in length
nLength = l_h*nHeight
-- Element size
dx = diameter/nDiameter 
-- Half of element size
dx_half = dx*0.5
-- Number of elements in bounding cube
-- = number of elements in channel length + west bnd + east bnd
nLength_bnd = nLength + 2
-- Level required to reach computed dx
level = math.ceil(math.log(nLength_bnd)/math.log(2))
-- Length of the bounding cube
length_bnd = (2^level)*dx
-- Refinement level for sphere boundary
sphere_level = level + 3
dx_sphere = length_bnd/2^sphere_level
-- Refinement radius
refine_radius1 = 0.1*diameter
refine_radius2 = 0.4*diameter
refine_radius3 = refine_radius2 + diameter
-- Smallest possible element size
dx_eps = length_bnd/2^20
-- Origin of west and east boundary
origin_west = {-westbnd_to_sphere-dx_eps, -height/2.0-dx_eps, -width/2.0-dx_eps}
origin_east = {sphere_to_eastbnd+dx_eps, -height/2.0-dx_eps, -width/2.0-dx_eps}
--! [Global variables]

---------------------------Seeder configurations -------------------------------
-- Directory to write mesh
folder = 'mesh/'

-- file to write measurement time of mesh generation each step
timing_file = 'sdr_timing.res'

-- How detailed should the output be?
-- The higher the level, the more output you'll get from seeder.
logging = { level = 3 }

-- Debug outputs
NOdebug = {
  debugMode = true,
  debugFiles = true
}

-- Bounding_cube: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
  origin = { -westbnd_to_sphere-dx, -height/2.0-dx, -width/2.0 - dx },
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
-- - ellipsoid
spatial_object = {
--! [Seed point]
  {
    attribute = {
      kind = 'seed'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -westbnd_to_sphere+2*dx, -height/2.0+2*dx, -width/2.0+2*dx },
        vec = {
          { length-4*dx, 0.0, 0.0 },
          { 0.0, height-4*dx, 0.0 }
        }
      }
    }
  },
--! [Seed point]

--! [West boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'west',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = origin_west,
        vec = {
          { 0.0, height+2*dx_eps, 0.0 },
          { 0.0, 0.0, width+2*dx_eps }
        }
      }
    }
  },
--! [West boundary]

--! [East boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'east',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = origin_east,
        vec = {
          { 0.0, height+2*dx_eps, 0.0 },
          { 0.0, 0.0, width+2*dx_eps }
        },
      },
    },
  },
--! [East boundary]

--! [South boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'south',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = origin_west,
        vec = {
          { length+2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, width+2*dx_eps }
        }
      }
    }
  },
--! [South boundary]

--! [North boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'north',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { origin_west[1], height/2.0+dx_eps, origin_west[3] },
        vec = {
          { length+2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, width+2*dx_eps }
        }
      }
    }
  },
--! [North boundary]

--! [Bottom boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'bottom',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = origin_west,
        vec = {
          { length+2*dx_eps, 0.0, 0.0 },
          { 0.0, height+2*dx_eps, 0.0 }
        }
      }
    }
  },
--! [Bottom boundary]

--! [Top boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'top',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { origin_west[1], origin_west[2], width/2.0+dx_eps },
        vec = {
          { length+2*dx_eps, 0.0, 0.0 },
          { 0.0, height+2*dx_eps, 0.0 }
        }
      }
    }
  },
--! [Top boundary]


--! [Sphere boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'sphere',
      calc_dist = true,
      flood_diagonal = true,
      level = sphere_level,
      distance_refine = {
        {
          level_offset =  0,
          radius = refine_radius1
        },
        {
          level_offset = -1,
          radius = refine_radius2
        },
        {
          level_offset = -2,
          radius = refine_radius3
        },
      }
    },
    geometry = {
      kind = 'stl',
      format = 'binary',
      object = {
        filename = stlFolder..'sphere.stl',
      }
    },
    NOtransformation = {
      deformation = 2.0
    --   deformation =  {          R,          R, R },
    --   translation =  { cylinder_x, cylinder_y, 0 },
    }
  },
--! [Sphere boundary]
--! [Refinement1 ellipsoid]
  {
    attribute = {
      kind = 'refinement',
      level = sphere_level,
    },
    geometry = {
      kind = 'ellipsoid',
      object = {
        origin = { 0.0, 0.0, 0.0 },
        radius = { diameter, radius+refine_radius1, radius+refine_radius1 }
      }
    },
    transformation = {
      translation = { radius, 0.0, 0.0 }
    }
  },
--! [Refinement1 ellipsoid]
--! [Refinement2 ellipsoid]
  {
    attribute = {
      kind = 'refinement',
      level = sphere_level-1,
    },
    geometry = {
      kind = 'ellipsoid',
      object = {
        origin = { 0.0, 0.0, 0.0 },
        radius = {
          diameter + radius,
          radius + refine_radius2,
          radius + refine_radius2
        }
      }
    },
    transformation = {
      translation = { 2.0*radius, 0.0, 0.0 }
    }
  },
--! [Refinement2 ellipsoid]
--! [Refinement3 ellipsoid]
  {
    attribute = {
      kind = 'refinement',
      level = sphere_level-2,
    },
    geometry = {
      kind = 'ellipsoid',
      object = {
        origin = { 0.0, 0.0, 0.0 },
        radius = {
          diameter + 6*radius,
          radius + refine_radius3,
          radius + refine_radius3
        }
      }
    },
    transformation = {
      translation = { 3*radius, 0.0, 0.0 }
    }
  },
--! [Refinement3 ellipsoid]
--! [Refinement4 box]
  {
    attribute = {
      kind = 'refinement',
      level = minlevel+1,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { sphere_center[1]-5*diameter, -height/4.0, -width/4.0 },
        vec = {
          { length-4*diameter, 0.0, 0.0 },
          { 0.0, height-height/2.0, 0.0 },
          { 0.0, 0.0, width-width/2.0 }
        }
      }
    },
  },
--! [Refinement4 box]
}
------------------- End of seeder configurations -------------------------------
