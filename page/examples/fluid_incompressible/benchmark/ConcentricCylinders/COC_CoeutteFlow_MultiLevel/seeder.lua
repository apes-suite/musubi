--! [Global variables]
shepherd = true
if shepherd then
  stlFolder = '$!stl_path!$'
else
  stlFolder = '../'
end
-- Ratio between inner and outer cylinder
ratio = 0.25
-- Diameter of the outer cylinder [m]
dia_outer = 1.0
-- Radius of the outer cylinder [m]
rad_outer = dia_outer/2.0
-- Diameter of the inner cylinder [m]
dia_inner = dia_outer * ratio
-- Radius of the inner cylinder [m]
rad_inner = dia_inner/2.0
-- Number of elements in inner diameter
nDia_inner = 16
-- Element size
dx = dia_inner/nDia_inner
-- Half of element size
dx_half = dx*0.5
-- Number of elements in outer diameter
nDia_outer = math.ceil(dia_outer/dx)
-- Number of elements in bounding cube
-- = number of elements in outer diameter + 2
nLength_bnd = nDia_outer + 2
-- Level required to reach computed dx
level = math.ceil(math.log(nLength_bnd)/math.log(2))
-- Length of the bounding cube
length_bnd = (2^level)*dx
-- Smallest possible element size
dx_eps = length_bnd/2^20
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

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
  origin = { -length_bnd/2.0,  -length_bnd/2.0,  -length_bnd/2.0 - dx/2.0 },
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
-- - periodic
--
-- Periodic boundaries are special spatial objects.
-- Attribute with kind = 'periodic' can only make use of geometric objects of
-- kind = 'periodic'.
spatial_object = {
--! [Seed point]
  {
    attribute = {
      kind = 'seed',   ----seed
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { rad_inner+dx, rad_inner+dx, 0.0 },
      }
    }
  },
--! [Seed point]

--! [Inner cylinder]
  {
    attribute = {
      kind = 'boundary',
      label = 'inner',
      calc_dist = true,
      flood_diagonal = true,
      level = level+1,
      distance_refine = {
        {
          radius = 3*dx,
          level_offset = 0
        },
        {
          radius = 6*dx,
          level_offset = -1
        }
      }
    },
    geometry = {
      kind = 'stl',
      object = {
        filename = stlFolder..'cylinder_Dia1m.stl'
      }
    },
    transformation = {
      -- Deform cylinder stil to ratio
      deformation = { ratio, ratio, 1.0 }
    }
  },
--! [Inner cylinder]

--! [Outer cylinder]
  {
    attribute = {
      kind = 'boundary',
      label = 'outer',
      calc_dist = true,
      flood_diagonal = true,
      level = level+1,
      distance_refine = {
        {
          radius = 3*dx,
          level_offset = 0
        },
        {
          radius = 6*dx,
          level_offset = -1
        }
      }
    },
    geometry = {
      kind = 'stl',
      object = {
        filename = stlFolder..'cylinder_Dia1m.stl'
      }
    }
  },
--! [Outer cylinder]

--! [Top and bottom boundary (periodic)]
  {
    attribute = {
      kind = 'periodic',
      label = 'periodic',
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = { -length_bnd/2.0, -length_bnd/2.0, -dx/2.0-dx_eps },
          vec = {
            { 0.0, length_bnd, 0.0 },
            { length_bnd, 0.0, 0.0 },
          }
        },
        plane2 = {
          origin = { -length_bnd/2.0, -length_bnd/2.0, dx/2.0+dx_eps },
          vec = {
            { length_bnd, 0.0, 0.0 },
            { 0.0, length_bnd, 0.0 }
          }
        }
      }
    }
  }
--! [Top and bottom boundary (periodic)]
}
---------------------------Seeder configurations -------------------------------

