----------------------- PLEASE READ THIS ---------------------------!!!
-- This input file is set up to run for regression check
-- Please make sure you DO NOT MODIFY AND PUSH it to the repository
--------------------------------------------------------------------!!!

--! [Global variables]
-- Length of the square cavity [m]
length = 1.0
-- Number of elements in height
nLength = 128
-- Element size
dx = length/nLength
-- Width of the cavity
width = dx
-- Number of elements in length + 1 element for either side of boundary
nLength_bnd = nLength + 2
-- Level required to reach defined dx
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

-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0.
minlevel = level

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
  origin = {-dx, -dx, -dx},
  length = length_bnd
}

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
      kind = 'seed'
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {length/2.0, length/2.0, width/2.}
      }
    }
  },
--! [Seed point]

--! [North boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'north',
      calc_dist = true
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = { -dx_eps, length + dx_eps, -dx_eps },
        vec = {
          { length + 2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, width + 2*dx_eps }
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
          { 0.0, 0.0, width + 2*dx_eps }
        }
      }
    }
  },
--! [South boundary]

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
          { 0.0, length + 2*dx_eps, 0.0 },
          {0.0, 0.0, width + 2*dx_eps }
        }
      }
    }
  },
--! [West boundary]

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
          { 0.0, length + 2*dx_eps, 0.0 },
          { 0.0, 0.0, width + 2*dx_eps }
        }
      }
    }
  },
--! [East boundary]

--! [Top and bottom boundary (periodic)]
  -- No label required for periodic boundary.
  -- Normal direction of each plane should point outward from fluid domain.
  {
    attribute = {
      kind = 'periodic'
    },
    geometry = {
      kind = 'periodic',
      object = {
        -- Top boundary
        plane1 = {
          origin = { -dx, -dx, dx+dx/2.0 },
          vec = {
            { length + 2*dx, 0.0, 0.0 },
            { 0.0, length + 2*dx, 0.0 }
          }
        },
        -- Bottom boundary
        plane2 = {
          origin = { -dx, -dx, -dx/2.0 },
          vec = {
            { 0.0, length + 2*dx, 0.0 },
            { length + 2*dx, 0.0, 0.0 }
          }
        }
      }
    }
  },
--! [Top and bottom boundary (periodic)]
}
------------------- End of seeder configurations -------------------------------
