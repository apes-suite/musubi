----------------------- PLEASE READ THIS ---------------------------!!!
-- This input file is set up to run for regression check
-- Please make sure you DO NOT MODIFY AND PUSH it to the repository
--------------------------------------------------------------------!!!

--! [Global variables]
-- Path to stl file
recheck = true
if recheck then
  stlFolder = '$!stl_path!$'
else
  stlFolder = './'
end
-- Geometry information of pipeSplit.stl
pipe_origin = { -6.10638, -4.0, -1.0 }
pipe_length = { 12.1064, 8.0, 2.0 }
pipe_center = { 0.0, 0.0, 0.0 }
-- Information for inlet and outlet boundaries
pipe_west_origin = { -5.0, -1.0, -1.0 }
pipe_northeast_origin = { 3.0, 2.0, -1.0 }
pipe_southeast_origin = { 5.0, -4.0, -1.0 }
-- Diameter of the pipe [m]
diameter =  2.0
-- Radius of the pipe
radius = diameter/2.0
-- Number of elements in diameter
nDiameter = 50
-- Element size
dx = diameter/nDiameter
-- length of the simulation domain
length = pipe_southeast_origin[1] - pipe_west_origin[1]
-- Number of elements in length
nLength = math.ceil(length / dx)
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
  origin = { pipe_west_origin[1]-dx, pipe_origin[2]-dx, pipe_origin[3]-dx },
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
spatial_object = {
--! [Seed point]
  {
    attribute = {
      kind = 'seed',
    },
    geometry = {
      kind = 'canoND',
      object = { origin = pipe_center }
    }
  },
--! [Seed point]

--! [Pipe STL]
  {
    attribute = {
      kind = 'boundary',
      label = 'pipe',
      calc_dist = true, -- calculate qValues to increase boundary approximation
      flood_diagonal = false,
      level = level,
    },
    geometry = {
      kind = 'stl',
      object = {
        filename = stlFolder..'pipeSplit.stl',
        format = 'binary'
      }
    }
  },
--! [Pipe STL]

--! [West boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'west',
      calc_dist = false, -- calculate qValues to increase boundary approximation
      flood_diagonal = false,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {
          pipe_west_origin[1] - dx_eps,
          pipe_west_origin[2] - dx,
          pipe_west_origin[3] - dx },
        vec = {
          { 0.0, diameter + 2.*dx, 0.0 },
          { 0.0, 0.0, diameter + 2.*dx }
        }
      }
    }
  },
--! [West boundary]

--! [North-East boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'north_east',
      calc_dist = false, -- calculate qValues to increase boundary approximation
      flood_diagonal = false,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {
          pipe_northeast_origin[1] + dx_eps,
          pipe_northeast_origin[2] - dx,
          pipe_northeast_origin[3] - dx
        },
        vec = {
          { 0.0, diameter + 2.*dx, 0.0 },
          { 0.0, 0.0, diameter + 2.*dx }
        }
      }
    }
  },
--! [North-East boundary]

--! [South-East boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'south_east',
      calc_dist = false, -- calculate qValues to increase boundary approximation
      flood_diagonal = false,
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {
          pipe_southeast_origin[1] + dx_eps,
          pipe_southeast_origin[2] - dx,
          pipe_southeast_origin[3] - dx
        },
        vec = {
          { 0.0, diameter + 2.*dx, 0.0 },
          { 0.0, 0.0, diameter + 2.*dx }
        }
      }
    }
  },
--! [North-East boundary]

}
------------------- End of seeder configurations -------------------------------
