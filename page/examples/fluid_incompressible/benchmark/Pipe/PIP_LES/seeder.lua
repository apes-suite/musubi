----------------------- PLEASE READ THIS ---------------------------!!!
-- This input file is set up to run for regression check
-- Please make sure you DO NOT MODIFY AND PUSH it to the repository
--------------------------------------------------------------------!!!

--! [Global variables]
shepherd = true
if shepherd then
  stlFolder = '$!stl_path!$'
  -- Number of elements in diameter
  nDiameter = 16
  -- diameter to length ratio
  d_l_ratio = 1.0
else
  stlFolder = './../PIP_Simple/'
  -- Number of elements in diameter
  nDiameter = 64
  -- diameter to length ratio
  d_l_ratio = 5.0
end

-- Diameter of the pipe [m]
diameter =  0.41
-- Radius of the pipe
radius = diameter/2.0
-- length of the pipe [m]
length = d_l_ratio * diameter
-- Geometry information of pipe.stl
pipe_origin = { -0.1, -0.205, -0.205 }
pipe_length = { length, 0.41, 0.41 }
pipe_center = { length/2.0, 0.0, 0.0 }
-- Number of elements in length
nLength = d_l_ratio * nDiameter
-- Element size
dx = diameter/nDiameter
-- Number of elements in length + 1 element for either side of boundary
nLength_bnd = nLength + 2
-- Level required to reach defined dx
level = math.ceil(math.log(nLength_bnd)/math.log(2))
-- Level to refine the pipe
if shepherd then
  level_pipe = level + 1
else
  level_pipe = level + 2
end
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
logging = { level = 10 }

-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0.
minlevel = level

-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
  origin = { -dx, pipe_origin[2]-dx, pipe_origin[3]-dx },
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
      object = {
        origin = {dx/2.0, pipe_center[2], pipe_center[3]},
        vec = {length-2*dx, 0.0, 0.0}
      }
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
      level = level_pipe,
      distance_refine = {
        {
          radius = 8 * dx * 2^(minlevel-(level+2)),
          level_offset = 0
        },
        {
          radius = 8 * dx * 2^(minlevel-(level+1)),
          level_offset = -1
        }
      }
    },
    geometry = {
      kind = 'stl',
      object = {
        filename = stlFolder..'pipe.stl',
        format = 'binary'
      }
    }
  },
--! [Pipe STL]

--! [West-East periodic boundary]
  {
    attribute = {
      kind = 'periodic',
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = { -dx_eps, pipe_origin[2] - dx, pipe_origin[3] - dx },
          vec = {
            { 0.0, 0.0, diameter + 2.*dx },
            { 0.0, diameter + 2.*dx, 0.0 }
          }
        },
        plane2 = { 
          origin = {
            length + dx_eps,
            pipe_origin[2] - dx,
            pipe_origin[3] - dx
          },
          vec = {
            { 0.0, diameter + 2.*dx, 0.0 },
            { 0.0, 0.0, diameter + 2.*dx }
          }
        }
      }
    }
  },
--! [West-East periodic boundary]
}
------------------- End of seeder configurations -------------------------------
