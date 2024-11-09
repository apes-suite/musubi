--! [Global variables]
-- actual height of the channel
height =  0.41 --m
-- Number of elements in channel height
nHeight = 40
-- logical to generate a symmetric testcase on north boundary
sym_case = true

-- seeder only use half channel height if we use symmetry condition
-- addapt number of elements and channel length acordingly
if ( sym_case == true ) then
  -- Number of elements in height_sdr
  height_sdr =  height * 0.5
  nHeight = math.floor(nHeight * 0.5)
  -- Length to height ratio
  l_h = 14.0
else
  -- Number of elements in height_sdr
  height_sdr =  height
  nHeight = nHeight
  -- Length to height ratio
  l_h = 7.0
end

-- Number of elements in length
nLength = nHeight*l_h 
-- Element size
dx = height_sdr/nHeight
-- Half of element size
dx_half = dx*0.5
-- Length of the channel
length = nLength*dx
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
          length/2.0 + dx_eps + x_offset, 
          height_sdr*0.5 + dx_eps + y_offset, 
          dx_half 
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
        origin = { -x_offset-dx_eps, -y_offset-dx_eps, -dx/2.0 },
        vec = {
          { 0.0, height_sdr + 2*y_offset + 2*dx_eps, 0.0 },
          { 0.0, 0.0, 2*dx } 
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
        origin = { length+x_offset+dx_eps, -y_offset-dx_eps, -dx/2.0 },
        vec = { 
          { 0.0, height_sdr + 2*y_offset + 2*dx_eps, 0.0 },
          { 0.0, 0.0, 2*dx }
        }
      }
    }
  },
--! [East boundary]

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
        origin = { -x_offset-dx_eps, -y_offset-dx_eps, -dx/2.0 },
        vec = {
          { length + 2*x_offset + 2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, 2*dx }
        }
      }
    }
  },
--! [South boundary]

--! [North boundary]
  {
    attribute = {
      kind = 'boundary',  
      calc_dist = false, -- to calculate q-values
      label = 'north',   -- label for the boundary at y = height_sdr + y_offset
    },
    geometry = {
      kind = 'canoND', 
      object = {
        origin = { -x_offset-dx_eps, height_sdr+y_offset+dx_eps, -dx/2 },
        vec = {
          { length + 2*x_offset + 2*dx_eps, 0.0, 0.0 },
          { 0.0, 0.0, 2*dx } 
        }
      }
    }
  },
--! [North boundary]

--! [Top and bottom boundary (periodic)]
  -- No label required for periodic boundary.
  -- Normal direction of each plane should point outward from fluid domain.
  {
    attribute = {
      kind = 'periodic', 
    },
    geometry = {
      kind = 'periodic',
      object = {
        -- Top boundary
        plane1 = {
          origin = { -x_offset-dx_eps, -y_offset-dx_eps, dx+dx_eps },
          vec = {
            { length + 2*x_offset + 2*dx_eps, 0.0, 0.0 },
            { 0.0, height_sdr + 2*y_offset + 2*dx_eps, 0.0} 
          }
        }, 
        -- Bottom boundary
        plane2 = {
          origin = { -x_offset-dx_eps, -y_offset-dx_eps, -dx_eps },
          vec = {
            { 0.0, height_sdr + 2*y_offset + 2*dx_eps, 0.0 },
            { length + 2*x_offset + 2*dx_eps, 0.0, 0.0}
          }
        }
      }  
    }
  },
--! [Top and bottom boundary (periodic)]
} -- spatial object
------------------- End of seeder configurations -------------------------------

