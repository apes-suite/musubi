--! [Global variables]
-- Length of the nozzle [m]
l_nozzle = 20e-2
-- Outer diameter of nozzle inlet [m]
outer_dia_nozzle = 2e-2
-- Outer radius of nozzle inlet [m]
outer_rad_nozzle = outer_dia_nozzle*0.5
-- Ratio of inner dia to outer dia
inner_to_outer_ratio = 0.4
-- Inner diameter of nozzle [m]
inner_dia_nozzle = outer_dia_nozzle*inner_to_outer_ratio
-- Inner radius of nozzle [m]
inner_rad_nozzle = inner_dia_nozzle*0.5
-- Origin of nozzle
origin_nozzle = {0.0,0.0,0.0}
-- Position of inner diameter area in nozzle
l_neck = -l_nozzle/3.0
nozzle_inner_dia_X = l_neck
-- Distance of inlet to nozzle center
inlet_2_nozzleCenter = l_nozzle/2.0
-- Distance of outlet to nozzle center
outlet_2_nozzleCenter = l_nozzle/2.0
-- Height of free flow area
h_ch = outer_dia_nozzle
-- Length of free flow area
l_ch = inlet_2_nozzleCenter + outlet_2_nozzleCenter
-- Number of elements in height of free flow area
nElems_h_ch = 64
-- Element size
dx = h_ch/nElems_h_ch
-- Number of elements in length
nLength = math.ceil(l_ch/dx)
-- Number of elements in bounding cube length
-- +2 for inlet and outlet plane
nLength_bnd = nLength+2
-- Refinement level required to achieve nElems_h_ch
level = math.ceil(math.log(nLength_bnd)/math.log(2))
-- Physical length of bounding cube
length_bnd = (2^level)*dx
-- Refinement level near nozzle
nozzleLevel = level + 0
-- Element size near nozzle
dx_nozzle = length_bnd/2^nozzleLevel
-- Smallest element size
dx_eps = length_bnd/2^(20)
-- Half of coarsest and finest element size
dx_c_half = dx/2.0
dx_f_half = dx_nozzle/2.0
z_pos = dx/2.0
-- nozzle inlet
inlet_nozzle = {-l_nozzle/2.0, -outer_rad_nozzle, z_pos}
-- nozzle outlet
outlet_nozzle = {l_nozzle/2.0, -outer_rad_nozzle, z_pos}
-- nozzle center
center_nozzle_l = {nozzle_inner_dia_X, -inner_rad_nozzle, z_pos}
center_nozzle_u = {nozzle_inner_dia_X, inner_rad_nozzle, z_pos}

-- Origin of inlet BC
origin_inletBC = {
  -inlet_2_nozzleCenter + dx_c_half,
  -h_ch/2.0-dx_c_half,
  -h_ch/2.0-dx_c_half
}

-- Origin of outlet BC
origin_outletBC = {
  outlet_2_nozzleCenter - dx_c_half,
  -h_ch/2.0- dx_c_half,
  -h_ch/2.0 - dx_c_half
}

-- z position of nozzle
nozzle_zPos = outer_rad_nozzle
--! [Global variables]

---------------------------Seeder configurations -------------------------------
simulation_name = 'nozzle'
-- Directory to write mesh
folder = 'mesh/'

-- file to write measurement time of mesh generation each step
timing_file = 'sdr_timing.res'

-- How detailed should the output be?
-- The higher the level, the more output you'll get from seeder.
logging = { level = 5 }

-- Debug outputs
NOdebug = {
  debugMode = true,
  debugFiles = true
}

-- bounding cube: two entries: origin and length in this
-- order, if no keys are used
bounding_cube = {
  origin = {
    -inlet_2_nozzleCenter - dx,
    -h_ch/2.0-dx,
    -dx-h_ch/2.0
  },
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
-- - triangle
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
        origin = { origin_inletBC[1]+dx, 0.0, dx_c_half },
        vec = { l_ch-2*dx, 0.0,0.0}
      }
    }
  },
--! [Seed point]


--! [Nozzle boundary]
  {
    -- Defining a domain boundary
    attribute = {
      kind = 'boundary', -- or seed, refinement
      label = 'nozzle',   -- some label to identify the boundary condition
      level = nozzleLevel,
      calc_dist = true,
      flood_diagonal = false,
      distance_refine = {
        radius = inner_rad_nozzle/2,
        level_offset = 0
      },
    },
    geometry = {
      kind = 'canoND',
      object = {
        -- Upper triangle
        -- plane 1
        {
          origin = { -l_nozzle/2.0, outer_rad_nozzle, -nozzle_zPos },
          vec = {
            { l_nozzle/2.0-math.abs(nozzle_inner_dia_X),
              -(outer_rad_nozzle-inner_rad_nozzle),
              0.0
            },
            { 0.0, 0.0, 2*nozzle_zPos }
          }
        },
        -- plane 2
        {
          origin = { -l_nozzle/2.0, outer_rad_nozzle, -nozzle_zPos },
          vec = {
            { l_nozzle, (outer_rad_nozzle-outer_rad_nozzle), 0.0 },
            { 0.0, 0.0, 2*nozzle_zPos }
          }
        },
        -- plane 3
        {
          origin = { nozzle_inner_dia_X, inner_rad_nozzle, -nozzle_zPos },
          vec = {
            {
              l_nozzle/2.0+math.abs(nozzle_inner_dia_X),
              outer_rad_nozzle-inner_rad_nozzle,
              0.0
            },
            { 0.0, 0.0, 2*nozzle_zPos }
          }
        },
        -- Lower triangle
        -- plane 1
        {
          origin = { -l_nozzle/2.0, -outer_rad_nozzle, -nozzle_zPos },
          vec = {
            { l_nozzle/2.0-math.abs(nozzle_inner_dia_X),
              outer_rad_nozzle-inner_rad_nozzle,
              0.0
            },
            { 0.0, 0.0, 2*nozzle_zPos }
          }
        },
        -- plane 2
        {
          origin = { -l_nozzle/2.0, -outer_rad_nozzle, -nozzle_zPos },
          vec = {
            {
              l_nozzle,
              -(outer_rad_nozzle-outer_rad_nozzle),
              0.0
            },
            { 0.0, 0.0, 2*nozzle_zPos }
          }
        },
        -- plane 3
        {
          origin = { nozzle_inner_dia_X, -inner_rad_nozzle, -nozzle_zPos },
          vec = {
            { l_nozzle/2.0+math.abs(nozzle_inner_dia_X),
              -(outer_rad_nozzle-inner_rad_nozzle),
              0.0
            },
            { 0.0, 0.0, 2*nozzle_zPos }
          }
        }
      }
    }
  },
--! [Nozzle boundary]

--! [West boundary]
  {
    attribute = {
      kind = 'boundary',
      label = 'west',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = origin_inletBC,
        vec = {
          { 0.0, h_ch+2*dx, 0.0 },
          { 0.0, 0.0, h_ch+2*dx }
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
        {
          origin = origin_outletBC,
          vec = {
            { 0.0, h_ch+2*dx, 0.0 },
            { 0.0, 0.0, h_ch+2*dx }
          }
        },
      }
    }
  },
--! [East boundary]


--! [Top and bottom boundary (periodic)]
  {
    -- Defining a periodic
    attribute = {
      kind = 'periodic',
      level = level
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = {origin_inletBC[1], -h_ch/2.0-dx_c_half, dx+dx_eps },
          vec = { {l_ch+2*dx, 0.0, 0.0},
                  {0.0, h_ch+2*dx, 0.0}}
        },
        plane2 = {
          origin = {origin_inletBC[1], -h_ch/2.0-dx_c_half, -dx_eps },
          vec = { {0.0, h_ch+2*dx, 0.0},
                  {l_ch+2*dx, 0.0, 0.0}}
        }
      }
    }
  },
--! [Top and bottom boundary (periodic)]
}
