-- Draw mesh for 2d cylinder with the diameter of 40

require 'args'

folder    = 'mesh/'
comment   = "cylinder2d"

minlevel  = level_usr
print("level of treelm: ", minlevel)

bc_origin = {-2.5, -2.5, -0.5 - 4}
length_bnd = length_usr

bounding_cube = { origin = bc_origin,
                  length = length_bnd }

NOdebug = {debugMode = true, debugFiles = false, debugMesh='debug/' }
spatial_object = {
  {
    attribute   = {
      kind      = 'boundary',
      label     = 'circle',
      level     = minlevel,
      calc_dist = true,
    },
    geometry  = {
      kind    = 'stl',
      object  = {
          filename = stl_path..'cylinder40.stl'
        },
     },
  },
  {
    attribute   = {
      kind      = 'periodic'
    },
    geometry  = {
      kind = 'periodic',
      object  = {
        plane1 = {
          origin = { -1.5, -1.5, 4.25},
          vec = { {length_bnd-10, 0.0, 0.0},
                {0.0, length_bnd-10, 0.0}
          },
        },
        plane2 = {
          origin = { -1.5, -1.5, 5.75},
          vec = { {length_bnd-10, 0.0, 0.0},
                {0.0, length_bnd-10, 0.0}
          },
        },
      }
    },
  },
  {
    attribute = {
      kind    = 'seed',
      label   = 'seed',
    },
    geometry  = {
      kind    = 'canoND',
      object  = { origin = {nelem,nelem,5.} }
    }
  }
}
