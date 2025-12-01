-- draw mesh for a 2d cube

require 'args'

folder    = 'mesh/'
comment   = "gaussian_pulse"

minlevel  = level_usr

bc_origin = {-nelem-2.5, -nelem-2.5, -0.5 - 4}
length_bnd = length_usr

bounding_cube = { origin = bc_origin,
                  length = length_bnd }

spatial_object = {
  {
    attribute   = {
      kind      = 'periodic',
    },
    geometry  = {
      kind    = 'periodic',
      object  = {
        plane1 = {
            origin = { -nelem-1, -nelem-1, -1},
            vec = { {0., 0.0, 2.0},
                  {0.0, 2.0*(nelem+1), 0.0}
          },
        },
        plane2 = {
          origin = { nelem+1, -nelem-1, -1},
          vec = { {0.0, 0., 2.0},
                {0., 2*(nelem+1), 0.0}
        },
      },
      }
    }
  },
  {
    attribute   = {
      kind      = 'periodic',
    },
    geometry  = {
      kind    = 'periodic',
      object  = {
        plane1 = {
            origin = { -nelem-1, -nelem-1, -1},
            vec = { {2*(nelem+1), 0.0, 0.0},
                  {0.0, 0., 2.}
          },
        },
        plane2 = {
          origin = { -nelem-1, nelem+1, -1},
          vec = { {2*(nelem+1), 0., 0.0},
                {0., 0.0, 2.}
        },
      },
      }
    }
  },
  {
    attribute   = {
      kind      = 'periodic',
    },
    geometry  = {
      kind    = 'periodic',
      object  = {
        plane1 = {
            origin = { -nelem-1, -nelem-1, -1},
            vec = { {2*(nelem+1), 0.0, 0.0},
                  {0.0, 2.0*(nelem+1), 0.0}
          },
        },
        plane2 = {
          origin = { -nelem-1, -nelem-1, 1},
          vec = { {2*(nelem+1), 0., 0.0},
                {0., 2*(nelem+1), 0.0}
        },
      },
      }
    }
  },
  {
    attribute = {
      kind    = 'seed',
      label   = 'seed',
    },
    geometry  = {
      kind    = 'canoND',
      object  = { origin = {0.,0.,0.} }
    }
  }
}
