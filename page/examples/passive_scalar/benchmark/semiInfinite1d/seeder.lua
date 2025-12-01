require 'params'

folder    = 'mesh/'
comment   = "semi1d_seeder"

minlevel  = level_usr

bc_origin = {0, -dx*3/2, -dx*3/2}

bounding_cube = { origin = bc_origin,
                  length = length_bnd }

spatial_object = {
  {
    attribute   = {
      kind      = 'periodic'
    },
    geometry  = { 
      kind = 'periodic',
      object  = {
        plane1 = {
          origin = { 0, -dx, -dx},
          vec = { {L, 0, 0.0}, 
                {0.0, 2*dx, 0.0}
          },
        },
        plane2 = {
          origin = { 0, -dx, dx},
          vec = { {L, 0.0, 0.0}, 
                {0.0, 2*dx, 0.0}
          },
        },
      }
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
          origin = { 0, -dx, -dx},
          vec = { {L, 0, 0.0}, 
                {0.0, 0.0, 2*dx}
          },
        },
        plane2 = {
          origin = { 0, dx, -dx},
          vec = { {L, 0.0, 0.0}, 
                {0.0, 0.0, 2*dx}
          },
        },
      }
    },
  },
  {
    attribute   = {
      kind      = 'boundary',
      label     = 'inlet',
      level     = minlevel,
      calc_dist = true,
    },
    geometry  = { 
      kind    = 'canoND',
      object  = {
        origin = { 0, -dx, -dx},
        vec = { {0, 2*dx, 0.0}, 
              {0.0, 0, 2*dx}
        },
        only_surface = true,
      }
    }
  },
  {
    attribute   = {
      kind      = 'boundary',
      label     = 'outlet',
      level     = minlevel,
      calc_dist = true,
    },
    geometry  = { 
      kind    = 'canoND',
      object  = {
        origin = { L, -dx, -dx},
        vec = { {0, 2*dx, 0.0}, 
              {0.0, 0, 2*dx}
        },
        only_surface = true,
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
      object  = { origin = {L/2, 0, 0} }
    }                
  }
}
