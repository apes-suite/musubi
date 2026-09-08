require 'params'

folder    = 'mesh/'
comment   = "BrinkmanChannel2D_seeder"

minlevel  = level

bounding_cube = { origin = bc_origin,
                  length = length_bnd }

spatial_object = {
  {
    attribute   = {
      kind      = 'boundary',
      label     = 'bottom',
      level     = minlevel,
      calc_dist = true,
    },
    geometry  = { 
      kind    = 'canoND',
      object  = {
        origin = { dx/2-dx_eps, dx/2-dx_eps, -dx-dx_eps},
        vec = { {domainlen-dx+2*dx_eps, 0.0, 0.0}, 
                {0.0, 0.0, 2*dx+2*dx_eps}
        },
        only_surface = true,
      }
    }
  },
  {
    attribute   = {
      kind      = 'boundary',
      label     = 'top',
      level     = minlevel,
      calc_dist = true,
    },
    geometry  = {  
      kind    = 'canoND',
      object  = {
        origin = { dx/2-dx_eps, Dia - dx/2 + dx_eps, -dx-dx_eps},
        vec = { {domainlen-dx+2*dx_eps, 0.0, 0.0}, 
                {0.0, 0.0, 2*dx+2*dx_eps}
        },
        only_surface = true,     
      }
    }
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
        origin = { dx/2-dx_eps, dx/2-dx_eps, -dx-dx_eps},
        vec = { {0.0, Dia-dx+2*dx_eps, 0.0}, 
                {0.0, 0.0, 2*dx+2*dx_eps}
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
        origin = { domainlen-dx/2+dx_eps, dx/2-dx_eps, -dx-dx_eps},
        vec = { {0.0, Dia-dx+2*dx_eps, 0.0}, 
                {0.0, 0.0, 2*dx+2*dx_eps}
        },
        only_surface = true,
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
            origin = { dx/2-dx_eps, dx/2-dx_eps, -dx-dx_eps},
            vec = { {domainlen-dx+2*dx_eps, 0.0, 0.0}, 
                  {0.0, Dia-dx+2*dx_eps, 0.0}
          },
        },
        plane2 = {
          origin = { dx/2-dx_eps, dx/2-dx_eps, dx+dx_eps},
          vec = { {0.0, Dia-dx+2*dx_eps, 0.0}, 
                {domainlen-dx+2*dx_eps, 0.0, 0.0}
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
      object  = { origin = seed_orig }
    }                
  }
}
