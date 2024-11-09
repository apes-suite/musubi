-- choose one of the following
-- xM xP yM yP zM zP
checkDir = 'xP'
level = 6

comment = 'cube refined'
folder = 'mesh/'
label = 'nrbc'
timing_file = 'sdr_timing.res'

minlevel = level
refinement = 0
usePeriodic = true
walls = true
length = 10.
dx= length/2^level
dxDash = dx*0.01
shift = dx-dxDash
origin = {-5., -5., -5. }
size = {length, length, length }
bounding_cube = {origin = origin, 
               length = size[1]}
-- refinebox: three entries: origin, length and refinementlevel
box2 = { 
                 origin = origin,
                   vec = {
                     {length/2.-dxDash, 0., 0. },
                     {0., length/1.-dxDash, 0. },
                     {0., 0., length/1.-dxDash }}
               }

box1 = { 
                 origin = {-1.0, -1.0, -1.4},
                   vec = {
                     {length/10., 0., 0. },
                     {0., length/10., 0. },
                     {0., 0., length/10. }}
               }

spatial_object = {
  {
    -- Defining a domain boundary
    attribute = { kind = 'refinement', 
                  label = 'refine', 
                  level = minlevel + refinement
    },
    geometry = { -- Example for a sphere definition
      kind = 'canoND',
      object = box2
    }},
    {
    -- Defining a seed to identify the part of the computational domain in
    -- the universe cube.
    attribute = { kind = 'seed' },
    geometry = { -- single point definition with a canoND object.
                 kind = 'canoND',
                 object = { origin = {origin[1]+length*0.5, origin[2]+length*0.5, origin[3]+length*0.5} }
               }
  }
}


if checkDir == 'xM' or checkDir == 'xP' 
or checkDir == 'yM' or checkDir == 'yP' then
  table.insert(spatial_object, { 
    attribute = { 
      kind = 'periodic', 
      label = 'periodic'
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = { origin[1], origin[2], dx+shift},
          vec = { { length, 0.0, 0.0},
                  { 0.0, length, 0.0},}
        }, -- plane 1
        plane2 = {
          origin = { origin[1], origin[2], -shift},
          vec = { { 0., length, 0.0},
                  { length, 0., 0.0},}
        }, -- plane 2
      } -- object
    } -- geometry

  })
end

if checkDir == 'xM' or checkDir == 'xP'  then
  table.insert(spatial_object, { 
    attribute = { 
      kind = 'periodic', 
      label = 'periodic'
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = { origin[1], dx+dxDash, origin[3]},
          vec = { { 0., 0.0, length },
                  { length, 0., 0.0},}
        }, -- plane 1
        plane2 = {
          origin = { origin[1], -dxDash  , origin[3]},
          vec = { { length, 0.0, 0.0},
                  { 0.0, 0.0, length},}
        }, -- plane 2
      } -- object
    } -- geometry

  })
end


if checkDir == 'yM' or checkDir == 'yP'  then
  table.insert(spatial_object, { 
    attribute = { 
      kind = 'periodic', 
      label = 'periodic'
    },
    geometry = {
      kind = 'periodic',
      object = {
        plane1 = {
          origin = { dx+dxDash, origin[2], origin[3]},
          vec = { { 0.0, length, 0.0},
                  { 0.0, 0.0, length},}
        }, -- plane 1
        plane2 = {
          origin = { -dxDash  , origin[2], origin[3]},
          vec = { { 0., 0.0, length },
                  { 0., length, 0.0},}
        }, -- plane 2
      } -- object
    } -- geometry

  })
end


table.insert( spatial_object, 
  { -- Channel Walls
    attribute = { kind = 'boundary',
                  label = 'wall_xM', 
                },
    geometry =  {
      kind = 'canoND',
      object = {
      {  -- x-
         origin = {origin[1]+shift, origin[2], origin[3]},
         vec = {{ 0., length,0.0},
              {0.0,0., size[3]}}
      }}}})
table.insert( spatial_object, 
  { -- Channel Walls
    attribute = { kind = 'boundary',
                  label = 'wall_xP', 
                },
    geometry =  {
      kind = 'canoND',
      object = {      {  -- x+
         origin = {origin[1]+size[1]-shift, origin[2], origin[3]},
         vec = {{ 0., length,0.0},
              {0.0,0., size[3]}}
      }}}})
--table.insert( spatial_object, 
--  { -- Channel Walls
--    attribute = { kind = 'boundary',
--                  label = 'wall_yM', 
--                },
--    geometry =  {
--      kind = 'canoND',
--      object = {      {  -- y-
--         origin = {origin[1], origin[2]+shift, origin[3]},
--         vec = {{ length,0.,0.0},
--              {0.0,0., size[3]}}
--      }}}})
--table.insert( spatial_object, 
--  { -- Channel Walls
--    attribute = { kind = 'boundary',
--                  label = 'wall_yP', 
--                },
--    geometry =  {
--      kind = 'canoND',
--      object = {      {  -- y+
--         origin = {origin[1], origin[2]+size[2]-shift, origin[3]},
--         vec = {{ length,0.,0.0},
--              {0.0,0., size[3]}}
--      }}}})
--table.insert( spatial_object, 
--  { -- Channel Walls
--    attribute = { kind = 'boundary',
--                  label = 'wall_zM', 
--                },
--    geometry =  {
--      kind = 'canoND',
--      object = {      {  -- z-
--         origin = {origin[1], origin[2], origin[3]+shift},
--         vec = {{ length,0.,0.0},
--              {0.0,size[2],0. }}
--      }}}})
--table.insert( spatial_object, 
--  { -- Channel Walls
--    attribute = { kind = 'boundary',
--                  label = 'wall_zP', 
--                },
--    geometry =  {
--      kind = 'canoND',
--      object = {      {  -- z+
--         origin = {origin[1], origin[2], origin[3]+size[3]-shift},
--         vec = {{ length,0.,0.0},
--              {0.0, size[2],0.}}
--      },
--    }}})

