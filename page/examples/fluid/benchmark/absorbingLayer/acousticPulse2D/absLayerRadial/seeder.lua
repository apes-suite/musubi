--! [Global variables]
-- Thickness of the absorbing layer
abs_thickness = 0.1 --m

-- Length of the domain
length = 1.0
-- Height of the domain
height = length
-- nElements in length
nLength = 256
-- Element size
dx = length/nLength
-- Level required to reach computed dx
level = math.ceil(math.log(nLength)/math.log(2))
-- Smallest possible element size
dx_eps = length/2.0^20
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
  origin = {-length/2.0, -length/2.0, -length/2.0},
  length = length
}

-- A minimum level, by which all parts in the computational domain should at
-- least be resolved with. Default is 0.
minlevel = level

-- spatial object is composed  by the attribute and geometry
spatial_object = {
  {
   attribute = {
      kind = 'seed',
    },
    geometry = {
      kind = 'canoND',
      object = {
        origin = {0.0, 0.0, dx/2.0}
        } ---object
      } --- geometry
  }, --- attribute
  {
   attribute = {
     kind = 'periodic',
   },               
   geometry = {
     kind = 'periodic',
     object = {
       plane1 = {
         origin = {-length*0.5, -length*0.5, dx+dx_eps/2 },
         vec = {
           {length, 0.0, 0.0},
           {0.0, length, 0.0}}  
         },        
        plane2 = {
          origin = {-length*0.5, -length*0.5, -dx_eps/2},
          vec = {
            {0.0, length, 0.0},
            {length, 0.0, 0.0}
          }  
        }        
      } --- object
    } -- geometry
  } --- attribute
}--- spatial object
