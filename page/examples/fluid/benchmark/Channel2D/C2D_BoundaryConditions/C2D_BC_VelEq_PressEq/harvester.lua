-- Example script showing the harvesting configuration for Seeder mesh files --
--
-- Use this as input to sdr_harvesting to obtain visualization files.
-- -----------------------------------------------------------------------------

-- Mesh to visualize, this is the common mesh definition, as found in
-- treelm/config.lua.
-- Either a predefined mesh, or the prefix for a mesh on disk.
mesh = 'mesh/'
--mesh = { predefined = 'cube',
--         origin = {0., 0., 0.},
--         length = 1.0,
--         refinementLevel = 2  }

-- If restart table and meah are present then priority is given to restart
norestart = {
  read = 'final100_restart_ProtoData_header_100.lua'
}


-- Subsampling (for subresoloved color information):
ply_sampling = { nlevels = 0,    -- maximal level to use in subsampling
                 --method  = 'fixed' -- method to use for subsampling
                                     -- currently only 'fixed' is implemented,
                                     -- which will refine all elements by
                                     -- nlevels
                                     -- 'fixed' is also the default, thus it
                                     -- is sufficient to only provide nlevels
}

-- Verbosity of logging:
logging = {level=2}

-- Define tracking objects to further restrict the created output data.
-- If no tracking is defined, all the information will be written 'as is' to
-- the output files as given by the output table below.
-- For example in the following, we only write out the treeIDs.
tracking = {
    label = 'allIDs',      -- a label for the tracking object,
                           -- default: 'unnamed'
    folder = 'output',      -- output folder name
                           -- if it ends in a '/', this is a directory,
                           -- which needs to exist!
    output = { 
      format = 'vtk',      -- Format of the output data (set to harvest!)
      write_pvd = false    -- write pvd file containing all time steps
                           -- Default is true
    },
    variable = {'treeid'}, -- variables to track
    --variable = {'treeid'}, -- variables to track
    shape = {kind='all'},  -- shape to track, 'all' or canoND are typical
                           -- other options see treelm/config.lua for examples.
}

-- To dump complete mesh or restart file with out any tracking table
-- just default output table with basename, format as shown below
-- Configuration for the output to produce
output_folder = 'output'    -- Prefix to use for the produced files,
                            -- if it ends in a '/', this is a directory,
                            -- which needs to exist!

NOoutput = {   
           format   = 'vtk',     -- Visualization file format, defaults to 'vtk'
                                 -- Currently only vtk is supported here.
           dataform = 'binary',   -- Form to write data in:
                                 --   - 'binary',  or
                                 --   - 'ascii'
                                 -- Default is 'binary'
           write_pvd = false     -- write pvd file containing all time steps
                                 -- Default is true
}


