-- Local variable
-- If true post-process debug output from seeder else post-process mesh files
debug = false

-- Simulation name to append to output file 
simulation_name = 'channel'

-- define the output
if debug then
  restart = { read = 'final100_restart_ProtoData_header_100.lua'  }
else
  mesh = 'mesh/'
end

-- define the output
tracking = {
    label = 'mesh',
    variable = { 'treeid', 'level' },
    shape = { kind='all' },
    folder = 'mesh/'
    output = { format = 'vtk', write_pvd = false },
}
