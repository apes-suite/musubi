simulation_name = 'channel'
debug = false

-- define the output

if debug then
  -- define the input
  restart = { read = 'final100_restart_ProtoData_header_100.lua'  }
else
  -- define the input
  mesh = 'mesh/'
end

-- define the output
output_folder = 'mesh/'     -- Output location
output = {format = 'vtk'}   -- Output format

