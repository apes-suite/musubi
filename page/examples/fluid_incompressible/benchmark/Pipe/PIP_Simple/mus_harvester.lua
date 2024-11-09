require 'musubi'

restart = {
  read = 'restart/pipe_lastHeader.lua'
}

-- define the output
--output_folder = 'tracking/'  -- Output location
--output = { format = 'vtk' }   -- Output format
tracking = {
  label = 'vtk_hvs',
  folder = tracking_fol,
  variable = { 'pressure_phy', 'velocity_phy' },
  shape = {
    kind = 'canoND',
    object = {
      origin = {0.0, pipe_origin[2], pipe_center[3]},
      vec = {
        {length, 0.0, 0.0},
        {0.0, diameter, 0.0}
      }
    }
  },
  time_control = {
    min = tmax_phy,
    max = tmax_phy,
    interval = interval_phy
  },
  output = {format = 'vtk'}
}

