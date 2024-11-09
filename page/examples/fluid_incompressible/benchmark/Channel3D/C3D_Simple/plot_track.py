## This is the user-script for plotting using gleaner tool.
import sys
import os

# Path to gleaner (Better use environment variable PYTHONPATH!)
if os.path.exists( os.getenv('HOME')+'/apes/gleaner'):
  glrPath = os.getenv('HOME')+'/apes/gleaner'
else:
  print('Gleaner library not found')
  exit

# Import gleaner module
sys.path.append(glrPath)
import gleaner

# Do not use X-server to create and save plot
import matplotlib
matplotlib.use('Agg')

## Import all required modules
import matplotlib.ticker as mtick
import matplotlib.pyplot as mplt
import numpy as np
import logging

# font setting
from matplotlib import rc
font_size = 12

#axis without scientific notation
y_formatter = mtick.ScalarFormatter(useOffset=False)

## -------------------------------------------------------------------------- ##
logging.basicConfig(level=logging.INFO)
logging.info('Started creating plots ...')

# data base filename
dbname = 'C3D_Simple.db'
# load database if exist else load tracking files and add to database
if os.path.isfile(dbname):
  logging.info('Processing data from existing database')
#  os.remove(dbname)
  import sqlite3
  sqlcon = sqlite3.connect(dbname)
else:
  logging.info('Processing data from tracking files')

  # Load tracking output with label probeAtCenter and store in tabname=probe
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*probeAtCenter*.res'], \
                                  dbname=dbname, tabname='probe')

  # Load tracking output with pressAlongLength and store in tabname=press_line
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*pressAlongLength*.res'], \
                                  dbname=dbname, tabname='press_line')

  # Load tracking output with velAlongHeight and store in tabname=vel_spatial
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*velAlongHeight*.res'], \
                                  dbname=dbname, tabname='vel_spatial')

  # Load tracking output with wssAlongHeight and store in tabname=wss_spatial
  #sqlcon = gleaner.tracking_to_db(fname = ['tracking/*wssAlongHeight*.res'], \
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*wss_spatial*.res'], \
                                  dbname=dbname, tabname='wss_spatial')

  # Reference from musubi default
  # Load tracking output with pressAlongLength and store in tabname=press_line
  sqlcon = gleaner.tracking_to_db(fname = ['reference/converged/*pressAlongLength*.res'], \
                                  dbname=dbname, tabname='ref_press_line')

  # Load tracking output with velAlongHeight and store in tabname=vel_spatial
  sqlcon = gleaner.tracking_to_db(fname = ['reference/converged/*velAlongHeight*.res'], \
                                  dbname=dbname, tabname='ref_vel_spatial')

## -------------------------------------------------------------------------- ##
logging.info('Pressure over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract simulation time and pressure
get_data_for_cols = ['time','pressure_phy']
[x, y] = gleaner.get_columns(sqlcon, tabname='probe', \
                             columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x,y))) # sort of needed
# Plot simulation result
mplt.plot(x, y, '-', color = 'b')

# plot setting
mplt.xlabel('time (s)')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,25.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureOverTime.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
### -------------------------------------------------------------------------- ##
logging.info('Velocity X over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract simulation time and velocity-X
get_data_for_cols = ['time','velocity_phy_01']
[x, y] = gleaner.get_columns(sqlcon, tabname='probe', \
                             columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x,y))) # sort of needed
mplt.plot(x, y, '-', color = 'b')

# plot setting
mplt.xlabel('time (s)')
mplt.ylabel('Velocity X ($m/s$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,25.0)
mplt.ylim(0.0,0.6)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelocityXOverTime.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Pressure drop across the channel length:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-X, pressure (from simulation)
get_data_for_cols = ['coordX','pressure_phy']
[x, y] = gleaner.get_columns(sqlcon, tabname='press_line', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'b', label = 'Simulation')

# Reference from Musubi
[x, y] = gleaner.get_columns(sqlcon, tabname='ref_press_line', \
                                  columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'r', label = 'Reference Musubi')

# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,2.5)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('Pressure_Profile.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Velocity Mag profile across the channel height:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-Y, velocity magnitude (from simulation)
get_data_for_cols = ['coordY','vel_mag_phy']
[x, y] = gleaner.get_columns(sqlcon, tabname='vel_spatial', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'b', label = 'Simulation')

# reference data
[x, y] = gleaner.get_columns(sqlcon, tabname='ref_vel_spatial', \
                                  columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, 'x', color = 'r', markevery=2, label = 'Reference Musubi')

## initial condition
height = 0.41
u_m = 0.45
y = np.arange(0, height+0.01, 0.01)
z = 0.5* height
x = 4 * u_m * y * ( height - y) / height**2
x = 16 * u_m * y * z * (height - y) * ( height - z ) / height**4
mplt.plot(y, x, 'o', color = 'm', markevery=1, label = 'Initial condition')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('y (m)')
mplt.ylabel('Velocity ($m/s$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,0.41)
#mplt.ylim(0.05,0.45)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('Velocity_Profile.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
#### -------------------------------------------------------------------------- ##
#logging.info('Wall shear stress profile across the channel height:')
#fig = mplt.figure()
#ax = fig.add_subplot(111)
## Plot x, y ... at certain time step
## Extract coordinate-Y, wall shear stress (from simulation)
## and analytical wall shear stress
#get_data_for_cols = ['coordY','wss_phy','wss_analy']
#[x, y,z] = gleaner.get_columns(sqlcon, tabname='wss_spatial', \
#                             columns=get_data_for_cols)
## sort loaded data according to x
#x, y,z = zip(*sorted(zip(x,y,z))) # sort of needed
## Plot simulation result
#mplt.plot(x, y, '-s', color = 'k', label = 'Simulation')
## Plot Analytical result
#mplt.plot(x, z, '-', color = 'r', label = 'Analytic')
#
## plot setting
#mplt.legend(loc=5, ncol=1,borderaxespad=0, \
#            prop={'size':font_size}).get_frame().set_lw(0.0)
#mplt.xlabel('x (m)')
#mplt.ylabel('WSS ($Pa$)')
#mplt.grid(True,which="major",ls="-")
#ax.yaxis.set_major_formatter(y_formatter)
#mplt.xlim(0.0,1.0)
#
## save fig
#figsize = [8,6]
#fig = mplt.gcf()
#fig.set_size_inches(figsize[0],figsize[1])
#mplt.savefig('WSS_Profile.png', dpi=100, format='png', \
#             bbox_inches="tight",interpolation=None)
### -------------------------------------------------------------------------- ##
logging.info('Plots created')
