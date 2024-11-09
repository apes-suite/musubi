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
dbname = 'pipe.db'
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

  # Load tracking output with planeCenter and store in tabname=plane_center
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*planeCenter*.res'], \
                                  dbname=dbname, tabname='plane_center')

  # Load tracking output with pressAlongLength and store in tabname=press_line
  sqlcon = gleaner.tracking_to_db(fname = ['reference/musubi_default_branch/*pressAlongLength*.res'], \
                                  dbname=dbname, tabname='ref_press_line')

  # Load tracking output with velAlongHeight and store in tabname=vel_spatial
  sqlcon = gleaner.tracking_to_db(fname = ['reference/musubi_default_branch/*velAlongHeight*.res'], \
                                  dbname=dbname, tabname='ref_vel_spatial')
  # reference probe at center
  sqlcon = gleaner.tracking_to_db(fname = ['reference/*probeAtCenter*.res'], \
                                  dbname=dbname, tabname='ref_probe')


## -------------------------------------------------------------------------- ##
logging.info('Velocity Mag profile across the channel height:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-Y, velocity magnitude (from simulation)
# and analytical velocity
get_data_for_cols = ['coordY','vel_mag_phy','vel_analy']
[x, y, z] = gleaner.get_columns(sqlcon, tabname='vel_spatial', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y, z = zip(*sorted(zip(x, y, z)))
# Plot simulation result
mplt.plot(x, y, '-s', color = 'k', label = 'Simulation')
# Plot Analytical result
mplt.plot(x, z, '-', color = 'r', label = 'Analytic')

# Plot reference from Musubi
#[x, y, z] = gleaner.get_columns(sqlcon, tabname='ref_vel_spatial', \
#                                columns=get_data_for_cols)
## sort loaded data according to x
#x, y, z = zip(*sorted(zip(x, y, z)))
## Plot simulation result
#mplt.plot(x, y, '-s', color = 'b', label = 'ReferenceMusubi')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('y (m)')
mplt.ylabel('Velocity ($m/s$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('Velocity_Profile.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Pressure drop across the channel length:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-X, pressure (from simulation) and analytical pressure
get_data_for_cols = ['coordX','pressure_phy','press_analy']
[x, y, z] = gleaner.get_columns(sqlcon, tabname='press_line', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y, z = zip(*sorted(zip(x, y, z)))
# Plot simulation result
mplt.plot(x, y, '-s', color = 'k', label = 'Simulation')
## Plot Analytical result
mplt.plot(x, z, '-', color = 'r', label = 'Analytic')

# Plot reference from Musubi
#[x, y, z] = gleaner.get_columns(sqlcon, tabname='ref_press_line', \
#                                columns=get_data_for_cols)
## sort loaded data according to x
#x, y, z = zip(*sorted(zip(x, y, z)))
## Plot simulation result
#mplt.plot(x, y, '-s', color = 'b', label = 'ReferenceMusubi')


# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('Pressure_Profile.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('3D Velocity profile in plane at center:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
ax = Axes3D(fig)
# Extract coordY, coordZ, vel_mag_phy
get_data_for_cols = ['coordY', 'coordZ', 'vel_mag_phy']
[x, y, z] = gleaner.get_columns(sqlcon, tabname='plane_center', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y, z = zip(*sorted(zip(x,y,z))) # sort of needed
import numpy as np
xx = np.array([])
yy = np.array([])
zz = np.array([])
for ii in range(len(x)):
  xx = np.append(xx, x[ii])
  yy = np.append(yy, y[ii])
  zz = np.append(zz, z[ii])
  
# Plot simulation result
from matplotlib import cm
#ax.plot_surface(xx, yy, zz, cmap='viridis', edgecolor='none')
surf = ax.plot_trisurf(xx, yy, zz, cmap=cm.jet, linewidth=0.1)
#surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

# plot setting
mplt.legend(loc=5, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('y (m)')
mplt.ylabel('z (m)')
#mplt.zlabel('$|u| (m/s)$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
#mplt.xlim(0.0,1.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelMag_3DProfile.png', dpi=100, format='png', \
             bbox_inches="tight")
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
mplt.plot(x, y, '-', color = 'r', label='Simulation')

[x, y] = gleaner.get_columns(sqlcon, tabname='ref_probe', \
                             columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x,y))) # sort of needed
# Plot simulation result
mplt.plot(x, y, '-', color = 'b', label='ReferenceMusubi')


# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('time (s)')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureOverTime.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Velocity X over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract simulation time and velocity-X
get_data_for_cols = ['time','velocity_phy_01']
[x, y] = gleaner.get_columns(sqlcon, tabname='probe', \
                             columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x,y))) # sort of needed
mplt.plot(x, y, '-', color = 'r', label='Simulation')
[x, y] = gleaner.get_columns(sqlcon, tabname='ref_probe', \
                             columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x,y))) # sort of needed
mplt.plot(x, y, '-', color = 'b', label='ReferenceMusubi')


# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('time (s)')
mplt.ylabel('Velocity X ($m/s$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelocityXOverTime.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##

logging.info('Plots created')
