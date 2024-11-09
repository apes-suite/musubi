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

# Axis without scientific notation
y_formatter = mtick.ScalarFormatter(useOffset=False)

# Inflow velocity
from subprocess import Popen, PIPE
mus_file = 'musubi.lua'
vel_phy = float(Popen(['lua', '-e',"dofile '"+mus_file+"'; \
                       print(string.format('%.2E',vel_phy))"],\
                       stdout=PIPE).communicate()[0].decode('ascii'))
## -------------------------------------------------------------------------- ##
logging.basicConfig(level=logging.INFO)
logging.info('Started creating plots ...') 

# Load reference data
import numpy as np
refFile_velX = 'referenceDatas/velocityX_GeometryCenter.txt'
y_ref, velX_ref = np.loadtxt(refFile_velX, usecols=(1,3), comments='#', \
                             unpack=True)
refFile_velY = 'referenceDatas/velocityY_GeometryCenter.txt'
x_ref, velY_ref = np.loadtxt(refFile_velY, usecols=(1,3), comments='#', \
                             unpack=True)

# Data base filename
dbname = 'lidDrivenCavity.db'
# Load database if exist else load tracking files and add to database
if os.path.isfile(dbname):
  logging.info('Processing data from existing database')
#  os.remove(dbname)
  import sqlite3
  sqlcon = sqlite3.connect(dbname)
else:
  logging.info('Processing data from tracking files')

  # Load tracking output with label probe and store in tabname=probe
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*probe*.res'], \
                                  dbname=dbname, tabname='probe')

  # Load tracking output with verticalLine and store in tabname=vertline
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*verticalLine*.res'], \
                                  dbname=dbname, tabname='vertline')

  # Load tracking output with horizondalLine and store in tabname=horiline
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*horizondalLine*.res'], \
                                  dbname=dbname, tabname='horiline')
## -------------------------------------------------------------------------- ##
logging.info('Velocity-X along the cavity height:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-Y, velocity magnitude (from simulation)
# and analytical velocity
get_data_for_cols = ['coordY','velocity_phy_01']
[x, y] = gleaner.get_columns(sqlcon, tabname='vertline', \
                             columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
y_norm = []
for ii in y:
  y_norm.append((ii[0]/vel_phy))
# Plot simulation result
mplt.plot(x, y_norm, '-', color = 'k', label = 'Simulation')
# Plot Reference data 
mplt.plot(y_ref, velX_ref, 's', color = 'r', label = 'Ghia et al')

# plot setting
mplt.legend(loc=9, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('y/L', fontsize=font_size)
mplt.ylabel('$u_x/u_L$', fontsize=font_size)
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,1.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelocityXAlongHeight.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Velocity-Y along the cavity length:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-Y, velocity magnitude (from simulation)
# and analytical velocity
get_data_for_cols = ['coordX','velocity_phy_02']
[x, y] = gleaner.get_columns(sqlcon, tabname='horiline', \
                             columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
y_norm = []
for ii in y:
  y_norm.append((ii[0]/vel_phy))

# Plot simulation result
mplt.plot(x, y_norm, '-', color = 'k', label = 'Simulation')
# Plot Reference data 
mplt.plot(x_ref, velY_ref, 's', color = 'r', label = 'Ghia et al')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x/L', fontsize=font_size)
mplt.ylabel('$u_y/u_L$', fontsize=font_size)
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,1.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelocityYAlongLength.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Velocity-X at center of cavity over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract simulation time and pressure
get_data_for_cols = ['time','velocity_phy_01']
[x, y] = gleaner.get_columns(sqlcon, tabname='probe', \
                             columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x,y))) # sort of needed
# Plot simulation result
y_norm = []
for ii in y:
  y_norm.append((ii[0]/vel_phy))

mplt.plot(x, y_norm, '-', color = 'r')

# plot setting
mplt.xlabel('time (s)', fontsize=font_size)
mplt.ylabel('u_x/u_L', fontsize=font_size)
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelocityXOverTime.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##

logging.info('Plots created')
