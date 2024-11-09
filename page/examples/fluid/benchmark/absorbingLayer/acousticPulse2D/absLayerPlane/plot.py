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
font_family = 'serif'
font_type = 'Times New Roman'
#rc('text',usetex=True)
#font = {'family':font_family,'%s'%font_family:font_type,'size':font_size}
#rc('font',**font)

#axis without scientific notation
y_formatter = mtick.ScalarFormatter(useOffset=False)

## -------------------------------------------------------------------------- ##
logging.basicConfig(level=logging.INFO)

# data base filename
dbname = 'acousticPulse2D-box.db'
# load database if exist else load tracking files and add to database
if os.path.isfile(dbname):
  print ('Processing data from existing database')
#  os.remove(dbname)
  import sqlite3
  sqlcon = sqlite3.connect(dbname)
else:
  print ('Processing data from tracking files')

  sqlcon = gleaner.tracking_to_db(fname = ['tracking/pulse2D_probe_p00000.res'], \
                                  dbname=dbname, tabname='probe')
## -------------------------------------------------------------------------- ##
#print ('Pressure along the center line:')
#fig = mplt.figure()
#ax = fig.add_subplot(111)
## Plot x, y ... at certain time step
## Simulation result
#get_data_for_cols = ['coordX','pressure_phy']
#[x, y] = gleaner.get_columns(sqlcon, tabname='line', \
#                             columns=get_data_for_cols)
#x, y = zip(*sorted(zip(x,y))) # sort of needed
#mplt.axvline(x=4.0, color = 'r', label = 'absLayer start')
#
#mplt.plot(x, y, '-', color = 'b', label = 'Musubi')
#
## plot setting
#mplt.legend(loc=8, ncol=1,borderaxespad=0, \
#            prop={'size':font_size}).get_frame().set_lw(0.0)
#mplt.xlabel('x (m)')
#mplt.ylabel('Pressure')
#mplt.grid(True,which="major",ls="-")
#ax.yaxis.set_major_formatter(y_formatter)
#mplt.xlim(0.0,4.8)
#
## save fig
#figsize = [8,6]
#fig = mplt.gcf()
#fig.set_size_inches(figsize[0],figsize[1])
#mplt.savefig('PressureAlongLine.png', dpi=100, format='png', \
#             bbox_inches="tight")
#
### -------------------------------------------------------------------------- ##
print ('Pressure over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Plot x, y ... at certain time step
# Simulation result
get_data_for_cols = ['time','pressure_phy']
[x, y] = gleaner.get_columns(sqlcon, tabname='probe', \
                             columns=get_data_for_cols)
x, y = zip(*sorted(zip(x,y))) # sort of needed
#mplt.axvline(x=4.0, color = 'r', label = 'absLayer start')

mplt.plot(x, y, '-', color = 'b', label = 'Musubi')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('Pressure')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
#mplt.xlim(0.0,4.8)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureOverTime.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Plots created')
