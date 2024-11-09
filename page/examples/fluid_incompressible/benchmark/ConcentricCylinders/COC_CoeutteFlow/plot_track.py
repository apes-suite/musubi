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
dbname = 'circularCyl.db'
# load database if exist else load tracking files and add to database
if os.path.isfile(dbname):
  logging.info('Processing data from existing database')
#  os.remove(dbname)
  import sqlite3
  sqlcon = sqlite3.connect(dbname)
else:
  logging.info('Processing data from tracking files')

  # Load tracking output with line and store in tabname=line
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*line*.res'], \
                                  dbname=dbname, tabname='line')


## -------------------------------------------------------------------------- ##
logging.info('Velocity Mag profile across the radius:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-Y, velocity magnitude (from simulation)
# and analytical velocity
get_data_for_cols = ['coordX','vel_mag_phy','vel_an']
[x, y, z] = gleaner.get_columns(sqlcon, tabname='line', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y, z = zip(*sorted(zip(x, y, z)))
# Plot simulation result
mplt.plot(x, y, '-s', color = 'k', label = 'Simulation')
# Plot Analytical result
mplt.plot(x, z, '-', color = 'r', label = 'Analytic')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('Velocity ($m/s$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
#mplt.xlim(0.0,1.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('media/Velocity_Profile.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##

logging.info('Plots created')
