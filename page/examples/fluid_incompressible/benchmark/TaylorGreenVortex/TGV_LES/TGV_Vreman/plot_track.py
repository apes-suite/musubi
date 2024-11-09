## This is the user-script for plotting using gleaner tool.
import sys
import os
import math

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
dbname = 'TGV_Vreman.db'
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

  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*kE_all*.res'], \
                                  dbname=dbname, tabname='kE')
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
mplt.plot(x, y, '-', color = 'b', label='Musubi')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('Time (s)')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,20.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureOverTime.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Kinetic energy over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract simulation time and pressure
get_data_for_cols = ['time','kinetic_energy_phy_red']
[time, ke] = gleaner.get_columns(sqlcon, tabname='kE', \
                             columns=get_data_for_cols)
# sort loaded data according to x
time, ke = zip(*sorted(zip(time,ke))) # sort of needed
volume = ( 2.0 * math.pi ) ** 3
ke_list = []
time_list = []
for i in range(len(ke)):
  time_list.append(time[i][0])
  ke_list.append( ke[i][0]/volume )

# Plot simulation result
mplt.plot(time_list, ke_list, '-', color = 'b', label = 'Musubi')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('Time (s)')
mplt.ylabel('Kinetic energy')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,10.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('KineticEnergy.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Dissipation rate over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Calculate dissipation rate from kinetic energy
# -dE/dt = ( E(t) - E(t+1) ) / dt
dr_list = []
for i in range(len(ke)-1):
  dr_list.append( (ke_list[i]-ke_list[i+1]) / (time_list[i+1]-time_list[i]) )

# Plot simulation result
mplt.plot(time_list[:len(ke)-1], dr_list, '-', color = 'b', label = 'Musubi')

data = np.loadtxt('./Re1600DNS_Brachet.inp')
mplt.plot(data[:,0], data[:,1], '-', color = 'k', label = 'Brachet')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('Time (s)')
mplt.ylabel('Dissipation rate')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,10.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('DissipationRate.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##


logging.info('Plots created')
