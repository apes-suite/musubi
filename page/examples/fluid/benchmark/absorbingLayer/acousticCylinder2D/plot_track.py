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
#matplotlib.use('Agg')

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

time_stamp = '42.272E-03'
# data base filename
dbname = 'channel.db'
# load database if exist else load tracking files and add to database
if os.path.isfile(dbname):
  logging.info('Processing data from existing database')
#  os.remove(dbname)
  import sqlite3
  sqlcon = sqlite3.connect(dbname)
else:
  logging.info('Processing data from tracking files')

  # Load tracking output with label Cp - pressure coefficient
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*Cp_p*_t'+time_stamp+'.res'], \
                                  dbname=dbname, tabname='Cp')

  # Load tracking output with label cyl_force - Lift and drag coefficient
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*cyl_force*.res'], \
                                  dbname=dbname, tabname='cLcDCoeff')
## -------------------------------------------------------------------------- ##
logging.info('Pressure coeff Cp over angle:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-Y, velocity magnitude (from simulation)
# and analytical velocity
get_data_for_cols = ['coordX','coordY','coordZ','press_coeff_avg']
[x, y, z, Cp] = gleaner.get_columns(sqlcon, tabname='Cp', \
                                columns=get_data_for_cols)

# Average Cp over z
import numpy as np
import math
Cp_avg= dict()
theta_list = []
ii_new = 0
for ii in range(len(Cp)):
  theta = 180-math.atan2(y[ii][0],x[ii][0])*180/math.pi
  if theta in Cp_avg:
    Cp_avg[theta].append(Cp[ii][0])
  else:
    Cp_avg[theta] = [Cp[ii][0]]
    theta_list.append(theta)

Cp_avg_list = []
for theta in theta_list:
  Cp_theta = np.array(Cp_avg[theta])
  Cp_avg_list.append(np.average(Cp_theta))

## sort loaded data according to theta
theta_list, Cp_avg_list = zip(*sorted(zip(theta_list, Cp_avg_list)))
## Plot simulation result
mplt.plot(theta_list, Cp_avg_list, '-s', color = 'k')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$\\theta$')
mplt.ylabel('$C_p$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,180)
mplt.ylim(-1.5,1.5)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureCoeffOverAngle.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Drag coeff C_d over time')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-Y, velocity magnitude (from simulation)
# and analytical velocity
get_data_for_cols = ['time','coeff_red_01']
[t, Cd] = gleaner.get_columns(sqlcon, tabname='cLcDCoeff', \
                                columns=get_data_for_cols)

## Plot simulation result
mplt.plot(t, Cd, '-', color = 'k')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$time [s]$')
mplt.ylabel('$C_D$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.035,0.04)
mplt.ylim(1.25,1.45)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('DragCoeffOverTime.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Lift coeff C_l over time')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-Y, velocity magnitude (from simulation)
# and analytical velocity
get_data_for_cols = ['time','coeff_red_02']
[t, Cl] = gleaner.get_columns(sqlcon, tabname='cLcDCoeff', \
                                columns=get_data_for_cols)

## Plot simulation result
mplt.plot(t, Cl, '-', color = 'k')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$time [s]$')
mplt.ylabel('$C_L$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.035,0.04)
mplt.ylim(-0.55,0.55)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('LiftCoeffOverTime.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##

mplt.show()
logging.info('Plots created')
