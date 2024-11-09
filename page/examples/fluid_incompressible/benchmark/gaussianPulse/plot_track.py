# This is the user-script for plotting using gleaner tool.
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
## -------------------------------------------------------------------------- ##
from subprocess import Popen, PIPE
# Simulation parameters for calculation of analytical solution of gaussian pulse
simfile = 'musubi.lua'
length = float(Popen(['lua', '-e',"dofile '"+simfile+"'; \
               print(string.format('%.2E', length))"],\
               stdout=PIPE).communicate()[0].decode('ascii'))

centerX = float(Popen(['lua', '-e',"dofile '"+simfile+"'; \
                print(string.format('%.2E', centerX))"],\
                stdout=PIPE).communicate()[0].decode('ascii'))

halfwidth = float(Popen(['lua', '-e',"dofile '"+simfile+"'; \
                  print(string.format('%.2E', halfwidth))"],\
                  stdout=PIPE).communicate()[0].decode('ascii'))

amplitude = float(Popen(['lua', '-e',"dofile '"+simfile+"'; \
                  print(string.format('%.2E', amplitude))"],\
                  stdout=PIPE).communicate()[0].decode('ascii'))

background = float(Popen(['lua', '-e',"dofile '"+simfile+"'; \
                   print(string.format('%.6E', background))"],\
                   stdout=PIPE).communicate()[0].decode('ascii'))

## -------------------------------------------------------------------------- ##

logging.info('Started creating plots ...')

# data base filename
dbname = 'gaussianPulse.db'
# load database if exist else load tracking files and add to database
if os.path.isfile(dbname):
  logging.info('Processing data from existing database')
#  os.remove(dbname)
  import sqlite3
  sqlcon = sqlite3.connect(dbname)
else:
  logging.info('Processing data from tracking files')

  # Load tracking output with label probeAt1 and store in tabname=probe1
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*probeAt1*.res'], \
                                  dbname=dbname, tabname='probe1')

  # Load tracking output with pressAlongLength at diffent points in time and
  # store it in tabname=press_line_t*

  # t = 0
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*pressAlongLength*t0*.res'], \
                                  dbname=dbname, tabname='press_line_t0')
  # t = 1s
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*pressAlongLength*t10*.res'], \
                                  dbname=dbname, tabname='press_line_t10')


  ##### Load reference data obtained with finer resolutions ####
  ##### (tracking output is with level 4. ####
  # Load tracking output with pressAlongLength at t=0s and t=10s for level 5
  sqlcon = gleaner.tracking_to_db(fname = ['./reference/*L5*pressAlongLength*t0*.res'], \
                                  dbname=dbname, tabname='press_line_t0_L5')

  sqlcon = gleaner.tracking_to_db(fname = ['./reference/*L5*pressAlongLength*t10*.res'], \
                                  dbname=dbname, tabname='press_line_t10_L5')

  # Load tracking output with pressAlongLength at t=0s and t=10s for level 6
  sqlcon = gleaner.tracking_to_db(fname = ['reference/*L6*pressAlongLength*t0*.res'], \
                                  dbname=dbname, tabname='press_line_t0_L6')

  sqlcon = gleaner.tracking_to_db(fname = ['reference/*L6*pressAlongLength*t10*.res'], \
                                  dbname=dbname, tabname='press_line_t10_L6')

## -------------------------------------------------------------------------- ##
logging.info('Pressure over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract simulation time and pressure
get_data_for_cols = ['time','pressure_phy']
[x, y] = gleaner.get_columns(sqlcon, tabname='probe1', \
                             columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x,y))) # sort of needed
# Plot simulation result
mplt.plot(x, y, '-', color = 'b', label='(1.0,5.0,5.0)')

# plot setting
mplt.xlabel('Time')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlim(0.0,1.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureOverTime.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Pressure profile across the channel length:')
fig = mplt.figure()
ax = fig.add_subplot(111)

# Analytic initial condition
x = np.arange(0, length, 0.01)
r = (x-centerX)**2
y = background + amplitude * np.exp( (-np.log(2.0) / (halfwidth**2)) * r)
## Plot analytic solution
mplt.plot(x, y, '-', color = 'm', label = 'Analytic')

# Extract coordinate-X, pressure (from simulation)
# t = 0
get_data_for_cols = ['coordX','pressure_phy']
[x, y] = gleaner.get_columns(sqlcon, tabname='press_line_t0', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'b', label = 'L4')

# For level 5 data
[x, y] = gleaner.get_columns(sqlcon, tabname='press_line_t0_L5', \
                            columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'c', label = 'L5')

# For level 6 data
[x, y] = gleaner.get_columns(sqlcon, tabname='press_line_t0_L6', \
                            columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'k', label = 'L6')

# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,10.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('Pressure_Profile_t0.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Pressure profile across the channel length over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)

# Extract coordinate-X, pressure (from simulation)
get_data_for_cols = ['coordX','pressure_phy']
# t = 10s
[x, y] = gleaner.get_columns(sqlcon, tabname='press_line_t10', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-x', color = 'b', label = 'L4')

# For level 5 data
[x, y] = gleaner.get_columns(sqlcon, tabname='press_line_t10_L5', \
                            columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-x', color = 'c', label = 'L5')

# For level 6 data
[x, y] = gleaner.get_columns(sqlcon, tabname='press_line_t10_L6', \
                            columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-x', color = 'k', label = 'L6')

# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,10.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('Pressure_Profile_t10.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Plots created')
