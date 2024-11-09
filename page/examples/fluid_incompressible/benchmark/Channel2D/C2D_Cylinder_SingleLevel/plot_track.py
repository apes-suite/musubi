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
import numpy as np

# font setting
from matplotlib import rc
font_size = 12

#axis without scientific notation
y_formatter = mtick.ScalarFormatter(useOffset=False)

## -------------------------------------------------------------------------- ##
logging.basicConfig(level=logging.INFO)
logging.info('Started creating plots ...') 

# data base filename
dbname = 'channelCyl.db'
# load database if exist else load tracking files and add to database
if os.path.isfile(dbname):
  logging.info('Processing data from existing database')
#  os.remove(dbname)
  import sqlite3
  sqlcon = sqlite3.connect(dbname)
else:
  logging.info('Processing data from tracking files')

  # Load tracking output with lift and drag coeff and store in tabname=cyl_force
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*cyl_force*.res'], \
                                  dbname=dbname, tabname='cyl_force')

  # Load tracking output along centerLine and store in tabname=press_line
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*centerLine*.res'], \
                                  dbname=dbname, tabname='center_line')

  # Load tracking output pressure before cylinder
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*probePressAtCylBack_*.res'], \
                                  dbname=dbname, tabname='probePressAtCylBack')

  # Load tracking output pressure after cylinder
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*probePressAtCylFront_*.res'], \
                                  dbname=dbname, tabname='probePressAtCylFront')

  # Load tracking output with lift and drag coeff and store in tabname=cyl_force
  sqlcon = gleaner.tracking_to_db(fname = ['tracking/*Cp*.res'], \
                                  dbname=dbname, tabname='cp')
## -------------------------------------------------------------------------- ##
logging.info('Pressure coeff Cp over angle:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-Y, velocity magnitude (from simulation)
# and analytical velocity
get_data_for_cols = ['coordX','coordY','coordZ','coeffPressureAvg']
[x, y, z, Cp] = gleaner.get_columns(sqlcon, tabname='Cp', \
                                columns=get_data_for_cols)

# Average Cp over z
import numpy as np
import math
radius = 0.05
diameter = radius*2.0
# Channel height [m]
height = 0.41
# cylinder offset from origin
cylinder_x = 0.2
cylinder_y = 0.2
Cp_avg= dict()
theta_list = []
ii_new = 0
for ii in range(len(Cp)):
  theta = 180-math.atan2(y[ii] - cylinder_y,x[ii] - cylinder_x)*180/math.pi
  if theta in Cp_avg:
    Cp_avg[theta].append(Cp[ii])
  else:
    Cp_avg[theta] = [Cp[ii]]
    theta_list.append(theta)
#print(theta_list)
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
#mplt.xlim(0.0,180)
#mplt.ylim(-1.5,1.5)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureCoeffOverAngle.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Drag coefficient over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract time, coeff_red_01
get_data_for_cols = ['time','coeff_red_01']
[x, y] = gleaner.get_columns(sqlcon, tabname='cyl_force', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'k')
mplt.axhline(y=3.22, color = 'b', linestyle = '--')
mplt.axhline(y=3.24, color = 'b', linestyle = '--')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('t (s)')
mplt.ylabel('$C_D$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(5.0,7.5)
mplt.ylim(3.0,3.5)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('DragCoeffOverTime.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Lift coefficient over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract time, coeff_red_01
get_data_for_cols = ['time','coeff_red_02']
[x, y] = gleaner.get_columns(sqlcon, tabname='cyl_force', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'k')
mplt.axhline(y=1.0, color = 'b', linestyle = '--')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('t (s)')
mplt.ylabel('$C_L$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(5.0,7.5)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('LiftCoeffOverTime.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Pressure profile across the channel length:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-X and pressure
get_data_for_cols = ['coordX','pressure_phy']
[x, y] = gleaner.get_columns(sqlcon, tabname='center_line', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y= zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'k')

# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,0.41*5)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureAlongLength.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Velocity X profile across the channel length:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-X and pressure
get_data_for_cols = ['coordX','velocity_phy_01']
[x, y] = gleaner.get_columns(sqlcon, tabname='center_line', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y= zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'k')

# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('$u_x$ ($m/s$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,0.41*5)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelocityXAlongLength.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Velocity Y profile across the channel length:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract coordinate-X and pressure
get_data_for_cols = ['coordX','velocity_phy_02']
[x, y] = gleaner.get_columns(sqlcon, tabname='center_line', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y= zip(*sorted(zip(x, y)))
# Plot simulation result
mplt.plot(x, y, '-', color = 'k')

# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('$u_y$ ($m/s$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,0.41*5)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelocityYAlongLength.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Pressure difference between the cylinder:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract time and pressure
get_data_for_cols = ['time','pressure_phy']
[x, p1] = gleaner.get_columns(sqlcon, tabname='probePressAtCylBack', \
                              columns=get_data_for_cols)
# sort loaded data according to x
x, p1 = zip(*sorted(zip(x, p1)))
[x, p2] = gleaner.get_columns(sqlcon, tabname='probePressAtCylFront', \
                              columns=get_data_for_cols)
# sort loaded data according to x
x, p2 = zip(*sorted(zip(x, p2)))

# Calculate pressure difference
dp = np.array([])
for ii in range(len(x)):
  dp = np.append(dp, p2[ii]-p1[ii])

# Plot simulation result
mplt.plot(x, dp, '-', color = 'k')
mplt.axhline(y=2.46, color = 'b', linestyle = '--')
mplt.axhline(y=2.50, color = 'b', linestyle = '--')

# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('t (s)')
mplt.ylabel('$\Delta P$ ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(5.0,7.5)
mplt.ylim(2.0,3.0)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureDiffOverTime.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Lift coefficient in Frequency:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Extract time, coeff_red_01
get_data_for_cols = ['time','coeff_red_02']
[x, y] = gleaner.get_columns(sqlcon, tabname='cyl_force', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))

# convert to numpy data
cL = np.array([])
for ii in range(len(x)):
  if x[ii] > 4.0:
    cL = np.append(cL, y[ii])
nData = len(cL)
# time difference between two outputs
interval = x[-2] - x[-1]
fft = 2.0*np.absolute(np.fft.fft(cL, n=nData)/nData)
freq = np.fft.fftfreq(nData, d=interval)
frequency = abs(freq[fft.argmax()])
logging.info('Frequency with maximum amplitude: '+str(frequency))
avg_Cl = cL.mean()
logging.info('Average Cl: '+str(avg_Cl))

# Extract time, coeff_red_01
get_data_for_cols = ['time','coeff_red_01']
[x, y] = gleaner.get_columns(sqlcon, tabname='cyl_force', \
                                columns=get_data_for_cols)
# sort loaded data according to x
x, y = zip(*sorted(zip(x, y)))

# convert to numpy data
cD = np.array([])
for ii in range(len(x)):
  if x[ii] > 4.0:
    cD = np.append(cD, y[ii])
nData = len(cD)

avg_CD = cD.mean()
logging.info('Average Cd: '+str(avg_CD))

# Plot simulation result
mplt.plot(freq, fft, '-', color = 'k')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('frequency (Hz)')
mplt.ylabel('$C_L$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,7.5)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('FFT.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##

logging.info('Plots created')
