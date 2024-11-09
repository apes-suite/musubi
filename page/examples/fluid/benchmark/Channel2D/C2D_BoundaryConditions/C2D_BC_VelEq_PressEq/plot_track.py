## This is the user-script for plotting using gleaner tool.
mus_path = ['./']
omega = [1.9]

import sys
import os
import numpy as np
# Path to gleaner (Better use environment variable PYTHONPATH!)
if os.path.exists('/work/ws-tmp/gk779-SimTec-II/apes_executables/'):
  glrPath = '/work/ws-tmp/gk779-SimTec-II/apes_executables/'
elif os.path.exists( os.getenv('HOME')+'/apes/gleaner'):
  glrPath = os.getenv('HOME')+'/apes/gleaner'
else:
  print('Gleaner library not found')
  exit

# converrt
omega_scaled = [ int(100*i) for i in omega ]
omega = np.array(omega_scaled,dtype=int)

# Import gleaner module
sys.path.append(glrPath)
import gleaner

# Do not use X-server to create and save plot
import matplotlib
matplotlib.use('Agg')

## Import all required modules
import matplotlib.ticker as mtick
import matplotlib.pyplot as mplt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import logging

# font setting
from matplotlib import rc
font_size = 12
#font_family = 'serif'
#font_type = 'Times New Roman'
#rc('text',usetex=True)
#font = {'family':font_family,'%s'%font_family:font_type,'size':font_size}
#rc('font',**font)

#axis without scientific notation
y_formatter = mtick.ScalarFormatter(useOffset=False)
color = ['r','y','g','m','b']
markerstyle = ['x','+','v','o','D']

## -------------------------------------------------------------------------- ##
logging.basicConfig(level=logging.INFO)
logging.info('Started creating plots ...')

# data base filename
dbname = 'convergence.db'
# load database if exist else load tracking files and add to database
if os.path.isfile(dbname):
  logging.info('Processing data from existing database')
#  os.remove(dbname)
  import sqlite3
  sqlcon = sqlite3.connect(dbname)
else:
  logging.info('Processing data from tracking files')

  for iFile in range(len(mus_path)):
    omega_l = str(omega[iFile])
    file_path = mus_path[iFile]+'/tracking/'
    sqlcon = gleaner.tracking_to_db(fname = [file_path+'*velAlongHeight*.res'], \
                                    dbname=dbname, tabname='vel_spatial'+ omega_l )

    sqlcon = gleaner.tracking_to_db(fname = [file_path+'*pressAlongLength*.res'], \
                                    dbname=dbname, tabname='press_line'+ omega_l )

    sqlcon = gleaner.tracking_to_db(fname = [file_path+'*wss_spatial*.res'], \
                                    dbname=dbname, tabname='wss_spatial'+ omega_l )

    sqlcon = gleaner.tracking_to_db(fname = [file_path+'*velAlongLength*.res'], \
                                    dbname=dbname, tabname='vel_centerline'+ omega_l )

    sqlcon = gleaner.tracking_to_db(fname = [file_path+'*probe*.res'], \
                                    dbname=dbname, tabname='probe'+ omega_l )

    reference_path = mus_path[iFile]+'/reference/musubi_default_branch/'
    sqlcon = gleaner.tracking_to_db(fname = [reference_path+'*velAlongHeight*.res'], \
                                      dbname=dbname, tabname='ref_vel_spatial'+ omega_l )

    sqlcon = gleaner.tracking_to_db(fname = [reference_path+'*pressAlongLength*.res'], \
                                    dbname=dbname, tabname='ref_press_line'+ omega_l )

## -------------------------------------------------------------------------- ##
logging.info('Velocity Mag profile across the channel height:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Simulation result
get_data_for_cols = ['coordY','vel_mag_phy']
for iFile in range(len(mus_path)):
  omega_l = omega[iFile]
  [x, y] = gleaner.get_columns(sqlcon, tabname='vel_spatial'+str(omega_l), \
                               columns=get_data_for_cols)
  x, y = zip(*sorted(zip(x,y))) # sort of needed

  mplt.plot(x, y, '-', color = color[iFile],    \
            label = '$Omega$='+str(omega_l/100) )

  # Reference result
  [x, y] = gleaner.get_columns(sqlcon, tabname='ref_vel_spatial'+str(omega_l), \
                               columns=get_data_for_cols)
  x, y = zip(*sorted(zip(x,y))) # sort of needed
      
  mplt.plot(x, y, '-', color = color[iFile+1],    \
            label = '$Ref Omega$='+str(omega_l/100) )

# Plot x, y ... at certain time step
# Analytic result
omega_l = str(omega[-1])
get_data_for_cols = ['coordY','vel_analy']
[x, y] = gleaner.get_columns(sqlcon, tabname='vel_spatial'+omega_l, \
                             columns=get_data_for_cols)
x, y = zip(*sorted(zip(x,y))) # sort of needed
mplt.plot(x, y, '--', color = 'k', label = 'Analytic')

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
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Pressure drop across the channel length:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Simulation result
get_data_for_cols = ['coordX','pressure_phy']
for iFile in range(len(mus_path)):
  omega_l = omega[iFile]
  [x, y] = gleaner.get_columns(sqlcon, tabname='press_line'+str(omega_l), \
                               columns=get_data_for_cols)
  x, y = zip(*sorted(zip(x,y))) # sort of needed

  mplt.plot(x, y, '-', color = color[iFile],    \
            label = '$Omega$='+str(omega_l/100) )
# Reference result
  [x, y] = gleaner.get_columns(sqlcon, tabname='ref_press_line'+str(omega_l), \
                               columns=get_data_for_cols)
  x, y = zip(*sorted(zip(x,y))) # sort of needed
  
  mplt.plot(x, y, '-', color = color[iFile+1],    \
            label = '$Ref Omega$='+str(omega_l/100) )

# Plot x, y ... at certain time step
# Analytic result
omega_l = str(omega[-1])
get_data_for_cols = ['coordX','press_analy']
[x, y] = gleaner.get_columns(sqlcon, tabname='press_line'+omega_l, \
                             columns=get_data_for_cols)
x, y = zip(*sorted(zip(x,y))) # sort of needed
mplt.plot(x, y, '--', color = 'k', label = 'Analytic')

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
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Wall shear stress profile across the channel height:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# Simulation result
get_data_for_cols = ['coordY','wss_phy']
for iFile in range(len(mus_path)):
  omega_l = omega[iFile]
  [x, y] = gleaner.get_columns(sqlcon, tabname='wss_spatial'+str(omega_l), \
                               columns=get_data_for_cols)
  x, y = zip(*sorted(zip(x,y))) # sort of needed

  mplt.plot(x, y, '-', color = color[iFile],    \
            label = '$Omega$='+str(omega_l/100) )

# Plot x, y ... at certain time step
# Analytic result
omega_l = omega[-1]
get_data_for_cols = ['coordY','wss_analy']
[x, y] = gleaner.get_columns(sqlcon, tabname='wss_spatial'+str(omega_l), \
                             columns=get_data_for_cols)
x, y = zip(*sorted(zip(x,y))) # sort of needed
mplt.plot(x, y, '--', color = 'k', label = 'Analytic')

# plot setting
mplt.legend(loc=9, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('WSS ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(0.0,0.41)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('WSS_Profile.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Velocity X across the channel length:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('x (m)')
mplt.ylabel('Velocity X ($m/s$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# Plot x, y ... at certain time step
# Simulation result
get_data_for_cols = ['coordX','velocity_phy_01']
for iFile in range(len(mus_path)):
  omega_l = omega[iFile]
  [x, y] = gleaner.get_columns(sqlcon, tabname='vel_centerline'+str(omega_l), \
                               columns=get_data_for_cols)
  x, y = zip(*sorted(zip(x,y))) # sort of needed

  mplt.plot(x, y, '-', color = color[iFile],    \
            label = '$Omega$='+str(omega_l/100) )

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelLine_Profile.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Pressure over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('time (s)')
mplt.ylabel('Pressure ($Pa$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# Plot x, y ... at certain time step
# Simulation result
get_data_for_cols = ['time','pressure_phy']
for iFile in range(len(mus_path)):
  omega_l = omega[iFile]
  [x, y] = gleaner.get_columns(sqlcon, tabname='probe'+str(omega_l), \
                               columns=get_data_for_cols)
  x, y = zip(*sorted(zip(x,y))) # sort of needed

  mplt.plot(x, y, '-', color = color[iFile],    \
            label = '$Omega$='+str(omega_l/100) )


# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('PressureOverTime.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)

mplt.xlim(xmin=0, xmax=1)
mplt.savefig('PressureOverTime_zoomed.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##
logging.info('Velocity X over time:')
fig = mplt.figure()
ax = fig.add_subplot(111)
# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('time (s)')
mplt.ylabel('Velocity X ($m/s$)')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# Plot x, y ... at certain time step
# Simulation result
get_data_for_cols = ['time','velocity_phy_01']
for iFile in range(len(mus_path)):
  omega_l = omega[iFile]
  [x, y] = gleaner.get_columns(sqlcon, tabname='probe'+str(omega_l), \
                               columns=get_data_for_cols)
  x, y = zip(*sorted(zip(x,y))) # sort of needed

  mplt.plot(x, y, '-', color = color[iFile],    \
            label = '$Omega$='+str(omega_l/100) )


# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('VelocityXOverTime.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)

mplt.xlim(xmin=0, xmax=1)
mplt.savefig('VelocityXOverTime_zoomed.png', dpi=100, format='png', \
             bbox_inches="tight",interpolation=None)
## -------------------------------------------------------------------------- ##

logging.info('Plots created')
