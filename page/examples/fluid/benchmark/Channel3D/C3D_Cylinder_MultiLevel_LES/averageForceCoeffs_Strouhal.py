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

# font setting
from matplotlib import rc
font_size = 12

#axis without scientific notation
y_formatter = mtick.ScalarFormatter(useOffset=False)

import numpy as np
import math

def to_str(var):
	return str(list(np.reshape(np.asarray(var), (1, np.size(var)))[0]))[1:-1]

import logging
logging.basicConfig(level=logging.INFO)
## -------------------------------------------------------------------------- ##
# data base filename
os.remove('channel.db')
dbname = 'channel.db'
# load database if exist else load tracking files and add to database
#if os.path.isfile(dbname):
#  logging.info('Processing data from existing database')
#  os.remove(dbname)
#  import sqlite3
# sqlcon = sqlite3.connect(dbname)
#else:
logging.info('Processing data from tracking files')

# Load tracking output with label cyl_force - Lift and drag coefficient
sqlcon = gleaner.tracking_to_db(fname = ['tracking/*_force*.res'], \
                                  dbname=dbname, tabname='cLcDCoeff')

sqlcon = gleaner.tracking_to_db(fname = ['tracking/*_probe*.res'], \
                                  dbname=dbname, tabname='probe')
## -------------------------------------------------------------------------- ##

# Evaluate Cl and Cd average
get_data_for_cols = ['time','coeff_red_01']
[t_Cd, Cd] = gleaner.get_columns(sqlcon, tabname='cLcDCoeff', \
                                columns=get_data_for_cols)

get_data_for_cols = ['time','coeff_red_02']
[t_Cl, Cl] = gleaner.get_columns(sqlcon, tabname='cLcDCoeff', \
                                columns=get_data_for_cols)

Cl_avg = np.array(Cl)
Cd_avg = np.array(Cd)
Cl_avg = np.average(Cl[len(Cl)//2:-1])
Cd_avg = np.average(Cd[len(Cl)//2:-1])
Cl_max = np.amax(Cl[len(Cl)//2:-1])
logging.info('Cl average = ' + to_str(Cl_avg))
logging.info('Cl max     = ' + to_str(Cl_max))
logging.info('Cd average = ' + to_str(Cd_avg))


# Strouhal number
from scipy import signal
from scipy.fft import fft, fftfreq, fftshift
get_data_for_cols = ['time','pressure_phy','velocity_phy_01','velocity_phy_02','velocity_phy_03']
[t_V, p, u, v, w] = gleaner.get_columns(sqlcon, tabname='probe', \
                                columns=get_data_for_cols)

fig = mplt.figure()
ax = fig.add_subplot(111)
mplt.plot(t_V, v, '-', color = 'k')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$s$')
mplt.ylabel('$V_y$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
#mplt.xlim(0.12,0.15)
#mplt.ylim(3.8,4.8)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('U_y.png', dpi=100, format='png', \
             bbox_inches="tight")

v_array = np.array(Cl)
v_array = np.concatenate(v_array)
vF = np.correlate(v_array,v_array,'full')
norm_fac = np.sum(v_array**2)
vF = vF/norm_fac
tF = np.array(t_Cl)
tF = tF.T
tF = tF - tF[0][0]
tF_flipped = np.fliplr(-tF)
tF_flipped = np.concatenate(tF_flipped)
tF_shift = np.concatenate([tF_flipped,tF[0][1:]])

fig = mplt.figure()
ax = fig.add_subplot(111)
mplt.plot(tF_shift, vF, '-', color = 'k')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('time')
mplt.ylabel('power')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
#mplt.xlim(0.12,0.15)
#mplt.ylim(3.8,4.8)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('XCorr.png', dpi=100, format='png', \
             bbox_inches="tight")

peaks, _ = signal.find_peaks(vF) #, height=0)
tF_peaks = tF_shift[peaks]
t = np.diff(tF_peaks)
period = np.mean(t)
freq = 1/period

radius = 0.00145275
Vel = 20.0772
L = 2 * radius
St = freq * L / Vel
logging.info('Strouhal (v) = ' + to_str(St))
