## This is the user-script for plotting using gleaner tool.
##
## Simulation parameters
import math
# Frictional Reynolds number
Re_tau = 1000
Re_bulk = (8.0/0.073)**(4.0/7.0)*Re_tau**(8.0/7.0)
half_height = 1.0
nhalf_height = 10
# Kinematic Viscosity
#nu = 5e-5
# mean bulk velocity m/s
#vel_bulk = Re_bulk * nu / (2*half_height)
vel_bulk = 1.0
nu = vel_bulk * (2*half_height) / Re_bulk
dx = half_height/nhalf_height
vel_tau = Re_tau*nu/half_height
#vel_bulk = 1.1287146295082
length = 2*math.pi*half_height
Ma = 0.1
cs_phy = vel_bulk / Ma
cs_L = math.sqrt(1.0/3.0)
dt = cs_L * dx / cs_phy
T_c = length/vel_bulk ## Turn over time
y_plus_mesh = dx*vel_tau/nu
y_plus_bc = dx/2.0*vel_tau/nu
samp_avg_iter = math.ceil(T_c/dt)
print('vel_tau ', vel_tau)

import numpy as np
# Load DNS data for 1000
vel_tau_dns = 0.05
y_plus_dns, u_plus_dns = np.loadtxt('DNS_data/DNS_Re1000.dat', usecols=(1,2), 
                                    dtype = 'float', comments='#', unpack=True)

y_dns, uu_plus_dns, vv_plus_dns, ww_plus_dns, uv_plus_dns = np.loadtxt(
       'DNS_data/DNS_Re1000_stress.dat', usecols=(0,2,3,4,5),
       dtype = 'float', comments='#', unpack=True)

# Function to convert array from gleaner to numpy array
def convert_array_to_nparray(x):
  x_new = np.array([])
  for i in x:
     if type(i) is tuple:
        x_new = np.append(x_new, sum(i)/len(i))
     else:
        x_new = np.append(x_new, i)
  return x_new

import sys
import os
# Path to gleaner (Better use environment variable PYTHONPATH!)
if os.path.exists('/scratch/ws/0/masi_ka-apes_repo/apes/gleaner/'):
  glrPath = '/scratch/ws/0/masi_ka-apes_repo/apes/gleaner'
elif os.path.exists( os.getenv('HOME')+'/apes/gleaner'): 
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
font_family = 'serif'
font_type = 'Times New Roman'
rc('text',usetex=True)
font = {'size':font_size}
rc('font',**font)

#axis without scientific notation
y_formatter = mtick.ScalarFormatter(useOffset=False)

## -------------------------------------------------------------------------- ##
logging.basicConfig(level=logging.INFO)
logging.info('Started creating plots ...') 

# data base filename
dbname = 'turbChannel.db'
# load database if exist else load tracking files and add to database
if os.path.isfile(dbname):
  logging.info('Processing data from existing database')
#  os.remove(dbname)
  import sqlite3
  sqlcon = sqlite3.connect(dbname)
else:
  logging.info('Processing data from tracking files')

  filename = 'tracking/Channel_half_p*_t3.770E+03.res'
  sqlcon = gleaner.tracking_to_db(fname = [filename], \
                                  dbname=dbname, tabname='yzPlane')
  fieldnames = gleaner.get_tracking_header(
                 'tracking/Channel_half_p00000_t3.770E+03.res')
  gleaner.spatial_reduction_in_db(sqlcon, tabname = 'yzPlane',
                                  columns = fieldnames,
                                  reduce_coord_column='coordY')
  sqlcon = gleaner.tracking_to_db(fname = 'tracking/Channel_bcTurbWall*.res', \
                                  dbname=dbname, tabname='bcTurbWall')
## -------------------------------------------------------------------------- ##

## -------------------------------------------------------------------------- ##
logging.info('Time average spatial averaged friction velocity at boundary:')
get_data_for_cols = ['time', 'bc_fric_velocity_phy_red']
[time, vel_tau_avg] = gleaner.get_columns(sqlcon, tabname='bcTurbWall', \
                          columns=get_data_for_cols, as_nparray=True)
vel_tau_sim = np.average(vel_tau_avg[-samp_avg_iter:])
print('vel_tau_sim: ', vel_tau_sim)
Re_tauN = vel_tau_sim*half_height/nu
print('Re_tauN:', Re_tauN)
vel_tau_sim = vel_tau_dns
#vel_tau_sim = 0.052

# Simulation result
logging.info('Loading all variables from yzPlane')
get_data_for_cols = ['coordY', 'vel_avg_01', 'vel_avg_02', 'vel_avg_03',
                     're_stress_avgX_01', 're_stress_avgY_02', 
                     're_stress_avgZ_03', 're_stress_avgX_02']
[y_sim, u_sim, v_sim, w_sim, uu_sim, vv_sim, ww_sim, uv_sim] = gleaner.get_columns(
        sqlcon, tabname='yzPlane_red', columns=get_data_for_cols, as_nparray=True)

y_plus_sim = y_sim*(vel_tau_sim/nu)
u_plus_sim = u_sim/vel_tau_sim
v_plus_sim = v_sim/vel_tau_sim
w_plus_sim = w_sim/vel_tau_sim
uu_plus_sim = uu_sim/vel_tau_sim**2 - u_plus_sim**2
vv_plus_sim = vv_sim/vel_tau_sim**2 - v_plus_sim**2
ww_plus_sim = ww_sim/vel_tau_sim**2 - w_plus_sim**2
uv_plus_sim = uv_sim/vel_tau_sim**2 - u_plus_sim*v_plus_sim

## -------------------------------------------------------------------------- ##
logging.info('Space-time averaged velocityX over height dim:')
fig = mplt.figure()
ax = fig.add_subplot(111)

mplt.plot(y_dns, u_plus_dns*vel_tau_dns, ls='-', color = 'k', label='DNS')
mplt.plot(y_sim, u_sim, ls=' ', marker='x', color = 'r', label='Musubi')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$y/H$')
mplt.ylabel('$\\langle{u}\\rangle$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('MeanVelocity_Re1000_dim.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Space-time averaged velocityX over height non-dim:')
fig = mplt.figure()
ax = fig.add_subplot(111)

mplt.plot(y_plus_dns, u_plus_dns, ls='-', color = 'k', label='DNS')
mplt.plot(y_plus_sim, u_plus_sim, ls=' ', marker='x', color = 'r', label='Musubi')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$y^+$')
mplt.ylabel('$u^+$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)
mplt.xlim(1,2000)
#mplt.ylim(0,25)
# Setting a logarithmic scale for x-axis
mplt.xscale('log')
mplt.xticks([1, 10, 100, 1000], [1, 10, 100, 1000])

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('MeanVelocity_Re1000_nondim.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Reynolds normal stress over height dim:')
fig = mplt.figure()
ax = fig.add_subplot(111)

mplt.plot(y_dns, uu_plus_dns*vel_tau_dns**2, ls='-', color = 'r', label='DNS ${u\'}^2$')
mplt.plot(y_dns, vv_plus_dns*vel_tau_dns**2, ls='-', color = 'b', label='DNS ${v\'}^2$')
mplt.plot(y_dns, ww_plus_dns*vel_tau_dns**2, ls='-', color = 'g', label='DNS ${w\'}^2$')
mplt.plot(y_sim, uu_plus_sim*vel_tau_sim**2, ls=' ', marker='x', color = 'r', label='Musubi ${u\'}^2$')
mplt.plot(y_sim, vv_plus_sim*vel_tau_sim**2, ls=' ', marker='x', color = 'b', label='Musubi ${v\'}^2$')
mplt.plot(y_sim, ww_plus_sim*vel_tau_sim**2, ls=' ', marker='x', color = 'g', label='Musubi ${w\'}^2$')

# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$y/H$')
#mplt.ylabel('$\\langle u\'u\' \\rangle$, $\\langle v\'v\' \\rangle$, $\\langle w\'w\' \\rangle$')
mplt.ylabel('$\\langle u\'u\' \\rangle$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('ReNormStress_Re1000_dim.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Reynolds normal stress over height non-dim:')
fig = mplt.figure()
ax = fig.add_subplot(111)

mplt.plot(y_dns, np.sqrt(uu_plus_dns), ls='-', color = 'r', label='DNS $u_{rms}$')
mplt.plot(y_dns, np.sqrt(vv_plus_dns), ls='-', color = 'b', label='DNS $v_{rms}$')
mplt.plot(y_dns, np.sqrt(ww_plus_dns), ls='-', color = 'g', label='DNS $w_{rms}$')
mplt.plot(y_sim, np.sqrt(uu_plus_sim), ls=' ', marker='x', color = 'r', label='Musubi $u_{rms}$')
mplt.plot(y_sim, np.sqrt(vv_plus_sim), ls=' ', marker='x', color = 'b', label='Musubi $v_{rms}$')
mplt.plot(y_sim, np.sqrt(ww_plus_sim), ls=' ', marker='x', color = 'g', label='Musubi $w_{rms}$')


# plot setting
mplt.legend(loc=1, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$y/H$')
#mplt.ylabel('$u_{rms}/u_\\tau$, $v_{rms}/u_\\tau$, $u_{rms}/u_\\tau$')
mplt.ylabel('$u_{rms}/u_\\tau$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('ReNormStress_Re1000_nondim.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Reynolds shear stress over height dim:')
fig = mplt.figure()
ax = fig.add_subplot(111)

mplt.plot(y_dns, uv_plus_dns*vel_tau_dns**2, ls='-', color = 'k', label='DNS')
mplt.plot(y_sim, uv_plus_sim*vel_tau_sim**2, ls=' ', marker='x', color = 'r', label='Musubi')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$y/H$')
mplt.ylabel('$\\langle u\'v\' \\rangle$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('ReShearStress_Re1000_dim.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##
logging.info('Reynolds shear stress over height nondim:')
fig = mplt.figure()
ax = fig.add_subplot(111)

mplt.plot(y_dns, uv_plus_dns, ls='-', color = 'k', label='DNS')
mplt.plot(y_sim, uv_plus_sim, ls=' ', marker='x', color = 'r', label='Musubi')

# plot setting
mplt.legend(loc=8, ncol=1,borderaxespad=0, \
            prop={'size':font_size}).get_frame().set_lw(0.0)
mplt.xlabel('$y/H$')
mplt.ylabel('$\\langle u\'v\' \\rangle/ {u_{\\tau}}^2$')
mplt.grid(True,which="major",ls="-")
ax.yaxis.set_major_formatter(y_formatter)

# save fig
figsize = [8,6]
fig = mplt.gcf()
fig.set_size_inches(figsize[0],figsize[1])
mplt.savefig('ReShearStress_Re1000_nondim.png', dpi=100, format='png', \
             bbox_inches="tight")
## -------------------------------------------------------------------------- ##

#mplt.show()
logging.info('Plots created')
