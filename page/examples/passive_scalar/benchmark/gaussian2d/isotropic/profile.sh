#! /bin/bash

# create a argument configuration named arg_given.lua
echo "tau = 5" > arg_given.lua
echo "u_field = 0." >> arg_given.lua
echo "order = 'first'" >> arg_given.lua
echo "collision = 'bgk'" >> arg_given.lua
echo "sigma0 = 40" >> arg_given.lua

# meshing and simulation
mkdir -p mesh tracking restart
seeder  # Attention: execution command might differ 
../../../../../build/musubi  # Attention: execution command might differ 

# plot with Gnuplot
# please check "tutorial/1_musubi_config/tut_1_mus_config.md" for more information
# reference file can be found in "./reference"
gnuplot -p -e "plot \"tracking/simulation_spc1_p00000_t200.000E+00.res\" \
using 1:4 title \"simulation\" with points, \
\"tracking/simulation_spc1_p00000_t200.000E+00.res\" \
using 1:5 title \"analytical\" with lines "
