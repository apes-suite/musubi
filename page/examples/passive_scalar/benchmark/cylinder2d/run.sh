#! /bin/bash

# build necessary folders and meshing
mkdir -p mesh tracking restart
seeder  # Attention: execution command might differ 

# do simulation
../../../../../build/musubi

# plot
# reference file can be found in "./reference"
gnuplot -p -e "plot \"tracking/simulation_spc1_p00000_t30.000E+03.res\" \
using 1:4 title \"computed\" with points,
\"tracking/simulation_spc1_p00000_t30.000E+03.res\"  \
using 1:5 title \"analytical\" with lines " 