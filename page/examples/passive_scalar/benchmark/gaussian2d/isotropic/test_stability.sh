#! /bin/bash

# remove stability.res if it already exists
rm -f stability.res

# build necessary folders and meshing
mkdir -p mesh tracking restart
seeder  # Attention: execution command might differ 

# do simulation with different relaxation time 'tau'
# then compute relative error of computed diffusion factor
for tau in $(seq 0.5001 0.0002 0.505)
do
    echo "tau = $tau" > arg_given.lua
    echo "u_field = 0." >> arg_given.lua
    echo "order = 'first'" >> arg_given.lua
    echo "collision = 'bgk'" >> arg_given.lua

    ../../../../../build/musubi
    lua calculateD.lua >> stability.res
done

for tau in $(seq 0.507 0.002 0.8)
do
    echo "tau = $tau" > arg_given.lua
    echo "u_field = 0." >> arg_given.lua
    echo "order = 'first'" >> arg_given.lua
    echo "collision = 'bgk'" >> arg_given.lua

    ../../../../../build/musubi
    lua calculateD.lua >> stability.res
done

for tau in $(seq 1 1 5)
do
    echo "tau = $tau" > arg_given.lua
    echo "u_field = 0." >> arg_given.lua
    echo "order = 'first'" >> arg_given.lua
    echo "collision = 'bgk'" >> arg_given.lua

    ../../../../../build/musubi
    lua calculateD.lua >> stability.res
done

# plot with Gnuplot
# please check "tutorial/1_musubi_config/tut_1_mus_config.md" for more information
gnuplot -p -e "plot \"stability.res\" \
using 2:6 title \"D_err\" with linespoints " 