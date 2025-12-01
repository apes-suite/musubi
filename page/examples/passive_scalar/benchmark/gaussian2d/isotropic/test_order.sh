#! /bin/bash

# remove stability.res if it already exists
rm -f order.res

# select relaxation time 'tau'
#   in 0.501 0.503 0.505 0.51 0.55 0.8 2 5
tau=5

# build necessary folders and meshing
mkdir -p mesh tracking restart
seeder  # Attention: execution command might differ 

# do simulations with different advection velocities
#   and order of equilibrium function in collison scheme
# compare the errors of diffusion factors between
#   first and second orders
for u in 0.001 0.005 0.01 0.05 0.1
do
    for order in 'first' 'second'
    do
        echo "tau = $tau" > arg_given.lua
        echo "u_field = $u" >> arg_given.lua
        echo "order = '$order'" >> arg_given.lua
        echo "collision = 'bgk'" >> arg_given.lua

        ../../../../../build/musubi
        lua calculateD.lua >> order.res
    done
done

# plot the error of diffusion factor between 
#   first and second order
grep '1st' order.res > 1st.res
grep '2nd' order.res > 2nd.res
gnuplot -p -e "plot \"1st.res\" \
using 3:7 title \"1st_order\" with linespoints,
\"2nd.res\" using 3:7 title \"2nd_order\" \
with linespoints " 
rm 1st.res 2nd.res