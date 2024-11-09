#!/bin/bash

# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

# Remove old directories
rm -rf tracking restart vtkfiles *.db

# Create directories for Seeder and Musubi output
mkdir tracking restart vtkfiles

# Run Musubi
mpirun --oversubscribe -np 12 $musubi_path musubi.lua

# Create 2D Plots using Gleaner
python plot_track.py

### kinetic energy and dissipation rate can be used to plot using gnuplot
### Compute disspation rate from kinetic energy
##export ke_file='tracking/TGV_kE_all_p00000.res'
##export dr_file='tracking/TGV_dr_all.res'
##lua calc_diss.lua
### Plot kinetic energy and disspiation rate
##export Re=1600
##export ke_pic='kineticEnergy.pdf'
##gnuplot plot_ke.gpl
##
##export dr_pic='dissipationRate.pdf'
##gnuplot plot_dr.gpl

# List the created plots
echo "List of created plots:"
ls *.png
