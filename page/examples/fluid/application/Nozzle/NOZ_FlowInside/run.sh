#!/bin/bash

# path to seeder executable
seeder_path=~/apes/seeder/build/seeder
# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

# Remove old directories
rm -rf tracking mesh restart nozzle.db

# Create directories for Seeder and Musubi output
mkdir mesh tracking restart 

# Run Seeder
$seeder_path seeder.lua

# Run Musubi
mpirun -np 12 $musubi_path musubi.lua

# Create 2D Plots using Gleaner
python plot.py

# Convert harvest output to vtk output for animation
python ~/apes/musubi/treelm/peons/harvest_series.py --config series_config.py
