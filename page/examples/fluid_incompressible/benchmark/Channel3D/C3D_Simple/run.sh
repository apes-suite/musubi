#!/bin/bash

# path to seeder executable
seeder_path=~/apes/seeder/build/seeder
# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

# Remove old directories and database
rm -rf tracking mesh restart *.db

# Create directories for Seeder and Musubi output
mkdir mesh tracking restart

# Run printParams.lua to print informations to screen
lua printParams.lua

# Run Seeder
$seeder_path seeder.lua

# Run Musubi
mpirun --oversubscribe -np 12 $musubi_path musubi.lua

# Create 3D Plots using Gleaner
python plot_track.py
# List the created plots
echo "List of created plots:"
ls *.png
