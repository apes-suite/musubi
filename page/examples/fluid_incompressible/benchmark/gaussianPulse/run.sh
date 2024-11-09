#!/bin/bash

# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

# Remove old directories
rm -rf tracking restart

# Create directories for Seeder and Musubi output
mkdir tracking restart

# Run Musubi
mpirun --oversubscribe -np 2 $musubi_path musubi.lua

# Remove database if existent
rm *.db
# Create 3D Plots using Gleaner
python plot_track.py
# List the created plots
echo "List of created plots:"
ls *.png
