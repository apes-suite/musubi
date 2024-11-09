#!/bin/bash

# path to seeder executable
seeder_path=~/apes/seeder/build/seeder
# path to musubi executable
musubi_path=~/apes/musubi/build/musubi

# Remove old directories
rm -rf tracking mesh restart

# Create directories for Seeder and Musubi output
mkdir mesh tracking restart 

# Run Seeder
$seeder_path seeder.lua

# Run Musubi
mpirun --oversubscribe -np 12 $musubi_path musubi.lua
