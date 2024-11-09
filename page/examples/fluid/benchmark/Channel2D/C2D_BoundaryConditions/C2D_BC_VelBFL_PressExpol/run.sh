#!/bin/bash
rm -rf tracking vtkfiles mesh restart

mkdir mesh tracking vtkfiles restart

~/apes/seeder/build/seeder seeder.lua

mpirun --oversubscribe -np 13 ~/apes/musubi/build/musubi musubi.lua

python plot_track.py
