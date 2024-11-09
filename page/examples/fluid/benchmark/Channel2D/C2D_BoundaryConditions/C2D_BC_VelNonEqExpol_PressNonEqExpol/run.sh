#!/bin/bash
rm -rf tracking vtkfiles mesh restart *.db

mkdir mesh tracking vtkfiles restart

~/apes/seeder/build/seeder seeder.lua

mpirun -np 5 ~/apes/musubi/build/musubi musubi.lua

python plot_track.py
