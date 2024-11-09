#!/bin/bash
rm -rf tracking mesh restart *.db

mkdir mesh tracking restart

~/apes/seeder/build/seeder seeder.lua

mpirun --oversubscribe -np 4 ~/apes/musubi/build/musubi musubi.lua

python plot_track.py
