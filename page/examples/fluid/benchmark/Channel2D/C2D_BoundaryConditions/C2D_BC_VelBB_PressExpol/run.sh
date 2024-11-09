#!/bin/bash
rm -rf tracking mesh restart channel2D.db

mkdir mesh tracking restart

~/apes/seeder/build/seeder seeder.lua

mpirun --oversubscribe -np 12 ~/apes/musubi/build/musubi musubi.lua

python plot_track.py
