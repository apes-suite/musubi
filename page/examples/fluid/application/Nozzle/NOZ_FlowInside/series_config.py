## input restart file which replaces the string
## $!file!$ in the harvest_series.template.
## Can be single file or multiple files
files: tracking/*hvs*header*.lua

## harvester input file template
## the filename for read and folder name for output are
## defined by harvest_series.py so they are fixed to
## $!file!$ and $!out!$
template: harvest_series.template

## output folder where replaces the string 
## ${folder}$ in the harvester.template.
out: harvest/

## harvester exectuable
#harvester: /work/ws-tmp/gk779-SimTec-II/apes_executables/mus_harvesting
harvester: /home/gk779/apes/musubi/build/mus_harvesting

## command to run in parallel
run: mpirun -np 8
