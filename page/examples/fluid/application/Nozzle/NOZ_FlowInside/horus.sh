#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --mail-type=END
#SBATCH --exclusive
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=12

module load PrgEnv/intel-openmpi

cd /work/ws-tmp/gk779-scratch/lehre/SimTecII/nozzle

date
mpiexec -n 12 /work/ws-tmp/gk779-SimTec-II/apes_executables/musubi musubi.lua
date
