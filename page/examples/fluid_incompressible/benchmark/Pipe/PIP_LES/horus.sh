#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --mail-type=END
#SBATCH --mail-user=kannan.masilamani@uni-siegen.de
#SBATCH --exclusive
#SBATCH --nodes 10
#SBATCH --dependency=singleton
#SBATCH --ntasks-per-node=12

# Load modules
module load PrgEnv/intel-openmpi

# Change to test case folder path
cd /home/gk779/apes/musubi/examples/fluid_incompressible/benchmark/Pipe/PIP_LES 
mpiexec /home/gk779/apes/musubi/build/musubi musubi.lua
