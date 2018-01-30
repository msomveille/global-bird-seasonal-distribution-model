#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=99:00:00
#SBATCH --job-name=juliaModel

module purge
module load julia
julia migraModel_parallel3.jl
