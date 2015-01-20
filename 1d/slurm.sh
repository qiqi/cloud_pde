#!/bin/bash
#SBATCH --job-name=diamond
#SBATCH --workdir=/master/home/qiqi/git/cloud_pde/1d
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
 
make run
