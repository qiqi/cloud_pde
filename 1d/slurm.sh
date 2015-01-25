#!/bin/bash
#SBATCH --job-name=kuramotoDiamondTest
#SBATCH --workdir=/master/home/qiqi/git/cloud_pde/1d
#SBATCH --output=kuramotoDiamondTest.out
#SBATCH --error=kuramotoDiamondTest.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
 
make run
