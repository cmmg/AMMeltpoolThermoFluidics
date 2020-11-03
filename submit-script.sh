#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=int # default "univ", if not specified
#SBATCH --time=0-00:29:00 # run time in days-hh:mm:ss
#SBATCH --nodes=1# require 2 nodes
#SBATCH --ntasks-per-node=16            # (by default, "ntasks"="cpus")
#SBATCH --mem-per-cpu=4000# RAM per CPU core, in MB (default 4 GB/core)
#SBATCH --job-name="cmake"
#SBATCH --job-name="SIERBM"
#SBATCH --error=file%j.err
#SBATCH --output=file%j.out
#SBATCH --mail-user=kbhagat2@wisc.edu
##SBATCH --mail-type=ALL
#Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output

#Now list your executable command (or a string of them).
# Example for non-SLURM-compiled code:
source ~/.bashrc
cd /home/kbhagat2/workspace/RBM/StagMultImpExp
rm -rf CMakeCache.txt cmake_install.cmake Makefile CMakeFiles 
cmake .
cmake -DCMAKE_BUILD_TYPE=Release . 
make
mpirun -np 10 ./main
