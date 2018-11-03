#!/bin/bash	
#SBATCH	-N	1	
#SBATCH	-p	RM	
#SBATCH	--ntasks-per-node	28	
#SBATCH	-t	00:05:00	
# echo commands to stdout	
set	-x	
cd /home/yipengm/eecs587/EECS587/yipengm_hw4
# run OpenMP program	
export	OMP_NUM_THREADS=4
./hw4 1 100 0.000001 12
