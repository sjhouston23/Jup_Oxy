#!/bin/bash                                                                                          
#SBATCH --partition=sixhour
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --time=06:00:00
#SBATCH --array=1-200

R=$RANDOM
./oxy.x $R
echo $R >> ../scratch/Jup_Oxy/Output/Trials.dat
exit 0
