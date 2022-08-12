#!/bin/bash
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Number of processes
#SBATCH --cpus-per-task=12 # bb2021.11.17 I think there is one task for each par for iteration, taking one full CPU. Might have to adjust this. 
#SBATCH --time=10:00:00   # bb2021.11.17 Might have to refine this, just to make sure a failed run will be terminated at a decent time.             # Time limit hrs:min:sec
##SBATCH --nice=1000 # bb2021.11.17 Be nice to the other researchers and let their scripts take some precedence. Not in use right now, will be once I really start things. 
#SBATCH --output=batch/run_info/slurm.%j.out
##SBATCH --oversubscribe

# Below things might not be needed. 
##SBATCH --mem=2G #bb2021.11.17 default option on Bellows is unlimited (try scontrol show config). So for Bellows, we do not need to specify this (I think)                     # Total memory limit
##SBATCH --output=multiprocess_%j.log # Standard output and error log
##SBATCH --job-name=parallel_job_test # Job name
##SBATCH --mail-type=END            # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=brennanbrunsvik@ucsb.edu  #bb2021.11.17 This doesn't seem to be working.   # Where to send mail	

module load MatLab/R2021b || echo "Did not load Matlab. Are we not using Knot?"

# Info above only applies if running with sbatch. Script still runs if executed from bash script. 
# First argument is station.network. 
# Script executes RUN_one_station.m on that station. 


# #nohup 
matlab -nodisplay -nosplash -nodesktop -r "RUN_example" # &> batch/run_info/nohup$network.$station.out 
# }
