#!/bin/bash
#PBS -l nodes=1:ppn12 # Run all processes on a single node	
#PBS -l --walltime=00:01:00   # bb2021.11.17 Might have to refine this, just to make sure a failed run will be terminated at a decent time.             # Time limit hrs:min:sec
#PBS -o batch/run_info/pbs.%j.out
export NCORES=12

# First argument is station.network. 
# Script executes RUN_one_station.m on that station. 

# line=$1
# line=$2
# network=$(echo $line | head -n1 | awk '{print $1;}')
# station=$(echo $line | head -n1 | awk '{print $2;}')
network=$1
station=$2
echo Running: $network $station
pwd
# #nohup 
# mpirun -np $NCORES -machinefile matlab -nodisplay -nosplash -nodesktop -r "network_manual='$network'; station_manual='$station'; RUN_one_station" # &> batch/run_info/nohup$network.$station.out 
# }
