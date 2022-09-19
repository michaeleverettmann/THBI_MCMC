#/bin/bash 
staList='batch/staList.txt'
stampList='batch/stampList.txt'
pid='batch/run_info/process_ID.txt'

cd .. 

# module load MatLab/R2021b || echo "Did not load Matlab. Are we not using Knot?"

# First, make sure we have the correct station files and paths on this computer. (otherwise paths might be generated on different computer than what you are using. )
echo "NOT RUNNING RUN_prep_data.m"
# matlab -nodisplay -nosplash -nodesktop -r "RUN_prep_data; exit"

# # The task to run for each parallel iteration. 
# # Get network and station name, then either run or submit request to run inversion. 
# task(){
# network=$(echo $line | head -n1 | awk '{print $1;}')
# station=$(echo $line | head -n1 | awk '{print $2;}')
# echo Running: $network $station
# pwd
# nohup matlab -nodisplay -nosplash -nodesktop -r "network_manual='$network'; station_manual='$station'; RUN_one_station" &> batch/run_info/nohup$network.$station.out 
# }

echo '************** New Run **************' !>> $pid
echo $(date) !>> $pid # Put the date in process_ID file so we know what's what. Also, don't want to remove old PIDs, in case one keeps running and we have to shut it down. 
while read -r lineSta; # Loop over each station. 
do
    # # Run "task" on network, station, and echo the process_ID to $pid file. 
    # ftest="~/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv/${lineSta//[ ]/_}_dat1/standard/allmodels_perchain_orig.mat"
    # # echo "Check if file there. $ftest"
    # if test -f "$ftest"; then 
    #     echo "$ftest exists."
    # else
    #     echo "sbatch batch/executeOneStation.bash '$lineSta'"
    #     # sbatch batch/executeOneStation.bash $lineSta
    #     # sleep 1.0
    # fi



    echo "sbatch batch/executeOneStation.bash '$lineSta'"
    sbatch batch/executeOneStation.bash $lineSta



done < $staList

wait 
echo "Done starting jobs."
