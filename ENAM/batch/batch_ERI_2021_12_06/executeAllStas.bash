#/bin/bash 
staList='batch/staList.txt'
pid='batch/run_info/process_ID.txt'

cd .. 

# First, make sure we have the correct station files and paths on this computer. (otherwise paths might be generated on different computer than what you are using. )
echo Temporarilly not running RUN_prep_data.m
matlab -nodisplay -nosplash -nodesktop -r "RUN_prep_data; exit"

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
while read -r line; # Loop over each station. 
do
    # Run "task" on network, station, and echo the process_ID to $pid file. 
#     ./batch/executeOneStation.bash $line & echo $! >> $pid
    sbatch batch/executeOneStation.bash $line
    # $(echo $line | head -n1 | awk '{print $1;}').$(echo $line | head -n1 | awk '{print $2;}')
done < $staList

wait 
echo "Done starting jobs."
