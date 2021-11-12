#/bin/bash 
staList="staList.txt"

# The task to run for each parallel iteration. 
# Get network and station name, then either run or submit request to run inversion. 
cd .. 

staList='batch/staList.txt'

task(){
network=$(echo $line | head -n1 | awk '{print $1;}')
station=$(echo $line | head -n1 | awk '{print $2;}')
echo Running: $network $station
matlab -nodisplay -nosplash -nodesktop -r "network='$network'; station='$station'; RUN_one_station"
}

while read -r line; 
do
    task $line &
done < $staList

wait 
echo "Done starting jobs."