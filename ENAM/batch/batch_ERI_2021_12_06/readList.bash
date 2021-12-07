#/bin/bash 
staList="staList.txt"

# The task to run for each parallel iteration. 
# Get network and station name, then either run or submit request to run inversion. 
task(){
network=$(echo $line | head -n1 | awk '{print $1;}')
station=$(echo $line | head -n1 | awk '{print $2;}')
echo Running: $network $station
}

while read -r line; 
do
    task $line
done < $staList