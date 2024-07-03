#/bin/bash 
staList='batch/staList.txt'
stampList='batch/stampList.txt'
pid='batch/run_info/process_ID.txt'

cd .. 

# module load MatLab/R2021b || echo "Did not load Matlab. Are we not using Knot?"

# First, make sure we have the correct station files and paths on this computer. (otherwise paths might be generated on different computer than what you are using. )
echo "NOT RUNNING RUN_prep_data.m"
# matlab -nodisplay -nosplash -nodesktop -r "RUN_prep_data; exit"

echo '************** New Run **************' !>> $pid
echo $(date) !>> $pid # Put the date in process_ID file so we know what's what. Also, don't want to remove old PIDs, in case one keeps running and we have to shut it down. 
while read -r lineSta; # Loop over each station. 
do

    echo "sbatch batch/executeOneStation.bash '$lineSta'"
    sbatch batch/executeOneStation.bash $lineSta

    # echo "NOT USING SLURM batch/executeOneStation.bash '$lineSta'"
    # ./batch/executeOneStation.bash $lineSta

done < $staList

wait 
echo "Done starting jobs."
