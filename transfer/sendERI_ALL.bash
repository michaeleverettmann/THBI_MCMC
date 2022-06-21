# This script should point to and send each of the scripts I use to transfer stuff from personal computer to HPC. 

computer=${1:-brunsvik@tong.eri.ucsb.edu} # This is the computer you are sending things to. 
# See the default value for an example: brunsvik@tong.eri.ucsb.edu
# It is used like this: $computer:path_to_send
echo Sending to: $computer

# Make sure everything in transferFolds is correct on your computer. 
# The folder structure (relative to ~) needs to be the same on the destination computer. 
transferFolds=\
'/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/transfer
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/CADMINEOS/transfer
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/PropMat/transfer
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/HV_ellipticity/transfer
/Users/brennanbrunsvik/MATLAB/transfer
/Users/brennanbrunsvik/Documents/repositories/data/models_seismic/transfer
/Users/brennanbrunsvik/Documents/repositories/Base_code/transfer
'

for fold in $transferFolds
do 
    cd $fold
    echo '**********'
    echo $fold 
    ./sendERI.bash $computer
done