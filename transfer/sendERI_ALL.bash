# This script should point to and run each of the scripts I use to transfer stuff from personal computer to HPC. 
transferFolds=\
'/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/transfer
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/CADMINEOS/transfer
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/PropMat/transfer
/Users/brennanbrunsvik/MATLAB/transfer
/Users/brennanbrunsvik/Documents/repositories/data/models_seismic/transfer
'

for fold in $transferFolds
do 
    cd $fold
    ./sendERI.bash
done