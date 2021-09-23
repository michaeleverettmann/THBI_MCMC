#!/bin/bash 
# This script loops through several folders, compiles the code, and tests the compilation success. 
# It checks which computer you are on (by looking at where your home directory is). 
# Then you might have different locations where code is stored depending on the different computers. 
homeDir=~ 

# TODO much easier if we just use ~/stuff. I didn't realize that my path structure is the same on the mac and linux computers after home. 
if [ $homeDir = '/Users/brennanbrunsvik' ]; then # My mac computer. Enter relevant folders here. 
    transferFolds=\
'/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/transfer
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/CADMINEOS/transfer
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/PropMat/transfer
/Users/brennanbrunsvik/MATLAB/transfer'
elif [ $homeDir = '/home/brunsvik' ]; then # ERI computers. Enter relevant folders here. 
    transferFolds=\
'/home/brunsvik/Documents/UCSB/ENAM/THBI_ENAM/transfer
/home/brunsvik/Documents/repositories/Peoples_codes/CADMINEOS/transfer
/home/brunsvik/Documents/repositories/Peoples_codes/PropMat/transfer
/home/brunsvik/MATLAB/transfer'
fi 

for fold in $transferFolds
do 
    cd $fold
    ./compileTest.bash
done