transferFolds=\
'/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/transfer
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/CADMINEOS/transfer
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/PropMat/transfer'

for fold in $transferFolds
do 
    cd $fold
    ./compileTest.bash
done