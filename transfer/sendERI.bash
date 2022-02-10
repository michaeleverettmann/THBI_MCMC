#/bin/bash 
# Needs to be ran from the main THBI_ENAM folder. 

computer=${1:-brunsvik@bellows.eri.ucsb.edu} # First argument is the computer you will send to. 
# excludeExtra=${2:---exclude={"**.mat","**.jpg","**.png","**.jpg"}} # Second argument is just more things to pass to rsync. 


rsync -ahv \
--exclude-from='exclude_THBI.txt' \
--exclude={'**.mat','**.jpg','**.png','**.jpg','**invState'} \
../ \
$computer:~/Documents/UCSB/ENAM/THBI_ENAM

 # Important! Comment out the .mat, .jpg, etc. stuff if you are sending to a new computer. 
#  --exclude={'**.mat','**.jpg','**.png','**.jpg','**invState'} \
