#/bin/bash 
# Needs to be ran from the main THBI_ENAM folder. 

computer=${1:-brunsvik@tong.eri.ucsb.edu} # First argument is the computer you will send to. 

rsync -ahv \
--exclude-from='exclude_THBI.txt' \
../ \
$computer:~/Documents/UCSB/ENAM/THBI_ENAM