#/bin/bash 
# Needs to be ran from the main THBI_ENAM folder. 
# Eventually replace the destination with a variable. 
rsync -ahv --dry-run \
--exclude-from='transfer/exclude_THBI.txt' \
/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ \
brunsvik@tong.eri.ucsb.edu:~/Documents/UCSB/ENAM/THBI_ENAM