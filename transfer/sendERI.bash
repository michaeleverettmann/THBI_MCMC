#/bin/bash 
# Needs to be ran from the main THBI_ENAM folder. 
# Eventually replace the destination with a variable. 
rsync -ahv \
--exclude-from='exclude_THBI.txt' \
../ \
brunsvik@tong.eri.ucsb.edu:~/Documents/UCSB/ENAM/THBI_ENAM