#/bin/bash 
rsync -amh \
--exclude-from='exclude_CADMINEOS.txt' \
/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/CADMINEOS/ \
brunsvik@tong.eri.ucsb.edu:~/repositories/Peoples_codes/CADMINEOS



### Some examples below. 
# #!/bin/bash
# logfile=/Users/brennanbrunsvik/Documents/for_backing_up/Records/rsync_cubase_ext_drive_to_mybook/$(date +%Y%m%d).log
# rm -f $logfile



# #!/bin/bash
# logfile=/Users/brennanbrunsvik/Documents/for_backing_up/Records/rsync_mac_to_mybook/$(date +%Y%m%d).log
# rm -f $logfile

# rsync -amh \
# --exclude-from='/Users/brennanbrunsvik/Documents/for_backing_up/rsync_exclude.txt' \
# --stats \
# /Users/brennanbrunsvik/Documents/ /Volumes/mybook/mac_rsync \
# >> $logfile

# # echo 'Backed up today I think' >> $logfile
