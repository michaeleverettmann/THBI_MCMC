This folder should contain various scripts to transfer code from a local machine to HPC. I intend also to include code to quickly run compilers for the various Fortran codes and C code and test that compilation was successful. bb2021.09.15

______________ Jupyter 
Sometimes, it's easiest to set up a Jupyter server on ERI/other computer and ssh into that one time on our remote computer. Some rough instructions: 

Open Jupyter on the remote machines terminal. 
jupyter lab --no-browser --port=8080

Connect local computer to the 8080 port at remote computer using local machines terminal (change computer name to whatever you need... e.g.) 
ssh -N -L 8080:localhost:8080 brunsvik@anvil.eri.ucsb.edu
ssh -N -L 8080:localhost:8080 brunsvik@tong.eri.ucsb.edu
ssh -N -L 8080:localhost:8080 brunsvik@bellows.eri.ucsb.edu



________________ Running scripts
Paste one of the returned links from terminal into my local machines browser. Should be good to go!

One good way to run matlab scripts over ssh
rm nohup.out; nohup matlab -nodisplay -nodesktop -nosplash -batch "run('RUN_ALL_STAS.m')"



________________ Checking on code 

See how much cpu in total you are using on a machine. 
top -b -n 1 -u brunsvik| awk 'NR>7 { sum += $9; } END { print sum; }' 