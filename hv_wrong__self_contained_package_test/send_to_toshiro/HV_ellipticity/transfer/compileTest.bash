#!/bin/bash

echo " ************* "
echo " HV code compilation "

cd ..

cd src 
make 

echo " *************** "
echo " Finished compiling HV code. "

cd .. 
echo "HV code test: this will be ran directly from the THBI folder (something like  THBI/matlab_to_hv_kernel/example/test_compare_to_MINEOS."