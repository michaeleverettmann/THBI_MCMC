#!/bin/bash

shopt -s nullglob

files=( *.run_HVker )
if (( ${#files[@]} )); then
    

for i in *.run_HVker
 do 
 j=`echo $i | sed s/.run_HVker//g`
 echo $j
 rm $j*
done


fi
