#!/bin/bash

#for i in {0..63}
#do
#  #./wrmsr -p $i 0x1a4 15 # disable
#  ./wrmsr -p $i 0x1a4 0 # enable prefetching
#done


echo "Reading back"

for i in {0..63}
do
  ./rdmsr -p $i 0x1a4
done
