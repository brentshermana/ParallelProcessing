#!/bin/bash
for ((CORES=3;CORES<16;CORES++))
do
    #echo $CORES
    ./prs 16000               .2               16                 3 $CORES
done
