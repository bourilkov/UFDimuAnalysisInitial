#!/bin/bash

# comparison plots, standard ranges, variable binning for low error in ratio plots
for i in {0..9}
do
echo "Running ./categorizeRun1 ${i} and saving output into logs/sig${i}.log ..."
   ./significanceRun1 $i > "logs/sig${i}.log" 2>&1 &
done
