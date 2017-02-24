#!/bin/bash

# comparison plots, standard ranges, variable binning for low error in ratio plots
for i in {0..9}
do
echo "Running ./categorizeRun1 ${i} 0 and saving output into logs/${i}.log0 ..."
   ./categorizeRun1 $i 1 0 > "logs/${i}.log10" 2>&1 &
done
# plots needed for limit setting
    # plot dimuon mass in 110 to 160 with constant binning
   ./categorizeRun1 0 0 1 > "logs/0.log01" 2>&1 &
    # plot dimuon mass in 110 to 160 with variable binning
   ./categorizeRun1 0 1 1 > "logs/0.log11" 2>&1 &
