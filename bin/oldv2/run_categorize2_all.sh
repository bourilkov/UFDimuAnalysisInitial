#!/bin/bash

# comparison plots, standard ranges, variable binning for low error in ratio plots
for i in {0..23}
do
echo "Running ./categorizeRun2 ${i} 0 0 and saving output into logs/${i}.log00 ..."
   ./categorizeRun2 $i 0 0 > "logs/${i}.log00" 2>&1 &
done
# plots needed for limit setting
echo "Running ./categorizeRun2 0 0 1 and saving output into logs/0.log01 ..."
    # plot dimuon mass in 110 to 160 with constant binning
   ./categorizeRun2 0 0 1 > "logs/0.log01" 2>&1 &
echo "Running ./categorizeRun2 0 1 1 and saving output into logs/0.log11 ..."
    # plot dimuon mass in 110 to 160 with variable binning
   ./categorizeRun2 0 1 1 > "logs/0.log11" 2>&1 &
