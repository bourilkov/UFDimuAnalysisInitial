#!/bin/bash

# -999 muon comparison
# ./fewzCategorize 0 1 1 1 1 0> "logs/log0-11110" 2>&1 & #get this in the smearing set, no need to duplicate it
# dimu mass, recomuCuts?, recomuPlot?, recojetCuts?, recojetsPlot?, cut no gens? 
./fewzCategorize 0 1 1 1 1 1 > "logs/fewzcat-log0-11111" 2>&1 &

# Smearing
./fewzCategorize 0 0 0 1 1 1 > "logs/fewzcat-log0-00111" 2>&1 &
./fewzCategorize 0 1 1 1 1 0 > "logs/fewzcat-log0-11110" 2>&1 &
