#!/bin/bash

for i in {0..9}
do
echo "Running ./exploreDY ${i} 0 and saving output into logs/${i}.log0 ..."
echo "Running ./exploreDY ${i} 1 and saving output into logs/${i}.log1 ..."
   ./exploreDY $i 0 > "logs/${i}.log0" 2>&1 &
   ./exploreDY $i 1 > "logs/${i}.log1" 2>&1 &
done
