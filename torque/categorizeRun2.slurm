#!/bin/bash
#SBATCH --job-name=categorizeRun2
#SBATCH -o categorizeRun2_%A_%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=acarnes@phys.ufl.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=10:00:00
#SBATCH --array=0-23

#Set the compiler architecture and CMSSW
export SCRAM_ARCH=slc6_amd64_gcc491
export CMS_PATH=/cvmfs/cms.cern.ch
source ${CMS_PATH}/cmsset_default.sh
export CVSROOT=:pserver:anonymous@cmssw.cvs.cern.ch:/local/reps/CMSSW

# source CMSSW libs to get root
cd /home/puno/cmsswinit/CMSSW_7_5_8/src/
eval `scram runtime -sh`
cd /home/puno/h2mumu/UFDimuAnalysis_v2/bin/

# run the executable
date
hostname
pwd
echo "JOB ID: ${SLURM_ARRAY_JOB_ID}"
echo "ARRAY ID: ${SLURM_ARRAY_TASK_ID}"
echo ""
#./categorizeRun2 varToPlot varRebin? plot110-160? nPartitions partition#
./categorizeRun2 ${SLURM_ARRAY_TASK_ID} 0 0 1 0
date
