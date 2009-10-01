#!/bin/sh

LocalVersion=/net/hisrv0001/home/yetkin/CMSSW_3_3_0_pre5

echo "Hello World - script working"
echo "Arguments" " 1="$1 " 2="$2 " 3="$3

echo "Setting up hacked CMSSW pre-release"
mit_workdir=`pwd`
export SCRAM_ARCH=slc4_ia32_gcc345
cd $LocalVersion/src
source /app/cms-soft/cmsset_default.sh
eval `scramv1 ru -sh`
cd $mit_workdir

echo "CMSSW Environment"
echo "CMSSW_BASE=" $CMSSW_BASE
echo "CMSSW_VERSION=" $CMSSW_VERSION
echo "CMSSW_BASE=" $CMSSW_BASE
echo "CMSSW_RELEASE_BASE=" $CMSSW_RELEASE_BASE
echo "CMSSW_SEARCH_PATH=" $CMSSW_SEARCH_PATH

echo "PY files in current directory are : "

ls *.py

echo "Using PARAMETER SET pset.py : "
cat pset.py

echo "Running the CMSSW job : pset.py"

#NJob=$1
cmsRun -j $RUNTIME_AREA/crab_fjr_$NJob.xml -p pset.py






