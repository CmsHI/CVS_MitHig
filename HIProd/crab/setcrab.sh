#!/bin/sh

CMSSW_BASE=$1
cd $CMSSW_BASE/src
source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh
grid-proxy-init -valid 720:00
eval `scramv1 runtime -sh`
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh
cd -
