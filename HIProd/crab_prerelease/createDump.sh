#!/bin/sh

v1=$CMSSW_BASE

cd $1/src
eval `scramv1 ru -sh`
cd -

edmConfigDump cfg.py > dumped.py

cd $v1/src
eval `scramv1 ru -sh`
cd -



