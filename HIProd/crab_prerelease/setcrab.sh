#!/bin/sh

source /osg/grid/setup.sh
grid-proxy-init

source /osg/app/crab/crab.sh
source /osg/app/glite/etc/profile.d/grid_env.sh

cd $1/src
. ~/setosg
cd -

#export EDG_WL_LOCATION=/osg/app/glite
