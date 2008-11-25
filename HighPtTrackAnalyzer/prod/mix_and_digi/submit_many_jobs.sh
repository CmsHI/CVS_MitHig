#!/bin/sh

if [ $# -ne 3 ]
then
	echo "This script takes three arguments (e.g. run_mixb0 mix_and_digi.py 50).  Exiting..."
	exit
fi

script=$1
cfg=$2

maxjobs=$3
job=0

while [ $job -lt $maxjobs ]
do
	run=`expr 1 + $job / 100`
	event=`expr  $job % 100 + 1`
	mixrun=`expr  $job % 100 + 1`
	./condor.sh $script $cfg $run $event $mixrun
	echo "Submitted to condor using script ($script), config ($cfg), run ($run), event ($event), and mixrun ($mixrun)";
	job=`expr $job + 1`
	sleep 2s
done



