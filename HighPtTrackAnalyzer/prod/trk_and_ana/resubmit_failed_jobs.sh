#!/bin/sh

if [ $# -ne 3 ]
then
	echo "This script takes three arguments (e.g. run_mixb0 trk_and_ana.py 50).  Exiting..."
	exit
fi

script=$1
cfg=$2

maxjobs=$3
job=0
resubmitted=0

case $script in
run_mixb0)
dir="hydjet_x2_b0_oldPL_d20081106";;
run_mixb6)
dir="hydjet_x2_b6_oldPL_d20081106";;
run_mixb12)
dir="hydjet_x2_b12_oldPL_d20081106";;
*)
echo "Don't know what to do with this script.  Exiting..."
exit;;
esac

while [ $job -lt $maxjobs ]
do
	run=`expr 1 + $job / 100`
	event=`expr  $job % 100 + 1`
	mixrun=`expr  $job % 100 + 1`
	
	runstr=`printf "%06d" $run`
	mixrunstr=`printf "%06d" $mixrun`
			
	if [ ! -f ~/dcache/trk_and_ana/${dir}/${dir}_r${mixrunstr}_embedded_with_pythia_dijet_pt100to9999_d20081021_r${runstr}_e${event}.root ]
	then
		echo "run ($run), event ($event), mixrun ($mixrun) needs to be re-submitted"
		resubmitted=`expr $resubmitted + 1`
		./condor.sh $script $cfg $run $event $mixrun
		sleep 2s
	fi
	
	job=`expr $job + 1`
done

echo "Resubmitted $resubmitted failed jobs out of $maxjobs."

