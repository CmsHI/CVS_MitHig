#!/bin/sh

tag=$1
run=$2
output_dest="/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin/reco/$tag/merged"

cmsRun merge_cfg.py 1> merge.out 2>merge.err

output=merged.root
if test -e $output; then
    edmEventSize $output
    check_file=$?
    if [ $check_file = "0" ]; then
	echo /opt/dcap/bin/dccp $output $output_dest/${tag}_$run.root
	/opt/dcap/bin/dccp $output $output_dest/${tag}_$run.root
	echo "File copied to Output Directory : $output" >&2
    else
	echo "Bad File : $output" >&2
	mv $output ${output}.bad
    fi
else
    echo "Failed to produce $output" >&2
    touch ${output}.none
fi

rm *.root



