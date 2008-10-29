#!/usr/bin/perl

# parameter: 1. event list directory.  The files has to named in numbers
#            2. background files
#
# for template configuration files, assume template_cfg.py

#chomp $ARGV[0];
chomp $ARGV[2];

# @filelist = `ls $ARGV[0]`;

$number_of_events = $ARGV[0];

$first_event = $ARGV[1];

@joblist = ();

$run = 13;
$event_per_run = 5;
$run_per_job = 10;
$do_edm = 1;

$tag="hydjet_x2_mb_oldPL_d20081021";

$last_event = $first_event + $number_of_events - 1;

@events = `seq $first_event $last_event`;

$work_dest = "DataFrom" . $first_event . "To" . $last_event;

`mkdir $work_dest`;

$dummy = 0;
$folder = 1;
$event = $first_event;

while($event <= $last_event){

#   $grandrunnum = `printf "%03d_%06d" $folder $event`;
    $grandrunnum = `printf "%06d" $event`;
    $random = int(rand(999999));
    
    `cat cfg1.py |sed "s/__MAXEVENTS__/$event_per_run/g" | sed "s/__SKIP__/$skip/g" | sed "s/__OUTPUT__/${grandrunnum}_1.root/g" | sed "s/__RANDOM__/$random/g" | sed "s/__MIX__/$backgroundlist[$background]/g" | sed "s/__INPUT__/$signalfile/g" | sed "s/__LIST__/$grandrunnum/g" | sed "s/__FIRSTEVENT__/$first_event/g" | sed "s/__RUN__/$run/g" >> ${grandrunnum}_cfg1.py`;
    `cat cfg2.py | sed "s/__OUTPUT__/${grandrunnum}_2.root/g" | sed "s/__RANDOM__/$random/g" | sed "s/__MIX__/$backgroundlist[$background]/g" | sed "s/__INPUT__/${grandrunnum}_1.root/g" >> ${grandrunnum}_cfg2.py`;
    `cat cfg3.py | sed "s/__OUTPUT__/${grandrunnum}_3.root/g" | sed "s/__RANDOM__/$random/g" | sed "s/__MIX__/$backgroundlist[$background]/g" | sed "s/__INPUT__/${grandrunnum}_2.root/g" >> ${grandrunnum}_cfg3.py`;
    `cat cfg4.py | sed "s/__OUTPUT__/${grandrunnum}_4.root/g" | sed "s/__RANDOM__/$random/g" | sed "s/__MIX__/$backgroundlist[$background]/g" | sed "s/__INPUT__/${grandrunnum}_3.root/g" >> ${grandrunnum}_cfg4.py`;
    `cat cfg5.py | sed "s/__OUTPUT__/${grandrunnum}_5.root/g" | sed "s/__RANDOM__/$random/g" | sed "s/__MIX__/$backgroundlist[$background]/g" | sed "s/__INPUT__/${grandrunnum}_4.root/g" >> ${grandrunnum}_cfg5.py`;

    $dummy++;
    push @joblist, $grandrunnum;
    
    if($dummy == $run_per_job)
    {
	$dummy = 0;
	$folder++;
    }

    $event = $event + $event_per_run;
}

`cp condor_backup condor`;
$size = scalar @joblist;

$index = 0;
for($index = 0; $index < scalar @joblist; $index = $index + $run_per_job)
{
    
    `mkdir $work_dest/$joblist[$index]`;
    `cat run_template.pl | sed "s/__TAG__/${tag}/g" > $work_dest/$joblist[$index]/run_$joblist[$index].pl`;
    `cp check.pl $work_dest/$joblist[$index]/`;
    `cp finalize.sh $work_dest/$joblist[$index]/`;

    @inputfilelist = ();
    
    $filelist = "";
    
    for($index2 = 0; $index2 < $run_per_job; $index2++)
    {
	if($index2 + $index < scalar @joblist)
	{
	    `mv $joblist[$index2+$index]_cfg1.py $work_dest/$joblist[$index]`;
	    `mv $joblist[$index2+$index]_cfg2.py $work_dest/$joblist[$index]`;
            `mv $joblist[$index2+$index]_cfg3.py $work_dest/$joblist[$index]`;
            `mv $joblist[$index2+$index]_cfg4.py $work_dest/$joblist[$index]`;
            `mv $joblist[$index2+$index]_cfg5.py $work_dest/$joblist[$index]`;
	    
	    push @inputfilelist, $joblist[$index+$index2];
	    
	    if($index2 > 0)
	    {
		$filelist = $filelist . ", ";
	    }
	    
	    $filelist = $filelist . "\\\"file:" . $joblist[$index+$index2] . "_4.root\\\"";
	}
    }
        
#    $merged_output = $tag . "_r" . $grandrunnum . ".root";
    $merged_output="merged.root";
    `cat merge_cfg.py | sed "s/__OUTPUT__/$merged_output/g" | sed "s/__INPUT__/$filelist/g" > $work_dest/$joblist[$index]/merge_cfg.py`;
    
    $filestring = "";
    $commafilestring = "";
    foreach $argument (@inputfilelist)
    {
	chomp $argument;
	if(length $commafilestring != 0)
	{
	    $filestring = $filestring . " ";
	    $commafilestring = $commafilestring . ", ";
	}

	$filestring = $filestring . $argument;
        $commafilestring = $commafilestring . $argument . "_cfg1.py, " . $argument . "_cfg2.py, " . $argument . "_cfg3.py," . $argument . "_cfg4.py," . $argument . "_cfg5.py,";
    }

    $commafilestring = $commafilestring . ", merge_cfg.py, run_$joblist[$index].pl, run_$joblist[$index].sh, finalize.sh";

    `cat run_template.sh | sed "s/__INTAG__/$signal/g" >> $work_dest/$joblist[$index]/run_$joblist[$index].sh`;
    `echo "perl check.pl $filestring 1> check.out 2> check.err &" >> $work_dest/$joblist[$index]/run_$joblist[$index].sh`;
    `echo "perl run_$joblist[$index].pl $filestring" >> $work_dest/$joblist[$index]/run_$joblist[$index].sh`;

    `echo "./finalize.sh $tag $joblist[$index]" >> $work_dest/$joblist[$index]/run_$joblist[$index].sh`;
#    `echo "cmsRun merge_cfg.py 1> merge.out 2>merge.err" >> $work_dest/$joblist[$index]/run_$joblist[$index].sh`;

    `echo "edm=$do_edm" >>$work_dest/$joblist[$index]/run_$joblist[$index].sh`;
    `echo "hidout=hid_$joblist[$index].root" >> $work_dest/$joblist[$index]/run_$joblist[$index].sh`;
    `echo "edmout=edm_$joblist[$index].root" >> $work_dest/$joblist[$index]/run_$joblist[$index].sh`;
    `echo "Initialdir           = $work_dest/$joblist[$index]" >> condor`;
    `echo "Executable           = $work_dest/$joblist[$index]/run_$joblist[$index].sh" >> condor`;
    `echo "Arguments            = " >> condor`;
    `echo "transfer_input_files = $commafilestring, check.pl" >> condor`;
    `echo "Queue" >> condor`;
    `echo >> condor`;

}

