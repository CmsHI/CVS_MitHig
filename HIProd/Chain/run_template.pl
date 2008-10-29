#!/usr/bin/perl

# $ARGC is the number of parameters
# @ARGV is the paramters (files needed to be run)

$outdir="/pnfs/cmsaf.mit.edu/hibat/cms/users/yetkin";
$tag="__TAG__";
@dir=("gen","sim","mix","digi","reco","time");

$maxtry=2;

chomp @ARGV;
$localdir=`pwd`;
chomp $localdir;
$done = 0;
$try = 0;
$input="";
$step=0;
$runnum="";

for($num = 0; $num < scalar @ARGV; $num++){
	
    $runnum = $ARGV[$num];

    for($step = 1; $step <= 5; $step++){ 
	$done = 0;
	$try = 0;
	while($done == 0 && $try < $maxtry){
	    $input=`ls $outdir/$dir[${step}]/${tag}/${tag}_r${runnum}.root`;
	    $nextstep = $step+1;
	    if($step == 2) $input=`ls $outdir/$dir[${nextstep}]/${tag}/${tag}_r${runnum}.root`;
	    chomp $input;
	    if($input eq "$outdir/$dir[${step}]/${tag}/${tag}_r${runnum}.root"){
		print "File for this run already exists. Retrieving from dcache... \n";
                print "opt/bin/dccp ${input} ${localdir}/${runnum}_${step}.root \n";
                `/opt/bin/dccp ${input} ${localdir}/${runnum}_${step}.root`;
	    }else{
		`echo "cmsRun -p ${runnum}_cfg$step.py 1> $runnum.$step.$try.out 2> $runnum.$step.$try.err"`;
		`cmsRun -p ${runnum}_cfg$step.py 1> $runnum.$step.$try.out 2> $runnum.$step.$try.err`;
	    }
	    &assert($runnum);
	}
    }
}

exit;

sub assert
{
   # check if the HID root file is valid and has 1 event

    $try++;
    $edm = $_[0] . "_$step.root"; 
    $inputstep=$step-1;
    $foundedm = 0;

    @edmsizes = `edmEventSize -n Events -v $edm`;
    if($? == 0){
        foreach $line2 (@edmsizes)
        {
            if($line2 =~ /edmHepMCProduct/)
            {
                $foundedm = 1;
                last;
            }
        }
    }

    if($step == 5){
	$foundedm = 1;
    }
    
    if($foundedm == 0)   # oops, file broken
    {
	`echo Bad Run : $_[0]`;
	`mv $edm $edm.bad`;

         # was this file retrieved form dcache? clean it
	if($input eq "$outdir/$dir[${step}]/${tag}/${tag}_r${runnum}.root"){
	    print "File for this run exists in dcache but corrupted. Deleting... \n";
	    print "rm $input \n";	
	    `rm $input`;
	}
	if($try == $maxtry){
	    $reco = $_[0] . "_4.root";
	    `echo ${runnum}.${step} >> BadRunList`;
            `cat merge_cfg.py | replace \"\\\"file:${reco}\\\",\"  \" \" | replace \", \\\"file:${reco}\\\"\"  \" \" > tm.py`;
	    `mv tm.py merge_cfg.py`;
	}
    }
    else   # put a file (tag) for resubmitting check
    {
	$done = 1;
	`md5sum $edm >> $edm.tag`;
	print $dir[$inputstep];
        if($dir[$inputstep] ne "reco"){
	    print "rm ${runnum}_${inputstep}.root";
	    `rm ${runnum}_${inputstep}.root`;
	}
	if($input eq "" && ($step == 1 || $step == 3 || $step == 4)){
	    `/opt/bin/dccp $edm $outdir/$dir[$step]/${tag}/${tag}_r$_[0].root`;
	}
    }
}


