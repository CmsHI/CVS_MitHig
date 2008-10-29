#!/usr/bin/perl

$user = yetkin;
$max_exe_minutes = 1800;
$max_file_up_sec = 18000;

use Time::Local;

@list = @ARGV;
foreach $item (@list)
{
   chomp $item;
}

$no_count = 0;

while(1 == 1)
{
   `sleep 5`;

   @found = ();
   @runnum = ();

   foreach $item (@list)
   {
      @answer = `ps aux | grep $user | grep "cmsRun $item.cfg"`;
      foreach $answer_item (@answer)
      {
         chomp $answer_item;
         if(!($answer_item =~ /ps aux/) && !($answer_item =~ /grep/) && !($answer_item =~ /sh -c/))
         {
            # print "answer = \"$answer_item\"\n";
            push @found, $answer_item;
            push @runnum, $item;
         }
      }
   }

   if(scalar @found == 0)
   {
      $time = `date`; chomp $time;
      print "$time: No file are found running.\n";
      $no_count = $no_count + 1;
   }
   else   # found a file running, check for stuck
   {
      $no_count = 0;

      $item = $found[0];   # there should be only one process

      $stuck = 0;

      # first cut: total execution time
      $time = `echo $item | awk '{print \$10}'`;
      chomp $time;
      $time =~ s/:/ /g;
      $minute = `echo $time | awk '{print \$1}'`;
      if($minute >= $max_exe_minutes)
      {
         $stuck = 1;
         $time = `date`; chomp $time;
         print "$time: total execution time too long for run number $runnum[0]\n";
      }

      # second cut: last update time of the output file
      $now = time();
      $filename = $runnum[0] . ".out";
      $update = (stat("$filename"))[9];
      if($now - $update > $max_file_up_sec)
      {
         $stuck = 1;
         $time = `date`; chomp $time;
         print "$time: last update of output file too long for run number $runnum[0]\n";
      }

      # kill the process if it is stuck
      if($stuck != 0)
      {
         $pid = `echo $item | awk '{print \$2}'`;
         chomp $pid;
         $time = `date`; chomp $time;
         print "$time: killing process with pid = $pid, run number = $runnum[0]\n";
         `kill -kill $pid`;
      }
      else
      {
         $time = `date`; chomp $time;
         print "$time: process with run number $runnum[0] is still running\n";
      }
   }

   if($no_count == 5)
   {
      $time = `date`; chomp $time;
      print "$time: Not running for 5 consecutive checks.  Terminating....\n";
      last;
   }
}

$time = `date`; chomp $time;
print "$time: leaving the script\n";
exit;

