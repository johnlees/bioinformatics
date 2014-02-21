#!/usr/bin/perl -w

#****************************************************************************************#
#* resubmit.pl                                                                          *#
#  Resubmits LSF jobs that have failed due to running out of memory                     *#
#****************************************************************************************#

use strict;

use Getopt::Long;

my $script_location = "~jl11/scripts";

my $help_message = <<END;
Usage: resubmit.pl --job <jobid>
Resubmits a job with a higher memory request if it fails due to running out 
of memory

Run with bsub -o /dev/null -w 'ended(jobid)' -q small resubmit.pl --job <jobid> 
for an job which is already pending or running. Run outside of bsub for jobs 
which have already failed
	
    -j, --job           The jobid to resubmit
    -i, --increment     The increase in memory (in GiB) for each resubmission
                        [default=2]
    -b, --options       A string of options to pass to bsub (except memory
                        requirement)
    -r, --max_retries   Maximum number of retries before giving up (will email
                        upon giving up) [default = 5]
    -h, --help          Displays this help message
END

#* Read command line options
my ($jobid, $increment, $bsub_options, $max_retries, $current_retries, $help);
GetOptions ("job|j=s"  => \$jobid,
            "increment|i=o"  => \$increment,
            "options|b=s" => \$bsub_options,
            "max_retries|r=o"  => \$max_retries,
            "retries|c=o" => \$current_retries,
            "help|h" => \$help
		   ) or die($help_message);
		   

if (!defined($jobid))
{
   print $help_message;
}
else
{		   
   # Give options defaults if not set
   if (!defined($increment))
   {
      $increment = 2;
   }
   if (!defined($max_retries))
   {
      $max_retries = 5;
   }
   if (!defined($current_retries))
   {
      $current_retries = 0;
   }

   # Find info about finished job
   my $bjobs_info_command = 'bjobs -o "command stat max_mem memlimit sub_cwd delimiter=' . "','" . '"' . " $jobid";
   my @bjobs_line = split(/\n/, `$bjobs_info_command`);
   my ($command, $stat, $max_mem, $memlimit, $sub_cwd) = split(',', $bjobs_line[1]);
   
   # Parse memory fields
   if ($max_mem =~ /^(\d+) (.+)bytes$/)
   {
     if ($2 eq "G")
     {
     	$max_mem = $1*1000;
     }
     else
     {
     	$max_mem = $1;
     }
   }
   if ($memlimit =~ /^(\d+) (.+)$/)
   {
   	 if ($2 eq "G")
     {
     	$memlimit = $1*1000;
     }
     else
     {
     	$memlimit = $1;
     }
   }
   
   if ($stat eq "EXIT")
   {
     if ($current_retries > $max_retries || $max_mem < $memlimit)
     {
     	#Send email saying failed
     	# Find user to send message to
        my $user_name = `whoami`;
        chomp($user_name);

        open(TMP, ">tmpemail.txt") || die("Couldn't open tmpemail.txt");
        print TMP "Job $jobid failed after $max_retries retries using mem limit $memlimit.\n Command used was $command";
        close TMP;

        # send the message
        system("mail -s 'Job $jobid failed after $max_retries" .'@sanger.ac.uk < tmpemail.txt');
        system("rm tmpemail.txt");
     }
     else
     {
     	$current_retries++;
     	
     	my $new_mem = $max_mem + $increment*1000;
     	my $mem_resource = '-R "select[mem>' . "$new_mem] rusage[mem=$new_mem]" .'" ' . "-M$new_mem";
     	
     	my $new_submission = "bsub $bsub_options $mem_resource $command";
     	print "Running command $new_submission from $sub_cwd\n";
     	chdir($sub_cwd);
     	my $array_job_num = `$new_submission`;
     	
     	# Set this script to run again once the new job finishes
     	if ($array_job_num =~ /<(\d+)>/)
        {
          system("bsub -o /dev/null -q small -w 'ended($1)' $script_location/resubmit.pl -j $1 -i $increment -b $bsub_options -r $max_retries -c $current_retries");
        }
     }
   }
   
}
    
exit(0);