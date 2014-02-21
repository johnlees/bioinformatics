#!/usr/bin/perl -w

#****************************************************************************************#
#* jlsub.pl                                                                             *#
#* Submits a job array to process the same command on multiple files as defined by an   *#
#* input file                                                                           *#
#****************************************************************************************#

#* NOTES
#  Run from directory files are in
#  Have filenames as relative, use something like ls *.bam >filelist to create input file
#  Currently memory resource request doesn't work, so use on pcs4 only
#  Put command in quotes
#  Use FILE in command to be replaced by filename in each array job
#  e.g. ~/scripts/jlsub.pl -f test_file_list -c "~sh16/scripts/bam_to_coverage_plot.py -b FILE -o ../coverage_plots/FILE_coverage.plot"


use strict;

use Getopt::Long;
use Cwd;

my $script_location = "~jl11/scripts";

my $help_message = <<END;
Usage: jlsub.pl <options> -f <file_list> -c "command"
Submits a job array to process the same command on multiple files as defined by an
input file of relative paths (easily created by e.g. ls *.bam > file_list)

Should be run from the directory containing the input files, specify output
directory in the command
	
	-f, --files    A newline separated list of files to process the command on
	-c, --command  The command to run, which should be enclosed in quotes.
	               Where the keyword 'FILE' appears in the command, this will
	               be replaced by the actual filename for each job in the array
	-l, --log      Prefix for the log files of stdout and stderr [default=logs]   
	-m, --mem      Memory required in GB [default=0.1]
	-j, --job      The name of the job [default=jlsub]
	-q, --queue    The name of the bsub queue to submit to [default=normal]
	-l, --limit    The maximum number of array jobs to process simultaneously
	               [default=number of files in file_list]
	-h, --help     Displays this help message
END

#* Read command line options
my ($logprefix, $memrequired, $jobname, $filelist, $command, $queuename, $joblimit, $help);
GetOptions ("log|l=s"  => \$logprefix,
            "mem|m=f"  => \$memrequired,
            "job|j=s"  => \$jobname,
            "files|f=s" => \$filelist,
            "command|c=s" => \$command,
            "queue|q=s" => \$queuename,
            "limit=o" => \$joblimit,
            "help|h" => \$help
		   ) or die($help_message);
		   

if (!defined($filelist) || !defined($command))
{
   print $help_message;
}
else
{		   
   my $dir = cwd;

   # Give options defaults if not set
   if (!defined($logprefix))
   {
      $logprefix = "logs";
   }

   # 100MB job default
   if (!defined($memrequired))
   {
      $memrequired = "100";
   }
   else
   {
      # Input supplied in GB
      $memrequired *= 1000;
   }

   if (!defined($queuename))
   {
      $queuename = "normal";
   }

   if (!defined($jobname))
   {
      $jobname = "jlsub";
   }

   #file length equals array length
   my $number_jobs = `wc -l $filelist`;
   if ($number_jobs =~ /^(\d+) .+/)
   {
      $number_jobs = $1;
   }

   #allow simultaneous execution of all jobs, unless otherwise specified
   if (!defined($joblimit))
   {
      $joblimit = $number_jobs;
   }	

   # Check which server being run on, for pcs4 don't bother with memory allocation
   my $bsub_command;
   my $lsid_return = `lsid`;
   if ($lsid_return =~ /My cluster name is (.+)\n/ && $1 eq "farm3")
   {
      #for farm3
      $bsub_command = "bsub -o $logprefix.%J.%I.out -e $logprefix.%J.%I.err -J\"$jobname" . "[1-$number_jobs]%" . "$joblimit\" -R\”select[mem>$memrequired] rusage[mem=$memrequired]\” -M$memrequired -q $queuename";
   }
   else
   {
      #for pcs4/5
      $bsub_command = "bsub -o $logprefix.%J.%I.out -e $logprefix.%J.%I.err -J\"$jobname" . "[1-$number_jobs]%" . "$joblimit\" -q $queuename";
   }
   
   # pass to wrapper script
   my $perl_command = "$script_location/jlsub_indexjob.pl $dir $filelist '$command'";

   # Run command, and show the LSF output
   print "Running command '$bsub_command $perl_command'\n";
   my $array_job_num = `$bsub_command $perl_command`;
   print "$array_job_num";
   
   # Send an email when we're done
   if ($array_job_num =~ /<(\d+)>/)
   {
      system("bsub -o /dev/null -q small -w 'done($1)' $script_location/jlmail.pl -i '$command' -j $1");
   }
}
    
exit(0);