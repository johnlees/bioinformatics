#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $dsm_location = "/nfs/users/nfs_j/jl11/installations/metaminer_dev";

sub check_status($)
{
   my ($jobid) = @_;

   my $status;

   # Check the status of the jobid specified, returning status and exit code
   # comma separated
   my $bjobs = `bjobs -a -noheader -o "stat exit_code delimiter=','" $jobid`;
   chomp($bjobs);

   my ($bjobs_stat, $exit_code) = split(/,/, $bjobs);

   if ($bjobs_stat eq "RUN")
   {
      # Job still in queue, or still running
      $status = "RUN";
   }
   elsif ($bjobs_stat eq "PEND")
   {
      $status = "PEND";
   }
   elsif ($bjobs_stat eq "EXIT" && $exit_code eq "130")
   {
      # Job terminated by LSF - check it's due to memory limit being exceeded
      # using another bjobs command
      my $mem_over = `bjobs -l $jobid | grep -l "TERM_MEMLIMIT" | wc -l`;
      chomp($mem_over);

      if ($mem_over eq "1")
      {
         $status = "MEMLIMIT";
      }
      else
      {
         $status = "Status: $bjobs_stat. Exit code: $exit_code. NOT memlimit exceeded";
      }
   }
   elsif ($bjobs_stat eq "DONE")
   {
      # Job done - check the LSF output to check it was actually completed
      my $successfully_complete = `bjobs -l $jobid | grep -l "Done successfully" | wc -l`;
      chomp($successfully_complete);

      if ($successfully_complete eq "1")
      {
         $status = "DONE";
      }
      else
      {
         $status = "Status: $bjobs_stat. Exit code: $exit_code";
      }
   }
   else
   {
      # Other status
      $status = "Status: $bjobs_stat. Exit code: $exit_code";
   }

   return $status;
}

my $usage_message = <<USAGE;
./run_client.pl -p 52000 --processes 4 --mindepth 4 --maxdepth 100 --emax 19 --pmin 1

Options
   -p, --port     Port range to communicate across
   --process_depth    Number of clients to run
   --mindepth     Min depth
   --maxdepth     Max depth
   --emax         Maxmimum entropy to report
   --pmin         Minimum number of samples a kmer must appear in

   -h, --help     This message
USAGE

#* gets input parameters
my ($start_port, $processes, $min_depth, $max_depth, $max_entropy, $pmin, $help);
GetOptions ("port|p=s" => \$start_port,
            "process_depth=s" =>\$processes,
            "mindepth=s" => \$min_depth,
            "maxdepth=s" => \$max_depth,
            "emax=s" => \$max_entropy,
            "pmin=s" => \$pmin,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help))
{
   print $help;
}
else
{
   system("cat dsm-tmp/meta-server_config_*.txt > dsm-tmp/hostinfo.txt");

   my $manager_id = `bsub -o manager.%J.o -e manager.%J.e mpirun $dsm_location/meta-manager -v -H$processes -p$start_port dsm-tmp/jobs.txt`;
   # Job <1592673> is submitted to default queue <normal>
   $manager_id =~ m/^Job <(\d+)>/;
   my $check_id = $1;

   sleep 5;
   while (check_status($check_id) eq "PEND")
   {
      sleep 10;
   }

   my $manager_host = `bjobs -a -noheader -o "exec_host" $check_id`;
   chomp $manager_host;

   my $num_jobs = 4 ** $processes;
   for (my $i = 1; $i <= $num_jobs; $i++)
   {
      system("bsub -o client.%J.$i.o -e client.%J.$i.e 'cat dsm-tmp/hostinfo.txt | mpirun $dsm_location/meta-client --pmin $pmin --emax $max_entropy --mindepth $min_depth --maxDepth $max_depth --printNames -v --debug -M $manager_host:$start_port | gzip - > output.$i.txt.gz'");
   }
}

exit(0);

