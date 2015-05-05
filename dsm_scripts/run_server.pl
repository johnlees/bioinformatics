#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $dsm_location = "/nfs/users/nfs_j/jl11/installations/metaminer_dev";

my $usage_message = <<USAGE;
./run_server.pl -f assembly_fasta_files -p 52000 -d 100 --fmin 1

Options
   -M, --memory   Memory to use, in Mb

   -f, --files    List of fasta files to operate on
   -p, --port     Port range to communicate across
   -d, --depth    Max depth
   --fmin         Minimum number of samples a kmer must appear in

   -h, --help     This message
USAGE

#* gets input parameters
my ($memory, $file_list, $start_port, $max_depth, $fmin, $help);
GetOptions ("memory|M=s" => \$memory,
            "files|f=s"  => \$file_list,
            "port|p=s" => \$start_port,
            "depth|d=s" => \$max_depth,
            "fmin=s" => \$fmin,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help) || !-e $file_list)
{
   print $usage_message;
}
else
{
   if (-e "dsm-tmp")
   {
      die("tmp dir dsm-tmp already exists\n");
   }
   else
   {
      mkdir "dsm-tmp";
   }

   my $num_jobs = `wc -l $file_list | cut -d " " -f 1`;
   chomp $num_jobs;

   if (!defined($memory))
   {
      $memory = 100;
   }

   system("bsub -J \"server[1-$num_jobs]\" -o server.%J.%I.o -e server.%J.%I.e -R \"select[mem>$memory] rusage[mem=$memory]\" -M$memory mpirun $dsm_location/wrapper-LSF/server_wrapper.pl $fmin $max_depth $start_port $file_list");
}

exit(0);

