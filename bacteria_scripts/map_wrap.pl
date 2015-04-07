#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use File::Spec;

my $usage = <<HELP;
./map_wrap.pl
Usage: bsub -J "map[1-20]" ./map_wrap.pl --reference ref.fa --lanes lanes.txt

Options:
   --reference     Fasta reference file to map to
   --lanes         A text file, each line a lane to map to ref

   --job           Run a specific job (line of above file)
   --cores         Number of cores to use (default 2)
HELP

sub get_fastq($$)
{
   my ($lane, $output_dir) = @_;

   my $fastq_dir = `pathfind -t lane -i $lane`;
   chomp($fastq_dir);

   symlink "$fastq_dir/$lane" . "_1.fastq.gz", "$output_dir/$lane" . "_1.fastq.gz";
   symlink "$fastq_dir/$lane" . "_2.fastq.gz", "$output_dir/$lane" . "_2.fastq.gz";
}

my ($reference_file, $lane_file, $job_in, $map_cores, $help);
GetOptions ("reference|a=s"  => \$reference_file,
            "lanes|l=s" => \$lane_file,
            "job=s" => \$job_in,
            "cores|c=s"  => \$map_cores,
            "help|h"     => \$help
		   ) or die($usage);

my $job = $ENV{'LSB_JOBINDEX'};
if (!defined($job))
{
   $job = $job_in;
}

if (!defined($map_cores))
{
   $map_cores = 2;
}

if (defined($help))
{
   print $usage;
}
else
{
   my $sed_command = "sed '$job" . "q;d' $lane_file";
   my $sample = `$sed_command`;
   chomp $sample;

   $sample =~ m/^(\d+)_(\d+)[#_](\d+)/;
   my $sample_name = "$1_$2_$3";
   my $sample_pathfind_name = "$1_$2#$3";

   if (!-d $sample_name)
   {
      mkdir $sample_name;
      get_fastq($sample_pathfind_name, $sample_name);

      chdir $sample_name;

      open(READS, ">reads.txt") || die("Could not write to $sample_name/reads.txt\n");
      print READS join("\t", $sample_name, File::Spec->rel2abs("$sample_pathfind_name" . "_1.fastq.gz"), File::Spec->rel2abs("$sample_pathfind_name" . "_2.fastq.gz")) . "\n";
      close READS;
   }
   else
   {
      chdir $sample_name;
   }

   my $map_command = "perl ~/bioinformatics/bacteria_scripts/map_snp_call.pl -m bwa -a $reference_file -r reads.txt -o $sample_name -p 1e-3 -t $map_cores";
   print STDERR $map_command . "\n";

   system($map_command);
}

exit(0);

