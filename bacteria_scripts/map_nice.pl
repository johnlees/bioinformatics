#!/usr/bin/perl -w

#
# Uses bwa to map, does initial checks on resulting bam
#

use strict;
use warnings;

use Getopt::Long;
use File::Spec;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use mapping;

#
# Help and usage
#
my $usage_message = <<USAGE;
Usage: ./map_nice.pl -r <reference.fasta> -i <read_prefix> -o <output_prefix>

Maps input reads to reference, marks duplicates and sorts

   Options

   Required
   -r, --reference     Fasta file of reference to map to.
   -i, --reads         Prefix of reads to map (ending _1.fastq.gz)
   -s, --sample        Sample name
   -o, --output        Output prefix

   Optional
   -t, --threads       Number of threads to use (default 1)

   -h, --help          Shows more detailed help.

USAGE

#* gets input parameters
my ($reference_file, $read_prefix, $output_prefix, $sample_name, $threads, $help);
GetOptions ("reference|r=s"  => \$reference_file,
            "reads|i=s" => \$read_prefix,
            "output|o=s" => \$output_prefix,
            "sample|s=s" => \$sample_name,
            "threads|t=i" => \$threads,
            "help|h"     => \$help
		   ) or die($usage_message);

# Parse input
if (defined($help))
{
   print $usage_message;
}
elsif (!defined($reference_file) || !defined($read_prefix) || !defined($output_prefix))
{
   print STDERR $usage_message;
}
else
{
   if (!defined($threads))
   {
      $threads = 1;
   }

   my $read_fw = "$read_prefix\_1.fastq.gz";
   my $read_rev = "$read_prefix\_2.fastq.gz";

   if (!-e $read_fw || !-e $read_rev)
   {
      die("Could not find forward ($read_fw) and/or reverse reads ($read_rev)\n")
   }

   my $mapped_bam = mapping::run_bwa($reference_file, $sample_name, $read_fw, $read_rev, $threads);
   mapping::sort_bam($mapped_bam, $threads);
   my $final_bam = mapping::mark_dups($mapped_bam);

   unlink($mapped_bam, glob("*.bai"));

   rename $final_bam, "$output_prefix.bam";
   system("samtools index $output_prefix.bam");
}

exit(0);

