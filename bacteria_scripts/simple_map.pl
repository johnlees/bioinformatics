#!/usr/bin/perl -w

#
# Uses htslib v1.2 to map haploid sample's fastq files to one of their
# assemblies
#

use strict;
use warnings;

use Getopt::Long;
use File::Spec;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 ) . "/../assembly_scripts";
use lib dirname( abs_path $0 );

use quake_wrapper;
use assembly_common;
use mapping;

#
# Help and usage
#
my $usage_message = <<USAGE;
Usage: ./simple_map.pl -a <reference.fasta> -r <reads_file> -o <output_prefix> <OPTIONS>

Maps input reads to reference. Marks duplicates

   Options

   Required
   -a, --assembly      Fasta file of reference to map to.
   -r, --reads         Tab delimited file of fastq or fastq.gz file locations.
                       One line per sample. First column sample name, second
                       column forward reads, third column reverse reads.
                       It is best to give absolute paths
   -o, --output        Prefix for output files

   Optional
   --dirty             Don't clean output
   -h, --help          Shows more detailed help.

Requires: samtools, bwa, picard
USAGE

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($reference_file, $read_file, $output_prefix, $no_picard, $dirty, $help);
GetOptions ("assembly|a=s"  => \$reference_file,
            "reads|r=s"  => \$read_file,
            "output|o=s" => \$output_prefix,
            "no_picard" => \$no_picard,
            "dirty" => \$dirty,
            "help|h"     => \$help
		   ) or die($usage_message);

# Parse input
if (defined($help))
{
   print $usage_message;
}
elsif (!defined($reference_file) || !defined($read_file) || !defined($output_prefix))
{
   print STDERR $usage_message;
}
else
{
   # Get array of samples and hash of read locations
   my ($samples, $reads) = quake_wrapper::parse_read_file($read_file, 1);

   # Map read pairs to reference with smalt (output sam files). Use threads
   # where possible
   print STDERR "Now mapping...\n\n";

   mapping::bwa_index($reference_file);

   assembly_common::add_tmp_file("$reference_file.amb");
   assembly_common::add_tmp_file("$reference_file.ann");
   assembly_common::add_tmp_file("$reference_file.bwt");
   assembly_common::add_tmp_file("$reference_file.pac");
   assembly_common::add_tmp_file("$reference_file.sa");


   my $sample = $$samples[0];
   my $forward_reads = $$reads{$sample}{"forward"};
   my $reverse_reads = $$reads{$sample}{"backward"};

   my $bam_file = mapping::run_bwa($reference_file, $sample, $forward_reads, $reverse_reads);

   assembly_common::add_tmp_file("$sample.mapping.log");

   # Use gatk to improve bams
   # This is done serially as processing time is short but mem use is large
   # GATK can use threads effectively itself
   unless (defined($no_picard))
   {
      print STDERR "bam improvement...\n\n";
      mapping::create_fa_dict($reference_file);

      # Use Picard to remove duplicates (cannot rename output)
      my $improved_bam = mapping::mark_dups($bam_file);

      assembly_common::add_tmp_file("$bam_file.picard.log");
      assembly_common::add_tmp_file("$bam_file.dups");

      assembly_common::add_tmp_file("$reference_file.dict");
      assembly_common::add_tmp_file("$reference_file.fai");
   }

   # Remove temporary files
   unless(defined($dirty))
   {
      print STDERR "removing temporary files...\n\n";
      assembly_common::clean_up();
   }
   print STDERR "Done.\n";
}

exit(0);

