#!/usr/bin/perl -w

#
# reference_free_variant_caller.pl
# arguments:   assembly.fa
#              annotation.gff
#              forward_reads.fastq
#              reverse_reads.fastq
#
# bsub this script with:
# bsub -o variation.%J.o -e variation.%J.e -R "select[mem>1000] rusage[mem=1000]" -M1000 -n6 -R "span[hosts=1]" ./reference_free_variant_caller.pl <arguments>

#
# Takes an assembly and two sets of paired reads and gets variation between
# samples using denovo assembly only. Output is a vcf where variants are called
# differently between the pair of samples, placed on the assembly provided, and
# annotated if they fall within a gene.
# Expected coverage for true variants is also output
#
# Requires cortex_var (>1.0.5.21), bcftools (>v1.0), quake (>v0.3) and
# jellyfish (>v2.0)
# Cortex binaries compiled:   cortex_var_31_c1
#                             cortex_var_63_c2
#
# To create an assembly for use with the script, try running (with parameters
# and jobids changed where appropriate):
# bsub -o spades.%J.o -e spades.%J.e -n4 -R "span[hosts=1]" -R "select[mem>6000] rusage[mem=6000]" -M6000 /nfs/users/nfs_j/jl11/software/bin/spades.py -o ./ -1 ../11822_8_30_1.fastq.gz -2 ../11822_8_30_2.fastq.gz --careful -t 4 -m 6 -k 21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85
# bsub -o SPAdes/filter.log -q small -w "done(1080103)" ~/bioinformatics/bacteria_scripts/filter_contigs.pl -i SPAdes/scaffolds.fasta -o SPAdes/scaffolds.filtered.fasta
# bsub -o improved_assemblies/SPAdes/improve_assembly.%J.o -e improved_assemblies/SPAdes/improve_assembly.%J.e -w "done(1081014)" -R "select[mem>3000] rusage[mem=3000]" -M3000 improve_assembly -a SPAdes/scaffolds.filtered.fasta -f 11822_8_30_1.fastq -r 11822_8_30_2.fastq -o improved_assemblies/SPAdes/
# and to annotate this
# bsub -w "done(1083756)" -o annotate.pipeline.%J.o -e annotate.pipeline.%J.e -M3000 -R "select[mem>3000] rusage[mem=3000]" -n4 -R "span[hosts=1]" annotate_bacteria -a scaffolds.scaffolded.gapfilled.length_filtered.sorted.fa --sample_name 2070227 --genus Streptococcus --cpus 4
#
# Note this is only use to place the called variants, and is not used during
# assembly of the pair
#

use strict;
use warnings;

use threads;
use Getopt::Long;

use vcf_to_gff;

# Global locations of software needed
my $bcftools_location = $vcf_to_gff::bcftools;
my $quake_location = "/nfs/users/nfs_j/jl11/software/bin/quake/quake.py";
my $cortex_wrapper = "/nfs/users/nfs_j/jl11/installations/CORTEX_release_v1.0.5.21/scripts/calling/run_calls.pl";
my $cortex_binaries = "/nfs/users/nfs_j/jl11/installations/CORTEX_release_v1.0.5.21/bin";

my @required_binaries = ("cortex_var_31_c1", "cortex_var_63_c2");

# Other global parameters
my $quake_threads = 4;

#
# Help and usage
#
my $usage_message = <<USAGE;
Usage: ./reference_free_variant_caller.pl -a <assembly.fasta> -g <annotation.gff> -r <reads_file>

Uses cortex_var to call variants between two samples without providing a reference
sequence

   Options
   -a, --assembly		Contigs of a de-novo assembly of one of the samples in
   	            	multi-fasta format. See help for more info.
   -g, --annotation	Annotation of the -a assembly in gff v3.
   -r, --reads	      Tab delimited file of fastq or fastq.gz file locations.
   	               One line per sample. First column forward reads, second
                     column reverse reads.
   -h, --help	   	Shows more detailed help.

Requires: cortex_var, bcftools, quake and jellyfish. See help for more details
USAGE

my $help_message = <<HELP;

 reference_free_variant_caller.pl
 arguments:   assembly.fa
              annotation.gff
              forward_reads.fastq
              reverse_reads.fastq

 bsub this script with:
 bsub -o variation.%J.o -e variation.%J.e -R "select[mem>1000] rusage[mem=1000]" -M1000 -n6 -R "span[hosts=1]" ./reference_free_variant_caller.pl <arguments>

 Takes an assembly and two sets of paired reads and gets variation between
 samples using denovo assembly only. Output is a vcf where variants are called
 differently between the pair of samples, placed on the assembly provided, and
 annotated if they fall within a gene.
 Expected coverage for true variants is also output

 Requires cortex_var (>1.0.5.21), bcftools (>v1.0), quake (>v0.3) and
 jellyfish (>v2.0)
 Cortex binaries compiled:   cortex_var_31_c1
                             cortex_var_63_c2

 To create an assembly for use with the script, try running (with parameters
 and jobids changed where appropriate):
 bsub -o spades.%J.o -e spades.%J.e -n4 -R "span[hosts=1]" -R "select[mem>6000] rusage[mem=6000]" -M6000 /nfs/users/nfs_j/jl11/software/bin/spades.py -o ./ -1 ../11822_8_30_1.fastq.gz -2 ../11822_8_30_2.fastq.gz --careful -t 4 -m 6 -k 21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85
 bsub -o SPAdes/filter.log -q small -w "done(1080103)" ~/bioinformatics/bacteria_scripts/filter_contigs.pl -i SPAdes/scaffolds.fasta -o SPAdes/scaffolds.filtered.fasta
 bsub -o improved_assemblies/SPAdes/improve_assembly.%J.o -e improved_assemblies/SPAdes/improve_assembly.%J.e -w "done(1081014)" -R "select[mem>3000] rusage[mem=3000]" -M3000 improve_assembly -a SPAdes/scaffolds.filtered.fasta -f 11822_8_30_1.fastq -r 11822_8_30_2.fastq -o improved_assemblies/SPAdes/
 and to annotate this
 bsub -w "done(1083756)" -o annotate.pipeline.%J.o -e annotate.pipeline.%J.e -M3000 -R "select[mem>3000] rusage[mem=3000]" -n4 -R "span[hosts=1]" annotate_bacteria -a scaffolds.scaffolded.gapfilled.length_filtered.sorted.fa --sample_name 2070227 --genus Streptococcus --cpus 4

 Note this is only use to place the called variants, and is not used during
 assembly of the pair

HELP

#****************************************************************************************#
#* Functions                                                                            *#
#****************************************************************************************#
sub check_binaries()
{
   my $missing = 0;

   foreach my $binary (@required_binaries)
   {
      if (!-e "$cortex_binaries/$binary")
      {
         print STDERR "Cortex binary $binary does not exist in $cortex_binaries. It must be built before this pipeline can be run";
         $missing = 1;
      }
   }

   return($missing);
}

sub quake_error_correct($)
{

}

sub prepare_reference($)
{

}

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($assembly_file, $annotation_file, $read_file, $help);
GetOptions ("assembly|a=s"  => \$assembly_file,
            "annotation|g=s" => \$annotation_file,
            "reads|r=s"  => \$read_file,
            "help|h"     => \$help
		   ) or die($usage_message);

# Parse input
if (!defined($assembly_file) || !defined($annotation_file) || !defined($read_file))
{
   print STDERR $usage_message;
}
elsif (!-e $assembly_file || !-e $annotation_file || !-e $read_file)
{
   print STDERR "One or more specified input files do not exist\n";
   print STDERR $usage_message;
}
elsif (defined($help))
{
   print $help_message;
}
elsif (check_binaries())
{
   print STDERR "See compile instructions in cortex documentation";
}
else
{
   my $reads = parse_read_file($read_file);

   # Thread to error correct reads
   my $quake_thread = threads->create(\&quake_error_correct, $reads);

   # Thread to prepare reference with cortex and stampy
   my $reference_thread = threads->create(\&prepare_reference, $assembly_file);

   # Run cortex once threads have finished


   # Output expected coverage for each variant type

   # Fix vcf output

   # Annotate vcf, and extract passed variant sites only
}

exit(0);

