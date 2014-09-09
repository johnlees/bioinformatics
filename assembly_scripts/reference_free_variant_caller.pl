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
# Requires cortex_var (>1.0.5.21), stampy, bcftools (>v1.0), quake (>v0.3) and
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
use Cwd;
use File::Spec;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

# Perl modules - assembly
use gff_to_vcf;
use quake_wrapper;
use assembly_common;

# Global locations of software needed
my $old_vcftools_location =  "/nfs/users/nfs_j/jl11/installations/vcftools_0.1.12a";
my $cortex_wrapper = "/nfs/users/nfs_j/jl11/installations/CORTEX_release_v1.0.5.21/scripts/calling/run_calls.pl";
my $cortex_binaries = "/nfs/users/nfs_j/jl11/installations/CORTEX_release_v1.0.5.21/bin";
my $stampy_location = "/nfs/users/nfs_j/jl11/installations/stampy-1.0.23/stampy.py";
my $bcftools_location = "/nfs/users/nfs_j/jl11/software/bin/bcftools";

my @required_binaries = ("cortex_var_31_c1", "cortex_var_63_c2");

# Other global parameters

# Quake
my $quake_threads = 4;
my $quake_kmer_size = 14;

# Run calls.pl
my $first_kmer = 31;
my $last_kmer = 61;
my $kmer_step = $last_kmer - $first_kmer;
my $auto_cleaning = "yes";
my $bc = "yes";
my $pd = "no";
my $out_dir = "ctx_out";
my $ploidy = "1";
my $ref_bin_dir = "ref_binaries";
my $mem_height = 18;
my $mem_width = 100;
my $do_union = "yes";
my $ref_usage = "CoordinatesOnly";
my $workflow = "joint";
my $ctx_logfile = "cortex.bc.log";
my $ref_se = "REF_se";

# bcftools
my $bcftools_stderr_file = "bcftools.err";

my $filter_lines = <<'FILTERS';
##FILTER=<ID=PF_FAIL_ERROR,Description="Failed population filter">
##FILTER=<ID=PF_FAIL_REPEAT,Description="Identified as a repeat by population filter">
FILTERS

#
# Help and usage
#
my $usage_message = <<USAGE;
Usage: ./reference_free_variant_caller.pl -a <assembly.fasta> -g <annotation.gff> -r <reads_file>

Uses cortex_var to call variants between two samples without providing a reference
sequence

   Options
   -a, --assembly      Contigs of a de-novo assembly of one of the samples in
                       multi-fasta format. See help for more info.
   -g, --annotation    Annotation of the -a assembly in gff v3.
   -r, --reads         Tab delimited file of fastq or fastq.gz file locations.
                       One line per sample. First column sample name, second
                       column forward reads, third column reverse reads.
                       It is best to give absolute paths
   --no-error-correct  Don't run quake on reads. Input reads are fed directly
                       into cortex
   -o, --output        Prefix for output vcfs
   -h, --help          Shows more detailed help.

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

 Requires cortex_var (>1.0.5.21), stampy, bcftools (>v1.0), quake (>v0.3) and
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

# Creates binaries and a stampy hash of a reference sequence
sub prepare_reference($)
{
   my ($reference_file) = @_;

   my @cortex_threads;
   my $ref_new = "reference.renamed.fa";

   # Standardise contig names
   assembly_common::standardise_contig_names($reference_file, $ref_new);

   # Put assembly.fa location into file for cortex input
   system("echo $ref_new > $ref_se");
   mkdir "$ref_bin_dir" || die("Could not make directory $ref_bin_dir: $!\n");

   # Dump a binary for each kmer being used
   my $i = 0;

   for (my $kmer = $first_kmer; $kmer<=$last_kmer; $kmer+= $kmer_step)
   {
      $cortex_threads[$i] = threads->create(\&build_binary, $kmer, 1, "$ref_se");
      $i++;
   }

   # Make a stampy hash of the reference
   my $stampy_command = "$stampy_location -G REF $ref_new";
   system("$stampy_command &> stampy.genome.log");

   $stampy_command = "$stampy_location -g REF -H REF";
   system("$stampy_command &> stampy.hash.log");

   # Wait for cortex jobs to finish
   foreach my $thread (@cortex_threads)
   {
      $thread->join();
   }

   # Give new reference name
   return($ref_new);

}

# Builds binaries of a reference using cortex
sub build_binary($$$)
{
   my ($kmer, $colours, $ref_file) = @_;

   # Find required binary
   my $correct_binary;
   foreach my $binary (@required_binaries)
   {
      if ($binary =~ /^cortex_var_(\d+)_c(\d+)$/)
      {
         if ($1 >= $kmer && (($1-30) <= $kmer) && $2 >= $colours)
         {
            $correct_binary = $binary;
         }
      }
   }

   # At the moment, should only encounter this due to an error in the script so
   # don't worry too much about giving a helpful error message
   if (!defined($correct_binary))
   {
      die("No suitable binary for kmer $kmer and $colours colour(s) exists.\n");
   }

   my $cortex_command = "$cortex_binaries/$correct_binary --kmer_size $kmer --mem_height $mem_height --mem_width $mem_width --se_list $ref_file --dump_binary $ref_bin_dir/ref.k$kmer.ctx --sample_id REF";
   system("$cortex_command &> $ref_bin_dir/ctx.ref.k$kmer.log");
}

# Gets the total sequence length of a multi-fasta file
sub reference_length($)
{
   my ($reference_file) = @_;

   my $length = `grep -v ">" $reference_file | wc -m`;
   chomp($length);

   return($length);
}

# Creates an index file for use by run_calls.pl
sub create_cortex_index($)
{
   # A reference to a hash of read locations
   # %$reads{sample}{direction}
   my ($reads) = @_;

   my $index_name = "INDEX";

   # INDEX file is one row per sample, columns: sample name, se reads, forward
   # pe reads, reverse pe reads. Empty fields marked by periods
   open(INDEX, ">$index_name") || die("Could not write to $index_name: $!\n");

   foreach my $sample (keys %$reads)
   {
      my $pe_1 = $sample . "_pe_1";
      my $pe_2 = $sample . "_pe_2";

      print INDEX join("\t", $sample, ".", $pe_1, $pe_2) . "\n";
      open(PE_F, ">$pe_1") || die("Could not open $pe_1 for writing: $!\n");
      open(PE_R, ">$pe_2") || die("Could not open $pe_2 for writing: $!\n");

      print PE_F $$reads{$sample}{"forward"} . "\n";
      print PE_R $$reads{$sample}{"backward"}. "\n";

      close PE_F;
      close PE_R;
   }

   close INDEX;
   return($index_name);
}

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($assembly_file, $annotation_file, $read_file, $no_quake, $output_prefix, $help);
GetOptions ("assembly|a=s"  => \$assembly_file,
            "annotation|g=s" => \$annotation_file,
            "reads|r=s"  => \$read_file,
            "no-error-correct" => \$no_quake,
            "output|o=s"  => \$output_prefix,
            "help|h"     => \$help
		   ) or die($usage_message);

# Parse input
if (defined($help))
{
   print $help_message;
}
elsif (!defined($assembly_file) || !defined($annotation_file) || !defined($read_file) || !defined($output_prefix))
{
   print STDERR $usage_message;
}
elsif (!-e $assembly_file || !-e $annotation_file || !-e $read_file)
{
   print STDERR "One or more specified input files do not exist\n";
   print STDERR $usage_message;
}
elsif (check_binaries())
{
   print STDERR "See compile instructions in cortex documentation";
}
else
{
   my $cwd = getcwd();
   my ($samples, $reads) = quake_wrapper::parse_read_file($read_file, 0);

   # Thread to error correct reads
   # Note this returns location of corrected reads
   my $quake_thread;
   unless ($no_quake)
   {
      print STDERR "Error correcting reads and preparing assembly\n";
      $quake_thread = threads->create(\&quake_wrapper::quake_error_correct, $reads, $quake_kmer_size, $quake_threads);
   }
   else
   {
      print STDERR "Preparing assembly\n";
   }

   # Thread to prepare reference with cortex and stampy
   my $reference_thread = threads->create(\&prepare_reference, $assembly_file);

   # Run cortex once threads have finished
   $assembly_file = $reference_thread->join();

   my $index_name;
   if ($no_quake)
   {
      $index_name = create_cortex_index($reads);
   }
   else
   {
      $index_name = create_cortex_index($quake_thread->join());
   }

   my $approx_length = reference_length($assembly_file);

   my $cortex_command = "perl $cortex_wrapper --first_kmer $first_kmer --last_kmer $last_kmer --kmer_step $kmer_step --fastaq_index $index_name --auto_cleaning $auto_cleaning --bc $bc --pd $pd --outdir $out_dir --outvcf $output_prefix --ploidy $ploidy --refbindir $ref_bin_dir --list_ref_fasta $ref_se --stampy_hash REF --stampy_bin $stampy_location --genome_size $approx_length --mem_height $mem_height --mem_width $mem_width --squeeze_mem --vcftools_dir $old_vcftools_location --do_union $do_union --ref $ref_usage --workflow $workflow --logfile $ctx_logfile --apply_pop_classifier";

   print STDERR "Running cortex\n";
   print STDERR "$cortex_command\n\n";
   system("$cortex_command &> cortex.err");

   # Output expected coverage for each variant type, but use correct
   # coverage and read length
   #
   # TODO: This assumes no further cleaning by cortex. It might be best either
   # just to take cortex's calculated value, or input the quality threhold to
   # quake rather than cortex
   my $coverage_fastq = $$reads{$$samples[0]}{"forward"};
   print STDERR "Using reads in $coverage_fastq for expected coverage calculations\n";

   # Expected coverage is 2 times this as we've only counted one end of each
   # read
   my $read_length = assembly_common::read_length($coverage_fastq);
   my $expected_coverage = 2*assembly_common::expected_coverage($coverage_fastq, $approx_length);

   print "Effective coverage at k=$first_kmer: " . assembly_common::effective_coverage($first_kmer, $read_length, $expected_coverage) . "\n";
   print "Effective coverage at k=$last_kmer: " . assembly_common::effective_coverage($last_kmer, $read_length, $expected_coverage) . "\n";

   # Fix vcf output
   # TODO: This isn't ideal as some of the name is hard coded when it can be inferred
   # from the run_calls.pl parameters
   my $output_vcf = "$cwd/$out_dir/vcfs/$output_prefix" . "_wk_flow_J_RefCO_FINALcombined_BC_calls_at_all_k.decomp.vcf";
   print STDERR "Cortex output vcf is: $output_vcf\n\n";
   print STDERR "Fixing and annotating vcf\n";

   # Reheader vcf with bcftools as pop filter FILTER fields not included (bug in cortex)
   my $fixed_vcf = "$cwd/$output_prefix.decomposed_calls.vcf";
   my $filtered_vcf = "$cwd/$output_prefix.filtered_calls.vcf.gz";

   system($bcftools_location . " view -h $output_vcf > vcf_header.tmp 2> $bcftools_stderr_file");
   system("head -n -1 vcf_header.tmp > vcf_header_reduced.tmp");
   system("tail -1 vcf_header.tmp > vcf_column_headings.tmp");

   open(FILTERS, ">filters.tmp") || die("Could not write to filters.tmp: $!\n");
   print FILTERS $filter_lines;
   close FILTERS;

   system("cat vcf_header_reduced.tmp filters.tmp vcf_column_headings.tmp > new_header.tmp");
   system($bcftools_location . " reheader -h new_header.tmp $output_vcf > $fixed_vcf 2>> $bcftools_stderr_file");
   unlink "vcf_header.tmp", "vcf_header_reduced.tmp", "vcf_column_headings.tmp", "filters.tmp", "new_header.tmp";

   # Fix error in filter fields introduced by population filter fields, then
   # bgzip and index
   system("sed -i -e 's/,PF/;PF/g' $fixed_vcf");
   system("bgzip $fixed_vcf");
   system($bcftools_location . " index $fixed_vcf.gz 2>> $bcftools_stderr_file");

   # Annotate vcf, and extract passed variant sites only
   vcf_to_gff::transfer_annotation($annotation_file, "$fixed_vcf.gz");

   system($bcftools_location . " view -C 2 -c 2 -f PASS $fixed_vcf.gz -o $filtered_vcf -O z 2>> $bcftools_stderr_file");
   system($bcftools_location . " index $filtered_vcf 2>> $bcftools_stderr_file");

   print "Final output:\n$filtered_vcf\n";
}

exit(0);

