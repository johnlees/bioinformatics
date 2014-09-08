#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use File::Spec;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );
use lib dirname( abs_path $0 ) . "/../assembly_scripts";

use quake_wrapper;
use gff_to_vcf;

use threads;

#
# Globals
#

# Smalt parameters
my $smalt_k = 13;
my $smalt_s = 6;

my $smalt_max_insert = 750;

# htslib parameters
my $max_depth = 1000;

#
# Help and usage
#
my $usage_message = <<USAGE;
Usage: ./map_snp_call.pl -r <reference.fasta> -a <reference_annotation.gff> -r <reads_file> -o <output_prefix>

Maps input reads to reference, and calls snps

   Options
   -a, --assembly      Fasta file of reference to map to.
   -g, --annotation    Annotation of the reference to annotate variants with.
   -r, --reads         Tab delimited file of fastq or fastq.gz file locations.
                       One line per sample. First column sample name, second
                       column forward reads, third column reverse reads.
                       It is best to give absolute paths
   -o, --output        Prefix for output files
   -t, --threads       Number of threads to use. Default 1
   -p, --prior         Prior for number of snps expected. Default 1e-3
   -h, --help          Shows more detailed help.

Requires: smalt, samtools, bcftools
USAGE

#****************************************************************************************#
#* Subs                                                                                 *#
#****************************************************************************************#

# Uses smalt to map paired end reads to an indexed reference. Returns the
# location of the sam file produced
sub run_smalt($$$$)
{
   my ($reference_name, $sample_name, $forward_reads, $reverse_reads) = @_;

   my $output_name = "$sample_name.mapping.sam";
   my $log_file = "$sample_name.mapping.log";

   my $smalt_command = "smalt map -f samsoft -i $smalt_max_insert -o $output_name $reference_name $forward_reads $reverse_reads &> $log_file";
   system($smalt_command);

   return($output_name);
}

# Converts sam to bam, sorts and indexes bam, deletes sam
# Returns name of bam
sub sort_sam($)
{
   my ($sam_file) = @_;

   # Set up file names
   my ($volume ,$directories, $file) = File::Spec->splitpath($sam_file);
   $file =~ m/^(.+)\.sam$/;
   my $file_prefix = $1;

   my $bam_file = $file_prefix . ".bam";

   # Do sorting and conversion (requires directory for temporary files
   my $sort_command = "samtools sort -O bam -o $bam_file -T /tmp/$file_prefix $sam_file";

   if (!-d "tmp")
   {
      mkdir "tmp" || die("Could not make directory tmp: $!\n");
   }

   print STDERR "$sort_command\n";
   system($sort_command);
   unlink $sam_file;

   # Index
   system("samtools index $bam_file");

   return($bam_file);
}

# Takes a list of bams and their sample ids and merges them into a single bam
sub merge_bams($$$)
{
   my ($samples, $bam_files, $output_bam) = @_;

   # Write a header of read groups and sample names to use in the merge
   my $rg_name = "rg.txt";
   open(RG, ">$rg_name") || die("Could not open $rg_name for writing: $!\n");

   my $i = 0;
   foreach my $sample (@$samples)
   {
      $$bam_files[$i] =~ m/^(.+)\.bam/;
      my $id = $1;

      print RG join("\t", '@RG', "ID:$id", "PL:ILLUMINA", "SM$sample") . "\n";
      $i++;
   }

   close RG;

   # Do merge
   my $merge_command = "samtools merge -rh $rg_name $output_bam " . join(" ", @$bam_files);
   system ($merge_command);

   system ("samtools index $output_bam");
}

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($reference_file, $annotation_file, $read_file, $output_prefix, $threads, $prior, $help);
GetOptions ("assembly|a=s"  => \$reference_file,
            "annotation|g=s" => \$annotation_file,
            "reads|r=s"  => \$read_file,
            "output|o=s" => \$output_prefix,
            "threads|t=i" => \$threads,
            "prior|p=s"  => \$prior,
            "help|h"     => \$help
		   ) or die($usage_message);

# Parse input
if (defined($help))
{
   print $usage_message;
}
elsif (!defined($reference_file) || !defined($annotation_file) || !defined($read_file) || !defined($output_prefix))
{
   print STDERR $usage_message;
}
else
{
   # Set default threads
   if (!defined($threads))
   {
      $threads = 1;
   }

   # Get array of samples and hash of read locations
   my ($samples, $reads) = quake_wrapper::parse_read_file($read_file, 1);

   # Index reference
   my ($volume ,$directories, $file) = File::Spec->splitpath($reference_file);
   $file =~ m/^(.+)\.(fasta|fa)$/;
   my $ref_name = $1;

   symlink $reference_file, "./$ref_name.fa";
   $reference_file = "./$ref_name.fa";

   system("smalt index -k $smalt_k -s $smalt_s $ref_name $reference_file");

   # Map read pairs to reference with smalt (output sam files). Use threads
   # where possible
   my @sam_files;

   for (my $i=0; $i<scalar(@$samples); $i+=$threads)
   {
      my @smalt_threads;
      for (my $thread = 1; $thread <= $threads; $thread++)
      {
         my $sample = $$samples[$i+$thread-1];
         my $forward_reads = $$reads{$sample}{"forward"};
         my $reverse_reads = $$reads{$sample}{"backward"};

         push(@smalt_threads, threads->create(\&run_smalt, $ref_name, $sample, $forward_reads, $reverse_reads));
      }

      # Wait for threads to rejoin before starting more
      foreach my $smalt_thread (@smalt_threads)
      {
         push (@sam_files, $smalt_thread->join());
      }

   }

   # Convert sam to bam, and sort. Again, threading where possible
   my @bam_files;
   for (my $i=0,; $i<scalar(@sam_files); $i+=$threads)
   {
      my @sort_threads;
      for (my $thread =1; $thread <= $threads; $threads++)
      {
         push(@sort_threads, threads->create(\&sort_sam, $sam_files[$i+$thread-1]));
      }

      foreach my $sort_thread (@sort_threads)
      {
         push(@bam_files, $sort_thread->join());
      }
   }

   # Merge bams
   # Sample array is in the same order as bam file name array
   my $merged_bam = "$output_prefix.merged.bam";
   merge_bams($merged_bam, $samples, \@bam_files);

   # Call variants, running mpileup and then calling through a pipe
   my $calling_command;
   my $output_vcf = "$output_prefix.vcf.gz";

   # Write a sample file, defining all as haploid
   my $sample_file = "samples.txt";
   open SAMPLES, ">$sample_file" || die("Could not write to $sample_file: $!\n");
   foreach my $sample (@$samples)
   {
      print SAMPLES "$sample 1\n";
   }
   close SAMPLES;

   # Use prior if given
   if (defined($prior))
   {
      $calling_command = "samtools mpileup -d $max_depth -t DP,SP -ug -f $reference_file $merged_bam | bcftools call -vm -P $prior -S $sample_file -O z -o $output_vcf";
   }
   else
   {
      $calling_command = "samtools mpileup -d $max_depth -t DP,SP -ug -f $reference_file $merged_bam | bcftools call -vm -S $sample_file -O z -o $output_vcf";
   }

   system($calling_command);

   # Annotate variants
   vcf_to_gff::transfer_annotation($annotation_file, $output_vcf);
   system("bcftools index $output_vcf");

   # Produced diff only vcf
   my $min_ac = 1;
   my $max_ac = scalar(@$samples) - 1;
   my $diff_vcf_name = "$output_prefix.diff.vcf.gz";
   system("bcftools view -c $min_ac -C $max_ac -o $diff_vcf_name -O z $output_vcf");
   system("bcftools index $diff_vcf_name");

   # TODO: bcftools stats, plot-vcfstats, bcftools filter
}

exit(0);
