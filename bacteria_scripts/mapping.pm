#!/usr/bin/perl -w

package mapping;

use strict;
use warnings;

# Smalt parameters
my $smalt_k = 13;
my $smalt_s = 6;

my $smalt_max_insert = 750;

my $snap_location = "/nfs/users/nfs_j/jl11/software/bin/snap";

# Random string of 8 alphanumeric characters
sub random_string()
{
   my $string = join'', map +(0..9,'a'..'z','A'..'Z')[rand(10+26*2)], 1..8;

   return $string;
}

# Sorts a given bam file to coordinate order
sub sort_bam($)
{
   my ($bam_file) = @_;

   my $bam_sort_prefix = "tmp" . random_string();
   my $tmp_bam = random_string() . "$bam_file";

   my $sort_command = "samtools sort -O bam -o $tmp_bam -T $bam_sort_prefix $bam_file";
   system($sort_command);

   rename $tmp_bam, $bam_file;
}

# Run bwa mem, producing sorted and indexed bam. Returns location of this bam
sub run_snap($$$$)
{
   my ($reference_name, $sample_name, $forward_reads, $reverse_reads) = @_;

   # Create index for reference if required
   if (!-d "snap_index")
   {
      snap_index($reference_name);
   }

   my $output_name = "$sample_name.mapping.bam";
   my $output_name_tmp = $output_name . random_string();
   my $log_file = "$sample_name.mapping.log";

   my $snap_command = "$snap_location paired snap_index $forward_reads $reverse_reads -o $output_name -so -sm 1 -R '" . join('\t', '@RG', "ID:$sample_name", "PL:ILLUMINA", "SM:$sample_name") . "'";
   system($snap_command);

   system("samtools fixmate -O bam $output_name $output_name_tmp");
   rename $output_name_tmp, $output_name;

   system("samtools index $output_name");

   return($output_name);
}

sub snap_index($)
{
   my ($reference_name) = @_;

   if (!-d "snap_index")
   {
      mkdir "snap_index";
   }

   system("$snap_location index $reference_name snap_index -large");
}

# Run bwa mem, producing sorted and indexed bam. Returns location of this bam
sub bwa_mem($$$$)
{
   my ($reference_name, $sample_name, $forward_reads, $reverse_reads) = @_;

   # Create index for reference if required
   if (!-e "$reference_name.bwt")
   {
      bwa_index($reference_name);
   }

   my $output_name = "$sample_name.mapping.bam";
   my $log_file = "$sample_name.mapping.log";

   my $bwa_command = "bwa mem -R '" . join('\t', '@RG', "ID:$sample_name", "PL:ILLUMINA", "SM:$sample_name") . "' $reference_name $forward_reads $reverse_reads 2>> $log_file | samtools fixmate -O bam - $output_name";
   system($bwa_command);

   sort_bam($output_name);
   system("samtools index $output_name");

   return($output_name);
}

sub bwa_index($)
{
   my ($reference_name) = @_;

   system("bwa index $reference_name");
}

# Uses smalt to map paired end reads to an indexed reference. Returns the
# location of the sorted, indexed bam file produced
sub run_smalt($$$$)
{
   my ($reference_name, $sample_name, $forward_reads, $reverse_reads) = @_;

   if (!-e "$reference_name.smi")
   {
      smalt_index($reference_name);
   }

   my $output_name = "$sample_name.mapping.bam";
   my $log_file = "$sample_name.mapping.log";

   my $smalt_command = "smalt map -f samsoft -p -i $smalt_max_insert $reference_name $forward_reads $reverse_reads 2> $log_file | samtools fixmate -O bam - $output_name";
   system($smalt_command);

   sort_bam($output_name);
   system("samtools index $output_name");

   return($output_name);
}

sub smalt_index($)
{
   my ($reference_name) = @_;

   system("smalt index -k $smalt_k -s $smalt_s $reference_name $reference_name");
}

1;
