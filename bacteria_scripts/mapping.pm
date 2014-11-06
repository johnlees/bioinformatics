#!/usr/bin/perl -w

package mapping;

use strict;
use warnings;

#
# Globals
#

# Software locations
my $java_location = "/software/bin/java";

# Set java proxy for sanger cluster
my $java_flags = "-DproxySet=true -Dhttp.proxyHost=wwwcache.sanger.ac.uk -Dhttp.proxyPort=3128";
my $picard_location = "/nfs/users/nfs_j/jl11/software/bin/picard.jar";
my $gatk_location = "/nfs/users/nfs_j/jl11/software/bin/GenomeAnalysisTK.jar";
my $gatk_key_location = "/nfs/users/nfs_j/jl11/installations/jl11_sanger.ac.uk.key";
my $gatk_key = "";

# Smalt parameters
my $smalt_k = 13;
my $smalt_s = 6;

my $smalt_max_insert = 750;

my $snap_location = "/nfs/users/nfs_j/jl11/software/bin/snap";

# Set java flags on module load. Use GATK key to disable phone home
INIT
{
   my $hostname = `hostname`;
   chomp($hostname);

   if ($hostname =~ /^farm3/ || $hostname =~ /pcs5/)
   {
      $java_location .= " $java_flags";
      $gatk_location .= " -et NO_ET -K $gatk_key_location";
   }
}

# Random string of 8 alphanumeric characters
sub random_string()
{
   my $string = join'', map +(0..9,'a'..'z','A'..'Z')[rand(10+26*2)], 1..8;

   return $string;
}

# Sorts a given bam file to coordinate order
sub sort_bam($;$)
{
   my ($bam_file, $threads) = @_;

   my $bam_sort_prefix = "tmp" . random_string();
   my $tmp_bam = random_string() . "$bam_file";

   # Do multithreaded if specified
   my $sort_command = "samtools sort -O bam -o $tmp_bam -T $bam_sort_prefix ";
   if (defined($threads) && $threads > 1)
   {
      $sort_command .= "-@ $threads ";
   }
   $sort_command .= "$bam_file";

   system($sort_command);

   rename $tmp_bam, $bam_file;
}

# Run snap, producing sorted and indexed bam. Returns location of this bam
sub run_snap($$$$$)
{
   my ($reference_name, $sample_name, $forward_reads, $reverse_reads, $threads) = @_;

   # Create index for reference if required
   if (!-d "snap_index")
   {
      snap_index($reference_name);
   }

   my $output_name = "$sample_name.mapping.bam";
   my $output_name_tmp = $output_name . random_string();
   my $log_file = "$sample_name.mapping.log";

   my $snap_command = "$snap_location paired snap_index $forward_reads $reverse_reads -o $output_name -R '" . join('\t', '@RG', "ID:$sample_name", "PL:ILLUMINA", "SM:$sample_name");

   # Use multiple cores if available
   if ($threads > 1)
   {
      $snap_command .= " -t $threads -b";
   }

   $snap_command .=  "' -= &>> $log_file";

   system($snap_command);

   system("samtools fixmate -O bam $output_name $output_name_tmp");
   rename $output_name_tmp, $output_name;

   sort_bam($output_name, $threads);

   system("samtools index $output_name");

   return($output_name);
}

# Make hash for snap (~15x genome size)
sub snap_index($)
{
   my ($reference_name) = @_;

   if (!-d "snap_index")
   {
      mkdir "snap_index";
   }

   system("$snap_location index $reference_name snap_index -large 1>&2");
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

# Make bwa index
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

# Make smalt index/hash
sub smalt_index($)
{
   my ($reference_name) = @_;

   system("smalt index -k $smalt_k -s $smalt_s $reference_name $reference_name");
}

# Creates a .dict file for a fasta. Required for any picard or GATK modules to
# work with a fasta
sub create_fa_dict($)
{
   my ($fasta_file) = @_;

   $fasta_file =~ m/^(.+)\.(fasta|fa)$/;
   my $ref_name = $1;

   my $dict_command = "$java_location -Xmx100M -jar $picard_location CreateSequenceDictionary R=$fasta_file O=$ref_name.dict";
   system($dict_command);
}

# Uses picard to mark PCR duplicates (if any)
sub mark_dups($)
{
   my ($bam_file) = @_;

   my $dup_file = "$bam_file.dups";
   my $marked_bam = random_string() . "marked.$bam_file";

   my $picard_command = "$java_location -Xmx2g -jar $picard_location MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=$bam_file OUTPUT=$marked_bam METRICS_FILE=$dup_file";

   rename $marked_bam, $bam_file;
}

# Realigns bam files around indels
sub indel_realign($$$)
{
   my ($reference_file, $bam_file, $threads) = @_;

   if(!defined($threads) || $threads eq 0)
   {
      $threads = 1;
   }

   $reference_file =~ m/^(.+)\.(fasta|fa)$/;
   my $ref_name = $1;

   # Create dict if it doesn't already exist
   if (!-e "$ref_name.dict")
   {
      create_fa_dict($reference_file);
   }
   # Create fasta index if it doesn't already exist
   if (!-e "$reference_file.fai")
   {
      system("samtools faidx $reference_file");
   }

   # Create targets for realignment
   my $interval_command = "$java_location -Xmx2g -jar $gatk_location -T RealignerTargetCreator -R $reference_file -I $bam_file -o $bam_file.intervals --num_threads $threads";
   system($interval_command);

   # Actually do realignment
   my $realigned_bam = random_string() . "realigned.$bam_file";
   my $realign_command = "$java_location -Xmx3g -jar $gatk_location -T IndelRealigner -R $reference_file -I $bam_file -targetIntervals $bam_file.intervals -o $realigned_bam";
   system($realign_command);

   # Overwrite with output files
   $realigned_bam =~ m/(.+)\.bam$/;
   my $realigned_bam_index = "$1.bai";

   rename $realigned_bam, $bam_file;
   rename $realigned_bam_index, "$bam_file.bai";
}

1;
