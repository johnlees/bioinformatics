#!/usr/bin/perl -w

package assembly_common;

use strict;
use warnings;

#
# Functions common to scripts dealing with assemblies and fastq files
#
#
# Note on fastq formats. 4 lines per read
# ReadID @machine_run:lane:(next three I think are coordinates of read in the
# library)#tag/direction
# Read sequence
# +
# Read quality (cigar string)
#

# Returns the effective coverage when reads have been split into kmers
sub effective_coverage($$$)
{
   my ($kmer, $read_length, $expected_coverage) = @_;

   if (!defined($read_length) || $read_length <= 0)
   {
      print STDERR "WARNING: Invalid read length, assume length of 100 instead\n";
      $read_length = 100;
   }

   # Formula from cortex supp mat
   my $effective_coverage = (($read_length - $kmer + 1)/$read_length) * $expected_coverage;

   return($effective_coverage);
}

# Extracts read length from a fastq file (assuming all are the same length)
sub read_length($)
{
   my ($fastq_file) = @_;

   open (FASTQ, "$fastq_file") || die("Could not open $fastq_file: $!\n");

   # Throw header away, store read
   my $header_line = <FASTQ>;
   my $read = <FASTQ>;

   close FASTQ;

   return(length($read));
}

# Gets the total number of reads in a fastq file
sub number_reads($)
{
   my ($fastq_file) = @_;

   # Quick and dirty
   my $command = 'grep -c "^+$" ' . $fastq_file;
   my $num_reads = `$command`;

   chomp($num_reads);

   return($num_reads);
}

# Gets expected coverage given number of reads, length and expected genome size
sub expected_coverage($$)
{
   my ($fastq_file, $genome_size) = @_;

   # Get total number of bases sequenced
   my $total_sequenced = number_reads($fastq_file) * read_length($fastq_file);

   if (!defined($genome_size) || $genome_size <= 0)
   {
      print STDERR "WARNING: Invalid genome size. Assuming 2Mbases\n";
      $genome_size = 2000000;
   }

   my $expected_coverage = $total_sequenced/$genome_size;

   return($expected_coverage);
}

1;
