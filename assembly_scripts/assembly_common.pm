#!/usr/bin/perl -w

package assembly_common;

use strict;
use warnings;

use Bio::SeqIO;
use POSIX;

use File::Spec;
use File::Path qw(remove_tree);

our @tmp_file_list;

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

# Renames contigs from the assembly pipeline to match the annotation pipeline
# Overwrites file in place if new file name not provided
sub standardise_contig_names($$)
{
   my ($contigs_file, $new_file) = @_;

   my $tmp_file = "tmp_contigs_rename.tmp";

   if(!defined($contigs_file) || !-e($contigs_file))
   {
      die("Contigs file $contigs_file doesn't exist or can't be read: $!\n");
   }

   # Open multi-fasta objects
   my $multi_fasta = Bio::SeqIO->new(-file => $contigs_file, -format => 'fasta') || die("Failed to open $contigs_file: $!\n");
   my $renamed_out;

   if(!defined($new_file))
   {
      $renamed_out = Bio::SeqIO->new(-file => ">$tmp_file", -format => 'fasta') || die("Couldn't write to $tmp_file: $!\n");
   }
   else
   {
      $renamed_out = Bio::SeqIO->new(-file => ">$new_file", -format => 'fasta') || die("Couldn't write to $tmp_file: $!\n");
   }

   while (my $sequence = $multi_fasta->next_seq())
   {
      # Contigs are named .run_lane_tag_contig
      # Rename as contig000001 (uses 6 spaces, so 2 digit contig numbers have
      # 4 zeros rather than 5 - hence use sprintf)
      my $new_id;

      if ($sequence->display_id =~ m/^\.\d+_\d+_\d+\.(\d+)$/)
      {
         # Velvet header
         $new_id = "contig" . sprintf("%06d", $1);
      }
      elsif ($sequence->display_id =~ m/^contig(\d+)$/)
      {
         # SPAdes header
         $new_id = "contig" . sprintf("%06d", $1);
      }

      # Write new name to tmp file
      $sequence->display_id($new_id);
      $renamed_out->write_seq($sequence);
   }

   # Overwrite original if necessary
   if (!defined($new_file))
   {
      rename $tmp_file, $contigs_file;
      unlink $tmp_file;
   }
}

# Takes a sequence as a string, and changes some letters to create SNPs
# Mutation rate is indepenedent of context, and the same rate between all bases
# i.e. JC69 model
#
# Returns sequence and an array (reference) of mutations and locations
sub create_snps($$$)
{
   my ($sequence, $rate, $seed) = @_;

   my @alphabet = ("a", "c", "t", "g");
   my @mutations;

   my $sequence_length = length($sequence);

   my $num_mutations = $sequence_length * $rate;
   print STDERR "Seq length: $sequence_length expected mutations: $num_mutations\n";

   # Use the seed passed if it exists to allow replication of results
   if (defined($seed))
   {
      srand($seed);
   }

   for (my $i = 0; $i< $sequence_length; $i++)
   {
      if (rand($sequence_length) <= $num_mutations)
      {
         # Randomly pick a new base ~ U(0-3)
         my $random_num = floor(rand(4));
         if ($random_num == 4)
         {
            $random_num = 3;
         }

         my $new_base = $alphabet[$random_num];
         my $base = substr($sequence, $i, 1, $new_base);

         # May be the same as old base, but if not record the mutation
         unless ($new_base eq $base)
         {
            push(@mutations, $i+1 . ": $base" . "->$new_base");
         }

      }
   }

   return($sequence, \@mutations);
}

# Adds given file to an array of files to delete
sub add_tmp_file($$)
{
   my ($tmp_file, $file_array) = @_;

   # Get absolute path
   my $abs_path = File::Spec->rel2abs($tmp_file);

   push(@$file_array, $abs_path);
}

# Adds an array of temporary files to be deleted
sub add_tmp_files($$)
{
   my ($tmp_files, $tmp_file_array) = @_;

   foreach my $tmp_file (@$tmp_files)
   {
      add_tmp_file($tmp_file, $tmp_file_array);
   }
}

# Cleans up temporary files
sub clean_up($)
{
   my ($tmp_files) = @_;

   foreach my $file (@$tmp_files)
   {
      if (-d $file)
      {
         # Remove recursively
         print STDERR "rm -rf $file\n";
         #remove_tree($file);
      }
      elsif (-e $file)
      {
         print STDERR "rm $file\n";
         #unlink $file;
      }
      else
      {
         print STDERR "cannot remove temporary file $file\n";
      }
   }
}

1;
