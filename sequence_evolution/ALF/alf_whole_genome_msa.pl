#!/usr/bin/perl -w

# Concatenates MSA output files which are by gene, into a full alignment
# usage: ./alf_msa_concat.pl msa_files.txt samples.txt

use strict;
use warnings;

use Bio::SeqIO;

# Read list of MSAs from input
my $seq_alignment_files = $ARGV[0];
my $sample_file = $ARGV[1];

# Read in sample list
open(SAMPLES, $sample_file) || die("Could not open $sample_file\n");
my @samples;
while (my $sample_line = <SAMPLES>)
{
   chomp $sample_line;
   push (@samples, $sample_line);
}
close SAMPLES;

# Read in all MSAs into a hash
open(LIST, $seq_alignment_files) || die("Could not open $seq_alignment_files\n");

my (%out_seqs, %seq_lengths);
my $largest_seq = 0;

while (my $msa_file = <LIST>)
{
   chomp $msa_file;

   my $msa_in = Bio::SeqIO->new(-file => "$msa_file",
                                 -format => 'Fasta');

   while (my $msa_seq = $msa_in->next_seq())
   {
      $msa_seq->id() =~ m/^(SE\d+)\/0+(\d+)/; # $1 sample, $2 sequence

      $out_seqs{$1}{$2} = $msa_seq->seq();

      # Keep track of number of sequences, and their lengths
      if ($2 > $largest_seq)
      {
         $largest_seq = $2;
      }

      if (!defined($seq_lengths{$2}))
      {
         $seq_lengths{$2} = length($msa_seq->seq());
      }
   }

   $msa_in->close();
}

close LIST;

# Write out. Loop over samples, then over sequences
my $fasta_out = Bio::SeqIO->new(-file => ">msa_all.fa",
                -format => 'Fasta');

foreach my $sample (@samples)
{
   my $sequence_out;
   for (my $sequence_number = 0; $sequence_number <= $largest_seq; $sequence_number++)
   {
      if (defined($out_seqs{$sample}{$sequence_number}))
      {
         $sequence_out .= $out_seqs{$sample}{$sequence_number};
      }
      else
      {
         $sequence_out .= "-" x $seq_lengths{$sequence_number};
      }
   }

   $fasta_out->write_seq(Bio::Seq->new(-display_id => $sample, -seq =>$sequence_out));
}

$fasta_out->close();

exit(0);

