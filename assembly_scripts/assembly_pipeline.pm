#!/usr/bin/perl -w

package assembly_pipeline;

use strict;
use warnings;

use File::Path qw(remove_tree);
use File::Copy;

sub run_spades($)
{
   
}

sub filter_contigs($$$$)
{
   my ($len_cutoff, $cov_cutoff, $input_file, $output_file) = @_;

   my $fasta_in = Bio::SeqIO->new(-file => $input_file,
                                  -format => 'fasta') || die($!);

   my $fasta_out = Bio::SeqIO->new(-file => ">$output_file",
                                  -format => 'fasta') || die($!);

   while (my $contig = $fasta_in->next_seq())
   {
      my ($length, $coverage);
      if ($contig->id =~ /^NODE_\d+_length_(\d+)_cov_(.+)_ID_\d+$/)
      {
         $length = $1;
         $coverage = $2;
      }

      if ($length >= $len_cutoff && $coverage >= $cov_cutoff)
      {
         $fasta_out->write_seq($contig);
      }
   }
}

sub run_improvement
{

}

sub run_annotation
{

}

1;
