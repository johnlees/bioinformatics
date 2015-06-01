#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::SeqIO;

my $seq_alignment_files = $ARGV[0];

open(LIST, $seq_alignment_files) || die("Could not open $seq_alignment_files\n");

my %out_seqs;
while (my $msa_file = <LIST>)
{
   chomp $msa_file;

   my $msa_in = Bio::SeqIO->new(-file => "$msa_file",
                                 -format => 'Fasta');

   while (my $msa_seq = $msa_in->next_seq())
   {
      $msa_seq->id() =~ m/^SE(\d+)\/(\d+)/;

      my $output_file = "organism_" . $1 . ".fa";
      if (!-e $output_file)
      {
         $out_seqs{$1} = Bio::SeqIO->new(-file => ">$output_file",
                                 -format => 'Fasta');
      }

      $out_seqs{$1}->write_seq($msa_seq);
   }

   $msa_in->close();
}

exit(0);

