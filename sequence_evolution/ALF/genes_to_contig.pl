#!/usr/bin/perl -w

use strict;
use warnings;

# Converts a multicontig file into one with just one contig, for use as input
# as an alignment file.
# Usage ./genes_to_contig.pl input_files.txt

use Bio::SeqIO;

my $input_files = $ARGV[0];

open(LIST, $input_files) || die("Could not open $input_files\n");

while (my $fasta_file = <LIST>)
{
   chomp $fasta_file;

   my $contigs_in = Bio::SeqIO->new(-file => "$fasta_file",
                                    -format => 'Fasta');

   my $contig = "";
   while (my $fasta_seq = $contigs_in->next_seq())
   {
      $contig .= $fasta_seq->seq();
   }

   my $contigs_out = Bio::SeqIO->new(-file => ">$fasta_file",
                                     -format => 'Fasta');
   my $seq_out = Bio::Seq->new(-display_id => "$fasta_file",
                              -seq => $contig);

   $contigs_out->write_seq($seq_out);

   $contigs_out->close();
}

close LIST;

exit(0);


