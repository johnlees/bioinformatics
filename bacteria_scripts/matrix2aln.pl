#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::SeqIO;

my $file_in = $ARGV[0];

open(IN, $file_in) || die("Couldn't open $file_in: $!\n");

my %sequences;
while (my $line_in = <IN>)
{
   chomp $line_in;

   my @fields = split("\t", $line_in);

   foreach my $field (@fields)
   {
      my ($sample, $gt) = split("=", $field);

      if ($gt =~ /^(.)\/.$/)
      {
         $gt = $1;
      }

      $sequences{$sample} .= $gt;
   }

}

close IN;

my $msa_out = Bio::SeqIO->new(-file => ">$file_in.aln",
                             -format => 'Fasta');

foreach my $sample (sort keys %sequences)
{
   my $seq_out = Bio::Seq->new(-id => $sample,
                               -seq => $sequences{$sample});
   $msa_out->write_seq($seq_out);
}

$msa_out->close();

exit(0);

