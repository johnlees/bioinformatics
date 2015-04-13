#!/usr/bin/perl -w

use strict;
use warnings;

my %all_consq;

my $sample_file = $ARGV[0];
open(SAMPLES, $sample_file) || die("Couldn't open $sample_file\n");

while (my $line_in = <SAMPLES>)
{
   chomp $line_in;

   my $dlt_vars = `grep dlt $line_in/$line_in.vep`;
   my @dlt_lines = split("\n", $dlt_vars);

   foreach my $variant (@dlt_lines)
   {
      my ($id, $pos, $allele, $gene_id, $trans_id, $loc, $consequences, $pos_start,
         $pos_end, $aa_pos, $aa, $codon, $existing, $extra) = split("\t", $variant);

      my @variant_consequences = split(",", $consequences);

      foreach my $consequence (@variant_consequences)
      {
         $all_consq{$consequence}++;
      }
   }
}

close SAMPLES;

foreach my $consequence (sort keys %all_consq)
{
   print join("\t", $consequence, $all_consq{$consequence}) . "\n";
}

exit(0);

