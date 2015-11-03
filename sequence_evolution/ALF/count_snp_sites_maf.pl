#!/usr/bin/perl -w

# Counts MAF from vcfs produced by snp_sites
# alignment -> snp_sites -> this script
# snp_sites -v MSA_88_dna.fa
# ./count_snp_sites_maf.pl MSA_88_dna.fa.vcf

use strict;
use warnings;

my $vcf_in = $ARGV[0];

open(VCF, $vcf_in) || die("Could not open $vcf_in\n");

while (my $vcf_line = <VCF>)
{
   if ($vcf_line !~ /^#/)
   {
      chomp $vcf_line;

      my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @vals) = split("\t", $vcf_line);

      my $mac = 0;
      foreach my $val (@vals)
      {
         if ($val ne "0")
         {
            $mac++
         }
      }

      my $maf = $mac/scalar(@vals);

      print join("\t", $pos, $ref, $alt, $maf) . "\n";
   }
}

exit(0);

