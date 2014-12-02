#!/usr/bin/perl -w

#
# Runs compare_variants.pl on a batch of call sets, and extracts power and
# number of false positives
#

use strict;
use warnings;

my @types = ("snps", "indels");

print join("\t", "Run", "Type", "Power", "False positives");

for (my $i = 1; $i <= 100; $i++)
{
   chdir "R6_mutant$i";

   foreach my $type (@types)
   {
      my $compare_command = "~/bioinformatics/bacteria_scripts/compare_variants.pl --vcf2 R6_mut$i.diff.vcf.gz --ref2 ~/pathogen_lustre/pairs_analysis/simulated_data/mapping/R6_mutant1/reference.renamed.fa --vcf1 ~/pathogen_lustre/simulated_bacteria/sim_pneumo/mut$i.vcf --ref1 ~/pathogen_lustre/simulated_bacteria/Streptococcus_pneumoniae_R6.fa --top-hit --type $type | grep -v \"MISMATCH\" | wc -l";

      my $true_positives = `$compare_command`;
      my $total_positives = `bcftools view -v $type -H ~/pathogen_lustre/simulated_bacteria/sim_pneumo/mut$i.vcf 2> /dev/null | wc -l`;
      my $calls = `bcftools view -v $type -H R6_mut$i.diff.vcf.gz 2> /dev/null | wc -l`;

      my $power = $true_positives/$total_positives;
      my $false_positives = $calls - $total_positives;

      print join("\t", $i, $type, $power, "$false_positives\n");
   }

   chdir "..";
}


exit(0);

