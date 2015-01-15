#!perl -w

use strict;
use warnings;

my $pairs_file = $ARGV[0];

open(PAIRS, $pairs_file) || die("Could not open $pairs_file\n");

print join("\t", "Sample", "SNPs", "INDELs", "Total", "Genes") . "\n";

my $header = <PAIRS>;
while (my $line_in = <PAIRS>)
{
   chomp $line_in;

   my ($sample, $csf, $blood) = split("\t", $line_in);

   if (-d $sample)
   {
      chdir $sample;

      my $vcf_name = "$sample.vcf.gz";

      if (-e $vcf_name)
      {
         my $indels = `bcftools view -H -v indels $vcf_name | wc -l`;
         chomp($indels);

         my $snps = `bcftools view -H -v snps $vcf_name | wc -l`;
         chomp($snps);

         my $total = $snps + $indels;

         my $genes = `bcftools query -f '%GENE\n' $vcf_name`;
         my @gene_list = split("\n", $genes);

         pop(@gene_list);

         print join("\t", $sample, $snps, $indels, $total, @gene_list) . "\n";
      }
      else
      {
         print STDERR "$sample has no vcf\n";
      }

      chdir "..";
   }
   else
   {
      print STDERR "$sample directory not found\n";
   }
}

close PAIRS;

exit(0);

