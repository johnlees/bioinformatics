#!perl -w

use strict;
use warnings;

my %gene_counts;
my $mode = $ARGV[0];
my $pairs_file = $ARGV[1];

open(PAIRS, $pairs_file) || die("Could not open $pairs_file\n");

unless ($mode = "genes")
{
   print join("\t", "Sample", "SNPs", "INDELs", "Total", "Genes") . "\n";
}

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

         unless ($mode = "genes")
         {
            print join("\t", $sample, $snps, $indels, $total, @gene_list) . "\n";
         }

         my %counted;
         foreach my $gene (@gene_list)
         {
            unless(defined($counted{$gene}))
            {
               $gene_counts{$gene}++;
               $counted{$gene} = 1;
            }
         }

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

if ($mode eq "genes")
{
   print join("\t", "Gene name", "Samples with variation") . "\n";

   foreach my $gene (sort keys %gene_counts)
   {
      print join("\t", $gene, $gene_counts{$gene}) . "\n";
   }
}


exit(0);

