#!perl -w

use strict;
use warnings;

use Getopt::Long;

my %gene_counts;

my $usage_message = <<USAGE;
Usage: ./summarise_output.pl <--genes|--orfs> --assemblies > result.txt 2> failed.txt

Gather the output of pairs_pipe.pl

   --samples         Sample, lane 1, lane 2. Same file as input to pairs_pipe.pl

   Mode options
   <none>            Number of mutations by sample
   --genes           Number of mutations by gene (count each sample once)
   --orfs            Names of all genes, to investigate unannotated proteins

   -h, --help        Shows this help.

USAGE

my ($pairs_file, $gene_mode, $orf_mode, $help);
GetOptions("samples=s" => \$pairs_file,
           "genes" => \$gene_mode,
           "orfs" => \$orf_mode,
           "help|h" => \$help) || die("$!\n$usage_message");

if (defined($help) || !defined($pairs_file) || (defined($gene_mode) && defined($orf_mode)))
{
   print STDERR $usage_message;
   exit 0;
}

open(PAIRS, $pairs_file) || die("Could not open $pairs_file\n");

if (!$gene_mode && !$orf_mode)
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

         my ($orfs, @orf_list);
         if ($orf_mode)
         {
            my $orfs = `bcftools query -f '%ANNOT_ID\n' $vcf_name`;
            my @orf_list = split("\n", $orfs);
         }

         for (my $i = 0; $i < scalar(@gene_list); $i++)
         {
            if ($gene_list[$i] =~ /^(.+)_(\d+)$/)
            {
               $gene_list[$i] = $1;
            }
            elsif ($orf_mode && $gene_list[$i] eq "1")
            {
               $gene_list[$i] = "ORF_" . $orf_list[$i];
            }
         }

         if (!$gene_mode && !$orf_mode)
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

if ($gene_mode || $orf_mode)
{
   print join("\t", "Gene name", "Samples with variation") . "\n";

   foreach my $gene (sort keys %gene_counts)
   {
      print join("\t", $gene, $gene_counts{$gene}) . "\n";
   }
}


exit(0);

