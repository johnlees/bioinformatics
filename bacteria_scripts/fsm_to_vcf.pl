#!/usr/bin/perl -w

use strict;
use warnings;

# Converts kmers from fsm to a vcf file for use in plink

my $pheno_file = $ARGV[0];
my $kmer_file = $ARGV[1];
my $maf = $ARGV[2];

my $kmer_regex = qr/^(.+):\d+$/;

if ($maf eq "")
{
   $maf = 0.01;
}

open(PHENO, $pheno_file) || die("Could not open pheno file $pheno_file: $!\n");
my @samples;
while (my $line_in = <PHENO>)
{
   chomp $line_in;
   my ($fid, $iid, $pheno) = split("\t", $line_in);

   push(@samples, $iid);
}
close PHENO;

my $vcf_header = <<VCF;
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20160524
##source=PLINKv1.90
##contig=<ID=26,length=198249>
##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=FM211187>
##bcftools_viewVersion=1.2-5-g7fa0d25+htslib-1.2.1-23-g47a2046
##bcftools_viewCommand=view MA_snps.vcf.gz
VCF

$vcf_header .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" . join("\t", @samples) . "\n";

print $vcf_header;

open(KMER, $kmer_file) || die("Could not open pheno file $kmer_file: $!\n");
my $i = 1;
while (my $line_in = <KMER>)
{
   chomp $line_in;
   my ($kmer, $bar, @found) = split(" ", $line_in);

   if (scalar(@found)/scalar(@samples) < $maf || scalar(@found)/scalar(@samples) > (1-$maf))
   {
      next;
   }

   my %found_h;
   foreach my $found_s (@found)
   {
      $found_s =~ $kmer_regex;
      $found_h{$1} = 1;
   }

   print join("\t", 26, $i, ".", "-", $kmer, ".", ".", ".", "GT");
   $i++;

   foreach my $sample (@samples)
   {
      if (defined($found_h{$sample}))
      {
         print "\t1/1";
      }
      else
      {
         print "\t0/0";
      }
   }
   print "\n"
}

exit(0);

