#!/usr/bin/perl -w

# This converts a MSA processed by snp sites into a vcf, into a list of samples
# and whether or not they have a specified snp
# Usage:
# ./snp_sites_to_pheno.pl MSA_78_dna.fa.vcf 875 0
# 875 snp position
# 0 reference coding

use warnings;
use strict;

my $file_in = $ARGV[0];
my $snp_position = $ARGV[1];
my $ref_coding = $ARGV[2];

open(VCF, $file_in) || die("Could not open $file_in\n");

my (@genes, @vals);
while (my $line_in = <VCF>)
{
   my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format);

   chomp $line_in;
   if ($line_in =~ /^#CHROM/)
   {
      ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genes) = split("\t", $line_in);
   }
   elsif ($line_in =~ /^##/)
   {
      next;
   }
   else
   {
      ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @vals) = split("\t", $line_in);
      if ($pos == $snp_position)
      {
         my $i = 0;
         my $last_sample = "";
         my $last_present = 0;
         foreach my $gene (@genes)
         {
            $gene =~ m/^(SE\d+)\/(\d+)$/;
            my $sample = $1;
            if ($last_sample ne $sample)
            {
               if ($last_sample ne "")
               {
                  print join("\t", $last_sample, $last_present) . "\n";
               }
               $last_present = 0;
            }

            if ($vals[$i] != $ref_coding)
            {
               $last_present = 1;
            }

            $last_sample = $sample;
            $i++;
         }
         print join("\t", $last_sample, $last_present) . "\n";

         last;
      }
   }
}

close VCF;

exit(0);

