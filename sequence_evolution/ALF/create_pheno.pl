#!/usr/bin/perl -w

# Creates a pheno file with the specified case ratio and OR, given a snp
# position
#
# ./create_pheno.pl msa.fa.vcf pos ratio OR

use strict;
use warnings;

my $vcf_in = $ARGV[0];
my $pos_in = $ARGV[1];
my $ratio = $ARGV[2];
my $or = $ARGV[3];

my ($ref_prob, $alt_prob);
if ($or > 0)
{
   $ref_prob = $ratio;
   $alt_prob = $ratio * $or;
}
else
{
   $ref_prob = 0;
   $alt_prob = 1;
}

open(VCF, $vcf_in) || die("Could not open $vcf_in\n");

my @sample_names;
while (my $vcf_line = <VCF>)
{
   if ($vcf_line !~ /^##/)
   {
      chomp $vcf_line;
      my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @vals) = split("\t", $vcf_line);

      if ($vcf_line =~ /^#/)
      {
         foreach my $sample (@vals)
         {
            $sample =~ m/^SE(\d+)\/(\d+)$/;
            push(@sample_names, "SE$1");
         }
      }
      elsif ($pos_in == $pos)
      {
         my $i = 0;
         foreach my $sample (@sample_names)
         {
            my $pheno = 0;
            if ($vals[$i] eq ".")
            {
               if (rand() < $ref_prob)
               {
                  $pheno = 1;
               }
            }
            else
            {
               if (rand() < $alt_prob)
               {
                  $pheno = 1;
               }
            }

            print join("\t", $sample, $sample, $pheno) . "\n";

            $i++;
         }

         last;
      }
   }
}

exit(0);

