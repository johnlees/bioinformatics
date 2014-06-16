#!/usr/bin/perl -w

use strict;
use warnings;

my $file_prefix = "cases-ALS-BPROOF";

my @het_range = (0.2526039, 0.3411034);
my $missing_limit = 0.03;

# First remove heterozygosity outliers
open (HET, "$file_prefix.het") || die ("Couldn't open $file_prefix.het: $!");
my $header = <HET>;

open (FAILHET, ">fail_het_qc.txt") || die ("Couldn't write to fail_het_qc.txt: $!");

while (my $het_line = <HET>)
{
   chomp($het_line);
   my ($blank, $FID, $IID, $OHOM, $EHOM, $NNM, $F) = split(/\s+/, $het_line);

   my $het_rate = ($NNM - $OHOM)/$NNM;

   if ($het_rate > $het_range[1] || $het_rate < $het_range[0])
   {
      print FAILHET "$FID\t$IID\n";
   }
}

close HET;
close FAILHET;

open (MISS, "$file_prefix.imiss") || die ("Couldn't open $file_prefix.imiss: $!");
$header = <MISS>;

open (FAILMISS, ">fail_imiss_qc.txt") || die ("Couldn't write to fail_imiss_qc.txt: $!");

while (my $miss_line = <MISS>)
{
   chomp($miss_line);
   my ($blank, $FID, $IID, $MISS_PHENO, $N_MISS, $N_GENO, $F_MISS) = split(/\s+/, $miss_line);

   if ($F_MISS > $missing_limit)
   {
      print FAILMISS "$FID\t$IID\n";
   }
}

close MISS;
close FAILMISS;

exit(0);

