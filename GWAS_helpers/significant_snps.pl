#!/usr/bin/perl -w
#
# Gives a list of all snps that pass genome wide significance
# Input is in qqman (R package) format
#
use strict;
use warnings;

# Set p value for significance
my $significance = 5e-8;

# Open input and output files
my $file = $ARGV[0];
my $out_file = $ARGV[0] . ".significants";

open (ASSOC, $file) || die("Could not open $file: $!\n");
open (OUT, ">$out_file") || die("Could not open $out_file for writing: $!\n");

my $header = <ASSOC>;
print OUT $header;

while (my $row = <ASSOC>)
{
   chomp($row);
   my ($rsid, $chr, $bp, $p) = split(/\s+/, $row);

   if ($p ne "NA" && $p <= $significance)
   {
      print OUT "$row\n";
   }
}

close ASSOC;
close OUT;

exit(0);

