#!/usr/bin/perl -w

use strict;
use warnings;

open(FCR, "omnix_extmeng_20160906.fcr.txt") || die("Could not open fcr file\n");
for (my $i = 0; $i < 10; $i++)
{
   my $header = <FCR>;
}

for (my $i = 1; $i <= 904; $i++)
{
   open(SAMPLE, ">$i.int") || die("Could not open sample file\n");

   my @fields;
   for (my $j = 1; $j <= 713014; $j++)
   {
      my $line_in = <FCR>;
      chomp $line_in;
      @fields = split("\t", $line_in);
      print SAMPLE join("\t", $fields[0], $fields[7], $fields[8]) . "\n";
   }
   print STDERR $fields[1] . "\n";
   close SAMPLE;
}

exit(0);

