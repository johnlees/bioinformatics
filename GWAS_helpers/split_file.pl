#!/usr/bin/perl -w

use strict;
use warnings;

open(FCR, "omnix_extmeng_20160906.fcr.txt") || die("Could not open fcr file\n");
for (my $i = 0; $i < 10; $i++)
{
   my $header = <FCR>;
}

my $line_in = <FCR>;
chomp $line_in;
my @fields = split("\t", $line_in);
my $last_sample = $fields[1];

open(SAMPLE, ">$last_sample.int") || die("Couldn't open sample file\n");
print SAMPLE join("\t", $fields[0], $fields[7], $fields[8]) . "\n";
print STDERR $last_sample . "\n";

while ($line_in = <FCR>)
{
   chomp $line_in;
   @fields = split("\t", $line_in);

   if ($fields[1] ne $last_sample)
   {
      $last_sample = $fields[1];
      print STDERR $last_sample . "\n";

      close SAMPLE;
      open(SAMPLE, ">$last_sample.int") || die("Couldn't open sample file\n");
   }
   print SAMPLE join("\t", $fields[0], $fields[7], $fields[8]) . "\n";
}

close SAMPLE;

exit(0);

