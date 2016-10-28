#!/usr/bin/perl -w

use strict;
use warnings;

# Converting sanger genotype data
# Usage: ./sanger_fcr_to_opticall.pl humanomniexpress-24-v1-1-a.bpm.csv
#

my $manifest_file = $ARGV[0];

# Read manifest
my $snp_regex = qr/^\[(.)\/(.)\]$/;

my @sample_names;
open(SAMPLES, "sample_names.txt") || die("Could not open sample_names.txt\n");
while (my $line_in = <SAMPLES>)
{
   chomp $line_in;
   push(@sample_names, $line_in);
}
close SAMPLES;

open(MANIFEST, $manifest_file) || die("Could not open manifest file $manifest_file: $!\n");
my $header = <MANIFEST>;

# Print headers
my @file_h;
for (my $i = 1; $i <= 22; $i++)
{
   open($file_h[$i], ">", "extmeng2.chr$i.int.txt") || die("Could not open output files $i: $!\n");
   print { $file_h[$i] } join("\t", "SNP", "Coor", "Alleles");
   foreach my $sample (@sample_names)
   {
      print { $file_h[$i] } "\t" . join("\t", $sample . "A", $sample . "B");
   }
   print { $file_h[$i] }  "\n";
}

# Open files
my @in_file_h;
for (my $i = 1; $i <= 904; $i++)
{
   open($in_file_h[$i], "$i.int") || die("Could not read input file $i: $!\n");
}

while (my $line_in = <MANIFEST>)
{
   chomp $line_in;

   my @fields = split(",", $line_in);
   $fields[5] =~ $snp_regex;

   unless ($fields[2] eq "X" || $fields[2] eq "Y" || $fields[2] eq "XY" || $fields[2] eq "0")
   {
      print { $file_h[$fields[2]] } join("\t", $fields[1], $fields[3], "$1$2") ;

      for (my $i = 1; $i <= 904; $i++)
      {
         my $sample_line = readline($in_file_h[$i]);
         chomp $sample_line;
         my @sample_fields = split("\t", $sample_line);

         print { $file_h[$fields[2]] } "\t" . join("\t", @sample_fields[1 .. 2]);
      }
      print { $file_h[$fields[2]] } "\n";
   }
}
close MANIFEST;

exit(0);

