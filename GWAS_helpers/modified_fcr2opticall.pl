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
open(SAMPLES, "../sample_names_new.txt") || die("Could not open ../sample_names_new.txt\n");
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
my %in_file_h;
foreach my $sample_name (@sample_names)
{
   open($in_file_h{$sample_name}, "$sample_name.int") || die("Could not read input file $sample_name.int: $!\n");
}

while (my $line_in = <MANIFEST>)
{
   chomp $line_in;

   my @fields = split(",", $line_in);
   $fields[5] =~ $snp_regex;
   my $chr = $fields[2];

   if ($chr eq "X" || $chr eq "Y" || $chr eq "XY" || $chr eq "0" || $chr > 22)
   {
      # read but don't do anything
      foreach my $sample_name (@sample_names)
      {
         my $sample_line = readline($in_file_h{$sample_name});
      }
   }
   else
   {
      print { $file_h[$chr] } join("\t", $fields[1], $fields[3], "$1$2") ;

      foreach my $sample_name (@sample_names)
      {
         my $sample_line = readline($in_file_h{$sample_name});
         chomp $sample_line;
         my @sample_fields = split("\t", $sample_line);

         if ($sample_fields[0] ne $fields[1])
         {
            print STDERR $sample_fields[0] . " ne $fields[1] in $sample_name\n";
            die;
         }

         print { $file_h[$chr] } "\t" . join("\t", @sample_fields[1 .. 2]);
      }
      print { $file_h[$chr] } "\n";
   }
}
close MANIFEST;

exit(0);

