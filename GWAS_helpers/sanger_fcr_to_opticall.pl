#!/usr/bin/perl -w

use strict;
use warnings;

# Converting sanger genotype data
# Usage: ./sanger_fcr_to_opticall.pl humanomniexpress-24-v1-1-a.bpm.csv omnix_extmeng_20160906.fcr.txt
#

my $manifest_file = $ARGV[0];
my $genotype_file = $ARGV[1];

# Read manifest
my $snp_regex = qr/^\[(.)\/(.)\]$/;

open(MANIFEST, $manifest_file) || die("Could not open manifest file $manifest_file: $!\n");
my $header = <MANIFEST>;

my %mapping;
while (my $line_in = <MANIFEST>)
{
   chomp $line_in;

   my @fields = split(",", $line_in);
   $mapping{$fields[1]}{chr} = $fields[2];
   $mapping{$fields[1]}{pos} = $fields[3];

   $fields[5] =~ $snp_regex;
   $mapping{$fields[1]}{A} = $1;
   $mapping{$fields[1]}{B} = $2;
}

close MANIFEST;

open(GENOTYPES, $genotype_file) || die("Could not open genotype file $genotype_file: $!\n");
for (my $i = 0; $i < 10; $i++)
{
   $header = <GENOTYPES>;
}

my $prev_sample = "";
my @sample_names;

# Read intensities
my %ints;
my $row = 1;
while (my $line_in = <GENOTYPES>)
{
   chomp $line_in;
   my (@fields) = split("\t", $line_in);

   # Check sample names
   if ($fields[1] ne $prev_sample)
   {
      $prev_sample = $fields[1];
      push(@sample_names, $prev_sample);
   }

   push(@{$ints{$fields[0]}}, $fields[7], $fields[8]);
   if ($row % 1000000 == 0)
   {
      print STDERR "Read $row rows\n";
   }
   $row++;
}

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

# Print to output
foreach my $rsid (keys %mapping)
{
   print { $file_h[$mapping{$rsid}{chr}] } join("\t", $rsid, $mapping{$rsid}{pos}, $mapping{$rsid}{A} . $mapping{$rsid}{B}, @{$ints{$rsid}}) . "\n";
}

exit(0);

