#!/usr/bin/perl -w

use strict;
use warnings;

my %names;
open (IDATS, "sample_array_names.txt") || die ("Couldn't open sample_array_names.txt");

# Get reference of array number to sample name from idat naming scheme
while (my $idat = <IDATS>)
{
   chomp ($idat);
   my ($array_name, $sample_name);

   if ($idat =~ /^(.+_.+)_(ALS.+)_Grn\.idat$/)
   {
      $array_name = $1;
      $sample_name = $2;
   }

   $names{$array_name} = $sample_name;
}

close IDATS;

# Get genders from info file
my %genders;
open (INFO, "case-ALS.info") || die("Couldn't open case-ALS.info");

while (my $info_line = <INFO>)
{
   chomp($info_line);

   my ($sample_name, $gender, $exclusion, $group) = split(/\s+/, $info_line);
   $genders{$sample_name} = $gender;
}

close INFO;

# Fix tfam
open (TFAM_IN, "ALS.0.tfam") || die("Couldn't open ALS.0.tfam");
open (TFAM_OUT, ">new.tfam") || die("Couldn't open new.tfam for writing");

while (my $tfam_line = <TFAM_IN>)
{
   chomp($tfam_line);
   my ($array_name1, $array_name2, $dad, $mum, $pheno) = split(/\s+/, $tfam_line);

   my $sample_name = $names{$array_name1};
   my $gender = $genders{$sample_name};
   my $new_line = join(" ", $sample_name, $sample_name, 0, 0, $gender, $pheno) . "\n";

   print TFAM_OUT $new_line;
}

close TFAM_OUT;
close TFAM_IN;

exit(0);

