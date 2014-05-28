#!/usr/bin/perl -w

use strict;
use warnings;

my $sample_file_name = "case_sample_names.txt";
my $info_file = "case-ALS.info";
my $gender_file = "case_gender.txt";
my $group = 0;

# Read sample IDs
my @samples;

open (SAMPLES, $sample_file_name) || die ("Could not open $sample_file_name");
while (my $sample_line = <SAMPLES>)
{
   chomp($sample_line);
   push @samples,$sample_line;
}
close SAMPLES;

# Get genders
my %genders;

open (GENDERS, $gender_file) || die ("Could not open $gender_file");
while (my $gender_line = <GENDERS>)
{
   chomp($gender_line);
   my ($sample_id, $gender) = split(/\s+/, $gender_line);

   $genders{$sample_id} = $gender;
}
close GENDERS;

# Write out info file
open (INFO, ">$info_file") || die ("Could not open $info_file");

foreach my $sample_id (@samples)
{
   my $gender;

   if (defined($genders{$sample_id}))
   {
      $gender = $genders{$sample_id};
   }
   else
   {
      $gender = 0;
   }

   print INFO "$sample_id $gender 0 $group\n";
}
close INFO;

exit(0);
