#!/usr/bin/perl -w
#
use strict;
use warnings;

my $chr = $ARGV[0];

my $output_file = "cases-ALS-BPROOF.impute2.chr$chr.gen";
my @files = glob("cases-ALS-BPROOF.impute2.$chr.*");

my %file_sort;

foreach my $file (@files)
{
   if ($file =~ /^cases-ALS-BPROOF\.impute2\..+\.(\d+)$/)
   {
      $file_sort{$1} = $file;
   }
   else
   {
      die("Could not parse $file\n");
   }
}

foreach my $file_number (sort { $a <=> $b } keys %file_sort)
{
   system("cat $file_sort{$file_number} >> $output_file");
   unlink($file_sort{$file_number});
}

system("gzip $output_file");

exit(0);

