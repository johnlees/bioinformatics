#!/software/bin/perl -w
#
use strict;
use warnings;

my $file_prefix = "cases-ALS-BPROOF.impute2";
my $chromosome = $ARGV[0];

my @files = glob("$file_prefix.$chromosome.*.gz");

my %file_hash;

foreach my $file (@files)
{
   if ($file =~ /\.(\d+)\.gz$/)
   {
      $file_hash{$1} = $file;
   }
}

my $i = 1;
my $j = 2;
foreach my $file_num (sort {$a <=> $b} keys %file_hash)
{
   my $new_name = "$file_prefix.$chromosome.part$j.$i.gz";
   rename $file_hash{$file_num}, $new_name;

   if ($i == 22)
   {
      $i = 1;
      $j++;
   }
   else
   {
      $i++;
   }
}

exit(0);

