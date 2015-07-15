#!/usr/bin/perl -w

use strict;
use warnings;

print STDERR "Usage: ./add_metadata_column.pl name hpgen_ids -9 phenotypes.txt > new_phenotypes.txt\n";

my $title = $ARGV[0];
my $in_file = $ARGV[1];
my $missing_val = $ARGV[2];
my $out_file = $ARGV[3];

open(IN, $in_file) || die("Could not open $in_file\n");

my @new_ids;
while (my $line_in = <IN>)
{
   chomp $line_in;

   push(@new_ids, $line_in);
}

close IN;

open(PHENO, $out_file) || die("Could not open $out_file\n");

my $header = <PHENO>;
chomp $header;

print "$header $title\n";

while (my $pheno_line = <PHENO>)
{
   chomp $pheno_line;

   my ($fid, $iid, @columns) = split(/\s+/, $pheno_line);

   if ($iid =~ m/^(.+)_(.+)_(hpgen\d+)$/)
   {
      my $found = 0;
      foreach my $new_id (@new_ids)
      {
         if ($new_id eq $3)
         {
            print join(" ", $fid, $iid, @columns, "1") . "\n";
            $found = 1;
            last;
         }
      }

      if (!$found)
      {
         print join(" ", $fid, $iid, @columns, $missing_val) . "\n";
      }
   }
   else
   {
      print join(" ", $pheno_line, "0") . "\n";
   }
}

close PHENO;

exit(0);

