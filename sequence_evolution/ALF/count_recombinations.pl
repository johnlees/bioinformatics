#!/usr/bin/perl -w

use strict;
use warnings;

my $input_file = $ARGV[0];
my $core_name = $ARGV[1];

my $core = 0;
my $recomb = 0;

open(INPUT, $input_file) || die("Could not open $input_file\n");
while (my $line_in = <INPUT>)
{
   chomp $line_in;

   if ($line_in =~ /^>SE\d+\/0*(\d+)$/)
   {
      if ($1 eq $core_name)
      {
         $core++;
      }
      else
      {
         $recomb++;
      }
   }
}

close INPUT;

print join("\t", $core_name, $core, $recomb) . "\n";

exit(0);

