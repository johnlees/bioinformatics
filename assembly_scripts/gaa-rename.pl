#!/usr/bin/perl -w
#
use strict;

my $input_file = $ARGV[0];
my $output_file = "gaa.$input_file";

open (IN, $input_file) || die("Couldn't open $input_file");
open (OUT, ">$output_file") || die("Couldn't open $output_file");

my $i = 1;

while (my $line_in = <IN>)
{
   if ($line_in =~ /^>/)
   {
      print OUT ">Contig1.$i\n";
      $i++;
   }
   else
   {
      print OUT $line_in;
   }
}

close IN;
close OUT;

exit(0);
