#!/usr/bin/perl -w
#
use strict;

# First run
# cut -d "," -f 2,21 humanomniexpress-12v1-1_b.csv > humanomniexpress-12v1-1_b.strands
#
# run this as ./strand_flip_list.pl > snplist.txt
#

open (STRANDS, "humanomniexpress-12v1-1_b.strands") || die("Could not open strand file");

while (my $strand_line = <STRANDS>)
{
   chomp($strand_line);

   my($rsid, $strand) = split(",", $strand_line);

   if ($strand eq "-")
   {
      print "$rsid\n"
   }
}

close STRANDS;

exit(0);
