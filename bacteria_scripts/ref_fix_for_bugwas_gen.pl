#!/usr/bin/perl -w

use strict;
use warnings;

my $file_in = $ARGV[0];

open(IN, $file_in) || die("Could not open $file_in\n");
while (my $line_in = <IN>)
{
   chomp $line_in;

   if ($line_in =~ /^1/)
   {
      my @fields = split("\t", $line_in);
      my @alts = split(",", $fields[4]);
      my $i = 0;
      my $rep = "NA";
      foreach my $alt (@alts)
      {
         if ($alt eq "*")
         {
            $rep = $i;
            last;
         }
         $i++;
      }

      if ($rep ne "NA")
      {
         my $gts = join("\t", @fields[10 .. $#fields]);
         $gts =~ s/$rep/0/g;
         $line_in = join("\t", @fields[0 .. 9], $gts);
      }
   }
   print $line_in . "\n";
}

exit(0);


