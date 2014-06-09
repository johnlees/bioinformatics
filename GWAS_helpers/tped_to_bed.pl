#!/usr/bin/perl -w

use strict;
use warnings;

my @tpeds = glob("*.tped");

foreach my $tped (@tpeds)
{
   my $prefix;

   if ($tped =~ /^(.+)\.tped$/)
   {
      $prefix = $1;
   }

   my $command = "plink --noweb --allow-no-sex --tfile $prefix --out $prefix --make-bed";
   system($command);
}

exit(0);

