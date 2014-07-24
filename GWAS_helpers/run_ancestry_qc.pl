#!/usr/bin/perl -w

use strict;
use warnings;

my $file_name = "cases-ALS-BPROOF.hapmap3r2.pruned.pca.evec";

my $pc2_lower_limit = 0.068;

open (PCA, "$file_name") || die ("Couldn't open $file_name: $!");
my $header = <PCA>;

open (FAILPCA, ">fail_ancestry_qc.txt") || die ("Couldn't write to fail_imiss_qc.txt: $!");

while (my $pca_line = <PCA>)
{
   chomp($pca_line);
   my @fields = split(/\s+/, $pca_line);

   my $group = pop(@fields);
   my $pc2 = pop(@fields);
   my $pc1 = pop(@fields);
   my $id = pop(@fields);

   if (($id =~ /^(.+):(.+)$/) && ($1 eq $2))
   {
      $id = $1;
   }

   if ($pc2 < $pc2_lower_limit)
   {
      if ($group eq "Control")
      {
         print FAILPCA "$id\n";
      }
      elsif ($group eq "Case" && $id =~ /^(.+_.+_)(.+)$/)
      {
         $id = "$1hpgen$2";
         print FAILPCA "$id\n";
      }
   }
}

close PCA;
close FAILPCA;

exit(0);

