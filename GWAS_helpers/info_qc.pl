#!/usr/bin/perl -w

# Performs QC based on the info field reported by snptest, concatenates all
# chromosome results into a format ready for input to R

use strict;

my $file_list = $ARGV[0];

open (FILES, "$file_list") || die("Could not open $file_list\n");

while (my $line = <FILES>)
{
   chomp $line;
   my ($path, $chr, $chunk, $suffix) = split("\t", $line);

   my ($info_field, $out_fields);
   if (($chr eq "X" || $chr eq "Y") && $chunk !~ /PAR/)
   {
      $info_field = 15;
      $out_fields = "2,3,4,15,64";
   }
   else
   {
      $info_field = 9;
      $out_fields = "2,3,4,9,42";
   }

   my $command = "awk '\$0!~\"^#\" && \$$info_field>0.3 {print \$0}' $line | sed '1d' | cut -d \" \" -f $out_fields >> meningitis_snptest.$chr.txt";
   system($command);
}

close FILES;

exit(0);
