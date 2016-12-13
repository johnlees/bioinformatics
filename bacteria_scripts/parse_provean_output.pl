#!/usr/bin/perl -w

use strict;
use warnings;

my $file_list = $ARGV[0];
open(FILES, $file_list) || die("Could not open $file_list\n");

while (my $line_in = <FILES>)
{
   chomp $line_in;

   open(PROVEAN, $line_in) || die("Could not open output file $line_in\n");
   my $info_file = "$line_in.info"; # To diagnose error
   my @pathogenic_vars;
   while (my $provean_line = <PROVEAN>)
   {
      chomp $provean_line;

      if ($provean_line =~ /^# Query sequence file:\s+(FM211187.+)\.fa$/)
      {
         $info_file = "$1.info";
      }
      elsif ($provean_line =~ /^# VARIATION/)
      {
         while (my $var_line = <PROVEAN>)
         {
            chomp $var_line;
            if ($var_line eq "")
            {
               last;
            }
            else
            {
               my ($var, $score) = split("\t", $var_line);
               if ($score < -2.5)
               {
                  push (@pathogenic_vars, 1);
               }
               else
               {
                  push (@pathogenic_vars, 0);
               }
            }
         }
         last;
      }
   }
   open(INFO, $info_file) || die("Could not open info file $info_file\n");
   for (my $i = 0; $i < scalar(@pathogenic_vars); $i++)
   {
      my $info_line = <INFO>;
      chomp $info_line;

      if ($pathogenic_vars[$i] == 1)
      {
         print "$info_line\n";
      }
   }
}

exit(0);

