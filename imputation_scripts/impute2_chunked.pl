#!/usr/bin/perl -w

use strict;

# Set up environment
system("source ~/.bashrc");

# Get index number of this process
my $array_num = $ENV{'LSB_JOBINDEX'};

my $start = ($array_num - 1)*5 . "e6";
my $end = ($array_num)*5 . "e6";

system("impute2 -m genetic_map_chr20_combined_b37.txt -h REF.haps -l REF.legend -known_haps_g TEST.frame.gen -int $start $end -Ne 20000 -o phased.impute2.$array_num -phase -allow_large_regions");

exit(0);