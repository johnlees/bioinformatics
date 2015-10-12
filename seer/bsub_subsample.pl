#!/usr/bin/perl -w

use strict;
use warnings;

# bsub -o logs/power.%J.%I.o -e logs/power.%J.%I.e -J "power[1-13]" power_subsample/bsub_subsample.pl
# number of jobs = 1+log(max_OR)/(log(start_or) + log(or_step))

my $start_or = 1.2;
my $or_step = 1.2;

my $job_nr = $ENV{'LSB_JOBINDEX'};
my $or = $start_or*($or_step**($job_nr-1));

system("power_subsample/subsample_seer G346_present.txt ../simulation_kmers/scaled_structure G346_kmers_in.txt.gz $or > G346_power.$job_nr.txt");

exit(0);

