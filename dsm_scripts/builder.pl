#!/usr/bin/perl -w

use strict;
use warnings;

my $job_num = $ENV{'LSB_JOBINDEX'};
if (!defined($job_num))
{
   $job_num = $ARGV[1];
}

my $file_list = $ARGV[0];
my $builder_path = "~/installations/metaminer_dev/builder";

my $sed_command = "sed '$job_num" . "q;d' $file_list";
my $file = `$sed_command`;
chomp $file;

my $command = "$builder_path -v $file";
system($command);

exit(0);

