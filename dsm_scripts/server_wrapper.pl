#!/usr/bin/perl -w

use strict;
use warnings;

my $dsm_location = "/nfs/users/nfs_j/jl11/installations/metaminer_dev";

my $process = $ENV{'LSB_JOBINDEX'};

my $sed_command = "sed '$process" . "q;d' $ARGV[3]";
my $file = `$sed_command`;
chomp $file;

my $process_port = $ARGV[2] + $process;

my $file_name = $file;
if ($file =~ /^(.+)\.(fasta|fa)$/)
{
   $file_name = $1;
}

my $hostname = `hostname`;
chomp $hostname;

open(SERVER, ">dsm-tmp/meta-server_config_$process.txt") || die("Couldn't open config file $process");
print SERVER join("\t", $file_name, $hostname, $process_port) . "\n";
close SERVER;

system("$dsm_location/meta-server --fmin $ARGV[0] --maxdepth $ARGV[1] -v -P $process_port $file");

exit(0);

