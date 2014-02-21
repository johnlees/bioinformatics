#!/usr/bin/perl -w

#****************************************************************************************#
#* strat_sample.pl                                                                      *#
#* Creates a stratified random samples from 2-column data                               *#
#****************************************************************************************#

use strict;

use Getopt::Long;
use List::Util qw(shuffle);

my $help_message = <<END;
Usage: strat_sample.pl --input <input_file> --num <number of samples>

    -h, --help          Displays this help message
END

#* Read command line options
my ($input_file, $sample_size, $help);
GetOptions ("input|i=s"  => \$input_file,
            "num|n=o"  => \$sample_size,
            "help|h" => \$help
		   ) or die($help_message);
		   
open (INPUT, "$input_file") || die("Couldn't open $input_file");

my %data;
my %size;
my @samples;
my $totalsize;

while (my $input_line = <INPUT>)
{
	chomp($input_line);
	my ($label, $value) = split(/,/, $input_line);
	
	push(@{$data{$label}}, $value);
	$size{$label}++;
	$totalsize++;
}

foreach my $label (keys %data)
{
	my $samples_required = sprintf("%.0f", $size{$label}/$totalsize * $sample_size);
	shuffle(@{$data{$label}});
	
	print "$label = $samples_required\n";
	
	push(@samples, @{$data{$label}}[0..$samples_required-1]);
	
}

foreach my $sample (@samples)
{
	print("$sample\n");
}
    
exit(0);