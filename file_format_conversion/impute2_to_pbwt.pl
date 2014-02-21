#!/usr/bin/perl -w

#****************************************************************************************#
#* impute2_to_pbwt.pl                                                                   *#
#* Runs a pbwt imputation on a range reference panels                                   *#
#****************************************************************************************#

use strict;

use Getopt::Long;
use List::Util qw( min max );
use Data::Dumper;

my $help_message = "Usage: impute2_to_pbwt.pl -i <input.gen> -o <output prefix> -c <chromosome>";

my ($input_file, $output_prefix, $chromosome);
GetOptions ("input|i=s"  => \$input_file,
            "output|o=s" => \$output_prefix,
            "chr|c=s" => \$chromosome
		   ) or die("Couldn't read command line input");

if (!defined($input_file) || !defined($output_prefix) || !defined($chromosome))
{
   die($help_message);
}

open(GEN, $input_file) || die ("Couldn't open $input_file");
open(NEWGEN, ">tmp.gen") || die ("Could write to tmp.gen");

while (my $gen_line = <GEN>)
{
	my @samples = split(" ", $gen_line);

	my @new_hap;
	push(@new_hap, @samples[0..4]);

	for (my $i = 5; $i < scalar(@samples); $i += 3)
	{
		my @triplet = @samples[$i..$i+2];

		my $max = max @triplet;

		my $max_found = 0;
		for (my $j = 0; $j < 3; $j++)
		{
			if (!$max_found && $triplet[$j] == $max)
			{
				push(@new_hap, "1");
				$max_found = 1;
			}
			else
			{
				push(@new_hap, "0");
			}
		}
	}

	my $new_gen_line = join(" ", @new_hap) . "\n";
	print NEWGEN $new_gen_line;
}

system("pbwt -readGen tmp.gen $chromosome -writeAll $output_prefix");
unlink("tmp.gen");

exit(0);
