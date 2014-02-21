#!/usr/bin/perl -w

#****************************************************************************************#
#* pbwt_corrupt.pl                                                                      *#
#* Runs a pbwt imputation on a range of corrupted references                            *#
#****************************************************************************************#

use strict;

my @range = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1);

open(PLOT, ">corrupt_performance.plot");

for (my $i = 0; $i < scalar(@range); $i++)
{
	my $p = $range[$i];
	print PLOT "$p";
	
	for (my $j = 0; $j < scalar(@range); $j++)
	{
		my $q = $range[$j];
		
		# Create corrupted reference
		system("pbwt -readAll ../REF -corruptSamples $p $q -writeAll REF.corrupt");
		
		# Run imputation, and extract performance at 1%
		my $impute_performance = `pbwt -readAll ../TEST.frame -referenceImpute REF.corrupt -removeSites ../TEST.frame.sites -genotypeCompare ../TEST.noframe | cut -f 14 | head -6 | tail -1`;
		
		chomp($impute_performance);		
		print PLOT " $impute_performance";
	}
	print PLOT "\n";
}

# Clean up temp files
unlink("REF.corrupt.pbwt");
unlink("REF.corrupt.sites");
unlink("REF.corrupt.samples");

close PLOT;

exit(0);