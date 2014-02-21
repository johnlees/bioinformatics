#!/usr/bin/perl -w

#****************************************************************************************#
#* pbwt_ref_explore.pl                                                                  *#
#* Runs a pbwt imputation on a range reference panels                                   *#
#****************************************************************************************#

# bsub command:
# bsub -J "pbwt[1-6]" -o pbwt.%J.%I.o -e pbwt.%J.%I.e -R "select[mem>250] rusage[mem=250]" -M250 -u jl11@sanger.ac.uk ./pbwt_ref_explore.pl

use strict;

use Math::Combinatorics;

my @refs = ("AMD.pbwt", "FINNS.pbwt", "GoNL.pbwt", "SARDINIA.pbwt", "1000GP.pbwt", "UK10K.pbwt");

my $array_num = $ENV{'LSB_JOBINDEX'};
my $combinat = Math::Combinatorics->new(count => $array_num,
                                          data => [@refs],
                                         );
                                         
open(PLOT, ">impute_performance.$array_num.plot") or die("Couldn't open impute_performance.$array_num.plot for writing\n");

while(my @combo = $combinat->next_combination)
{
	# Create reference
	my $merge_command;
	if ($array_num == 1)
	{
		if ($combo[0] =~ /^(.+)\.pbwt$/)
		{
			$merge_command = "pbwt -readAll $1 -writeAll noORCADES.$array_num";
		}
	}
	else
	{
		$merge_command = "pbwt -merge " . join(" ", @combo) . " -writeAll noORCADES.$array_num";
	}
	print("merge_command=$merge_command\n");	
	system($merge_command);

	# Run imputation, and extract performance at 1%
	my $impute_command = "pbwt -readAll ORCADES.frame -referenceImpute noORCADES.$array_num -removeSites ORCADES.frame.sites -genotypeCompare ORCADES.noframe | cut -f 14 | head -6 | tail -1";
	print("impute_command=$impute_command\n");
	my $impute_performance = `$impute_command`;
		
	chomp($impute_performance);		
	print PLOT join(",", @combo) . " $impute_performance\n";
}

# Clean up temp files
unlink("noORCADES.$array_num.pbwt");
unlink("noORCADES.$array_num.sites");
unlink("noORCADES.$array_num.samples");

close PLOT;

exit(0);