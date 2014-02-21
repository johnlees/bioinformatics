#!/usr/bin/perl -w

#****************************************************************************************#
#* run_abacas.pl                                                                        *#
#* Runs abacas on all assemblies in the pwd                                             *#
#****************************************************************************************#

# Required modules
use strict;

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

# Get list of files to process
my @file_list = glob "*_automated_velvet.fa";

# Submit each file to abacas using bsub
# NB reference file has been hard coded
foreach my $assembly_file (@file_list)
{
	if ($assembly_file =~ /^(\d+_\d+#\d+)_automated_velvet\.fa$/)
	{	
		my $command = "bsub -o abacas/logs/$1.abacas.log abacas.pl -r ../2011_data/N16961.fa -q $assembly_file -p nucmer -b -d -a -o abacas.$1";
		system($command);
	}
}

exit(0);

