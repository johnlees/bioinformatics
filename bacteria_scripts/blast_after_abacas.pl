#!/usr/bin/perl -w

#****************************************************************************************#
#* blast_after_abacas.pl                                                                *#
#* Blasts all abacas orderings against their reference                                  *#
#* Everything is hard coded...                                                          *#
#****************************************************************************************#

# Required modules
use strict;

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

# Get list of files to process
my @file_list = glob "*.fasta";

# Submit each file to abacas using bsub
# NB reference file has been hard coded
foreach my $assembly_file (@file_list)
{	
	my $command = "~sh16/scripts/better_blast.py -q ../N16961.fa -s $assembly_file";
	system($command);
}

exit(0);
