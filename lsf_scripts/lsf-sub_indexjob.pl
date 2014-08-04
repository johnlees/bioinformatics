#!/usr/bin/perl -w

#****************************************************************************************#
#* jlsub_indexjob.pl                                                                    *#
#* Runs each array job as specified by the jlsub master script                          *#
#****************************************************************************************#

use strict;

# Takes in execution directory, list of files and command to run
my ($dir, $file_list, $parsed_command) = @ARGV;

# If the options are missing, assume jlsub has not been used so don't do anything
if (!defined($dir) || !defined($file_list) || !defined($parsed_command))
{
	print "This script should not be run directly, run through jlsub.pl instead\n"
}
else
{
	# Set up environment
	chdir($dir);
	system("source ~/.bashrc");

	# Get index number of this process
	my $array_num = $ENV{'LSB_JOBINDEX'};

	# Extract file name for this index number
	my $file_name_command = "sed '$array_num" . "!d' $file_list";
	my $file_name = `$file_name_command`;
	chomp($file_name);

	# Replace FILE keyword with actual file name in command
	$parsed_command =~ s/FILE/$file_name/g;

	# Run command
    system($parsed_command);
    #print "$parsed_command\n";
}

exit(0);