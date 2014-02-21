#!/usr/bin/perl -w

#****************************************************************************************#
#* gff_chromosome_switcher.pl															*#
#* Switches the order of two chromosomes in a gff file e.g. Chr I and Chr II in         *# 
#* V. Cholerae                   														*#
#****************************************************************************************#

# Required modules
use strict;

# Other modules used
use Getopt::Long;

# Globals
my $helpmessage = <<END;
Usage: gff_chromosome_switcher.pl -i [FILE] -o [FILE]
Switches the order of two chromosomes in a gff file e.g. Chr II and Chr I in the 
V. Cholerae N16961 reference
	
	-i, --input    The input gff file
	-o, --output   The output gff file
	-h, --help     Displays this help message
END

#****************************************************************************************#
#* 																						*#
#* Subroutines																            *# 
#*                   			 														*#
#****************************************************************************************#
sub GrepLineNumbers ($$)
{
	my ($grep_query, $inputfile) = @_;
	
	# Returns an array of line numbers where the query is found in the input file
	my @line_numbers;
	
	my $grep_command = "grep -n \"$grep_query\" $inputfile";
	my $query_return = `$grep_command`;

	my @grep_lines = split(/\n/, $query_return);
	
	foreach my $grep_line (@grep_lines) 
	{
		if ($grep_line =~ /^(\d+):/) 
		{
			push(@line_numbers, $1);
		}
		else
		{
			die("File in wrong format!\n");
		}
	}
	
	return (\@line_numbers);
}

#****************************************************************************************#
#* 																						*#
#* Main																		            *# 
#*                   			 														*#
#****************************************************************************************#

my ($inputfile, $outputfile, $helprequested);
my @headerboundaries;

# Get options and error check first
GetOptions ("input|i=s"  => \$inputfile,
			"output|o=s" => \$outputfile,
			"help|h"       => \$helprequested
		   ) or die($helpmessage);
		   
if ($helprequested)
{
	print $helpmessage;
}
elsif (!-e $inputfile)
{
	print "Input file $inputfile doesn't exist!\n";
	print $helpmessage;
}
else
{
	# First grep for relevant file numbers
	my $annotation_header_lines = GrepLineNumbers("chromosome=", $inputfile);
	my $fasta_header_lines = GrepLineNumbers(">", $inputfile);
	
	# Count lines in the file
	my $line_count = `wc -l $inputfile`;
	my @wc_return = split(/\s+/, $line_count);
	my $file_end = $wc_return[0];
	
	# Now a big hacky sed command that does all the swapping, then cats it all together in
	# the right order
	# NB Probably really slow as it will have to read through the entire file multiple
	# times
	# Also no check for stderr or success - so for now better hope it works...
	my $fasta_boundary = @$fasta_header_lines[0];
	my $fasta_title = $fasta_boundary - 1;
	my $second_annotation_end = $fasta_boundary - 2;
	
	my $first_annotation_start = @$annotation_header_lines[0];
	my $second_annotation_start = @$annotation_header_lines[1];
	my $first_annotation_end = $second_annotation_start - 1;
	
	my $second_fasta_start = @$fasta_header_lines[1];
	my $first_fasta_end = $second_fasta_start - 1;
	
	my $sed_command = "( sed -n '1p' $inputfile ; sed -n '3p' $inputfile ; sed -n '2p' $inputfile ; sed -n '$second_annotation_start," 
	. $second_annotation_end . "p' $inputfile ; sed -n '$first_annotation_start," . $first_annotation_end
	. "p' $inputfile ; sed -n '" . $fasta_title . "p' $inputfile ; sed -n '$second_fasta_start,"
	. $file_end . "p' $inputfile ; sed -n '$fasta_boundary," . $first_fasta_end . "p' $inputfile ) | cat > $outputfile";
	print "About to run '$sed_command' \n\n";
	system($sed_command);
}


exit(0);
		