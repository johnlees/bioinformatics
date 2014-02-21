#!/usr/bin/perl -w

use strict;

# CPAN modules
use Text::CSV;

#
# Main
#

my $csv = Text::CSV->new ( { binary => 1 } );
open my $csv_in "<:encoding(utf8)", "Tree_annotation.csv" or die "Tree_annotation.csv: $!";

open my $csv_out "<:encoding(utf8)", "Tree_annotation_out.csv" or die "Tree_annotation_out.csv: $!";

while ( my $row = $csv_obj->getline( $annotation_csv ) ) 
{
	if ($row->[0] =~ m/(.*)_(.*)_(.*)/)
	{
		my $ref_name = "$1_$2#$3";
		my $new_name = $row->[1] . "_" . $row->[2] . $row ->[3];
		
		my $print_out = $csv-> print ($csv_out, $ref_name, $new_name);
	}
	
}

close $csv_in;
close $csv_out;