#!/usr/bin/perl -w

use strict;

# Assuming everything in the map has been flipped so everything refers to the
# top strand

use Text::CSV;

my $csv = Text::CSV->new ( { binary => 1, sep_char => ",", auto_diag => 1, } )  # should set binary attribute.
                 or die "Cannot use CSV: ".Text::CSV->error_diag ();

my $manifest_file = "humanomniexpress-12v1-1_b.csv";
my $flip_file = "flip_list.txt";

open my $manifest_in, "<:encoding(utf8)", "$manifest_file" or die "$manifest_file: $!";

while (my $row = $csv->getline($manifest_in))
{
   if ($row->[0] eq "[Assay]")
   {
      print "Header row found\n";
      # eat header row
      $csv->column_names($csv->getline($manifest_in));
      last;
   }
}

open (FLIP, ">$flip_file") || die ("Could not open $flip_file: $!");

while (my $row = $csv->getline_hr($manifest_in))
{
   if (($row->{IlmnStrand} eq "TOP" && $row->{RefStrand} eq "-") || ($row->{IlmnStrand} eq "BOT" && $row->{RefStrand} eq "+"))
   {
      print FLIP $row->{Name} . "\n";
   }
}

close FLIP;
close $manifest_in;

exit(0);
