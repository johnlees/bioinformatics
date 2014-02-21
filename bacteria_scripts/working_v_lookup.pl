#!/usr/bin/perl -w

#****************************************************************************************#
#* working_v_lookup.pl                                                                  *#
#*                                                                                      *#
#****************************************************************************************#

use strict;

use Getopt::Long;

#* Gets input directory supplied by -i and reads the names of all contained files
my ($inputfile);
GetOptions ("input|i=s"  => \$inputfile
		   ) or die($!);
		   
open (CSVIN, $inputfile);

my %data;
my @lookup;

while (my $linein = <CSVIN>)
{
   # Line endings are /r/n
   chop($linein);
   chop($linein);
   my @fields = split(/,/, $linein);
   
   #$data{$fields[0]} = join(',', $fields[2], $fields[1]);
   $data{$fields[0]} = $fields[1];
   
   #print $fields[0] . " " . $data{$fields[0]} . "\n";
   push(@lookup, $fields[2]);
   # print $fields[3] . "\n";
}

close CSVIN;

foreach my $lookupval (@lookup)
{
   print $data{$lookupval} . "\n";
}

exit(0);