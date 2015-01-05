#!/usr/bin/perl -w
#
use strict;
use warnings;

my $usage = "./ivr_to_pairs.pl ivr_mapped_alleles.txt pairs.txt\n";

open(IVR, $ARGV[0]) || die($usage);

my $i = 1;
my %input;

while (my $line_in = <IVR>)
{
   chomp $line_in;

   my @fields = split("\t", $line_in);

   $input{$fields[0]} = $i;

   $i++;
}

close IVR;

open(PAIRS, $ARGV[1]) || die ($usage);

print "Sample\tCSF_index\tBlood_index\n";
while (my $line_in = <PAIRS>)
{
   chomp $line_in;

   my ($sample, $csf, $blood) = split("\t", $line_in);

   if (defined($input{$csf}) && defined($input{$blood}))
   {
      print join("\t", $sample, $input{$csf}, $input{$blood} . "\n");
   }
}

exit(0);

