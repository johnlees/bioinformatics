#!/usr/bin/perl

use strict;

my %imiss;
my %removed;

open IMISS, '<', $ARGV[0].".imiss"
        or die "Cannot open genotypes file (".$ARGV[0].".imiss): $!\n";
print "Reading PLINK .imiss file ".$ARGV[0].".imiss\n";
while(<IMISS>){
	s/^\s+//;
    my @fields = split /\s+/, $_;
    $imiss{$fields[0]}{$fields[1]} = $fields[5];
}

open GENOME, '<', $ARGV[0].".kin0"
        or die "Cannot open genotypes file (".$ARGV[0].".kin0): $!\n";
open OUT, '>', "fail-kinship-QC.txt";

my @kin_out;
for (my $i = 0; $i <= 3; $i++)
{
   open $kin_out[$i], ">", "kinship-$i.txt";
}

print "Reading KING .ibs0 file ".$ARGV[0].".ibs0\n";
while(<GENOME>){
    s/^\s+//;
    my @fields = split /\s+/, $_;
 	if($fields[7] > 0.0442){

      my $relatedness;
      if ($fields[7] < 0.0884)
      {
         $relatedness = 3;
      }
      elsif ($fields[7] < 0.177)
      {
         $relatedness = 2;
      }
      elsif ($fields[7] < 0.354)
      {
         $relatedness = 1;
      }
      else
      {
         $relatedness = 0;
      }

 		if($imiss{$fields[0]}{$fields[1]}>$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[0]}{$fields[1]}){
 				print OUT "$fields[0] $fields[1]\n";
            $kin_out[$relatedness]->print("$fields[0] $fields[2] $fields[6]\n");
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 		elsif($imiss{$fields[0]}{$fields[1]}<$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[2]}{$fields[3]}){
 				print OUT "$fields[2] $fields[3]\n";
            $kin_out[$relatedness]->print("$fields[0] $fields[2] $fields[6]\n");
 				$removed{$fields[2]}{$fields[3]} = 1;
 			}
 		}
 		else{
 			unless($removed{$fields[0]}{$fields[1]}){
 				print OUT "$fields[0] $fields[1]\n";
            $kin_out[$relatedness]->print("$fields[0] $fields[2] $fields[6]\n");
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 	}
}



