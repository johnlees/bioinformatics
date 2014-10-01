#!/usr/bin/perl -w

use strict;
use warnings;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 );

use Getopt::Long;

use compare_variants;

my $usage_message = "./compare_variants.pl --vcf1 1.vcf.gz --vcf2 2.vcf.gz --ref1 ref1.fa --ref2 ref2.fa\n";

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($vcf1, $ref1, $vcf2, $ref2, $help);
GetOptions( "vcf1=s" => \$vcf1,
            "vcf2=s" => \$vcf2,
            "ref1=s" => \$ref1,
            "ref2=s" => \$ref2,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help))
{
   print "$usage_message";
}
else
{
   my @vcfs = ($vcf1, $vcf2);
   my @refs = ($ref1, $ref2);

   #TODO: it might actually be better to use BLAT. It probably won't make much
   # difference but this is kind of what it's designed for from what I can tell,
   # and the input creation may be easier. It will run faster, but on 300bp
   # this doesn't make a difference
   compare_variants::create_blastn_input(\@vcfs, \@refs, "blast_windows");
   my $blast_scores = compare_variants::blastn_pairwise("blast_windows.1.fa", "blast_windows.2.fa");

   foreach my $q_id (sort keys %$blast_scores)
   {
      foreach my $s_id (sort keys %{$$blast_scores{$q_id}})
      {
         print "$q_id\t$s_id\t$$blast_scores{$q_id}{$s_id}\n";
      }
   }

   unlink "blast_windows.1.fa", "blast_windows.2.fa";
}

exit(0);

