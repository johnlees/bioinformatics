#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::SeqIO;

my $roary_directory = $ARGV[0];
my $number_samples = $ARGV[1];

sub mat_val($$$)
{
   my ($i, $j, $array) = @_;

   my $val = 0;
   if ($j!=$i || defined($$array[$i][$j]))
   {
      if ($j<$i)
      {
         $val = $$array[$j][$i];
      }
      else
      {
         $val = $$array[$i][$j];
      }
   }
   return $val;
}

if(!defined($roary_directory) || !defined($number_samples) || !-d $roary_directory)
{
   print "Usage: ./roary_bigsdb_distances.pl <roary_alignment_directory> <number_of_samples>\n";
}
else
{
   my @alignment_files = glob("$roary_directory/*.aln");

   my @dist_mat;
   foreach my $alignment_file (@alignment_files)
   {
      my $alignment = Bio::SeqIO->new(-file => $alignment_file, -format => 'fasta')
      || die("Failed to open $alignment_file: $!\n");

      my @sequences;
      while (my $base_sequence = $alignment->next_seq())
      {
         $base_sequence->display_id() =~ m/^SE0+(\d+)_\d+$/;
         my $base_sample_nr = $1;

         # Takes only the first sequence if there are multiple copies of a gene
         if(!defined($sequences[$base_sample_nr]))
         {
            $sequences[$base_sample_nr] = $base_sequence->seq();
         }
      }

      for (my $i = 1; $i <= $number_samples; $i++)
      {
         for (my $j = $i+1; $j <= $number_samples; $j++)
         {
            if (!defined($sequences[$i]) || !defined($sequences[$j]))
            {
               print STDERR "Not found: pair $i, $j in $alignment_file\n";
               $dist_mat[$i-1][$j-1]++;
            }
            elsif ($sequences[$i] ne $sequences[$j])
            {
               $dist_mat[$i-1][$j-1]++;
            }
         }
      }

      $alignment->close();
   }

   for (my $i = 0; $i< $number_samples; $i++)
   {
      print mat_val($i, 0, \@dist_mat);

      for (my $j = 1; $j < $number_samples; $j++)
      {
         print "," . mat_val($i, $j, \@dist_mat);
      }
      print "\n";
   }
}

exit(0);

