#!/usr/bin/perl -w

# Gives an estimate of coverage over a certain region
# See simple_serotype.pl

use strict;
use warnings;

my $read_count_lim = 10;

# Returns an 8 character string of random alphanumeric characters
sub random_string()
{
   my $string = join'', map +(0..9,'a'..'z','A'..'Z')[rand(10+26*2)], 1..8;

   return $string;
}

my $bam_in = $ARGV[0];
my $region = $ARGV[1]; # CP006046:1679226-1719774
my $ref = $ARGV[2];  # /lustre/scratch108/pathogen/pathpipe/refs/Listeria/monocytogenes_J1816/Listeria_monocytogenes_J1816_v1.fa

$region =~ m/^(.+):(\d+)-(\d+)$/;
my $region_start = $2;
my $region_end= $3;

my $tmp_file = random_string() . ".mpileup.tmp";

system("/nfs/users/nfs_j/jl11/software/bin/samtools view -h $bam_in $region | /nfs/users/nfs_j/jl11/software/bin/samtools mpileup -q 20 -f $ref - > $tmp_file");

open(MPILEUP, $tmp_file) || die("Could not open mpileup result $tmp_file\n");

my $total = 0;
my $reads_mapped = 0;
while (my $pileup_line = <MPILEUP>)
{
   chomp $pileup_line;
   my ($chr, $pos, $base, $reads, @junk) = split("\t", $pileup_line);

   if ($pos >= $region_start && $pos <= $region_end)
   {
      $total++;

      if ($reads >= $read_count_lim)
      {
         $reads_mapped++;
      }
   }
}

close MPILEUP;

my $coverage = $reads_mapped/$total;
print $coverage . "\n";

unlink $tmp_file;

exit(0);

