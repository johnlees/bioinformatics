#!/usr/bin/perl -w

use strict;
use warnings;

sub max ($$) { $_[$_[0] < $_[1]] }

sub compressed_size($)
{
   my ($file) = @_;

   system("7z a -m0=lzma $file.7z $file 1>&2");
   my $size = -s "$file.7z";

   unlink "$file.7z";

   return $size;
}

sub ncd($$$$)
{
   my ($file1_size, $file1, $file2_size, $file2) = @_;

   my $f12 = "cat_1_2.tmp";
   my $f21 = "cat_2_1.tmp";

   system("cat $file1 $file2 > $f12");
   system("cat $file2 $file1 > $f21");

   my $size_12 = compressed_size($f12);
   my $size_21 = compressed_size($f21);
   unlink $f12, $f21;

   my $ncd = max($size_12-$file1_size, $size_21-$file2_size) / max($file1_size, $file2_size);
   return $ncd;
}

sub compressed_sizes($)
{
   my ($files) = @_;

   my @compressed_sizes;
   foreach my $file (@$files)
   {
      push(@compressed_sizes, compressed_size($file));
   }

   return(@compressed_sizes);
}

my $input_genomes = $ARGV[0];
open(INPUT, $input_genomes) || die("Could not open $input_genomes\n");

my @input_files;
while (my $input_file = <INPUT>)
{
   chomp $input_file;
   push(@input_files, $input_file);
}

close INPUT;

my @compressed_sizes = compressed_sizes(\@input_files);

my @distances;
for (my $i=0; $i<scalar(@input_files); $i++)
{
   for (my $j=$i; $j<scalar(@input_files); $j++)
   {
      if ($i == $j)
      {
         $distances[$i][$j] = 0;
      }
      else
      {
         my $ncd_pair = ncd($compressed_sizes[$i], $input_files[$i], $compressed_sizes[$j], $input_files[$j]);
         $distances[$i][$j] = $ncd_pair;
         $distances[$j][$i] = $ncd_pair;
      }
   }
   print join("\t", @{$distances[$i]}) . "\n";
}

exit(0);

