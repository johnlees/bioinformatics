#!/usr/bin/perl -w

#****************************************************************************************#
#* impute2_genotype_compare.pl                                                          *#
#* Calculates r^2 between imputed and actual genotype (in .gen format)                  *#
#****************************************************************************************#

use strict;

use Getopt::Long;
use Statistics::LineFit;

my ($input_gen, ,$ref_gen, $input_frequencies, $output_file);
GetOptions ("input|i=s"  => \$input_gen,
            "ref|r=s" => \$ref_gen,
            "output|o=s" => \$output_file,
            "frequencies|f=s" => \$input_frequencies,
		   ) or die("Couldn't read command line input");

sub hash_freq($)
{
   my ($freq_file) = @_;

   my %freq_hash;
   my @bins = (0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0);

   open (FREQ, $freq_file) || die ("Couldn't open $freq_file");

   my $header = <FREQ>;

   while (my $freq_line = <FREQ>)
   {
      chomp($freq_line);

      # Parse each line
      my ($chr, $pos, $n_allele, $n_chr, $a0, $a1) = split(/\s+/, $freq_line);
      #print("$chr, $pos, $n_allele, $n_chr, $a0, $a1\n");
      my ($a1_gt, $maf) = split(/:/, $a1);
      $maf *= 100;
      #print("$a1_gt, $maf\n");

      # Bin values
      for (my $i=0; $i<scalar(@bins); $i++)
      {
         if ($maf <= $bins[$i])
         {
            $freq_hash{$pos} = $bins[$i];
            #print("$pos, $bins[$i]\n");
            last;
         }
      }
   }

   close FREQ;

   return \%freq_hash
}

sub parse_gen_gt ($)
{
   my ($gt_ref) = @_;

   my $gt_length = scalar(@$gt_ref);

   # Get minor allele count for each sample
   my @dosage;
   for (my $i=0; $i<$gt_length; $i+= 3)
   {
      # Allele dosage as in impute2 paper doi: 10.1534/g3.111.001198
      $dosage[$i/3] = 0;
      for (my $x=0; $x<=2; $x++)
      {
         $dosage[$i/3] += $$gt_ref[$i+$x] * $x;
      }

   }

   return \@dosage;
}

# Build freq hash for binning
print("Building site-maf hash table\n");
my $freq_hash_ref = hash_freq($input_frequencies);

# Set up hashes for regression calculations
my (%x, %y);

print("Reading input gen files\n");
# Read from the gen files
open(IN, $input_gen) || die("Couldn't open $input_gen");
open(REF, $ref_gen) || die("Couldn't open $ref_gen");

# Set up ref gen
my ($ref_chr, $ref_site_ref, $ref_pos, $ref_a0, $ref_a1, @ref_gt);
my $ref_line = <REF>;
chomp($ref_line);

while (my $gen_line = <IN>)
{
   chomp($gen_line);

   my ($chr, $site_ref, $pos, $a0, $a1, @gt) = split(/\s+/, $gen_line);
   my $imputed_dosage = parse_gen_gt(\@gt);

   ($ref_chr, $ref_site_ref, $ref_pos, $ref_a0, $ref_a1, @ref_gt) = split(/\s+/, $ref_line);
   my $ref_dosage = parse_gen_gt(\@ref_gt);

   # Only take intersection of sites
   if ($pos eq $ref_pos)
   {
      # Read forward in REF
      $ref_line = <REF>;
      chomp($ref_line);

      # Find correct bin, and put in frequencies
      my $maf_bin = $$freq_hash_ref{$pos};

      push(@{$x{$maf_bin}}, @$ref_dosage);
      push(@{$y{$maf_bin}}, @$imputed_dosage);
   }
}

close IN;
close REF;

print("Performing regression on each bin\n\n");
if ($output_file)
{
   open (OUT, ">$output_file");
   print OUT "sites maf r2\n"
}
print("sites maf r2\n");

foreach my $maf (sort{$a<=>$b} keys %x)
{
   my $linearFit = Statistics::LineFit->new();
   $linearFit->setData (\@{$x{$maf}}, \@{$y{$maf}}) or die "Invalid data";

   my $r2 = $linearFit->rSquared;
   print(scalar(@{$x{$maf}}) . " $maf $r2\n");
   if ($output_file)
   {
      print OUT " $maf $r2\n"
   }
}

close OUT;

exit(0);
