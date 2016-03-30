#!/usr/bin/perl -w

use strict;
use warnings;

my $max_log = 709;

sub safe_log($)
{
   my ($val) = @_;

   if ($val == 0)
   {
      $val = $max_log;
   }
   else
   {
      $val = -log($val);
   }

}

sub read_srst2($)
{
   my ($file) = @_;

   open(SRST2, $file) || die("Couldn't open $file: $!\n");
   my $header = <SRST2>;

   my %info;
   while (my $line_in = <SRST2>)
   {
      chomp $line_in;
      my ($allele, $score, $depth1, $depth2, $depth3, $coverage, $size, $mismatches, $indels, $truncated, $depth, $maf, $lr, $lm, $ld, $lp) = split("\t", $line_in);

      # Get allele
      $allele =~ m/^0__pspC__(pspC-\d+\.\d+)__\d+$/;
      $allele = $1;

      # Some results appear twice
      if (defined($info{$allele}{coverage}))
      {
         next;
      }

      # Transform variables
      $coverage = safe_log((100 - $coverage)/100);

      # Save variables
      $info{$allele}{coverage} = $coverage;
      $info{$allele}{mismatches} = $mismatches;
      $info{$allele}{indels} = $indels;
      $info{$allele}{truncated} = $truncated;
   }

   return(\%info);
}

sub read_blastp($)
{
   my ($file) = @_;

   open(BLASTP, $file) || die("Couldn't open $file: $!\n");

   my %info;
   while (my $line_in = <BLASTP>)
   {
      chomp $line_in;
      my ($gene, $allele, $p_id, $length, $mismatches, $gaps, $qstart, $qend, $tstart, $tend, $evalue, $bitscore) = split("\t", $line_in);

      # Transform variables
      $p_id = safe_log((100 - $p_id)/100);
      $evalue = safe_log($evalue);

      $info{$allele}{p_id} = $p_id;
      $info{$allele}{length} = $length;
      $info{$allele}{mismatches} = $mismatches;
      $info{$allele}{gaps} = $gaps;
      $info{$allele}{evalue} = $evalue;
      $info{$allele}{bitscore} = $bitscore;
   }

   return(\%info);
}

sub print_results($$$)
{
   my ($allele,$hash, $fields) = @_;

   if (defined($$hash{$allele}{mismatches}))
   {
      foreach my $field (sort keys %{ $$hash{$allele}})
      {
         print "\t" . $$hash{$allele}{$field};
      }
   }
   else
   {
      print "\t" . join("\t", ("NA") x $fields);
   }
}

# Main

my @cbpA_alleles = ("pspC-1.1","pspC-2.2","pspC-2.3","pspC-2.4","pspC-2.5","pspC-3.10","pspC-3.11","pspC-3.12","pspC-3.13","pspC-3.1","pspC-3.5","pspC-3.6","pspC-3.7","pspC-3.8","pspC-3.9","pspC-4.1","pspC-5.1","pspC-5.2","pspC-6.10","pspC-6.11","pspC-6.12","pspC-6.13","pspC-6.14","pspC-6.1","pspC-6.2","pspC-6.3","pspC-6.4","pspC-6.5","pspC-6.6","pspC-6.7","pspC-6.8","pspC-6.9");
my @pspC_alleles = ("pspC-10.1","pspC-11.1","pspC-11.2","pspC-11.3", "pspC-7.1","pspC-7.2","pspC-7.3","pspC-7.4","pspC-8.1","pspC-8.2","pspC-8.3","pspC-8.4","pspC-9.1","pspC-9.2","pspC-9.3","pspC-9.4");

my @all_alleles = @cbpA_alleles;
push(@all_alleles, @pspC_alleles);

# Print the header
print "sample.name";

my @fields = ("srst2.coverage", "srst2.mismatches", "srst2.indels", "srst2.truncated", "velvet.p.id", "velvet.length", "velvet.mismatches", "velvet.gaps", "velvet.evalue", "velvet.bitscore", "spades.p.id", "spades.length", "spades.mismatches", "spades.gaps", "spades.evalue", "spades.bitscore");
foreach my $allele (@all_alleles)
{
   foreach my $field (@fields)
   {
      print "\t$allele.$field";
   }
}

print "\n";

# Print the matrix
my @sample_list = glob("blastp_typing/*.velvet.blastp.top");

foreach my $sample (@sample_list)
{
   my ($lane, $sample_name);
   if ($sample =~ /^blastp_typing\/(\d+_\d+)_(\d+)\.velvet\.blastp\.top$/)
   {
      $lane = "$1#$2";
      $sample_name = "$1_$2";
   }
   elsif ($sample =~ /^blastp_typing\/(pspC-\d+\.\d+)\.velvet\.blastp\.top$/)
   {
      $lane = $1;
      $sample_name = $1;
   }

   my $srst2_results = read_srst2("srst2_typing/$sample_name\__$lane.pspC_typing.dna.scores");
   my $blastp_velvet_results = read_blastp("blastp_typing/$sample_name.velvet.blastp.top");
   my $blastp_spades_results = read_blastp("blastp_typing/$sample_name.spades.blastp.top");

   # print here
   print "$sample_name";

   foreach my $allele (@all_alleles)
   {
      # print srst2 results
      print_results($allele,$srst2_results, 4);
      print_results($allele,$blastp_velvet_results, 6);
      print_results($allele,$blastp_spades_results, 6);
   }

   print "\n";
}

exit(0);

