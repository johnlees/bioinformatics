#!/usr/bin/perl -w

use strict;
use warnings;

use Text::CSV;

my $snp_ref_file = "HumanOmniExpress-12v1_J.bpm.csv";
my $intensity_summary_file = "/lustre/scratch113/teams/barrett/meningitis/BPROOF_intensity_report.txt";
my $header_file_name = "BPROOF_header_row";
my $out_suffix = "BPROOF.int";

my %sex_map = ("23" => "X",
               "24" => "Y",
               "25" => "XY");

# Get index number of this process
#my $array_num = $ENV{'LSB_JOBINDEX'};
#my $chr_name;

#if ($array_num <= 23)
#{
#   $chr_name = $array_num-1;
#}
#else
#{
#   $chr_name = $sex_map{$array_num};
#}
#
#my $out_file = "$chr_name.$out_suffix";

# Look up SNPs: ref/alt alleles by SNP ID
my %snps;
open (SNPS, $snp_ref_file) || die("Could not open $snp_ref_file");

while (my $snp_line = <SNPS>)
{
   chomp $snp_line;
   my ($id, $name, $chr, $pos, $gentrain, $allele, $ill_str, $cus_str, $norm_id) = split(/,/, $snp_line);

   if ($allele =~ /^\[(.+)\/(.+)\]$/)
   {
      $snps{"allele"}{$name} = $1 . $2;
      $snps{"pos"}{$name} = $pos;
   }
}

close SNPS;
print "SNPS parsed\n";

# TODO IT WOULD BE WISE to first split the intensity report into chromosomes
my $csv = Text::CSV->new ( { binary => 1, sep_char => "\t", auto_diag => 1, eol => $/ } )  # should set binary attribute.
                 or die "Cannot use CSV: ".Text::CSV->error_diag ();

open my $fh_in, "<:encoding(utf8)", "$intensity_summary_file" or die "$intensity_summary_file: $!";

my %fh_out;
for (my $i = 0; $i <= 25; $i++)
{
   my $chr_name;

   if ($i < 23)
   {
      $chr_name = $i;
   }
   else
   {
      $chr_name = $sex_map{$i};
   }

   my $out_file = "$chr_name.$out_suffix";
   open $fh_out{$chr_name}, ">:encoding(utf8)", "$out_file" or die "$out_file: $!";
}

# Skip header
my $row = $csv->getline( $fh_in );

while ($row = $csv->getline( $fh_in ))
{
   my $chr_name = $row->[1];

   if (defined($snps{"pos"}{$row->[0]}))
   {
      # Splice removes fields 1 and 2 (zero indexed) and adds in position
      # and alleles instead
      my @new_fields = ($snps{"pos"}{$row->[0]},$snps{"allele"}{$row->[0]});
      splice @$row, 1, 2, @new_fields;
   }
   else
   {
      print "Error: " . $row->[0] . " not found\n"
   }

   $csv->print($fh_out{$chr_name}, $row);
}

close $fh_in;
foreach my $chr (keys %fh_out)
{
   close $fh_out{$chr};
}

print "File spliced\n";

#system("cat $header_file_name $out_file > tmpout.$chr_name.txt");
#rename "tmpout.$chr_name.txt", $out_file or die $!;

exit(0);
