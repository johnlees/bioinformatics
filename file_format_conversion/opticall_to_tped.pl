#!/usr/bin/perl -w

use strict;
use warnings;

use Text::CSV;
use Getopt::Long;

my $help_message = "Usage: opticall_to_tped -i in.calls -o out_prefix --chr <chromosome> --phen <phenotype> --info <opticall info file>\n";

sub allele_convert($$$)
{
   my ($opti_val, $a0, $a1) = @_;

   my @ped_al;

   if ($opti_val == 1)
   {
      @ped_al = ($a0, $a0);
   }
   elsif ($opti_val == 2)
   {
      @ped_al = ($a0, $a1);
   }
   elsif ($opti_val == 3)
   {
      @ped_al = ($a1, $a1);
   }
   else
   {
      @ped_al = (0, 0);
   }

   return(\@ped_al);
}

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

my ($input_file, $output_prefix, $help, $info, $phenotype, $chromosome, $info_file);
GetOptions ("input|i=s"  => \$input_file,
            "output|o=s" => \$output_prefix,
            "chr=s" => \$chromosome,
            "phen=s" => \$phenotype,
            "info=s" => \$info_file,
            "help|h"     => \$help
		   ) or die($help_message);

if (defined($help) || !defined($input_file) || !defined($output_prefix))
{
   print $help_message;
}
else
{
   my $out_tped = "$output_prefix.$chromosome.tped";
   my $out_tfam = "$output_prefix.$chromosome.tfam";

   my $csv = Text::CSV->new ( { binary => 1, sep_char => " ", auto_diag => 1, eol => $/ } )  # should set binary attribute.
                    or die "Cannot use CSV: ".Text::CSV->error_diag ();

   # First parse .info file to get genders
   open my $info_in, "<:encoding(utf8)", "$info_file" or die "$input_file: $!";

   my %genders;
   while (my $row = $csv->getline($info_in))
   {
      my $sample = $row->[0];
      my $gender = $row->[1];
      $genders{$sample} = $gender;
   }

   close $info_in;

   # Work through .call input file
   open my $calls_in, "<:encoding(utf8)", "$input_file" or die "$input_file: $!";
   open my $tped_out, ">:encoding(utf8)", "$out_tped" or die "$out_tped: $!";
   open my $tfam_out, ">:encoding(utf8)", "$out_tfam" or die "$out_tfam: $!";

   # Header line goes to tfam (this will be the same for all chromosome files)
   my $header_row = $csv->getline($calls_in);
   my $row_length = scalar(@$header_row);

   for (my $i = 4; $i < $row_length; $i++)
   {
      my $sample = $header_row->[$i];

      my $gender;
      if (defined($genders{$sample}) && $genders{$sample} == 1 || $genders{$sample} == 2)
      {
         $gender =  $genders{$sample};
      }
      else
      {
         $gender = 0;
      }

      my @tfam_line = ($sample, $sample, 0, 0, $gender, $phenotype);

      # write tfam line
      $csv->print($tfam_out, \@tfam_line);
   }

   # Remaining lines need editing slightly
   while (my $row = $csv->getline($calls_in))
   {
      my $snp_name = $row->[0];
      my $snp_pos = $row->[1];

      my ($a0, $a1);
      if ($row->[2] =~ /^(\w{1})(\w{1})$/)
      {
         $a0 = $1;
         $a1 = $2;
      }

      my @tped_line = ($chromosome, $snp_name, 0, $snp_pos);

      for (my $i = 4; $i <$row_length; $i++)
      {
         my $alleles = allele_convert($row->[$i], $a0, $a1);
         push(@tped_line, @$alleles);
      }
      # write tped line
      $csv->print($tped_out, \@tped_line);
   }

   close $calls_in;
   close $tped_out;
   close $tfam_out;
}

exit(0);

