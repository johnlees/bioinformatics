#!/usr/bin/env perl -w

# Converts a bacterial GFF and FASTA file to cache for use with VEP
# An interface to the ensembl gtf2vep script, which requires some extra info in
# the GFF description fields to work
#
use strict;
use warnings;

my $gtf2vep_location = "~/installations/ensembl-tools-release-78/scripts/variant_effect_predictor/gtf2vep.pl";
my $tmp_gff = "tmp.gff";
my $species = "Streptoccocus pneumoniae R6";
my $database_version = "25";

my $usage = "perl gff2vepcache.pl genes.gff sequence.fa";

# Add to necessary attribute fields
sub add_trans_exon($$)
{
   my ($attributes, $cds_nr) = @_;

   my (%new_attributes, @order);
   foreach my $attribute (@$attributes)
   {
      my ($name, $value) = split("=", $attribute);
      $new_attributes{$name} = $value;
      push(@order, $name);
   }

   if (!defined($new_attributes{transcript_id}))
   {
      if(defined($new_attributes{locus_tag}))
      {
         $new_attributes{transcript_id} = $new_attributes{locus_tag};
      }
      else
      {
         $new_attributes{transcript_id} = "cds.$cds_nr";
      }
   }

   # These are bacteria
   if(!defined($new_attributes{exon_number}))
   {
      $new_attributes{exon_number} = 1;
   }

   # Reformat as a string to return
   my $attr_string;
   foreach my $name (@order)
   {
      $attr_string .= "$name=" . $new_attributes{$name} . ";";
   }
   chop($attr_string); # Remove final semi-colon

   return $attr_string;
}

# Inputs
my $gff = $ARGV[0];
my $fasta = $ARGV[1];

if (!-e $gff || !-e $fasta)
{
   print STDERR "Could not find input files\n";
   print STDERR $usage;
}
else
{
   open(GFF, $gff) || die("Could not open gff file $gff: $!\n");
   open(TMP, ">$tmp_gff") || die("Could not write to tmp gff file $tmp_gff: $!\n");

   # Parse gff fields. Add necessary description
   my $cds_nr = 0;
   while (my $gff_line = <GFF>)
   {
      if ($gff_line =~ /^##/)
      {
         print TMP $gff_line;
      }
      else
      {
         my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split("\t", $gff_line);

         chomp $attribute;
         my @attributes = split(";", $attribute);

         if ($feature eq "CDS")
         {
            $cds_nr++;
            my $new_attributes = add_trans_exon(\@attributes, $cds_nr);
            print TMP join("\t", $seqname, $source, $feature, $start, $end, $score, $strand, $frame, $new_attributes) . "\n";
         }
         else
         {
            print TMP $gff_line;
         }
      }
   }

   close GFF;
   close TMP;

   # Run script
   print STDERR "Running gtf2vep\n";
   system("perl $gtf2vep_location -i $tmp_gff -f $fasta -s $species -d $database_version");

   unlink $tmp_gff;
}

exit(0);

