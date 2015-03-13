#!/usr/bin/env perl

# Converts a bacterial GFF and FASTA file to cache for use with VEP
# An interface to the ensembl gtf2vep script, which requires some extra info in
# the GFF description fields to work
#
use strict;
use warnings;

my $gtf2vep_location = "~/installations/ensembl-tools-release-78/scripts/variant_effect_predictor/gtf2vep.pl";
my $tmp_gff = "tmp.gtf";
my $database_version = "25";

my $usage = "perl gff2vepcache.pl species_name genes.gff sequence.fa";

# Add to necessary attribute fields
sub add_trans_exon($$$)
{
   my ($attributes, $gene_id, $trans_id) = @_;

   my (%new_attributes, @order);
   foreach my $attribute (@$attributes)
   {
      my ($name, $value) = split("=", $attribute);
      $new_attributes{$name} = $value;
      push(@order, $name);
   }

   if (!defined($new_attributes{transcript_id}))
   {
      push(@order, "transcript_id");
      $new_attributes{transcript_id} = $trans_id;
   }

   # These are bacteria
   if(!defined($new_attributes{exon_number}))
   {
      push(@order, "exon_number");
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

sub reformat_attribute($)
{
   my ($attribute_string) = @_;

   my @attributes = split(";", $attribute_string);

   my $new_string;
   foreach my $attribute (@attributes)
   {
      my ($key, $value) = split("=", $attribute);
      $new_string .= "$key \"$value\"; ";
   }

   return $new_string;
}

# Inputs
my $species = $ARGV[0];
my $gff = $ARGV[1];
my $fasta = $ARGV[2];

if (!-e $gff || !-e $fasta)
{
   print STDERR "Could not find input files\n";
   print STDERR $usage;
}
else
{
   open(GFF, $gff) || die("Could not open gff file $gff: $!\n");
   open(TMP, ">$tmp_gff") || die("Could not write to tmp gff file $tmp_gff: $!\n");

   # Persist between rows
   my ($gene_id, $trans_id, $gene_last, %copy, @genes);
   my $cds_nr = 0;

   # Parse gff fields. Add necessary description
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

         if ($feature eq "gene")
         {
            $gene_last = 1;

            # rename ID -> gene_id
            $attribute = "";

            foreach my $attr_pair (@attributes)
            {
               my ($key, $value) = split("=", $attr_pair);
               if ($key eq "ID")
               {
                  $gene_id = $value;
                  $attribute .= "gene_id=$gene_id;";
                  $trans_id = "trans_$gene_id";
               }
               elsif ($key eq "gene")
               {
                  # Same name as last gene?
                  my $repeat = 0;
                  foreach my $gene_name (@genes)
                  {
                     $copy{$gene_name}++;

                     if ($value eq $gene_name)
                     {
                        $trans_id = "$gene_name\_" . ($copy{$gene_name} - 1);
                        $repeat = 1;
                        last;
                     }
                  }

                  if (!$repeat)
                  {
                     $trans_id = $value;
                     push (@genes, $value);
                  }

                  $attribute .= "gene=$trans_id;";
               }
               else
               {
                  $attribute .= "$attr_pair;";
               }

            }

            print TMP join("\t", $seqname, $source, $feature, $start, $end, $score, $strand, $frame, reformat_attribute($attribute)) . "\n";

            # Need to print a transcript and an exon
            print TMP join("\t", $seqname, $source, "transcript", $start, $end, $score, $strand, $frame,
               reformat_attribute("gene_id=$gene_id;transcript_id=$trans_id")) . "\n";
            print TMP join("\t", $seqname, "protein_coding", "exon", $start, $end, $score, $strand, $frame,
               reformat_attribute("gene_id=$gene_id;transcript_id=$trans_id;exon_number=1")) . "\n";
         }
         elsif ($feature eq "CDS")
         {
            $cds_nr++;
            if ($gene_last)
            {
               $attribute = add_trans_exon(\@attributes, $gene_id, $trans_id);
            }
            else
            {
               $attribute = add_trans_exon(\@attributes, $cds_nr, "trans_$cds_nr");
            }

            $source = "protein_coding";

            print TMP join("\t", $seqname, $source, $feature, $start, $end, $score, $strand, $frame, reformat_attribute($attribute)) . "\n";

            $gene_last = 0;
         }
         else
         {
            $gene_last = 0;

            print TMP join("\t", $seqname, $source, $feature, $start, $end, $score, $strand, $frame, reformat_attribute($attribute)) . "\n";
         }

      }
   }

   close GFF;
   close TMP;

   # Run script
   print STDERR "Running gtf2vep\n";
   system("perl $gtf2vep_location -i $tmp_gff -f $fasta -s $species -d $database_version");

   #unlink $tmp_gff;
}

exit(0);

