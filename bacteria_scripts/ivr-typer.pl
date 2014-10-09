#!/usr/bin/perl -w
#
use strict;
use warnings;

use Getopt::Long;

use Bio::SeqIO;
use Bio::Tools::GFF;

#
# Globals
#
my $query_fasta = "hsdS_query.fa";

my $N_term_ref = "Ntermini.fasta";
my $C_term_ref = "Ctermini.fasta";

my $N_blast_out = "Ntermini.blast.out";
my $C_blast_out = "Ctermini.blast.out";
my $blast_err = "blast.err";

my %allele_map = ("Aa" => "A",
                  "Ab" => "B",
                  "Ac" => "E",
                  "Ba" => "D",
                  "Bb" => "C",
                  "Bc" => "F");

my $help_message = <<HELP;
Usage: ./ivr_typer.pl [--map|--assembly] <options>

Given S. pneumo read (via mapping mode) or annotated assembly (via assembly mode),
returns information on the likely allele type for the ivr/hsd R-M system locus

   Options
   --map                 Infer by mapping reads to reference alleles
   --assembly            Infer by annotation of gene, then a blast with reference
                         genes

   --ref_dir             Directory containing sequences of reference alleles
                         (default ./)

   <map mode>
   -f, --forward_reads   Location of forward reads (fastq(.gz))
   -r, --reverse_reads   Location of reverse reads (fastq(.gz))

   <assembly mode>
   -a, --annotation      Annotation of the -a assembly in gff v3.


   -h, --help            Displays this message

HELP

# Need to remove contigs with no features due to failings of Bio::Tools::GFF
# used with attach_seqs
# Should all be sorted by contig:position
sub pre_process_gff($$)
{
   my ($file_in, $file_out) = @_;

   my $gff_tmp = "tmp.gff";
   my $fasta_tmp = "tmp.fa";

   if (-e $gff_tmp || -e $fasta_tmp)
   {
      die("$gff_tmp or $fasta_tmp already exist. Will not overwrite\n");
   }

   my $gff_in = Bio::Tools::GFF->new(-file => $file_in,
                                     -gff_version => 3) || die ("Could not open $file_in as gff v3: $!");
   my $gff_out = Bio::Tools::GFF->new(-file => ">$gff_tmp",
                                     -gff_version => 3) || die ("Could not write to $gff_tmp as gff v3: $!");

   my (@contigs, $last_contig);
   foreach my $feature ($gff_in->next_feature())
   {
      $gff_out->write_feature($feature);

      if ($feature->seq_id() ne $last_contig)
      {
         $last_contig = $feature->seq_id();
         push(@contigs, $last_contig);
      }
   }

   my @contig_sequences = $gff_in->get_seqs();

   $gff_in->close();

   my $fasta_out = Bio::SeqIO->new( -file   => ">$fasta_tmp",
                                    -format => "fasta") || die("Could not write to $fasta_tmp: $!");

   my $next_contig = shift(@contigs);
   foreach my $contig (@contig_sequences)
   {
      if ($contig->display_id() eq $next_contig)
      {
         $next_contig = shift(@contigs);
         $fasta_out->write_seq($contig);
      }
   }

   #system("echo \"##FASTA\" | cat $gff_tmp - $fasta_tmp > $file_out");
   open(GFF, "$gff_tmp") || die("Could not open $gff_tmp\n");
   open(FASTA, "$fasta_tmp") || die("Could not open $fasta_tmp\n");
   open(OUT, ">$file_out") || die("Could not write to $file_out\n");

   while (my $gff_line = <GFF>)
   {
      print OUT $gff_line;
   }
   print OUT "##FASTA";
   while (my $fasta_line = <FASTA>)
   {
      print OUT $fasta_line;
   }

   unlink $gff_tmp, $fasta_tmp;

}

sub make_query_fasta($$$)
{
   my ($genes, $sequences, $f_out) = @_;

   my $sequence_out = Bio::SeqIO->new( -file   => ">$f_out",
                                       -format => "fasta") || die ($!);

   foreach my $sequence (@$sequences)
   {
      foreach my $feature ($sequence->get_SeqFeatures())
      {
         my @feature_ids = $feature->get_tag_values("ID");
         foreach my $gene_id (@$genes)
         {
            my @hsds_ids = $gene_id->get_tag_values("ID");
            if ($hsds_ids[0] eq $feature_ids[0])
            {
               my $hsds_out = Bio::Seq->new( -seq => $feature->seq->seq(),
                                             -display_id => $hsds_ids[0]);

               $sequence_out->write_seq($hsds_out);
            }
         }
      }
   }

   $sequence_out->close();
}

sub extract_hsds($)
{
   my ($annotation_file) = @_;

   my $tmp_annot_file = "annotation.tmp";
   pre_process_gff($annotation_file, $tmp_annot_file);

   my $gff_in = Bio::Tools::GFF->new(-file => $tmp_annot_file,
                                     -gff_version => 3) || die ("Could not open $tmp_annot_file as gff v3: $!");
   $gff_in->features_attached_to_seqs(1);

   # Hash of arrays, whose reference will be returned
   my %hsd_genes;
   my (@gene_order, @gene_strands);

   # Look through serially, find any hsd relevant genes
   while (my $feature = $gff_in->next_feature())
   {
      # Looking at products is more reliable than gene tags
      if ($feature->has_tag("product"))
      {
         my $strand = $feature->strand();
         push(@gene_strands, $strand);

         my $product = join(",", $feature->get_tag_values("product"));

         # This annotation is part of a type I R-M system
         if ($product =~ /type I restriction/)
         {
            # Specificity subunit
            if ($product =~ /hsds/i || $product =~ /S protein/ || $product =~ /S subunit/)
            {
               push(@{ $hsd_genes{hsdS} }, $feature);
               push(@gene_order, "hsdS");
            }
            # Restriction subunit
            elsif ($product =~ /hsdr/i || $product =~ /R protein/ || $product =~ /R subunit/)
            {
               push(@{ $hsd_genes{hsdR} }, $feature);
               push(@gene_order, "hsdR");
            }
            # Methylation subunit
            elsif ($product =~ /hsdm/i || $product =~ /M protein/ || $product =~ /M subunit/)
            {
               push(@{ $hsd_genes{hsdM} }, $feature);
               push(@gene_order, "hsdM");

            }
         }
         # Possible creX/xerC gene, not too important to be specific here
         elsif($product =~ /recombinase/i)
         {
            push(@gene_order, "recombinase");
         }
         else
         {
            push(@gene_order, "-");
         }
      }
   }

   my @sequences = $gff_in->get_seqs();

   $gff_in->close();

   # Now look through gene order, and try and find hsdS flanked by hsdM and
   # creX/xerC
   my (@confident_hsds_genes, @likely_hsds_genes, @possible_hsds_genes, $hsds_genes);
   my $j = 0;
   for (my $i = 0; $i < scalar(@gene_order); $i++)
   {
      if ($gene_order[$i] eq "hsdS")
      {
         if ($gene_order[$i-(1*$gene_strands[$i])] eq "hsdM" && $gene_order[$i+(1*$gene_strands[$i])] eq "recombinase")
         {
            push(@confident_hsds_genes, ${$hsd_genes{hsdS}}[$j]);
         }
         elsif ($gene_order[$i-(1*$gene_strands[$i])] eq "hsdM" || $gene_order[$i+(1*$gene_strands[$i])] eq "recombinase")
         {
            push(@likely_hsds_genes, ${$hsd_genes{hsdS}}[$j]);
         }
         else
         {
            push(@possible_hsds_genes, ${$hsd_genes{hsdS}}[$j]);
         }

         $j++;
      }
   }

   # Pick most confident existing set
   if (scalar(@confident_hsds_genes) > 0)
   {
      $hsds_genes = \@confident_hsds_genes;
   }
   elsif (scalar(@likely_hsds_genes) > 0)
   {
      $hsds_genes = \@likely_hsds_genes;
   }
   elsif (scalar(@possible_hsds_genes) > 0)
   {
      $hsds_genes = \@possible_hsds_genes;
   }
   else
   {
      die("No plausible hsdS genes!\n");
   }

   return($hsds_genes, \@sequences);
}

sub tblastx($$$)
{
   my ($subject, $query, $output_file) = @_;

   my $blast_command = "tblastx -subject $subject -query $query -outfmt \"6 qseqid sallseqid evalue score\" > $output_file 2> $blast_err";
   system($blast_command);

   # Parse output to extract best hit
   open(BLAST, "$output_file") || die("Could not open $output_file: $!\n");

   my $high_score = 0;
   my $top_hit;

   while (my $blast_line = <BLAST>)
   {
      chomp($blast_line);

      my ($query_id, $subject_id, $e_value, $score) = split("\t", $blast_line);
      if ($score > $high_score)
      {
         $top_hit = $subject_id;
         $high_score = $score;
      }
   }

   return($top_hit);
}

sub print_allele($$)
{
   my ($N_term, $C_term, $naming) = @_;

   print STDERR "Most likely (highest scoring) allele:\n";

   print join("\t", "$N_term$C_term", $allele_map{"$N_term$C_term"}) . "\n";
}


#
# Main
#
my ($map, $assembly, $ref_dir, $forward_reads, $reverse_reads, $annotation_file, $help);
GetOptions ("map"       => \$map,
            "assembly"    => \$assembly,
            "ref_dir=s" => \$ref_dir,
            "forward_reads|f=s" => \$forward_reads,
            "reverse_reads|r=s" => \$reverse_reads,
            "annotation|a=s" => \$annotation_file,
            "help|h"     => \$help
		   ) or die($help_message);

if (!defined($ref_dir))
{
   $ref_dir = ".";
}

if (defined($help))
{
   print $help_message;
}
elsif (defined($map))
{

}
elsif (defined($assembly))
{
   if (!defined($annotation_file) || !-e $annotation_file)
   {
      die("Must set $annotation_file for assembly mode\n");
   }

   print STDERR "Using assembly\n\n";

   # This returns gff features, with attached sequence
   my ($hsds_genes, $sequences) = extract_hsds($annotation_file);

   # Process these objects into a multifasta to use as a blast query
   make_query_fasta($hsds_genes, $sequences, $query_fasta);

   # Run a protein blast for C and N termini
   my $N_term = tblastx("$ref_dir/$N_term_ref", $query_fasta, $N_blast_out);
   if ($N_term =~ /^segment(.)$/)
   {
      $N_term = $1;
   }

   my $C_term = tblastx("$ref_dir/$C_term_ref", $query_fasta, $C_blast_out);
   if ($C_term =~ /^segment(.)$/)
   {
      $C_term = $1;
   }

   print_allele($N_term, $C_term);

}
else
{
   print STDERR "Mode must be set as one of map or assembly\n";
   print STDERR $help_message;
}

exit(0);

