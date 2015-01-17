#!perl -w

#
# Extracts genes which are not present between set of highly
# similar sequences
#

use strict;
use warnings;

# Allows use of perl modules in ./
use Cwd 'abs_path';
use File::Basename;
use lib dirname( abs_path $0 ) . "/..";

use Bio::SeqIO;
use compare_variants;

my $min_ident = 90;
my $blat_log = "blat.log";
my $lane_regex = qr/^(\d+_\d+)#(\d+)$/;

# Rename lane
sub convert_lane($)
{
   my ($lane_in) = @_;

   $lane_in =~ $lane_regex;

   return("$1_$2");
}

# Extract a list of gene ids from a fasta
sub gene_list($)
{
   my ($fasta_file) = @_;

   my @genes;

   my $fasta_in = Bio::SeqIO->new(-file => "$fasta_file",
                                -format => 'Fasta');

   while (my $seq = $fasta_in->next_seq())
   {
      push(@genes,$seq->id);
   }

   $fasta_in->close();
   return(\@genes);
}

# from perl4faq
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

my $pairs_file = $ARGV[0];
my $proteome_dir = $ARGV[1];

open(PAIRS, $pairs_file) || die("Could not open pairs file: $pairs_file\n");

my $header = <PAIRS>;

# Output format to stdout
print join("\t", "Sample", "Tissue", "Missing", "Genes") . "\n";

# For every pair
while (my $line_in = <PAIRS>)
{
   chomp($line_in);
   my ($sample, $csf_lane, $blood_lane) = split("\t", $line_in);

   $csf_lane = convert_lane($csf_lane);
   $blood_lane = convert_lane($blood_lane);

   my %all_genes;
   $all_genes{csf} = gene_list("$proteome_dir/$csf_lane.gff.proteome.faa");
   $all_genes{blood} = gene_list("$proteome_dir/$blood_lane.gff.proteome.faa");

   # Protein BLAT with high identifity
   my $blat_command = "blat -prot out=blast8 -minIdentity=$min_ident $proteome_dir/$csf_lane.gff.proteome.faa $proteome_dir/$blood_lane.gff.proteome.faa $sample.blat.psl >> $blat_log";
   system($blat_command);

   open(BLAT, "$sample.blat.psl") || die("Could not open $sample.blat.psl\n");

   # List of genes matched in samples
   my (%hits, %filtered_hits);
   while (my $blat_line = <BLAT>)
   {
      my ($blood_gene, $csf_gene, @fields) = split("\t", $blat_line);
      push(@{ $hits{blood} }, $blood_gene);
      push(@{ $hits{csf} }, $csf_gene);
   }

   close BLAT;

   # Find those not matched by blat, count and print them
   foreach my $tissue (keys %hits)
   {
      @{ $filtered_hits{$tissue} } = uniq(@{ $hits{$tissue}});

      my $diff_genes = compare_variants::list_compare("diff", \@{ $filtered_hits{$tissue} }, $all_genes{$tissue});

      print join("\t", $sample, $tissue, scalar(@$diff_genes), @$diff_genes) . "\n";
   }

}

close PAIRS;

exit(0);

