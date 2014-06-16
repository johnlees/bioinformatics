#!/usr/bin/perl -w

use strict;

# Run with bsub -J "gen[1-26]" -o gen.%J.%I.o -e gen.%J.%I.e -R "select[mem>1000] rusage[mem=1000]" -M1000 ./bed_to_gen.pl

my $input_bed = "../QC/cases-ALS-BPROOF-clean";

# Get index number of this process
my $array_num = $ENV{'LSB_JOBINDEX'};
my $chr = $array_num - 1;

my $plink_command = "plink --noweb --bfile $input_bed --chr $chr --recode --out tmp-ped.$chr";
my $gtool_command = "gtool -P --ped tmp-ped.$chr.ped --map tmp-ped.$chr.map --og cases-ALS-BPROOF-clean.$chr.gen --os cases-ALS-BPROOF-clean.$chr.sample --binary_phenotype --log gtool.$chr.log";

print "Running plink command: $plink_command\n";
system("$plink_command");

print "Running gtool command: $gtool_command\n";
system("$gtool_command");

print "Removing temporary peds and maps\n";
unlink("tmp-ped.$chr.ped");
unlink("tmp-ped.$chr.map");
unlink("tmp-ped.$chr.log");
unlink("gtool.$chr.log");

print "Zipping gen\n";
system("gzip cases-ALS-BPROOF-clean.$chr.gen");

exit(0);

