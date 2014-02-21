#!/usr/bin/perl -w

use strict;

# Set up environment
chdir("/lustre/scratch109/sanger/jl11/rotation_2/omni20/beagle");
system("source ~/.bashrc");

my $b4 = "/nfs/users/nfs_j/jl11/bin/b4.r1230.jar";
my $bcftools = "/usr/bin/bcftools";

my $gt = "TEST.frame";
my $ref = "REF.vcf.gz";
my $chrom = "20";

# Get index number of this process
my $array_num = $ENV{'LSB_JOBINDEX'};

my $input_file = "$gt.$array_num.vcf.gz";

# before the discovery of splitvcf.jar
#my $start = ($array_num - 1)*5*1000000;
#my $end = ($array_num)*5*1000000;

my $start = `$bcftools view -S $input_file | grep -v "#" | head -1 | cut -f 2`;
chomp($start);

my $end = `$bcftools view -S $input_file | tail -1 | cut -f 2`;
chomp($end);

my $out = "beagle.$array_num.vcf.gz";

my $command = "/software/bin/java -Xmx3000m -Xms3000m -jar $b4 gt=$input_file ref=$ref chrom=$chrom:$start-$end impute=true out=$out";
print("Running command: $command\n");
system($command);

exit(0);
