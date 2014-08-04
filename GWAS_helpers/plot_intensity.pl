#!/usr/bin/perl -w
#
use strict;
use warnings;

my $intensity_file = $ARGV[0];
my $calls_file = $ARGV[1];
my $rsid = $ARGV[2];

my $R_tmp = "plot.Rscript";
my $plot_tmp = "intensity_plot.tmp";

my $tmp_intensity = "int.tmp";
my $tmp_calls = "calls.tmp";

# R script to do the plot
my $R_script = <<'RSCRIPT';
data = read.table("intensity_plot.tmp",h=T)
AA=which(data$call=="1")
AB=which(data$call=="2")
BB=which(data$call=="3")
NN=which(data$call=="4")
pdf("intensity_plot.pdf")
plot(0,0,pch="",xlim = c(-0.1, 1.5),ylim = c(-0.1, 1.5), xlab="red int", ylab="green int")
points(data$red[AA], data$green[AA],pch=20, col="RED")
points(data$red[AB], data$green[AB],pch=20, col="BLUE")
points(data$red[BB], data$green[BB],pch=20, col="GREEN")
points(data$red[NN], data$green[NN],pch=20, col="GRAY")
dev.off()
quit()
RSCRIPT

if (!defined($intensity_file) || (!-e ($intensity_file)) || !defined($calls_file) || (!-e ($calls_file)))
{
   print "Usage: ./plot_intensity.pl <intensity_file_line.txt> <calls.txt>\n";
}
else
{
   if (defined($rsid))
   {
      system("grep $rsid $calls_file > $tmp_calls");
      system("grep $rsid $intensity_file > $tmp_intensity");

      $intensity_file = $tmp_intensity;
      $calls_file = $tmp_calls;
   }

   open(INT, "$intensity_file") || die("Could not open $intensity_file");
   open(CALLS, "$calls_file") || die("Could not open $calls_file");

   # Read int file in
   my $int_line = <INT>;
   chomp($int_line);

   # Get header line from int file
   my @intensities = split(/\s+/, $int_line);
   $rsid = shift(@intensities);
   my $pos = shift(@intensities);
   my $alleles = shift(@intensities);

   # Read call files in
   my $calls_line = <CALLS>;
   chomp ($calls_line);

   # Remove first four columns containing rsid etc
   my @calls = split(/\s+/, $calls_line);
   splice(@calls, 0, 4);

   close INT;
   close CALLS;

   # Make file for plotting
   open (DATA, ">$plot_tmp") || die("Could not open $plot_tmp for writing");
   print DATA join("\t", "red", "green", "call") . "\n";

   for (my $i=0; $i<scalar(@calls); $i++)
   {
      print DATA join("\t", $intensities[$i*2], $intensities[$i*2+1], $calls[$i]) . "\n";
   }
   close DATA;

   # Make and run R script to make the plot
   open (RSCRIPT, ">$R_tmp") || die ("Could not write to $R_tmp");
   print RSCRIPT $R_script;
   close RSCRIPT;

   system("R CMD BATCH $R_tmp");
   unlink $R_tmp, $plot_tmp, "$R_tmp.Rout", $tmp_calls, $tmp_intensity;

   rename "intensity_plot.pdf", "$rsid" . "_intensity_plot.pdf";
}

exit(0);

