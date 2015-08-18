#!/usr/bin/perl -w
#

use strict;
use warnings;

use POSIX;

# list diff, union, intersection
sub list_compare($$$)
{
   my ($mode, $list1, $list2) = @_;

   my %count;
   my (@union, @intersection, @diff);

   foreach my $element (@$list1, @$list2)
   {
      $count{$element}++;
   }
   foreach my $element (keys %count)
   {
      push(@union, $element);

      if ($count{$element} > 1)
      {
         push(@intersection, $element);
      }
      else
      {
         push(@diff, $element);
      }
   }

   my $return_ref;
   if ($mode eq "diff")
   {
      $return_ref = \@diff;
   }
   elsif ($mode eq "intersection")
   {
      $return_ref = \@intersection;
   }
   elsif ($mode eq "union")
   {
      $return_ref = \@union;
   }
   else
   {
      die("Invalid list mode $mode\n");
   }

   return($return_ref);
}

sub in_list($$)
{
   my ($search, $array_ref) = @_;

   my $in = 0;
   my $i = 0;
   foreach my $ref_el (@$array_ref)
   {
      if ($search eq $ref_el)
      {
         $in = $i+1;
         last;
      }
      $i++;
   }

   return $in;
}

open(TIME, ">raxml_time.txt") || die("Couldn't write to raxml_time.txt\n");

my $total_time = 0;
my $total_jobs = 0;
my $i = 0;

my @jobids;
while(1)
{
   $i++;
   my @new_jobs;

   my $jobs = `bjobs -G team81 -o "jobid job_name run_time cmd delimiter=','" | grep raxml`;
   chomp $jobs;

   my @job_lines = split("\n", $jobs);

   foreach my $job_line (@job_lines)
   {
      my ($jobid, $job_name, $run_time, $crapola) = split(",", $job_line);

      if ($job_name =~ /\[(\d+)\]$/)
      {
         $jobid = "$jobid\[$1\]";
      }

      if (!in_list($jobid, \@jobids))
      {
         push(@jobids, $jobid);
      }
      push(@new_jobs, $jobid);
   }

   my $done_jobs = list_compare("diff", \@jobids, \@new_jobs);
   foreach my $done_job (@$done_jobs)
   {
      print STDERR "$done_job finished\n";

      # Remove job from checked list
      splice(@jobids,in_list($done_job, \@jobids)-1, 1);

      my $time_return = `bjobs -noheader -o "run_time" $done_job`;
      chomp $time_return;

      $time_return =~ m/^(\d+) /;

      $total_time += $1;
      $total_jobs++;
   }

   print TIME join("\t", strftime("%Y-%m-%d_%H-%M-%S", localtime), $total_time, $total_jobs, scalar(@jobids)) . "\n";

   print STDERR "loop $i done\n";
   sleep(600);

}

