#!/usr/bin/perl -w

#****************************************************************************************#
#* jlmail.pl                                                                            *#
#* Sends an email saying the job has finished                                           *#
#****************************************************************************************#

use strict;

use Getopt::Long;

my $help_message = <<END;
This script should not be run directly, run through jlsub.pl instead
END

#* Read command line options
my ($command_run, $job_number);
GetOptions ("i=s" => \$command_run,
            "j=s" => \$job_number
		   ) or die($help_message);
		   
if (!defined($command_run) || !defined($job_number))
{
	print $help_message;
}
else
{
   # Set up environment
   system("source ~/.bashrc");

   # Find user to send message to
   my $user_name = `whoami`;
   chomp($user_name);

   open(TMP, ">tmpemail.txt") || die("Couldn't open tmpemail.txt");
   print TMP "Command: '$command_run' was completed successfully!\n";
   close TMP;

   # send the message
   system("mail -s 'jlsub job $job_number complete!' $user_name" .'@sanger.ac.uk < tmpemail.txt');
   system("rm tmpemail.txt");
}

exit(0);
