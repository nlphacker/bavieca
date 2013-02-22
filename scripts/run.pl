#! /usr/bin/perl

use warnings;
use strict;

# input parameters
my $times = $ARGV[0];
my $script = $ARGV[1];	     			# perl script 
my $parameters = $ARGV[2];				# parameters passed to the script
if ((scalar @ARGV) != 3) {
	die("wrong number of parameters");
}

# execute command n times
my @args = ($script , $parameters);
for (my $i = 0 ; $i < $times ; ++$i) {
   my $pid = fork();
   if ($pid == -1) {
       die;
   } elsif ($pid == 0) {
      system("perl $script $parameters");
      exit;
   }
}
while (wait() != -1) {}
