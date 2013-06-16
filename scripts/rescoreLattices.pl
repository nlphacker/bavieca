#! /usr/bin/perl

use warnings;
use strict;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use Fcntl;
use File::Basename;

use FindBin;                                  
use lib "$FindBin::Bin/../scripts";  # use the scripts directory
use BaviecaLattices;

#------------------------------------------------------------------- 
# Daniel Bolanos 2011
# Boulder Language Technologies / University of Colorado at Boulder
#
# description: performs lattice rescoring based on likelihood or 
#	posterior probabilities
# 
#-------------------------------------------------------------------

# input parameters
my ($filePhoneSet,$fileLexicon,$ams,$lms,$ip,$ipf,$method,$dirLattices,$dirHypotheses) = @ARGV;

# create the directory to keep aligned lattices
system("mkdir -p \"$dirHypotheses\"");

# create auxiliar directories
system("mkdir -p \"$dirHypotheses/hypotheses\"");
system("mkdir -p \"$dirHypotheses/temp\"");

# process lattices session by session
opendir(FOLDER_SESSIONS, $dirLattices) || die("Could not open folder: $dirLattices");
my @data = grep(/[^\.]+$/,readdir(FOLDER_SESSIONS));
foreach my $dirSession (sort @data) {

	# semaphore for exclusive processing of sessions 
	my $fileTemp = "$dirHypotheses/temp/$dirSession\.tmp";
  	if (! defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
  		next;
  	}
  	close(FILE);

	BaviecaLattices::rescoring("$dirLattices/$dirSession/","$dirHypotheses/hypotheses/$dirSession",
		$filePhoneSet,$fileLexicon,$ams,$lms,$ip,$ipf,$method);
}
closedir(FOLDER_SESSIONS);
