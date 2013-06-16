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
# description: generate n-best lists from lattices
#
#-------------------------------------------------------------------

# input parameters
my ($filePhoneSet,$fileLexicon,$n,$ams,$lms,$ip,$ipf,$method,$dirLattices,$dirNBest) = @ARGV;

# create the directory to keep aligned lattices
system("mkdir -p \"$dirNBest\"");

# create auxiliar directories
system("mkdir -p \"$dirNBest/n-bests\"");
system("mkdir -p \"$dirNBest/temp\"");

# process lattices session by session
opendir(FOLDER_SESSIONS, $dirLattices) || die("Could not open folder: $dirLattices");
my @data = grep(/[^\.]+$/,readdir(FOLDER_SESSIONS));
foreach my $dirSession (sort @data) {

	# semaphore for exclusive processing of sessions 
	my $fileTemp = "$dirNBest/temp/$dirSession\.tmp";
  	if (! defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
  		next;
  	}
  	close(FILE);

	BaviecaLattices::nbest("$dirLattices/$dirSession/","$dirNBest/n-bests/$dirSession",
		$filePhoneSet,$fileLexicon,$n,$ams,$lms,$ip,$ipf,$method);
}
closedir(FOLDER_SESSIONS);
