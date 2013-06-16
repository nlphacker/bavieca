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
# description: compute posterior probabilities over lattices
#
#-------------------------------------------------------------------

# input parameters
my ($filePhoneSet,$fileLexicon,$ams,$lms,$ip,$ipf,$conf,$dirLatticesInput,$dirLatticesOutput) = @ARGV;

# create the directory to keep aligned lattices
system("mkdir -p \"$dirLatticesOutput\"");

# create auxiliar directories
system("mkdir -p \"$dirLatticesOutput/lattices\"");
system("mkdir -p \"$dirLatticesOutput/temp\"");

# process lattices session by session
opendir(FOLDER_SESSIONS, $dirLatticesInput) || die("Could not open folder: $dirLatticesInput");
my @data = grep(/[^\.]+$/,readdir(FOLDER_SESSIONS));
foreach my $dirSession (sort @data) {

	# semaphore for exclusive processing of sessions 
	my $fileTemp = "$dirLatticesOutput/temp/$dirSession\.tmp";
  	if (! defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
  		next;
  	}
  	close(FILE);

	BaviecaLattices::pp("$dirLatticesInput/$dirSession/","$dirLatticesOutput/lattices/$dirSession",
		$filePhoneSet,$fileLexicon,$ams,$lms,$ip,$ipf,$conf);
	
}
closedir(FOLDER_SESSIONS);
