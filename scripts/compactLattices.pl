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
# description: compact lattices (forward backward determinization) 
# 
#-------------------------------------------------------------------

# input parameters
my ($filePhoneSet,$fileLexicon,$dirLattices,$dirLatticesCompacted) = @ARGV;

# create the directory to keep aligned lattices
system("mkdir -p \"$dirLatticesCompacted\"");

# create auxiliar directories
system("mkdir -p \"$dirLatticesCompacted/lattices\"");
system("mkdir -p \"$dirLatticesCompacted/temp\"");

# process lattices session by session
opendir(FOLDER_SESSIONS, $dirLattices) || die("Could not open directory: $dirLattices");
my @data = grep(/[^\.]+$/,readdir(FOLDER_SESSIONS));
foreach my $dirSession (sort @data) {

	# semaphore for exclusive processing of sessions 
	my $fileTemp = "$dirLatticesCompacted/temp/$dirSession\.tmp";
  	if (! defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
  		next;
  	}
  	close(FILE);

	BaviecaLattices::compact("$dirLattices/$dirSession",
		"$dirLatticesCompacted/lattices/$dirSession",$filePhoneSet,$fileLexicon);
}
closedir(FOLDER_SESSIONS);


	

