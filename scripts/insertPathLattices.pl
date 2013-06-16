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
# description: insert a path in the lattice (typically an alignment
#              containing the transcription)
# 
#-------------------------------------------------------------------

# input parameters
my ($filePhoneSet,$fileLexicon,$dirLattices,$dirAlignments,$dirLatticesOutput) = @ARGV;

# create the directory to keep aligned lattices
system("mkdir -p \"$dirLatticesOutput\"");

# create auxiliar directories
system("mkdir -p \"$dirLatticesOutput/lattices\"");
system("mkdir -p \"$dirLatticesOutput/temp\"");

# process lattices session by session
opendir(FOLDER_SESSIONS, $dirLattices) || die("Could not open folder: $dirLattices");
my @data = grep(/[^\.]+$/,readdir(FOLDER_SESSIONS));
foreach my $dirSession (sort @data) {

	# semaphore for exclusive processing of sessions 
	my $fileTemp = "$dirLatticesOutput/temp/$dirSession\.tmp";
  	if (! defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
  		next;
  	}
  	close(FILE);

	BaviecaLattices::insertPath("$dirLattices/$dirSession","$dirAlignments/$dirSession",
		"$dirLatticesOutput/lattices/$dirSession",$filePhoneSet,$fileLexicon);
}
closedir(FOLDER_SESSIONS);

	

