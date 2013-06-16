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
# description: align lattices in order to get phon-boundaries and
#              acoustic log-likelihood  
# 
#-------------------------------------------------------------------

# input parameters
my ($filePhoneSet,$fileLexicon,$fileHMM,$dirFeatures,$dirLattices,$dirLatticesAM) = @ARGV;

# create the directory to keep aligned lattices
system("mkdir -p \"$dirLatticesAM\"");

# create auxiliar directories
system("mkdir -p \"$dirLatticesAM/lattices\"");
system("mkdir -p \"$dirLatticesAM/temp\"");

# process lattices session by session
opendir(FOLDER_SESSIONS, $dirLattices) || die("Could not open folder: $dirLattices");
my @data = grep(/[^\.]+$/,readdir(FOLDER_SESSIONS));
foreach my $dirSession (sort @data) {

	# semaphore for exclusive processing of sessions 
	my $fileTemp = "$dirLatticesAM/temp/$dirSession\.tmp";
  	if (! defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
  		next;
  	}
  	close(FILE);

	BaviecaLattices::amMarking("$dirLattices/$dirSession",
		"$dirLatticesAM/lattices/$dirSession","$dirFeatures/$dirSession/",$filePhoneSet,$fileLexicon,$fileHMM);
}
closedir(FOLDER_SESSIONS);

