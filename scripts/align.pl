#! /usr/bin/perl

use warnings;
use strict;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use Fcntl;
use File::Basename;

#------------------------------------------------------------------- 
# Daniel Bolanos 2011
# Boulder Language Technologies / University of Colorado at Boulder
#
# description: align data in a master label file
# 
# parameters:
#
#	$ARGV[0] = acoustic models (used for the alignment)
#	$ARGV[1] = folder containing mlf segments to be processed in parallel
#	$ARGV[2] = output folder
#
#-------------------------------------------------------------------

if ((scalar @ARGV) != 7) { die("Error: wrong number of parameters: "); }

# input parameters
my ($filePhoneSet,$fileLexicon,$fileHMMInput,$dirFeatures,$dirMLF,$format,$dirOutput) = @ARGV;

# auxiliar folders
my $dirLog = "$dirOutput/log";
my $dirTemp = "$dirOutput/tmp";
my $dirAlignments = "$dirOutput/ali";
my $dirAlignmentOutput = "$dirOutput/out";

# create the auxiliar folders
system("mkdir -p $dirLog");
system("mkdir -p $dirTemp");
system("mkdir -p $dirAlignments");
system("mkdir -p $dirAlignmentOutput");

# configuration parameters for the estimation

# process each of the segments
opendir(DIR_MLF,"$dirMLF") || die("unable to open the folder: $dirMLF");
my @filesMLF = grep(/\.txt/,readdir(DIR_MLF));
foreach my $fileMLF (sort @filesMLF) {

	# semaphore for exclusive processing of mlf segments 
	my $fileTemp = "$dirTemp/$fileMLF";
  	if (! defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
  		next;
  	}
  	close(FILE); 
  	
	my $cpu = ($fileMLF =~ m/(\d+)/)[0];
	my $base = ($fileMLF =~ m/(.+)\./)[0];
	my $fileMLFPath = "$dirMLF/$fileMLF";
	my $fileOutput = "$dirAlignmentOutput/$base\.out\.txt";   	
	my $fileError = "$dirAlignmentOutput/$base\.err\.txt";
	
	system("aligner -pho $filePhoneSet -mod $fileHMMInput -lex $fileLexicon -for $format -fof $dirFeatures -mlf $fileMLFPath -dir $dirAlignments 1> $fileOutput 2> $fileError");
}
closedir(DIR_MLF);

