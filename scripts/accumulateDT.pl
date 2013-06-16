#! /usr/bin/perl

use warnings;
use strict;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use Fcntl;
use File::Basename;

#------------------------------------------------------------------- 
# Daniel Bolanos 2012
# Boulder Language Technologies / University of Colorado at Boulder
#
# description: parallel accumulation of statistics for discriminative
#              training
# 
# parameters:
#
#-------------------------------------------------------------------

# input parameters
if ((scalar @ARGV) != 13) {
	die("wrong number of parameters");
}
my ($filePhoneSet,$fileLexicon,$fileHMM,$nphones,$dirFeatures,$fileFeaturesConfig,
$dirLattices,$objectiveFunction,$boostingFactor,$cancelation,$amScaling,$dirMLFSegments,$dirOutput) = @ARGV;

# auxiliar folders
my $dirLog = "$dirOutput/log";
my $dirTemp = "$dirOutput/tmp";
my $dirMLF = $dirMLFSegments;
my $dirAccumulators = "$dirOutput/acc";
my $dirTrainingOutput = "$dirOutput/out";

# create the auxiliar folders
system("mkdir -p $dirTemp");
system("mkdir -p $dirAccumulators/num");
system("mkdir -p $dirAccumulators/den");
system("mkdir -p $dirTrainingOutput");

# configuration parameters for the estimation

# process each of the segments
opendir(DIR_MLF,"$dirMLF") || die("unable to open the folder: $dirMLF");
my @filesMLF = grep(/\.txt/,readdir(DIR_MLF));
foreach my $fileMLF (sort @filesMLF) {

	#print "processing mlf: $fileMLF $nphones\n";

	# semaphore for exclusive processing of mlf segments 
	my $fileTemp = "$dirTemp/$fileMLF";
  	if (! defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
  		next;
  	}
  	close(FILE);
  	
	my $base = ($fileMLF =~ m/(.+)\./)[0];
	my $fileOutput = "$dirTrainingOutput/$base\.out\.txt";   	
	my $fileError = "$dirTrainingOutput/$base\.err\.txt";   	
	my $fileMLFPath = "$dirMLF/$fileMLF";
	my $fileAccNum = "$dirAccumulators/num/$base\.num\.bin";
	my $fileAccDen = "$dirAccumulators/den/$base\.den\.bin";
	my $context = "";
	if ($nphones ne "physical") {
		$context = "-ww $nphones -cw $nphones";
	}
	
	system("dtaccumulator -pho \"$filePhoneSet\" -mod \"$fileHMM\" -lex \"$fileLexicon\" -ams \"$amScaling\" -fea \"$dirFeatures\" -cfg \"$fileFeaturesConfig\" -mlf \"$fileMLFPath\" -lat \"$dirLattices\" $context -dAccNum \"$fileAccNum\" -dAccDen \"$fileAccDen\" -obj $objectiveFunction -bst $boostingFactor -can $cancelation 1> $fileOutput 2> $fileError");
}
closedir(DIR_MLF);

