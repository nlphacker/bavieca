#! /usr/bin/perl

use warnings;
use strict;
use Fcntl;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use File::Basename;

#------------------------------------------------------------------- 
# Daniel Bolanos 2012
# Boulder Language Technologies
#
# description: extract features from the training data
# 
# parameters:
#
#	$ARGV[0] = configuration file
#	$ARGV[1] = corpus directory (raw audio)
#	$ARGV[2] = master label file
#	$ARGV[3] = feature directory
#	$ARGV[4] = normalization mode (utterance, session, etc)
#	$ARGV[5] = normalization method (CMN, CMVN, etc)
#	$ARGV[6] = warping factors (optional)
#
#-------------------------------------------------------------------

# input parameters
my ($fileConfiguration,$dirCorpus,$dirMLF,$dirFeatures,$normMode,$normMethod) = @ARGV;
my $fileWarpingFactors = "";
if ((scalar @ARGV) > 6) {
	$fileWarpingFactors = $ARGV[6];
}

# create the feature directory
if (-d $dirFeatures) {
	#die("feature directory \"$dirFeatures\" already exists!");
}

my $dirTemp = "$dirFeatures/tmp";
my $dirBatch = "$dirFeatures/batch";

system("mkdir -p $dirFeatures");
system("mkdir -p $dirTemp");
system("mkdir -p $dirBatch");

# load the warping factors if necessary
my %warpingFactors = ();
if ($fileWarpingFactors ne "") {
	open(FILE,$fileWarpingFactors) || die("unable to open the file: $fileWarpingFactors");
	foreach my $line (<FILE>) {
		if ($line =~ m/^([^\s]+)\s+([^\s]+)/) {
			my $speaker = $1;
			my $warpingFactor = $2;
			$warpingFactors{$speaker} = $warpingFactor;
		}	
	}	
	close(FILE);
}

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
  	
  	# create the batch file
	open(FILE,"$dirMLF/$fileMLF") || die("unable to open the file: $fileMLF");
	my $iLines = 0;
	my $speaker;
	my $warpingFactor = 1.0;
	my %speakerData = ();
	my $fileBatch = "$dirBatch/$fileMLF";
	open(FILE_BATCH,">$fileBatch") || die("unable to create the batch file: $fileBatch");
	foreach my $line (<FILE>) {
		if ($line =~ m/^\"\/([^\/]+)\/.+\"/) {
			$speaker = $1;
			my $fileFea = ($line =~ m/\"(.+)\"/)[0];
			my $fileFeaPath = "$dirFeatures/$fileFea";
			my $fileRawPath = "$dirCorpus/$fileFea";
			$fileRawPath =~ s/\.fea/\.raw/g;
			
			my $dirFea = ($line =~ m/\"(.+)\/[^\/]+/)[0];
			$dirFea = "$dirFeatures/$dirFea";
			system("mkdir -p $dirFea");			
			
			if (defined($warpingFactors{$speaker})) {
				$warpingFactor = $warpingFactors{$speaker};
			}
			
			print FILE_BATCH "$fileRawPath $fileFeaPath\n";			
		}	
	}
	close(FILE_BATCH);
	close(FILE);
	
	# actual parameterization
	system("param -cfg \"$fileConfiguration\" -bat \"$fileBatch\" -wrp $warpingFactor -nrm $normMode -met $normMethod -hlt no");		
}
closedir(DIR_MLF);







