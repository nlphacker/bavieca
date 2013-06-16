#! /usr/bin/perl

use warnings;
use strict;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use Fcntl;
use File::Basename;

#------------------------------------------------------------------- 
# Daniel Bolanos 2012
# Boulder Language Technologies 
#
# description: accumulate statistics for Maximum Likelihood Estimation  
# 
# parameters:
#
#-------------------------------------------------------------------

# input parameters
my $mode = $ARGV[0];

my $filePhoneSet;
my $fileLexicon;
my $dirMLF;
my $fileHMMInput;
my $nphones;
my $optionalSymbols;
my $pronOpt;
my $dirFeatures;
my $dirFeaturesAli;
my $dirFeaturesAcc;
my $fileFeaturesConfig;
my $fileFeaturesConfigAli;
my $fileFeaturesConfigAcc;
my $covarianceAcc;
my $dirOutput;

# single-stream accumulation
if ($ARGV[0] eq "single") {
	if ((scalar @ARGV) != 11) {
		die("wrong number of parameters");
	}
	($mode,$filePhoneSet,$fileLexicon,$dirMLF,$fileHMMInput,$nphones,$optionalSymbols,$pronOpt,
		$dirFeatures,$fileFeaturesConfig,$dirOutput) = @ARGV;
}
# double-stream accumulation (useful for single pass retraining)
elsif ($mode eq "double") {
	if ((scalar @ARGV) != 14) {
		die("wrong number of parameters");
	}
	($mode,$filePhoneSet,$fileLexicon,$dirMLF,$fileHMMInput,$nphones,$optionalSymbols,$pronOpt,
		$dirFeaturesAli,$dirFeaturesAcc,$fileFeaturesConfigAli,$fileFeaturesConfigAcc,
		$covarianceAcc,$dirOutput) = @ARGV;
} else {
	die("unrecognized mode: \"$mode\"");
}

# auxiliar folders
my $dirTemp = "$dirOutput/tmp";
my $dirAccumulators = "$dirOutput/acc";
my $dirTrainingOutput = "$dirOutput/out";

# create the auxiliar folders
system("mkdir -p $dirTemp");
system("mkdir -p $dirAccumulators");
system("mkdir -p $dirTrainingOutput");

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
  	
	my $base = ($fileMLF =~ m/(.+)\./)[0];
	my $fileOutput = "$dirTrainingOutput/$base\.out\.txt"; 
	my $fileError = "$dirTrainingOutput/$base\.err\.txt"; 
	my $fileMLFPath = "$dirMLF/$fileMLF";
	my $fileAcc = "$dirAccumulators/$base\.bin";

	my $context = "";
	if ($nphones ne "physical") {
		$context = "-ww $nphones -cw $nphones";
	}
	
	my $fillersPron = "";
	if ($optionalSymbols ne "no") {
		$fillersPron .= "-opt $optionalSymbols";
	}
	if ($pronOpt eq "yes") {
		$fillersPron .= " -pro yes";
	}	
	
	# single-stream accumulation
	if ($mode eq "single") {	
		system("mlaccumulator -pho \"$filePhoneSet\" -mod \"$fileHMMInput\" -lex \"$fileLexicon\" -fea \"$dirFeatures\" -cfg \"$fileFeaturesConfig\" -mlf \"$fileMLFPath\" $context $fillersPron -dAcc \"$fileAcc\" 1> $fileOutput 2> $fileError");		
	}
	# double-stream accumulation 
	elsif ($mode eq "double") {
		system("mlaccumulator -pho \"$filePhoneSet\" -mod \"$fileHMMInput\" -lex \"$fileLexicon\" -fea \"$dirFeaturesAli\" -feaA \"$dirFeaturesAcc\" -cfg \"$fileFeaturesConfigAli\" -cfgA \"$fileFeaturesConfigAcc\" -covA $covarianceAcc -mlf \"$fileMLFPath\" $context $fillersPron -dAcc \"$fileAcc\" 1> $fileOutput 2> $fileError");	
	}
}
closedir(DIR_MLF);

