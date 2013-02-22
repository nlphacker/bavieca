#! /usr/bin/perl

use warnings;
use strict;
use POSIX;

use FindBin;                                  
use lib "$FindBin::Bin/../../../../scripts";  # use the scripts directory
use BaviecaTrain;

use TrainingConfig;

#------------------------------------------------------------------- 
# Daniel Bolanos 2012
# Boulder Language Technologies
#
# description: simple script for Maximum Likelihood based acoustic
#              model training 
# 
# parameters:
#
#	$ARGV[0] = acoustic model folder
#
#-------------------------------------------------------------------

use File::Basename;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use Fcntl;

# check parameters
if (@ARGV != 1) {
	die("wrong number of parameters, please use: train.pl [acoustic model directory]");
}

# create the training folder
my $dirTraining = "$TrainingConfig::dirModels/$ARGV[0]";
system("mkdir -p $dirTraining");

# (0) extract features
BaviecaTrain::extractFeatures($TrainingConfig::fileConfigFeatures,$TrainingConfig::dirCorpus,$TrainingConfig::dirMLFSegments,$TrainingConfig::dirFeatures,$TrainingConfig::normMode,$TrainingConfig::normMethod);
BaviecaTrain::extractFeatures($TrainingConfig::fileConfigFeatures,$TrainingConfig::dirCorpus,$TrainingConfig::dirMLFSpeakers,$TrainingConfig::dirFeatures,$TrainingConfig::normMode,$TrainingConfig::normMethod);

# (1) initialize the HMM parameters to the global distribution of the data (flat-start)
system("mkdir $dirTraining/init");
my $fileHMMInitial = "$dirTraining/init/models00.bin";
BaviecaTrain::initialize($TrainingConfig::dirFeatures,$TrainingConfig::fileConfigFeatures,$TrainingConfig::filePhoneSet,$TrainingConfig::fileLexicon,$TrainingConfig::fileMLF,$fileHMMInitial);

# (2) HMM parameter estimation (iterative reestimation)

# - context independent (CI) single Gaussian estimation
my $dirCI = "$dirTraining/ci";
system("mkdir $dirCI");
BaviecaTrain::reestimate($TrainingConfig::filePhoneSet,$TrainingConfig::fileLexicon,
	$TrainingConfig::dirMLFSegments,$fileHMMInitial,4,"physical",$TrainingConfig::optionalSymbols,
	$TrainingConfig::pronOpt,$TrainingConfig::dirFeatures,$TrainingConfig::fileConfigFeatures,"$dirCI");
system("mkdir $dirCI/5");
my $output = BaviecaTrain::accumulate($TrainingConfig::filePhoneSet,$TrainingConfig::fileLexicon,
	$TrainingConfig::dirMLFSegments,"$dirCI/4/models.bin",$TrainingConfig::nphones,
	$TrainingConfig::optionalSymbols,$TrainingConfig::pronOpt,$TrainingConfig::dirFeatures,
	$TrainingConfig::fileConfigFeatures,"$dirCI/5");
printf "iteration: %3d %s\n",5,$output;

# - context modeling
BaviecaTrain::clusterContext($TrainingConfig::filePhoneSet,$TrainingConfig::filePhoneticRules,
	$TrainingConfig::nphones,$TrainingConfig::minimumOccupationCount,
	$TrainingConfig::minimumLikelihoodGain,"$dirCI/4/models.bin","$dirCI/5/acc","$dirCI/5/models.bin");

# - context dependent (CD) single Gaussian estimation
my $dirCD = "$dirTraining/cd";
system("mkdir $dirCD");
BaviecaTrain::reestimate($TrainingConfig::filePhoneSet,$TrainingConfig::fileLexicon,
	$TrainingConfig::dirMLFSegments,"$dirCI/5/models.bin",3,"physical",$TrainingConfig::optionalSymbols,
	$TrainingConfig::pronOpt,$TrainingConfig::dirFeatures,$TrainingConfig::fileConfigFeatures,"$dirCD");

# - multiple Gaussian estimation
BaviecaTrain::multipleGaussianTrainingIncrement($TrainingConfig::filePhoneSet,
	$TrainingConfig::gaussianIncrement,$TrainingConfig::splittingCriterion,
	$TrainingConfig::minimumGaussianOccupation,$TrainingConfig::minimumGaussianWeight,$TrainingConfig::epsilon,
	$TrainingConfig::optionalSymbols,$TrainingConfig::pronOpt,$TrainingConfig::fileLexicon,
	$TrainingConfig::dirMLFSegments,"$dirCD/3/models.bin","$dirCD/3/acc",$TrainingConfig::dirFeatures,
	$TrainingConfig::fileConfigFeatures,26,"$dirTraining/cd_m_increment");



