#!/usr/local/bin/perl

use strict;
use warnings;
use Cwd;

package TrainingConfig;

use FindBin;                 # locate this script
use lib "$FindBin::Bin/..";  # use the parent directory

#------------------------------------------------------------------- 
# Daniel Bolanos 2012
# Boulder Language Technologies
# 
# description: training definitions for the WSJ task
#
#-------------------------------------------------------------------

# initialization

our $dirBase = Cwd::getcwd();
our $dirCorpus = "<corpus directory, one subdir per speaker containing a raw and a txt file (transcription) per utterance>";

# feature extraction
our $dirFeatures = "<directory where features will be extracted>";
our $fileConfigFeatures = "$dirBase/config/features.cfg";
our $normMode = "utterance";
our $normMethod = "CMN";

# general
our $dirModels = "<base directory to store acoustic models>";
our $filePhoneSet = "$dirBase/config/lexiconCMU/phoneset.txt";
our $fileLexicon = "$dirBase/config/lexiconCMU/lexicon.txt";

# master label files
our $fileMLF = "$dirBase/mlf/mlfAll.txt";			# MLF containing all training data (paths are relative to $dirFeatures)
our $dirMLFSegments = "$dirBase/mlf/segments64";		# directory containing MLF segments (resulting from splitting $fileMLF), at least one segment per cpu u
our $dirMLFSpeakers = "$dirBase/mlf/speakers";			# same, but one MLF per speaker	

# HMM initialization
our $initializationMethod = "flatStart";

# HMM parameters update
our $covariaceUpdate = "diagonal";

# context modeling
our $filePhoneticRules = "$dirBase/config/cmu.rules";		# file containing phonetic rules
our $minimumOccupationCount = 200;
our $minimumLikelihoodGain = 1800;
our $nphones = "triphones";

# accumulation
our $optionalSymbols = "no";		# "no" means only silence insertion
our $pronOpt = "no";			# "yes" means alternative pronunciations will be considered

# mixture splitting
our $gaussianIncrement = 2;
our $epsilon = 0.05;
our $splittingCriterion = "covariance";
our $minimumGaussianOccupation = 100.0;
our $minimumGaussianWeight = 0.00001;


















1;