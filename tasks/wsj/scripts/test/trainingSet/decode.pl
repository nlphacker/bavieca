#!/usr/local/bin/perl

use strict;
use warnings;

use FindBin;                                  
use lib "$FindBin::Bin/../../../../../scripts";  # use the scripts directory
use BaviecaDecode;


#------------------------------------------------------------------- 
# Daniel Bolanos 2011
# Boulder Language Technologies / University of Colorado at Boulder
# 
# description: decoding
#
# parameters:
#
#  $ARGV[0] = configuration file
#  $ARGV[1] = experiment directory
#  $ARGV[2] = acoustic models (typically resulting from MLE)
#
#-------------------------------------------------------------------

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use Fcntl;
use File::Basename;

# check parameters
if ((scalar @ARGV) != 3) { die("Error: wrong number of parameters"); }

# load parameters
my $fileConfiguration = $ARGV[0];			# configuration file
my $dirExperiment = $ARGV[1];					# experiment folder
my $fileHMM = $ARGV[2];							# speaker independent models

# make sure the experiment folder exists
if (! (-d $dirExperiment)) { die("Error: experiment folder $dirExperiment does not exist"); }

# configuration parameters
my $fileLexicon = "../../set-up/lm64/lexicon64.txt";
my $fileLM = "../../set-up/lm64/lm.arpa";
#my $fileLM = "../../set-up/lm64/lm-unigram.arpa";

my $fileFeatureConfig = "./config/features.cfg";
my $filePhoneSet = "./languageModel/CMU/phoneset.txt";

# data partition
my $dirControl = "../../set-up/control";
my $fileList = "$dirControl/list.txt";

# configuration parameters
my %paramValue = ();
$paramValue{"phoneticSymbolSet.file"} = $filePhoneSet;
$paramValue{"acousticModels.file"} = $fileHMM;
$paramValue{"feature.configurationFile"} = $fileFeatureConfig;			
$paramValue{"languageModel.file"} = $fileLM;
$paramValue{"lexicon.file"} = $fileLexicon;

# speaker independent
&BaviecaDecode::si("$dirExperiment/si",$fileConfiguration,\%paramValue,$dirControl,$fileList);

