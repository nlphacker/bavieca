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
#	$ARGV[0] = configuration file
#  $ARGV[1] = experiment folder
#  $ARGV[2] = acoustic models
#
#-------------------------------------------------------------------

# check parameters
if ((scalar @ARGV) != 3) {
	die("Error: wrong number of parameters");
}

# load parameters
my $fileConfiguration = $ARGV[0];			# configuration file
my $dirExperiment = $ARGV[1];					# experiment folder
my $fileHMM = $ARGV[2];							# speaker independent models

# make sure the experiment folder exists
if (! (-d $dirExperiment)) {
	die("Error: experiment folder $dirExperiment does not exists");
}

my $filePhoneSet = "./languageModel/CMU/phoneset.txt";
my $fileLexicon = "./languageModel/CMU/lexicon20k.txt";
#my $fileLM = "./languageModel/wsj-20k-onp-trigram.arpa";
my $fileLM = "./languageModel/wsj-20k-onp-bigram.arpa";
my $fileFeatureConfig = "./config/features.cfg";

# data partition
my $dirControl = "./control";
my $fileList = "$dirControl/list.txt";
#my $fileList = "$dirControl/listSingle.txt";

# configuration parameters
my %paramValue = ();
$paramValue{"phoneticSymbolSet.file"} = $filePhoneSet;
$paramValue{"acousticModels.file"} = $fileHMM;
$paramValue{"feature.configurationFile"} = $fileFeatureConfig;			
$paramValue{"languageModel.file"} = $fileLM;
$paramValue{"lexicon.file"} = $fileLexicon;

# speaker independent
&BaviecaDecode::si("$dirExperiment/si",$fileConfiguration,\%paramValue,$dirControl,$fileList);
#&BaviecaDecode::fmllr("$dirExperiment/si","$dirExperiment/fmllr1",$fileHMM,$fileConfiguration,\%paramValue,$dirControl,$fileList);
#&BaviecaDecode::mllr("$dirExperiment/fmllr1","$dirExperiment/mllr1",$fileHMM,$fileConfiguration,\%paramValue,$dirControl,$fileList);
