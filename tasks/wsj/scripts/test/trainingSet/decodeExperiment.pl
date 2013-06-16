#!/usr/local/bin/perl

use strict;
use warnings;

use FindBin;                 # locate this script
use lib "$FindBin::Bin/../../";  # use the parent directory
use WSJ;
use Configuration;

use Cwd 'abs_path';

#------------------------------------------------------------------- 
# Daniel Bolanos 2010
# Boulder Language Technologies / University of Colorado at Boulder
# 
# description: decoding (multiple processors) 
#
# parameters:
#
#  $ARGV[0] = folder containing acoustic models
#  $ARGV[1] = experiment id
#
#-------------------------------------------------------------------

my $cpuCores = 8;

# load parameters
my $dirModels = $ARGV[0];						# dir containing acoustic models
my $experimentId = $ARGV[1];					# experiment id

# check parameters
if ((scalar @ARGV) != 2) {
	die("Error: wrong number of parameters");
}

my $iterationStart = 17;
my $iterationEnd = 17;
my $dirExperiment = "$WSJ::dirData/experiments/$experimentId";
#my $fileConfiguration = "./config/configurationBigram.txt";   
my $fileConfiguration = "./config/configurationUnigram.txt";   

# create the experiment directory
system("mkdir -p $dirExperiment");

# copy the configuration file
my $fileConfigurationExperiment = "$dirExperiment/configuration.txt";
system("cp \"$fileConfiguration\" \"$fileConfigurationExperiment\"");

# run the decoding for each acoustic model
opendir(INPUT_DIR, $dirModels) || die("unable to open the directory: $dirModels");
my @subdirs = grep(/^\d+$/,readdir(INPUT_DIR));
foreach my $subdir (sort {$a <=> $b} @subdirs) {

	# process subdir? 
	my $iteration = $subdir;
	if (($iteration < $iterationStart) || ($iteration > $iterationEnd)) {
		next;
	}

	my $fileHMM = "$dirModels/$subdir/models.bin";
	if (!(-e $fileHMM)) {
		die("acoustic model file: \"$fileHMM\" was not found");
	}

	# create the decoding folder
	my $dirDecoding = "$dirExperiment/$subdir";
	system("mkdir -p $dirDecoding");		
		
	# run the decoding in multiple processors
	system("run.pl $cpuCores ./decode.pl \"\"$fileConfigurationExperiment\" \"$dirDecoding\" \"$fileHMM\"\"");

}
closedir(INPUT_DIR);

