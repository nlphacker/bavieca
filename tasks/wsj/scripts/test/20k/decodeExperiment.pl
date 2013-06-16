#!/usr/local/bin/perl

use strict;
use warnings;

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
#	$ARGV[0] = acoustic models directory
#  $ARGV[1] = experiment id
#
#-------------------------------------------------------------------

my $cpuCores = 4;

# load parameters
my $dirModels = $ARGV[0];						# dir containing acoustic models
my $experimentId = $ARGV[1];					# experiment id

# check parameters
if ((scalar @ARGV) != 2) {
	die("Error: wrong number of parameters");
}

my $iterationStart = 14;
my $iterationEnd = 26;
my $dirExperiment = "$WSJ::dirData/experiments/$experimentId";
my $fileConfiguration = "./config/configuration.txt";   

# create the experiment directory
system("mkdir -p $dirExperiment");

# copy the configuration file
my $fileConfigurationExperiment = "$dirExperiment/configuration.txt";
system("cp \"$fileConfiguration\" \"$fileConfigurationExperiment\"");

# run the decoding for acoustic model
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










