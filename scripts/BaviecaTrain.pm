#! /usr/bin/perl

use strict;
use warnings;
use Fcntl;
use File::Basename;


package BaviecaTrain;

my $cpuCores = 8;
my $dirScripts = File::Basename::dirname(__FILE__);

#------------------------------------------------------------------- 
# Daniel Bolanos 2012
# Boulder Language Technologies
# 
# description: generic functions to perform acoustic model training
#
#-------------------------------------------------------------------

# extract features from raw audio
sub extractFeatures {

	if ((scalar @_) ne 6) { die("wrong number of parameters"); }
	my ($fileConfigFeatures,$dirCorpus,$dirMLFSegments,$dirFeatures,$normMode,$normMethod) = @_;

	system("run.pl $cpuCores $dirScripts/extractFeatures.pl \"$fileConfigFeatures $dirCorpus $dirMLFSegments $dirFeatures $normMode $normMethod\"");
}

# initialize acoustic model parameters
sub initialize {

	if ((scalar @_) ne 6) { die("wrong number of parameters"); }
	my ($dirFeatures,$fileConfigFeatures,$filePhoneSet,$fileLexicon,$fileMLF,$fileHMMInitial) = @_;

	system("hmminitializer -fea $dirFeatures -cfg $fileConfigFeatures -pho $filePhoneSet -lex $fileLexicon -mlf $fileMLF -met flatStart -mod $fileHMMInitial");	
}


# reestimate acoustic model parameters (iterative reestimation)
sub reestimate {

	if ((scalar @_) ne 11) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$dirMLF,$fileHMMBase,$iIterations,$nphones,
		$optionalSymbols,$pronOpt,$dirFeatures,$dirConfigurationFeatures,$dirOutput) = @_;

	# process each of the segments
	my $fileHMMInput = $fileHMMBase;

	for(my $iIteration = 1 ; $iIteration <= $iIterations ; ++$iIteration) {

		# create the output directory for the iteration
		my $dirOutputIteration = "$dirOutput/$iIteration";
		system("mkdir -p $dirOutputIteration");
	
		# accumulate statistics
		my $output = &accumulate($filePhoneSet,$fileLexicon,$dirMLF,$fileHMMInput,$nphones,
			$optionalSymbols,$pronOpt,$dirFeatures,$dirConfigurationFeatures,$dirOutputIteration);
	
		# estimate model parameters from the accumulators
		my $fileHMMOutput = "$dirOutputIteration/models.bin";
		&estimate($filePhoneSet,$fileHMMInput,"$dirOutputIteration/acc",$fileHMMOutput);

		# get ready for next iteration		
		$fileHMMInput = $fileHMMOutput;
		
		printf "iteration: %3d %s\n",$iIteration,$output;
	}
}

# accumulate statistics (single feature stream)
sub accumulate {

	if ((scalar @_) ne 10) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$dirMLF,$fileHMM,$nphones,$optionalSymbols,$pronOpt,
		$dirFeatures,$fileConfigurationFeatures,$dirOutput) = @_;

	my $fileOutput = "$dirOutput/accumulation.txt";	
	system("run.pl $cpuCores $dirScripts/accumulate.pl \"single \"$filePhoneSet\" \"$fileLexicon\" \"$dirMLF\" \"$fileHMM\" $nphones \"$optionalSymbols\" $pronOpt \"$dirFeatures\" \"$fileConfigurationFeatures\" \"$dirOutput\"\"");
	system("cat $dirOutput/out/* > $fileOutput");
		
	return parseOutput($fileOutput);	
}

# estimate HMM parameters
sub estimate {
	
	if ((scalar @_) ne 4) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileHMMInput,$dirAccumulators,$fileHMMOutput) = @_;

	# list accumulators
	my $fileAccumulatorList = "$dirAccumulators/list.txt";
	&listAccumulators($dirAccumulators,$fileAccumulatorList);	

	# estimate model parameters from the accumulators
	system("mlestimator -pho $filePhoneSet -mod $fileHMMInput -acc $fileAccumulatorList -cov 0.05 -out $fileHMMOutput");
}

# cluster logical HMM-states using phonetic contexts
sub clusterContext {

	if ((scalar @_) ne 8) { die("wrong number of parameters"); }
	my ($filePhoneSet,$filePhoneticRules,$nphones,$minimumOccupationCount,$minimumLikelihoodGain,
		$fileHMMInput,$dirAccumulators,$fileHMMOutput) = @_;

	# list accumulators	
	my $fileAccumulatorList = "$dirAccumulators/list.txt";
	&listAccumulators($dirAccumulators,$fileAccumulatorList);
	
	# context clustering
	system("contextclustering -pho $filePhoneSet -mod $fileHMMInput -rul $filePhoneticRules -ww $nphones -cw $nphones -met local -mrg yes -acc $fileAccumulatorList -occ $minimumOccupationCount -gan $minimumLikelihoodGain -out $fileHMMOutput");
}

# mixture increment
sub mixtureIncrement() {

	if ((scalar @_) ne 9) { die("wrong number of parameters"); }
	my ($filePhoneSet,$gaussianIncrement,$splittingCriterion,$minimumGaussianOccupation,
		$minimumGaussianWeight,$epsilon,$fileHMMInput,$dirAccumulators,$fileHMMOutput) = @_;
	
	# list accumulators	
	my $fileAccumulatorList = "$dirAccumulators/list.txt";
	&listAccumulators($dirAccumulators,$fileAccumulatorList); 

	# mixture increment	
	system("gmmeditor -pho $filePhoneSet -mod $fileHMMInput -acc $fileAccumulatorList -inc $TrainingConfig::gaussianIncrement -crt $TrainingConfig::splittingCriterion -occ $TrainingConfig::minimumGaussianOccupation -wgh $TrainingConfig::minimumGaussianWeight -eps $TrainingConfig::epsilon -mrg yes -out $fileHMMOutput");
}

# multiple Gaussian training (Gaussian increment)
sub multipleGaussianTrainingIncrement {

	if ((scalar @_) ne 16) { die("wrong number of parameters"); }
	my ($filePhoneSet,$gaussianIncrement,$splittingCriterion,$minimumGaussianOccupation,
		$minimumGaussianWeight,$epsilon,$optionalSymbols,$pronOpt,$fileLexicon,$dirMLF,
		$fileHMMInputBase,$dirAccBase,$dirFeatures,$fileConfigurationFeatures,$iIterations,$dirOutput) = @_;

	system("mkdir -p $dirOutput");
	my $fileHMMInput = $fileHMMInputBase;
	my $dirAcc = $dirAccBase;
	for(my $i = 1 ; $i <= $iIterations ; ++$i) {
		system("mkdir -p $dirOutput/$i");
		&mixtureIncrement($filePhoneSet,$gaussianIncrement,$splittingCriterion,$minimumGaussianOccupation,
			$minimumGaussianWeight,$epsilon,$fileHMMInput,$dirAcc,"$dirOutput/$i/models_b.bin");
		if ($dirAcc ne $dirAccBase) {
			&removeAccumulators("$dirAcc/list.txt");
		}	
		my $output = &accumulate($filePhoneSet,$fileLexicon,$dirMLF,"$dirOutput/$i/models_b.bin","physical",
			$optionalSymbols,$pronOpt,$dirFeatures,$fileConfigurationFeatures,"$dirOutput/$i");
		&estimate($filePhoneSet,"$dirOutput/$i/models_b.bin","$dirOutput/$i/acc","$dirOutput/$i/models.bin");
		$fileHMMInput = "$dirOutput/$i/models.bin";
		$dirAcc = "$dirOutput/$i/acc";
		
		printf "iteration: %3d %s\n",$i,$output;
	}
	if ($dirAcc ne $dirAccBase) {
		&removeAccumulators("$dirAcc/list.txt");
	}	
}

# create a list with the accumulators in the given folder
sub listAccumulators {

	if ((scalar @_) ne 2) { die("wrong number of parameters"); }
	my ($dirAccumulators,$fileAccumulatorList) = @_;

	# get the accumulators
	opendir(DIR_ACCUMULATORS,$dirAccumulators) || die("unable to open the folder: $dirAccumulators");
	my @filesAccumulators = grep(/\.bin/,readdir(DIR_ACCUMULATORS));
	closedir(DIR_ACCUMULATORS);
	
	# create a list with all the accumulators
	open(FILE_LIST,">$fileAccumulatorList") || die("unable to create file: $fileAccumulatorList");
	foreach my $fileAccumulator (@filesAccumulators) {
		print FILE_LIST "$dirAccumulators/$fileAccumulator\n";
	}
	close(FILE_LIST);
}

# remove all the accumulators in the given list
sub removeAccumulators {

	if ((scalar @_) ne 1) { die("wrong number of parameters"); }
	my ($fileAccumulatorList) = @_;

	open(FILE_LIST,"$fileAccumulatorList") || die("unable to open the file: $fileAccumulatorList");
	foreach my $fileAccumulator(<FILE_LIST>) {
		chomp($fileAccumulator);
		system("rm $fileAccumulator");
	}
	close(FILE_LIST);	
}

# parse the accumulation results 
sub parseOutput() {

	my ($fileOutput) = @_;
	
	if ((-s $fileOutput) == 0) {
		return "empty";
	}

	my $likelihood = 0.0;
	my $gaussian = 0;
	my $secondsTotal = 0;
	my $secondsEstimation = 0;
	my $lines = 0;
	my $rtfTotal = 0;
	open(FILE,"$fileOutput") || die("unable to open the file: $fileOutput");
	foreach my $line (<FILE>) {
		if ($line =~ m/likelihood=\s+([-\d\.]+)/) {
			$likelihood += $1;			
		}
		if ($line =~ m/\s+(\d+)\s+Gauss/) {
			$gaussian = $1;
		}
		if ($line =~ m/\[(\d+):(\d+)'(\d+)''\]\[(\d+):(\d+)'(\d+)''\]/) {
			$secondsTotal += $1*3600+$2*60+$3;
			$secondsEstimation += $4*3600+$5*60+$6;
		}
		# RTF
		if ($line =~ m/\[RTF=([^\]]+)\]/) {
			++$lines;
			$rtfTotal += $1;
		}
	}
	close(FILE);	
	
	my $likelihoodFrame = $likelihood/($secondsEstimation*100);
	my $rtf = $rtfTotal/$lines;

	my $hoursTotal = int($secondsTotal/3600);
	my $minutesTotal = int(($secondsTotal-($hoursTotal*3600))/60);
	$secondsTotal = int(($secondsTotal-($hoursTotal*3600)-($minutesTotal*60)));

	my $hoursEstimation = int($secondsEstimation/3600);
	my $minutesEstimation = int(($secondsEstimation-($hoursEstimation*3600))/60);
	$secondsEstimation = int(($secondsEstimation-($hoursEstimation*3600)-($minutesEstimation*60)));	
	
	return sprintf("likelihood: %18.4f (%5.2f) Gauss: %7d \[RTF=%.4f\]\[%d:%02d'%2d''\]\[%d:%02d'%2d''\]",$likelihood,$likelihoodFrame,$gaussian,$rtf,$hoursTotal,$minutesTotal,$secondsTotal,$hoursEstimation,$minutesEstimation,$secondsEstimation);	
}





1;