#! /usr/bin/perl

package BaviecaDecode;

use strict;
use warnings;
use Fcntl;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use File::Basename;

use FindBin;                                  
use lib "$FindBin::Bin/../scripts";  # use the scripts directory
use Configuration;

my $dirScripts = File::Basename::dirname(__FILE__);

#------------------------------------------------------------------- 
# Daniel Bolanos 2012
# Boulder Language Technologies
# 
# description: generic functions to perform speech recognition
#
#-------------------------------------------------------------------

# speaker independent decoding
sub si {

	if ((scalar @_) ne 5) { die("wrong number of parameters"); }
	my ($dirSI,$fileConfiguration,$paVa,$dirControl,$fileList) = @_;
	my %paramValue = %$paVa;	 

	my $dirConfig = "$dirSI/config";				  # configuration folder
	my $dirHypotheses = "$dirSI/hypotheses";    # hypotheses folder 
	my $dirTemp = "$dirSI/temp";				     # to keep track of the control files already processed
	my $dirOutput = "$dirSI/output";            # output folder
	
	system("mkdir -p $dirConfig");
	system("mkdir -p $dirHypotheses");
	system("mkdir -p $dirTemp");
	system("mkdir -p $dirOutput");
	
	# process each session
  	my %speakerSessions = getSpeakerSessionsMap($fileList);  # get the session <-> speaker mapping	
	foreach my $speaker (keys %speakerSessions) {
		foreach my $fileControl (@{$speakerSessions{$speaker}}) {

			my	$session = $fileControl;
			$session =~ s/\.ctl$//;
			
			# semaphore for exclusive processing of the given control file 
			my $fileTemp = "$dirTemp/$session\.txt";
		  	if (!defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
		  		next;
		  	}
		  	close(FILE);
		   	
	      # create speaker dependent directories
			my $dirLattices = "$dirSI/lattices/$session";         # folder to store the lattices
			my $dirFeatures = "$dirSI/features/$session";         # folder to store the features
			my $dirAlignments = "$dirSI/alignments/$session";     # folder to store the alignments
			system("mkdir -p $dirLattices");
			system("mkdir -p $dirFeatures");
			system("mkdir -p $dirAlignments");      
	
			my $fileConfigurationSpecific = "$dirConfig/$session\.txt";	
	
			# create a suitable copy of the configuration file (with the appropiate language model/lexicon)
			$paramValue{"output.features.folder"} = $dirFeatures;
			$paramValue{"output.alignment.folder"} = $dirAlignments;
			$paramValue{"output.lattice.folder"} = $dirLattices;
					
			Configuration::setConfigurationParameters($fileConfiguration,$fileConfigurationSpecific,\%paramValue);   
		   
			my $fileHypothesis = "$dirHypotheses/$session\.trn";  
		  	my $fileControlPath = "$dirControl/$fileControl";   
			my $fileOutput = "$dirOutput/$session\.out\.txt";
			my $fileError = "$dirOutput/$session\.err\.txt";
		
		  	# actual decoding
			system("dynamicdecoder -cfg $fileConfigurationSpecific -hyp $fileHypothesis -bat $fileControlPath 1> $fileOutput 2> $fileError");
		}
	}  	
	close(FILE_LIST);
}

# fmllr speaker adaptation (feature-space adaptation)
sub fmllr {

	if ((scalar @_) ne 7) { die("wrong number of parameters"); }
	my ($dirBase,$dirFMLLR,$fileHMM,$fileConfiguration,$paVa,$dirControl,$fileList) = @_;
	my %paramValue = %$paVa;	
	
	my $dirWarping = "";
	
	# create auxiliar directories
	my $dirConfig = "$dirFMLLR/config";				   # configuration folder
	my $dirHypotheses = "$dirFMLLR/hypotheses";     # hypotheses folder 
	my $dirTemp = "$dirFMLLR/temp";				      # to keep track of the control files already processed
	my $dirOutput = "$dirFMLLR/output";             # output folder
	my $dirNBest = "$dirFMLLR/nbest";               # folder to store the n-best files
	my $dirFeatures = "$dirFMLLR/features";         # folder to store the features
	my $dirAlignments = "$dirFMLLR/alignments";     # folder to store the alignments
	my $dirTransforms = "$dirFMLLR/transforms";     # model transforms
	my $dirModels = "$dirFMLLR/models";             # adapted models
	my $dirBatch = "$dirFMLLR/batch";               # folder to store the batch files
	
	system("mkdir -p $dirConfig");
	system("mkdir -p $dirHypotheses");
	system("mkdir -p $dirTemp");
	system("mkdir -p $dirOutput");
	system("mkdir -p $dirTransforms");	
	system("mkdir -p $dirModels");	
	system("mkdir -p $dirBatch");	 	
	
	# process each session
  	my %speakerSessions = getSpeakerSessionsMap($fileList);  # get the session <-> speaker mapping	
	foreach my $speaker (keys %speakerSessions) {
		foreach my $fileControl (@{$speakerSessions{$speaker}}) {

			my	$session = $fileControl;
			$session =~ s/\.ctl$//;				  	

			# semaphore for exclusive processing of the given control file 
			# note: the session's speaker is extracted and all sessions for that speaker are processed
			my $fileTemp = "$dirTemp/$session\.txt";
		  	if (!defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
		  		next;
		  	}
		  	close(FILE);		  	
	  	
	  		# create the batch file for fMLLR
	  		my $fileBatchFMLLR = "$dirBatch/$session\.fmllr.txt";
		  	open(FILE_BATCH,">$fileBatchFMLLR") || die("unable to create the batch file for the speaker: $fileBatchFMLLR");
		  	my $dirAlignmentsInput = "$dirBase/alignments";
		  	my $dirFeaturesInput = "$dirBase/features";	
			opendir(SESSION_DIR, "$dirFeaturesInput/$session") || die("unable to open the directory: $dirFeaturesInput/$session");
			my @filesFea = grep(/\.fea$/,readdir(SESSION_DIR));
			foreach my $fileFea (@filesFea) {
				my $fileAli = "$dirAlignmentsInput/$session/$fileFea";
				$fileAli =~ s/\.fea$/\.ali/g;
				$fileFea = "$dirFeaturesInput/$session/$fileFea";
				if (-e $fileAli) {
					print FILE_BATCH "$fileFea $fileAli\n";
				}
			}			
			closedir(SESSION_DIR);			
			close(FILE_BATCH);	
				
			if ((-s $fileBatchFMLLR) == 0) {
				die("empty batch file!");
			}		
	
			# compute the fMLLR transform
			my $fileOutputFMLLR = "$dirOutput/$session\.fmllr\.txt";   	
			my $fileTransformFMLLR = "$dirTransforms/$session\.fmllr.bin";
			my $filePhoneSet = $paramValue{"phoneticSymbolSet.file"};
			system("fmllrestimator -pho $filePhoneSet -for text -mod $fileHMM -tra $fileTransformFMLLR -bat $fileBatchFMLLR > $fileOutputFMLLR");
			my $fileFeatureTransformList = "$dirTransforms/$session\.fmllr.list.txt";
			system("echo \"$fileTransformFMLLR\" > $fileFeatureTransformList");		
			
			# get the warping factor if needed
			my $warpingFactor;
			if ($dirWarping ne "") {
				my $fileWarpingFactor = "$dirWarping/$session\.warpingFactor\.txt";		
				$warpingFactor = &getWarpingFactor($fileWarpingFactor);		
			}
			   	
	      # create session directories
			my $dirLattices = "$dirFMLLR/lattices/$session";         # folder to store the lattices
			my $dirFeatures = "$dirFMLLR/features/$session";         # folder to store the features
			my $dirAlignments = "$dirFMLLR/alignments/$session";     # folder to store the alignments
			system("mkdir -p $dirLattices");
			system("mkdir -p $dirFeatures");
			system("mkdir -p $dirAlignments");      
	
			my $fileConfigurationSpecific = "$dirConfig/$session\.txt";	
	
			# create a suitable copy of the configuration file (with the appropiate language model/lexicon)
			$paramValue{"feature.transformFile"} = $fileFeatureTransformList;			
			$paramValue{"output.features.folder"} = $dirFeatures;
			$paramValue{"output.alignment.folder"} = $dirAlignments;
			#$paramValue{"output.lattice.folder"} = $dirLattices;
			
			Configuration::setConfigurationParameters($fileConfiguration,$fileConfigurationSpecific,\%paramValue);   
		   
			my $fileHypothesis = "$dirHypotheses/$session\.trn";  
		  	my $fileControlPath = "$dirControl/$fileControl";   
			my $fileOutput = "$dirOutput/$session\.out\.txt";
			my $fileError = "$dirOutput/$session\.err\.txt";
		
		  	# actual decoding
			system("dynamicdecoder -cfg $fileConfigurationSpecific -hyp $fileHypothesis -bat $fileControlPath 1> $fileOutput 2> $fileError");
		}	
	}  	
}

# speaker independent decoding
sub mllr {

	if ((scalar @_) ne 7) { die("wrong number of parameters"); }
	my ($dirBase,$dirMLLR,$fileHMM,$fileConfiguration,$paVa,$dirControl,$fileList) = @_;
	my %paramValue = %$paVa;	
	
	my $dirWarping = "";	

	my $dirConfig = "$dirMLLR/config";				 # configuration folder
	my $dirHypotheses = "$dirMLLR/hypotheses";    # hypotheses folder 
	my $dirTemp = "$dirMLLR/temp";				    # to keep track of the control files already processed
	my $dirOutput = "$dirMLLR/output";            # output folder
	# mllr specific
	my $dirRegTree = "$dirMLLR/regTree";           # model transforms
	my $dirTransforms = "$dirMLLR/transforms";     # model transforms
	my $dirModels = "$dirMLLR/models";             # adapted models	
	my $dirBatch = "$dirMLLR/batch";               # folder to store the batch files
	
	system("mkdir -p $dirConfig");
	system("mkdir -p $dirHypotheses");
	system("mkdir -p $dirTemp");
	system("mkdir -p $dirOutput");
	# mllr specific
	system("mkdir -p $dirRegTree");	
	system("mkdir -p $dirTransforms");	
	system("mkdir -p $dirModels");
	system("mkdir -p $dirBatch");		
	
	# create the regression tree (only once)
	my $fileBefore = "$dirRegTree/creatingRT.txt";
	my $fileAfter = "$dirRegTree/createdRT.txt";
	my $fileRegTree = "$dirRegTree/regTree.bin";
	my $filePhoneSet = $paramValue{"phoneticSymbolSet.file"};
  	if (!defined(sysopen(FILE, $fileBefore, O_WRONLY|O_CREAT|O_EXCL))) {
		until(-e $fileAfter) {
			sleep(1);
		}
  	} else {
		my $fileRegTreeOutput = "$dirRegTree/output.txt";
		system("regtree -pho $filePhoneSet -mod $fileHMM -met kMeans -rgc 50 -gau 50 -out \"$fileRegTree\" > $fileRegTreeOutput");			
		system("echo \"created\" > $fileAfter");
  	}
  	close(FILE);
  	
	# process each session
  	my %speakerSessions = getSpeakerSessionsMap($fileList);  # get the session <-> speaker mapping	
	foreach my $speaker (keys %speakerSessions) {
		foreach my $fileControl (@{$speakerSessions{$speaker}}) {

			my	$session = $fileControl;
			$session =~ s/\.ctl$//;
  	
			# semaphore for exclusive processing of the given control file 
			my $fileTemp = "$dirTemp/$session\.txt";
		  	if (!defined(sysopen(FILE, $fileTemp, O_WRONLY|O_CREAT|O_EXCL))) {
		  		next;
		  	}
		  	close(FILE);
		  	
	  		# create the batch file 
	  		my $fileBatchMLLR = "$dirBatch/$session\.mllr.txt";
		  	open(FILE_BATCH,">$fileBatchMLLR") || die("unable to create the batch file for the speaker: $fileBatchMLLR");
		  	my $dirAlignmentsInput = "$dirBase/alignments";
		  	my $dirFeaturesInput = "$dirBase/features";	  	
			opendir(SESSION_DIR, "$dirFeaturesInput/$session") || die("unable to open the directory: $dirFeaturesInput/$session");
			my @filesFea = grep(/\.fea$/,readdir(SESSION_DIR));
			foreach my $fileFea (@filesFea) {
				my $fileAli = "$dirAlignmentsInput/$session/$fileFea";
				$fileAli =~ s/\.fea$/\.ali/g;
				$fileFea = "$dirFeaturesInput/$session/$fileFea";
				if (-e $fileAli) {
					print FILE_BATCH "$fileFea $fileAli\n";
				}
			}			
			closedir(SESSION_DIR);
			close(FILE_BATCH);	
			
			if ((-s $fileBatchMLLR) == 0) {
				die("empty batch file!");
			}
	
			# get the fMLLR transform
			my $fileTransformFMLLR = "$dirBase/transforms/$session\.fmllr.bin";
			my $fileFeatureTransformList = "$dirTransforms/$session\.fmllr.list.txt";
			system("echo \"$fileTransformFMLLR\" > $fileFeatureTransformList");
			
			my $fileOutputMLLR = "$dirOutput/$speaker\.mllr\.txt";   	
			my $fileTransformMLLR = "$dirTransforms/$session\.mllr.bin";
			my $fileModelsOutput = "$dirTransforms/$session\.models.bin";
			
			system("mllrestimator -pho $filePhoneSet -mod $fileHMM -tra $fileTransformMLLR -rgt $fileRegTree -bat $fileBatchMLLR -for text -occ 3500 -gau 50 -cov no -bst no > $fileOutputMLLR");		
			
			# get the likelihood before and after adaptation		
			my $likelihoodIncrease = &getLikelihoodIncrease($fileOutputMLLR);
			my $str = sprintf("%-14s %s",$speaker,$likelihoodIncrease);		
			system("echo \"$str\" >> $dirOutput/increase.txt");
			
			# create the transformed HMMs
			my $fileHMM_MLLR = "$dirModels/$session\.models\.bin";
			system("hmmx -pho $filePhoneSet -tra $fileTransformMLLR -rgt $fileRegTree -in $fileHMM -out $fileHMM_MLLR");		
			
			# get the warping factor if needed
			my $warpingFactor;
			if ($dirWarping ne "") {
				my $fileWarpingFactor = "$dirWarping/$session\.warpingFactor\.txt";		
				$warpingFactor = &getWarpingFactor($fileWarpingFactor);		
			}	  	
			
			# decode the session		
		   	
	      # create speaker dependent directories
			my $dirLattices = "$dirMLLR/lattices/$session";         # folder to store the lattices
			my $dirFeatures = "$dirMLLR/features/$session";         # folder to store the features
			my $dirAlignments = "$dirMLLR/alignments/$session";     # folder to store the alignments
			system("mkdir -p $dirLattices");
			system("mkdir -p $dirFeatures");
			system("mkdir -p $dirAlignments");      
	
			my $fileConfigurationSpecific = "$dirConfig/$session\.txt";	
	
			# create a suitable copy of the configuration file (with the appropiate language model/lexicon)
			$paramValue{"acousticModels.file"} = $fileHMM_MLLR;			
			if (-e $fileTransformFMLLR) {
				$paramValue{"feature.transformFile"} = $fileFeatureTransformList;			
			}
			$paramValue{"output.features.folder"} = $dirFeatures;
			$paramValue{"output.features.folder"} = $dirFeatures;
			$paramValue{"output.alignment.folder"} = $dirAlignments;
			#$paramValue{"output.lattice.folder"} = $dirLattices;
			
			Configuration::setConfigurationParameters($fileConfiguration,$fileConfigurationSpecific,\%paramValue);   
		   
			my $fileHypothesis = "$dirHypotheses/$session\.trn";  
		  	my $fileControlPath = "$dirControl/$fileControl";   
			my $fileOutput = "$dirOutput/$session\.out\.txt";
			my $fileError = "$dirOutput/$session\.err\.txt";
		
		  	# actual decoding
			system("dynamicdecoder -cfg $fileConfigurationSpecific -hyp $fileHypothesis -bat $fileControlPath 1> $fileOutput 2> $fileError");
			
			system("rm $fileHMM_MLLR");
		}	
	}  	
	close(FILE_LIST);
}

# return the likelihood increase on the adaptation data
sub getLikelihoodIncrease {

	if ((scalar @_) ne 1) { die("wrong number of parameters"); }
	my ($file) = @_;
	
	open(FILE,$file) || die("unable to open the file: $file");
	my @like = ();
	foreach my $line (<FILE>) {
		if ($line =~ m/total likelihood\:\s+([\d\.\-]+)/) {
			push(@like,$1);
		}
	}
	close(FILE);
	
	my $data = "";
	if ((scalar @like) == 2) {
		my $increase = 1.0-($like[1]/$like[0]);
		$data .= sprintf "%12.2f %12.2f %6.2f",$like[0],$like[1],$increase;
	}
	
	return $data;
}

# get the warping factor from the file
sub getWarpingFactor {

	if ((scalar @_) ne 1) { die("wrong number of parameters"); }
	my ($file) = @_;

	# get the warp factor from the output
	my $warpFactor = 1.0;
	open(READFILE, "$file");
  	my @lines = <READFILE>;
  	foreach my $line (@lines) {
		if ($line =~ m/^global warp factor: ([\d\.]+)\s.*$/) {
			$warpFactor = $1;
		}
	}
	close(READFILE);
	
	return $warpFactor;
}

# score a recognition experiment using sclite
sub score {

	if ((scalar @_) ne 3) { die("wrong number of parameters"); }
	my ($dirExperiment,$stage,$fileReference) = @_;

	# open the experiments folder
	opendir(INPUT_DIR, $dirExperiment) || die("unable to open the directory: $dirExperiment");
	my @dirs = grep(/\d+/,readdir(INPUT_DIR));
	
	# compute the WER for each experiment subfolder 
	foreach my $dir (sort @dirs) {
	
		# only model folders
		my $iteration = $dir;
		printf("iteration: %3s",$iteration);		
	
		my $dirHypotheses = "$dirExperiment/$dir/$stage/hypotheses";
		if (! -d $dirHypotheses) {
			next;
		}	
		
		my $dirScoring = "$dirExperiment/$dir/$stage/scoring";
		system("mkdir -p \"$dirScoring\"");
		
		my $fileOutput = "$dirScoring/$stage\_hyp.trn";
		my $fileAlignmentOutput = "$dirScoring/$stage\_sclite.out";
		
		# remove the contents of the scoring folder (just in case)
		if (-d $dirScoring) {
	      system("rm -rf $fileOutput");
	      system("rm -rf $fileAlignmentOutput");
		}
		
		# create a file with the global hypothesis
		opendir(FOLDER_SPEAKERS, $dirHypotheses) || die("Could not open folder: $dirHypotheses");
		my @data = grep(/[^\.]+$/,readdir(FOLDER_SPEAKERS));
		foreach my $file (sort @data) {
	
			chomp($file);			
		   my $fileHypothesis = "$dirHypotheses/$file";		  
		  	if (-e $fileHypothesis) {
		   	system("cat \"$fileHypothesis\" >> \"$fileOutput\"");
		   }
		}
		closedir(FOLDER_SPEAKERS);
	
		# score the hypothesis
		my $fileDetailedScoring = "${fileAlignmentOutput}.detailed.txt";
		system("sclite -h $fileOutput trn -r $fileReference trn -i wsj -o dtl > \"$fileDetailedScoring\"");				
		system("sclite -h $fileOutput trn -r $fileReference trn -i wsj > \"$fileAlignmentOutput\"");		
	
		&extractWERSclite($fileAlignmentOutput);
	}   
	closedir(INPUT_DIR);
}

# parse sclite output and get WER
sub extractWERSclite {

	my($fileSclite) = @_;

	# open the clite output
	open(INFO, $fileSclite) || die("Could not open file: $fileSclite");
	my @lines = <INFO>;
	foreach my $line (@lines) {
     	if ($line =~ m/\s*\|\s*Sum\/Avg.*/) {
			print $line;     	
     	} 
	}	 	
	close(INFO);
}

# return the session speaker map
sub getSpeakerSessionsMap {

	if ((scalar @_) ne 1) { die("wrong number of parameters"); }
	my ($file) = @_;
	
	my %speakerSessions = ();
	open(FILE,$file) || die("unable to open the file: $file");
	foreach my $line (<FILE>) {
		if ($line =~ m/^([^\s]+)\s+([^\s]+)\s*/) {
			my $fileControl = $1; 
			my $speaker = $2;
			push(@{$speakerSessions{$speaker}},$fileControl);
		}
	}
	close(FILE);
	
	return %speakerSessions;
}


1;



