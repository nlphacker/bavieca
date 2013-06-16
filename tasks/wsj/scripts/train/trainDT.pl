use warnings;
use strict;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use Fcntl;
use File::Basename;

use FindBin;                                  
use lib "$FindBin::Bin/../../../../scripts";  # use the scripts directory
use BaviecaTrain;
use BaviecaLattices;

use TrainingConfig;

#------------------------------------------------------------------- 
# Daniel Bolanos 2012
# Boulder Language Technologies / University of Colorado at Boulder
#
# description: discriminative training of acoustic models for the 
#              WSJ task
# 
# parameters:
#
#-------------------------------------------------------------------

my ($dirLattices,$fileLexicon,$fileHMM,$fileLM,$dirFeatures,$dirDT) = @ARGV;
my $filePhoneSet = $TrainingConfig::filePhoneSet;

# create the discriminative training directory
system("mkdir -p $dirDT");

# align transcriptions (a path containing the transcription might need to be added to denominator lattices)
my $dirAlignments = "$dirDT/alignments";
BaviecaTrain::alignMLF($filePhoneSet,$fileLexicon,$fileHMM,$dirFeatures,
	$TrainingConfig::dirMLFSegments,"binary",$dirAlignments);
	
my $dirLatticesCompacted = "$dirDT/lattices/latticesCompacted";
BaviecaLattices::compactParallel($filePhoneSet,$fileLexicon,$dirLattices,$dirLatticesCompacted);

my $dirLatticesBestPath = "$dirDT/lattices/latticesBestPath";
BaviecaLattices::insertPathParallel($filePhoneSet,$fileLexicon,"$dirLatticesCompacted/lattices",
	"$dirAlignments/ali",$dirLatticesBestPath);

my $ip = -18.0;
my $ipFillers = "$TrainingConfig::fillerIP";
my $dirLatticesLM = "$dirDT/lattices/latticesLM";
BaviecaLattices::lmMarkingParallel($filePhoneSet,$fileLexicon,$fileLM,$ip,$ipFillers,
	"$dirLatticesBestPath/lattices","$dirLatticesLM");

my $dirLatticesAM = "$dirDT/lattices/latticesAM";
BaviecaLattices::amMarkingParallel($filePhoneSet,$fileLexicon,$fileHMM,$dirFeatures,
	"$dirLatticesLM/lattices",$dirLatticesAM);
	
# accumulation parameters
my $amScaling = 1.0/25.0;
my $objectiveFunction = "bMMI";
my $boostingFactor = 0.5;
my $cancelation = "yes";
	
# estimation parameters
my $ismoothing = "prev";
my $learningRate = 3.0;
my $tau = 200.0;
my $iterations = 4;

BaviecaTrain::discriminativeTraining($filePhoneSet,$fileLexicon,$fileHMM,$dirFeatures,
	$TrainingConfig::fileConfigFeatures,$TrainingConfig::dirMLFSegments,"$dirLatticesAM/lattices",
	$objectiveFunction,$boostingFactor,$cancelation,$ismoothing,$learningRate,$tau,$amScaling,$iterations,
	"$dirDT/bMMI_${tau}_${learningRate}_${boostingFactor}");
	





