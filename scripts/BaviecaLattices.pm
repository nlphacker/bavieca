#! /usr/bin/perl

package BaviecaLattices;

use strict;
use warnings;
use Fcntl;
use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );
use File::Basename;

use FindBin;                                  
use lib "$FindBin::Bin/../scripts";  # use the scripts directory
use Configuration;

my $dirScripts = File::Basename::dirname(__FILE__);
my $cpuCores = 8;

#------------------------------------------------------------------- 
# Daniel Bolanos 2012
# Boulder Language Technologies
# 
# description: generic functions to perform lattice processing
#
#-------------------------------------------------------------------

my $profiler = "";
#my $profiler = "valgrind -v --leak-check=yes";

# language model marking
sub lmMarking {

	if ((scalar @_) ne 7) { die("wrong number of parameters"); }
	my ($dirIn,$dirOut,$filePhoneSet,$fileLexicon,$fileLM,$ip,$ipFiller) = @_;
	
	system("mkdir -p \"$dirOut\"");
	my $fileBatch = "$dirOut/list.bat";
	open(FILE,">$fileBatch") || die("unable to create the batch file"); 
	opendir(DIR,"$dirIn") || die("unable to open the input lattice directory: $dirIn");
	my @files = grep(/\.bin$/,readdir(DIR));
	foreach my $file (@files) {
		chomp($file);
		my $fileLatticeInput = "$dirIn/$file";
		my $fileLatticeOutput = "$dirOut/$file";
		print FILE "$fileLatticeInput\t$fileLatticeOutput\n";
	}	
	closedir(DIR);
	close(FILE);
	
	my $fileOutput = "$dirOut/output.txt";		
	my $fileError = "$dirOut/output.err.txt";		
	system("$profiler latticeeditor -pho \"$filePhoneSet\" -lex \"$fileLexicon\" -lm \"$fileLM\" -ip $ip -ipf $ipFiller -bat \"$fileBatch\" -act lm 1> $fileOutput 2> $fileError");	
}	

# acoustic model marking
sub amMarking {

	if ((scalar @_) ne 6) { die("wrong number of parameters"); }
	my ($dirIn,$dirOut,$dirFeatures,$filePhoneSet,$fileLexicon,$fileHMM) = @_;
	
	system("mkdir -p \"$dirOut\"");
	my $fileBatch = "$dirOut/list.bat";
	open(FILE,">$fileBatch") || die("unable to create the batch file"); 
	opendir(DIR,"$dirIn") || die("unable to open the input lattice directory: $dirIn");
	my @files = grep(/\.bin$/,readdir(DIR));
	foreach my $file (@files) {
		chomp($file);
		my $utteranceID = $file;
		$utteranceID =~ s/\.bin//g;
		my $fileLatticeInput = "$dirIn/$file";
		my $fileLatticeOutput = "$dirOut/$file";
		my $fileFeatures = "$dirFeatures/$utteranceID\.fea";
		if (! (-e $fileFeatures)) { next; }
		print FILE "$fileLatticeInput\t$fileFeatures\t$fileLatticeOutput\n";
	}	
	closedir(DIR);
	close(FILE);
	
	my $fileOutput = "$dirOut/output.txt";		
	my $fileError = "$dirOut/output.err.txt";		

	system("$profiler latticeeditor -pho \"$filePhoneSet\" -lex \"$fileLexicon\" -mod \"$fileHMM\" -bat \"$fileBatch\" -act align 1> $fileOutput 2> $fileError");
}	

# posterior probability + confidence estimates
sub pp {

	if ((scalar @_) ne 9) { die("wrong number of parameters"); }
	my ($dirIn,$dirOut,$filePhoneSet,$fileLexicon,$ams,$lms,$ip,$ipf,$conf) = @_;
	
	system("mkdir -p \"$dirOut\"");
	my $fileBatch = "$dirOut/list.bat";
	open(FILE,">$fileBatch") || die("unable to create the batch file"); 
	opendir(DIR,"$dirIn") || die("unable to open the input lattice directory: $dirIn");
	my @files = grep(/\.bin$/,readdir(DIR));
	foreach my $file (@files) {
		chomp($file);
		my $fileLatticeInput = "$dirIn/$file";
		my $fileLatticeOutput = "$dirOut/$file";
		print FILE "$fileLatticeInput\t$fileLatticeOutput\n";
	}	
	closedir(DIR);
	close(FILE);

	my $fileHypothesis = "$dirOut/hypothesis.hyp";	
	my $fileOutput = "$dirOut/output.txt";		
	my $fileError = "$dirOut/output.err.txt";		
	system("$profiler latticeeditor -pho \"$filePhoneSet\" -lex \"$fileLexicon\" -ams $ams -lms $lms -ip $ip -ipf $ipf -conf $conf -bat \"$fileBatch\" -act pp 1> $fileOutput 2> $fileError");
}
 
# rescoring
sub rescoring {

	if ((scalar @_) ne 9) { die("wrong number of parameters"); }
	my ($dirIn,$dirOut,$filePhoneSet,$fileLexicon,$ams,$lms,$ip,$ipf,$method) = @_;
	
	system("mkdir -p \"$dirOut\"");
	my $fileBatch = "$dirOut/list.bat";
	open(FILE,">$fileBatch") || die("unable to create the batch file"); 
	opendir(DIR,"$dirIn") || die("unable to open the input lattice directory: $dirIn");
	my @files = grep(/\.bin$/,readdir(DIR));
	foreach my $file (@files) {
		chomp($file);
		my $utteranceID = $file;
		$utteranceID =~ s/\.bin//g;
		my $fileLattice = "$dirIn/$file";
		my $fileNBest = "$dirOut/$utteranceID\.txt";
		print FILE "$fileLattice\t$utteranceID\n";
	}	
	closedir(DIR);
	close(FILE);

	my $fileHypothesis = "$dirOut/hypothesis.hyp";	
	my $fileOutput = "$dirOut/output.txt";		
	my $fileError = "$dirOut/output.err.txt";
	system("$profiler latticeeditor -pho \"$filePhoneSet\" -lex \"$fileLexicon\" -ams $ams -lms $lms -ip $ip -ipf $ipf -bat \"$fileBatch\" -res $method -act rescore -hyp \"$fileHypothesis\" 1> $fileOutput 2> $fileError");
}

# n-best generation
sub nbest {

	if ((scalar @_) ne 10) { die("wrong number of parameters"); }
	my ($dirIn,$dirOut,$filePhoneSet,$fileLexicon,$n,$ams,$lms,$ip,$ipf,$method) = @_;
	
	system("mkdir -p \"$dirOut\"");
	my $fileBatch = "$dirOut/list.bat";
	open(FILE,">$fileBatch") || die("unable to create the batch file"); 
	opendir(DIR,"$dirIn") || die("unable to open the input lattice directory: $dirIn");
	my @files = grep(/\.bin$/,readdir(DIR));
	foreach my $file (@files) {
		chomp($file);
		my $utteranceID = $file;
		$utteranceID =~ s/\.bin//g;
		my $fileLattice = "$dirIn/$file";
		my $fileNBest = "$dirOut/$utteranceID\.txt";
		print FILE "$fileLattice\t$fileNBest\n";
	}	
	closedir(DIR);
	close(FILE);
	
	my $fileOutput = "$dirOut/output.txt";
	my $fileError = "$dirOut/output.err.txt";		
	system("$profiler latticeeditor -pho \"$filePhoneSet\" -lex \"$fileLexicon\" -nbest $n -ams $ams -lms $lms -ip $ip -ipf $ipf -bat \"$fileBatch\" -res $method -act nbest 1> $fileOutput 2> $fileError");
}	


# lattice compacting (forward backward determinization)
sub compact {

	if ((scalar @_) ne 4) { die("wrong number of parameters"); }
	my ($dirIn,$dirOut,$filePhoneSet,$fileLexicon) = @_;
	
	system("mkdir -p \"$dirOut\"");
	my $fileBatch = "$dirOut/list.bat";
	open(FILE,">$fileBatch") || die("unable to create the batch file"); 
	opendir(DIR,"$dirIn") || die("unable to open the input lattice directory: $dirIn");
	my @files = grep(/\.bin$/,readdir(DIR));
	foreach my $file (@files) {
		chomp($file);
		my $fileLatticeInput = "$dirIn/$file";
		my $fileLatticeOutput = "$dirOut/$file";
		print FILE "$fileLatticeInput\t$fileLatticeOutput\n";
	}	
	closedir(DIR);
	close(FILE);
	
	my $fileOutput = "$dirOut/output.txt";		
	my $fileError = "$dirOut/output.err.txt";
	system("$profiler latticeeditor -pho \"$filePhoneSet\" -lex \"$fileLexicon\" -bat \"$fileBatch\" -act compact 1> $fileOutput 2> $fileError");
}	

# path insertion
sub insertPath {

	if ((scalar @_) ne 5) { die("wrong number of parameters"); }
	my ($dirIn,$dirAlignments,$dirOut,$filePhoneSet,$fileLexicon) = @_;
	
	system("mkdir -p \"$dirOut\"");
	my $fileBatch = "$dirOut/list.bat";
	open(FILE,">$fileBatch") || die("unable to create the batch file"); 
	opendir(DIR,"$dirIn") || die("unable to open the input lattice directory: $dirIn");
	my @files = grep(/\.bin$/,readdir(DIR));
	foreach my $file (@files) {
		chomp($file);
		my $utteranceID = $file;
		$utteranceID =~ s/\.bin//g;		
		my $fileLatticeInput = "$dirIn/$file";
		my $fileAlignment = "$dirAlignments/$utteranceID\.ali";
		my $fileLatticeOutput = "$dirOut/$file";
		print FILE "$fileLatticeInput\t$fileAlignment\t$fileLatticeOutput\n";
	}	
	closedir(DIR);
	close(FILE);
	
	my $fileOutput = "$dirOut/output.txt";
	my $fileError = "$dirOut/output.err.txt";		
	
	system("$profiler latticeeditor -pho \"$filePhoneSet\" -lex \"$fileLexicon\" -bat \"$fileBatch\" -act addpath 1> $fileOutput 2> $fileError");	
}


# compute Word Error Rate
sub wer {

	if ((scalar @_) ne 5) { die("wrong number of parameters"); }
	my ($dirLattices,$dirOut,$filePhoneSet,$fileLexicon,$fileReference) = @_;
	
	system("mkdir -p \"$dirOut\"");
	my $fileBatch = "$dirOut/list.bat";
	open(FILE,">$fileBatch") || die("unable to create the batch file"); 
	opendir(DIR,"$dirLattices") || die("unable to open the input lattice directory: $dirLattices");
	my @files = grep(/\.bin$/,readdir(DIR));
	foreach my $file (@files) {
		chomp($file);
		my $utteranceID = $file;
		$utteranceID =~ s/\.bin//g;		
		my $fileLattice = "$dirLattices/$file";
		print FILE "$fileLattice\t$utteranceID\n";
	}	
	closedir(DIR);
	close(FILE);
	
	my $fileHypothesis = "$dirOut/hypothesis.hyp";	
	my $fileOutput = "$dirOut/output.txt";
	my $fileError = "$dirOut/output.err.txt";		
	
	system("$profiler latticeeditor -pho \"$filePhoneSet\" -lex \"$fileLexicon\" -bat \"$fileBatch\" -trn \"$fileReference\" -hyp \"$fileHypothesis\" -act wer 1> \"$fileOutput\" 2> \"$fileError\"");	
}



# compact lattices (forward backward determinization)
sub lmMarkingParallel {

	if ((scalar @_) ne 7) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$fileLM,$ip,$ipf,$dirLatticesIn,$dirLatticesOut) = @_;

	system("run.pl $cpuCores $dirScripts/lmMarkingLattices.pl \"$filePhoneSet $fileLexicon $fileLM $ip $ipf $dirLatticesIn $dirLatticesOut\"");
}

# compact lattices (forward backward determinization)
sub amMarkingParallel {

	if ((scalar @_) ne 6) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$fileHMM,$dirFeatures,$dirLatticesIn,$dirLatticesOut) = @_;

	system("run.pl $cpuCores $dirScripts/amMarkingLattices.pl \"$filePhoneSet $fileLexicon $fileHMM $dirFeatures $dirLatticesIn $dirLatticesOut\"");
}

# compact lattices (forward backward determinization)
sub compactParallel {

	if ((scalar @_) ne 4) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$dirLatticesIn,$dirLatticesOut) = @_;

	system("run.pl $cpuCores $dirScripts/compactLattices.pl \"$filePhoneSet $fileLexicon $dirLatticesIn $dirLatticesOut\"");
}

# insert a path in the lattice
sub insertPathParallel {

	if ((scalar @_) ne 5) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$dirLatticesIn,$dirAlignments,$dirLatticesOut) = @_;		

	system("run.pl $cpuCores $dirScripts/insertPathLattices.pl \"$filePhoneSet $fileLexicon $dirLatticesIn $dirAlignments $dirLatticesOut\"");
}

# compute posterior probabilities
sub ppParallel {

	if ((scalar @_) ne 9) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$ams,$lms,$ip,$ipf,$conf,$dirLatticesIn,$dirLatticesOut) = @_;

	system("run.pl $cpuCores $dirScripts/ppLattices.pl \"$filePhoneSet $fileLexicon $ams $lms $ip $ipf $conf $dirLatticesIn $dirLatticesOut\"");
}

# rescore lattices and generate hypothesis
sub rescoreParallel {

	if ((scalar @_) ne 9) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$ams,$lms,$ip,$ipf,$method,$dirLatticesIn,$dirOut) = @_;

	system("run.pl $cpuCores $dirScripts/rescoreLattices.pl \"$filePhoneSet $fileLexicon $ams $lms $ip $ipf $method $dirLatticesIn $dirOut\"");
}

# generate n-best from lattices
sub nbestParallel {

	if ((scalar @_) ne 10) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$n,$ams,$lms,$ip,$ipf,$method,$dirLattices,$dirNBest) = @_;

	system("run.pl $cpuCores $dirScripts/nbestLattices.pl \"$filePhoneSet $fileLexicon $n $ams $lms $ip $ipf $method $dirLattices $dirNBest\"");
}

# compute the lattice Word Error Rate
sub werParallel {

	if ((scalar @_) ne 5) { die("wrong number of parameters"); }
	my ($filePhoneSet,$fileLexicon,$fileReference,$dirLattices,$dirOut) = @_;

	system("run.pl $cpuCores $dirScripts/werLattices.pl \"$filePhoneSet $fileLexicon $fileReference $dirLattices $dirOut\"");
}
1;
