#!/usr/local/bin/perl

package Configuration;

use strict;
use warnings;

#------------------------------------------------------------------- 
# Daniel Bolanos 2011
# Boulder Language Technologies
# 
# description: generic routines to manipulate configuration files
#
#-------------------------------------------------------------------

# creates a suitable copy of a configuration file
sub setConfigurationParameters {

	if ((scalar @_) ne 3) { die("wrong number of parameters"); }
	my($fileConfigurationInput,$fileConfigurationOutput,$paramValue) = @_;
	my %original = %$paramValue;
	my %mapParameterValue = %original; #make a copy
	
	open(READFILE, "<$fileConfigurationInput") || die("unable to open the file: $fileConfigurationInput");
   my @lines = <READFILE>;
	close READFILE;
	
	# override parameters
   open(WRITEFILE,">$fileConfigurationOutput") || die("unable to create the file: $fileConfigurationOutput");
   foreach(@lines) {
   	foreach my $parameterName (keys %mapParameterValue) {
			if (s/^$parameterName =.*/$parameterName = $mapParameterValue{$parameterName}/) {				
				delete $mapParameterValue{$parameterName};
			}
		}		
   	print WRITEFILE $_;
	}
   close WRITEFILE;
   
   # add remaining parameters (if any)
   open(WRITEFILE,">>$fileConfigurationOutput") || die("unable to append to the file: $fileConfigurationOutput");
   print WRITEFILE"\n# dynamically added parameters:"; 
  	foreach my $parameterName (keys %mapParameterValue) {
   	print WRITEFILE "\n$parameterName \= $mapParameterValue{$parameterName}";
	}		
   close WRITEFILE;
}

# return the value of a parameter in the configuration file
sub getParameterValue($$) {

	my($parameterName,$fileConfiguration) = @_;
	
	open(READFILE, "<$fileConfiguration");
   my @lines = <READFILE>;
   foreach my $line (@lines) {
   	if ($line =~ m/^$parameterName\s+=\s+([^\s]+)\s*$/) {
   		return $1;
   	}
	}
	close(READFILE);
	
	die("parameter \"$parameterName\" was not found\n");	
}



1;