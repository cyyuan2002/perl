#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Find;

my $version="1.0 Alvin Chen 2011-07-31";

my %files={};
my %dirs;
my $lastdir="";

my %backfiles={};
my @backdirs;

my %opts;
GetOptions(\%opts,"s=s","t=s","help");
if((!defined $opts{s}) || (!defined $opts{t})){
    &Usage();
}

my $parSourceDir=$opts{s};
my $parTargetDir=$opts{t};

$parSourceDir=&foldercheckfull($parSourceDir);
$parTargetDir=&foldercheckfull($parTargetDir);

#print "$parSourceDir\t$parTargetDir\n";


my @sourcedirs=split("\/",$parSourceDir);
my $sourceDir=pop(@sourcedirs);
my $sourceDirlen=length($sourceDir);
my $tarDir="$parTargetDir$sourceDir";
#print "$tarDir\n";
if(!(-e $tarDir)){
	mkdir($tarDir);
	print "Folder: $tarDir created\n";
}

###create links for exist files;
find(\&searchfiles,($parSourceDir));
foreach my $dir(keys %dirs){
	my $dirN=&transfilename($dir,$parTargetDir);
	#print "New Directory: $dirN\n";
	if(!(-e $dirN)){
		mkdir ($dirN);
		print "Folder: $dir created\n";
	}
	if(exists($files{$dir})){
		my %dirfiles=%{$files{$dir}};
		foreach my $file(keys %dirfiles){
			my $fileN=&transfilename($file,$parTargetDir);
			link ($file,$fileN) if(!(-e $fileN));
			print "File: $file linked\n"
		}
	}
}

find(\&searchbackfiles,($tarDir));
my $sourceFD=join("\/",@sourcedirs);
$sourceFD=$sourceFD."\/";
#print "sourceFD: $sourceFD\n";
foreach my $dir(@backdirs){
	my $dirN=&transfilename($dir,$sourceFD);
	## IF directory is deleted, delete all the files and remove directory
	if(!exists($dirs{$dirN})){ 
		if(exists($backfiles{$dir})){
			my %dirfiles=%{$backfiles{$dir}};
			foreach my $file(keys %dirfiles){
				unlink($file);
				print "File: $file unlinked\n";
			}
		}
		rmdir($dir);
		print "Folder: $dir removed\n";
	}
	##IF directory exists
	else{
		if(exists($backfiles{$dir})){
			my %dirfiles=%{$backfiles{$dir}};
			my %sourcefiles=%{$files{$dirN}};
			foreach my $file(keys %dirfiles){
				my $fileN=&transfilename($file,$sourceFD);
				if(!exists($sourcefiles{$fileN})){
					unlink($file);
					print "File: $file unlinked\n";
				}
			}
		}
	}
}

exit(0);

sub foldercheckfull{
	my $foldername=shift;
	my @folderN=split("",$foldername);
	my $lastc=pop(@folderN);
	if($lastc eq "\/"){
		return $foldername;
	}
	else{
		$foldername="$foldername"."\/";
		return $foldername;
	}
}

sub transfilename{
	my ($filename,$tarD)=@_;
	my $subname;
	if($filename=~/$sourceDir/g){
		my $dirpos=pos($filename)-$sourceDirlen;
		$subname=substr($filename,$dirpos);
	}
	my $subfile="$tarD$subname";
	return "$subfile";
}

sub searchbackfiles{
	if(-d){
		#print "Directory: ",$File::Find::name,"\n";
		$lastdir=$File::Find::name;
		push(@backdirs,$lastdir);
	}
	if(-f){
		my $filename=$File::Find::name;
		$backfiles{$lastdir}->{$filename}=1;
	}
}

sub searchfiles{
	if(-d){
		#print "Directory: ",$File::Find::name,"\n";
		$lastdir=$File::Find::name;
		$dirs{$lastdir}=1;
	}
	if(-f){
		my $filename=$File::Find::name;
		$files{$lastdir}->{$filename}=1;
	}
}


sub Usage(){
  print << "    Usage";
	
	This program is used to create/sync the backup folder form source folder. 
	All the files will be created by hard links, which will not spend extra spaces.
	
	Usage:  $0 (version $version)

	<options>
		-s     Source Folder (eg:/home/XXX/Data)
		-t     Target Folder (eg:/home/XXX/backup)
		
    Usage
	exit(0);
}