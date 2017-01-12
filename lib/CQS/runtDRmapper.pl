use strict;
use warnings;

use File::Basename;

my $outFile=$ARGV[0];

#my $perlFile="/scratch/cqs/zhaos/vickers/otherPipeline/tDRmapper/Scripts/TdrMappingScripts.pl";
#my $faDatabaseFile="hg19_mature_and_pre";

my $perlFile=$ARGV[1];
my $faDatabaseFile=$ARGV[2];
my $fastqFile=$ARGV[3];

my $fastqFileInCurrentFolder=basename($fastqFile);
if ($fastqFile=~/.gz$/) {
	print "Now Unzip $fastqFile \n";
	$fastqFileInCurrentFolder=~s/.gz$//;
	`zcat $fastqFile >$fastqFileInCurrentFolder`
} else {
	symlink ( $fastqFile, $fastqFileInCurrentFolder );
}

print "Now Running tDRmapper \n";
`perl $perlFile $faDatabaseFile $fastqFileInCurrentFolder` if (! -e "$fastqFileInCurrentFolder.hq_cs");
#perl Scripts/TdrMappingScripts.pl hg19_mature_and_pre.fa trimmed_small_RNA-seq.fastq

if ( -l $fastqFileInCurrentFolder or -e $fastqFileInCurrentFolder ) {
    unlink $fastqFileInCurrentFolder
        or die "Failed to remove file $fastqFileInCurrentFolder: $!\n";
}
