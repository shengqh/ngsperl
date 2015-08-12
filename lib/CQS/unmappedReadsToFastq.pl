use strict;
use warnings;

use XML::Simple;
use File::Basename;

#my $resultFolder="/scratch/cqs/zhaos/vickers/20150728_3018-KCV-39-40-3220-KCV-1/";
#my $sampleName="Ctr_261";

my $outFile=$ARGV[0];
#my $resultFolder=$ARGV[0];
#my $sampleName=$ARGV[1];

my $smallRNAreadsFile=$ARGV[1];
my $perfectmatchReadsFile=$ARGV[2];
my $identicalFastqFile=$ARGV[3];

my %readsDel;

#match to small RNA reads
#my $smallRNAreadsFile=$resultFolder.'/bowtie1_genome_1mm_NTA_smallRNA_count/result/'.$sampleName."/$sampleName.bam.count.mapped.xml";
my $smallRNAreadsFileContent = XMLin($smallRNAreadsFile);
my @mappedReads=keys %{${$smallRNAreadsFileContent}{'queries'}{'query'}};
foreach my $temp (@mappedReads) {
	$temp=~s/:CLIP_(\w+)?$//;
	$temp='@'.$temp;
	$readsDel{$temp}='';
}

#perfect match reads
#my $perfectmatchReadsFile=$resultFolder.'/bowtie1_genome_1mm_NTA_pmnames/result/'.$sampleName.'.pmnames';
open READ, "<$perfectmatchReadsFile" or die $!;
while (<READ>) {
	chomp;
	if (/_$/) { #
		s/:CLIP_$//;
		$readsDel{'@'.$_}='';
#		push @mappedReads, $_;
	}
}
my $readsDelCount= scalar (keys %readsDel);
print"$readsDelCount reads labeled as maped\n";

#go through identical folder to make new fastq files
#my $identicalFastqFile=$resultFolder.'/identical/result/'.$sampleName.'_clipped_identical.fastq.gz';
if ( $identicalFastqFile =~ /\.gz$/ ) {
	open( FASTQ, "zcat $identicalFastqFile|" ) or die $!;
}
else {
	open( FASTQ, $identicalFastqFile ) or die $!;
}
#my $identicalFastqFileBase=basename($identicalFastqFile);
open RESULT, "| gzip -c > $outFile" or die "error writing result: $!";

while ( my $line1 = <FASTQ> ) {
	my $line2 = <FASTQ>;
	my $line3 = <FASTQ>;
	my $line4 = <FASTQ>;
	my $readKey=(split(" ",$line1))[0];
	if (exists $readsDel{$readKey}) {
		delete $readsDel{$readKey};
	} else {
		print RESULT $line1.$line2.$line3.$line4;
	}
}
close RESULT;

if (%readsDel) {
	$readsDelCount= scalar (keys %readsDel);
	print"WARNING: $readsDelCount reads labeled as unmapped but not found in fastq\n";
	my $temp=(keys %readsDel )[0];
	print ("First unmapped reads but not in fastq: ".$temp);
}

print ("Success!\n");


