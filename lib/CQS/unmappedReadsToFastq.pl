use strict;
use warnings;

#use XML::Simple;
use File::Basename;

#my $resultFolder="/scratch/cqs/zhaos/vickers/20150728_3018-KCV-39-40-3220-KCV-1/";
#my $sampleName="Ctr_261";

my $outFile = $ARGV[0];

#my $resultFolder=$ARGV[0];
#my $sampleName=$ARGV[1];

my $identicalFastqFile    = $ARGV[1];
my $smallRNAreadsFile     = $ARGV[2];
my $perfectmatchReadsFile = $ARGV[3];

my %readsDel;

#match to small RNA reads
#my $smallRNAreadsFile=$resultFolder.'/bowtie1_genome_1mm_NTA_smallRNA_count/result/'.$sampleName."/$sampleName.bam.count.mapped.xml";
foreach my $smallRNAreadsFileEach ( split( ",", $smallRNAreadsFile ) ) {
  open XML, "grep 'query name=' $smallRNAreadsFileEach|" or die "Can't read $smallRNAreadsFileEach\n";
  print "Reading $smallRNAreadsFileEach\n";
  while (<XML>) {
    my $temp = ( split '"', $_ )[1];
    $temp =~ s/:CLIP_(\w+)?$//;
    $temp = '@' . $temp;
    $readsDel{$temp} = '';
  }

  #	my $smallRNAreadsFileContent = XMLin($smallRNAreadsFileEach);
  #	my @mappedReads              = keys %{ ${$smallRNAreadsFileContent}{'queries'}{'query'} };
  #	foreach my $temp (@mappedReads) {
  #		$temp =~ s/:CLIP_(\w+)?$//;
  #		$temp = '@' . $temp;
  #		$readsDel{$temp} = '';
  #	}
}

if ( defined $perfectmatchReadsFile ) {
  foreach my $perfectmatchReadsFileEach ( split( ",", $perfectmatchReadsFile ) ) {

    #perfect match reads
    #my $perfectmatchReadsFile=$resultFolder.'/bowtie1_genome_1mm_NTA_pmnames/result/'.$sampleName.'.pmnames';
    open READ, "<$perfectmatchReadsFileEach" or die $!;
    print "Reading $perfectmatchReadsFileEach\n";
    while (<READ>) {
      chomp;
      if (/_$/) {    #
        s/:CLIP_$//;
        $readsDel{ '@' . $_ } = '';

        #		push @mappedReads, $_;
      }
    }
  }
}
my $readsDelCount = scalar( keys %readsDel );
print "$readsDelCount reads labeled as mapped\n";

#go through identical folder to make new fastq files
#my $identicalFastqFile=$resultFolder.'/identical/result/'.$sampleName.'_clipped_identical.fastq.gz';
if ( $identicalFastqFile =~ /\.gz$/ ) {
	open( FASTQ, "zcat $identicalFastqFile|" ) or die $!;
	my $identicalFastqCountFile=$identicalFastqFile;
	$identicalFastqCountFile=~s/\.gz$//;
	open (DUPCOUNT, $identicalFastqCountFile.'.dupcount') or die $!;
	
} else {
	open( FASTQ, $identicalFastqFile ) or die $!;
	open (DUPCOUNT, $identicalFastqFile.'.dupcount') or die $!;
}

my %fastq2Count;
while(<DUPCOUNT>) {
	my @lines=( split '\t', $_ );
	$fastq2Count{$lines[0]}=$_;
}
my $dupCount = scalar( keys %fastq2Count );
print "$dupCount reads recorded in DupCount file\n";

#my $identicalFastqFileBase=basename($identicalFastqFile);
open RESULT, "| gzip -c > $outFile" or die "error writing result: $!";
my $outDupcountFile=$outFile;
$outDupcountFile=~s/\.gz$/.dupcount/;
#$outDupcountFile=$outDupcountFile.".dupcount";
open RESULTCOUNT, ">$outDupcountFile" or die "error writing result: $!";
print RESULTCOUNT "Query\tCount\tSequence\n";
my $outMappedDupcountFile=$outFile;
$outMappedDupcountFile=~s/\.unmapped\.fastq\.gz$/.mapped.dupcount/;
open MAPPEDCOUNT, ">$outMappedDupcountFile" or die "error writing result: $!";
print MAPPEDCOUNT "Query\tCount\tSequence\n";

while ( my $line1 = <FASTQ> ) {
  my $line2   = <FASTQ>;
  my $line3   = <FASTQ>;
  my $line4   = <FASTQ>;
  my $readKey = ( split( " ", $line1 ) )[0];
  if ( exists $readsDel{$readKey} ) {
    delete $readsDel{$readKey};
    $readKey=~s/^@//;
    print MAPPEDCOUNT $fastq2Count{$readKey};
  }
  else {
    print RESULT $line1 . $line2 . $line3 . $line4;
    $readKey=~s/^@//;
    print RESULTCOUNT $fastq2Count{$readKey};
  }
}
close RESULT;
close RESULTCOUNT;
close MAPPEDCOUNT;

if (%readsDel) {
  $readsDelCount = scalar( keys %readsDel );
  print "WARNING: $readsDelCount reads labeled as mapped but not found in fastq\n";
  my $temp = ( keys %readsDel )[0];
  print( "First mapped reads but not in fastq: " . $temp );
}

print("Success!\n");

