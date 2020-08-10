use strict;
use warnings;

#use XML::Simple;
use File::Basename;

my $DEBUG = 0;

my $outFile = $ARGV[0];
my $identicalFastqFile    = $ARGV[1];
my $smallRNAreadsFile     = $ARGV[2];
my $perfectmatchReadsFile = $ARGV[3];
my $shortReads            = $ARGV[4];
if ( !defined $shortReads ) {
  $shortReads = 20;
}

if($DEBUG){
  $outFile = "/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_unmapped_reads/pbs/test/test.unmapped.fastq.gz";
  $identicalFastqFile    = "/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/preprocessing/identical/result/Urine_WT_14_clipped_identical.fastq.gz";
  $smallRNAreadsFile     = "/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_count/result/Urine_WT_14/Urine_WT_14.count.mapped.xml";
  $perfectmatchReadsFile = "/scratch/cqs/shengq2/vickers/20170628_smallRNA_3018-KCV-77_78_79_mouse_v3/host_genome/bowtie1_genome_1mm_NTA_pmnames/result/Urine_WT_14.pmnames";
}

die "no identical fastq file: " . $identicalFastqFile if ! -e $identicalFastqFile;
#die "no smallRNA reads file:" . $smallRNAreadsFile if ! -e $smallRNAreadsFile;

sub getDupcountFile {
  my $result = shift;
  if ( $result =~ /\.gz$/ ) {
    $result =~ s/\.gz$/.dupcount/;
  }
  else {
    $result = $result . '.dupcount';
  }
  return $result;
}

my %readsDel;

my $totalReads    = 0;
my $featureReads  = 0;
my $genomeReads   = 0;
my $tooShortReads = 0;
my $unmappedReads = 0;

#match to small RNA reads
#my $smallRNAreadsFile=$resultFolder.'/bowtie1_genome_1mm_NTA_smallRNA_count/result/'.$sampleName."/$sampleName.bam.count.mapped.xml";
foreach my $smallRNAreadsFileEach ( split( ",", $smallRNAreadsFile ) ) {
  die "no smallRNA reads file:" . $smallRNAreadsFileEach if ! -e $smallRNAreadsFileEach;
  open XML, "grep 'query name=' $smallRNAreadsFileEach|" or die "Can't read $smallRNAreadsFileEach\n";
  print "Reading $smallRNAreadsFileEach\n";
  while (<XML>) {
    my $temp = ( split '"', $_ )[1];
    $temp =~ s/:CLIP_(\w+)?$//;
    $temp = '@' . $temp;
    $readsDel{$temp} = '';
  }
}

my %perfectmatchOnlyReads;
if ( defined $perfectmatchReadsFile ) {
  foreach my $perfectmatchReadsFileEach ( split( ",", $perfectmatchReadsFile ) ) {
    die "no perfectmatchReadsFile:" . $perfectmatchReadsFileEach if ! -e $perfectmatchReadsFileEach;
    #perfect match reads
    #my $perfectmatchReadsFile=$resultFolder.'/bowtie1_genome_1mm_NTA_pmnames/result/'.$sampleName.'.pmnames';
    open READ, "<$perfectmatchReadsFileEach" or die $!;
    print "Reading $perfectmatchReadsFileEach\n";
    while (<READ>) {
      chomp;
      if (/_$/) {    #
        s/:CLIP_$//;
        my $readKey = '@' . $_;
        if ( !exists $readsDel{$readKey} ) {    #Not in small RNA mapping result
          $readsDel{$readKey}              = '';
          $perfectmatchOnlyReads{$readKey} = '';
        }
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
}
else {
  open( FASTQ, $identicalFastqFile ) or die $!;
}

my $identicalFastqCountFile = getDupcountFile($identicalFastqFile);
open( DUPCOUNT, $identicalFastqCountFile) or die $!;

my %fastq2Count;
while (<DUPCOUNT>) {
  my @lines = ( split '\t', $_ );
  $fastq2Count{'@' . $lines[0] } = [$_, $lines[1]];
}
my $dupCount = scalar( keys %fastq2Count ) - 1;    #-1 becasue of the title in dupcount file
print "$dupCount reads recorded in DupCount file\n";

#my $identicalFastqFileBase=basename($identicalFastqFile);
open RESULT, "| gzip -c > $outFile" or die "error writing result: $!";

my $outDupcountFile = getDupcountFile($outFile);
open RESULTCOUNT, ">$outDupcountFile" or die "error writing result: $!";
print RESULTCOUNT "Query\tCount\tSequence\n";

if (%perfectmatchOnlyReads) {                      #Has perfect matched file
  my $outMappedFastqFile = $outFile;
  $outMappedFastqFile =~ s/\.unmapped\.fastq\.gz$/\.mappedToHostGenome\.fastq\.gz/;
  open RESULTGENOME, "| gzip -c > $outMappedFastqFile" or die "error writing result: $!";

  my $outMappedDupcountFile = getDupcountFile($outMappedFastqFile);
  open MAPPEDCOUNT, ">$outMappedDupcountFile" or die "error writing result: $!";
  print MAPPEDCOUNT "Query\tCount\tSequence\n";

  my $outShortFastqFile = $outFile;
  $outShortFastqFile =~ s/\.unmapped\.fastq\.gz$/\.short\.fastq\.gz/;
  open RESULTSHORT, "| gzip -c > $outShortFastqFile" or die "error writing result: $!";

  my $outShortDupcountFile = getDupcountFile($outShortFastqFile);
  open SHORTCOUNT, ">$outShortDupcountFile" or die "error writing result: $!";
  print SHORTCOUNT "Query\tCount\tSequence\n";

  while ( my $line1 = <FASTQ> ) {
    my $line2   = <FASTQ>;
    my $line3   = <FASTQ>;
    my $line4   = <FASTQ>;
    my $readKey = ( split( " ", $line1 ) )[0];
    my $count   = $fastq2Count{$readKey};
    $totalReads = $totalReads + $count->[1];

    if ( exists $readsDel{$readKey} ) {
      delete $readsDel{$readKey};

      if ( exists $perfectmatchOnlyReads{$readKey} ) {    #Reads from perfect match genoem, not small RNA
        $readKey =~ s/^@//;
        if ( ( length($line2) - 1 ) < $shortReads ) {     # -1 because of next line sign, < means short Reads>
          print RESULTSHORT $line1 . $line2 . $line3 . $line4;
          print SHORTCOUNT $count->[0];
          $tooShortReads = $tooShortReads + $count->[1];
        }
        else {                                            #long reads
          print RESULTGENOME $line1 . $line2 . $line3 . $line4;
          print MAPPEDCOUNT $count->[0];
          $genomeReads = $genomeReads + $count->[1];
        }
      }
      else {
        $featureReads = $featureReads + $count->[1];
      }
    }
    else {
      $readKey =~ s/^@//;
      if ( ( length($line2) - 1 ) < $shortReads ) {       # -1 because of next line sign, < means short Reads>
        print RESULTSHORT $line1 . $line2 . $line3 . $line4;
        print SHORTCOUNT $count->[0];
        $tooShortReads = $tooShortReads + $count->[1];
      }
      else {
        print RESULT $line1 . $line2 . $line3 . $line4;
        print RESULTCOUNT $count->[0];
        $unmappedReads = $unmappedReads + $count->[1];
      }
    }
  }
  close RESULTGENOME;
  close MAPPEDCOUNT;
  close RESULTSHORT;
  close SHORTCOUNT;
}
else {    #Don't have perfect matched file
  while ( my $line1 = <FASTQ> ) {
    my $line2   = <FASTQ>;
    my $line3   = <FASTQ>;
    my $line4   = <FASTQ>;
    my $readKey = ( split( " ", $line1 ) )[0];
    my $count   = $fastq2Count{$readKey};
    $totalReads = $totalReads + $count->[1];

    if ( exists $readsDel{$readKey} ) {
      delete $readsDel{$readKey};
      $featureReads = $featureReads + $count->[1];
    }
    else {
      print RESULT $line1 . $line2 . $line3 . $line4;
      $readKey =~ s/^@//;
      print RESULTCOUNT $count->[0];
      $unmappedReads = $unmappedReads + $count->[1];
    }
  }

}
close RESULT;
close RESULTCOUNT;

if (%readsDel) {
  $readsDelCount = scalar( keys %readsDel );
  print "WARNING: $readsDelCount reads labeled as mapped but not found in fastq\n";
  my $temp = ( keys %readsDel )[0];
  print( "First mapped reads but not in fastq: " . $temp );
}

my $infofile = $outFile . ".info";
open RESULTINFO, ">$infofile" or die "error writing result: $!";
my $mappedReads = $featureReads + $genomeReads;
print RESULTINFO "Category\tCount
TotalReads\t$totalReads
MappedReads\t$mappedReads
FeatureReads\t$featureReads
GenomeReads\t$genomeReads
TooShortReads\t$tooShortReads
UnmappedReads\t$unmappedReads
";

print("Success!\n");

