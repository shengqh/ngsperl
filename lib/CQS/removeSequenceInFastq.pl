use strict;
use warnings;

use File::Basename;

my $finalFile=$ARGV[0];
my $outFile=$finalFile . ".tmp";
my $sequenceToDel=$ARGV[1];
my $fastqFile=$ARGV[2];

#my $fastqFile = "/scratch/cqs/zhaos/temp/test.fastq";
#my $sequenceToDel  = "aaaaaaaaataaaaaccccccc;gggtgggtgggggggg";

my @sequencesToDel=( split /;/, $sequenceToDel );

if ( ! -e $fastqFile) {
    die "the file does not existi: " . $fastqFile;
}

if ( $fastqFile =~ /\.gz$/ ) {
	open( FASTQ, "zcat $fastqFile|" ) or die $!;
}
else {
	open( FASTQ, $fastqFile ) or die $!;
}
#my $fastqFileBase=basename($fastqFile);
open RESULT, "| gzip -c > $outFile" or die "error writing result: $!";

my $delCount=0;

while ( my $line1 = <FASTQ> ) {
	my $line2 = <FASTQ>;
	chomp $line2;
	my $line3 = <FASTQ>;
	my $line4 = <FASTQ>;
	my $keepSign=1;
	foreach my $sequence (@sequencesToDel) {
		if (length($line2)>=length($sequence)) {
			if ($line2=~/$sequence/) {
				$keepSign=0;
				$delCount++;
				last;
			}
		} else { #length($line2)<length($sequence)
			if ($sequence=~/$line2/) {
				$keepSign=0;
				$delCount++;
				last;
			}
		}
	}
	if ($keepSign) {
			print RESULT $line1.$line2."\n".$line3.$line4
	}
}

if (-e $finalFile){
  unlink($finalFile);
}
rename($outFile, $finalFile);

print "Success: $delCount reads were deleted\n";
