#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

sub run_command {
  my $command = shift;
  print "$command \n";
  `$command `;
}

my $usage = "

Synopsis:

Depth -b bedFile -c configFile -r totalReadsFile -v cnvrFile [-s]

Options:

  -b|--bedFile            Input bed file, the forth column should be result file prefix
  -c|--configFile         Input file list with two columns: sample name and bam file
  -v|--cnvrFile           Input file from cn.mops module which includes loci and CNx categories cross samples
  -s|--singlePdf          Output as single pdf (for small dataset)
  -h|--help               This page.
";

Getopt::Long::Configure('bundling');

my $help;
my $bedFile;
my $configFile;
my $cnvrFile;
my $singlePdf;

GetOptions(
  'h|help'             => \$help,
  'b|bedFile=s'        => \$bedFile,
  'c|configFile=s'     => \$configFile,
  'v|cnvrFile=s'       => \$cnvrFile,
  's|siglepdf'         => \$singlePdf,
);

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined($bedFile) ) {
  print $usage;
  exit(1);
}

if ( !defined($configFile) ) {
  print $usage;
  exit(1);
}

my @bamNames    = ();
my @bamFiles    = ();
my $cutindecies = "1,2";
my $curcutindex = 1;
open( my $config, $configFile ) or die "Cannot open file $configFile";
while (<$config>) {
  s/\r|\n//g;
  my ( $name, $file ) = split "\t";
  if ( defined $name && defined $file ) {
    push( @bamNames, $name );
    push( @bamFiles, $file );
    $curcutindex += 3;
    $cutindecies = $cutindecies . "," . $curcutindex;
  }
}
close $config;

my $bamNamesStr = join( '\\t', @bamNames );
my $bamFilesStr = join( ' ',   @bamFiles );

my $r = dirname(__FILE__) . "/Depth.r";

my $readsFile = basename($bedFile) . ".reads";
open( my $reads, "> $readsFile" ) or die "Cannot open file $readsFile";
print $reads "Sample\tMappedReads\n";
for ( my $index = 0 ; $index < scalar(@bamNames) ; $index++ ) {
  my $curName        = $bamNames[$index];
  my $curBamFile     = $bamFiles[$index];
  my $curBamStatFile = $curBamFile . ".stat";
  if ( !-e $curBamStatFile ) {
    `samtools flagstat $curBamFile > $curBamStatFile`;
  }

  my $curMappedReads = `grep " mapped (" $curBamStatFile | cut -d ' ' -f1`;
  chomp($curMappedReads);
  print $reads $curName, "\t", $curMappedReads, "\n";
}
close($reads);

my $depthFile = basename($bedFile) . ".depth";
if ( !-e $depthFile ) {
  open( my $bed, $bedFile ) or die "Cannot open file $bedFile";
  `printf "Chr\tPosition\t${bamNamesStr}\tFile\n" > $depthFile`;

  my $keys = {};
  while (<$bed>) {
    s/\r|\n//g;
    my ( $chr, $start, $end, $fileprefix ) = split "\t";
    if ( defined $start && defined $end ) {
      if ( !defined $fileprefix ) {
        $fileprefix = "${chr}_${start}_${end}";
      }

      if ( exists $keys->{$fileprefix} ) {
        next;
      }

      print( $fileprefix . "\n" );
      my $cmd        = "samtools depth -d 1000000 -r ${chr}:${start}-${end} $bamFilesStr | sed -e \"s/\$/\t$fileprefix/g \" >> $depthFile";
      my $returnCode = system($cmd);
      if ( $returnCode != 0 ) {
        die("Error return code = $returnCode, exit.");
      }

      $keys->{$fileprefix} = 1;
    }
  }
  close $bed;
}

my $singlePdfStr = ( defined $singlePdf ) ? 1 : 0;
my $outputFile = ( defined $singlePdf ) ? "${depthFile}.pdf" : "";

if(-e "Depth.r"){
  unlink("Depth.r");
}

`echo "readFile<-\\"$readsFile\\"" >> Depth.r`;

if ( defined($cnvrFile) ) {
  `echo "cnvrFile<-\\"$cnvrFile\\"" >> Depth.r`;
}

`cat $r >> Depth.r`;
system("R --vanilla -f Depth.r --args $singlePdfStr $depthFile $bedFile $outputFile");

