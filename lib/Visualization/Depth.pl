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

Depth -b bedFile -c configFile [-s]

Options:

  -b|--bedFile            Input bed file, the forth column should be result file prefix
  -c|--configFile         Input file list with two columns: sample name and bam files
  -s|--singlePdf          Output as single pdf (for small dataset)
  -h|--help               This page.
";

Getopt::Long::Configure('bundling');

my $help;
my $bedFile;
my $configFile;
my $singlePdf;

GetOptions(
  'h|help'         => \$help,
  'b|bedFile=s'    => \$bedFile,
  'c|configFile=s' => \$configFile,
  's|siglepdf'     => \$singlePdf,
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
open( CON, $configFile ) or die "Cannot open file $configFile";
while (<CON>) {
  s/\r|\n//g;
  my ( $name, $file ) = split "\t";
  if ( defined $name && defined $file ) {
    push( @bamNames, $name );
    push( @bamFiles, $file );
    $curcutindex += 3;
    $cutindecies = $cutindecies . "," . $curcutindex;
  }
}
close CON;

my $bamNamesStr = join( '\\t', @bamNames );
my $bamFilesStr = join( ' ',   @bamFiles );

my $r = dirname(__FILE__) . "/Depth.r";

my $dataFile = basename($bedFile) . ".depth";
open( BED, $bedFile ) or die "Cannot open file $bedFile";
`printf "Chr\tPosition\t${bamNamesStr}\tFile\n" > $dataFile`;

my $keys = {};
while (<BED>) {
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
    my $cmd = "samtools depth -r ${chr}:${start}-${end} $bamFilesStr | sed -e \"s/\$/\t$fileprefix/g \" >> $dataFile";
    system($cmd);

    $keys->{$fileprefix} = 1;
  }
}
close BED;

my $singlePdfStr = ( defined $singlePdf ) ? 1 : 0;
my $outputFile = ( defined $singlePdf ) ? "${dataFile}.pdf" : basename($bedFile) . ".";
system("R --vanilla -f $r --args $singlePdfStr $dataFile $outputFile");

