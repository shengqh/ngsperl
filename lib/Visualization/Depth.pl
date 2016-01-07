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

Depth -b bedFile -c configFile

Options:

  -b|--bedFile            Input bed file, the forth column should be result file prefix
  -c|--configFile         Input file list with two columns: sample name and bam files
  -h|--help               This page.
";

Getopt::Long::Configure('bundling');

my $help;
my $bedFile;
my $configFile;

GetOptions(
  'h|help'         => \$help,
  'b|bedFile=s'    => \$bedFile,
  'c|configFile=s' => \$configFile,
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

my $bamFilesStr = join( ' ', @bamFiles );

open( BED, $bedFile ) or die "Cannot open file $bedFile";
while (<BED>) {
  s/\r|\n//g;
  my ( $chr, $start, $end, $fileprefix ) = split "\t";
  if ( defined $start && defined $end && defined $fileprefix ) {
    print "samtools depth -r ${chr}:${start}-${end} $bamFilesStr | cut -f${cutindecies} > ${fileprefix}.depth \n";
  }
}
close BED;

