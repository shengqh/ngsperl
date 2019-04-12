#!/usr/bin/perl
package GWAS::GwasUtils;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(readChromosomesFromBedFile
    readGwasDataFile)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub readChromosomesFromBedFile {
  my ( $bedFile ) = @_;
  my @chroms = ();
  open( my $fin, '<', $bedFile ) or die "Could not open $bedFile\n";
  while ( my $line = <$fin> ) {
    chomp $line;
    my @parts = split( "\t", $line );
    my $chrom = $parts[0];
    if ( not grep( /^$chrom$/, @chroms ) ) {
      push @chroms, $chrom;
    }
  }
  return @chroms;
}

sub readGwasDataFile {
  my ( $config, $section, $key, $rawFiles ) = @_;
  my $result = get_raw_files( $config, $section, $key );
  my @resultKeys = (sort keys %$result);
  
  for my $sample (keys %$rawFiles){
    if (not defined $result->{$sample}){
      my ( $chrKey ) = $sample =~ /_(chr\S+)$/;
      if (not defined $result->{$chrKey}){
        die "$key in $section doesn't have key $sample or $chrKey";
      }
      $result->{$sample} = $result->{$chrKey}
    }
  }
  
  return $result;
}


1;
