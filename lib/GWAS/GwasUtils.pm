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
    qw(readChromosomesFromBedFile)
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


1;
