#!/usr/bin/env perl
use strict;
use warnings;

use List::Util qw(first);

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(readDictionaryByIndex)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

sub readDictionaryByIndex {
  my ( $fileName, $keyIndex, $valueIndex, $haveHeader ) = @_;
  if ( !defined $keyIndex || !defined $valueIndex ) {
    $keyIndex   = 0;
    $valueIndex = 1;
  }

  my $maxIndex = $keyIndex > $valueIndex ? $keyIndex : $valueIndex;

  if ( !defined $haveHeader ) {
    $haveHeader = 0;
  }

  my $result = {};
  open( my $file, "<$fileName" ) or die "Cannot open file: " . $fileName;
  if ($haveHeader) {
    my $header = <$file>;
  }
  while ( my $line = (<$file>) ) {
    chomp $line;
    my @parts = split( '\t', $line );
    next unless scalar(@parts) > $maxIndex;
    $result->{ $parts[$keyIndex] } = $parts[$valueIndex];
  }
  close($file);

  return ($result);
}

sub readDictionaryByColumnName {
  my ( $fileName, $keyName, $valueName ) = @_;

  my $result = {};
  open( my $file, "<$fileName" ) or die "Cannot open file: " . $fileName;
  my $header = <$file>;
  chomp $header;
  my @headers = split( '\t', $header );
  my $keyIndex   = first { $headers[$_] eq $keyName } 0 .. $#headers;
  my $valueIndex = first { $headers[$_] eq $valueName } 0 .. $#headers;
  close($file);

  return readDictionaryByIndex( $fileName, $keyIndex, $valueIndex, 1 );
}

1;
