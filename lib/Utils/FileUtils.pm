#!/usr/bin/env perl
use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(changeExtension)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

sub changeExtension {
  my ( $result, $newExtension ) = @_;
  $result =~ s/\.[^.]+$//;
  return ( $result . $newExtension );
}

1;
