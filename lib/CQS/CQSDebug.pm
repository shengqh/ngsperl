#!/usr/bin/perl
package CQS::CQSDebug;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(is_debug)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

my $DEBUG = 0;

sub is_debug() {
  return $DEBUG;
}

1;