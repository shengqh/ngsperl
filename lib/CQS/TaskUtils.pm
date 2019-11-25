package CQS::TaskUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( get_key_name )] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub get_key_name {
  my ($sample_name, $scatter_name) = @_;
  return ($sample_name . "." . $scatter_name);
}

1;
