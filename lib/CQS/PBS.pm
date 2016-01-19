#!/usr/bin/perl
package CQS::PBS;

use strict;
use warnings;
use File::Basename;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(init_dir)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub init_dir {
  my ( $rootDir, $create ) = @_;

  if ( !defined $create ) {
    $create = 1;
  }

  #defined several folders
  my $pbs_dir    = "$rootDir/pbs";
  my $result_dir = "$rootDir/result";
  my $log_dir    = "$rootDir/log";

  if ($create) {
    create_directory_or_die($rootDir);
    create_directory_or_die($pbs_dir);
    create_directory_or_die($result_dir);
    create_directory_or_die($log_dir);
  }

  return ( $log_dir, $pbs_dir, $result_dir );
}

1;

