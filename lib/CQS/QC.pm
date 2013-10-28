#!/usr/bin/perl
package CQS::QC;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::ConfigUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(fastqc_by_pbs RNASeQC)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub fastqc_by_pbs {
  my ( $config, $section ) = @_;
  my $obj = instantiate("FastQC");
  $obj->perform( $config, $section );
}

sub RNASeQC {
  my ( $config, $section ) = @_;
  my $obj = instantiate("RNASeQC");
  $obj->perform( $config, $section );
}

1;
