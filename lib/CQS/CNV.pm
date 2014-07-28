#!/usr/bin/perl
package CQS::CNV;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::DNASeq;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::ClassFactory;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(cnvnator conifer cnmops freec)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub cnvnator {
  my ( $config, $section ) = @_;
  my $obj = instantiate("CNV::CNVnator");
  $obj->perform( $config, $section );
}

sub conifer {
  my ( $config, $section ) = @_;
  my $obj = instantiate("CNV::Conifer");
  $obj->perform( $config, $section );
}

sub cnmops {
  my ( $config, $section ) = @_;
  my $obj = instantiate("CNV::cnMops");
  $obj->perform( $config, $section );
}

sub freec {
  my ( $config, $section ) = @_;
  my $obj = instantiate("CNV::Freec");
  $obj->perform( $config, $section );
}

1;
