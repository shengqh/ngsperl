#!/usr/bin/perl
package CQS::Task;

use strict;
use warnings;
use CQS::ConfigUtils;

sub new {
  my ($class) = @_;
  my $self = { _name => undef, _suffix => "" };
  bless $self, $class;
  return $self;
}

sub name {
  my ($self) = @_;
  return $self->{_name};
}

sub perform {
}

sub result {
  my $result = {};
  return $result;
}

sub require {
  my $result = [];
  return $result;
}

sub pbsname {
  my ( $self, $sampleName ) = @_;
  return $sampleName . $self->{_suffix} . ".pbs";
}

sub logname {
  my ( $self, $dir, $sampleName ) = @_;
  return $dir . "/" . $sampleName . $self->{_suffix} . ".log";
}

sub taskfile {
  my ( $self, $dir, $task_name ) = @_;
  my $shfile = "${dir}/${task_name}" . $self->{_suffix} . ".sh";
  return $shfile;
}

sub pbsfiles {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %fqFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( sort keys %fqFiles ) {
    $result->{$sampleName} = $pbsDir . "/" . $self->pbsname($sampleName);
  }

  return $result;
}

sub is_individual {
  return 1;
}

1;
