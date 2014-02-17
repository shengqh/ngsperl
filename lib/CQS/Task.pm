#!/usr/bin/perl
package CQS::Task;

use strict;
use warnings;
use CQS::ConfigUtils;

sub new {
  my ($class) = @_;
  my $self = { _name => undef, _suffix => "", _task_prefix => "", _task_suffix => "", _pbskey => "source" };
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

sub getname {
  my ( $self, $name, $extension, $hassuffix ) = @_;
  if ( !defined $extension ) {
    $extension = "";
  }
  if ( !defined $hassuffix ) {
    $hassuffix = 1;
  }

  if ($hassuffix) {
    return $self->{_task_prefix} . $name . $self->{_task_suffix} . $self->{_suffix} . $extension;
  }
  else {
    return $self->{_task_prefix} . $name . $self->{_task_suffix} . $extension;
  }
}

sub getfile {
  my ( $self, $dir, $name, $extension, $hassuffix ) = @_;
  if ( !defined $extension ) {
    $extension = "";
  }
  if ( !defined $hassuffix ) {
    $hassuffix = 1;
  }

  return $dir . "/" . $self->getname( $name, $extension, $hassuffix );
}

sub pbsname {
  my ( $self, $sampleName ) = @_;
  return $self->getname( $sampleName, ".pbs" );
}

sub pbsfile {
  my ( $self, $dir, $sampleName ) = @_;
  return $self->getfile( $dir, $sampleName, ".pbs" );
}

sub logfile {
  my ( $self, $dir, $sampleName ) = @_;
  return $self->getfile( $dir, $sampleName, ".log" );
}

sub taskfile {
  my ( $self, $dir, $task_name ) = @_;
  return $self->getfile( $dir, $task_name, ".sh" );
}

sub pbsfiles {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};
  if ( $self->{_pbskey} eq "" ) {
    $result->{$task_name} = $self->pbsfile( $pbsDir, $task_name );
  }
  else {
    my %fqFiles = %{ get_raw_files( $config, $section, $self->{_pbskey} ) };

    for my $sampleName ( sort keys %fqFiles ) {
      $result->{$sampleName} = $self->pbsfile( $pbsDir, $sampleName );
    }
  }

  return $result;
}

1;
