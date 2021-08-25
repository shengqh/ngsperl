#!/usr/bin/perl
package CQS::UniqueTask;

use strict;
use warnings;
use CQS::Task;
use CQS::ConfigUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_pbskey} = "";
  $self->{_depend_all} = 1;
  $self->{_forbid_tmp_folder} = 1;
  bless $self, $class;
  return $self;
}

sub get_clear_map {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $comparison_result = $self->SUPER::get_clear_map( $config, $section, $pattern );
  my $result            = {};
  my @result_files      = ();
  for my $crs ( values %$comparison_result ) {
    for my $cr (@$crs) {
      push( @result_files, $cr );
    }
  }
  $result->{$task_name} = \@result_files;

  return $result;
}

sub get_pbs_source {
  my ( $self, $config, $section ) = @_;

  my $task_section = $config->{$section};
  my $sources = {};
  for my $key ( keys %$task_section ) {
    my $mapname = $key;
    if ( $mapname =~ /_ref$/ ) {
      $mapname =~ s/_config_ref//g;
      $mapname =~ s/_ref//g;
      my $refpbsmap = get_ref_section_pbs( $config, $section, $mapname );
      for my $refkey (keys %$refpbsmap){
        $sources->{$refkey} = 1;
      }
    }
  }
  my $sourceNames = [sort keys %$sources];

  my $pbsFiles = $self->get_pbs_files( $config, $section );
  my $result   = {};
  for my $resKey ( keys %$pbsFiles ) {
    $result->{ $pbsFiles->{$resKey} } = [$resKey, @$sourceNames];
  }
  return $result;
}

1;
