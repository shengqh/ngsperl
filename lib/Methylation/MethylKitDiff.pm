#!/usr/bin/perl
package Methylation::MethylKitDiff;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::GroupTask;
use CQS::IndividualR;
use CQS::NGSCommon;
use CQS::StringUtils;
use Data::Dumper;

our @ISA = qw(CQS::IndividualR);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "";
  #$self->{_group_keys} = ["source"];
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern, $removeEmpty ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $comparisons = get_raw_files( $config, $section );

  my $result = {};
  for my $group_name ( keys %{$comparisons} ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$group_name";
    my $filtered;
    my @sampleNames = @{ $comparisons->{$group_name}; };

    my $dmcpgsFile1=${group_name}."_".$sampleNames[0].".dmcpgs";
    my $dmcpgsFile2=${group_name}."_".$sampleNames[1].".dmcpgs";

    @result_files = ();
    push( @result_files, "$cur_dir/${dmcpgsFile1}" );
    $filtered = filter_array( \@result_files, $pattern, $removeEmpty );
    if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
      $result->{$group_name . "_" . $sampleNames[0]} = $filtered;
    }
    
    @result_files = ();
    push( @result_files, "$cur_dir/${dmcpgsFile2}" );
    $filtered = filter_array( \@result_files, $pattern, $removeEmpty );
    if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
      $result->{$group_name . "_" . $sampleNames[1]} = $filtered;
    }
  }

  return $result;
}

1;
