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
  for my $comparison_name ( keys %{$comparisons} ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$comparison_name";
    my $filtered;
    my @group_names = @{ $comparisons->{$comparison_name}; };

    my $dmcpgsFile1=${comparison_name}."_".$group_names[0].".dmcpgs";
    my $dmcpgsFile2=${comparison_name}."_".$group_names[1].".dmcpgs";

    @result_files = ();
    push( @result_files, "$cur_dir/${dmcpgsFile1}" );
    $filtered = filter_array( \@result_files, $pattern, $removeEmpty );
    if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
      $result->{$comparison_name . "_" . $group_names[0]} = $filtered;
    }
    
    @result_files = ();
    push( @result_files, "$cur_dir/${dmcpgsFile2}" );
    $filtered = filter_array( \@result_files, $pattern, $removeEmpty );
    if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
      $result->{$comparison_name . "_" . $group_names[1]} = $filtered;
    }
  }

  return $result;
}

sub get_absolute_final_file {
  my ( $self, $config, $section, $comparison_name ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $comparisons = get_raw_files( $config, $section );

  my $cur_dir      = $result_dir . "/$comparison_name";

  my @group_names = @{ $comparisons->{$comparison_name} };

  my $dmcpgsFile1=${comparison_name}."_".$group_names[0].".dmcpgs";

  return( "$cur_dir/${dmcpgsFile1}" );
}

sub get_result_pbs {
  my ( $self, $config, $section ) = @_;
  
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  
  my $comparisons = get_raw_files( $config, $section );

  my $result = {};
  
  for my $comparison_name ( sort keys %$comparisons ) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $comparison_name );

    my @group_names = @{ $comparisons->{$comparison_name} };

    $result->{$comparison_name . "_" . $group_names[0]} = $pbs_file;
    $result->{$comparison_name . "_" . $group_names[1]} = $pbs_file;
  }

  #print(Dumper($result));
  
  return ($result);
}

1;
