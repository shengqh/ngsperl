#!/usr/bin/perl
package Chipseq::MACS2Callpeak;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;
use Data::Dumper;

our @ISA = qw(CQS::AbstractMACSCallpeak);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %raw_files = %{ $self->get_current_raw_files( $config, $section, "groups" ) };

  my $has_control = has_raw_files( $config, $section, "controls" );
  print "has_control = " . $has_control . "\n";
  my %control_files = {};
  if ($has_control) {
    %control_files = %{ $self->get_current_raw_files( $config, $section, "controls" ) };
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my @sample_files = @{ $raw_files{$sample_name} };
    my $treatment = "-t " . join( " ", @sample_files );

    my $control = "";
    if ($has_control) {
      my @control_files = @{ $control_files{$sample_name} };
      $control = "-c " . join( " ", @control_files );
    }
    my $final_file = "${sample_name}_peaks.bed";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "   
macs2 callpeak $option $treatment $control -n $sample_name
cut -f1-6 ${sample_name}_peaks.narrowPeak > ${sample_name}_peaks.narrowPeak.bed
";

    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;
  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ $self->get_current_raw_files( $config, $section ) };
  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $cur_dir      = $result_dir . "/$sample_name";
    my @result_files = ();
    push( @result_files, $cur_dir . "/${sample_name}_treat_pileup.bdg" );
    push( @result_files, $cur_dir . "/${sample_name}_control_lambda.bdg" );
    push( @result_files, $cur_dir . "/${sample_name}_peaks.narrowPeak" );
    push( @result_files, $cur_dir . "/${sample_name}_peaks.narrowPeak.bed" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
