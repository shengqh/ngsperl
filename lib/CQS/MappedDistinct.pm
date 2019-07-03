#!/usr/bin/perl
package CQS::MappedDistinct;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_dt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %firstFiles = %{ get_raw_files( $config, $section ) };
  my %secondFiles = %{ get_raw_files( $config, $section, "second" ) };
  my $firstSuffix  = get_option( $config, $section, "first_suffix" );
  my $secondSuffix = get_option( $config, $section, "second_suffix" );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $log = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  for my $sample_name ( sort keys %firstFiles ) {
    my $firstFile    = $firstFiles{$sample_name}->[0];
    my $secondFile   = $secondFiles{$sample_name}->[0];
    my $firstoutput  = $firstSuffix . $sample_name . ".distinct.count";
    my $secondoutput = $secondSuffix . $sample_name . ".distinct.count";

    print $pbs "cqstools mapped_distinct $option --inputfile1 $firstFile --outputfile1 $firstoutput --inputfile2 $secondFile --outputfile2 $secondoutput
";
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $firstSuffix  = $config->{$section}{first_suffix}  or die "define ${section}::first_suffix first";
  my $secondSuffix = $config->{$section}{second_suffix} or die "define ${section}::second_suffix first";

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $firstoutput  = $result_dir . "/" . $firstSuffix . $sample_name . ".distinct.count";
    my $secondoutput = $result_dir . "/" . $secondSuffix . $sample_name . ".distinct.count";

    my @result_files = ();

    push( @result_files, $firstoutput );
    push( @result_files, $secondoutput );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
