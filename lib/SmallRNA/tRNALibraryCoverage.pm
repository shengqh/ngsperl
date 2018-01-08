#!/usr/bin/perl
package SmallRNA::tRNALibraryCoverage;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_nlc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $py_script = dirname(__FILE__) . "/tRNALibraryCoverage.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }

  my $fastaFile      = get_option_file( $config, $section, "fasta_file" );
  my $speciesMapFile = get_option_file( $config, $section, "species_map_file" );
  my $species = get_option( $config, $section, "species" );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $final_file = $task_name . ".position";
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];
    my $final_file   = $sample_name . ".position";

    print $pbs "python $py_script -i $sampleFile -o $final_file $option --fasta $fastaFile --speciesMap $speciesMapFile --species $species \n";
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $final_file = $task_name . ".position";

  my $result = {};
  $result->{$task_name} = filter_array( [ $result_dir . "/" . $final_file ], $pattern );
  return $result;
}

1;
