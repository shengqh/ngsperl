#!/usr/bin/perl
package Annotation::CBioPortal;

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
  $self->{_suffix} = "_cbp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $sampleNamePattern = get_option( $config, $section, "sample_name_pattern" );
  my $sampleNameSuffix  = get_option( $config, $section, "sample_name_suffix", "" );
  my $geneNames         = get_option( $config, $section, "gene_names" );

  my $script = dirname(__FILE__) . "/CBioPortal.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $log = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  for my $sample_name ( sort keys %raw_files ) {
    my @inputFiles = @{ $raw_files{$sample_name} };
    for my $inputFile (@inputFiles) {
      my $filename  = basename($inputFile);
      my $prefix    = change_extension( $filename, "${sampleNameSuffix}.cBioPortal" );
      my $finalFile = $prefix . ".oncoprinter.txt";
      print $pbs "if [ ! -e $finalFile ]; then 
  R --vanilla -f $script --args $inputFile $prefix $sampleNamePattern $geneNames
fi
";
    }
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $sampleNameSuffix = get_option( $config, $section, "sample_name_suffix" );

  my $result = {};

  for my $sample_name ( sort keys %raw_files ) {
    my @result_files = ();
    my @inputFiles   = @{ $raw_files{$sample_name} };
    for my $inputFile (@inputFiles) {
      my $filename = basename($inputFile);
      my $prefix = change_extension( $filename, "${sampleNameSuffix}.cBioPortal" );
      push( @result_files, "$result_dir/${prefix}.oncoprinter.txt" );
      push( @result_files, "$result_dir/${prefix}.mutationmapper.txt" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
