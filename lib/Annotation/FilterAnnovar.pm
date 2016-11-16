#!/usr/bin/perl
package Annotation::FilterAnnovar;

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
  $self->{_suffix} = "_fa";
  bless $self, $class;
  return $self;
}

sub get_maximum_exac_values {
  my ($config, $section ) = @_;
  my $exac_values = get_option($config, $section, "maximum_exac_values", "0.01,0.001");
  my @result = split(',', $exac_values);
  return @result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my @exac_values = get_maximum_exac_values($config, $section);
  my $sampleNamePattern = get_option( $config, $section, "sample_name_pattern", "" );
  my $sampleNameSuffix = "";
  if($sampleNamePattern ne ""){
    $sampleNameSuffix = get_option( $config, $section, "sample_name_suffix" );
    $sampleNamePattern = "-r $sampleNamePattern"
  }

  my $script = dirname(__FILE__) . "/filterAnnovar.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }
  
  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $log = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );
  
  for my $sample_name ( sort keys %raw_files ) {
    my @annovar_files = @{ $raw_files{$sample_name} };
    my $annovar_file  = $annovar_files[0];
    
    if(scalar(@exac_values) > 0){
      for my $exac_value (@exac_values){
        my $finalFilePrefix = "${sample_name}${sampleNameSuffix}.exac${exac_value}";
        my $finalFile = $finalFilePrefix . ".gene.tsv";
        print $pbs "if [ ! -e $finalFile ]; then 
  python $script $option -i $annovar_file -e $exac_value -o $finalFilePrefix $sampleNamePattern
fi
";
      }
    }else{
        my $finalFilePrefix = "${sample_name}${sampleNameSuffix}";
        my $finalFile = $finalFilePrefix . ".gene.tsv";
        print $pbs "if [ ! -e $finalFile ]; then 
  python $script $option -i $annovar_file -o $finalFilePrefix $sampleNamePattern
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
  my @exac_values = get_maximum_exac_values($config, $section);
  my $sampleNamePattern = get_option( $config, $section, "sample_name_pattern", "" );
  my $sampleNameSuffix = "";
  if($sampleNamePattern ne ""){
    $sampleNameSuffix = get_option( $config, $section, "sample_name_suffix" );
  }

  my $result = {};
  
  for my $sample_name ( sort keys %raw_files ) {
    my @result_files = ();
    if(scalar(@exac_values) > 0){
      for my $exac_value (@exac_values){
        push(@result_files, "$result_dir/${task_name}${sampleNameSuffix}.exac${exac_value}.snv.tsv");
        push(@result_files, "$result_dir/${task_name}${sampleNameSuffix}.exac${exac_value}.gene.tsv");
      }
    }else{
      push(@result_files, "$result_dir/${task_name}${sampleNameSuffix}.snv.tsv");
      push(@result_files, "$result_dir/${task_name}${sampleNameSuffix}.gene.tsv");
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
