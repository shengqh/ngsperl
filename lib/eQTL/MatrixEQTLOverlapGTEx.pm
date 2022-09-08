#!/usr/bin/perl
package eQTL::MatrixEQTLOverlapGTEx;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::Task;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mg";
  bless $self, $class;
  return $self;
}

sub acceptSample {
  my ( $self, $config, $section, $sample_name ) = @_;
  my $gtex_files = get_raw_files( $config, $section );
  my $matrixEQTL_files = get_raw_files( $config, $section, "matrix_eqtl" );
  return (exists($gtex_files->{$sample_name}) & exists($matrixEQTL_files->{$sample_name}));
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $gtex_files = get_raw_files( $config, $section );
  my $matrixEQTL_files = get_raw_files( $config, $section, "matrix_eqtl" );
  my $snpbim_files = get_raw_files( $config, $section, "matrix_eqtl_bim" );
 
  my $script = dirname(__FILE__) . "/MatrixEQTLOverlapGTEx.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $expectFiles = $self->result( $config, $section );

  for my $sampleName ( keys %$gtex_files ) {
    if (!$self->acceptSample($config, $section, $sampleName)){
      next;
    }
    
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sampleName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sampleName );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_file = basename( $expectFiles->{$sampleName}->[0] );
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    my $gtex_file    = $gtex_files->{$sampleName}[0];
    my $matrixeqtl_file    = $matrixEQTL_files->{$sampleName}[0];
    my $snpbim_file    = $snpbim_files->{$sampleName}[0];

    print $pbs "python3 $script -i $matrixeqtl_file -b $snpbim_file -g $gtex_file -o $final_file\n";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section, 0 );
  my $prefix = get_option( $config, $section, "prefix", "" );

  my $snp_genotype_files = get_raw_files( $config, $section );
  my $result = {};

  for my $sampleName (keys %$snp_genotype_files) {
    if (!$self->acceptSample($config, $section, $sampleName)){
      next;
    }
    
    my @result_files = ();
    push( @result_files, "$result_dir/${prefix}${sampleName}.matrixeqtl_gtex.tsv");
    $result->{$sampleName} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
