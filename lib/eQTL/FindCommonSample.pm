#!/usr/bin/perl
package eQTL::FindCommonSample;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_dp";
  bless $self, $class;
  return $self;
}

sub acceptSample {
  my ( $self, $config, $section, $sample_name ) = @_;
  my $sample_set = get_option( $config, $section, "sample_set", [] );
  if ( ( scalar(@$sample_set) > 0 ) & !grep( /^$sample_name$/, @$sample_set ) ) {
    return (0);
  }
  return (1);
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );
  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $prefix = get_option( $config, $section, "prefix", "" );

  my $rnaseq_files = get_raw_files( $config, $section );
  my $rnaseq_names = get_option( $config, $section, "rnaseq_names" );

  my $fam_files = get_raw_files( $config, $section, "plink_fam_files" );
  my $fam_names = get_option( $config, $section, "plink_fam_names" );

  my $pattern = get_option( $config, $section, "pattern", 1 );
  my $find_common_sample_script = dirname(__FILE__) . "/findCommonSample.py";
  if ( !-e $find_common_sample_script ) {
    die "File not found : " . $find_common_sample_script;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $expect_files = $self->result($config, $section);
  for my $sample_name ( sort keys %$rnaseq_files ) {
    if ( !$self->acceptSample( $config, $section, $sample_name ) ) {
      next;
    }

    my $sampleRnaseqFiles    = $rnaseq_files->{$sample_name};
    my $sampleRnaseqFilesStr = join( ',', @$sampleRnaseqFiles );
    my $sampleFamFiles       = $fam_files->{$sample_name};
    my $sampleFamFilesStr    = join( ',', @$sampleFamFiles );

    my $file_name          = $prefix . $sample_name;
    my $common_prefix      = $file_name . "_common";
    my $common_snp_file    = $common_prefix . ".snp.samples";
    my $common_rnaseq_file = $common_prefix . ".rnaseq.samples";

    my $common_bed_file = $common_prefix . ".bed";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";
    
    my $final_file = $expect_files->{$sample_name}->[0];
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $pbs "
echo findCommonSample=`date`
python3 $find_common_sample_script -r $sampleRnaseqFilesStr --rnaseqNames $rnaseq_names -f $sampleFamFilesStr --famNames $fam_names -o $common_prefix -p $pattern
";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all fastx_trimmer tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $prefix = get_option( $config, $section, "prefix", "" );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $rnaseq_files = get_raw_files( $config, $section );
  my $rnaseq_names = get_option( $config, $section, "rnaseq_names" );

  my $fam_files = get_raw_files( $config, $section, "plink_fam_files" );
  my $fam_names = get_option( $config, $section, "plink_fam_names" );

  my @rnaseqNames = split( ',', $rnaseq_names );
  my @famNames    = split( ',', $fam_names );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    if ( !$self->acceptSample( $config, $section, $sample_name ) ) {
      next;
    }

    my $file_name     = $prefix . $sample_name;
    my $common_prefix = $file_name . "_common";

    my @result_files = ();
    for my $rnaseqName (@rnaseqNames) {
      push( @result_files, "${result_dir}/${common_prefix}.rnaseq.${rnaseqName}.samples" );
    }
    for my $famName (@famNames) {
      push( @result_files, "${result_dir}/${common_prefix}.snp.${famName}.samples" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
