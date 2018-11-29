#!/usr/bin/perl
package QC::FastQCSummary;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_fqs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name   = basename($pbs_file);
  my $log        = $self->get_log_filename( $log_dir, $task_name );
  my $cqstools   = get_param_file( $config->{$section}{cqstools}, "cqstools", 1 );
  my $fastqc_dir = get_directory( $config, $section, "fastqc_dir", 0 );
  if ( !defined $fastqc_dir ) {
    $fastqc_dir = $result_dir;
  }
  my $log_desc = $cluster->get_log_description($log);

  my $final_file = "${task_name}.FastQC.summary.tsv";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  print $pbs "
qcimg2pdf.sh -o $task_name
mono $cqstools fastqc_summary -i $fastqc_dir -o $final_file 
";
  $self->close_pbs( $pbs, $pbs_file );

  if ( is_linux() ) {
    chmod 0755, $pbs_file;
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  push( @result_files, "${result_dir}/${task_name}.FastQC.summary.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.summary.reads.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.summary.overrepresented.tsv" );
#  push( @result_files, "${result_dir}/${task_name}.FastQC.pdf" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
