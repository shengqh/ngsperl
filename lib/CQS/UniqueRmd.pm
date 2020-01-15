#!/usr/bin/perl
package CQS::UniqueRmd;

use strict;
use warnings;
use File::Basename;
use File::Copy;
use Data::Dumper;
use String::Util ':all';
use CQS::UniqueWrapper;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use CQS::RmdUtils;

our @ISA = qw(CQS::UniqueWrapper);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ur";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $final_pbs = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $final_log = $self->get_log_filename( $log_dir, $task_name );
  my $final_log_desp = $cluster->get_log_description($final_log);

  my $output_file_ext = get_option( $config, $section, "output_file_ext", ".html" );

  my $rfilename = ${task_name}. $self->{_suffix} . ".Rmd";
  my $rfile     = build_rmd_file($config, $section, $result_dir, $rfilename);

  my $removeEmpty = get_option( $config, $section, "remove_empty_parameter", 0 );
  my $parametersample_files1 = writeParameterSampleFile( $config, $section, $result_dir, 1, $removeEmpty );
  my $parametersample_files2 = writeParameterSampleFile( $config, $section, $result_dir, 2, $removeEmpty );
  my $parametersample_files3 = writeParameterSampleFile( $config, $section, $result_dir, 3, $removeEmpty );
  my $parametersample_files4 = writeParameterSampleFile( $config, $section, $result_dir, 4, $removeEmpty );

  
  my $final_file = ${task_name}. $self->{_suffix} . ".html";
  my $final = $self->open_pbs( $final_pbs, $pbs_desc, $final_log_desp, $path_file, $result_dir, $final_file );
  print $final "
R -e \"library(knitr);rmarkdown::render('$rfilename');\"
";

  $self->close_pbs( $final, $final_pbs );
}

1;
