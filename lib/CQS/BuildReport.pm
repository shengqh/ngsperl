#!/usr/bin/perl
package CQS::BuildReport;

use strict;
use warnings;
use File::Basename;
use File::Copy;
use Data::Dumper;
use String::Util ':all';
use CQS::UniqueTask;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_br";
  bless $self, $class;
  return $self;
}


sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $final_file     = $task_name . ".html";
  my $final_pbs      = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $final_log      = $self->get_log_filename( $log_dir, $task_name );
  my $final_log_desp = $cluster->get_log_description($final_log);

  my $rtemplate = dirname(__FILE__) . "/" . $config->{$section}{report_rmd_file};
  my $rfile     = $result_dir . "/${task_name}.Rmd";
  copy( $rtemplate, $rfile ) or die "Copy failed: $!";

  my @additional_rtemplates = split(';', $config->{$section}{additional_rmd_files});
  for my $additional (@additional_rtemplates){
    my $additional_rtempalte = dirname(__FILE__) . "/" . trim($additional);
    my $additional_rfile = $result_dir . "/" . basename($additional);
    copy( $additional_rtempalte, $additional_rfile ) or die "Copy failed: $!";
  }

  my $parametersample_files1 = writeParameterSampleFile($config, $section, $result_dir, 1);
  my $parametersample_files2 = writeParameterSampleFile($config, $section, $result_dir, 2);

  my $final = $self->open_pbs( $final_pbs, $pbs_desc, $final_log_desp, $path_file, $result_dir );
  
  print $final "R --slave -e \"library(knitr);rmarkdown::render('" . $final_file . "');\" \n";

  $self->close_pbs( $final, $final_pbs );
}

1;
