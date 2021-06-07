#!/usr/bin/perl
package CQS::UniqueR;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::UniqueWrapper;
use CQS::StringUtils;
use File::Spec;

our @ISA = qw(CQS::UniqueWrapper);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_uniqueR";
  bless $self, $class;
  return $self;
}

sub getRcode {
  my ( $self, $config, $section ) = @_;
  return get_option( $config, $section, "rCode", "" );
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );
  my $task_suffix = $self->{_task_suffix};

  my $rCode = $self->getRcode( $config, $section );
  my $output_file     = get_option( $config, $section, "output_file",     "" );
  my $output_to_result_directory = get_option( $config, $section, "output_to_result_directory", 0 );

  my $removeEmpty = get_option( $config, $section, "remove_empty_parameter", 0 );
  my $paramSampleFiles = {};
  for my $myidx (1..10) {
    my $paramSampleFile = writeParameterSampleFile( $config, $section, $result_dir, $myidx, $removeEmpty );
    if ($paramSampleFile ne ""){
      $paramSampleFiles->{"parSampleFile" . $myidx} = $paramSampleFile;
    }
  }

  my $parameterFile1 = parse_param_file( $config, $section, "parameterFile1", 0 );
  my $parameterFile2 = parse_param_file( $config, $section, "parameterFile2", 0 );
  my $parameterFile3 = parse_param_file( $config, $section, "parameterFile3", 0 );
  if ( !defined($parameterFile1) ) {
    $parameterFile1 = "";
  }
  if ( !defined($parameterFile2) ) {
    $parameterFile2 = "";
  }
  if ( !defined($parameterFile3) ) {
    $parameterFile3 = "";
  }

  my $rfile = $result_dir . "/${task_name}${task_suffix}.r";
  open( my $rf, ">$rfile" ) or die "Cannot create $rfile";
  print $rf "rm(list=ls()) \n";

  my $result_files = $self->result( $config, $section )->{$task_name};
  my $final_file = $result_files->[-1];
  my $output_file_r;
  if ( $output_file eq "parameterSampleFile1" or $output_file eq "parameterSampleFile2" or $output_file eq "parameterSampleFile3" ) {
    $output_file_r = "";
  }
  else {
    $output_file_r = $task_name . $output_file;
  }

  my $rParameter = "outFile='$output_file_r'\n";
  for my $file_key (sort keys %$paramSampleFiles) {
    $rParameter = $rParameter . $file_key . "='" . $paramSampleFiles->{$file_key} . "'\n";
  }
  if ( defined($parameterFile1) ) {
    $rParameter = $rParameter . "parFile1='$parameterFile1'\n";
  }
  if ( defined($parameterFile2) ) {
    $rParameter = $rParameter . "parFile2='$parameterFile2'\n";
  }
  if ( defined($parameterFile3) ) {
    $rParameter = $rParameter . "parFile3='$parameterFile3'\n";
  }
  if ($output_to_result_directory) {
    $rParameter = $rParameter . "outputDirectory='.'\n";
  }

  if ( $rParameter ne "" ) {
    print $rf $rParameter;
  }
  if ( defined($rCode) ) {
    print $rf $rCode . "\n";
  }
  
  print $rf "\nsetwd('$result_dir')\n\n";

  my $rtemplates = get_option( $config, $section, "rtemplate" );
  my @rtemplates = split( ",|;", $rtemplates );
  foreach my $rtemplate (@rtemplates) {
    my $is_absolute = File::Spec->file_name_is_absolute($rtemplate);
    if ( !$is_absolute ) {
      $rtemplate = dirname(__FILE__) . "/$rtemplate";
    }
    if ( !( -e $rtemplate ) ) {
      die("rtemplate $rtemplate defined but not exists!");
    }
    open( my $rt, "<$rtemplate" ) or die $!;
    while ( my $row = <$rt> ) {
      chomp($row);
      $row =~ s/\r//g;
      print $rf "$row\n";
    }
    close($rt);
  }
  close($rf);

  my $rReportTemplates = get_option( $config, $section, "rReportTemplate", "" );
  if ( $rReportTemplates ne "" ) {
    #my $rReportFile = $result_dir . "/" . basename($rReportTemplates);
    my $rReportFile = $result_dir . "/" . basename($rReportTemplates);
    open( my $rf, ">$rReportFile" ) or die "Cannot create $rReportFile";

    my @rReportTemplates = split( ",|;", $rReportTemplates );
    foreach my $rReportTemplate (@rReportTemplates) {
      my $is_absolute = File::Spec->file_name_is_absolute($rReportTemplate);
      if ( !$is_absolute ) {
        $rReportTemplate = dirname(__FILE__) . "/$rReportTemplate";
      }
      if ( !( -e $rReportTemplate ) ) {
        die("rReportTemplate $rReportTemplate defined but not exists!");
      }
      open( my $rt, "<$rReportTemplate" ) or die $!;
      while ( my $row = <$rt> ) {
        chomp($row);
        $row =~ s/\r//g;
        print $rf "$row\n";
      }
      close($rt);
    }
    close($rf);
  }

  my $additional_rmd_files = get_option( $config, $section, "additional_rmd_files","" );
  if ( $additional_rmd_files ne "" ) {
    my @additional_rmd_files = split( ",|;", $additional_rmd_files );
    foreach my $additional_rmd_files (@additional_rmd_files) {
      my $is_absolute = File::Spec->file_name_is_absolute($additional_rmd_files);
      if ( !$is_absolute ) {
        $additional_rmd_files = dirname(__FILE__) . "/$additional_rmd_files";
      }
      if ( !( -e $additional_rmd_files ) ) {
        die("additional_rmd_files $additional_rmd_files defined but not exists!");
      }
      my $additional_rmd_files_inResult = $result_dir . "/" . basename($additional_rmd_files);
      open( my $rf, ">$additional_rmd_files_inResult" ) or die "Cannot create $additional_rmd_files_inResult";
      open( my $rt, "<$additional_rmd_files" ) or die $!;
      while ( my $row = <$rt> ) {
        chomp($row);
        $row =~ s/\r//g;
        print $rf "$row\n";
      }
      close($rt);
      close($rf);
    }
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  my $rlibs = get_option_include_general($config, $section, "R_LIBS", "");
  if($rlibs ne ""){
    print $pbs "export R_LIBS=$rlibs \n\n";
  }

  my $rscript = get_option_include_general($config, $section, "Rscript", "Rscript");
  if ( defined($option) and $option ne "" ) {
    print $pbs "$rscript --vanilla " . basename($rfile) . " $option";
  }
  else {
    print $pbs "$rscript --vanilla " . basename($rfile);
  }
  $self->close_pbs( $pbs, $pbs_file );
}

1;
