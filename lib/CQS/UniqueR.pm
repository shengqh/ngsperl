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
use CQS::UniqueTask;
use CQS::StringUtils;
use File::Spec;

our @ISA = qw(CQS::UniqueTask);

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  my $rCode = $self->getRcode( $config, $section );
  my $output_file     = get_option( $config, $section, "output_file",     "" );
  my $output_file_ext = get_option( $config, $section, "output_file_ext", "" );
  my @output_file_exts = split( ";", $output_file_ext );
  if ( scalar(@output_file_exts) == 0 ) {
    push( @output_file_exts, "" );
  }
  my $output_to_result_directory = get_option( $config, $section, "output_to_result_directory", 0 );

  my $removeEmpty = get_option( $config, $section, "remove_empty_parameter", 0 );
  my $parametersample_files1 = writeParameterSampleFile( $config, $section, $result_dir, 1, $removeEmpty );
  my $parametersample_files2 = writeParameterSampleFile( $config, $section, $result_dir, 2, $removeEmpty );
  my $parametersample_files3 = writeParameterSampleFile( $config, $section, $result_dir, 3, $removeEmpty );
  my $parametersample_files4 = writeParameterSampleFile( $config, $section, $result_dir, 4, $removeEmpty );

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
  if ( defined($parametersample_files1) ) {
    $rParameter = $rParameter . "parSampleFile1='$parametersample_files1'\n";
  }
  if ( defined($parametersample_files2) ) {
    $rParameter = $rParameter . "parSampleFile2='$parametersample_files2'\n";
  }
  if ( defined($parametersample_files3) ) {
    $rParameter = $rParameter . "parSampleFile3='$parametersample_files3'\n";
  }
  if ( defined($parametersample_files4) ) {
    $rParameter = $rParameter . "parSampleFile4='$parametersample_files4'\n";
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
    my $rfile = $result_dir . "/" . basename($rReportTemplates);
    open( my $rf, ">$rfile" ) or die "Cannot create $rfile";

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

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  if ( defined($option) and $option ne "" ) {
    print $pbs "R --vanilla --slave -f $rfile --args $option";

    #    "R --vanilla --slave -f $rfile --args $task_name$output_file $option $parametersample_files1 $parametersample_files2 $parametersample_files3 $parameterFile1 $parameterFile2 $parameterFile3";
  }
  else {
    print $pbs "R --vanilla --slave -f $rfile";
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $output_file                = get_option( $config, $section, "output_file",                "" );
  my $output_file_ext            = get_option( $config, $section, "output_file_ext",            "" );
  my $output_to_result_directory = get_option( $config, $section, "output_to_result_directory", 0 );
  
  my @output_file_exts;
  if ( $output_file_ext =~ /;/ ) {
    @output_file_exts = split( ";", $output_file_ext );
  }else{
    @output_file_exts = ($output_file_ext );
  }

  my $result       = {};
  my @result_files = ();

  if ( $output_file eq "parameterSampleFile1" or $output_file eq "parameterSampleFile2" or $output_file eq "parameterSampleFile3" ) {
    if ( has_raw_files( $config, $section, $output_file ) ) {
      my %temp = %{ get_raw_files( $config, $section, $output_file ) };
      foreach my $sample_name ( keys %temp ) {
        if ( ref( $temp{$sample_name} ) eq "HASH" ) {
          foreach my $output_file_ext_one (@output_file_exts) {
            push( @result_files, "${result_dir}/${sample_name}${output_file_ext_one}" );
          }
        }
        else {
          foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
            my $subSampleName = $output_to_result_directory ? $result_dir . "/" . basename($subSampleFile) : $subSampleFile;
            foreach my $output_file_ext_one (@output_file_exts) {
              push( @result_files, "${subSampleName}${output_file_ext_one}" );
            }
          }
        }
      }
    }
    else {
      die "output_file defined as " . $output_file . ", but " . $output_file . " not in configure\n";
    }
  }
  else {
    foreach my $output_file_ext_one (@output_file_exts) {
      push( @result_files, "${result_dir}/${task_name}${output_file}${output_file_ext_one}" );
    }
  }

  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
