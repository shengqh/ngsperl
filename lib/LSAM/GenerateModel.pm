#!/usr/bin/perl
package LSAM::GenerateModel;

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

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $methods = get_raw_files( $config, $section );
  my $timeRanges = get_option( $config, $section, "time_ranges" );
  my $dsaFiles   = $config->{$section}{dsa_files};
  
  my $defOption = "";
  if(defined $config->{$section}{"def_data_def"}){
    $defOption = $defOption . " --defdatadef " . $config->{$section}{"def_data_def"};
  }
  if(defined $config->{$section}{"opt_data_def"}){
    $defOption = $defOption . " --optdatadef " . $config->{$section}{"opt_data_def"};
  }

  my $pbs_file = $self->get_file( $pbs_dir, $task_name, ".bat" );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $py_script = dirname(__FILE__) . "/generateModel.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }

  my $pbs = $self->open_pbs( $pbs_file, "", "", $path_file, $result_dir );

  for my $timeName ( keys %$timeRanges ) {
    my $start    = $timeRanges->{$timeName}->[0];
    my $end      = $timeRanges->{$timeName}->[1];
    my $template = $timeRanges->{$timeName}->[2];
    for my $methodName ( keys %$methods ) {
      my $methodFile = $methods->{$methodName}[0];
      if ( defined $dsaFiles ) {
        for my $dsaName ( keys %$dsaFiles ) {
          my $dsaFile   = $dsaFiles->{$dsaName}[0];
          my $finalFile = $timeName . "_" . $methodName . "_" . $dsaName . ".inp";

          print $pbs "
python $py_script -i $template -o $finalFile $defOption --name ${timeName}_${methodName} --method $methodFile --dsa $dsaFile --start \"$start\" --end \"$end\" 
";
        }
      }
      else {
        my $finalFile = $timeName . "_" . $methodName . ".inp";
        print $pbs "
python $py_script -i $template -o $finalFile $defOption --name ${timeName}_${methodName} --method $methodFile --start \"$start\" --end \"$end\"
";
      }
    }
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  $result_dir =~ s/\//\\/g;

  my $methods = get_raw_files( $config, $section );
  my $timeRanges = get_option( $config, $section, "time_ranges" );
  my $dsaFiles = $config->{$section}{dsa_files};

  my $result = {};
  for my $timeName ( keys %$timeRanges ) {
    for my $methodName ( keys %$methods ) {
      if ( defined $dsaFiles ) {
        for my $dsaName ( keys %$dsaFiles ) {
          my $key          = $timeName . "_" . $methodName . "_" . $dsaName;
          my $finalFile    = $key . ".inp";
          my @result_files = ( $result_dir . "\\" . $finalFile );
          $result->{$key} = filter_array( \@result_files, $pattern );
        }
      }
      else {
        my $key          = $timeName . "_" . $methodName;
        my $finalFile    = $key . ".inp";
        my @result_files = ( $result_dir . "\\" . $finalFile );
        $result->{$key} = filter_array( \@result_files, $pattern );
      }
    }
  }

  return $result;
}
1;
