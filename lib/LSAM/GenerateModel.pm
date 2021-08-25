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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $generate_iteration = get_option($config, $section, "generate_iteration", 0);
  if ($generate_iteration) {
    $option = $option . " --generate_iteration";
  }
  my $methods = get_raw_files( $config, $section );
  my $timeRanges = get_option( $config, $section, "time_ranges" );
  my $dsaFiles = $config->{$section}{dsa_files};

  my $mortalityFiles = $config->{$section}{motarlity_files};

  my $defOption = "";
  if ( defined $config->{$section}{"def_data_def"} ) {
    $defOption = $defOption . " --defdatadef " . $config->{$section}{"def_data_def"};
  }
  if ( defined $config->{$section}{"opt_data_def"} ) {
    $defOption = $defOption . " --optdatadef " . $config->{$section}{"opt_data_def"};
  }
  if ( defined $config->{$section}{"border_state_file"} ) {
    $defOption = $defOption . " --border " . $config->{$section}{"border_state_file"};
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

    die "template file not exists: " . $template if !-e $template;

    for my $methodName ( keys %$methods ) {
      my $methodFile = $methods->{$methodName}[0];
      die "method file not exists: " . $methodFile if !-e $methodFile;

      if ( defined $dsaFiles ) {
        my @dsaNames = keys %$dsaFiles;
        for my $dsaName (@dsaNames) {
          my $dsaFile = $dsaFiles->{$dsaName}[0];
          die "dsa file not exists: " . $dsaFile if !-e $dsaFile;

          my $key = $timeName . "_" . $methodName;
          if ( scalar(@dsaNames) > 1 ) {
            $key = $key . "_" . $dsaName;
          }
          my $finalFile = $key . ".inp";

          print $pbs "
python $py_script $option -i $template -o $finalFile $defOption --name $key --method $methodFile --dsa $dsaFile --start \"$start\" --end \"$end\" 
";
        }
      }
      elsif ( defined $mortalityFiles ) {
        my @mNames = keys %$mortalityFiles;
        for my $mName (@mNames) {
          my $mFile = $mortalityFiles->{$mName};
          die "motarlity file not exists: " . $mFile if !-e $mFile;

          my $key = $timeName . "_" . $methodName;
          if ( scalar(@mNames) > 1 ) {
            $key = $key . "_" . $mName;
          }
          my $finalFile = $key . ".inp";

          print $pbs "
python $py_script $option -i $template -o $finalFile $defOption --name $key --method $methodFile --mortality $mFile --startTime \"$start\" --endTime \"$end\" 
";
        }
      }
      else {
        my $finalFile = $timeName . "_" . $methodName . ".inp";
        print $pbs "
python $py_script $option -i $template -o $finalFile $defOption --name ${timeName}_${methodName} --method $methodFile --startTime \"$start\" --endTime \"$end\"
";
      }
    }
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub get_result_files {
  my ($result_dir, $key, $generate_iteration) = @_;
  if ($generate_iteration) {
    my $result_files = [ 
      "$result_dir\\${key}.01.inp",
      "$result_dir\\${key}.02.inp",
      "$result_dir\\${key}.03.inp",
      "$result_dir\\${key}.04.inp",
      "$result_dir\\${key}.05.inp",
      "$result_dir\\${key}.06.inp",
      "$result_dir\\${key}.07.inp",
      "$result_dir\\${key}.08.inp",
      "$result_dir\\${key}.09.inp",
      "$result_dir\\${key}.10.inp" ];
    return ($result_files);
  }else{
    my $finalFile    = $key . ".inp";
    my $result_files = [ $result_dir . "\\" . $finalFile ];
    return ($result_files);
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $result_dir =~ s/\//\\/g;

  my $methods = get_raw_files( $config, $section );
  my $timeRanges     = get_option( $config, $section, "time_ranges" );
  my $dsaFiles       = $config->{$section}{dsa_files};
  my $mortalityFiles = $config->{$section}{motarlity_files};
  my $generate_iteration = get_option($config, $section, "generate_iteration", 0);

  my $result = {};
  for my $timeName ( keys %$timeRanges ) {
    for my $methodName ( keys %$methods ) {
      if ( defined $dsaFiles ) {
        my @dsaNames = keys %$dsaFiles;
        for my $dsaName (@dsaNames) {
          my $key = $timeName . "_" . $methodName;
          if ( scalar(@dsaNames) > 1 ) {
            $key = $key . "_" . $dsaName;
          }
          $result->{$key} = filter_array( get_result_files($result_dir, $key, $generate_iteration), $pattern );
        }
      }
      elsif ( defined $mortalityFiles ) {
        my @mNames = keys %$mortalityFiles;
        for my $mName (@mNames) {
          my $key = $timeName . "_" . $methodName;
          if ( scalar(@mNames) > 1 ) {
            $key = $key . "_" . $mName;
          }
          $result->{$key} = filter_array( get_result_files($result_dir, $key, $generate_iteration), $pattern );
        }
      }
      else {
        my $key          = $timeName . "_" . $methodName;
        $result->{$key} = filter_array( get_result_files($result_dir, $key, $generate_iteration), $pattern );
      }
    }
  }

  return $result;
}
1;
