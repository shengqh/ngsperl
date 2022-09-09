#!/usr/bin/perl
package CQS::CombineTask;

use strict;
use warnings;
use File::Basename;
use File::Copy;
use Data::Dumper;
use CQS::Task;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ct";
  bless $self, $class;
  return $self;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section, 0 );

  my $curSection = get_config_section( $config, $section );
  my $tasks = $curSection->{source};

  my $samples = {};
  for my $taskSectionName (@$tasks) {
    my $taskSection = get_config_section( $config, $taskSectionName );
    my $classname = $taskSection->{class};
    if ( !defined $classname ) {
      die "$taskSectionName is not a valid task section.";
    }
    my $myclass = instantiate($classname);
    my $pbs_file_map = $myclass->get_pbs_files( $config, $taskSectionName );
    for my $sample ( sort keys %{$pbs_file_map} ) {
      $samples->{$sample} = 1;
    }
  }

  my $result = {};
  for my $sample ( sort keys %$samples ) {
    $result->{$sample} = $self->get_pbs_filename( $pbs_dir, $sample );
  }
  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $curSection = get_config_section( $config, $section );
  my $tasks = $curSection->{source};

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my $samples     = {};
  my $taskpbs     = {};
  my $expectFiles = {};
  my $clears      = {};

  my $lastSectionName = $tasks->[ scalar(@$tasks) - 1 ];
  for my $taskSectionName (@$tasks) {
    my $taskSection = get_config_section( $config, $taskSectionName );
    my $classname = $taskSection->{class};
    if ( !defined $classname ) {
      die "$taskSectionName is not a valid task section.";
    }
    my $myclass = instantiate($classname);
    $myclass->perform($config, $taskSectionName);
    $expectFiles->{$taskSectionName} = $myclass->result( $config, $taskSectionName );
    $clears->{$taskSectionName} = $myclass->get_clear_map( $config, $taskSectionName );
    $taskpbs->{$taskSectionName} = $myclass->get_pbs_files( $config, $taskSectionName );
    for my $sample (keys %{$taskpbs->{$taskSectionName}}){
      $samples->{$sample}=1;
    }
  }

  for my $sample ( sort keys %$samples ) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample );
    my $logdesp  = $cluster->get_log_description($log);

    my $lastExpects= $expectFiles->{$lastSectionName};
    my $lastExpectResult = $lastExpects->{$sample};
    if(!defined $lastExpectResult){
      die "No expect result from $lastSectionName"
    }
    my $finalFile = $expectFiles->{$lastSectionName}->{$sample}[0];
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $logdesp, $path_file, $result_dir, $finalFile );

    for my $taskSectionName (@$tasks) {
      my $pbs_map = $taskpbs->{$taskSectionName};
      if ( defined $pbs_map->{$sample} ) {
        my $samplepbs = $pbs_map->{$sample};
        if ( is_array($samplepbs) ) {
          for my $subpbs ( @{$samplepbs} ) {
            print $pbs "bash " . $subpbs . "\n";
          }
        }
        else {
          print $pbs "bash " . $samplepbs . "\n";
        }
      }
    }

    print $pbs "
if [[ ( -s $finalFile ) || ( -d $finalFile ) ]]; then
";
    for my $taskSectionName (@$tasks) {
      if ( $taskSectionName eq $lastSectionName ) {
        last;
      }
      my $clear_map = $clears->{$taskSectionName};
      if ( defined $clear_map->{$sample} ) {
        my $clearFiles = $clear_map->{$sample};
        if ( is_array($clearFiles) ) {
          for my $clearFile (@$clearFiles) {
            print $pbs "  rm " . $clearFile . "\n";
          }
        }
        else {
          print $pbs "  rm " . $clearFiles . "\n";
        }
      }
    }
    print $pbs "
fi
";
    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }
  print $sh "exit 0\n";
  close $sh;
  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
}

1;
