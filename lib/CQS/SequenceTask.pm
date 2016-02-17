#!/usr/bin/perl
package CQS::SequenceTask;

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::ClassFactory;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_st";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %step_map = %{ get_raw_files( $config, $section ) };

  #print Dumper(\%step_map);

  my $cluster = get_cluster( $config, $section );

  my $final_name     = $task_name . "_pipeline";
  my $final_pbs      = $self->get_pbs_filename( $pbs_dir, $final_name );
  my $final_log      = $self->get_log_filename( $log_dir, $final_name );
  my $final_log_desp = $cluster->get_log_description($final_log);

  my $summary_name     = $task_name . "_summary";
  my $summary_pbs      = $self->get_pbs_filename( $pbs_dir, $summary_name );
  my $summary_log      = $self->get_log_filename( $log_dir, $summary_name );
  my $summary_log_desp = $cluster->get_log_description($summary_log);

  my $result_list_file = $self->get_file( $result_dir, $task_name, "_expect_result.tsv" );
  
  my $summary = $self->open_pbs( $summary_pbs, $pbs_desc, $summary_log_desp, $path_file, $result_dir );
  my $rtemplate = dirname(__FILE__) . "/summaryResultFiles.R";
  print $summary "R --vanilla --slave -f $rtemplate --args $result_list_file \n";
  $self->close_pbs( $summary, $summary_pbs );
  
  my $final   = $self->open_pbs( $final_pbs,   $pbs_desc, $final_log_desp,   $path_file, $pbs_dir );
  
  open( my $result_list, ">$result_list_file" ) or die $!;
  print $result_list "StepName\tTaskName\tSampleName\tFileList\n";

  for my $step_name ( sort keys %step_map ) {
    my $shfile = $self->get_task_filename( $pbs_dir, $step_name );
    open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
    print $sh get_run_command($sh_direct);

    my @tasks = @{ $step_map{$step_name} };

    my $samples = {};
    my $taskpbs = {};
    for my $task_section (@tasks) {
      my $classname = $config->{$task_section}{class};
      if ( !defined $classname ) {
        die "$task_section is not a valid task section.";
      }
      my $myclass         = instantiate( $classname );
      my %expect_file_map;
      eval {
        %expect_file_map = %{ $myclass->result( $config, $task_section ) };
      } or do {
        my $e = $@;
        die ("Something went wrong to get result of section $task_section : $e\n");
      };
      
      my $pbs_file_map    = $myclass->get_pbs_files( $config, $task_section );
      for my $expect_name ( sort keys %expect_file_map ) {
        my $expect_files = $expect_file_map{$expect_name};
        my $expect_file_list = join( ",", @{$expect_files} );
        print $result_list $step_name, "\t", $task_section, "\t", $expect_name, "\t", $expect_file_list, "\n";
      }

      #print "task " . $task_section . " ...\n";
      $taskpbs->{$task_section} = $pbs_file_map;

      for my $sample ( sort keys %{$pbs_file_map} ) {
        $samples->{$sample} = 1;
      }
    }

    for my $sample ( sort keys %{$samples} ) {
      my $pbs_file = $self->get_step_sample_pbs( $pbs_dir, $step_name, $sample );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_step_sample_log( $log_dir, $step_name, $sample );
      my $logdesp  = $cluster->get_log_description($log);

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $logdesp, $path_file, $result_dir );

      for my $task_section (@tasks) {

        #print "task " . $task_section . " ...\n";
        my $pbs_files = $taskpbs->{$task_section};

        if ( exists $pbs_files->{$sample} ) {
          my $samplepbs = $pbs_files->{$sample};
          if ( ref($samplepbs) eq 'ARRAY' ) {
            for my $pbs ( @{$samplepbs} ) {
              if ( !-e $pbs ) {
                die "Task " . $task_section . ", file not exists " . $pbs . "\n";
              }
              print $pbs "bash " . $pbs . "\n";
            }
          }
          else {
            if ( !-e $samplepbs ) {
              die "Task " . $task_section . ", file not exists " . $samplepbs . "\n";
            }
            print $pbs "bash " . $samplepbs . "\n";
          }
        }
      }

      $self->close_pbs( $pbs, $pbs_file );

      print $sh "\$MYCMD ./$pbs_name \n";
      print $final "bash ./$pbs_name \n";

    }
    print $sh "exit 0\n";
    close $sh;

    if ( is_linux() ) {
      chmod 0755, $shfile;
    }
  }

  close($result_list);
  
  my $summary_pbs_name = basename($summary_pbs);
  print $final "\nbash ./$summary_pbs_name \n";
  $self->close_pbs( $final, $final_pbs );
}

sub get_step_sample_name {
  my ( $self, $step_name, $sample ) = @_;
  return $sample . "_" . $step_name;
}

sub get_step_sample_pbs {
  my ( $self, $pbs_dir, $step_name, $sample ) = @_;
  my $task_sample = $self->get_step_sample_name( $step_name, $sample );
  return $self->get_pbs_filename( $pbs_dir, $task_sample );
}

sub get_step_sample_log {
  my ( $self, $log_dir, $step_name, $sample ) = @_;
  my $task_sample = $self->get_step_sample_name( $step_name, $sample );
  return $self->get_log_filename( $log_dir, $task_sample );
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = {};

  my %step_map = %{ get_raw_files( $config, $section ) };

  for my $step_name ( sort keys %step_map ) {
    my @tasks = @{ $step_map{$step_name} };

    my $samples = {};
    my $taskpbs = {};
    for my $task_section (@tasks) {
      my $classname = $config->{$task_section}{class};
      if ( !defined $classname ) {
        die "$task_section is not a valid task section.";
      }
      my $myclass = instantiate( $classname );
      my $pbs_file_map = $myclass->get_pbs_files( $config, $task_section );
      for my $sample ( sort keys %{$pbs_file_map} ) {
        $samples->{$sample} = 1;
      }
    }

    for my $sample ( sort keys %{$samples} ) {
      my $step_sample = $self->get_step_sample_name( $step_name, $sample );
      $result->{$step_sample} = $self->get_step_sample_pbs( $pbs_dir, $step_name, $sample );
    }
  }
  return $result;
}

1;
