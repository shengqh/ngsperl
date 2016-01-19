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

  my %fqFiles = %{ get_raw_files( $config, $section ) };

  #print Dumper(\%fqFiles);

  my $cluster = get_cluster( $config, $section );

  my $finalName    = $task_name . "_pipeline";
  my $finalpbs     = $self->get_pbs_filename( $pbs_dir, $finalName );
  my $finallog     = $self->get_log_filename( $log_dir, $finalName );
  my $finallogdesp = $cluster->get_log_description($finallog);

  my $final = $self->open_pbs( $finalpbs, $pbs_desc, $finallogdesp, $path_file, $result_dir );

  for my $taskName ( sort keys %fqFiles ) {
    my $shfile = $self->get_task_filename( $pbs_dir, $taskName );
    open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
    print $sh get_run_command($sh_direct);

    my @tasks = @{ $fqFiles{$taskName} };

    my $samples = {};
    my $taskpbs = {};
    for my $tasksection (@tasks) {

      #print "task " . $tasksection . " ...\n";
      my $get_pbs_files = getget_pbs_files( $config, $tasksection );
      for my $sample ( sort keys %{$get_pbs_files} ) {

        #print "\t", $sample, " => ", $get_pbs_files->{$sample}, "\n";
        $samples->{$sample} = 1;
      }

      $taskpbs->{$tasksection} = $get_pbs_files;
      for my $sample ( sort keys %{$get_pbs_files} ) {
        $samples->{$sample} = 1;
      }
    }

    for my $sample ( sort keys %{$samples} ) {
      my $taskSample = $sample . "_" . $taskName;

      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $taskSample );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_log_filename( $log_dir, $taskSample );
      my $logdesp  = $cluster->get_log_description($log);

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $logdesp, $path_file, $result_dir );

      for my $tasksection (@tasks) {

        #print "task " . $tasksection . " ...\n";
        my $get_pbs_files = $taskpbs->{$tasksection};

        if ( exists $get_pbs_files->{$sample} ) {
          my $samplepbs = $get_pbs_files->{$sample};
          if ( ref($samplepbs) eq 'ARRAY' ) {
            for my $pbs ( @{$samplepbs} ) {
              if ( !-e $pbs ) {
                die "Task " . $tasksection . ", file not exists " . $pbs . "\n";
              }
              print $pbs "bash " . $pbs . "\n";
            }
          }
          else {
            if ( !-e $samplepbs ) {
              die "Task " . $tasksection . ", file not exists " . $samplepbs . "\n";
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

  $self->close_pbs( $final, $finalpbs );

}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = {};

  my %fqFiles = %{ get_raw_files( $config, $section ) };

  for my $taskName ( sort keys %fqFiles ) {
    my @tasks = @{ $fqFiles{$taskName} };

    my $samples = {};
    my $taskpbs = {};
    for my $tasksection (@tasks) {

      #print "task " . $tasksection . " ...\n";
      my $get_pbs_files = getget_pbs_files( $config, $tasksection );
      for my $sample ( sort keys %{$get_pbs_files} ) {

        #print "\t", $sample, " => ", $get_pbs_files->{$sample}, "\n";
        $samples->{$sample} = 1;
      }

      $taskpbs->{$tasksection} = $get_pbs_files;
      for my $sample ( sort keys %{$get_pbs_files} ) {
        $samples->{$sample} = 1;
      }
    }

    for my $sample ( sort keys %{$samples} ) {
      my $taskSample = $sample . "_" . $taskName;
      $result->{$taskSample} = $self->get_pbs_filename( $pbs_dir, $taskSample );
    }
  }

  return $result;
}
1;
