#!/usr/bin/perl
package CQS::SequenceTask;

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
use CQS::ClassFactory;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "SequenceTask";
  $self->{_suffix} = "_st";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %fqFiles = %{ get_raw_files( $config, $section ) };
  my $cluster = get_cluster($config, $section);

  for my $taskName ( sort keys %fqFiles ) {
    my $shfile = $self->taskfile( $pbsDir, $taskName );
    open( SH, ">$shfile" ) or die "Cannot create $shfile";
    print SH get_run_command($sh_direct);

    my @tasks = @{ $fqFiles{$taskName} };

    my $samples = {};
    my $taskpbs = {};
    for my $tasksection (@tasks) {
      #print "task " . $tasksection . " ...\n";
      my $pbsfiles = getPbsFiles( $config, $tasksection );
      for my $sample ( sort keys %{$pbsfiles} ) {
        #print "\t", $sample, " => ", $pbsfiles->{$sample}, "\n";
        $samples->{$sample} = 1;
      }

      $taskpbs->{$tasksection} = $pbsfiles;
      for my $sample ( sort keys %{$pbsfiles} ) {
        $samples->{$sample} = 1;
      }
    }

    for my $sample ( sort keys %{$samples} ) {
      my $taskSample = $taskName . "_" . $sample;

      my $pbsFile = $self->pbsfile( $pbsDir, $taskSample );
      my $pbsName = basename($pbsFile);
      my $log     = $self->logfile( $logDir, $taskSample );
      my $logdesp = $cluster->get_log_desc($log);

      open( OUT, ">$pbsFile" ) or die $!;

      print OUT "$pbsDesc
$logdesp

$path_file

echo sequenceTaskStart=`date` 
";
      for my $tasksection (@tasks) {

        #print "task " . $tasksection . " ...\n";
        my $pbsfiles = $taskpbs->{$tasksection};

        if ( exists $pbsfiles->{$sample} ) {
          my $samplepbs = $pbsfiles->{$sample};
          if ( ref($samplepbs) eq 'ARRAY' ) {
            for my $pbs ( @{$samplepbs} ) {
              if ( !-e $pbs ) {
                die "Task " . $tasksection . ", file not exists " . $pbs . "\n";
              }
              print OUT "bash " . $pbs . "\n";
            }
          }
          else {
            if ( !-e $samplepbs ) {
              die "Task " . $tasksection . ", file not exists " . $samplepbs . "\n";
            }
            print OUT "bash " . $samplepbs . "\n";
          }
        }
      }

      print OUT "
echo sequenceTaskEnd=`date` 
";
      close(OUT);

      print SH "\$MYCMD ./$pbsName \n";
      print "$pbsFile created\n";
    }
    print SH "exit 0\n";
    close(SH);

    if ( is_linux() ) {
      chmod 0755, $shfile;
    }
  }
}

sub pbsfiles {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = {};

  my %fqFiles = %{ get_raw_files( $config, $section ) };

  for my $taskName ( sort keys %fqFiles ) {
    my @tasks = @{ $fqFiles{$taskName} };

    my $samples = {};
    my $taskpbs = {};
    for my $tasksection (@tasks) {
      #print "task " . $tasksection . " ...\n";
      my $pbsfiles = getPbsFiles( $config, $tasksection );
      for my $sample ( sort keys %{$pbsfiles} ) {
        #print "\t", $sample, " => ", $pbsfiles->{$sample}, "\n";
        $samples->{$sample} = 1;
      }

      $taskpbs->{$tasksection} = $pbsfiles;
      for my $sample ( sort keys %{$pbsfiles} ) {
        $samples->{$sample} = 1;
      }
    }

    for my $sample ( sort keys %{$samples} ) {
      my $taskSample = $taskName . "_" . $sample;
      $result->{$taskSample} = $self->pbsfile( $pbsDir, $taskSample );
    }
  }

  return $result;
}
1;
