#!/usr/bin/perl
package CQS::SequenceTaskSlurm;

use strict;
use warnings;
use File::Basename;
use File::Copy;
use Data::Dumper;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use CQS::SequenceTask;

our @ISA = qw(CQS::SequenceTask);


sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %step_map = %{ get_raw_files( $config, $section ) };
  
  my $task_shell = get_option($config, $section, "task_shell", "bash");

  #print Dumper(\%step_map);

  my $cluster = get_cluster( $config, $section );

  my $final_name     = $task_name . "_pipeline";
  my $final_pbs      = $self->get_pbs_filename( $pbs_dir, $final_name );
  my $final_log      = $self->get_log_filename( $log_dir, $final_name );
  my $final_log_desp = $cluster->get_log_description($final_log);

  my $summary_name = $task_name . "_summary";

  my $report_name     = $task_name . "_report";
  my $report_pbs      = $self->get_pbs_filename( $pbs_dir, $report_name );
  my $report_log      = $self->get_log_filename( $log_dir, $report_name );
  my $report_log_desp = $cluster->get_log_description($report_log);
  my $summary_pbs     = $self->get_pbs_filename( $pbs_dir, $summary_name );

  my $result_list_file = $self->get_file( $result_dir, $task_name, "_expect_result.tsv" );

  my $report  = $self->open_pbs( $report_pbs,  $pbs_desc, $report_log_desp, $path_file, $result_dir );
  my $summary = $self->open_pbs( $summary_pbs, $pbs_desc, $report_log_desp, $path_file, $result_dir );

  #Make Summary Figure
  my $rtemplate = dirname(__FILE__) . "/summaryResultFiles.R";
  my $rfile     = $result_dir . "/${summary_name}.r";
  open( my $rf, ">$rfile" )     or die "Cannot create $rfile";
  open( my $rt, "<$rtemplate" ) or die $!;
  print $rf "parFile1='$result_list_file'\n";
  while (<$rt>) {
    print $rf $_;
  }
  close($rt);
  close($rf);

  print $report "R --vanilla --slave -f $rfile \n";
  print $summary "R --vanilla --slave -f $rfile \n";
  $self->close_pbs( $summary, $summary_pbs );

  #Make Report
  my $rtemplateReport = dirname(__FILE__) . "/MakeReport.R";
  my $projectRmd      = dirname(__FILE__) . "/ProjectReport.Rmd";
  my $taskRmd         = dirname(__FILE__) . "/TaskReport.Rmd";
  copy( $projectRmd, $result_dir . "/ProjectReport.Rmd" ) or die "Copy failed: $!";
  copy( $taskRmd,    $result_dir . "/TaskReport.Rmd" )    or die "Copy failed: $!";

  my $rfileReport = $result_dir . "/${report_name}.r";
  my $task_dir    = $target_dir;
  $task_dir =~ s/\/sequencetask$//;
  open( my $rfReport, ">$rfileReport" )     or die "Cannot create $rfileReport";
  open( my $rtReport, "<$rtemplateReport" ) or die $!;
  print $rfReport "parFile1='$task_name'\n";
  print $rfReport "parFile2='$task_dir'\n";

  while (<$rtReport>) {
    print $rfReport $_;
  }
  close($rtReport);
  close($rfReport);

  print $report "R --vanilla --slave -f $rfileReport \n";
  $self->close_pbs( $report, $report_pbs );

  my $final = $self->open_pbs( $final_pbs, $pbs_desc, $final_log_desp, $path_file, $pbs_dir );

  open( my $result_list, ">$result_list_file" ) or die $!;
  print $result_list "StepName\tTaskName\tSampleName\tFileList\tCanFileEmpty\n";

  for my $step_name ( sort keys %step_map ) {
    my $shfile = $self->get_task_filename( $pbs_dir, $step_name );
    open( my $sh, ">$shfile" ) or die "Cannot create $shfile";

    my @tasks = @{ $step_map{$step_name} };

    my $samples    = {};
    my $taskpbs    = {};
    my $clears     = {};
    my $clear_keys = {};

    my $count = 0;
    for my $task_section (@tasks) {
      print $final "# " . $task_section . "\n";
      $count = $count + 1;
      my $classname = $config->{$task_section}{class};
      if ( !defined $classname ) {
        die "$task_section is not a valid task section.";
      }
      my $myclass = instantiate($classname);

      my $can_result_be_empty_file = $myclass->can_result_be_empty_file( $config, $task_section ) ? "True" : "False";
      my $expect_file_map;
      eval { $expect_file_map = $myclass->result( $config, $task_section ); } or do {
        my $e = $@;
        die("Something went wrong to get result of section $task_section : $e\n");
      };

      if ( !$config->{$task_section}{not_clean} ) {
        my $clear_map = $myclass->get_clear_map( $config, $task_section );
        $clears->{$task_section} = $clear_map;
        for my $sample ( sort keys %$clear_map ) {
          $clear_keys->{$sample} = 1;
        }
      }
      my $pbs_file_map = $myclass->get_pbs_files( $config, $task_section );
      for my $expect_name ( sort keys %$expect_file_map ) {
        my $expect_files = $expect_file_map->{$expect_name};
        my $expect_file_list = join( ",", @{$expect_files} );
        print $result_list $step_name, "\t", $task_section, "\t", $expect_name, "\t", $expect_file_list, "\t", $can_result_be_empty_file, "\n";
      }

      for my $sample ( sort keys %{$pbs_file_map} ) {
        my $samplepbs = $pbs_file_map->{$sample};
        if ( ref($samplepbs) eq 'ARRAY' ) {
          for my $subpbs ( @{$samplepbs} ) {
            if ( !-e $subpbs ) {
              die "Task " . $task_section . ", file not exists " . $subpbs . "\n";
            }
            print $final "bash $subpbs \n";
          }
        }
        else {
          if ( !-e $samplepbs ) {
            die "Task " . $task_section . ", file not exists " . $samplepbs . "\n";
          }
          print $final "bash $samplepbs \n";
        }
      }

      #print "task " . $task_section . " ...\n";
      $taskpbs->{$task_section} = $pbs_file_map;

      for my $sample ( sort keys %$pbs_file_map ) {
        $samples->{$sample} = 1;
      }
    }

    #print dumper($expects), "\n";
    #print dumper($samples), "\n";

    for my $sample ( sort keys %{$samples} ) {
      my $pbs_file = $self->get_step_sample_pbs( $pbs_dir, $step_name, $sample );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_step_sample_log( $log_dir, $step_name, $sample );
      my $logdesp  = $cluster->get_log_description($log);

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $logdesp, $path_file, $result_dir );
      
      my $index = 0;
      for my $task_section (@tasks) {

        #print "task " . $task_section . " ...\n";
        my $pbs_map = $taskpbs->{$task_section};
        if ( defined $pbs_map->{$sample} ) {
          my $samplepbs = $pbs_map->{$sample};
          if ( ref($samplepbs) eq 'ARRAY' ) {
            for my $subpbs ( @{$samplepbs} ) {
              my $lastjid = $index == 0?"":"--dependency=afterany:\$jid" . $index;
              $index = $index + 1;
              print $pbs "jid" . $index . "=\$(sbatch $lastjid " . $subpbs . " | awk '{print \$NF}')\n";
            }
          }
          else {
            my $lastjid = $index == 0?"":"--dependency=afterany:\$jid" . $index;
            $index = $index + 1;
            print $pbs "jid" . $index . "=\$(sbatch $lastjid " . $samplepbs . " | awk '{print \$NF}')\n";
          }
        }
      }
      $self->close_pbs( $pbs, $pbs_file );

      print $sh "bash ./$pbs_name \n";
    }
    print $sh "exit 0\n";
    close $sh;
    if ( is_linux() ) {
      chmod 0755, $shfile;
    }
    for my $clear_key ( sort keys %$clear_keys ) {
      my $pbs_file = $self->get_step_sample_pbs( $pbs_dir, $step_name, $clear_key );
      my $clear_file = change_extension( $pbs_file, "_clear.sh" );
      open my $clear, ">$clear_file" or die "Cannot create file $clear_file";
      print $clear "
read -p \"Are you sure you want to clear all result for $clear_key? \" -n 1 -r
echo
if [[ \$REPLY =~ ^[Yy]\$ ]]
then
";
      for my $task_section (@tasks) {

        my $clear_files = $clears->{$task_section}{$clear_key};
        if ( defined $clear_files ) {
          for my $expect_file (@$clear_files) {
            print $clear "  rm -rf $expect_file \n";
          }
        }
      }

      print $clear "fi \n";
      close($clear);
    }
  }
  close($result_list);

  print $final "\nbash $report_pbs \n";
  $self->close_pbs( $final, $final_pbs );
}

1;
