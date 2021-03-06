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

sub getAllDependentJobids {
  my ($final_pbs_id_map) = @_;

  my $result = "";
  if ( keys %$final_pbs_id_map ) {
    $result = "--dependency=afterany";
    for my $each_dep_pbs ( keys %$final_pbs_id_map ) {
      if ($final_pbs_id_map->{$each_dep_pbs}[1] == 0){
        $result = $result . ":\$" . $final_pbs_id_map->{$each_dep_pbs}[0];
      }
    }
  }

  return $result;
}

sub getDependentJobids {
  my ( $task_dep_pbs_map, $pbs_id_map, $task_section, $task_name, $sample ) = @_;
  my $dep_pbs_map = $task_dep_pbs_map->{$task_section};

  my $result = "";
  if ( $task_name eq $sample ) {
    if ( keys %$dep_pbs_map ) {
      $result = "--dependency=afterany";
      for my $sample ( keys %$dep_pbs_map ) {
        my $sample_pbs_map = $dep_pbs_map->{$sample};
        if ( keys %$sample_pbs_map ) {
          for my $each_dep_pbs ( keys %$sample_pbs_map ) {
            my $pid = $pbs_id_map->{$each_dep_pbs};
            if(!defined $pid){
              die "Undefined $each_dep_pbs";
            }
            $result = $result . ":\$" . $pid->[0];
            $pbs_id_map->{$each_dep_pbs}[1] = 1;
          }
        }
      }
    }
  }
  else {
    my $taskname_pbs_map = $dep_pbs_map->{$task_name};
    my $dep_pbs          = $dep_pbs_map->{$sample};

    if ( defined $dep_pbs ) {
      $result = "--dependency=afterany";
      for my $each_dep_pbs ( keys %$dep_pbs ) {
        $result = $result . ":\$" . $pbs_id_map->{$each_dep_pbs}[0];
        $pbs_id_map->{$each_dep_pbs}[1] = 1;
      }
    }
    if ( defined $taskname_pbs_map ) {
      if ( $result eq "" ) {
        $result = "--dependency=afterany";
      }
      for my $each_dep_pbs ( keys %$taskname_pbs_map ) {
        if (!defined $pbs_id_map->{$each_dep_pbs} ){
          die "Didn't find $each_dep_pbs in $pbs_id_map";
        }
        $result = $result . ":\$" . $pbs_id_map->{$each_dep_pbs}[0];
        $pbs_id_map->{$each_dep_pbs}[1] = 1;
      }
    }
  }

  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section );

  my $task_dep_pbs_map = $self->get_all_dependent_pbs_map( $config, $section );

  my %step_map = %{ get_raw_files( $config, $section ) };

  my $task_shell = get_option( $config, $section, "task_shell", "bash" );

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

  open( my $final_submit, ">${final_pbs}.submit" ) or die $!;
  my $final = $self->open_pbs( $final_pbs, $pbs_desc, $final_log_desp, $path_file, $pbs_dir );

  open( my $result_list, ">$result_list_file" ) or die $!;
  print $result_list "StepName\tTaskName\tSampleName\tFileList\tCanFileEmpty\n";

  my $pbs_id_map = {};

  my $final_pbs_id_map = {};
  my $final_index      = 0;
  for my $step_name ( sort keys %step_map ) {
    my $shfile = $self->get_file( $pbs_dir, $step_name, ".submit" );
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
      eval { $expect_file_map = $myclass->result( $config, $task_section, "(?<!version)\$", 1 ); } or do {
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
        if ( !-e $samplepbs ) {
          die "Task " . $task_section . ", file not exists " . $samplepbs . "\n";
        }

        my $depjid = getDependentJobids( $task_dep_pbs_map, $final_pbs_id_map, $task_section, $task_name, $sample );

        $final_index = $final_index + 1;
        
        if(not defined $expect_file_map->{$sample}){
          print $final_submit "jid" . $final_index . "=\$(sbatch $depjid " . $samplepbs . " | awk '{print \$NF}') \n";
          $final_pbs_id_map->{$samplepbs} = ["jid" . $final_index, 0];
        }else{
          my $expect_file = $expect_file_map->{$sample}[0];
          print $final_submit "if [[ (1 -eq \$1) || ((! -s $expect_file) && (! -d $expect_file)) ]]; then \n";
          print $final_submit "  jid" . $final_index . "=\$(sbatch $depjid " . $samplepbs . " | awk '{print \$NF}') \n";
          print $final_submit "else \n";
          print $final_submit "  jid" . $final_index . "=1000000000 \n";
          print $final_submit "fi \n";
          $final_pbs_id_map->{$samplepbs} = ["jid" . $final_index, 0];
        }
        print $final "bash $samplepbs \n";
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
      my $submit_file = "${pbs_file}.submit";
      open( my $pbs_submit, ">$submit_file" ) or die $!;

      my $index = 0;
      for my $task_section (@tasks) {

        #print "task " . $task_section . " ...\n";
        my $pbs_map          = $taskpbs->{$task_section};
        my $dep_pbs_map      = $task_dep_pbs_map->{$task_section};
        my $taskname_pbs_map = $dep_pbs_map->{$task_name};
        if ( defined $pbs_map->{$sample} ) {
          my $depjid = getDependentJobids( $task_dep_pbs_map, $pbs_id_map, $task_section, $task_name, $sample );
          my $samplepbs = $pbs_map->{$sample};

          $index = $index + 1;
          print $pbs_submit "jid" . $index . "=\$(sbatch $depjid " . $samplepbs . " | awk '{print \$NF}')\n";
          $pbs_id_map->{$samplepbs} = ["jid" . $index, 0];

          print $pbs "bash " . $samplepbs . " \n";
        }
      }
      close($pbs_submit);
      if ( is_linux() ) {
        chmod 0755, $submit_file;
      }
      
      $self->close_pbs( $pbs, $pbs_file );

      print $sh "bash ./${pbs_name}.submit \n";
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
      if ( is_linux() ) {
        chmod 0755, $clear_file;
      }
    }
  }
  close($result_list);

  my $alldepids = getAllDependentJobids($final_pbs_id_map);
  $final_index = $final_index + 1;
  print $final_submit "jid" . $final_index . "=\$(sbatch $alldepids " . $report_pbs . " | awk '{print \$NF}')\n";
  close($final_submit);
  if ( is_linux() ) {
    chmod 0755, "${final_pbs}.submit";
  }

  print $final "\nbash $report_pbs \n";
  $self->close_pbs( $final, $final_pbs );
}

1;
