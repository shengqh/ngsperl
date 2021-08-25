#!/usr/bin/perl
package CQS::SequenceTaskSlurmSlim;

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
      if ( $final_pbs_id_map->{$each_dep_pbs}[1] == 0 ) {
        $result = $result . ":\$" . $final_pbs_id_map->{$each_dep_pbs}[0];
      }
    }
  }

  return $result;
}

sub get_dependent_job_ids {
  my ( $task_dep_pbs_map, $pbs_id_map, $task_section, $task_name, $pbs_file ) = @_;
  my $dep_pbs_map = $task_dep_pbs_map->{$task_section}{$pbs_file};

  my $result = "";
  if(not defined $dep_pbs_map){
    return($result);
  }

  for my $dep_pbs (keys %$dep_pbs_map){
    if (!defined $pbs_id_map->{$dep_pbs}){
      my $str = Dumper($pbs_id_map);
      open(my $idmap, ">id.map");
      print ($idmap $str);
      close($idmap);
      #die "$dep_pbs is not in pbs_id_map: " . Dumper($pbs_id_map);
      die "processing $task_section, $dep_pbs is not in pbs_id_map: id.map " ;
    }
    
    if ($result eq ""){
      $result = "--dependency=afterany";
    }
    
    $result = $result . ":\$" . $pbs_id_map->{$dep_pbs}[0];
    $pbs_id_map->{$dep_pbs}[1] = 1;
  }

  return $result;
}

sub getDependentJobids {
  my ( $task_dep_pbs_map, $pbs_id_map, $depend_all, $task_section, $task_name, $dep_sample_names ) = @_;
  my $dep_pbs_map = $task_dep_pbs_map->{$task_section};

  #if ($task_section eq "fastqc_post_trim"){
  #  print Dumper($dep_pbs_map);
  #  print @$dep_sample_names;
  #}

  my $result = "";
  if (not keys %$dep_pbs_map){
    return($result);
  }

  if ($depend_all){
    $result = "--dependency=afterany";
    for my $sample ( keys %$dep_pbs_map ) {
      my $sample_pbs_map = $dep_pbs_map->{$sample};
      if ( keys %$sample_pbs_map ) {
        for my $each_dep_pbs ( keys %$sample_pbs_map ) {
          my $pid = $pbs_id_map->{$each_dep_pbs};
          if ( !defined $pid ) {
            die "Undefined $each_dep_pbs";
          }
          $result = $result . ":\$" . $pid->[0];
          $pbs_id_map->{$each_dep_pbs}[1] = 1;
        }
      }
    }
    return($result);
  }

  my $taskname_pbs_map = $dep_pbs_map->{$task_name};
  if ( defined $taskname_pbs_map ) {
    $result = "--dependency=afterany";
    for my $each_dep_pbs ( keys %$taskname_pbs_map ) {
      if ( !defined $pbs_id_map->{$each_dep_pbs} ) {
        die "Didn't find $each_dep_pbs in $pbs_id_map";
      }
      $result = $result . ":\$" . $pbs_id_map->{$each_dep_pbs}[0];
      $pbs_id_map->{$each_dep_pbs}[1] = 1;
    }
  }

  for my $dep_sample_name (@$dep_sample_names){
    my $dep_pbs = $dep_pbs_map->{$dep_sample_name};
    if ( defined $dep_pbs ) {
      if ($result eq ""){
        $result = "--dependency=afterany";
      }
      for my $each_dep_pbs ( keys %$dep_pbs ) {
        if (!defined $pbs_id_map->{$each_dep_pbs}){
          die "$each_dep_pbs is not in pbs_id_map: " . Dumper($pbs_id_map);
        }
        $result = $result . ":\$" . $pbs_id_map->{$each_dep_pbs}[0];
        $pbs_id_map->{$each_dep_pbs}[1] = 1;
      }
    }
  }

  #print(@$dep_sample_names);
  #print($result);
  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section );

  my $task_dep_pbs_map = $self->get_all_dependent_pbs_map( $config, $section );

  # if ($config->{general}{debug}){
  #   print Dumper($task_dep_pbs_map->{bowtie1_contamination_all_02_align});
  # }

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

  my $result_list_file = $self->get_file( $result_dir, $task_name, "_expect_result.tsv" );

  my $report = $self->open_pbs( $report_pbs, $pbs_desc, $report_log_desp, $path_file, $result_dir );

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
        print(Dumper(@tasks));
        die "$task_section is not a valid task section.";
      }
      
      my $cur_task_dep_pbs_map = $task_dep_pbs_map->{$task_section};
      
      my $myclass = instantiate($classname);

      my $depend_all = $myclass->{_depend_all};

      #print($classname . " = " . $depend_all . "\n");

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

        my @canEmpty = ();
        for my $expect_file (@$expect_files){
          push(@canEmpty, $myclass->can_result_be_empty_file( $config, $task_section, $expect_file) ? "True" : "False");
        }
        my $canEmptyStr = join(",", @canEmpty);
        print $result_list $step_name, "\t", $task_section, "\t", $expect_name, "\t", $expect_file_list, "\t", $canEmptyStr, "\n";
      }

      for my $sample ( sort keys %{$pbs_file_map} ) {
        my $samplepbs = $pbs_file_map->{$sample};
        if ( !-e $samplepbs ) {
          die "Task " . $task_section . ", file not exists " . $samplepbs . "\n";
        }

        my $depjid = get_dependent_job_ids( $task_dep_pbs_map, $final_pbs_id_map, $task_section, $task_name, $samplepbs );

        $final_index = $final_index + 1;

        if ( not defined $expect_file_map->{$sample} ) {
          print $final_submit "jid" . $final_index . "=\$(sbatch $depjid " . $samplepbs . " | awk '{print \$NF}') \n";
          $final_pbs_id_map->{$samplepbs} = [ "jid" . $final_index, 0 ];
        }
        else {
          my $expect_file = $expect_file_map->{$sample}[0];
          
          if(not defined $expect_file){
            print "Found error";
          }
          print $final_submit "if [[ (1 -eq \$1) || ((! -s $expect_file) && (! -d $expect_file)) ]]; then \n";
          print $final_submit "  jid" . $final_index . "=\$(sbatch $depjid " . $samplepbs . " | awk '{print \$NF}') \n";
          print $final_submit "else \n";
          print $final_submit "  jid" . $final_index . "=1000000000 \n";
          print $final_submit "fi \n";
          $final_pbs_id_map->{$samplepbs} = [ "jid" . $final_index, 0 ];
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

  }
  close($result_list);

  #print(Dumper(%$final_pbs_id_map));
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
