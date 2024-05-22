#!/usr/bin/perl
package Pipeline::SLAMseq;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performSLAMseq)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing",   1 );
  initDefaultValue( $def, "DE_export_significant_gene_name",   1 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [@$individual, @$summary];

  my $target_dir      = $def->{target_dir};

  my $reference_fasta = getValue($def, "fasta_file");
  my $utr3_bed = getValue($def, "utr3_bed");
  my $max_thread = getValue($def, "max_thread", 8);

  #https://github.com/t-neumann/slamdunk/issues/25#issuecomment-385730642
  #--max-polya 4 is required to generate XA tag, otherwise, count will raise error
  #since filter command will sort the bam file, we don't need to sort bam file in ngm mapping
  my $nextgenmap_task = "T01_nextgenmap";
  $config->{$nextgenmap_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$nextgenmap_task",
    init_command          => "",
    option                => "
echo ngm=`date`
ngm -r $reference_fasta \\
  -1 __FILE__ \\
  -t $max_thread \\
  --no-progress \\
  --slam-seq 2 \\
  -5 12 \\
  -l \\
  --max-polya 4 \\
  --rg-id 1 \\
  --rg-sm __NAME__:NA:-1 | samtools view -bS -o __NAME__.tmp.bam -
  
status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.ngm.failed
  rm -f __NAME__.tmp.bam __NAME__.ngm.succeed
else
  touch __NAME__.ngm.succeed
  rm -f __NAME__.ngm.failed
  mv __NAME__.tmp.bam __NAME__.bam
fi
",
    docker_prefix         => "nextgenmap_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $source_ref,
    source_arg            => "",
    source_join_delimiter => " -2 ",
    output_to_same_folder => 0,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => ".bam",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $nextgenmap_task);

  my $slamdunk_filter_task = "T02_slamdunk_filter";
  $config->{$slamdunk_filter_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$slamdunk_filter_task",
    init_command          => "",
    option                => "
slamdunk filter -o . \\
  --threads $max_thread \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.filter.failed
  rm -f __NAME___filtered.bam __NAME__.filter.succeed
else
  touch __NAME__.filter.succeed
  rm -f __NAME__.filter.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $nextgenmap_task,
    no_bai => 1,
    source_arg            => "",
    source_join_delimiter => "",
    no_output => 1,
    output_to_same_folder => 1,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => "_filtered.bam",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $slamdunk_filter_task);

  my $slamdunk_snp_task = "T03_slamdunk_snp";
  $config->{$slamdunk_snp_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$slamdunk_snp_task",
    init_command          => "",
    option                => "
slamdunk snp -o . \\
  -t $max_thread \\
  -r $reference_fasta \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.snp.failed
  rm -f __NAME__.snp.succeed
else
  touch __NAME__.snp.succeed
  rm -f __NAME__.snp.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_filter_task,
    source_arg            => "",
    source_join_delimiter => "",
    no_output => 1,
    output_to_same_folder => 1,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => "_filtered_snp.vcf",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $slamdunk_snp_task);

  my $read_length = getValue($def, "read_length");
  my $slamdunk_count_task = "T04_slamdunk_count";
  $config->{$slamdunk_count_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$slamdunk_count_task",
    init_command          => "",
    option                => "
slamdunk count -o . \\
  -t $max_thread \\
  -l $read_length \\
  -v __FILE2__ \\
  -b $utr3_bed \\
  -r $reference_fasta \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.count.failed
  rm -f __NAME__.count.succeed
else
  touch __NAME__.count.succeed
  rm -f __NAME__.count.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_filter_task,
    source_arg            => "",
    source_join_delimiter => "",
    parameterSampleFile2_ref => $slamdunk_snp_task,
    no_output => 1,
    output_to_same_folder => 1,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => "_filtered_tcount.tsv",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $slamdunk_count_task);

  #https://github.com/nf-core/slamseq/blob/1.0.0/main.nf
  my $alleyoop_collapse_task = "T05_alleyoop_collapse";
  $config->{$alleyoop_collapse_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$alleyoop_collapse_task",
    init_command          => "",
    option                => "
alleyoop collapse \\
    -o . \\
    -t $max_thread \\
    __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.collapse.failed
  rm -f __NAME__.collapse.succeed
else
  touch __NAME__.collapse.succeed
  rm -f __NAME__.collapse.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_count_task,
    source_arg            => "",
    source_join_delimiter => "",
    no_output => 1,
    output_to_same_folder => 1,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => "_filtered_tcount_collapsed.csv",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $alleyoop_collapse_task);

  my $alleyoop_rates_task = "T06_alleyoop_rates";
  $config->{$alleyoop_rates_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$alleyoop_rates_task",
    init_command          => "",
    option                => "
alleyoop rates \\
  -o . \\
  -r $reference_fasta \\
  -t $max_thread \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.rates.failed
  rm -f __NAME__.rates.succeed
else
  touch __NAME__.rates.succeed
  rm -f __NAME__.rates.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_filter_task,
    source_arg            => "",
    source_join_delimiter => "",
    no_output => 1,
    output_to_same_folder => 1,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => "_filtered_overallrates.csv",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $alleyoop_rates_task);

  my $alleyoop_utrrates_task = "T07_alleyoop_utrrates";
  $config->{$alleyoop_utrrates_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$alleyoop_utrrates_task",
    init_command          => "",
    option                => "
alleyoop utrrates \\
  -o . \\
  -r $reference_fasta \\
  -b $utr3_bed \\
  -l $read_length \\
  -t $max_thread \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.utrrates.failed
  rm -f __NAME__.utrrates.succeed
else
  touch __NAME__.utrrates.succeed
  rm -f __NAME__.utrrates.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_filter_task,
    source_arg            => "",
    source_join_delimiter => "",
    no_output => 1,
    output_to_same_folder => 1,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => "_filtered_mutationrates_utr.csv",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $alleyoop_utrrates_task);

  my $alleyoop_tcperreadpos_task = "T08_alleyoop_tcperreadpos";
  $config->{$alleyoop_tcperreadpos_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$alleyoop_tcperreadpos_task",
    init_command          => "",
    option                => "
alleyoop tcperreadpos \\
  -o . \\
  -r $reference_fasta \\
  -v __FILE2__ \\
  -l $read_length \\
  -t $max_thread \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.tcperreadpos.failed
  rm -f __NAME__.tcperreadpos.succeed
else
  touch __NAME__.tcperreadpos.succeed
  rm -f __NAME__.tcperreadpos.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_filter_task,
    source_arg            => "",
    source_join_delimiter => "",
    parameterSampleFile2_ref => $slamdunk_snp_task,
    no_output => 1,
    output_to_same_folder => 1,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => "_filtered_tcperreadpos.csv",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $alleyoop_tcperreadpos_task);

  my $alleyoop_tcperutrpos_task = "T09_alleyoop_tcperutrpos";
  $config->{$alleyoop_tcperutrpos_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$alleyoop_tcperutrpos_task",
    init_command          => "",
    option                => "
alleyoop tcperutrpos \\
  -o . \\
  -r $reference_fasta \\
  -b $utr3_bed \\
  -v __FILE2__ \\
  -l $read_length \\
  -t $max_thread \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.tcperutrpos.failed
  rm -f __NAME__.tcperutrpos.succeed
else
  touch __NAME__.tcperutrpos.succeed
  rm -f __NAME__.tcperutrpos.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_filter_task,
    source_arg            => "",
    source_join_delimiter => "",
    parameterSampleFile2_ref => $slamdunk_snp_task,
    no_output => 1,
    output_to_same_folder => 1,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => "_filtered_tcperutr.csv",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $alleyoop_tcperutrpos_task);

  my $alleyoop_summary_task = "T10_alleyoop_summary";
  my $count_folder = $config->{$slamdunk_count_task}{target_dir} . "/result";
  $config->{$alleyoop_summary_task} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "$target_dir/$alleyoop_summary_task",
    init_command          => "",
    option                => "
alleyoop summary \\
  -o __NAME__.summary.txt \\
  -t $count_folder \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.summary.failed
  rm -f __NAME__.summary.succeed
else
  touch __NAME__.summary.succeed
  rm -f __NAME__.summary.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_filter_task,
    source_arg            => "",
    source_type           => "array",
    source_join_delimiter => " ",
    parameterSampleFile2_ref => $slamdunk_count_task,
    no_output => 1,
    output_to_same_folder => 1,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext => ".summary.txt",
    output_other_ext => ".summary_PCA.txt",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $alleyoop_summary_task);

  my $count_table_task = "T11_count_table";
  $config->{$count_table_task} = {
    class                 => "CQS::UniqueR",
    perform               => 1,
    target_dir            => "$target_dir/$count_table_task",
    init_command          => "",
    rtemplate => "../CQS/countTableVisFunctions.R,../SLAMseq/slamseq_count_table.r",
    parameterSampleFile1_ref => $alleyoop_collapse_task,
    parameterFile1 => getValue($def, "transcript_gtf"),
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext => ".tc_read.csv",
    output_other_ext => ".all_read.csv,.sizeFactor.csv",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $count_table_task);

  if(defined $def->{pairs}){
    my $deseq2taskname = addDEseq2( $config, $def, $tasks, "tc_read", [ $count_table_task, ".tc_read.csv" ], $target_dir, 1, [$count_table_task, ".sizeFactor.csv"], "sizeFactor"  );

    if ( getValue( $def, "perform_webgestalt" ) ) {
      my $webgestaltTaskName = addWebgestalt($config, $def, $tasks, $target_dir, $deseq2taskname, [ $deseq2taskname, "sig_genename.txt\$" ]);

      #if ( defined $def->{perform_link_webgestalt_deseq2} ) {
      my $linkTaskName = $webgestaltTaskName . "_link_deseq2";
      $config->{$linkTaskName} = {
        class                      => "CQS::UniqueR",
        perform                    => 1,
        target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $linkTaskName,
        rtemplate                  => "../Annotation/WebGestaltReportFunctions.r;../Annotation/WebGestaltDeseq2.r",
        rReportTemplate            => "../Annotation/WebGestaltDeseq2.rmd",
        output_to_result_directory => 1,
        output_perSample_file      => "parameterSampleFile1",
        output_perSample_file_regex => "enrichment_results_(.+).txt",
        output_perSample_file_ext  => ".html;.html.rds",
        parameterSampleFile1_ref   => [ $webgestaltTaskName, ".txt\$" ],
        parameterSampleFile2_ref   => [ $deseq2taskname, "sig.csv\$" ],
        sh_direct                  => 1,
        rCode                      => "",
        pbs                        => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "23",
          "mem"       => "10gb"
        },
      };
      push( @$tasks, $linkTaskName );

      if (getValue( $def, "perform_webgestaltHeatmap", 0 )) {
        my $webgestaltHeatmapTaskName = $webgestaltTaskName . "_heatmap_deseq2";
        $config->{$webgestaltHeatmapTaskName} = {
          class                      => "CQS::UniqueR",
          perform                    => 1,
          target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $webgestaltHeatmapTaskName,
          rtemplate                  => "../Annotation/WebGestaltReportFunctions.r;../Annotation/WebGestaltHeatmap.r",
          output_to_result_directory => 1,
          output_perSample_file      => "parameterSampleFile1",
          output_perSample_file_ext  => ".heatmap.png",
          parameterSampleFile1_ref   => [ $webgestaltTaskName, ".txt\$" ],
          parameterSampleFile2_ref   => [ $deseq2taskname, "sig.csv\$" ],
          parameterSampleFile3_ref   => [ $deseq2taskname, "vsd.csv\$" ],
          parameterSampleFile4_ref   => [ $deseq2taskname, ".design\$" ],
          sh_direct                  => 1,
          rCode                      => "",
          pbs                        => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "23",
            "mem"       => "10gb"
          },
        };
        push( @$tasks, $webgestaltHeatmapTaskName );
      }
      #}
    }

    #print(Dumper($def));

    if ( getValue( $def, "perform_gsea" ) ) {
      my $gseaTaskName = $deseq2taskname . "_GSEA";

      if(getValue($def, "use_mouse_gsea_db", 0)){
        $gseaTaskName = $gseaTaskName . "_Mm";
      }else{
        $gseaTaskName = $gseaTaskName . "_Hs";
      }

      my $pairs = $config->{pairs};
      my $keys = [keys %$pairs];

      add_gsea($config, $def, $tasks, $target_dir, $gseaTaskName, [ $deseq2taskname, "_GSEA.rnk\$" ], $keys, "" );
    }
  }

  $config->{report} = {
    class                      => "CQS::BuildReport",
    perform                    => 1,
    target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . "report",
    report_rmd_file            => "../Pipeline/SLAMseq.Rmd",
    additional_rmd_files       => "../Pipeline/Pipeline.R;reportFunctions.R",
    parameterSampleFile1 => {
      task_name => $task_name,
      email     => $email,
      affiliation => getValue($def, "affiliation", ""),
    },
    parameterSampleFile2 => $def->{groups},
    parameterSampleFile4_ref => $alleyoop_summary_task,
    parameterSampleFile5_ref => $alleyoop_rates_task,
    parameterSampleFile6_ref => $alleyoop_utrrates_task,
    parameterSampleFile7_ref => $alleyoop_tcperreadpos_task,
    parameterSampleFile8_ref => $alleyoop_tcperutrpos_task,
    sh_direct => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$tasks, "report" );

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      tasks => $tasks,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
    },
  };

  return ($config);
}

sub performSLAMseq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

1;
