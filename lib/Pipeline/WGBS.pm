#!/usr/bin/perl
package Pipeline::WGBS;

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
use Pipeline::MethylationUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performWGBS)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType",              "FAIL" );
  initDefaultValue( $def, "cluster",                "slurm" );
  initDefaultValue( $def, "perform_preprocessing",  1 );
  initDefaultValue( $def, "perform_age_estimation", 0 );
  initDefaultValue( $def, "methylation_mincov",     "20" );
  initDefaultValue( $def, "perform_multiqc",        1 );

  return $def;
} ## end sub initializeDefaultOptions


sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  #$def->{perform_cutadapt} = 0;

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);
  my $tasks = [ @$individual, @$summary ];

  my $target_dir = $def->{target_dir};

  if ( $def->{perform_preprocessing_only} ) {

    $config->{sequencetask} = {
      class      => getSequenceTaskClassname($cluster),
      perform    => 1,
      target_dir => "${target_dir}/sequencetask",
      option     => "",
      source     => { tasks => $tasks, },
      sh_direct  => 0,
      cluster    => $cluster,
      pbs        => {
        "nodes"    => "1:ppn=" . $def->{max_thread},
        "walltime" => $def->{sequencetask_run_time},
        "mem"      => "40gb"
      },
    };
    return $config;
  } ## end if ( $def->{perform_preprocessing_only...})

  my $thread             = getValue( $def, "thread", 4 );
  my $abismal_index      = getValue( $def, "abismal_index" );
  my $chr_fasta          = getValue( $def, "chr_fasta" );
  my $chr_size_file      = getValue( $def, "chr_size_file" );
  my $annovar_buildver   = getValue( $def, "annovar_buildver" );
  my $annovar_db         = getValue( $def, "annovar_db" );
  my $annovar_param      = getValue( $def, "annovar_param" );
  my $HOMER_perlFile     = getValue( $def, "HOMER_perlFile" );
  my $addqual_pythonFile = dirname(__FILE__) . "/../Methylation/add_qual.py";
  my $picard             = getValue( $def, "picard" );
  my $gatk               = getValue( $def, "gatk" );
  my $interval_list      = getValue( $def, "interval_list" );

  # my $trimgalore_task = "trimgalore";
  # $config->{$trimgalore_task} = {
  #   class => "Trimmer::TrimGalore",
  #   perform => 1,
  #   target_dir => "${target_dir}/" . getNextFolderIndex($def) . "trimgalore",
  #   option => getValue($def, "trimgalore_option", ""),
  #   source_ref => $untrimed_ref,
  #   extension => "_val.fq.gz",
  #   pairend => is_paired_end($def),
  #   do_fastqc => getValue($def, "trimgalore_do_fastqc", 1),
  #   init_command => "",
  #   pbs        => {
  #     "nodes"     => "1:ppn=1",
  #     "walltime"  => "48",
  #     "mem"       => "10gb"
  #   },
  # };
  # push(@$tasks, $trimgalore_task);
  my $methylation_key = "methylation_key";

  my $abismal_task  = "abismal";
  my $abismal_class = $def->{abismal_not_filter_intervals} ? "Alignment::Abismal" : "Alignment::AbismalFilterIntervals";
  $config->{$abismal_task} = {
    class              => $abismal_class,
    perform            => 1,
    option             => "",
    thread             => $thread,
    target_dir         => "${target_dir}/" . getNextFolderIndex($def) . "abismal",
    chr_fasta          => $chr_fasta,
    picard             => $picard,
    gatk               => $gatk,
    interval_list      => $interval_list,
    addqual_pythonFile => $addqual_pythonFile,
    abismal_index      => $abismal_index,
    source_ref         => $source_ref,
    dnmtools_command   => getValue( $def, "dnmtools_command", "dnmtools" ),
    preseq_command     => getValue( $def, "preseq_command",   "preseq" ),
    docker_prefix      => "dnmtools_",
    pbs                => {
      "nodes"    => "1:ppn=8",
      "walltime" => getValue( $def, "abismal_walltime", "128" ),
      "mem"      => "80gb"
    },
  };
  push( @$tasks, $abismal_task );

  my $abismal_summary_task = "abismal_summary";
  $config->{$abismal_summary_task} = {
    class                    => "CQS::UniqueRmd",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "$abismal_summary_task",
    report_rmd_file          => "../Alignment/AbismalSummary.rmd",
    additional_rmd_files     => "../CQS/reportFunctions.R;../CQS/countTableVisFunctions.R",
    option                   => "",
    parameterSampleFile1_ref => [$abismal_task],
    parameterSampleFile2     => {
      task_name   => getValue( $def, "task_name" ),
      email       => getValue( $def, "email" ),
      affiliation => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
    },
    suffix                   => ".abismal",
    output_file_ext          => ".abismal.html;.reads_mapped.bar.png;.reads_off_bait.bar.png;.final_unique_reads.bar.png;.duplication_rate.png;.rrbs_rate.png",
    can_result_be_empty_file => 0,
    sh_direct                => 1,
    no_docker                => 0,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "40gb"
    },
  };
  push( @$tasks, $abismal_summary_task );

  my $dnmtools_task = "DNMTools";
  $config->{$dnmtools_task} = {
    class            => "Methylation::DNMTools",
    perform          => 1,
    target_dir       => "${target_dir}/" . getNextFolderIndex($def) . "dnmtools",
    option           => "",
    chr_fasta        => $chr_fasta,
    chr_size_file    => $chr_size_file,
    source_ref       => [ "abismal", ".intervals.dnmtools_format.uniq.addqual.bam\$" ],
    dnmtools_command => getValue( $def, "dnmtools_command", "dnmtools" ),
    docker_prefix    => "dnmtools_",
    pbs              => {
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "80gb"
    },
  };
  push( @$tasks, $dnmtools_task );

  my $dnmtoolsannovar_task = "dnmtoolsAnnovar";
  $config->{$dnmtoolsannovar_task} = {
    class            => "Annotation::Annovar",
    perform          => 1,
    target_dir       => "${target_dir}/" . getNextFolderIndex($def) . "dnmtoolsAnnovar",
    option           => $annovar_param,
    annovar_db       => $annovar_db,
    buildver         => $annovar_buildver,
    perform_splicing => 0,
    docker_prefix    => "annovar_",
    isBed            => 1,
    source_ref       => [ "DNMTools", ".hmr\$|.amr\$|.pmd\$" ],
    pbs              => {
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  };
  push( @$tasks, $dnmtoolsannovar_task );

  my $methylkitprep_task = "MethylKitPreparation";
  $config->{$methylkitprep_task} = {
    class                    => "CQS::IndividualR",
    perform                  => 1,
    target_dir               => "${target_dir}/" . getNextFolderIndex($def) . "MethylKitPreparation",
    option                   => "",
    docker_prefix            => "wgbs_r_",
    rtemplate                => "../Methylation/prepare_CpG_input.R",
    parameterSampleFile1_ref => [ "DNMTools", ".cpg.meth.gz\$" ],
    output_file_ext          => ".CpG.txt.gz",
    sh_direct                => 1,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "100gb"
    },
  };
  push( @$tasks, $methylkitprep_task );

  my $methylkitcorr_task = add_MethylKitCorr( $config, $def, $tasks, $target_dir, "MethylKitCorr", [ $methylkitprep_task, ".CpG.txt.gz" ], "amp" );

  if ( $def->{perform_age_estimation} ) {
    my $methy_age_task = add_MethylAgeEstimation( $config, $def, $tasks, $target_dir, "dnaMethyAge", $methylkitcorr_task );
  }

  my $methylkitdiff_task             = undef;
  my $methylkitdiffannovar_task      = undef;
  my $MethylKitDiffAnnovarGenes_task = undef;
  my $webgestalt_task                = undef;

  if ( defined $def->{pairs} ) {
    my $task_map = add_MethylDiffAnalysis( $config, $def, $tasks, $target_dir, $methylkitcorr_task );
    $methylkitdiff_task             = $task_map->{methylkitdiff_task};
    $methylkitdiffannovar_task      = $task_map->{methylkitdiffannovar_task};
    $MethylKitDiffAnnovarGenes_task = $task_map->{MethylKitDiffAnnovarGenes_task};
    $webgestalt_task                = $task_map->{webgestalt_task};

    #  my $homer_task = "HOMER_DMR";
    #  $config->{$homer_task} = {
    #    class        => "Homer::FindMotifs",
    #    perform      => 1,
    #    target_dir   => "${target_dir}/HOMER_DMR",
    #    option       => "-nomotif",
    #    homer_genome => "hg19",
    #    source_ref   => [ "DNMToolsDiff", ".DMR.filtered\$" ],
    #    remove_empty_source => 1,
    #    sh_direct    => 0,
    #    pbs          => {
    #      "nodes"     => "1:ppn=1",
    #      "walltime"  => "2",
    #      "mem"       => "40gb"
    #    },
    #  };
    #  push(@$tasks, "HOMER_DMR");
  } ## end if ( defined $def->{pairs...})

  my $version_files = get_version_files($config);

  my @report_files = ();
  my @report_names = ();
  my @copy_files   = ();
  if ( defined $config->{fastqc_raw_summary} ) {
    push( @report_files, "fastqc_raw_summary", ".FastQC.baseQuality.tsv.png" );
    push( @report_files, "fastqc_raw_summary", ".FastQC.sequenceGC.tsv.png" );
    push( @report_files, "fastqc_raw_summary", ".FastQC.adapter.tsv.png" );
    push( @report_names, "fastqc_raw_per_base_sequence_quality", "fastqc_raw_per_sequence_gc_content", "fastqc_raw_adapter_content" );
  } ## end if ( defined $config->...)

  if ( defined $config->{fastqc_post_trim_summary} ) {
    push( @report_files, "fastqc_post_trim_summary", ".FastQC.baseQuality.tsv.png" );
    push( @report_files, "fastqc_post_trim_summary", ".FastQC.sequenceGC.tsv.png" );
    push( @report_files, "fastqc_post_trim_summary", ".FastQC.adapter.tsv.png" );
    push( @report_names, "fastqc_post_trim_per_base_sequence_quality", "fastqc_post_trim_per_sequence_gc_content", "fastqc_post_trim_adapter_content" );
  } ## end if ( defined $config->...)

  if ( ( defined $dnmtools_task ) && ( defined $config->{$dnmtools_task} ) ) {
    push( @copy_files, $dnmtools_task, ".amr\$" );
    push( @copy_files, $dnmtools_task, ".hmr\$" );
    push( @copy_files, $dnmtools_task, ".pmd\$" );
  }
  if ( ( defined $dnmtoolsannovar_task ) && ( defined $config->{$dnmtoolsannovar_task} ) ) {
    push( @copy_files, $dnmtoolsannovar_task, ".annovar.final.tsv\$" );
  }
  if ( ( defined $methylkitcorr_task ) && ( defined $config->{$methylkitcorr_task} ) ) {
    push( @copy_files, $methylkitcorr_task, ".png\$|.rds\$" );
  }
  if ( ( defined $methylkitdiff_task ) && ( defined $config->{$methylkitdiff_task} ) ) {
    push( @copy_files, $methylkitdiff_task, ".dmcpgs.tsv\$" );
  }
  if ( ( defined $methylkitdiffannovar_task ) && ( defined $config->{$methylkitdiffannovar_task} ) ) {
    push( @copy_files, $methylkitdiffannovar_task, ".annovar.final.tsv\$" );
  }
  # if ( defined $webgestalt_task ) {
  #   push( @copy_files, $webgestalt_task, "_geneontology_Biological_Process\$" );
  #   push( @copy_files, $webgestalt_task, "_geneontology_Cellular_Component\$" );
  #   push( @copy_files, $webgestalt_task, "_geneontology_Molecular_Function\$" );
  #   push( @copy_files, $webgestalt_task, "_pathway_KEGG\$" );
  # } ## end if ( defined $webgestalt_task)

  if ( $def->{perform_multiqc} ) {
    addMultiQC( $config, $def, $tasks, $target_dir, $target_dir, $dnmtools_task );
  }

  $config->{report} = {
    class      => "CQS::BuildReport",
    perform    => 1,
    target_dir => "$target_dir/report",
    #docker_command             => "singularity exec -c -e -B /home,/gpfs51,/gpfs52,/panfs,/data,/dors,/nobackup,/tmp -H `pwd` /data/cqs/softwares/singularity/wgbs_r.1.1.sif ",
    docker_prefix              => "wgbs_r_",
    report_rmd_file            => "../Pipeline/WGBS.Rmd",
    additional_rmd_files       => "../Pipeline/Pipeline.R;reportFunctions.R",
    parameterSampleFile1_ref   => \@report_files,
    parameterSampleFile1_names => \@report_names,
    parameterSampleFile2       => {
      task_name          => $task_name,
      meta_data          => "../../" . $task_name . "_meta.tsv",
      abismal_path       => $config->{abismal}{target_dir} . "/result/",
      dnmtools_path      => $config->{DNMTools}{target_dir} . "/result/",
      MethylKitCorr_path => $config->{MethylKitCorr}{target_dir} . "/result/",
      MethylKitDiff_path => defined $config->{MethylKitDiff} ? $config->{MethylKitDiff}{target_dir} . "/result/" : undef,
    },
    parameterSampleFile3_ref => \@copy_files,
    parameterSampleFile4_ref => $webgestalt_task,
    parameterSampleFile5_ref => $abismal_task,
    parameterSampleFile6_ref => $dnmtools_task,
    parameterSampleFile7_ref => [ $abismal_summary_task, ".png\$", $methylkitcorr_task, ".png\$" ],
    parameterSampleFile8_ref => [ $methylkitdiff_task,   $methylkitdiffannovar_task ],
    parameterSampleFile9     => $version_files,
    sh_direct                => 1,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "40gb"
    },
  };
  push( @$tasks, "report" );

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => { tasks => $tasks, },
    sh_direct  => 0,
    pbs        => {
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return ($config);
} ## end sub getConfig


sub performWGBS {
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
} ## end sub performWGBS

1;

