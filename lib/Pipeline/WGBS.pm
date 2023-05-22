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

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing",   1 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  $def->{perform_cutadapt} = 0;

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);
  my $tasks = [@$individual, @$summary];

  my $targetDir = $def->{target_dir};

  my $thread = getValue($def, "thread", 4);
  my $abismal_index = getValue($def, "abismal_index");
  my $chr_dir = getValue($def, "chr_dir");
  my $chr_size_file = getValue($def, "chr_size_file");
  my $annovar_buildver = getValue($def, "annovar_buildver");
  my $annovar_db = getValue($def, "annovar_db");
  my $annovar_param = getValue($def, "annovar_param");
  my $HOMER_perlFile = getValue($def, "HOMER_perlFile");
  my $addqual_perlFile = dirname(__FILE__) . "/../Methylation/add_qual.pl";
  my $picard = getValue($def, "picard");
  my $interval_list = getValue($def, "interval_list");

  my $trimgalore_task = "trimgalore";
  $config->{$trimgalore_task} = {
    class      => "Trimmer::TrimGalore",
    perform    => 1,
    target_dir => "${targetDir}/trimgalore",
    option     => "",
    source_ref => $untrimed_ref,
    extension  => "_val.fq.gz",
    pairend    => 1,
    init_command => "",
    pbs        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "48",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, "trimgalore");

  my $abismal_task = "abismal";
  $config->{$abismal_task} = {
    class      => "Alignment::abismal",
    perform    => 1,
    option     => "",
    thread     => $thread,
    target_dir => "${targetDir}/abismal",
    chr_dir    => $chr_dir,
    picard     => $picard,
    interval_list => $interval_list,
    addqual_perlFile => $addqual_perlFile,
    abismal_index => $abismal_index,
    source_ref => "trimgalore",
    pbs        => {
      "nodes"     => "1:ppn=4",
      "walltime"  => "128",
      "mem"       => "80gb"
    },
  };
  push(@$tasks, "abismal");

  my $dnmtools_task = "DNMTools";
  $config->{$dnmtools_task} = {
    class         => "Methylation::DNMTools",
    perform       => 1,
    target_dir    => "${targetDir}/dnmtools",
    option        => "",
    chr_dir       => $chr_dir,
    chr_size_file => $chr_size_file,
    source_ref    => "abismal",
    pbs           => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "72",
      "mem"       => "80gb"
    },
  };
  push(@$tasks, "DNMTools");

  my $dnmtoolsannovar_task = "annovar";
  $config->{$dnmtoolsannovar_task} = {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${targetDir}/dnmtoolsAnnovar",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => $annovar_buildver,
    isBed      => 1,
    source_ref => [ "DNMTools", ".hmr\$|.amr\$|.pmr\$|.pmd\$" ],
    pbs        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "40gb"
    },
  };
  push(@$tasks, "annovar");

  my $methylkitprep_task = "MethylKitPreparation";
  $config->{$methylkitprep_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${targetDir}/methylkitprep",
    option                   => "",
    docker_command           => "singularity exec -c -e -B /home,/gpfs23,/gpfs51,/gpfs52,/panfs,/data,/dors,/nobackup,/tmp -H `pwd` /data/cqs/softwares/singularity/wgbs_r.1.1.sif ",
    rtemplate                => "../Methylation/prepare_CpG_input.R",
    output_file_ext          => ".CpG.txt",
    source_ref => [ "DNMTools", ".meth\$" ],
    #parameterSampleFile2_ref => [ "featurecount", ".count.summary" ],
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, "MethylKitPreparation");

  my $methylkitcorr_task = "MethylKitCorr";
  $config->{$methylkitcorr_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${targetDir}/methylkitcorr",
    #option                   => " --args ${task_name} hg19 group 4 ",
    docker_command           => "singularity exec -c -e -B /home,/gpfs23,/gpfs51,/gpfs52,/panfs,/data,/dors,/nobackup,/tmp -H `pwd` /data/cqs/softwares/singularity/wgbs_r.1.1.sif ",
    rtemplate                => "../Methylation/methylkit_corr.R",
    output_file_ext          => "_methyl_CpG_bvalue_corr_MDS_plot.png;_methyl_CpG_bvalue_corr_MDS_plot.pdf;.filtered.cpg.meth.rds",
    parameterSampleFile2 => {
      task_name => $task_name,
      org       => "hg19",
      var       => "group",
      mincov     => 4
    },
    #parameterSampleFile2_ref => [ "featurecount", ".count.summary" ],
    source_ref => [ "MethylKitPreparation", ".CpG.txt\$" ],
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "80gb"
    },
  };
  push(@$tasks, "MethylKitCorr");

  my $methylkitdiff_task = "MethylKitDiff";
  $config->{$methylkitdiff_task} = {
    class                    => "Methylation::MethylKitDiff",
    target_dir               => "${targetDir}/methylkitdiff",
    docker_command           => "singularity exec -c -e -B /home,/gpfs23,/gpfs51,/gpfs52,/panfs,/data,/dors,/nobackup,/tmp -H `pwd` /data/cqs/softwares/singularity/wgbs_r.1.1.sif ",
    rtemplate                => "../Methylation/methylkit_diff.R",
    #option                   => " --args ${task_name} 25 0.01 16 ",
    source_ref => "pairs",
    #parameterSampleFile1_ref => "pairs",
    parameterSampleFile2_ref => "groups",
    parameterSampleFile3 => {
      task_name => $task_name,
      diff      => 25,
      qvalue    => 0.01,
      ncore     => 16
    },
    parameterSampleFile4_ref => [ "MethylKitCorr", ".filtered.cpg.meth.rds\$" ],
    output_file_ext          => ".dmcpgs",
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "80gb"
    },
  };
  push(@$tasks, "MethylKitDiff");


#  my $dnmtoolsdiff_task = "DNMToolsDiff";
#  $config->{$dnmtoolsdiff_task} = {
#    class      => "Methylation::DNMToolsDiff",
#    perform    => 1,
#    target_dir => "${targetDir}/DNMToolsDiff",
#    option     => "",
#    source_ref    => "pairs",
#    methfile_ref  => [ "DNMTools", "^(?!.*?all).*\.meth\$" ],
#    hmrfile_ref   => [ "DNMTools", ".hmr\$" ],
#    parameterSampleFile1_ref => "pairs",
#    parameterSampleFile2_ref => "groups",
#    minCpG        => 10,
#    minSigCpG     => 5,
#    perc_cut      => 0.25,
#    FDR           => 0.05,
#    mincov        => 4,
#    chr_size_file => $chr_size_file,
#    pbs => {
#      "nodes"     => "1:ppn=1",
#      "walltime"  => "12",
#      "mem"       => "60gb"
#    },
#  };
#  push(@$tasks, "DNMToolsDiff");

  my $methylkitdiffannovar_task = "MethylKitDiffAnnovar";
  $config->{$methylkitdiffannovar_task} = {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${targetDir}/MethylKitDiffAnnovar",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => $annovar_buildver,
    perform_splicing => 0,
    remove_empty_source => 1,
    isBed      => 1,
    source_ref => [ "MethylKitDiff", ".dmcpgs\$" ],
    pbs        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "40gb"
    },
  };
  push(@$tasks, "MethylKitDiffAnnovar");

  $config->{MethylKitDiffAnnovarGenes} = {
    class              => "CQS::ProgramWrapperOneToOne",
    perform            => 1,
    target_dir         => "$targetDir/MethylKitDiffAnnovarGenes",
    interpretor => "perl",
    program => "../Methylation/get_gene_names.pl",
    #parameterSampleFile1_ref => [ "DNMToolsDiffAnnovar" ],
    source_ref => ["MethylKitDiffAnnovar", ".dmcpgs.annovar.final.tsv\$" ],
    output_file_prefix => ".dmcpgs.annovar.final.tsv.genename.txt",
    output_ext => ".dmcpgs.annovar.final.tsv.genename.txt",
    output_by_file => 0,
    sh_direct          => 1,
    pbs                => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, "MethylKitDiffAnnovarGenes");

  my $webgestalt_task = addWebgestalt($config, $def, $tasks, $targetDir, "MethylKitDiffAnnovarGenes",  [ "MethylKitDiffAnnovarGenes", ".genename.txt\$" ]);

#  my $homer_task = "HOMER_DMR";
#  $config->{$homer_task} = {
#    class        => "Homer::FindMotifs",
#    perform      => 1,
#    target_dir   => "${targetDir}/HOMER_DMR",
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

  my @report_files = ();
  my @report_names = ();
  my @copy_files   = ();
  if ( defined $config->{fastqc_raw_summary} ) {
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.baseQuality.tsv.png" );
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.sequenceGC.tsv.png" );
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.adapter.tsv.png" );
    push( @report_names, "fastqc_raw_per_base_sequence_quality", "fastqc_raw_per_sequence_gc_content", "fastqc_raw_adapter_content" );
  }
  if(( defined $dnmtools_task ) && ( defined $config->{$dnmtools_task} )) {
    push( @copy_files, $dnmtools_task, ".amr\$" );
    push( @copy_files, $dnmtools_task, ".hmr\$" );
    push( @copy_files, $dnmtools_task, ".pmr\$" );
    push( @copy_files, $dnmtools_task, ".pmd\$" );
  }
  if(( defined $dnmtoolsannovar_task ) && ( defined $config->{$dnmtoolsannovar_task} )) {
    push( @copy_files, $dnmtoolsannovar_task, ".annovar.final.tsv\$" );
  }
  if (( defined $methylkitcorr_task ) && ( defined $config->{$methylkitcorr_task} )) {
    push( @copy_files, $methylkitcorr_task, ".pdf\$|.png\$|.rds\$" );
  }
  if (( defined $methylkitdiff_task ) && ( defined $config->{$methylkitdiff_task} )) {
    push( @copy_files, $methylkitdiff_task, ".dmcpgs\$" );
  }
  if (( defined $methylkitdiffannovar_task ) && ( defined $config->{$methylkitdiffannovar_task} )) {
    push( @copy_files, $methylkitdiffannovar_task, ".annovar.final.tsv\$" );
  }
  if ( defined $webgestalt_task ) {
    push( @copy_files, $webgestalt_task, "_geneontology_Biological_Process\$" );
    push( @copy_files, $webgestalt_task, "_geneontology_Cellular_Component\$" );
    push( @copy_files, $webgestalt_task, "_geneontology_Molecular_Function\$" );
    push( @copy_files, $webgestalt_task, "_pathway_KEGG\$" );
  }

  $config->{report} = {
    class                      => "CQS::BuildReport",
    perform                    => 1,
    target_dir                 => "$targetDir/report",
    docker_command             => "singularity exec -c -e -B /home,/gpfs23,/gpfs51,/gpfs52,/panfs,/data,/dors,/nobackup,/tmp -H `pwd` /data/cqs/softwares/singularity/wgbs_r.1.1.sif ",
    report_rmd_file            => "../Methylation/dnmtools_report.Rmd",
    additional_rmd_files       => "../Pipeline/Pipeline.R;reportFunctions.R",
    parameterSampleFile1_ref   => \@report_files,
    parameterSampleFile1_names => \@report_names,
    parameterSampleFile2 => {
      task_name => $task_name,
      meta_data => "../../" . $task_name . "_meta.tsv",
      abismal_path  => $config->{abismal}{target_dir} . "/result/",
      dnmtools_path => $config->{DNMTools}{target_dir} . "/result/",
      MethylKitCorr_path => $config->{MethylKitCorr}{target_dir} . "/result/",
      MethylKitDiff_path => $config->{MethylKitDiff}{target_dir} . "/result/",
    },
    parameterSampleFile3_ref   => \@copy_files,
    parameterSampleFile4_ref   => $webgestalt_task,
    parameterSampleFile5_ref   => $abismal_task ,
    parameterSampleFile6_ref   => $dnmtools_task ,
    parameterSampleFile7_ref   => [ $methylkitcorr_task, ".pdf\$|.png\$" ],
    parameterSampleFile8_ref   => $methylkitdiff_task,
    sh_direct                  => 1,
    pbs                        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "40gb"
    },
  };
  push( @$tasks, "report" );

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${targetDir}/sequencetask",
    option     => "",
    source     => {
      tasks => $tasks,
    },
    sh_direct => 0,
    pbs       => {
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return($config);
};


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
}

sub performWGBSTask {
  my ( $def, $task ) = @_;

  my $config = performScRNASeq($def, 0);

  performTask( $config, $task );

  return $config;
}


1;

