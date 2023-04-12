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

  $config->{trimgalore} = {
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

  $config->{abismal} = {
    class      => "Alignment::abismal",
    perform    => 1,
    option     => "",
    target_dir => "${targetDir}/abismal",
    chr_dir    => $chr_dir,
    picard     => $picard,
    interval_list => $interval_list,
    addqual_perlFile => $addqual_perlFile,
    abismal_index => $abismal_index,
    source_ref => "trimgalore",
    pbs        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "128",
      "mem"       => "80gb"
    },
  };
  push(@$tasks, "abismal");

  $config->{DNMTools} = {
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

  $config->{annovar} = {
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

  $config->{DNMToolsDiff} = {
    class      => "Methylation::DNMToolsDiff",
    perform    => 1,
    target_dir => "${targetDir}/DNMToolsDiff",
    option     => "",
    source_ref    => "pairs",
    methfile_ref  => [ "DNMTools", "^(?!.*?all).*\.meth\$" ],
    hmrfile_ref   => [ "DNMTools", ".hmr\$" ],
    minCpG        => 10,
    minSigCpG     => 5,
    perc_cut      => 0.25,
    FDR           => 0.5,
    mincov        => 2,
    chr_size_file => $chr_size_file,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => "60gb"
    },
  };
  push(@$tasks, "DNMToolsDiff");

  $config->{DNMToolsDiffAnnovar} = {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${targetDir}/DNMToolsDiffAnnovar",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => $annovar_buildver,
    remove_empty_source => 1,
    isBed      => 1,
    source_ref => [ "DNMToolsDiff", ".DMR.filtered\$|\.dmcpgs\$" ],
    pbs        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "40gb"
    },
  };
  push(@$tasks, "DNMToolsDiffAnnovar");

  $config->{DNMToolsDiffAnnovarGenes} = {
    class              => "CQS::ProgramWrapperOneToOne",
    perform            => 1,
    target_dir         => "$targetDir/DNMToolsDiffAnnovarGenes",
    interpretor => "perl",
    program => "../Methylation/get_gene_names.pl",
    #parameterSampleFile1_ref => [ "DNMToolsDiffAnnovar" ],
    source_ref => ["DNMToolsDiffAnnovar", ".dmcpgs.annovar.final.tsv\$" ],
    output_file_prefix => ".dmcpgs.annovar.final.tsv.genename.txt",
    output_ext => ".dmcpgs.annovar.final.tsv.genename.txt",
    output_by_file => 1,
    sh_direct          => 1,
    pbs                => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, "DNMToolsDiffAnnovarGenes");

  my $webgestalt_task = addWebgestalt($config, $def, $tasks, $targetDir, "DNMToolsDiffAnnovarGenes",  [ "DNMToolsDiffAnnovarGenes", ".genename.txt\$" ]);

  $config->{HOMER_DMR} = {
    class        => "Homer::FindMotifs",
    perform      => 1,
    target_dir   => "${targetDir}/HOMER_DMR",
    option       => "-nomotif",
    homer_genome => "hg19",
    source_ref   => [ "DNMToolsDiff", ".DMR.filtered\$" ],
    remove_empty_source => 1,
    sh_direct    => 0,
    pbs          => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "40gb"
    },
  };
  push(@$tasks, "HOMER_DMR");

  my @report_files = ();
  my @report_names = ();
  if ( defined $config->{fastqc_raw_summary} ) {
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.baseQuality.tsv.png" );
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.sequenceGC.tsv.png" );
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.adapter.tsv.png" );
    push( @report_names, "fastqc_raw_per_base_sequence_quality", "fastqc_raw_per_sequence_gc_content", "fastqc_raw_adapter_content" );
  }

  $config->{report} = {
    class                      => "CQS::BuildReport",
    perform                    => 1,
    target_dir                 => "$targetDir/report",
    report_rmd_file            => "../Methylation/dnmtools_report.Rmd",
    additional_rmd_files       => "../Pipeline/Pipeline.R;reportFunctions.R",
    parameterSampleFile1_ref   => \@report_files,
    parameterSampleFile1_names => \@report_names,
    parameterSampleFile2 => {
      task_name => $task_name,
      dnmtools_path => $config->{DNMTools}{target_dir} . "/result/",
    },
    parameterSampleFile3 => [],
    parameterSampleFile4_ref => [ $webgestalt_task, ".txt\$" ],
    sh_direct                  => 1,
    pbs                        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
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

