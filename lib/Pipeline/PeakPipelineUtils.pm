#!/usr/bin/perl
package Pipeline::PeakPipelineUtils;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Data::Dumper;
use Text::CSV;
use List::MoreUtils qw(first_index);
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(
      add_chipqc
      addPeakPipelineReport
    )
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub add_chipqc {
  my ($config, $def, $tasks, $target_dir, $task_name, $bed_ref, $bam_ref) = @_;

  if(!defined $def->{"design_table"}){
    if(defined $def->{design_table_condition_pattern}){
      $def = init_design_table_by_pattern($def);
    }
  }

  $config->{$task_name} = {
    class          => "QC::ChipseqQC",
    perform        => 1,
    target_dir     => "${target_dir}/" . getNextFolderIndex($def) . $task_name,
    option         => "",
    source_ref     => $bam_ref,
    groups         => $def->{"treatments"},
    controls       => $def->{"controls"},
    qctable        => $def->{"design_table"},
    peaks_ref      => $bed_ref,
    peak_software  => "bed",
    genome         => getValue( $def, "chipqc_genome" ),
    combined       => getValue( $def, "chipqc_combined", 1 ),
    blacklist_file => $def->{"blacklist_file"},
    chromosomes    => $def->{"chipqc_chromosomes"},
    is_paired_end => getValue($def, "is_paired_end"),
    sh_direct      => 0,
    pbs            => {
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };
  push (@$tasks, $task_name);
}

sub addPeakPipelineReport {
  my ($config, $def, $tasks, $target_dir, $task_dic) = @_;

  my $options = {
    "perform_cutadapt"                     => [ getValue( $def, "perform_cutadapt") ],
    "cutadapt_option"                          => [ getValue( $def, "cutadapt_option",  "" ) ],
    "cutadapt_min_read_length"                  => [ getValue( $def, "min_read_length", 0 ) ],
    "is_paired_end" => [getValue( $def, "is_paired_end", 0 ) ? "TRUE" : "FALSE"],
    "aligner" => [ getValue( $def, "aligner") ],
    "peak_caller" => [ getValue( $def, "peak_caller") ]
  };

  my @report_files = ();
  my @report_names = ();
  my @copy_files   = ();

  my $version_files = get_version_files($config);

  if ( defined $config->{fastqc_raw_summary} ) {
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.baseQuality.tsv.png" );
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.sequenceGC.tsv.png" );
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.adapter.tsv.png" );
    push( @report_names, "fastqc_raw_per_base_sequence_quality", "fastqc_raw_per_sequence_gc_content", "fastqc_raw_adapter_content" );
  }

  if ( defined $config->{fastqc_post_trim_summary} ) {
    push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.baseQuality.tsv.png" );
    push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.sequenceGC.tsv.png" );
    push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.adapter.tsv.png" );
    push( @report_names, "fastqc_post_trim_per_base_sequence_quality", "fastqc_post_trim_per_sequence_gc_content", "fastqc_post_trim_adapter_content" );
  }

  if ( defined $config->{bowtie2_summary} ) {
    push( @report_files, "bowtie2_summary", ".reads.png" );
    push( @report_files, "bowtie2_summary", ".chromosome.png" );
    push( @report_names, "bowtie2_summary_reads", "bowtie2_summary_chromosome" );
  }

  if ( defined $config->{bowtie2_cleanbam_summary} ) {
    push( @report_files, "bowtie2_cleanbam_summary", ".chromosome.png" );
    push( @report_names, "bowtie2_cleanbam_summary_chromosome" );
  }

  my $peakCallerTask = $task_dic->{peak_caller};
  push( @copy_files, $peakCallerTask, ".bed\$" );

  push( @report_files, $task_dic->{peak_count}, ".txt" );
  push( @report_names, "peak_count" );

  if ( $def->{perform_chipqc} ) {
    my $chipqc_taskname = $task_dic->{chipqc};
    push( @report_files, $chipqc_taskname, ".rdata" );
    push( @report_files, $chipqc_taskname, ".html" );
    push( @report_names, "chipqc_data", "chipqc_html" );
  }

  if( $def->{"perform_homer"}){
    my $homer_name = $task_dic->{homer};

    $options->{homer_result} = $config->{$homer_name}{target_dir};

    my $treatments = $def->{treatments};

    for my $key ( keys %$treatments ) {
      push( @report_files, $homer_name, "$key/knownResults.txt" );
      push( @report_names, "homer_" . $key );
    }
    push( @copy_files, $homer_name, "homerMotifs.all.motifs" );
  }
  if ($peakCallerTask =~ /macs2/){
    $options->{macs2_result} = $config->{$peakCallerTask}{target_dir};
  }

  if($def->{perform_diffbind}) {
    my $bindName = $task_dic->{diffbind};
    $options->{diffbind_result} = $config->{$bindName}{target_dir};
    if($def->{"perform_homer"}) {
      my $bind_homer_name = $task_dic->{diffbind_homer};
      $options->{diffbind_homer_result} = $config->{$bind_homer_name}{target_dir};
    }
  }

  $config->{report} = {
    class                      => "CQS::BuildReport",
    perform                    => 1,
    target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . "report",
    report_rmd_file            => "../Pipeline/ChIPSeq.rmd",
    additional_rmd_files       => "Functions.R;../Pipeline/Pipeline.R",
    docker_prefix              => "report_",
    parameterSampleFile1_ref   => \@report_files,
    parameterSampleFile1_names => \@report_names,
    parameterSampleFile2       => $options,
    parameterSampleFile3_ref   => \@copy_files,
    parameterSampleFile4       => $version_files,
    # parameterSampleFile5       => $def->{software_version},
    # parameterSampleFile6       => $def->{groups},
    sh_direct                  => 1,
    pbs                        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$tasks, "report" );
}
