#!/usr/bin/perl
package Pipeline::PeakPipelineUtils;

use strict;
use warnings;
use CQS::StringUtils;
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
      init_treatments_design_table
      add_chipqc
      addPeakPipelineReport
      add_bamplot_by_gff
      add_nearest_gene
    )
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub init_treatments_design_table {
  my ($def) = @_;

  if (getValue($def, "treatments_auto")){
    my $files = getValue($def, "files");
    my $treatments = {};
    for my $sample (sort keys %$files){
      $treatments->{$sample} = [$sample];
    }
    $def->{treatments} = $treatments;
  }

  my $treatments_pattern = $def->{"treatments_pattern"};
  if ($treatments_pattern){
    my $files = getValue($def, "files");
    my $treatments = {};
    for my $sample (sort keys %$files){
      if ($sample =~ /$treatments_pattern/){
        $treatments->{$sample} = [$sample];
      }
    }
    $def->{treatments} = $treatments;
  }

  my $inputs_pattern = $def->{"inputs_pattern"};
  if ($inputs_pattern){
    my $treatments_inputs_match_pattern = getValue($def, "treatments_inputs_match_pattern");
    my $files = getValue($def, "files");
    my $input_files = {};
    for my $sample (sort keys %$files){
      if ($sample =~ /$inputs_pattern/){
        my $match_sample = capture_regex_groups($sample, $treatments_inputs_match_pattern);
        $input_files->{$match_sample} = [$sample];
      }
    }
  
    my $inputs = {};
    my $treatments = getValue($def, "treatments");
    for my $sample (sort keys %$treatments){
      my $match_sample = capture_regex_groups($sample, $treatments_inputs_match_pattern);
      my $match_input = $input_files->{$match_sample};
      die "cannot find input file of $sample, check treatments_inputs_match_pattern" if (!defined $match_input);
      $inputs->{$sample} = $match_input;
    }

    $def->{controls} = $inputs;
  }

  checkFileGroupPairNames($def, ["treatments", "controls"]);

  my $design_table_condition_pattern = $def->{design_table_condition_pattern};
  if(getValue($def, "design_table_auto", 0)){
    my $task_name = getValue($def, "task_name");
    my $treatments = getValue($def, "treatments");
    my $design_table = {};

    my $rep_index = {};

    for my $sample (sort keys %$treatments){
      my $condition = (defined $design_table_condition_pattern) ? capture_regex_groups($sample, $design_table_condition_pattern) : $sample;
      
      $design_table->{$sample} = {
        Condition => $condition,
        Replicate => get_next_index($rep_index, $condition, 0)
      };
    }

    $def->{design_table} = {
      $task_name => $design_table
    };
  }
}

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
    consensus => getValue($def, "chipqc_consensus", 1),
    pcaAttributes => getValue($def, "chipqc_pcaAttributes", "Tissue,Factor"),
    pcaLabels => getValue($def, "chipqc_pcaLabels", "Replicate"),
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
    additional_rmd_files       => "reportFunctions.R;../Pipeline/Pipeline.R",
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

sub add_bamplot_by_gff {
  my ($config, $def, $tasks, $target_dir, $task_name, $bam_ref) = @_;

  my $plotgroups = $def->{plotgroups};
  if ( !defined $plotgroups ) {
    my $files         = getValue($config, "files");
    my @sortedSamples = sort keys %$files;
    $plotgroups = { getValue( $def, "task_name" ) => \@sortedSamples };
  }
  $config->{plotgroups} = $plotgroups;

  my $gffFile;
  if (defined $def->{"bamplot_gff"}){
    $gffFile = $def->{"bamplot_gff"};
  }elsif (defined $def->{annotation_locus}){
    $gffFile = $target_dir . "/annotation_locus.gff";
    open(my $fh, '>', $gffFile) or die "Could not open file '$gffFile' $!";
    my $locusList = $def->{annotation_locus};
    my $count = 0;
    for my $locus (@$locusList){
      $count = $count + 1;
      $locus =~ s/,//g;

      my $locusName = $locus;
      $locusName =~ s/:/_/g; 
      $locusName =~ s/-/_/g; 

      my @parts = split /:/, $locus;
      my $chr = $parts[0];
      my $positions = $parts[1];
      my @pos = split /-/, $positions;
      my $start = $pos[0];
      my $end = $pos[1];
      my $strand = scalar(@parts) >= 3 ? $parts[2] : "+";
      print $fh $chr . "\t" . $locusName . "\tLOCUS\t" . $start . "\t" . $end . "\t.\t" . $strand . "\t.\n";
    }
    close($fh);
  }else{
    getValue( $def, "bamplot_gff" );
  }
  $config->{"bamplot"} = {
    class              => "Visualization::Bamplot",
    perform            => 1,
    target_dir         => "${target_dir}/" . getNextFolderIndex($def) . "bamplot",
    option             => getValue( $def, "bamplot_option" ),
    source_ref         => $bam_ref,
    groups_ref         => "plotgroups",
    gff_file           => $gffFile,
    is_rainbow_color   => 0,
    is_draw_individual => 0,
    is_single_pdf      => 1,
    draw_by_r => getValue($def, "bamplot_draw_by_r", 1),
    draw_by_r_width => getValue($def, "bamplot_draw_by_r_width", 10),
    draw_by_r_height => getValue($def, "bamplot_draw_by_r_height", 10),
    sh_direct          => 1,
    pbs                => {
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  };
  push @$tasks, ("bamplot");
}

sub add_nearest_gene {
  my ($config, $def, $summary, $target_dir, $callName, $callFilePattern, $gene_bed) = @_;
  my $geneName = $callName . "_gene";
  $config->{$geneName} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/${callName}",
    rtemplate                => "../Annotation/findNearestGene.r",
    output_file              => "",
    output_file_ext          => ".files.csv",
    parameterSampleFile1_ref => [ $callName, $callFilePattern ],
    parameterFile1           => $gene_bed,
    rCode                    => '',
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push @$summary, ($geneName);
}
