#!/usr/bin/perl
package Pipeline::PeakPipelineUtils;

use strict;
use warnings;
use File::Basename;
use CQS::StringUtils;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Data::Dumper;
use Text::CSV;
use List::MoreUtils qw(first_index);
use Hash::Merge     qw( merge );

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
        add_homer_makeTagDirectory
        add_homer_mergePeaks
        add_homer_annotatePeaks
    )
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub init_treatments_design_table {
  my ($def) = @_;

  if ( getValue( $def, "treatments_auto" ) ) {
    my $files      = getValue( $def, "files" );
    my $treatments = {};
    for my $sample ( sort keys %$files ) {
      $treatments->{$sample} = [$sample];
    }
    $def->{treatments} = $treatments;
  } ## end if ( getValue( $def, "treatments_auto"...))

  my $treatments_pattern = $def->{"treatments_pattern"};
  if ($treatments_pattern) {
    my $files      = getValue( $def, "files" );
    my $treatments = {};
    for my $sample ( sort keys %$files ) {
      if ( $sample =~ /$treatments_pattern/ ) {
        $treatments->{$sample} = [$sample];
      }
    }
    $def->{treatments} = $treatments;
  } ## end if ($treatments_pattern)

  my $inputs_pattern = $def->{"inputs_pattern"};
  if ($inputs_pattern) {
    my $treatments_inputs_match_pattern = getValue( $def, "treatments_inputs_match_pattern" );
    my $files                           = getValue( $def, "files" );
    my $input_files                     = {};
    for my $sample ( sort keys %$files ) {
      if ( $sample =~ /$inputs_pattern/ ) {
        my $match_sample = capture_regex_groups( $sample, $treatments_inputs_match_pattern );
        $input_files->{$match_sample} = [$sample];
      }
    } ## end for my $sample ( sort keys...)

    my $inputs     = {};
    my $treatments = getValue( $def, "treatments" );
    for my $sample ( sort keys %$treatments ) {
      my $match_sample = capture_regex_groups( $sample, $treatments_inputs_match_pattern );
      my $match_input  = $input_files->{$match_sample};
      die "cannot find input file of $sample, check treatments_inputs_match_pattern" if ( !defined $match_input );
      $inputs->{$sample} = $match_input;
    } ## end for my $sample ( sort keys...)

    $def->{controls} = $inputs;
  } ## end if ($inputs_pattern)

  checkFileGroupPairNames( $def, [ "treatments", "controls" ] );

  my $design_table_condition_pattern = $def->{design_table_condition_pattern};
  if ( getValue( $def, "design_table_auto", 0 ) ) {
    my $task_name    = getValue( $def, "task_name" );
    my $treatments   = getValue( $def, "treatments" );
    my $design_table = {};

    my $rep_index = {};

    for my $sample ( sort keys %$treatments ) {
      my $condition = ( defined $design_table_condition_pattern ) ? capture_regex_groups( $sample, $design_table_condition_pattern ) : $sample;

      $design_table->{$sample} = {
        Condition => $condition,
        Replicate => get_next_index( $rep_index, $condition, 0 )
      };
    } ## end for my $sample ( sort keys...)

    $def->{design_table} = { $task_name => $design_table };
  } ## end if ( getValue( $def, "design_table_auto"...))
} ## end sub init_treatments_design_table


sub add_chipqc {
  my ( $config, $def, $tasks, $target_dir, $task_name, $bed_ref, $bam_ref ) = @_;

  if ( !defined $def->{"design_table"} ) {
    if ( defined $def->{design_table_condition_pattern} ) {
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
    is_paired_end  => getValue( $def, "is_paired_end" ),
    consensus      => getValue( $def, "chipqc_consensus",     1 ),
    pcaAttributes  => getValue( $def, "chipqc_pcaAttributes", "Tissue,Factor" ),
    pcaLabels      => getValue( $def, "chipqc_pcaLabels",     "Replicate" ),
    sh_direct      => 0,
    pbs            => {
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };
  push( @$tasks, $task_name );
} ## end sub add_chipqc


sub addPeakPipelineReport {
  my ( $config, $def, $tasks, $target_dir, $task_dic ) = @_;

  my $options = {
    "perform_cutadapt"         => [ getValue( $def, "perform_cutadapt" ) ],
    "cutadapt_option"          => [ getValue( $def, "cutadapt_option", "" ) ],
    "cutadapt_min_read_length" => [ getValue( $def, "min_read_length", 0 ) ],
    "is_paired_end"            => [ getValue( $def, "is_paired_end",   0 ) ? "TRUE" : "FALSE" ],
    "aligner"                  => [ getValue( $def, "aligner" ) ],
    "peak_caller"              => [ getValue( $def, "peak_caller" ) ]
  };

  my @report_files = ();
  my @report_names = ();
  my @copy_files   = ();

  my $version_files = get_version_files($config);

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

  if ( defined $config->{bowtie2_summary} ) {
    push( @report_files, "bowtie2_summary",       ".reads.png" );
    push( @report_files, "bowtie2_summary",       ".chromosome.png" );
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
    push( @report_names, "chipqc_data",    "chipqc_html" );
  } ## end if ( $def->{perform_chipqc...})

  if ( $def->{"perform_homer"} ) {
    my $homer_name = $task_dic->{homer};

    $options->{homer_result} = $config->{$homer_name}{target_dir};

    my $treatments = $def->{treatments};

    for my $key ( keys %$treatments ) {
      push( @report_files, $homer_name, "$key/knownResults.txt" );
      push( @report_names, "homer_" . $key );
    }
    push( @copy_files, $homer_name, "homerMotifs.all.motifs" );
  } ## end if ( $def->{"perform_homer"...})
  if ( $peakCallerTask =~ /macs2/ ) {
    $options->{macs2_result} = $config->{$peakCallerTask}{target_dir};
  }

  if ( $def->{perform_diffbind} ) {
    my $bindName = $task_dic->{diffbind};
    $options->{diffbind_result} = $config->{$bindName}{target_dir};
    if ( $def->{"perform_homer"} ) {
      my $bind_homer_name = $task_dic->{diffbind_homer};
      $options->{diffbind_homer_result} = $config->{$bind_homer_name}{target_dir};
    }
  } ## end if ( $def->{perform_diffbind...})

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
    sh_direct => 1,
    pbs       => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  };
  push( @$tasks, "report" );
} ## end sub addPeakPipelineReport


sub add_bamplot_by_gff {
  my ( $config, $def, $tasks, $target_dir, $task_name, $bam_ref ) = @_;

  my $plotgroups = $def->{plotgroups};
  if ( !defined $plotgroups ) {
    my $files         = getValue( $config, "files" );
    my @sortedSamples = sort keys %$files;
    $plotgroups = { getValue( $def, "task_name" ) => \@sortedSamples };
  }
  $config->{plotgroups} = $plotgroups;

  my $gffFile;
  if ( defined $def->{"bamplot_gff"} ) {
    $gffFile = $def->{"bamplot_gff"};
  }
  elsif ( defined $def->{annotation_locus} ) {
    $gffFile = $target_dir . "/annotation_locus.gff";
    open( my $fh, '>', $gffFile ) or die "Could not open file '$gffFile' $!";
    my $locusList = $def->{annotation_locus};
    my $count     = 0;
    for my $locus (@$locusList) {
      $count = $count + 1;
      $locus =~ s/,//g;

      my $locusName = $locus;
      $locusName =~ s/:/_/g;
      $locusName =~ s/-/_/g;

      my @parts     = split /:/, $locus;
      my $chr       = $parts[0];
      my $positions = $parts[1];
      my @pos       = split /-/, $positions;
      my $start     = $pos[0];
      my $end       = $pos[1];
      my $strand    = scalar(@parts) >= 3 ? $parts[2] : "+";
      print $fh $chr . "\t" . $locusName . "\tLOCUS\t" . $start . "\t" . $end . "\t.\t" . $strand . "\t.\n";
    } ## end for my $locus (@$locusList)
    close($fh);
  } ## end elsif ( defined $def->{annotation_locus...})
  else {
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
    draw_by_r          => getValue( $def, "bamplot_draw_by_r",        1 ),
    draw_by_r_width    => getValue( $def, "bamplot_draw_by_r_width",  10 ),
    draw_by_r_height   => getValue( $def, "bamplot_draw_by_r_height", 10 ),
    sh_direct          => 1,
    pbs                => {
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  };
  push @$tasks, ("bamplot");
} ## end sub add_bamplot_by_gff


sub add_nearest_gene {
  my ( $config, $def, $summary, $target_dir, $callName, $callFilePattern, $gene_bed ) = @_;
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
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  };
  push @$summary, ($geneName);
} ## end sub add_nearest_gene


sub add_homer_makeTagDirectory {
  my ( $config, $def, $tasks, $target_dir, $homer_makeTagDirectory_task, $bam_ref, $dependent_task ) = @_;

  my $homer_makeTagDirectory_option = getValue( $def, "homer_makeTagDirectory_option", "" );

  $config->{$homer_makeTagDirectory_task} = {
    class         => "CQS::ProgramWrapperOneToOne",
    perform       => 1,
    target_dir    => "${target_dir}/" . getNextFolderIndex($def) . "$homer_makeTagDirectory_task",
    interpretor   => "",
    program       => "",
    check_program => 0,
    option        => "
rm -rf __NAME__.makeTagDirectory.failed __NAME__.makeTagDirectory.succeed __NAME__.failed

# makeTagDirectory will create a file named without .bam under the folder of the original bam file, so we need to create a link with the same name to make homer work
ln -s __FILE__ __NAME__.tag.bam

makeTagDirectory __NAME__ $homer_makeTagDirectory_option __NAME__.tag.bam

status=\$?
if [[ \$status -eq 0 && -s __NAME__/tagAutocorrelation.txt && -s __NAME__/chrY.tags.tsv ]]; then
  touch __NAME__.makeTagDirectory.succeed
else
  echo \$status > __NAME__.makeTagDirectory.failed
  mv __NAME__ __NAME__.failed
fi

rm -f __NAME__.tag.bam

# __FILE2__

",
    source_ref               => $bam_ref,
    parameterSampleFile2_ref => $dependent_task,
    check_file_ext           => "__NAME__/chrY.tags.tsv",
    sh_direct                => 0,
    no_output                => 1,
    output_file_ext          => "__NAME__",
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "10gb"
    },
  };
  push @$tasks, ($homer_makeTagDirectory_task);
} ## end sub add_homer_makeTagDirectory


sub add_homer_mergePeaks {
  my ( $config, $def, $tasks, $target_dir, $homer_mergePeaks_task, $peak_ref, $dependent_task ) = @_;

  my $rename_1_py             = dirname(__FILE__) . "/../Homer/rename_find_peaks_path.py";
  my $rename_2_py             = dirname(__FILE__) . "/../Homer/rename_annotate_peaks_path.py";
  my $homer_genome            = getValue( $def, "homer_genome" );
  my $homer_mergePeaks_option = getValue( $def, "homer_mergePeaks_option", "-d given" );

  $config->{$homer_mergePeaks_task} = {
    class         => "CQS::ProgramWrapper",
    perform       => 1,
    target_dir    => $target_dir . "/" . $homer_mergePeaks_task,
    program       => "",
    check_program => 0,
    option        => "

names=(__SAMPLE_NAMES__)

files=(__FILE__)

for i in \${!files[\@]}; do
  zcat \${files[\$i]} > \${names[\$i]}.narrowPeaks
done

rm -f __NAME__.mergePeaks.failed __NAME__.mergePeaks.succeed

mergePeaks *.narrowPeaks \\
  $homer_mergePeaks_option -matrix __NAME__ \\
  -venn __NAME__.venn.txt > tmp.__NAME__.all_dGiven.peaks.txt

if [[ \$(wc -l <tmp.__NAME__.all_dGiven.peaks.txt) -ge 2 ]]; then
  export status=0
else
  export status=1
fi

if [[ \$status -ne 0 ]]; then
  mv tmp.__NAME__.all_dGiven.peaks.txt __NAME__.mergePeaks.failed
else
  touch __NAME__.mergePeaks.succeed

  mv tmp.__NAME__.all_dGiven.peaks.txt __NAME__.all_dGiven.peaks.txt

  echo rename_names=`date`
  python $rename_1_py \\
    __NAME__.all_dGiven.peaks.txt \\
    __NAME__.venn.txt \\
    __NAME__.count.matrix.txt \\
    __NAME__.logPvalue.matrix.txt \\
    __NAME__.logRatio.matrix.txt

  awk -v OFS='\\t' -F'\\t' '{ print \$2,\$3,\$4,\$1,\$6,\$5}' __NAME__.all_dGiven.peaks.txt > __NAME__.all_dGiven.peaks.bed
fi

for i in \${!files[\@]}; do
  rm -f \${names[\$i]}.narrowPeaks
done

# __parameterSampleFile2__

",
    source_ref                 => $peak_ref,
    source_join_delimiter      => " \\\n  ",
    source_name_join_delimiter => " \\\n  ",
    source_type                => "array",
    parameterSampleFile2_ref   => $dependent_task,
    output_file_ext            => ".all_dGiven.peaks.txt",
    output_to_result_dir       => 1,
    output_to_same_folder      => 1,
    no_output                  => 1,
    sh_direct                  => 1,
    pbs                        => {
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "40G"
    }
  };
  push( @$tasks, $homer_mergePeaks_task );

} ## end sub add_homer_mergePeaks


sub add_homer_annotatePeaks {
  my ( $config, $def, $tasks, $target_dir, $homer_annotatePeaks_task, $tagDirectory_task, $mergePeaks_ref ) = @_;

  my $rename_2_py  = dirname(__FILE__) . "/../Homer/rename_annotate_peaks_path.py";
  my $homer_genome = getValue( $def, "homer_genome" );

  $config->{$homer_annotatePeaks_task} = {
    class         => "CQS::ProgramWrapper",
    perform       => 1,
    target_dir    => $target_dir . "/" . $homer_annotatePeaks_task,
    program       => "",
    check_program => 0,
    option        => "

rm -f __NAME__.annotatePeaks.failed __NAME__.annotatePeaks.succeed

echo annotatePeaks_raw=`date`
annotatePeaks.pl __parameterFile1__ \\
  $homer_genome -raw -d \\
  __FILE__ > __NAME__.raw.count

status=\$?
if [[ \$status -ne 0 ]]; then
  mv __NAME__.raw.count __NAME__.annotatePeaks.failed
else
  touch __NAME__.annotatePeaks.succeed
  python $rename_2_py __NAME__.raw.count
fi

",
    source_ref                 => $tagDirectory_task,
    source_join_delimiter      => " \\\n  ",
    source_name_join_delimiter => " \\\n  ",
    source_type                => "array",
    parameterFile1_ref         => $mergePeaks_ref,
    output_file_ext            => ".raw.count",
    output_to_result_dir       => 1,
    output_to_same_folder      => 1,
    no_output                  => 1,
    sh_direct                  => 1,
    pbs                        => {
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "40G"
    }
  };
  push( @$tasks, $homer_annotatePeaks_task );

} ## end sub add_homer_annotatePeaks

1;
