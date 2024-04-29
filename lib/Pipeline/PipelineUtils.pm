#!/usr/bin/perl
package Pipeline::PipelineUtils;

use strict;
use warnings;
use File::Basename;
use CQS::StringUtils;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Data::Dumper;
use Text::CSV;
#use Pipeline::WdlPipeline qw(addCollectAllelicCounts);
use List::MoreUtils qw(first_index);
use Hash::Merge qw( merge );
use POSIX qw/ceil/;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw( 
    do_add_gene_locus
    getIntermidiateDir
    getIndexName
    initPipelineOptions 
    readChromosomeFromDictFile
    addPreprocess 
    addPairendFastqValidation
    addFastQC 
    addBlastn 
    addBowtie 
    addPARalyzer
    addBowtie1PARalyzer 
    addBamStat 
    addOutputOption 
    getOutputFormat
    add_table_correlation
    getDEseq2TaskName 
    addDEseq2 
    addDeseq2Visualization 
    addDeseq2SignificantSequenceBlastn
    getBatchGroups 
    addHomerMotif 
    addDiffbind    
    addHomerAnnotation 
    addEnhancer 
    init_design_table_by_pattern    
    writeDesignTable 
    addMultiQC
    getNextIndex
    getNextFolderIndex 
    addCleanBAM 
    getReportDir 
    getSequenceTaskClassname
    addAnnovar 
    addAnnovarFilter 
    addAnnovarFilterGeneannotation
    addAnnovarMafReport
    addFilterMafAndReport
    addMafToIntervals
    addGATK4CNVGermlineCohortAnalysis
    addAllelicCountsForClonalAnalysis
    addSciCloneAndClonevol
    addPyCloneVIAndClonevol
    addApePhylogeneticTree
    AddMothurPipeline
    AddMothurPipelineVis
    addXHMM
    addGeneLocus
    annotateNearestGene
    checkFileGroupPairNames
    addStarFeaturecount
    add_salmon
    add_pairedend_fastq_gather
    add_pairedend_fastq_to_ubam
    add_BWA
    add_BWA_WGS
    add_BWA_summary
    add_BWA_and_summary
    get_expect_trunk
    add_BWA_and_summary_scatter
    addMarkduplicates
    addMarkduplicates_merge    
    addSequenceTask
    addFilesFromSraRunTable
    addWebgestalt
    addBamsnap
    addBamsnapLocus
    addPlotGene
    addSizeFactor
    addAnnotationLocus
    addAnnotationGenes
    add_peak_count
    add_alignment_summary
    add_bam_validation
    add_gsea
    add_unique_r
    add_maf_filter
    add_bowtie_index
    addLocusCoverage    
    addGeneCoverage
    get_next_index
    add_extract_bam_locus
    has_comparison
    add_featurecount
    add_md5
    add_bamplot
    add_fastq_screen
    writeAnnotationLocus_bed
    writeAnnotationLocus_gff
    add_split_fastq_dynamic
    )
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getIntermidiateDir {
  my ($defaultDir, $def) = @_;

  my $result = $defaultDir;
  if($def->{use_intermediate_dir} ){
    $result = $def->{target_dir} . "/intermediate_data";
    if (-e $result) {
      return($result);
    }
    if ($def->{use_intermediate_dir_aside}){
      $result = $def->{target_dir} . ".intermediate_data";
    }
    $result = create_directory_or_die( $result );
  }

  return $result;
}

sub readChromosomeFromDictFile {
  my ($dictFile, $primaryOnly) = shift;
  if(! defined $primaryOnly){
    $primaryOnly = 1;
  }

  my $result = [];
  open( my $fin, "<$dictFile" ) or die "Cannot open $dictFile";
  while ( my $line = (<$fin>) ) {
    chomp $line;
    my @parts = split( '\t', $line );
    if ($parts[0] eq "\@SQ"){
      my $chrName = substr($parts[1], 3);
      if ($primaryOnly){
        my $maxLength = rindex($chrName, "chr", 0) == 0 ? 5 : 2;
        if (length($chrName) <= $maxLength){
          push @$result, $chrName;
        }
      }else{
        push @$result, $chrName;
      }
    }
  }
  close($fin);

  return($result);
}

sub getNextIndex {
  my ($def, $key, $digital) = @_;

  if (! defined $digital){
    $digital = 2;
  }

  my $result = "";
  my $index = getValue( $def, $key, 1 );
  $result = sprintf( "%0" . $digital . "d", $index );
  $def->{$key} = $index + 1;

  return $result;
}

sub get_next_index {
  my ($def, $key, $format_index) = @_;
  if (not defined $def->{$key}){
    $def->{$key} = 1;
  }else{
    $def->{$key} = $def->{$key} + 1;
  }

  if(not defined $format_index){
    $format_index = 0;
  }

  #print("current_index = " . $def->{$key} . "\n");
  my $res = $format_index ? sprintf("_%02d", $def->{$key}) : $def->{$key};
  return($res);
}

sub getNextFolderIndex {
  my ($def) = @_;

  my $result = "";
  my $add_folder_index = getValue( $def, "add_folder_index", 0 );
  if ($add_folder_index) {
    my $folder_index = getValue( $def, "folder_index", 1 );
    $result = sprintf( "T%03d_", $folder_index );
    $def->{folder_index} = $folder_index + 1;
  }

  return $result;
}

sub addPairendFastqValidation {
  my ($config, $def, $individual, $parent_dir, $task_name, $source_ref) = @_;
  $config->{"$task_name"} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => $parent_dir . "/" . getNextFolderIndex($def) . "$task_name",
    option => "",
    use_tmp_folder => 1,
    suffix  => "_qc",
    interpretor => "python3",
    program => "../QC/validatePairendFastq.py",
    source_arg => "-i",
    source_ref => $source_ref,
    output_arg => "-o",
    output_file_prefix => ".txt",
    output_file_ext => ".txt",
    output_to_same_folder => 1,
    can_result_be_empty_file => 1,
    use_tmp_folder => getValue($def, "use_tmp_folder_paired_end_validation", 0),
    sh_direct   => 0,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "6",
      "mem"       => "10gb"
    }
  };
  push(@$individual, $task_name);
}

sub addFastQC {
  my ( $config, $def, $individual, $summary, $fastqcTask, $source_ref, $parentDir ) = @_;

  my $intermediateDir = getIntermidiateDir($parentDir, $def);

  my $pairend = is_paired_end( $def );
  my $curThread = $pairend? 2: 1;

  $config->{$fastqcTask} = {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => $intermediateDir . "/" . getNextFolderIndex($def) . $fastqcTask,
    option     => "",
    source_ref => $source_ref,
    cluster    => $def->{"cluster"},
    fastqc     => $def->{"fastqc"},
    use_tmp_folder => $def->{use_tmp_folder_fastqc},
    sh_direct  => 0,
    pbs        => {
      "nodes"    => "1:ppn=" . $curThread,
      "walltime" => "12",
      "mem"      => "10gb"
    },
  };
  push @$individual, $fastqcTask;

  my $summaryTask = $fastqcTask . "_summary";
  my $python_script = dirname(dirname(__FILE__)) . "/QC/fastQCSummary.py";
  my $task_name = getValue($def, "task_name");
  $config->{$summaryTask} = {
    class => "CQS::UniqueR",
    target_dir => $parentDir . "/" . getNextFolderIndex($def) . $summaryTask,
    init_command => "
python3 $python_script -i fileList1.txt -o ${task_name}.FastQC      
",
    rtemplate => "../QC/fastQCSummary.r",
    rReportTemplate => "../QC/fastQCSummary.Rmd;reportFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => ".FastQC.html",
    option => "",
    output_file_ext => ".FastQC.summary.tsv,.FastQC.summary.tsv.png,.FastQC.reads.tsv,.FastQC.reads.tsv.png,.FastQC.baseQuality.tsv,.FastQC.baseQuality.tsv.png,.FastQC.sequenceGC.tsv,.FastQC.sequenceGC.tsv.png,.FastQC.adapter.tsv,.FastQC.adapter.tsv.png,.FastQC.overrepresented.tsv",
    docker_prefix => "report_",
    can_result_be_empty_file => 0,
    parameterSampleFile1_ref => [$fastqcTask, "data.txt"],
    parameterSampleFile2 => {
      task_name => $task_name,
      email => getValue($def, "email"),
    },
    sh_direct  => 1,
    pbs        => {
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  };
  push @$summary, $summaryTask;

  return($summaryTask);
}

sub addBlastn {
  my ( $config, $def, $summary, $blastTask, $fastaTask, $filePattern, $parentDir ) = @_;

  $config->{$blastTask} = {
    class      => "Blast::Blastn",
    perform    => 1,
    target_dir => $parentDir . "/" . getNextFolderIndex($def) . "$blastTask",
    option     => "",
    source_ref => [ $fastaTask, $filePattern ],
    sh_direct  => 0,
    localdb    => $def->{blast_localdb},
    cluster    => $def->{cluster},
    pbs        => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => "10",
      "mem"       => "10gb"
    }
  };

  push @$summary, $blastTask;
}

sub addBowtie {
  my ( $config, $def, $individual, $taskName, $parentDir, $bowtieIndex, $sourceRef, $bowtieOption, $hours ) = @_;

  if ( !defined $hours ) {
    $hours = "23";
  }

  my $intermediateDir = getIntermidiateDir($parentDir, $def);

  $config->{$taskName} = {
    class                 => "Alignment::Bowtie1",
    perform               => 1,
    target_dir            => $intermediateDir . "/" . getNextFolderIndex($def) . $taskName,
    option                => $bowtieOption,
    source_ref            => $sourceRef,
    samonly               => 0,
    sh_direct             => getValue($def, "bowtie1_direct", 0),
    mappedonly            => 1,
    export_max_mapped     => $def->{export_max_mapped},
    cluster               => $def->{cluster},
    output_to_same_folder => $def->{bowtie1_output_to_same_folder},
    pbs                   => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $hours,
      "mem"       => "40gb"
    },
  };

  if(is_array($bowtieIndex) or defined $config->{$bowtieIndex}){
    $config->{$taskName}{bowtie1_index_ref} = $bowtieIndex;
  }else{
    $config->{$taskName}{bowtie1_index} = $bowtieIndex;
  }


  push @$individual, $taskName;

  return ($taskName);
}

sub addPARalyzer {
  my ( $config, $def, $individual, $taskName, $parentDir, $sourceRef ) = @_;

  $config->{$taskName} = {
    class          => "ParClip::PARalyzer",
    perform        => 1,
    target_dir     => $parentDir . "/" . getNextFolderIndex($def) . $taskName,
    option         => "",
    source_ref     => $sourceRef,
    genome2bit     => getValue( $def, "genome_2bit" ),
    mirna_db       => getValue( $def, "mirna_db" ),
    sorted_by_name => 0,
    sh_direct      => 0,
    pbs            => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "20gb"
    },
  };

  return ($taskName);
}

sub addBowtie1PARalyzer {
  my ( $config, $def, $individual, $taskName, $parentDir, $bowtieIndex, $sourceRef, $bowtieOption ) = @_;

  $config->{$taskName} = {
    class         => "ParClip::Bowtie1PARalyzer",
    perform       => 1,
    target_dir    => $parentDir . "/" . getNextFolderIndex($def) . $taskName,
    option        => $bowtieOption,
    source_ref    => $sourceRef,
    bowtie1_index => $bowtieIndex,
    genome2bit    => getValue( $def, "paralyzer_genome_2bit" ),
    mirna_db      => getValue( $def, "paralyzer_mirna_db" ),
    sh_direct     => 0,
    cluster       => $def->{cluster},
    pbs           => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "40gb"
    },
  };

  push @$individual, $taskName;

  return ($taskName);
}

sub addBamStat {
  my ( $config, $def, $summary, $taskName, $targetDir, $sourceRef ) = @_;

  $config->{$taskName} = {
    class                    => "CQS::UniqueR",
    target_dir               => $targetDir,
    perform                  => 1,
    suffix                   => "_bs",
    rtemplate                => "../Samtools/BamStat.r",
    output_file              => ".bamstat.csv",
    sh_direct                => 1,
    cluster                  => $def->{cluster},
    parameterSampleFile1_ref => $sourceRef,
    pbs                      => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$summary, $taskName;
}

sub getDEseq2TaskName {
  my ( $taskKey, $libraryKey, $def ) = @_;
  my $result = "deseq2_" . $taskKey;
  if ( defined $libraryKey && $libraryKey ne 'None' ) {
    $result = $result . "_" . $libraryKey;
    $result =~ s/\s/_/g;
  }
  if ( defined $def->{DE_task_suffix} ) {
    $result = $result . $def->{DE_task_suffix};
  }
  return $result;
}

sub getReportDir {
  my $def        = shift;
  my $report_dir = undef;
  if ( defined $def->{"output_to_report_dir"} && $def->{"output_to_report_dir"} ) {
    $report_dir = $def->{target_dir} . "/report";
  }
  return ($report_dir);
}

sub addOutputOption {
  my ( $def, $rcode, $key, $defaultValue, $alternativeKey ) = @_;
  my $result = $rcode;
  my $newkey = ( defined $alternativeKey ) ? $alternativeKey : $key;
  if ( $result !~ /$key/ ) {
    my $value = getValue( $def, $key, $defaultValue );
    if($value eq "FALSE"){
      $result = $result . "$newkey<-FALSE;";
    }elsif( $value ) {
      $result = $result . "$newkey<-TRUE;";
    }
    else {
      $result = $result . "$newkey<-FALSE;";
    }
  }
  return ($result);
}

sub getOutputFormat {
  my ( $def, $rcode, $is_DE ) = @_;
  my $result = $rcode;

  if(!defined $is_DE){
    $is_DE = 0;
  }

  $result = addOutputOption( $def, $result, "DE_outputPdf",         0,                          "outputPdf" );
  $result = addOutputOption( $def, $result, "DE_outputPng",         1,                          "outputPng" );
  $result = addOutputOption( $def, $result, "DE_outputTIFF",        0,                          "outputTIFF" );
  $result = addOutputOption( $def, $result, "DE_showVolcanoLegend", 0,                          "showVolcanoLegend" );

  return ($result);
}

sub get_output_hash {
  my ( $def, $result, $is_DE ) = @_;

  if(!defined $is_DE){
    $is_DE = 0;
  }

  $result->{outputPdf} = getValue($def, "DE_outputPdf", 0);
  $result->{outputPng} = getValue($def, "DE_outputPng", 1);
  $result->{outputTIFF} = getValue($def, "DE_outputTIFF", 0);
  $result->{showVolcanoLegend} = getValue($def, "DE_showVolcanoLegend", 0);
  $result->{usePearsonInHCA} = getValue($def, "use_pearson_in_hca", 0);

  if($is_DE){
    $result->{showLabelInPCA} = getValue($def, [ "DE_show_label_PCA", "showLabelInPCA", "show_label_PCA"], 1);
  }else{
    $result->{showLabelInPCA} = getValue($def, [ "showLabelInPCA", "show_label_PCA"], 1);
  }

  return ($result);
}

sub has_contrast {
  my $def = shift;
  my $pairs = $def->{pairs};
  for my $cname (keys %$pairs){
    my $cdef = $pairs->{$cname};
    if(ref($cdef) eq "HASH"){
      for my $cname (keys %$cdef){
        if($cname eq "contrast"){
          return(1);
        }
      }
    }
  }
  return(0);
}

sub addDEseq2 {
  my ( $config, $def, $summary, $taskKey, $countfileRef, $deseq2Dir, $DE_min_median_read, $libraryFile, $libraryKey, $feature_name_regex, $n_first ) = @_;

  my $taskName = getDEseq2TaskName( $taskKey, $libraryKey, $def );

  my $libraryFileKey = "library_file";
  if ( is_array($libraryFile) ) {
    $libraryFileKey = "library_file_ref";
  }

  my $rCode = getOutputFormat( $def, getValue($def, "DE_rCode", ""), 1 );
  $rCode = addOutputOption( $def, $rCode, "top25cv_in_hca", $def->{top25cv_in_hca}, "top25cvInHCA" );

  my $raw_rCode = $rCode;
  my $filterBaseMean = (defined $def->{filterBaseMean}) && $def->{filterBaseMean};
  if($filterBaseMean){
    #die("filterBaseMean");
    $rCode = $rCode . "filterBaseMean=1;filterBaseMeanValue=" . getValue($def, "filterBaseMeanValue", 30) . ";";
  }

  if(not defined $feature_name_regex){
    $feature_name_regex = getValue($def, "DE_feature_name_regex", "");
  }

  if(not defined $n_first){
    $n_first = getValue($def, "DE_n_first", -1);
  }

  $config->{$taskName} = {
    perform                      => 1,
    target_dir                   => $deseq2Dir . "/" . getNextFolderIndex($def) . "$taskName",
    output_to_dir                => getReportDir($def),
    option                       => "",
    sh_direct                    => 1,
    show_label_PCA               => getValue($def, ["DE_showLabelInPCA", "showLabelInPCA", "show_label_PCA"], 1),
    use_pearson_in_hca           => $def->{use_pearson_in_hca},
    show_DE_gene_cluster         => $def->{DE_show_gene_cluster},
    pvalue                       => $def->{DE_pvalue},
    fold_change                  => $def->{DE_fold_change},
    min_median_read              => $DE_min_median_read,
    add_count_one                => $def->{DE_add_count_one},
    top25only                    => $def->{DE_top25only},
    detected_in_both_group       => $def->{DE_detected_in_both_group},
    use_raw_p_value              => $def->{DE_use_raw_pvalue},
    text_size                    => $def->{DE_text_size},
    cluster                      => $def->{cluster},
    export_significant_gene_name => $def->{DE_export_significant_gene_name},
    cooksCutoff                  => $def->{DE_cooksCutoff},
    independentFiltering => $def->{DE_independentFiltering},
    covariance_name_index        => getValue($def, "covariance_name_index", 0),
    $libraryFileKey              => $libraryFile,
    library_key                  => $libraryKey,
    rCode                        => $rCode,
    docker_prefix => "deseq2_",
    parameterSampleFile1 => {
      feature_name_regex => $feature_name_regex,
      feature_name_filter => getValue($def, "DE_feature_name_filter", ""),
      enhanced_volcano_red_blue_only => getValue($def, "DE_enhanced_volcano_red_blue_only", 0),
      title_in_volcano => getValue($def, "DE_title_in_volcano", 1),
      caption_in_volcano => getValue($def, "DE_caption_in_volcano", 1),
      n_first => $n_first,
    },
    pbs                          => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => getValue($def, "DESeq2_walltime", "23"),
      "mem"       => getValue($def, "DESeq2_mem", "40gb"),
    },
  };

  if (defined $def->{pairs_config}) {
    $config->{$taskName}{class} = "Comparison::DESeq2config";
    $config->{$taskName}{source} = $def->{pairs_config};
  }else{
    my $groupNames = defined $def->{deseq2_groups} ? "deseq2_groups" : "groups";
    $config->{$taskName}{source_ref} = "pairs";
    $config->{$taskName}{groups_ref} = $groupNames;
    if(has_contrast($def)){
      $config->{$taskName}{class} = "Comparison::DESeq2contrast";
    }elsif(defined $def->{covariance_file}){
      $config->{$taskName}{class} = "Comparison::DESeq2covariance";
      $config->{$taskName}{covariance_file} = $def->{covariance_file};
    }else{
      $config->{$taskName}{class} = "Comparison::DESeq2";
    }
  }

  if ( is_array($countfileRef) ) {
    $config->{$taskName}{countfile_ref} = $countfileRef;
  }
  else {
    $config->{$taskName}{countfile} = $countfileRef;
  }
  push @$summary, $taskName;

  if($filterBaseMean){
    my $raw_section = $taskName . "_noFilterBaseMean";
    my $base_mean_task = {};
    my $filterSection = $config->{$taskName};
    for my $task_key (sort keys %$filterSection){
      $base_mean_task->{$task_key} = $filterSection->{$task_key};
    }
    $base_mean_task->{rCode} = $raw_rCode;
    $base_mean_task->{target_dir} = $base_mean_task->{target_dir} . "_noFilterBaseMean";
    $config->{$raw_section} = $base_mean_task;
    push @$summary, $raw_section;
  }

  return $taskName;
}

sub addDeseq2Visualization {
  my ( $config, $def, $summary, $taskKey, $deseq2Tasks, $dataVisualizationDir, $layoutName, $libraryKey ) = @_;

  my $taskName = getDEseq2TaskName( $taskKey, $libraryKey, $def ) . "_vis";

  my $deseq2FileRef = [];
  for my $deseq2Task (@$deseq2Tasks) {
    push @$deseq2FileRef, ( getDEseq2TaskName( $deseq2Task, $libraryKey, $def ), "_DESeq2.csv\$" );
  }

  $config->{$taskName} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $dataVisualizationDir . "/" . getNextFolderIndex($def) . "$taskName",
    output_to_dir            => getReportDir($def),
    rtemplate                => "DESeq2_all_vis.R",
    output_file              => ".${taskKey}.DESeq2.Matrix",
    output_file_ext          => ".png",
    remove_empty_parameter   => 1,
    parameterSampleFile1_ref => $deseq2FileRef,
    parameterSampleFile2     => $def->{$layoutName},
    rCode                    => 'useRawPvalue=' . $def->{DE_use_raw_pvalue} . ";",
    sh_direct                => 1,
    cluster                  => $def->{cluster},
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };
  push @$summary, $taskName;
  return $taskName;
}

sub addDeseq2SignificantSequenceBlastn {
  my ( $config, $def, $summary, $deseq2Task, $parentDir ) = @_;

  my $fastaTask = $deseq2Task . "_sequences";
  $config->{$fastaTask} = {
    class                  => "Blast::DESeq2SignificantReadToFasta",
    perform                => 1,
    target_dir             => $parentDir . "/" . getNextFolderIndex($def) . "$fastaTask",
    option                 => "",
    remove_empty_parameter => 1,
    source_ref             => [ $deseq2Task, "_DESeq2_sig.csv\$" ],
    sh_direct              => 1,
    cluster                => $def->{cluster},
    pbs                    => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    }
  };

  push @$summary, ($fastaTask);

  addBlastn( $config, $def, $summary, $fastaTask . "_blastn", $fastaTask, ".fasta\$", $parentDir );
}

sub getBatchGroups {
  my ($def) = @_;
  my $files = $def->{files};
  my $result;
  my $layout;

  if ( defined $def->{batch_groups_name_file_regex} ) {
    for my $regexName ( keys %{ $def->{batch_groups_name_file_regex} } ) {
      my $regexHash = $def->{batch_groups_name_file_regex}->{$regexName};
      my $fileRegex = $regexHash->{file};
      my $nameRegex = $regexHash->{name};

      my $lCols   = [];
      my $lRows   = [];
      my $lGroups = [];
      $layout->{$regexName} = {
        "Col_Group" => $lCols,
        "Row_Group" => $lRows,
        "Groups"    => $lGroups
      };
      my $pushed    = {};
      my $curGroups = {};
      $result->{$regexName} = $curGroups;

      for my $sample ( keys %$files ) {
        my $sampleFile = $files->{$sample}[0];
        my $group;
        if ( $sampleFile =~ /$fileRegex/igs ) {
          my $fileGroupName = $1;
          if ( $sample =~ /$nameRegex/igs ) {
            my $sampleGroupName = $1;
            my $groupName       = $sampleGroupName . "_" . $fileGroupName;
            if ( !exists( $pushed->{$groupName} ) ) {
              push( @$lCols,   $fileGroupName );
              push( @$lRows,   $sampleGroupName );
              push( @$lGroups, $groupName );
              $pushed->{$groupName} = "";
            }

            if ( !defined $curGroups->{$groupName} ) {
              $curGroups->{$groupName} = [];
            }
            my $groups = $curGroups->{$groupName};
            push @$groups, $sample;
          }
          else {
            die( $sample . " didn't match with name regex " . $nameRegex );
          }
        }
        else {
          die( $sample . " didn't match with file regex " . $fileRegex . " : " . $sampleFile );
        }
      }
    }
  }
  elsif ( defined $def->{batch_groups_file_regex} ) {
    for my $regexName ( keys %{ $def->{batch_groups_file_regex} } ) {
      my $regex     = $def->{batch_groups_file_regex}->{$regexName};
      my $curGroups = {};
      $result->{$regexName} = $curGroups;
      for my $sample ( keys %$files ) {
        my $sampleFile = $files->{$sample}[0];
        my $group;
        if ( $sampleFile =~ /$regex/igs ) {
          my $groupName = $1;
          if ( !defined $curGroups->{$groupName} ) {
            $curGroups->{$groupName} = [];
          }
          my $groups = $curGroups->{$groupName};
          push @$groups, $sample;
        }
        else {
          die( $sample . " didn't match with regex " . $regex . " : " . $sampleFile );
        }
      }
    }
  }
  elsif ( defined $def->{batch_groups_name_regex} ) {
    for my $regexName ( keys %{ $def->{batch_groups_file_regex} } ) {
      my $regex     = $def->{batch_groups_name_regex}->{$regexName};
      my $curGroups = {};
      $result->{$regexName} = $curGroups;
      for my $sample ( keys %$files ) {
        my $group;
        if ( $sample =~ /$regex/igs ) {
          my $groupName = $1;
          if ( !defined $curGroups->{$groupName} ) {
            $curGroups->{$groupName} = [];
          }
          my $groups = $curGroups->{$groupName};
          push @$groups, $sample;
        }
        else {
          die( $sample . " didn't match with regex " . $regex );
        }
      }
    }
  }
  elsif ( defined $def->{batch_groups} ) {
    $result = $def->{batch_groups};
  }

  return ( $result, $layout );
}

sub addHomerMotif {
  my ( $config, $def, $summary, $target_dir, $callName, $callFilePattern ) = @_;
  my $homerName = $callName . "_homer_motifs";
  $config->{$homerName} = {
    class        => "Homer::FindMotifs",
    option       => getValue( $def, "homer_option" ),
    perform      => 1,
    target_dir   => $target_dir . "/" . getNextFolderIndex($def) . $homerName,
    source_ref   => [ $callName, $callFilePattern ],
    homer_genome => getValue( $def, "homer_genome" ),
    sh_direct    => 1,
    pbs          => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };
  push @$summary, ($homerName);
  return $homerName;
}

sub addDiffbind {
  my ( $config, $def, $tasks, $target_dir, $task_name, $bam_ref, $peaks_ref ) = @_;

  $config->{$task_name} = {
    class                   => "Comparison::DiffBind",
    perform                 => 1,
    target_dir              => "${target_dir}/" . getNextFolderIndex($def) . "${task_name}",
    option                  => "",
    source_ref              => $bam_ref,
    groups                  => getValue($def, "treatments"),
    controls                => $def->{"controls"},
    design_table            => getValue($def, "design_table"),
    peaks_ref               => $peaks_ref,
    peak_software => getValue($def, "peak_software","bed"),
    use_version2 => getValue($def, "use_version2", 1),
    homer_annotation_genome => $def->{homer_annotation_genome},
    parameterSampleFile1 => {
      summits => getValue($def, "diffbind_summits", 200),
      consensus_minOverlap => getValue($def, "diffbind_consensus_minOverlap", 2),
    },
    can_result_be_empty_file => 1,
    sh_direct               => 0,
    pbs                     => {
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => getValue($def, "diffbind_mem", "80gb")
    },
  };
  push @$tasks, $task_name;
}

sub add_table_correlation {
  my ($config, $def, $tasks, $task_name, $task_dir, $count_file_ref, $parameterFile3_ref) = @_;

  my $rCode = getValue( $def, "correlation_rcode", "" );
  my $project_name = getValue($def, "task_name");

  my $options = {     
    "task_name" => $project_name,
    "email" => getValue($def, "email"),
    "affiliation" => getValue($def, "affiliation", ""),
    "draw_all_groups_in_HCA" => getValue($def, "draw_all_groups_in_HCA", 0),
    "draw_umap" => getValue($def, "draw_umap", 0),
    "use_green_red_color_in_hca" => getValue($def, "use_green_red_color_in_hca", 0),
    "top25cv_in_hca" => getValue($def, "top25cv_in_hca", 0),
    "showLabelInPCA" => getValue($def, ["showLabelInPCA", "show_label_PCA"], 0),
    "outputPdf" => getValue($def, "outputPdf", 0),
    "outputPng" => getValue($def, "outputPng", 1),
    "outputTIFF" => getValue($def, "outputTIFF", 0),
    "usePearsonInHCA" => getValue($def, "use_pearson_in_hca", 0),
    "n_first" => getValue($def, "correlation_n_first", -1),
    "onlySamplesInGroup" => getValue($def, "correlation_onlySamplesInGroup", 0),
    "useLeastGroups" => getValue($def, "useLeastGroups", 0),
    "minMedian" => getValue($def, "minMedian", 0),
    "minMedianInGroup" => getValue($def, "minMedianInGroup", 1),
    "totalCountKey" => getValue($def, "totalCountKey", "None"),
  };

  my $corr_output_file_ext = ".density.png;.density.individual.png.Correlation.png;.heatmap.png;.PCA.png;";
  my $corr_output_file_task_ext = "";
  if ( defined $def->{groups} ){
    my $groups = $def->{groups};
    if(scalar( keys %$groups ) >= 3 ) {
      $corr_output_file_task_ext = ".Group.heatmap.png;.Group.Correlation.png";
    } 
  }

  my $correlation_script="countTableGroupCorrelation.v2.R";

  $config->{$task_name} = {
    class           => "CQS::CountTableGroupCorrelation",
    perform         => 1,
    rCode           => $rCode,
    target_dir      => $task_dir,
    parameterSampleFile4 => $options,
    parameterFile3_ref        => $parameterFile3_ref,
    parameterFile4 => $def->{covariance_file},

    rtemplate       => "countTableVisFunctions.R,$correlation_script",
    rReportTemplate => "CountTableGroupCorrelation.Rmd;reportFunctions.R",
    run_rmd_independent => 1,

    output_file     => "parameterSampleFile1",
    output_file_ext           => $corr_output_file_ext,
    output_file_task_ext      => $corr_output_file_task_ext,
    output_to_result_dir      => getValue($def, "correlation_output_to_result_dir", 0),
    output_include_folder_name => getValue($def, "correlation_output_include_folder_name", 1),

    can_result_be_empty_file => 1,
    docker_prefix => "correlation_",
    sh_direct       => 1,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => getValue($def, "correlation_mem", "40gb")
    },
  };
  if ( is_array($count_file_ref) ) {
    $config->{$task_name}{parameterSampleFile1_ref} = $count_file_ref;
  }
  else {
    $config->{$task_name}{parameterSampleFile1} = { $project_name => [$count_file_ref] };
    $config->{$task_name}{output_to_result_dir} = 1;
  }

  if ( defined $def->{groups} ) {
    $config->{$task_name}{parameterSampleFile2} = { all => $def->{groups} };
  }

  if ( defined $def->{correlation_groups} ) {
    $config->{$task_name}{parameterSampleFile2} = $def->{correlation_groups};
  }

  if ( defined $def->{groups_colors} ) {
    $config->{$task_name}{parameterSampleFile3} = $def->{groups_colors};
  }

  if ( defined $def->{correlation_groups_colors} ) {
    $config->{$task_name}{parameterSampleFile3} = $def->{correlation_groups_colors};
  }

  if ( defined $config->{$task_name}{parameterSampleFile3} ) {
    my $colorGroups = $config->{$task_name}{parameterSampleFile3};
    my $corGroups   = $config->{$task_name}{parameterSampleFile2};
    for my $title ( keys %$corGroups ) {
      my $titleGroups = $corGroups->{$title};
      for my $subGroup ( keys %$titleGroups ) {
        if ( !defined $colorGroups->{$subGroup} ) {
          my %cgroups = %$colorGroups;
          die "Color of group '$subGroup' was not define in " . Dumper($colorGroups);
        }
      }
    }
  }
  push @$tasks, $task_name;
}

sub addHomerAnnotation {
  my ( $config, $def, $summary, $target_dir, $callName, $callFilePattern, $homerName ) = @_;
  if(!defined $homerName){
    $homerName = $callName . "_homer_annotation";
  }
  $config->{$homerName} = {
    class        => "Homer::Annotation",
    option       => getValue( $def, "homer_option" ),
    perform      => 1,
    target_dir   => $target_dir . "/" . getNextFolderIndex($def) . $homerName,
    source_ref   => [ $callName, $callFilePattern ],
    homer_genome => getValue( $def, "homer_genome" ),
    can_result_be_empty_file => 1,
    sh_direct    => 1,
    pbs          => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };
  push @$summary, ($homerName);
  return $homerName;
}

sub addEnhancer {
  my ( $config, $def, $individual, $summary, $target_dir, $enhancerName, $bam_ref, $peak_ref ) = @_;
  $config->{$enhancerName} = {
    class         => "Chipseq::Enhancer",
    perform       => 1,
    target_dir    => "${target_dir}/" . getNextFolderIndex($def) . $enhancerName,
    option        => "",
    source_ref    => $bam_ref,
    treatments    => $def->{treatments},
    peaks_ref     => $peak_ref,
    pipeline_dir  => getValue( $def, "enhancer_folder" ),
    genome        => getValue( $def, "enhancer_genome" ),
    genome_path   => getValue( $def, "enhancer_genome_path" ),
    gsea_path     => getValue( $def, "enhancer_gsea_path" ),
    gmx_path      => getValue( $def, "enhancer_gmx_path" ),
    cpg_path      => getValue( $def, "enhancer_cpg_path" ),
    activity_file => $def->{enhancer_activity_file},
    sh_direct     => 1,
    pbs           => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "40gb"
    },
  };
  push @$individual, ($enhancerName);

  my $enhancerVis = $enhancerName . "_vis";
  $config->{$enhancerVis} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/" . getNextFolderIndex($def) . $enhancerVis,
    option                   => "",
    rtemplate                => "../Chipseq/enhancerVis.R",
    output_file              => ".enhancer",
    output_file_ext          => ".log.tss.png",
    output_other_ext          => ".log.distal.png;.tss.tsv;.distal.tsv",
    sh_direct                => 1,
    parameterSampleFile1_ref => [ "$enhancerName", ".txt\$" ],
    sh_direct                => 1,
    pbs                      => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "40gb"
    },
  };
  push @$summary, $enhancerVis;

  my $enhancerVisCor = $enhancerVis . "_correlation";
  $config->{$enhancerVisCor} = {
    class                    => "CQS::CountTableGroupCorrelation",
    perform                  => 1,
    suffix                   => "_corr",
    target_dir               => $config->{$enhancerVis}->{target_dir},
    rtemplate                => "countTableVisFunctions.R,countTableGroupCorrelation.v2.R",
    output_file              => "parameterSampleFile1",
    output_file_ext          => ".Correlation.png",
    parameterSampleFile1_ref => [ $enhancerVis, ".tsv\$" ],
    sh_direct                => 1,
    pbs                      => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push @$summary, $enhancerVisCor;
}

sub addMultiQC {
  my ( $config, $def, $summary, $target_dir, $root_dir, $multiqc_depedents ) = @_;
  $config->{multiqc} = {
    class         => "QC::MultiQC",
    option        => getValue( $def, "multiqc_option", "" ),
    perform       => 1,
    target_dir    => $target_dir . "/" . getNextFolderIndex($def) . "multiqc",
    docker_prefix => "multiqc_",
    output_to_dir => getReportDir($def),
    source_ref    => $multiqc_depedents,
    root_dir      => $root_dir,
    sh_direct     => 1,
    pbs           => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };

  push @$summary, ("multiqc");
  return "multiqc";
}

sub addCleanBAM {
  my ( $config, $def, $individual, $task_name, $target_dir, $bam_ref ) = @_;

  my $pairend = is_paired_end( $def );

  my $cleanbam_option;
  my $minimum_insert_size;
  my $maximum_insert_size;
  if ($pairend) {
    $cleanbam_option     = "-f 3 -F 3852";
    $minimum_insert_size = getValue( $def, "minimum_insert_size" );
    $maximum_insert_size = getValue( $def, "maximum_insert_size" );
  }
  else {
    $cleanbam_option = "-F 3844";

  }

  $config->{$task_name} = {
    class                   => "ATACseq::CleanBam",
    perform                 => 1,
    target_dir              => $target_dir,
    option                  => $cleanbam_option,
    source_ref              => $bam_ref,
    picard_jar              => getValue( $def, "picard_jar" ),
    remove_chromosome       => $def->{remove_chromosome},
    keep_chromosome         => $def->{keep_chromosome},
    minimum_maq             => getValue( $def, "minimum_maq" ),
    minimum_insert_size     => $minimum_insert_size,
    maximum_insert_size     => $maximum_insert_size,
    blacklist_file          => $def->{blacklist_file},
    mark_duplicates         => getValue($def, "mark_duplicates", 1),
    is_paired_end           => $pairend,
    is_sorted_by_coordinate => 1,
    sh_direct               => 0,
    pbs                     => {
      "nodes"    => "1:ppn=1",
      "walltime" => "48",
      "mem"      => "40gb"
    },
  };
  push @$individual, $task_name;
}

sub writeDesignTable {
  my ( $target_dir, $section, $designtable, $bamfiles, $peaksfiles, $peakSoftware, $merged, $task_name, $treatments, $controls ) = @_;

  #print(Dumper($treatments));

  my $defaultTissue = getValue( $designtable, "Tissue", "" );
  my $defaultFactor = getValue( $designtable, "Factor", "" );

  my $result = {};
  my $condition_map = {};

  if ($merged) {
    my $conditions = {};
    my $mapFileName = "${task_name}.config.txt";
    my $mapfile     = $target_dir . "/" . $mapFileName;
    open( my $map, ">$mapfile" ) or die "Cannot create $mapfile";
    print $map "SampleID\tTissue\tFactor\tCondition\tReplicate\tbamReads\tControlID\tbamControl\tPeaks\tPeakCaller\n";
    for my $name ( sort keys %$designtable ) {
      if ( $name eq "Tissue" || $name eq "Factor" ) {
        next;
      }

      my $sampleList        = $designtable->{$name};
      my $defaultNameTissue = getValue( $sampleList, "Tissue", $defaultTissue );
      my $defaultNameFactor = getValue( $sampleList, "Factor", $defaultFactor );

      for my $sampleName ( sort keys %$sampleList ) {
        if ( $sampleName eq "Tissue" || $sampleName eq "Factor" || $sampleName eq "Comparison" || $sampleName eq "MinOverlap" ) {
          next;
        }

        my $entryMap = getValue( $sampleList, $sampleName );
        my $tissue   = getValue( $entryMap,   "Tissue", $defaultNameTissue );
        my $factor   = getValue( $entryMap,   "Factor", $defaultNameFactor );
        my $condition = $entryMap->{Condition} or die "Define Condition for $sampleName in designtable of section $section";
        my $replicate = $entryMap->{Replicate} or die "Define Replicate for $sampleName in designtable of section $section";
        my $peakFile  = $peaksfiles->{$sampleName};
        if(is_array($peakFile)){
          $peakFile = $peakFile->[0];
        }

        my $sampleId   = $treatments->{$sampleName};
        if(is_array($sampleId)){
          $sampleId = $sampleId->[0];
        }

        my $bamReads   = $bamfiles->{$sampleId}[0];
        my $controlId  = "";
        my $bamControl = "";

        if ( defined $controls && defined $controls->{$sampleName} ) {
          $controlId  = $controls->{$sampleName};
          if(is_array($controlId)){
            $controlId = $controlId->[0];
          }

          $bamControl = $bamfiles->{$controlId}[0];
        }

        print $map $sampleId . "\t"
          . $tissue . "\t"
          . $factor . "\t"
          . $condition . "\t"
          . $replicate . "\t"
          . $bamReads . "\t"
          . $controlId . "\t"
          . $bamControl . "\t"
          . $peakFile . "\t"
          . $peakSoftware . "\n";

        $conditions->{$condition} = 1;
      }
    }
    close($map);

    $result->{$task_name} = $mapfile;
    $condition_map->{$task_name} = $conditions;
  }
  else {
    for my $name ( sort keys %$designtable ) {
      if ( $name eq "Tissue" || $name eq "Factor" ) {
        next;
      }

      my $conditions = {};

      my $sampleList        = $designtable->{$name};
      my $defaultNameTissue = getValue( $sampleList, "Tissue", $defaultTissue );
      my $defaultNameFactor = getValue( $sampleList, "Factor", $defaultFactor );

      my $hasMinOverlap = 0;
      my $curdir      = create_directory_or_die( $target_dir . "/" . $name );
      my $mapFileName = "${name}.config.txt";
      my $mapfile     = $curdir . "/" . $mapFileName;
      open( my $map, ">$mapfile" ) or die "Cannot create $mapfile";
      print $map "SampleID\tTissue\tFactor\tCondition\tReplicate\tbamReads\tControlID\tbamControl\tPeaks\tPeakCaller\n";
      for my $sampleName ( sort keys %$sampleList ) {
        if ( $sampleName eq "MinOverlap" ) {
          $hasMinOverlap = 1;
          next;
        }

        if ( $sampleName eq "Tissue" || $sampleName eq "Factor" || $sampleName eq "Comparison" ) {
          next;
        }

        my $entryMap = getValue( $sampleList, $sampleName );
        my $tissue   = getValue( $entryMap,   "Tissue", $defaultNameTissue );
        my $factor   = getValue( $entryMap,   "Factor", $defaultNameTissue );
        my $condition = $entryMap->{Condition} or die "Define Condition for $sampleName in designtable of section $section";
        my $replicate = $entryMap->{Replicate} or die "Define Replicate for $sampleName in designtable of section $section";
        my $peakFile  = $peaksfiles->{$sampleName}->[0];

        my $sampleId   = $treatments->{$sampleName}->[0];
        my $bamReads   = $bamfiles->{$sampleId}[0];
        my $controlId  = "";
        my $bamControl = "";

        if ( defined $controls ) {
          $controlId  = $controls->{$sampleName}->[0];
          $bamControl = $bamfiles->{$controlId}[0];
        }

        print $map $sampleId . "\t"
          . $tissue . "\t"
          . $factor . "\t"
          . $condition . "\t"
          . $replicate . "\t"
          . $bamReads . "\t"
          . $controlId . "\t"
          . $bamControl . "\t"
          . $peakFile . "\t"
          . $peakSoftware . "\n";

        $conditions->{$condition} = 1;
      }
      close($map);

      if ($hasMinOverlap){
        my $overlapFileName = "${curdir}/${name}.minoverlap.txt";
        open( my $overlap, ">$overlapFileName" ) or die "Cannot create $overlapFileName";
        print $overlap "Condition\tminoverlap\n";
        my $overlapDef =  $sampleList->{"MinOverlap"};
        for my $oKey (sort keys %$overlapDef){
          print $overlap $oKey . "\t" . $overlapDef->{$oKey} . "\n";
        }
      }

      $result->{$name} = $mapfile;
      $condition_map->{$name} = $conditions;
    }
  }

  return $result, $condition_map;
}

sub init_design_table_by_pattern {
  my ( $def ) = @_;

  my $task_name = getValue($def, "task_name");
  my $treatments = getValue($def, "treatments");

  my $condition_pattern = getValue($def, "design_table_condition_pattern");
  my $factor_pattern = getValue($def, "design_table_factor_pattern", "");
  my $replicate_pattern = getValue($def, "design_table_replicate_pattern", "");

  my $samples = {};

  for my $sampleName ( sort keys %$treatments ) {
    my $condition = capture_regex_groups($sampleName, $condition_pattern);
    my $factor = ($factor_pattern eq "") ? undef : capture_regex_groups($sampleName, $factor_pattern);
    my $replicate = ($replicate_pattern eq "") ? "1" : capture_regex_groups($sampleName, $replicate_pattern);
    $samples->{$sampleName} = {
      'Condition' => $condition,
      'Factor' => $factor,
      'Replicate' => $replicate
    };
  }

  $def->{"design_table"} = {
    $task_name => $samples
  };

  return $def;
}

sub getSequenceTaskClassname {
  my $cluster = shift;
  my $result = $cluster eq "slurm" ? "CQS::SequenceTaskSlurmSlim" : "CQS::SequenceTask";
  return ($result);
}

sub getIndexName{
  my ($prefix, $suffix, $indexDic, $indexKey) = @_;
  my $index = defined $indexDic ? getNextIndex($indexDic, $indexKey) : "";
  return($prefix . $index . $suffix);
}

sub addAnnovar {
  my ( $config, $def, $summary, $target_dir, $source_name, $source_pattern, $prefix, $indexDic, $indexKey, $clean_folder, $perform_splicing, $output_to_same_folder ) = @_;
  if (not defined $prefix){
    $prefix = $source_name;
  }

  my $annovar_name;
  if(defined $indexDic){
    $annovar_name = getIndexName($prefix, "_annovar", $indexDic, $indexKey);
  }else{
    $annovar_name = $prefix . "_annovar";
  }

  my $isBed=0;
  if(defined $source_pattern){
    $isBed = $source_pattern =~ /.bed$/;
  } 

  my $source_ref = ( defined($source_pattern) and ( $source_pattern ne "" ) ) ? [ $source_name, $source_pattern ] : $source_name;
  $config->{$annovar_name} = {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/$annovar_name",
    source_ref => $source_ref,
    option     => getValue( $def, "annovar_param" ),
    annovar_db => getValue( $def, "annovar_db" ),
    buildver   => getValue( $def, "annovar_buildver" ),
    annovar_gzipped => getValue( $def, "annovar_gzipped", 0 ),
    docker_prefix => "annovar_",
    clean_folder => $clean_folder,
    perform_splicing => $perform_splicing,
    output_to_same_folder => $output_to_same_folder,
    sh_direct  => getValue($def, "annovar_sh_direct", 1),
    isBed => $isBed,
    isvcf      => $isBed?0:1,
    pbs        => {
      "nodes"    => "1:ppn=8",
      "walltime" => "10",
      "mem"      => "40gb"
    },
  };
  push @$summary, $annovar_name;

  return ($annovar_name);
}

sub addAnnovarFilter {
  my ( $config, $def, $summary, $target_dir, $annovar_name, $prefix, $indexDic, $indexKey ) = @_;
  if (not defined $prefix){
    $prefix = $annovar_name;
  }

  my $annovar_filter_name = getIndexName($prefix, "_filter", $indexDic, $indexKey);
  
  $config->{$annovar_filter_name} = {
    class               => "Annotation::FilterAnnovar",
    perform             => 1,
    target_dir          => "${target_dir}/$annovar_filter_name",
    source_ref          => $annovar_name,
    option              => $def->{annovar_filter},
    sh_direct           => 1,
    maximum_freq_values => getValue($def, "maximum_freq_values", "0.001,0.01,0.1,1.0"),
    filter_fq_equal_1   => $def->{filter_variants_fq_equal_1},
    pbs                 => {
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  };
  push @$summary, $annovar_filter_name;
  return ($annovar_filter_name);
}

sub addAnnovarFilterGeneannotation {
  my ( $config, $def, $summary, $target_dir, $annovar_filter_name ) = @_;
  my $annovar_filter_geneannotation_name = $annovar_filter_name . "_geneannotation";

  $config->{$annovar_filter_geneannotation_name} = {
    class                   => "Annotation::GenotypeAnnotation",
    perform                 => 1,
    target_dir              => $config->{$annovar_filter_name}{target_dir},
    source_ref              => [ $annovar_filter_name, ".snv.missense.tsv" ],
    detail_file_ref         => [ $annovar_filter_name, ".filtered.missense.tsv" ],
    detail_genes            => $def->{annotation_genes},
    option                  => "",
    sample_name_pattern     => ".",
    sample_name_suffix      => "",
    gene_names              => $def->{annotation_genes},
    draw_onco_print         => 1,
    onco_options            => $def->{onco_options},
    prepare_cbioportal_data => 1,
    sh_direct               => 1,
    pbs                     => {
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  };

  push @$summary, $annovar_filter_geneannotation_name;
  return ($annovar_filter_geneannotation_name);
}

sub addAnnovarMafReport {
  my ( $config, $def, $summary, $target_dir, $annovar_filter_name, $prefix, $indexDic, $indexKey ) = @_;

  my $annovar_to_maf = $prefix . getNextIndex($indexDic, $indexKey) . "_toMAF";
  $config->{$annovar_to_maf} = {
    class      => "Annotation::Annovar2Maf",
    perform    => 1,
    target_dir => $target_dir . "/" . $annovar_to_maf,
    source_ref => [ $annovar_filter_name, "\\.freq0\\..*.filtered.tsv" ],
    refBuild   => getValue( $def, "annovar_buildver" ),
    docker_prefix => "mafreport_",
    sh_direct  => 1,
    pbs        => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push @$summary, $annovar_to_maf;

  my $annovar_to_maf_report = $prefix . getNextIndex($indexDic, $indexKey) . "_report";

  my $re_rcode = ( defined $def->{family_info_file} ? "clinicalFeatures='" . $def->{family_info_feature} . "';" : "" );
  if(defined $def->{annovar_buildver}){
    $re_rcode = $re_rcode . "genome='" . $def->{annovar_buildver} . "'";
  } 

  $config->{$annovar_to_maf_report} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . $annovar_to_maf_report,
    rtemplate                => "../Annotation/mafReport.r",
    output_file              => "parameterSampleFile1",
    output_file_ext          => ".report.html",
    docker_prefix            => "mafreport_",
    parameterSampleFile1_ref => [ $annovar_to_maf, ".tsv.maf\$" ],
    parameterFile1           => $def->{family_info_file},
    sh_direct                => 1,
    rCode                    => $re_rcode,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "24",
      "mem"       => "10gb"
    },
  };
  push (@$summary, $annovar_to_maf_report);
  return($annovar_to_maf, $annovar_to_maf_report);
}

sub addFilterMafAndReport {
  my ( $config, $def, $summary, $target_dir, $mutect2call ) = @_;

#  my $mutect2_index_dic = {};
#  my $mutect2_index_key = "mafReport_Index";
#  my $taskName = $mutect2call . getNextIndex($mutect2_index_dic, $mutect2_index_key) . "_mergeAndMafreport";
  $mutect2call = (is_array($mutect2call) ? $mutect2call->[0] : $mutect2call);
  my $taskName = $mutect2call . "_mergeAndMafreport";

  my $rCode=( defined $def->{family_info_file} ? "clinicalFeatures=" . $def->{family_info_feature} . ";" : "" );
  $rCode=$rCode."genome=\"" . getValue($def, "annovar_buildver", "hg38") . "\";";

  $config->{$taskName}={
    class      => "CQS::UniqueR",
    perform    => 1,
    target_dir => "${target_dir}/${taskName}",
    rtemplate                  => "../CQS/muTect2MergeAndMafreport.R",
    parameterSampleFile1_ref=> [$mutect2call, ".maf"],
    parameterFile1           => $def->{family_info_file},
    rCode                    => $rCode,
    output_file_ext            => ".filter.allSamples.maf;.filter.allSamples.maf.report.html.RData",
    docker_prefix => "mafreport_",
    sh_direct  => 0,
    pbs        => {
      "nodes"    => "1:ppn=8",
      "walltime" => "4",
      "mem"      => "30gb"
    }
  };
  push @$summary, $taskName;
  return($taskName);
}

sub addMafToIntervals {
  my ( $config, $def, $target_dir,$summary, $prefix, $indexDic, $indexKey,$mafResults ) = @_;

  my $task = $prefix . getNextIndex($indexDic, $indexKey) . "_mafToIntervals";
  $config->{$task} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . '/' . $task,
      rtemplate                  => "../Variants/mafToIntervals.R",
      parameterSampleFile1_ref   => [ $mafResults, ".maf\$" ],
      output_to_result_directory => 1,
      output_file                => "",
      output_file_ext            => ".maf.intervals",
      sh_direct                  => 1,
      'pbs'                      => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '10'
      },
    };
  push @$summary, $task;
  return($task);
}


sub addAllelicCountsForClonalAnalysis {
  my ( $config, $def, $summary,$target_dir,$AllelicCountsFiles,$mafResults,$cnvFile ) = @_;

  if (!defined($def->{family_info_file}) | !defined($def->{patient_info_feature}) ) {
    die("Need define family_info_file and patient_info_feature in def so that samples from the same patients can be extracted")
  }
  my $rCode=( defined $def->{family_info_file} ? "patientFeature='" . $def->{patient_info_feature} . "';" : "" );

  my $task = "${AllelicCountsFiles}_PrepareClonalAnalysis";
  $config->{$task} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . '/' . $task,
      rtemplate                  => "../Variants/AllelicCountsForClonalAnalysis.R",
      parameterSampleFile1_ref   => [ $AllelicCountsFiles, ".allelicCounts.tsv\$" ],
      parameterFile1             => getValue( $def, "family_info_file" ),
      parameterFile2_ref         => [ $mafResults, ".maf\$" ],
      parameterFile3_ref         => [ $cnvFile, ".seg\$" ],
      output_to_result_directory => 1,
      output_file                => "",
      output_file_ext            => ".pyCloneVIInputList.txt;.sciCloneInputList.txt",
      rCode                      => $rCode,
      sh_direct                  => 1,
      'pbs'                      => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '10'
      },
  };
  push @$summary, $task;
  return($task);
}

sub addSciCloneAndClonevol {
  my ( $config, $def, $summary,$target_dir,$prepareClonalAnalysis ) = @_;

  my $task = "${prepareClonalAnalysis}_SciCloneAndClonevol";
  $config->{$task} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . '/' . $task,
      rtemplate                  => "../Variants/sciCloneAndClonEov.R",
      parameterFile1_ref         => [ $prepareClonalAnalysis, ".sciCloneInputList.txt\$" ],
      output_to_result_directory => 1,
      output_file                => "",
      output_file_ext            => "",
      sh_direct                  => 1,
      'pbs'                      => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '10'
      },
  };
  push @$summary, $task;
  return($task);
}

sub addPyCloneVIAndClonevol {
  my ( $config, $def, $summary,$target_dir,$prepareClonalAnalysis ) = @_;

  my $task = "${prepareClonalAnalysis}_PyCloneVIAndClonevol";
  $config->{$task} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . '/' . $task,
      rtemplate                  => "../Variants/pyCloneVIAndClonEov.R",
      parameterFile1_ref         => [ $prepareClonalAnalysis, ".pyCloneVIInputList.txt\$" ],
      output_to_result_directory => 1,
      output_file                => "",
      output_file_ext            => "",
      sh_direct                  => 1,
      'pbs'                      => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '10'
      },
  };
  push @$summary, $task;
  return($task);
}

sub addApePhylogeneticTree {
  my ( $config, $def, $summary,$target_dir,$mafResults ) = @_;

  if (!defined($def->{family_info_file}) | !defined($def->{patient_info_feature}) ) {
    die("Need define family_info_file and patient_info_feature in def so that samples from the same patients can be extracted")
  }
  my $rCode=( defined $def->{family_info_file} ? "patientFeature='" . $def->{patient_info_feature} . "';" : "" );

  my $task = "${mafResults}_ApePhylogeneticTree";
  $config->{$task} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . '/' . $task,
      rtemplate                  => "../Variants/ApePhylogeneticTree.R.R",
      parameterFile1             => getValue( $def, "family_info_file" ),
      parameterFile2_ref         => [ $mafResults, ".maf\$" ],
      output_to_result_directory => 1,
      output_file                => "",
      output_file_ext            => ".ApePhylogeneticTree.txt",
      rCode                      => $rCode,
      sh_direct                  => 1,
      'pbs'                      => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '10'
      },
  };
  push @$summary, $task;
  return($task);
}


sub AddMothurPipeline {
  my ( $config, $def, $summary,$target_dir, $source_def ) = @_;

my $taskName="mothur_pipeline";

my $projectName=getValue($def, "task_name");
my $mothurPipelineCodeFile = getValue($def, "mothurPipelineCodeFile");
#copy mothur code to current folder and add $batch file name to the beginning 
my $mothurPipelineCodeFileOut = "$target_dir/$taskName/result/$projectName.mothurPipeline.code";

#softlink of
#ln -s /scratch/cqs/zhaos/reference/mothur/silva/silva.seed_v132.pcr.align silva.v4.fasta
#ln -s /scratch/cqs/zhaos/reference/mothur/trainset16_022016.pds/* .
my $mothurSilvaFile = getValue($def, "mothurSilvaFile");
my $mothurTrainsetFastaFile = getValue($def, "mothurTrainsetFastaFile");
my $mothurTrainsetTaxFile = getValue($def, "mothurTrainsetTaxFile");

symlink ( $mothurSilvaFile, "$target_dir/$taskName/result/silva.v4.fasta" );
symlink ( $mothurTrainsetFastaFile, "$target_dir/$taskName/result/trainset16_022016.pds.fasta" );
symlink ( $mothurTrainsetTaxFile, "$target_dir/$taskName/result/trainset16_022016.pds.tax" );

#my $stability_batch = "$target_dir/$taskName/result/stability.files";
my $sampleToFiles = "$target_dir/$taskName/result/${projectName}__fileList1.list";

open CODESIN,"<$mothurPipelineCodeFile";
open CODESOUT,">$mothurPipelineCodeFileOut";
print CODESOUT "make.contigs(file='$sampleToFiles', processors=8) #make .trim.contigs";
while(<CODESIN>) {
  print CODESOUT $_;
}
close(CODESIN);
close(CODESOUT);

$config->{$taskName} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "$target_dir/$taskName",
    #init_command          => "",
    option                => "
mothur $mothurPipelineCodeFileOut

#__FILE__           
#__OUTPUT__
",
    interpretor           => "",
    check_program         => 0,
    program               => "mothur",
    source_ref            => $source_def,
    source_fileFirst => 0,
    source_join_delimiter => "\t",
    source_arg            => "",
    output_arg            => "",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => ".trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared",
    output_other_ext      => ".trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "40gb"
    },
  };
  return($taskName)
}

sub AddMothurPipelineVis {
  my ( $config, $def, $tasks, $target_dir, $mothurPipelineTask, $groups, $fastqQC_summary ) = @_;

  die "mothurPipelineTask is required." if not defined $mothurPipelineTask;
  die "fastqQC_summary is required." if not defined $fastqQC_summary;

  my $task = "${mothurPipelineTask}_vis";
  $config->{$task} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . '/' . $task,
      rtemplate                  => "../Microbiome/mothurPipeline.vis.R",
      parameterSampleFile1_ref   => $groups,
      parameterFile1_ref         => [ $mothurPipelineTask, ".shared\$" ],
      parameterFile2_ref         => [ $mothurPipelineTask, ".taxonomy\$" ],
      parameterFile3_ref         => [ $fastqQC_summary, ".reads.tsv\$" ],
      output_to_result_directory => 1,
      output_file                => "",
      output_file_ext            => "",
      sh_direct                  => 1,
      'pbs'                      => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '10'
      },
  };
  push @$tasks, $task;
  return($task);
}


sub addGATK4PreprocessIntervals {
  my ( $config, $def, $target_dir, $bam_ref, $prefix, $step1, $step2, $step3, $step4, $step5, $step6, $index ) = @_;
  if (!defined $index){
    $index = "";
  }

  my $result = $prefix . "gatk4_CNV_Germline${index}_PreprocessIntervals";
  if ( !defined $config->{$result} ) {
    my $interval_file;
    my $bin_option = "";
    if (defined $def->{is_wgs}) {
      if($def->{is_wgs}){
        $interval_file = getValue($def, "wgs_calling_regions_file");
        my $cnv_bin_length = getValue($def, "cnv-bin-length", 1000);
        $bin_option = "--bin-length $cnv_bin_length";
      }else{
        $interval_file = getValue($def, "covered_bed");
      }
    }else{
      $interval_file = getValue($def, "covered_bed");
    }

    #PreprocessIntervals at summary level
    $config->{$result} = {
      class             => "GATK4::PreprocessIntervals",
      option            => $bin_option,
      interval_file     => $interval_file,
      ref_fasta_dict    => getValue( $def, "ref_fasta_dict" ),
      ref_fasta         => getValue( $def, "ref_fasta" ),
      'sh_direct'       => 1,
      'perform'         => 1,
      'target_dir'      => $target_dir . "/" . $result,
      'pbs'             => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '2'
      },
    };
    push( @$step2, $result );
  }
  return ($result);
}

sub addGATK4CNVGermlineCohortAnalysis {
  my ( $config, $def, $target_dir, $bam_ref, $prefix, $step1, $step2, $step3, $step4, $step5, $step6 ) = @_;

  my $result = {};

  my $preprocessIntervalsTask = addGATK4PreprocessIntervals( $config, $def, $target_dir, $bam_ref, $prefix, $step1, $step2, $step3, $step4, $step5, $step6, "_01" );
  $result->{interval} = $preprocessIntervalsTask;

  my $chrCode = getValue($def, "has_chr_in_chromosome_name") ? ";addChr=1" : "";

  #CollectReadCounts at sample level
  my $CollectReadCounts = $prefix . "gatk4_CNV_Germline_02_CollectReadCounts";
  $result->{CollectReadCounts} = $CollectReadCounts;
  $config->{$CollectReadCounts} = {
    class                      => "GATK4::CollectReadCounts",
    source_ref                 => $bam_ref,
    option                     => getValue($def, "gatk4_CollectReadCounts_option", ""),
    preprocessed_intervals_ref => $preprocessIntervalsTask,
    ref_fasta_dict             => getValue( $def, "ref_fasta_dict" ),
    ref_fasta                  => getValue( $def, "ref_fasta" ),
    'sh_direct'                => 0,
    'perform'                  => 1,
    'target_dir'               => $target_dir . '/' . $CollectReadCounts,
    'pbs'                      => {
      'nodes'    => '1:ppn=1',
      'mem'      => '40gb',
      'walltime' => '10'
    },
  };
  push( @$step3, $CollectReadCounts );

  #FilterIntervals at summary level
  my $FilterIntervals = $prefix . "gatk4_CNV_Germline_03_FilterIntervals";
  $result->{FilterIntervals} = $FilterIntervals;
  $config->{$FilterIntervals} = {
    class                      => "GATK4::FilterIntervals",
    source_ref                 => $CollectReadCounts,
    option                     => "",
    preprocessed_intervals_ref => $preprocessIntervalsTask,
    ref_fasta_dict             => getValue( $def, "ref_fasta_dict" ),
    ref_fasta                  => getValue( $def, "ref_fasta" ),
    blacklist_file             => $def->{blacklist_file},
    contig_ploidy_priors       => getValue( $def, "contig_ploidy_priors_file" ),
    'sh_direct'                => 0,
    'perform'                  => 1,
    'target_dir'               => $target_dir . '/' . $FilterIntervals,
    'pbs'                      => {
      'nodes'    => '1:ppn=1',
      'mem'      => '40gb',
      'walltime' => '10'
    },
  };
  push( @$step4, $FilterIntervals );

  #DetermineGermlineContigPloidy at summary level
  my $DetermineGermlineContigPloidyCohortMode = $prefix . "gatk4_CNV_Germline_04_DetermineGermlineContigPloidyCohortMode";
  $result->{DetermineGermlineContigPloidyCohortMode} = $DetermineGermlineContigPloidyCohortMode;
  $config->{$DetermineGermlineContigPloidyCohortMode} = {
    class                  => "GATK4::DetermineGermlineContigPloidy",
    source_ref             => $CollectReadCounts,
    option                 => "",
    filtered_intervals_ref => $FilterIntervals,
    contig_ploidy_priors   => getValue( $def, "contig_ploidy_priors_file" ),
    'sh_direct'            => 0,
    'perform'              => 1,
    'target_dir'           => $target_dir . '/' . $DetermineGermlineContigPloidyCohortMode,
    'pbs'                  => {
      'nodes'    => '1:ppn=1',
      'mem'      => '40gb',
      'walltime' => '10'
    },
  };
  push( @$step4, $DetermineGermlineContigPloidyCohortMode );

  #GermlineCNVCaller at summary level
  my $GermlineCNVCaller = $prefix . "gatk4_CNV_Germline_05_GermlineCNVCaller";
  if ($def->{gatk4_cnv_by_scatter}){
    #scatter filter intervals
    my $ScatterIntervals = $GermlineCNVCaller . "_1_scatterIntervals";
    $config->{$ScatterIntervals} = {
      class       => "CQS::ProgramWrapperOneToMany",
      perform     => 1,
      target_dir  => "$target_dir/$ScatterIntervals",
      option      => "",
      interpretor => "python3",
      program     => "../GATK4/scatterIntervals.py",
      source_ref      => $FilterIntervals,
      source_arg            => "-i",
      output_to_same_folder => 1,
      output_arg            => "-o",
      output_file_prefix    => "",
      output_file_ext       => "._ITER_.interval_list",
      iteration_arg         => "-n",
      iteration             => $def->{gatk4_cnv_scatter_count},
      sh_direct             => 1,
      pbs                   => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "4gb"
      },
    };
    push( @$step4, $ScatterIntervals );

    my $GermlineCNVCaller_scatter = $GermlineCNVCaller . "_2_scatterCall";
    $config->{$GermlineCNVCaller_scatter} = {
      class                       => "GATK4::GermlineCNVCallerScatter",
      counts_ref                  => $CollectReadCounts,
      option                      => "",
      source_ref                  => $ScatterIntervals,
      contig_ploidy_calls_dir_ref => [ $DetermineGermlineContigPloidyCohortMode, "calls" ],
      'sh_direct'                 => 0,
      'perform'                   => 1,
      'target_dir'                => $target_dir . '/' . $GermlineCNVCaller_scatter,
      'pbs'                       => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '48'
      },
    };
    push( @$step4, $GermlineCNVCaller_scatter );
    $GermlineCNVCaller = $GermlineCNVCaller_scatter;
  }else{
    $config->{$GermlineCNVCaller} = {
      class                       => "GATK4::GermlineCNVCaller",
      source_ref                  => $CollectReadCounts,
      option                      => "",
      filtered_intervals_ref      => $FilterIntervals,
      contig_ploidy_calls_dir_ref => [ $DetermineGermlineContigPloidyCohortMode, "calls" ],
      'sh_direct'                 => 0,
      'perform'                   => 1,
      'target_dir'                => $target_dir . '/' . $GermlineCNVCaller,
      'pbs'                       => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '48'
      },
    };
    push( @$step4, $GermlineCNVCaller );
  }
  $result->{GermlineCNVCaller} = $GermlineCNVCaller;

  #PostprocessGermlineCNVCalls at sample level
  my $PostprocessGermlineCNVCalls = $prefix . "gatk4_CNV_Germline_06_PostprocessGermlineCNVCalls";
  $result->{PostprocessGermlineCNVCalls} = $PostprocessGermlineCNVCalls;
  $config->{$PostprocessGermlineCNVCalls} = {
    class                       => "GATK4::PostprocessGermlineCNVCalls",
    source_ref                  => $CollectReadCounts,
    calls_shard_path_ref        => [ $GermlineCNVCaller, "calls\$" ],
    model_shard_path_ref        => [ $GermlineCNVCaller, "model\$" ],
    option                      => "",
    contig_ploidy_calls_dir_ref => [ $DetermineGermlineContigPloidyCohortMode, "calls" ],
    has_chr_in_chromosome_name  => getValue($def, "has_chr_in_chromosome_name"),
    'sh_direct'                 => 0,
    'perform'                   => 1,
    'target_dir'                => $target_dir . '/' . $PostprocessGermlineCNVCalls,
    'pbs'                       => {
      'nodes'    => '1:ppn=1',
      'mem'      => '40gb',
      'walltime' => '10'
    },
  };
  push( @$step5, $PostprocessGermlineCNVCalls );

  if(not $def->{is_wgs}) {
    #CombineGCNV at summary level
    my $CombineGCNV = $prefix . "gatk4_CNV_Germline_07_CombineGCNV";
    $result->{CombineGCNV} = $CombineGCNV;
    $config->{$CombineGCNV} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               => $target_dir . '/' . $CombineGCNV,
      interpretor              => "python3",
      program                  => "../GATK4/combineGCNV.py",
      parameterSampleFile1_arg => "-i",
      parameterSampleFile1_ref => [ $PostprocessGermlineCNVCalls, ".genotyped_intervals.vcf.gz" ],
      parameterFile1_arg       => "-b",
      parameterFile1           => getValue( $def, "covered_bed" ),
      parameterFile2_arg       => "--annovar_db",
      parameterFile2           => $def->{perform_annovar} ? $def->{annovar_db} : undef,
      parameterFile3_arg       => "--annovar_buildver",
      parameterFile3           => $def->{perform_annovar} ? $def->{annovar_buildver} : undef,
      output_arg               => "-o",
      output_file_ext          => ".txt",
      sh_direct                => 1,
      'pbs'                    => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '10'
      },
    };
    push( @$step6, $CombineGCNV );
    
    my $sizeFactorTask = $prefix . "gatk4_CNV_Germline_08_SizeFactor";
    $result->{sizeFactor} = $sizeFactorTask;
    $config->{$sizeFactorTask} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               => $target_dir . '/' . $sizeFactorTask,
      interpretor              => "python3",
      program                  => "../GATK4/getBackgroundCount.py",
      parameterSampleFile1_arg => "-b",
      parameterSampleFile1_ref => $bam_ref,
      parameterFile1_arg       => "-i",
      parameterFile1           => getValue( $def, "covered_bed" ),
      parameterFile2_arg       => "-c",
      parameterFile2_ref       => [ $CombineGCNV ],
      output_arg               => "-o",
      output_file_ext          => ".txt",
      output_other_ext         => ".txt.sizefactor",
      sh_direct                => 1,
      'pbs'                    => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '10'
      },
    };
    push( @$step6, $sizeFactorTask );

    my $cnvIndex;
    if($def->{plotCNVGenes} && $def->{annotation_genes}){
      my $cnvGenes = $prefix . "gatk4_CNV_Germline_09_CNVGenesLocus";
      $result->{cnvGenes} = $cnvGenes;
      $config->{$cnvGenes} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $target_dir . '/' . $cnvGenes,
        rtemplate                => "../Annotation/getCNVGeneLocus.r",
        rCode                    => "host=\"" . getValue($def, "biomart_host") . "\";dataset=\"" . getValue($def, "biomart_dataset") . "\";symbolKey=\"" . getValue($def, "biomart_symbolKey") . "\"" . $chrCode,
        parameterFile1_ref       => [$CombineGCNV, ".txt"],
        output_file_ext          => ".bed",
        sh_direct                => 1,
        'pbs'                    => {
          'nodes'    => '1:ppn=1',
          'mem'      => '40gb',
          'walltime' => '10'
        },
      };
      push( @$step6, $cnvGenes );
      
      my $plotCNVgenes = $prefix . "gatk4_CNV_Germline_10_CNVGenesPlot";
      $result->{plotCNVgenes} = $plotCNVgenes;
      $config->{$plotCNVgenes} = {
        class                 => "CQS::ProgramWrapper",
        perform               => 1,
        target_dir            => $def->{target_dir} . "/$plotCNVgenes",
        option                => "",
        interpretor           => "python3",
        program               => "../Visualization/plotCNV.py",
        parameterFile1_arg => "-i",
        parameterFile1_ref => [ $cnvGenes, ".bed" ],
        parameterFile2_arg => "-c",
        parameterFile2_ref => [$CombineGCNV, ".txt"],
        parameterFile3_arg => "-s",
        parameterFile3_ref => [$sizeFactorTask, ".sizefactor"],
        parameterSampleFile1_arg => "-b",
        parameterSampleFile1_ref => $bam_ref,
        output_to_result_directory => 1,
        output_arg            => "-o",
        output_file_ext       => ".position.txt",
        sh_direct             => 1,
        pbs                   => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
        },
      };

      push( @$step6, $plotCNVgenes );
      $cnvIndex = "11"; 
    }else{
      $cnvIndex = "07"; 
    }

    my $annotationGenesPlot = undef;
    if(defined $def->{annotation_genes} && defined $config->{"annotation_genes_locus"}){
      $annotationGenesPlot = $prefix . "gatk4_CNV_Germline_" . $cnvIndex . "_AnnotationGenesPlot";
      $result->{annotationGenesPlot} = $annotationGenesPlot;
      $config->{$annotationGenesPlot} = {
        class                 => "CQS::ProgramWrapper",
        perform               => 1,
        target_dir            => $def->{target_dir} . "/$annotationGenesPlot",
        option                => "",
        interpretor           => "python3",
        program               => "../Visualization/plotCNV.py",
        parameterFile1_arg => "-i",
        parameterFile1_ref => [ "annotation_genes_locus", ".bed" ],
        parameterFile2_arg => "-c",
        parameterFile2_ref => [$CombineGCNV, ".txt"],
        parameterFile3_arg => "-s",
        parameterFile3_ref => [$sizeFactorTask, ".sizefactor"],
        parameterSampleFile1_arg => "-b",
        parameterSampleFile1_ref => $bam_ref,
        output_to_result_directory => 1,
        output_arg            => "-o",
        output_file_ext       => ".position.txt",
        output_other_ext      => ".position.txt.slim",
        sh_direct             => 1,
        pbs                   => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
        },
      };

      push( @$step6, $annotationGenesPlot );
    }
  }elsif(getValue($def, "perform_gistic2", 0)){
    my $Gistic2SegFile = $prefix . "gatk4_CNV_Germline_07_Gistic2SegFile";
    $result->{Gistic2SegFile} = $Gistic2SegFile;
    $config->{$Gistic2SegFile} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               => $target_dir . '/' . $Gistic2SegFile,
      interpretor              => "python3",
      program                  => "../Format/gatk2gistic_segment_file.py",
      option                   => getValue($def, "Gistic2SegFile_option", "--no_chr --no_y"),
      parameterSampleFile1_arg => "-i",
      parameterSampleFile1_ref => [ $PostprocessGermlineCNVCalls, ".genotyped_segments.vcf.gz" ],
      parameterFile1_arg       => "-b",
      output_arg               => "-o",
      output_file_ext          => ".segmentation.txt",
      sh_direct                => 1,
      'pbs'                    => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '10'
      },
    };
    push( @$step6, $Gistic2SegFile );

    my $gistic2_sif = getValue($def, "gistic2_docker");
    my $gistic2_refgene = getValue($def, "gistic2_refgene");
    my $Gistic2 = $prefix . "gatk4_CNV_Germline_08_Gistic2";
    $result->{Gistic2} = $Gistic2;
    $config->{$Gistic2} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               => $target_dir . '/' . $Gistic2,
      interpretor              => "",
      program                  => "",
      check_program            => 0,
      option                   => "
GISTIC_LOC=/opt/GISTIC
DOCKER_OUTDIR=\${GISTIC_LOC}/run_result

singularity run -e --bind `pwd`:\${DOCKER_OUTDIR}  $gistic2_sif \\
  -b \${DOCKER_OUTDIR} -seg __parameterFile1__ -refgene \${GISTIC_LOC}/refgenefiles/$gistic2_refgene \\
  -rx 0 -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -twosize 1 \\
  -armpeel 1 -savegene 1 -maxseg 10000 -conf 0.99

#__OUTPUT__
",
      parameterFile1_arg => "-seg",
      parameterFile1_ref => $Gistic2SegFile,
      output_arg               => "-o",
      output_file_ext          => ".segmentation.txt",
      sh_direct                => 1,
      no_docker                => 1,
      'pbs'                    => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '10'
      },
    };
    push( @$step6, $Gistic2 );
  }
  
  return($result);
}

sub addXHMM {
  my ( $config, $def, $target_dir, $bam_ref, $step1, $step2, $step3, $step4, $step5, $step6 ) = @_;

  my $gatk3_jar = getValue( $def, "gatk3_jar" );

  my $interval_key = "interval_file";
  my $interval_value = getValue( $def, "covered_bed" );

  my $individual = $step1;
  my $summary    = $step2;
  my $docTask    = "XHMM_GATK3_DepthOfCoverage";
  my $cnvTask    = "XHMM_CNV";

  if ( $def->{cnv_xhmm_preprocess_intervals} ) {
    $interval_key   = "interval_file_ref";
    $interval_value = addGATK4PreprocessIntervals( $config, $def, $target_dir, $bam_ref, $step1, $step2, $step3, $step4, $step5, $step6 ), $individual = $step3;
    $summary        = $step4;
    $docTask        = "XHMM_gatk4_Intervals_GATK3_DepthOfCoverage";
    $cnvTask        = "XHMM_gatk4_Intervals_CNV";
  }

  $config->{$docTask} = {
    class         => "GATK::DepthOfCoverage",
    gatk_jar      => $gatk3_jar,
    source_ref    => $bam_ref,
    option        => "",
    $interval_key => $interval_value,
    ref_fasta     => getValue( $def, "ref_fasta" ),
    perform       => 1,
    sh_direct     => 0,
    target_dir    => $target_dir . "/$docTask",
    pbs           => {
      nodes    => '1:ppn=1',
      mem      => '40gb',
      walltime => '10'
    },
  };
  push( @$individual, $docTask );

  $config->{$cnvTask} = {
    class      => "CNV::XHMM",
    source_ref => [ $docTask, ".sample_interval_summary" ],
    ref_fasta  => getValue( $def, "ref_fasta" ),
    option     => "",
    perform    => 1,
    sh_direct  => 1,
    target_dir => $target_dir . "/$cnvTask",
    pbs        => {
      nodes    => '1:ppn=1',
      mem      => '40gb',
      walltime => '10'
    },
  };
  push( @$summary, $cnvTask );

  return ($cnvTask);
}

sub do_add_gene_locus {
  my ($config, $def, $tasks, $target_dir, $task_name, $genes_str) = @_;
  $genes_str =~ s/\s+/,/g;
  $config->{$task_name} = {
    class      => "CQS::UniqueR",
    perform    => 1,
    target_dir => $target_dir . '/' . $task_name,
    rtemplate  => "../Annotation/getGeneLocus.r",
    parameterSampleFile1 => {
      host=>getValue( $def, "biomart_host" ),
      dataset=>getValue( $def, "biomart_dataset" ),
      symbolKey=>getValue( $def, "biomart_symbolKey" ),
      genesStr=> $genes_str,
      add_chr => getValue($def, "annotation_genes_add_chr", getValue($def, "has_chr_in_chromosome_name", 0))
    },
    docker_prefix => "biomart_",
    rCode => "",
    output_file => "",
    output_file_ext => ".bed",
    sh_direct       => 1,
    'pbs'           => {
      'nodes'    => '1:ppn=1',
      'mem'      => '40gb',
      'walltime' => '10'
    },
  };
  push( @$tasks, $task_name );

  return ($task_name);
}

sub addGeneLocus {
  my ($config, $def, $summary, $target_dir) = @_;
  my $geneLocus = undef;
  if ( defined $def->{annotation_genes} ) {
    $geneLocus = "annotation_genes_locus";
    my $genesStr = getValue( $def, "annotation_genes" );
    do_add_gene_locus($config, $def, $summary, $target_dir, $geneLocus, $genesStr);
  }
  return ($geneLocus);
}

sub annotateNearestGene {
  my ($config, $def, $summary, $target_dir, $source_file_ref) = @_;

  my $task_name = "nearest_gene";

  $config->{$task_name} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/${task_name}",
    rtemplate                => "../Annotation/findNearestGene.r",
    output_file              => "",
    output_file_ext          => ".Category.Table.csv",
    parameterSampleFile1_ref => $source_file_ref,
    parameterFile1           => getValue($def, "gene_bed"),
    rCode                    => '',
    sh_direct                => 1,
    pbs                      => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push @$summary, ($task_name);
}

sub checkFileGroupPairNames {
  my ($def, $groupKeys, $pairKeys, $fileKey, $remove_missing) = @_;

  if(!defined $remove_missing){
    $remove_missing = 0;
  }

  my $files = {};
  if(!defined $fileKey){
    $fileKey = "files";
  }
  if(is_array($fileKey)){
    for my $fname (@$fileKey){
      $files->{$fname} = 1;
    }
  }else{
    $files = getValue($def, $fileKey);
  }
  my $bFailed = 0;

  if (!defined $groupKeys){
    $groupKeys = ["groups"];
  }

  if (!defined $pairKeys){
    $pairKeys = ["pairs"];
  }

  my $allGroupNames = {};
  for my $groupKey (@$groupKeys){
    if(defined $def->{$groupKey}){
      my $groups = $def->{$groupKey};
      for my $groupName (keys %$groups){
        if ($groupName =~ /^[.]/) {
          next;
        }
        my $sampleNames = $groups->{$groupName};

        my $newSampleNames = [];
        for my $sampleName (@$sampleNames){
          if (!defined $files->{$sampleName}){
            if(!$remove_missing){
              print STDERR "Sample $sampleName in $groupKey $groupName is not defined in files.\n";
              $bFailed = 1;
            }
          }else{
            push(@$newSampleNames, $sampleName);
          }
        }

        $groups->{$groupName} = $newSampleNames;

        $allGroupNames->{$groupName} = 1;
      }
    }
  }

  if(scalar(keys %$allGroupNames) > 0){
    for my $pairKey (@$pairKeys){
      if (defined $def->{$pairKey}){
        my $pairs = $def->{$pairKey};
        for my $pairName (keys %$pairs){
          my $groupNames = $pairs->{$pairName};
          if (is_array($groupNames)) {
            for my $groupName (@$groupNames){
              if (!defined $allGroupNames->{$groupName}){
                print STDERR "Group $groupName in $pairKey $pairName is not defined in groups.\n";
                $bFailed = 1;
              }
            }
          }
        }
      }
    }
  }

  if ($bFailed){
    print("Wrong definition detected, please fix it and run again.\n");
    exit;
  }
}

sub addStarFeaturecount {
  my ($config, $def, $individual, $summary, $target_dir, $source_ref, $suffix) = @_;
  my $aligner_index              = $def->{star_index} or die "Define star_index at definition first";
  my $transcript_gtf             = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
  my $star_featurecount_walltime = getValue( $def, "star_featurecount_walltime", 23 );
  my $star_memory = getValue( $def, "star_memory", 40 );
  my $star_option     = $def->{star_option};

  if(not defined $suffix) {
    $suffix = "";
  }
  my $starFolder                 = $target_dir . "/" . getNextFolderIndex($def) . "star_featurecount" . $suffix;

  my $star_task = "star_featurecount" . $suffix;
  
  $config->{$star_task} = {
    class                     => "Alignment::STARFeatureCount",
    perform                   => 1,
    target_dir                => $starFolder,
    option                    => $star_option,
    source_ref                => $source_ref,
    genome_dir                => $aligner_index,
    output_sort_by_coordinate => 1,
    use_tmp_folder            => $def->{star_use_tmp_folder},
    output_to_same_folder     => $def->{output_bam_to_same_folder},
    featureCount_option       => getValue( $def, "featureCount_option" ),
    star_location             => $def->{star_location},
    gff_file                  => $transcript_gtf,
    is_paired_end             => is_paired_end($def),
    delete_star_featureCount_bam => $def->{delete_star_featureCount_bam},
    sh_direct                 => 0,
    pbs                       => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => "$star_featurecount_walltime",
      "mem"       => "${star_memory}gb"
    },
  };
  push @$individual, ($star_task);

  my $bam_validation_ref = undef;
  if(getValue($def, "perform_bam_validation", 0)){
    add_bam_validation($config, $def, $individual, $target_dir, $star_task . "_bam_validation", [ $star_task, '_Aligned.sortedByCoord.out.bam$' ] );
    $bam_validation_ref = $star_task . "_bam_validation";
  }
  
  add_alignment_summary($config, $def, $summary, $target_dir, "${star_task}_summary", "../Alignment/AlignmentUtils.r;../Alignment/STARFeatureCount.r", ".FeatureCountSummary.csv;.FeatureCountSummary.csv.png;.STARSummary.csv;.STARSummary.csv.png;.chromosome.csv;.chromosome.png", [ $star_task, "_Log.final.out" ], [ $star_task, ".count.summary" ], [$star_task, ".chromosome.count"], $bam_validation_ref, [$star_task, '^(?!.*\.chromosome\.count).*\.count$'] );

  return($star_task);
};

sub add_salmon {
  my ($config, $def, $individual, $summary, $target_dir, $source_ref, $suffix) = @_;
  my $aligner_index = $def->{salmon_index} or die "Define salmon_index at definition first";
  my $transcript_gtf = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
  my $salmon_walltime = getValue( $def, "salmon_walltime", 23 );
  my $salmon_memory = getValue( $def, "salmon_memory", 40 );
  my $salmon_option = getValue( $def, "salmon_option", "");
  my $thread = getValue( $def, "max_thread", 8);

  if(not defined $suffix) {
    $suffix = "";
  }
  my $salmon_folder = $target_dir . "/" . getNextFolderIndex($def) . "salmon" . $suffix;

  my $salmon_task = "salmon" . $suffix;
  
  $config->{$salmon_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => $salmon_folder,
    option                => "
salmon quant -i $aligner_index \\
  --libType A \\
  -1 __FILE__ \\
  -g $transcript_gtf \\
  -p $thread \\
  -z __NAME__.sam \\
  -o .

status=\$?
if [[ \$status -ne 0 ]]; then
  rm __NAME__.sam
  touch __NAME__.failed
else
  echo sort __NAME__.sam
  samtools sort -O BAM -o __NAME__.bam -@ $thread __NAME__.sam
  echo index __NAME__.bam
  samtools index __NAME__.bam
  rm __NAME__.sam
  touch __NAME__.succeed
fi

#__OUTPUT__
",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_arg            => "-1",
    source_ref            => $source_ref,
    source_join_delimiter => " -2 ",
    output_to_same_folder => 0,
    output_arg            => "",
    output_file_ext       => "quant.genes.sf",
    output_other_ext => "__NAME__.bam",
    samplename_in_result => 0,
    sh_direct => 0,
    pbs                       => {
      "nodes"     => "1:ppn=" . $thread,
      "walltime"  => "$salmon_walltime",
      "mem"       => "${salmon_memory}gb"
    },
  };
  push @$individual, ($salmon_task);

  return($salmon_task);
};

sub add_BWA {
  my ($config, $def, $tasks, $target_dir, $bwa_name, $source_ref, $rg_name_regex) = @_;

  $config->{ $bwa_name } = {
    class                 => "Alignment::BWA",
    perform               => 1,
    target_dir            => "${target_dir}/${bwa_name}",
    option                => getValue( $def, "bwa_option" ),
    bwa_index             => getValue( $def, "bwa_fasta" ),
    output_unmapped_fastq   => getValue( $def, "bwa_output_unmapped_fastq", 0),
    source_ref            => $source_ref,
    do_bam_stat => getValue($def, "bwa_do_bam_stat", 1),
    output_to_same_folder => 1,
    rg_name_regex         => $rg_name_regex,
    sh_direct             => 0,
    pbs                   => {
      "nodes"    => "1:ppn=" . getValue( $def, "max_thread"),
      "walltime" => "24",
      "mem"      => getValue($def, "bwa_memory", "40gb")
    },
  };

  push @$tasks, ( $bwa_name );
}

sub add_BWA_WGS {
  my ($config, $def, $tasks, $target_dir, $bwa_name, $source_ref, $rg_name_regex) = @_;

  $config->{ $bwa_name } = {
    class                 => "Alignment::BWA_WGS",
    perform               => 1,
    target_dir            => "${target_dir}/${bwa_name}",
    option                => getValue( $def, "bwa_option" ),
    bwa_index             => getValue( $def, "bwa_fasta" ),
    source_ref            => $source_ref,
    output_to_same_folder => 1,
    rg_name_regex         => $rg_name_regex,
    output_to_same_folder => 1,
    sh_direct             => 0,
    pbs                   => {
      "nodes"    => "1:ppn=" . getValue( $def, "max_thread"),
      "walltime" => getValue($def, "bwa_walltime", "24"),
      "mem"      => getValue($def, "bwa_memory", "40gb")
    },
  };

  push @$tasks, ( $bwa_name );
}

sub add_BWA_summary {
  my ($config, $def, $tasks, $target_dir, $bwa_summary, $bwa, $rg_name_regex) = @_;

  my $bwa_do_bam_stat = getValue($def, "bwa_do_bam_stat", 1);

  my $bamstat_ref = $bwa_do_bam_stat ? [ $bwa, ".bamstat" ] : undef;

  $config->{ $bwa_summary } = {
    class                 => "CQS::UniqueR",
    perform               => 1,
    target_dir            => "${target_dir}/${bwa_summary}",
    option                => "",
    rtemplate             => "../Alignment/AlignmentUtils.r;../Alignment/BWASummary.r",
    parameterSampleFile1_ref    => $bamstat_ref,
    parameterSampleFile2_ref    => [$bwa, ".chromosome.count"],
    output_file           => "",
    output_file_ext       => ".BWASummary.csv",
    output_other_ext      => ".reads.png;.chromosome.png",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };

  if (defined $rg_name_regex) {
    $config->{ $bwa_summary }{rCode} = "rg_name_regex='$rg_name_regex'";
  }

  push @$tasks, $bwa_summary;
}

sub add_BWA_and_summary {
  my ($config, $def, $tasks, $target_dir, $source_ref, $rg_name_regex) = @_;

  my $bwa = "bwa";
  add_BWA($config, $def, $tasks, $target_dir, $bwa, $source_ref, $rg_name_regex); 

  my $bwa_summary = "bwa_summary";
  add_BWA_summary($config, $def, $tasks, $target_dir, $bwa_summary, $bwa, $rg_name_regex); 

  return( [$bwa, ".bam\$"], $bwa);
}

sub add_merge_bam {
  my ($config, $def, $tasks, $target_dir, $merge_name, $source_ref) = @_;

  $config->{$merge_name} = {
    class                 => "CQS::ProgramWrapperManyToOne",
    perform               => 1,
    target_dir            => "$target_dir/$merge_name",
    option                => "
sambamba merge -t 8 __OUTPUT__ \\
  __INPUT__ 
",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $source_ref,
    source_arg            => "",
    source_join_delimiter => " \\\n  ",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => ".bam",
    output_file_ext       => ".bam,.bam.bai",
    use_tmp_folder => getValue($def, "merge_bam_use_tmp_folder", 0),
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=" . getValue($def, "merge_bam_threads", 8),
      "walltime"  => "48",
      "mem"       => "40gb"
    },
  };

  push @$tasks, (  $merge_name );
}

sub get_expect_trunk{
  my ($file_map, $max_file_size_gb, $trunk_file_size_gb) = @_;
  my $result = {};

  for my $sample (sort keys %$file_map){
    my $files = $file_map->{$sample};
    my $idx = 0;
    my $total = 0;
    while ($idx < scalar(@$files)){
      my $read1 = $files->[$idx];
      my $read1_filesize = -s $read1;
      $idx = $idx + 2;
      my $read1_gb = $read1_filesize / (1024 * 1024 * 1024);
      if ($read1_gb <= $max_file_size_gb){
        $total = $total + 1;
        next
      }

      my $trunk = ceil($read1_gb / $trunk_file_size_gb);
      $total = $total + $trunk;
    }
    $result->{$sample} = $total;
  }

  return($result);
}

sub add_split_fastq_dynamic {
  my ($config, $def, $tasks, $target_dir, $split_fastq, $source_ref) = @_;

  my $min_file_size_gb = getValue($def, "split_fastq_min_file_size_gb", 10);
  my $trunk_file_size_gb = getValue($def, "split_fastq_trunk_file_size_gb", 5);
  
  my $options = getValue($def, "call_fastqsplitter", 1) ? "--call_fastqsplitter ":"";
  my $no_docker = getValue($def, "fastqsplitter_by_docker", 1) ? 0 : 1;

  $config->{ $split_fastq } = {
    class       => "CQS::ProgramWrapperOneToMany",
    perform     => 1,
    target_dir  => "$target_dir/$split_fastq",
    interpretor => "python3",
    init_command => "set -o pipefail",
    program     => "../Format/splitFastqDynamic.py",
    option      => $options . " \\
  --min_file_size_gb $min_file_size_gb \\
  --trunk_file_size_gb $trunk_file_size_gb \\
  --fill_length 3 \\
  -i __FILE__ \\
  -o __NAME__

status=\$?
if [[ \$status -ne 0 ]]; then
  rm -f __NAME__.*
  touch __NAME__.failed
else
  touch __NAME__.succeed
  rm -f __NAME__.failed
fi
",
    source_ref      => $source_ref,
    source_arg            => "-i",
    source_join_delimiter => ",",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => "._ITER_.1.fastq.gz",
    output_other_ext      => "._ITER_.2.fastq.gz",
    zfill_iter_in_result  => 1,
    iteration_arg         => "",
    iteration_fill_length => 3,
    sh_direct             => 1,
    no_output => 1,
    no_docker => $no_docker,
    use_tmp_folder => 0, #for split fastq, we don't need tmp folder
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };

  my $file_map = get_raw_files($config, $split_fastq, "source");

  my $trunk_map = get_expect_trunk($file_map, $min_file_size_gb, $trunk_file_size_gb);

  $config->{$split_fastq}{iteration} = $trunk_map;

  push @$tasks, (  $split_fastq );
}

sub add_split_fastq {
  my ($config, $def, $tasks, $target_dir, $split_fastq, $source_ref) = @_;

  $config->{ $split_fastq } = {
    class       => "CQS::ProgramWrapperOneToMany",
    perform     => 1,
    target_dir  => "$target_dir/$split_fastq",
    option      => "",
    interpretor => "python3",
    program     => "../Format/splitFastq.py",
    source_ref      => $source_ref,
    source_arg            => "-i",
    source_join_delimiter => ",",
    parameterSampleFile2_arg => "--total_reads",
    parameterSampleFile2 => $def->{"total_reads"},
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => "._ITER_.1.fastq.gz",
    output_other_ext      => "._ITER_.2.fastq.gz",
    iteration_arg         => "--trunk",
    iteration             => getValue($def, "aligner_scatter_count"),
    iteration_fill_length => 3,
    zfill_iter_in_result  => 1,
    use_tmp_folder => getValue($def, "split_fastq_use_tmp_folder", 0),
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };

  push @$tasks, (  $split_fastq );
}

sub add_pairedend_fastq_gather {
  my ($config, $def, $tasks, $target_dir, $merge_fastq, $read1_ref, $read2_ref) = @_;
  $config->{$merge_fastq} = {
    class                 => "CQS::ProgramWrapperManyToOne",
    perform               => 1,
    target_dir            => "$target_dir/$merge_fastq",
    option                => "
rm -f __NAME__.failed __NAME__.succeed

cat __INPUT__ > __NAME__.unmapped.1.fq.gz

status=\$?
if [[ \$status -ne 0 ]]; then
  rm __NAME__.unmapped.1.fq.gz
  touch __NAME__.failed
else
  cat __INPUT2__ > __NAME__.unmapped.2.fq.gz

  status=\$?
  if [[ \$status -ne 0 ]]; then
    rm __NAME__.unmapped.1.fq.gz __NAME__.unmapped.2.fq.gz
    touch __NAME__.failed
  else
    touch __NAME__.succeed
  fi
fi

#__OUTPUT__
",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $read1_ref,
    source_arg            => "",
    source_join_delimiter => " ",
    parameterSampleFile2_ref            => $read2_ref,
    parameterSampleFile2_arg            => "",
    parameterSampleFile2_join_delimiter => " ",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => ".unmapped.1.fq.gz",
    output_file_ext       => ".unmapped.1.fq.gz,.unmapped.2.fq.gz",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };

  push(@$tasks, $merge_fastq);
}

sub add_pairedend_fastq_to_ubam {
  my ($config, $def, $tasks, $target_dir, $ubam_fastq, $fastq_ref) = @_;
  $config->{$ubam_fastq} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$ubam_fastq",
    docker_prefix         => "gatk4_",
    option                => "
rm -f __NAME__.failed __NAME__.succeed

gatk FastqToSam \\
  -FASTQ __FILE__ \\
  -OUTPUT __NAME__.tmp.bam \\
  -READ_GROUP_NAME __NAME__ \\
  -SAMPLE_NAME __NAME__ \\
  -LIBRARY_NAME __NAME__ 

gatk RevertSam \\
  -I __NAME__.tmp.bam \\
  -O __NAME__.unmapped.bam \\
  --SANITIZE true \\
  --MAX_DISCARD_FRACTION 0.005 \\
  --ATTRIBUTE_TO_CLEAR XT \\
  --ATTRIBUTE_TO_CLEAR XN \\
  --ATTRIBUTE_TO_CLEAR AS \\
  --ATTRIBUTE_TO_CLEAR OC \\
  --ATTRIBUTE_TO_CLEAR OP \\
  --SORT_ORDER queryname \\
  --RESTORE_ORIGINAL_QUALITIES true \\
  --REMOVE_DUPLICATE_INFORMATION true \\
  --REMOVE_ALIGNMENT_INFORMATION true

status=\$?
if [[ \$status -ne 0 ]]; then
  rm __NAME__.tmp.bam __NAME__.unmapped.bam
  touch __NAME__.failed
else
  rm __NAME__.tmp.bam
  touch __NAME__.succeed
fi

#__OUTPUT__
",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $fastq_ref,
    source_arg            => "",
    source_join_delimiter => " -FASTQ2 ",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => ".unmapped.bam",
    output_file_ext       => ".unmapped.bam",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };

  push(@$tasks, $ubam_fastq);
}

sub add_BWA_and_summary_scatter {
  my ($config, $def, $tasks, $target_dir, $source_ref) = @_;

  my $bwa_key = "bwa";

  my $pairend = is_paired_end( $def );
  my $step = $pairend ? 2 : 1;

  my $splitFastq = "bwa_00_splitFastq";

  my $perform_split_fastq = getValue($def, "perform_split_fastq", 0);
  if($perform_split_fastq eq "by_dynamic") {
    if (is_array($source_ref)){
      if ($source_ref->[0] ne "files"){
        die "split fastq files can only work when source_ref is files directly, now is [" . join(',', @$source_ref) . "]";
      }
    }else{
      if ($source_ref ne "files"){
        die "split fastq files can only work when source_ref is files directly, now is $source_ref";
      }
    }
    add_split_fastq_dynamic($config, $def, $tasks, $target_dir, $splitFastq, $source_ref);
  }elsif($perform_split_fastq eq "by_file"){
    $config->{ $splitFastq } = {
      class       => "CQS::FileScatterTask",
      source_ref  => $source_ref,
      step => $step,
    };
  }elsif($perform_split_fastq eq "by_scatter"){
    add_split_fastq($config, $def, $tasks, $target_dir, $splitFastq, $source_ref);
  }else{
    die 'wrong perform_split_fastq value, it should be one of [ 0, "by_dynamic", "by_file", "by_scatter"]';
  }

  if(getValue($def, "perform_paired_end_validation", 1)){
    my $fastq_validator = $splitFastq . "_validation";
    addPairendFastqValidation($config, $def, $tasks, $target_dir, $fastq_validator, $splitFastq);
  }

  my $rg_name_regex = "(.+)_ITER_";

  my $bwa = "bwa_01_alignment";
  if(0){
    add_BWA_WGS($config, $def, $tasks, $target_dir, $bwa, $splitFastq, $rg_name_regex);
  }else{
    add_BWA($config, $def, $tasks, $target_dir, $bwa, $splitFastq, $rg_name_regex); 
    my $bwa_summary = "bwa_03_summary";
    add_BWA_summary($config, $def, $tasks, $target_dir, $bwa_summary, $bwa, $rg_name_regex);
  }

  my $bam_section = undef;
  my $perform_mark_duplicates = getValue($def, "perform_mark_duplicates", 1);
  if ($perform_mark_duplicates) {
    $bam_section = "bwa_02_markduplicates";
    addMarkduplicates_merge($config, $def, $tasks, $target_dir, $bam_section, $bwa);
  }else{
    $bam_section = "bwa_02_merge";
    add_merge_bam($config, $def, $tasks, $target_dir, $bam_section, $bwa); 
  }

  if (getValue($def, "bwa_output_unmapped_fastq", 0)){
    my $merge_fastq = "bwa_04_unmapped_fastq";
    add_pairedend_fastq_gather($config, $def, $tasks, $target_dir, $merge_fastq, [$bwa, ".1.fq.gz"], [$bwa, ".2.fq.gz"]);

    if (getValue($def, "bwa_output_unmapped_bam", 0)){
      my $ubam_fastq = "bwa_05_unmapped_bam";
      add_pairedend_fastq_to_ubam($config, $def, $tasks, $target_dir, $ubam_fastq, [$merge_fastq, ".fq.gz"]);
    }
  }
  return( [ $bam_section, ".bam\$" ], $bam_section);
}


sub addMarkduplicates_merge {
  my ($config, $def, $tasks, $target_dir, $task_name, $source_ref) = @_;

  my $mem = getValue($def, "mark_duplicates_mem", "40gb");
  #replace gb with empty and minus 2
  my $java_memory_size = $mem;
  $java_memory_size =~ s/gb//g;
  $java_memory_size = $java_memory_size - 2;

  my $ref_fasta = getValue($def, "ref_fasta");

  #default is the location in gatk docker
  my $picard_jar = getValue($def, "picard_jar", "/opt/picard.jar");

  my $sort_order = getValue($def, "sort_order", "coordinate");

  $config->{$task_name} = {
    class                 => "CQS::ProgramWrapperManyToOne",
    perform               => 1,
    target_dir            => "$target_dir/$task_name",
    option                => "
gatk --java-options \"-Dsamjdk.compression_level=2 -Xms${java_memory_size}g\" \\
  MarkDuplicates \\
  --INPUT __INPUT__ \\
  --OUTPUT __OUTPUT__ \\
  --METRICS_FILE __NAME__.duplicates_metrics.txt \\
  --VALIDATION_STRINGENCY SILENT \\
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
  --ASSUME_SORT_ORDER $sort_order \\
  --CLEAR_DT false \\
  --ADD_PG_TAG_TO_READS false \\
  --CREATE_INDEX true

status=\$?
if [[ \$status -ne 0 ]]; then
  touch __OUTPUT__.failed
  rm -f __OUTPUT__.succeed __NAME__.duplicates_marked.bam  __NAME__.duplicates_marked.bai
else
  rm -f __OUTPUT__.failed
  touch __OUTPUT__.succeed
fi
",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $source_ref,
    source_arg            => "",
    source_join_delimiter => " \\\n  --INPUT ",
    output_to_same_folder => 1,
    output_arg            => "--OUTPUT",
    output_file_prefix    => ".duplicates_marked.bam",
    output_file_ext       => ".duplicates_marked.bam,.duplicates_marked.bai",
    use_tmp_folder => getValue($def, "mark_duplicates_use_tmp_folder", 0),
    sh_direct             => 0,
    docker_prefix        => "gatk4_",
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "mark_duplicates_walltime", "48"),
      "mem"       => getValue($def, "mark_duplicates_mem", "40gb")
    },
  };

  push(@$tasks, $task_name);
}

sub addMarkduplicates {
  my ($config, $def, $tasks, $target_dir, $task_name, $source_ref) = @_;

  $config->{$task_name} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$task_name",
    option                => "
sambamba markdup -p -t __THREAD__ --tmpdir tmp __FILE__ __NAME__.tmp.bam

status=\$?
if [[ \$status -ne 0 ]]; then
  rm -f __OUTPUT__.succeed
  touch __OUTPUT__.failed
  rm __NAME__.tmp.bam __NAME__.tmp.bam.bai
else
  rm -f __OUTPUT__.failed
  touch __OUTPUT__.succeed
  mv __NAME__.tmp.bam __OUTPUT__
  mv __NAME__.tmp.bam.bai __OUTPUT__.bai
fi
",
    interpretor           => "",
    program               => "",
    check_program         => 0,
    source_arg            => "",
    source_ref            => $source_ref,
    source_join_delimiter => " ",
    docker_prefix         => "cqs_",
    output_to_same_folder => 1,
    output_arg            => "",
    output_file_prefix    => ".duplicates_marked.bam",
    output_file_ext       => ".duplicates_marked.bam, .duplicates_marked.bam.bai",
    sh_direct             => 0,
    use_tmp_folder => getValue($def, "mark_duplicates_use_tmp_folder", 1),
    pbs                   => {
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => getValue($def, "mark_duplicates_memory", "40gb")
    },
  };

  push(@$tasks, $task_name);
}

sub addSequenceTask {
  my ($config, $def, $tasks, $target_dir, $summary_tasks) = @_;
  $config->{sequencetask} =  {
    class      => getSequenceTaskClassname(getValue($def, "cluster")),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => $tasks
    },
    sh_direct => 0,
    pbs       => {
      "nodes"     => "1:ppn=" . getValue($def, "max_thread"),
      "walltime"  => getValue($def, "sequencetask_run_time", 48),
      "mem"       => "40gb"
    },
  };

  if(defined $summary_tasks){
    $config->{sequencetask}{source}{step_2} = $summary_tasks;
  }
}

sub addFilesFromSraRunTable {
  my ($def, $filename) = @_;

  my $sample_name_column = getValue($def, "SraRunTable_sample_name_column", 'Sample Name');

  my $csv = Text::CSV->new ({
    binary    => 1,
    auto_diag => 1,
    sep_char  => ','    # not really needed as this is the default
  });

  my $files = {};
  open(my $fh, '<:encoding(utf8)', $filename) or die $!;
  my $headers = $csv->getline( $fh );
  my $sample_name_index = first_index { $_ eq $sample_name_column } @$headers;
  if ($sample_name_index < 0){
    die "cannot find '$sample_name_column' in '@$headers'"
  }

  #print("$sample_name_index = " . $headers->[$sample_name_index]);
  while (my $fields = $csv->getline( $fh )) {
    my $srr = $fields->[0];
    my $name = $fields->[$sample_name_index];
    $files->{$name} = [$srr];
  }
  close($fh);

  $def->{files} = $files;
}

sub addWebgestalt {
  my ($config, $def, $tasks, $target_dir, $prefix, $source_ref) = @_;

  my $webgestaltTaskName = $prefix . "_WebGestalt";
  $config->{$webgestaltTaskName} = {
    class            => "Annotation::WebGestaltR",
    perform          => 1,
    target_dir       => $target_dir . "/" . getNextFolderIndex($def) . $webgestaltTaskName,
    option           => "",
    source_ref       => $source_ref,
    organism         => getValue( $def, "webgestalt_organism" ),
    interestGeneType => $def->{interestGeneType},
    referenceSet     => $def->{referenceSet},
    sh_direct        => 1,
    pbs              => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "10gb"
    },
  };
  push @$tasks, "$webgestaltTaskName";

  return($webgestaltTaskName);
}

sub addBamsnap {
  my ($config, $def, $tasks, $target_dir, $task_name, $bed_ref, $bam_ref, $params) = @_;

  my $parameterFile1_key = "parameterFile1_ref";
  if( -e $bed_ref) {
    $parameterFile1_key = "parameterFile1";
  }

  if(not defined $params){
    $params = $def;
  }

  my $gene_track = ($def->{bamsnap_option} =~ /no_gene_track/) ? "" : "gene";
  my $width = getValue($def, "bamsnap_width", 2000);
  my $height = getValue($def, "bamsnap_height", 3000);

  $config->{$task_name} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "$target_dir/$task_name",
    docker_prefix         => "bamsnap_",
    option                => getValue($params, "bamsnap_option", ""),
    interpretor           => "python3",
    check_program         => 1,
    program               => "../Visualization/bamsnap.py",
    parameterSampleFile1_arg => "-b",
    parameterSampleFile1_ref => $bam_ref,
    parameterSampleFile2_arg => "-c",
    parameterSampleFile2  => getValue($params, "bamsnap_raw_option", {
      "-width" => $width,
      "-height" => $height
    }),
    parameterFile1_arg    => "-i",
    $parameterFile1_key   => $bed_ref,
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => ".txt",
    output_other_ext      => "",
    sh_direct             => 1,
    can_result_be_empty_file => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };
  push( @$tasks, $task_name );
}

sub addGeneCoverage {
  my ($config, $def, $tasks, $target_dir, $task_name, $gene_ref, $bam_ref, $locus_ref) = @_;

  my $python_script = dirname(__FILE__) . "/../Visualization/gene_coverage.py";
  my $width = getValue($def, "coverage_width", 3000);
  my $height = getValue($def, "coverage_height", 1500);

  $config->{$task_name} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$task_name",
    docker_prefix         => "",
    #init_command          => "ln -s __FILE__ __NAME__.bam",
    option                => "python3 $python_script -n __NAME__ --width $width --height $height",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $gene_ref,
    source_arg            => "-n",
    parameterSampleFile2_ref => $bam_ref,
    parameterSampleFile2_arg => "-b",
    parameterFile1_ref => $locus_ref,
    parameterFile1_arg => "-l",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => ".coverage.png",
    output_other_ext      => "",
    use_tmp_folder        => 0,
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push( @$tasks, $task_name );
}


sub addLocusCoverage {
  my ($config, $def, $tasks, $target_dir, $task_name, $locus_ref, $bam_ref) = @_;

  my $python_script = dirname(__FILE__) . "/../Visualization/gene_coverage.py";
  my $width = getValue($def, "coverage_width", 3000);
  my $height = getValue($def, "coverage_height", 1500);

  $config->{$task_name} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$task_name",
    docker_prefix         => "",
    #init_command          => "ln -s __FILE__ __NAME__.bam",
    option                => "python3 $python_script -n __NAME__ --width $width --height $height",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $locus_ref,
    source_arg            => "-l",
    parameterSampleFile2_ref => $bam_ref,
    parameterSampleFile2_arg => "-b",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => ".coverage.png",
    output_other_ext      => "",
    use_tmp_folder        => 0,
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push( @$tasks, $task_name );
}

sub addBamsnapLocus {
  my ($config, $def, $tasks, $target_dir, $task_name, $bam_ref) = @_;

  my $option = getValue($def, "bamsnap_option", "-coverage_color 000000 -coverage_vaf 2.0 -draw coordinates bamplot gene -bamplot coverage base read -width 2000 -height 3000 -margin 40");

  $config->{$task_name} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$task_name",
    docker_prefix         => "bamsnap_",
    option                => $option . " -out __NAME__.png",
    interpretor           => "",
    check_program         => 0,
    program               => "bamsnap",
    source                => getValue($def, "bamsnap_locus"),
    source_arg            => "-pos",
    parameterSampleFile2_ref => $bam_ref,
    parameterSampleFile2_arg => "\\\n  -bam",
    parameterSampleFile2_type => "array",
    parameterSampleFile2_join_delimiter => " \\\n       ",
    parameterSampleFile2_name_arg => "\\\n  -title",
    parameterSampleFile2_name_join_delimiter => '" ' . "\\\n         " . '"',
    parameterSampleFile2_name_has_comma => 1,
    output_to_same_folder => 1,
    output_arg            => "-out",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => ".png",
    output_other_ext      => "",
    use_tmp_folder        => 0,
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push( @$tasks, $task_name );

  if (getValue($def, "bamsnap_coverage", 1)){
    $config->{bamsnap_locus_list} = getValue($def, "bamsnap_locus");
    my $coverage_task = "bamsnap_coverage";
    addLocusCoverage($config, $def, $tasks, $target_dir, $coverage_task, "bamsnap_locus_list", $bam_ref);
  }
}

sub addPlotGene {
  my ($config, $def, $tasks, $target_dir, $task_name, $sizeFactorTask, $bed_ref, $bam_ref) = @_;

  my $parameterFile1_key = "parameterFile1_ref";
  if( -e $bed_ref) {
    $parameterFile1_key = "parameterFile1";
  }

  $config->{$task_name} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => $target_dir . "/$task_name",
    option                => getValue($def, "plot_gene_option", ""),
    interpretor           => "python3",
    program               => "../Visualization/plotGene.py",
    parameterFile1_arg    => "-i",
    $parameterFile1_key   => $bed_ref,
    parameterFile3_arg    => "-s",
    parameterFile3_ref    => [$sizeFactorTask, ".sizefactor"],
    parameterSampleFile1_arg => "-b",
    parameterSampleFile1_ref => $bam_ref,
    output_to_result_directory => 1,
    output_arg            => "-o",
    output_file_ext       => ".position.txt",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };

  push( @$tasks, $task_name );
}

sub addSizeFactor {
  my ($config, $def, $tasks, $target_dir, $sizeFactorTask, $bam_ref) = @_;
  $config->{$sizeFactorTask} = {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    #docker_prefix            => "report_",
    target_dir               => $target_dir . '/' . $sizeFactorTask,
    interpretor              => "python3",
    program                  => "../Annotation/getBackgroundCount.py",
    parameterSampleFile1_arg => "-b",
    parameterSampleFile1_ref => $bam_ref,
    output_arg               => "-o",
    output_file_ext          => ".count",
    output_other_ext         => ".count.sizefactor",
    sh_direct                => 1,
    'pbs'                    => {
      'nodes'    => '1:ppn=1',
      'mem'      => '40gb',
      'walltime' => '10'
    },
  };
  push( @$tasks, $sizeFactorTask );
}

sub writeAnnotationLocus_bed {
  my ($locusList, $locusFile) = @_;
  open(my $fh, '>', $locusFile) or die "Could not open file '$locusFile' $!";

  my $locusName;
  if(is_array($locusList)){
    my $count = 0;
    for my $locus (@$locusList){
      $count = $count + 1;
      $locus =~ s/,//g;

      $locusName = $locus;
      $locusName =~ s/:/_/g; 
      $locusName =~ s/-/_/g; 

      my @parts = split /:/, $locus;
      my $chr = $parts[0];
      my $positions = $parts[1];
      my @pos = split /-/, $positions;
      my $start = $pos[0];
      my $end = $pos[1];
      print $fh $chr . "\t" . $start . "\t" . $end . "\t1000\t" . $locusName . "\t+\t" . $locusName . "\n";
    }
    close($fh);
  }elsif(is_hash($locusList)){
    for $locusName (sort keys %$locusList){
      my $locus = $locusList->{$locusName};
      my @parts = split /:/, $locus;
      my $chr = $parts[0];
      my $positions = $parts[1];
      my @pos = split /-/, $positions;
      my $start = $pos[0];
      my $end = $pos[1];
      print $fh $chr . "\t" . $start . "\t" . $end . "\t1000\t" . $locusName . "\t+\t" . $locusName . "\n";
    }
    close($fh);
  }else{
    close($fh);
    die("locus should be either array or hash");
  }
}

sub writeAnnotationLocus_gff {
  my ($locusList, $locusFile) = @_;
  open(my $fh, '>', $locusFile) or die "Could not open file '$locusFile' $!";

  my $locusName;
  if(is_array($locusList)){
    my $count = 0;
    for my $locus (@$locusList){
      $count = $count + 1;
      $locus =~ s/,//g;

      $locusName = $locus;
      $locusName =~ s/:/_/g; 
      $locusName =~ s/-/_/g; 

      my @parts = split /:/, $locus;
      my $chr = $parts[0];
      my $positions = $parts[1];
      my @pos = split /-/, $positions;
      my $start = $pos[0];
      my $end = $pos[1];
      print $fh $chr . "\t" . $locusName . "\tGENE\t".  $start . "\t" . $end . "\t.\t+\t.\n";
    }
    close($fh);
  }elsif(is_hash($locusList)){
    for $locusName (sort keys %$locusList){
      my $locus = $locusList->{$locusName};
      my @parts = split /:/, $locus;
      my $chr = $parts[0];
      my $positions = $parts[1];
      my @pos = split /-/, $positions;
      my $start = $pos[0];
      my $end = $pos[1];
      print $fh $chr . "\t" . $locusName . "\tGENE\t".  $start . "\t" . $end . "\t.\t+\t.\n";
    }
    close($fh);
  }else{
    close($fh);
    die("locus should be either array or hash");
  }
}

sub addAnnotationLocus {
  my ($config, $def, $tasks, $target_dir, $task_name, $sizeFactorTask, $bam_ref) = @_;

  my $locusFile = $target_dir . "/annotation_locus.bed";
  writeAnnotationLocus_bed($def->{annotation_locus}, $locusFile);

  if($def->{perform_bamsnap}){
    my $bamsnap_task = $task_name . "_bamsnap";
    addBamsnap($config, $def, $tasks, $target_dir, $bamsnap_task, $locusFile, $bam_ref);
  }

  addPlotGene($config, $def, $tasks, $target_dir, $task_name, $sizeFactorTask, $locusFile, $bam_ref);
}

sub addAnnotationGenes {
  my ($config, $def, $tasks, $target_dir, $task_name, $sizeFactorTask, $bam_ref) = @_;

  my $genes_str = $def->{annotation_genes};
  my @genes = split /[;, ]+/, $genes_str;
  my %gene_map = map { $_ => 1 } @genes;
  $config->{annotation_genes} = \%gene_map;
  #print(Dumper($config->{annotation_genes}));

  my $geneLocus = addGeneLocus($config, $def, $tasks, $target_dir);

  if($def->{perform_bamsnap}){
    my $bamsnap_task = $task_name . "_bamsnap";
    addBamsnap($config, $def, $tasks, $target_dir, $bamsnap_task, [$geneLocus, "bed"], $bam_ref);
  }

  addPlotGene($config, $def, $tasks, $target_dir, $task_name, $sizeFactorTask, [ $geneLocus, ".bed" ], $bam_ref);
}

sub add_peak_count {
  my ($config, $def, $tasks, $target_dir, $task_name, $callName) = @_;
  $config->{$task_name} = {
    class      => "CQS::ProgramWrapper",
    perform    => 1,
    suffix     => "_pc",
    target_dir => "${target_dir}/" . $task_name,
    interpretor => "python3",
    program    => "../Count/bedCount.py",
    option     => "",
    source_arg => "-i",
    source_ref => [ $callName, ".bed\$" ],
    output_arg => "-o",
    output_prefix => ".txt",
    output_ext => ".txt",
    can_result_be_empty_file => 1,
    sh_direct  => 1,
    pbs        => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "2gb"
    },
  };
  push @$tasks, ($task_name);
  return($task_name);
}

sub add_alignment_summary {
  my ($config, $def, $tasks, $target_dir, $task_name, $rtemplate, $output_file_ext, $read_1_ref, $read_2_ref, $read_3_ref, $read_4_ref, $read_5_ref ) = @_;

  $config->{$task_name} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    rCode                    => "",
    target_dir               => "${target_dir}/" . getNextFolderIndex($def) . ${task_name},
    option                   => "",
    parameterSampleFile1_ref => $read_1_ref,
    parameterSampleFile2_ref => $read_2_ref,
    parameterSampleFile3_ref => $read_3_ref,
    parameterSampleFile4_ref => $read_4_ref,
    parameterSampleFile5_ref => $read_5_ref,
    rtemplate                => $rtemplate,
    output_file              => "",
    output_file_ext          => $output_file_ext,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "5gb"
    },
  };
  push(@$tasks, $task_name);
  return($task_name);
}

sub add_bam_validation {
  my ($config, $def, $tasks, $target_dir, $task_name, $source_ref) = @_;

  my $bam_validation_option = getValue($def, "bam_validation_option", "--IGNORE_WARNINGS --SKIP_MATE_VALIDATION --VALIDATE_INDEX false --INDEX_VALIDATION_STRINGENCY NONE");

  $config->{$task_name} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$task_name",
    option                => "
gatk ValidateSamFile -I __FILE__ -O __NAME__.txt $bam_validation_option
    
status=\$?
if [[ \$status -ne 0 ]]; then
  if [[ -e __NAME__.txt ]]; then
    mv __NAME__.txt __NAME__.failed
  else
    touch __NAME__.failed
  fi
fi

",
    interpretor           => "",
    program               => "",
    check_program         => 0,
    source_arg            => "",
    source_ref            => $source_ref,
    docker_prefix         => "gatk4_",
    output_to_same_folder => 1,
    output_arg            => "-O",
    output_file_prefix    => ".txt",
    output_file_ext       => ".txt",
    sh_direct             => 0,
    can_result_be_empty_file => 1,
    pbs                   => {
      "nodes"    => "1:ppn=1",
      "walltime" => getValue($def, "bam_validation_walltime", "24"),
      "mem"      => getValue($def, "bam_validation_mem", "5gb")
    },
  };

  push(@$tasks, $task_name);
}

sub add_unique_r {
  my ($config, $def, $tasks, $target_dir, $task_name, $rtemplate, $output_file_ext, $ref_array ) = @_;

  $config->{$task_name} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    rCode                    => "",
    target_dir               => "${target_dir}/" . getNextFolderIndex($def) . ${task_name},
    option                   => "",
    rtemplate                => $rtemplate,
    output_file              => "",
    output_file_ext          => $output_file_ext,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "5gb"
    },
  };

  for my $idx (0..scalar(@$ref_array)){
    my $idx2 = $idx+1;
    $config->{$task_name}{"parameterSampleFile${idx2}_ref"} = $ref_array->[$idx];
  }
  push(@$tasks, $task_name);
  return($task_name);
}

sub add_gsea {
  my ($config, $def, $tasks, $target_dir, $gseaTaskName, $rnk_file_ref, $keys, $suffix ) = @_;

  my $gsea_jar        = getValue($def, "gsea_jar");
  my $gsea_db         = getValue($def, "gsea_db");
  my $gsea_categories = getValue($def, "gsea_categories");
  my $rCode = "gseaDb='" . $gsea_db . "'; gseaJar='" . $gsea_jar . "'; gseaCategories=c(" . $gsea_categories . "); makeReport=0;";

  if(defined $def->{"gsea_chip"}){
    my $gsea_chip = getValue($def, "gsea_chip");
    if( ! -e $gsea_chip ){
      $gsea_chip = $gsea_db . "/" . $gsea_chip;
    }
    $rCode = $rCode . "gseaChip='" . $gsea_chip . "';";
  }

  #my $gseaCategories = "'h.all.v6.1.symbols.gmt','c2.all.v6.1.symbols.gmt','c5.all.v6.1.symbols.gmt','c6.all.v6.1.symbols.gmt','c7.all.v6.1.symbols.gmt'";
  $config->{$gseaTaskName} = {
    class                      => "CQS::UniqueR",
    perform                    => 1,
    target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gseaTaskName,
    rtemplate                  => "GSEAPerform.R",
    parameterSampleFile1_ref   => $rnk_file_ref,
    parameterSampleFile2   => {
      task_name => getValue($def, "task_name")
    },
    output_to_result_directory => 1,
    output_perSample_file      => "parameterSampleFile1",
    output_perSample_file_byName => 1,
    output_perSample_file_regex => "(.+?)[_min5_fdr|_GSEA.rnk]",
    output_perSample_file_ext  => ".gsea.csv;.gsea",
    no_docker                  => getValue($def, "gsea_no_docker", 0),
    #has_empty_ext              => 1,
    sh_direct                  => 1,
    rCode                      => $rCode,
    pbs                        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "10gb"
    },
  };
  push( @$tasks, $gseaTaskName );

  my @gsea_report_files = ();
  my @gsea_report_names = ();
  my $pairs = $config->{pairs};
  for my $key (sort keys %$pairs ) {
    push( @gsea_report_files, $gseaTaskName, $key . $suffix . ".gsea.csv" );
    push( @gsea_report_names, "gsea_" . $key );
  }

  my $gsea_report = $gseaTaskName . "_report";
  $config->{$gsea_report} = {
    class                      => "CQS::UniqueR",
    perform                    => 1,
    target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gsea_report,
    rtemplate                  => "GSEAReport.R",
    rReportTemplate            => "GSEAReport.Rmd;../Pipeline/Pipeline.R;Functions.R,reportFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => ".gsea.html",
    parameterSampleFile1_ref   => \@gsea_report_files,
    parameterSampleFile1_names => \@gsea_report_names,
    parameterSampleFile2 => {
      task_name => getValue($def, "task_name")
    },
    remove_empty_parameter => 1,
    output_ext => "gsea_files.csv",
    samplename_in_result => 0,
    sh_direct                  => 1,
    pbs                        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$tasks, $gsea_report );
}

sub add_maf_filter {
  my ($config, $def, $tasks, $target_dir, $maf_filter_name, $source_ref ) = @_;

  my $allele_frequency_percentage = getValue($def, "filter_variants_by_allele_frequency_percentage");
  my $allele_frequency_maf = getValue($def, "filter_variants_by_allele_frequency_maf");
  my $min_inbreeding_coeff = getValue($def, "filter_variants_by_min_inbreeding_coeff");

  $config->{$maf_filter_name} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "${target_dir}/${maf_filter_name}",
    option                => "-p $allele_frequency_percentage -f $allele_frequency_maf --min_inbreeding_coeff $min_inbreeding_coeff",
    interpretor           => "python3",
    program               => "../Annotation/filterVcf.py",
    parameterFile1_arg    => "-i",
    parameterFile1_ref    => $source_ref,
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_ext       => ".maf_filtered.vcf.gz",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$tasks, $maf_filter_name;
}

sub add_bowtie_index {
  my ($config, $def, $tasks, $parent_dir, $bowtie_index_task, $fasta) = @_;
  
  $config->{$bowtie_index_task} = {
    class        => "CQS::ProgramWrapper",
    perform      => 1,
    target_dir   => $parent_dir . "/" . $bowtie_index_task,
    option       => "
bowtie-build $fasta __NAME__

#__FILE__
#__OUTPUT__
",
    source => {
      $def->{task_name} => [$fasta]
    },
    interpretor           => "",
    check_program         => 0,
    program               => "",
    output_result_folder  => 0,
    output_arg            => "",
    output_file_ext       => "",
    output_to_same_folder => 1,
    check_file_ext        => ".4.ebwt",
    sh_direct             => 1,
    pbs          => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "10gb"
    },
  };

  push(@$tasks, $bowtie_index_task);
}

sub add_extract_bam_locus {
  my ($config, $def, $tasks, $target_dir, $task_name, $locus, $bam_ref) = @_;

  $config->{$task_name} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => $target_dir . "/" . getNextFolderIndex($def) . $task_name,
    interpretor => "",
    check_program => 0,
    option => "
samtools view -b -o __OUTPUT__ __FILE__ $locus
samtools index __OUTPUT__
",
    program => "",
    source_ref => $bam_ref,
    output_arg => "",
    output_file_prefix => ".bam",
    output_file_ext => ".bam",
    output_to_same_folder => 1,
    sh_direct   => 1,
    use_tmp_folder => 0,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    }
  };
  push @$tasks, ($task_name);
}

sub has_comparison {
  my $def = shift;
  return ((defined $def->{pairs}) || (defined $def->{pairs_config}));
}

sub add_featurecount {
  my ($config, $def, $tasks, $target_dir, $task_name, $bam_ref, $transcript_gtf, $no_docker) = @_;

  my $featureCountFolder = $target_dir . "/" . getNextFolderIndex($def) . "$task_name";
  $config->{$task_name} = {
    class         => "Count::FeatureCounts",
    perform       => 1,
    target_dir    => $featureCountFolder,
    option        => getValue($def, "featureCounts_option", "-g gene_id -t exon"),
    source_ref    => $bam_ref,
    gff_file      => $transcript_gtf,
    is_paired_end => is_paired_end($def),
    no_docker => $no_docker,
    sh_direct     => 0,
    pbs           => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "40gb"
    },
  };
  $config->{"${task_name}_summary"} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${featureCountFolder}_summary",
    option                   => "",
    rtemplate                => "../Alignment/AlignmentUtils.r,../Alignment/STARFeatureCount.r",
    output_file_ext          => ".FeatureCountSummary.csv;.FeatureCountSummary.csv.png",
    parameterSampleFile2_ref => [ $task_name, ".count.summary" ],
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };

  push @$tasks, $task_name;
  push @$tasks, "${task_name}_summary";
}

sub add_md5 {
  my ($config, $def, $tasks, $target_dir, $source_ref) = @_;
  $config->{md5} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => $target_dir . "/" . getNextFolderIndex($def) . "md5",
    check_program => 0,
    program => "",
    option => "
rm -f tmp.md5

md5sum __FILE__ >> tmp.md5

mv tmp.md5 __OUTPUT__
",
    source_arg => "",
    source_join_delimiter => " >> tmp.md5 \nmd5sum ",
    source_ref => $source_ref,
    output_arg => "",
    output_file_prefix => ".md5",
    output_file_ext => ".md5",
    output_to_same_folder => 1,
    can_result_be_empty_file => 0,
    sh_direct   => 0,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    }
  };
  push @$tasks, ("md5");

  $config->{md5_merge} = {
    class => "CQS::UniqueR",
    target_dir => $target_dir . "/" . getNextFolderIndex($def) . "md5_merge",
    check_program => 0,
    rtemplate                => "../QC/validateMD5merge.r",
    parameterSampleFile1_ref => "md5",
    output_file_ext => ".md5.txt",
    sh_direct   => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    }
  };
  push @$tasks, ("md5_merge");
  return($config);
}

sub add_bamplot {
  my ($config, $def, $tasks, $target_dir, $bam_ref) = @_;
  my $gff_key = "";
  my $gff_value = "";

  if ( not defined $def->{bamplot_gff} ) {
    if($def->{annotation_locus}){
      $gff_key = "gff_file";
      $gff_value = $target_dir . "/annotation_locus.gff";
      print("writing annotation_locus to " . $gff_value . "\n");
      writeAnnotationLocus_gff($def->{annotation_locus}, $gff_value);
    }else{
      $config->{bamplot_gene_gff} = {
        class => "CQS::UniqueR",
        target_dir => $target_dir . "/" . getNextFolderIndex($def) . "bamplot_gene_gff",
        check_program => 0,
        rtemplate => "../CQS/countTableVisFunctions.R,../Annotation/get_gene_locus.r",
        output_file_ext => ".gff",
        parameterSampleFile1 => {
          task_name => getValue($def, "task_name"),
          email => getValue($def, "email"),

          biomart_host => getValue($def, "biomart_host"),
          biomart_dataset   => getValue($def, "biomart_dataset"),
          biomart_symbolKey => getValue($def, "biomart_symbolKey"),
          biomart_add_chr => getValue($def, "biomart_add_chr"),
          biomart_add_prefix => getValue($def, "biomart_add_prefix", ""),

          gene_names => getValue($def, "gene_names"),
          gene_shift => getValue($def, "gene_shift", 0),

          output_gff => 1,
        },
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "2",
          "mem"       => "10gb"
        }
      };
      $gff_key = "gff_file_ref";
      $gff_value = "bamplot_gene_gff";
      push( @$tasks, "bamplot_gene_gff" );
    }
  }else{
    $gff_key = "gff_file";
    $gff_value = getValue($def, "bamplot_gff");
  }

  $config->{bamplot} = {
    class              => "Visualization::Bamplot",
    perform            => 1,
    target_dir         => $target_dir . "/" . getNextFolderIndex($def) . "bamplot",
    option             => "-g " . getValue($def, "dataset_name") . " -y uniform -r --save-temp",
    source_ref         => $bam_ref,
    $gff_key           => $gff_value,
    is_rainbow_color   => 0,
    is_single_pdf      => 0,
    is_draw_individual => 0,
    groups             => $def->{"plotgroups"},
    colors             => $def->{"colormaps"},
    docker_prefix => "bamplot_",
    sh_direct          => 1,
    pbs                => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "10gb"
    },
  };
  push( @$tasks, "bamplot" );

  $config->{bamplot_html} = {
    class              => "CQS::UniqueRmd",
    perform            => 1,
    target_dir         => $target_dir . "/" . getNextFolderIndex($def) . "bamplot",
    report_rmd_file => "../Visualization/bamplot.Rmd",
    additional_rmd_files => "../CQS/reportFunctions.R",
    option => "",
    parameterSampleFile1 => {
      task_name => getValue($def, "task_name"),
      email => getValue($def, "email"),
      affiliation => getValue($def, "affiliation", ""),
      gene_shift => getValue($def, "gene_shift", 0),
    },
    suffix => ".report",
    output_file_ext => ".report.html",
    can_result_be_empty_file => 0,
    sh_direct   => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };
  push( @$tasks, "bamplot_html" );
};

sub add_fastq_screen {
  my ($config, $def, $tasks, $parent_dir, $task_name, $source_ref) = @_;

  my $conf = "";
  if($def->{fastq_screen_configuration_file}){
     $conf = "--conf " . $def->{fastq_screen_configuration_file};
  }

  $config->{$task_name} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => $parent_dir . "/" . getNextFolderIndex($def) . $task_name,
    use_tmp_folder => 0,
    suffix  => "_fc",
    interpretor => "",
    program => "",
    check_program => 0,
    option => "
fastq_screen $conf --outdir . --threads 8 __FILE__

#__OUTPUT__
",
    source_arg => "",
    source_join_delimiter => " ",
    source_ref => $source_ref,
    output_arg => "",
    output_by_file => 1,
    output_by_file_remove_pattern => ".fastq.gz",
    output_file_prefix => "_screen.txt",
    output_file_ext => "_screen.txt",
    output_to_same_folder => 1,
    can_result_be_empty_file => 1,
    sh_direct   => 0,
    no_docker => 1,
    pbs => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "20gb"
    }
  };
  push(@$tasks, $task_name);
}

1;

