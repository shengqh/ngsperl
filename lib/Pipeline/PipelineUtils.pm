#!/usr/bin/perl
package Pipeline::PipelineUtils;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(getValue 
    getIntermidiateDir
    getIndexName
    initPipelineOptions 
    readChromosomeFromDictFile
    addPreprocess 
    addFastQC 
    addBlastn 
    addBowtie 
    addPARalyzer
    addBowtie1PARalyzer 
    addBamStat 
    addOutputOption 
    getOutputFormat
    getDEseq2TaskName 
    addDEseq2 
    addDeseq2Visualization 
    addDeseq2SignificantSequenceBlastn
    getBatchGroups 
    addHomerMotif 
    addHomerAnnotation 
    addEnhancer 
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
    addGATK4CNVGermlineCohortAnalysis 
    addXHMM
    addGeneLocus
    annotateNearestGene
    checkFileGroupPairNames
    addStarFeaturecount)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getValue {
  my ( $def, $name, $defaultValue ) = @_;
  if ( defined $def->{$name} ) {
    return $def->{$name};
  }
  elsif ( defined $defaultValue ) {
    return $defaultValue;
  }
  else {
    die "Define $name in user definition first.";
  }
}

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
    cluster    => $def->{cluster},
    sh_direct  => 1,
    pbs        => {
#      "email"    => $def->{email},
      "nodes"    => "1:ppn=" . $curThread,
      "walltime" => "4",
      "mem"      => "4gb"
    },
  };

  my $summaryTask = $fastqcTask . "_summary";

  $config->{$summaryTask} = {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => $parentDir . "/" . getNextFolderIndex($def) . $summaryTask,
    option     => "",
    cluster    => $def->{cluster},
    source_ref => [$fastqcTask, "data.txt"],
    sh_direct  => 1,
    can_result_be_empty_file => 1,
    pbs        => {
#      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  };
  push @$individual, $fastqcTask;
  push @$summary,    $summaryTask;
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
    bowtie1_index         => $bowtieIndex,
    samonly               => 0,
    sh_direct             => 0,
    mappedonly            => 1,
    export_max_mapped     => $def->{export_max_mapped},
    cluster               => $def->{cluster},
    output_to_same_folder => $def->{bowtie1_output_to_same_folder},
    pbs                   => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $hours,
      "mem"       => "40gb"
    },
  };

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
  if ( defined $libraryKey ) {
    $result = $result . "_" . $libraryKey;
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
    if ( getValue( $def, $key, $defaultValue ) ) {
      $result = $result . "$newkey<-TRUE;";
    }
    else {
      $result = $result . "$newkey<-FALSE;";
    }
  }
  return ($result);
}

sub getOutputFormat {
  my ( $def, $rcode ) = @_;
  my $result = $rcode;

  $result = addOutputOption( $def, $result, "DE_outputPdf",         0,                          "outputPdf" );
  $result = addOutputOption( $def, $result, "DE_outputPng",         1,                          "outputPng" );
  $result = addOutputOption( $def, $result, "DE_outputTIFF",        0,                          "outputTIFF" );
  $result = addOutputOption( $def, $result, "DE_showVolcanoLegend", 1,                          "showVolcanoLegend" );
  $result = addOutputOption( $def, $result, "use_pearson_in_hca",   $def->{use_pearson_in_hca}, "usePearsonInHCA" );
  $result = addOutputOption( $def, $result, "showLabelInPCA",       1 );

  return ($result);
}

sub addDEseq2 {
  my ( $config, $def, $summary, $taskKey, $countfileRef, $deseq2Dir, $DE_min_median_read, $libraryFile, $libraryKey ) = @_;

  my $taskName = getDEseq2TaskName( $taskKey, $libraryKey, $def );

  my $libraryFileKey = "library_file";
  if ( ref($libraryFile) eq 'ARRAY' ) {
    $libraryFileKey = "library_file_ref";
  }

  my $rCode = getOutputFormat( $def, getValue($def, "DE_rCode", "") );
  $rCode = addOutputOption( $def, $rCode, "top25cv_in_hca", $def->{top25cv_in_hca}, "top25cvInHCA" );

  my $raw_rCode = $rCode;
  my $filterBaseMean = (defined $def->{filterBaseMean}) && $def->{filterBaseMean};
  if($filterBaseMean){
    #die("filterBaseMean");
    $rCode = $rCode . "filterBaseMean=1;filterBaseMeanValue=" . getValue($def, "filterBaseMeanValue", 30) . ";";
  }

  $config->{$taskName} = {
    perform                      => 1,
    target_dir                   => $deseq2Dir . "/" . getNextFolderIndex($def) . "$taskName",
    output_to_dir                => getReportDir($def),
    option                       => "",
    sh_direct                    => 1,
    show_label_PCA               => $def->{show_label_PCA},
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
    $libraryFileKey              => $libraryFile,
    library_key                  => $libraryKey,
    rCode                        => $rCode,
    pbs                          => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => "10",
      "mem"       => "20gb"
    },
  };

  if (defined $def->{pairs_config}) {
    $config->{$taskName}{class} = "Comparison::DESeq2config";
    $config->{$taskName}{source} = $def->{pairs_config};
  }else{
    my $groupNames = defined $def->{deseq2_groups} ? "deseq2_groups" : "groups";
    $config->{$taskName}{source_ref} = "pairs";
    $config->{$taskName}{groups_ref} = $groupNames;
    if(defined $def->{covariance_file}){
      $config->{$taskName}{class} = "Comparison::DESeq2covariance";
      $config->{$taskName}{covariance_file} = $def->{covariance_file};
    }else{
      $config->{$taskName}{class} = "Comparison::DESeq2";
    }
  }

  if ( ref($countfileRef) eq "ARRAY" ) {
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
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
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

sub addHomerAnnotation {
  my ( $config, $def, $summary, $target_dir, $callName, $callFilePattern ) = @_;
  my $homerName = $callName . "_homer_annotation";
  $config->{$homerName} = {
    class        => "Homer::Annotation",
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
    output_file_ext          => ".log.tss.png;.log.distal.png;.tss.tsv;.distal.tsv",
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
    rtemplate                => "countTableVisFunctions.R,countTableGroupCorrelation.R",
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
    output_to_dir => getReportDir($def),
    source_ref    => $multiqc_depedents,
    root_dir      => $root_dir,
    sh_direct     => 1,
    pbs           => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
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
    is_paired_end           => $pairend,
    is_sorted_by_coordinate => 1,
    sh_direct               => 0,
    pbs                     => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "48",
      "mem"      => "40gb"
    },
  };
  push @$individual, $task_name;
}

sub writeDesignTable {
  my ( $target_dir, $section, $designtable, $bamfiles, $peaksfiles, $peakSoftware, $merged, $task_name, $treatments, $controls ) = @_;

  my $defaultTissue = getValue( $designtable, "Tissue", "" );
  my $defaultFactor = getValue( $designtable, "Factor", "" );

  my $result = {};

  if ($merged) {
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
        my $peakFile  = $peaksfiles->{$sampleName}->[0];

        my $sampleId   = $treatments->{$sampleName}->[0];
        my $bamReads   = $bamfiles->{$sampleId}[0];
        my $controlId  = "";
        my $bamControl = "";

        if ( defined $controls && defined $controls->{$sampleName} ) {
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
      }
    }
    close($map);

    $result->{$task_name} = $mapfile;
  }
  else {
    for my $name ( sort keys %$designtable ) {
      if ( $name eq "Tissue" || $name eq "Factor" ) {
        next;
      }

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
    }
  }

  return $result;
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
  my ( $config, $def, $summary, $target_dir, $source_name, $source_pattern, $prefix, $indexDic, $indexKey ) = @_;
  if (not defined $prefix){
    $prefix = $source_name;
  }

  my $annovar_name = getIndexName($prefix, "_annovar", $indexDic, $indexKey);
  my $source_ref = ( defined($source_pattern) and ( $source_pattern ne "" ) ) ? [ $source_name, $source_pattern ] : $source_name;
  $config->{$annovar_name} = {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/$annovar_name",
    source_ref => $source_ref,
    option     => getValue( $def, "annovar_param" ),
    annovar_db => getValue( $def, "annovar_db" ),
    buildver   => getValue( $def, "annovar_buildver" ),
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
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
  $config->{$annovar_to_maf_report} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . $annovar_to_maf_report,
    rtemplate                => "../Annotation/mafReport.r",
    output_file              => "parameterSampleFile1",
    output_file_ext          => ".report.html",
    parameterSampleFile1_ref => [ $annovar_to_maf, ".tsv.maf\$" ],
    parameterFile1           => $def->{family_info_file},
    sh_direct                => 1,
    rCode                    => ( defined $def->{family_info_file} ? "clinicalFeatures=\"" . $def->{family_info_feature} . "\";" : "" ),
#            rCode                    => ( defined $def->{family_info_file} ? "clinicalFeatures=\"" . $def->{family_info_feature} . "\";" : "" )
#              . ( defined $def->{annotation_genes} ? "interestedGeneStr=\"" . $def->{annotation_genes} . "\"" : "" ),
    pbs => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "24",
      "mem"       => "10gb"
    },
  };
  push @$summary, $annovar_to_maf_report;
}

sub addFilterMafAndReport {
  my ( $config, $def, $summary, $target_dir, $mutect2call ) = @_;

#  my $mutect2_index_dic = {};
#  my $mutect2_index_key = "mafReport_Index";
#  my $taskName = $mutect2call . getNextIndex($mutect2_index_dic, $mutect2_index_key) . "_mergeAndMafreport";
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


sub addGATK4PreprocessIntervals {
  my ( $config, $def, $target_dir, $bam_ref, $prefix, $step1, $step2, $step3, $step4, $step5, $step6, $index ) = @_;
  if (!defined $index){
    $index = "";
  }

  my $result = $prefix . "_gatk4_CNV_Germline${index}_PreprocessIntervals";
  if ( !defined $config->{$result} ) {
    #PreprocessIntervals at summary level
    $config->{$result} = {
      class             => "GATK4::PreprocessIntervals",
      option            => "",
      interval_file     => getValue( $def, "covered_bed" ),
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

  my $preprocessIntervalsTask = addGATK4PreprocessIntervals( $config, $def, $target_dir, $bam_ref, $prefix, $step1, $step2, $step3, $step4, $step5, $step6, "_01" );

  my $chrCode = getValue($def, "has_chr_in_chromosome_name") ? ";addChr=1" : "";

  #CollectReadCounts at sample level
  my $CollectReadCounts = $prefix . "_gatk4_CNV_Germline_02_CollectReadCounts";
  $config->{$CollectReadCounts} = {
    class                      => "GATK4::CollectReadCounts",
    source_ref                 => $bam_ref,
    option                     => "",
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
  my $FilterIntervals = $prefix . "_gatk4_CNV_Germline_03_FilterIntervals";
  $config->{$FilterIntervals} = {
    class                      => "GATK4::FilterIntervals",
    source_ref                 => $CollectReadCounts,
    option                     => "",
    preprocessed_intervals_ref => $preprocessIntervalsTask,
    ref_fasta_dict             => getValue( $def, "ref_fasta_dict" ),
    ref_fasta                  => getValue( $def, "ref_fasta" ),
    blacklist_file             => $def->{blacklist_file},
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
  my $DetermineGermlineContigPloidyCohortMode = $prefix . "_gatk4_CNV_Germline_04_DetermineGermlineContigPloidyCohortMode";
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
  my $GermlineCNVCaller = $prefix . "_gatk4_CNV_Germline_05_GermlineCNVCaller";
  if ($def->{gatk4_cnv_by_scatter}){
    #scatter filter intervals
    my $ScatterIntervals = $GermlineCNVCaller . "_1_scatterIntervals";
    $config->{$ScatterIntervals} = {
      class       => "CQS::ProgramWrapperOneToMany",
      perform     => 1,
      target_dir  => "$target_dir/$ScatterIntervals",
      option      => "",
      interpretor => "python",
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

  #PostprocessGermlineCNVCalls at sample level
  my $PostprocessGermlineCNVCalls = $prefix . "_gatk4_CNV_Germline_06_PostprocessGermlineCNVCalls";
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

  #CombineGCNV at summary level
  my $CombineGCNV = $prefix . "_gatk4_CNV_Germline_07_CombineGCNV";
  $config->{$CombineGCNV} = {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    target_dir               => $target_dir . '/' . $CombineGCNV,
    interpretor              => "python",
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
  
  my $sizeFactorTask = $prefix . "_gatk4_CNV_Germline_08_SizeFactor";
  $config->{$sizeFactorTask} = {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    target_dir               => $target_dir . '/' . $sizeFactorTask,
    interpretor              => "python",
    program                  => "../GATK4/getBackgroundCount.py",
    parameterSampleFile1_arg => "-b",
    parameterSampleFile1_ref => $bam_ref,
    parameterFile1_arg       => "-i",
    parameterFile1           => getValue( $def, "covered_bed" ),
    parameterFile2_arg       => "-c",
    parameterFile2_ref       => [ $CombineGCNV ],
    output_arg               => "-o",
    output_file_ext          => ".txt.sizefactor;.txt",
    sh_direct                => 1,
    'pbs'                    => {
      'nodes'    => '1:ppn=1',
      'mem'      => '40gb',
      'walltime' => '10'
    },
  };
  push( @$step6, $sizeFactorTask );

  my $cnvIndex = "09";
  if($def->{plotCNVGenes} && $def->{annotation_genes}){
    my $cnvGenes = $prefix . "_gatk4_CNV_Germline_09_CNVGenesLocus";
    $config->{$cnvGenes} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $target_dir . '/' . $cnvGenes,
      rtemplate                => "../Annotation/getCNVGeneLocus.r",
      rCode                    => "host=\"" . getValue($def, "biomart_host") . "\";dataset=\"" . getValue($def, "biomart_dataset") . "\";symbolKey=\"" . getValue($def, "biomart_symbolKey") . "\"" . $chrCode,
      parameterFile1_ref       => [$CombineGCNV, ".txt"],
      output_file_ext          => ".bed;.missing",
      sh_direct                => 1,
      'pbs'                    => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '10'
      },
    };
    push( @$step6, $cnvGenes );
    
    my $plotCNVgenes = $prefix . "_gatk4_CNV_Germline_10_CNVGenesPlot";
    $config->{$plotCNVgenes} = {
      class                 => "CQS::ProgramWrapper",
      perform               => 1,
      target_dir            => $def->{target_dir} . "/$plotCNVgenes",
      option                => "",
      interpretor           => "python",
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
  }
  
  my $annotationGenesPlot = undef;
  if(defined $def->{annotation_genes} && defined $config->{"annotation_genes_locus"}){
    $annotationGenesPlot = $prefix . "_gatk4_CNV_Germline_" . $cnvIndex . "_AnnotationGenesPlot";
    $config->{$annotationGenesPlot} = {
      class                 => "CQS::ProgramWrapper",
      perform               => 1,
      target_dir            => $def->{target_dir} . "/$annotationGenesPlot",
      option                => "",
      interpretor           => "python",
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
      output_file_ext       => ".position.txt.slim;.position.txt",
      sh_direct             => 1,
      pbs                   => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };

    push( @$step6, $annotationGenesPlot );
  }
  
  return($annotationGenesPlot);
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

sub addGeneLocus {
  my ($config, $def, $summary, $target_dir) = @_;
  my $geneLocus = undef;
  if ( defined $def->{annotation_genes} ) {
    $geneLocus = "annotation_genes_locus";
    $config->{$geneLocus} = {
      class      => "CQS::UniqueR",
      perform    => 1,
      target_dir => $target_dir . '/' . $geneLocus,
      rtemplate  => "../Annotation/getGeneLocus.r",
      rCode      => "host=\""
        . getValue( $def, "biomart_host" )
        . "\";dataset=\""
        . getValue( $def, "biomart_dataset" )
        . "\";symbolKey=\""
        . getValue( $def, "biomart_symbolKey" )
        . "\";genesStr=\""
        . getValue( $def, "annotation_genes" ) 
        . "\";shift=" . getValue( $def, "annotation_genes_shift", 0),
      output_file_ext => ".bed;.missing",
      sh_direct       => 1,
      'pbs'           => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '10'
      },
    };
    push( @$summary, $geneLocus );
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
  my ($def, $groupKeys, $pairKeys, $fileKey) = @_;
  if(!defined $fileKey){
    $fileKey = "files";
  }
  my $files = getValue($def, $fileKey);
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
        if ($groupName =~ /^./) {
          next;
        }
        my $sampleNames = $groups->{$groupName};
        for my $sampleName (@$sampleNames){
          if (!defined $files->{$sampleName}){
            print STDERR "Sample $sampleName in $groupKey $groupName is not defined in files.\n";
            $bFailed = 1;
          }
        }

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
          if (ref($groupNames) eq 'ARRAY') {
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
    output_to_same_folder     => $def->{output_bam_to_same_folder},
    featureCount_option       => getValue( $def, "featureCount_option" ),
    star_location             => $def->{star_location},
    gff_file                  => $transcript_gtf,
    is_paired_end             => is_paired_end($def),
    delete_star_featureCount_bam => $def->{delete_star_featureCount_bam},
    sh_direct                 => 0,
    pbs                       => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => "$star_featurecount_walltime",
      "mem"       => "${star_memory}gb"
    },
  };
  push @$individual, ($star_task);
  
  my $summary_task = "star_featurecount${suffix}_summary";
  $config->{$summary_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $starFolder . "_summary",
    option                   => "",
    rtemplate                => "../Alignment/STARFeatureCount.r",
    output_file_ext          => ".FeatureCountSummary.csv;.FeatureCountSummary.csv.png;.STARSummary.csv;.STARSummary.csv.png",
    parameterSampleFile1_ref => [ $star_task, "_Log.final.out" ],
    parameterSampleFile2_ref => [ $star_task, ".count.summary" ],
    sh_direct                => 1,
    pbs                      => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };
  push @$summary, ($summary_task);

  return($star_task);
};

1;
