#!/usr/bin/perl
package Pipeline::ATACSeq;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::Preprocession;
use Pipeline::PipelineUtils;
use Pipeline::PeakPipelineUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performATACSeq performATACSeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "cluster",               "slurm" );
  initDefaultValue( $def, "max_thread",            8 );
  initDefaultValue( $def, "sequencetask_run_time", 12 );

  #Tasks
  initDefaultValue( $def, "sra_to_fastq", 0 );

  initDefaultValue( $def, "fastq_remove_N", 0 );

  initDefaultValue( $def, "minimum_maq",         30 );
  initDefaultValue( $def, "minimum_insert_size", 30 );
  initDefaultValue( $def, "maximum_insert_size", 1000 );

  initDefaultValue( $def, "perform_cutadapt", 0 );
  if ( getValue( $def, "perform_cutadapt" ) ) {
    #initDefaultValue( $def, "adapter", "CTGTCTCTTATA" );
    initDefaultValue( $def, "min_read_length", 30 );
    initDefaultValue( $def, "cutadapt_option", "-m " . $def->{min_read_length} );
    initDefaultValue( $def, "trim_polyA", 0 );
  }

  initDefaultValue( $def, "perform_rose",          0 );
  initDefaultValue( $def, "perform_coltron",       0 );
  initDefaultValue( $def, "perform_diffbind",      0 );
  initDefaultValue( $def, "annotate_nearest_gene", 0 );
  initDefaultValue( $def, "perform_enhancer",      0 );

  initDefaultValue( $def, "caller_type", "macs2" );
  if ( ( getValue( $def, "caller_type" ) eq "macs2" ) || ( getValue( $def, "caller_type" ) eq "both" ) ) {
    initDefaultValue( $def, "macs2_peak_type",             "board" );
    initDefaultValue( $def, "macs2_callpeak_as_singleend", 0 );
  }

  initDefaultValue( $def, "perform_homer_motifs", 0 );
  if ( getValue( $def, "perform_homer_motifs" ) ) {
    initDefaultValue( $def, "homer_option", "" );
  }

  initDefaultValue( $def, "perform_bamplot", 0 );

  if ( !defined $def->{"treatments"} ) {
    my $files = getValue( $def, "files" );
    my $groups = {};
    for my $sample_name ( sort keys %$files ) {
      $groups->{$sample_name} = [$sample_name];
    }
    $def->{"treatments"} = $groups;
  }

  initDefaultValue( $def, "perform_report", 1 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig( $def, 1 );

  $config->{"treatments"} = $def->{"treatments"};

  my $task    = $def->{task_name};

  my $email   = getValue( $def, "email" );
  my $pairend = is_paired_end( $def );

  my $perform_rose = getValue( $def, "perform_rose" );
  my $perform_coltron = 0;
  if ($perform_rose) {
    $perform_coltron = getValue( $def, "perform_coltron" );
  }

  my $gene_bed;
  my $annotate_nearest_gene = getValue( $def, "annotate_nearest_gene" );
  if ($annotate_nearest_gene) {
    $gene_bed = getValue( $def, "gene_bed" );
  }

  my $perform_diffbind = getValue( $def, "perform_diffbind" );
  if ($perform_diffbind) {
    getValue( $def, "design_table" );
  }

  $def->{aligner} = "bwa";
  $config->{"bwa"} = {
    class              => "Alignment::BWA",
    perform            => 1,
    target_dir         => "${target_dir}/bwa",
    option             => "",
    bwa_index          => getValue( $def, "bwa_fasta" ),
    picard_jar         => getValue( $def, "picard_jar" ),
    source_ref         => $source_ref,
    sort_by_coordinate => 1,
    mark_duplicates    => 0,
    sh_direct          => 0,
    pbs                => {
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  };
  push @$individual, "bwa";
  add_alignment_summary($config, $def, $summary, $target_dir, "bwa_summary", "../Alignment/AlignmentUtils.r;../Samtools/BamStat.r", ".reads.csv;.reads.png;.chromosome.csv;.chromosome.png", [ "bwa", "bam.stat" ], ["bwa", ".chromosome.count"] );

  $config->{bwa_insertsize} = {
    class                    => "CQS::UniqueR",
    target_dir               => $target_dir . "/bwa_insertsize",
    perform                  => 1,
    rtemplate                => "../Visualization/insertSize.r",
    output_file              => ".insertsize.png",
    sh_direct                => 1,
    parameterSampleFile1_ref => [ "bwa", ".bam\$" ],
    cluster                  => $def->{cluster},
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$summary, "bwa_insertsize";

  my $cleanbam_option = $pairend ? "-f 3 -F 3852" : "-F 3844";
  my $cleanbam_walltime = getValue( $def, "bwa_cleanbam_walltime", "24");
  $config->{"bwa_cleanbam"} = {
    class                   => "ATACseq::CleanBam",
    perform                 => 1,
    target_dir              => "${target_dir}/bwa_cleanbam",
    option                  => $cleanbam_option,
    source_ref              => "bwa",
    picard_jar              => getValue( $def, "picard_jar" ),
    remove_chromosome       => "M",
    keep_chromosome         => "chr",
    minimum_maq             => getValue( $def, "minimum_maq" ),
    minimum_insert_size     => getValue( $def, "minimum_insert_size" ),
    maximum_insert_size     => getValue( $def, "maximum_insert_size" ),
    blacklist_file          => $def->{blacklist_file},
    is_sorted_by_coordinate => 1,
    mark_duplicates         => 1,
    sh_direct               => 0,
    pbs                     => {
      "nodes"    => "1:ppn=1",
      "walltime" => "$cleanbam_walltime",
      "mem"      => "40gb"
    },
  };
  push @$individual, "bwa_cleanbam";
  add_alignment_summary($config, $def, $summary, $target_dir, "bwa_cleanbam_summary", "../Alignment/AlignmentUtils.r;../Samtools/BamStat.r", ".chromosome.csv;.chromosome.png", undef, ["bwa_cleanbam", ".chromosome.count"] );

  my $bam_ref = ["bwa_cleanbam", ".bam\$"];
  if($def->{perform_bamsnap} && $def->{"bamsnap_locus"}){
    addBamsnapLocus($config, $def, $summary, $target_dir, "bamsnap_locus", $bam_ref);
  }

  if(defined $def->{annotation_locus} or defined $def->{annotation_genes}){
    my $sizeFactorTask = "size_factor";
    addSizeFactor($config, $def, $summary, $target_dir, $sizeFactorTask, $bam_ref);

    if(defined $def->{annotation_locus}){
      addAnnotationLocus($config, $def, $summary, $target_dir, "annotation_locus_plot", $sizeFactorTask, $bam_ref);
    }

    if(defined $def->{annotation_genes}){
      addAnnotationGenes($config, $def, $summary, $target_dir, "annotation_genes_plot", $sizeFactorTask, $bam_ref);
    }
  }

  my $callerType = getValue( $def, "caller_type" );
  my @callers = ();
  if ( $callerType eq "both" ) {
    push( @callers, ( "macs", "macs2" ) );
  }
  elsif ( $callerType eq "macs" ) {
    push( @callers, ("macs") );
  }
  else {
    push( @callers, ("macs2") );
  }

  my $task_dic = {
  };

  for my $caller (@callers) {
    my @peakTypes = ();
    my $callOption;
    my $callClass;
    my $callerName;

    if ( $caller eq "macs" ) {
      $callClass  = "Chipseq::MACS";
      $callerName = "macs1callpeak";
      push( @peakTypes, "macs" );
      if ( defined $def->{"macs_option"} ) {
        $callOption = $def->{"macs_option"};
      }
      else {
        $callOption = "-p 1e-9 -w -S --space=50";
      }
    }
    else {
      $callClass  = "Chipseq::MACS2Callpeak";
      $callerName = "macs2callpeak";
      if ( defined $def->{"macs2call_option"} ) {
        $callOption = $def->{"macs2call_option"};
      }
      else {
        my $default_macs2call_option = $def->{macs2_callpeak_as_singleend} ? "-f BAM" : $pairend ? "-f BAMPE" : "-f BAM";
        $callOption = $default_macs2call_option . " -g " . getValue( $def, "macs2_genome" ) . " -B -q 0.01 --nomodel --slocal 20000 --llocal 20000 --keep-dup all";
      }

      my $macs2_peak_type = getValue( $def, "macs2_peak_type" );
      if ( $macs2_peak_type eq "both" ) {
        push( @peakTypes, ( "broad", "narrow" ) );
      }
      elsif ( $macs2_peak_type eq "narrow" ) {
        push( @peakTypes, ("narrow") );
      }
      else {
        push( @peakTypes, ("broad") );
      }
    }

    for my $peakType (@peakTypes) {
      my $curCallOption = $callOption;
      my $callFilePattern;
      my $callName = "bwa_${callerName}_${peakType}";
      if ( $peakType eq "macs" ) {
        $callName        = "bwa_${callerName}";
        $callFilePattern = ".name.bed\$";
      }
      elsif ( $peakType eq "broad" ) {
        $curCallOption   = " --broad --broad-cutoff 0.01 " . $curCallOption;
        $callFilePattern = "broadPeak.bed\$";
      }
      else {
        $callFilePattern = "narrowPeak.bed\$";
        if ( $def->{macs2_callpeak_as_singleend} and $pairend ) {
          $callName = $callName . "_CallAsSingleEnd";
        }
      }

      $config->{$callName} = {
        class      => $callClass,
        perform    => 1,
        target_dir => "${target_dir}/" . $callName,
        option     => $curCallOption,
        source_ref => [ "bwa_cleanbam", ".bam\$" ],
        groups_ref  => "treatments",
        sh_direct  => 0,
        pbs        => {
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
      push @$individual, ($callName);

      add_peak_count($config, $def, $summary, $target_dir, $callName . "_count", $callName);

      if ($annotate_nearest_gene) {
        my $geneName = $callName . "_gene";
        $config->{$geneName} = {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => "${target_dir}/${callName}",
          rtemplate                => "../Annotation/findNearestGene.r",
          output_file              => "",
          output_file_ext          => ".Category.Table.csv",
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

      $task_dic->{peak_caller} = $callName;
      $task_dic->{peak_count} = $callName . "_count";

      my $chipqc_taskname = $callName . "_chipqc";
      if ($def->{perform_chipqc}) {
        add_chipqc($config, $def, $summary, $target_dir, $chipqc_taskname, [ $callName, $callFilePattern ], [ "bwa_cleanbam", ".bam\$" ]);
        $task_dic->{chipqc} = $chipqc_taskname;
      }

      if ( getValue( $def, "perform_homer" ) ) {
        my $homer = addHomerAnnotation( $config, $def, $summary, $target_dir, $callName, $callFilePattern );
        $task_dic->{homer} = $homer;
      }

      if ( getValue( $def, "perform_enhancer" ) ) {
        addEnhancer( $config, $def, $individual, $summary, $target_dir, $callName . "_enhancer", [ "bwa_cleanbam", ".bam\$" ], [$callName, $callFilePattern] );
      }

      if ($perform_diffbind) {
        my $bindName = $callName . "_diffbind";
        addDiffbind($config, $def, $summary, $target_dir, $bindName, "bwa_cleanbam", [ $callName, $callFilePattern ]);
        $task_dic->{diffbind} = $bindName;
        if ( getValue( $def, "perform_homer" ) ) {
          my $diffbind_homer = addHomerAnnotation( $config, $def, $summary, $target_dir, $bindName, ".sig.bed" );
          $task_dic->{diffbind_homer} = $diffbind_homer;
        }
      }

      if ( getValue( $def, "perform_homer_motifs" ) ) {
        addHomerMotif( $config, $def, $summary, $target_dir, $callName, $callFilePattern );
      }

      if ($perform_rose) {
        my $roseName = $callName . "_rose";
        $config->{$roseName} = {
          class                => "Chipseq::Rose2",
          perform              => 1,
          target_dir           => "${target_dir}/${roseName}",
          option               => "",
          source_ref           => "bwa_cleanbam",
          groups_ref           => "treatments",
          pipeline_dir         => getValue( $def, "rose_folder" ),
          genome               => getValue( $def, "rose_genome" ),
          binding_site_bed_ref => [ $callName, $callFilePattern ],
          sh_direct            => 1,
          pbs                  => {
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "40gb"
          },
        };
        push @$summary, ($roseName);

        if ($perform_coltron) {
          my $coltronName = $roseName . "_coltron";
          $config->{$coltronName} = {
            class              => "Chipseq::Coltron",
            perform            => 1,
            target_dir         => "${target_dir}/${coltronName}",
            option             => "",
            source_ref         => "bwa_cleanbam",
            groups_ref         => "treatments",
            enhancer_files_ref => [ $roseName, "_AllEnhancers.table.txt" ],
            genome             => getValue( $def, "coltron_genome" ),
            pipeline_dir       => getValue( $def, "rose_folder" ),
            sh_direct          => 1,
            pbs                => {
              "nodes"    => "1:ppn=1",
              "walltime" => "72",
              "mem"      => "40gb"
            },
          };
          push @$summary, ($coltronName);
        }
      }

      if ( ( $caller eq "macs2" ) && ( defined $def->{comparison} ) ) {
        my $peakTask = $callName;
        if ( defined $def->{replicates} ) {
          my $repCallName = $callName . "_replicates";
          $config->{$repCallName} = {
            class      => $callName,
            perform    => 1,
            target_dir => "${target_dir}/${repCallName}",
            option     => $curCallOption,
            source_ref => "bwa_cleanbam",
            groups     => $def->{"replicates"},
            sh_direct  => 0,
            pbs        => {
              "nodes"    => "1:ppn=1",
              "walltime" => "72",
              "mem"      => "40gb"
            },
          };
          push @$summary, ($repCallName);
          $peakTask = $repCallName;

          if ( getValue( $def, "perform_homer_motifs" ) ) {
            addHomerMotif( $config, $def, $summary, $target_dir, $repCallName, $callFilePattern );
          }
        }

        my $diffName = $peakTask . "_bdgdiff";
        $config->{$diffName} = {
          class      => "Chipseq::MACS2Bdgdiff",
          perform    => 1,
          target_dir => "${target_dir}/$diffName",
          option     => "",
          source_ref => $peakTask,
          groups     => $def->{comparison},
          sh_direct  => 0,
          pbs        => {
            "nodes"    => "1:ppn=1",
            "walltime" => "72",
            "mem"      => "40gb"
          },
        };
        push @$summary, ($diffName);

        if ($annotate_nearest_gene) {
          my $geneName = $diffName . "_gene";
          $config->{$geneName} = {
            class                    => "CQS::UniqueR",
            perform                  => 1,
            target_dir               => "${target_dir}/${geneName}",
            rtemplate                => "../Annotation/findNearestGene.r",
            output_file              => "",
            output_file_ext          => ".Category.Table.csv",
            parameterSampleFile1_ref => $diffName,
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

        if ( getValue( $def, "perform_homer_motifs" ) ) {
          addHomerMotif( $config, $def, $summary, $target_dir, $diffName, ".bed\$" );
        }
      }
    }
  }

  if ( getValue( $def, "perform_bamplot" ) ) {
    my $plot_gff = $def->{bamplot_gff};

    # "-g HG19 -y uniform -r"
    my $bamplot_option = getValue( $def, "bamplot_option" );
    my $plotgroups = $def->{plotgroups};
    if ( !defined $plotgroups ) {
      my $files         = $def->{files};
      my @sortedSamples = sort keys %$files;
      $plotgroups = { $task => \@sortedSamples };
    }
    $config->{"plotgroups"} = $plotgroups;
    $config->{"bamplot"}    = {
      class              => "Visualization::Bamplot",
      perform            => 1,
      target_dir         => "${target_dir}/bamplot",
      option             => $bamplot_option,
      source_ref         => "bwa_cleanbam",
      groups_ref         => "plotgroups",
      gff_file           => $plot_gff,
      is_rainbow_color   => 0,
      is_draw_individual => 0,
      is_single_pdf      => 1,
      sh_direct          => 1,
      pbs                => {
        "nodes"    => "1:ppn=1",
        "walltime" => "1",
        "mem"      => "10gb"
      },
    };
    push @$summary, ("bamplot");
  }

  if ( getValue( $def, "perform_report" ) ) {
    addPeakPipelineReport($config, $def, $summary, $target_dir, $task_dic);
  }

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => $individual,
      step_2 => $summary,
    },
    sh_direct => 0,
    pbs       => {
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  #print(Dumper($def->{treatments}));
  #print(Dumper($def));

  return ($config);
}

sub performATACSeq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);
  #print(Dumper($def));

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }
  return $config;
}

sub performATACSeqTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
