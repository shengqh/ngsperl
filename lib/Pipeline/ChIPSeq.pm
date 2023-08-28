#!/usr/bin/perl
package Pipeline::ChIPSeq;

use strict;
use warnings;
use CQS::StringUtils;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Pipeline::PeakPipelineUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performChIPSeq performChIPSeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;
  initDefaultValue( $def, "subdir", 0 );

  initDefaultValue( $def, "sra_to_fastq", 0 );
  initDefaultValue( $def, "mark_duplicates", 1 );

  my $is_paired = is_paired_end($def);
  
  if ($is_paired) {
    initDefaultValue( $def, "aligner", "bwa" );
  }
  else {
    initDefaultValue( $def, "aligner", "bowtie" );
  }
  if ( $def->{aligner} eq "bowtie1" ) {
    initDefaultValue( $def, "bowtie1_option", "-v 1 -m 1 --best --strata" );
  }
  elsif ( $def->{aligner} eq "bwa" ) {
    initDefaultValue( $def, "bwa_option", "" );
  }
  elsif ( $def->{aligner} eq "bowtie2" ) {
    initDefaultValue( $def, "bowtie2_option", "" );
  }

  initDefaultValue( $def, "peak_caller", "macs" );
  if ( $def->{"peak_caller"} eq "macs" ) {
    initDefaultValue( $def, "macs_option", "-p 1e-9 -w -S --space=50" );
  }
  elsif ( $def->{peak_caller} eq "macs2" ) {
    initDefaultValue( $def, "macs2_peak_type", "narrow" );
#    initDefaultValue( $def, "macs2_output_bigwig", 0 );
    if ( not defined $def->{"macs2_option"} ) {
      my $macs2_genome    = getValue( $def, "macs2_genome" );      #hs
      my $macs2_peak_type = getValue( $def, "macs2_peak_type" );
      my $pairend         = is_paired_end( $def );
      if($macs2_peak_type eq "broad"){
        if(defined $def->{"macs2_broad_option"}){
          $def->{"macs2_option"} = $def->{"macs2_broad_option"};
        }
      }
      if($macs2_peak_type eq "narrow"){
        if(defined $def->{"macs2_narrow_option"}){
          $def->{"macs2_option"} = $def->{"macs2_narrow_option"};
        }
      }

      if(!defined $def->{"macs2_option"}){
        my $defaultOption   = "-B -q 0.01 -g " . $macs2_genome;
        if ( $macs2_peak_type ne "narrow" ) {
          $defaultOption = "--broad " . $defaultOption;
        }
        if ($pairend) {
          $defaultOption = "-f BAMPE " . $defaultOption;
        }

        $def->{"macs2_option"} = $defaultOption;
      }
    }
  }

  if(defined $def->{treatments}){
    $def->{treatments_auto} = 0;
  }else{
    initDefaultValue( $def, "treatments_auto",     1 );
  }

  if(defined $def->{design_table}){
    $def->{design_table_auto} = 0;
  }else{
    initDefaultValue( $def, "design_table_auto",     0 );
  }

  initDefaultValue( $def, "perform_rose",     0 );
  initDefaultValue( $def, "perform_bamplot",  0 );
  initDefaultValue( $def, "perform_cleanbam", 0 );

  if ( $def->{perform_cleanbam} ) {
    initDefaultValue( $def, "minimum_maq",         10 );
    initDefaultValue( $def, "minimum_insert_size", 30 );
    initDefaultValue( $def, "maximum_insert_size", 1000 );
  }

  initDefaultValue( $def, "perform_chipqc",   1 );
  initDefaultValue( $def, "perform_diffbind", 0 );
  initDefaultValue( $def, "perform_enhancer", 0 );
  initDefaultValue( $def, "perform_multiqc",  1 );

  initDefaultValue( $def, "perform_homer", 1 );
  initDefaultValue( $def, "perform_merge_peaks", 0 );

  initDefaultValue( $def, "perform_report", 1 );
  
  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster, $run_cutadapt_test ) = getPreprocessionConfig($def);

  init_treatments_design_table($def);

  my $perform_chipqc = getValue( $def, "perform_chipqc" );

  my $perform_diffbind = getValue( $def, "perform_diffbind" );
  if ($perform_diffbind) {
    getValue( $def, "design_table" );
  }

  my $step2 = [];

  my $email    = getValue( $def, "email" );

  if (! $run_cutadapt_test){
    if ( $def->{aligner} eq "bowtie1" ) {
      $config->{bowtie1} = {
        class                 => "Alignment::Bowtie1",
        perform               => 1,
        target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "bowtie1",
        option                => getValue( $def, "bowtie1_option" ),
        fasta_file            => getValue( $def, "bowtie1_fasta" ),
        bowtie1_index         => getValue( $def, "bowtie1_index" ),
        source_ref            => $source_ref,
        output_to_same_folder => 1,
        picard_jar            => getValue( $def, "picard_jar" ),
        mark_duplicates       => getValue( $def, "mark_duplicates" ),
        sh_direct             => 0,
        pbs                   => {
          "email"    => $email,
          "nodes"    => "1:ppn=8",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
      $config->{bowtie1_summary} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        rCode                    => "",
        target_dir               => "${target_dir}/" . getNextFolderIndex($def) . "bowtie1_summary",
        option                   => "",
        parameterSampleFile1_ref => [ "bowtie1", ".log" ],
        rtemplate                => "../Alignment/Bowtie1Summary.r",
        output_file              => "",
        output_file_ext          => ".csv",
        pbs                      => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "10gb"
        },
      };

      push @$summary, ("bowtie1_summary");
    }
    elsif ( $def->{aligner} eq "bwa" ) {
      $config->{ $def->{aligner} } = {
        class                 => "Alignment::BWA",
        perform               => 1,
        target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "bwa",
        option                => getValue( $def, "bwa_option" ),
        bwa_index             => getValue( $def, "bwa_fasta" ),
        source_ref            => $source_ref,
        output_to_same_folder => 1,
        picard_jar            => getValue( $def, "picard_jar" ),
        mark_duplicates       => 1,
        sh_direct             => 0,
        pbs                   => {
          "email"    => $email,
          "nodes"    => "1:ppn=8",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
    }
    elsif ( $def->{aligner} eq "bowtie2" ) {
      $config->{ bowtie2 } = {
        class                 => "Alignment::Bowtie2",
        perform               => 1,
        target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "bowtie2",
        option                => getValue( $def, "bowtie2_option" ),
        bowtie2_index         => getValue( $def, "bowtie2_index" ),
        source_ref            => $source_ref,
        output_to_same_folder => 1,
        sh_direct             => 0,
        pbs                   => {
          "email"    => $email,
          "nodes"    => "1:ppn=8",
          "walltime" => getValue($def, "bowtie2_walltime", "72"),
          "mem"      => getValue($def, "bowtie2_mem", "40gb")
        },
      };
      add_alignment_summary($config, $def, $summary, $target_dir, "bowtie2_summary", "../Alignment/AlignmentUtils.r;../Alignment/Bowtie2Summary.r", ".reads.csv;.reads.png;.chromosome.csv;.chromosome.png", [ "bowtie2", ".log" ], ["bowtie2", ".chromosome.count"] );
    }
    else {
      die "Unknown alinger " . $def->{aligner};
    }
    my $bam_ref = [ $def->{aligner}, ".bam\$" ];
    push @$individual, ( $def->{aligner} );

    add_bam_validation($config, $def, $individual, $target_dir, $def->{aligner} . "_bam_validation", $bam_ref );

    if ( $def->{perform_cleanbam} ) {
      my $taskName = $def->{aligner} . "_cleanbam";
      addCleanBAM( $config, $def, $individual, $taskName, "${target_dir}/" . getNextFolderIndex($def) . $taskName, $bam_ref);
      $bam_ref = [ $taskName, ".bam\$" ];

      add_alignment_summary($config, $def, $summary, $target_dir, "${taskName}_summary", "../Alignment/AlignmentUtils.r;../Alignment/Bowtie2Summary.r", ".chromosome.csv;.chromosome.png", undef, [$taskName, ".chromosome.count"] );
      add_bam_validation($config, $def, $individual, $target_dir, "${taskName}_bam_validation", [$taskName, ".bam\$"] );
    }

    if ( getValue( $def, "perform_bamplot" ) ) {
      my $plotgroups = $def->{plotgroups};
      if ( !defined $plotgroups ) {
        my $files         = $def->{files};
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
        is_multi_page => getValue($def, "bamplot_multi_page", 1),
        draw_by_r => getValue($def, "bamplot_draw_by_r", 1),
        draw_by_r_width => getValue($def, "bamplot_draw_by_r_width", 10),
        draw_by_r_height => getValue($def, "bamplot_draw_by_r_height", 10),
        sh_direct          => 1,
        pbs                => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      };
      push @$summary, ("bamplot");
    }

    if($def->{perform_bamsnap} && $def->{"bamsnap_locus"}){
      addBamsnapLocus($config, $def, $summary, $target_dir, "bamsnap_locus", $bam_ref);
    }

    if(defined $def->{annotation_locus} or defined $def->{annotation_genes}){
      my $sizeFactorTask = "size_factor";
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
      push( @$summary, $sizeFactorTask );

      if(defined $def->{annotation_locus}){
        my $locusFile = $target_dir . "/annotation_locus.bed";
        open(my $fh, '>', $locusFile) or die "Could not open file '$locusFile' $!";
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
          print $fh $chr . "\t" . $start . "\t" . $end . "\t1000\t" . $locusName . "\t+\t" . $locusName . "\n";
        }
        close($fh);

        addPlotGene($config, $def, $summary, $target_dir, "annotation_locus_plot", $sizeFactorTask, $locusFile, $bam_ref);
      }

      if(defined $def->{annotation_genes}){
        my $genes_str = $def->{annotation_genes};
        my @genes = split /[;, ]+/, $genes_str;
        my %gene_map = map { $_ => 1 } @genes;
        $config->{annotation_genes} = \%gene_map;
        #print(Dumper($config->{annotation_genes}));

        my $geneLocus = addGeneLocus($config, $def, $summary, $target_dir);

        if($def->{perform_bamsnap}){
          my $bamsnap_task = "annotation_genes_bamsnap";
          my $annotation_genes_shift = $def->{"annotation_genes_shift"};
          my $bamsnap_option = defined $annotation_genes_shift ? "-e " . $annotation_genes_shift : "";
          my $params = {
            bamsnap_option => $bamsnap_option,
            bamsnap_raw_option => $def->{bamsnap_raw_option}
          };

          addBamsnap($config, $def, $summary, $target_dir, $bamsnap_task, [$geneLocus, "bed"], $bam_ref, $params);
        }

        addPlotGene($config, $def, $summary, $target_dir, "annotation_genes_plot", $sizeFactorTask, [ $geneLocus, ".bed" ], $bam_ref);
      }
    }

    my $peakCallerTask;
    my $callFilePattern;

    if ( $def->{peak_caller} eq "macs" ) {
      $peakCallerTask            = "macs1callpeak";
      $callFilePattern           = ".name.bed\$";
      $config->{$peakCallerTask} = {
        class      => "Chipseq::MACS",
        perform    => 1,
        target_dir => "${target_dir}/" . getNextFolderIndex($def) . "${peakCallerTask}",
        option     => getValue( $def, "macs_option" ),
        source_ref => $bam_ref,
        groups     => $def->{"treatments"},
        controls   => $def->{"controls"},
        sh_direct  => 0,
        pbs        => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
      $callFilePattern = ".name.bed\$";
    }
    elsif ( $def->{peak_caller} eq "macs2" ) {
      $peakCallerTask = "macs2callpeak";
      my $macs2option = getValue( $def, "macs2_option" );
      if ( $macs2option =~ /broad/ ) {
        $callFilePattern = "broadPeak.bed\$";
      }
      else {
        $callFilePattern = "narrowPeak.bed\$";
      }
      my $macs2_output_bigwig = getValue( $def, "macs2_output_bigwig" ,0);
      #my $chr_size_file = "/scratch/cqs_share/references/gencode/GRCh38.p13/bowtie2_index_2.4.3/GRCh38.primary_assembly.genome.len";
      my $chr_size_file = getValue( $def, "chr_size_file" ,"");
      $peakCallerTask = $peakCallerTask . "_" . $def->{macs2_peak_type};
      $config->{$peakCallerTask} = {
        class      => "Chipseq::MACS2Callpeak",
        perform    => 1,
        target_dir => "${target_dir}/" . getNextFolderIndex($def) . "$peakCallerTask",
        option     => $macs2option,
        source_ref => $bam_ref,
        groups     => $def->{"treatments"},
        controls   => $def->{"controls"},
        output_bigwig=>$macs2_output_bigwig,
        chr_size_file   => $chr_size_file,
        can_result_be_empty_file => 1,
        sh_direct  => 0,
        pbs        => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
    }
    else {
      die "Unknown peak caller " . $def->{"peak_caller"};
    }
    push (@$step2, $peakCallerTask);

    my $peak_count_task = add_peak_count($config, $def, $summary, $target_dir, $peakCallerTask . "_count", $peakCallerTask);

    # if(getValue($def, "perform_annotateNearestGene", 0)){
    #   annotateNearestGene($config, $def, $summary, $target_dir,  [$peakCallerTask, ".bed"]);
    # }

    if(getValue($def, "perform_annovar", 1)){
      my $annovar_task = addAnnovar( $config, $def, $summary, $target_dir, $peakCallerTask, ".bed", undef, undef, undef, 1, 0, 1 );

      if(0){
        my $annovar_vis_task = $annovar_task . "_vis";
        $config->{$annovar_vis_task} = {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          rCode                    => "",
          target_dir               => "${target_dir}/" . getNextFolderIndex($def) . $annovar_vis_task,
          option                   => "",
          parameterSampleFile1_ref => $annovar_task,
          rtemplate                => "../Visualization/peakAnnovarVis.r",
          output_file              => "",
          output_file_ext          => ".png",
          pbs                      => {
            "nodes"    => "1:ppn=1",
            "walltime" => "1",
            "mem"      => "10gb"
          },
        };
        push(@$summary, $annovar_vis_task);
      }
    }

    if ( getValue( $def, "perform_merge_peaks" )) {
      my $mergePeakTask = $peakCallerTask . "_mergePeaks";
      $config->{$mergePeakTask} = {
        class      => "Homer::MergePeaks",
        perform    => 1,
        target_dir => "${target_dir}/" . getNextFolderIndex($def) . "$mergePeakTask",
        option     => "",
        source_ref => [ $peakCallerTask, $callFilePattern ],
        groups     => getValue($def, "merge_peaks_groups"),
        sh_direct  => 1,
        pbs        => {
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "5gb"
        },
      };
      push @$summary, $mergePeakTask;
    }

    my $homer_name;
    if ( getValue( $def, "perform_homer" ) ) {
      $homer_name = addHomerAnnotation( $config, $def, $summary, $target_dir, $peakCallerTask, $callFilePattern );
    }

    if(getValue($def, "perform_crc", 0)){
      $def->{perform_activeGene} = 1;
      $def->{perform_rose} = 1;
    }

    my $active_gene_task;
    if (getValue($def, "perform_activeGene", 0)) {
      $active_gene_task = $peakCallerTask . "_active_gene";
      $config->{$active_gene_task} = {
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => "${target_dir}/" . getNextFolderIndex($def) . $active_gene_task,
        interpretor => "python3",
        program => "../Chipseq/activeGene.py",
        option => "-g " . getValue($def, "active_gene_genome"),
        source_arg => "-i",
        source_ref => [$peakCallerTask, ".bed"],
        output_arg => "-o",
        output_file_prefix => "",
        output_file_ext => ".TSS_ACTIVE_-1000_1000.txt",
        output_to_same_folder => 1,
        can_result_be_empty_file => 0,
        docker_prefix => "crc_",
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push (@$step2, $active_gene_task);
    }

    my $roseTask;
    if ( $def->{perform_rose} ) {
      $config->{bam_treatments} = {
        class => "CQS::GroupPickTask",
        source_ref => $bam_ref,
        groups => getValue($def, "treatments"),
      };

      my $output_file_ext;
      if ($def->{peak_caller} eq "macs"){
        $output_file_ext = "__NAME__/__NAME___peaks_SuperEnhancers_ENHANCER_TO_GENE.txt",
      }else{
        die "contact quanhu.sheng.1\@vumc.org to fix the problem.";
        #"__NAME__/overlap_SuperEnhancers_ENHANCER_TO_GENE.txt"
      }

      $roseTask = $peakCallerTask . "_rose2";
      $config->{$roseTask} = {
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => "${target_dir}/$roseTask",
        interpretor => "",
        program => "ROSE2",
        check_program => 0,
        docker_prefix => "crc_",
        option => "-o __NAME__/ -g " . getValue($def, "rose_genome"),
        source_arg => "-i",
        source_ref => [ $peakCallerTask, ".bed\$" ],
        parameterSampleFile2_arg => "-r",
        parameterSampleFile2_ref => "bam_treatments",
        output_arg => "-o",
        output_file_prefix => "",
        output_file_ext => $output_file_ext,
        output_to_same_folder => 1,
        can_result_be_empty_file => 0,
        sh_direct => 0,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      },

      my $has_control = defined $def->{controls};
      if($has_control){
        $config->{bam_controls} = {
          class => "CQS::GroupPickTask",
          source_ref => $bam_ref,
          groups => getValue($def, "controls"),
        };
        $config->{$roseTask}{parameterSampleFile3_arg} = "-c";
        $config->{$roseTask}{parameterSampleFile3_ref} = "bam_controls";
      }

      push @$summary, ($roseTask);
    }

    if(getValue($def, "perform_crc", 0)){
      my $crc_task = $peakCallerTask . "_crc";
      $config->{$crc_task} = { 
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => "${target_dir}/$crc_task",
        interpretor => "",
        program => "crc",
        check_program => 0,
        docker_prefix => "crc_",
        option => "-n __NAME__ -c " . getValue($def, "enhancer_genome_path") . " -g " . getValue($def, "rose_genome"),
        source_arg => "-s",
        source_ref => [ $peakCallerTask, ".bed\$" ],
        parameterSampleFile2_arg => "-e",
        parameterSampleFile2_ref => [$roseTask, "SuperEnhancers_ENHANCER_TO_GENE.txt"],
        parameterSampleFile3_arg => "-a",
        parameterSampleFile3_ref => [$active_gene_task, ".txt"],
        output_arg => "-o",
        output_file_prefix => "",
        output_file_ext => "__NAME__/peak/overlap_reproducibility/overlap.optimal_peak.narrowPeak.gz.bed,__NAME__/qc/qc.html,__NAME__/align/rep1/__NAME__.nodup.bam,__NAME__/align/ctl1/input.nodup.bam",
        output_to_same_folder => 1,
        can_result_be_empty_file => 0,
        sh_direct => 0,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary, ($crc_task);
    }

    my $chipqc_taskname = $peakCallerTask . "_chipqc";
    if ($perform_chipqc) {
      add_chipqc($config, $def, $summary, $target_dir, $chipqc_taskname, [ $peakCallerTask, ".bed\$" ], $bam_ref);
    }

    my $bindName = $peakCallerTask . "_diffbind";
    my $bind_homer_name;
    if ($perform_diffbind) {
      addDiffbind($config, $def, $summary, $target_dir, $bindName, $bam_ref, [ $peakCallerTask, ".bed\$" ]);

      if ( getValue( $def, "perform_homer" ) ) {
        $bind_homer_name = addHomerAnnotation( $config, $def, $summary, $target_dir, $bindName, ".sig.bed" );
      }
    }

    if (getValue($def, "perform_bdgdiff", 0)) {
    #if (1) {
      my $bindName = $peakCallerTask . "_bdgdiff";
      $config->{$bindName} = {
        class                   => "Chipseq::MACS2Bdgdiff",
        perform                 => 1,
        target_dir              => "${target_dir}/" . getNextFolderIndex($def) . "${bindName}",
        option                  => "",
        source_ref              => $peakCallerTask,
        groups                  => $def->{"treatments_bdgdiff"},
        sh_direct               => 0,
        pbs                     => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
    }

    if ( getValue( $def, "perform_enhancer" ) ) {
      addEnhancer( $config, $def, $individual, $summary, $target_dir, $peakCallerTask . "_enhancer", $bam_ref, [ $peakCallerTask, ".bed\$" ] );
    }

    if ( getValue( $def, "perform_multiqc" ) ) {
      addMultiQC( $config, $def, $summary, $target_dir, $target_dir );
    }

    if ( getValue( $def, "perform_report" ) ) {
      my $task_dic = {
        peak_caller => $peakCallerTask,
        peak_count => $peak_count_task,
      };
      
      if ( $perform_chipqc ) {
        $task_dic->{chipqc} = $chipqc_taskname;
      }

      if( $def->{perform_homer}) {
        $task_dic->{homer} = $homer_name;
      }

      if($perform_diffbind) {
        $task_dic->{diffbind} = $bindName;
        if( $def->{"perform_homer"}) {
          $task_dic->{diffbind_homer} = $bind_homer_name;
        }
      }

      addPeakPipelineReport($config, $def, $summary, $target_dir, $task_dic);
    }
  }

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => $individual,
      step_2 => $step2,
      step_3 => $summary,
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return ($config);
}

sub performChIPSeq {
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

sub performChIPSeqTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
