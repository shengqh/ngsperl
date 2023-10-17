#!/usr/bin/perl
package Pipeline::Preprocession;

use strict;
use warnings;

use Text::CSV qw( csv );
use Data::Dumper;
use Hash::Merge qw( merge );

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(
  getPreprocessionConfig
  addCutadapt
  addFastqLen
  addExtractSingleEndFastqFromPairend
  add_fastq_join
  )
  ] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "perform_preprocessing",     1 );
  initDefaultValue( $def, "generate_md5",     0 );
  initDefaultValue( $def, "perform_check_fastq_duplicate", 1);
  initDefaultValue( $def, "cluster",                   'slurm' );
  initDefaultValue( $def, "sra_to_fastq",              0 );
  initDefaultValue( $def, "check_file_exists",         1 );
  initDefaultValue( $def, "merge_fastq",               0 );
  initDefaultValue( $def, "perform_dedup_fastq",       0 );
  initDefaultValue( $def, "fastq_remove_N",            0 );
  initDefaultValue( $def, "perform_fastqc",            1 );
  initDefaultValue( $def, "perform_cutadapt_test",     0 );
  initDefaultValue( $def, "perform_cutadapt",          0 );
  initDefaultValue( $def, "remove_sequences",          "" );
  initDefaultValue( $def, "table_vis_group_text_size", '10' );
  initDefaultValue( $def, "max_thread",                '8' );
  initDefaultValue( $def, "sequencetask_run_time",     '12' );
  initDefaultValue( $def, "emailType",                 "FAIL" );
  initDefaultValue( $def, "extractSingleEndFastqFromPairend", 0 );
  initDefaultValue( $def, "sra_to_fastq_sh_direct", 0 );
  
  # initDefaultValue( $def, "rmdformats_theme", "rmdformats::readthedown");
  # initDefaultValue( $def, "toc_depth", "3");

  return $def;
}

sub addCutadapt {
  my ($config, $def, $individual, $summary, $cutadapt_task, $fastqcName, $intermediate_dir, $preprocessing_dir, $source_ref, $is_pairend, $cluster) = @_;
  my $default_thread = $is_pairend?2:1;
  my $cutadapt_thread = getValue($def, "cutadapt_thread", $default_thread);
  if($def->{extractSingleEndFastqFromPairend}){
    $def->{cutadapt_option} = $def->{cutadapt_option} . " -O " . getValue($def, "cutadapt_overlap", 8);
  }
  print("cutadapt_option=" . $def->{cutadapt_option} . "\n");
  my $cutadapt_class = ( defined $def->{cutadapt_config} ) ? "Trimmer::CutadaptByConfig" : "Trimmer::Cutadapt";
  my $cutadapt = {
    "$cutadapt_task" => {
      class                            => $cutadapt_class,
      perform                          => 1,
      target_dir                       => $intermediate_dir . "/" . getNextFolderIndex($def) . "$cutadapt_task",
      option                           => $def->{cutadapt_option},
      use_option_only                  => $def->{use_cutadapt_option_only},
      source_ref                       => $source_ref,
      config                           => $def->{cutadapt_config},
      adapter                          => $def->{adapter},
      adapter_5                        => $def->{adapter_5},
      adapter_3                        => $def->{adapter_3},
      random_bases_remove_after_trim   => $def->{"fastq_remove_random"},
      random_bases_remove_after_trim_5 => $def->{"fastq_remove_random_5"},
      random_bases_remove_after_trim_3 => $def->{"fastq_remove_random_3"},
      fastq_remove_random              => $def->{"fastq_remove_random"},
      fastq_remove_random_5            => $def->{"fastq_remove_random_5"},
      fastq_remove_random_3            => $def->{"fastq_remove_random_3"},
      trim_polyA                   => $def->{trim_polyA},
      trim_base_quality_after_adapter_trim => $def->{trim_base_quality_after_adapter_trim},
      hard_trim                        => $def->{"hard_trim"},
      extension                        => "_clipped.fastq",
      is_paired_end                    => $is_pairend,
      sh_direct                        => 0,
      cluster                          => $cluster,
      pbs                              => {
#        "email"     => $def->{email},
#        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=" . $cutadapt_thread,
        "walltime"  => "14",
        "mem"       => "20gb"
      },
    }
  };
  if ( defined $def->{cutadapt} ) {
    $cutadapt->{cutadapt} = merge_hash_right_precedent( $def->{cutadapt}, $cutadapt->{cutadapt} );
  }

  for my $key (keys %$cutadapt){
    $config->{$key} = $cutadapt->{$key};
  }
  push @$individual, ($cutadapt_task);

  if ($is_pairend) {
    if (getValue($def, "perform_paired_end_validation", 1)){
      my $fastq_validator = $cutadapt_task . "_validation";
      addPairendFastqValidation($config, $def, $individual, $intermediate_dir, $fastq_validator, [$cutadapt_task, ".fastq.gz"]);
    }
  }

  if ( $def->{perform_fastqc} ) {
    addFastQC( $config, $def, $individual, $summary, $fastqcName, [ $cutadapt_task, ".fastq.gz" ], $preprocessing_dir );
  }
  $source_ref = [ $cutadapt_task, ".fastq.gz" ];
  return ($source_ref, $cutadapt_task);
}

sub addFastqLen {
  my ($config, $def, $individual, $summary, $fastqLenName, $parentFolder, $source_ref, $cluster) = @_;

  my $fastq_len_dir = $parentFolder . "/" . getNextFolderIndex($def) . $fastqLenName;
  my $fastq_len     = {
    "$fastqLenName" => {
      class => "CQS::ProgramWrapperOneToOne",
      target_dir => $fastq_len_dir,
      #option => "-c",
      option => "",
      use_tmp_folder => 1,
      suffix  => "_flen",
      interpretor => "python3",
      program => "../QC/fastq_len.py",
      source_arg => "-i",
      source_ref => $source_ref,
      output_arg => "-o",
      output_file_prefix => ".len",
      output_file_ext => ".len",
      output_to_same_folder => 1,
      can_result_be_empty_file => 0,
      sh_direct   => 0,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "40gb"
      }
    },
    "${fastqLenName}_vis" => {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $parentFolder . "/" . getNextFolderIndex($def) . "${fastqLenName}_summary",
      rtemplate                => "countTableVisFunctions.R,fastqLengthVis.R",
      output_file              => ".lengthDistribution",
      output_file_ext          => ".png",
      output_other_ext         => ".csv",
      parameterSampleFile1_ref => [ $fastqLenName, ".len\$" ],
      parameterSampleFile2     => $def->{groups},
      sh_direct                => 1,
      pbs                      => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    }
  };

  for my $key (keys %$fastq_len){
    $config->{$key} = $fastq_len->{$key};
  }

  push @$individual, ($fastqLenName);
  push @$summary,    ("${fastqLenName}_vis");

  return($config);
}

sub initCutadaptOption {
  my ( $config, $def, $fastq_remove_N ) = @_;

  if ( $def != $config ) {
    initDefaultValue( $config, "min_read_length", getValue($def, "min_read_length", 0) );
  }
  my $cutadapt_option = getValue( $config, "cutadapt_option", getValue( $def, "cutadapt_option", "" ) );
  if(! getValue($def, "use_cutadapt_option_only", 0)){
    if ( ( $cutadapt_option !~ /-a/ ) && ( $cutadapt_option !~ /-g/ ) ) {
      defined $config->{"adapter_5"} or defined $config->{"adapter_3"} or getValue( $config, "adapter" );
    }

    if ( $cutadapt_option !~ /\-m/ ) {
      my $min_read_length = getValue( $config, "min_read_length" );
      $cutadapt_option = $cutadapt_option . " -m " . $min_read_length;
    }

    if ( $cutadapt_option !~ /\-n/ ) {
      my $max_adapter_count = getValue( $config, "max_adapter_count", 1 );
      $cutadapt_option = $cutadapt_option . " -n " . $max_adapter_count;
    }

    $config->{cutadapt_option} = $cutadapt_option;
  }
}

sub getPreprocessionConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = create_directory_or_die( getValue( $def, "target_dir" ) );
  my $task    = getValue( $def, "task_name" );

  if ($def->{sra_to_fastq}) {
    if (not defined $def->{files}) {
      if (defined $def->{SraRunTable_file}){
        addFilesFromSraRunTable($def, $def->{SraRunTable_file});
      }
    }
  }

  if (defined $def->{groups_pattern}) {
    $def->{groups} = get_groups_by_pattern($def);
  }

  if(defined $def->{covariance_file}){
    if(!defined $def->{groups} && defined $def->{sample_column} && defined $def->{group_column}){
      $def->{groups} = get_groups_from_covariance_file($def->{covariance_file}, $def->{sample_column}, $def->{group_column});
    }
  }

  if (defined $def->{files}) {
    my $all_sample_names = get_all_sample_names($def);

    if(not defined $def->{groups}){
      my $files = $def->{files};
      my $sampleNames = [keys %$files];
      $def->{groups} = {"All" => $sampleNames};
    }

    if(!defined $def->{HTO_samples}){
      checkFileGroupPairNames($def, ["groups"], ["pairs"], $all_sample_names, getValue($def, "remove_missing_samples_in_group", 0));
    }
  }

  if(defined $def->{ignore_samples}){
    my $ignore_samples = $def->{ignore_samples};
    
    if(! is_array($ignore_samples)){
      die("ignore_samples should be array.")
    }

    my %ignore_sample_map = map { $_ => 1 } @$ignore_samples;
    
    # my $files = $def->{files};
    # for my $ignore_sample (@$ignore_samples){
    #   delete $files->{$ignore_sample};
    # }
    # $def->{files} = $files;

    my $groups = $def->{groups};
    if(defined $groups){
      my $group_names = [keys %$groups];
      for my $group_name (@$group_names){
        my $filter_samples = [];
        my $cur_samples = $groups->{$group_name};
        for my $sample (@$cur_samples){
          if (not defined $ignore_sample_map{$sample}) {
            push @$filter_samples, $sample;
          }
        }
        $groups->{$group_name} = $filter_samples;
      }
      $def->{groups} = $groups;
    }
  }

  if(defined $def->{groups}){
    my $groups = $def->{groups};
    if(defined $def->{unique_group_names}){
      my $unique_group_names = $def->{unique_group_names};
      my $unique_groups = {};
      for my $ugname (@$unique_group_names){
        my $ugroup = $groups->{$ugname};
        if(defined $ugroup){
          $unique_groups->{$ugname} = $ugroup;
        }
      }
      $def->{unique_groups} = $unique_groups;
    }else{
      $def->{unique_groups} = get_unique_groups($groups);
    }

    if(not defined $def->{correlation_groups}){
      if(defined $def->{correlation_groups_dic}){
        my $correlationGroups = get_pair_group_sample_map( $def->{correlation_groups_dic}, $def->{groups} );
        if ( getValue( $def, "correlation_all", 1 ) and (not defined $correlationGroups->{all}) ) {
          $correlationGroups->{all} = $def->{unique_groups};
        }
        $def->{correlation_groups} = $correlationGroups;
      }else{
        $def->{correlation_groups} = get_correlation_groups_by_pattern($def);
      }
    }
  }
  if (defined $def->{covariant_patterns}){
    $def->{covariance_patterns} = $def->{covariant_patterns};
  }

  if (defined $def->{covariance_patterns}){
    $def->{covariance_file} = create_covariance_file_by_pattern($def);
  }

  $def = initializeDefaultOptions($def);

  my $preprocessing_dir = $target_dir;
  if ( $def->{perform_preprocessing} && $def->{subdir} ) {
    $preprocessing_dir = $target_dir . "/preprocessing";
    if(!$def->{donot_create_directory_or_die}){
      create_directory_or_die($preprocessing_dir);
    }
  }

  my $intermediate_dir = getIntermidiateDir($preprocessing_dir, $def);

  my $is_pairend = is_paired_end($def);
  if (defined $is_pairend){
    print ("is_pairend=" . $is_pairend . "\n");
  }

  #general
  my $cluster = getValue( $def, "cluster" );
  my $email   = getValue( $def, "email" );

  if ((! $def->{sra_to_fastq}) && (defined $def->{files}) && $def->{check_file_exists}){
    #all file defined in $files should be hard-coding with absolute path, check file
    my $sourcefiles   = getValue( $def, "files" );
    foreach my $filename (keys %$sourcefiles){
      my $fileref = $sourcefiles->{$filename};
      foreach my $eachfile (@$fileref) {
        if (! -e $eachfile){
          die "file not exists (" . $filename . "): " . $eachfile;
        }
      }
    }
  }

  #data
  my $config = {
    general => {
      task_name  => $task,
      cluster    => $cluster,
      email      => $email,
      emailType  => getValue( $def, "emailType", "FAIL" ),
      constraint => $def->{constraint},
      account    => $def->{account},
      debug      => $def->{debug},
      sratoolkit_setting_file => $def->{sratoolkit_setting_file},
      interval_list_file => $def->{interval_list_file},
      R          => getValue($def, "R", "R"),
      Rscript    => getValue($def, "Rscript", "Rscript"),
      R_LIBS     => $def->{"R_LIBS"},
      localize_to_local_folder => getValue($def, "localize_to_local_folder", 0),
      use_tmp_folder => getValue($def, "use_tmp_folder", 0),
    },
    ignore_samples => $def->{ignore_samples},
    files                => $def->{files},
    raw_files => $def->{raw_files},
    groups               => $def->{groups},
    deseq2_groups        => $def->{deseq2_groups},
    pairs                => $def->{pairs},
    additional_bam_files => $def->{additional_bam_files},
    singularity_image_files => $def->{singularity_image_files},
  };
  
  if(not defined $def->{ignore_docker} or not $def->{ignore_docker}){
    foreach my $key (keys %$def){
      if ($key =~ /docker_/){
        $config->{general}{$key} = $def->{$key};
      }
    }
  }
  
  my $source_ref = ["files"];
  my $individual = [];
  my $summary    = [];

  if ( !$def->{perform_preprocessing} ) {
    return ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $source_ref, $cluster, $intermediate_dir );
  }

  #task
  if ( $def->{sra_to_fastq} ) {
    defined $is_pairend or die "Define is_paired_end first!";
    defined $def->{is_restricted_data} or die "Define is_restricted_data first!";
    #defined $def->{sra_table} or die "Define sra_table first, can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab";
  }

  if($def->{generate_md5}){
    $config = add_md5($config, $def, $individual, $preprocessing_dir, $source_ref);
  }

  if (defined $def->{files}) {
    if(getValue($def, "perform_md5_validation", 0)){
      die "Define md5_files first if you want to validate file m5." if not defined $def->{md5_files};

      if(!defined $config->{md5_merge}){
        $config = add_md5($config, $def, $individual, $preprocessing_dir, $source_ref);
      }

      $config->{md5_validation} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $preprocessing_dir . "/" . getNextFolderIndex($def) . "md5_validation",
        rtemplate                => "countTableVisFunctions.R,../QC/validateMD5.r",
        output_file              => ".csv",
        output_file_ext          => ".csv",
        output_other_ext         => ".png",
        parameterFile1_ref => "md5_merge",
        parameterSampleFile1_ref => $source_ref,
        parameterSampleFile2     => {
          "md5files" => $def->{md5_files},
        },
        sh_direct                => 1,
        pbs                      => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      }      
    }
  }

  if ( $def->{merge_fastq} ) {
    defined $is_pairend or die "Define is_paired_end first!";
  }

  my $fastq_remove_N   = getValue( $def, "fastq_remove_N" );
  my $remove_sequences = getValue( $def, "remove_sequences" );    #remove contamination sequences from sequence kit before adapter trimming

  if ( (defined $remove_sequences) and ($remove_sequences ne "") ) {
    defined $is_pairend or die "Define is_paired_end first!";
  }

  my $run_cutadapt     = getValue( $def, "perform_cutadapt" );
  my $run_cutadapt_test  = getValue( $def, "perform_cutadapt_test" );
  
  if ($run_cutadapt or $run_cutadapt_test) {
    if ( defined $def->{cutadapt} ) {
      my $cconfig = $def->{cutadapt};
      initCutadaptOption( $cconfig, $def, $fastq_remove_N );
    }
    elsif ( defined $def->{cutadapt_config} ) {
      my $cconfig = $def->{cutadapt_config};
      for my $key ( keys %$cconfig ) {
        my $kconfig = $cconfig->{$key};
        initCutadaptOption( $kconfig, $def, $fastq_remove_N );
      }
    }
    else {
      initCutadaptOption( $def, $def, $fastq_remove_N );
    }
    $fastq_remove_N = 0;
  }

  if ( $def->{sra_to_fastq} ) {
    my $class = getValue($def, "sra_to_fastq_wget", 0)? "SRA::Wget" :"SRA::FastqDump";
    #my $class = getValue($def, "sra_to_fastq_wget", 0)? "SRA::Wget" :"SRA::FasterqDump";
    $config->{sra2fastq} = {
      class      => $class,
      perform    => 1,
      is_paired_end   => $is_pairend,
      target_dir => $intermediate_dir . "/" . getNextFolderIndex($def) . "sra2fastq",
      option     => "",
      source_ref => $source_ref,
      sra_table  => $def->{sra_table},
      sh_direct  => getValue($def, "sra_to_fastq_sh_direct", 0),
      cluster    => $def->{cluster},
      not_clean  => getValue( $def, "sra_not_clean", 1 ),
      is_restricted_data => getValue($def, "is_restricted_data"),
      docker_prefix => "sratools_",
      no_docker => $def->{"no_docker"},
      pbs        => {
        "nodes"     => "1:ppn=1",
        "walltime"  => getValue($def, "sra_to_fastq_walltime", "24"),
        "mem"       => "10gb"
      },
    };
    $source_ref = "sra2fastq";
    push @$individual, ("sra2fastq");
  }

  if (getValue($def, "perform_paired_end_validation", 1)){
    defined $is_pairend or die "Define is_paired_end first!";
    if($is_pairend){
      my $fastq_validator = "paired_end_validation";
      addPairendFastqValidation($config, $def, $individual, $intermediate_dir, $fastq_validator, $source_ref);
    }
  }

  if ( $def->{merge_fastq} ) {
    if ( $def->{perform_fastqc} ) {
      if(getValue($def, "perform_fastqc_before_merge", 1)){
        addFastQC( $config, $def, $individual, $summary, "fastqc_raw_before_merge", $source_ref, $preprocessing_dir );
      }
    }

    my $class = getValue($def, "is_fastq_gzipped", 0) ? "Format::MergeFastqGzipped" : "Format::MergeFastq";
    $config->{merge_fastq} = {
      class       => $class,
      perform     => 1,
      target_dir  => $intermediate_dir . "/" . getNextFolderIndex($def) . "merge_fastq",
      option      => "",
      source_ref  => $source_ref,
      sh_direct   => 0,
      is_paired_end   => $is_pairend,
      is_bzipped  => $def->{is_bzipped},
      is_collated => $def->{is_collated},
      cluster     => $def->{cluster},
      pbs         => {
        "nodes"     => "1:ppn=1",
        "walltime"  => getValue($def, "merge_fastq_time", 4),
        "mem"       => "10gb"
      }
    };
    $source_ref = "merge_fastq";
    push @$individual, ("merge_fastq");
  }

  if ($def->{perform_qc_check_fastq_duplicate}){
    $config->{qc_check_fastq_duplicate} = {
      class => "CQS::ProgramWrapperOneToOne",
      target_dir => $intermediate_dir . "/" . getNextFolderIndex($def) . "qc_check_fastq_duplicate",
      interpretor => "python3",
      program => "../QC/checkFastqDuplicate.py",
      source_arg => "-i",
      source_ref => $source_ref,
      output_arg => "-o",
      output_file_prefix => "",
      output_file_ext => ".txt",
      output_to_same_folder => 1,
      can_result_be_empty_file => 1,
      sh_direct   => 1,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "40gb"
      }
    };
    push @$individual, ("qc_check_fastq_duplicate");
  }

  if ($def->{perform_dedup_fastq}){
    $config->{dedup_fastq} = {
      class => "CQS::ProgramWrapperOneToOne",
      target_dir => $intermediate_dir . "/" . getNextFolderIndex($def) . "dedup_fastq",
      interpretor => "python3",
      program => "../Format/dedupFastq.py",
      source_arg => "-i",
      source_ref => $source_ref,
      output_arg => "-o",
      output_file_prefix => ".dedup",
      output_file_ext => ".dedup.1.fastq.gz",
      output_other_ext => ".dedup.2.fastq.gz",
      output_to_same_folder => 1,
      sh_direct   => 0,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      }
    };
    $source_ref = "dedup_fastq";
    push @$individual, ("dedup_fastq");
  }

  if ($fastq_remove_N) {
    $config->{fastq_remove_N} = {
      class      => "CQS::FastqTrimmer",
      perform    => 1,
      target_dir => $intermediate_dir . "/" . getNextFolderIndex($def) . "fastq_remove_N",
      option     => "",
      extension  => "_trim.fastq.gz",
      source_ref => $source_ref,
      sh_direct  => 1,
      cluster    => $def->{cluster},
      pbs        => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "2",
        "mem"       => "10gb"
      }
    };
    $source_ref = "fastq_remove_N";
    push @$individual, ("fastq_remove_N");
  }

  if ( ( $source_ref ne "files" ) and ( defined $def->{fastqs} ) ) {
    $config->{fastqs} = $def->{fastqs};
    $source_ref = [ $source_ref, "fastq.gz\$", "fastqs" ];
  }

  if ( $def->{perform_fastqc} ) {
    addFastQC( $config, $def, $individual, $summary, "fastqc_raw", $source_ref, $preprocessing_dir );
  }

  if ( $remove_sequences ne "" ) {
    if ($is_pairend){
      $config->{"remove_contamination_sequences"} = {
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => $intermediate_dir . "/" . getNextFolderIndex($def) . "remove_contamination_sequences",
        option => "-s " . $remove_sequences,
        interpretor => "python3",
        program => "../Format/removeSequence.py",
        source_arg => "-i",
        source_ref => $source_ref,
        source_join => ",",
        output_arg => "-o",
        output_file_prefix => "_removeSeq",
        output_file_ext => "_removeSeq.1.fastq.gz",
        output_other_ext => "_removeSeq.2.fastq.gz",
        output_to_same_folder => 1,
        sh_direct   => 0,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
        }
      };
    }else{
      $config->{"remove_contamination_sequences"} = {
        class      => "CQS::Perl",
        perform    => 1,
        target_dir => $intermediate_dir . "/" . getNextFolderIndex($def) . "remove_contamination_sequences",
        option     => $remove_sequences,
        output_ext => "_removeSeq.fastq.gz",
        perlFile   => "removeSequenceInFastq.pl",
        source_ref => $source_ref,
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "2",
          "mem"       => "20gb"
        },
      };
    }
    push @$individual, ("remove_contamination_sequences");
    $source_ref = [ "remove_contamination_sequences", ".fastq.gz" ];

    if ( $def->{perform_fastqc} ) {
      addFastQC( $config, $def, $individual, $summary, "fastqc_post_remove", $source_ref, $preprocessing_dir );
    }
  }

  my $untrimed_ref = $source_ref;

  if ($run_cutadapt_test){
    my $tasks = [];
    my $test_dir = create_directory_or_die( $target_dir . "/test_cutadapt" );
  
    my $extractTask = "test_extract";
    my $output_file_ext = $is_pairend? ".1.fastq.gz":".fastq.gz";
    my $output_other_ext = $is_pairend? ".2.fastq.gz":undef;
    my $extract_option = $is_pairend? "-p":"";
    $config->{$extractTask} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               =>  $test_dir . "/" . $extractTask,
      option                   => $extract_option,
      interpretor              => "python3",
      program                  => "../Pipeline/extractFirstNFastq.py",
      parameterSampleFile1_arg => "-i",
      parameterSampleFile1_ref => $source_ref,
      output_arg               => "-o",
      output_file              => "",
      output_file_ext          => $output_file_ext,
      output_other_ext         => $output_other_ext,
      sh_direct                => 1,
      pbs                      => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$tasks, $extractTask );

    my ($test_ref, $cutadapt_task) = addCutadapt($config, $def, $tasks, $tasks, "test_cutadapt", "test_fastqc_post_trim", $test_dir, $test_dir, $extractTask, $is_pairend, $cluster);
    addFastqLen($config, $def, $tasks, $tasks, "test_fastq_len", $test_dir, [$cutadapt_task, ".gz"], $cluster );

    print(join(',', @$tasks));
    $config->{"test_sequencetask"} = {
      class      => getSequenceTaskClassname($cluster),
      perform    => 1,
      target_dir => "${test_dir}/test_sequencetask",
      option     => "",
      source     => {
        step_1 => $tasks,
      },
      sh_direct => 0,
      pbs       => {
        "nodes"    => "1:ppn=8",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    performConfig($config);
    delete $config->{test_sequencetask};
  }

  if ($run_cutadapt) {
    my $cutadapt_task;
    ($source_ref, $cutadapt_task)  = addCutadapt($config, $def, $individual, $summary, "cutadapt", "fastqc_post_trim", $intermediate_dir, $preprocessing_dir, $source_ref, $is_pairend, $cluster);
    addFastqLen($config, $def, $individual, $summary, "fastq_len", $preprocessing_dir, [$cutadapt_task, ".gz"], $cluster );
  }elsif ($def->{fastq_len} ) {
    addFastqLen($config, $def, $individual, $summary, "fastq_len", $preprocessing_dir, $source_ref, $cluster );
  }

  if($def->{extractSingleEndFastqFromPairend}){
    my $extract_task = "extract_singleend_fastq";
    my $fastqc_task = "fastqc_extract_singleend";
    addExtractSingleEndFastqFromPairend($config, $def, $individual, $summary, $extract_task, $fastqc_task, $intermediate_dir, $preprocessing_dir, $source_ref, $cluster );
    $source_ref = [$extract_task, ".se.fastq.gz"];
    addFastqLen($config, $def, $individual, $summary, "fastq_len_extract", $preprocessing_dir, $source_ref, $cluster );
  }elsif($def->{perform_fastq_join}){
    my $extract_task = "fastq_join";
    my $fastqc_task = "fastqc_fastq_join";
    add_fastq_join($config, $def, $individual, $summary, $extract_task, $fastqc_task, $intermediate_dir, $preprocessing_dir, $source_ref, $cluster );
    $source_ref = [$extract_task, ".join"];
    addFastqLen($config, $def, $individual, $summary, "fastq_len_extract", $preprocessing_dir, $source_ref, $cluster );
  }

  if ( $def->{perform_fastqc} ) {
    my $fastqc_count_vis_files = undef;
    if ( length($remove_sequences) && $run_cutadapt ) {
      $fastqc_count_vis_files = {
        target_dir         => $config->{fastqc_post_trim_summary}->{target_dir},
        parameterFile2_ref => [ "fastqc_post_remove_summary", ".FastQC.reads.tsv\$" ],
        parameterFile3_ref => [ "fastqc_post_trim_summary", ".FastQC.reads.tsv\$" ],
      };
    }
    elsif ( length($remove_sequences) ) {
      $fastqc_count_vis_files = {
        target_dir         => $config->{fastqc_post_remove_summary}->{target_dir},
        parameterFile2_ref => [ "fastqc_post_remove_summary", ".FastQC.reads.tsv\$" ],
      };
    }
    elsif ($run_cutadapt) {
      $fastqc_count_vis_files = {
        target_dir         => $config->{fastqc_post_trim_summary}->{target_dir},
        parameterFile2_ref => [ "fastqc_post_trim_summary", ".FastQC.reads.tsv\$" ],
      };
    }

    if ( defined $fastqc_count_vis_files ) {
      $config->{"fastqc_count_vis"} = merge_hash_right_precedent(
        {
          class              => "CQS::UniqueR",
          perform            => 1,
          suffix => "_fqcv",
          rtemplate          => "countTableVisFunctions.R,countInFastQcVis.R",
          output_file        => ".countInFastQcVis.Result",
          output_file_ext    => ".Reads.csv",
          output_other_ext   => ".pdf,.png",
          sh_direct          => 1,
          parameterFile1_ref => [ "fastqc_raw_summary", ".FastQC.reads.tsv\$" ],
          can_result_be_empty_file => 1,
          pbs                => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        },
        $fastqc_count_vis_files
      );
      push @$summary, ("fastqc_count_vis");
    }
  }

  $config->{def} = $def;
  
  return ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster, $run_cutadapt_test );
}

sub addExtractSingleEndFastqFromPairend {
  my ($config, $def, $individual, $summary, $extract_task, $fastqc_task, $intermediate_dir, $preprocessing_dir, $source_ref) = @_;
  my $minReadLength = getValue($def, "minReadLength", 16);
  my $minSimilarityRatio = getValue($def, "minSimilarityRatio", 0.8);
  $config->{$extract_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$intermediate_dir/$extract_task",
    #init_command          => "source /data/cqs/softwares/cqsperl/scripts/path_conda.txt",
    option                => "-o __NAME__.se.fastq.gz --minReadLength $minReadLength --minSimilarityRatio 0.8",
    interpretor           => "python3",
    check_program         => 1,
    program               => "../SmallRNA/getFastqFromTrimmedPairendReads.py",
    source_ref            => $source_ref,
    source_arg            => "--read1",
    source_join_delimiter => " --read2 ",
    no_output             => 0,
    output_to_same_folder => 1, 
    output_arg            => "-o",
    output_ext            => ".se.fastq.gz",
    sh_direct             => 0,
    use_tmp_folder        => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push @$individual, ($extract_task);

  get_result_file($config, $extract_task, "");

  if ( $def->{perform_fastqc} ) {
    addFastQC( $config, $def, $individual, $summary, $fastqc_task, [ $extract_task, ".se.fastq.gz" ], $preprocessing_dir );
  }
}

sub add_fastq_join {
  my ($config, $def, $individual, $summary, $fastq_join_task, $fastqc_task, $intermediate_dir, $preprocessing_dir, $source_ref) = @_;
  $config->{$fastq_join_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$intermediate_dir/$fastq_join_task",
    option                => "fastq-join __FILE__ -o __NAME__ ",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $source_ref,
    source_arg            => "",
    source_join_delimiter => " ",
    no_output             => 0,
    output_to_same_folder => 1, 
    output_arg            => "-o",
    output_ext            => ".join",
    sh_direct             => 0,
    use_tmp_folder        => 1,
    no_docker => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push @$individual, ($fastq_join_task);

  get_result_file($config, $fastq_join_task, "");

  if ( $def->{perform_fastqc} ) {
    addFastQC( $config, $def, $individual, $summary, $fastqc_task, $fastq_join_task, $preprocessing_dir );
  }
}

1;
