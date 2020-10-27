#!/usr/bin/perl
package Pipeline::Preprocession;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(
  getPreprocessionConfig
  addCutadapt
  addFastqLen)
  ] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "perform_preprocessing",     1 );
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

  return $def;
}

sub addCutadapt {
  my ($config, $def, $individual, $summary, $cutadaptName, $fastqcName, $intermediate_dir, $preprocessing_dir, $source_ref, $is_pairend, $cluster) = @_;
  my $default_thread = $is_pairend?2:1;
  my $cutadapt_thread = getValue($def, "cutadapt_thread", $default_thread);
  print("cutadapt_thread=" . $cutadapt_thread . "\n");
  my $cutadapt_class = ( defined $def->{cutadapt_config} ) ? "Trimmer::CutadaptByConfig" : "Trimmer::Cutadapt";
  my $cutadapt = {
    "$cutadaptName" => {
      class                            => $cutadapt_class,
      perform                          => 1,
      target_dir                       => $intermediate_dir . "/" . getNextFolderIndex($def) . "$cutadaptName",
      option                           => $def->{cutadapt_option},
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
      trim_poly_atgc                   => $def->{trim_poly_atgc},
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
  push @$individual, ($cutadaptName);

  if ( $def->{perform_fastqc} ) {
    addFastQC( $config, $def, $individual, $summary, $fastqcName, [ $cutadaptName, ".fastq.gz" ], $preprocessing_dir );
  }
  $source_ref = [ $cutadaptName, ".fastq.gz" ];
  return ($source_ref, $cutadaptName);
}

sub addFastqLen {
  my ($config, $def, $individual, $summary, $fastqLenName, $parentFolder, $source_ref, $cluster) = @_;

  my $fastq_len_dir = $parentFolder . "/" . getNextFolderIndex($def) . $fastqLenName;
  my $fastq_len     = {
    "$fastqLenName" => {
      class      => "CQS::FastqLen",
      perform    => 1,
      target_dir => $fastq_len_dir,
      option     => "",
      source_ref => $source_ref,
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "23",
        "mem"       => "20gb"
      },
    },
    "${fastqLenName}_vis" => {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $fastq_len_dir,
      rtemplate                => "countTableVisFunctions.R,fastqLengthVis.R",
      output_file              => ".lengthDistribution",
      output_file_ext          => ".png;.csv",
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
    initDefaultValue( $config, "min_read_length", $def->{"min_read_length"} );
  }
  my $cutadapt_option = getValue( $config, "cutadapt_option", getValue( $def, "cutadapt_option", "" ) );

  if ( ( $cutadapt_option !~ /-a/ ) && ( $cutadapt_option !~ /-g/ ) ) {
    defined $config->{"adapter_5"} or defined $config->{"adapter_3"} or getValue( $config, "adapter" );
  }

  if ( $cutadapt_option !~ /\-m/ ) {
    my $min_read_length = getValue( $config, "min_read_length" );
    $cutadapt_option = $cutadapt_option . " -m " . $min_read_length;
  }

  if ( $cutadapt_option !~ /\-n/ ) {
    my $max_adapter_count = getValue( $config, "max_adapter_count", 2 );
    $cutadapt_option = $cutadapt_option . " -n " . $max_adapter_count;
  }

  $config->{cutadapt_option} = $cutadapt_option;
}

sub getPreprocessionConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = create_directory_or_die( getValue( $def, "target_dir" ) );
  my $task    = getValue( $def, "task_name" );

  if (defined $def->{groups_pattern}) {
    my $gpattern = $def->{groups_pattern};
    my $files = $def->{files};
    my $groups = {};
    for my $samplename (sort keys %$files) {
      $samplename =~ /$gpattern/;
      my $groupname = $1;
      #print($groupname . " : " . $samplename . "\n");
      if (not defined $groups->{$groupname}){
        $groups->{$groupname} = [$samplename];
      }
      else{
        my $samples = $groups->{$groupname};
        push (@$samples, $samplename);
      }
    }

    #print("groups=" . Dumper($groups));
    $def->{groups} = $groups;
  }

  if (defined $def->{covariance_patterns}){
    my $files = $def->{files};
    my $covariance_patterns = $def->{covariance_patterns};
    my @covariances = (sort keys %$covariance_patterns);
    my $cov_map = {};
    for my $covariance (@covariances) {
      my $cov_pattern_def = $covariance_patterns->{$covariance};

      my $cov_pattern;
      my $cov_prefix;
      if (ref $cov_pattern_def eq 'HASH'){
        $cov_pattern = $cov_pattern_def->{pattern};
        $cov_prefix = $cov_pattern_def->{prefix};
      }
      else{
        $cov_pattern = $cov_pattern_def;
        $cov_prefix = "";
      }

      $cov_map->{$covariance} = {};
      for my $samplename (sort keys %$files) {
        $samplename =~ /$cov_pattern/;
        my $cov_value = $1;
        $cov_map->{$covariance}{$samplename} = $cov_prefix . $cov_value;
      }
    }

    #print("covariances=" . Dumper($cov_map));

    my $cov_file = $target_dir . "/covariance.txt";
    open( my $cov, ">$cov_file" ) or die "Cannot create $cov_file";
    print $cov "Sample";
    for my $covariance (@covariances) {
      print $cov "\t" . $covariance;
    }
    print $cov "\n";
    for my $samplename (sort keys %$files) {
      print $cov "$samplename";
      for my $covariance (@covariances) {
        print $cov "\t" . $cov_map->{$covariance}{$samplename};
      }
      print $cov "\n";
    }
    close($cov);

    $def->{covariance_file} = $cov_file;
  }

  if(defined $def->{ignore_samples}){
    my $ignore_samples = $def->{ignore_samples};

    my %ignore_sample_map = map { $_ => 1 } @$ignore_samples;
    
    my $files = $def->{files};
    for my $ignore_sample (@$ignore_samples){
      delete $files->{$ignore_sample};
    }
    $def->{files} = $files;

    my $groups = $def->{groups};
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

  if (defined $def->{pool_sample}){
    checkFileGroupPairNames($def, ["pool_sample_groups"], ["pairs"], "files");
    checkFileGroupPairNames($def, ["groups"], ["pairs"], "pool_sample_groups");
  }else{
    checkFileGroupPairNames($def, ["groups"], ["pairs"], "files");

    if(not defined $def->{groups}){
      my $files = $def->{files};
      my $sampleNames = [keys %$files];
      $def->{groups} = {"All" => $sampleNames};
    }
  }

  $def = initializeDefaultOptions($def);

  my $preprocessing_dir = $target_dir;
  if ( $def->{perform_preprocessing} && $def->{subdir} ) {
    $preprocessing_dir = create_directory_or_die( $target_dir . "/preprocessing" );
  }

  my $intermediate_dir = getIntermidiateDir($preprocessing_dir, $def);

  my $is_pairend = is_paired_end($def);
  if (defined $is_pairend){
    print ("is_pairend=" . $is_pairend . "\n");
  }

  #general
  my $cluster = getValue( $def, "cluster" );
  my $email   = getValue( $def, "email" );

  if ((! $def->{sra_to_fastq}) && $def->{check_file_exists}){
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
    },
    files                => $def->{files},
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

  if ( $def->{merge_fastq} ) {
    defined $is_pairend or die "Define is_paired_end first!";
  }

  my $fastq_remove_N   = getValue( $def, "fastq_remove_N" );
  my $remove_sequences = getValue( $def, "remove_sequences" );    #remove contamination sequences from sequence kit before adapter trimming
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
    $config->{sra2fastq} = {
      class      => "SRA::FastqDump",
      perform    => 1,
      is_paired_end   => $is_pairend,
      target_dir => $intermediate_dir . "/" . getNextFolderIndex($def) . "sra2fastq",
      option     => "",
      source_ref => $source_ref,
      sra_table  => $def->{sra_table},
      sh_direct  => 1,
      cluster    => $def->{cluster},
      not_clean  => getValue( $def, "sra_not_clean", 1 ),
      is_restricted_data => getValue($def, "is_restricted_data"),
      pbs        => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    $source_ref = "sra2fastq";
    push @$individual, ("sra2fastq");
  }

  if ( $def->{merge_fastq} ) {
    $config->{merge_fastq} = {
      class       => "Format::MergeFastq",
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
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
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
      interpretor => "python",
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
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
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
      interpretor => "python",
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
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
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
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
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

  if ( length($remove_sequences) ) {
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
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "2",
        "mem"       => "20gb"
      },
    };
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
    my $output_file_ext = $is_pairend? ".1.fastq.gz;.2.fastq.gz":".fastq.gz";
    my $extract_option = $is_pairend? "-p":"";
    $config->{$extractTask} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               =>  $test_dir . "/" . $extractTask,
      option                   => $extract_option,
      interpretor              => "python",
      program                  => "../Pipeline/extractFirstNFastq.py",
      parameterSampleFile1_arg => "-i",
      parameterSampleFile1_ref => $source_ref,
      output_arg               => "-o",
      output_file              => "",
      output_file_ext          => $output_file_ext,
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

    my ($test_ref, $cutadapt_name) = addCutadapt($config, $def, $tasks, $tasks, "test_cutadapt", "test_fastqc_post_trim", $test_dir, $test_dir, $extractTask, $is_pairend, $cluster);
    addFastqLen($config, $def, $tasks, $tasks, "test_fastq_len", $test_dir, [$cutadapt_name, ".gz"], $cluster );

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
    my $cutadapt_name;
    ($source_ref, $cutadapt_name)  = addCutadapt($config, $def, $individual, $summary, "cutadapt", "fastqc_post_trim", $intermediate_dir, $preprocessing_dir, $source_ref, $is_pairend, $cluster);
    addFastqLen($config, $def, $individual, $summary, "fastq_len", $preprocessing_dir, [$cutadapt_name, ".gz"], $cluster );
  }elsif ($def->{fastq_len} ) {
    addFastqLen($config, $def, $individual, $summary, "fastq_len", $preprocessing_dir, $source_ref, $cluster );
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
          rtemplate          => "countInFastQcVis.R",
          output_file        => ".countInFastQcVis.Result",
          output_file_ext    => ".Reads.csv;.pdf",
          sh_direct          => 1,
          parameterFile1_ref => [ "fastqc_raw_summary", ".FastQC.reads.tsv\$" ],
          pbs                => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
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

  return ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster, $run_cutadapt_test );
}

1;
