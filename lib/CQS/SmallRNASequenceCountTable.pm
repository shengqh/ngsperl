#!/usr/bin/perl
package CQS::SmallRNASequenceCountTable;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;
use Data::Dumper;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_srst";
  bless $self, $class;
  return $self;
}

sub get_sample_map {
  my ( $groups, $rawfiles, $group_names ) = @_;
  my $result = {};
  for my $group_name (@$group_names) {
    my @sample_names = @{ $groups->{$group_name} };
    for my $sample_name (@sample_names) {
      die "Unknown sample $sample_name in group $group_name" if !$rawfiles->{$sample_name};
      $result->{$sample_name} = $rawfiles->{$sample_name};
    }
  }
  return $result;
}

sub get_name_files_map {
  my ( $config, $section, $task_name ) = @_;

  my $raw_files = get_raw_files( $config, $section );

  my $name_files_map = {};

  if ( has_raw_files( $config, $section, "pairs" ) ) {
    my $comparisons      = get_raw_files( $config, $section, "pairs" );
    my @comparison_names = sort keys %{$comparisons};
    my $groupsName       = defined $config->{$section}{deseq2_groups} ? "deseq2_groups" : "groups";
    my $groups           = get_raw_files( $config, $section, $groupsName );

    for my $comparison_name (@comparison_names) {
      my $gNames = $comparisons->{$comparison_name};
      my @group_names;

      if ( ref $gNames eq ref {} ) {
        @group_names = @{ $gNames->{groups} };
      }
      else {
        @group_names = @{$gNames};
      }

      $name_files_map->{$comparison_name} = get_sample_map( $groups, $raw_files, \@group_names );
    }
  }
  
  if ( has_raw_files( $config, $section, "groups" ) ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $group_name ( sort keys %{$groups} ) {
      my @group_names = ($group_name);
      $name_files_map->{$group_name} = get_sample_map( $groups, $raw_files, \@group_names );
    }
  }

  if ( scalar( keys %$name_files_map ) > 1 ) {
    $name_files_map->{$task_name} = $raw_files;
  }
  else {
    my $is_total = 0;
    for my $name ( sort keys %{$name_files_map} ) {
      my $samples = $name_files_map->{$name};

      #print Dumper($samples);
      #print Dumper($raw_files);
      if ( scalar( keys %$samples ) == scalar( keys %$raw_files ) ) {
        $is_total = 1;
        last;
      }
    }
    if ($is_total) {
      $name_files_map = { $task_name => $raw_files };
    }
    else {
      $name_files_map->{$task_name} = $raw_files;
    }
  }

  #print Dumper($name_files_map);

  return $name_files_map;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  $option = $option . " --exportFasta";

  my $fastq_files;
  if ( has_raw_files( $config, $section, "fastq_files" ) ) {
    $fastq_files = get_raw_files( $config, $section, "fastq_files" );
  }

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $name_files_map = get_name_files_map( $config, $section, $task_name );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);
  my $pbs      = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );
  for my $name ( sort keys %{$name_files_map} ) {
    my $filelist   = $self->get_file( $pbs_dir,    $name, ".filelist", 0 );
    my $outputfile = $self->get_file( $result_dir, $name, ".count",    0 );
    my $outputname = basename($outputfile);
    my $finalfile = basename( $self->get_file( $result_dir, $name, ".minicontig.count", 0 ) );

    open( my $fl, ">$filelist" ) or die "Cannot create $filelist";
    my $samples = $name_files_map->{$name};
    for my $sample_name ( sort keys %$samples ) {
      my @count_files = @{ $samples->{$sample_name} };
      my $countFile   = $count_files[0];
      print $fl $sample_name, "\t", $countFile;
      if ( defined $fastq_files ) {
        print $fl "\t", $fastq_files->{$sample_name}->[0];
      }
      print $fl "\n";
    }
    close($fl);

    print $pbs "
if [ ! -s $finalfile ]; then
  cqstools smallrna_sequence_count_table $option -o $outputname -l $filelist
fi

";
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern, $canReturnEmpty ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  my $name_files_map = get_name_files_map( $config, $section, $task_name );

  for my $name ( sort keys %$name_files_map ) {
    my @result_files   = ();
    my $outputfile     = $self->get_file( $result_dir, $name, ".count", 0 );
    my $miniContigFile = $self->get_file( $result_dir, $name, ".minicontig.count", 0 );
    my $filelist       = $self->get_file( $pbs_dir, $name, ".filelist", 0 );
    my $readoutputfile = $self->get_file( $result_dir, $name, ".read.count", 0 );
    push( @result_files, $outputfile );
    push( @result_files, $outputfile . ".fasta" );
    push( @result_files, $outputfile . ".details" );
    push( @result_files, $filelist );
    push( @result_files, $miniContigFile );
    push( @result_files, $miniContigFile . ".fasta" );
    push( @result_files, $miniContigFile . ".details" );
    push( @result_files, $readoutputfile );
    push( @result_files, $readoutputfile . ".fasta" );
    $result->{$name} = filter_array( \@result_files, $pattern, $canReturnEmpty );
  }

  return $result;
}

1;
