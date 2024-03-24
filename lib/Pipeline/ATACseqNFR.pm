#!/usr/bin/perl
package Pipeline::ATACseqNFR;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Pipeline::WdlPipeline;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performATACseqNFR)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "cluster", "slurm" );
  initDefaultValue( $def, "max_thread", 8 );
  initDefaultValue( $def, "sequencetask_run_time", 12 );

  if ( !defined $def->{"treatments"} ) {
    my $files = getValue( $def, "files" );
    my $groups = {};
    for my $sample_name ( sort keys %$files ) {
      $groups->{$sample_name} = [$sample_name];
    }
    $def->{"treatments"} = $groups;
  }

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [ @$individual, @$summary ];

  $def->{replicates} = $config->{groups};

  $def->{adapter} = "";

  my $encode_atac_folder = getValue($def, "encode_atac_folder");
  my $macs2_genome = getValue($def, "macs2_genome");
  
  my $nfr_task = "NFR_filter";
  $config->{$nfr_task} = {
    class                    => "CQS::ProgramWrapperOneToOne",
    perform                  => 1,
    target_dir               => $target_dir . "/" . $nfr_task,
    program => "",
    check_program => 0,
    option                   => "
echo samtools view -h -b -e 'tlen < 150 && tlen > -150' -o __NAME__.NFR_filtered.bam __FILE__ 
samtools view -h -b -e 'tlen < 150 && tlen > -150' -o __NAME__.NFR_filtered.bam __FILE__ 
samtools index __NAME__.NFR_filtered.bam

echo python3 $encode_atac_folder/src/encode_task_bam2ta.py __NAME__.NFR_filtered.bam --paired-end --mito-chr-name chrM --subsample 0 --mem-gb 6.008817602694035
python3 $encode_atac_folder/src/encode_task_bam2ta.py __NAME__.NFR_filtered.bam --paired-end --mito-chr-name chrM --subsample 0 --mem-gb 6.008817602694035

echo macs2 callpeak -t __NAME__.NFR_filtered.tn5.tagAlign.gz -f BED -n __NAME__.NFR_filtered -g $macs2_genome -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits
macs2 callpeak -t __NAME__.NFR_filtered.tn5.tagAlign.gz -f BED -n __NAME__.NFR_filtered -g $macs2_genome -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits

echo makeTagDirectory __NAME__.NFR_filtered.tagdir __NAME__.NFR_filtered.bam -format sam
makeTagDirectory __NAME__.NFR_filtered.tagdir __NAME__.NFR_filtered.bam -format sam

#__OUTPUT__
",
    source_ref               => [ "files" ],
    output_file_ext          => ".NFR_filtered.bam,.NFR_filtered.tn5.tagAlign.gz",
    output_to_result_dir     => 1,
    output_to_same_folder => 0,
    sh_direct                => 0,
    pbs                      => {
      "nodes"     => "1:ppn=8",
      "walltime" => "10",
      "mem"      => "10G"
    }
  };
  push (@$tasks, $nfr_task);

  addSequenceTask($config, $def, $tasks, $target_dir, $summary);

  return($config);
}

sub performATACseqNFR {
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

1;
