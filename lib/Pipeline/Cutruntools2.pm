#!/usr/bin/perl
package Pipeline::Cutruntools2;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performCutruntools2)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing",   0 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [@$individual, @$summary];

  my $target_dir = $def->{target_dir};

  my $config_json = getValue($def, "cutruntools2-bulk-config", dirname(__FILE__) . "/../Chipseq/cutruntools2-bulk-config.json");
  my $adaptor_type = getValue($def, "cutruntools2_adaptor_type");

  my $cutruntools2 = "cutruntools2";
  $config->{$cutruntools2} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_dir}/${cutruntools2}",
    option                => "
IN='__FILE__'
arrIN=(\${IN//,/ })
r1=\${arrIN[0]}
r2=\${arrIN[1]}

if [[ ! -s \$r1 ]]; then
  echo file not exists: \$r1
  exit 1
fi

if [[ ! -s \$r2 ]]; then
  echo file not exists: \$r2
  exit 1
fi

rm -f __NAME___R1_001.fastq.gz __NAME___R2_001.fastq.gz

ln -s \$r1 __NAME___R1_001.fastq.gz
ln -s \$r2 __NAME___R2_001.fastq.gz

cur_dir=`pwd`
#echo cur_dir = \$cur_dir

cat $config_json | sed \"s#INPUT_FASTQ_DIRECTORY#\${cur_dir}#g\" | sed \"s#INPUT_WORKDIR#\${cur_dir}#g\" | sed \"s#Truseq#${adaptor_type}#g\" > __NAME__.config.json

bash /opt/CUT-RUNTools-2.0/run_bulkModule.sh __NAME__.config.json __NAME__

rm -f __NAME___R1_001.fastq.gz __NAME___R2_001.fastq.gz

#__FILE__ __OUTPUT__
",
    interpretor           => "",
    program               => "",
    check_program         => 0,
    source_arg            => "-i",
    source_ref            => $untrimed_ref,
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_ext    => "__NAME__/peakcalling/seacr/__NAME__.spikein_normalized.bw,__NAME__/aligned/dup.marked.120bp/__NAME__.bam",
    use_tmp_folder        => 0,
    sh_direct             => 0,
    docker_prefix => "cutruntools2_",
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "40gb"
    },
  };
  push @$tasks, ($cutruntools2);

  if ( $def->{perform_bamplot} ) {
    add_bamplot($config, $def, $tasks, $target_dir, [$cutruntools2, ".bam"]);
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      tasks => $tasks,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
    },
  };

  return ($config);
}

sub performCutruntools2 {
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

1;
