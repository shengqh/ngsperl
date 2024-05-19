#!/usr/bin/perl
package Pipeline::SLAMseq;

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

our %EXPORT_TAGS = ( 'all' => [qw(performSLAMseq)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing",   1 );

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

  my $target_dir      = $def->{target_dir};

  my $reference_fasta = getValue($def, "fasta_file");
  my $utr3_bed = getValue($def, "utr3_bed");
  my $max_thread = getValue($def, "max_thread", 8);

  my $nextgenmap_task = "T01_nextgenmap";
  $config->{$nextgenmap_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$nextgenmap_task",
    init_command          => "",
    option                => "
echo ngm=`date`
ngm -r $reference_fasta \\
  -1 __FILE__ \\
  -o __NAME__.sam \\
  --rg-id 1 \\
  --rg-sm __NAME__ \\
  --rg-lb __NAME__ \\
  --rg-pu __NAME__ \\
  --rg-pl ILLUMINA \\
  -t $max_thread
  
status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.ngm.failed
  rm -f __NAME__.sam __NAME__.ngm.succeed
else
  touch __NAME__.ngm.succeed
  rm -f __NAME__.ngm.failed

  echo sort=`date`
  samtools sort -@ $max_thread -o __NAME__.sorted.tmp.bam __NAME__.sam
  
  status=\$?
  if [ \$status -ne 0 ]; then
    touch __NAME__.sort.failed
    rm -f __NAME__.sorted.tmp.bam __NAME__.sort.succeed
  else
    touch __NAME__.sort.succeed
    rm -f __NAME__.sam __NAME__.sort.failed
    mv __NAME__.sorted.tmp.bam __NAME__.sorted.bam
  fi
fi
",
    docker_prefix         => "nextgenmap_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => "files",
    source_arg            => "",
    source_join_delimiter => " -2 ",
    output_to_same_folder => 0,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => ".sorted.bam",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $nextgenmap_task);

  my $slamdunk_filter_task = "T02_slamdunk_filter";
  $config->{$slamdunk_filter_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$slamdunk_filter_task",
    init_command          => "",
    option                => "
slamdunk filter -o . \\
  --bed $utr3_bed \\
  --threads $max_thread \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.filter.failed
  rm -f __NAME__.bam __NAME__.filter.succeed
else
  touch __NAME__.filter.succeed
  rm -f __NAME__.filter.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $nextgenmap_task,
    source_arg            => "",
    source_join_delimiter => "",
    no_output => 1,
    output_to_same_folder => 0,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => ".bam",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $slamdunk_filter_task);

  my $slamdunk_snp_task = "T03_slamdunk_snp";
  $config->{$slamdunk_snp_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$slamdunk_snp_task",
    init_command          => "",
    option                => "
slamdunk snp -o . \\
  -t $max_thread \\
  -r $reference_fasta \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.snp.failed
  rm -f __NAME__.snp.succeed
else
  touch __NAME__.snp.succeed
  rm -f __NAME__.snp.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_filter_task,
    source_arg            => "",
    source_join_delimiter => "",
    no_output => 1,
    output_to_same_folder => 0,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => ".bam",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $slamdunk_snp_task);

  my $slamdunk_count_task = "T04_slamdunk_count";
  $config->{$slamdunk_count_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$slamdunk_count_task",
    init_command          => "",
    option                => "
slamdunk count -o . \\
  -t $max_thread \\
  -r $reference_fasta \\
  __FILE__

status=\$?
if [ \$status -ne 0 ]; then
  touch __NAME__.count.failed
  rm -f __NAME__.count.succeed
else
  touch __NAME__.count.succeed
  rm -f __NAME__.count.failed
fi
",
    docker_prefix         => "slamdunk_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $slamdunk_snp_task,
    source_arg            => "",
    source_join_delimiter => "",
    no_output => 1,
    output_to_same_folder => 0,
    output_arg            => "",
    output_to_folder      => 0,
    output_file_prefix    => "",
    output_file_ext       => ".txt",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=$max_thread",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push (@$tasks, $slamdunk_count_task);

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

sub performSLAMseq {
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
