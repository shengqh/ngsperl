#!/usr/bin/perl
package Pipeline::Qiime2;

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

our %EXPORT_TAGS = ( 'all' => [qw(initializeQiime2DefaultOptions performQiime2 performQiime2)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeQiime2DefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "ALL" );
  initDefaultValue( $def, "cluster",   "slurm" );

  initDefaultValue( $def, "perform_preprocessing", 1 );
  initDefaultValue( $def, "perform_mapping",       1 );
  initDefaultValue( $def, "perform_counting",      1 );
  initDefaultValue( $def, "perform_count_table",   1 );
  initDefaultValue( $def, "perform_report",        1 );

  initDefaultValue( $def, "perform_cutadapt", 0 );

  initDefaultValue( $def, "max_thread", 8 );
  initDefaultValue( $def, "sequencetask_run_time", '24' );

  initDefaultValue( $def, "is_paired_end", 1 );

  initDefaultValue( $def, "DE_outputPdf",  getValue( $def, "outputPdf",  0 ) );
  initDefaultValue( $def, "DE_outputPng",  getValue( $def, "outputPng",  1 ) );
  initDefaultValue( $def, "DE_outputTIFF", getValue( $def, "outputTIFF", 0 ) );

  return $def;
}

sub getQiime2Config {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeQiime2DefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my $thread = getValue($def, "max_thread", 8);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $target_dir = $def->{target_dir};

  my $cache_dir = create_directory_or_die($target_dir . "/qiime2_cache");

  my $files = $def->{files};
  my $groups = $def->{groups};
  my $file_group_map = {};
  my $group_header = "";
  if(defined $groups){
    $group_header = "\tGroup";
    for my $gname (keys %$groups){
      my $gfiles = $groups->{$gname};
      for my $gfile (@$gfiles){
        $file_group_map->{$gfile} = $gname;
      }
    }
  }

  my $link_dir = create_directory_or_die($target_dir . "/linked_data");
  my $manifest_file = $target_dir . "/MANIFEST";
  my $meta_file = $target_dir . "/metafile.txt";
  open(my $mani_fout, '>', $manifest_file) or die "Cannot create file $manifest_file";
  open(my $meta_fout, '>', $meta_file) or die "Cannot create file $meta_file";
  $mani_fout->print("sample-id,filename,direction\n");
  $meta_fout->print("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tgroup\n");
  $meta_fout->print("#q2:types\tcategorical\tcategorical\tcategorical\n");
  
  for my $name (sort keys %$files){
    my $cfiles = $files->{$name};
    my $cfile1 = $cfiles->[0];
    my $cname = basename($cfile1);
    $cname =~ s/_S\d+_L\d+_R\d+.+//g;
    $meta_fout->print("$cname\t\t\t" . ($group_header eq "" ? "" : $file_group_map->{$name}) . "\n");

    my $idx = 0;
    for my $cfile (@$cfiles){
      my $tfile = $link_dir. "/" . basename($cfile);
      if ( ! -e $tfile ) {
        `ln -s $cfile $tfile`;
      }
      $mani_fout->print("$cname,$cfile," . ($idx == 0?"forward":"reverse") . "\n");
      $idx = $idx + 1;
    }
  }
  close($mani_fout);
  close($meta_fout);

  $def->{folder_index} = 1;

  my $qiime2_artifact = "qiime2_artifact";
  $config->{$qiime2_artifact} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => $target_dir . "/" . getNextFolderIndex($def) .  $qiime2_artifact,
    option => "
export CONDA_PREFIX=$cache_dir
export MPLCONFIGDIR=$cache_dir

if [[ -e __NAME__.failed ]]; then
  rm __NAME__.failed
fi

qiime tools import \\
  --type 'SampleData[PairedEndSequencesWithQuality]' \\
  --input-path __FILE__ \\
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \\
  --output-path __NAME__.qza | tee __NAME__.import.log

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.succeed

  qiime demux summarize \\
    --i-data __NAME__.qza \\
    --o-visualization __NAME__.qzv  
else
  touch __NAME__.failed
  rm __NAME__.qza
fi

#__OUTPUT__
",
    use_tmp_folder => 0,
    suffix  => "_qa",
    docker_prefix => "qiime2_",
    interpretor => "",
    program => "",
    check_program => 0,
    source_arg => "--input-path",
    source => {
      $task_name => [$link_dir],
    },
    output_arg => "--output-path",
    output_file_prefix => ".qza",
    output_file_ext => ".qza",
    output_to_same_folder => 1,
    sh_direct   => getValue($def, "sh_direct", 0),
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    }
  };
  push (@$individual, $qiime2_artifact);

  my $source = $qiime2_artifact;
  if(getValue($def, "perform_qiime2_cutadapt", 0)){
    my $qiime2_cutadapt = "qiime2_cutadapt";
    my $qiime2_cutadapt_front_f = getValue($def, "qiime2_cutadapt_front_f", "GTGCCAGCMGCCGCGGTAA");
    my $qiime2_cutadapt_front_r = getValue($def, "qiime2_cutadapt_front_r", "GACTACNVGGGTATCTAATCC");
    $config->{$qiime2_cutadapt} = {
      class => "CQS::ProgramWrapperOneToOne",
      target_dir => $target_dir . "/" . getNextFolderIndex($def) .$qiime2_cutadapt,
      option => "
export CONDA_PREFIX=$cache_dir
export MPLCONFIGDIR=$cache_dir

if [[ -e __NAME__.failed ]]; then
  rm __NAME__.failed
fi

echo qiime cutadapt trim-paired

qiime cutadapt trim-paired \\
  --i-demultiplexed-sequences __FILE__ \\
  --p-cores $thread \\
  --p-front-f $qiime2_cutadapt_front_f \\
  --p-front-r $qiime2_cutadapt_front_r \\
  --o-trimmed-sequences __NAME__.trimmed.qza \\
  --verbose | tee __NAME__.trim.log


status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.succeed

  qiime demux summarize \\
    --i-data __NAME__.trimmed.qza \\
    --o-visualization __NAME__.trimmed.qzv
else
  touch __NAME__.failed
  rm __NAME__.trimmed.qza
fi

#__OUTPUT__
  ",
      use_tmp_folder => 0,
      suffix  => "_qc",
      docker_prefix => "qiime2_",
      interpretor => "",
      program => "",
      check_program => 0,
      source_arg => "–i-demultiplexed-sequences",
      source_ref => $qiime2_artifact,
      output_arg => "–o-trimmed-sequences",
      output_file_prefix => ".trimmed.qza",
      output_file_ext => ".trimmed.qza",
      output_to_same_folder => 1,
      sh_direct   => getValue($def, "sh_direct", 0),
      pbs => {
        "nodes"     => "1:ppn=$thread",
        "walltime"  => "20",
        "mem"       => "40gb"
      }
    };
    push (@$individual, $qiime2_cutadapt);
    $source = $qiime2_cutadapt;
  }

  my $qiime2_dada2 = "qiime2_dada2";
  my $trim_left_f = getValue($def, "p-trim-left-f", 13);
  my $trim_left_r = getValue($def, "p-trim-left-r", 13);
  my $trunc_len_f = getValue($def, "p-trunc-len-f", 230);
  my $trunc_len_r = getValue($def, "p-trunc-len-r", 160);
  
  $config->{$qiime2_dada2} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => $target_dir . "/" . getNextFolderIndex($def) . $qiime2_dada2,
    option => "
export CONDA_PREFIX=$cache_dir
export MPLCONFIGDIR=$cache_dir

if [[ -e __NAME__.failed ]]; then
  rm __NAME__.failed
fi

echo qiime dada2 denoise-paired

qiime dada2 denoise-paired \\
  --i-demultiplexed-seqs __FILE__ \\
  --p-trim-left-f $trim_left_f \\
  --p-trim-left-r $trim_left_r \\
  --p-trunc-len-f $trunc_len_f \\
  --p-trunc-len-r $trunc_len_r \\
  --o-representative-sequences __NAME__.dada2.qza \\
  --o-table __NAME__.dada2.table.qza \\
  --o-denoising-stats __NAME__.dada2.stats.qza \\
  --p-n-threads $thread \\
  --verbose | tee __NAME__.dada2.denoising.log

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.succeed
  qiime feature-table summarize \\
    --i-table __NAME__.dada2.table.qza \\
    --o-visualization __NAME__.dada2.table.qzv

  qiime feature-table tabulate-seqs \\
    --i-data __NAME__.dada2.qza \\
    --o-visualization __NAME__.dada2.qzv  

  qiime metadata tabulate \\
    --m-input-file __NAME__.dada2.stats.qza \\
    --o-visualization __NAME__.dada2.stats.qzv    
else
  touch __NAME__.failed
  rm __NAME__.dada2.qza
fi

#__OUTPUT__
",
    use_tmp_folder => 0,
    suffix  => "_qd",
    docker_prefix => "qiime2_",
    interpretor => "",
    program => "",
    check_program => 0,
    source_arg => "--i-demultiplexed-seqs",
    source_ref => $source,
    output_arg => "--o-representative-sequences",
    output_file_prefix => ".dada2.qza",
    output_file_ext => ".dada2.qza,.dada2.table.qza",
    output_to_same_folder => 1,
    sh_direct   => getValue($def, "sh_direct", 0),
    pbs => {
      "nodes"     => "1:ppn=$thread" ,
      "walltime"  => "20",
      "mem"       => "40gb"
    }
  };
  push (@$individual, $qiime2_dada2);

  my $qiime2_phylogenetic = "qiime2_phylogenetic";
  $config->{$qiime2_phylogenetic} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => $target_dir . "/" . getNextFolderIndex($def) . $qiime2_phylogenetic,
    option => "
export CONDA_PREFIX=$cache_dir
export MPLCONFIGDIR=$cache_dir

if [[ -e __NAME__.failed ]]; then
  rm __NAME__.failed
fi

qiime phylogeny align-to-tree-mafft-fasttree \\
  --i-sequences __FILE__ \\
  --o-alignment __NAME__.aligned.qza \\
  --o-masked-alignment __NAME__.aligned.masked.qza \\
  --o-tree __NAME__.unrooted-tree.qza \\
  --o-rooted-tree __NAME__.rooted-tree.qza \\
  --verbose | tee __NAME__.rooted-tree.log

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.succeed

else
  touch __NAME__.failed
  rm __NAME__.aligned.qza __NAME__.aligned.masked.qza __NAME__.unrooted-tree.qza __NAME__.rooted-tree.qza
fi

#__OUTPUT__
",
    suffix  => "_qp",
    docker_prefix => "qiime2_",
    interpretor => "",
    program => "",
    check_program => 0,
    source_arg => "--i-sequences",
    source_ref => [ $qiime2_dada2, ".dada2.qza" ],
    output_arg => "--o-alignment",
    output_file_prefix => ".rooted-tree.qza",
    output_file_ext => ".rooted-tree.qza",
    output_to_same_folder => 1,
    sh_direct   => getValue($def, "sh_direct", 0),
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    }
  };
  push (@$individual, $qiime2_phylogenetic);

  my $qiime2_diversity = "qiime2_diversity";
  $config->{$qiime2_diversity} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => $target_dir . "/" . getNextFolderIndex($def) . $qiime2_diversity,
    option => "
export CONDA_PREFIX=$cache_dir
export MPLCONFIGDIR=$cache_dir

if [[ -e __NAME__.failed ]]; then
  rm __NAME__.failed
fi

qiime diversity core-metrics-phylogenetic \\
  --i-phylogeny __FILE__ \\
  --i-table __FILE2__ \\
  --m-metadata-file $meta_file \\
  --p-sampling-depth 1103 \\
  --output-dir __NAME__ \\
  --verbose | tee __NAME__.phylogenetic.log

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.succeed

  qiime metadata tabulate \\
    --m-input-file __NAME__/faith_pd_vector.qza \\
    --o-visualization __NAME__/faith_pd_vector.qzv

  qiime diversity alpha-group-significance \\
    --i-alpha-diversity __NAME__/faith_pd_vector.qza \\
    --m-metadata-file $meta_file \\
    --o-visualization __NAME__/faith-pd-group-significance.qzv

  qiime metadata tabulate \\
    --m-input-file __NAME__/shannon_vector.qza \\
    --o-visualization __NAME__/shannon_vector.qzv

  qiime diversity alpha-group-significance \\
    --i-alpha-diversity __NAME__/shannon_vector.qza \\
    --m-metadata-file $meta_file \\
    --o-visualization __NAME__/shannon-group-significance.qzv

  qiime metadata tabulate \\
    --m-input-file __NAME__/evenness_vector.qza \\
    --o-visualization __NAME__/evenness_vector.qzv

  qiime diversity alpha-group-significance \\
    --i-alpha-diversity __NAME__/evenness_vector.qza \\
    --m-metadata-file $meta_file \\
    --o-visualization __NAME__/evenness-group-significance.qzv

  qiime diversity beta-group-significance \\
    --i-distance-matrix __NAME__/unweighted_unifrac_distance_matrix.qza \\
    --m-metadata-file $meta_file \\
    --m-metadata-column group \\
    --o-visualization __NAME__/unweighted-unifrac-group-significance.qzv \\
    --p-pairwise    

  qiime diversity alpha-rarefaction \\
    --i-table __FILE2__ \\
    --i-phylogeny __FILE__ \\
    --p-max-depth 4000 \\
    --m-metadata-file $meta_file \\
    --o-visualization __NAME__/alpha-rarefaction.qzv
else
  touch __NAME__.failed
fi

#__OUTPUT__
",
    use_tmp_folder => 0,
    suffix  => "_qd",
    docker_prefix => "qiime2_",
    interpretor => "",
    program => "",
    check_program => 0,
    source_arg => "--i-phylogeny",
    source_ref => [ $qiime2_phylogenetic, ".rooted-tree.qza" ],
    parameterSampleFile2_arg => "--i-table",
    parameterSampleFile2_ref => [ $qiime2_dada2, ".table.qza" ],
    output_arg => "--output-dir",
    output_file_prefix => "shannon_vector.qza",
    output_file_ext => "__NAME__/shannon_vector.qza",
    output_to_same_folder => 1,
    sh_direct   => getValue($def, "sh_direct", 0),
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    }
  };
  push (@$individual, $qiime2_diversity);

  my $qiime2_taxonomic = "qiime2_taxonomic";
  my $classifier_file = getValue($def, "classifier_file");
  $config->{$qiime2_taxonomic} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => $target_dir . "/" . getNextFolderIndex($def) . $qiime2_taxonomic,
    option => "
export CONDA_PREFIX=$cache_dir
export MPLCONFIGDIR=$cache_dir

if [[ -e __NAME__.failed ]]; then
  rm __NAME__.failed
fi

qiime feature-classifier classify-sklearn \\
  --i-classifier $classifier_file \\
  --p-n-jobs $thread \\
  --i-reads __FILE__ \\
  --o-classification __NAME__.taxonomy.qza \\
  --verbose | tee __NAME__.taxonomy.log

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.succeed

  qiime metadata tabulate \\
    --m-input-file __NAME__.taxonomy.qza \\
    --o-visualization __NAME__.taxonomy.qzv

  qiime taxa barplot \\
    --i-table __FILE2__ \\
    --i-taxonomy __NAME__.taxonomy.qza \\
    --m-metadata-file $meta_file \\
    --o-visualization __NAME__.taxa-bar-plots.qzv

else
  touch __NAME__.failed
fi

#__OUTPUT__
",
    use_tmp_folder => 0,
    suffix  => "_qt",
    docker_prefix => "qiime2_",
    interpretor => "",
    program => "",
    check_program => 0,
    source_arg => "--i-reads",
    source_ref => [ $qiime2_dada2, ".dada2.qza" ],
    parameterSampleFile2_arg => "--i-table",
    parameterSampleFile2_ref => [ $qiime2_dada2, ".table.qza" ],
    output_arg => "-o-classification",
    output_file_prefix => ".taxonomy.qza",
    output_file_ext => ".taxonomy.qza",
    output_to_same_folder => 1,
    sh_direct   => getValue($def, "sh_direct", 0),
    pbs => {
      "nodes"     => "1:ppn=$thread",
      "walltime"  => "40",
      "mem"       => "40gb"
    }
  };
  push (@$individual, $qiime2_taxonomic);

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $individual,
      step2 => $summary,
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

sub performQiime2 {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getQiime2Config($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performQiime2Task {
  my ( $def, $task ) = @_;

  my $config = getQiime2Config($def);

  performTask( $config, $task );

  return $config;
}

1;
