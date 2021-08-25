#!/usr/bin/perl
package CNV::XHMM;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "CNV::XHMM";
  $self->{_suffix} = "_xhmm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $raw_files = get_raw_files( $config, $section );
  my $ref_fasta = get_param_file( $config->{$section}{ref_fasta},     "ref_fasta",     1 );

  my $inputOption = get_rawfiles_option( $raw_files, "--GATKdepths" );
  
  my $c2r_script = dirname(__FILE__) . "/ColumnToRow.py";
  if ( !-e $c2r_script ) {
    die "File not found : " . $c2r_script;
  }
  
  my $params_file = dirname(__FILE__) . "/XHMM.params";
  if ( !-e $params_file ) {
    die "File not found : " . $params_file;
  }

  `cp $params_file $result_dir`;
  
  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc   = $cluster->get_log_description($log);
  my $final_file = "${task_name}.vcf";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  
  my $rd_file = "${task_name}.RD.txt";
  my $rd_c2r_file = "${task_name}.RD.c2r.txt";
  my $rd_centered_file = "${task_name}.filtered_centered.RD.txt";
  my $rd_centered_filtered_targets_file = "${task_name}.filtered_centered.RD.txt.filtered_targets.txt";
  my $rd_centered_filtered_samples_file = "${task_name}.filtered_centered.RD.txt.filtered_samples.txt";
  my $rd_centered_pca_file = "${task_name}.RD_PCA";
  my $rd_centered_pca_normalized_file = "${task_name}.PCA_normalized.txt";
  my $rd_centered_pca_normalized_zscore_file = "${task_name}.PCA_normalized.filtered.sample_zscores.RD.txt";
  my $rd_centered_pca_normalized_zscore_filtered_targets_file = "${task_name}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt";
  my $rd_centered_pca_normalized_zscore_filtered_samples_file = "${task_name}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt";
  my $rd_filtered_file = "${task_name}.same_filtered.RD.txt";
  my $rd_filtered_c2r_file = "${task_name}.same_filtered.RD.c2r.txt";
  my $rd_xcnv_file = "${task_name}.xcnv";
  my $rd_aux_xcnv_file = "${task_name}.aux_xcnv";
  
  print $pbs "
if [[ ! -s $rd_file ]]; then
  xhmm --mergeGATKdepths -o $rd_file $inputOption
  python3 $c2r_script $rd_file $rd_c2r_file
fi

if [[ ! -s $rd_centered_file ]]; then
  xhmm --matrix -r $rd_file --centerData --centerType target \\
    -o $rd_centered_file \\
    --outputExcludedTargets $rd_centered_filtered_targets_file \\
    --outputExcludedSamples $rd_centered_filtered_samples_file \\
    --minTargetSize 10 --maxTargetSize 10000 \\
    --minMeanTargetRD 10 --maxMeanTargetRD 500 \\
    --minMeanSampleRD 25 --maxMeanSampleRD 200 \\
    --maxSdSampleRD 150
fi

if [[ ! -s $rd_centered_pca_file ]]; then
  xhmm --PCA -r $rd_centered_file --PCAfiles $rd_centered_pca_file
fi

if [[ ! -s $rd_centered_pca_normalized_file ]]; then
  xhmm --normalize -r $rd_centered_file --PCAfiles $rd_centered_pca_file \\
    --normalizeOutput $rd_centered_pca_normalized_file \\
    --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
fi
  
if [[ ! -s $rd_centered_pca_normalized_zscore_file ]]; then
  xhmm --matrix -r $rd_centered_pca_normalized_file --centerData --centerType sample --zScoreData \\
    -o $rd_centered_pca_normalized_zscore_file \\
    --outputExcludedTargets $rd_centered_pca_normalized_zscore_filtered_targets_file \\
    --outputExcludedSamples $rd_centered_pca_normalized_zscore_filtered_samples_file \\
    --maxSdTargetRD 30
fi

if [[ ! -s $rd_filtered_file ]]; then
  xhmm --matrix -r $rd_file \\
    --excludeTargets $rd_centered_filtered_targets_file \\
    --excludeTargets $rd_centered_pca_normalized_zscore_filtered_targets_file \\
    --excludeSamples $rd_centered_filtered_samples_file \\
    --excludeSamples $rd_centered_pca_normalized_zscore_filtered_samples_file \\
    -o $rd_filtered_file
  python3 $c2r_script $rd_filtered_file $rd_filtered_c2r_file
fi
  
if [[ ! -s $rd_xcnv_file ]]; then
  xhmm --discover -p XHMM.params \\
    -r $rd_centered_pca_normalized_zscore_file \\
    -R $rd_filtered_file \\
    -c $rd_xcnv_file \\
    -a $rd_aux_xcnv_file \\
    -s ${task_name}.depth
fi

xhmm --genotype -p XHMM.params \\
  -r $rd_centered_pca_normalized_zscore_file \\
  -R $rd_filtered_file \\
  -g $rd_xcnv_file \\
  -F $ref_fasta \\
  -v $final_file
  
#if [[ -s $final_file ]]; then
#  rm $rd_centered_file $rd_centered_pca_normalized_file $rd_centered_pca_normalized_zscore_file $rd_filtered_file $rd_xcnv_file $rd_aux_xcnv_file 
#fi

";

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result = { $task_name => [ $result_dir . "/${task_name}.vcf" ] };

  return $result;
}

1;
