#!/usr/bin/perl
package GATK4::GenotypeGVCFs;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::GATK4UniqueTask;
use Data::Dumper;
use Tie::IxHash;

our @ISA = qw(GATK4::GATK4UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gg";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $ref_fasta = get_option_file( $config, $section, "fasta_file", 1 );
  my $dbsnp_vcf = get_option_file( $config, $section, "dbsnp_vcf", 1 );
  my $target_intervals_file = get_option($config, $section, "target_intervals_file"),

  my $import_folders = get_raw_files( $config, $section );

  my $WORKSPACE = $import_folders->{$task_name}[0];

  my $output_vcf_filename    = $task_name . ".g.vcf.gz";
  my $output_vcf_filename_index = $output_vcf_filename . ".tbi";

  my $tmp_file    = $task_name . ".tmp.g.vcf.gz";
  my $tmp_index = $tmp_file . ".tbi";

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $output_vcf_filename_index );

  print $pbs "
gatk --java-options -Xms8g \\
  GenotypeGVCFs $option \\
  -R ${ref_fasta} \\
  -O ${tmp_file} \\
  -D ${dbsnp_vcf} \\
  -G StandardAnnotation -G AS_StandardAnnotation \\
  --only-output-calls-starting-in-intervals \\
  -V gendb://$WORKSPACE \\
  -L $target_intervals_file \\
  --merge-input-intervals

if [[ -s $tmp_index ]]; then
  mv $tmp_file $output_vcf_filename
  mv $tmp_index $output_vcf_filename_index
fi
";

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $result = {};
  $result->{$task_name} = filter_array( ["${result_dir}/${task_name}.g.vcf.gz"], $pattern );
  return $result;
}

1;
