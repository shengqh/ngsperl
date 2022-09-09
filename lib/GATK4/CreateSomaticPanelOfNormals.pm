#!/usr/bin/perl
package GATK4::CreateSomaticPanelOfNormals;

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

our @ISA = qw(GATK4::GATK4UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pon";
  $self->{_docker_prefix}   = "gatk4_";
  bless $self, $class;

  return $self;
}

sub get_sample_names {
  my ($self, $config, $section) = @_;
  my ( $task_name ) = $self->init_parameter( $config, $section, 0 );
  return [$task_name];
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);

  my $gvcf_files = get_raw_files( $config, $section );
  my $target_intervals_file = get_option($config, $section, "intervals");
  my $ref_fasta = get_option_file($config, $section, "fasta_file");

  my $germline_resource = get_param_file( $config->{$section}{germline_resource}, "germline_resource", 0 );
  if(defined $germline_resource and ($germline_resource ne "")){
    $option = $option . " --germline-resource $germline_resource ";
  }

  my $sample_name_map_file = $self->get_file( $result_dir, $task_name, ".sample_name_map.txt" );
  open( my $sf, ">$sample_name_map_file" ) or die "Cannot create $sample_name_map_file";
  for my $sample_name ( sort keys %$gvcf_files ) {
    my $gvcf_file = $gvcf_files->{$sample_name}[0];
    print $sf $sample_name . "\t" . $gvcf_file . "\n";
  }
  close($sf);

  my $prefix = $task_name;
  
  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $prefix );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $prefix );

  my $final_folder = $result_dir . "/$prefix";

  my $log_desc = $cluster->get_log_description($log);
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_folder );
  print $pbs "
rm -rf $prefix

gatk --java-options $java_option \\
    GenomicsDBImport \\
    --genomicsdb-workspace-path $prefix \\
    --batch-size 50 \\
    --sample-name-map $sample_name_map_file \\
    --reader-threads 5 \\
    --merge-input-intervals \\
    -L $target_intervals_file \\
    --consolidate

gdb_exit_code=\$?
if [[ \$gdb_exit_code -eq 0 ]]; then
  gatk --java-options $java_option \\
    CreateSomaticPanelOfNormals $option \\
    -R ${ref_fasta} \\
    -V gendb://$prefix \\
    -O ${task_name}.pon.vcf.gz
  pon_exit_code=\$?
  if [[ \$pon_exit_code -ne 0 ]]; then
    touch $prefix.pon.failed
  fi
else
  touch $prefix.GenomicsDBImport.failed
fi

rm -rf $prefix
";
    
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $files = [ "${result_dir}/${task_name}.pon.vcf.gz", "${result_dir}/${task_name}.pon.vcf.gz.tbi" ];
  my $result = { $task_name => filter_array( $files, $pattern ) };

  return $result;
}

1;
