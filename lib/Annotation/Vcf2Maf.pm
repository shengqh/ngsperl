#!/usr/bin/perl
package Annotation::Vcf2Maf;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_v2m";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $vcf2mgf = get_option_file($config, $section, "vcf2maf_pl");
  my $vep_path = get_option($config, $section, "vep_path");
  my $vep_data = get_option($config, $section, "vep_data");
  my $species = get_option($config, $section, "species");
  my $ncbi_build = get_option($config, $section, "ncbi_build");
  my $ref_fasta = get_option_file($config, $section, "ref_fasta");
  my $filter_vcf;
  if (defined $config->{$section}{filter_vcf}){
    $filter_vcf = "--filter-vcf " . get_option_file($config, $section, "filter_vcf");
  }else{
    $filter_vcf = "--filter-vcf 0"
  }

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $vcf_file     = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $final_file = basename($vcf_file) . ".maf";
    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

    print $pbs "
if [ ! -s $final_file ]; then
  perl $vcf2mgf --vep-forks $thread --input-vcf $vcf_file --output-maf $final_file --vep-path $vep_path --vep-data $vep_data --species $species --ncbi-build $ncbi_build --ref-fasta $ref_fasta $filter_vcf
fi
";
    $self->close_pbs( $pbs, $pbs_file );
  }

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $vcf_file     = $sample_files[0];
    my $final_file = basename($vcf_file) . ".maf";

    my @result_files = ();
    push( @result_files, $result_dir . "/" . $final_file );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
