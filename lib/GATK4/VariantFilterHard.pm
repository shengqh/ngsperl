#!/usr/bin/perl
package GATK4::VariantFilterHard;

use strict;
use warnings;
use File::Basename;
use List::Util qw[min];
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::GATK4Task;
use GATK4::VariantFilterUtils;

our @ISA = qw(GATK4::GATK4Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_vfh";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);
  $self->get_docker_value(1);

  #https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4-local.wdl
  my $excess_het_threshold = get_option($config, $section, "excess_het_threshold", 54.69);

  my $is_sample_size_small = get_option($config, $section, "is_sample_size_small", 0);

  my $ExcessHet_filter = "";
  my $InbreedingCoeff_filter = "";
  if(! $is_sample_size_small){
    $ExcessHet_filter = "-filter \"ExcessHet > ${excess_het_threshold}\" --filter-name \"ExcessHet${excess_het_threshold}\" ";
    $InbreedingCoeff_filter = "-filter \"InbreedingCoeff < -0.8\" --filter-name \"InbreedingCoeff-0.8\" ";
  }

  my $vcf_files = get_raw_files( $config, $section );

  my $script = dirname(__FILE__) . "/fixLeftTrimDeletion.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $interval_name (sort keys %$vcf_files) {

    my $cur_dir = create_directory_or_die( $result_dir . "/$interval_name" );

    my $vcf = $vcf_files->{$interval_name}[0];
    my $final_file = $interval_name . ".indels.snp.hardfilter.pass.vcf.gz";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $interval_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $interval_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file, $init_command );
    print $pbs "  

gatk SelectVariants \\
  -V ${vcf} \\
  -select-type SNP \\
  -O snps.vcf.gz

gatk SelectVariants \\
  -V ${vcf} \\
  -select-type INDEL \\
  -select-type MIXED \\
  -O indels.vcf.gz

gatk VariantFiltration \\
  -V snps.vcf.gz \\
  -filter \"QD < 2.0\" --filter-name \"QD2\" \\
  -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \\
  -filter \"SOR > 3.0\" --filter-name \"SOR3\" \\
  -filter \"FS > 60.0\" --filter-name \"FS60\" \\
  -filter \"MQ < 40.0\" --filter-name \"MQ40\" \\
  -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \\
  -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" $ExcessHet_filter $ExcessHet_filter \\
  -O snps_filtered.vcf.gz

gatk VariantFiltration \\
  -V indels.vcf.gz \\
  -filter \"QD < 2.0\" --filter-name \"QD2\" \\
  -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \\
  -filter \"FS > 200.0\" --filter-name \"FS200\" \\
  -filter \"ReadPosRankSum < -20.0\" --filter-name \"ReadPosRankSum-20\" $InbreedingCoeff_filter  \\
  -O indels_filtered.vcf.gz

gatk MergeVcfs \\
  -I snps_filtered.vcf.gz \\
  -I indels_filtered.vcf.gz \\
  -O merged_filtered.vcf.gz

gatk SelectVariants \\
  -O $final_file \\
  -V merged_filtered.vcf.gz \\
  --exclude-filtered
";

    print $pbs "
if [[ -s $final_file ]]; then
  rm snps.vcf.gz snps.vcf.gz.tbi \\
      indels.vcf.gz indels.vcf.gz.tbi \\
      snps_filtered.vcf.gz snps_filtered.vcf.gz.tbi \\
      indels_filtered.vcf.gz indels_filtered.vcf.gz.tbi \\
      merged_filtered.vcf.gz merged_filtered.vcf.gz.tbi
fi

  ";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks. \n ";

}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result = {};
  my $vcf_files = get_raw_files( $config, $section );
  for my $interval_name (sort keys %$vcf_files) {
    my $final_file = $interval_name . ".indels.snp.hardfilter.pass.vcf.gz";
    $result->{$interval_name} = filter_array( ["$result_dir/$interval_name/$final_file"], $pattern );
  }

  return $result;
}

1;
