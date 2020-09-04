#!/usr/bin/perl
package GATK4::VariantFilterGather;

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
use CQS::GatherTask;
use CQS::TaskUtils;
use GATK4::VariantFilterUtils;

our @ISA = qw(CQS::GatherTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_vr";
  $self->{_docker_prefix} = "gatk4_";
  $self->{_export_home} = 1;
  bless $self, $class;
  return $self;
}

#return a gather_name/scatter_name_list map
sub get_gather_map {
  my ($self, $config, $section) = @_;
  my ($task_name) = get_task_name($config, $section);
  
  my $vcf_files = get_raw_files($config, $section);
  
  my @sample_names = ();
  if (has_option($config, $section, "chromosome_names")){
    my $chromosomeStr = get_option($config, $section, "chromosome_names");
    my @chromosomes = split /,/, $chromosomeStr;
    for my $chr (@chromosomes) {
      my $chrTaskName = get_key_name($task_name, $chr);
      push @sample_names, $chrTaskName;
    }
  }else{
    @sample_names = sort keys %$vcf_files;
  }
  
  return {$task_name => \@sample_names};  
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $gather_name ) = @_;
  my $final_file = $result_dir . "/" . $gather_name . ".indels.snp.recal.pass.norm.nospan.vcf.gz";
  return [$final_file];
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);
  $self->get_docker_value(1);

  my %vcfFiles = %{ get_raw_files( $config, $section ) };

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $pass_file = $task_name . ".indels.snp.recal.pass.vcf.gz";
  my $final_file = $task_name . ".indels.snp.recal.pass.norm.nospan.vcf.gz";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );
  print $pbs "  
if [ ! -s $pass_file ]; then
  echo GatherVcfsCloud=`date` 
  gatk --java-options \"$java_option\" \\
    GatherVcfsCloud \\
    --ignore-safety-checks \\
    --gather-type BLOCK \\
";

  my $sample_names = $self->get_gather_map($config, $section)->{$task_name}; 
  for my $sample_name ( @$sample_names ) {
    my @sample_files = @{ $vcfFiles{$sample_name} };
    my $individualPassFile     = $sample_files[0];
    print $pbs "    -I $individualPassFile \\\n";
  }

  print $pbs "    -O $pass_file

  tabix -p vcf $pass_file
fi
";

  my ($rmlist) = add_left_trim_pbs($self, $config, $section, $pbs, $task_name, $pass_file, $final_file);

print $pbs "
if [[ -s $final_file ]]; then
  rm -rf $rmlist \\
    .conda
fi

";
  $self->close_pbs( $pbs, $pbs_file );
}

1;
