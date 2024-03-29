#!/usr/bin/perl
package QC::MultiQC;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;
use Pipeline::PipelineUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mqc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $root_dir = get_option( $config, $section, "root_dir" );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);
  
  my $final = $task_name . ".html";
  my $pbs      = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  print $pbs "
if [[ -s ${task_name}.html ]]; then
  rm -rf ${task_name}_data ${task_name}_plots ${task_name}.html
fi

multiqc $option -f -p -n $final $root_dir 
";

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $result = {};

  my @result_files = ("$result_dir/${task_name}.html");
  if(defined $config->{fastqc_raw}){
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_fastqc_adapter_content_plot_1.png");    
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_fastqc_overrepresented_sequences_plot_1.png");    
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_fastqc_per_base_n_content_plot_1.png");    
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_fastqc_per_base_sequence_quality_plot_1.png");    
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_fastqc_per_sequence_gc_content_plot_Counts.png");    
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_fastqc_per_sequence_gc_content_plot_Percentages.png");    
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_fastqc_per_sequence_quality_scores_plot_1.png");    
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_fastqc_sequence_duplication_levels_plot_1.png");    
  }
  if(defined $config->{star} || defined $config->{star_featurecount}){
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_star_alignment_plot_1_pc.png");    
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_star_alignment_plot_1.png");    
  }
  if(defined $config->{featurecount} || defined $config->{star_featurecount}){
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_featureCounts_assignment_plot_1_pc.png");    
    push(@result_files, "${result_dir}/${task_name}_plots/png/mqc_featureCounts_assignment_plot_1.png");    
    push(@result_files, "${result_dir}/${task_name}_data/multiqc_featureCounts.txt");    
  }
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
