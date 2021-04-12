#!/usr/bin/perl
package scRNA::Indrops;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_id";
  bless $self, $class;
  return $self;
}

#Based on https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/germline/cnv_germline_cohort_workflow.wdl
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  #parameter files
  $self->get_docker_value(1);

  my $indrops_script = get_param_file( $config->{$section}{indrops_script}, "indrops_script", 1, not $self->using_docker() );
  my $indrops_version = get_option( $config, $section, "indrops_version", "v2" );
  my $bowtie_index = get_option($config, $section, "bowtie_index");

  my $raw_files = get_raw_files( $config, $section );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name (sort keys %$raw_files){
    my $samples = $raw_files->{$sample_name};
    my $sample_file = $samples->[0];

    my $sample_folder = dirname($sample_file);
    $sample_file = basename($sample_file);
    $sample_file =~ s/R1/{read}/g;

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my $cdfile = $cur_dir . "/${sample_name}.yaml";
    open( my $cd, ">$cdfile" ) or die "Cannot create $cdfile";
    print $cd "
project_name : \"$task_name\"
project_dir : \"$cur_dir\"
paths :
  bowtie_index : \"$bowtie_index\"
sequencing_runs :
  - name : \"$sample_name\"
    version : \"$indrops_version\"
    dir : \"$sample_folder\"
    fastq_path : \"$sample_file\"
    library_name : \"$sample_name\"
parameters :
  umi_quantification_arguments :
    m : 10
    u : 1
    d : 600
    split-ambigs : False
    min_non_polyA : 15
  output_arguments :
    output_unaligned_reads_to_other_fastq : False
    low_complexity_mask : False
  bowtie_arguments :
    m : 200
    n : 1
    l : 15
    e : 1000
  trimmomatic_arguments :
    LEADING : \"28\"
    SLIDINGWINDOW : \"4:20\"
    MINLEN : \"16\"
";
    close($cd);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_file = $cur_dir . "/${sample_name}.txt";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );
    print $pbs "  
echo python3 $indrops_script ${sample_name}.yaml filter -l $sample_name
python3 $indrops_script ${sample_name}.yaml filter -l $sample_name

echo python3 $indrops_script ${sample_name}.yaml identify_abundant_barcodes -l $sample_name
python3 $indrops_script ${sample_name}.yaml identify_abundant_barcodes -l $sample_name

echo python3 $indrops_script ${sample_name}.yaml sort -l $sample_name
python3 $indrops_script ${sample_name}.yaml sort -l $sample_name

echo python3 $indrops_script ${sample_name}.yaml quantify -l $sample_name
python3 $indrops_script ${sample_name}.yaml quantify -l $sample_name
";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all GATK refine tasks. \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $raw_files = get_raw_files( $config, $section );

  my $result = {};
  for my $sample_name (sort keys %$raw_files){
    my $cur_dir = $result_dir . "/$sample_name";
    $result->{$sample_name} = filter_array( [ "${cur_dir}/${sample_name}.txt" ], $pattern );
  }

  return $result;
}

1;
