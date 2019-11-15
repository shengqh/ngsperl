#!/usr/bin/perl
package GATK4::VariantGather;

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
use GATK4::GATK4UniqueTask;

our @ISA = qw(GATK4::GATK4UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_vr";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);
  $self->get_docker_value(1);

  my %vcfFiles = %{ get_raw_files( $config, $section ) };
  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $script = dirname(__FILE__) . "/fixLeftTrimDeletion.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $pass_file = $task_name . ".indels.snp.recal.pass.vcf.gz";
  my $split_file = $task_name . ".indels.snp.recal.pass.split.vcf";
  my $left_trim_file = $task_name . ".indels.snp.recal.pass.norm.vcf";
  my $fix_file = $task_name . ".indels.snp.recal.pass.norm.nospan.vcf";
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

  my @sample_names = ();
  if (has_option($config, $section, "chromosome_names")){
    my $chromosomeStr = get_option($config, $section, "chromosome_names");
    my @chromosomes = split /,/, $chromosomeStr;
    for my $chr (@chromosomes) {
      my $chrTaskName = $task_name . "." . $chr;
      push @sample_names, $chrTaskName;
    }
  }else{
    @sample_names = sort keys %vcfFiles;
  }

  for my $sample_name ( @sample_names ) {
    my @sample_files = @{ $vcfFiles{$sample_name} };
    my $individualPassFile     = $sample_files[0];
    print $pbs "    -I $individualPassFile \\\n";
  }

  print $pbs "    -O $pass_file

  tabix -p vcf $pass_file
fi

if [[ -s $pass_file && ! -s $left_trim_file ]]; then
  echo LeftAlignAndNorm=`date`
  bcftools norm -m- -o $split_file $pass_file 
  bcftools norm -f $faFile -o $left_trim_file $split_file 
fi

if [[ -s $left_trim_file && ! -s $final_file ]]; then
  echo noSpanDeletion=`date`
  python $script -i $left_trim_file -o $fix_file
  bgzip $fix_file
  tabix -p vcf $final_file
fi

if [[ -s $final_file ]]; then
  rm -rf $split_file ${split_file}.idx $left_trim_file ${left_trim_file}.idx \\
    .conda
fi

";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $final_file = $result_dir . "/" . $task_name . ".indels.snp.recal.pass.norm.nospan.vcf.gz";

  my $result = {};
  $result->{$task_name} = filter_array( [$final_file], $pattern );

  return $result;
}

1;
