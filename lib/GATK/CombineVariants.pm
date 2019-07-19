#!/usr/bin/perl
package GATK::CombineVariants;

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

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cv";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $faFile   = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1, not $self->using_docker() );
  my $extension = get_option( $config, $section, "extension", ".vcf" );

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option || $java_option eq "" ) {
    $java_option = "-Xmx${memory}";
  }

  my %vcfFiles = %{ get_raw_files( $config, $section ) };

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $final_file = $task_name . $extension;
  my $merged_file = $task_name . ".tmp.vcf";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  print $pbs "
java $java_option -jar $gatk_jar \\
  -T CombineVariants \\
  -R $faFile \\
";

  for my $sample_name ( sort keys %vcfFiles ) {
    my @sample_files = @{ $vcfFiles{$sample_name} };
    my $vcfFile      = $sample_files[0];
    print $pbs "  --variant $vcfFile \\\n";
  }
  
  print $pbs "  -o $merged_file \\
  -genotypeMergeOptions UNIQUIFY
  
  if [[ -s $merged_file ]]; then
    grep \"^##\" $merged_file > $final_file
    grep -v \"^##\" $merged_file | grep \"^#\" | sed -e \"s/.variant\\S*//g\" >> $final_file
    grep -v \"^#\" $merged_file |sed -e \"s/;set=variant\\S*//g\" >> $final_file
    rm $merged_file ${merged_file}.idx
  fi 
";


  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $extension = get_option( $config, $section, "extension", ".vcf" );

  my @result_files = ();
  my $merged_file  = $task_name . $extension;

  push( @result_files, $result_dir . "/" . $merged_file );
  my $result = {};
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
