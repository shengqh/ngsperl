#!/usr/bin/perl
package GATK::UnifiedGenotyperCombine;

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
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ug";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $extension = get_option( $config, $section, "extension", ".vcf" );

  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar}, "gatk_jar", 1 );

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option || $java_option eq "" ) {
    $java_option = "-Xmx${memory}";
  }

  my $bedFile = get_param_file( $config->{$section}{bed_file}, "bed_file", 0 );
  my $interval_padding = get_option( $config, $section, "interval_padding", 0 );
  my $restrict_intervals = "";
  if ( defined $bedFile and $bedFile ne "" ) {
    if ( defined $interval_padding and $interval_padding != 0 ) {
      $restrict_intervals = "-L $bedFile -ip $interval_padding";
    }
    else {
      $restrict_intervals = "-L $bedFile";
    }
  }

  my %bam_files = %{ get_raw_files( $config, $section ) };

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $snvOut = $task_name . $extension;

  #if the program throw exception, the idx file will not be generated.
  my $snvOutIndex = $snvOut . ".idx";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $snvOutIndex );
  print $pbs "java $java_option -jar $gatk_jar -T HaplotypeCaller $option -R $faFile -nct $thread --out $snvOut $restrict_intervals \\\n";
  for my $sample_name ( sort keys %bam_files ) {
    my @sample_files = @{ $bam_files{$sample_name} };
    my $bam_file     = $sample_files[0];
    print $pbs "  -I $bam_file \\\n";

  }
  $self->close_pbs( $pbs, $pbs_file );

  print "!!!shell file $pbs_file created, you can run this pbs file to do task.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  my $result = {};

  my $extension = get_option( $config, $section, "extension", ".vcf" );
  my @result_files = ();
  push( @result_files, "${result_dir}/${task_name}.${extension}" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
