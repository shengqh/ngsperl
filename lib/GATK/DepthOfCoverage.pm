#!/usr/bin/perl
package GATK::DepthOfCoverage;

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
  $self->{_suffix} = "_doc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $intervals = parse_param_file( $config, $section, "interval_file", 1 );
  my $ref_fasta = get_param_file( $config->{$section}{ref_fasta},     "ref_fasta",     1 );
  my $gatk_jar  = get_param_file( $config->{$section}{gatk_jar},      "gatk_jar",      1, $self->using_docker() );

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option || $java_option eq "" ) {
    $java_option = "-Xmx${memory}";
  }

  my %bam_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %bam_files ) {
    my @sample_files = @{ $bam_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $final_file = $sample_name . ".depth";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $pbs "echo Processing $sample_name \n";

    print $pbs "java $java_option -jar $gatk_jar -T DepthOfCoverage $option \\
  -I $bam_file \\
  -L $intervals \\
  -R $ref_fasta \\
  -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \\
  --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \\
  --includeRefNSites \\
  --countType COUNT_FRAGMENTS \\
  -o $final_file
";

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result = {};
  my %bam_files = %{ get_raw_files( $config, $section ) };
  for my $sample_name ( sort keys %bam_files ) {
    my $result_files = ["${result_dir}/${sample_name}.depth.sample_interval_summary",
    "${result_dir}/${sample_name}.depth.sample_interval_statistics",
    "${result_dir}/${sample_name}.depth.sample_statistics",
    "${result_dir}/${sample_name}.depth.sample_summary"];
    $result->{$sample_name} = filter_array( $result_files, $pattern );
  }
  return $result;
}

1;
