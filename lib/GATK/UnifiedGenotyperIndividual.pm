#!/usr/bin/perl
package GATK::UnifiedGenotyperIndividual;

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
  $self->{_suffix} = "_ug";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $extension = get_option( $config, $section, "extension", ".vcf" );

  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar}, "gatk_jar", 1, not $self->using_docker() );

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

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %bam_files ) {
    my @sample_files = @{ $bam_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $snvOut = $sample_name . $extension;

    #if the program throw exception, the idx file will not be generated.
    my $snvOutIndex = $snvOut . ".idx";
    my $snvOutTmp   = $sample_name . ".tmp" . $extension;
    my $snvStat     = $sample_name . ".stat";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $snvOutIndex );

    print $pbs "echo Processing $sample_name \n";

    print $pbs "java $java_option -jar $gatk_jar -T HaplotypeCaller $option -R $faFile -I $bam_file -nct $thread --out $snvOut $restrict_intervals
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
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $result = {};

  my $extension = get_option( $config, $section, "extension", ".vcf" );

  my %bam_files = %{ get_raw_files( $config, $section ) };
  for my $sample_name ( sort keys %bam_files ) {
    my $snvOut       = $sample_name . $extension;
    my @result_files = ();
    push( @result_files, "${result_dir}/${snvOut}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
