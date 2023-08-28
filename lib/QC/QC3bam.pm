#!/usr/bin/perl
package QC::QC3bam;

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
  $self->{_suffix} = "_qc3";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $r_script = dirname(__FILE__) . "/QC3bam.r";

  my $target_region_bed = get_param_file( $config->{$section}{target_region_bed}, "target_region_bed", 0 );
  if ( defined $target_region_bed ) {
    $option = $option . " -r $target_region_bed";
    my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );
    if(defined $transcript_gtf){
      $option = $option . " -g $transcript_gtf";
    }
  }
  else {
    my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 1 );
    $option = $option . " -g $transcript_gtf";
  }
  my $qc3_perl = get_param_file( $config->{$section}{qc3_perl}, "qc3_perl", 1, not $self->using_docker() );

  my $raw_files = get_raw_files( $config, $section );

  my $mapfile = $result_dir . "/${task_name}_sample.list";
  open( MAP, ">$mapfile" ) or die "Cannot create $mapfile";
  for my $sample_name ( sort keys %{$raw_files} ) {
    my @bam_files = @{ $raw_files->{$sample_name} };
    for my $bam (@bam_files) {
      print MAP $bam, "\t", $sample_name, "\n";
    }
  }
  close(MAP);

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $final_file = "bamResult/bamSummary.txt";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  print $pbs "
perl $qc3_perl $option -t $thread -m b -i $mapfile -o $result_dir

R --vanilla -f $r_script --args $result_dir $task_name
";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  push( @result_files, $result_dir . "/bamResult/bamSummary.txt" );
  push( @result_files, $result_dir . "/bamFigure/${task_name}.readsByCategory.png" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
