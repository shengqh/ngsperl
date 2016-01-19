#!/usr/bin/perl
package QC::QC3fastq;

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $qc3_perl = get_param_file( $config->{$section}{qc3_perl}, "qc3_perl", 1 );

  my $raw_files = get_raw_files( $config, $section );

  my $mapfile = $result_dir . "/${task_name}_sample.list";
  open( MAP, ">$mapfile" ) or die "Cannot create $mapfile";
  for my $sample_name ( sort keys %{$raw_files} ) {
    my @fastqFiles = @{ $raw_files->{$sample_name} };
    for my $fastq (@fastqFiles) {
      print MAP $fastq, "\t", $sample_name, "\n";
    }
  }
  close(MAP);

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log     = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  open( my $out, ">$pbs_file" ) or die $!;
  print $out "$pbs_desc
$log_desc

$path_file

cd $result_dir

echo QC3bam=`date`
 
perl $qc3_perl $option -m f -i $mapfile -o $result_dir -t $thread

echo finished=`date`

exit 0
";
  close $out;

  print "$pbs_file created. \n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result      = {};
  my @result_files = ();
  push( @result_files, $result_dir );
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
