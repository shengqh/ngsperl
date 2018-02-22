#!/usr/bin/perl
package Bedtools::Merge;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_merge";
  bless $self, $class;
  return $self;
}

sub get_final_file {
  my ( $sample_name, $blacklistfile, $shiftPosition ) = @_;
  my $result = $sample_name;
  if ( defined $blacklistfile ) {
    $result = $sample_name . ".confident";
  }

  if ($shiftPosition) {
    $result = $sample_name . ".shifted";
  }
  $result = $result . ".bed";
  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $raw_files = get_raw_files( $config, $section );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $final_file = $task_name . ".bed";
  my $gff3_file = $task_name . ".gff3";
  
  my $rscript = dirname(__FILE__) . "/bed2gff3.r";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  my @bedfiles = ();
  for my $sample_name ( sort keys %$raw_files ) {
    my $sample_files = $raw_files->{$sample_name};
    push @bedfiles, @$sample_files;
  }

  my $fileOption = join(" ", @bedfiles);
  my $tempFile = $task_name . ".tmp.bed";
  print $pbs "
cat $fileOption | sort -k1,1 -k2,2n > $tempFile
bedtools merge -i $tempFile -c 4 -o collapse > $final_file
R --vanilla -f $rscript --args $final_file $gff3_file
rm $tempFile 
";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  my $final_file   = $task_name . ".bed";
  my $gff3_file = $task_name . ".gff3";
  push( @result_files, "${result_dir}/${final_file}" );
  push( @result_files, "${result_dir}/${gff3_file}" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
