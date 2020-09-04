#!/usr/bin/perl
package Oasis::Compress;

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
  $self->{_suffix} = "_oc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my %raw_files     = %{ get_raw_files( $config, $section ) };
  my $compressor    = dirname(__FILE__) . "/OasisCompressor.jar";

  my $final_file = $task_name . ".zip";
  my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name   = basename($pbs_file);
  my $log        = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc   = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  print $pbs "
if [ -d $task_name ]; then
  rm -rf $task_name
fi
mkdir $task_name

if [ -d ${task_name}_tmp ]; then
  rm -rf ${task_name}_tmp
fi
mkdir ${task_name}_tmp

cd ${task_name}_tmp
";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample       = $sample_files[0];
    print $pbs "
cp $sample ${sample_name}.fastq.gz
";
  }
  print $pbs "
cd ..
java -jar $compressor -inputDir ${task_name}_tmp -outputDir $task_name -N $thread
rm -rf ${task_name}_tmp
";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $final_file = $task_name . ".zip";

  my $result       = {};
  my @result_files = ();

  push( @result_files, "${result_dir}/$final_file" );

  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}
1;
