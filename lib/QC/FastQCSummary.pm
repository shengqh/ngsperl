#!/usr/bin/perl
package QC::FastQCSummary;

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
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_fqs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $pbs_file    = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name    = basename($pbs_file);
  my $log        = $self->get_log_filename( $log_dir, $task_name );
  my $cqstools   = get_param_file( $config->{$section}{cqstools}, "cqstools", 1 );
  my $fastqc_dir = get_directory( $config, $section, "fastqc_dir", 0 );
  if ( !defined $fastqc_dir ) {
    $fastqc_dir = $result_dir;
  }
  my $log_desc = $cluster->get_log_description($log);

  open( my $out, ">$pbs_file" ) or die $!;
  print $out "$pbs_desc
$log_desc

$path_file

cd $result_dir

if [ -s ${task_name}.FastQC.summary.tsv ]; then
  echo job has already been done. if you want to do again, delete ${result_dir}/${task_name}.FastQC.summary.tsv and submit job again.
  exit 0;
fi

qcimg2pdf.sh -o $task_name

mono $cqstools fastqc_summary -i $fastqc_dir -o ${task_name}.FastQC.summary.tsv 
";
  close $out;

  if ( is_linux() ) {
    chmod 0755, $pbs_file;
  }

  print "!!!shell file $pbs_file created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result      = {};
  my @result_files = ();
  push( @result_files, "${result_dir}/${task_name}.FastQC.reads.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.summary.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.overrepresented.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.pdf" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
