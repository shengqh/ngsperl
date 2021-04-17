#!/usr/bin/perl
package QC::FastQCSummary;

use strict;
use warnings;
use File::Basename;
use File::Copy;
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name   = basename($pbs_file);
  my $log        = $self->get_log_filename( $log_dir, $task_name );

  my $python_script = dirname(__FILE__) . "/fastQCSummary.py";
  if ( !-e $python_script ) {
    die "File not found : " . $python_script;
  }

  my $r_script = dirname(__FILE__) . "/fastQCSummary.r";
  if ( !-e $r_script ) {
    die "File not found : " . $r_script;
  }

  my $rmd_script = dirname(__FILE__) . "/fastQCSummary.Rmd";
  if ( !-e $rmd_script ) {
    die "File not found : " . $rmd_script;
  }

  my $rmdfile = $task_name . ".FastQC.Rmd";
  copy($rmd_script, $result_dir . "/" . $rmdfile) or die "Copy failed: $!";
  
  my $fastqc_file_list = "fileList1.txt";
  save_parameter_sample_file( $config, $section, "source", "${result_dir}/$fastqc_file_list" );

  my $fastqc_dir = get_directory( $config, $section, "fastqc_dir", 0 );
  if ( !defined $fastqc_dir ) {
    $fastqc_dir = $result_dir;
  }
  my $log_desc = $cluster->get_log_description($log);

  my $final_file = "${task_name}.FastQC.html";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  print $pbs "
#qcimg2pdf.sh -o $task_name

python3 $python_script -i $fastqc_file_list -o ${task_name}.FastQC

R --vanilla -f $r_script --args ${task_name}.FastQC $rmdfile
";
  $self->close_pbs( $pbs, $pbs_file );

  if ( is_linux() ) {
    chmod 0755, $pbs_file;
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  push( @result_files, "${result_dir}/${task_name}.FastQC.summary.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.summary.tsv.png" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.reads.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.reads.tsv.png" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.baseQuality.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.baseQuality.tsv.png" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.sequenceGC.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.sequenceGC.tsv.png" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.adapter.tsv" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.adapter.tsv.png" );
  push( @result_files, "${result_dir}/${task_name}.FastQC.overrepresented.tsv" );
#  push( @result_files, "${result_dir}/${task_name}.FastQC.pdf" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
