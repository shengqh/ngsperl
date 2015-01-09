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

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "FastQCSummary";
  $self->{_suffix} = "_fqs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $summaryfile = $self->taskfile( $pbsDir, $task_name . "_summary" );
  open( SH, ">$summaryfile" ) or die "Cannot create $summaryfile";
  print SH "cd $resultDir
qcimg2pdf.sh -o $task_name
";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $summaryfile;
  }

  print "!!!shell file $summaryfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result      = {};
  my @resultFiles = ();
  push( @resultFiles, "${resultDir}/${task_name}.pdf" );
  $result->{$task_name} = filter_array( \@resultFiles, $pattern );
  return $result;
}

1;
