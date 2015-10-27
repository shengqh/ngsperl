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
  $self->{_name}   = "FastQCSummary";
  $self->{_suffix} = "_fqs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $pbsFile    = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName    = basename($pbsFile);
  my $log        = $self->logfile( $logDir, $task_name );
  my $cqstools   = get_param_file( $config->{$section}{cqstools}, "cqstools", 1 );
  my $fastqc_dir = get_directory( $config, $section, "fastqc_dir", 0 );
  if ( !defined $fastqc_dir ) {
    $fastqc_dir = $resultDir;
  }
  my $log_desc = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

qcimg2pdf.sh -o $task_name

mono $cqstools fastqc_summary -i $fastqc_dir -o ${task_name}.FastQC.summary.tsv 
";
  close(OUT);

  if ( is_linux() ) {
    chmod 0755, $pbsFile;
  }

  print "!!!shell file $pbsFile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result      = {};
  my @resultFiles = ();
  push( @resultFiles, "${resultDir}/${task_name}.FastQC.summary.reads.tsv" );
  push( @resultFiles, "${resultDir}/${task_name}.FastQC.summary.tsv" );
  push( @resultFiles, "${resultDir}/${task_name}.FastQC.overrepresented.tsv" );
  push( @resultFiles, "${resultDir}/${task_name}.FastQC.pdf" );
  $result->{$task_name} = filter_array( \@resultFiles, $pattern );
  return $result;
}

1;
