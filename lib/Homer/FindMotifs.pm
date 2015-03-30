#!/usr/bin/perl
package Homer::FindMotifs;

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
  $self->{_name}   = "Homer::FindMotifs";
  $self->{_suffix} = "_fm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my %tagDirectories = %{ get_raw_files( $config, $section ) };

  my $pairs = get_raw_files( $config, $section );
  my $fasta  =get_param_file($config, $section, "genome_fasta", 1);

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $pairName ( sort keys %{$pairs} ) {
    my $pairFile = $pairs->{$pairName};
    
    my $pbsFile = $self->pbsfile( $pbsDir, $pairName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $pairName );

    my $finalFile = $pairName . ".motif";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir 

echo homer_FindMotifs_start=`date` 

if [ -s $finalFile ];then
  echo job has already been done. if you want to do again, delete ${resultDir}/${finalFile} and submit job again.
  exit 0;
fi

findMotifsGenome.pl $pairFile $fasta peakAnalysis -size 200 -len 8

echo homer_FindMotifs_finished=`date` 

exit 0

";

    close(OUT);

    print SH "\$MYCMD ./$pbsName \n";
    print "$pbsFile created\n";
  }
  print SH "exit 0\n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $pairs = get_raw_files( $config, $section );

  my $result = {};
  for my $pairName ( sort keys %{$pairs} ) {
    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${pairName}.motif" );
    $result->{$pairName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
