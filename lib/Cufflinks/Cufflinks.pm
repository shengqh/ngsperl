#!/usr/bin/perl
package Cufflinks::Cufflinks;

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
  $self->{_name}   = "Cufflinks";
  $self->{_suffix} = "_clinks";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );
  my $gtf = "";
  if ( defined $transcript_gtf ) {
    $gtf = "-g $transcript_gtf";
  }

  my %tophat2map = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";

  for my $sampleName ( sort keys %tophat2map ) {
    my @tophat2Files = @{ $tophat2map{$sampleName} };
    my $tophat2File  = $tophat2Files[0];

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );
    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    my $cluster = get_cluster( $config, $section );
    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir

if [ -s transcripts.gtf ];then
  echo job has already been done. if you want to do again, delete ${curDir}/transcripts.gtf and submit job again.
  exit 0;
fi

echo cufflinks=`date`
 
cufflinks $option $gtf -o . $tophat2File

echo finished=`date`

exit 0
";

    close(OUT);

    print "$pbsFile created. \n";

    print SH "\$MYCMD ./$pbsName \n";
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

  my %tophat2map = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( sort keys %tophat2map ) {
    my $curDir      = $resultDir . "/$sampleName";
    my @resultFiles = ();
    push( @resultFiles, $curDir . "/transcripts.gtf" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
