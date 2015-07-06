#!/usr/bin/perl
package Alignment::Gsnap;

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
  $self->{_name}   = "Alignment::Gsnap";
  $self->{_suffix} = "_gs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $selfname = $self->{_name};

    print $option;
  if ( !( $option =~ /\-t/ ) ) {
    print "no thread defined";
    if ( $thread > 1 ) {
      $option = $option . " -t " . $thread;
    }
  }

  my $gsnap_index_directory = get_directory($config, $section, "gsnap_index_directory", 1);
  my $gsnap_index_name = get_option($config, $section, "gsnap_index_name");

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $resultFile     = $sampleName . ".txt";

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );
    
    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc
$path_file

cd $curDir

echo gsnap_start=`date`

if [ -s $resultFile ]; then
  echo job $selfname has already been done. if you want to do again, delete $resultFile and submit job again.
  exit 0
fi

gsnap $option -D $gsnap_index_directory -d $gsnap_index_name -o $resultFile";

  for my $sampleFile (@sampleFiles){
    print OUT " $sampleFile";
  }
  print OUT "
  
echo finished=`date`

exit 0;
";

    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $bamFile     = "${resultDir}/${sampleName}/${sampleName}.txt";
    my @resultFiles = ();
    push( @resultFiles, $bamFile );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
