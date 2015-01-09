#!/usr/bin/perl
package CQS::FastqTrimmer;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "CQS::FastqTrimmer";
  $self->{_suffix} = "_ft";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );
  my $extension = get_option( $config, $section, "extension" );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );
    
    print SH "\$MYCMD ./$pbsName \n";

    my $cluster = get_cluster( $config, $section );
    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

";
    for my $sampleFile (@sampleFiles) {
      my $trimFile = get_trim_file( $sampleFile, $extension );
      print OUT "if [ ! -s $trimFile ]; then
  mono-sgen $cqstools fastq_trimmer $option -i $sampleFile -o $trimFile 
fi
";
    }
    print OUT "
echo finished=`date`

exit 0 
";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->name() . " tasks.\n";

  #`qsub $pbsFile`;
}

sub get_trim_file {
  my ( $sampleFile, $extension ) = @_;
  my $fileName = basename($sampleFile);
  if ( $fileName =~ /.gz$/ ) {
    $fileName = change_extension( $fileName, "" );
  }
  $fileName = change_extension( $fileName, "" );
  my $result = $fileName . $extension;
  return ($result);
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $extension = get_option( $config, $section, "extension" );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my @resultFiles = ();

    for my $sampleFile (@sampleFiles) {
      my $trimFile = get_trim_file( $sampleFile, $extension );
      push( @resultFiles, "${resultDir}/${trimFile}" );
    }
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
