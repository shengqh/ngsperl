#!/usr/bin/perl
package QC::QC3vcf;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "QC3vcf";
  $self->{_suffix} = "_qc3";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $qc3_perl = get_param_file( $config->{$section}{qc3_perl}, "qc3_perl", 1 );
  my $annovarDB = get_directory( $config, $section, "annovar_db", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $vcffile     = $sampleFiles[0];

    my $resultFile = $resultDir . "/" . $sampleName;

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);
    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

echo QC3bam=`date`
 
perl $qc3_perl $option -m v -i $vcffile -a $annovarDB -o $resultFile -rp

echo finished=`date`

exit 0
";
    close(OUT);

    print "$pbsFile created. \n";
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

  my $rawFiles = get_raw_files( $config, $section );

  my $result = {};
  for my $sampleName ( sort keys %{$rawFiles} ) {

    my $resultFile  = $resultDir . "/" . $sampleName;
    my @resultFiles = ();
    push( @resultFiles, $resultFile );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;

  return $result;
}

1;
