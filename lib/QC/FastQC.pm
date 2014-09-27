#!/usr/bin/perl
package QC::FastQC;

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
  $self->{_name}   = "FastQC";
  $self->{_suffix} = "_fq";
  bless $self, $class;
  return $self;
}

sub parsePairedSamples {
  my ($samples)   = @_;
  my @sampleFiles = @{$samples};
  my @result      = ();
  for my $sample (@sampleFiles) {
    if ( $sample =~ /,/ ) {
      my @files = split( ',', $sample );
      for my $file (@files) {
        push( @result, $file );
      }
    }
    else {
      push( @result, $sample );
    }
  }

  return @result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $summaryfile = $self->taskfile( $pbsDir, $task_name . "_summary" );
  open( SH, ">$summaryfile" ) or die "Cannot create $summaryfile";
  print SH "cd $resultDir
qcimg2pdf.sh -o $task_name
";
  close(SH);

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command(1);

  my $result = $self->result( $config, $section );

  for my $sampleName ( sort keys %rawFiles ) {
    my @originalFiles = @{ $rawFiles{$sampleName} };
    my @sampleFiles   = parsePairedSamples( \@originalFiles );

    my $sampleCount = scalar(@sampleFiles);
    my $samples     = join( ' ', @sampleFiles );
    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    my @expectresult = @{ $result->{$sampleName} };
    my $expectname   = $expectresult[0];

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo fastqc=`date`

if [ -e $expectname ]; then
  echo job has already been done. if you want to do again, delete ${expectname} and submit job again.
  exit 0;
fi

fastqc $option --extract -t $sampleCount -o $curDir $samples

echo finished=`date`
";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all FastQC tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @originalFiles = @{ $rawFiles{$sampleName} };
    my @sampleFiles   = parsePairedSamples( \@originalFiles );
    my @resultFiles   = ();
    for my $sampleFile (@sampleFiles) {
      my $name = basename($sampleFile);
      if ( $name =~ /gz$/ ) {
        $name = change_extension( $name, "" );
      }
      $name = change_extension( $name, "_fastqc" );
      push( @resultFiles, "${resultDir}/${sampleName}/${name}" );
    }
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
