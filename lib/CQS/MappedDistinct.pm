#!/usr/bin/perl
package CQS::CQSMappedCount;

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
  $self->{_name} = "MappedDistinct";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_param_file( $config->{$section}{cqs_tools}, "cqs_tools", 1 );

  my %firstFiles = %{ get_raw_files( $config, $section ) };
  my %secondFiles = %{ get_raw_files( $config, $section, "second" ) };
  my $firstSuffix  = $config->{$section}{first_suffix}  or die "define ${section}::first_suffix first";
  my $secondSuffix = $config->{$section}{second_suffix} or die "define ${section}::second_suffix first";

  my $shfile = $pbsDir . "/${task_name}_dt.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "cd $resultDir
echo CQSMappedDistinct=`date`
";

  for my $sampleName ( sort keys %firstFiles ) {
    my $firstFile  = $firstFiles{$sampleName}->[0];
    my $secondFile = $secondFiles{$sampleName}->[0];
    print SH "mono-sgen $cqsFile mapped_distinct $option --file1 $firstFile --name1 $firstSuffix --file2 $secondFile --name2 $secondSuffix -o .
";
  }
  print SH "
echo finished=`date`
exit 1
";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all trna_count tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $fasta_format = $config->{$section}{fasta_format};
  if ( !defined $fasta_format ) {
    $fasta_format = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $curDir = $resultDir . "/$sampleName";

    my @bamFiles = @{ $rawFiles{$sampleName} };
    my $bamFile  = $bamFiles[0];
    my $fileName = basename($bamFile);

    my @resultFiles = ();
    my $countFile   = "${curDir}/${fileName}.count";
    push( @resultFiles, $countFile );
    push( @resultFiles, "${countFile}.mapped.xml" );
    push( @resultFiles, "${curDir}/${fileName}.info" );

    my $unmapped;
    if ($fasta_format) {
      $unmapped = change_extension( $countFile, ".unmapped.fasta" );
    }
    else {
      $unmapped = change_extension( $countFile, ".unmapped.fastq" );
    }
    push( @resultFiles, $unmapped );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
