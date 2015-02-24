#!/usr/bin/perl
package CQS::SmallRNACategory;

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
  $self->{_name}   = "SmallRNACategory";
  $self->{_suffix} = "_cat";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqsFile = get_cqstools( $config, $section, 1 );

  my %rawFiles   = %{ get_raw_files( $config, $section ) };
  my $hasMiRNA   = 0;
  my %miRNAFiles = ();
  if ( defined $config->{$section}{"mirna_count"} || defined $config->{$section}{"mirna_count_ref"} ) {
    %miRNAFiles = %{ get_raw_files( $config, $section, "mirna_count" ) };
    $hasMiRNA = 1;
  }

  my $ispdf = $config->{$section}{pdfgraph};
  if ( !defined $ispdf ) {
    $ispdf = 0;
  }
  my $pdfoption = "";
  if ($ispdf) {
    $pdfoption = "--pdf";
  }

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my $filelist = $self->getfile( $pbsDir, $task_name, ".filelist" );
  open( FL, ">$filelist" ) or die "Cannot create $filelist";

  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $groupName ( sort keys %{$groups} ) {
      my @samples = @{ $groups->{$groupName} };
      for my $sampleName ( sort @samples ) {
        my @smallRNAFiles = @{ $rawFiles{$sampleName} };
        my $smallRNAFile  = $smallRNAFiles[0];

        my $mirnaFileOption = "";
        if ($hasMiRNA) {
          my @miRNAs = @{ $miRNAFiles{$sampleName} };
          my $miRNA  = $miRNAs[0];

          $mirnaFileOption = "\t" . $miRNA;
        }

        print FL $groupName, "\t", $sampleName, "\t", $smallRNAFile, $mirnaFileOption, "\n";
      }
    }
  }
  else {
    for my $sampleName ( sort keys %rawFiles ) {
      my @smallRNAFiles = @{ $rawFiles{$sampleName} };
      my $smallRNAFile  = $smallRNAFiles[0];

      my $mirnaFileOption = "";
      if ($hasMiRNA) {
        my @miRNAs = @{ $miRNAFiles{$sampleName} };
        my $miRNA  = $miRNAs[0];

        $mirnaFileOption = "\t" . $miRNA;
      }

      print FL $task_name, "\t", $sampleName, "\t", $smallRNAFile, $mirnaFileOption, "\n";
    }
  }
  close(FL);

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  my $log_desc = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

echo SmallRNACategory=`date` 

mono-sgen $cqsFile smallrna_group -i $filelist -o $resultDir

echo finished=`date`

exit 0 
";
  close(OUT);
  print SH "\$MYCMD ./$pbsName \n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all smallRNA category tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $countFile = "${resultDir}/${sampleName}.catcount";

    my @resultFiles = ();
    push( @resultFiles, $countFile );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
