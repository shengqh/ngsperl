#!/usr/bin/perl
package Cufflinks::Cuffmerge;

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
  $self->{_name}   = "Cuffmerge";
  $self->{_suffix} = "_cmerge";
  bless $self, $class;
  return $self;
}

sub get_assemblies_file {
  my ( $config, $section, $target_dir ) = @_;

  my $result = get_param_file( $config->{$section}{source}, "${section}::source", 0 );

  if ( defined $result ) {
    return $result;
  }

  my $cufflinks_gtf = get_raw_files( $config, $section, "source", ".gtf\$" );

  #print_hash($cufflinks_gtf);

  $result = $target_dir . "/assemblies.txt";
  open( OUT, ">$result" ) or die $!;

  foreach my $k ( sort keys %{$cufflinks_gtf} ) {
    my @gtfs = @{ $cufflinks_gtf->{$k} };

    #print "$gtfs[0]\n";
    print OUT "$gtfs[0]\n";
  }
  close OUT;

  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );
  my $gtfparam = "";
  if ( defined $transcript_gtf ) {
    $gtfparam = "-g $transcript_gtf";
  }

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $assembliesfile = get_assemblies_file( $config, $section, $resultDir );

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  my $cluster = get_cluster( $config, $section );
  my $log_desc = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

echo cuffmerge=`date` 

cuffmerge $option $gtfparam -s $faFile -o . $assembliesfile 

echo finished=`date`

exit 0
";

  close(OUT);

  print "$pbsFile created\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = { $task_name => [ $resultDir . "/merged.gtf" ] };

  return $result;
}

1;
