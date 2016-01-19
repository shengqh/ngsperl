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
  $self->{_name}   = __PACKAGE__;
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
  open( my $out, ">$result" ) or die $!;

  foreach my $k ( sort keys %{$cufflinks_gtf} ) {
    my @gtfs = @{ $cufflinks_gtf->{$k} };

    #print "$gtfs[0]\n";
    print $out "$gtfs[0]\n";
  }
  close $out;

  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );
  my $gtfparam = "";
  if ( defined $transcript_gtf ) {
    $gtfparam = "-g $transcript_gtf";
  }

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $assembliesfile = get_assemblies_file( $config, $section, $result_dir );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $final_file = "merged.gtf";
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  print $pbs "cuffmerge $option $gtfparam -s $faFile -o . $assembliesfile";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = { $task_name => [ $result_dir . "/merged.gtf" ] };

  return $result;
}

1;
