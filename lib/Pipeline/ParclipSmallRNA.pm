#!/usr/bin/perl
package Pipeline::ParclipSmallRNA;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::SmallRNAUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(getParclipSmallRNAConfig performParclipSmallRNA performParclipSmallRNATask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getParclipSmallRNAConfig {
  my ($def) = @_;

  my ( $config, $individual_ref, $summary_ref, $cluster ) = getPrepareConfig($def);
  my @individual = @{$individual_ref};
  my @summary    = @{$summary_ref};

  my $gsnap = {
    gsnap => {
      class                 => 'Alignment::Gsnap',
      perform               => 1,
      target_dir            => $def->{target_dir} . '/gsnap',
      option                => '-y 0 -z 0 -Y 0 -Z 0 -m 1 -Q --trim-mismatch-score 0 --trim-indel-score 0 --mode ttoc-nonstranded --gunzip',
      gsnap_index_directory => $def->{gsnap_index_directory},
      gsnap_index_name      => $def->{gsnap_index_name},
      source_ref            => [ 'identical_NTA', '.fastq.gz\$' ],
      sh_direct             => 0,
      cluster               => $cluster,
      pbs                   => {
        'email'    => $def->{email},
        'nodes'    => '1:ppn=' . $def->{max_thread},
        'walltime' => '72',
        'mem'      => '40gb'
      }
    },
    gsnap_smallRNA_count => {
      class           => 'CQS::SmallRNACount',
      perform         => 1,
      target_dir      => $def->{target_dir} . '/gsnap_smallRNA_count',
      option          => '-s -e 4',
      source_ref      => 'gsnap',
      seqcount_ref    => [ 'identical_NTA', '.dupcount\$' ],
      coordinate_file => $def->{coordinate},
      fasta_file      => $def->{coordinate_fasta},
      cqs_tools       => $def->{cqstools},
      sh_direct       => 0,
      cluster         => $cluster,
      pbs             => {
        'email'    => $def->{email},
        'walltime' => '72',
        'mem'      => '40gb',
        'nodes'    => '1:ppn=1'
      },
    },
    gsnap_smallRNA_t2c_summary => {
      class      => 'SmallRNA::T2CSummary',
      perform    => 1,
      target_dir => $def->{target_dir} . '/gsnap_smallRNA_t2c',
      option     => '',
      source_ref => [ 'gsnap_smallRNA_count', '.mapped.xml$' ],
      cqs_tools  => $def->{cqstools},
      sh_direct  => 0,
      cluster    => $cluster,
      pbs        => {
        'email'    => $def->{email},
        'walltime' => '72',
        'mem'      => '40gb',
        'nodes'    => '1:ppn=1'
      },
    },
  };

  push @individual, ( 'gsnap', 'gsnap_smallRNA_count' );
  push @summary, ('gsnap_smallRNA_t2c_summary');

  $config = merge( $config, $gsnap );
  $config->{sequencetask} = {
    class      => 'CQS::SequenceTask',
    perform    => 1,
    target_dir => $def->{target_dir} . '/sequencetask',
    option     => '',
    source     => {
      step1 => \@individual,
      step2 => \@summary,
    },
    sh_direct => 0,
    cluster   => 'slurm',
    pbs       => {
      'email'    => $def->{email},
      'nodes'    => '1:ppn=' . $def->{max_thread},
      'walltime' => '72',
      'mem'      => '40gb'
    },
  };

  return ($config);
}

sub performParclipSmallRNA {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getParclipSmallRNAConfig($def);

  if ($perform) {

    my $configFile = $def->{target_dir} . '/' . $def->{task_name} . '.config';
    open( SH, '>$configFile' ) or die 'Cannot create $configFile';
    print SH Dumper($config);
    close(SH);

    performConfig($config);
  }

  return $config;
}

sub performParclipSmallRNATask {
  my ( $def, $task ) = @_;

  my $config = getParclipSmallRNAConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
