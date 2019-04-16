#!/usr/bin/perl
package Blast::FastaScatter;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;
use Bio::SeqIO;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_fs";
  bless $self, $class;
  return $self;
}

sub getSeq{
  my ( $self, $config, $section ) = @_;
  my $fasta        = get_option_file( $config, $section, "source" );
  my $seqio = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
  my $seqs = [];
  while(my $seq = $seqio->next_seq) {
    push @$seqs, $seq;
  }
  return $seqs;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $seqs = $self->getSeq();
  
  for my $seq (@$seqs) {
    my $sample_name = $seq->accession_number;
    my $sample_fasta = $result_dir . "/" . $sample_name . ".fa";
    my $seq_out = Bio::SeqIO->new( -file   => ">$sample_fasta", -format => "Fasta" );
    $seq_out->write_seq($seq);
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $seqs = $self->getSeq();

  my $result = {};
  for my $seq (@$seqs) {
    my $sample_name = $seq->accession_number;
    my $sample_fasta = $result_dir . "/" . $sample_name . ".fa";
    $result->{$sample_name} = filter_array( [$sample_fasta], $pattern );
  }

  return $result;
}
1;
