#!/usr/bin/perl
package Alignment::STARIndex;

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
  $self->{_suffix} = "_si";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my %sjdbFiles = %{ get_raw_files( $config, $section ) };

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $final = $result_dir . "/" . $task_name . ".tab";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final );

  for my $sample_name ( sort keys %sjdbFiles ) {
    my @sjdbs = @{ $sjdbFiles{$sample_name} };
    for my $sjdb (@sjdbs) {
      print $pbs "awk 'BEGIN {OFS=\"\\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' $sjdb >> $final \n";
    }
  }

  print $pbs "STAR $option --runThreadN $thread --runMode genomeGenerate --genomeDir . --genomeFastaFiles $faFile --sjdbGTFfile $transcript_gtf --sjdbFileChrStartEnd $final
  ";

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result = { $task_name => [$result_dir] };

  return $result;
}

1;
