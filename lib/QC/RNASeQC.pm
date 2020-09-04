#!/usr/bin/perl
package QC::RNASeQC;

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
  $self->{_suffix} = "_qc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $faFile         = get_param_file( $config->{$section}{fasta_file},     "fasta_file",     1 );
  my $jar            = get_param_file( $config->{$section}{jar},            "jar",            1, not $self->using_docker() );
  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 1 );
  my $sorted = get_option_value( $config->{$section}{sorted}, 1 );

  my $rrna_fasta = get_param_file( $config->{$section}{rrna_fasta}, "rrna_fasta", 0 );
  if ($rrna_fasta) {
    $option = $option . " -BWArRNA " . $rrna_fasta;
  }

  my $raw_files = get_raw_files( $config, $section );

  my $mapfile = $result_dir . "/${task_name}_sample.list";
  open( MAP, ">$mapfile" ) or die "Cannot create $mapfile";
  print MAP "SampleID\tbam_file\tNotes\n";
  for my $sample_name ( sort keys %{$raw_files} ) {
    my @bam_files = @{ $raw_files->{$sample_name} };
    for my $bam (@bam_files) {
      print MAP $sample_name, "\t", $bam, "\t", $sample_name, "\n";
    }
  }
  close(MAP);

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $final_file = "metrics.tsv";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );
  print $pbs "
java -jar $jar $option -s $mapfile -t $transcript_gtf -ttype 2 -r $faFile -o .

rm refGene.txt*
rm exons.rpkm.gct
rm */*.tmp.txt*
rm */*/perBaseDoC.out

";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  push( @result_files, $result_dir . "/metrics.tsv" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
