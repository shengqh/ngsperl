#!/usr/bin/perl
package CQS::ParclipTarget;

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
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pmt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );
  my $refgeneFile = get_param_file( $config->{$section}{refgene_file}, "refgene_file (Refgene file downloaded from UCSC Genome Browser)", 1 );
  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file (genome fasta file)", 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my %target_files = %{ get_raw_files( $config, $section, "target" ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my $xmlFile    = $raw_files{$sample_name}->[0];
    my $final_file = $sample_name . ".target.tsv";

    my $targetXml = $target_files{$sample_name}->[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );
    print $pbs "mono $cqstools parclip_mirna_target $option -i $xmlFile -t $targetXml -r $refgeneFile -g $faFile -o $final_file";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $final_file = $sample_name . ".target.tsv";

    my @result_files = ();
    push( @result_files, $result_dir . "/" . $final_file );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
