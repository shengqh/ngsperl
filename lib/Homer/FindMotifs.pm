#!/usr/bin/perl
package Homer::FindMotifs;

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
  $self->{_suffix} = "_fm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my %tagDirectories = %{ get_raw_files( $config, $section ) };

  my $pairs = get_raw_files( $config, $section );
  my $fasta = get_param_file( $config, $section, "genome_fasta", 1 );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $pair_name ( sort keys %{$pairs} ) {
    my $pairFile = $pairs->{$pair_name};

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $pair_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $pair_name );

    my $final_file = $pair_name . ".motif";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    print $pbs "findMotifsGenome.pl $pairFile $fasta peakAnalysis -size 200 -len 8";
    $self->close_pbs( $pbs, $pbs_file );
  
    print $sh "\$MYCMD ./$pbs_name \n";
  }
  print $sh "exit 0\n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $pairs = get_raw_files( $config, $section );

  my $result = {};
  for my $pair_name ( sort keys %{$pairs} ) {
    my @result_files = ();
    push( @result_files, "${result_dir}/${pair_name}.motif" );
    $result->{$pair_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
