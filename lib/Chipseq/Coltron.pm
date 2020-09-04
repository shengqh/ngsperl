#!/usr/bin/perl
package Chipseq::Coltron;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cl";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $genome = get_option( $config, $section, "genome" );
  my %treatments_files = %{ get_grouped_raw_files( $config, $section, "groups" ) };
  my $pipeline_dir = get_directory( $config, $section, "pipeline_dir", 1 );

  #print Dumper(%treatments_files);

  my %enhancer_files = %{ get_raw_files( $config, $section, "enhancer_files" ) };

  #print Dumper(%enhancer_files);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %enhancer_files ) {
    my @e_files  = @{ $enhancer_files{$sample_name} };
    my $enhancer = "-e " . $e_files[0];

    my @b_files = @{ $treatments_files{$sample_name} };
    my $bam     = "-b " . $b_files[0];

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name/" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = "${cur_dir}/${sample_name}_CLIQUES_RANKED.txt";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );
    print $pbs "cd $pipeline_dir
coltron $enhancer $bam -g $genome -o $cur_dir -n $sample_name $option
";
    $self->close_pbs( $pbs, $pbs_file );
    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %enhancer_files = %{ get_raw_files( $config, $section, "enhancer_files" ) };

  my $result = {};
  for my $sample_name ( sort keys %enhancer_files ) {
    my $cur_dir      = $result_dir . "/$sample_name/";
    my @result_files = ();
    push( @result_files, "${cur_dir}/${sample_name}_CANIDATE_TF_AND_SUPER_TABLE.txt" );
    push( @result_files, "${cur_dir}/${sample_name}_CLIQUES_ALL.txt" );
    push( @result_files, "${cur_dir}/${sample_name}_CLIQUES_RANKED.txt" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
