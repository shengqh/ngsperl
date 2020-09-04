#!/usr/bin/perl
package Samtools::Mpileup;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my %group_sample_map = %{ get_group_sample_map( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $group_name ( sort keys %group_sample_map ) {
    my @sample_files = @{ $group_sample_map{$group_name} };
    my $sampleCount  = scalar(@sample_files);
    my $samples      = "";
    for ( my $index = 0 ; $index < $sampleCount ; $index++ ) {
      $samples = $samples . " " . $sample_files[$index][1];
    }

    my $cur_dir = create_directory_or_die( $result_dir . "/$group_name" );

    my $mpileup = "${group_name}.mpileup";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $mpileup );
    print $pbs "samtools mpileup $option -f $faFile $samples > $mpileup";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $group_name ( keys %{$groups} ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$group_name";
    my $mpileup      = "${group_name}.mpileup";
    push( @result_files, "$cur_dir/$mpileup" );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
