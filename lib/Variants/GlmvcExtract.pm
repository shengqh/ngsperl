#!/usr/bin/perl
package Variants::GlmvcExtract;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ge";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $glmvcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1, not $self->using_docker() );

  my $raw_files = get_raw_files( $config, $section );
  my $bam_files = get_raw_files( $config, $section, "bam_files" );

  #print Dumper($raw_files);
  #print Dumper($bam_files);

  my %group_sample_map = ();
  my %group_name_map   = ();

  my $fafile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $groups = get_raw_files( $config, $section, "groups" );
  for my $group_name ( sort keys %{$raw_files} ) {
    my @samples = @{ $groups->{$group_name} };
    my @gfiles  = ();
    my @names   = ();
    my $index   = 0;
    foreach my $sample_name (@samples) {
      my @samplebam_files = @{ $bam_files->{$sample_name} };
      push( @gfiles, $samplebam_files[0] );
      push( @names,  $sample_name );
    }
    $group_sample_map{$group_name} = \@gfiles;
    $group_name_map{$group_name}   = \@names;
  }

  #print Dumper(%group_name_map);
  #print Dumper(%group_sample_map);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $group_name ( sort keys %{$raw_files} ) {
    my $validateFile = $raw_files->{$group_name}[0];
    my @sample_files = @{ $group_sample_map{$group_name} };
    my @sample_names = @{ $group_name_map{$group_name} };

    my $samples = join( ',', @sample_files );
    my $names   = join( ',', @sample_names );

    my $cur_dir = create_directory_or_die( $result_dir . "/$group_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc   = $cluster->get_log_description($log);
    my $final_file = "${group_name}.tsv";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );
    print $pbs "mono $glmvcfile extract $option --bam_files $samples --bam_names $names -f $fafile -o ${cur_dir}/$final_file -v $validateFile";
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $raw_files = get_raw_files( $config, $section );

  my $result = {};
  for my $group_name ( keys %{$raw_files} ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$group_name";

    push( @result_files, "$cur_dir/${group_name}.tsv" );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
