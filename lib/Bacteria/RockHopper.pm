#!/usr/bin/perl
package Bacteria::RockHopper;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::PairTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::PairTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_rh";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $raw_files = get_raw_files( $config, $section );
  my $groups    = get_raw_files( $config, $section, "groups" );
  my $pairs     = get_raw_files( $config, $section, "pairs" );
  my $rockhopper_jar = get_param_file( $config->{$section}{rockhopper_jar}, "rockhopper_jar", 1 );
  my $genome_dir  = get_option( $config, $section, "genome_dir",  1 );
  my $java_option = get_option( $config, $section, "java_option", "" );

  my %tpgroups = ();
  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    my @gfiles  = ();
    foreach my $sample_name ( sort @samples ) {
      my @fastqFiles = @{ $raw_files->{$sample_name} };
      push( @gfiles, join( '%', @fastqFiles ) );
    }
    $tpgroups{$group_name} = join( ",", @gfiles );
  }

  my $mapfile = $self->get_file( $result_dir, $task_name, ".map" );
  open( my $map, ">$mapfile" ) or die "Cannot create $mapfile";
  print $map "Fastq\tSample\tGroup\n";
  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    foreach my $sample_name ( sort @samples ) {
      my @fastqFiles = @{ $raw_files->{$sample_name} };
      foreach my $fastq ( sort @fastqFiles ) {
        my $fname = basename($fastq);
        $fname = change_extension( $fname, "" );
        print $map "$fname\t$sample_name\t$group_name\n";
      }
    }
  }
  close $map;

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct), "\n";

  for my $pair_name ( sort keys %{$pairs} ) {
    my @group_names = @{ $pairs->{$pair_name} };
    if ( scalar(@group_names) != 2 ) {
      die $pair_name . " should include and only include two groups, currently is [" . join( ", ", @group_names );
    }
    my @fastqs = ();
    push( @fastqs, $tpgroups{ $group_names[1] } );
    push( @fastqs, $tpgroups{ $group_names[0] } );
    my $fastqstrs = join( " ", @fastqs );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $pair_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $pair_name );

    my $cur_dir = create_directory_or_die( $result_dir . "/$pair_name" );

    my $labels = $group_names[1] . "," . $group_names[0];

    my $log_desc   = $cluster->get_log_description($log);
    my $final_file = "summary.txt";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "
java $java_option -cp $rockhopper_jar Rockhopper $option -g $genome_dir -L $labels -o . $fastqstrs

if [ -e intermediary ]; then
  rm -rf intermediary
fi

";

    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all FastQC tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};

  my $raw_files = get_raw_files( $config, $section );
  my $groups    = get_raw_files( $config, $section, "groups" );
  my $pairs     = get_raw_files( $config, $section, "pairs" );
  my $genome_name = get_option( $config, $section, "genome_name", 1 );

  for my $pair_name ( sort keys %{$pairs} ) {
    my $cur_dir      = $result_dir . "/$pair_name";
    my @result_files = ();
    push( @result_files, $cur_dir . "/" . $genome_name . "_operons.txt" );
    push( @result_files, $cur_dir . "/" . $genome_name . "_transcripts.txt" );
    $result->{$pair_name} = filter_array( \@result_files, $pattern );
  }
  my $mapfile = $self->get_file( $result_dir, $task_name, ".map" );
  my @result_files = ();
  push( @result_files, $mapfile );
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
