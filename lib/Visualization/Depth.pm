#!/usr/bin/perl
package Visualization::Depth;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::GroupTask;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_vd";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $bedFiles  = get_raw_files( $config, $section );
  my $groups    = get_raw_files( $config, $section, "groups" );
  my $bam_files = get_raw_files( $config, $section, "bam_files" );
  my $singlepdf   = get_option( $config, $section, "single_pdf",   0 ) ? "-s" : "";
  my $facetSample = get_option( $config, $section, "facet_sample", 1 ) ? "-f" : "";
  my $drawLine    = get_option( $config, $section, "draw_line",    0 ) ? "-l" : "";

  my $cnvr_files;
  if ( has_raw_files( $config, $section, "cnvr_files" ) ) {
    $cnvr_files = get_raw_files( $config, $section, "cnvr_files" );
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my $perl = dirname(__FILE__) . "/Depth.pl";

  for my $name ( sort keys %{$bedFiles} ) {
    my $cur_dir = create_directory_or_die( $result_dir . "/$name" );

    my @curBedFiles  = @{ $bedFiles->{$name} };
    my @curBamNames  = @{ $groups->{$name} };
    my @curbam_files = ();
    for my $bamName (@curBamNames) {
      push( @curbam_files, $bam_files->{$bamName}->[0] );
    }

    my $cnvr_option = "";
    if ( defined $cnvr_files ) {
      $cnvr_option = "--cnvrFile " . $cnvr_files->{$name}->[0];
    }
    my $curbam_fileStr = join( ' ', @curbam_files );

    my $bamCount = scalar(@curBamNames);

    my $configFileName = "${name}.filelist";
    my $configFile     = $pbs_dir . "/${configFileName}";
    open( my $con, ">$configFile" ) or die "Cannot create $configFile";
    for ( my $index = 0 ; $index < $bamCount ; $index++ ) {
      print $con $curBamNames[$index], "\t", $curbam_files[$index], "\n";
    }
    close $con;

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $name );

    print $sh "\$MYCMD ./$pbs_name \n";
    my $log_desc = $cluster->get_log_description($log);

    my $depthFile = $cur_dir . "/" . basename( $curBedFiles[0] ) . ".depth";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $depthFile );
    for my $bedFile (@curBedFiles) {
      print $pbs "perl $perl -b $bedFile -c $configFile $cnvr_option $singlepdf $facetSample $drawLine \n";
    }
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

  my $bedFiles = get_raw_files( $config, $section );

  my $result = {};
  for my $name ( sort keys %{$bedFiles} ) {
    my $cur_dir      = create_directory_or_die( $result_dir . "/$name" );
    my @result_files = ();
    push( @result_files, $pbs_dir . "/" . $name . ".filelist" );

    my @curBedFiles = @{ $bedFiles->{$name} };
    for my $bedFile (@curBedFiles) {
      push( @result_files, $cur_dir . "/" . basename($bedFile) . ".depth" );
      push( @result_files, $cur_dir . "/" . basename($bedFile) . ".reads" );
    }

    $result->{$name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
