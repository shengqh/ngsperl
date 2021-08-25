#!/usr/bin/perl
package Visualization::Heatmap;

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
  $self->{_suffix} = "_hm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $bamFiles = get_raw_files( $config, $section );
  my $groups   = get_raw_files( $config, $section, "groups" );
  my $bedFiles = get_raw_files( $config, $section, "bed_files" );
  my $bamliquidator = get_option( $config, $section, "bamliquidator" );
  my $window        = get_option( $config, $section, "window", 5000 );
  my $extension     = get_option( $config, $section, "extension" );
  my $color         = get_option( $config, $section, "color", "blue" );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my $python_script = dirname(__FILE__) . "/heatmap.py";
  my $r_script      = dirname(__FILE__) . "/heatmap.r";

  for my $name ( sort keys %$groups ) {
    my $curBedFile  = $bedFiles->{$name}->[0];
    my @curBamNames = @{ $groups->{$name} };

    my $outputPrefix = $name;
    my $finalFile    = $outputPrefix . ".unscaled_heat.png";

    my $bamListFile = "$name.bam.filelist";
    open( my $list, ">$result_dir/$bamListFile" ) or die "Cannot create $result_dir/$bamListFile";
    my $lastGff = "";
    foreach my $sample_name (@curBamNames) {
      print $list $bamFiles->{$sample_name}->[0] . "\t$sample_name\n";
      $lastGff = $name . "." . $sample_name . ".gff";
    }
    close($list);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $gffListFile = "$name.gff.filelist";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $finalFile );
    print $pbs "if [ ! -s $lastGff ]; then
  python3 $python_script --bamliquidator $bamliquidator -i $curBedFile -d $bamListFile -o $outputPrefix -w $window -e $extension
fi

if [ -s $lastGff ]; then
  R --vanilla -f $r_script --args $gffListFile $color $name TRUE $result_dir
fi
";
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
  for my $name ( sort keys %$groups ) {
    my @result_files = ();
    push( @result_files, $result_dir . "/" . $name . ".unscaled_heat.png" );
    $result->{$name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
