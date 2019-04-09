#!/usr/bin/perl
package LSAM::Lsam;

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
  $self->{_suffix} = "_lsam";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my ( $raw_files, $samples ) = get_raw_files_and_keys( $config, $section );

  my $lsamSoftware = get_option( $config, $section, "liver_model_console" );
  my $randomSeed   = get_option( $config, $section, "random_seed" );
  my $iteration    = get_option( $config, $section, "iteration" );

  my $shfile = $self->get_file( $pbs_dir, $task_name, ".bat" );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";

  for my $name (@$samples) {
    my @sample_files = @{ $raw_files->{$name} };
    my $sample       = $sample_files[0];

    my $pbs_file = $self->get_file( $pbs_dir, $name, ".bat" );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $name );

    my $log_desc = $cluster->get_log_description($log);

    my $sampleFolder = create_directory_or_die( $result_dir . "/" . $name );

    my $pbs = $self->open_pbs( $pbs_file, "", "", $path_file, $sampleFolder );
    my $curRandomSeed = $randomSeed;

    print $pbs "
\@echo off

set minbytesize=200

";
    for my $iter ( 1 .. $iteration ) {
      my $curName = sprintf( "${name}_iter%02d_", $iter );
      my $curFinalFile = $curName . "Summary.out";
      my $command = "$lsamSoftware $option -i:$sample -o:$sampleFolder -n:$curName -s:$curRandomSeed";
      
      print $pbs "
echo Checking $curFinalFile ...
set file=$curFinalFile
if EXIST \%file\% (
  echo \%file\% exists, check file size ...
  set curSize=0
  FOR /F \"usebackq\" \%\%A IN ('\%file\%') DO (
    echo filesize=\%\%~zA
    if \%\%~zA LSS \%minbytesize\% (
      echo failed, do simulation, start at %date% %time%
      $command
    ) ELSE (
      echo simulation has been done, ignored.
    )
  )
) ELSE (
  echo no result, do simulation, start at %date% %time%
  $command
)
";
      $curRandomSeed = $curRandomSeed + 1;
    }
    $self->close_pbs( $pbs, $pbs_file );

    open( $pbs, "<", $pbs_file ) or die "Could not open file $pbs_file: $!";
    while ( my $line = <$pbs> ) {
      if ( $line !~ /exit/ ) {
        print $sh $line;
      }
    }
  }
  close($sh);
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $iteration    = get_option( $config, $section, "iteration" );

  my $result = {};

  for my $name ( sort keys %raw_files ) {
    my @result_files = ();

    my $sampleFolder = create_directory_or_die( $result_dir . "/" . $name );
    for my $iter ( 1 .. $iteration ) {
      my $curName = sprintf( "${name}_iter%02d_", $iter );
      push( @result_files, $sampleFolder . "/" . $curName . "Death.out" );
    }

    $result->{$name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}
1;
