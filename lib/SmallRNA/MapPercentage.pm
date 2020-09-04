#!/usr/bin/perl
package SmallRNA::MapPercentage;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

my $directory;

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

  my $rawFiles = get_raw_files( $config, $section );

  my $rtemplate = dirname(__FILE__) . "/mapPercentage.R";
  if ( !-e $rtemplate ) {
    die "File not found : " . $rtemplate;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $fileName ( keys %$rawFiles ) {
    my $fileListFile = "${result_dir}/${fileName}.filelist";
    open( my $df, ">$fileListFile" ) or die "Cannot create $fileListFile";
    my $files = $rawFiles->{$fileName};
    for my $file (@$files) {
      print $df "$file\n";
    }
    close($df);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $fileName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $fileName );
    my $log_desc = $cluster->get_log_description($log);

    my $final_file = $fileName . ".png";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $pbs "R --vanilla -f $rtemplate --args $fileListFile $fileName\n";

    $self->close_pbs( $pbs, $pbs_file );
    print $sh "\$MYCMD ./$pbs_name \n";
  }
  print $sh "exit 0\n";
  close $sh;
  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $rawFiles = get_raw_files( $config, $section );
  my $result = {};
  for my $fileName ( keys %$rawFiles ) {
    my @result_files = ();
    push( @result_files, $result_dir . "/${fileName}.png" );
    push( @result_files, $result_dir . "/${fileName}.csv" );
    $result->{$fileName} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
