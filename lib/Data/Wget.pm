#!/usr/bin/perl
package Data::Wget;

use strict;
use warnings;
use File::Basename;
use List::Util;
use File::Slurp;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::Task;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_wget";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $batch = get_option( $config, $section, "batch", 20 );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @listFiles = @{ $raw_files{$sample_name} };
    my $listFile  = $listFiles[0];
    my @urls      = read_file( $listFile, chomp => 1 );

    my $urlCount = scalar(@urls);

    for ( my $i = 0 ; $i < $urlCount ; $i += $batch ) {
      my $iend = $i + $batch;
      if ( $iend >= $urlCount ) {
        $iend = $urlCount;
      }

      my $name     = "${sample_name}_${i}_${iend}";
      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $name );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_log_filename( $log_dir, $name );

      my $log_desc = $cluster->get_log_description($log);

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

      for ( my $j = $i ; $j < $iend ; $j++ ) {
        my $url      = $urls[$j];
        my $filename = basename($url);

        print $pbs "if [ ! -s $filename ]; then
  echo 'wget $url ...\n'
  wget $url
fi

";
      }

      $self->close_pbs( $pbs, $pbs_file );

      print $sh "\$MYCMD ./$pbs_name \n";
    }

    close $sh;

    if ( is_linux() ) {
      chmod 0755, $shfile;
    }

    print "!!!shell file $shfile created, you can run this shell file to submit all ", $self->{_name}, " tasks.\n";
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = {};

  my %raw_files = %{ get_raw_files( $config, $section ) };
  for my $sample_name ( sort keys %raw_files ) {
    my @listFiles = @{ $raw_files{$sample_name} };
    my $listFile  = $listFiles[0];
    my @urls      = read_file( $listFile, chomp => 1 );
    for my $url (@urls) {
      my $filename    = basename($url);
      my $result_file = $result_dir . "/" . $filename;

      my @result_files = ();
      push( @result_files, $result_dir . "/${filename}" );
      $result->{$filename} = filter_array( \@result_files, $pattern );
    }
  }

  return $result;
}

1;
