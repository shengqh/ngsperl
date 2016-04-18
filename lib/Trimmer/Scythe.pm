#!/usr/bin/perl
package Trimmer::Scythe;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_scy";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $selfname = $self->{_name};

  my $faFile = get_param_file( $config->{$section}{adapter_file}, "adapter_file", 1 );
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

    if ( scalar(@sample_files) == 2 ) {
      my $sample1 = $sample_files[0];
      my $sample2 = $sample_files[1];

      my $trim1 = change_extension_gzipped( basename($sample1), "_scythe.fastq" );
      my $trim2 = change_extension_gzipped( basename($sample2), "_scythe.fastq" );

      my $final_file1 = $trim1 . ".gz";
      my $final_file2 = $trim2 . ".gz";

      print $pbs "
if [ ! -s $final_file1 ]; then
  if [ ! -s $trim1 ]; then
    scythe $option -a $faFile -o $trim1 $sample1
  fi

  if [ ! -s $trim2 ]; then
    scythe $option -a $faFile -o $trim2 $sample2
  fi

  gzip $trim1
  gzip $trim2
fi
";
    }
    else {
      my $sample1 = $sample_files[0];

      my $trim1 = change_extension_gzipped( basename($sample1), "_scythe.fastq" );

      my $final_file1 = $trim1 . ".gz";

      print $pbs "
if [ ! -s $final_file1 ]; then
  if [ ! -s $trim1 ]; then
    scythe $option -a $faFile -o $trim1 $sample1
  fi

  gzip $trim1
fi
";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my @result_files = ();

    if ( scalar(@sample_files) == 2 ) {
      my $sample1 = $sample_files[0];
      my $sample2 = $sample_files[1];

      my $trim1 = change_extension_gzipped( basename($sample1), "_scythe.fastq" );
      my $trim2 = change_extension_gzipped( basename($sample2), "_scythe.fastq" );

      my $final_file1 = $trim1 . ".gz";
      my $final_file2 = $trim2 . ".gz";

      push( @result_files, "${result_dir}/${final_file1}" );
      push( @result_files, "${result_dir}/${final_file2}" );
    }
    else {
      my $sample1 = $sample_files[0];

      my $trim1 = change_extension_gzipped( basename($sample1), "_scythe.fastq" );

      my $final_file1 = $trim1 . ".gz";

      push( @result_files, "${result_dir}/${final_file1}" );
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
