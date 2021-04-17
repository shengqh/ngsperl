#!/usr/bin/perl
package Trimmer::RemoveFastqTerminalN;

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
  $self->{_suffix} = "_ft";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $extension = get_option( $config, $section, "extension" );
  
  my $py_script = dirname(__FILE__) . "/removeFastqTerminalN.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }
  
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );
    if ( scalar(@sample_files) == 1 ) {
      my $sampleFile = $sample_files[0];
      my $trimFile = get_trim_file( $sampleFile, $extension );
      print $pbs "if [ ! -s $trimFile ]; then
  python3 $py_script $option -i $sampleFile -o $trimFile 
fi
";
    }
    else {
      my $read1file = $sample_files[0];
      my $read2file = $sample_files[1];
      my $trim1file = get_trim_file( $read1file, $extension );
      my $trim2file = get_trim_file( $read2file, $extension );
      print $pbs "if [ ! -s $trim1file ]; then
  python3 $py_script $option -i $read1file,$read2file -o $trim1file,$trim2file
fi
";
    }

    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->name() . " tasks.\n";

  #`qsub $pbs_file`;
}

sub get_trim_file {
  my ( $sampleFile, $extension ) = @_;
  my $fileName = basename($sampleFile);
  if ( $fileName =~ /.gz$/ ) {
    $fileName = change_extension( $fileName, "" );
  }
  $fileName = change_extension( $fileName, "" );
  my $result = $fileName . $extension;
  return ($result);
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0);

  my $extension = get_option( $config, $section, "extension" );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my @result_files = ();

    if ( scalar(@sample_files) == 1 ) {
      my $trimFile = get_trim_file( $sample_files[0], $extension );
      push( @result_files, "${result_dir}/${trimFile}" );
    }
    else {
      my $trim1file = get_trim_file( $sample_files[0], $extension );
      my $trim2file = get_trim_file( $sample_files[1], $extension );
      push( @result_files, "${result_dir}/${trim1file}" );
      push( @result_files, "${result_dir}/${trim2file}" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
