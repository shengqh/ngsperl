#!/usr/bin/perl
package Proteomics::Format::PrepareShiftMgf;

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
  $self->{_suffix} = "_psm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $shift_dalton      = get_option( $config, $section, "shift_dalton", "7" );
  my $shift_software    = get_option( $config, $section, "shift_software" );
  my $shift_option_file = get_option( $config, $section, "shift_option_file" );

  my $python_script = dirname(__FILE__) . "/prepareShiftMgf.py";
  if ( !-e $python_script ) {
    die "File not found : " . $python_script;
  }

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

    for my $sampleFile (@sample_files) {
      my $sampleBasename = basename($sampleFile);
      my $result_prefix  = change_extension( $sampleBasename, ".shifted" . round($shift_dalton) . "daltons" );
      my $result_file    = $result_prefix . ".optimal.mgf";

      print $pbs "if [ ! -s $result_file ]; then
  python3 $python_script -i $sampleFile -o $result_prefix -d $shift_dalton -p $shift_software -c $shift_option_file
fi
";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all ", $self->{_name}, " tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $shift_dalton = get_option( $config, $section, "shift_dalton", "7" );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my @result_files = ();

    for my $sampleFile (@sample_files) {
      my $sampleBasename = basename($sampleFile);
      my $result_prefix = change_extension( $sampleBasename, ".shifted" . round($shift_dalton) . "daltons" );

      #push( @result_files, "$result_dir/${result_prefix}.original.mgf" );
      #push( @result_files, "$result_dir/${result_prefix}.center.mgf" );
      push( @result_files, "$result_dir/${result_prefix}.optimal.mgf" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
