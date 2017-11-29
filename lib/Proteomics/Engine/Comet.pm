#!/usr/bin/perl
package Proteomics::Engine::Comet;

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
  $self->{_suffix} = "_comet";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $param_file = get_param_file( $config->{$section}{param_file}, "param_file", 1 );
  my $database   = get_param_file( $config->{$section}{database},   "database",   1 );

  open my $handle, '<', $param_file;
  chomp( my @lines = <$handle> );
  close $handle;

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
      my $sname = basename($sampleFile);
      my $result_prefix = change_extension( $sname, "" );
      my $result_file = $result_prefix . ".pep.xml";

      print $pbs "if [ ! -s $result_file ]; then\n";

      print $pbs "  comet -P$param_file -D$database -N$result_prefix $sampleFile
  if [ -s $result_file ]; then
    RefreshParser $result_file $database
  fi
";

      print $pbs "
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $param_file = get_param_file( $config->{$section}{param_file}, "param_file", 1 );
  open my $handle, '<', $param_file;
  chomp( my @lines = <$handle> );
  close $handle;
  my $outputPercolatorPin = grep( /output_percolatorfile\s*=\s*1/, @lines );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my @result_files = ();

    for my $sampleFile (@sample_files) {
      my $sname = basename($sampleFile);
      my $result_file = change_extension( $sname, ".pep.xml" );
      push( @result_files, "${result_dir}/${result_file}" );
      if ($outputPercolatorPin) {
        push( @result_files, "${result_dir}/" . change_extension( $sname, ".pin" ) );
      }
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
