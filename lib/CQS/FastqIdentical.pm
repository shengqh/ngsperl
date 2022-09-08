#!/usr/bin/perl
package CQS::FastqIdentical;

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
  $self->{_suffix} = "_IQB";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $extension = get_option( $config, $section, "extension" );
  my $use_first_read_after_trim = get_option( $config, $section, "use_first_read_after_trim");

  my $minlen = $config->{$section}{minlen};
  if ( defined $minlen ) {
    $minlen = "-l $minlen";
  }
  else {
    $minlen = "";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $final_file   = $sample_name . $extension;

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    if ( scalar(@sample_files) == 1 or $use_first_read_after_trim ) {
      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
      print $pbs "cqstools fastq_identical $option -i $sample_files[0] $minlen -o $final_file \n";
      $self->close_pbs( $pbs, $pbs_file );
    }
    else {
      my $firstOutputFile = change_extension_gzipped( basename($sample_files[0]), $extension );
      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $firstOutputFile );
      for my $sampleFile (@sample_files) {
        my $outputFile = change_extension_gzipped( basename($sampleFile), $extension );
        print $pbs "cqstools fastq_identical -i $sampleFile $minlen -o $outputFile \n";
      }
      $self->close_pbs( $pbs, $pbs_file );
    }
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $extension = get_option( $config, $section, "extension" );
  my $merge_result = get_option( $config, $section, "merge_result", 0 );
  my $use_first_read_after_trim = get_option( $config, $section, "use_first_read_after_trim");

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $final_file   = $result_dir . "/" . $sample_name . $extension;
    my @result_files = ();

    if ( scalar(@sample_files) == 1 or $use_first_read_after_trim) {
      push( @result_files, $final_file );
      push( @result_files, change_extension( $final_file, ".dupcount" ) );
    }
    else {
      for my $sampleFile (@sample_files) {
        my $outputFile = $result_dir . "/" . change_extension_gzipped( basename($sampleFile), $extension );
        push( @result_files, $outputFile );
        push( @result_files, change_extension( $outputFile, ".dupcount" ) );
      }
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
