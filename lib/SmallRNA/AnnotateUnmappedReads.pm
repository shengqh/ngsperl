#!/usr/bin/perl
package SmallRNA::AnnotateUnmappedReads;

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
  $self->{_suffix} = "_fm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $raw_files = get_raw_files( $config, $section );
  my $mapped_files = get_raw_files( $config, $section, "mapped_files" );
  my $mapped_names = get_option( $config, $section, "mapped_names" );
  my $min_count = get_option($config, $section, "min_count");
  my $extension = get_option( $config, $section, "extension", ".unmapped.tsv" );

  my $py_script = dirname(__FILE__) . "/AnnotateUnmappedReads.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( keys %$raw_files ) {
    my @sample_files = @{ $raw_files->{$sample_name} };
    my $samples = join(',', @sample_files);
    
    my @mapped_files = @{ $mapped_files->{$sample_name} };
    my $mappeds = join(',', @mapped_files);

    my $final_file   = $sample_name . $extension;

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    print $pbs "python3 $py_script -i $samples -b $mappeds -n $mapped_names -m $min_count -o $final_file $option";
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $extension = get_option( $config, $section, "extension", ".unmapped.tsv" );

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $final_file = $result_dir . "/" . $sample_name . $extension;

    my @result_files = ();
    push( @result_files, $final_file );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
