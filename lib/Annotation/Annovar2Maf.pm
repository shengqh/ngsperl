#!/usr/bin/perl
package Annotation::Annovar2Maf;

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
  $self->{_suffix} = "_a2m";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $combine_result = get_option( $config, $section, "combine_result", 0 );

  my $r_script = dirname(__FILE__) . "/annovar2maf.r";
  if ( !-e $r_script ) {
    die "File not found : " . $r_script;
  }

  my $py_script = dirname(__FILE__) . "/annovar2maf.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }

  my $refBuild = get_option( $config, $section, "refBuild" );

  my $sampleCount = scalar keys %raw_files;

  my $shfile;
  my $sh;
  if ( $sampleCount > 1 ) {
    $shfile = $self->get_task_filename( $pbs_dir, $task_name );
    open( $sh, ">$shfile" ) or die "Cannot create $shfile";
    print $sh get_run_command($sh_direct);
  }

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $vcf_file     = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    if ( $sampleCount > 1 ) {
      print $sh "\$MYCMD ./$pbs_name \n";
    }

    my $final_file = get_final_file($config, $section, $result_dir, $sample_name);
    my $log_desc   = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    my @maf_files = ();
    for my $sample_file (@sample_files) {
      my $cur_maf_file   = basename($sample_file) . ".tmp";
      my $cur_final_file = basename($sample_file) . ".maf";
      print $pbs "
R --vanilla -f $r_script --args $sample_file $cur_maf_file $refBuild
  
if [[ -s $cur_maf_file ]]; then
  python $py_script -i $cur_maf_file -o $cur_final_file
fi

if [[ -s $cur_final_file ]]; then
  rm $cur_maf_file
fi

";
      push @maf_files, $cur_final_file;
    }

    $self->close_pbs( $pbs, $pbs_file );
  }

  if ( $sampleCount > 1 ) {
    close $sh;

    if ( is_linux() ) {
      chmod 0755, $shfile;
    }

    print " !!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks . \n ";
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my @result_files = ();
    for my $sample_file (@sample_files) {
      my $cur_final_file = basename($sample_file) . ".maf";
      push @result_files, $result_dir . "/" . $cur_final_file;
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
