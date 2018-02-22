#!/usr/bin/perl
package CQS::CountTableGroupCorrelation;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::UniqueR;
use CQS::StringUtils;
use File::Spec;

our @ISA = qw(CQS::UniqueR);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cr";
  bless $self, $class;
  return $self;
}

sub getRcode {
  my ( $self, $config, $section ) = @_;
  my $result = get_option( $config, $section, "rCode", "" );

  my $suffix = get_option( $config, $section, "suffix", "" );
  if ( $suffix ne "" ) {
    $result = $result . 'suffix="' . $suffix . '";';
  }

  my $output_to_result_dir = get_option( $config, $section, "output_to_result_dir", 0 );
  if ($output_to_result_dir) {
    $result = $result . 'outputDirectory=".";';
  }

  return $result;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $output_to_result_dir = get_option( $config, $section, "output_to_result_dir", 0 );
  my $output_file_ext      = get_option( $config, $section, "output_file_ext",      "" );
  my $suffix               = get_option( $config, $section, "suffix",               "" );
  my $result               = {};
  my @result_files         = ();

  my %temp = %{ get_raw_files( $config, $section, "parameterSampleFile1" ) };
  foreach my $sample_name ( keys %temp ) {
    foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
      my $prefix = $subSampleFile;
      if ($output_to_result_dir) {
        my $file = basename($subSampleFile);
        my $pdir = dirname($subSampleFile);
        while ( basename($pdir) eq "result" ) {
          $pdir = dirname($pdir);
        }
        $prefix = $result_dir . "/" . basename($pdir) . "." . $file;
      }
      
      $prefix = $prefix . $suffix;

      if ( $output_file_ext =~ /;/ ) {
        my @output_file_exts = split( ";", $output_file_ext );
        foreach my $output_file_ext_one (@output_file_exts) {
          push( @result_files, "${prefix}${output_file_ext_one}" );
        }
      }
      else {
        push( @result_files, "${prefix}${output_file_ext}" );
      }
    }
  }
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
