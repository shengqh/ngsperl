#!/usr/bin/perl
package CQS::UniqueWrapper;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::UniqueTask;
use CQS::StringUtils;
use File::Spec;

our @ISA = qw(CQS::UniqueTask);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_uw";
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $output_file                = get_option( $config, $section, "output_file",                "" );
  my $output_file_ext            = get_option( $config, $section, "output_file_ext",            "" );
  my $output_to_result_directory = get_option( $config, $section, "output_to_result_directory", 0 );
  my $output_perSample_file                = get_option( $config, $section, "output_perSample_file",                "" );
  my $output_perSample_file_ext            = get_option( $config, $section, "output_perSample_file_ext",            "" );

  my @output_file_exts;
  if ( $output_file_ext =~ /;/ ) {
    @output_file_exts = split( ";", $output_file_ext );
  }else{
    @output_file_exts = ($output_file_ext );
  }
  my @output_perSample_file_exts;
  if ( $output_perSample_file_ext =~ /;/ ) {
    @output_perSample_file_exts = split( ";", $output_perSample_file_ext );
  }else{
    @output_perSample_file_exts = ($output_perSample_file_ext );
  }

  my $result       = {};
  my @result_files = ();

  if ( $output_file eq "parameterSampleFile1" or $output_file eq "parameterSampleFile2" or $output_file eq "parameterSampleFile3" ) {
    if ( has_raw_files( $config, $section, $output_file ) ) {
      my %temp = %{ get_raw_files( $config, $section, $output_file ) };
      foreach my $sample_name ( keys %temp ) {
#        print "SampleName: $sample_name\n";
        if ( ref( $temp{$sample_name} ) eq "HASH" ) {
          foreach my $output_file_ext_one (@output_file_exts) {
            push( @result_files, "${result_dir}/${sample_name}${output_file_ext_one}" );
          }
        }
        else {
          foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
            my $subSampleName = $output_to_result_directory ? $result_dir . "/" . basename($subSampleFile) : $subSampleFile;
            foreach my $output_file_ext_one (@output_file_exts) {
              push( @result_files, "${subSampleName}${output_file_ext_one}" );
            }
          }
        }
      }
    }
    else {
      die "output_file defined as " . $output_file . ", but " . $output_file . " not in configure\n";
    }
  }
  else {
    foreach my $output_file_ext_one (@output_file_exts) {
      push( @result_files, "${result_dir}/${task_name}${output_file}${output_file_ext_one}" );
    }
  }
  $result->{$task_name} = filter_array( \@result_files, $pattern );

#per sample result
  if ( $output_perSample_file eq "parameterSampleFile1" or $output_perSample_file eq "parameterSampleFile2" or $output_perSample_file eq "parameterSampleFile3" ) {
    if ( has_raw_files( $config, $section, $output_perSample_file ) ) {
      my %temp = %{ get_raw_files( $config, $section, $output_perSample_file ) };
      foreach my $sample_name ( keys %temp ) {
#        print "SampleName: $sample_name\n";
        my @result_perSample_files = ();
        if ( ref( $temp{$sample_name} ) eq "HASH" ) {
          foreach my $output_file_ext_one (@output_perSample_file_exts) {
            push( @result_perSample_files, "${result_dir}/${sample_name}${output_file_ext_one}" );
          }
        }
        else {
          foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
            my $subSampleName = $output_to_result_directory ? $result_dir . "/" . basename($subSampleFile) : $subSampleFile;
            foreach my $output_file_ext_one (@output_perSample_file_exts) {
              push( @result_perSample_files, "${subSampleName}${output_file_ext_one}" );
            }
          }
        }
        $result->{$sample_name} = filter_array( \@result_perSample_files, $pattern );
      }
    }
    else {
      die "output_file defined as " . $output_file . ", but " . $output_file . " not in configure\n";
    }
  }

  return $result;
}

1;
