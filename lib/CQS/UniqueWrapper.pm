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
  $self->{_forbid_tmp_folder} = 1;
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern,$removeEmpty ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $output_no_name             = get_option( $config, $section, "output_no_name",                0 );
  my $output_taskname = $output_no_name ? "" : $task_name;

  my $output_file                = get_option( $config, $section, "output_file",                "" );
  my $output_file_exts = get_output_ext_list( $config, $section );

  #print($section .  Dumper($output_file_exts));

  my $output_to_result_directory = get_option( $config, $section, "output_to_result_directory", 0 );
  my $output_perSample_file                = get_option( $config, $section, "output_perSample_file",                "" );
  my $output_perSample_file_byName         = get_option( $config, $section, "output_perSample_file_byName",                0 );
  my $output_perSample_file_ext            = get_option( $config, $section, "output_perSample_file_ext",            "" );

  my @output_perSample_file_exts;
  if ( $output_perSample_file_ext =~ /;/ ) {
    @output_perSample_file_exts = split( ";", $output_perSample_file_ext );
  }else{
    @output_perSample_file_exts = ($output_perSample_file_ext );
  }
  @output_perSample_file_exts = grep { $_ ne '' } @output_perSample_file_exts; #remove empty elements

  my $result       = {};

  if ( $output_file eq "parameterSampleFile1" or $output_file eq "parameterSampleFile2" or $output_file eq "parameterSampleFile3" ) {
    if ( has_raw_files( $config, $section, $output_file ) ) {
      my @result_files = ();
      my %temp = %{ get_raw_files( $config, $section, $output_file ) };
      foreach my $sample_name ( keys %temp ) {
#        print "SampleName: $sample_name\n";
        if ( is_hash( $temp{$sample_name} ) ) {
          foreach my $output_file_ext_one (@$output_file_exts) {
            push( @result_files, "${result_dir}/${sample_name}${output_file_ext_one}" );
          }
        }
        else {
          foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
            my $subSampleName = $output_to_result_directory ? $result_dir . "/" . basename($subSampleFile) : $subSampleFile;
            foreach my $output_file_ext_one (@$output_file_exts) {
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
    my @result_files = ();
    foreach my $output_file_ext_one (@$output_file_exts) {
      if (($output_file ne "") or ($output_file_ext_one ne "")) {
        push( @result_files, "${result_dir}/${output_taskname}${output_file}${output_file_ext_one}" );
      }
    }
    
    if(scalar(@result_files) > 0){
      my $filtered = filter_array( \@result_files, $pattern, 1 );
      if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
        $result->{$task_name} = $filtered;
      }
    }
  }

  if ( $output_perSample_file eq "parameterSampleFile1" or $output_perSample_file eq "parameterSampleFile2" or $output_perSample_file eq "parameterSampleFile3" ) {
    #per sample result
    if ( has_raw_files( $config, $section, $output_perSample_file ) ) {
      my %temp = %{ get_raw_files( $config, $section, $output_perSample_file ) };
      foreach my $sample_name ( keys %temp ) {
#        print "SampleName: $sample_name\n";
        my @result_perSample_files = ();
        if ( is_hash( $temp{$sample_name}) or $output_perSample_file_byName ) {
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
        #$result->{$sample_name} = filter_array( \@result_perSample_files, $pattern );
        my $filtered = filter_array( \@result_perSample_files, $pattern, $removeEmpty );
        if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
          $result->{$sample_name} = $filtered;
        }
      }
    }
    else {
      die "output_file defined as " . $output_file . ", but " . $output_file . " not in configure\n";
    }
  }

  #print(Dumper($result));

  return $result;
}

1;
