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
    my $output_include_folder_name = get_option( $config, $section, "output_include_folder_name", 1 );
    $result = $result . "output_include_folder_name=" . ($output_include_folder_name?"TRUE":"FALSE") . ";";
  }

  return $result;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $output_to_result_dir = get_option( $config, $section, "output_to_result_dir", 0 );
  my $output_file_exts     = get_output_ext_list( $config, $section );
  my $output_file_task_ext = get_option( $config, $section, "output_file_task_ext", "" );
  my $suffix               = get_option( $config, $section, "suffix",               "" );
  my $result               = {};
  my $output_include_folder_name = get_option( $config, $section, "output_include_folder_name", 1 );
  #print("output_include_folder_name =". $output_include_folder_name."\n");

  my %temp = %{ get_raw_files( $config, $section, "parameterSampleFile1" ) };
  my @names = keys %temp;
  if ( exists( $temp{$task_name} ) ) {
    @names = grep { $_ ne $task_name } @names;
    push( @names, $task_name );
  }

  my $titles = [$task_name];
  if ( defined $config->{$section}{parameterSampleFile2} ) {
    my $correlationGroups = get_raw_files( $config, $section, "parameterSampleFile2" );
    for my $correlationTitle ( keys %$correlationGroups ) {
      my $groups = $correlationGroups->{$correlationTitle};
      if ( is_hash($groups) ) {
        if ($correlationTitle ne "all"){
          push( @$titles, $correlationTitle );
        }
      }
    }
  }

  my @result_files = ();
  for my $title (@$titles) {
    foreach my $sample_name ( keys %temp ) {
      my $curSuffix = $suffix;
      if ( $title ne $task_name ) {
        $curSuffix = $suffix . "." . $title;
      }
      $curSuffix =~ s/\\s+/_/g;

      foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
        my $prefix = $subSampleFile;
        if ($output_to_result_dir) {
          my $file = basename($subSampleFile);
          if($output_include_folder_name){
            my $pdir = dirname($subSampleFile);
            while ( basename($pdir) eq "result" ) {
              $pdir = dirname($pdir);
            }
            $prefix = $result_dir . "/" . basename($pdir) . "." . $file;
          }else{
            $prefix = $result_dir . "/" . $file;
          }
        }

        $prefix = $prefix . $curSuffix;

        foreach my $output_file_ext_one (@$output_file_exts) {
          push( @result_files, "${prefix}${output_file_ext_one}" );
        }

        if ( ( $output_file_task_ext ne "" ) && ( $sample_name eq $task_name ) ) {
          if ( $output_file_task_ext =~ /;/ ) {
            my @output_file_task_exts = split( ";", $output_file_task_ext );
            foreach my $output_file_ext_one (@output_file_task_exts) {
              push( @result_files, "${prefix}${output_file_ext_one}" );
            }
          }
          else {
            push( @result_files, "${prefix}${output_file_task_ext}" );
          }
        }

      }
    }
  }
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  #print(Dumper($result));
  return $result;
}

1;
