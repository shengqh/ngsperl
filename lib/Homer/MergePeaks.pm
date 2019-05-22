#!/usr/bin/perl
package Homer::MergePeaks;

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $overlapRate = get_option($config, $section, "overlap_rate", 0.5);
  my %raw_files = %{ get_grouped_raw_files( $config, $section, "groups" ) };

  my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name   = basename($pbs_file);
  my $log        = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc   = $cluster->get_log_description($log);

  my @groupNames = sort keys %raw_files;
  my $final_file = $groupNames[-1]. ".filtered.sorted.bed";
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  for my $sample_name ( @groupNames ) {
    my $merged_file = "${sample_name}.merged.bed";
    my $filtered_file = "${sample_name}.filtered.bed";
    my $final_file = "${sample_name}.filtered.sorted.bed";

    print $pbs "mergePeaks \\";
    my @peak_files = @{ $raw_files{$sample_name} };
    my $fileCount = scalar(@peak_files);
    for my $peak_file (sort @peak_files) {
      print $pbs "  \"$peak_file\" \\\n";
    }
    print $pbs " > $merged_file \n\n";
    print $pbs "awk '{n=split(\$7, array, \"|\"); if(n / $fileCount > $overlapRate) print \$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$1\"\\t\"n\"\\t\"\$5}' $merged_file > $filtered_file \n\n";
    print $pbs "bedtools sort -i $filtered_file > $final_file \n\n";
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result = {};

  my %raw_files = %{ get_grouped_raw_files( $config, $section, "groups" ) };
  for my $sample_name ( sort keys %raw_files ) {
    $result->{$sample_name} = filter_array( ["$result_dir/${sample_name}.filtered.sorted.bed"], $pattern );
  }    
  return $result;
}

1;
