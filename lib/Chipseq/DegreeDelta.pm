#!/usr/bin/perl
package Chipseq::DegreeDelta;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use Data::Dumper;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_dd";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $edge_files = get_raw_files( $config, $section);
  my $bam_files = get_raw_files($config, $section, "bam_files");

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  print $pbs "out_degree_delta \\\n";

  my @samples = sort keys %$edge_files;
  #print(Dumper(@samples));
  #print(Dumper($edge_files));
  #print(Dumper($bam_files));

  my $sample_count = scalar(@samples) - 1;
  for my $idx (0..$sample_count){
    my $sample = $samples[$idx];
    #print($idx . " : " . $sample . "\n");
    my $edge_files = $edge_files->{$sample};
    my $bam_files = $bam_files->{$sample};
    print $pbs "  --edge_table_" . ($idx+1) . " " . $edge_files->[0] . " \\\n";
    print $pbs "  --bams_" . ($idx+1) . " ";
    for my $bam_file (@$bam_files) {
      print $pbs "    " . $bam_file . " \\\n";
    }
  }
  print $pbs "  -o $result_dir \\
  -n $task_name
";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result = {
    $task_name => filter_array( ["$result_dir/placehoder.txt"], $pattern )
  };

  return $result;
}

1;
