#!/usr/bin/perl
package RNAediting::ParseMutation;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $python_script = dirname(__FILE__) . "/parseMutation.py";
  if ( !-e $python_script ) {
    die "File not found : " . $python_script;
  }
  
  #my $rTemplate = dirname(__FILE__) . "/parseMutation.r";
  #if ( !-e $rTemplate ) {
  #  die "File not found : " . $rTemplate;
  #}
  
  #my $picture_width = get_option($config, $section, "picture_width", 3000);
  #my $picture_height = get_option($config, $section, "picture_height", 3000);
  #my $percentageThresholds = get_option($config, $section, "percentage_thresholds", "1 2 5");
  
  my $positions = "\"" . get_option($config, $section, "positions") . "\"";
  my $sequence_fasta = get_param_file($config->{$section}{"sequence_fasta"}, "sequence_fasta", 1);
  
  my $raw_files = get_raw_files( $config, $section );
  my @bamfiles = ();
  for my $key (sort keys %$raw_files){
    push(@bamfiles, $raw_files->{$key}[0])
  }
  my $bamfile = join(',', @bamfiles);

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $final_file = $task_name . ".tsv";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  print $pbs "
python3 $python_script $option -i $bamfile -o $final_file -f $sequence_fasta -p $positions

";
#R --vanilla -f $rTemplate --args ${task_name}.SNV.tsv ${task_name}.SNV $picture_width $picture_height $percentageThresholds
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  
  push( @result_files, $result_dir . "/${task_name}.tsv" );
  push( @result_files, $result_dir . "/${task_name}.summary.tsv" );
  push( @result_files, $result_dir . "/${task_name}.reads.tsv" );
  #push( @result_files, $result_dir . "/${task_name}.SNV.tsv" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
