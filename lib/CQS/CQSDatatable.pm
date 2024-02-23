#!/usr/bin/perl
package CQS::CQSDatatable;

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
  $self->{_suffix} = "_tb";
  bless $self, $class;
  return $self;
}

sub get_result {
  my ( $self, $result_dir, $task_name, $option ) = @_;

  my $result;
  if ( $option =~ /-o\s+(\S+)/ ) {
    $result = $1;
  }
  else {
    $result = $self->get_file( $result_dir, $task_name, ".count", 0 );
    $option = $option . " -o " . $result;
  }
  return ( $result, $option );
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $mapFile = get_param_file( $config->{$section}{name_map_file}, "name_map_file", 0 );

  my $cpm_r = dirname(__FILE__) . "/../Count/count2cpm.r";
  if ( !-e $cpm_r ) {
    die "File not found : " . $cpm_r;
  }

  my $mapoption = "";
  if ( defined $mapFile ) {
    $mapoption = "-m $mapFile";
  }
  my $output_proteincoding_gene=get_option( $config, $section, "output_proteincoding_gene", 0 );
  if ($output_proteincoding_gene eq '' && $mapFile ne "") {
    $output_proteincoding_gene=1;
  }

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $filelist = $self->get_file( $pbs_dir, $task_name, ".filelist" );
  open( FL, ">$filelist" ) or die "Cannot create $filelist";
  for my $sample_name ( sort keys %raw_files ) {
    my @bam_files = @{ $raw_files{$sample_name} };
    my $bam_file  = $bam_files[0];
    print FL $sample_name, "\t", $bam_file, "\n";
  }
  close(FL);

  my $result_file = $self->get_file( ".", $task_name, ".count", 0 );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $final_file = $result_file;  
  if ($output_proteincoding_gene) {
    $final_file = $self->get_file( ".", $task_name, ".proteincoding.count", 0 );
  }

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  print $pbs "
cqstools data_table $option -o $result_file -l $filelist $mapoption 

R --vanilla -f $cpm_r --args $result_file $result_file
";
  if ($output_proteincoding_gene) {
  print $pbs "
R --vanilla -f $cpm_r --args $final_file $final_file
";
  }

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $mapFile = $config->{$section}{name_map_file};
  my $output_proteincoding_gene=get_option( $config, $section, "output_proteincoding_gene", 0 );
  if ($output_proteincoding_gene eq '' && $mapFile ne "") {
    $output_proteincoding_gene=1;
  }

  my $result       = {};
  my @result_files = ();
  push( @result_files, $self->get_file( $result_dir, $task_name, ".count", 0 ) );
  push( @result_files, $self->get_file( $pbs_dir, $task_name, ".filelist" ) );

  if ( defined $config->{$section}{name_map_file} ) {
    push( @result_files, $self->get_file( $result_dir, $task_name, ".fpkm.tsv", 0 ) );
    if ($output_proteincoding_gene) {
      push( @result_files, $self->get_file( $result_dir, $task_name, ".fpkm.proteincoding.tsv", 0 ) );
    }
  }
  push( @result_files, $self->get_file( $result_dir, $task_name, ".cpm.csv", 0 ) );
  if ($output_proteincoding_gene) {
    push( @result_files, $self->get_file( $result_dir, $task_name, ".proteincoding.count", 0 ) );
    push( @result_files, $self->get_file( $result_dir, $task_name, ".proteincoding.cpm.csv", 0 ) );
  }
  
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
