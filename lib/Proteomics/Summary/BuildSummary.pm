#!/usr/bin/perl
package Proteomics::Summary::BuildSummary;

use strict;
use warnings;
use File::Basename;
use File::Slurp;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::StringUtils;
use Proteomics::Summary::AbstractBuildSummary;

our @ISA = qw(Proteomics::Summary::AbstractBuildSummary);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_bs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $proteomicstools = get_param_file( $config->{$section}{proteomicstools}, "proteomicstools", 1, not $self->using_docker() );

  my ( $datasets, $lines, $dataset ) = $self->get_datasets( $config, $section );
  my %datasets = %{$datasets};
  my @lines    = @{$lines};
  my @dataset  = @{$dataset};

  my $currentParamFile = $result_dir . "/" . $task_name . ".param";
  open( my $out, ">$currentParamFile" ) or die $!;

  for ( my $index = 0 ; $index < scalar(@lines) ; $index++ ) {
    print $out $lines[$index] . "\n";
    if ( $lines[$index] =~ "<Datasets>" ) {
      last;
    }
  }

  #print @dataset;

  for my $sample_name ( sort keys %datasets ) {
    print $out "    <Dataset>\n";
    print $out "      <Name>$sample_name</Name>\n";
    foreach my $dsline (@dataset) {
      print $out $dsline . "\n";
    }
    print $out "      <PathNames>\n";
    my @sample_files = @{ $datasets{$sample_name} };
    for my $sampleFile (@sample_files) {
      print $out "        <PathName>$sampleFile</PathName>\n";
    }
    print $out "      </PathNames>\n";
    print $out "    </Dataset>\n";
  }

  for ( my $index = 0 ; $index < scalar(@lines) ; $index++ ) {
    if ( $lines[$index] =~ "</Datasets>" ) {
      for ( my $nextindex = $index ; $nextindex < scalar(@lines) ; $nextindex++ ) {
        print $out $lines[$nextindex] . "\n";
      }
      last;
    }
  }

  close $out;

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $result_file = $result_dir . "/" . $task_name . ".noredundant";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $result_file );
  print $pbs "mono $proteomicstools buildsummary -i $currentParamFile";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result       = {};
  my $result_file  = $result_dir . "/" . $task_name . ".noredundant";
  my @result_files = ();
  push( @result_files, $result_file );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
