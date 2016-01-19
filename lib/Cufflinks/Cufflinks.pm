#!/usr/bin/perl
package Cufflinks::Cufflinks;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_clinks";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );
  my $gtf = "";
  if ( defined $transcript_gtf ) {
    $gtf = "-g $transcript_gtf";
  }

  my %tophat2map = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";

  for my $sample_name ( sort keys %tophat2map ) {
    my @tophat2Files = @{ $tophat2map{$sample_name} };
    my $tophat2File  = $tophat2Files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $cur_dir  = create_directory_or_die( $result_dir . "/$sample_name" );

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = "transcripts.gtf";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    print $pbs "cufflinks $option $gtf -o . $tophat2File";
    $self->close_pbs( $pbs, $pbs_file );
    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;
  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %tophat2map = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %tophat2map ) {
    my $cur_dir      = $result_dir . "/$sample_name";
    my @result_files = ();
    push( @result_files, $cur_dir . "/transcripts.gtf" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
