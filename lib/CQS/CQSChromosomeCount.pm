#!/usr/bin/perl
package CQS::CQSChromosomeCount;

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
  $self->{_suffix} = "_cc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my %seqcount_files = ();
  if ( has_raw_files( $config, $section, "seqcount" ) ) {
    %seqcount_files = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my %pmNameFiles = ();
  if ( has_raw_files( $config, $section, "perfect_mapped_name" ) ) {
    %pmNameFiles = %{ get_raw_files( $config, $section, "perfect_mapped_name" ) };
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @bam_files  = @{ $raw_files{$sample_name} };
    my $bam_file   = $bam_files[0];
    my $fileName   = basename($bam_file);
    my $final_file = $fileName . ".count";

    my $seqcountFile = "";
    if ( defined $seqcount_files{$sample_name} ) {
      my @seqcounts = @{ $seqcount_files{$sample_name} };
      my $seqcount  = $seqcounts[0];
      $seqcountFile = " -c $seqcount";
    }

    my $pmNameFile = "";
    if ( defined $pmNameFiles{$sample_name} ) {
      my @pmNames = @{ $pmNameFiles{$sample_name} };
      my $pmname  = $pmNames[0];
      $pmNameFile = " -n $pmname";
    }

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "mono $cqstools chromosome_count $option $pmNameFile -i $bam_file -o $final_file $seqcountFile";

    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all trna_count tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $cur_dir = $result_dir . "/$sample_name";

    my @bam_files = @{ $raw_files{$sample_name} };
    my $bam_file  = $bam_files[0];
    my $fileName  = basename($bam_file);

    my @result_files = ();
    my $countFile    = "${cur_dir}/${fileName}.count";
    push( @result_files, $countFile );
    push( @result_files, "${countFile}.mapped.xml" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
