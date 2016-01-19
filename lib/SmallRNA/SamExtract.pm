#!/usr/bin/perl
package SmallRNA::SamExtract;

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
  $self->{_suffix} = "_se";
  $self->{_name} = __PACKAGE__;
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );
  my %count_files = %{ get_raw_files( $config, $section ) };
  my %bam_files = %{ get_raw_files( $config, $section, "bam_files" ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %count_files ) {
    my @sample_counts = @{ $count_files{$sample_name} };
    my $sample_count  = $sample_counts[0];

    my @sample_bams = @{ $bam_files{$sample_name} };
    my $sample_bam  = $sample_bams[0];

    my $final_sam  = $sample_name . ".mapped.sam";
    my $final_file = $sample_name . ".mapped.bam";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs($pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $pbs "
mono $cqstools sam_extract $option --bam $sample_bam --count $sample_count -o $final_sam

if [ -s $final_sam ]; then
  samtools view -S -b $final_sam > $final_file
  if [ -s $final_file ]; then
    samtools index $final_file
    rm $final_sam
  fi
fi
";

    $self->close_pbs( $pbs );

    print "$pbs_file created \n";
  }
  close($sh);
  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );
  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my @result_files = ();
    push( @result_files, $result_dir . "/" . $sample_name . ".mapped.bam" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
