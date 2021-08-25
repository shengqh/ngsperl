#!/usr/bin/perl
package Bedtools::Bamtobed;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_b2b";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $blacklistfile = get_param_file( $config->{$section}{"blacklist_file"}, "blacklist_file", 0 );
  my $isPairedEnd = get_is_paired_end_option( $config, $section);
  my $is_sorted_by_name = get_option( $config, $section, "is_sorted_by_name" );
  my $sort_memory = $thread == 1 ? $memory : "4G";

  if ($isPairedEnd) {
    $option = $option . " -bedpe";
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_file = get_final_file($config, $section, $result_dir, $sample_name);
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    my $rmlist = "";
    if ( $isPairedEnd and not $is_sorted_by_name ) {
      my $sortedFile = $sample_name . ".sortedByName.bam";
      print $pbs "
if [ ! -s $sortedFile ]; then
  echo sort=`date`
  samtools sort -n -m $sort_memory -o $sortedFile $bam_file
fi
";
      $rmlist   = $rmlist . " " . $sortedFile;
      $bam_file = $sortedFile;
    }

    my $bed_file = $sample_name . ".converted.bed";
    $rmlist = $rmlist . " $bed_file";
    print $pbs "
if [[ -s $bam_file && ! -s $bed_file ]]; then
  echo bamtobed=`date` 
  bedtools bamtobed $option -i $bam_file > $bed_file 
fi
";
    if ( defined $blacklistfile ) {
      my $confident_file = $sample_name . ".confident.bed";
      print $pbs "
if [[ -s $bed_file && ! -s $confident_file ]]; then
  echo remove_read_in_blacklist=`date` 
  bedtools intersect -v -a $bed_file -b $blacklistfile > $confident_file
fi
";
      $rmlist   = $rmlist . " " . $confident_file;
      $bed_file = $confident_file;
    }

    print $pbs "
if [ -s $bed_file ]; then
  echo sort_bed=`date` 
  sort -k1,1V -k2,2n -k3,3n $bed_file > $final_file
  #rm $rmlist;
fi
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();

    my $final_file = $sample_name . ".bed";
    push( @result_files, "${result_dir}/${final_file}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
