#!/usr/bin/perl
package ATACseq::BamToBed;

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

sub get_final_file {
  my ( $sample_name, $blacklistfile, $shiftPosition ) = @_;
  my $result = $sample_name;
  if ( defined $blacklistfile ) {
    $result = $sample_name . ".confident";
  }

  if ($shiftPosition) {
    $result = $sample_name . ".shifted";
  }
  $result = $result . ".bed";
  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );
  
  my $sort_memory = $thread == 1? $memory:"4G";

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $blacklistfile = get_param_file( $config->{$section}{"blacklist_file"}, "blacklist_file", 0 );
  my $isPairedEnd       = get_option( $config, $section, "is_paired_end" );
  my $isSortedByName    = 0;
  my $maxFragmentLength = 0;
  my $minFragmentLength = 0;
  if ($isPairedEnd) {
    $option            = $option . " -bedpe";
    $isSortedByName    = get_option( $config, $section, "is_sorted_by_name" );
    $maxFragmentLength = get_option( $config, $section, "maximum_fragment_length" );
    $minFragmentLength = get_option( $config, $section, "minimum_fragment_length" );
  }
  my $shiftPosition = get_option( $config, $section, "shift_position", 0 );

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

    my $final_file = $sample_name . ".bed";
    
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    my $rmlist = "";

    my $sortCmd = "";
    if ( !$isSortedByName ) {
      my $presortedFile = $sample_name . ".sortedByName.bam";
      print $pbs "
if [ ! -s $presortedFile ]; then
  echo SortBamByName=`date` 
  samtools sort -n  -@ $thread -m $sort_memory -o $presortedFile $bam_file 
fi
";
      $rmlist   = $presortedFile;
      $bam_file = $presortedFile;
    }

    my $bed_file = $isPairedEnd ? $sample_name . ".raw.bedpe" : $sample_name . ".raw.bed";
    print $pbs "
if [ ! -s $bed_file ]; then
  echo bamtobed=`date` 
  bedtools bamtobed $option -i $bam_file > $bed_file
fi
";
    if ($isPairedEnd) {
      my $slim_file = $sample_name . ".slim.bed";
      print $pbs "
if [[ -s $bed_file && ! -s $slim_file ]]; then
  echo convert_paired_end_bed=`date`
  awk 'BEGIN {OFS = \"\\t\"} ; {if (\$1 == \$4 && \$6 - \$2 <= $maxFragmentLength && \$6 - \$2 >= $minFragmentLength ) print \$1, \$2, \$6, \$7, \$8, \$9}' $bed_file > $slim_file 
fi
";
      $rmlist   = $rmlist . " " . $bed_file;
      $bed_file = $slim_file;
    }

    if ( defined $blacklistfile ) {
      my $confident_file = $sample_name . ".confident.bed";
      print $pbs "
if [[ -s $bed_file && ! -s $confident_file ]]; then
  echo remove_read_in_blacklist=`date` 
  bedtools intersect -v -a $bed_file -b $blacklistfile > $confident_file
fi
";
      $rmlist   = $rmlist . " " . $bed_file;
      $bed_file = $confident_file;
    }

    if ($shiftPosition) {
      my $shiftPositive = $shiftPosition;
      my $shiftNegative = $shiftPosition + 1;
      my $shiftFile = $sample_name + ".shifted.bed";
      print $pbs "
if [[ -s $bed_file && ! -s $shiftFile ]]; then
  echo shift_reads=`date` 
  awk 'BEGIN {OFS = \"\\t\"} ; {if (\$6 == \"+\") print \$1, \$2 + $shiftPositive, \$3 + $shiftPositive, \$4, \$5, \$6; else print \$1, \$2 - $shiftNegative, \$3 - $shiftNegative, \$4, \$5, \$6}' $bed_file > $shiftFile
fi
";
      $rmlist = $rmlist . " " . $bed_file . " " . $shiftFile;
      $bed_file = $shiftFile;
    }

    print $pbs "
if [ -s $bed_file ]; then
  sort -k1,1V -k2,2n -k3,3n $bed_file > $final_file
  rm $rmlist;
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

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
