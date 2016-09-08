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

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $blacklistfile = get_param_file( $config->{$section}{"blacklist_file"}, "blacklist_file", 0 );
  my $isPairedEnd    = get_option( $config, $section, "is_paired_end" );
  my $isSortedByName = 0;
  my $maxFragmentLength = 0;
  my $minFragmentLength = 0;
  if($isPairedEnd){
    $option = $option . " -bedpe";
    $isSortedByName = get_option( $config, $section, "is_sorted_by_name" );
    $maxFragmentLength = get_option( $config, $section, "maximum_fragment_length" );
    $minFragmentLength = get_option( $config, $section, "minimum_fragment_length" );
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

    my $bed_file       = $sample_name . ".bed";
    my $confident_file = $sample_name . ".confident.bed";
    my $final_file     = defined $blacklistfile ? $sample_name . ".confident.shifted.bed" : $sample_name . ".shifted.bed";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    my $rmlist = "";

    my $sortCmd = "";
    if ( !$isSortedByName ) {
      my $presortedFile   = $sample_name . ".sortedByName.bam";
      print $pbs "
if [ ! -s $presortedFile ]; then
  echo SortBamByName=`date` 
  samtools sort -n  -@ $thread -m $memory -o $presortedFile $bam_file 
fi
";
      $rmlist   = $presortedFile;
      $bam_file = $presortedFile;
    }

    print $pbs "
if [ ! -s $bed_file ]; then
  echo bamtobed=`date` 
  bedtools bamtobed $option -i $bam_file > $bed_file
fi
";
    $rmlist = $rmlist . " " . $bed_file;
    
    if($isPairedEnd){
      my $slim_file = $sample_name . ".paired_end.bed";
      print $pbs "
if [[ -s $bed_file && ! -s $slim_file ]]; then
  echo convert_paired_end_bed=`date`
  awk 'BEGIN {OFS = \"\\t\"} ; {if (\$1 == \$4 && \$6 - \$2 <= $maxFragmentLength && \$6 - \$2 >= $minFragmentLength ) print \$1, \$2, \$6, \$7, \$8, \$9}' $bed_file > $slim_file 
fi
";
      $bed_file = $slim_file;
      $rmlist = $rmlist . " " . $slim_file;
    }

    if ( defined $blacklistfile ) {
      print $pbs "
if [[ -s $bed_file && ! -s $confident_file ]]; then
  echo remove_read_in_blacklist=`date` 
  bedtools intersect -v -a $bed_file -b $blacklistfile > $confident_file
fi
";
      $bed_file = $confident_file;
      $rmlist   = $rmlist . " " . $confident_file;
    }

    print $pbs "
if [[ -s $bed_file && ! -s $final_file ]]; then
  echo shift_reads=`date` 
  awk 'BEGIN {OFS = \"\\t\"} ; {if (\$6 == \"+\") print \$1, \$2 + 4, \$3 + 4, \$4, \$5, \$6; else print \$1, \$2 - 5, \$3 - 5, \$4, \$5, \$6}' $bed_file > $final_file
fi

if [ -s $final_file ]; then
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
  my $blacklistfile = get_param_file( $config->{$section}{"blacklist_file"}, "blacklist_file", 0 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();

    my $final_file = defined $blacklistfile ? $sample_name . ".confident.shifted.bed" : $sample_name . ".shifted.bed";
    push( @result_files, "${result_dir}/${final_file}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
