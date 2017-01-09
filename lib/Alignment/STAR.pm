#!/usr/bin/perl
package Alignment::STAR;

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
  $self->{_suffix} = "_star";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  if ( $option !~ /outSAMprimaryFlag/ ) {
    $option = $option . "  --outSAMprimaryFlag AllBestScore";
  }

  my $output_sort_by_coordinate = get_option( $config, $section, "output_sort_by_coordinate", 0 );
  my $output_unsorted           = get_option( $config, $section, "output_unsorted",           0 );
  if ( !$output_sort_by_coordinate && !$output_unsorted ) {
    $output_unsorted = 1;
  }
  my $output_format = "--outSAMtype BAM";
  if ($output_sort_by_coordinate) {
    $output_format = $output_format . " SortedByCoordinate";
  }
  if ($output_unsorted) {
    $output_format = $output_format . " Unsorted";
  }
  my $genome_dir = parse_param_file( $config, $section, "genome_dir", 1 );

  my %fqFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %fqFiles ) {
    my @sample_files = @{ $fqFiles{$sample_name} };

    my $uncompress = ( $sample_files[0] =~ /.gz$/ ) ? " --readFilesCommand zcat" : "";

    my $samples = join( " ", @sample_files );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );
    my $rgline  = "ID:$sample_name SM:$sample_name LB:$sample_name PL:ILLUMINA PU:ILLUMINA";

    my $final = $output_sort_by_coordinate ? $sample_name . "_Aligned.sortedByCoord.out.bam" : $sample_name . "_Aligned.out.bam";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final );

    print $pbs "
STAR $option --outSAMattrRGline $rgline --runThreadN $thread --genomeDir $genome_dir --readFilesIn $samples $uncompress --outFileNamePrefix ${sample_name}_ $output_format  

if [ -s $final ]; then
  samtools index $final
  samtools flagstat $final > ${final}.stat
  rm -rf ${sample_name}__STARgenome ${sample_name}__STARpass1 ${sample_name}_SJ.out.tab ${sample_name}_Log.progress.out
fi

";

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $sort_by_coordinate = get_option_value( $config->{$section}{sort_by_coordinate}, 0 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files              = ();
    my $output_sort_by_coordinate = get_option( $config, $section, "output_sort_by_coordinate", 0 );
    my $output_unsorted           = get_option( $config, $section, "output_unsorted", 0 );
    if ( !$output_sort_by_coordinate && !$output_unsorted ) {
      $output_unsorted = 1;
    }
    if ($output_sort_by_coordinate) {
      push( @result_files, "${result_dir}/${sample_name}/${sample_name}_Aligned.sortedByCoord.out.bam" );

    }
    if ($output_unsorted) {
      push( @result_files, "${result_dir}/${sample_name}/${sample_name}_Aligned.out.bam" );
    }
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}_Log.final.out" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
