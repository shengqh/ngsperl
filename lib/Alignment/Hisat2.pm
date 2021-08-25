#!/usr/bin/perl
package Alignment::Hisat2;

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
use Alignment::AlignmentUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_h2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  if ( $option !~ /outSAMprimaryFlag/ ) {
    $option = $option . "  --outSAMprimaryFlag AllBestScore";
  }

  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );

  my $output_to_same_folder     = get_option( $config, $section, "output_to_same_folder",     1 );
  my $genome_dir = parse_param_file( $config, $section, "genome_dir", 1 );

  my %fqFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %fqFiles ) {
    my @sample_files = @{ $fqFiles{$sample_name} };
    my $samples;
    if(scalar(@sample_files) == 2){
      $samples = "-1 " . $sample_files[0] . " -2 " . $sample_files[1];
    }else{
      $samples = $sample_files[0];
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );
    my $finalSam = $sample_name . ".sam";
    my $final = $sample_name . ".bam";
    my $finalStat = $final . ".stat";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final );

    my $chromosome_grep_command = getChromosomeFilterCommand( $final, $chromosome_grep_pattern );

    print $pbs "
hisat2 -p $thread --dta -x $genome_dir $samples --novel-splicesite-outfile ${sample_name}.splicing  -S $finalSam 2> $finalStat

if [ -s $finalSam ]; then
  echo Samtools::Sort_start=`date`
  samtools sort -o ${final}.tmp $finalSam
  mv ${final}.tmp $final
  rm $finalSam
fi

if [ -s $final ]; then
  echo Samtools::Index_start=`date`
  samtools index $final
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder", 1 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $cur_dir = $output_to_same_folder ? $result_dir : $result_dir . "/$sample_name";

    my @result_files              = ();
    push( @result_files, "$cur_dir/${sample_name}.bam" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
