#!/usr/bin/perl
package Picard::MarkDuplicates;

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
  $self->{_suffix} = "_rd";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];

    my $rmdupFile    = $sample_name . ".rmdup.bam";
    my $sortedPrefix = $sample_name . ".rmdup_sorted";
    my $sortedFile   = $sortedPrefix . ".bam";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $sortedFile );

    print $pbs "

if [ ! -s $rmdupFile ]; then
  echo RemoveDuplicate=`date` 
  java $option -jar $picard_jar MarkDuplicates I=$sampleFile O=$rmdupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT M=${rmdupFile}.metrics
fi

if [[ -s $rmdupFile && ! -s $sortedFile ]]; then
  echo BamSort=`date` 
  samtools sort $rmdupFile $sortedPrefix 
fi

if [[ -s $sortedFile && ! -s ${sortedFile}.bai ]]; then
  echo BamIndex=`date` 
  samtools index $sortedFile
  samtools flagstat $sortedFile > ${sortedFile}.stat
  rm $rmdupFile
fi
  
";

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Picard::MarkDuplicates tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $sortedFile = $sample_name . ".rmdup_sorted.bam";

    my @result_files = ();
    push( @result_files, "${result_dir}/${sample_name}/${sortedFile}" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
