#!/usr/bin/perl
package UMI::UmiReadsProcessing;

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
  $self->{_suffix} = "_ur";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );


  my $bwa_index = $config->{$section}{bwa_index};
  if ( !defined $bwa_index ) {
    $bwa_index = $config->{$section}{fasta_file} or die "define ${section}::bwa_index first";
  }
  my $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1 );
  my %unmapped_bam_files = %{ get_raw_files( $config, $section ) };
  
  my $extension = get_option( $config, $section, "extension", "_umi" );
  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option || $java_option eq "" ) {
    $java_option = "-Xmx${memory}";
  }
  
  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  
  for my $sample_name ( sort keys %unmapped_bam_files ) {
    my @sample_files = @{ $unmapped_bam_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

    print $pbs "echo Processing $sample_name \n"; 

    print $pbs "
if [[ ! -s $bam_file.fastq ]]; then
  java $java_option -jar $picard_jar SamToFastq I=$bam_file F=$bam_file.fastq INTERLEAVE=true
fi
";
    print $pbs "
if [[ ! -s $bam_file$extension.sam ]]; then
  bwa mem -p -t 8 $bwa_index $bam_file.fastq > $bam_file$extension.sam
fi
";

    print $pbs "
if [[ ! -s $bam_file.fastq ]]; then
  java $java_option -jar $picard_jar MergeBamAlignment UNMAPPED=1356-YZ-71.bam ALIGNED=1356-YZ-71_umi.sam O=1356-YZ-71_umim.bam R=/scratch/cqs/zhaos/reference/hg38/hg38bundle/bwa_index_0.7.12/Homo_sapiens_assembly38.fasta SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true 
fi
";

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  
  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  my $result = {};

  my $extension = get_option( $config, $section, "extension", "_umi" );

  my %unmapped_bam_files = %{ get_raw_files( $config, $section ) };
  for my $sample_name ( sort keys %unmapped_bam_files ) {
    my $bamOut       = $sample_name .$extension."_consensus.bam";
    my @result_files = ();
    push( @result_files, "${result_dir}/${bamOut}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
